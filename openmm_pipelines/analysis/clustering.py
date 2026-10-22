"""Conformational clustering for the default-overview analysis, anchored to the design.

The question this answers for a *designed* structure is conformational heterogeneity: does
the sampled ensemble stay as one tight cluster near the design pose, or split into several
conformations? So alongside cluster populations we report each representative's RMSD to the
design and the fraction of the ensemble in the design-like basin — the stability /
variant-ranking signal.

Two methods (`AnalysisConfig.cluster_method`):
- **hierarchical** (default): scipy average-linkage on the OPTIMAL pairwise Cα-RMSD matrix
  (`md.rmsd` re-aligns each pair). An interpretable Å cutoff ("within X Å = same
  conformation") that is comparable across any system. O(N^2); above `cluster_max_frames`
  it uniformly strides down, recording the applied stride (never silent).
- **dbscan** (optional, needs scikit-learn): density clustering on the flattened,
  already-design-superposed Cα coordinates, with `eps = cutoff * sqrt(n_atoms)`. NOTE it
  therefore clusters on COMMON-FRAME RMSD (distance to the shared reference frame), an upper
  bound on the optimal pairwise RMSD the hierarchical path uses — the same cutoff gives a
  slightly different partition. Valid in the small-drift regime; use for very large ensembles.

No hand-rolled clustering — scipy / scikit-learn do the work.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import mdtraj as md
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from . import plots
from .config import AnalysisConfig, ClusterMethod
from .metrics import rmsd_to_reference

_NM_TO_ANG = 10.0


@dataclass(frozen=True)
class ClusterResult:
    """Outcome of clustering an ensemble, ranked by population and anchored to the design.

    `labels` / `frame_indices` are per CLUSTERED frame (after any striding), aligned to each
    other; `labels` are 0-based ranks by descending population (DBSCAN noise kept as -1).
    `frame_indices` index the ORIGINAL trajectory. The per-cluster lists are in rank order
    and EXCLUDE noise. `design_like_cluster` is the rank whose representative is closest to
    the design pose; `design_like_fraction` is that cluster's share of the clustered frames.
    """

    labels: np.ndarray
    frame_indices: np.ndarray
    populations: list[int]
    representative_frames: list[int]
    representative_rmsd_to_design_ang: list[float]
    design_like_cluster: int
    design_like_fraction: float
    n_clusters: int
    n_frames_clustered: int
    cutoff_ang: float
    method: str
    stride: int
    min_samples: int | None


def pairwise_rmsd(traj: md.Trajectory, selection: str = "name CA") -> np.ndarray:
    """Symmetric (n_frames, n_frames) optimal pairwise Cα RMSD matrix in Å."""
    idx = traj.topology.select(selection)
    assert len(idx) > 0, f"selection {selection!r} matched no atoms"
    n = traj.n_frames
    dist = np.empty((n, n))
    for i in range(n):
        dist[i] = md.rmsd(traj, traj, i, atom_indices=idx)
    dist *= _NM_TO_ANG
    dist = 0.5 * (dist + dist.T)  # enforce symmetry against FP asymmetry
    np.fill_diagonal(dist, 0.0)
    return dist


def rank_by_population(labels: np.ndarray) -> tuple[np.ndarray, list[int]]:
    """Relabel clusters 0,1,… by DESCENDING population (noise -1 kept as -1).

    Returns ``(ranked_labels, populations)`` where ``populations[r]`` is the size of rank r.
    """
    present = [lab for lab in set(labels.tolist()) if lab != -1]
    order = sorted(present, key=lambda lab: (-(labels == lab).sum(), lab))  # size desc, det.
    remap = {old: new for new, old in enumerate(order)}
    ranked = np.array([remap[lab] if lab != -1 else -1 for lab in labels], dtype=int)
    populations = [int((labels == old).sum()) for old in order]
    return ranked, populations


def _subsample(traj: md.Trajectory, max_frames: int) -> tuple[md.Trajectory, np.ndarray, int]:
    """Uniformly stride `traj` down to at most `max_frames`. Returns (sub, original_indices, stride)."""
    n = traj.n_frames
    if n <= max_frames:
        return traj, np.arange(n), 1
    stride = int(np.ceil(n / max_frames))
    idx = np.arange(0, n, stride)
    return traj[idx], idx, stride


def _medoid(dist: np.ndarray, member_positions: np.ndarray) -> int:
    """Position (into the clustered frames) of the medoid: min summed in-cluster RMSD."""
    sub = dist[np.ix_(member_positions, member_positions)]
    return int(member_positions[int(sub.sum(axis=1).argmin())])


def _centroid_representative(coords: np.ndarray, member_positions: np.ndarray) -> int:
    """Position of the frame nearest the cluster's mean (flattened) coordinate."""
    centroid = coords[member_positions].mean(axis=0)
    d = ((coords[member_positions] - centroid) ** 2).sum(axis=1)
    return int(member_positions[int(d.argmin())])


def _dbscan_labels(sub: md.Trajectory, selection: str, cutoff_ang: float, min_samples: int):
    """(labels, flattened coords) from DBSCAN on design-superposed Cα coordinates."""
    try:
        from sklearn.cluster import DBSCAN
    except ImportError as exc:  # optional dependency — fail loud with the fix
        raise ImportError(
            "cluster_method='dbscan' needs scikit-learn (optional). "
            "Install it: pip install 'openmm_pipelines[cluster]'."
        ) from exc
    idx = sub.topology.select(selection)
    assert len(idx) > 0, f"selection {selection!r} matched no atoms"
    coords = sub.xyz[:, idx, :].reshape(sub.n_frames, -1)  # nm, flattened; frames pre-superposed
    # For frames superposed to one common frame, Euclidean distance = sqrt(n_atoms) * RMSD,
    # so an Å RMSD cutoff maps to eps in nm as below (coords are nm).
    eps = (cutoff_ang / _NM_TO_ANG) * np.sqrt(len(idx))
    labels = DBSCAN(eps=eps, min_samples=min_samples, metric="euclidean").fit_predict(coords)
    return labels, coords


def cluster_conformations(
    traj: md.Trajectory, reference: md.Trajectory, config: AnalysisConfig
) -> ClusterResult:
    """Cluster `traj`'s conformations and anchor the result to the `reference` design pose."""
    assert traj.n_frames >= 1, "need at least 1 frame to cluster"
    rmsd_to_design = rmsd_to_reference(traj, reference, config.cluster_selection)  # Å, per frame

    sub, frame_indices, stride = _subsample(traj, config.cluster_max_frames)
    m = sub.n_frames

    if m == 1:  # degenerate: a single conformation
        ranked = np.array([0], dtype=int)
        populations = [1]
        rep_positions = [0]
        min_samples = None
    elif config.cluster_method is ClusterMethod.HIERARCHICAL:
        dist = pairwise_rmsd(sub, config.cluster_selection)
        raw = fcluster(
            linkage(squareform(dist, checks=False), method="average"),
            t=config.cluster_cutoff_ang, criterion="distance",
        )
        ranked, populations = rank_by_population(raw)
        rep_positions = [_medoid(dist, np.where(ranked == r)[0]) for r in range(len(populations))]
        min_samples = None
    elif config.cluster_method is ClusterMethod.DBSCAN:
        raw, coords = _dbscan_labels(
            sub, config.cluster_selection, config.cluster_cutoff_ang, config.cluster_min_samples
        )
        ranked, populations = rank_by_population(raw)
        rep_positions = [
            _centroid_representative(coords, np.where(ranked == r)[0])
            for r in range(len(populations))
        ]
        min_samples = config.cluster_min_samples
    else:  # unreachable — enum is exhaustive
        raise ValueError(f"unknown cluster_method: {config.cluster_method}")

    representative_frames = [int(frame_indices[p]) for p in rep_positions]
    assert not representative_frames or max(representative_frames) < traj.n_frames
    representative_rmsd = [float(rmsd_to_design[f]) for f in representative_frames]

    # Design-like basin: the cluster whose representative is closest to the design pose.
    design_like = int(np.argmin(representative_rmsd)) if representative_rmsd else -1
    design_like_fraction = (populations[design_like] / m) if design_like >= 0 else 0.0

    return ClusterResult(
        labels=ranked,
        frame_indices=frame_indices,
        populations=populations,
        representative_frames=representative_frames,
        representative_rmsd_to_design_ang=representative_rmsd,
        design_like_cluster=design_like,
        design_like_fraction=float(design_like_fraction),
        n_clusters=len(populations),
        n_frames_clustered=m,
        cutoff_ang=config.cluster_cutoff_ang,
        method=config.cluster_method.value,
        stride=stride,
        min_samples=min_samples,
    )


def write_representative_pdbs(traj: md.Trajectory, result: ClusterResult, adir) -> list[Path]:
    """Write one PDB per cluster (rank order; noise excluded) — the distinct conformations."""
    adir = Path(adir)
    paths = []
    for rank, frame in enumerate(result.representative_frames):
        path = adir / f"cluster_rep_{rank:02d}.pdb"
        traj[frame].save_pdb(str(path))
        paths.append(path)
    return paths


def cluster_and_report(
    traj: md.Trajectory, reference: md.Trajectory, config: AnalysisConfig, adir, x, xlabel: str
) -> ClusterResult:
    """Cluster the ensemble; write the representative PDBs and the population + assignment plots
    into an `analysis/clustering/` subdirectory. Shared by the plain-MD and multistate composers."""
    cdir = Path(adir) / "clustering"
    cdir.mkdir(parents=True, exist_ok=True)
    result = cluster_conformations(traj, reference, config)
    write_representative_pdbs(traj, result, cdir)
    plots.plot_cluster_populations(
        result.populations, cdir / "cluster_populations.png", result.design_like_cluster
    )
    plots.plot_cluster_timeseries(
        np.asarray(x)[result.frame_indices], result.labels,
        cdir / "cluster_timeseries.png", xlabel,
    )
    return result


def cluster_summary(result: ClusterResult) -> dict:
    """JSON-ready cluster block for the analysis report."""
    m = result.n_frames_clustered
    return {
        "method": result.method,
        "cutoff_ang": result.cutoff_ang,
        "stride": result.stride,
        "min_samples": result.min_samples,
        "n_frames_clustered": m,
        "n_clusters": result.n_clusters,
        "n_noise": int((result.labels == -1).sum()),
        "populations": result.populations,
        "population_fractions": [p / m for p in result.populations],
        "representative_frames": result.representative_frames,
        "representative_rmsd_to_design_ang": result.representative_rmsd_to_design_ang,
        "design_like_cluster": result.design_like_cluster,
        "design_like_fraction": result.design_like_fraction,
    }
