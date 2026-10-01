"""Conformational clustering by pairwise Cα RMSD (average-linkage hierarchical).

Uses scipy for the linkage/flat-cluster — no hand-rolled clustering.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

_NM_TO_ANG = 10.0


def pairwise_rmsd(traj: md.Trajectory, selection: str = "name CA") -> np.ndarray:
    """Symmetric (n_frames, n_frames) Cα RMSD matrix in Å."""
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


def cluster_conformations(
    traj: md.Trajectory, cutoff_ang: float = 2.0, selection: str = "name CA"
) -> np.ndarray:
    """Flat cluster labels (1..k) joining frames within `cutoff_ang` Å (average link)."""
    assert traj.n_frames >= 2, "need at least 2 frames to cluster"
    dist = pairwise_rmsd(traj, selection)
    linkage_matrix = linkage(squareform(dist, checks=False), method="average")
    return fcluster(linkage_matrix, t=cutoff_ang, criterion="distance")
