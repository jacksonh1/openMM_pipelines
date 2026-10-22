"""Multi-chain / interchain observables for complexes (objective #4: bound-state sampling).

Auto-run by the composers when the protein has >=2 chains. All distance observables are
minimum-image-based, so they are correct regardless of how whole molecules were imaged:

- **interchain minimum distance** (closest heavy-atom distance per chain pair) — the
  binding/unbinding observable (contact ~1.5 Å; a sustained rise = separation).
- **interchain contact count** (inter-chain residue pairs within `contact_cutoff_nm`,
  closest-heavy) — interface size / stability.
- **interchain center-of-mass distance** (minimum-image) — overall chain separation.
- **per-chain internal RMSD** — each chain's Cα RMSD to the design, fit to that chain alone,
  so it reports fold retention independent of docking drift.

Chain correspondence between the trajectory and the design `reference` is **positional**:
OpenMM prep preserves chain order (PDBFixer edits within chains; `Modeller.addSolvent`
appends solvent after the protein), and both are stripped to protein-only, so protein chain
i in the trajectory is protein chain i in the reference. We assert the counts agree and
fail loud otherwise — no content-based remap (unlike GROMACS, which renumbers).

Pure mdtraj; no OpenMM import.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import mdtraj as md
import numpy as np

from . import plots
from .config import AnalysisConfig

_NM_TO_ANG = 10.0


def protein_chain_ids(traj: md.Trajectory) -> list[int]:
    """Indices of chains containing protein atoms (mdtraj's `protein` selection), ascending."""
    sel = traj.topology.select("protein")
    assert len(sel) > 0, "no protein atoms — cannot identify protein chains"
    return sorted({traj.topology.atom(int(i)).residue.chain.index for i in sel})


def chain_pairs(chain_ids: list[int]) -> list[tuple[int, int]]:
    """All unordered pairs of chains."""
    return list(combinations(chain_ids, 2))


def assert_chain_correspondence(traj: md.Trajectory, reference: md.Trajectory) -> list[int]:
    """Positional traj<->reference chain match: equal protein-chain count and equal per-chain
    Cα count. Returns the shared protein chain ids. Fails loud on any mismatch."""
    t_chains = protein_chain_ids(traj)
    r_chains = protein_chain_ids(reference)
    assert t_chains == r_chains, (
        f"protein chains differ between trajectory {t_chains} and reference {r_chains} — "
        "chain order/count must match (both are protein-only, prep preserves order)"
    )
    for c in t_chains:
        nt = len(traj.topology.select(f"chainid {c} and name CA"))
        nr = len(reference.topology.select(f"chainid {c} and name CA"))
        assert nt == nr > 0, (
            f"chain {c} Cα count mismatch: trajectory {nt} vs reference {nr}"
        )
    return t_chains


def _chain_residue_indices(traj: md.Trajectory, chain_id: int) -> list[int]:
    return [res.index for res in traj.topology.chain(chain_id).residues]  # .residues is an iterator property


def interchain_residue_distances(
    traj: md.Trajectory, chain_a: int, chain_b: int
) -> np.ndarray:
    """Per-frame closest-heavy distances (nm) for every cross-chain residue pair, shape
    (n_frames, n_pairs). Minimum-image (`periodic=True`). Compute once, derive min & count."""
    res_a = _chain_residue_indices(traj, chain_a)
    res_b = _chain_residue_indices(traj, chain_b)
    pairs = [(ra, rb) for ra in res_a for rb in res_b]
    assert pairs, f"no residue pairs between chains {chain_a} and {chain_b}"
    dist, _ = md.compute_contacts(traj, contacts=pairs, scheme="closest-heavy", periodic=True)
    return dist  # nm


def interchain_com_distance(
    traj: md.Trajectory, chain_a: int, chain_b: int
) -> np.ndarray:
    """Per-frame minimum-image center-of-mass distance (Å) between two chains.

    COM has no built-in minimum image, so we reduce the per-frame COM displacement into the
    (triclinic) unit cell via fractional coordinates — correct across a periodic boundary.
    """
    assert traj.unitcell_vectors is not None, "interchain COM distance needs box vectors"
    com_a = md.compute_center_of_mass(traj.atom_slice(traj.topology.select(f"chainid {chain_a}")))
    com_b = md.compute_center_of_mass(traj.atom_slice(traj.topology.select(f"chainid {chain_b}")))
    d = com_a - com_b  # (n_frames, 3), nm
    box = traj.unitcell_vectors  # (n_frames, 3, 3); rows are the box vectors a, b, c
    frac = np.einsum("fj,fjk->fk", d, np.linalg.inv(box))  # Cartesian -> fractional
    frac -= np.round(frac)  # minimum image
    d_min = np.einsum("fk,fkl->fl", frac, box)  # back to Cartesian
    return np.linalg.norm(d_min, axis=1) * _NM_TO_ANG


def per_chain_rmsd(
    traj: md.Trajectory, reference: md.Trajectory, chain_id: int, selection: str = "name CA"
) -> np.ndarray:
    """Per-frame Cα RMSD (Å) of one chain to the design, fit to THAT chain alone.

    `md.rmsd` aligns using only this chain's atoms and does not mutate `traj`, so the global
    superposition is untouched — this isolates intra-chain fold drift from docking motion.
    """
    t_idx = traj.topology.select(f"chainid {chain_id} and {selection}")
    r_idx = reference.topology.select(f"chainid {chain_id} and {selection}")
    assert len(t_idx) == len(r_idx) > 0, (
        f"chain {chain_id} selection {selection!r} mismatch: {len(t_idx)} vs {len(r_idx)}"
    )
    rmsd_nm = md.rmsd(traj, reference, frame=0, atom_indices=t_idx, ref_atom_indices=r_idx)
    return rmsd_nm * _NM_TO_ANG


def _stats(arr: np.ndarray) -> dict:
    return {"mean": float(arr.mean()), "min": float(arr.min()),
            "max": float(arr.max()), "final": float(arr[-1])}


def interchain_report(
    traj: md.Trajectory, reference: md.Trajectory, adir, x, xlabel: str, config: AnalysisConfig
) -> tuple[dict, dict]:
    """Interchain + per-chain metrics vs `reference`, write the PNGs into `adir`, and return
    ``(summary, arrays)`` — the same shape as `drift_report`, so composers fold them in.
    """
    adir = Path(adir)
    chains = assert_chain_correspondence(traj, reference)
    pairs = chain_pairs(chains)

    min_dist, contacts, com_dist = {}, {}, {}
    pair_summary = {}
    for a, b in pairs:
        key = f"{a}-{b}"
        d_nm = interchain_residue_distances(traj, a, b)
        min_dist[key] = d_nm.min(axis=1) * _NM_TO_ANG  # Å
        contacts[key] = (d_nm < config.contact_cutoff_nm).sum(axis=1).astype(float)
        com_dist[key] = interchain_com_distance(traj, a, b)
        pair_summary[key] = {
            "min_distance_ang": _stats(min_dist[key]),
            "contact_count": _stats(contacts[key]),
            "com_distance_ang": _stats(com_dist[key]),
        }

    per_chain = {str(c): per_chain_rmsd(traj, reference, c, config.cluster_selection)
                 for c in chains}
    chain_summary = {c: _stats(r) for c, r in per_chain.items()}

    plots.plot_series(x, min_dist, "interchain min distance (Å)",
                      adir / "interchain_mindist.png", "interchain minimum distance", xlabel)
    plots.plot_series(x, contacts, "interchain contacts",
                      adir / "interchain_contacts.png",
                      f"interchain contacts (< {config.contact_cutoff_nm} nm)", xlabel)
    plots.plot_series(x, com_dist, "interchain COM distance (Å)",
                      adir / "interchain_com.png", "interchain centre-of-mass distance", xlabel)
    plots.plot_series(x, per_chain, "Cα RMSD to design (Å)",
                      adir / "per_chain_rmsd.png", "per-chain fold retention", xlabel)

    summary = {"n_protein_chains": len(chains),
               "chain_pairs": pair_summary,
               "per_chain_rmsd_ang": chain_summary}
    arrays = {}
    for key in min_dist:
        a_b = key.replace("-", "_")
        arrays[f"interchain_mindist_ang_{a_b}"] = min_dist[key]
        arrays[f"interchain_contacts_{a_b}"] = contacts[key]
        arrays[f"interchain_com_ang_{a_b}"] = com_dist[key]
    for c, r in per_chain.items():
        arrays[f"chain{c}_rmsd_ang"] = r
    return summary, arrays
