"""Per-frame / per-residue drift metrics, all in angstroms.

- `rmsd_to_reference` — backbone-Cα RMSD of each frame to the **input design pose**
  (the reference), mapping Cα across the two topologies explicitly (the trajectory
  carries hydrogens the input PDB may not).
- `radius_of_gyration` — protein Rg per frame.
- `rmsf_per_residue` — Cα fluctuation about the trajectory mean (intrinsic
  flexibility; does not need the external reference).
"""

from __future__ import annotations

import mdtraj as md
import numpy as np

_NM_TO_ANG = 10.0


def rmsd_to_reference(
    traj: md.Trajectory, reference: md.Trajectory, selection: str = "name CA"
) -> np.ndarray:
    """Cα RMSD (Å) of each frame to the reference's frame 0."""
    t_idx = traj.topology.select(selection)
    r_idx = reference.topology.select(selection)
    assert len(t_idx) == len(r_idx) and len(t_idx) > 0, (
        f"selection {selection!r} mismatch: {len(t_idx)} vs {len(r_idx)} atoms"
    )
    rmsd_nm = md.rmsd(traj, reference, frame=0, atom_indices=t_idx, ref_atom_indices=r_idx)
    return rmsd_nm * _NM_TO_ANG


def radius_of_gyration(traj: md.Trajectory) -> np.ndarray:
    """Protein radius of gyration (Å) per frame."""
    return md.compute_rg(traj) * _NM_TO_ANG


def rmsf_per_residue(
    traj: md.Trajectory, selection: str = "name CA"
) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue Cα RMSF (Å) about the trajectory mean. Returns ``(resSeq, rmsf)``.

    ASSUMES the trajectory is already superposed onto the reference (so the fit frame
    is consistent); fluctuation is then measured about the per-atom mean position.
    """
    ca = traj.topology.select(selection)
    assert len(ca) > 0, f"selection {selection!r} matched no atoms"
    xyz = traj.xyz[:, ca, :]  # (frames, n_ca, 3) in nm
    mean = xyz.mean(axis=0)
    rmsf_nm = np.sqrt(((xyz - mean) ** 2).sum(axis=2).mean(axis=0))
    resseq = np.array([traj.topology.atom(i).residue.resSeq for i in ca])
    return resseq, rmsf_nm * _NM_TO_ANG
