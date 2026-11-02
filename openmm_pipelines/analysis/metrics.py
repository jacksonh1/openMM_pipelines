"""Per-frame / per-residue drift metrics, all in angstroms.

- `rmsd_to_reference` — backbone-Cα RMSD of each frame to the **input design pose**
  (the reference), mapping Cα across the two topologies explicitly (the trajectory
  carries hydrogens the input PDB may not).
- `radius_of_gyration` — protein Rg per frame.
- `rmsf_per_residue` — Cα fluctuation about the trajectory mean (intrinsic
  flexibility; does not need the external reference).

Three **alignment-free** fold-retention metrics vs the design pose (`native_reference`, which
must share the trajectory's atom layout — use `trajectory.load_prepared_reference`):

- `native_contact_fraction` — Q (Best, Hummer & Eaton, PNAS 2013): soft fraction of native
  all-heavy-atom contacts retained.
- `ca_drmsd` — distance-RMSD: the superposition-free continuous analog of RMSD.
- `lddt_to_reference` — lDDT (Mariani et al., Bioinformatics 2013) over a chosen atom set, global
  and per-residue, via `biotite.structure.lddt`.

All three use only intra-structure distances, so they need no superposition and are invariant to it.
"""

from __future__ import annotations

import biotite.structure as struc
import mdtraj as md
import numpy as np

from .config import LddtAtomSet

_NM_TO_ANG = 10.0

# lDDT atom set -> (mdtraj selection, human label). Single source of truth for selection + display.
LDDT_SPEC: dict[LddtAtomSet, tuple[str, str]] = {
    LddtAtomSet.CA_CB: ("name CA CB", "Cα+Cβ"),       # Gly has no CB -> contributes Cα only
    LddtAtomSet.CA: ("name CA", "Cα"),
    LddtAtomSet.BACKBONE: ("name N CA C O", "backbone"),
    LddtAtomSet.HEAVY: ("not element H", "all-heavy"),
}


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


# --- alignment-free fold-retention vs the design pose ---------------------------------------
# `native_reference` must share `traj`'s atom layout (same topology.cif origin); these ASSUME it.


def _native_contact_pairs(
    reference: md.Trajectory, cutoff_nm: float, min_seq_sep: int
) -> tuple[np.ndarray, np.ndarray]:
    """Heavy-atom native contacts of the design: ``(pairs (n, 2) atom indices, r0 (n,) nm)``.

    A native contact is a pair of heavy atoms from residues at least `min_seq_sep` apart in
    sequence (ordinal) whose distance in `reference` frame 0 is below `cutoff_nm`.
    """
    heavy = reference.topology.select("not element H")
    assert len(heavy) > 1, "need >= 2 heavy atoms to define native contacts"
    res = np.array([reference.topology.atom(int(i)).residue.index for i in heavy])
    ii, jj = np.triu_indices(len(heavy), k=1)  # all unordered heavy-atom pairs (positions)
    keep = np.abs(res[ii] - res[jj]) >= min_seq_sep
    pairs = np.column_stack([heavy[ii[keep]], heavy[jj[keep]]])
    r0 = md.compute_distances(reference, pairs)[0]  # nm, frame 0
    native = r0 < cutoff_nm
    return pairs[native], r0[native]


def native_contact_fraction(
    traj: md.Trajectory,
    native_reference: md.Trajectory,
    *,
    cutoff_ang: float = 4.5,
    min_seq_sep: int = 4,
    beta_per_nm: float = 50.0,
    lam: float = 1.8,
) -> np.ndarray:
    """Per-frame fraction of native contacts Q (Best, Hummer & Eaton, PNAS 2013).

    Native contacts are fixed on `native_reference` (the design): all-heavy-atom pairs from
    residues ``>= min_seq_sep`` apart in sequence within `cutoff_ang`. Each contact's on/off is
    softened by the sigmoid ``1 / (1 + exp[beta (r - lam r0)])`` so Q is smooth and tolerant of
    thermal breathing; the defaults are the paper's all-atom values (beta = 5/Å = 50/nm,
    lambda = 1.8). Returns Q per frame in [0, 1].
    """
    pairs, r0_nm = _native_contact_pairs(native_reference, cutoff_ang / _NM_TO_ANG, min_seq_sep)
    assert len(pairs) > 0, (
        f"no native contacts (cutoff {cutoff_ang} Å, seq-sep {min_seq_sep}) — check the reference"
    )
    r_nm = md.compute_distances(traj, pairs)  # (n_frames, n_contacts)
    q = 1.0 / (1.0 + np.exp(beta_per_nm * (r_nm - lam * r0_nm)))
    return q.mean(axis=1)


def ca_drmsd(
    traj: md.Trajectory, native_reference: md.Trajectory, *, selection: str = "name CA"
) -> np.ndarray:
    """Per-frame distance-RMSD (Å) to the design: ``sqrt(mean_{i<j} (d_ij - d_ij^0)^2)`` over all
    `selection` pairs. Superposition-free (internal distances only), continuous, unbounded."""
    t_idx = traj.topology.select(selection)
    r_idx = native_reference.topology.select(selection)
    assert len(t_idx) == len(r_idx) > 1, (
        f"selection {selection!r} mismatch or too few atoms: {len(t_idx)} vs {len(r_idx)}"
    )
    ii, jj = np.triu_indices(len(t_idx), k=1)
    d0 = md.compute_distances(native_reference, np.column_stack([r_idx[ii], r_idx[jj]]))[0]  # nm
    d = md.compute_distances(traj, np.column_stack([t_idx[ii], t_idx[jj]]))  # (n_frames, n_pairs) nm
    return np.sqrt(((d - d0) ** 2).mean(axis=1)) * _NM_TO_ANG


def _biotite_reference(reference: md.Trajectory, atom_indices: np.ndarray) -> struc.AtomArray:
    """Build a biotite `AtomArray` (coords in Å) for `atom_indices` of a single-frame mdtraj
    `reference`, carrying the annotations lDDT needs (chain, residue, atom name, element)."""
    top = reference.topology
    atoms = [top.atom(int(i)) for i in atom_indices]
    arr = struc.AtomArray(len(atoms))
    arr.coord = (reference.xyz[0, atom_indices, :] * _NM_TO_ANG).astype(np.float32)  # nm -> Å
    arr.chain_id = np.array([str(a.residue.chain.index) for a in atoms])
    arr.res_id = np.array([a.residue.resSeq for a in atoms], dtype=int)
    arr.res_name = np.array([a.residue.name for a in atoms])
    arr.atom_name = np.array([a.name for a in atoms])
    arr.element = np.array([a.element.symbol if a.element is not None else "" for a in atoms])
    return arr


def lddt_to_reference(
    traj: md.Trajectory,
    native_reference: md.Trajectory,
    selection: str,
    *,
    inclusion_radius_ang: float = 15.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Per-frame lDDT (Mariani et al., Bioinformatics 2013) of `traj` to the design over atoms
    matching `selection`, via `biotite.structure.lddt`.

    Whole-system (all local pairs within `inclusion_radius_ang`, only same-residue pairs excluded);
    thresholds fixed at the GDT-HA set 0.5/1/2/4 Å. `native_reference` must share `traj`'s atom
    layout. Returns ``(global (n_frames,), per_residue (n_frames, n_res), resseq (n_res,))``, where
    the residue axis is in trajectory order.

    NOTE: biotite does not resolve the symmetric-residue atom-naming ambiguity (Asp/Glu/Val/Tyr/
    Leu/Phe/Arg), so an all-heavy lDDT reads a physical 180° side-chain flip as a small dip; coarser
    atom sets (Cα, Cα+Cβ, backbone) are immune.
    """
    t_idx = traj.topology.select(selection)
    r_idx = native_reference.topology.select(selection)
    assert len(t_idx) == len(r_idx) > 0, (
        f"selection {selection!r} mismatch: trajectory {len(t_idx)} vs reference {len(r_idx)}"
    )
    assert np.array_equal(t_idx, r_idx), (
        "trajectory and reference atom layouts differ — pass load_prepared_reference()"
    )
    ref_aa = _biotite_reference(native_reference, t_idx)
    subject = (traj.xyz[:, t_idx, :] * _NM_TO_ANG).astype(np.float32)  # (n_frames, n_sel, 3) Å
    glob = np.asarray(struc.lddt(ref_aa, subject, inclusion_radius=inclusion_radius_ang))
    per_res = np.asarray(
        struc.lddt(ref_aa, subject, aggregation="residue", inclusion_radius=inclusion_radius_ang)
    )
    # Distinct residues in first-appearance (trajectory) order — the order biotite aggregates in.
    ordered: list[int] = []
    seen: set[int] = set()
    for i in t_idx:
        ri = native_reference.topology.atom(int(i)).residue.index
        if ri not in seen:
            seen.add(ri)
            ordered.append(ri)
    resseq = np.array([native_reference.topology.residue(ri).resSeq for ri in ordered])
    assert per_res.shape[-1] == len(resseq), (
        f"lDDT per-residue width {per_res.shape[-1]} != residue count {len(resseq)}"
    )
    return glob, per_res, resseq
