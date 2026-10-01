"""Load and condition a production trajectory for analysis.

Everything measures drift from the **input pose** (the design), so the flow is: load
the raw trajectory, image molecules whole across PBC, strip to protein, and superpose
onto the input structure. The topology is parsed from the mmCIF via OpenMM (mdtraj's
own CIF support is unreliable) and converted to an mdtraj topology.
"""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
from openmm.app import PDBxFile


def load_system_trajectory(outdir) -> md.Trajectory:
    """Load `production/trajectory.xtc` with the `prepared/topology.cif` topology, and
    image molecules whole (undo PBC splits). Full solvated system."""
    outdir = Path(outdir)
    omm_topology = PDBxFile(str(outdir / "prepared" / "topology.cif")).topology
    topology = md.Topology.from_openmm(omm_topology)
    traj = md.load(str(outdir / "production" / "trajectory.xtc"), top=topology)
    traj.image_molecules(inplace=True)
    return traj


def to_protein(traj: md.Trajectory) -> md.Trajectory:
    """Strip solvent/ions — keep protein atoms only."""
    sel = traj.topology.select("protein")
    assert len(sel) > 0, "no protein atoms in trajectory"
    return traj.atom_slice(sel)


def load_input_reference(outdir) -> md.Trajectory:
    """The design pose: `input.pdb`, protein atoms only (single frame reference)."""
    ref = md.load(str(Path(outdir) / "input.pdb"))
    return to_protein(ref)


def superpose_to_reference(
    traj: md.Trajectory, reference: md.Trajectory, selection: str = "name CA"
) -> md.Trajectory:
    """Least-squares-fit every frame onto the reference using `selection` (default Cα).

    Fails loud if the selection counts differ between trajectory and reference (atom
    mapping is wrong — e.g. residue mismatch after structure fixing).
    """
    t_idx = traj.topology.select(selection)
    r_idx = reference.topology.select(selection)
    assert len(t_idx) == len(r_idx) and len(t_idx) > 0, (
        f"selection {selection!r} mismatch: trajectory has {len(t_idx)} atoms, "
        f"reference has {len(r_idx)} — residues do not correspond"
    )
    traj.superpose(reference, atom_indices=t_idx, ref_atom_indices=r_idx)
    return traj


def load_protein_trajectory(outdir, selection: str = "name CA"):
    """Convenience: load → image → strip to protein → superpose onto the input pose.

    Returns ``(protein_trajectory, protein_reference)``, both ready for metrics.
    """
    traj = to_protein(load_system_trajectory(outdir))
    reference = load_input_reference(outdir)
    superpose_to_reference(traj, reference, selection)
    return traj, reference
