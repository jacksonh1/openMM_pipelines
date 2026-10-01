"""Harmonic position restraints on a chosen set of atoms.

The standard OpenMM cookbook idiom: a `CustomExternalForce` whose stiffness `k` is a
*global* parameter, so the restraint is released with a single
`context.setParameter("k", 0.0)` — no force removal, no context rebuild. Reference
positions `x0,y0,z0` are per-particle parameters; `periodicdistance` keeps the
restraint well-defined across periodic boundaries.
"""

from __future__ import annotations

import openmm
from openmm import unit


def add_restraint(system, positions, atom_indices, k):
    """Add a harmonic position restraint to `system` and return the force.

    `k` is in kJ/mol/nm^2. Reference positions are taken from `positions` at each
    index in `atom_indices` (add-order == iteration-order, which `reanchor` relies
    on). Fails loudly on an empty selection — a restraint on nothing is a bug.
    """
    atom_indices = list(atom_indices)
    assert atom_indices, "add_restraint got an empty atom_indices selection"

    force = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    for i in atom_indices:
        xyz = positions[i].value_in_unit(unit.nanometer)
        force.addParticle(i, [xyz[0], xyz[1], xyz[2]])
    system.addForce(force)
    return force


def reanchor(force, context, positions, atom_indices):
    """Move the restraint reference positions to the current coordinates.

    Call once after minimization so the restraint reference is the minimized pose and
    never drifts during equilibration. `atom_indices` must be the *same sequence*
    passed to `add_restraint` (same order), so each restraint slot maps back to its
    atom.
    """
    atom_indices = list(atom_indices)
    assert force.getNumParticles() == len(atom_indices), (
        f"reanchor atom_indices ({len(atom_indices)}) does not match the force's "
        f"particle count ({force.getNumParticles()}) — wrong selection or order"
    )
    for slot, atom_idx in enumerate(atom_indices):
        xyz = positions[atom_idx].value_in_unit(unit.nanometer)
        force.setParticleParameters(slot, atom_idx, [xyz[0], xyz[1], xyz[2]])
    force.updateParametersInContext(context)
