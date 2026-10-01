"""Restraint force behaviour: energy zero at reference, released by k=0, reanchorable.

Uses a trivial 3-particle System with a Reference platform — no force field, no PDB.
"""

import openmm
import pytest
from openmm import unit

from openmm_pipelines.lib.restraints import add_restraint, reanchor

K = 1000.0  # kJ/mol/nm^2


def _system(n):
    system = openmm.System()
    for _ in range(n):
        system.addParticle(1.0 * unit.amu)
    return system


def _context(system):
    integ = openmm.VerletIntegrator(1.0 * unit.femtosecond)
    platform = openmm.Platform.getPlatformByName("Reference")
    return openmm.Context(system, integ, platform)


def _restraint_energy(context):
    return (
        context.getState(getEnergy=True, groups={_GROUP})
        .getPotentialEnergy()
        .value_in_unit(unit.kilojoule_per_mole)
    )


_GROUP = 1  # isolate the restraint's energy from any other force


def _positions(coords):
    return [openmm.Vec3(*c) for c in coords] * unit.nanometer


def test_energy_zero_at_reference():
    system = _system(2)
    pos = _positions([(0, 0, 0), (1, 0, 0)])
    force = add_restraint(system, pos, [0, 1], K)
    force.setForceGroup(_GROUP)
    ctx = _context(system)
    ctx.setPositions(pos)
    assert _restraint_energy(ctx) == pytest.approx(0.0, abs=1e-6)


def test_energy_quadratic_when_displaced():
    system = _system(1)
    pos = _positions([(0, 0, 0)])
    force = add_restraint(system, pos, [0], K)
    force.setForceGroup(_GROUP)
    ctx = _context(system)
    ctx.setPositions(_positions([(0.1, 0, 0)]))  # 0.1 nm off reference
    # k * d^2 = 1000 * 0.01 = 10 kJ/mol
    assert _restraint_energy(ctx) == pytest.approx(10.0, rel=1e-5)


def test_release_with_k_zero():
    system = _system(1)
    pos = _positions([(0, 0, 0)])
    force = add_restraint(system, pos, [0], K)
    force.setForceGroup(_GROUP)
    ctx = _context(system)
    ctx.setPositions(_positions([(0.1, 0, 0)]))
    assert _restraint_energy(ctx) > 1.0
    ctx.setParameter("k", 0.0)
    assert _restraint_energy(ctx) == pytest.approx(0.0, abs=1e-6)


def test_reanchor_moves_reference():
    system = _system(1)
    pos = _positions([(0, 0, 0)])
    force = add_restraint(system, pos, [0], K)
    force.setForceGroup(_GROUP)
    ctx = _context(system)
    moved = _positions([(0.1, 0, 0)])
    ctx.setPositions(moved)
    assert _restraint_energy(ctx) > 1.0
    reanchor(force, ctx, moved, [0])  # reference is now the displaced pose
    assert _restraint_energy(ctx) == pytest.approx(0.0, abs=1e-6)


def test_empty_selection_raises():
    system = _system(1)
    pos = _positions([(0, 0, 0)])
    with pytest.raises(AssertionError):
        add_restraint(system, pos, [], K)


def test_reanchor_wrong_length_raises():
    system = _system(2)
    pos = _positions([(0, 0, 0), (1, 0, 0)])
    force = add_restraint(system, pos, [0, 1], K)
    ctx = _context(system)
    ctx.setPositions(pos)
    with pytest.raises(AssertionError):
        reanchor(force, ctx, pos, [0])  # only one index for a two-particle force
