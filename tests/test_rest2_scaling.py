"""REST2 scale_solute correctness on toy Systems — no force field, Reference platform.

Pins the scaling contract directly on force parameters (and an energy identity at λ=1):
solute charge -> √λ·q, solute ε -> λ·ε (σ unchanged), solute-only bonded k -> λ·k, CMAP
maps -> λ, NonbondedForce exceptions rescaled with the same combining rules, solvent
untouched, and a loud refusal on force types whose scaling is not implemented.
"""

import math

import openmm
import pytest
from openmm import unit

from openmm_pipelines.rest2.scaling import scale_solute

LAM = 0.49  # sqrt = 0.7, a clean non-trivial factor
SOLUTE = [0, 1]  # particles 0,1 are solute; 2,3 are solvent


def _charge(force, i):
    return force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)


def _epsilon(force, i):
    return force.getParticleParameters(i)[2].value_in_unit(unit.kilojoule_per_mole)


def _sigma(force, i):
    return force.getParticleParameters(i)[1].value_in_unit(unit.nanometer)


def _nonbonded_system():
    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0 * unit.amu)
    nb = openmm.NonbondedForce()
    for i in range(4):
        nb.addParticle((i + 1) * 0.1, 0.3, 0.5)  # charge, sigma, epsilon
    nb.addException(0, 1, 0.08, 0.3, 0.4)  # solute-solute 1-4-style exception
    system.addForce(nb)
    return system, nb


def test_lambda_one_is_identity_parameters():
    system, nb = _nonbonded_system()
    scaled = scale_solute(system, SOLUTE, 1.0)
    snb = scaled.getForce(0)
    for i in range(4):
        assert _charge(snb, i) == pytest.approx(_charge(nb, i))
        assert _epsilon(snb, i) == pytest.approx(_epsilon(nb, i))


def test_solute_charge_and_epsilon_scaled():
    system, _ = _nonbonded_system()
    scaled = scale_solute(system, SOLUTE, LAM)  # keep System alive (force proxy references it)
    snb = scaled.getForce(0)
    sqrt_lam = math.sqrt(LAM)
    # solute scaled, solvent untouched
    assert _charge(snb, 0) == pytest.approx(0.1 * sqrt_lam)
    assert _epsilon(snb, 0) == pytest.approx(0.5 * LAM)
    assert _sigma(snb, 0) == pytest.approx(0.3)  # sigma never scaled
    assert _charge(snb, 2) == pytest.approx(0.3)  # solvent charge untouched
    assert _epsilon(snb, 2) == pytest.approx(0.5)


def test_solute_solute_exception_scaled_by_lambda():
    system, _ = _nonbonded_system()
    scaled = scale_solute(system, SOLUTE, LAM)
    snb = scaled.getForce(0)
    p1, p2, cprod, sigma, eps = snb.getExceptionParameters(0)
    # both atoms solute: chargeProd ∝ (√λ)^2 = λ, ε ∝ √(λ·λ) = λ
    assert cprod.value_in_unit(unit.elementary_charge**2) == pytest.approx(0.08 * LAM)
    assert eps.value_in_unit(unit.kilojoule_per_mole) == pytest.approx(0.4 * LAM)


def test_lambda_one_energy_identity():
    system, _ = _nonbonded_system()
    scaled = scale_solute(system, SOLUTE, 1.0)
    pos = [openmm.Vec3(0, 0, 0), openmm.Vec3(0.2, 0, 0),
           openmm.Vec3(0, 0.5, 0), openmm.Vec3(0.5, 0.5, 0)] * unit.nanometer

    def energy(s):
        ctx = openmm.Context(s, openmm.VerletIntegrator(1 * unit.femtosecond),
                             openmm.Platform.getPlatformByName("Reference"))
        ctx.setPositions(pos)
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)

    assert energy(scaled) == pytest.approx(energy(system), rel=1e-9)


def test_original_system_untouched():
    system, nb = _nonbonded_system()
    scale_solute(system, SOLUTE, LAM)
    assert _charge(nb, 0) == pytest.approx(0.1)  # caller's system not mutated


def test_harmonic_bond_scaled_only_within_solute():
    system = openmm.System()
    for _ in range(4):
        system.addParticle(12.0 * unit.amu)
    bond = openmm.HarmonicBondForce()
    bond.addBond(0, 1, 0.15, 1000.0)  # solute-solute -> scaled
    bond.addBond(2, 3, 0.15, 2000.0)  # solvent-solvent -> untouched
    system.addForce(bond)
    scaled = scale_solute(system, SOLUTE, LAM)
    sb = scaled.getForce(0)
    assert sb.getBondParameters(0)[3].value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer**2) == pytest.approx(1000.0 * LAM)
    assert sb.getBondParameters(1)[3].value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer**2) == pytest.approx(2000.0)


def test_periodic_torsion_scaled_only_within_solute():
    system = openmm.System()
    for _ in range(5):
        system.addParticle(12.0 * unit.amu)
    tor = openmm.PeriodicTorsionForce()
    tor.addTorsion(0, 1, 2, 3, 2, 0.0, 5.0)   # all solute? 2,3 not in SOLUTE -> no
    system.addForce(tor)
    # with SOLUTE={0,1,2,3}, the torsion is fully inside solute and scales
    s_all = scale_solute(system, [0, 1, 2, 3], LAM)
    assert s_all.getForce(0).getTorsionParameters(0)[6].value_in_unit(
        unit.kilojoule_per_mole) == pytest.approx(5.0 * LAM)
    # with SOLUTE={0,1} the torsion crosses the boundary and is left alone
    s_part = scale_solute(system, [0, 1], LAM)
    assert s_part.getForce(0).getTorsionParameters(0)[6].value_in_unit(
        unit.kilojoule_per_mole) == pytest.approx(5.0)


def _cmap_system(n_particles=8):
    system = openmm.System()
    for _ in range(n_particles):
        system.addParticle(12.0 * unit.amu)
    cmap = openmm.CMAPTorsionForce()
    size = 2
    energy = [1.0, 2.0, 3.0, 4.0]  # size*size grid
    m = cmap.addMap(size, energy)
    cmap.addTorsion(m, 0, 1, 2, 3, 4, 5, 6, 7)
    system.addForce(cmap)
    return system, cmap


def test_cmap_maps_scaled_when_torsions_inside_solute():
    system, _ = _cmap_system()
    scaled = scale_solute(system, list(range(8)), LAM)
    sc = scaled.getForce(0)
    _, energy = sc.getMapParameters(0)
    got = [e.value_in_unit(unit.kilojoule_per_mole) for e in energy]
    assert got == pytest.approx([e * LAM for e in (1.0, 2.0, 3.0, 4.0)])


def test_cmap_refuses_when_torsion_reaches_outside_solute():
    system, _ = _cmap_system()
    with pytest.raises(AssertionError):
        scale_solute(system, list(range(7)), LAM)  # atom 7 left out of the solute


def test_refuses_unknown_force_type():
    system = openmm.System()
    for _ in range(2):
        system.addParticle(12.0 * unit.amu)
    cf = openmm.CustomBondForce("k*r^2")  # stand-in for CHARMM NBFIX LJ bonds
    cf.addPerBondParameter("k")
    cf.addBond(0, 1, [1.0])
    system.addForce(cf)
    with pytest.raises(NotImplementedError):
        scale_solute(system, [0, 1], LAM)


def _external_restraint_system(k_default):
    system = openmm.System()
    for _ in range(2):
        system.addParticle(12.0 * unit.amu)
    f = openmm.CustomExternalForce("k*(x-x0)^2")
    f.addGlobalParameter("k", k_default)
    f.addPerParticleParameter("x0")
    f.addParticle(0, [0.0])
    system.addForce(f)
    return system


def test_released_external_restraint_ignored():
    system = _external_restraint_system(0.0)  # released (prepare bakes k=0)
    scaled = scale_solute(system, [0, 1], LAM)  # no raise
    assert scaled.getNumForces() == 1


def test_active_external_restraint_refused():
    system = _external_restraint_system(1000.0)  # still-active restraint
    with pytest.raises(AssertionError):
        scale_solute(system, [0, 1], LAM)


def test_rejects_bad_lambda():
    system, _ = _nonbonded_system()
    for bad in (0.0, -0.1, 1.5):
        with pytest.raises(AssertionError):
            scale_solute(system, SOLUTE, bad)
