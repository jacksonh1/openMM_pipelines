"""scale_solute — REST2 solute Hamiltonian scaling (the one genuinely new REST2 piece).

REST2 (Wang-Friesner-Berne 2011) tempers only the *solute*: every replica runs at the
same physical temperature T_0, but replica m's solute interactions are scaled by
λ_m = T_0 / T_m ∈ (0, 1] so only solute degrees of freedom enter the exchange acceptance.
This primitive takes a System + the solute atom indices + λ and returns an independent
copy with the solute pre-scaled — "path A" from the kickoff note (N pre-scaled Systems,
one plain ThermodynamicState each). Single purpose, energy-testable in isolation.

The scaling (scale the force-field terms; all replicas share T_0):

- **Charges:** q_solute -> √λ · q. Coulomb solute-solute ∝ λ, solute-solvent ∝ √λ
  (falls out of the product q_i·q_j).
- **LJ ε:** ε_solute -> λ · ε (σ unchanged). Lorentz-Berthelot ε_ij = √(ε_i ε_j) gives
  solute-solute ∝ λ, solute-solvent ∝ √λ automatically.
- **NonbondedForce exceptions (1-4, exclusions):** OpenMM stores the *combined* chargeProd
  and ε per exception and does NOT re-derive them from per-particle params, so they must be
  rescaled explicitly with the same combining factors (solute-solute ∝ λ; a solute-solvent
  exception, which bonded topology should never produce, ∝ √λ).
- **Solute-only bonded terms** (bonds / angles / torsions fully inside the solute): k · λ.
- **CMAP (CHARMM36m backbone correction):** a solute intramolecular term — the map
  energies are scaled by λ. Forgetting it silently gives a wrong REST2, so this is gated:
  if any CMAP torsion reaches outside the solute, we refuse.

Solvent-solvent is left untouched. At λ = 1 every factor is 1, so the scaled System is
energetically identical to the input (the identity the tests pin).

**Force-field support:** the AMBER force set (HarmonicBond/Angle, PeriodicTorsion,
NonbondedForce) plus CMAP is handled. CHARMM36(m) additionally builds the LJ as
`CustomNonbondedForce` + `CustomBondForce` (NBFIX) and impropers as `CustomTorsionForce`;
scaling those correctly is not yet implemented, so a System carrying any unrecognized
force with solute atoms is **refused** rather than silently mis-scaled. See TODO.
"""

from __future__ import annotations

import math

from openmm import (
    CMAPTorsionForce,
    CMMotionRemover,
    CustomExternalForce,
    HarmonicAngleForce,
    HarmonicBondForce,
    MonteCarloAnisotropicBarostat,
    MonteCarloBarostat,
    MonteCarloMembraneBarostat,
    NonbondedForce,
    PeriodicTorsionForce,
    System,
)

from ..lib.system_edits import clone_system

# Forces that carry no solute-dependent potential energy — safe to leave as-is.
_IGNORED_FORCES = (
    CMMotionRemover,
    MonteCarloBarostat,
    MonteCarloAnisotropicBarostat,
    MonteCarloMembraneBarostat,
)


def scale_solute(system: System, solute_indices, lam: float) -> System:
    """Return an independent copy of `system` with the solute scaled by REST2 factor `lam`.

    `solute_indices` is the set of atoms to temper (whole protein, or a chain subset).
    `lam` = λ_m = T_0/T_m ∈ (0, 1]; λ=1 is the physical, unscaled Hamiltonian.
    Refuses (raises) on any force type whose correct REST2 scaling is not implemented.
    """
    assert 0.0 < lam <= 1.0, f"lambda must be in (0, 1], got {lam}"
    solute = frozenset(int(i) for i in solute_indices)
    assert solute, "empty solute selection — nothing to scale"

    scaled = clone_system(system)
    sqrt_lam = math.sqrt(lam)

    for force in (scaled.getForce(i) for i in range(scaled.getNumForces())):
        if isinstance(force, NonbondedForce):
            _scale_nonbonded(force, solute, lam, sqrt_lam)
        elif isinstance(force, HarmonicBondForce):
            _scale_harmonic_bond(force, solute, lam)
        elif isinstance(force, HarmonicAngleForce):
            _scale_harmonic_angle(force, solute, lam)
        elif isinstance(force, PeriodicTorsionForce):
            _scale_periodic_torsion(force, solute, lam)
        elif isinstance(force, CMAPTorsionForce):
            _scale_cmap(force, solute, lam)
        elif isinstance(force, CustomExternalForce):
            # The pipeline's position restraint (prepare() bakes it released, k=0). It is
            # an external restraint, not a force-field solute term, so it is NOT λ-scaled —
            # but it MUST be released, or we'd be running a restrained REST2.
            _assert_released_external(force)
        elif isinstance(force, _IGNORED_FORCES):
            continue
        else:
            raise NotImplementedError(
                f"REST2 scaling not implemented for {type(force).__name__}. The AMBER "
                "force set + CMAP is supported; CHARMM's Custom{Nonbonded,Bond,Torsion}Force "
                "(NBFIX LJ / impropers) is not yet — refusing rather than mis-scaling the "
                "solute. See rest2/scaling.py and knowledgebase/TODO.md."
            )
    return scaled


def _assert_released_external(force: CustomExternalForce) -> None:
    """A CustomExternalForce must be fully released (every global parameter default 0) to be
    carried into REST2 — an active external restraint would bias every replica. Refuse otherwise.
    """
    for i in range(force.getNumGlobalParameters()):
        name = force.getGlobalParameterName(i)
        val = force.getGlobalParameterDefaultValue(i)
        assert val == 0.0, (
            f"CustomExternalForce global parameter '{name}' = {val} != 0 — an ACTIVE "
            "external restraint under REST2. prepare() should bake it released (k=0) before "
            "production; refusing to scale a restrained system."
        )


def _scale_nonbonded(force: NonbondedForce, solute, lam: float, sqrt_lam: float) -> None:
    # Particle parameter offsets would make the per-particle charge/ε λ-dependent in ways
    # this static rescale can't track — refuse rather than silently ignore them.
    assert force.getNumParticleParameterOffsets() == 0, (
        "NonbondedForce has particle parameter offsets — unsupported under REST2 scaling"
    )
    assert force.getNumExceptionParameterOffsets() == 0, (
        "NonbondedForce has exception parameter offsets — unsupported under REST2 scaling"
    )

    for i in range(force.getNumParticles()):
        if i not in solute:
            continue
        charge, sigma, epsilon = force.getParticleParameters(i)
        force.setParticleParameters(i, charge * sqrt_lam, sigma, epsilon * lam)

    # Exceptions store the combined chargeProd / ε; rescale with the same combining rules.
    # charge factor per atom = √λ (solute) else 1; ε factor per atom = λ (solute) else 1,
    # combined as √(ε_i·ε_j). solute-solute -> (λ, λ); solute-solvent -> (√λ, √λ).
    for k in range(force.getNumExceptions()):
        p1, p2, charge_prod, sigma, epsilon = force.getExceptionParameters(k)
        cf = (sqrt_lam if p1 in solute else 1.0) * (sqrt_lam if p2 in solute else 1.0)
        ef = math.sqrt((lam if p1 in solute else 1.0) * (lam if p2 in solute else 1.0))
        force.setExceptionParameters(k, p1, p2, charge_prod * cf, sigma, epsilon * ef)


def _scale_harmonic_bond(force: HarmonicBondForce, solute, lam: float) -> None:
    for i in range(force.getNumBonds()):
        p1, p2, length, k = force.getBondParameters(i)
        if p1 in solute and p2 in solute:
            force.setBondParameters(i, p1, p2, length, k * lam)


def _scale_harmonic_angle(force: HarmonicAngleForce, solute, lam: float) -> None:
    for i in range(force.getNumAngles()):
        p1, p2, p3, angle, k = force.getAngleParameters(i)
        if p1 in solute and p2 in solute and p3 in solute:
            force.setAngleParameters(i, p1, p2, p3, angle, k * lam)


def _scale_periodic_torsion(force: PeriodicTorsionForce, solute, lam: float) -> None:
    for i in range(force.getNumTorsions()):
        p1, p2, p3, p4, periodicity, phase, k = force.getTorsionParameters(i)
        if p1 in solute and p2 in solute and p3 in solute and p4 in solute:
            force.setTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, k * lam)


def _scale_cmap(force: CMAPTorsionForce, solute, lam: float) -> None:
    """Scale every CMAP energy map by λ (CMAP is a solute backbone term).

    Maps are shared across torsions, so gate first: every CMAP torsion must lie entirely
    within the solute, else scaling a shared map would wrongly temper a non-solute residue
    (e.g. tempering only one chain of a CHARMM complex). Refuse in that case.
    """
    for t in range(force.getNumTorsions()):
        atoms = force.getTorsionParameters(t)[1:]  # (map, a1..a4, b1..b4) -> atom indices
        assert all(a in solute for a in atoms), (
            f"CMAP torsion {t} reaches outside the solute selection — cannot scale the "
            "shared CMAP map without mis-tempering a non-solute residue. Refusing."
        )
    for m in range(force.getNumMaps()):
        size, energy = force.getMapParameters(m)
        force.setMapParameters(m, size, [e * lam for e in energy])
