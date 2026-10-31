"""Cation-π NBFIX: the pinned Table-1 subset, unit conversion, and the rendered fragment.

The pure tests (table/conversion/render) need no OpenMM. The two tagged tests load the
fragment into a real CHARMM36m ForceField to prove the pairs actually register.
"""

import io
import math
from pathlib import Path

import pytest

from openmm_pipelines.cation_pi import (
    COMPATIBLE_FF,
    nbfix_pairs,
    render_nbfix_xml,
)
from openmm_pipelines.forcefield import ForceFieldError, ProteinFF, WaterModel, resolve

_INPUTS = Path(__file__).parents[1] / "examples" / "input_structures"
_CATION_PI_RICH = _INPUTS / "1a22-fixed.pdb"   # HGH: many Arg/Lys + Trp/Tyr/Phe (real cation-π)
_NO_CATION_NO_AROMATIC = _INPUTS / "helix_fusion.pdb"  # poly-Ala: neither partner present

_KJ_PER_KCAL = 4.184
_RMIN_OVER_SIGMA = 2.0 ** (1.0 / 6.0)


def test_pairs_are_the_paper_protein_subset():
    classes = {(p.cation_class, p.aromatic_class) for p in nbfix_pairs()}
    assert classes == {
        ("NC2", "CA"), ("NC2", "CAI"), ("NC2", "CPT"),   # Arg guanidinium
        ("NH3", "CA"), ("NH3", "CAI"), ("NH3", "CPT"),   # Lys ammonium
    }
    # Non-protein rows are excluded: CAB (para-cresol toy atom type) and the TMA+ classes.
    aromatic = {p.aromatic_class for p in nbfix_pairs()}
    cation = {p.cation_class for p in nbfix_pairs()}
    assert "CAB" not in aromatic
    assert cation == {"NC2", "NH3"}  # no NG3P0 / CG334 (TMA+)


def test_pairs_carry_paper_values_verbatim():
    # Spot-check two rows against Table 1 (the pinned source must stay in paper units).
    nc2_ca = next(p for p in nbfix_pairs() if p[:2] == ("NC2", "CA"))
    assert (nc2_ca.rmin_angstrom, nc2_ca.epsilon_kcal_per_mol) == (3.90, 0.42500)
    nh3_cpt = next(p for p in nbfix_pairs() if p[:2] == ("NH3", "CPT"))
    assert (nh3_cpt.rmin_angstrom, nh3_cpt.epsilon_kcal_per_mol) == (3.20, 0.80000)


def test_unit_conversion_to_openmm():
    p = next(p for p in nbfix_pairs() if p[:2] == ("NC2", "CA"))  # 3.90 Å, 0.425 kcal/mol
    assert p.epsilon_kj_per_mol == pytest.approx(0.425 * _KJ_PER_KCAL)
    assert p.sigma_nm == pytest.approx(3.90 / _RMIN_OVER_SIGMA / 10.0)
    assert 0.33 < p.sigma_nm < 0.37  # physical range for an N–C(aromatic) pair


def test_render_fragment_structure():
    xml = render_nbfix_xml(ProteinFF.CHARMM36M)
    assert xml.count("<NBFixPair") == 6
    assert 'lj14scale="1.0"' in xml  # must match charmm36_2024.xml's LennardJonesForce header
    assert 'class1="NC2" class2="CA"' in xml
    assert 'class1="NH3" class2="CPT"' in xml


def test_render_refuses_non_charmm36m():
    for ff in (ProteinFF.AMBER19SB, ProteinFF.CHARMM36):  # C36 ≠ C36m; paper is C36m-specific
        with pytest.raises(ForceFieldError):
            render_nbfix_xml(ff)
    assert COMPATIBLE_FF is ProteinFF.CHARMM36M


def _lennard_jones_generator(ff):
    gens = [g for g in ff._forces if type(g).__name__ == "LennardJonesGenerator"]
    assert len(gens) == 1, "CHARMM force field must have exactly one LennardJonesGenerator"
    return gens[0]


def _nbfix_table_size(ff):
    g = _lennard_jones_generator(ff)
    return len(g.nbfixTypes1) + len(g.nbfixTypes2)


def test_fragment_registers_extra_nbfix_pairs():
    """Loading the fragment onto CHARMM36m must enlarge the generator's NBFIX pair table."""
    import openmm.app as app

    xmls, _ = resolve(ProteinFF.CHARMM36M, WaterModel.TIP3P)
    base = app.ForceField(*xmls)
    with_catpi = app.ForceField(*xmls, io.StringIO(render_nbfix_xml(ProteinFF.CHARMM36M)))
    assert _nbfix_table_size(with_catpi) > _nbfix_table_size(base)


def _vacuum_energy_kj(pdb_path: Path, *, with_cation_pi: bool) -> float:
    """Potential energy (kJ/mol) of a fixed+protonated structure built in vacuum (NoCutoff).

    Mirrors `build_system`'s structure prep (PDBFixer → add atoms/H) but skips solvation
    and PME: vacuum + NoCutoff so the base and with-fragment builds are byte-identical
    EXCEPT the overridden cation↔aromatic LJ pairs. The energy DIFFERENCE therefore equals
    the override's effect on real atoms, which `createSystem` either wired in or did not.
    """
    import openmm
    import openmm.app as app
    from openmm import unit
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    xmls, _ = resolve(ProteinFF.CHARMM36M, WaterModel.TIP3P)
    sources: list = list(xmls)
    if with_cation_pi:
        sources.append(io.StringIO(render_nbfix_xml(ProteinFF.CHARMM36M)))
    ff = app.ForceField(*sources)
    system = ff.createSystem(fixer.topology, nonbondedMethod=app.NoCutoff, constraints=None)

    context = openmm.Context(
        system,
        openmm.VerletIntegrator(1.0 * unit.femtosecond),
        openmm.Platform.getPlatformByName("CPU"),
    )
    context.setPositions(fixer.positions)
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )
    del context
    return energy


def test_override_changes_energy_of_cation_pi_rich_protein():
    """Positive control: on a protein with real Arg/Lys↔aromatic contacts, enabling the
    fragment must change the vacuum energy — proof `createSystem` applied the overridden LJ
    to real atoms (the end-to-end step the generator-table test cannot see)."""
    e_base = _vacuum_energy_kj(_CATION_PI_RICH, with_cation_pi=False)
    e_catpi = _vacuum_energy_kj(_CATION_PI_RICH, with_cation_pi=True)
    assert math.isfinite(e_base) and math.isfinite(e_catpi)
    delta = e_catpi - e_base
    # Sign is not asserted: Δ is a sum over many pairs at many distances (the override both
    # deepens ε and shifts r_min outward). What must hold is that it is real and bounded —
    # the override reached the force (|Δ| well above numerical noise) without a blow-up.
    assert 1.0 < abs(delta) < 1.0e4, f"ΔE = {delta} kJ/mol is noise or a blow-up"


def test_override_is_inert_without_cation_or_aromatic():
    """Negative control: on poly-alanine (no NC2/NH3 cation, no aromatic C), the fragment
    must change nothing — the override touches ONLY cation↔aromatic pairs, never global LJ."""
    e_base = _vacuum_energy_kj(_NO_CATION_NO_AROMATIC, with_cation_pi=False)
    e_catpi = _vacuum_energy_kj(_NO_CATION_NO_AROMATIC, with_cation_pi=True)
    assert e_catpi == pytest.approx(e_base)
