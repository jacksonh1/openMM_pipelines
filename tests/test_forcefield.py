"""resolve() pairing + packing-by-site-count, and the per-FF nonbonded specs."""

import pytest

from openmm_pipelines.forcefield import (
    ForceFieldError,
    ProteinFF,
    WaterModel,
    _PACKING_BOX,
    _RECOMMENDED_WATER,
    resolve,
)


def test_amber_water_in_same_family_tree():
    xmls, packing = resolve(ProteinFF.AMBER14SB, WaterModel.TIP3P)
    assert xmls == ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"]
    assert packing == "tip3p"


def test_amber19_opc_is_four_site_packing():
    xmls, packing = resolve(ProteinFF.AMBER19SB, WaterModel.OPC)
    assert xmls == ["amber19/protein.ff19SB.xml", "amber19/opc.xml"]
    assert packing == "tip4pew"  # 4-site model packs in the tip4pew box


def test_packing_chosen_by_site_count_not_name():
    # every 4-site water packs in the same box regardless of its own name
    for w in (WaterModel.OPC, WaterModel.TIP4PFB, WaterModel.TIP4PEW):
        _, packing = resolve(ProteinFF.AMBER14SB, w)
        assert packing == _PACKING_BOX[4] == "tip4pew"


def test_charmm_uses_bundled_water_xml():
    xmls, packing = resolve(ProteinFF.CHARMM36M, WaterModel.TIP3P)
    assert xmls == ["charmm36_2024.xml", "charmm36_2024/water.xml"]
    assert packing == "tip3p"


def test_charmm_with_amber_water_raises():
    with pytest.raises(ForceFieldError):
        resolve(ProteinFF.CHARMM36, WaterModel.OPC)
    with pytest.raises(ForceFieldError):
        resolve(ProteinFF.CHARMM36M, WaterModel.TIP4PEW)


def test_every_protein_ff_has_family_and_nonbonded():
    # exhaustiveness: no member falls through the match statements
    for ff in ProteinFF:
        assert ff.family
        spec = ff.nonbonded
        assert spec.cutoff_nm > 0


def test_charmm_is_force_switched():
    for ff in (ProteinFF.CHARMM36, ProteinFF.CHARMM36M):
        spec = ff.nonbonded
        assert spec.cutoff_nm == 1.2
        assert spec.switch_nm == 1.0


def test_amber_has_no_switch():
    for ff in (ProteinFF.AMBER14SB, ProteinFF.AMBER19SB, ProteinFF.AMBER99SBILDN):
        assert ff.nonbonded.switch_nm is None


def test_every_ff_has_a_recommended_water():
    for ff in ProteinFF:
        assert ff in _RECOMMENDED_WATER


def test_water_site_counts():
    assert WaterModel.TIP3P.n_sites == 3
    assert WaterModel.OPC.n_sites == 4
    assert WaterModel.OPC.stem == "opc"
