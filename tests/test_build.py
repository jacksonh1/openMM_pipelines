"""build_system end-to-end on a small real structure + the structure-prep guards.

Runs real PDBFixer + solvation + createSystem, so it needs the openmm env. Uses the
small single-chain helix_fusion.pdb to stay fast.
"""

from pathlib import Path

import pytest
from openmm.app import Element, Topology

from openmm_pipelines.preparation.build import _count_disulfides, build_system
from openmm_pipelines.forcefield import ProteinFF, WaterModel
from openmm_pipelines.md.config import MDConfig

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"


@pytest.fixture(scope="module")
def tip3p_build(tmp_path_factory):
    cfg = MDConfig(pdb_in=STRUCT, outdir=tmp_path_factory.mktemp("tip3p"),
                   padding_nm=0.8)  # small box -> faster test
    return build_system(cfg), cfg


def _water_residues(topology):
    return [r for r in topology.residues() if r.name == "HOH"]


def test_build_returns_solvated_system(tip3p_build):
    built, _ = tip3p_build
    waters = _water_residues(built.topology)
    assert waters, "no water added — solvation failed"
    assert built.system.getNumParticles() == built.topology.getNumAtoms()


def test_tip3p_is_three_site(tip3p_build):
    built, _ = tip3p_build
    for r in _water_residues(built.topology):
        assert sum(1 for _ in r.atoms()) == 3


def test_tip3p_has_no_virtual_sites(tip3p_build):
    built, _ = tip3p_build
    sys = built.system
    assert not any(sys.isVirtualSite(i) for i in range(sys.getNumParticles()))


def test_positions_match_particle_count(tip3p_build):
    built, _ = tip3p_build
    assert len(built.positions) == built.system.getNumParticles()


@pytest.mark.slow
def test_opc_is_four_site_with_virtual_sites(tmp_path):
    # ff19SB + OPC: packing-by-site-count must yield real 4-site water with M-sites.
    cfg = MDConfig(pdb_in=STRUCT, outdir=tmp_path / "opc",
                   protein_ff=ProteinFF.AMBER19SB, water=WaterModel.OPC, padding_nm=0.8)
    built = build_system(cfg)
    for r in (r for r in built.topology.residues() if r.name == "HOH"):
        assert sum(1 for _ in r.atoms()) == 4
    sys = built.system
    assert any(sys.isVirtualSite(i) for i in range(sys.getNumParticles()))


# --- disulfide counter (guard input) ------------------------------------------

def _two_cys_topology(bonded: bool):
    top = Topology()
    chain = top.addChain()
    atoms = []
    for _ in range(2):
        res = top.addResidue("CYS", chain)
        atoms.append(top.addAtom("SG", Element.getBySymbol("S"), res))
    if bonded:
        top.addBond(atoms[0], atoms[1])
    return top


def test_count_disulfides_counts_sg_sg_bond():
    assert _count_disulfides(_two_cys_topology(bonded=True)) == 1


def test_count_disulfides_zero_when_unbonded():
    assert _count_disulfides(_two_cys_topology(bonded=False)) == 0
