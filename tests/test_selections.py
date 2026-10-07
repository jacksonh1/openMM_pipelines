"""Atom selections on a hand-built Topology — pure (no force field, no PDB).

Covers the REST2 solute selection `protein_atoms` (whole protein + per-chain subset, H
included) and that non-protein residues are excluded.
"""

import pytest
from openmm.app import Topology, element

from openmm_pipelines.lib.selections import protein_atoms, protein_heavy_atoms


def _two_chain_topology():
    """Chain 0: ALA (CA + HA). Chain 1: GLY (CA). Plus a HOH water (O)."""
    top = Topology()
    c0 = top.addChain()
    r0 = top.addResidue("ALA", c0)
    top.addAtom("CA", element.carbon, r0)
    top.addAtom("HA", element.hydrogen, r0)
    c1 = top.addChain()
    r1 = top.addResidue("GLY", c1)
    top.addAtom("CA", element.carbon, r1)
    cw = top.addChain()
    rw = top.addResidue("HOH", cw)
    top.addAtom("O", element.oxygen, rw)
    return top


def test_protein_atoms_whole_protein_includes_hydrogens_excludes_water():
    top = _two_chain_topology()
    idx = protein_atoms(top)
    assert len(idx) == 3  # ALA CA+HA, GLY CA — water O excluded
    names = {list(top.atoms())[i].name for i in idx}
    assert "HA" in names  # hydrogens kept (REST2 scales them)
    assert "O" not in names  # water excluded


def test_protein_atoms_chain_subset():
    top = _two_chain_topology()
    assert len(protein_atoms(top, chain_index=0)) == 2  # ALA CA + HA
    assert len(protein_atoms(top, chain_index=1)) == 1  # GLY CA


def test_protein_atoms_empty_selection_raises():
    top = _two_chain_topology()
    with pytest.raises(AssertionError):
        protein_atoms(top, chain_index=2)  # the water chain — no protein residues


def test_protein_heavy_atoms_excludes_hydrogen():
    top = _two_chain_topology()
    idx = protein_heavy_atoms(top)
    assert len(idx) == 2  # two CAs; HA (hydrogen) and water O excluded
