"""Topology atom selections.

Operates on an OpenMM `Topology`. `protein_heavy_atoms` is the restraint selection:
protein heavy atoms (no hydrogens, no solvent/ions), used to hold the designed pose
through equilibration.
"""

from __future__ import annotations

# Standard amino acids + common protonation/termini variants (AMBER + CHARMM naming).
_PROTEIN_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        # protonation / disulfide / terminal variants
        "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH",
        "LYN", "ACE", "NME", "NMA",
    }
)


def protein_heavy_atoms(topology) -> list[int]:
    """Indices of protein heavy atoms (element present and not hydrogen).

    Fails loud on an empty selection — restraining nothing is a bug (wrong residue
    naming, or a non-protein input).
    """
    indices = [
        atom.index
        for atom in topology.atoms()
        if atom.residue.name in _PROTEIN_RESIDUES
        and atom.element is not None
        and atom.element.symbol != "H"
    ]
    assert indices, (
        "protein_heavy_atoms selected no atoms — no standard protein residues found "
        "in the topology (check residue naming)"
    )
    return indices
