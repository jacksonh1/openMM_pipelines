"""PDB -> solvated, parameterized OpenMM System (the first stage of prepare()).

No force introspection — all nonbonded settings flow through `createSystem` kwargs,
so OpenMM applies them correctly whether the force field uses a `NonbondedForce`
(AMBER) or a `CustomNonbondedForce` (CHARMM). Three fail-loud structure-prep guards run
here: no backbone chain break between consecutive residues (a SEQRES-independent catch
for a missing internal loop), disulfide bonds must survive PDBFixer, and solvation must
produce the water model's expected per-residue geometry.
"""

from __future__ import annotations

import io
from typing import TYPE_CHECKING

import numpy as np
from openmm import app, unit
from pdbfixer import PDBFixer

from ..cation_pi import render_nbfix_xml
from ..forcefield import resolve
from ..lib.selections import _PROTEIN_RESIDUES
from .handoff import BuiltSystem

if TYPE_CHECKING:
    from ..config import PrepConfig


class PrepareError(RuntimeError):
    pass


def _count_disulfides(topology) -> int:
    """SG-SG bonds in the topology."""
    return sum(
        1 for a, b in topology.bonds() if a.name == "SG" and b.name == "SG"
    )


def _assert_no_internal_gaps(fixer) -> None:
    """Fail loud if the input has an INTERNAL missing-residue gap (chain break).

    A complete designed structure should have no internal chain breaks. An internal gap
    would make PDBFixer silently rebuild a span AND can make one logical chain read as
    several chains in later multi-chain analysis — so we catch it here, at the boundary,
    rather than guarding for it downstream. Terminal missing residues (at a chain's very
    start/end) are allowed. ASSUMES `fixer.findMissingResidues()` has already run.
    """
    # PDBFixer keys missingResidues by (chain_index, insertion_index), where chain_index
    # indexes list(topology.chains()) and insertion_index is the position within that
    # chain's residues before which residues would be inserted (0 = N-term, len = C-term).
    residues_per_chain = [sum(1 for _ in chain.residues()) for chain in fixer.topology.chains()]
    internal = sorted(
        (chain_index, ins_idx)
        for (chain_index, ins_idx) in fixer.missingResidues
        if 0 < ins_idx < residues_per_chain[chain_index]
    )
    if internal:
        raise PrepareError(
            f"input has internal missing-residue gap(s) at (chain, position) {internal}; "
            "a complete design should have no internal chain breaks. Model the missing "
            "residues in the input structure before preparing."
        )


_PEPTIDE_BOND_MAX_NM = 0.2  # a real C(i)-N(i+1) peptide bond is ~0.133 nm; a consecutive
# protein pair farther apart than this is a backbone break, not a bond.


def _assert_no_backbone_breaks(topology, positions, max_cn_nm: float = _PEPTIDE_BOND_MAX_NM) -> None:
    """Fail loud if consecutive protein residues are too far apart to be peptide-bonded.

    OpenMM bonds consecutive residues within a chain by their order, not their distance,
    so a missing internal loop that carries no SEQRES record — nothing for
    `_assert_no_internal_gaps` to see — yields a spurious long C(i)->N(i+1) bond that
    silently corrupts the system (e.g. a relaxed structure whose SEQRES was stripped).
    This is the SEQRES-independent backstop: it measures the actual C-N distance between
    consecutive protein residues. ASSUMES backbone atoms are present — call after
    `addMissingAtoms`. ACE/NME caps legitimately lack an N or C, so they are skipped on
    the side they lack; a non-protein residue resets the run (a real chain boundary).
    """
    raw = positions.value_in_unit(unit.nanometer)
    breaks = []
    for chain in topology.chains():
        prev = None  # (residue, index of backbone C) of the previous protein residue
        for residue in chain.residues():
            if residue.name not in _PROTEIN_RESIDUES:
                prev = None
                continue
            backbone = {a.name: a.index for a in residue.atoms() if a.name in ("N", "C")}
            if prev is not None and "N" in backbone:
                prev_res, c_idx = prev
                d = float(np.linalg.norm(np.array(raw[c_idx]) - np.array(raw[backbone["N"]])))
                if d > max_cn_nm:
                    breaks.append((chain.id, prev_res.id, residue.id, d))
            prev = (residue, backbone["C"]) if "C" in backbone else None
    if breaks:
        detail = ", ".join(
            f"chain {ch} {r0}->{r1} (C-N {d * 10:.1f} A)" for ch, r0, r1, d in breaks
        )
        raise PrepareError(
            f"backbone chain break(s) between consecutive residues: {detail}. A peptide "
            "bond is ~1.3 A, so a gap this large means missing residues that OpenMM would "
            "bond across, silently corrupting the system. Model the missing loop before "
            "preparing."
        )


def build_system(cfg: PrepConfig) -> BuiltSystem:
    fixer = PDBFixer(filename=str(cfg.pdb_in))
    ss_before = _count_disulfides(fixer.topology)  # capture the input pose's disulfides
    fixer.findMissingResidues()
    _assert_no_internal_gaps(fixer)  # GUARD: no internal chain breaks in a designed input
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(cfg.ph)

    # GUARD: PDBFixer must not add/remove disulfides — critical for
    # disulfide-constrained designs (XBL peptides). A changed count means the input
    # SG-SG bonding was not preserved.
    ss_after = _count_disulfides(fixer.topology)
    if ss_after != ss_before:
        raise PrepareError(
            f"disulfide count changed during structure fixing: {ss_before} -> "
            f"{ss_after}. The input pose's CYS bonding was not preserved."
        )

    # GUARD: no spatial backbone break between consecutive residues. SEQRES-independent,
    # so it catches a missing internal loop that _assert_no_internal_gaps cannot see.
    _assert_no_backbone_breaks(fixer.topology, fixer.positions)

    xmls, packing_model = resolve(cfg.protein_ff, cfg.water)  # FF xmls + packing box
    ff_sources: list = list(xmls)
    if cfg.cation_pi_nbfix:
        # Cation-π NBFIX overrides merge into CHARMM36m's LennardJonesForce, so the fragment
        # is loaded LAST. PrepConfig has already enforced protein_ff == CHARMM36M; render
        # re-checks it (fail loud). One ForceField feeds both addSolvent and createSystem.
        ff_sources.append(io.StringIO(render_nbfix_xml(cfg.protein_ff)))
    ff = app.ForceField(*ff_sources)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        ff,
        model=packing_model,  # packing GEOMETRY (site count), not the param name
        boxShape=cfg.box_shape.value,
        padding=cfg.padding_nm * unit.nanometer,
        ionicStrength=cfg.ionic_strength_molar * unit.molar,
        neutralize=cfg.neutralize,
    )

    spec = cfg.protein_ff.nonbonded
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=spec.cutoff_nm * unit.nanometer,
        switchDistance=(spec.switch_nm * unit.nanometer) if spec.switch_nm else None,
        constraints=app.HBonds,
        hydrogenMass=cfg.hydrogen_mass_amu * unit.amu,
    )

    _assert_water_topology(modeller.topology, system, cfg.water)  # catch box-reuse change
    return BuiltSystem(system, modeller.topology, modeller.positions)


def _assert_water_topology(topology, system, water) -> None:
    """Fail loud if solvation did not produce the water model's expected geometry.

    Guards the site-count / packing contract: a wrong packing box would otherwise
    reach production. Every `HOH` must have the model's particle count, and 4-site
    models must carry virtual sites.
    """
    for res in topology.residues():
        if res.name == "HOH":
            n = sum(1 for _ in res.atoms())
            if n != water.n_sites:
                raise PrepareError(
                    f"{water.name} expects {water.n_sites} particles/water but "
                    f"solvation produced {n}. The packing-box reuse for this model "
                    f"may have changed."
                )
    if water.n_sites >= 4 and not any(
        system.isVirtualSite(i) for i in range(system.getNumParticles())
    ):
        raise PrepareError(
            f"{water.name} is a {water.n_sites}-site model but the System has no "
            f"virtual sites; the M-site was not applied."
        )
