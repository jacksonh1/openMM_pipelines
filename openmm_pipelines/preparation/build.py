"""PDB -> solvated, parameterized OpenMM System (the first stage of prepare()).

No force introspection — all nonbonded settings flow through `createSystem` kwargs,
so OpenMM applies them correctly whether the force field uses a `NonbondedForce`
(AMBER) or a `CustomNonbondedForce` (CHARMM). Two fail-loud structure-prep guards run
here: disulfide bonds must survive PDBFixer, and solvation must produce the water
model's expected per-residue geometry.
"""

from __future__ import annotations

from openmm import app, unit
from pdbfixer import PDBFixer

from ..forcefield import resolve
from .handoff import BuiltSystem


class PrepareError(RuntimeError):
    pass


def _count_disulfides(topology) -> int:
    """SG-SG bonds in the topology."""
    return sum(
        1 for a, b in topology.bonds() if a.name == "SG" and b.name == "SG"
    )


def build_system(cfg) -> BuiltSystem:  # cfg: PrepConfig
    fixer = PDBFixer(filename=str(cfg.pdb_in))
    ss_before = _count_disulfides(fixer.topology)  # capture the input pose's disulfides
    fixer.findMissingResidues()
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

    xmls, packing_model = resolve(cfg.protein_ff, cfg.water)  # FF xmls + packing box
    ff = app.ForceField(*xmls)
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
