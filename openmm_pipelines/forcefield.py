"""Force-field and water-model registry, and the single ``resolve()`` that pairs them.

Pure Python — no OpenMM import. Every nonbonded setting is a *property of the force
field* (`NonbondedSpec`), so it flows through `createSystem` kwargs and is never
special-cased by name downstream. Water packing is chosen by *site count*, water
parameters by *force field*, and the FF/water pairing is validated in exactly one
place (`resolve`) so the solute and solvent can never silently disagree.
"""

from __future__ import annotations

from enum import Enum
from typing import NamedTuple


class ForceFieldError(Exception):
    """An illegal force-field / water pairing (a correctness error, not a warning)."""


class NonbondedSpec(NamedTuple):
    cutoff_nm: float
    switch_nm: float | None  # None -> no switching function
    validated: bool  # False -> config refuses it until energy-matched


class ProteinFF(Enum):
    AMBER99SBILDN = "amber99sbildn.xml"  # legacy; for reproducing older work
    AMBER14SB = "amber14/protein.ff14SB.xml"
    AMBER19SB = "amber19/protein.ff19SB.xml"  # QM-trained backbone; recommended with OPC
    CHARMM36 = "charmm36.xml"  # original C36  (par_all36_prot)
    CHARMM36M = "charmm36_2024.xml"  # C36m      (par_all36m_prot)
    # CRITICAL: CHARMM36 != CHARMM36M. 36m is the Huang-2017 backbone-CMAP + glycine
    # refinement (usually the right choice for folding/stability systems such as WW
    # domains). All verified to load in OpenMM 8.6.1; modular paths preferred over the
    # -all bundles, which also pull nucleic-acid params. No openmmforcefields needed.

    @property
    def family(self) -> str:
        """Solvent-parameter tree this FF draws water + ions from."""
        match self:
            case ProteinFF.AMBER99SBILDN | ProteinFF.AMBER14SB:
                return "amber14"
            case ProteinFF.AMBER19SB:
                return "amber19"
            case ProteinFF.CHARMM36:
                return "charmm36"
            case ProteinFF.CHARMM36M:
                return "charmm36_2024"
        # no default: a new member with no family fails loudly here (exhaustiveness)

    @property
    def nonbonded(self) -> NonbondedSpec:
        match self:
            case ProteinFF.AMBER99SBILDN:
                return NonbondedSpec(1.0, None, True)
            case ProteinFF.AMBER14SB:
                return NonbondedSpec(1.0, None, True)
            case ProteinFF.AMBER19SB:
                return NonbondedSpec(1.0, None, True)
            case ProteinFF.CHARMM36:
                return NonbondedSpec(1.2, 1.0, True)  # CHARMM: force-switched vdW
            case ProteinFF.CHARMM36M:
                return NonbondedSpec(1.2, 1.0, True)
        # no default: a new member with no spec fails loudly here (exhaustiveness)


class WaterModel(Enum):
    # (xml stem within the FF family tree, particle count per water)
    TIP3P = ("tip3p", 3)
    SPCE = ("spce", 3)
    TIP3PFB = ("tip3pfb", 3)
    OPC3 = ("opc3", 3)
    TIP4PEW = ("tip4pew", 4)
    TIP4PFB = ("tip4pfb", 4)
    OPC = ("opc", 4)  # 4-site; recommended with ff19SB

    @property
    def stem(self) -> str:
        return self.value[0]

    @property
    def n_sites(self) -> int:
        return self.value[1]


# GENERAL water-packing rule. Modeller.addSolvent ships pre-equilibrated boxes ONLY
# for a few geometries and rejects any other name (model='opc' -> ValueError). The
# packing box only needs to match the water's SITE COUNT; parameters come from the FF
# xml, and the minimization step (already in the pipeline) corrects the residual
# geometry. OpenMM's own documented idiom ("a box of TIP4P-Ew water can be used for
# most four-site water models"). Adding a new non-polarizable water = one WaterModel
# entry with its site count; packing is derived, never special-cased. (Drude /
# polarizable water such as swm4ndp is out of scope — needs a Drude integrator.)
_PACKING_BOX = {3: "tip3p", 4: "tip4pew", 5: "tip5p"}  # site count -> addSolvent model=

_RECOMMENDED_WATER = {
    ProteinFF.AMBER99SBILDN: WaterModel.TIP3P,
    ProteinFF.AMBER14SB: WaterModel.TIP3P,
    ProteinFF.AMBER19SB: WaterModel.OPC,  # ff19SB trained for OPC; TIP3P degrades it
    ProteinFF.CHARMM36: WaterModel.TIP3P,  # CHARMM-modified TIP3P (see resolve)
    ProteinFF.CHARMM36M: WaterModel.TIP3P,
}


def resolve(protein: ProteinFF, water: WaterModel) -> tuple[list[str], str]:
    """``(ForceField xml list, addSolvent packing model)``.

    Raises `ForceFieldError` on an illegal pairing — so the solute FF and the water
    can never disagree — and the packing geometry is chosen by site count in ONE
    place.
    """
    fam = protein.family
    if fam in ("charmm36", "charmm36_2024"):
        # CHARMM water is version-locked to its protein FF; only its bundled
        # modified-TIP3P is valid. Mixing an AMBER-tree water with CHARMM is a
        # correctness error.
        if water is not WaterModel.TIP3P:
            raise ForceFieldError(
                f"{protein.name} must use its bundled CHARMM-modified TIP3P "
                f"({fam}/water.xml); got {water.name}. Do not mix AMBER water params "
                f"with CHARMM."
            )
        return [protein.value, f"{fam}/water.xml"], _PACKING_BOX[3]
    # AMBER: water parameters + ions come from the same family tree as the protein FF.
    return [protein.value, f"{fam}/{water.stem}.xml"], _PACKING_BOX[water.n_sites]
