"""Cation-π NBFIX correction for CHARMM36m (Liu, Fu, Shao, Cai & Chipot,
J. Chem. Theory Comput. 2020, 16, 6397 — https://doi.org/10.1021/acs.jctc.0c00637).

Pairwise-additive force fields systematically under-bind cation-π contacts because they
average out the cation-induced polarization of the aromatic ring. The paper's remedy is a
handful of *pair-specific* Lennard-Jones overrides (CHARMM's NBFIX feature) that deepen the
LJ well for cation-N / aromatic-C atom-CLASS pairs, fit to reproduce QM/MM potentials of mean
force. It changes no functional form and costs nothing at runtime — OpenMM already routes all
CHARMM LJ (these pairs included) through a ``CustomNonbondedForce``; we only add rows to the
pair table it builds. This module is the pinned source of truth for those rows; `build_system`
renders them to a ForceField fragment and loads it after the base CHARMM36m XML.

We do NOT reimplement the paper's QM/MM + meta-eABF parameterization — we apply its *output*
(Table 1), which is the whole point of the method.

**Scope — the protein-relevant rows of Table 1.** Two basic side chains carry the cation:

- **Arg** guanidinium N → atom class ``NC2``
- **Lys** ammonium N   → atom class ``NH3``  (also the paper's free-NH4+ model)

against the aromatic carbons of **Phe/Tyr** (class ``CA``) and **Trp** (classes ``CA``,
``CAI``, ``CPT``). Phe and Tyr are covered: their ring carbons are class ``CA`` in CHARMM36m
(verified against the shipped residue templates), the same class the paper fit to toluene, so
they ride on the ``*-CA`` rows. Deliberately NOT included:

- **His** — the paper never parameterized imidazole; no cation-π parameters exist for it. It
  stays at stock CHARMM36m (the paper's omission, not a scoping choice here).
- **CAB** (para-cresol) — a *new atom type the paper introduced only for its para-cresol free
  molecule* so that toy's PMF could be fit independently of toluene. In standard CHARMM36 both
  toluene and para-cresol carbons are class ``CA`` (supp Table S2), and real Tyr is ``CA`` too,
  so ``CAB`` is never present in a protein. The ``NC2/NH3-CAB`` rows are toy-model-only.
- **TMA+** (tetramethylammonium; classes ``NG3P0``/``CG334``) — a *quaternary* ammonium. No
  standard amino acid has one; it appears only in ligands (choline, acetylcholine) or modified
  residues (trimethyl-Lys), which are typed by CGenFF, not the protein FF.

**NBFIX is global by atom CLASS, not per residue pair.** Enabling this deepens the well for
*every* Arg/Lys ↔ aromatic contact in the system (Arg–Phe and Arg–Tyr as well as the Arg–Trp
of interest), because these classes are shared across those residues. That is the paper's
intent — a uniform FF correction validated on whole proteins — not a single-contact patch, and
there is no way to restrict an NBFIX to one residue pair.

**CHARMM36m only.** The values are atom-class LJ parameters specific to the CHARMM36m protein
force field and are meaningless against AMBER's atom typing. `COMPATIBLE_FF` pins this and
`PrepConfig` refuses any other FF. REST2 is additionally incompatible: the correction lives in
CHARMM's ``CustomNonbondedForce``, whose solute-scaling is unimplemented — `REST2Config` refuses
it (and `rest2.scaling.scale_solute` refuses CHARMM as a backstop). Use plain MD or T-REMD.

**Extending to non-natural / CGenFF-typed residues.** ncAAs (e.g. trimethyl-Lys, naphthyl-Ala)
are typed by CGenFF, which names the same groups differently (guanidinium ``NG2P1``, TMA+
``NG3P0``/``CG334``, benzene-like aromatic ``CG2R61``, fused-ring ``CG2RC0``). The paper's supp
Table S1 is the CGenFF twin of Table 1 for exactly this case. Adding support is literally
appending rows below (the render/merge mechanism is class-name-agnostic) — not done now because
no such residue is in use, and note the paper leaves genuine gaps (mixed natural-cation ×
CGenFF-aromatic pairs, and partially methylated Lys, are tabulated nowhere).
"""

from __future__ import annotations

from typing import NamedTuple

from .forcefield import ForceFieldError, ProteinFF

# The paper optimized CHARMM36m specifically; its NC2/NH3/CA/CAI/CPT class LJ overrides do not
# transfer to another force field's atom typing.
COMPATIBLE_FF = ProteinFF.CHARMM36M

# Conversions from the paper's CHARMM-convention units (Table 1) to OpenMM's NBFixPair units.
_KJ_PER_KCAL = 4.184  # thermochemical calorie (the CHARMM/OpenMM convention)
_RMIN_OVER_SIGMA = 2.0 ** (1.0 / 6.0)  # pair r_min = 2^(1/6) · sigma (the LJ minimum)
_ANGSTROM_PER_NM = 10.0


class CationPiPair(NamedTuple):
    """One pair-specific LJ override, in the paper's own units (Table 1).

    `cation_class`/`aromatic_class` are CHARMM36m atom *classes*. `rmin_angstrom` is the pair
    minimum-energy distance and `epsilon_kcal_per_mol` the well-depth MAGNITUDE, both exactly
    as tabulated (the paper writes epsilon negative; we store the magnitude). `cation` names the
    residue the cation class belongs to, for provenance/readability only.
    """

    cation_class: str
    aromatic_class: str
    rmin_angstrom: float
    epsilon_kcal_per_mol: float
    cation: str

    @property
    def sigma_nm(self) -> float:
        return self.rmin_angstrom / _RMIN_OVER_SIGMA / _ANGSTROM_PER_NM

    @property
    def epsilon_kj_per_mol(self) -> float:
        return self.epsilon_kcal_per_mol * _KJ_PER_KCAL


# Paper Table 1, protein-FF rows only (CAB and TMA+ rows excluded — see module docstring).
# These are the OPTIMIZED NBFIX parameters; the standard-FF originals they replace are supp
# Table S2. Arg = guanidinium (NC2); Lys = ammonium (NH3).
_CATION_PI_NBFIX: tuple[CationPiPair, ...] = (
    # Arg guanidinium N (class NC2) ↔ aromatic C
    CationPiPair("NC2", "CA", 3.90, 0.42500, "ARG"),
    CationPiPair("NC2", "CAI", 4.00, 0.50000, "ARG"),
    CationPiPair("NC2", "CPT", 4.00, 0.54000, "ARG"),
    # Lys ammonium N (class NH3) ↔ aromatic C
    CationPiPair("NH3", "CA", 3.30, 0.70000, "LYS"),
    CationPiPair("NH3", "CAI", 3.20, 0.70000, "LYS"),
    CationPiPair("NH3", "CPT", 3.20, 0.80000, "LYS"),
)


def nbfix_pairs() -> tuple[CationPiPair, ...]:
    """The pinned cation-π NBFIX pairs (paper Table 1, protein subset)."""
    return _CATION_PI_NBFIX


def render_nbfix_xml(protein_ff: ProteinFF) -> str:
    """An OpenMM ForceField fragment (XML string) adding the cation-π NBFIX overrides.

    Loaded *after* the base CHARMM36m XML, its ``<NBFixPair>`` rows merge into the
    ``<LennardJonesForce>`` that force field already defines (``lj14scale`` must match, hence
    "1.0"). Refuses any force field other than `COMPATIBLE_FF` — defense in depth; `PrepConfig`
    enforces the same pairing, so a mismatch here is a programming error, not user input.
    """
    if protein_ff is not COMPATIBLE_FF:
        raise ForceFieldError(
            f"cation-π NBFIX is parameterized for {COMPATIBLE_FF.name} only (Liu et al. 2020); "
            f"its atom-class LJ overrides have no equivalent under {protein_ff.name}."
        )
    lines = [
        "<ForceField>",
        " <!-- Cation-π NBFIX overrides for CHARMM36m — Liu, Fu, Shao, Cai & Chipot,",
        "      J. Chem. Theory Comput. 2020, 16, 6397, Table 1 (protein subset).",
        "      Generated from openmm_pipelines.cation_pi; do not edit by hand. -->",
        ' <LennardJonesForce lj14scale="1.0">',
    ]
    for p in _CATION_PI_NBFIX:
        lines.append(
            f"  <!-- {p.cation} {p.cation_class}-{p.aromatic_class}: "
            f"rmin={p.rmin_angstrom} Å, eps={p.epsilon_kcal_per_mol} kcal/mol -->"
        )
        lines.append(
            f'  <NBFixPair class1="{p.cation_class}" class2="{p.aromatic_class}" '
            f'sigma="{p.sigma_nm}" epsilon="{p.epsilon_kj_per_mol}"/>'
        )
    lines += [" </LennardJonesForce>", "</ForceField>", ""]
    return "\n".join(lines)
