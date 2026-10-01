# openmm_pipelines — engine skeleton (design sketch)

Status: **design sketch, not yet implemented.** This file records the intended shape of
the engine so implementation can start from an agreed plan. No `.py` files exist yet.

## What this is

A small, Pythonic openMM engine that ports the **general science flow** of the existing
GROMACS plain-MD pipeline (`RELE_simulations/gromacs_REMD/scripts/simulation/MD-gromacs.sbatch`),
using **openMM conventions** rather than mirroring the GROMACS engine. It is a `tools/`
library — same tier as `snekwrap` and `FragForge` — consumed by campaigns and projects
(e.g. `Ellen_WW_simulations`, XBL bound-state sampling) via `import openmm_md`.

Scope now: **plain production MD.** But the whole design is organized around one seam so
that **T-REMD and REST2 can be added later without reworking any of the shared code.**
See "Architecture: the prepare/produce seam" and "Extending to REMD/REST2" below.
(`tools/openmm_enh_samp/` is an older stub for the broader scope and is expected to be
retired; this package is the foundation those engines would import.)

Design rules inherited from the workspace: **fail loudly** (crash on violated assumptions,
no silent fallbacks), **make illegal states unrepresentable**, **simple single-purpose
primitives composed into behavior**.

## Architecture: the prepare/produce seam

Everything up to production is **identical across MD, T-REMD and REST2**: build the solvated
system, restrain, minimize, NVT-equilibrate, density-equilibrate at the reference (lowest)
temperature, release restraints. Only **production** differs:

- **plain MD** — one trajectory via `app.Simulation.step()`.
- **T-REMD / REST2** — many replicas via an `openmmtools.multistate` sampler, which drives its
  own `Context`s from `ThermodynamicState` + `SamplerState` objects. It does **not** use
  `app.Simulation`.

So the engine is split at that seam:

```
prepare(cfg) ---------------------------> EquilibratedSystem      (SHARED — ~90% of the code)
                                          (System, topology, positions, box_vectors)
                                                   |
                    +------------------------------+------------------------------+
                    v                              v                              v
           produce_md(eq, cfg)          produce_remd(eq, cfg)          produce_rest2(eq, cfg)
           app.Simulation.step()        multistate.ReplicaExchange     multistate + REST region
           (this package, now)          (later; imports prepare)       (later; imports prepare)
```

`prepare()` returns **plain openMM objects** (a serialized `System` + a `State` with positions
and box vectors). That is exactly what both `app.Simulation` and `openmmtools.multistate`
consume, so the production backends share nothing but this handoff. `prepare()` also writes the
equilibrated system to disk, so equilibration can be run once and several production backends
(or a separate SLURM job / separate package) can launch from the same equilibrated box.

## Science flow (ported) and openMM mapping

| GROMACS stage | openMM equivalent | prepare / produce |
|---|---|---|
| build: pdb2gmx → editconf → solvate → genion | `pdbfixer` + `ForceField` + `Modeller.addSolvent` | prepare |
| EM (steepest descent, restrained) | `simulation.minimizeEnergy()` with restraints | prepare |
| heat (NVT, restrained) | `setVelocitiesToTemperature` + NVT segment | prepare |
| density (NPT, restrained, iterative plateau) | `MonteCarloBarostat` + segment loop | prepare |
| relax (unrestrained NPT, optional) | set restraint `k=0`, keep running NPT | prepare |
| production (NPT/NVT) | MD: `Simulation.step()`; REMD/REST2: `multistate` sampler | **produce** |
| export PDB + analysis | `PDBFile` / mdtraj PBC-wrap + a new in-package mdtraj analysis module | produce |

## Dependencies

Core openMM covers every protein force field needed — **no `openmmforcefields` dependency**
(verified against openMM 8.6.1 in `openmm_env`). `openmmtools` is required only for the future
REMD/REST2 backends, so it is an **optional extra**, not a core dependency of plain MD.

```
# plain MD (now):
conda create -n openmm_md -c conda-forge openmm pdbfixer mdtraj
# + REMD/REST2 (later):
conda install -c conda-forge openmmtools
```

## Module layout

```
tools/openmm_pipelines/
  pyproject.toml         # extras: [remd] -> openmmtools
  CLAUDE.md
  README.md
  openmm_md/
    __init__.py          # exports PrepConfig, MDConfig, ProteinFF, WaterModel, prepare, run_md
    config.py            # PrepConfig (shared base) + MDConfig; enums; derived step counts
    build.py             # pdb -> solvated System + Topology + positions
    restraints.py        # harmonic position restraint, k as global param
    density.py           # volume-plateau convergence (port density_converged.py verbatim)
    stages.py            # minimize / equilibrate_nvt / equilibrate_density / relax primitives
    prepare.py           # prepare(cfg: PrepConfig) -> EquilibratedSystem  (SHARED pipeline)
    produce.py           # produce_md(eq, cfg) + run_md(cfg)  (plain-MD backend)
    equilibrated.py      # EquilibratedSystem dataclass + to_disk / from_disk
    analysis/            # new mdtraj-based analysis subpackage (NOT gromd_analysis)
      __init__.py        #   re-exports the analysis entry points
      trajectory.py      #   load xtc+topology, PBC image_molecules, strip solvent, align to input pose
      metrics.py         #   rmsd_to_reference / radius_of_gyration / rmsf_per_residue (vs input pose)
      dssp.py            #   secondary structure via mdtraj.compute_dssp
      clustering.py      #   conformational clustering (scipy)
      plots.py           #   RMSD / Rg / RMSF / DSSP figures
      # REMD-specific metrics (exchange acceptance, round-trip mixing) land here with that backend
    # future (here or in openmm_enh_samp), all importing prepare():
    #   remd.py          # REMDConfig + run_remd  (openmmtools multistate)
    #   rest2.py         # REST2Config + run_rest2
    #   ladder.py        # geometric temperature ladder, acceptance helpers
  tests/
    test_config.py
    test_density.py       # port gromacs density_converged tests verbatim
```

Package import name is `openmm_md` (underscore) — deliberately **not** `openmm`, to avoid
shadowing the real package.

## config.py

Config is split by the same seam: a shared `PrepConfig` base (everything through
equilibration) and a thin `MDConfig` subclass adding production fields. Future `REMDConfig` /
`REST2Config` subclass the **same** `PrepConfig`. Inheritance (not nesting) keeps attribute
access flat — `cfg.temperature_k` works on every engine's config — while `prepare()` only ever
depends on `PrepConfig` fields.

```python
from enum import Enum
from pathlib import Path
from typing import NamedTuple
from pydantic import BaseModel, model_validator


class BoxShape(str, Enum):
    CUBE = "cube"
    DODECAHEDRON = "dodecahedron"
    OCTAHEDRON = "octahedron"


class NonbondedSpec(NamedTuple):
    cutoff_nm: float
    switch_nm: float | None          # None -> no switching function
    validated: bool                  # False -> gate behind an energy-match check


class ForceFieldError(Exception):
    pass


class ProteinFF(Enum):
    AMBER99SBILDN = "amber99sbildn.xml"              # == gromacs default; legacy/reproduction
    AMBER14SB     = "amber14/protein.ff14SB.xml"
    AMBER19SB     = "amber19/protein.ff19SB.xml"     # QM-trained backbone; recommended with OPC
    CHARMM36      = "charmm36.xml"                    # ORIGINAL C36 (par_all36_prot)
    CHARMM36M     = "charmm36_2024.xml"              # 36m (par_all36m_prot)
    # CRITICAL: CHARMM36 != CHARMM36M. 36m is the Huang-2017 backbone-CMAP + glycine
    # refinement (usually the right choice for folding/stability systems like WW domains).
    # All verified to load in openMM 8.6.1 (modular paths preferred over the -all bundles,
    # which also pull nucleic-acid params). No openmmforcefields needed.

    @property
    def family(self) -> str:
        """Solvent-parameter tree this FF draws water + ions from."""
        match self:
            case ProteinFF.AMBER99SBILDN | ProteinFF.AMBER14SB: return "amber14"
            case ProteinFF.AMBER19SB:                           return "amber19"
            case ProteinFF.CHARMM36:                            return "charmm36"
            case ProteinFF.CHARMM36M:                           return "charmm36_2024"

    @property
    def nonbonded(self) -> NonbondedSpec:
        match self:
            case ProteinFF.AMBER99SBILDN: return NonbondedSpec(1.0, None, True)
            case ProteinFF.AMBER14SB:     return NonbondedSpec(1.0, None, True)
            case ProteinFF.AMBER19SB:     return NonbondedSpec(1.0, None, True)
            case ProteinFF.CHARMM36:      return NonbondedSpec(1.2, 1.0, True)
            case ProteinFF.CHARMM36M:     return NonbondedSpec(1.2, 1.0, True)
        # no default: a new member with no spec fails loudly here (exhaustiveness)


class WaterModel(Enum):
    # (xml stem within the FF family tree, particle count per water)
    TIP3P   = ("tip3p",   3)
    SPCE    = ("spce",    3)
    TIP3PFB = ("tip3pfb", 3)
    OPC3    = ("opc3",    3)
    TIP4PEW = ("tip4pew", 4)
    TIP4PFB = ("tip4pfb", 4)
    OPC     = ("opc",     4)          # 4-site; recommended with ff19SB

    @property
    def stem(self) -> str:   return self.value[0]
    @property
    def n_sites(self) -> int: return self.value[1]


# GENERAL water-packing rule. Modeller.addSolvent ships pre-equilibrated boxes ONLY for a few
# geometries and rejects any other name (model='opc' -> ValueError). The packing box only needs
# to match the water's SITE COUNT; parameters come from the FF xml, and the minimization step
# (already in the pipeline) corrects the residual geometry. This is OpenMM's own documented
# idiom ("a box of TIP4P-Ew water can be used for most four-site water models"). Verified in
# openmm_env: model='tip4pew' + opc.xml yields real 4-site OPC (charges O=0, H=+0.6791,
# M=-1.3583 — NOT TIP4P-Ew's +0.5242). Adding a new non-polarizable water = one WaterModel
# entry with its site count; packing is derived, never special-cased. (Drude/polarizable water
# such as swm4ndp is out of scope — needs a Drude integrator.)
_PACKING_BOX = {3: "tip3p", 4: "tip4pew", 5: "tip5p"}      # site count -> addSolvent model=

_RECOMMENDED_WATER = {
    ProteinFF.AMBER99SBILDN: WaterModel.TIP3P,
    ProteinFF.AMBER14SB:     WaterModel.TIP3P,
    ProteinFF.AMBER19SB:     WaterModel.OPC,               # ff19SB trained for OPC; TIP3P degrades it
    ProteinFF.CHARMM36:      WaterModel.TIP3P,             # CHARMM-modified TIP3P (see resolve)
    ProteinFF.CHARMM36M:     WaterModel.TIP3P,
}


def resolve(protein: ProteinFF, water: WaterModel) -> tuple[list[str], str]:
    """(ForceField xml list, addSolvent packing model). Raises on an illegal pairing — so the
    solute FF and the water can never disagree, and the packing geometry is chosen by site count
    in ONE place."""
    fam = protein.family
    if fam in ("charmm36", "charmm36_2024"):
        # CHARMM water is version-locked to its protein FF; only its bundled modified-TIP3P is
        # valid. Mixing an AMBER-tree water with CHARMM is a correctness error.
        if water is not WaterModel.TIP3P:
            raise ForceFieldError(
                f"{protein.name} must use its bundled CHARMM-modified TIP3P ({fam}/water.xml); "
                f"got {water.name}. Do not mix AMBER water params with CHARMM.")
        return [protein.value, f"{fam}/water.xml"], _PACKING_BOX[3]
    # AMBER: water parameters + ions come from the same family tree as the protein FF.
    return [protein.value, f"{fam}/{water.stem}.xml"], _PACKING_BOX[water.n_sites]


class PrepConfig(BaseModel, frozen=True):
    """Everything through equilibration — shared by MD, REMD and REST2."""
    # --- system ---
    pdb_in: Path
    outdir: Path
    protein_ff: ProteinFF = ProteinFF.AMBER14SB
    water: WaterModel = WaterModel.TIP3P
    ph: float = 7.4                         # physiological (blood); PDBFixer protonation heuristic
    box_shape: BoxShape = BoxShape.DODECAHEDRON
    padding_nm: float = 1.0                 # openMM Modeller term (was gromacs box_buffer)
    ionic_strength_molar: float = 0.15      # openMM addSolvent term (was gromacs salt_molar)
    neutralize: bool = True

    # --- dynamics ---
    temperature_k: float = 300.0            # reference/equilibration temperature (= T_min for REMD)
    dt_ps: float = 0.004                    # 4 fs -> requires HMR (default on; see below)
    hydrogen_mass_amu: float = 4.0          # HMR; set 1.008 + dt_ps=0.002 to disable
    friction_ps: float = 1.0               # Langevin friction (openMM convention; was gromacs gamma_ln 2.0)
    ref_p_bar: float = 1.0
    restraint_k: float = 1000.0            # kJ/mol/nm^2 on protein heavy atoms

    # --- equilibration lengths ---
    nvt_equil_ns: float = 0.2              # NVT equilibration (was gromacs "heat")
    relax_ns: float = 0.0                  # optional unrestrained NPT before production

    # --- density equilibration (ported adaptive protocol; see audit note) ---
    density_seg_steps: int = 10_000
    density_min_seg: int = 8
    density_max_seg: int = 20
    density_tol_rel: float = 0.005

    # --- runtime ---
    platform: str = "CUDA"
    seed: int | None = None               # None = random (pythonic; was gromacs -1 sentinel)

    @model_validator(mode="after")
    def _check(self) -> "PrepConfig":
        assert self.pdb_in.exists(), f"pdb_in not found: {self.pdb_in}"
        assert self.density_min_seg <= self.density_max_seg, "density_min_seg > density_max_seg"
        resolve(self.protein_ff, self.water)     # raises ForceFieldError on an illegal pairing
        assert self.protein_ff.nonbonded.validated, (
            f"{self.protein_ff.name} nonbonded spec is unvalidated — "
            "run an energy-match check before production"
        )
        if self.water is not _RECOMMENDED_WATER[self.protein_ff]:
            import warnings
            warnings.warn(
                f"{self.protein_ff.name} is recommended with "
                f"{_RECOMMENDED_WATER[self.protein_ff].name}; you chose {self.water.name}. "
                "Proceeding — verify this is intentional (e.g. ff19SB+TIP3P is a known-degraded "
                "pairing)."
            )
        if self.dt_ps > 0.0025:            # HMR sanity: 4 fs without repartitioned H is unstable
            assert self.hydrogen_mass_amu >= 3.0, (
                f"dt_ps={self.dt_ps} needs HMR (hydrogen_mass_amu >= 3); "
                f"got {self.hydrogen_mass_amu}"
            )
        return self

    # Derived equilibration step counts.
    @property
    def nvt_equil_steps(self) -> int:  return round(self.nvt_equil_ns * 1000 / self.dt_ps)
    @property
    def relax_steps(self) -> int:      return round(self.relax_ns * 1000 / self.dt_ps)


class MDConfig(PrepConfig):
    """Plain production MD."""
    total_ns: float = 100.0
    traj_ps: float = 10.0

    @property
    def prod_steps(self) -> int:       return round(self.total_ns * 1000 / self.dt_ps)
    @property
    def traj_steps(self) -> int:       return round(self.traj_ps / self.dt_ps)
    @property
    def checkpoint_steps(self) -> int: return max(self.traj_steps * 10, 5000)


# --- future, sketched here so the shared base is designed for them ---
# class REMDConfig(PrepConfig):
#     t_max_k: float = 400.0            # temperature_k is T_min
#     n_replicas: int = 24
#     temps_list: tuple[float, ...] | None = None   # explicit ladder overrides geometric
#     exchange_ps: float = 1.0
#     ensemble: Ensemble = Ensemble.NVT             # production ensemble (NVT default)
#     iterations: int = ...            # from total_ns / exchange_ps
#
# class REST2Config(PrepConfig):
#     t_max_eff_k: float = 450.0       # max effective SOLUTE temperature
#     n_replicas: int = 16
#     rest_region: str = "protein"     # selection scaled by lambda
#     exchange_ps: float = 1.0
```

The `validated` flag on `NonbondedSpec` is the structural hook for the CHARMM caution: flip it
to `False` for any force field whose nonbonded handling has not been energy-matched, and config
validation refuses it until checked.

## equilibrated.py — the handoff object

```python
from dataclasses import dataclass
from openmm import System, XmlSerializer, unit
from openmm.app import Topology


@dataclass(frozen=True)
class EquilibratedSystem:
    """Output of prepare(): plain openMM objects that BOTH app.Simulation and
    openmmtools.multistate can consume. Serializable, so equilibration runs once and any
    production backend (MD / REMD / REST2, same or separate job) starts from it."""
    system: System
    topology: Topology
    positions: object            # list[Vec3] with units
    box_vectors: object          # 3x3 with units

    def to_disk(self, outdir):
        (outdir / "system.xml").write_text(XmlSerializer.serialize(self.system))
        # positions + box go in a State xml (write via a throwaway Context or Simulation)
        ...

    @classmethod
    def from_disk(cls, outdir, topology):
        system = XmlSerializer.deserialize((outdir / "system.xml").read_text())
        ...
        return cls(system, topology, positions, box_vectors)
```

## build.py

No force introspection — all nonbonded settings flow through `createSystem` kwargs, so openMM
applies them correctly whether the force field uses a `NonbondedForce` (AMBER) or a
`CustomNonbondedForce` (CHARMM).

```python
from openmm import app, unit
from pdbfixer import PDBFixer
from .config import resolve
from .equilibrated import BuiltSystem     # (system, topology, positions) frozen dataclass


class PrepareError(RuntimeError):
    pass


def _count_disulfides(topology) -> int:
    """SG-SG bonds in the topology."""
    return sum(1 for a, b in topology.bonds()
               if a.name == "SG" and b.name == "SG")


def build_system(cfg) -> BuiltSystem:     # cfg: PrepConfig
    fixer = PDBFixer(filename=str(cfg.pdb_in))
    ss_before = _count_disulfides(fixer.topology)     # capture the input pose's disulfides
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(cfg.ph)

    # GUARD: PDBFixer must not add/remove disulfides — critical for disulfide-constrained
    # designs (XBL peptides). A changed count means the input SG-SG bonding was not preserved.
    ss_after = _count_disulfides(fixer.topology)
    if ss_after != ss_before:
        raise PrepareError(
            f"disulfide count changed during structure fixing: {ss_before} -> {ss_after}. "
            "The input pose's CYS bonding was not preserved.")

    xmls, packing_model = resolve(cfg.protein_ff, cfg.water)   # one place: FF xmls + packing box
    ff = app.ForceField(*xmls)
    modeller = app.Modeller(fixer.topology, fixer.positions)
    modeller.addSolvent(
        ff,
        model=packing_model,                    # packing GEOMETRY (site count), not the param name
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

    _assert_water_topology(modeller.topology, system, cfg.water)   # catch a future box-reuse change
    return BuiltSystem(system, modeller.topology, modeller.positions)


def _assert_water_topology(topology, system, water) -> None:
    """Fail loud if solvation did not produce the water model's expected geometry. Guards the
    site-count / packing contract: a wrong packing box would otherwise reach production."""
    for res in topology.residues():
        if res.name == "HOH":
            n = sum(1 for _ in res.atoms())
            if n != water.n_sites:
                raise PrepareError(
                    f"{water.name} expects {water.n_sites} particles/water but solvation "
                    f"produced {n}. The packing-box reuse for this model may have changed.")
    if water.n_sites >= 4 and not any(
            system.isVirtualSite(i) for i in range(system.getNumParticles())):
        raise PrepareError(
            f"{water.name} is a {water.n_sites}-site model but the System has no virtual sites; "
            "the M-site was not applied.")
```

**Verified end to end in `openmm_env`** (1a22, ff19SB + OPC): `model='opc'` raises
`ValueError: Unknown water model: opc`; `model='tip4pew'` + `opc.xml` yields 4-site waters
`[O, H1, H2, M]` with charges `O=0.0, H=+0.6791, M=-1.3583` — real OPC, one virtual site per
water — confirming the packing-by-site-count rule produces true parameters, not TIP4P-Ew.

## restraints.py

Exactly the openMM cookbook idiom. `k` is a global parameter, so release is a single
`context.setParameter("k", 0.0)` — no force removal, no context rebuild.

```python
import openmm
from openmm import unit


def add_restraint(system, positions, atom_indices, k):
    force = openmm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    for i in atom_indices:
        force.addParticle(i, positions[i].value_in_unit(unit.nanometer))   # add-order == iter-order
    system.addForce(force)
    return force


def reanchor(force, context, positions, atom_indices):
    """Move restraint references to current coords (call once after minimization, so
    equilibration restrains to the MINIMIZED pose — matching the GROMACS non-drifting reference)."""
    for add_idx, atom_idx in enumerate(atom_indices):
        force.setParticleParameters(add_idx, atom_idx, positions[atom_idx].value_in_unit(unit.nanometer))
    force.updateParametersInContext(context)
```

## density.py

Port `RELE_simulations/gromacs_REMD/scripts/simulation/density_converged.py` **verbatim** —
pure numpy (least-squares slope/plateau test over the trailing window). Its tests port over
too, giving a passing suite on day one.

```python
def density_converged(volumes: list[float], tol_rel: float, min_seg: int) -> bool:
    """True when a least-squares line through the trailing `min_seg` volumes accounts for a
    fractional drift <= tol_rel. Plateau test, not a consecutive-segment difference."""
    ...  # verbatim port
```

## stages.py — equilibration primitives (shared)

`tau_p` is gone (MC barostat has no coupling time — uses a volume-move frequency, default 25).
Single integrator/thermostat for the whole system (no GROMACS Protein/Non-Protein split).

```python
from openmm import unit, MonteCarloBarostat
from .density import density_converged


class DensityNotConverged(RuntimeError):
    pass


def minimize(sim, cfg):
    sim.minimizeEnergy()
    # fail loud: a non-finite potential energy after minimization means a bad topology/params
    # (bad template, clashing added atoms) — do not let it reach dynamics.
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy()
    import math
    assert math.isfinite(pe.value_in_unit(pe.unit)), f"non-finite energy after minimization: {pe}"


def equilibrate_nvt(sim, cfg):                   # was GROMACS "heat"
    seed = cfg.seed if cfg.seed is not None else 0
    sim.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin, seed)
    sim.step(cfg.nvt_equil_steps)


def add_barostat(sim, cfg):
    sim.system.addForce(MonteCarloBarostat(cfg.ref_p_bar * unit.bar,
                                           cfg.temperature_k * unit.kelvin))
    sim.context.reinitialize(preserveState=True)  # required after addForce


def equilibrate_density(sim, cfg) -> int:
    volumes = []
    for seg in range(1, cfg.density_max_seg + 1):
        sim.step(cfg.density_seg_steps)
        v = sim.context.getState().getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
        volumes.append(v)
        if seg >= cfg.density_min_seg and \
           density_converged(volumes, cfg.density_tol_rel, cfg.density_min_seg):
            return seg
    raise DensityNotConverged(f"no volume plateau in {cfg.density_max_seg} segments: {volumes}")


def relax(sim, cfg):
    sim.step(cfg.relax_steps)
```

## prepare.py — the shared pipeline

Runs build → restrain → minimize → NVT → density → (relax) → release restraints on a throwaway
equilibration `Simulation`, then hands back plain openMM objects. Takes a `PrepConfig`, so it is
literally the same call for MD, REMD and REST2.

```python
from openmm import app, unit, LangevinMiddleIntegrator, Platform
from . import stages
from .build import build_system
from .restraints import add_restraint, reanchor
from .equilibrated import EquilibratedSystem


def prepare(cfg) -> EquilibratedSystem:          # cfg: PrepConfig (or any subclass)
    cfg.outdir.mkdir(parents=True, exist_ok=True)
    built = build_system(cfg)
    heavy = protein_heavy_atoms(built.topology)  # assert non-empty
    restraint = add_restraint(built.system, built.positions, heavy, cfg.restraint_k)

    # Equilibration integrator (separate from any production integrator).
    integ = LangevinMiddleIntegrator(cfg.temperature_k * unit.kelvin,
                                     cfg.friction_ps / unit.picosecond,
                                     cfg.dt_ps * unit.picosecond)
    if cfg.seed is not None:
        integ.setRandomNumberSeed(cfg.seed)

    platform = Platform.getPlatform(cfg.platform)
    props = {"Precision": "mixed"} if cfg.platform in ("CUDA", "OpenCL", "HIP") else {}
    sim = app.Simulation(built.topology, built.system, integ, platform, props)
    sim.context.setPositions(built.positions)

    stages.minimize(sim, cfg)
    pos = sim.context.getState(getPositions=True).getPositions()
    reanchor(restraint, sim.context, pos, heavy)     # restrain to minimized pose
    stages.equilibrate_nvt(sim, cfg)                 # NVT
    stages.add_barostat(sim, cfg)                    # -> NPT
    stages.equilibrate_density(sim, cfg)             # restrained, iterative

    sim.context.setParameter("k", 0.0)               # release restraints
    if cfg.relax_ns > 0:
        stages.relax(sim, cfg)

    state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    eq = EquilibratedSystem(built.system, built.topology,
                            state.getPositions(), state.getPeriodicBoxVectors())
    eq.to_disk(cfg.outdir)                            # equilibrate once, reuse for any backend
    return eq
```

Note: the restraint force stays in the System with `k=0`. For REMD/REST2 that is harmless
(zero energy/force), but the backends may strip it for cleanliness before building replica
states — a one-liner on the serialized system.

## produce.py — plain-MD backend + driver

```python
from pathlib import Path
from openmm import app, unit, LangevinMiddleIntegrator, Platform
from .prepare import prepare


def produce_md(eq, cfg) -> Path:                 # eq: EquilibratedSystem, cfg: MDConfig
    integ = LangevinMiddleIntegrator(cfg.temperature_k * unit.kelvin,
                                     cfg.friction_ps / unit.picosecond,
                                     cfg.dt_ps * unit.picosecond)
    if cfg.seed is not None:
        integ.setRandomNumberSeed(cfg.seed)
    platform = Platform.getPlatform(cfg.platform)
    props = {"Precision": "mixed"} if cfg.platform in ("CUDA", "OpenCL", "HIP") else {}

    sim = app.Simulation(eq.topology, eq.system, integ, platform, props)
    sim.context.setPeriodicBoxVectors(*eq.box_vectors)
    sim.context.setPositions(eq.positions)
    sim.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin,
                                           cfg.seed if cfg.seed is not None else 0)

    sim.reporters += [
        app.XTCReporter(str(cfg.outdir / "prod.xtc"), cfg.traj_steps),   # compact; mdtraj reads it
        app.StateDataReporter(str(cfg.outdir / "prod.log"), cfg.traj_steps,
            step=True, time=True, potentialEnergy=True, temperature=True, volume=True,
            density=True, progress=True, remainingTime=True, speed=True, totalSteps=cfg.prod_steps),
        app.CheckpointReporter(str(cfg.outdir / "prod.chk"), cfg.checkpoint_steps),
    ]
    sim.step(cfg.prod_steps)

    export_final_pdb(sim, cfg.outdir)
    sim.saveState(str(cfg.outdir / "final.xml"))
    sim.saveCheckpoint(str(cfg.outdir / "final.chk"))
    (cfg.outdir / "config.json").write_text(cfg.model_dump_json(indent=2))
    return cfg.outdir


def run_md(cfg) -> Path:                          # cfg: MDConfig
    return produce_md(prepare(cfg), cfg)
```

## Usage

The interface is the config object; YAML/CLI are optional sugar. Scripting is a plain loop:

```python
from openmm_md import MDConfig, run_md, ProteinFF, WaterModel

base = MDConfig(pdb_in="ww.pdb", outdir="out",
                protein_ff=ProteinFF.CHARMM36M, water=WaterModel.TIP3P, total_ns=200)
for t in (300, 310, 320, 330):
    run_md(base.model_copy(update={"temperature_k": t, "outdir": f"out/T{t}"}))
```

---

## Extending to REMD/REST2 (design target, not yet built)

The point of the prepare/produce seam: adding T-REMD is a new `produce_*` backend plus a
`PrepConfig` subclass — **zero changes to build/restraints/density/stages/prepare.** Sketch of
the T-REMD backend (lives here under `[remd]` extra, or in `openmm_enh_samp`, importing
`prepare`):

```python
from openmm import unit
from openmmtools import states, mcmc, multistate
from .prepare import prepare
from .ladder import geometric_ladder


def run_remd(cfg) -> Path:                        # cfg: REMDConfig
    eq = prepare(cfg)                             # SHARED — identical to MD
    temps = cfg.temps_list or geometric_ladder(cfg.temperature_k, cfg.t_max_k, cfg.n_replicas)

    system = strip_zero_restraint(eq.system)     # optional cleanup of the k=0 force
    thermo_states = [
        states.ThermodynamicState(
            system, temperature=T * unit.kelvin,
            pressure=(cfg.ref_p_bar * unit.bar if cfg.ensemble is Ensemble.NPT else None),
        )
        for T in temps
    ]
    sampler_state = states.SamplerState(positions=eq.positions, box_vectors=eq.box_vectors)

    exchange_steps = round(cfg.exchange_ps / cfg.dt_ps)
    move = mcmc.LangevinDynamicsMove(
        timestep=cfg.dt_ps * unit.picosecond,
        collision_rate=cfg.friction_ps / unit.picosecond,
        n_steps=exchange_steps,
    )
    sampler = multistate.ReplicaExchangeSampler(
        mcmc_moves=move, number_of_iterations=cfg.iterations)
    reporter = multistate.MultiStateReporter(
        str(cfg.outdir / "remd.nc"), checkpoint_interval=...)
    sampler.create(thermo_states, sampler_state, reporter)
    sampler.run()                                # handles exchanges, checkpointing, .nc output
    return cfg.outdir
```

What each engine adds on top of the shared base:

- **T-REMD** — a temperature ladder + `ReplicaExchangeSampler` (or `ParallelTemperingSampler`,
  which builds the ladder itself). Production ensemble NVT or NPT via the `ThermodynamicState`
  pressure argument.
- **REST2** — the same `multistate` machinery, but each replica's `ThermodynamicState` scales
  only the **solute** terms by λ (effective-temperature ladder at one physical temperature).
  openMM/openmmtools does not ship a turnkey REST2; the solute-scaling is the real work
  (a `CompoundThermodynamicState` with a REST-region state, or a per-λ modified system).
  **This is the one genuinely involved piece** — budget for it — but it still consumes
  `prepare()` unchanged. It also inherits the GROMACS lesson that not every force field is
  REST2-safe (CMAP), so the same energy-match check applies before trusting it.

Analysis: `multistate` writes a NetCDF (`.nc`) reporter with all replicas + exchange stats;
`openmmtools.multistate.MultiStateReporter` / `MultiStateSamplerAnalyzer` read it. Acceptance
rates and round-trip mixing (the GROMACS `gromd-acceptance` / `gromd-roundtrip` metrics) come
from that file — a separate analysis concern, not part of this engine.

**Design implications already baked in for REMD/REST2:**
1. `prepare()` returns serializable openMM objects (System xml + State), exactly what
   `multistate` consumes — no `app.Simulation` assumption leaks past the seam.
2. `temperature_k` is the *reference/equilibration* temperature (= T_min), so equilibration is
   already engine-agnostic; REMD just adds `t_max_k` + `n_replicas`.
3. Config uses inheritance from `PrepConfig`, so a `REMDConfig` reuses every shared field and
   validator with no duplication.
4. `openmmtools` is an optional extra, keeping the plain-MD install lean.

---

## GROMACS-ism audit

What was deliberately dropped, renamed, or kept when porting — recorded so the port does not
silently drag GROMACS conventions into openMM.

### Dropped (openMM handles it, or it was GROMACS-specific machinery)

- **Two-group thermostat** (`tc-grps = Protein Non-Protein`): openMM applies one Langevin
  friction to the whole system. The solute/solvent split is a discouraged GROMACS practice
  ("hot solvent / cold solute"); a single thermostat is standard and correct.
- **`nstcomm` / COM-motion removal**: `createSystem` adds a `CMMotionRemover` by default.
- **`constraint-algorithm = LINCS`**: openMM chooses CCMA/SETTLE internally.
- **Scratch folder-symlink model** (`SYMLINK_BULK`, `PRESERVE_SCRATCH_FROM`, `SCRATCH_DIR`,
  cross-fs `rename()` gotchas): reporters write wherever `outdir` points; set it to scratch.
- **`tau_p`**: the MC barostat has no coupling time (uses a move frequency; default kept).
- **Global `CUTOFF_NM` knob + `charmm*` string branch + `DispCorr`**: cutoff/switching are now
  a property of the force field (`NonbondedSpec`); dispersion correction is the `createSystem`
  default. No `if FF == charmm*` special-casing.
- **`SEED = -1` sentinel**: replaced by `seed: int | None = None`.

### Renamed (values were GROMACS/AMBER mental-model leaks)

- `gamma_ln` (2.0) -> `friction_ps` (1.0, openMM convention)
- `heat` / `HEAT_NS` -> `equilibrate_nvt` / `nvt_equil_ns` ("heat" implied a temperature ramp
  that never happened — it is NVT equilibration with velocities set once)
- `box_buffer` -> `padding_nm` (openMM Modeller term)
- `salt_molar` -> `ionic_strength_molar` (openMM addSolvent term)
- `T_SIM` / `T_MIN` -> `temperature_k` (the reference/equilibration temperature, engine-agnostic)

### Kept deliberately (ported science, not an openMM convention)

- **Iterative volume-plateau density equilibration** (`density_*` knobs + `density_converged`):
  openMM's own convention is a single fixed-length NPT equilibration (e.g. 1 ns). The adaptive
  plateau loop is carried over because it is genuinely more robust (catches a slowly
  contracting box) and reuses already-tested code. Simpler alternative if it feels heavy: a
  fixed-length NPT block.
- **Restrained EM + restraints held through NVT/density, reference = minimized pose**: the
  project's science goal (preserve the designed pose during equilibration), not an artifact.
- **`restraint_k = 1000 kJ/mol/nm^2`** on protein heavy atoms: matches the GROMACS POSRES
  default (the cookbook example uses 100 on CA only); a stronger, pose-preserving choice.

### Adopted from openMM (were not in the GROMACS engine)

- **4 fs timestep + hydrogen mass repartitioning** (`dt_ps=0.004`, `hydrogen_mass_amu=4`):
  standard modern openMM practice, ~2x throughput. Drop to 2 fs via `dt_ps=0.002`,
  `hydrogen_mass_amu=1.008` (validation enforces HMR when dt > 2.5 fs).
- **Mixed precision** on GPU platforms (the recommended MD default).
- **prepare/produce split + serialized equilibrated system**: enables equilibrate-once /
  produce-many and the REMD/REST2 seam — no GROMACS analogue.

---

## Structure-prep guards (implemented in build.py)

- **Disulfide preservation** — SG-SG bond count is captured before fixing and asserted unchanged
  after; a change raises `PrepareError`. Critical for disulfide-constrained designs.
- **Water topology** — after solvation, every `HOH` must have the water model's expected particle
  count, and 4-site models must have virtual sites, else `PrepareError`. Backstops the
  packing-by-site-count rule against a future OpenMM box-reuse change.

## Science-review flags (resolve before production, not port-mechanical)

1. **Restraint reference vs MC barostat.** Restraint `x0` are absolute coordinates; the MC
   barostat rescales molecule centers of mass. GROMACS used `refcoord-scaling = com`. Confirm
   the openMM interaction is benign (restraints on solute, barostat scales by molecule).
2. **Final-PDB PBC wrapping.** GROMACS used `trjconv -pbc mol -center -ur compact`; do the
   equivalent with `mdtraj` image_molecules on export.
3. **`protein_heavy_atoms` selection.** Must match GROMACS `-DPOSRES` (protein heavy atoms, no
   hydrogens). Assert non-empty.
4. **CHARMM energy-match.** Before trusting any CHARMM run, energy-match a known system to
   confirm switching / `CustomNonbondedForce` behavior. AMBER needs no such check.
5. **REST2 force-field safety (future).** CMAP-bearing force fields (CHARMM, AMBER ff19SB) are
   not safely solute-scalable without handling the cross-term — the GROMACS engine rejected
   `charmm*` for REST2 for exactly this. Carry that guard into the REST2 backend.

## Settled defaults

- **4 fs + HMR** (`dt_ps=0.004`, `hydrogen_mass_amu=4.0`). Drop to 2 fs via `dt_ps=0.002`,
  `hydrogen_mass_amu=1.008` per run if needed.
- **Adaptive density loop kept** — NPT in segments until the box volume plateaus, fail-loud if it
  never does. The one deliberately-ported GROMACS protocol; correct + fail-loud, so retained.
- **Default pH 7.4** (physiological).
- **Trajectory: XTC** via `app.XTCReporter` (verified present in openMM 8.6) — compact, mdtraj reads it.
- **mdtraj is the only trajectory library** (no MDAnalysis); analysis is a new in-package module.

## Analysis

A **new, in-package analysis subpackage** (`openmm_md/analysis/`, `mdtraj`-based) — the GROMACS
`gromd_analysis` package is **not** imported. `mdtraj` is the single trajectory library (no
MDAnalysis). Production writes `prod.xtc`; analysis reads it + the topology and measures drift from
the **input pose** (the reference): RMSD, Rg, per-residue RMSF, DSSP, and conformational clustering.

It is a subpackage of small single-purpose modules (not one file) because it will grow:
`trajectory` (load + PBC image + solvent strip + align), `metrics`, `dssp`, `clustering`, `plots`,
with REMD-specific metrics (exchange acceptance, round-trip mixing) added alongside the REMD backend
later. Each module is independently callable; `analysis/__init__.py` re-exports the entry points.
- REMD production ensemble default: NVT (GROMACS default) or NPT?
- Where the REMD/REST2 backends ultimately live: in this package under a `[remd]` extra, or in
  a separate package that imports `openmm_md.prepare`. (Either works given the seam.)
