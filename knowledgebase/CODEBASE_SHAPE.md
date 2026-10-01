# CODEBASE_SHAPE — openmm_pipelines

> Status (2026-10-01): **design sketch only.** No `.py` files exist yet. This
> document describes the *intended* architecture from `SKELETON_v2.md`; update it
> to reflect reality as modules land. The authoritative current design is always
> `SKELETON_v2.md` at the repo root until implementation catches up.

## What this package is

A small, Pythonic OpenMM engine for **explicit-solvent protein MD**. Given a
folded/designed input PDB, it runs the standard flow — structure prep + solvation,
energy minimization, equilibration (NVT then NPT density), production — and the
analysis measures how far the structure drifts from the **input pose** (RMSD =
drift from input, RMSF = local flexibility). It never folds from an unfolded
state; the input pose is the analysis reference.

It is a `tools/` library — same tier as `snekwrap` and `FragForge` — imported by
campaigns via `import openmm_pipelines`. Scope now is **plain production MD**;
T-REMD and REST2 are designed-for but not built.

## The central idea: the prepare/produce seam

Everything up to production is **identical across MD, T-REMD and REST2**. Only
production differs. So the engine splits at exactly that line:

```
prepare(cfg) ---------------------------> EquilibratedSystem      (SHARED — ~90% of code)
                                          (System, topology, positions, box_vectors)
                                                   |
                    +------------------------------+------------------------------+
                    v                              v                              v
           md.run(eq, cfg)              remd.run(eq, cfg)              rest2.run(eq, cfg)
           app.Simulation.step()        multistate.ReplicaExchange     multistate + REST region
           (this package, now)          (later; imports prepare)       (later; imports prepare)
```

`prepare()` returns **plain OpenMM objects** (serialized `System` + a `State` with
positions and box vectors) — exactly what both `app.Simulation` and
`openmmtools.multistate` consume. The production backends share nothing but this
handoff. `prepare()` also writes the equilibrated system to disk, so equilibration
runs once and any production backend (or a separate job) launches from the same box.

Why this matters: adding T-REMD is a new engine folder (`remd/`) + a `PrepConfig`
subclass — **zero changes** to `lib/` or `preparation/`.

## Module layout (tiers: types · lib · pipeline · engines)

Reorganized from SKELETON's flat list (see `DECISIONS.md`, 2026-10-01). **Top level
keeps only `__init__.py` + `config.py` + `forcefield.py`** (shared public types);
everything else is a folder. Three tiers: shared **types/registries** at root, the
reusable **toolkit** in `lib/`, the shared **pipeline** in `preparation/`, and one
folder per **engine**. `[BUILT]` = exists; everything else is planned.

```
openmm_pipelines/
  __init__.py              # public API
  config.py                # PrepConfig (base), BoxShape, derived step counts   [BUILT]
  forcefield.py            # ProteinFF, WaterModel, NonbondedSpec, resolve(),
                           #   ForceFieldError, registries                      [BUILT]

  lib/                     # the toolkit: reusable building blocks (not stages)
    restraints.py          #   add_restraint (global k), reanchor               [BUILT]
    density_convergence.py #   assess_plateau, density_converged (pure math)    [BUILT]

  preparation/             # the prepare side of the seam (shared by all engines)
    __init__.py
    prepare.py             #   prepare(cfg) -> EquilibratedSystem  (composer)   [soon]
    build.py               #   build_system -> BuiltSystem; disulfide+water guards [BUILT]
    minimize.py            #   minimize (static energy min + finite-E guard)    [BUILT]
    equilibration_steps.py #   equilibrate_nvt, add_barostat, equilibrate_density [BUILT]
    handoff.py             #   BuiltSystem, EquilibratedSystem (to_disk/from_disk(outdir)) [stub]

  md/                      # plain production MD  [NOW]
    __init__.py            #   -> run, MDConfig
    config.py              #   MDConfig(PrepConfig)                             [BUILT]
    production.py          #   run(eq, cfg) -> Path                             [soon]

  remd/                    # T-REMD  [SOON; needs [remd] extra]
    config.py              #   REMDConfig(PrepConfig)  (+ per-replica equilibration)
    ladder.py              #   geometric_ladder, acceptance helpers
    production.py          #   run(eq, cfg) -> Path  (openmmtools multistate)

  rest2/                   # REST2  [SOON; needs [remd] extra]
    config.py              #   REST2Config(PrepConfig)
    production.py          #   run(eq, cfg) -> Path  (multistate + solute scaling; CMAP gate)

  analysis/                # mdtraj-based, in-package  [BUILT]
    __init__.py            #   analyze(outdir) composer + re-exports
    trajectory.py          #   load (cif topo via openmm) + image_molecules + strip + superpose to input
    metrics.py             #   rmsd_to_reference / radius_of_gyration / rmsf_per_residue (Å)
    dssp.py                #   secondary_structure, ss_fractions (mdtraj compute_dssp)
    clustering.py          #   pairwise_rmsd + cluster_conformations (scipy hierarchical)
    plots.py               #   rmsd / rg / rmsf / dssp PNGs (matplotlib Agg)
    remd.py                #   [SOON] exchange acceptance, round-trip mixing (reads .nc)
tests/
  test_forcefield.py (10) · test_config.py (13) · test_density.py (10)   [BUILT, pass]
  test_restraints.py (6) · test_build.py (7) · test_equilibration.py (5, slow) [BUILT, pass]
  gpu_tests.sbatch         # runs the suite on CUDA (-w node3620)
```

**Tier rule — what goes in `lib/`:** a *substantial, reusable* building block —
used across multiple stages (`restraints`: add → reanchor → release span the
pipeline) or a pure, independently-testable unit (`density_convergence`: numpy, no
OpenMM). Trivial glue for one stage stays with that stage (`add_barostat` lives in
`equilibration_steps.py`, used only at the NVT→NPT transition).

**Honest stage names:** `minimize` (static) and the dropped `relax` are *not*
equilibration, so `equilibration_steps.py` holds only the restrained equilibration
dynamics (NVT + NPT density). `minimize` is its own file.

**Naming (seam verbs).** `prepare` (setup, `core/`) + `production` (run side). Each
engine exposes one public function `run(eq, cfg) -> Path` — `md.run`, `remd.run`,
`rest2.run`, same signature. No `produce_*`; no `cfg`-only wrapper. One-shot is
explicit: `eq = prepare(cfg); md.run(eq, cfg)` (one config object, since `MDConfig`
IS-A `PrepConfig`).

Package import name is `openmm_pipelines` (underscore) — deliberately **not**
`openmm`, to avoid shadowing the real package.

## Config shape

Config is split by the same seam: a shared `PrepConfig` base (everything through
equilibration) and a thin `MDConfig` subclass adding production fields. Future
`REMDConfig`/`REST2Config` subclass the **same** `PrepConfig`. Inheritance (not
nesting) keeps attribute access flat — `cfg.temperature_k` works on every engine's
config — while `prepare()` depends only on `PrepConfig` fields. Pydantic
`BaseModel, frozen=True`; step counts are derived `@property`s off ns + dt.

`temperature_k` is the *reference/equilibration* temperature (= T_min for REMD), so
equilibration is already engine-agnostic.

## Pipeline phases

| Phase | OpenMM mechanism | side |
|---|---|---|
| structure prep + solvation | `pdbfixer` + `ForceField` + `Modeller.addSolvent` | prepare |
| minimization (restrained) | `simulation.minimizeEnergy()` | prepare |
| NVT equilibration (restrained) | `setVelocitiesToTemperature` + Langevin steps | prepare |
| NPT density equil (restrained, adaptive plateau) | `MonteCarloBarostat` + segment loop | prepare |
| restraint release (+ optional unrestrained NPT relax) | `context.setParameter("k", 0)` | prepare |
| production (NPT or NVT) | MD: `Simulation.step()`; REMD/REST2: `multistate` | **produce** |
| export + analysis | `PDBFile` + `mdtraj` PBC wrap; in-package `mdtraj` analysis | produce |

Restraints are held at full strength through the **entire** equilibration and
released only at the **start of production**, so the production trajectory captures
drift from the input pose — not drift accumulated during equilibration.

## Dependencies

Core OpenMM bundles every protein force field needed — **no `openmmforcefields`
dependency** (verified against OpenMM 8.6.1 in `openmm_env`). `mdtraj` is the single
trajectory library (no MDAnalysis). `openmmtools` is required only for the future
REMD/REST2 backends — an **optional `[remd]` extra**, not a core dep of plain MD.

## Design conventions (quick reference)

- One `LangevinMiddleIntegrator` = integrator + thermostat for the whole system; no
  per-group thermostats.
- `MonteCarloBarostat` at default volume-move frequency.
- Nonbonded cutoff + vdW switching are **properties of the force field**
  (`NonbondedSpec`), passed through `createSystem` kwargs — never special-cased by
  FF name, never read back off the built System.
- `seed=None` = fresh RNG per run; an int pins setup + RNG streams (GPU still not
  bit-for-bit reproducible).
- 4 fs timestep with HMR (`hydrogenMass=4 amu`, `constraints=HBonds`); drop to 2 fs
  via physical H mass + `dt_ps=0.002`.
- Mixed precision on GPU.
- Water packing chosen by **site count**, parameters by **force field** (the
  `_PACKING_BOX` rule); FF/water pairing enforced by a single `resolve()`.

See `DECISIONS.md` for the *why* behind each.
