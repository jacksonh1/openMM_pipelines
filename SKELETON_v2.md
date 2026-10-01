# openmm_pipelines — design & architecture

Status: **plain MD implemented** (prepare → run → analyze, end-to-end, tested on GPU).
T-REMD / REST2 are designed-for but not built. This file is the design overview; the
living architecture map is `knowledgebase/CODEBASE_SHAPE.md`, the rationale is
`knowledgebase/DECISIONS.md`, and the code is the source of truth for signatures.

## What this is

A small, Pythonic OpenMM engine for **explicit-solvent protein MD**. Given a folded /
designed input PDB it runs the standard flow — structure prep + solvation, energy
minimization, equilibration (NVT then NPT density), production — and the analysis
measures how far the structure drifts from the **input pose**, which is the reference
(RMSD = drift from design, RMSF = local flexibility). It never folds from an unfolded
state.

It is a `tools/` library (same tier as `snekwrap`, `FragForge`), imported by campaigns
via `import openmm_pipelines`.

## The prepare/produce seam

Everything up to production is **identical across MD, T-REMD and REST2**: build the
solvated system, restrain, minimize, equilibrate NVT, equilibrate NPT density, release
restraints. Only **production** differs:

- **plain MD** — one trajectory via `app.Simulation.step()`.
- **T-REMD / REST2** — many replicas via an `openmmtools.multistate` sampler (its own
  `Context`s from `ThermodynamicState` + `SamplerState`); does **not** use `app.Simulation`.

```
prepare(cfg) ---------------------------> EquilibratedSystem        (SHARED — ~90% of code)
                                          (System, topology, positions, box_vectors)
                                                   |
                    +------------------------------+------------------------------+
                    v                              v                              v
            md.run(eq, cfg)              remd.run(eq, cfg)              rest2.run(eq, cfg)
            app.Simulation.step()        multistate.ReplicaExchange     multistate + REST region
            (now)                        (later; imports prepare)       (later; imports prepare)
```

`prepare()` returns plain OpenMM objects and serializes them (portable System + State
XML + mmCIF topology), so equilibration runs once and any production backend — or a
separate job on another node — launches from the same box. Validated against OpenMM's
docs (seam placement, the handoff, the restraint idiom). Adding T-REMD is a new engine
folder + a `PrepConfig` subclass — zero changes to `lib/` or `preparation/`.

## Module layout

Four tiers: shared **types** at root, the reusable **toolkit** in `lib/`, the shared
**pipeline** in `preparation/`, one folder per **engine**. Only `config.py` +
`forcefield.py` sit at the top level; everything else is a folder.

```
openmm_pipelines/
  __init__.py
  config.py                # PrepConfig (base), BoxShape, derived step counts
  forcefield.py            # ProteinFF, WaterModel, NonbondedSpec, resolve, ForceFieldError
  lib/                     # toolkit (reusable building blocks, not stages)
    restraints.py          #   add_restraint (global k), reanchor
    density_convergence.py #   assess_plateau, density_converged (pure numpy)
    selections.py          #   protein_heavy_atoms
  preparation/             # the prepare side of the seam (shared)
    prepare.py             #   prepare(cfg) -> EquilibratedSystem   (composer; the order lives here)
    build.py               #   build_system -> BuiltSystem; disulfide + water guards
    minimize.py            #   minimize (static energy min + finite-E guard)
    equilibration_steps.py #   equilibrate_nvt, add_barostat, equilibrate_density
    handoff.py             #   BuiltSystem, EquilibratedSystem (to_disk / from_disk)
  md/                      # plain production MD  [NOW]
    config.py              #   MDConfig(PrepConfig)
    production.py          #   run(eq, cfg) -> Path
  remd/  rest2/            # [SOON] config.py + production.py (+ ladder.py); need [remd] extra
  analysis/                # mdtraj-based; trajectory / metrics / dssp / clustering / plots + analyze()
  reports.py               # RunInfo / PrepReport / RunReport (typed output records)
```

Package import name is `openmm_pipelines` (underscore) — deliberately not `openmm`.

## Pipeline stages (what `prepare()` composes, in order)

| # | stage | module | restraints |
|---|---|---|---|
| 1 | build (PDBFixer → solvate → createSystem) | `preparation/build.py` | — |
| 2 | add restraints (protein heavy atoms, global `k`) | `lib/restraints.py` | on |
| 3 | minimize (static) | `preparation/minimize.py` | on |
| 4 | reanchor restraints to minimized pose | `lib/restraints.py` | on |
| 5 | equilibrate NVT | `preparation/equilibration_steps.py` | on |
| 6 | add barostat (→ NPT) | `preparation/equilibration_steps.py` | on |
| 7 | equilibrate density (adaptive plateau) | `preparation/equilibration_steps.py` | on |
| 8 | release restraints (`k`→0, baked into System default) | `preparation/prepare.py` | → off |
| 9 | hand off (serialize EquilibratedSystem) | `preparation/handoff.py` | off |

Restraints are held full-strength through all of equilibration and released only at the
end, so production (which starts at release) records all drift from the input pose.
There is **no relax step** — it would let early drift happen off-camera, contradicting
the analysis contract (enhanced-sampling pre-equilibration is a per-replica REMD/REST2
backend concern instead).

## Config

`PrepConfig` (pydantic, frozen) holds everything through equilibration; `MDConfig`
subclasses it with production fields (`total_ns`, `traj_ps`). Future `REMDConfig` /
`REST2Config` subclass the same base. `temperature_k` is the reference/equilibration
temperature (= T_min for REMD). Derived step counts are `@property`s off ns + `dt_ps`.

Force field + water are enums; `resolve(protein_ff, water)` returns the FF xml list +
`addSolvent` packing box and raises on an illegal pairing. Water packing is chosen by
**site count**, parameters by **force field** (so any non-polarizable water works with
one registry entry). Nonbonded cutoff/switch are a **property of the FF**
(`NonbondedSpec`) passed through `createSystem` kwargs — OpenMM does not pick these for
you, and CHARMM needs its switch (see GOTCHAS).

Defaults: AMBER14SB + TIP3P, dodecahedron box, 1.0 nm padding, 0.15 M NaCl, 300 K,
**4 fs + HMR**, mixed precision, CUDA, adaptive density plateau.

## Output contract (first draft — see DECISIONS.md)

```
outdir/
  config.json              # all config fields + derived step counts
  run_info.json            # versions, platform/precision, resolved FF xmls, host, SLURM id
  input.pdb                # exact input, copied
  prepared/  system.xml  state.xml (portable saveState)  topology.cif  prep_report.json
  production/ trajectory.xtc  production.log  final_state.xml  production.chk  run_report.json
  analysis/  analysis_report.json  analysis_data.npz  rmsd.png rg.png rmsf.png dssp.png
```

`state.xml` / `final_state.xml` are portable (`saveState`); `production.chk` is an exact
same-hardware restart only (checkpoints are not portable). Reports are typed dataclasses
(`reports.py`).

## Analysis

`analysis/` (mdtraj): load xtc with the mmCIF topology (parsed via OpenMM), image
molecules, strip to protein, superpose onto the input pose; then RMSD-to-design, Rg,
per-residue RMSF, DSSP, conformational clustering (scipy). `analyze(outdir)` runs them
all and writes `outdir/analysis/`.

## Dependencies

Core OpenMM bundles every protein FF needed (no `openmmforcefields`). `mdtraj` is the
only trajectory library; `pandas` parses the production log; `matplotlib` for plots;
`scipy` for clustering. `openmmtools` is required only for the future REMD/REST2
backends — an optional `[remd]` extra. Everything runs in the `openmm_env` conda env.

## Extending to REMD / REST2 (design target)

New engine folder (`remd/`, `rest2/`) + a `PrepConfig` subclass; consume `prepare()`
unchanged. T-REMD = temperature ladder + `ReplicaExchangeSampler`. REST2 = the same
multistate machinery with per-replica solute λ-scaling (the genuinely involved piece;
CMAP force fields need a safety check). Per-replica equilibration lives in these
backends (not shared `prepare()`). Open: production ensemble default (NVT/NPT);
job-submission story (CLI/sbatch generator vs left to campaigns).
