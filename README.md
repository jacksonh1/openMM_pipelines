# openmm_pipelines

A small, Pythonic OpenMM engine for **explicit-solvent MD of designed protein
structures** on single-node GPU clusters. Three protocols — **plain production MD**,
**T-REMD** (temperature replica exchange), and **REST2** (replica exchange with solute
tempering) — all built behind one shared *prepare → produce* seam.

Every run starts from a **folded/designed input structure**, never from an unfolded or
extended state. The input pose is the reference the analysis is measured against: RMSD is
drift from the design, RMSF is local flexibility. This is a characterization tool, not a
folding tool.

It is a `tools/` library in the `09-fragfold` workspace — same tier as `snekwrap` and
`FragForge` — imported by design campaigns via `import openmm_pipelines`.

---

## What it's for

| # | Objective | Protocol |
|---|---|---|
| 1 | Stability / rigidity of a designed structure | plain MD · T-REMD · REST2 |
| 2 | Identify flexible regions (per-residue RMSF) | plain MD · T-REMD · REST2 |
| 3 | Variant comparison — which design best retains its fold | T-REMD · REST2 |
| 4 | Bound-state sampling — simulate a complex in its bound pose | plain MD |

Plain MD is mainly for **#4** (and optionally #1/#2); T-REMD and REST2 are primarily for
enhanced sampling of **#1–3**. In explicit solvent, **REST2 is usually the right enhanced-
sampling tool** — it tempers only the solute, so it reaches a given effective temperature
with a handful of replicas where T-REMD needs many (replica count scales with √(box DOF)).

---

## Install

The package has its **own isolated conda env, `openmm_env`** — do not use the sibling
repos' envs.

```bash
conda env create -f environment.yml     # creates openmm_env (openmm 8.6, openmmtools, mdtraj, …)
conda activate openmm_env
pip install -e .                         # editable install (done by environment.yml's pip section)
```

- `openmm` + `pdbfixer` ship from **conda-forge** (not pip).
- `openmmtools` powers the T-REMD / REST2 backends; it's the optional `[remd]` extra for a
  pip-only install but is already in `environment.yml`.
- Core OpenMM bundles every protein force field used here — **no `openmmforcefields`
  dependency**.

Makefile shortcuts: `make create_environment`, `make requirements`, `make test`, `make format`.

---

## Quickstart

### Plain production MD

```python
from openmm_pipelines import MDConfig, ProteinFF, WaterModel, prepare, md, analysis

cfg = MDConfig(
    pdb_in="design.pdb",
    outdir="run1",
    protein_ff=ProteinFF.AMBER14SB,
    water=WaterModel.TIP3P,
    temperature_k=300.0,
    total_ns=100.0,
    traj_ps=10.0,
    platform="CUDA",
    seed=12345,                 # None = fresh RNG per run
)

eq = prepare(cfg)               # build + solvate + minimize + equilibrate (shared by all engines)
md.run(eq, cfg)                 # production
analysis.analyze(cfg.outdir)    # RMSD/Rg/RMSF/DSSP + clustering vs the input pose
```

One-shot: `md.run(prepare(cfg), cfg)`. `prepare()` also serializes the equilibrated system
to disk, so production (or a separate job) can relaunch from the same box.

### T-REMD

```python
from openmm_pipelines import ProteinFF, WaterModel, prepare, remd
from openmm_pipelines.remd.config import REMDConfig
from openmm_pipelines.analysis.remd import analyze_remd

cfg = REMDConfig(
    pdb_in="design.pdb", outdir="run_remd",
    temperature_k=300.0,            # T_min = reference / design temperature (ladder minimum)
    max_temperature_k=400.0,
    n_replicas=12,                  # geometric ladder from min/max/n; or pin temperatures_k=(...)
    total_ns=20.0,                  # MD time collected per replica
    exchange_attempt_ps=1.0,        # swap attempt every 1 ps
    platform="CUDA",
)
remd.run(prepare(cfg), cfg)
analyze_remd(cfg.outdir, reference_state=0)   # state 0 = T_min, the Boltzmann sample at design T
```

### REST2

```python
from openmm_pipelines import ProteinFF, WaterModel, prepare, rest2
from openmm_pipelines.rest2.config import REST2Config
from openmm_pipelines.analysis.rest2 import analyze_rest2

cfg = REST2Config(
    pdb_in="design.pdb", outdir="run_rest2",
    temperature_k=300.0,                   # physical T_0 — every replica runs here
    max_effective_temperature_k=450.0,     # top rung's EFFECTIVE temperature
    n_replicas=6,                          # λ ladder from T_0/max/n; or pin lambdas=(1.0, …)
    solute_chain_index=None,               # None = temper whole protein; else one chain
    total_ns=20.0,
    exchange_attempt_ps=1.0,
    platform="CUDA",
)
rest2.run(prepare(cfg), cfg)
analyze_rest2(cfg.outdir, reference_state=0)  # state 0 = λ=1, the unscaled physical ensemble
```

Each engine exposes the **same signature** — `run(eq, cfg) -> Path` — and `MDConfig` /
`REMDConfig` / `REST2Config` all subclass the shared `PrepConfig`, so `prepare()` depends
only on the base fields.

---

## The prepare/produce seam

Everything up to production is identical across the three protocols; only production
differs, so the engine splits exactly there:

```
prepare(cfg) ─────────────────────────────► EquilibratedSystem   (SHARED — ~90% of the code)
  pdbfixer prep · solvate · minimize ·        (System + positions + box, on disk)
  NVT equil · NPT density equil (adaptive)            │
  (protein heavy-atom restraints held full-            │
   strength through ALL of equilibration)    ┌─────────┼─────────┐
                                              ▼         ▼         ▼
                                        md.run     remd.run   rest2.run
                                        plain MD   T ladder   λ-scaled solute
                                        (NPT/NVT)  (NVT)      (NVT, all at T_0)
```

Restraints release only at the **start of production**, so the trajectory captures drift
from the input pose, not drift accumulated during equilibration. Density equilibration uses
an adaptive plateau (least-squares slope) test rather than a consecutive-segment difference.

---

## Force fields & water

Pick both deliberately and **pin them** — defaults can change as the codebase evolves. A
single `resolve()` enforces legal FF/water pairings (and warns on degraded ones like
ff19SB+TIP3P).

**Protein** (`ProteinFF`): `AMBER99SBILDN` (legacy), `AMBER14SB`, `AMBER19SB` (QM-trained,
recommended with OPC), `CHARMM36`, `CHARMM36M` (Huang-2017 CMAP; usually right for
folding/stability — and *different* from CHARMM36).

**Water** (`WaterModel`): `TIP3P`, `SPCE`, `TIP3PFB`, `OPC3` (3-site); `TIP4PEW`, `TIP4PFB`,
`OPC` (4-site). Packing box is chosen by **site count**; parameters come from the FF xml.

Recommended pairings: AMBER14SB→TIP3P, AMBER19SB→OPC, CHARMM36/36M→(CHARMM) TIP3P.

> **Note:** REST2 solute scaling currently supports **AMBER + CMAP**. CHARMM36m's
> `CustomNonbonded`/`CustomBond`/`CustomTorsion` forces are refused loudly (not silently
> mishandled) — see `knowledgebase/TODO.md`.

---

## Output layout

```
outdir/
  input.pdb                       # copy of the input pose (analysis reference)
  config.json, run_info.json      # pinned config + run metadata
  prepared/                       # the shared handoff: system.xml, state.xml, topology.cif, prep_report.json
  production/                     # plain MD: trajectory.xtc, final_state.xml, production.log, run_report.json
  remd/ (or rest2/)               # <protocol>.nc, ladder.json, run_report.json
    analysis/                     # exchange_report.json, demuxed state00_ensemble.xtc, plots, report JSON
  analysis/                       # plain MD: processed trajectory, analysis_report.json, data npz, RMSD/Rg/RMSF/DSSP PNGs
```

### REMD/REST2 demultiplexing — read this before trusting a replica-exchange trajectory

openmmtools stores **walkers** (it calls them "replicas") — continuous coordinate sets that
move *through* state space as swaps reassign their thermodynamic state. The `.nc` trajectories
are **not** fixed-temperature ensembles. To get a constant-T (REMD) or fixed-λ (REST2)
ensemble you **demux**: per iteration, pick the walker occupying the target slot. Only slot 0
(T_min for REMD, λ=1 for REST2) is a real Boltzmann sample at the design condition — the only
one that feeds RMSD/Rg/RMSF/clustering. `analyze_remd` / `analyze_rest2` do this for you (see
`analysis/multistate.py`, `demux_state`).

---

## Running on SLURM

Each demo ships a `.sbatch`. The pattern (see `examples/MD_demo/demo_short_run.sbatch`):

```bash
#SBATCH --gres=gpu:l40s:1
eval "$(conda shell.bash hook)"
conda activate openmm_env
cd "$SLURM_SUBMIT_DIR"       # sbatch copies the script to a temp path — BASH_SOURCE won't point at the repo
python ./demo_short_run.py
```

Submit from the demo folder: `cd examples/MD_demo && sbatch demo_short_run.sbatch`.

---

## Examples

Worked, ps-scale end-to-end smoke runs (with committed outputs) live in `examples/`:

- `examples/MD_demo/` — `prepare()` + `md.run()` + `analysis.analyze()`
- `examples/REMD_demo/` — T-REMD + `analyze_remd()` (shows why explicit solvent needs many rungs)
- `examples/REST2_demo/` — REST2 + `analyze_rest2()` (same reach, a handful of rungs)
- `examples/input_structures/` — shared test structures

---

## Tests

```bash
conda activate openmm_env
python -m pytest -q                 # full suite
python -m pytest -q -m "not slow"   # skip the real-OpenMM solvation/dynamics builds
```

The pure units (`forcefield.py`, `config.py`, density math) import no OpenMM and run on any
env; run the full suite in `openmm_env` so the solvation/dynamics modules are covered.
`tests/gpu_tests.sbatch` runs the suite on CUDA.

---

## Code style

**Fail loudly** — crashes are cheap, silent wrong results are expensive. No bare
`except`, no `.get(key, default)` or `x or fallback` to paper over uncertainty, `zip(strict=True)`,
dataclasses/pydantic over raw dicts, parse-don't-validate at boundaries. Details in `CLAUDE.md`.

---

## Docs

- `knowledgebase/CODEBASE_SHAPE.md` — architecture map of record (module-by-module)
- `knowledgebase/DECISIONS.md` — the *why* behind each design choice
- `knowledgebase/GOTCHAS.md` — discovered pitfalls (OpenMM/openmmtools/SLURM/CHARMM quirks)
- `knowledgebase/HISTORY.md` — dated session log
- `knowledgebase/TODO.md` — known gaps and planned features
- `CLAUDE.md` — working rules; `SKELETON_v2.md` — original design overview
