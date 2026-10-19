# openmm_pipelines

A small, Pythonic OpenMM engine for **explicit-solvent MD of folded protein structures** on
single-node GPU clusters. Three protocols — **plain production MD**, **T-REMD** (temperature
replica exchange), and **REST2** (replica exchange with solute tempering) — all built behind
one shared *prepare → produce* seam. It is system-agnostic: anything OpenMM can parameterize
from a folded PDB works.

Every run starts from a **folded input structure**, never from an unfolded or extended state.
The input pose is the reference every analysis metric is measured against (RMSD is drift from
it, RMSF is local flexibility). This is a characterization tool, not a folding tool.

---

## Choosing a protocol

- **plain MD** — a single trajectory of the system as given; sample a specific state (e.g. a
  bound complex) or the local relaxation/flexibility around the input pose.
- **T-REMD / REST2** — enhanced conformational sampling of the input structure. In explicit
  solvent **REST2 is usually the better choice** — it tempers only the solute, reaching a
  given effective temperature with far fewer replicas than T-REMD needs for the same range.
  Reach for **T-REMD** when REST2 can't scale your force field: REST2 solute scaling supports
  AMBER + CMAP only, so a CHARMM36m run uses T-REMD (see *Force fields*).

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

**Hardware:** production runs expect a **CUDA GPU** (`platform="CUDA"`, the default). With no
GPU, set `platform="CPU"` — fine for a tiny smoke test, far too slow for a real run.

---

## Preparing your input structure

`prepare()` runs PDBFixer on `pdb_in` before solvation, so a reasonably clean PDB is enough
— you do **not** need to pre-protonate or pre-fix it. PDBFixer adds missing heavy atoms and
short missing loops, adds hydrogens at the configured `ph` (default 7.4), and the build
**preserves existing disulfides** (an SG–SG count change during fixing is a hard error). What
it does **not** do:

- **Only standard protein residues are recognized** (the 20 AAs + common protonation/termini
  variants). Ligands, cofactors, nucleic acids, or non-standard residues have no parameters
  and will fail at system build — strip or separately parameterize them first.
- It does not choose protonation states beyond the pH heuristic, assign a fold, or dock
  anything. A multi-chain system (e.g. a complex) must already be assembled in the input PDB
  in the pose you want to simulate.

---

## Quickstart

### Plain production MD

```python
from openmm_pipelines import MDConfig, ProteinFF, WaterModel, prepare, md, analysis

cfg = MDConfig(
    pdb_in="structure.pdb",
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
    pdb_in="structure.pdb", outdir="run_remd",
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
    pdb_in="structure.pdb", outdir="run_rest2",
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
`REMDConfig` / `REST2Config` all subclass the shared `PrepConfig`. The quickstarts show the
fields you'll usually set; **every field, its default, and its units is in the config class**
(`openmm_pipelines/config.py` for the shared ones, each engine's `config.py` for the rest).
Configs are frozen — pin the parameters you care about rather than relying on defaults,
which can change as the tool evolves.

Because `prepare()` writes the equilibrated system to disk, you can equilibrate once and
launch several productions from it without re-running prepare:

```python
from openmm_pipelines import EquilibratedSystem
eq = EquilibratedSystem.from_disk(cfg.outdir)   # reload a prior prepare()
md.run(eq, cfg)
```

**Key defaults** you inherit if you don't set them (full list + units in the config classes):
`ph=7.4`, `box_shape=dodecahedron`, `padding_nm=1.0`, `ionic_strength_molar=0.15`,
`neutralize=True`, `temperature_k=300`, `dt_ps=0.004` (4 fs, needs the default HMR
`hydrogen_mass_amu=4.0`; for 2 fs set `hydrogen_mass_amu=1.008, dt_ps=0.002`),
`friction_ps=1.0`, `restraint_k=1000` (kJ/mol/nm² on protein heavy atoms during
equilibration), `nvt_equil_ns=0.2`. Plain MD adds `total_ns=100`, `traj_ps=10`.

### Choosing a replica ladder (T-REMD / REST2)

For both exchange protocols you either give a **replica count** (`n_replicas`) and the code
builds a geometric ladder between the endpoints, or pin the full ladder explicitly
(`temperatures_k=(...)` / `lambdas=(...)`). Per-adjacent-rung exchange acceptance is reported
in `<protocol>/analysis/exchange_report.json`. `solute_chain_index` (REST2) is **0-based**
(`1` = the second chain; `None` = whole protein).

---

## How a run is structured (the prepare/produce seam)

Everything up to production — structure prep, solvation, minimization, equilibration — is
**identical across all three protocols**. Only production differs, so a run always has two
steps: `prepare(cfg)` builds and equilibrates the solvated system, then one engine's
`run(eq, cfg)` does production. This is why you pass the same `cfg` to both.

```
prepare(cfg) ─────────────────────────────► EquilibratedSystem   (shared by every protocol)
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

**Protein** (`ProteinFF`): `AMBER99SBILDN`, `AMBER14SB`, `AMBER19SB` (QM-trained, pairs with
OPC), `CHARMM36`, `CHARMM36M` (Huang-2017 CMAP; distinct from CHARMM36).

**Water** (`WaterModel`): `TIP3P`, `SPCE`, `TIP3PFB`, `OPC3` (3-site); `TIP4PEW`, `TIP4PFB`,
`OPC` (4-site). Packing box is chosen by **site count**; parameters come from the FF xml.

Pairings `resolve()` recommends (warns otherwise): AMBER14SB→TIP3P, AMBER19SB→OPC,
CHARMM36/36M→(CHARMM) TIP3P.

> **Note:** REST2 solute scaling currently supports **AMBER + CMAP** only; CHARMM36m's
> `CustomNonbonded`/`CustomBond`/`CustomTorsion` forces are refused (not silently mishandled).

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

The `analysis/` folder holds the computed metrics, all measured against `input.pdb`:
`analysis_report.json` (RMSD drift from the input pose, radius of gyration, per-residue RMSF,
secondary-structure fractions, conformational-cluster count), the matching PNG plots vs time,
and `analysis_data.npz` with the raw arrays.

### REMD/REST2 demultiplexing

openmmtools stores **walkers** (it calls them "replicas") — continuous coordinate sets that
move *through* state space as swaps reassign their thermodynamic state. The `.nc` trajectories
are **not** fixed-temperature ensembles. To get a constant-T (REMD) or fixed-λ (REST2)
ensemble you **demux**: per iteration, pick the walker occupying the target slot. Only slot 0
(T_min for REMD, λ=1 for REST2) is a real Boltzmann sample at the reference condition — the
only one that feeds RMSD/Rg/RMSF/clustering. `analyze_remd` / `analyze_rest2` do this for you
(see `analysis/multistate.py`, `demux_state`).

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
- `examples/REMD_demo/` — T-REMD + `analyze_remd()`
- `examples/REST2_demo/` — REST2 + `analyze_rest2()`
- `examples/input_structures/` — shared test structures

Throughput depends on system size and GPU; the demo SLURM logs (`examples/*/*.out`) record
real wall-times for these boxes.

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
