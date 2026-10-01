# HISTORY — openmm_pipelines

Dated session log, newest on top. One entry per working session/day. Terse and
factual.

---

## 2026-10-01 — knowledgebase seeded

- Started the package. Repo currently holds only `CLAUDE.md`, `README.md`,
  `SKELETON_v2.md` (design sketch), and `archive/SKELETON.md` (prior sketch). **No
  `.py` files exist yet.**
- Read `SKELETON_v2.md` end to end and seeded `knowledgebase/`: `README.md`
  (maintenance protocol), `CODEBASE_SHAPE.md` (intended architecture), `DECISIONS.md`
  (the settled design decisions extracted from the sketch), `GOTCHAS.md`, `TODO.md`.
- Design intent as of today: plain-MD engine organized around the **prepare/produce
  seam** so T-REMD and REST2 drop in later as new `produce_*` backends with zero
  changes to the shared equilibration code. Dependencies kept minimal (core OpenMM +
  pdbfixer + mdtraj; `openmmtools` an optional `[remd]` extra).
- Implemented `density.py` first (pure, no GPU): `assess_plateau` + `density_converged`.
  Ported the *science* from the GROMACS engine's `scripts/simulation/density_converged.py`
  (trailing-window slope test) but **not** its packaging — the GROMACS version is a
  stdlib CLI script with hand-rolled OLS and stdout verdict-parsing, built to run under
  an sbatch's system python. The OpenMM version is an in-process library function using
  `numpy.polyfit` for the fit (per the global "don't hand-roll a standard algorithm"
  rule). Added `tests/test_density.py` (10 tests, all pass under `fragforge` env, numpy
  1.26.4) — the key case asserts a 0.4 %/segment drift that *passes* a pairwise
  difference still *fails* the slope test. Minimal `openmm_pipelines/__init__.py` stub
  (full exports land as their modules do). Fixed the stale GROMACS path in the CLAUDE.md
  density gotcha.
- **Layout reorg + naming (settled).** Moved from SKELETON's flat module list to
  per-engine folders: top level is just `__init__.py` + `forcefield.py`; shared
  equilibration in `core/`; `md/` `remd/` `rest2/` each self-contained
  (`config.py` + `production.py`). Seam verbs are `prepare` / `production`; dropped
  `produce_*`. Each engine exposes one `run(eq, cfg) -> Path` — no `cfg`-only
  wrapper (production always consumes an `EquilibratedSystem`; one-shot is explicit
  `md.run(prepare(cfg), cfg)`). `forcefield.py` split out of `config.py`. Relocated
  the built `density.py` -> `core/density.py`; tests still pass. See `DECISIONS.md`.
- **Validated the design against OpenMM docs** (User Guide, cookbook,
  openmmtools): seam placement, the `EquilibratedSystem` handoff, the global-`k`
  restraint idiom, 4 fs HMR, and that REST2 has no turnkey helper (manual custom
  forces + CMAP gap) all check out. See `DECISIONS.md`.
- **Topology output fix:** `prepare()` will write a topology file so
  `EquilibratedSystem.from_disk(outdir)` is one-arg (System xml lacks topology);
  also needed for analysis. Full output contract deferred to a later discussion.
- **Built the pure (openmm-free) units:** `pyproject.toml` (core deps + `[remd]`/
  `[test]` extras), `forcefield.py` (`ProteinFF`/`WaterModel`/`NonbondedSpec`/
  `resolve` + registries, split out of config per the layout decision),
  `core/config.py` (`PrepConfig` + `BoxShape`, pydantic v2 frozen, validators,
  derived step counts), `md/config.py` (`MDConfig`). Wired top-level `__init__.py`
  to the live surface. Tests: `test_forcefield.py` (10), `test_config.py` (13) —
  **33 pass total** under `fragforge` env (pydantic 2.12.5, py 3.12). These import
  no OpenMM, so config/FF validation runs on any env.
- **Env settled: `openmm_env`** (this is an isolated tool with its own env — do not
  use the sibling repos' envs). Has openmm 8.6.1, pdbfixer, mdtraj, numpy, scipy,
  openmmtools 0.26.0; added pydantic 2.13.5 + pytest 9.1.1; installed the package
  editable. Full suite (33) passes there. Recorded in the project `CLAUDE.md`
  Environment section.
- **Track A (compute core) built + tested in `openmm_env`:** `core/restraints.py`
  (cookbook global-`k` idiom; 6 tests on a trivial System), `core/build.py`
  (`build_system` + disulfide guard + water-topology guard; 7 tests on
  `helix_fusion.pdb`, incl. a real ff19SB+OPC 4-site/virtual-site build confirming
  packing-by-site-count), `core/stages.py` (`minimize`/`equilibrate_nvt`/
  `add_barostat`/`equilibrate_density`/`relax`; 5 `slow` tests). Added
  `core/handoff.py` (renamed from `equilibrated.py`) with the `BuiltSystem` +
  `EquilibratedSystem` dataclasses; `to_disk`/`from_disk` stubbed (NotImplementedError)
  pending the output contract. Non-slow suite: **45 pass**.
- **GPU test harness:** `tests/gpu_tests.sbatch` (pi_keating, l40s, `openmm_env`),
  platform env-driven via `OPENMM_TEST_PLATFORM` (default CPU; CUDA on GPU). CPU
  dynamics tests were too slow interactively; run them with
  `sbatch -w node3620 tests/gpu_tests.sbatch`.
- Open: **output contract discussion** before `core/prepare.py` / `md/production.py`
  and before implementing `handoff.to_disk/from_disk`. Still to write:
  `protein_heavy_atoms` selection (for `prepare`), `analysis/`.

### Naming / layout working session (later, 2026-10-01)

Reworked the module layout with the user into four tiers (see `DECISIONS.md`):
- `core/` dissolved: reusable toolkit → **`lib/`** (`restraints.py`,
  `density_convergence.py`); the shared pipeline → **`preparation/`** (composer +
  stages + `handoff.py`).
- **`config.py` moved to the package root** (base config for every engine, not part
  of the pipeline).
- Renames for honesty: `density.py` → `lib/density_convergence.py`; `stages.py` →
  `preparation/equilibration_steps.py` (now *only* nvt + barostat + density);
  `minimize` split into its own `preparation/minimize.py` (static, not
  equilibration). `add_barostat` stays with the density step (trivial glue, not
  toolkit).
- **`relax` removed** (function + config fields) — it contradicts the drift-from-
  input analysis contract and isn't needed; enhanced-sampling pre-equilibration is a
  per-replica REMD/REST2 backend concern (TODO), not shared `prepare()`.
- Full suite green after the refactor: **45 non-slow local, 51 on GPU** (node3620).
- Open unchanged: output contract before `prepare.py` / `md/production.py`.

### Plain MD wired end to end + first run (later, 2026-10-01)

- **Output contract (FIRST DRAFT)** recorded in `DECISIONS.md` and reviewed against
  OpenMM docs: `outdir/{config.json, run_info.json, input.pdb, prepared/{system.xml,
  state.xml, topology.cif, prep_report.json}, production/{trajectory.xtc,
  production.log, final_state.xml, production.chk, run_report.json}}`. `state.xml` =
  portable `saveState`; `.chk` = exact same-hardware restart. Typed `RunInfo`/
  `PrepReport`/`RunReport`. No input hash (not content-addressed).
- Built: `lib/selections.py` (`protein_heavy_atoms`), `reports.py`,
  `preparation/handoff.py` (to_disk/from_disk), `preparation/prepare.py` (composer +
  QC), `md/production.py` (`run`). `equilibrate_density` now returns per-segment
  volumes.
- **Bug caught before running:** a new `Context` resets global params to the force
  *default*, so production would re-enable restraints — fixed by
  `setGlobalParameterDefaultValue(0, 0.0)` at release. Logged in `GOTCHAS.md` (+ a
  second gotcha: checkpoint non-portability).
- **End-to-end demo** (`examples/demo_short_run.py`, GPU node3620, helix_fusion,
  ~20 ps): prepare + production succeeded, all output files written. QC sane —
  minimize PE −19.9k→−86.2k kJ/mol; NVT T 295.7 K; density converged in 5 segments
  (drift 0.9% < 2% tol); production T 300.9±2.9 K, **density 1.02±0.005 g/mL**, PE
  stable, totalEnergy drift small, `stable=true`; ~1455 ns/day on an L40S. **51 unit
  tests pass on GPU.**
- Open: `analysis/` (trajectory load + PBC image + RMSD/Rg/RMSF/DSSP vs input pose);
  output contract refinement; then REMD/REST2 backends. `relax`→per-replica
  equilibration still noted for REMD.

### analysis/ subpackage built (later, 2026-10-01)

- Built `analysis/`: `trajectory` (load xtc with cif topology parsed via OpenMM
  `PDBxFile` → `md.Topology.from_openmm`, since mdtraj's cif reader is unreliable;
  `image_molecules` → strip to protein → superpose onto the input pose), `metrics`
  (RMSD-to-design, Rg, per-residue RMSF, all in Å — RMSD maps Cα across the two
  topologies via `ref_atom_indices`; RMSF is about the trajectory mean, no external
  ref), `dssp` (mdtraj native `compute_dssp`, no mkdssp), `clustering` (scipy
  average-linkage on pairwise Cα RMSD), `plots` (matplotlib Agg PNGs), and an
  `analyze(outdir)` composer writing `outdir/analysis/` (report JSON + npz + 4 PNGs).
- Added **matplotlib** (+ confirmed pandas) to `environment.yml` and `pyproject.toml`;
  installed matplotlib-base via conda-forge.
- Validated `analyze()` on the demo output (helix_fusion, 20 ps, 43 res): **79.5%
  helix**, RMSD-to-design mean 0.73 Å / max 1.17 Å, Rg 10.6 Å, peak RMSF at residue
  26, 1 cluster — all sane. Tests: `test_analysis` (5 synthetic + 1 guarded demo
  smoke). **50 non-slow tests pass.** Wired `analyze()` into the demo driver.
- Open: output-contract refinement; REMD/REST2 backends (+ `analysis/remd.py` for
  exchange/round-trip metrics).
- Docs: removed the `validated` flag (YAGNI) + added CHARMM potential-vs-force-switch
  provenance/caveat (verified vs OpenMM/ParmEd docs); force-field energy-match tool
  deferred to TODO with the no-install (OpenMM-only) design recorded. Rewrote
  `SKELETON_v2.md` to the as-built state (implemented, trimmed — no duplicated code)
  and reduced `README.md` to a barebones install + 5-line usage + docs pointer.
