# HISTORY — openmm_pipelines

Dated session log, newest on top. One entry per working session/day. Terse and
factual.

---

## 2026-10-07 — REST2 review pass: decompose/parameterize cleanup + re-entrancy fix

Reviewed the REST2 work against the project principles (simple single-purpose primitives;
decompose eagerly, parameterize lazily; don't hand-roll standard algorithms):
- **Removed a dead parameter** — `cluster_cutoff_ang` was threaded through
  `analyze_multistate`/`analyze_remd`/`analyze_rest2` but never used there (pre-existing in
  remd; I'd propagated it). Dropped from that chain (`analysis.analyze` still uses it).
- **De-duplicated neighbour-acceptance** — identical `_mean_neighbor_acceptance` in both
  `remd`/`rest2` production + the same per-iteration-sum in `analysis/multistate` → one
  primitive `exchange.neighbor_acceptance(reporter)` (the `read_mixing_statistics` axis-0 sum
  gotcha now lives in exactly one place).
- **Inlined** the needless `effective_temperature_ladder` one-line wrapper.
- **Tested the solute-chain knob** — kept `REST2Config.solute_chain_index` (core to binder
  campaigns, per the kickoff note — not speculative) but it was unexercised; added
  `test_selections.py` (whole protein / per-chain / H-included / water-excluded).
- **Fixed `run()` re-entrancy** (latent bug the 2nd multistate slow test exposed):
  openmmtools' global context cache refuses `set_platform` once populated, so a 2nd `run()`
  in one process crashed (also broke `python -m pytest -q`). Factored the platform setup into
  `context_cache.configure_global_platform`, which `empty()`s the cache first. New GOTCHAS
  entry. No new dependencies (reviewed: `np.geomspace` + numpy sums, nothing hand-rolled).
- Suite: 84 fast + 5 slow (both multistate tests now coexist in one process).
- **Added conformational clustering** of the demuxed slot-0 ensemble (T_min for REMD, λ=1
  for REST2) to `analyze_multistate` — `cluster_cutoff_ang` is now live (Cα-RMSD average-link
  via `analysis/clustering`, mirroring plain-MD `analyze`); `n_clusters`/`cluster_cutoff_ang`
  in the summary, `cluster_labels` in the npz. User wants clusters of the lowest slot, which
  is exactly the demuxed reference ensemble (a real Boltzmann sample at the design condition —
  the only slot that is, since the raw `.nc` stores walkers that move through state space).
- **Demuxed `state00_ensemble.xtc` is now fully processed**: PBC-imaged (demux) + desolvated +
  aligned-to-design + **origin-centered** (frame-0 centroid → origin), matching plain-MD's
  processed copy. Centering is translation-only; `md.rmsd` superposes internally and Rg/RMSF/
  clustering are translation-invariant, so metrics are unaffected. Folded into the one ensemble
  file (not a separate one) — the "keep raw" concern is the primary sim output (the `.nc`),
  which analysis only ever opens `open_mode="r"`.

---

## 2026-10-06 — REST2 backend built end to end

- Built the **REST2** engine behind the prepare/produce seam (`rest2/`) + its analysis,
  following `plans/rest2-kickoff.md` path A (N pre-scaled Systems, one plain
  `ThermodynamicState` at the shared physical T_0 per replica). `prepare()` reused
  unchanged.
- **`rest2/scaling.py::scale_solute(system, solute_indices, lam)`** — the one genuinely new
  piece. Returns an independent scaled copy: solute charge→√λ·q, LJ ε→λ·ε (σ fixed),
  solute-only bonded k→λ·k, CMAP maps→λ. NonbondedForce **exceptions** rescaled explicitly
  (OpenMM stores combined chargeProd/ε; solute-solute→λ, cross→√λ via the same combining
  rules). λ=1 is a verified energy identity.
- **Fail-loud force-field support.** Handles the AMBER force set + CMAP. **Refuses**
  (`NotImplementedError`) any unrecognized force with solute atoms — so CHARMM's
  `Custom{Nonbonded,Bond,Torsion}Force` (NBFIX LJ / impropers) is NOT silently mis-scaled;
  CHARMM REST2 is deferred (see TODO). The pipeline's released position restraint
  (`CustomExternalForce`, k=0) is carried through but **asserted released** first — an
  active restraint under REST2 would bias every replica.
- **λ ladder** (`rest2/ladder.py`): geometric in *effective* temperature [T_0, T_max_eff]
  (reuses `remd.ladder.geometric_ladder`), then λ_m = T_0/T_m. State 0 = λ=1 = the unscaled
  physical ensemble. Config `REST2Config` has `max_effective_temperature_k` / `n_replicas`
  (or explicit `lambdas`) + `solute_chain_index` (temper a single chain for binder work).
- **Shared base `MultiStateProductionConfig(PrepConfig)`** in `config.py` now holds the
  production/exchange fields + derived iteration-count properties; both `REMDConfig` and
  `REST2Config` subclass it (killed the duplication). `remd/production._nvt_system` factored
  into `lib/system_edits.strip_barostats` (+ `clone_system`), reused by rest2.
- **Analysis factored**: `analysis/multistate.py` holds the shared exchange-diagnostics +
  de-multiplex core; `analysis/remd.py` and new `analysis/rest2.py` are thin wrappers pinning
  the subdir + `.nc` name. Demux of state 0 gives the fixed-λ (=1) reference ensemble.
- **Solute scaling verified in a real run** (the user's explicit ask): `test_rest2.py`
  prepares a real solvated AMBER system and asserts a protein atom gets √λ·q / λ·ε while a
  water atom is byte-identical — not just inferred. Full suite green: 93 fast + 3 rest2 slow
  + 2 remd slow (CPU). New gotcha recorded: a Force proxy from a temporary System reads freed
  memory (caused garbage/`std::bad_alloc` in the first scaling tests).
- Added `examples/REST2_demo/demo_rest2.py` — 6 rungs over an effective 300→450 K (contrast
  the 14-rung T-REMD demo over 300–360 K; that contrast is the showcase).

---

## 2026-10-06 — T-REMD backend built end to end

- Built the **T-REMD** engine behind the prepare/produce seam (`remd/`) + its analysis
  (`analysis/remd.py`). Zero changes to `lib/`/`preparation/` — new engine folder +
  `REMDConfig(PrepConfig)` only, as the seam design promised. `prepare()` was reused
  unchanged (REMDConfig IS-A PrepConfig).
- **Production is NVT** (user decision 2026-10-06): equilibrated box already at target
  density, replicas hold volume fixed, so `remd.production._nvt_system` strips the
  `MonteCarloBarostat` from a serialized copy and rungs differ only in temperature.
- Backend: `openmmtools.multistate.ReplicaExchangeSampler` + geometric ladder
  (`remd/ladder.py`, `np.geomspace`) + `LangevinDynamicsMove(reassign_velocities=True)`.
  Per-replica equilibration via `sampler.equilibrate()` (discarded) — the replacement
  for the dropped shared-reference `relax`. Platform set on
  `cache.global_context_cache` (openmmtools reads it there, not per-call).
- **Analysis resolves the "single coordinate set vs fixed-T ensemble" question**
  (the CLAUDE.md critical concept): `exchange_diagnostics` reports per-neighbour
  acceptance + the replica->state permutation, plus mixing from openmmtools'
  `ReplicaExchangeAnalyzer.generate_mixing_statistics` (transition-matrix subdominant
  eigenvalue + state statistical inefficiency — the library's canonical mixing metric);
  `demux_state(k)` rebuilds the fixed-temperature NVT ensemble for rung k by picking, at
  each checkpointed iteration, the replica then occupying state k. Only the demuxed T_min
  ensemble feeds RMSD/Rg/RMSF/SS (vs the input design pose).
- **Principle review (same session).** (1) Dropped a hand-rolled round-trip counter in
  favour of openmmtools' `generate_mixing_statistics` (don't hand-roll a standard
  algorithm). (2) Extracted `analysis/drift.py::drift_report` — the RMSD/Rg/RMSF/SS +
  4-plots + summary block shared by `analyze` (plain MD) and `analyze_remd`; both now
  compose it (was duplicated). Fixed a mislabel: demuxed-ensemble x-axis is "demuxed
  frame", not "time (ps)" (`plot_timeseries`/`plot_ss_fractions` took an `xlabel` param).
- **GPU showcase + examples reorg.** Split `examples/` into per-engine `MD_demo/` +
  `REMD_demo/` (shared `input_structures/`); sbatch partition → `mit_normal_gpu`. First
  REMD demo (4 rungs, 300–360 K, 5.4k-atom box) showed **zero exchange** — correctly: it
  is the explicit-solvent √DOF scaling wall (acceptance ~1e-7 at ΔT≈19 K), and the QC
  (mean acceptance 0, subdominant eigenvalue 1.0, statistical inefficiency ≈ n_iter)
  flagged it. Re-ran with 14 rungs (~6 K gaps): **mean neighbour acceptance 0.33**, all
  pairs 0.23–0.41 (uniform → geometric ladder validated), eigenvalue 1.0→0.985. Engine +
  QC + example all validated on GPU.
- **REST2 handoff written** — `plans/rest2-kickoff.md`: what carries over from T-REMD, the
  solute λ-scaling math (charges·√λ, ε·λ, bonded+CMAP·λ; CMAP safety gate), the two
  implementation paths (pre-scaled Systems first; `CompoundThermodynamicState`+
  `GlobalParameterState` later — openmmtools ships no turnkey REST, `alchemy` is for
  free-energy not REST), config/ladder, open decisions, first steps. REST2 to be built in
  a fresh session starting from that note.
- Output contract: `outdir/remd/` → remd.nc (+ checkpoint), ladder.json, run_report.json
  (`REMDRunReport`: mean neighbour acceptance etc.); `outdir/remd/analysis/` →
  exchange_report.json, acceptance.png, replica_walk.png, demuxed ensemble xtc/pdb,
  drift report + npz.
- `remd/` and `analysis/remd.py` are **not** imported by the top-level/`analysis`
  `__init__` (they need the `[remd]` extra); import them explicitly. `config.py` +
  `ladder.py` stay OpenMM-free.
- Tests: `test_remd_ladder.py` (3, pure), `test_remd_config.py` (7, pure),
  `test_remd.py` (2, slow: prepare+run+analyze+demux on helix_fusion, CPU). Full suite
  now **60 non-slow + 10 slow**; all green (slow on CPU). Two new gotchas found
  (`read_mixing_statistics` per-iteration shape; mdtraj float32 box) — see GOTCHAS.

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
