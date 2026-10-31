# TODO — openmm_pipelines

Concise backlog. One terse bullet per item. Promote a multi-session effort to a
`plans/` file if one grows.

## Implementation (plain MD — per-engine layout; see CODEBASE_SHAPE.md)

- [x] `pyproject.toml` — core deps + `[remd]`/`[test]` extras. (2026-10-01)
- [x] `forcefield.py` — `ProteinFF`/`WaterModel`/`NonbondedSpec`, `resolve()`, `ForceFieldError`, registries. (2026-10-01)
- [x] `test_forcefield.py` — pairings, packing-by-site-count, exhaustiveness (10 tests). (2026-10-01)
- [x] `core/config.py` — `PrepConfig` + `BoxShape`, validators, derived counts; `md/config.py` `MDConfig`. (2026-10-01)
- [x] `test_config.py` — validators + derived counts + MDConfig (13 tests). (2026-10-01)
- [x] `core/density.py` — `density_converged` + `assess_plateau` (numpy OLS, `PlateauTest`). (2026-10-01)
- [x] `test_density.py` — plateau vs drift vs scatter + arg-guards (10 tests, pass). (2026-10-01)
- [x] `preparation/handoff.py` — `BuiltSystem` + `EquilibratedSystem`; `to_disk`/`from_disk(prepared_dir)` implemented (portable System+State XML + mmCIF topology). (2026-10-01)
- [x] `lib/selections.py` — `protein_heavy_atoms`. (2026-10-01)
- [x] `reports.py` — `RunInfo` / `PrepReport` / `RunReport` dataclasses + `write_json`. (2026-10-01)
- [x] `preparation/prepare.py` — composer + QC capture + output-contract writes. (2026-10-01)
- [x] `md/production.py` — `run(eq, cfg)`; XTC/StateData/Checkpoint reporters + run report. (2026-10-01)
- [x] End-to-end demo on `helix_fusion` (GPU): prepare + production green, QC sane. (2026-10-01)
- [x] `preparation/build.py` — `build_system`, disulfide guard, `_assert_water_topology`. (2026-10-01)
- [x] `lib/restraints.py` — `add_restraint` (global `k`), `reanchor`. (2026-10-01)
- [x] `lib/density_convergence.py` — `assess_plateau`, `density_converged` (moved from density.py). (2026-10-01)
- [x] `preparation/minimize.py` — `minimize` (static; finite-energy guard). (2026-10-01)
- [x] `preparation/equilibration_steps.py` — `equilibrate_nvt`, `add_barostat`, `equilibrate_density`. `relax` dropped. (2026-10-01)
- [x] tests: `test_restraints` (6), `test_build` (7), `test_equilibration` (5, `slow`); `tests/gpu_tests.sbatch` runs the suite on CUDA (`-w node3620`). 51 pass on GPU. (2026-10-01)
- [ ] `core/prepare.py` — shared pipeline; `protein_heavy_atoms` selection (assert non-empty); writes topology.
- [ ] `md/config.py` — `MDConfig(PrepConfig)`.
- [ ] `md/production.py` — `run(eq, cfg)`; XTC/StateData/Checkpoint reporters; `export_final_pdb`.
- [x] `analysis/` — `trajectory` (load + PBC image + strip + superpose to input), `metrics` (RMSD/Rg/RMSF, Å), `dssp`, `clustering` (scipy), `plots`, `analyze(outdir)` composer. Validated on demo_out (5 synthetic + 1 smoke test). (2026-10-01)
- [ ] `__init__.py` top-level + `core/`/`md/` re-exports; `README.md` usage.

## Analysis output upgrades (2026-10-21; see DECISIONS/HISTORY)

- [x] `analysis/config.py` — `AnalysisConfig` (pure pydantic), serialized into every report; composers take `config`.
- [x] Clustering overhaul — `ClusterResult`, hierarchical default (auto-stride, recorded) + optional DBSCAN (`[cluster]` extra, scikit-learn), representative PDBs + populations + plots, design-anchored (rep RMSD-to-design, design-like fraction). Fixed noise-in-`n_clusters` + frame-index-through-stride.
- [x] `analysis/interchain.py` — min distance / contact count / COM distance (min-image) + per-chain RMSD; auto-run at ≥2 chains; positional chain correspondence; computed pre-superposition.
- [x] `build._assert_no_internal_gaps` — prep-phase chain-integrity guard (so analysis trusts chain counts).
- [x] `lib/run_layout.py` `RunLayout` + `storage_directory` on `PrepConfig` — offload/record/mirror/clean, deterministic collision-proof folder, O_EXCL marker; wired into md/remd/rest2 production + analysis composers.
- [x] `analysis/report.py` — self-contained interactive Plotly `report.html` (optional `[report]` extra), on by default (`AnalysisConfig.html_report`), mirrored into storage.
- [x] Per-residue DSSP timeline (`dssp.dssp_timeline` + `plots.plot_dssp_timeline`) → `analysis/dssp_timeline.png` + an interactive report heatmap; figure size bounded on both axes (GROMACS-style). (mdtraj DSSP lacks DSSP-4 PP-II/`=`; mkdssp-4 backend is a possible follow-on.)
- [x] Output layout: clustering artifacts → `analysis/clustering/`; multistate analysis → top-level `outdir/analysis/` (consistent with plain MD).
- [ ] `storage status`/gc helper — list orphaned storage folders (marker points to a gone/relinked outdir). Marker enables it; not yet built.
- [ ] Optional: content/config hash recorded for dedup (run identity stays outdir; layout unchanged).

## Force-field corrections

- [x] **Cation-π NBFIX (Liu et al. 2020)** — `cation_pi.py` + `PrepConfig.cation_pi_nbfix` (opt-in,
  CHARMM36m-only); Arg `NC2` + Lys `NH3` × aromatic `CA`/`CAI`/`CPT`; merged in `build_system`;
  REST2 refuses it. (2026-10-29) Runtime wiring validated by vacuum-energy ΔE controls
  (positive 1a22 = −20.8 kJ/mol attractive; negative poly-Ala = unchanged). (2026-10-30)
- [ ] **Cation-π CGenFF twin (supp Table S1)** for non-natural / CGenFF-typed residues
  (`NG2P1`/`NG3P0`/`CG334` × `CG2R61`/`CG2RC0`, incl. TMA⁺ for trimethyl-Lys/choline). Drop-in row
  append to `_CATION_PI_NBFIX`; deferred until an ncAA is actually in use. Note the paper leaves
  mixed natural-cation × CGenFF-aromatic pairs and partially-methylated Lys untabulated.

## Science-review items to resolve before trusting production

- [ ] Restraint `x0` (absolute) vs MC barostat COM rescaling — confirm benign.
- [x] PBC wrap on trajectory (`image_molecules` in analysis). Final-structure wrap still TBD.
- [ ] **Force-field energy-match validation tool** (deferred — decided 2026-10-01). A real, runnable per-force-group energy check to run once per FF before a big campaign; replaces the removed fake `validated` bool. Decisions reached:
  - **AMBER (the default) is low-risk** — standard, switch-free, well-trodden; no validation needed. CHARMM is the one to check.
  - **The real risk is GROSS errors** (missing/wrong switch, cutoff, method, FF file), not the subtle potential-switch-vs-force-switch delta (small, widely-accepted approximation; matters for free-energy, not stability/drift MD).
  - **Preferred no-install design (tier 1):** fully-local, OpenMM-only — build one `ForceField('charmm36*.xml')` system and compare per-force-group energies under our built-in *potential*-switch vs a `CustomNonbondedForce` *force*-switch (CHARMM-GUI's published formula). Turns the switch question into a kJ/mol number, no CHARMM/GROMACS install. Likely `scripts/validate_charmm_switch.py`.
  - **Tier 2 (gross-error cross-check):** `ForceField(xml)` route vs `CharmmPsfFile`+`CharmmParameterSet` route, both in OpenMM (needs a `.psf` from CHARMM-GUI web / VMD psfgen, not the CHARMM program).
  - NOTE: choderalab/OpenMMEnergyComparisons ships *code to generate* benchmarks, not cached reference energies — it requires installing the other engines, so it is NOT a no-install path.
  - Docs: method https://docs.openmm.org/latest/userguide/library/07_testing_validation.html ; per-force-group how-to https://openmm.github.io/openmm-cookbook/dev/notebooks/analysis_inspection/Analyzing%20Energy%20Contributions.html ; `parmed.openmm.energy_decomposition` https://parmed.github.io/ParmEd/html/openmmobj/parmed.openmm.energy_decomposition.html

## Open design decisions

- [ ] **Independent replicates & seeding** — `seed` drives reproducible replicates for plain MD only; it is **inert in REMD/REST2 production** (openmmtools RNG unseeded), and there is no replicate-management layer. No decision yet — see `replicates_and_seeding.md` (findings + best-practice standard + open questions).
- [~] **Output contract** — layout stable (prepare / production / remd|rest2 / analysis subdirs); run identity = **outdir** (path-based, decided 2026-10-21, not content-addressed); `storage_directory` offload/mirror added via `RunLayout`. Exact file naming/formats still open for a dedicated pass.
- [x] REMD production ensemble default: **NVT** (decided 2026-10-06; see DECISIONS).
- [x] REMD/REST2 backend home: per-engine folders; `[remd]` extra gates the openmmtools import and `remd/`+`analysis/remd.py` are not eagerly imported. (2026-10-06)
- [x] `analysis/remd.py` single module (not nested) — kept one module. (2026-10-06)

## Future backends (designed-for, not built)

- [x] `remd.py` + `ladder.py` — T-REMD via openmmtools multistate `ReplicaExchangeSampler`, NVT, geometric ladder. (2026-10-06)
- [x] **Per-replica equilibration** in REMD — `sampler.equilibrate()` (discarded) at each replica's own T; `REMDConfig.equilibration_ns`. (2026-10-06)
- [x] REMD analysis — exchange acceptance, round-trip mixing, de-multiplex to fixed-T ensembles (`analysis/remd.py`, reads the `.nc`). (2026-10-06)
- [x] **REST2 backend** — `rest2/` (scaling + ladder + config + production) + `analysis/rest2.py`, path A. Solute scaling verified on a real system. (2026-10-06; plan `plans/rest2-kickoff.md`)
- [x] REST2 production ensemble default: **NVT** (by analogy to REMD; barostat stripped). (2026-10-06)
- [x] `rest2/scaling.py` — solute λ-scaling (√λ charge, λ·ε, λ bonded) + exception rescale + CMAP safety gate + released-restraint gate. (2026-10-06)
- [ ] **CHARMM REST2** — `scale_solute` refuses CHARMM's `Custom{Nonbonded,Bond,Torsion}Force` (NBFIX LJ / impropers). Implement their scaling (and tie to the FF energy-match tool) before any CHARMM REST2 campaign. AMBER works today.
- [ ] REST2: 1-4 (scaled exception) handling under charge/ε scaling is done via combining rules — cross-check against a per-force-group energy decomposition on a real system (open item from the kickoff note).
- [ ] REMD/REST2: expose demux of **all** rungs (not just the reference state) if multi-state ensembles are wanted; currently `demux_state(k)` per call, `analyze_*` does the reference rung. (shared `analysis/multistate.py`)
- [ ] REMD: per-neighbour-acceptance-driven ladder auto-tuning (iterate `n_replicas` from diagnostics) — currently manual.
- [x] **Resume / extend** — `resume(outdir, extend_to_ns=...)` continues (or extends) an interrupted run via the documented OpenMM/openmmtools restart mechanisms. MD: ping-pong portable State + truncate-to-checkpoint + `append=True`. REMD/REST2: openmmtools `from_storage` + `run`/`extend`. Manifest-driven engine dispatch; library-only (no SLURM coupling). (2026-10-30; DECISIONS/GOTCHAS; `tests/test_resume.py`)
- [ ] Resume: exact same-hardware binary `.chk` fast path, gated by an environment fingerprint — portable State is the default + sufficient; this only saves the RNG-stream divergence on a same-node restart. Deferred 2026-10-30.
- [ ] Resume: a `should_stop` hook on the run loops + a non-exception `is_complete(outdir)` — only if a self-requeuing SLURM driver is adopted (MD's monolithic `sim.step` would need chunking; REMD/REST2 are already iteration-chunked). Deferred 2026-10-30.

## Analysis — clustering

- [x] Clustering over-fragmentation fix — DBSCAN default + adaptive `min_samples`, HDBSCAN opt-in, explicit `cluster_stride` (no silent 2000-frame cap), hierarchical fails loud over `MAX_HIERARCHICAL_FRAMES`, sklearn → core dep. (2026-10-28; DECISIONS)
- [ ] `design_like_fraction` is read off cluster populations, so a genuine design-like basin holding <1.5% of frames collapses into DBSCAN noise and under-reports retention. Consider computing retention directly as the fraction of frames within `cluster_cutoff_ang` RMSD of the design pose, independent of clustering.
