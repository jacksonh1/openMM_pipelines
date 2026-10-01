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
- [ ] `analysis/` — `trajectory` (load + PBC image + strip + align), `metrics`, `dssp`, `clustering`, `plots`.
- [ ] `__init__.py` top-level + `core/`/`md/` re-exports; `README.md` usage.

## Science-review items to resolve before trusting production

- [ ] Restraint `x0` (absolute) vs MC barostat COM rescaling — confirm benign.
- [ ] PBC wrap on trajectory + final PDB export (`mdtraj image_molecules`).
- [ ] CHARMM energy-match before flipping any CHARMM `NonbondedSpec.validated`.

## Open design decisions

- [ ] **Output contract** — what files `prepare()` / `run()` write, layout, naming. Topology output settled (needed for reload + analysis); the rest deferred to a dedicated discussion.
- [ ] REMD production ensemble default: NVT or NPT?
- [ ] REMD/REST2 backend home: now per-engine folders inside this package; confirm `[remd]` extra gates the openmmtools import.
- [ ] `analysis/remd.py` single module vs nested `analysis/remd/` — keep one module until it grows.

## Future backends (designed-for, not built)

- [ ] `remd.py` + `ladder.py` — T-REMD via openmmtools multistate.
- [ ] **Per-replica equilibration** in REMD/REST2 backends (replaces the dropped `relax`) — each replica equilibrates at its own T/λ before collecting, via the multistate sampler (e.g. `REMDConfig.equilibration_iterations` or run-and-discard), NOT a shared reference-T relax.
- [ ] `rest2.py` — solute λ-scaling (the one genuinely involved piece); CMAP safety gate.
- [ ] REMD analysis — exchange acceptance, round-trip mixing (read the `.nc` reporter).
