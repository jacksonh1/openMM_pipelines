# TODO — openmm_pipelines

Concise backlog. One terse bullet per item. Promote a multi-session effort to a
`plans/` file if one grows.

## Implementation (plain MD — per-engine layout; see CODEBASE_SHAPE.md)

- [ ] `pyproject.toml` — core deps (openmm, pdbfixer, mdtraj, numpy, pydantic), `[remd]` extra -> openmmtools.
- [ ] `forcefield.py` — `ProteinFF`/`WaterModel`/`NonbondedSpec`, `resolve()`, `ForceFieldError`, `_PACKING_BOX`, `_RECOMMENDED_WATER`.
- [ ] `test_forcefield.py` — resolve() legal/illegal pairings, packing-by-site-count.
- [ ] `core/config.py` — `PrepConfig` (base), `BoxShape`, validators (HMR assert, recommended-water warning), derived step counts.
- [ ] `test_config.py` — validators + derived counts (+ `md/config.py` `MDConfig`).
- [x] `core/density.py` — `density_converged` + `assess_plateau` (numpy OLS, `PlateauTest`). (2026-10-01)
- [x] `test_density.py` — plateau vs drift vs scatter + arg-guards (10 tests, pass). (2026-10-01)
- [ ] `core/equilibrated.py` — `BuiltSystem` + `EquilibratedSystem`; `to_disk` (system.xml + state.xml + topology PDB/CIF) + `from_disk(outdir)` **one-arg**.
- [ ] `core/build.py` — `build_system`, disulfide-preservation guard, `_assert_water_topology`.
- [ ] `core/restraints.py` — `add_restraint` (global `k`), `reanchor`.
- [ ] `core/stages.py` — `minimize` (finite-energy assert), `equilibrate_nvt`, `add_barostat`, `equilibrate_density`, `relax`.
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
- [ ] `rest2.py` — solute λ-scaling (the one genuinely involved piece); CMAP safety gate.
- [ ] REMD analysis — exchange acceptance, round-trip mixing (read the `.nc` reporter).
