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
- Open: `pyproject.toml`, then `forcefield.py` + `test_forcefield.py` and
  `core/config.py` + `test_config.py` (pure, no GPU) next. Output contract to
  discuss before `equilibrated.py` / `prepare.py` / `md/production.py`.
