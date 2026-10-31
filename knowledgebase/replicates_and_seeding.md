# Note: independent replicates & RNG seeding (open — no decision yet)

Status: **investigation note, 2026-10-28.** Captures (a) the field's best-practice
standard for replicate simulations and (b) what the codebase currently does with
`seed`. No design decision has been made — do not treat this as a decision record.

## Terminology (these get conflated)

- **Replica** — a parallel copy inside *one* REMD/REST2 run (the T / λ ladder).
  Replicas exchange, so they are **not** statistically independent; one REMD run is
  still one sample of the sampling process.
- **Replicate** — an independent *repeat* of the whole simulation from different
  initial velocities (different RNG), used to estimate statistical uncertainty.

## Industry best-practice standard (for literature review / later decision)

- **Never report from a single trajectory.** Run **≥3 independent replicates**;
  **3–5** is the common working range, more for quantitative observables
  (Δstability between variants, free energies).
- **Error bars = standard error across replicates**, not within-run block averaging
  alone (within-run blocks underestimate error due to autocorrelation).
- Most relevant to **variant comparison (objective #3)**: an A-vs-B claim needs the
  gap to exceed the replicate-to-replicate spread — a single run per variant can't
  support it.
- **REMD/REST2**: because ladder replicas aren't independent, run **2–3 independent
  REMD/REST2 runs** as well (expensive, so often fewer than for plain MD). Within one
  run, standard convergence checks: **neighbour exchange acceptance ≈ 20–40%**
  (theoretical optimum ~20–23%), and **each walker should make several full ladder
  round trips** (≥~5–10) — the round-trip test matters more than raw acceptance.
- Key refs: Grossfield et al. 2019 (*LiveCoMS* 1:5067, the citable best-practices
  article); Knapp, Ospina & Deane 2018 (*JCTC* 14:6127, "The Importance of
  Replicas"); Grossfield & Zuckerman 2009 (*Annu Rep Comput Chem* 5:23); Wang,
  Friesner & Berne 2011 (*JPCB* 115:9431, REST2); Kone & Kofke 2005 / Rathore et al.
  2005 (optimal exchange acceptance).

## What the codebase does today (findings, 2026-10-28)

Foundation is right: `EquilibratedSystem` (preparation/handoff.py) carries
positions + box only — **no velocities** — so every production `run()` draws fresh
velocities. One `prepare()` can feed many velocity-reseeded replicates.

- **Plain MD (`md/production.py`): seeded replicates work, but manually.** Each run
  sets `setVelocitiesToTemperature(T, seed if not None else 0)` and
  `integ.setRandomNumberSeed(seed)` (only if `seed is not None`). Distinct nonzero
  seeds → distinct, individually-reproducible replicates; `seed=None` → fresh-random
  each run (independent, not reproducible).
- **REMD / REST2 (`remd/production.py`, `rest2/production.py`): `seed` is inert.**
  Never read in production — `mcmc.LangevinDynamicsMove(reassign_velocities=True)`
  and `ReplicaExchangeSampler` are built with no seed, so openmmtools seeds its own
  RNG randomly. Re-running gives independent replicates, but **not controllable or
  reproducible via `seed`**. The field is inherited onto `REMDConfig`/`REST2Config`
  but only affects the equilibration that produced `eq`.
- **No replicate abstraction.** No `run_replicates(n=…)`, no per-replicate
  seed/outdir management. A replicate = calling `run(eq, cfg)` N times with distinct
  `outdir` (and, for MD, distinct `seed`).
- **GPU caveat:** even a fixed seed is not bit-reproducible (nondeterministic
  reduction order) — it pins RNG streams + setup, not the trajectory
  (reports.py:49-50).
- **`seed=0` footgun:** code gates on `is not None`, but OpenMM treats a seed of 0 as
  "pick random," so explicit `seed=0` passes the check yet still randomizes. Use a
  nonzero int.

## How openmmtools seeds REMD/REST2 (0.26.0 internals, 2026-10-28)

Verified against the installed source
(`.conda/envs/openmm_env/.../openmmtools/{mcmc,multistate/replicaexchange}.py`).
Three independent RNG sources, **none exposed as a seed argument**:

1. **Per-replica Langevin integrator.** `LangevinDynamicsMove._get_integrator`
   (mcmc.py:1167) builds `openmm.LangevinMiddleIntegrator(...)` and never calls
   `setRandomNumberSeed` (zero hits in mcmc.py). Default seed 0 → OpenMM assigns a
   **unique seed per Context**. The integrator is built lazily *per thermodynamic
   state* inside `apply()`, so **each slot gets its own random seed.**
2. **Velocity reassignment.** With `reassign_velocities=True`, `apply()` calls
   `context.setVelocitiesToTemperature(T)` (mcmc.py:711) with no seed → fresh random
   draw every iteration, per replica.
3. **Swap proposals + Metropolis acceptance.** `ReplicaExchangeSampler._mix_replicas`
   uses the **global** `np.random.randint` / `np.random.rand`
   (replicaexchange.py:324-400) — unseeded.

Consequences:

- **Today each slot gets a *different*, randomly-chosen seed** — no shared-seed
  correlation; replicas are genuinely independent, just not reproducible.
- **Naively "reading the seed" would be a bug.** We pass one shared `move`; giving it
  a single integrator seed makes every replica share one RNG stream → correlated
  thermal noise, partly defeating independent sampling. A correct version must derive
  **distinct per-slot seeds** from one master seed (e.g. `seed+k`, or
  `np.random.RandomState(seed)` drawing N child seeds).
- **Difficulty is split:** swaps (source 3) are trivial to pin — `np.random.seed(seed)`
  once before `create()`/`run()` controls every `np.random.*` in openmmtools.
  Integrator + velocities (sources 1-2) are moderate + invasive: cleanest is a **list
  of moves** (`mcmc_moves=[move_0..move_{N-1}]`), each a `LangevinDynamicsMove`
  subclass overriding `_get_integrator` to `setRandomNumberSeed(seed+k)` and seeding
  the velocity draw — touches openmmtools internals, version-fragile.
- **CUDA caveat dominates the cost/benefit:** on the GPU platform this pipeline
  targets, even full seeding gives **no reproducible trajectory** (non-deterministic
  reduction order; reports.py:49-50), and once dynamics diverge the seeded swap stream
  is moot too. Full seeding yields reproducibility only on CPU/Reference (tests,
  debugging); on GPU, at most a reproducible initial-velocity draw + seed provenance.
  Independent replicates — the actual science goal — already hold for free.

## Open design questions (for a later decision)

1. Wire `cfg.seed` into REMD/REST2. Likely split: do the one-line `np.random.seed`
   (cheap; CPU-test determinism + swap provenance) and skip the invasive per-slot
   integrator/velocity seeding unless CPU-reproducible REMD is specifically needed —
   given the CUDA caveat above, per-slot seeding buys little on GPU.
2. Add a replicate-management layer (fan out N replicates with paired seeds + labelled
   outdirs, aggregate cross-replicate error) vs. leave it to the campaign driver?
3. Default replicate count + the convergence/round-trip diagnostics to surface in the
   run reports.
