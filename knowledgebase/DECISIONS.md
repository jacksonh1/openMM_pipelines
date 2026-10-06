# DECISIONS — openmm_pipelines

Design decision records. One entry per significant decision, newest on top. Capture
the *why* and the rejected alternatives, not just the choice.

The decisions below are extracted from `SKELETON_v2.md` (the pre-implementation
design sketch). They are **settled design intent**, not yet validated against
running code except where "Verified" is noted.

---

## 2026-10-06 — T-REMD: NVT production + de-multiplex to fixed-T ensembles

**Decision.** The T-REMD backend runs **NVT** production, uses a **geometric**
temperature ladder, and the analysis **de-multiplexes** replica trajectories into
fixed-temperature ensembles before measuring anything.

**Why NVT.** `prepare()` already equilibrates the box to the target density under NPT;
production then holds that volume fixed. Replica exchange over temperature only (one
`ThermodynamicState` per rung) is the clean T-REMD formulation — a barostat per replica
would make the exchange acceptance depend on volume work and complicate the weights for
no benefit to the stability/variant-comparison objectives. So production strips the
barostat (`_nvt_system`). *Rejected:* NPT-REMD (unnecessary here; volume held is fine
once density is equilibrated).

**Why geometric ladder.** For roughly constant heat capacity, equal temperature ratios
give roughly equal neighbour exchange-acceptance — the right target, since one weak
("broken") link stalls replica diffusion through temperature space. It is a *starting*
ladder; `n_replicas` is tuned from the measured per-neighbour acceptance. *Rejected:*
hand-tuned/optimized ladders up front (premature; `np.geomspace` + a diagnostics plot
is enough to iterate).

**The critical concept (CLAUDE.md), resolved.** A NetCDF reporter stores **replica**
trajectories — continuous coordinate sets that *walk through temperature space* as
swaps reassign their state. A single replica's frames are NOT an ensemble at one
temperature. To characterize the design at its reference temperature you must
**de-multiplex**: at each iteration take the replica currently in the target state.
`analysis/remd.demux_state(k)` does this (from checkpointed frames); only the demuxed
T_min ensemble feeds RMSD/Rg/RMSF (via the shared `analysis/drift.drift_report`
primitive, reused by plain MD too). `exchange_diagnostics` reports the permutation +
neighbour acceptance + openmmtools' mixing statistics (transition-matrix subdominant
eigenvalue + statistical inefficiency — the library metric, not a hand-rolled round-trip
counter) so the mixing is auditable. This is the distinction the CLAUDE.md "critical
REMD/REST2 concept" demands be made explicit.

**Per-replica equilibration.** `sampler.equilibrate(n)` runs each replica at its own
temperature and is discarded — the designed replacement for the dropped shared-reference
`relax` (TODO item), so no rung starts production out of equilibrium at its own T.

**Gating.** `remd/` + `analysis/remd.py` require the `[remd]` extra (`openmmtools`) and
are not eagerly imported by the package `__init__`s; `config.py`/`ladder.py` stay
OpenMM-free so configs/ladders build on any node.

---

## 2026-10-01 — design validated against OpenMM docs/cookbook

**Decision.** Keep the architecture as sketched; confirmed it matches OpenMM's own
documented idioms before writing more code.

**Why / what was checked.**
- **Seam placement.** The canonical app-layer flow (User Guide: *Running
  Simulations*) is linear through `PDBFixer → ForceField → Modeller.addSolvent →
  createSystem(PME, HBonds, hydrogenMass) → LangevinMiddleIntegrator → Simulation →
  minimizeEnergy → setVelocitiesToTemperature → MonteCarloBarostat → reporters →
  step`. The prepare/produce seam cuts exactly at `step`: `core/` owns everything
  before it, each engine owns `step`.
- **Handoff completeness.** `app.Simulation` needs `(system, positions, box)`;
  `openmmtools.multistate` needs `ThermodynamicState(system, T)` +
  `SamplerState(positions, box_vectors)`. `EquilibratedSystem` is the union — no
  backend assumption leaks past the seam. Confirmed against the openmmtools
  multistate docs.
- **Restraint idiom.** The cookbook's *Restraining Atom Positions* uses exactly
  `CustomExternalForce('k*periodicdistance(...)^2')` with a global `k` released via
  `k=0`, no rebuild — verbatim our `restraints.py`.
- **4 fs HMR default** matches the docs' own `LangevinMiddleIntegrator(..., 0.004*ps)`.
- **REST2 is genuinely the hard part.** OpenMM's REST cookbook builds it *manually*
  (custom Bond/Angle/Torsion forces + per-term `lambda_rest_*` globals +
  `CompoundThermodynamicState`); there is no turnkey helper, and the tutorial
  silently ignores CMAP. This confirms the `rest2/` CMAP safety gate covers a real
  gap, not a hypothetical one.

Sources (2026-10-01): OpenMM User Guide *Running Simulations*; OpenMM Cookbook
*Restraining Atom Positions* and *Running a REST simulation*; openmmtools
*multistate* docs.

---

## 2026-10-01 — package layout: per-engine folders; naming: prepare/production, run(eq, cfg)

**Decision.** Reorganize from SKELETON's flat module list into per-engine
subpackages. Top level keeps only `__init__.py` + `forcefield.py`. Shared
equilibration lives in `core/`; each engine (`md/`, `remd/`, `rest2/`) is a
self-contained folder with its own `config.py` + `production.py` (+ `ladder.py` for
remd). See `CODEBASE_SHAPE.md` for the tree.

**Naming.** The seam verbs are `prepare` (setup side, in `core/`) and `production`
(run side). The `produce_*` nomenclature from SKELETON is dropped — "production" is
the literal MD term. Each engine exposes **one** public production function,
`run(eq, cfg) -> Path`, namespaced by package (`md.run`, `remd.run`, `rest2.run`).

**`run` always takes an `EquilibratedSystem`; there is no `cfg`-only wrapper.**
Production only ever runs on an equilibrated system, and the serialized
`EquilibratedSystem` is the point of the seam (equilibrate once — often a short
separate job — then launch a long production run, possibly across a temp sweep or
REMD ladder). A `run(cfg)` convenience would hide the expensive equilibration and
invite accidental re-prep. The one-shot is explicit composition:
`eq = prepare(cfg); md.run(eq, cfg)` — two lines, reads like an OpenMM script,
and because `MDConfig` IS-A `PrepConfig` the *same* config object flows to both.

**Also:** `forcefield.py` is pulled out of `config.py` (the `ProteinFF`/`WaterModel`/
`NonbondedSpec`/`resolve` machinery is ~100 lines and single-purpose);
`core/config.py` holds only `PrepConfig` + `BoxShape`; each engine's `config.py`
holds just its subclass. Config inheritance crosses package boundaries by design
(`md/config.py` imports `PrepConfig` from `core`) — a one-line import, flat
attribute access preserved.

**Rejected.** Flat layout (SKELETON) — ~11 top-level files once REMD/REST2 land.
"Split by the seam" (`prep/` + `backends/`) — considered; per-engine folders won
for self-containment (each engine's config + backend + helpers in one place).
Keeping `produce`/`simulate`+`run` dual functions — collapsed to one `run(eq, cfg)`.

---

## 2026-10-01 — final module tiers; `lib/`; honest stage names; relax dropped

A naming/layout working session replaced the `core/` grab-bag with four tiers:
shared **types** at root (`config.py`, `forcefield.py`), the reusable **toolkit**
in `lib/`, the shared **pipeline** in `preparation/`, one folder per **engine**.

**`core/` → split into `lib/` + `preparation/`.** "core" was too generic and mixed
two kinds of thing: reusable helpers and pipeline stages. Now `preparation/` holds
*only* the pipeline (composer + stages + the handoff output type), and `lib/` holds
the toolkit.

**`config.py` moved to the package root.** `PrepConfig` is the base config for
*every* engine, so it is not part of the preparation pipeline — it belongs at the
root with the other shared types. Each engine's config (`MDConfig`, …) stays in its
own folder and subclasses it.

**The `lib/` tier rule.** A module goes in `lib/` only if it is a *substantial,
reusable* building block: used across multiple stages (`restraints`: add → reanchor
→ release span the whole pipeline) OR pure and independently testable
(`density_convergence`: numpy, no OpenMM — a kernel). Trivial glue that serves one
stage stays with that stage — so `add_barostat` (2 calls, used only at the NVT→NPT
transition) lives in `equilibration_steps.py`, not `lib/`. (`lib/` matches the
sibling FragForge convention; verified no `lib/` pattern in `.gitignore`/
`.mutagenignore` so the source tracks and syncs.)

**Honest stage names.** `density.py` → `lib/density_convergence.py` (it is the
convergence logic, not a verb). `stages.py` → `equilibration.py` →
`equilibration_steps.py`, and it now holds *only* the restrained equilibration
dynamics (`equilibrate_nvt`, `add_barostat`, `equilibrate_density`). `minimize` is
static energy minimization — *not* equilibration — so it is its own `minimize.py`.
The earlier grab-bag (minimize + nvt + barostat + density + relax under
"equilibration") was the core dishonesty; splitting the non-equilibration ops out
fixes it. Line counts drove the granularity: the op bodies are tiny (nvt 3 lines,
relax 1), so one-file-*per-op* would be mostly boilerplate — group the coherent
equilibration dynamics, split only the genuinely different `minimize`.

**`relax` dropped entirely** (function + `relax_ns` + `relax_steps`). It was an
optional unrestrained NPT at the reference T after releasing restraints, default
0.0. It *contradicts* the tool's analysis contract: restraints are released only at
production start so the production trajectory captures *all* drift from the input
pose; a relax window lets early drift happen off-camera, so production under-reports
it. There is also no restraint-release "shock" to absorb (`k=0` injects no energy;
the thermostat handles it). If initial-settling exclusion is ever wanted it belongs
in **analysis** (discard initial frames, or `pymbar.timeseries.detect_equilibration`
— pymbar is in `openmm_env`), not baked into every prepare run. The REMD/REST2
argument for a pre-production settle is real but points elsewhere: those need
*per-replica* equilibration at each replica's own T/λ (a backend concern, via the
multistate sampler), not a single reference-T relax in shared `prepare()`. Deferred
to the REMD/REST2 backends (see TODO).

No `cfg`-only `run` wrapper, one `run(eq, cfg)` per engine — unchanged from the
earlier layout decision.

---

## 2026-10-01 — output contract (FIRST DRAFT — will be refined)

**Status: first draft.** Captured so implementation can proceed; expected to change as
we learn what analysis/QC actually needs. Goal: capture every parameter + piece of
information needed for reproducibility and QC.

**Layout** (one `outdir`, split by the seam):
```
outdir/
  config.json          # all PrepConfig/MDConfig fields + derived step counts (intent)
  run_info.json        # openmm/python/dep versions, platform+precision, GPU, timestamp,
                       #   host, SLURM id, RESOLVED ff xml list + water + packing box
  input.pdb            # exact input, copied (no hash — the copy is the record)
  prepared/
    system.xml         # serialized System — authoritative record of forces/params
    state.xml          # saveState XML (PORTABLE): positions + velocities + box
    topology.cif       # mmCIF (handles >99,999 atoms; PDB serial overflows)
    prep_report.json   # QC: build + minimize + density-equilibration metrics
  production/
    trajectory.xtc
    production.log     # StateDataReporter: step,time,PE,KE,totalEnergy,T,V,density,speed
    final_state.xml    # saveState (portable restart)
    production.chk      # saveCheckpoint (exact restart; same hardware+openmm ONLY)
    run_report.json    # QC: timing, T/density/energy stats, stability flag
```

**Grounded in OpenMM docs (2026-10-01):**
- `prepared/state.xml` is a **`saveState` XML** (portable: positions+velocities+box),
  *not* a checkpoint — that is what lets any node load the equilibrated box (the
  seam's "equilibrate once, produce anywhere"). Docs: a checkpoint "can only be loaded
  into another Simulation that has an identical System, uses the same Platform and
  OpenMM version, and is running on identical hardware."
- Production keeps **both** `final_state.xml` (portable resume) and `production.chk`
  (exact same-hardware resume, keeps RNG state).
- `production.log` includes **KE + totalEnergy** (energy-drift QC), confirmed
  reportable by `StateDataReporter`.
- `System` XML excludes the `Topology`, so `topology.cif` is written separately
  (`PDBxFile.writeFile`).

**Typed, not raw dicts:** `run_info.json` / `prep_report.json` / `run_report.json` are
serialized from `RunInfo` / `PrepReport` / `RunReport` dataclasses (fail-loud style).

**Not content-addressed:** no input hash, no FragForge-style `deterministic_id` — runs
are identified by `outdir`. (Earlier draft had a sha256; dropped as unjustified here.)

**Caveat (documented in run_info):** GPU runs are not bit-reproducible even with
`seed` (nondeterministic reduction order); `seed` pins setup + RNG streams, not the
trajectory. `saveState` also omits RNG state, so a state-resumed trajectory differs.

---

## 2026-10-01 — handoff dataclasses live in `preparation/handoff.py` (not `equilibrated.py`)

**Decision.** The module holding `BuiltSystem` + `EquilibratedSystem` is
`core/handoff.py`.

**Why.** The module holds *data classes* (the objects that cross the prepare/produce
seam), not the *action* of equilibrating — the verb lives in `stages.py` /
`prepare.py`, so `equilibrate.py`/`equilibration.py` would mislead. It also holds
`BuiltSystem`, which is *pre*-equilibration, so `equilibrated.py` (the original name)
was too narrow. `handoff` names the role we already use in prose for these objects.

---

## 2026-10-01 — GPU tests via SLURM on the lab's pinned nodes

**Decision.** `slow` tests (real OpenMM dynamics) run on a GPU through
`tests/gpu_tests.sbatch` (partition `pi_keating`, `gpu:l40s:1`, env `openmm_env`),
pinned to the user's usual nodes with `-w node3620` / `-w node3619`. Test platform is
env-driven: `OPENMM_TEST_PLATFORM` (default `CPU`, set to `CUDA` by the sbatch). CPU
runs of the solvated-system dynamics tests are too slow for an interactive loop
(>120 s) — GPU is the right harness.

**Why.** Keeps the fast pure/`build` tests runnable anywhere while the dynamics tests
get a GPU without hardcoding a platform. Mirrors the lab's existing
`Ellen_WW_simulations/openmm_env_test` sbatch header.

---

## 2026-10-01 — prepare() writes topology; from_disk(outdir) is one-arg

**Decision.** `prepare()` writes a topology file (PDB/CIF) alongside `system.xml` +
the State xml. `EquilibratedSystem.from_disk(outdir)` reloads everything from
`outdir` with no external arguments.

**Why.** A `System` xml does **not** carry the `Topology`. SKELETON's
`from_disk(outdir, topology)` required the caller to supply the topology from
elsewhere — which breaks the seam's core promise that a separate production job
starts from `outdir` alone. Writing the topology also serves analysis (we need it to
read trajectories and view results) and PBC imaging on export, so it is output we
want regardless. One-arg reload is the ease-of-use win.

**Deferred.** The full output contract (what files, what layout, naming) is a
separate discussion — see TODO.

---

## 2026-10-01 — prepare/produce seam as the one architectural axis

**Decision.** Split the engine at production: `prepare(cfg) -> EquilibratedSystem`
is shared by all backends; `produce_md` / `produce_remd` / `produce_rest2` differ.
`prepare()` returns plain OpenMM objects (serialized `System` + `State`).

**Why.** Everything up to production (build, restrain, minimize, NVT, NPT density,
release) is identical across plain MD, T-REMD and REST2 — ~90% of the code. Plain
MD drives one `app.Simulation.step()`; REMD/REST2 drive many `Context`s via an
`openmmtools.multistate` sampler and do **not** use `app.Simulation`. Returning
plain OpenMM objects is the only handoff both consume, so no `app.Simulation`
assumption leaks past the seam. Serializing it means equilibrate-once / produce-many
(and a separate production job can start from the same box).

**Rejected.** A monolithic `run()` per engine (would duplicate equilibration three
ways); passing an `app.Simulation` across the seam (multistate can't consume it).

---

## 2026-10-01 — config by inheritance, not nesting

**Decision.** `PrepConfig` (pydantic `BaseModel, frozen=True`) holds everything
through equilibration. `MDConfig(PrepConfig)` adds production fields. Future
`REMDConfig`/`REST2Config` subclass the same `PrepConfig`.

**Why.** Keeps attribute access flat (`cfg.temperature_k` on every engine's config)
and lets `prepare()` depend only on `PrepConfig` fields with no duplication of
fields or validators. `temperature_k` is defined as the *reference/equilibration*
temperature (= T_min for REMD), so the shared base is already engine-agnostic.

**Rejected.** A nested `cfg.prep.temperature_k` shape (noisier access, no real
benefit since the fields are genuinely shared).

---

## 2026-10-01 — water packing by site count, parameters by force field

**Decision.** `Modeller.addSolvent` is called with a packing box chosen only by the
water's **site count** (`_PACKING_BOX = {3: "tip3p", 4: "tip4pew", 5: "tip5p"}`);
the actual water *parameters* come from the FF xml. A single `resolve(protein, water)`
returns both the FF xml list and the packing-box name, and raises on an illegal
pairing.

**Why.** `Modeller.addSolvent` ships pre-equilibrated boxes for only a few geometry
names and rejects any other (`model='opc'` -> `ValueError`). The packing box only
needs the right site count; parameters come from the xml and minimization fixes the
residual geometry. This is OpenMM's own documented idiom ("a box of TIP4P-Ew water
can be used for most four-site water models"). Adding a new non-polarizable water is
then one `WaterModel` entry with its site count — packing is derived, never
special-cased.

**Verified** (openmm_env, OpenMM 8.6.1): `model='opc'` raises
`ValueError: Unknown water model: opc`; `model='tip4pew'` + `opc.xml` yields real
4-site OPC (charges O=0.0, H=+0.6791, M=−1.3583 — *not* TIP4P-Ew's +0.5242), one
virtual site per water. Backstopped at runtime by `_assert_water_topology` in
`build.py`.

**Scope limit.** Drude/polarizable water (e.g. swm4ndp) is out — needs a Drude
integrator.

---

## 2026-10-01 — FF/water pairing is enforced, not advisory

**Decision.** `resolve()` raises `ForceFieldError` on a correctness-breaking
pairing (any non-TIP3P water with a CHARMM FF, since CHARMM water is version-locked
to its bundled modified-TIP3P). A *degraded-but-legal* pairing (e.g. ff19SB+TIP3P)
only warns, via `_RECOMMENDED_WATER`.

**Why.** Mixing an AMBER-tree water with CHARMM is a silent correctness error;
make the illegal state unrepresentable. ff19SB+TIP3P is scientifically degraded but
not *wrong*, so it warns rather than blocks.

---

## 2026-10-01 — nonbonded settings are properties of the force field

**Decision.** `NonbondedSpec(cutoff_nm, switch_nm)` is a property of each `ProteinFF`
member. All nonbonded settings flow through `createSystem` kwargs; nothing is
special-cased by FF name or read back off the built System. The `match` in
`nonbonded` has no default — a new FF member with no spec fails loudly
(exhaustiveness).

**Why this file must exist (verified against OpenMM docs, 2026-10-01).** OpenMM does
**not** pick nonbonded settings from the force field: `createSystem` defaults to a
1.0 nm cutoff and *no* switch, and loading `charmm36.xml` does not auto-apply CHARMM's
force-switch. But CHARMM's LJ params were parameterized *with* a 10–12 Å force-switch,
so a plain cutoff with CHARMM is a silent correctness error. Making cutoff+switch a
property of the `ProteinFF` enum means you cannot select CHARMM and forget its switch
(fail-loud / make-illegal-states-unrepresentable). The kwargs apply correctly whether
the FF uses a `NonbondedForce` (AMBER) or a `CustomNonbondedForce` (CHARMM) — no force
introspection, no FF-family branching in `build.py`.

**Removed the `validated` flag (was in the first draft).** It was a bool that was
always `True`, gating a config assertion for an energy-match check never performed —
speculative machinery (YAGNI). Cut it and its config assert. The CHARMM energy-match
remains a **manual recommendation** (science-review flag #4), not a code gate.

---

## 2026-10-01 — restraints via a global `k`, released at production start

**Decision.** Harmonic position restraint on protein heavy atoms via a
`CustomExternalForce` with `k` as a **global parameter**. Reference coordinates are
re-anchored to the minimized pose once (`reanchor`), held at full strength through
*all* of equilibration, and released with a single
`context.setParameter("k", 0.0)` at production start. The force stays in the System
with `k=0`.

**Why.** Release is one call — no force removal, no context rebuild. Holding
restraints through equilibration and releasing only at production means the
production trajectory records drift from the *input pose* (the analysis reference),
not drift accumulated during equilibration steps that are never analyzed. The
leftover `k=0` force is harmless for MD; REMD/REST2 may strip it for cleanliness
before building replica states (a one-liner on the serialized system).

---

## 2026-10-01 — adaptive density plateau (slope) test, fail-loud

**Decision.** NPT density equilibration runs in fixed-length segments; after each,
mean box volume is recorded and a least-squares line is fit through the trailing
`density_min_seg` segments. Stop when the fractional drift that line accounts for
across the window is `<= density_tol_rel`. If no plateau in `density_max_seg`
segments, raise `DensityNotConverged`.

**Why.** A **plateau (slope) test, not a consecutive-segment difference.** A
sustained slow drift just under tolerance passes every pairwise segment comparison
while the box contracts several percent overall; the slope over the window catches
it, random scatter has ~zero slope and passes. This exact mistake was made in the
GROMACS pipeline (see `GOTCHAS.md`). Fail-loud over silently shipping an
unconverged box.

---

## 2026-10-01 — 4 fs + HMR as the default timestep

**Decision.** Default `dt_ps=0.004` with `hydrogen_mass_amu=4.0` and
`constraints=HBonds`. Config asserts: `dt_ps > 0.0025` requires
`hydrogen_mass_amu >= 3.0`. Drop to 2 fs via physical H mass + `dt_ps=0.002`.

**Why.** HMR repartitions hydrogen mass so a 4 fs step is stable, ~2× throughput.
4 fs without repartitioned H is unstable — the assertion makes that
unrepresentable.

---

## 2026-10-01 — dependency minimalism: no openmmforcefields, mdtraj only, openmmtools optional

**Decision.** Core OpenMM bundles every protein FF needed — no `openmmforcefields`.
`mdtraj` is the single trajectory library (no MDAnalysis). `openmmtools` is an
optional `[remd]` extra, not a core dependency of plain MD. Trajectories are XTC via
`app.XTCReporter`.

**Why.** Verified that OpenMM 8.6.1 loads all needed AMBER/CHARMM FFs directly
(modular paths preferred over `-all` bundles, which also pull nucleic-acid params).
Keeping `openmmtools` optional keeps the plain-MD install lean; the seam makes it a
clean extra. XTC is compact and read natively by mdtraj; verified present in
OpenMM 8.6.

---

## Open decisions (unresolved)

- REMD production ensemble default: NVT or NPT?
- Where the REMD/REST2 backends ultimately live: in this package under a `[remd]`
  extra, or a separate package that imports `openmm_pipelines.prepare`. (Either
  works given the seam.)

## Science-review flags to resolve before trusting production

1. **Restraint reference vs MC barostat.** Restraint `x0` are absolute coords; the
   barostat rescales molecule COMs. Confirm benign (restraints on solute; barostat
   scales by molecule).
2. **Trajectory / final-PDB PBC wrapping.** Wrap with `mdtraj` `image_molecules` on
   export; never ship a raw box.
3. **`protein_heavy_atoms` selection.** Protein heavy atoms only, no H, no
   solvent/ions; assert non-empty.
4. **CHARMM energy-match.** Energy-match a known system before trusting a CHARMM run
   (manual recommendation — no longer a code gate; the `validated` flag was removed).
   AMBER needs no such check.
5. **REST2 force-field safety (future).** CMAP-bearing FFs (CHARMM, ff19SB) are not
   safely solute-scalable without handling the cross-term; gate REST2 behind the
   check.
