# DECISIONS — openmm_pipelines

Design decision records. One entry per significant decision, newest on top. Capture
the *why* and the rejected alternatives, not just the choice.

The decisions below are extracted from `SKELETON_v2.md` (the pre-implementation
design sketch). They are **settled design intent**, not yet validated against
running code except where "Verified" is noted.

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

**Decision.** `NonbondedSpec(cutoff_nm, switch_nm, validated)` is a property of each
`ProteinFF` member. All nonbonded settings flow through `createSystem` kwargs;
nothing is special-cased by FF name or read back off the built System. The `match`
in `nonbonded` has no default — a new FF member with no spec fails loudly
(exhaustiveness).

**Why.** OpenMM applies the kwargs correctly whether the FF uses a `NonbondedForce`
(AMBER) or a `CustomNonbondedForce` (CHARMM, force-switched vdW at 1.2/1.0 nm). No
force introspection = no FF-family branching in `build.py`.

**The `validated` flag** is the structural hook for the CHARMM caution: set it
`False` for any FF whose nonbonded handling has not been energy-matched, and config
validation refuses it until checked. AMBER needs no such check.

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
   (the `validated` gate). AMBER needs no such check.
5. **REST2 force-field safety (future).** CMAP-bearing FFs (CHARMM, ff19SB) are not
   safely solute-scalable without handling the cross-term; gate REST2 behind the
   check.
