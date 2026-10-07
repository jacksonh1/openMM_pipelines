# GOTCHAS — openmm_pipelines

Specific pitfalls discovered in this environment. The repo-root `CLAUDE.md` keeps a
one-line index (phrased as the trigger); this file keeps the **full explanation +
fix**. When a new pitfall is discovered, add its detail here AND add a one-line
entry to the `CLAUDE.md` index — both, every time. Goal: never make the same
mistake twice.

---

## openmmtools `read_mixing_statistics()` returns PER-ITERATION matrices, not cumulative

**Trigger.** Computing exchange acceptance from
`MultiStateReporter.read_mixing_statistics()` (REMD run report / analysis).

**What goes wrong.** The docstring says shape `(n_states, n_states)`, but with the
default `iteration=slice(None)` it returns **per-iteration** stacks of shape
`(n_iterations, n_states, n_states)`. Indexing `n_accepted[i, i+1]` then yields a 1-D
array, not a scalar — `ValueError: The truth value of an array with more than one
element is ambiguous` the moment you compare or divide it.

**Fix.** Sum over the iteration axis first:
`n_accepted = np.asarray(accepted).sum(axis=0)` (same for `proposed`), giving the
cumulative `(n_states, n_states)` matrices. Done in `remd/production._mean_neighbor_acceptance`
and `analysis/remd.exchange_diagnostics`.

---

## mdtraj Trajectory needs float32 xyz/box, and takes no `unitcell_vectors` kwarg

**Trigger.** Building an `md.Trajectory` by hand from de-multiplexed REMD positions
(`analysis/remd.demux_state`) — positions/box come out of openmmtools as float64 Vec3
lists.

**What goes wrong.** Two separate traps: (1) `md.Trajectory(xyz, top,
unitcell_vectors=...)` raises `TypeError` — the ctor takes `unitcell_lengths` /
`unitcell_angles`, not vectors; set the full triclinic box via the
`traj.unitcell_vectors` *attribute* after construction. (2) Passing float64 arrays
makes `image_molecules` crash with `ValueError: Buffer dtype mismatch, expected 'float'
but got 'double'` — mdtraj's C core is hard-wired to float32.

**Fix.** `md.Trajectory(np.asarray(frames, dtype=np.float32), top)` then
`traj.unitcell_vectors = np.asarray(boxes, dtype=np.float32)`.

---

## A new Context resets global parameters to the force's DEFAULT (restraint release)

**Trigger.** Releasing a restraint (or any `addGlobalParameter`-controlled force) with
`context.setParameter(...)`, then building a *new* `Context` from the same `System` —
e.g. production from an `EquilibratedSystem`, or an eq reloaded from disk.

**What goes wrong.** `context.setParameter("k", 0)` changes the value only in *that*
context. The `System`'s force still carries the **default** set by
`addGlobalParameter("k", restraint_k)`. A new `Context` initializes every global to the
force's default — so production would silently re-enable restraints at full strength,
pinning the structure and invalidating the entire run.

**Fix.** Bake the release into the force default before handoff:
`force.setGlobalParameterDefaultValue(index, 0.0)` (then the serialized System and any
new context start at 0). Done in `prepare()` at restraint release; belt-and-suspenders
alternative is to strip the force entirely. Setting only the context value is not
enough.

---

## OpenMM checkpoint (.chk) is NOT portable — use a saveState XML to move between nodes

**Trigger.** Reloading a production restart on a different node, GPU, or after an
OpenMM upgrade; moving the equilibrated system from the prepare job to a production
job on another node.

**What goes wrong.** `simulation.saveCheckpoint()` writes a binary checkpoint that
"can only be loaded into another Simulation that has an identical System, uses the
same Platform and OpenMM version, and is running on identical hardware" (OpenMM
docs). Load it anywhere else and it fails.

**Fix.** For anything that must cross machines/versions, use the **portable** XML:
`simulation.saveState()` / `loadState()` (or `XmlSerializer` of a `State`) — it holds
positions, velocities, and box vectors. We use `saveState` for `prepared/state.xml`
(so any node can launch production from it) and for `production/final_state.xml`.
`production/production.chk` (checkpoint) is kept only for *exact* same-hardware resume
(it also preserves RNG state); never rely on it off the original node. Trade-off: a
state-resumed trajectory is not bit-identical (RNG state is not in the XML).

---

## Density convergence must be a plateau (slope) test, not a consecutive-segment difference

**Trigger.** Checking whether NPT density/box volume has equilibrated.

**What goes wrong.** Segment-to-segment volume noise is well under
`density_tol_rel`, so a *sustained slow drift* just under tolerance passes every
pairwise (segment N vs segment N−1) comparison while the box contracts several
percent overall. The run "converges" on a box that is still shrinking.

**Fix.** Fit a least-squares line through the trailing `density_min_seg` segment
volumes and test the fractional drift that line accounts for across the window
(`density.py: density_converged`). A real plateau has ~zero slope and random scatter
passes; a steady drift is caught. This is why the test is a slope over a window, not
a difference between adjacent segments.

**Origin.** Inherited lesson from the sibling GROMACS pipeline; baked into the
`openmm_pipelines` density design from the start.

---

## Modeller.addSolvent rejects most water-model names — pack by site count, not name

**Trigger.** Solvating with a water model other than TIP3P/TIP4P-Ew (e.g. OPC,
OPC3, TIP3P-FB).

**What goes wrong.** `Modeller.addSolvent(model=...)` ships pre-equilibrated boxes
for only a handful of names and raises on any other:
`model='opc'` -> `ValueError: Unknown water model: opc`.

**Fix.** Choose the packing box by the water's **site count**
(`_PACKING_BOX = {3: "tip3p", 4: "tip4pew", 5: "tip5p"}`) and let the actual
parameters come from the FF xml; minimization fixes residual geometry. Verified
(OpenMM 8.6.1): `model='tip4pew'` + `opc.xml` yields real 4-site OPC (charges
O=0.0, H=+0.6791, M=−1.3583 — not TIP4P-Ew's +0.5242), one virtual site per water.
`build.py: _assert_water_topology` backstops this by checking every `HOH` has the
expected particle count and 4-site models have virtual sites. See the water-packing
decision in `DECISIONS.md`.

---

## CHARMM nonbonded: OpenMM's switch is a POTENTIAL switch, not CHARMM's force-switch

**Trigger.** Running a CHARMM force field in OpenMM and trusting the energetics for a
large/long campaign.

**What's subtle (and could silently waste compute).** OpenMM does NOT pick nonbonded
settings from the FF — `createSystem` defaults to a 1.0 nm cutoff and *no* switch
([docs](https://docs.openmm.org/latest/userguide/application/02_running_sims.html)).
We set CHARMM to `1.2 nm cutoff / 1.0 nm switch / PME / HBonds`, which matches the
canonical ParmEd CHARMM/OpenMM example (12 Å / 10 Å). BUT CHARMM was parameterized
with **force-switching** (vfswitch), while OpenMM's built-in `NonbondedForce`
`switchDistance` applies a **potential switch** (multiplies the LJ energy by
`S = 1−6x⁵+15x⁴−10x³`) — a *different* function that leaves a small force "bump"
between switch and cutoff. OpenMM's built-in force cannot do force-switch; CHARMM-GUI
achieves the exact form via a `CustomNonbondedForce`.

**Fix / safeguard.** The 1.2/1.0 potential-switch is the standard documented
approximation and is fine for conformational stability/drift MD. For *exact* CHARMM
energetics, implement the CHARMM-GUI force-switch (`CustomNonbondedForce`) — out of
current scope. Either way, the real protection against a wrong setting is a **one-time
single-point energy match** vs a reference
([choderalab/OpenMMEnergyComparisons](https://github.com/choderalab/OpenMMEnergyComparisons))
before committing a production campaign. AMBER (1.0 nm, no switch, PME) needs no such
caution. See `DECISIONS.md` (nonbonded settings) and TODO (energy-match tool).

---

## CHARMM36 != CHARMM36M

**Trigger.** Picking a CHARMM protein force field.

**What goes wrong.** They are different force fields. CHARMM36M (`charmm36_2024.xml`,
`par_all36m_prot`) is the Huang-2017 backbone-CMAP + glycine refinement, usually the
right choice for folding/stability systems (e.g. WW domains). CHARMM36
(`charmm36.xml`, `par_all36_prot`) is the original. Silently using the wrong one
degrades results in a way that looks plausible.

**Fix.** Separate `ProteinFF` enum members (`CHARMM36` vs `CHARMM36M`) with distinct
`family` solvent trees. Choose deliberately.

**Also:** before trusting any CHARMM run, energy-match a known system to confirm the
switching / `CustomNonbondedForce` handling (the `NonbondedSpec.validated` gate).
AMBER needs no such check.

---

## SLURM: sbatch scripts are copied to a temp path (for future drivers)

**Trigger.** Writing an sbatch driver that references repo files relative to itself.

**What goes wrong.** `BASH_SOURCE[0]` inside an sbatch script does *not* point to the
repo — SLURM copies the script to a spool/temp path before running it. Paths
resolved from `BASH_SOURCE` break.

**Fix.** Use `$SLURM_SUBMIT_DIR`, or pass the repo path explicitly via `--export`.

**Status.** Inherited from the sibling GROMACS pipeline. No sbatch drivers exist in
`openmm_pipelines` yet; recorded here so the first one avoids it.

---

## Output model: folder-symlink stages (for future drivers)

**Trigger.** Pointing bulk stage output dirs at fast scratch via symlinks.

**What goes wrong.** If bulk stage dirs (`prod/ equil/ density/ heat/`) are symlinks
into scratch: never symlink `build/`/`em/` (a cross-filesystem rename fails);
pre-create the stage symlinks before `mdrun`/production; and a `finalize_outputs.sh`
that walks the tree must skip symlinks (`[[ -L ]]`) or it recurses into scratch.

**Status.** Inherited verbatim from the sibling GROMACS pipeline. The
`openmm_pipelines` design writes plain files to `outdir` (`prod.xtc`, `final.xml`,
`final.chk`, `config.json`) with **no** scratch-symlink scheme, so this does not
apply today. Recorded in case a scratch-staging output model is added later.

---

## Bash: a `[[ … ]] && echo` as the LAST line of a `set -e` script exits 1 (for future drivers)

**Trigger.** Ending a `set -e` shell script (e.g. an sbatch) with a bare test-and-act
one-liner.

**What goes wrong.** A false test (`[[ … ]] && echo`, `grep -q`, `((n++))`) as the
*final* statement makes the script exit 1 — so SLURM reports a successful job as
FAILED.

**Fix.** Use a real `if` block for the final statement, or add an explicit
`exit 0` / trailing `true`.

**Status.** Inherited from the sibling GROMACS pipeline; no shell drivers here yet.

---

## An OpenMM Force proxy from a temporary System reads freed memory

**Trigger.** `scaled = scale_solute(system, ...)` then inspecting a force — written as the
one-liner `scale_solute(system, ...).getForce(0)` (or any `make_system().getForce(i)`).

**What goes wrong.** `System.getForce(i)` returns a SWIG proxy that only *references* the
C++ System; it does not keep it alive. If the System has no surviving Python reference, it
is garbage-collected while the proxy is still in use, so reads hit freed memory —
`getParticleParameters` returns garbage denormals (e.g. `4.67e-310`) or `getMapParameters`
raises `OpenMMException: std::bad_alloc`. Bit it silently *looks* like a scaling bug.

**Fix.** Bind the System to a variable that outlives every proxy taken from it:
`scaled = scale_solute(...); nb = scaled.getForce(0)`. Not REST2-specific — any transient
OpenMM container (System/Force/Context) whose child proxy outlives it. Surfaced writing the
`scale_solute` tests (2026-10-06).

---

## openmmtools' global context cache cannot be re-platformed once populated — `run()` must empty it

**Trigger.** Calling a multistate `run()` (`remd`/`rest2`) twice in one process — a driver
that runs two backends, or the test suite running both slow multistate tests.

**What goes wrong.** openmmtools reads the compute platform off its *global* context cache,
and `cache.global_context_cache.set_platform(...)` raises `RuntimeError: Cannot change
platform of a Context cache already in use` once the cache holds any context — even when
setting the *same* platform. So a second `run()` in the process crashes, regardless of
platform.

**Fix.** `cache.global_context_cache.empty()` before `set_platform` makes the setup
idempotent (and frees the prior run's cached contexts). Centralized in
`openmm_pipelines/context_cache.configure_global_platform`, used by both backends.
Surfaced 2026-10-07 when a second multistate slow test (REST2) was added alongside REMD.
