# DECISIONS — openmm_pipelines

Design decision records. One entry per significant decision, newest on top. Capture
the *why* and the rejected alternatives, not just the choice.

The decisions below are extracted from `SKELETON_v2.md` (the pre-implementation
design sketch). They are **settled design intent**, not yet validated against
running code except where "Verified" is noted.

---

## 2026-11-02 — fold-retention metrics: Q + dRMSD + lDDT, alignment-free, vs a topology-matched design

**Decision.** Add three alignment-free drift metrics to `drift_report`, all referenced to the design
pose: native-contact fraction **Q** (Best, Hummer & Eaton 2013 — all-heavy-atom soft cutoff: 4.5 Å
native, β=50/nm, λ=1.8, |i−j|≥4), **dRMSD** (Cα distance-RMSD), and **lDDT** (Mariani 2013, via
`biotite.structure.lddt`) over one or more atom sets (`LddtAtomSet`: `ca_cb`/`ca`/`backbone`/`heavy`),
default `(ca_cb, ca)`, whole-system, thresholds fixed at the GDT-HA set 0.5/1/2/4 Å, R0=15 Å. The
per-atom metrics reference a **topology-matched** design pose (`load_prepared_reference`, the protein
slice of `prepared/topology.cif`), not raw `input.pdb`.

**Why.**
- *Alignment-free is the point.* These designs are sampled in the folded basin (not folding/unfolding —
  too slow to sample), so the question is "how far from the design." Distance-based metrics need no fit
  region, so a floppy tail or hinge doesn't corrupt the measure the way it biases RMSD's superposition.
- *Q alone is blunt near the folded state.* λ=1.8 keeps a contact "intact" until ~1.8× its native
  distance, so Q saturates near 1 and poorly ranks stable variants. lDDT's graded 0.5/1/2/4 Å thresholds
  give the near-native sensitivity Q lacks, plus a per-residue profile (serves flexible-region ID).
- *Atom-set knob, default Cα+Cβ.* All-heavy lDDT carries fast side-chain motion (ring flips, rotamers) —
  noise for fold retention. Cα+Cβ captures backbone + coarse side-chain orientation while immune to
  distal wiggle; Cα is the pLDDT-style named standard. A free list of sets (not a boolean), each plotted
  and labeled `lDDT (<label>)`, with the set slug in the filename/npz key.
- *Topology-matched reference (not input.pdb).* Q/lDDT compare specific atoms, so reference and
  trajectory must share the exact atom layout. `input.pdb` is the raw design (no H, PDBFixer may rebuild
  atoms, PDB vs FF naming) — not index-aligned to the simulated system; the prepared topology slice is,
  and is geometrically identical to the design (prep moves no existing heavy atom). dRMSD (Cα-only) could
  use input.pdb but shares the prepared reference for consistency. RMSD left on input.pdb (unchanged).

**Rejected / not done.**
- *Residue-level closest-heavy hard-cutoff Q* (FragForge `md_quant.contact_fraction` style) — robust but
  coarser and a step function; the paper's soft per-atom Q is exact here thanks to the matched reference.
- *`exclude_same_chain` lDDT knob* — dropped as confusing (biotite's `True` = inter-chain-**only**, not
  intra-chain; there is no intra-chain-only switch). Whole-system lDDT is the default; an explicit
  interface-lDDT / per-chain-lDDT can be added later as named outputs.
- *Unifying all metrics on the prepared reference* (dropping input.pdb for RMSD) — deferred; RMSD would
  be numerically identical (same Cα coords) but it changes which file RMSD reads.

**Library.** biotite added (lDDT only; mdtraj has none). Reinforces the md_quant principle as "no
MDAnalysis / MD-specific trajectory libraries," NOT "mdtraj-only" — general structure libraries are fine.

## 2026-10-31 — prep never models missing residues implicitly

**Decision.** `build_system` builds **no** missing residues. Internal missing-residue gaps hard-fail
(`_assert_no_internal_gaps` -> model them externally); terminal tails (SEQRES longer than the
coordinates) are skipped with a warning (`_drop_missing_residues`) and the structure is simulated as
given. Missing atoms/hydrogens within existing residues are still completed.

**Why.** The tool's premise is "start from the given designed structure." Silently completing a chain
invents residues — terminal runs become an idealized extended tail that inflates the solvent box
(~143 -> ~57 Å on a WW domain), adds a floppy segment, and distorts the very structure being assessed,
with no error raised. Measuring drift against geometry PDBFixer invented is scientifically misleading.
Chain completion, if ever wanted, must be an explicit opt-in and may be out of scope for this tool.

**Rejected.** (1) The prior implicit behavior — `addMissingAtoms` building terminal residues (the bug).
(2) Auto-filling small internal gaps via PDBFixer — kept as a hard fail; PDBFixer's modeled loops are
crude and unrefined, a liability for a structure under assessment.

## 2026-10-30 — resume / extend: continue an interrupted run (library-only, documented-mechanism)

**Decision.** A top-level `resume(outdir, *, extend_to_ns=None, platform=None)` continues an
interrupted (walltime/crash) production run — or extends it to a larger target — using the
**documented OpenMM / openmmtools restart mechanisms**, not a hand-rolled scheme. All three engines;
MD is first-class. Library-only: re-submitting the SLURM job (so it calls `resume` not `run`) is the
driver's job — no SLURM coupling in the library (`resume` reads disk state only). Reviewed by the
architect agent; correctness fixes folded in. Verified end-to-end on CPU (`tests/test_resume.py`).

**Mechanisms (verified on openmm 8.6.1 / openmmtools 0.26.0).**
- **MD**: portable State XML as the resume point (`simulation.saveState`/`loadState` — any node;
  drops the RNG stream, so the continued trajectory differs but is physically valid — GPU runs are
  not bit-reproducible anyway). Trajectory/log continue via reporter `append=True`. The State is
  written **ping-pong to two alternating files** so a SIGKILL mid-write never destroys the prior good
  one (OpenMM's `CheckpointReporter` overwrites a single file in place — no torn-write protection).
- **REMD/REST2**: openmmtools `ReplicaExchangeSampler.from_storage(reporter)` + `run()` (finish to the
  stored target) / `extend(n_iterations)` (extend past it). openmmtools rebuilds the whole sampler
  from the `.nc`, so resume needs no `eq`/ladder/re-scaling.

**Correctness (the MD sharp edges).** The resume checkpoint is written on a frame boundary
(`checkpoint_steps = checkpoint_every_n_frames * traj_steps`, a multiple by construction — replaces
`max(traj_steps*10, 5000)` which could violate it) with the checkpoint reporter registered **before**
the trajectory/log reporters; on resume the trajectory + log are **truncated back to the checkpoint
step** and `sim.currentStep` restored before appending. This guarantees no duplicate/divergent frames
and keeps the XTC and log aligned and monotonic. Full rationale: GOTCHAS.md.

**Engine discriminator = a `resume_manifest.json` at the outdir root**, written by `run` at the start
of production (survives a kill), carrying `{engine, production_config}`. `resume` reconstructs the
exact frozen config from the manifest's `production_config` (NOT `config.json`, which is the
prepare-time `PrepConfig` and would inject silent engine defaults). Missing manifest ⇒ fail loud
(never infer the engine from which subdir exists); a different-engine `run` into the same outdir ⇒
fail loud (one outdir = one run/engine).

**Semantics.** `extend_to_ns` is an absolute target (re-passing it after a crashed resume continues,
not double-extends). Resuming an already-complete run with no `extend_to_ns` ⇒ fail loud; completion
is a **derived** query (manifest target vs checkpoint/`.nc` progress), **not** a written sentinel
file. Each resume segment writes its own `run_info.resume.<utc>.json` (provenance may differ from the
prepare-time environment) and a per-segment `run_report.resume.<utc>.json` (does not overwrite the
original report).

**Rejected / deferred.** A `run(..., resume=True)` kwarg (different entry shape — not single-purpose);
segmented `trajectory.partNNNN.xtc` files (truncate-to-checkpoint heals the append, so single-file is
safe). **Deferred (not built):** the exact same-hardware binary `.chk` fast path gated by an
environment fingerprint (portable State is the sole, sufficient resume point); a `should_stop` hook +
non-exception `is_complete(outdir)` for a future self-requeuing SLURM driver (one-line lifts if
auto-requeue is adopted). The output-contract block below changes accordingly: MD `production.chk`
→ `checkpoint_state.{A,B}.xml`.

## 2026-10-30 — configs round-trip to JSON; enums serialized by NAME

**Decision.** `PrepConfig` gains `to_json(path)` / `from_json(path)` (inherited by all engine
configs), and `protein_ff` / `water` serialize by enum **name**, not value. A saved config is both
the run's provenance and a reloadable input: `submit → MDConfig.to_json → sbatch → MDConfig.from_json
→ prepare/run`. `from_json` re-runs validation, so a bad field or FF pairing fails at load.

**Why.** `prepare()` already wrote `config.json`, but it could not be read back:
`WaterModel`'s value is the tuple `("tip3p", 3)`, which JSON stores as a list and pydantic refuses
to coerce back to the enum — so `model_validate_json` raised. config.json was therefore write-only.
Serializing enums by name (`"TIP3P"`, `"CHARMM36M"`) makes the file human-readable, editable, and
round-trippable. `before` validators still accept the legacy value forms (xml string, `[name, n]`
list) so older config.json files load.

**Rejected.** (1) Passing run parameters through SLURM `--export` env vars (the GROMACS idiom) —
stringly-typed, unvalidated until the job starts, and not the pipeline's native mechanism. (2)
Per-run argparse flags on the driver — still re-specifies the typed config as CLI strings. (3)
Changing `WaterModel`'s value off the tuple — invasive; the `(name, n_sites)` value is used across
solvation/packing. A localized field (de)serializer is the smallest correct fix.

**Scope note.** Nothing consumes config.json's enum fields programmatically
(`read_storage_directory` reads only `storage_directory`; `RunInfo` reads the live cfg object), so
the on-disk format change is safe.

## 2026-10-29 — cation-π NBFIX as an opt-in, CHARMM36m-only FF overlay (not a source-file edit)

**Decision.** The cation-π correction (Liu et al., JCTC 2020, 16, 6397) is a `bool` flag on
`PrepConfig` (`cation_pi_nbfix`, default off) that makes `build_system` load a rendered
`<LennardJonesForce>` fragment *after* the base CHARMM36m XML. The 6 pairs (Arg `NC2` + Lys
`NH3`, each × `CA`/`CAI`/`CPT`) live as a typed table in `cation_pi.py` in the **paper's units**,
converted to OpenMM units by a unit-tested function. Enabled only for `CHARMM36M`.

**Why a table-in-paper-units + render, not a static XML.** Keeps the paper's actual Table 1
numbers visible as the source of truth and makes the kcal→kJ / rmin→sigma conversion
tamper-evident (a unit test pins it) rather than freezing hand-computed floats. *Rejected:* a
static XML fragment (simpler, but the converted floats are unchecked); editing OpenMM's
installed `charmm36_2024.xml` (not version-controlled, clobbered on reinstall).

**Why opt-in + CHARMM36m-only.** It changes the force field, so it must be explicit and pinned
(reproducibility). The values are CHARMM36m atom-class LJ parameters with no AMBER equivalent;
`PrepConfig._check` raises `ForceFieldError` for any other FF (incl. plain C36 — the paper fit
C36m). *Note:* NBFIX is **global by atom class** — it deepens every Arg/Lys↔aromatic contact
(Arg–Phe/Tyr too, not just the target Arg–Trp). That is the paper's intent (uniform correction
validated on whole proteins); a per-residue-pair restriction is not expressible as NBFIX.

**Scope = the protein subset of Table 1.** Phe/Tyr are covered via the shared `CA` class (supp
Table S2: para-cresol C = `CA`; real Tyr = `CA`). *Excluded:* `CAB` (a toy atom type the paper
introduced only for its para-cresol free molecule — never present in a protein), TMA⁺
(`NG3P0`/`CG334` — quaternary ammonium, CGenFF ligands/ncAAs only), His (the paper never
parameterized imidazole). The CGenFF twin (supp Table S1) is the documented extension path for
non-natural residues; deferred (no ncAA in use).

**REST2 incompatibility made explicit.** `REST2Config._check_rest2` refuses `cation_pi_nbfix`
at config time: the correction is a CHARMM `CustomNonbondedForce`, whose REST2 solute-scaling is
unimplemented (`scale_solute` refuses CHARMM; see the 2026-10-06 entry + TODO). An early,
cation-π-specific error beats the deep runtime refusal, and it stays correct even if plain-CHARMM
REST2 lands later without cation-π-pair scaling. MD + T-REMD are unaffected (no force scaling).

## 2026-10-28 — conformational clustering default: DBSCAN (noise-aware) + HDBSCAN option; no frame cap

**The default `cluster_method` is now `dbscan` (was `hierarchical`), with an adaptive `min_samples`
density floor; `hdbscan` is added as an opt-in method; `cluster_max_frames` (silent auto-stride to
2000) is replaced by an explicit `cluster_stride` (default 1 = all frames).**

**Why.** Hierarchical average-linkage with a fixed distance cut (`fcluster(criterion="distance")`)
has no noise category: on a continuous/diffuse ensemble every pair separated by more than the cutoff
becomes a new cluster, so it over-fragments. A real unbound-peptide run
(`a83660d9…_openmm_md_unbound_300K_400ns_seed1`) produced **444 clusters** at the 2.0 Å cutoff
(design-like basin rank 167, 0.2%) — the signature of a floppy continuum, not discrete states.
DBSCAN on the design-superposed Cα coordinates keeps the *same* interpretable Å cutoff
(`eps = cutoff·√n_atoms`, valid because callers superpose before clustering) but sends the diffuse
tail to a **noise** label instead of singletons, and scales to all frames (no O(N²) matrix).

**The real fix is the density floor, not just the method.** DBSCAN with a *fixed* `min_samples`
(the old `=5`) would also fragment, because 5 frames is a vanishing fraction of a 40k-frame run.
So `cluster_min_samples=None` now resolves to the adaptive `max(10, 1.5%×frames)` — a region must
hold ~1.5% of the trajectory to be a state. (The 1.5%/floor-10 values are re-derived from the sibling
GROMACS `gromd-cluster`, which hit and solved this same problem; the constant is **copied, not
imported** — the two packages stay decoupled, and we do not commit to bit-identical output.)

**HDBSCAN as the opt-in deep-dive.** HDBSCAN extracts each cluster at its own density scale, so it
keeps a broad-but-populated state that a single DBSCAN `eps` would drop into noise (varying-density
basins). Rejected as the *default* because it has no single physical Å cutoff — its partition adapts
per run, which undercuts this module's whole job (a system-agnostic overview *comparable across
variants*, per `config.py`). It is the right tool for enumerating one design's true metastable states.
**Rejected QT/Daura–GROMOS** (the classic `gmx cluster` method) outright: it needs the full O(N²)
pairwise-RMSD matrix and has no noise category — the exact cost DBSCAN avoids.

**No silent frame cap.** The old path silently strided to ≤2000 frames (hierarchical's O(N²) limit).
DBSCAN/HDBSCAN cluster all frames by default. Hierarchical stays available but fails loud above
`MAX_HIERARCHICAL_FRAMES` (10k; a dense 40k² matrix ≈ 13 GB) telling the caller to raise
`cluster_stride` or use dbscan — never silent decimation.

**Consequences.** scikit-learn moves from the optional `[cluster]` extra to a **core dependency**
(`>=1.3`, which ships `sklearn.cluster.HDBSCAN`); `pyproject.toml` + `environment.yml` updated. The
population plot gained a noise bar and now takes percentages over all clustered frames (incl. noise),
matching the report JSON. **Open:** `design_like_fraction` is read off cluster populations, so a
design-like basin holding <1.5% of frames collapses into noise and under-reports retention — see TODO.

---

## 2026-10-22 — chain-break detection trusts geometry; SEQRES/REMARK are hints only

**Added `_assert_no_backbone_breaks` (`preparation/build.py`): a geometry-only guard that measures the
actual C(i)–N(i+1) distance between consecutive protein residues and fails loud above
`_PEPTIDE_BOND_MAX_NM` (0.2 nm; a real peptide bond is ~0.133 nm).** It runs after `addMissingAtoms`,
skips ACE/NME caps on the side they lack a backbone atom, and resets across non-protein residues. This
is now **the** chain-break check; the pre-existing `_assert_no_internal_gaps` (SEQRES-based, via
`PDBFixer.findMissingResidues()`) is kept only as a secondary early hint with residue-identity in its
message.

**Why geometry is the authority.** For "is residue i actually bonded to i+1," the coordinates are the
experimental observable; everything else is annotation about them, and all the annotations fail:
- **Author residue numbers are a label, not connectivity.** A renumbered construct can read 10,11,15,16
  with 11 peptide-bonded straight to 15 (cosmetic gap); continuous numbers can hide a real break;
  insertion codes (52,52A,52B) add residues with no number change. A numbering-gap heuristic is wrong in
  both directions — rejected.
- **SEQRES is often missing (stripped/processed files, computational models) or genuinely wrong**
  (author/construct sequence ≠ coordinates). A hint, not truth — rejected as the authority.
- **REMARK 465 is free-text metadata decoupled from the coordinates**: usually right on a pristine wwPDB
  deposition, but dropped or left stale by any tool that rewrites the file (relax, renumber, chain
  extraction), and never updated downstream. Must never be trusted over geometry.
- A C–N distance test has a huge safe margin (bond ≤ ~1.5 Å, gap ≥ ~2.5 Å), and two non-bonded backbone
  atoms can't sit ~1.3 Å apart without clashing — so geometry can't false-pass a real break.

**Motivating failure (see GOTCHAS).** A relaxed input with SEQRES stripped defeated `_assert_no_internal_gaps`
(findMissingResidues → `{}`), so OpenMM bonded chain B 46→50 across a 10 Å gap; initial PE still looked
normal (torn-bond penalty masked by solvent) and minimization silently corrupted the pose. The geometric
guard catches it (`chain B 46->50 (C-N 10.0 A)`), SEQRES or not.

**Scope / caveat.** Geometry decides *whether/where* a loop is missing; it cannot give the *identity or
length* of the missing residues — that comes from the **design record** (designed constructs) or
**UniProt cross-checked against the observed flanks** (native), never from SEQRES/REMARK on faith. The
check also needs the backbone atoms to exist: a CA-only/coarse model would pass silently and would need a
CA(i)–CA(i+1) fallback (~3.8 Å bonded) — out of scope, since a CA-only input fails `createSystem` anyway.

---

## 2026-10-22 — multistate per-replica equilibration default = 0 (off)

**`MultiStateProductionConfig.equilibration_ns` default flipped 0.5 → 0.0** (REMD + REST2); the
discarded per-replica equilibration (`sampler.equilibrate`, `remd/production.py` / `rest2/production.py`,
gated on `equilibration_iterations > 0`) is now **off by default**, so production starts straight from
the prepared pose (minimized + restrained-NVT + NPT-density-equilibrated, restraints released).

**Why.** This tool measures **drift from the starting designed pose** (RMSD/RMSF vs the input; see
objectives #1–3). `sampler.equilibrate` discards frames *in the simulation* — irreversibly throwing
away the early relaxation away from the design, which is exactly the signal of interest, and
contradicting the prepare/produce seam's own stated intent ("restraints release at start of production
so the trajectory captures drift from the input pose"). A burn-in window is better chosen **in analysis**
(one line, revisable) than baked in by dropping frames. The only thing equilibration buys is slightly
cleaner early exchange-acceptance statistics, also recoverable by discarding early iterations post-hoc.
The knob stays — set `equilibration_ns > 0` when you want that and don't care about the initial transient.

**Also removed the duplicate `equilibration_ns` re-declaration in `REMDConfig`** — it restated the
parent value and was why the default lived in two places; REMD now inherits it like REST2 already did.

**Cross-engine note.** The sibling GROMACS engine (`RELE_simulations/gromacs_REMD`) defaults
`EQUIL_NS=0.2` (shell fallback in `REMD-/REST2-gromacs.sbatch`), not 0 — it permits 0 but doesn't
default to it. So the two engines now differ on this default (openmm 0.0 vs GROMACS 0.2); pin
explicitly on both for a matched cross-engine comparison. (`HEAT_NS`/`nvt_equil_ns` still agree at 0.2.)

---

## 2026-10-21 — default FF = ff19SB/OPC; OPC is HARD-required for ff19SB (not just recommended)

**Default force field is now ff19SB + OPC** (`PrepConfig.protein_ff`/`water`), reversing the prior
ff14SB + TIP3P default. ff19SB (Tian et al. 2020, JCTC 16:528) is the modern Amber backbone model —
amino-acid-specific φ/ψ CMAP trained against QM surfaces with a solution correction — and it is the
right default for this tool's job (stability/flexibility of *designed* structures). The examples
(`examples/*/demo_*.py`) still pin ff14SB/TIP3P explicitly for fast demos; only the default moved.

**OPC is enforced as a hard requirement for ff19SB, not a soft warning.** ff19SB was jointly
parameterized with OPC; the paper tested TIP3P and found it degrades the model (the CMAP + side-chain
balance assumes OPC's dipole response). So ff19SB + any non-OPC water is a *correctness error*, not a
tunable choice. New `_REQUIRED_WATER = {AMBER19SB: OPC}` in `forcefield.py`; `resolve()` raises
`ForceFieldError` on a mismatch — the single pairing-validation point, same place that already
hard-locks CHARMM to its bundled TIP3P. `_RECOMMENDED_WATER` + the `config._check` soft warning stay
for FFs whose water is a genuine choice (e.g. ff14SB with TIP3P *or* OPC). Rejected keeping it a soft
warning: a mandatory pairing that silently proceeds is exactly the "looks plausible but wrong" failure
this project's fail-loud rule targets. (OPC ≠ OPC3 — OPC3 is a different 3-point model, also rejected.)

**REST2 needs no change for the new default.** `scale_solute` scales any `CMAPTorsionForce` by λ
(keyed on force *type*, not FF name), so ff19SB's QM-CMAP is tempered exactly like CHARMM36m's and
like the proper torsions. Verified: on a solvated ff19SB peptide, `scale_solute(..., lam=0.5)` scales
every CMAP map energy by 0.5 and does not trip the outside-solute refusal gate.

---

## 2026-10-21 — analysis = default-overview; clustering method policy; storage_directory model; run identity

**Analysis is an autogenerated default overview, not the final analysis.** System-agnostic defaults
so runs/variants compare directly; every parameter + method serialized into the report. Consolidated
into one pure `AnalysisConfig` (not a `PrepConfig` subclass — different axis).

**Clustering: hierarchical is the DEFAULT, DBSCAN is optional.** A fixed Å-RMSD cutoff is
interpretable ("within X Å = same conformation") and comparable across any system; DBSCAN's
`eps`/`min_samples` are density knobs that shift with sampling density (not comparable run-to-run)
and it clusters on common-frame RMSD (an upper bound on optimal pairwise), so it is kept only as an
opt-in method (scikit-learn, lazy import, `[cluster]` extra) for very large ensembles. Large
hierarchical runs **auto-stride** (recorded, not silent) rather than fail loud — an autogenerated
overview must not crash on any system. Clustering is anchored to the design (rep RMSD-to-design +
design-like-basin fraction), since the question for a designed structure is conformational
heterogeneity vs the design, not folding basins. *Reversed* an earlier "drop DBSCAN entirely" and an
earlier "add DBSCAN as a co-default".

**Run identity = the outdir** (path-based), made explicit (was implicit in the output-contract
draft). Rejected content-addressing (the sibling FragForge model): MD is stochastic with intended
replicates and the workflow is human-organized variant comparison; a recorded config hash can be
added later for dedup without changing the layout.

**storage_directory model — inverted from the sibling GROMACS pipeline (no code shared or
imported).** storage = the complete record; outdir = a lightweight working view. Huge raw files are
written **directly onto storage** (never *through* a symlink → no cross-fs `rename()`/EXDEV, no
reliance on how a reporter tests file existence) and symlinked into outdir; small files + all
analysis outputs + the processed trajectory are copied into both. Folder name =
`<outdir name>-<hash(abs outdir)>` — deterministic (re-findable; reruns reuse it) and collision-proof
across same-basename runs; `.source_outdir` `O_EXCL` marker records ownership (enables orphan
detection); `offload` fails loud on a pre-existing output (never silently deletes a trajectory);
`clean()` is marker-interlocked. `storage_directory` lives on `PrepConfig` (set once at produce
time) and analysis reads it back from `config.json`; producing without storage then analyzing with
it is a non-goal (avoids a split-brain record). Rejected the GROMACS choices: timestamped per-job
folder names (orphan accumulation; non-deterministic → analysis can't re-find), folder-level
stage symlinks (would hide small files from outdir), and a `finalize_outputs.sh`-style whole-tree
walk (replaced by a scoped `RunLayout.mirror`).

**Chain-integrity guarded at prep, not analysis.** `build._assert_no_internal_gaps` fails loud on an
internal missing-residue gap in a designed input, so analysis trusts topology chain counts and uses
**positional** chain correspondence (OpenMM prep preserves chain order; PDBFixer edits within chains
and `Modeller.addSolvent` appends solvent after the protein — unlike GROMACS, which renumbers).

## 2026-10-15 — REST2 scales torsions + nonbonded, NOT bonds/angles (the standard); straddle guard

**Decision.** The default REST2 solute scaling (`rest2/scaling.py`) scales: charges √λ, LJ ε·λ,
NonbondedForce exceptions (1-4), proper + improper torsions ·λ, and CMAP ·λ. It does **not** scale
solute **bonds or angles** (factor 1). Bonds/angles are recognized and passed through unscaled, not
dropped — so the unhandled-force refusal still protects against genuinely unknown forces.

**Why.** The prior `scale_solute` (Claude-generated, never validated against the field) scaled bonds
and angles. That is non-standard. Scaling bonds/angles appears only in NAMD's opt-in `soluteScalingAll`
and a teaching tutorial; it is absent from the original REST2 paper (Wang-Friesner-Berne 2011, which
reports scaling them does not aid sampling), plumed `partial_tempering`, NAMD's default, and femto
(Psivant's production OpenMM REST2, which `raise`s on bonds/angles). Verified three ways, including an
empirical diff of the sibling `gromacs_REMD` pipeline's generated per-replica topologies: only LJ,
charges, and dihedrals change across λ; bonds and angles are identical. Physical rationale: bonds and
angles are stiff coordinates near equilibrium; softening them does not lower the conformational
barriers REST2 exists to cross. *Rejected:* the full-E_pp "literal solute temperature" scheme (scale
everything) — legitimate but non-standard, no sampling benefit, and it would have forked AMBER and
CHARMM into different methods.

**Straddle guard.** `_scale_periodic_torsion` **refuses** (NotImplementedError) a torsion with some —
but not all — atoms in the solute. Such a cross-boundary bonded term should be √λ-scaled (the E_pw
cross factor), which is not implemented. This pipeline only tempers whole molecules (whole protein or a
whole chain), where no covalent torsion crosses the boundary, so the guard is inert for supported
selections. It fires only if the selection cut a covalent chain mid-molecule — a future sub-chain
selection, or an **inter-chain disulfide** linking the hot chain to a cold one (a real case for the
disulfide-constrained designs). *Rejected:* silently leaving straddling torsions unscaled (the prior
all-or-nothing behavior) — it would produce a quietly-wrong Hamiltonian the day sub-molecule selection
is added. Torsions are the only term needing the guard: bonds/angles are no longer scaled, NB
exceptions already apply per-atom √λ to mixed pairs, and CMAP has its own all-inside gate.

**Mechanism note.** femto implements this scheme with global-parameter offsets (one System, λ per
context). We keep **path A** (N pre-scaled Systems). The *scheme* (which terms, which powers) is
adopted from femto/the standard; the *mechanism* stays path A — independent choices (see 2026-10-06).

## 2026-10-10 — Fail-loud must not fire on valid input (λ round-trip assert removed)

**Decision.** Removed the `_ladder` assertion in `rest2/production.py` that re-derived an
explicit `lambdas` ladder via `lambdas_from_effective_temperatures` and compared with exact
float `==`. Validation of the λ ladder lives **once**, at the config boundary
(`REST2Config._check_rest2`: λ[0]=1, strictly descending, every λ in (0,1]); production trusts it.

**Why.** "Fail loud" means crash on *violated assumptions*, not on *floating-point identity of
a redundant round-trip*. `T_0 / (T_0 / λ)` reintroduces rounding (0.43 → 0.43000000000000005),
so the assert raised `AssertionError` on legitimate ladders — a fail-loud check that fails on
valid input is worse than none: it makes the documented `lambdas=` path unusable and surfaces
only at job-submit. General rule going forward: validate a value **once at the boundary** where
it enters, parse-don't-validate; don't re-assert a derived quantity's exact equality downstream.
If a downstream sanity check is genuinely wanted, use a tolerance (`np.allclose`), never `==` on
floats.

## 2026-10-06 — REST2: pre-scaled Systems (path A), fail-loud FF support, released-restraint gate

**Decision.** REST2 scales the solute force-field terms into **N independent pre-scaled
Systems** (one per λ), each a plain `ThermodynamicState` at the shared physical T_0
(kickoff-note "path A"). `scale_solute` handles the **AMBER force set + CMAP** and
**refuses** any other force touching the solute. Production is **NVT**, same as T-REMD.

**Why path A over CompoundThermodynamicState.** Path A matches the existing `_nvt_system`
shape and makes scaling a pure, energy-testable primitive (`scale_solute(system, indices,
lam) -> System`), validated by a λ=1 energy identity + direct parameter checks. Cost is N
Systems in memory, fine at our replica counts (REST2's whole point is *few* replicas).
*Rejected (for now):* `GlobalParameterState` + `CompoundThermodynamicState` (the scalable
"openmmtools way", one System with a λ global parameter) — more upfront force-authoring; a
later refactor if memory ever bites. openmmtools ships **no** turnkey REST2; `alchemy` is
for free-energy, not this.

**Why fail loud on non-AMBER forces.** The combining-rule trick (charge→√λ, ε→λ) gives the
correct solute-solute (∝λ) / solute-solvent (∝√λ) scaling *only* for a plain
`NonbondedForce`. CHARMM36(m) builds LJ as `CustomNonbondedForce`+`CustomBondForce` (NBFIX)
and impropers as `CustomTorsionForce`; scaling those correctly is real work. Rather than
silently mis-scale (the exact "looks plausible but wrong" failure the project forbids),
`scale_solute` raises `NotImplementedError` on any unrecognized solute force. AMBER REST2
works today; CHARMM REST2 is an explicit TODO. The **CMAP safety gate** is the same spirit:
CMAP maps are shared across torsions, so scaling them is allowed only after asserting every
CMAP torsion lies inside the solute (else a one-chain solute would wrongly temper the rest).

**Released-restraint gate.** `prepare()` bakes the position restraint released (global
`k`=0) into the handed-off System, so a `CustomExternalForce` is present in production. It
is an *external* restraint, not a force-field solute term, so it is **not** λ-scaled — but
`scale_solute` asserts it is released (all global-parameter defaults 0) before carrying it,
because an active restraint under REST2 would bias every replica identically. *Rejected:*
blindly ignoring any `CustomExternalForce` (would hide an accidentally-live restraint).

**Shared config base + factored primitives.** `REMDConfig`/`REST2Config` now subclass
`MultiStateProductionConfig(PrepConfig)` (shared production/exchange fields + derived
iteration counts); the barostat strip moved to `lib/system_edits.strip_barostats`; the
exchange/demux analysis moved to `analysis/multistate.py` with `remd.py`/`rest2.py` as thin
subdir+`.nc` wrappers. Decompose eagerly — one primitive each, composed by both engines.

**λ ladder = T_0/T_eff on a geometric effective-T ladder.** State 0 is λ=1 (the unscaled,
physical ensemble — the only REST2 rung that is a real Boltzmann sample at the design
temperature, and the one analysis demuxes). Reuses the T-REMD geometric ladder on effective
temperatures. *Verified:* solute scaling checked on a real solvated AMBER system in
`test_rest2.py` (protein atom √λ·q/λ·ε; water byte-identical).

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
    final_state.xml    # saveState (portable end state)
    checkpoint_state.{A,B}.xml  # ping-pong portable resume checkpoints (see 2026-10-30 resume entry)
    run_report.json    # QC: timing, T/density/energy stats, stability flag
    run_report.resume.<utc>.json  # per-resume-segment QC (does not overwrite run_report.json)
  run_info.resume.<utc>.json      # per-resume-segment env provenance (outdir root)
  resume_manifest.json            # {engine, production_config} — authoritative, for resume (outdir root)
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
