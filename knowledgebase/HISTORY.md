# HISTORY — openmm_pipelines

Dated session log, newest on top. One entry per working session/day. Terse and
factual.

---

## 2026-10-31 — prep: never build missing residues (silent terminal-extension bug)

- **Found via the Ellen WW-domain MD runs:** PDBFixer was building 21 **terminal** residues into an
  idealized extended tail (`Refine43_chainA_only_WW34`: SEQRES 98 vs 77 coordinate residues) -> ~143 Å
  solvent box instead of ~57 Å. `_assert_no_internal_gaps` only catches internal gaps, so terminal runs
  passed silently. Also explains the large unrestrained-minimization drift (un-superposed Cα RMSD
  ~2.75 Å, one residue ~10 Å) and plausibly the rare NVT NaN seen earlier.
- **Fix:** `build_system` calls `_drop_missing_residues(fixer)` after the internal-gap guard — clears
  `fixer.missingResidues` (no residues built; internal gaps still fail loud) and warns with the count.
  Missing atoms/H within existing residues are still completed. Validated: both WW34 inputs build 77
  residues, 56.8 Å dodecahedron box, warning naming the 21 skipped. New `test_build.py` unit tests;
  full suite green (148 passed).
- **Box shape verified, not a bug:** the pipeline produces a correct rhombic dodecahedron
  (a=b=c=56.8 Å, angles 60/60/90). A viewer showing a "rectangular prism" renders the triclinic unit
  cell, whose a–b face is square (γ=90°); no box code changed.
- **NVT NaN (context):** the one-off R448A + cation-π NaN did **not** reproduce (0/17 builds survived
  50k NVT steps, incl. pe₀ up to 1.5×10¹¹ kJ/mol). `build_system` is nondeterministic even with stdlib
  `random` seeded in-process (not a seedable RNG — likely atom-object set iteration in `addSolvent`), so
  no clean seeded-solvation fix; rare (<~6%), re-run clears. The terminal tail was the real, fixable bug.
- Extended-structure jobs cancelled; their outdirs + storage scrubbed.

## 2026-10-30 — resume / extend an interrupted production run (MD, REMD, REST2)

- **Added `resume(outdir, *, extend_to_ns=None, platform=None)`** (top-level, `resume.py`) — continues
  an interrupted (walltime/crash) run, or extends it to a larger target, via the **documented**
  OpenMM / openmmtools restart mechanisms (not hand-rolled). Library-only: re-submitting the SLURM
  job is the driver's job; the library reads disk state only (no SLURM coupling). Planned with the
  architect agent (its correctness fixes folded in before coding).
- **Engine dispatch via `resume_manifest.json`** at the outdir root (written by `run` at production
  start; `lib/run_layout.py` `write_manifest`/`read_manifest`). Carries `{engine, production_config}`;
  `resume` reconstructs the exact frozen config from it (NOT `config.json`, which is the prepare-time
  `PrepConfig` → would inject silent engine defaults). Missing manifest / different-engine rerun ⇒
  fail loud.
- **MD (`md/production.py`)** — split `run`/`resume` over shared `_make_simulation`/`_add_reporters`.
  Resume point is a **portable State XML written ping-pong to two files** (`_PingPongStateReporter`;
  torn-write safe — OpenMM's `CheckpointReporter` overwrites one file in place). On resume: load the
  newest valid State, set `sim.currentStep`, **truncate `trajectory.xtc` + `production.log` back to the
  checkpoint step**, then `append=True`. Checkpoint cadence is now frame-based
  (`checkpoint_every_n_frames=10` → `checkpoint_steps = n*traj_steps`, a frame multiple by
  construction; fixes the latent `max(traj_steps*10, 5000)` misalignment). Dropped the binary
  `production.chk` (replaced by the ping-pong State; exact `.chk` fast path deferred).
- **REMD/REST2 (`multistate_resume.py` + each `production.py`)** — reopen the existing `.nc` via the
  new `RunLayout.existing()` (never `offload`/`create`, which refuse existing files), then openmmtools
  `ReplicaExchangeSampler.from_storage` + `run()` (finish) / `extend()` (extend). `continue_to` returns
  the sampler so the caller keeps it alive through the QC report — **the sampler closes its reporter on
  GC** (`__del__`), which otherwise makes `_run_report` read a closed `.nc` (`NetCDF: Not a valid ID`).
- **Semantics**: `extend_to_ns` absolute (idempotent re-resume); already-complete + no extend ⇒ fail
  loud (completion is a derived query, no sentinel file); per-segment `run_info.resume.<utc>.json` +
  `run_report.resume.<utc>.json`.
- **Tests** `tests/test_resume.py` (9): config round-trip (all engines), `existing()` + manifest
  guards (fast); MD extend continuity (frame count == analytic, log↔xtc aligned, monotonic steps, no
  dups), MD resume-to-planned via manifest, MD fail-loud guards, MD resume with `storage_directory`
  (symlink resolution in truncate+append), REMD extend via `from_storage`, REMD resume-to-planned via
  the `run()` branch (built a partial 5-of-10 `.nc`) (slow).
- **Verified**: resume 9/9; full non-slow suite 125/125; `test_md.py`+`test_remd.py`+`test_rest2.py`
  6+3 (fresh `run()` + analyze + storage offload unregressed). Docs: README (resume quickstart +
  `run_or_resume` resubmit idiom + output layout; fixed a stale `from_disk(outdir)` → `/prepared`) and
  `examples/MD_demo/resume_demo.{py,sbatch}`; DECISIONS (resume contract + output-contract update),
  GOTCHAS (append/truncate ordering, `existing()`+`from_storage`, sampler-GC-closes-reporter), TODO,
  CLAUDE.md index.
- **Deferred** (noted in TODO): exact same-hardware `.chk` fast path (fingerprint-gated); a
  `should_stop` hook + non-exception `is_complete(outdir)` if a self-requeuing SLURM driver is adopted.

## 2026-10-29 — cation-π NBFIX correction for CHARMM36m (Liu et al. 2020)

- **Added `cation_pi.py`** — pinned the paper's Table 1 protein subset as a typed table in the
  paper's own units (6 pairs: Arg `NC2` and Lys `NH3`, each × aromatic `CA`/`CAI`/`CPT`) and a
  tested conversion to OpenMM NBFixPair units (`sigma=rmin/2^(1/6)/10`, `eps=kcal·4.184`).
  `render_nbfix_xml()` emits a `<LennardJonesForce lj14scale="1.0">` fragment.
- **Not** a reimplementation of the paper — we consume its *output* (the NBFIX LJ values). The
  paper's method is a pure lookup table; the scientific core is 6 rows.
- **Verified the atom-class mapping** against OpenMM's shipped `charmm36_2024.xml` templates
  (not assumed): Arg NE/NH1/NH2→`NC2`; Trp CD1/CZ3/CH2→`CA`, CE3/CZ2→`CAI`, CE2/CD2→`CPT`;
  Phe/Tyr rings→`CA` (so both ride the `*-CA` rows); **Lys NZ→`NH3`** (so the paper's NH4+ rows
  apply to Lys). Excluded: `CAB` (para-cresol toy atom type, never in proteins — supp Table S2
  shows para-cresol C = `CA`), TMA⁺ `NG3P0`/`CG334` (CGenFF ligand/ncAA only), His (unparameterized).
- **Wiring:** `PrepConfig.cation_pi_nbfix: bool = False`; `_check` refuses it unless
  `protein_ff is CHARMM36M`. `build_system` appends the fragment (via `io.StringIO`) *after* the
  base FF so the pairs merge into CHARMM's `CustomNonbondedForce`. One `ForceField` feeds both
  `addSolvent` and `createSystem`; switchDistance already flows through (no extra work).
- **REST2 guard:** `REST2Config._check_rest2` refuses `cation_pi_nbfix` at config time (the
  correction is a CHARMM `CustomNonbondedForce`, unscalable by `scale_solute`; see DECISIONS/TODO).
  REMD + MD are unaffected — T-REMD runs it today (no force scaling).
- **Tests:** `test_cation_pi.py` (table/conversion/render + a real-ForceField check that the
  fragment enlarges the NBFIX pair table) + guard tests in `test_config.py` / `test_rest2_config.py`.
- **Open:** CGenFF twin (supp Table S1 — `NG2P1`/`NG3P0`/`CG334`/`CG2R61`/`CG2RC0`) for
  non-natural/CGenFF-typed residues; a drop-in row append, deferred until an ncAA is in use.

## 2026-10-30 — cation-π runtime-wiring validation (energy controls)

- **Closed the end-to-end gap** the generator-table test could not see: that `createSystem`
  actually applies the overridden LJ to real atoms. Added two `test_cation_pi.py` tests building
  a structure in **vacuum + NoCutoff** with and without the fragment (base vs. with-fragment are
  byte-identical except the overridden pairs, so ΔE isolates the override):
  - **Positive control** (`1a22-fixed.pdb`, HGH — many real Arg/Lys↔aromatic contacts): energy
    changes, bounded (`1 < |ΔE| < 1e4` kJ/mol), not noise, no blow-up. Measured **ΔE = −20.8
    kJ/mol (attractive**, the physically expected direction — deeper cation-π well).
  - **Negative control** (`helix_fusion.pdb`, poly-Ala — no cation, no aromatic): energy
    **unchanged**, proving the override touches only cation↔aromatic pairs, never global LJ.
- Un-minimized single-point absolute energies are ~1e11 kJ/mol (crystal-pose + added-H clashes);
  irrelevant — they cancel in the difference. Full suite still green.
- **Not done (user's call):** a full comparative production run (with/without correction on a real
  Trp–Arg design) — that validates the *science* on a given protein, not the implementation; the
  paper already validated the physics. Note the override is global-by-class, so such a run shifts
  *every* cation↔aromatic contact, not only the one of interest.
- **Comparative runs now being set up** (external to this repo, in `Ellen_WW_simulations/submit_jobs`):
  WW34 WT + R448A on CHARMM36m, each with cation-π ON and OFF (4 jobs/seed), via the config-file
  loader below. R448A removes the Arg (built-in control); ON-vs-OFF isolates the correction at
  fixed FF.

## 2026-10-30 — config-file loader for configs (to_json / from_json)

- **Added `PrepConfig.to_json(path)` / `PrepConfig.from_json(path)`** (inherited by MD/REMD/REST2
  configs). The pipeline had **no** way to reload a saved config: `prepare()` wrote `config.json`
  (`model_dump(mode="json")`) but `model_validate_json` on it **failed** — `WaterModel`'s value is a
  tuple `("tip3p", 3)`, which JSON stores as a list and pydantic cannot coerce back. So config.json
  was write-only provenance; now it round-trips.
- **Enums serialize by NAME** (`field_serializer` + `before` `field_validator` on `protein_ff` and
  `water`): config.json now reads `"protein_ff": "CHARMM36M"`, `"water": "TIP3P"` (human-readable +
  reloadable). The validators also accept the legacy value forms (xml string / `[name, n]` list), so
  older config.json files still load. `from_json` re-runs full validation (bad FF pairing fails loud).
- **Why:** setting up the WW-domain cation-π MD sweep, the clean way to specify per-run parameters
  is a saved typed `MDConfig`, not env-var `--export` smuggling (the GROMACS idiom). The submitter
  builds + validates + serializes one config per run; the driver does `MDConfig.from_json(path)`.
- **Tests:** `test_config.py` — round-trip preserves config + enums-by-name; legacy value-form loads;
  `from_json` revalidates (illegal combo raises). Nothing else reads config.json's enum fields
  (`read_storage_directory` takes only `storage_directory`; `RunInfo` uses the cfg object), so the
  format change is safe.

## 2026-10-28 — clustering overhaul: DBSCAN default, adaptive min_samples, HDBSCAN option, no frame cap

- **Diagnosed** the over-clustering report: an unbound-peptide run clustered into **444** states at
  2.0 Å (`…/a83660d9…_openmm_md_unbound_300K_400ns_seed1/analysis/clustering/`). Root cause: the
  default hierarchical distance-cut has no noise category, so a continuous ensemble fragments.
- **Verified on that run (all 40k frames, min_samples=600=1.5%):** DBSCAN → **8 clusters**, 61% noise,
  design-like frac 0.144; HDBSCAN → **3 clusters**, 13% noise, one dominant basin of 27511 frames
  (69%), design-like frac 0.688. The DBSCAN/HDBSCAN gap is the varying-density case in the flesh —
  HDBSCAN keeps the broad near-design basin that DBSCAN's tight 2 Å eps scatters into noise.
- **`analysis/config.py`**: default `cluster_method` `hierarchical → dbscan`; added
  `ClusterMethod.HDBSCAN`; `cluster_min_samples` now `int | None` (None = adaptive); replaced
  `cluster_max_frames=2000` (silent auto-stride) with `cluster_stride=1` (explicit, all frames).
- **`analysis/clustering.py`**: rewrote. Module constants `MIN_SAMPLES_FRAC=0.015`,
  `MIN_SAMPLES_FLOOR=10` (re-derived from GROMACS `gromd-cluster`, not imported),
  `MAX_HIERARCHICAL_FRAMES=10_000`. `_resolve_min_count` (adaptive density floor, used as DBSCAN
  `min_samples` / HDBSCAN `min_cluster_size`), `_apply_stride`, `_hdbscan_labels` (copy=True so it
  can't mutate `traj.xyz`). Hierarchical now fails loud above the frame cap. sklearn imported at top
  (no longer lazy/optional).
- **`analysis/plots.py`**: `plot_cluster_populations` gained a noise bar; percentages over all
  clustered frames (incl. noise) to match the JSON.
- **Deps**: scikit-learn moved from `[cluster]` extra to core (`>=1.3`) in `pyproject.toml` +
  `environment.yml`; installed into `openmm_env` (1.9.1).
- **Tests** (`tests/test_analysis.py`): replaced the 3 clustering tests with 6 — DBSCAN-default/all-
  frames, adaptive-min_samples-sends-tail-to-noise, HDBSCAN splits blobs, hierarchical opt-in +
  fails-loud-over-cap, explicit stride mapping. Added `_two_blob_traj` / `_blob_plus_outliers`
  helpers. Full non-slow suite: **107 passed**.
- Rationale + rejected alternatives (HDBSCAN-as-default, QT/Daura) in DECISIONS. Open item
  (`design_like_fraction` vs sub-1.5% basins → noise) in TODO.

## 2026-10-22 — per-replica equilibration default off (equilibration_ns 0.5 → 0.0)

- **`MultiStateProductionConfig.equilibration_ns` default 0.5 → 0.0** (`config.py`), so the discarded
  per-replica `sampler.equilibrate` is off by default for both REMD and REST2; production starts from
  the prepared pose and the trajectory captures the full drift from the designed input (the analysis
  signal). Burn-in is now an analysis-time choice, not baked in. Rationale in DECISIONS.md.
- **Removed the duplicate `equilibration_ns = 0.5` in `REMDConfig`** (`remd/config.py`) — it restated
  the parent default; REMD now inherits it like REST2 does.
- Verified both engines resolve `equilibration_ns=0.0` / `equilibration_iterations=0`. Config tests
  (`test_remd_config`, `test_rest2_config`, `test_config`) still pass (28) — they pin values explicitly
  and test the derivation, not the default.
- Indexed this repo with CodeGraph (`codegraph init`) at the user's request.
- Noted: sibling GROMACS engine defaults `EQUIL_NS=0.2` (not 0) — the two now differ on this default.

---

## 2026-10-21 — ff19SB/OPC is now the default; OPC hard-enforced for ff19SB

- **Default FF flipped to ff19SB + OPC** (`config.py` `PrepConfig.protein_ff`/`water`); was
  ff14SB/TIP3P. Flows down to MD/REMD/REST2 configs. Examples still pin ff14SB/TIP3P explicitly.
- **OPC hard-required for ff19SB** (`forcefield.py`): new `_REQUIRED_WATER = {AMBER19SB: OPC}`;
  `resolve()` now raises `ForceFieldError` for ff19SB + any non-OPC water (promoted from the
  `config._check` soft warning, which stays for FFs whose water is a real choice).
- **Confirmed (no code change) REST2 handles ff19SB CMAP**: `scale_solute._scale_cmap` scales any
  `CMAPTorsionForce` by λ — verified a solvated ff19SB peptide's CMAP maps scale by λ. Generalized the
  CMAP comments in `rest2/scaling.py` (was phrased CHARMM36m-only).
- Verified `amber19/protein.ff19SB.xml` + `amber19/opc.xml` ship in OpenMM 8.6.1 and
  `build_system` produces the CMAPTorsionForce + OPC M-site virtual sites.
- Tests: updated `test_config` defaults + soft-warning case, added `test_required_water_enforced`;
  added `test_ff19sb_requires_opc`; pinned the `tip3p_build` fixture to ff14SB/TIP3P. Fast suite:
  100 passed, 1 skipped, 15 slow deselected.

---

## 2026-10-21 — analysis output upgrades: clustering overhaul, interchain metrics, storage_directory

Reworked the analysis/output layer. Framing decided with the user: the analysis path is an
**autogenerated, re-runnable default overview** (system-agnostic defaults, comparable across runs,
every parameter + method serialized into the report), not the final analysis.

- **AnalysisConfig** (`analysis/config.py`, pure pydantic) carries every analysis knob, serialized
  into each report. Composers now `analyze(outdir, config)` / `analyze_remd/rest2(outdir, config)`.
- **Clustering** (`analysis/clustering.py`): `ClusterResult`; hierarchical default (scipy, Å-RMSD
  cutoff, auto-strides above `cluster_max_frames`, stride recorded) + optional DBSCAN (scikit-learn,
  lazy import, `[cluster]` extra). Writes representative PDBs + populations + cluster plots; anchored
  to the design (rep RMSD-to-design + design-like-basin fraction). Fixed noise-in-`n_clusters` and
  frame-index-through-stride bugs.
- **Interchain** (`analysis/interchain.py`, new): min distance / contact count / COM distance (all
  minimum-image) + per-chain fold-retention RMSD; auto-run at ≥2 protein chains, positional chain
  correspondence. Computed BEFORE superposition (which rotates coords but not the box).
- **Chain-integrity guard** moved to prep (`build.py _assert_no_internal_gaps`) — fail loud on an
  internal missing-residue gap so analysis can trust chain counts.
- **storage_directory** (`PrepConfig` + `lib/run_layout.RunLayout`): offload huge raw files
  (trajectory.xtc, .nc + checkpoint) to a deterministic, collision-proof storage folder with
  symlinks in outdir; mirror small files + the analysis overview + processed traj into storage
  (storage = complete record). `O_EXCL` ownership marker, fail-loud-on-existing, interlocked
  `clean()`. Wired into md/remd/rest2 production + both analysis composers.
- **Interactive HTML report** (`analysis/report.py`, new): `write_html_report` builds a
  self-contained Plotly `report.html` (plotly.js embedded once → offline-viewable, archivable) from
  the composer's `summary` + `arrays` — drift / RMSF / SS / clustering / interchain / exchange
  figures + a key-metrics table. On by default (`AnalysisConfig.html_report`); plotly is an optional
  `[report]` extra (lazy import, fail-loud with install hint). Written per run (plain MD
  `report.html`; multistate `<prefix>_report.html`) and mirrored into storage.
- **Analysis output layout:** clustering artifacts (representative PDBs + cluster plots) now go in
  an `analysis/clustering/` subdirectory; multistate analysis writes to top-level `outdir/analysis/`
  (was the nested `outdir/<engine>/analysis/`) — consistent with plain MD; the `.nc` stays under
  `outdir/<engine>/`.
- **Per-residue DSSP timeline** (`dssp.dssp_timeline` + `plots.plot_dssp_timeline`): per-residue ×
  per-frame SS map (full mdtraj H/G/I/E/B/T/S/coil/NA codes, `DSSP_CATEGORIES` palette) written as
  `analysis/dssp_timeline.png` + an interactive heatmap in the report. Figure size bounded on both
  axes (`min(max(8,n_frames·0.02),16) × min(max(4,n_res·0.05),16)`, `aspect="auto"`) — borrowed
  from the sibling GROMACS pipeline's dssp plot so it stays legible for any system/frame count.
  (mdtraj's native DSSP does NOT emit DSSP-4 PP-II `P` / `=`; would need mkdssp 4.)
- Deps: scikit-learn (optional `[cluster]` extra) + plotly (optional `[report]` extra), both in
  environment.yml; `RunInfo` records scipy + scikit-learn versions.
- Ran all three demos' analysis on CPU (REMD/REST2 analysis-only on existing `.nc`): interactive
  reports generated; REST2 λ=1 ensemble resolved 4 clusters.
- Tests: new `test_run_layout.py` (9) + `test_md.py` (slow: MD + analyze + storage); extended
  `test_analysis.py` (interchain, clustering, dbscan); storage assertions in
  `test_remd.py`/`test_rest2.py`. Green on CPU: 97 non-slow + MD/REMD/REST2 slow.

## 2026-10-21 — a99SB-disp OpenMM port: answered the external validation review (tools/convert_a99sbdisp)

`claude science` reviewed `validation_report.md` (see `plans/response_to_validation_report.md`) and
raised corrections + 2 harness questions + a reference-value table + 10 follow-up items. All closed;
report rewritten (new §13 follow-up, §14 artifacts; §5/§6/§8/§11/§12 corrected).

- **Torsion-suppression semantics proven directly** (`prove_torsion_suppression.py`): a specific k=0
  override makes OpenMM's most-specific matcher emit ZERO terms on that quartet (generic suppressed);
  the fix no longer rests on an inferred energy match.
- **7 C-terminal Tier-1 flags fully explained + dual-ref confirmed:** all the SAME one-term artifact —
  Amber `tleap` stacks the generic `X-CT-CT-X` (n=3, k=0.6508) on the C-terminal Cγ-Cβ-Cα-C quartet;
  GROMACS tpr (grompp, `verify_terminal_flag_gmx.py`) carries only the specific term → candidate is
  correct. Mirror image of the original internal-χ bug. (Earlier "Amber adds a terminal term" hunch
  was right in spirit; now pinned to the exact quartet + term, and it's NOT a term-set change for CYS
  mid/N-term — only C-term.)
- **Harness hardened (Q1):** `tier1_compare.torsion_diff` splits term-SET/count mismatch (never
  benign) from pure value diff; `tier1_run` routing can no longer excuse a nonzero term-set mismatch
  via the energy heuristic. **Symmetric exclusion (Q2)** proven.
- **Ions added:** `a99SBdisp_ions.xml` (11 ions, LJ in `<LennardJonesForce>`). Na⁺/Cl⁻ dual-confirmed
  (Amber tleap == GROMACS, JC/TIP4P-Ew); K/Mg/etc. follow the validated GROMACS port and DIFFER from
  the shipped Amber leaprc's legacy parm99 values (upstream inconsistency — see GOTCHAS).
- **New standing checks:** `reference_values.py` (16 source-derived asserts), `tier6_regression_guards.py`
  (strict term-set guard + k=0-shadow survival + negative control), `test_fix_zero_periodicity.py`
  (4 unit tests), `tier7_ions.py`, `tier8_coverage.py` (programmatic), `tier9_boundary_disulfide.py`,
  `tier4b_discriminators.py` (NVE drift with uncertainty + discriminators; the unprovenanced "5× bar"
  removed). Still open: LYN Amber-vs-GROMACS disagreement; REST2 scaling layer (next).

---

## 2026-10-15 — REST2 solute scaling corrected to the standard (bonds/angles unscaled) + straddle guard

Investigated what REST2 should scale, because the existing `scale_solute` (Claude-generated,
never validated) scaled solute **bonds and angles** — and found that is non-standard.

- **Evidence (three independent lines, all agree): do NOT scale bonds/angles.**
  1. Literature + implementations: original REST2 (Wang-Friesner-Berne 2011), plumed
     `partial_tempering`, NAMD default (bonds/angles only under `soluteScalingAll`), and
     femto (Psivant, production OpenMM REST2 — `raise NotImplementedError` on bonds/angles).
  2. The sibling `RELE_simulations/gromacs_REMD` repo's `docs/rest2-validation-brief.md`
     states bonds/angles = factor 1, citing Wang 2011 ("scaling them does not help sampling").
  3. **Empirical:** diffed that pipeline's generated per-replica topologies (λ=1 vs most-scaled)
     — only `[atomtypes]` (LJ), `[atoms]` (charges, exactly the 433 protein atoms), and
     `[dihedrals]` (torsions) change; `[bonds]` and `[angles]` are byte-identical.
- **Change:** `rest2/scaling.py` — bonds/angles no longer scaled; now recognized-but-pass-through
  (`HarmonicBond/AngleForce` → `continue`) so the unhandled-force refusal still fires for genuinely
  unknown forces. Removed `_scale_harmonic_bond`/`_scale_harmonic_angle`. Still scaled (all standard):
  charges √λ, LJ ε·λ, NB exceptions (1-4), proper+improper torsions ·λ, CMAP ·λ.
- **Straddle guard:** `_scale_periodic_torsion` now **refuses** (NotImplementedError) a torsion with
  some-but-not-all atoms in the solute. Such a cross-boundary term would need √λ (E_pw) scaling,
  which is unimplemented; whole-molecule selection (whole protein / whole chain — all this pipeline
  does) guarantees none. Catches a future sub-chain selection AND an inter-chain disulfide linking a
  hot chain to a cold one (relevant to the disulfide-constrained designs). Torsions are the only spot
  that needed it: bonds/angles no longer scaled, NB exceptions already do per-atom √λ, CMAP has its
  own all-inside gate.
- Tests updated (`test_bonds_and_angles_not_scaled`, `test_straddling_torsion_refused`); **85/85
  non-slow pass**. Slow GPU tests not run (no GPU on dev node) — they assert only on NonbondedForce,
  so unaffected.
- Wrote `knowledgebase/REST2_scaling_open_question.md` (neutral, bias-stripped) to feed to Claude
  Science for a literature cross-check.

**Still deferred / follow-ups:** (a) CHARMM36m REST2 — the `CustomNonbonded`/`CustomBond`/
`CustomTorsion` (NBFIX LJ / 1-4 LJ / impropers) handlers are still refused loudly; design is worked
out (per-particle √λ multiplier on the NBFIX force) and grounded in a real force-dump, not yet
implemented. (b) The AMBER REST2 demo ladder (11 rungs, ~0.42 acceptance) was tuned under the OLD
scheme — fewer terms scaled now → acceptance will rise → needs a re-tune on a GPU node. (c) a99SB-disp
needs a trustworthy OpenMM conversion first (openmmforcefields does NOT ship it; authors' repo has only
Amber/Gromacs/Desmond). (d) Adopt the gromacs_REMD validation brief's gates 0-7 as the REST2 validation
harness.

## 2026-10-10 — Architecture review + fixes (HIGH/MEDIUM/LOW)

Ran a principal-architect pass over the whole package. Verdict: sound — prepare/produce
seam clean, pure/OpenMM separation honored, fail-loud strong (zero bare except), demux
walker/slot indexing correct. Fixed the findings:

- **HIGH — spurious λ round-trip assert (`rest2/production.py`).** `_ladder` re-checked an
  explicit `lambdas` ladder with exact-float `==` against `lambdas_from_effective_temperatures`;
  the round-trip injects float error (e.g. 0.43 → 0.43000000000000005) so it `AssertionError`ed
  on *valid* input, discoverable only at job submit. Deleted — `REST2Config._check_rest2`
  already enforces λ[0]=1, strictly descending, (0,1]. (See DECISIONS, fail-loud overreach.)
- **MEDIUM — reporter leak.** `remd`+`rest2` `run()` now wrap create/equilibrate/run in
  `try/finally: reporter.close()`, mirroring `analysis/multistate`. A crashed run no longer
  orphans the `.nc` (which then trips openmmtools' "cowardly refusing to overwrite" on rerun).
- **MEDIUM — DSSP `NA` dropped.** `ss_fractions` silently omitted mdtraj's `NA` code (chain
  termini, ACE/NME caps — common here), so H+E+C ≠ 1 and coil was under-reported. Now reports
  `na` as its own fraction + asserts codes ⊆ {H,E,C,NA}; `plots.py` plots it.
- **LOW — platform/precision dedup.** New `lib/platforms.py` (`is_gpu`/`gpu_props`/
  `precision_label`); replaced 4 inline copies (prepare, md/production, context_cache, reports).
  `context_cache` keeps `None`-for-CPU (set_platform's validator short-circuits on None, else
  builds a throwaway probe Context on `{}`).
- **LOW — fail-loud tightening.** `reports.write_json` `default=str` → `_json_default` (Path
  only; `TypeError` otherwise). `md._col` now asserts exactly one column match (was silent
  `matches[0]`). Real type hints across the seam (`prepare`/`build_system`/`*.run`) via
  `TYPE_CHECKING`.
- Refreshed `CODEBASE_SHAPE.md`: dropped the stale "design sketch only, no .py yet" header;
  `[soon]`/`[stub]` → `[BUILT]`; added `platforms.py`.

Verified: 84/84 non-slow tests pass; all edited modules import. The 11 slow tests are live
GPU runs (need a GPU node) — not run this session (no GPU on the dev node).

## 2026-10-07 — REST2 demo sbatch + ladder sizing + solute-mask verification

Added `examples/REST2_demo/demo_rest2.sbatch` (mirrors the REMD demo submitter:
`mit_normal_gpu`, one L40S, 30 min, `$SLURM_SUBMIT_DIR`) and ran it on a GPU node.

- **Ladder sizing.** The 6-rung demo (300→450 K effective) gave ~0.11 mean neighbour
  acceptance — confirmed real (not statistics) by a 500 ps / 501-iteration rerun that
  reproduced it (0.107, flat). Too sluggish (`subdominant_eigenvalue` ~0.98). A √p-scaling
  estimate (halving spacing takes p→√p) predicted ~0.33 at 11 rungs; the 11-rung run gave
  **0.42 mean, flat across the ladder** (`subdominant_eigenvalue` 0.956). Bumped the demo to
  11 rungs. 9–10 would also land in the 0.2–0.4 band; 6 is too few.
- **Why REST2's rung count looked close to REMD's** (user flagged it as suspicious): pure
  range artifact. REMD demo = 14 rungs over a *narrow* 300→360 K (4.6 K/gap, whole box);
  REST2 = 6 rungs over a *wide* 300→450 K effective (30 K/gap). REST2 tolerates 6.5× coarser
  spacing precisely because it tempers only the solute — the expected solute-tempering payoff,
  not a broken mask.
- **Solute mask verified on the REAL prepared system** (`demo_rest2_out/prepared/system.xml`
  + `topology.cif`, no re-prepare). Of 5414 atoms, exactly **433 (8%) are solute and get
  scaled** (q→√λ·q, ε→λ·ε, σ fixed), and all **4981 solvent atoms are byte-identical** —
  particles-changed == solute-count, so the mask neither leaks into solvent nor misses solute.
  All **217 solute hydrogens** are included. The low 10.1 atoms/residue was just because
  `helix_fusion.pdb` is a 43-residue **poly-alanine** helix, not a selection bug. (Scratch
  check script, not committed; the synthetic-system equivalent is already in
  `tests/test_rest2_scaling.py` + `tests/test_rest2.py::test_solute_is_actually_scaled_on_real_system`.)
- **New gotcha:** openmmtools refuses to overwrite an existing `.nc` ("cowardly refusing" —
  the git/Perl idiom, from `multistatesampler.py:589`); a rerun needs a clean outdir. Logged
  in GOTCHAS + CLAUDE.md index.

Strengthened `test_solute_is_actually_scaled_on_real_system`: besides the per-atom spot-check
it now asserts the fixture is 433 poly-Ala solute atoms and that **exactly** that many
particles' nonbonded params change (none outside the solute) — catches a mask that is too
broad or too narrow, which single-atom checks cannot. Demo left at 11 rungs (0.42 acceptance).

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
