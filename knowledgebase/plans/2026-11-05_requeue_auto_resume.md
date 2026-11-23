# Requeue / preemption auto-resume (run-until-complete on short or preemptable nodes)

**Status:** 2026-11-19 — the `--requeue` (preemption) path is IMPLEMENTED (`lifecycle.py`:
`run_or_resume` + `is_complete`; example `examples/MD_demo/requeue_demo.{py,sbatch}`). The
short-walltime `afterany` self-resubmit tail and its crash-loop / forward-progress watchdog guards
(items 3b, 4) were deliberately NOT built: pure `--requeue` does not loop on a crash (SLURM does not
requeue a non-zero exit), so those guards only matter if/when the `afterany` path is added. See
HISTORY 2026-11-19. Original design outline below. Reviewed by the architect agent. 2026-11-05.
Verdict: **sound, needs changes.** Lower-risk and higher-value than the multi-GPU plan; do
this first. The resume substrate already exists — this is driver glue only.

## Existing substrate (already built)

`resume(outdir, extend_to_ns, platform)` dispatches by engine via `resume_manifest.json`;
multistate reopens the `.nc` via `from_storage`; MD uses ping-pong State XML + xtc/log
truncate; `offload`/`existing` are EXDEV-safe; periodic checkpoints bound data loss.
`DECISIONS.md` already notes a planned non-exception `is_complete(outdir)`. `resume.py`'s
docstring explicitly scopes OUT the re-submission driver — that's what this plan adds.

## Plan

1. **`is_complete(outdir) -> bool`** — read-only AND context-free. Do NOT reuse
   `continue_to`'s `from_storage` (heavy; rebuilds systems) or call
   `configure_global_platform` (that would populate/poison the global context cache
   `context_cache.py` exists to protect). Multistate: `MultiStateReporter(open_mode='r')`
   → last iteration vs the manifest's derived `cfg.n_iterations`. MD: `XmlSerializer`-
   deserialize the ping-pong State and read `.getTime()` vs `cfg.prod_steps` — no Context,
   no Simulation. Source of truth = the real artifacts, **not** a separate progress file
   (a parallel counter is a second source of truth that can desync). Fail loud if the
   manifest is missing. Placement: beside `read_manifest` in `lib/run_layout.py`.
2. **`run_or_resume(outdir/cfg)`** — idempotent entrypoint: fresh outdir → prepare+run;
   incomplete → resume; complete → no-op. The one call the sbatch invokes every time.
   Crosses into fresh-run territory (which `resume.py` disclaims) → put it in a new
   `lifecycle.py`, one layer above `resume.py`.
3. **sbatch — two triggers, two idioms:**
   - **Preemptable partition** → `#SBATCH --requeue` + `scontrol requeue $SLURM_JOB_ID`.
     SLURM re-runs the SAME job on preemption, so write NO self-resubmit code — the driver
     just always calls `run_or_resume`.
   - **Short walltime** → trailing real-`if` `sbatch --dependency=afterany:$SLURM_JOB_ID this.sbatch`.
   - **Reject** job arrays (those are for independent parallel tasks, not sequential
     continuation). No SIGTERM trap needed — openmmtools checkpoints leave the `.nc`
     consistent at the last completed interval, MD ping-pong is torn-write safe; a hard kill
     loses ≤ 1 checkpoint interval. Guard the `set -e` final-line gotcha (use a real `if`).
4. **Two guards (the real failure modes):**
   - **Forward-progress guard** — if `checkpoint_ps` / `checkpoint_interval_iterations`
     exceeds the segment walltime, no checkpoint ever lands → every segment restarts from 0
     and is killed again → infinite resubmit burning allocation. Assert checkpoint interval
     ≪ expected segment walltime; detect "resumed but iteration/step did not advance" → fail
     loud.
   - **Crash-loop guard** — a run that errors deterministically at startup is also
     "incomplete" + nonzero exit; `afterany` + is_complete would resubmit it forever.
     Resubmit only when the segment made forward progress AND segments < a cap. Distinguish
     a walltime-kill from a startup crash.
5. **Checkpoint tuning** — set `checkpoint_ps` / `checkpoint_steps` relative to expected
   segment walltime; document the I/O-vs-rework tradeoff (worse when `storage_directory`
   offloads over a network fs).

## Tradeoffs

No speedup — total compute is unchanged (slightly more, from re-equilibration/restart
overhead); it buys *access* to cheaper/preemptable/short-walltime nodes. More restart surface
(the ping-pong truncate + `currentStep` healing is already a known silent-desync gotcha). A
resumed GPU trajectory differs from an uninterrupted one (fine for ensembles, not exact
reproduction).

## Verdict

Factoring correct (`is_complete` + `run_or_resume`). Blockers to clear: (a) `is_complete`
read-only + context-free; (b) add the forward-progress and crash-loop guards; (c) `--requeue`
for preempt, `afterany` tail for walltime, no array; (d) derive completion from real artifacts,
not a side file; (e) `run_or_resume` in a new lifecycle layer, `is_complete` next to
`read_manifest`.
