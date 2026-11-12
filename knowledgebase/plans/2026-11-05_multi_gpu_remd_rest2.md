# Multi-GPU for REMD / REST2 (MPI, one rank per GPU)

**Status:** design outline, no code written. Reviewed by the architect agent + web research
of openmmtools' supported approach. 2026-11-05.

## Supported approach (confirmed)

openmmtools has **no single-process multi-GPU**; its parallelism *is* MPI via `mpiplus`
(since 0.18.0). You run **one MPI rank per GPU** and `mpiplus.distribute` scatters the
replicas across ranks. Ranks < replicas is normal and intended: each rank propagates its
several replicas **serially** on its one GPU (8 replicas, 2 GPUs → 2 ranks → 4 replicas/GPU).
GPU assignment is by **CUDA_VISIBLE_DEVICES masking per rank** (each rank sees one GPU as
index 0, library stays device-unaware) — **not** by setting an OpenMM DeviceIndex list.

Sources: openmmtools release history; issues #738 (launch recipes), #648 (replicas↔GPUs),
#449 (multi-GPU mixing bug). OpenMM-native multi-GPU (`CudaDeviceIndex='0,1,2'`) is the WRONG
tool — it splits ONE System's force eval across GPUs (one huge system; scales poorly past 2),
whereas REMD wants one replica per GPU (data-parallel). Optional: oversubscribe ranks per GPU
via CUDA MPS; the plain rank-per-GPU path already gives multiple-replicas-per-GPU.

## Plan

1. **Platform** — do NOT set DeviceIndex. Mask one GPU per rank via SLURM
   `--gpus-per-task=1 --gpu-bind=closest` (or an OpenMPI appfile with `-x CUDA_VISIBLE_DEVICES=N`
   per line, #738). `configure_global_platform` keeps `{'Precision':'mixed'}`, unchanged per
   rank. If a multi-rank launch exposes >1 or 0 GPUs to a rank → crash loud, never default.
2. **MPI-gate the non-openmmtools I/O** in `remd.run`/`rest2.run` AND `*.resume`: mkdir,
   `write_manifest`, `ladder.json`, **`layout.offload` (FS-mutating: O_EXCL marker + symlink —
   the most non-MPI-safe call, easy to miss)**, `_run_report` + `write_json`, `layout.mirror`
   → rank 0 only. Rank 0 resolves the `.nc`/checkpoint path, broadcasts it, barrier; ALL ranks
   build the reporter and call `create`/`equilibrate`/`run` (collective). `layout.existing` is
   read-only → safe on all ranks. (openmmtools already gates its own `.nc` writes to rank 0.)
3. **One extracted `_run_multistate` skeleton** shared by REMD + REST2 (they are ~90%
   identical), carrying the rank-0/collective split once, for run AND resume — not scattered
   `is_rank_zero()` guards. Add `lib/mpi.py` (`rank`/`barrier`/`bcast`; mpi4py imported lazily
   only under MPI — single-GPU path never imports it; `ntasks>1` but import fails → crash).
4. **`prepare()` stays single-GPU** (separate step) — short vs production; parallelizing buys
   nothing. Only run/resume go multi-rank.
5. **sbatch** — `--ntasks-per-node=K --gpus-per-task=1`, launch `mpirun python driver.py`
   (SLURM handles host discovery); K ≤ n_replicas. Analysis unchanged (post-hoc, single proc).
6. **Correctness gate (mandatory).** openmmtools multi-GPU mixing has a troubled history:
   0.19.0 "poor mixing of odd states, use 1 GPU"; 0.21.3–0.21.5 fixed several multi-GPU bugs;
   the thread was reopened with no maintainer "resolved". We are on 0.26.0 (past the fixes) but
   MUST validate: run a short job 1-GPU vs K-GPU and confirm the exchange-mixing matrix matches
   before any production campaign (our pipeline already computes mixing stats). Keep
   `group_size=1` for `mpiplus.distribute` as a known fallback knob if mixing looks wrong.

   **Build this as its own small deliverable (run it when we build multi-GPU — don't skip).**
   A reusable check (a `tests/`-style GPU validation like `gpu_tests.sbatch`, or a one-shot
   script), NOT a per-campaign gate: take one config, run the same short problem twice (same
   seed/ladder, a few hundred iterations) on 1 rank and on K ranks; diff the exchange-mixing
   statistics our analysis already computes (`analysis/multistate.py` mixing matrix + mean
   neighbor acceptance + all-states-visited) and assert they agree within tolerance → PASS/FAIL
   with the offending state pairs, so a bad launch or a resurfaced `group_size` bug crashes loud
   instead of silently producing a wrong ensemble. Comparison is STATISTICAL (GPU runs aren't
   bit-reproducible), so calibrate the tolerance + iteration count so 1-GPU-vs-1-GPU with
   different seeds already passes. Run once per cluster / openmmtools-version change, then trust
   multi-GPU production until the environment changes.

## Resume interaction

`from_storage` + run/extend works multi-rank. The **replica count** (ladder size) is fixed by
the `.nc` and cannot change on resume, but the **rank/GPU count is free** — checkpoint on 4
GPUs, resume on 2. The fresh-vs-resume decision must be made on rank 0 and broadcast so a race
can't split ranks down different paths.

## Env status (openmm_env, checked 2026-11-05)

- `mpi4py` 4.1.2 built against Open MPI 5.0.11; in-env `mpirun` is Open MPI 5.0.11 — **they
  match** (the failure mode #738 warns about is absent).
- `clusterutils` **NOT installed** (choderalab's hostfile/appfile generator; ~7yr stale — hand-
  roll the SLURM launch instead).
- `mpiplus` is a dirty untagged dev build; `openmmtools` 0.26.0; `openmm` 8.6.1.
- **Open risk:** the env ships its own conda Open MPI 5.0.11 — confirm it integrates with the
  cluster's SLURM/PMIx launcher (srun vs the bundled mpirun) before trusting multi-node.

## Verdict

Direction sound (= the supported approach). Minimal-correct version: SLURM GPU masking +
rank-0 FS setup with broadcast/barrier in one `_run_multistate` skeleton (run + resume) +
single-GPU prepare + the 1-vs-K mixing-validation gate. Heavier lift than the requeue plan;
builds on that plan's shared resume seam, so do requeue first.
