"""Requeue auto-resume demo: ONE idempotent entrypoint for a preemptable / requeued SLURM job.

`run_or_resume(cfg)` starts fresh on the first launch and continues from the last checkpoint on
every requeue. Submit `requeue_demo.sbatch` once with `#SBATCH --requeue`; on preemption or node
failure SLURM reruns the SAME job and this script resumes it — no self-resubmit code, no SIGTERM
trap. A hard kill loses at most one checkpoint interval.

Checkpoint cadence is the one knob that matters: `checkpoint_every_n_frames * traj_ps` (the sim
time between checkpoints) must be MUCH shorter than how long a node typically survives before
preemption, or a segment can be killed before it ever writes a checkpoint and every rerun restarts
from zero. Pure `--requeue` does not loop forever on that (a crash is not requeued) — it just fails
loud with "nothing to resume", which is your signal to checkpoint more often.
"""

from pathlib import Path

from openmm_pipelines import MDConfig, ProteinFF, WaterModel, run_or_resume

HERE = Path(__file__).parent
STRUCT = HERE.parent / "input_structures" / "helix_fusion.pdb"  # shared across demos


def main() -> None:
    cfg = MDConfig(
        pdb_in=STRUCT,
        outdir=HERE / "requeue_out",
        protein_ff=ProteinFF.AMBER14SB,
        water=WaterModel.TIP3P,
        padding_nm=0.8,
        nvt_equil_ns=0.01,  # 10 ps
        density_seg_steps=500,
        density_min_seg=3,
        density_max_seg=12,
        density_tol_rel=0.02,
        total_ns=0.2,  # 200 ps target
        traj_ps=2.0,
        checkpoint_every_n_frames=2,  # checkpoint every 4 ps (tiny demo cadence)
        platform="CUDA",
        seed=12345,
    )
    run_or_resume(cfg)  # fresh on first launch; continues from the last checkpoint on every requeue


if __name__ == "__main__":
    main()
