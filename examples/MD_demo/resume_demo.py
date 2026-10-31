"""Resume/extend demo: run a short MD, then CONTINUE it to a longer target with resume().

Shows the walltime-recovery / extend pattern end to end. Deliberately tiny (ps-scale) — a smoke
test of resume, not a scientific run.

In a real walltime-limited SLURM job you resubmit the SAME script after a kill and let it decide
fresh-vs-continue (the `run_or_resume` idiom below); here we call both in one process so a single
submission demonstrates the whole cycle.
"""

from pathlib import Path

from openmm_pipelines import MDConfig, ProteinFF, WaterModel, md, prepare, resume

HERE = Path(__file__).parent
STRUCT = HERE.parent / "input_structures" / "helix_fusion.pdb"  # shared across demos


def run_or_resume(cfg: MDConfig) -> None:
    """The resubmit pattern: continue a prior (partial) run if one exists, else start fresh.

    Resubmitting the same sbatch after a walltime kill re-enters this and continues from the last
    checkpoint. The library never talks to SLURM — resubmission is the driver's job.
    """
    if (Path(cfg.outdir) / "resume_manifest.json").exists():
        resume(cfg.outdir)              # a prior run is present -> continue it to its planned length
    else:
        md.run(prepare(cfg), cfg)       # first submission -> fresh run


def main() -> None:
    outdir = HERE / "resume_out"
    cfg = MDConfig(
        pdb_in=STRUCT,
        outdir=outdir,
        protein_ff=ProteinFF.AMBER14SB,
        water=WaterModel.TIP3P,
        padding_nm=0.8,
        nvt_equil_ns=0.01,  # 10 ps
        density_seg_steps=500,
        density_min_seg=3,
        density_max_seg=12,
        density_tol_rel=0.02,
        total_ns=0.1,  # 100 ps "first segment"
        traj_ps=2.0,
        checkpoint_every_n_frames=2,  # resume checkpoint every 2 frames (tiny demo cadence)
        platform="CUDA",
        seed=12345,
    )

    # First segment (a fresh run). In a real walltime job this is all that runs before the kill.
    run_or_resume(cfg)

    # Continue the finished segment to a longer target — the same call a resubmit would make, with
    # an explicit larger length. (extend_to_ns is an absolute target.)
    resume(outdir, extend_to_ns=0.2)  # -> 200 ps

    seg = sorted(outdir.glob("production/run_report.resume.*.json"))[-1]
    print(f"\n----- {seg.name} (resume-segment QC) -----")
    print(seg.read_text())


if __name__ == "__main__":
    main()
