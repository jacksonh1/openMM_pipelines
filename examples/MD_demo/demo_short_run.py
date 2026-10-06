"""Short end-to-end demo: prepare() + md.run() on a small structure, then print the
QC/metadata the output contract produced. Deliberately tiny (ps-scale) — a smoke test
of the full pipeline, not a scientific run."""

from pathlib import Path

from openmm_pipelines import MDConfig, ProteinFF, WaterModel, analysis, md, prepare

HERE = Path(__file__).parent
STRUCT = HERE.parent / "input_structures" / "helix_fusion.pdb"  # shared across demos


def main() -> None:
    cfg = MDConfig(
        pdb_in=STRUCT,
        outdir=HERE / "demo_out",
        protein_ff=ProteinFF.AMBER14SB,
        water=WaterModel.TIP3P,
        padding_nm=0.8,
        nvt_equil_ns=0.01,  # 10 ps
        density_seg_steps=500,
        density_min_seg=3,
        density_max_seg=12,
        density_tol_rel=0.02,
        total_ns=0.2,  # 20 ps production
        traj_ps=2.0,
        platform="CUDA",
        seed=12345,
    )

    eq = prepare(cfg)
    md.run(eq, cfg)
    analysis.analyze(cfg.outdir)

    print(f"\n=== outputs in {cfg.outdir} ===")
    for rel in (
        "run_info.json",
        "prepared/prep_report.json",
        "production/run_report.json",
        "analysis/analysis_report.json",
    ):
        print(f"\n----- {rel} -----")
        print((cfg.outdir / rel).read_text())


if __name__ == "__main__":
    main()
