"""Short end-to-end demo: prepare() + md.run() on a small structure, then print the
QC/metadata the output contract produced. Deliberately tiny (ps-scale) — a smoke test
of the full pipeline, not a scientific run."""

from pathlib import Path

from openmm_pipelines import MDConfig, ProteinFF, WaterModel, md, prepare

HERE = Path(__file__).parent
STRUCT = HERE / "input_structures" / "helix_fusion.pdb"


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
        total_ns=0.02,  # 20 ps production
        traj_ps=2.0,
        platform="CUDA",
        seed=12345,
    )

    eq = prepare(cfg)
    prod = md.run(eq, cfg)

    print(f"\n=== outputs in {cfg.outdir} ===")
    for rel in (
        "run_info.json",
        "prepared/prep_report.json",
        "production/run_report.json",
    ):
        print(f"\n----- {rel} -----")
        print((cfg.outdir / rel).read_text())
    print(f"\nproduction dir: {prod}")
    print("files:", sorted(p.name for p in prod.iterdir()))


if __name__ == "__main__":
    main()
