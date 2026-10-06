"""Short end-to-end T-REMD demo: prepare() + remd.run() + analyze_remd() on a small
structure, then print the QC the output contract produced. Deliberately tiny (ps-scale,
4 rungs) — a smoke test of the replica-exchange pipeline and the de-multiplex, NOT a
scientific run. Needs the `[remd]` extra (openmmtools).

The point to notice: a replica trajectory walks THROUGH temperature space, so the
analysis de-multiplexes the reference rung (state 0 = T_min) into a fixed-temperature
ensemble before measuring drift from the design pose. See `replica_walk.png`.
"""

from pathlib import Path

from openmm_pipelines import ProteinFF, WaterModel, prepare, remd
from openmm_pipelines.analysis.remd import analyze_remd
from openmm_pipelines.remd.config import REMDConfig

HERE = Path(__file__).parent
STRUCT = HERE.parent / "input_structures" / "helix_fusion.pdb"  # shared across demos


def main() -> None:
    cfg = REMDConfig(
        pdb_in=STRUCT,
        outdir=HERE / "demo_remd_out",
        protein_ff=ProteinFF.AMBER14SB,
        water=WaterModel.TIP3P,
        padding_nm=0.8,
        # fast prepare (shared with plain MD)
        nvt_equil_ns=0.01,  # 10 ps
        density_seg_steps=500,
        density_min_seg=3,
        density_max_seg=12,
        density_tol_rel=0.02,
        # ladder + short REMD (NVT production)
        temperature_k=300.0,  # T_min = reference / design temperature
        max_temperature_k=360.0,
        # Explicit-solvent T-REMD needs neighbour gaps small enough for the replicas'
        # potential-energy histograms to OVERLAP: adjacent acceptance ~ exp(-(Δβ·σ_E)²/2)
        # and σ_E grows with the box's heat capacity (∝ degrees of freedom). For this
        # ~5.4k-atom box, 4 rungs over 300-360 K gives ΔT≈19 K → acceptance ~1e-7 (no
        # exchange at all). ~14 rungs here drops neighbour gaps to ~6 K for real mixing.
        # (This √DOF replica-count scaling is exactly why REST2 — which scales only the
        # SOLUTE — is the right tool for explicit solvent; it is the next backend.)
        n_replicas=14,
        total_ns=0.05,  # 50 ps collected per replica
        exchange_attempt_ps=1.0,  # swap attempt every 1 ps
        equilibration_ns=0.005,  # 5 ps per-replica equilibration (discarded)
        checkpoint_ps=2.0,  # demux frame resolution
        platform="CUDA",
        seed=12345,
    )

    eq = prepare(cfg)
    remd.run(eq, cfg)
    analyze_remd(cfg.outdir, reference_state=0)

    print(f"\n=== outputs in {cfg.outdir} ===")
    for rel in (
        "run_info.json",
        "prepared/prep_report.json",
        "remd/ladder.json",
        "remd/run_report.json",
        "remd/analysis/exchange_report.json",
        "remd/analysis/remd_analysis_report.json",
    ):
        print(f"\n----- {rel} -----")
        print((cfg.outdir / rel).read_text())


if __name__ == "__main__":
    main()
