"""Short end-to-end REST2 demo: prepare() + rest2.run() + analyze_rest2() on a small
structure, then print the QC the output contract produced. Deliberately tiny (ps-scale,
few rungs) — a smoke test of the solute-tempering pipeline and the de-multiplex, NOT a
scientific run. Needs the `[remd]` extra (openmmtools).

The showcase vs the T-REMD demo: that one needed ~14 rungs over 300-360 K for the
explicit-solvent box's potential-energy histograms to overlap (acceptance scales with
√DOF of the WHOLE box). REST2 tempers only the SOLUTE, so the same effective-temperature
reach — here 300 K up to an effective 450 K — mixes with a HANDFUL of replicas. Every
replica runs at the physical 300 K; replica m's solute is scaled by λ_m = 300/T_eff_m.
State 0 is λ=1, the unscaled physical ensemble — the only one that is a real Boltzmann
sample at the design temperature, and the one analysis de-multiplexes. See `replica_walk.png`.
"""

from pathlib import Path

from openmm_pipelines import ProteinFF, WaterModel, prepare, rest2
from openmm_pipelines.analysis.rest2 import analyze_rest2
from openmm_pipelines.rest2.config import REST2Config

HERE = Path(__file__).parent
STRUCT = HERE.parent / "input_structures" / "helix_fusion.pdb"  # shared across demos


def main() -> None:
    cfg = REST2Config(
        pdb_in=STRUCT,
        outdir=HERE / "demo_rest2_out",
        protein_ff=ProteinFF.AMBER14SB,
        water=WaterModel.TIP3P,
        padding_nm=0.8,
        # fast prepare (shared with plain MD / T-REMD)
        nvt_equil_ns=0.01,  # 10 ps
        density_seg_steps=500,
        density_min_seg=3,
        density_max_seg=12,
        density_tol_rel=0.02,
        # λ ladder + short REST2 (NVT production, all replicas at temperature_k)
        temperature_k=300.0,  # physical T_0 = reference / design temperature (λ=1)
        max_effective_temperature_k=450.0,  # top rung's EFFECTIVE temperature
        # Only the solute is tempered, so this 300->450 K effective reach needs just a few
        # rungs (contrast the 14 the T-REMD demo needed for a narrower 300-360 K range).
        n_replicas=6,
        total_ns=0.05,  # 50 ps collected per replica
        exchange_attempt_ps=1.0,  # swap attempt every 1 ps
        equilibration_ns=0.005,  # 5 ps per-replica equilibration (discarded)
        checkpoint_ps=2.0,  # demux frame resolution
        # solute_chain_index=None,  # whole protein tempered; set an index to heat one chain
        platform="CUDA",
        seed=12345,
    )

    eq = prepare(cfg)
    rest2.run(eq, cfg)
    analyze_rest2(cfg.outdir, reference_state=0)  # state 0 = λ=1 physical ensemble

    print(f"\n=== outputs in {cfg.outdir} ===")
    for rel in (
        "run_info.json",
        "prepared/prep_report.json",
        "rest2/ladder.json",
        "rest2/run_report.json",
        "rest2/analysis/exchange_report.json",
        "rest2/analysis/rest2_analysis_report.json",
    ):
        print(f"\n----- {rel} -----")
        print((cfg.outdir / rel).read_text())


if __name__ == "__main__":
    main()
