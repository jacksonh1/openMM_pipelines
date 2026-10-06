"""End-to-end T-REMD on a small real system (CPU by default) — the prepare/produce seam
reused, replica exchange run, and the de-multiplex into a fixed-T ensemble.

Marked slow (real dynamics). Tiny ladder + tiny lengths — this checks wiring and the
demux, not sampling quality. Set OPENMM_TEST_PLATFORM=CUDA on a GPU node.
"""

import json
import os
from pathlib import Path

import pytest

from openmm_pipelines import prepare
from openmm_pipelines import remd
from openmm_pipelines.analysis.remd import analyze_remd
from openmm_pipelines.remd.config import REMDConfig

pytestmark = pytest.mark.slow

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"
PLATFORM = os.environ.get("OPENMM_TEST_PLATFORM", "CPU")


@pytest.fixture(scope="module")
def remd_run(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("remd_run")
    cfg = REMDConfig(
        pdb_in=STRUCT, outdir=outdir, platform=PLATFORM, padding_nm=0.6,
        # fast prepare
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        # tiny REMD: 3 rungs, 10 iterations of 25 steps, checkpoint every iteration
        n_replicas=3, max_temperature_k=330.0,
        total_ns=0.001, exchange_attempt_ps=0.1, equilibration_ns=0.0002,
        checkpoint_ps=0.1,
    )
    eq = prepare(cfg)
    remd.run(eq, cfg)
    return outdir, cfg


def test_run_writes_contract(remd_run):
    outdir, cfg = remd_run
    rdir = outdir / "remd"
    assert (rdir / "remd.nc").exists()
    ladder = json.loads((rdir / "ladder.json").read_text())["temperatures_k"]
    assert len(ladder) == cfg.n_replicas
    assert ladder[0] == pytest.approx(cfg.temperature_k)
    assert ladder[-1] == pytest.approx(cfg.max_temperature_k)

    report = json.loads((rdir / "run_report.json").read_text())
    assert report["n_replicas"] == 3
    assert report["n_iterations"] == cfg.n_iterations
    assert 0.0 <= report["mean_neighbor_acceptance"] <= 1.0


def test_analyze_demuxes_reference_state(remd_run):
    outdir, cfg = remd_run
    adir = analyze_remd(outdir, reference_state=0)

    exchange = json.loads((adir / "exchange_report.json").read_text())
    assert exchange["n_states"] == 3
    assert len(exchange["neighbor_acceptance"]) == 2
    assert "subdominant_eigenvalue" in exchange
    assert "state_statistical_inefficiency" in exchange

    summary = json.loads((adir / "remd_analysis_report.json").read_text())
    assert summary["reference_state"] == 0
    assert summary["n_frames_demuxed"] >= 1
    assert summary["rmsd_to_design_ang"]["mean"] >= 0.0

    assert (adir / "state00_ensemble.xtc").exists()
    assert (adir / "acceptance.png").exists()
    assert (adir / "replica_walk.png").exists()
