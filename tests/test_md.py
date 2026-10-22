"""End-to-end plain MD on a small real system (CPU by default) — produce + the default-overview
analysis + the storage_directory offload/mirror.

Marked slow (real dynamics). Tiny lengths — checks wiring, not sampling. Set
OPENMM_TEST_PLATFORM=CUDA on a GPU node.
"""

import json
import os
from pathlib import Path

import pytest

from openmm_pipelines import md, prepare
from openmm_pipelines.analysis import analyze
from openmm_pipelines.lib.run_layout import RunLayout, read_storage_directory
from openmm_pipelines.md.config import MDConfig

pytestmark = pytest.mark.slow

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"
PLATFORM = os.environ.get("OPENMM_TEST_PLATFORM", "CPU")


def _tiny_cfg(outdir, **extra) -> MDConfig:
    return MDConfig(
        pdb_in=STRUCT, outdir=outdir, platform=PLATFORM, padding_nm=0.6,
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        total_ns=0.002, traj_ps=0.1, **extra,
    )


@pytest.fixture(scope="module")
def md_run(tmp_path_factory):
    outdir = tmp_path_factory.mktemp("md_run")
    cfg = _tiny_cfg(outdir)
    md.run(prepare(cfg), cfg)
    return outdir, cfg


def test_run_writes_contract(md_run):
    outdir, _ = md_run
    prod = outdir / "production"
    for name in ("trajectory.xtc", "production.log", "final_state.xml", "run_report.json"):
        assert (prod / name).exists(), f"missing {name}"


def test_analyze_writes_default_overview(md_run):
    outdir, _ = md_run
    adir = analyze(outdir)
    for name in ("analysis_report.json", "analysis_data.npz", "rmsd.png", "rmsf.png",
                 "dssp_timeline.png", "clustering/cluster_populations.png",
                 "clustering/cluster_timeseries.png", "report.html", "processed.xtc",
                 "clustering/cluster_rep_00.pdb"):
        assert (adir / name).exists(), f"missing {name}"
    summary = json.loads((adir / "analysis_report.json").read_text())
    assert summary["clustering"]["n_clusters"] >= 1
    assert "analysis_config" in summary
    assert "interchain" not in summary  # helix_fusion is single-chain


def test_storage_directory_offloads_trajectory_and_completes_record(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    cfg = _tiny_cfg(outdir, storage_directory=storage)
    md.run(prepare(cfg), cfg)

    layout = RunLayout(outdir, read_storage_directory(outdir))
    home = layout.storage_home
    # the trajectory is offloaded: symlink in outdir, real bytes in storage
    assert (outdir / "production" / "trajectory.xtc").is_symlink()
    assert (home / "production" / "trajectory.xtc").exists()
    assert not (home / "production" / "trajectory.xtc").is_symlink()
    # small files mirrored into storage (storage = complete record)
    assert (home / "production" / "run_report.json").exists()

    # analysis works through the symlink and mirrors its overview into storage
    adir = analyze(outdir)
    assert (adir / "analysis_report.json").exists()
    assert (home / "analysis" / "analysis_report.json").exists()
    assert (home / "analysis" / "report.html").exists()  # interactive report in the record
    assert (home / "analysis" / "processed.xtc").exists()  # processed trajectory copied too
