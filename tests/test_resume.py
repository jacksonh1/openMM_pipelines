"""resume / extend — the continuation contract.

Fast tests (no dynamics): config round-trip through the manifest, RunLayout.existing, and the
manifest collision/absence guards. Slow tests (real dynamics, CPU by default): MD extend is a
*continuous* trajectory (frame count, log<->trajectory alignment, monotonic steps, no duplicates),
resume-to-planned via the manifest, the fail-loud guards, and REMD extend via openmmtools
from_storage. Set OPENMM_TEST_PLATFORM=CUDA on a GPU node.
"""

import json
import os
from pathlib import Path

import pytest

from openmm_pipelines import md, prepare, resume
from openmm_pipelines.lib.run_layout import RunLayout, read_manifest, write_manifest
from openmm_pipelines.md.config import MDConfig

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"
PLATFORM = os.environ.get("OPENMM_TEST_PLATFORM", "CPU")


def _md_cfg(outdir, **extra) -> MDConfig:
    return MDConfig(
        pdb_in=STRUCT, outdir=outdir, platform=PLATFORM, padding_nm=0.6,
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        traj_ps=0.1, **extra,
    )


def _xtc_frame_count(path) -> int:
    from mdtraj.formats import XTCTrajectoryFile

    with XTCTrajectoryFile(str(path), "r") as fh:
        xyz, _times, _steps, _box = fh.read()
    return len(xyz)


def _log_steps(path) -> list[int]:
    lines = Path(path).read_text().splitlines()
    out = []
    for line in lines[1:]:  # skip the single header row
        if line.strip():
            out.append(int(float(line.split(",")[0].strip().strip('"'))))
    return out


# --------------------------------------------------------------------------- fast (no dynamics)

def test_config_roundtrips_through_manifest():
    """A reconstructed config equals the original (no silent defaults) for every engine."""
    from openmm_pipelines.remd.config import REMDConfig
    from openmm_pipelines.rest2.config import REST2Config
    from openmm_pipelines.resume import _reconstruct_config

    cases = {
        "md": MDConfig(pdb_in=STRUCT, outdir=Path("/tmp/x"), total_ns=7.5, checkpoint_every_n_frames=4),
        "remd": REMDConfig(pdb_in=STRUCT, outdir=Path("/tmp/x"), total_ns=3.0, n_replicas=4),
        "rest2": REST2Config(pdb_in=STRUCT, outdir=Path("/tmp/x"), total_ns=3.0, n_replicas=4),
    }
    for engine, cfg in cases.items():
        again = _reconstruct_config(engine, cfg.model_dump(mode="json"))
        assert again == cfg, engine


def test_existing_resolves_symlink_and_fails_loud(tmp_path):
    outdir, storage = tmp_path / "run", tmp_path / "store"
    layout = RunLayout(outdir, storage)
    # offloaded file: symlink in outdir, real bytes in storage -> existing() returns the real path
    real = layout.offload("remd/remd.nc")
    real.write_text("x")
    got = layout.existing("remd/remd.nc")
    assert got == real and not got.is_symlink()
    # recorded file: real in outdir -> existing() returns it
    rec = layout.record("remd/ladder.json")
    rec.write_text("{}")
    assert layout.existing("remd/ladder.json") == rec
    # missing -> fail loud
    with pytest.raises(FileNotFoundError):
        layout.existing("remd/nope.nc")


def test_manifest_collision_and_absence(tmp_path):
    write_manifest(tmp_path, {"engine": "md", "production_config": {}})
    assert read_manifest(tmp_path)["engine"] == "md"
    write_manifest(tmp_path, {"engine": "md", "production_config": {"total_ns": 9}})  # same engine OK
    with pytest.raises(RuntimeError):  # different engine into the same outdir
        write_manifest(tmp_path, {"engine": "remd", "production_config": {}})
    with pytest.raises(FileNotFoundError):  # no manifest -> cannot resume
        read_manifest(tmp_path / "empty")


# --------------------------------------------------------------------------- slow (real dynamics)

@pytest.mark.slow
def test_md_extend_is_continuous(tmp_path):
    out = tmp_path / "run"
    cfg = _md_cfg(out, total_ns=0.002)
    md.run(prepare(cfg), cfg)

    resume(out, extend_to_ns=0.004)

    prod = out / "production"
    extended_steps = round(0.004 * 1000 / cfg.dt_ps)
    expected_frames = extended_steps // cfg.traj_steps

    frames = _xtc_frame_count(prod / "trajectory.xtc")
    steps = _log_steps(prod / "production.log")
    assert frames == expected_frames, (frames, expected_frames)
    assert len(steps) == frames, "log and trajectory desynced"
    assert steps == sorted(steps) and len(set(steps)) == len(steps), "non-monotonic/duplicate steps"
    assert steps[-1] == extended_steps, "did not reach the extended target"

    assert list(out.glob("run_info.resume.*.json")), "missing per-segment provenance"
    rep = list(prod.glob("run_report.resume.*.json"))
    assert rep
    data = json.loads(rep[0].read_text())
    assert data["resumed"] is True
    assert data["resumed_from_step"] == cfg.prod_steps and data["target_step"] == extended_steps


@pytest.mark.slow
def test_md_resume_to_planned_via_manifest(tmp_path):
    out = tmp_path / "run"
    cfg = _md_cfg(out, total_ns=0.002)
    md.run(prepare(cfg), cfg)

    # Simulate a run that was *supposed* to reach 0.004 ns but stopped at 0.002: bump the target
    # in the authoritative manifest, then plain resume (no extend_to_ns).
    mpath = out / "resume_manifest.json"
    m = json.loads(mpath.read_text())
    m["production_config"]["total_ns"] = 0.004
    mpath.write_text(json.dumps(m))

    resume(out)

    extended_steps = round(0.004 * 1000 / cfg.dt_ps)
    assert _log_steps(out / "production" / "production.log")[-1] == extended_steps


@pytest.mark.slow
def test_md_resume_guards_fail_loud(tmp_path):
    out = tmp_path / "run"
    cfg = _md_cfg(out, total_ns=0.002)
    md.run(prepare(cfg), cfg)

    with pytest.raises(RuntimeError):  # already complete, no extend target
        resume(out)
    with pytest.raises(RuntimeError):  # extend target below what is already done
        resume(out, extend_to_ns=0.001)


@pytest.mark.slow
def test_remd_extend_continues_the_nc(tmp_path):
    from openmm_pipelines import remd
    from openmm_pipelines.remd.config import REMDConfig
    from openmmtools.multistate import MultiStateReporter

    out = tmp_path / "run"
    cfg = REMDConfig(
        pdb_in=STRUCT, outdir=out, platform=PLATFORM, padding_nm=0.6,
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        n_replicas=3, max_temperature_k=330.0,
        total_ns=0.001, exchange_attempt_ps=0.1, equilibration_ns=0.0, checkpoint_ps=0.1,
    )
    remd.run(prepare(cfg), cfg)
    assert cfg.n_iterations == 10

    resume(out, extend_to_ns=0.002)  # -> 20 iterations

    reporter = MultiStateReporter(str(out / "remd" / "remd.nc"), open_mode="r")
    try:
        last = reporter.read_last_iteration(last_checkpoint=False)
    finally:
        reporter.close()
    assert last >= 20, last
    assert list((out / "remd").glob("run_report.resume.*.json"))


@pytest.mark.slow
def test_md_resume_with_storage_directory(tmp_path):
    """Resume through the storage symlink: existing() must resolve it for truncate + append."""
    out, store = tmp_path / "run", tmp_path / "store"
    cfg = _md_cfg(out, total_ns=0.002, storage_directory=store)
    md.run(prepare(cfg), cfg)

    resume(out, extend_to_ns=0.004)

    traj = out / "production" / "trajectory.xtc"
    assert traj.is_symlink(), "trajectory should be offloaded to storage"
    extended_steps = round(0.004 * 1000 / cfg.dt_ps)
    frames = _xtc_frame_count(traj)  # follows the symlink into storage
    steps = _log_steps(out / "production" / "production.log")
    assert frames == extended_steps // cfg.traj_steps
    assert len(steps) == frames and steps[-1] == extended_steps


@pytest.mark.slow
def test_remd_resume_to_planned_run_branch(tmp_path):
    """The common walltime case: an interrupted .nc (5 of 10 iterations) resumed to the PLANNED
    length via openmmtools run() (not extend). A partial .nc is built with sampler.run(n_iterations=5),
    mirroring remd.run's construction so resume's existing()+from_storage path is exercised exactly."""
    from openmm import unit
    from openmmtools import mcmc, states
    from openmmtools.multistate import MultiStateReporter, ReplicaExchangeSampler

    from openmm_pipelines.context_cache import configure_global_platform
    from openmm_pipelines.lib.run_layout import RunLayout, write_manifest
    from openmm_pipelines.lib.system_edits import strip_barostats
    from openmm_pipelines.remd.config import REMDConfig
    from openmm_pipelines.remd.production import _ladder, _manifest

    out = tmp_path / "run"
    cfg = REMDConfig(
        pdb_in=STRUCT, outdir=out, platform=PLATFORM, padding_nm=0.6,
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        n_replicas=3, max_temperature_k=330.0,
        total_ns=0.001, exchange_attempt_ps=0.1, equilibration_ns=0.0, checkpoint_ps=0.1,
    )
    eq = prepare(cfg)
    assert cfg.n_iterations == 10

    # Build a PARTIAL run: created with number_of_iterations=10 but stopped at 5 (a kill mid-run).
    layout = RunLayout(out, cfg.storage_directory)
    (out / "remd").mkdir(parents=True, exist_ok=True)
    write_manifest(out, _manifest(cfg))
    configure_global_platform(cfg.platform)
    temperatures = _ladder(cfg)
    system = strip_barostats(eq.system)
    thermo = [states.ThermodynamicState(system=system, temperature=t * unit.kelvin) for t in temperatures]
    ss = states.SamplerState(positions=eq.positions, box_vectors=eq.box_vectors)
    move = mcmc.LangevinDynamicsMove(
        timestep=cfg.dt_ps * unit.picosecond, collision_rate=cfg.friction_ps / unit.picosecond,
        n_steps=cfg.steps_per_iteration, reassign_velocities=True,
    )
    sampler = ReplicaExchangeSampler(
        mcmc_moves=move, number_of_iterations=cfg.n_iterations,
        replica_mixing_scheme=cfg.replica_mixing_scheme, online_analysis_interval=None,
    )
    reporter = MultiStateReporter(
        str(layout.offload("remd/remd.nc")), checkpoint_interval=cfg.checkpoint_interval_iterations,
        checkpoint_storage=str(layout.offload("remd/remd_checkpoint.nc")),
    )
    sampler.create(thermo, [ss] * len(temperatures), storage=reporter)
    sampler.run(n_iterations=5)  # stop at 5 of 10
    reporter.close()
    del sampler

    resume(out)  # no extend -> target == stored 10, current 5 -> run() branch finishes to 10

    reporter = MultiStateReporter(str(out / "remd" / "remd.nc"), open_mode="r")
    try:
        last = reporter.read_last_iteration(last_checkpoint=False)
    finally:
        reporter.close()
    assert last >= 10, last
