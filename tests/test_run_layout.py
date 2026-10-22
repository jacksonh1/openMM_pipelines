"""RunLayout — pure filesystem behaviour (no OpenMM needed)."""

import json

import pytest

from openmm_pipelines.lib.run_layout import RunLayout, read_storage_directory


def test_no_storage_is_plain_outdir(tmp_path):
    outdir = tmp_path / "run"
    lay = RunLayout(outdir, None)
    assert lay.storage_home is None

    off = lay.offload("production/trajectory.xtc")
    rec = lay.record("production/run_report.json")
    assert off == outdir / "production" / "trajectory.xtc"
    assert rec == outdir / "production" / "run_report.json"
    off.write_text("traj")
    assert not off.is_symlink()
    lay.mirror("production")  # no-op, must not raise
    lay.clean()  # no-op, must not raise


def test_offload_writes_to_storage_and_symlinks_outdir(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)

    target = lay.offload("remd/remd.nc")
    assert target == lay.storage_home / "remd" / "remd.nc"
    target.write_text("ncdata")  # writer writes directly onto storage

    link = outdir / "remd" / "remd.nc"
    assert link.is_symlink()
    assert link.read_text() == "ncdata"  # resolves through the symlink
    assert (lay.storage_home / ".source_outdir").exists()  # ownership marker written


def test_record_is_real_and_mirror_completes_storage(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)

    lay.offload("remd/remd.nc").write_text("nc")  # big file -> storage + symlink
    lay.record("remd/run_report.json").write_text("{}")  # small file -> real in outdir

    lay.mirror("remd")
    # small file copied into storage; big file already there (written via offload), symlink skipped
    assert (lay.storage_home / "remd" / "run_report.json").read_text() == "{}"
    assert (lay.storage_home / "remd" / "remd.nc").read_text() == "nc"
    assert not (lay.storage_home / "remd" / "remd.nc").is_symlink()


def test_deterministic_and_collision_proof_name(tmp_path):
    storage = tmp_path / "store"
    a = RunLayout(tmp_path / "projA" / "out", storage)
    b = RunLayout(tmp_path / "projB" / "out", storage)
    # deterministic: same outdir -> same home
    assert a.storage_home == RunLayout(tmp_path / "projA" / "out", storage).storage_home
    # collision-proof: same basename, different path -> different home
    assert a.storage_home != b.storage_home
    assert a.storage_home.name.startswith("out-") and b.storage_home.name.startswith("out-")


def test_foreign_owner_fails_loud(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)
    home = lay.storage_home
    home.mkdir(parents=True)
    (home / ".source_outdir").write_text(json.dumps({"abs_outdir": "/some/other/run"}))
    with pytest.raises(RuntimeError, match="owned by a different run"):
        lay.offload("x.nc")


def test_offload_refuses_existing_output(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)
    lay.offload("remd/remd.nc").write_text("first")
    with pytest.raises(RuntimeError, match="refusing to overwrite"):
        lay.offload("remd/remd.nc")  # target already exists (a rerun) -> fail loud


def test_clean_removes_folder_and_symlinks_after_marker_check(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)
    lay.offload("remd/remd.nc").write_text("nc")
    lay.record("remd/run_report.json").write_text("{}")
    lay.mirror("remd")
    home = lay.storage_home
    assert home.exists()

    lay.clean()
    assert not home.exists()
    assert not (outdir / "remd" / "remd.nc").exists()  # dangling symlink removed


def test_clean_refuses_foreign_folder(tmp_path):
    outdir = tmp_path / "run"
    storage = tmp_path / "store"
    lay = RunLayout(outdir, storage)
    home = lay.storage_home
    home.mkdir(parents=True)
    (home / ".source_outdir").write_text(json.dumps({"abs_outdir": "/some/other/run"}))
    with pytest.raises(RuntimeError, match="refusing to clean"):
        lay.clean()


def test_read_storage_directory(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    assert read_storage_directory(outdir) is None  # no config.json
    (outdir / "config.json").write_text(json.dumps({"storage_directory": None}))
    assert read_storage_directory(outdir) is None  # explicit null
    (outdir / "config.json").write_text(json.dumps({"storage_directory": "/data/store"}))
    assert read_storage_directory(outdir) == __import__("pathlib").Path("/data/store")
