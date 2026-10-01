"""PrepConfig / MDConfig validators and derived step counts.

Pure — no OpenMM. A dummy file stands in for `pdb_in` (validator only checks it
exists, not that it parses).
"""

import warnings

import pytest

from openmm_pipelines.config import BoxShape, PrepConfig
from openmm_pipelines.forcefield import ForceFieldError, ProteinFF, WaterModel
from openmm_pipelines.md.config import MDConfig


@pytest.fixture
def pdb(tmp_path):
    p = tmp_path / "in.pdb"
    p.write_text("REMARK dummy\nEND\n")
    return p


def test_defaults_construct(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out")
    assert cfg.protein_ff is ProteinFF.AMBER14SB
    assert cfg.water is WaterModel.TIP3P
    assert cfg.dt_ps == 0.004


def test_missing_pdb_raises(tmp_path):
    with pytest.raises(Exception):  # pydantic wraps the AssertionError
        PrepConfig(pdb_in=tmp_path / "nope.pdb", outdir=tmp_path / "out")


def test_frozen(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out")
    with pytest.raises(Exception):
        cfg.temperature_k = 310.0


def test_illegal_ff_water_rejected(pdb, tmp_path):
    with pytest.raises(ForceFieldError):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.CHARMM36M, water=WaterModel.OPC)


def test_hmr_required_for_4fs(pdb, tmp_path):
    with pytest.raises(Exception):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   dt_ps=0.004, hydrogen_mass_amu=1.008)


def test_2fs_physical_h_ok(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                     dt_ps=0.002, hydrogen_mass_amu=1.008)
    assert cfg.dt_ps == 0.002


def test_density_seg_bounds_checked(pdb, tmp_path):
    with pytest.raises(Exception):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   density_min_seg=20, density_max_seg=8)


def test_recommended_water_mismatch_warns(pdb, tmp_path):
    with pytest.warns(UserWarning, match="recommended with"):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.AMBER19SB, water=WaterModel.TIP3P)


def test_recommended_pairing_no_warn(pdb, tmp_path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning fails the test
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.AMBER19SB, water=WaterModel.OPC)


def test_derived_equil_steps(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                     nvt_equil_ns=0.2, dt_ps=0.004)
    assert cfg.nvt_equil_steps == 50_000  # 0.2 ns / 4 fs


def test_mdconfig_is_prepconfig_and_derives_prod_steps(pdb, tmp_path):
    cfg = MDConfig(pdb_in=pdb, outdir=tmp_path / "out", total_ns=100.0, traj_ps=10.0)
    assert isinstance(cfg, PrepConfig)
    assert cfg.prod_steps == 25_000_000  # 100 ns / 4 fs
    assert cfg.traj_steps == 2_500  # 10 ps / 4 fs
    assert cfg.checkpoint_steps == 25_000  # traj_steps * 10


def test_mdconfig_model_copy_update(pdb, tmp_path):
    base = MDConfig(pdb_in=pdb, outdir=tmp_path / "out", total_ns=50.0)
    hot = base.model_copy(update={"temperature_k": 330.0})
    assert hot.temperature_k == 330.0
    assert hot.total_ns == 50.0


def test_box_shape_enum(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out", box_shape=BoxShape.OCTAHEDRON)
    assert cfg.box_shape.value == "octahedron"
