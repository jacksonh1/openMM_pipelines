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
    assert cfg.protein_ff is ProteinFF.AMBER19SB  # modern default (QM-trained backbone-CMAP)
    assert cfg.water is WaterModel.OPC  # ff19SB's parameterization water
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


def test_cation_pi_nbfix_requires_charmm36m(pdb, tmp_path):
    # Default FF is ff19SB -> cation-π must be refused (parameters are CHARMM36m-specific).
    with pytest.raises(ForceFieldError):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out", cation_pi_nbfix=True)
    # Plain CHARMM36 is also refused: the paper parameterized C36m, not C36.
    with pytest.raises(ForceFieldError):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out", cation_pi_nbfix=True,
                   protein_ff=ProteinFF.CHARMM36, water=WaterModel.TIP3P)


def test_cation_pi_nbfix_ok_on_charmm36m(pdb, tmp_path):
    cfg = PrepConfig(pdb_in=pdb, outdir=tmp_path / "out", cation_pi_nbfix=True,
                     protein_ff=ProteinFF.CHARMM36M, water=WaterModel.TIP3P)
    assert cfg.cation_pi_nbfix is True


def test_cation_pi_nbfix_off_by_default(pdb, tmp_path):
    assert PrepConfig(pdb_in=pdb, outdir=tmp_path / "out").cation_pi_nbfix is False


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
    # ff14SB with a 4-site water is legal but not recommended -> soft warning (no hard lock).
    with pytest.warns(UserWarning, match="recommended with"):
        PrepConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.AMBER14SB, water=WaterModel.TIP4PEW)


def test_required_water_enforced(pdb, tmp_path):
    # ff19SB's OPC is mandatory, not merely recommended: a mismatch is rejected, not warned.
    with pytest.raises(ForceFieldError):
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


# --- config-file round-trip (to_json / from_json) ---

def test_json_round_trip_preserves_config(pdb, tmp_path):
    cfg = MDConfig(
        pdb_in=pdb, outdir=tmp_path / "out", storage_directory=tmp_path / "store",
        protein_ff=ProteinFF.CHARMM36M, water=WaterModel.TIP3P, cation_pi_nbfix=True,
        total_ns=250.0, traj_ps=10.0, seed=7,
    )
    path = cfg.to_json(tmp_path / "cfg" / "run.json")
    text = path.read_text()
    # Enums are written by NAME (human-readable + reloadable, unlike the tuple WaterModel value).
    assert '"water": "TIP3P"' in text
    assert '"protein_ff": "CHARMM36M"' in text
    back = MDConfig.from_json(path)
    assert back == cfg
    assert back.water is WaterModel.TIP3P and back.protein_ff is ProteinFF.CHARMM36M
    assert back.cation_pi_nbfix is True and back.seed == 7


def test_from_json_accepts_legacy_value_form(pdb, tmp_path):
    """A config.json written before name-serialization (water as a [name, n] list, protein_ff
    as its xml value) must still load."""
    import json
    cfg = MDConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.CHARMM36M, water=WaterModel.TIP3P)
    raw = cfg.model_dump(mode="json")
    raw["water"] = ["tip3p", 3]                 # legacy tuple-value form
    raw["protein_ff"] = "charmm36_2024.xml"     # legacy by-value form
    path = tmp_path / "legacy.json"
    path.write_text(json.dumps(raw))
    back = MDConfig.from_json(path)
    assert back.water is WaterModel.TIP3P and back.protein_ff is ProteinFF.CHARMM36M


def test_from_json_revalidates(pdb, tmp_path):
    """from_json runs the same validation as construction — an illegal combo fails loud."""
    import json
    cfg = MDConfig(pdb_in=pdb, outdir=tmp_path / "out",
                   protein_ff=ProteinFF.CHARMM36M, water=WaterModel.TIP3P, cation_pi_nbfix=True)
    raw = cfg.model_dump(mode="json")
    raw["protein_ff"] = "CHARMM36"  # non-36m + cation_pi_nbfix -> must raise on load
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(raw))
    with pytest.raises(ForceFieldError):
        MDConfig.from_json(path)
