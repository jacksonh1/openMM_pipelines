"""REMDConfig validation + derived iteration counts — pure (no OpenMM)."""

from pathlib import Path

import pytest

from openmm_pipelines.remd.config import REMDConfig

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"


def _cfg(**kw):
    return REMDConfig(pdb_in=STRUCT, outdir=Path("/tmp/_remd_cfg_unused"), **kw)


def test_derived_iteration_counts():
    cfg = _cfg(total_ns=20.0, exchange_attempt_ps=1.0, dt_ps=0.004,
               equilibration_ns=0.5, checkpoint_ps=100.0)
    assert cfg.steps_per_iteration == 250        # 1.0 ps / 0.004 ps
    assert cfg.n_iterations == 20_000            # 20 ns / 1 ps
    assert cfg.equilibration_iterations == 500   # 0.5 ns / 1 ps
    assert cfg.checkpoint_interval_iterations == 100


def test_min_is_base_temperature():
    cfg = _cfg(temperature_k=300.0, max_temperature_k=400.0, n_replicas=10)
    assert cfg.temperature_k == 300.0            # ladder minimum = reference T


def test_rejects_max_not_above_min():
    with pytest.raises(ValueError):
        _cfg(temperature_k=350.0, max_temperature_k=350.0)


def test_rejects_single_replica():
    with pytest.raises(ValueError):
        _cfg(n_replicas=1)


def test_explicit_ladder_must_start_at_base():
    with pytest.raises(ValueError):
        _cfg(temperature_k=300.0, temperatures_k=(310.0, 350.0, 400.0))


def test_explicit_ladder_must_ascend():
    with pytest.raises(ValueError):
        _cfg(temperature_k=300.0, temperatures_k=(300.0, 400.0, 350.0))


def test_explicit_ladder_accepted():
    cfg = _cfg(temperature_k=300.0, temperatures_k=(300.0, 330.0, 365.0, 400.0))
    assert cfg.temperatures_k == (300.0, 330.0, 365.0, 400.0)
