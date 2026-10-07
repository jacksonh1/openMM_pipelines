"""REST2Config validation + λ ladder — pure (no OpenMM)."""

from pathlib import Path

import pytest

from openmm_pipelines.rest2.config import REST2Config
from openmm_pipelines.rest2.ladder import lambda_ladder

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"


def _cfg(**kw):
    return REST2Config(pdb_in=STRUCT, outdir=Path("/tmp/_rest2_cfg_unused"), **kw)


def test_derived_iteration_counts_shared_with_remd():
    cfg = _cfg(total_ns=20.0, exchange_attempt_ps=1.0, dt_ps=0.004,
               equilibration_ns=0.5, checkpoint_ps=100.0)
    assert cfg.steps_per_iteration == 250
    assert cfg.n_iterations == 20_000
    assert cfg.equilibration_iterations == 500
    assert cfg.checkpoint_interval_iterations == 100


def test_rejects_max_eff_not_above_physical():
    with pytest.raises(ValueError):
        _cfg(temperature_k=300.0, max_effective_temperature_k=300.0)


def test_rejects_single_replica():
    with pytest.raises(ValueError):
        _cfg(n_replicas=1)


def test_explicit_lambdas_must_start_at_one():
    with pytest.raises(ValueError):
        _cfg(lambdas=(0.9, 0.7, 0.5))


def test_explicit_lambdas_must_descend_strictly():
    with pytest.raises(ValueError):
        _cfg(lambdas=(1.0, 0.7, 0.7))
    with pytest.raises(ValueError):
        _cfg(lambdas=(1.0, 0.5, 0.7))


def test_explicit_lambdas_accepted():
    cfg = _cfg(lambdas=(1.0, 0.8, 0.6, 0.5))
    assert cfg.lambdas == (1.0, 0.8, 0.6, 0.5)


def test_lambda_ladder_reference_is_one_and_descends():
    lambdas, temps = lambda_ladder(300.0, 450.0, 6)
    assert len(lambdas) == 6
    assert lambdas[0] == pytest.approx(1.0)
    assert temps[0] == pytest.approx(300.0)
    assert temps[-1] == pytest.approx(450.0)
    # λ = T0 / T_eff, strictly descending
    assert all(a > b for a, b in zip(lambdas, lambdas[1:]))
    assert all(lam == pytest.approx(300.0 / t) for lam, t in zip(lambdas, temps, strict=True))
