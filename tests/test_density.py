"""Plateau-detection unit tests for openmm_pipelines.density.

The science under test: a trailing-window slope test must PASS a noisy plateau and
FAIL a slow sustained drift that a consecutive-segment difference would miss.
"""

import numpy as np
import pytest

from openmm_pipelines.core.density import assess_plateau, density_converged

TOL = 0.005  # 0.5 % fractional drift across the window
MIN_SEG = 8


def test_flat_plateau_converges():
    volumes = [100.0] * MIN_SEG
    assert density_converged(volumes, TOL, MIN_SEG)


def test_noisy_plateau_converges():
    # Physical NPT scatter (~0.2 % std) about a constant mean, zero trend.
    rng = np.random.default_rng(0)
    volumes = (100.0 + rng.normal(0, 0.2, MIN_SEG)).tolist()
    assert density_converged(volumes, TOL, MIN_SEG)


def test_slow_sustained_drift_fails():
    # 0.4 %/segment contraction: clears a 0.5 % PAIRWISE threshold every step, but
    # the slope accounts for ~2.8 % drift across the 8-segment window -> not converged.
    volumes = [100.0 * (1 - 0.004) ** i for i in range(MIN_SEG)]
    res = assess_plateau(volumes, TOL, MIN_SEG)
    assert not res.converged
    assert res.drift_rel > TOL
    # the pairwise difference the old test used would have passed:
    assert abs(volumes[-1] - volumes[-2]) / volumes[-2] < TOL


def test_only_trailing_window_used():
    # A big early transient then a flat tail: convergence judged on the tail only.
    volumes = [130.0, 120.0, 110.0] + [100.0] * MIN_SEG
    assert density_converged(volumes, TOL, MIN_SEG)


def test_drift_direction_agnostic():
    # Expansion drifts just like contraction.
    expand = [100.0 * (1 + 0.004) ** i for i in range(MIN_SEG)]
    assert not density_converged(expand, TOL, MIN_SEG)


def test_window_boundary_drift_rel():
    # Exact linear ramp: slope*(window-1)/mean must equal the end-to-end fractional span.
    volumes = [100.0 + i for i in range(MIN_SEG)]  # 100..107
    res = assess_plateau(volumes, TOL, MIN_SEG)
    expected = (volumes[-1] - volumes[0]) / res.mean_volume
    assert res.drift_rel == pytest.approx(expected)


def test_too_few_volumes_raises():
    with pytest.raises(AssertionError):
        assess_plateau([100.0, 100.0], TOL, MIN_SEG)


def test_nonpositive_tol_raises():
    with pytest.raises(AssertionError):
        assess_plateau([100.0] * MIN_SEG, 0.0, MIN_SEG)


def test_min_seg_below_two_raises():
    with pytest.raises(AssertionError):
        assess_plateau([100.0], TOL, 1)


def test_nonpositive_volume_raises():
    with pytest.raises(AssertionError):
        assess_plateau([100.0] * (MIN_SEG - 1) + [-1.0], TOL, MIN_SEG)
