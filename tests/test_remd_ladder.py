"""Temperature-ladder math — pure, no OpenMM/openmmtools."""

import numpy as np
import pytest

from openmm_pipelines.remd.ladder import geometric_ladder


def test_geometric_ladder_endpoints_and_length():
    ladder = geometric_ladder(300.0, 400.0, 8)
    assert len(ladder) == 8
    assert ladder[0] == pytest.approx(300.0)
    assert ladder[-1] == pytest.approx(400.0)
    assert ladder == sorted(ladder)  # strictly ascending


def test_geometric_ladder_constant_ratio():
    ladder = geometric_ladder(300.0, 450.0, 6)
    ratios = [ladder[i + 1] / ladder[i] for i in range(len(ladder) - 1)]
    assert np.allclose(ratios, ratios[0])  # equal temperature ratios = geometric


def test_geometric_ladder_rejects_bad_args():
    with pytest.raises(AssertionError):
        geometric_ladder(300.0, 400.0, 1)  # need >= 2
    with pytest.raises(AssertionError):
        geometric_ladder(400.0, 300.0, 4)  # t_max must exceed t_min
