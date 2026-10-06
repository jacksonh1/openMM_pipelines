"""Temperature ladder for T-REMD.

Geometric (exponential) spacing is the standard default: for a system whose heat
capacity is roughly constant across the range, equal temperature *ratios* give roughly
equal exchange-acceptance between neighbours — the right target, since a ladder with a
weak link (a pair that never swaps) breaks replica diffusion through temperature space.
It is a starting ladder, not an optimized one; tune `n_replicas` from the measured
per-neighbour acceptance (see `analysis/remd.py`) — raise it where acceptance is low.

Pure numpy; no OpenMM / openmmtools import.
"""

from __future__ import annotations

import numpy as np


def geometric_ladder(t_min: float, t_max: float, n: int) -> list[float]:
    """`n` temperatures from `t_min` to `t_max`, geometrically spaced (endpoints incl.)."""
    assert n >= 2, f"need >= 2 replicas, got {n}"
    assert t_max > t_min > 0, f"require t_max ({t_max}) > t_min ({t_min}) > 0"
    return np.geomspace(t_min, t_max, n).tolist()
