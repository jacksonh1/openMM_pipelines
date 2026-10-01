"""NPT density-equilibration convergence: a trailing-window slope (plateau) test.

Why a slope test, not a consecutive-segment difference
------------------------------------------------------
For a solvated box the segment-to-segment volume noise is small
(sigma_V / V ~ sqrt(kT * kappa / V) ~ 2e-3, and the standard error of a segment
mean is smaller still). A consecutive-difference test -- stop once
|V_n - V_{n-1}| / V_{n-1} <= tol -- therefore passes readily *even while the box
is still shrinking steadily*: a sustained 0.4 %/segment drift clears a 0.5 %
pairwise threshold on every comparison while the box contracts several percent
overall.

This fits a least-squares line through the trailing window and asks how much
systematic drift that slope accounts for *across the whole window*. Random NPT
fluctuation has ~zero slope and passes; a steady drift does not, however small
each individual step is. The scatter itself is the physical volume fluctuation of
an NPT box, not a sign of non-convergence, so it is deliberately not part of the
test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PlateauTest:
    """Result of the trailing-window drift test."""

    converged: bool
    window: int  # segments used (== min_seg)
    drift_rel: float  # |slope| * (window - 1) / mean : fractional drift across the window
    mean_volume: float  # nm^3, over the window


def assess_plateau(volumes: list[float], tol_rel: float, min_seg: int) -> PlateauTest:
    """Fit a line through the last ``min_seg`` segment volumes and judge the plateau.

    ASSUMES: ``volumes`` are per-segment mean volumes (nm^3) in chronological order
    and the segments are of equal duration, so the segment index is a valid time
    axis. Fails loudly on bad arguments rather than returning a misleading verdict.
    """
    assert tol_rel > 0, f"tol_rel must be > 0, got {tol_rel}"
    assert min_seg >= 2, f"a slope needs at least two points, got min_seg={min_seg}"
    assert len(volumes) >= min_seg, (
        f"need at least min_seg={min_seg} volumes, got {len(volumes)}"
    )

    tail = np.asarray(volumes[-min_seg:], dtype=float)
    assert np.all(tail > 0), f"volumes must be positive, got {tail.tolist()}"

    mean_v = float(tail.mean())
    x = np.arange(min_seg, dtype=float)
    slope = float(np.polyfit(x, tail, 1)[0])  # OLS line; degree-1 polynomial fit

    drift_rel = abs(slope) * (min_seg - 1) / mean_v
    return PlateauTest(
        converged=drift_rel <= tol_rel,
        window=min_seg,
        drift_rel=drift_rel,
        mean_volume=mean_v,
    )


def density_converged(volumes: list[float], tol_rel: float, min_seg: int) -> bool:
    """True when a least-squares line through the trailing ``min_seg`` volumes
    accounts for a fractional drift <= ``tol_rel``. Plateau test, not a
    consecutive-segment difference. Thin wrapper over :func:`assess_plateau` for
    callers that only need the verdict."""
    return assess_plateau(volumes, tol_rel, min_seg).converged
