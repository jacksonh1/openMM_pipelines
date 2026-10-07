"""REST2 λ ladder.

REST2 spans an *effective-temperature* range [T_0, T_max_eff] while every replica runs at
the physical T_0; replica m's solute is scaled by λ_m = T_0 / T_m. We build the ladder in
effective temperature with the same geometric spacing as T-REMD (roughly equal acceptance
between neighbours when the tempered heat capacity is ~constant), then map to λ.

State 0 is the reference: T_0 -> λ = 1 (unscaled, the physical ensemble to analyze). λ
therefore descends from 1.0 as the effective temperature climbs. Returns the λ ladder and
the matching effective temperatures (kept for the run report / plots).

Pure numpy; no OpenMM / openmmtools import.
"""

from __future__ import annotations

from ..remd.ladder import geometric_ladder


def lambdas_from_effective_temperatures(t0: float, effective_temperatures_k) -> list[float]:
    """λ_m = t0 / T_m for each effective temperature (state 0 = t0 -> λ = 1)."""
    temps = list(effective_temperatures_k)
    assert temps and temps[0] == t0, (
        f"effective-temperature ladder must start at t0={t0}, got {temps[:1]}"
    )
    return [t0 / t for t in temps]


def lambda_ladder(t0: float, t_max_eff: float, n: int) -> tuple[list[float], list[float]]:
    """`(lambdas, effective_temperatures_k)` for a geometric effective-T ladder.

    `lambdas[0] == 1.0` (reference); λ descends as effective temperature climbs.
    """
    temps = geometric_ladder(t0, t_max_eff, n)  # geometric in effective temperature
    return lambdas_from_effective_temperatures(t0, temps), temps
