"""REST2Config — Replica Exchange with Solute Tempering (REST2), NVT.

Every replica runs at the physical `temperature_k` (T_0); replica m's solute is scaled by
λ_m = T_0 / T_m, spanning an *effective*-temperature range [T_0, max_effective_temperature_k].
Because only the solute is tempered, far fewer replicas span a given range than T-REMD needs
in explicit solvent. Pure Python (no OpenMM / openmmtools import). Production / exchange fields
and the derived iteration counts come from `MultiStateProductionConfig`.
"""

from __future__ import annotations

from pydantic import model_validator

from ..config import MultiStateProductionConfig


class REST2Config(MultiStateProductionConfig):
    """REST2 production over a geometric effective-temperature (λ) ladder (NVT)."""

    # --- λ ladder (effective-temperature range; physical T stays temperature_k) ---
    # Either give a replica count (geometric effective-T ladder built from T_0/max/n) or
    # pin the full λ ladder explicitly. λ[0] is always 1.0 (the reference, unscaled).
    max_effective_temperature_k: float = 450.0
    n_replicas: int = 6
    lambdas: tuple[float, ...] | None = None  # explicit λ ladder overrides min/max/n

    # --- solute selection (which atoms are "hot") ---
    solute_chain_index: int | None = None  # None = whole protein; else temper one chain

    @model_validator(mode="after")
    def _check_rest2(self) -> "REST2Config":
        if self.lambdas is not None:
            assert len(self.lambdas) >= 2, "need >= 2 lambdas"
            assert self.lambdas[0] == 1.0, (
                f"lambdas[0] must be 1.0 (the reference, unscaled state); got {self.lambdas[0]}"
            )
            assert all(0.0 < x <= 1.0 for x in self.lambdas), (
                f"every lambda must be in (0, 1]: {self.lambdas}"
            )
            assert list(self.lambdas) == sorted(self.lambdas, reverse=True) and len(
                set(self.lambdas)
            ) == len(self.lambdas), (
                f"lambdas must be strictly descending from 1.0: {self.lambdas}"
            )
        else:
            assert self.n_replicas >= 2, f"need >= 2 replicas, got {self.n_replicas}"
            assert self.max_effective_temperature_k > self.temperature_k, (
                f"max_effective_temperature_k ({self.max_effective_temperature_k}) must "
                f"exceed temperature_k ({self.temperature_k})"
            )
        return self
