"""REMDConfig — T-REMD production over a geometric temperature ladder (NVT).

Pure Python (no OpenMM / openmmtools import) so construction + validation run anywhere.
The base `temperature_k` is the ladder's **minimum** (= the reference / equilibration
temperature); the ladder runs from there up to `max_temperature_k`. Production / exchange
fields and the derived iteration counts come from `MultiStateProductionConfig`.
"""

from __future__ import annotations

from pydantic import model_validator

from ..config import MultiStateProductionConfig


class REMDConfig(MultiStateProductionConfig):
    """T-REMD production over a geometric temperature ladder (NVT)."""

    # --- temperature ladder ---
    # Minimum temperature is the base `temperature_k`. Either give a replica count
    # (geometric ladder is built from min/max/n) or pin the full ladder explicitly.
    max_temperature_k: float = 400.0
    n_replicas: int = 12
    temperatures_k: tuple[float, ...] | None = None  # explicit ladder overrides min/max/n

    # Per-replica equilibration default (at each rung's own temperature, discarded).
    equilibration_ns: float = 0.5

    @model_validator(mode="after")
    def _check_remd(self) -> "REMDConfig":
        if self.temperatures_k is not None:
            assert len(self.temperatures_k) >= 2, "need >= 2 temperatures"
            assert list(self.temperatures_k) == sorted(self.temperatures_k), (
                f"temperatures_k must be ascending: {self.temperatures_k}"
            )
            assert self.temperatures_k[0] == self.temperature_k, (
                f"temperatures_k[0] ({self.temperatures_k[0]}) must equal temperature_k "
                f"({self.temperature_k}) — the base is the ladder minimum"
            )
        else:
            assert self.n_replicas >= 2, f"need >= 2 replicas, got {self.n_replicas}"
            assert self.max_temperature_k > self.temperature_k, (
                f"max_temperature_k ({self.max_temperature_k}) must exceed temperature_k "
                f"({self.temperature_k})"
            )
        return self
