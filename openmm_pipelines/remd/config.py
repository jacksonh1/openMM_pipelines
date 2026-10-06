"""REMDConfig — T-REMD production. Adds replica-exchange fields to the shared PrepConfig.

Pure Python (no OpenMM / openmmtools import) so construction + validation run anywhere.
The base `temperature_k` is the ladder's **minimum** (= the reference / equilibration
temperature); the ladder runs from there up to `max_temperature_k`.

openmmtools runs in **iterations**: one iteration = `steps_per_iteration` MD steps at
each replica's own temperature, then one exchange-attempt sweep. So production length
is expressed as a per-replica wall of MD time (`total_ns`) plus the MD time between
exchange attempts (`exchange_attempt_ps`), and the iteration counts derive from those.
"""

from __future__ import annotations

from pydantic import model_validator

from ..config import PrepConfig


class REMDConfig(PrepConfig):
    """T-REMD production over a geometric temperature ladder (NVT)."""

    # --- temperature ladder ---
    # Minimum temperature is the base `temperature_k`. Either give a replica count
    # (geometric ladder is built from min/max/n) or pin the full ladder explicitly.
    max_temperature_k: float = 400.0
    n_replicas: int = 12
    temperatures_k: tuple[float, ...] | None = None  # explicit ladder overrides min/max/n

    # --- production (per replica) ---
    total_ns: float = 20.0  # MD time collected per replica
    exchange_attempt_ps: float = 1.0  # MD time per iteration (between exchange attempts)
    equilibration_ns: float = 0.5  # per-replica equilibration at its own T, discarded

    # --- replica exchange ---
    replica_mixing_scheme: str = "swap-all"  # openmmtools: 'swap-all' | 'swap-neighbors'
    checkpoint_ps: float = 100.0  # full-state checkpoint interval (restart granularity)

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
        assert self.exchange_attempt_ps > 0, "exchange_attempt_ps must be > 0"
        return self

    @property
    def steps_per_iteration(self) -> int:
        return round(self.exchange_attempt_ps / self.dt_ps)

    @property
    def n_iterations(self) -> int:
        return round(self.total_ns * 1000 / self.exchange_attempt_ps)

    @property
    def equilibration_iterations(self) -> int:
        return round(self.equilibration_ns * 1000 / self.exchange_attempt_ps)

    @property
    def checkpoint_interval_iterations(self) -> int:
        return max(round(self.checkpoint_ps / self.exchange_attempt_ps), 1)
