"""MDConfig — plain production MD. Adds production fields to the shared PrepConfig."""

from __future__ import annotations

from ..config import PrepConfig


class MDConfig(PrepConfig):
    """Plain production MD."""

    total_ns: float = 100.0
    traj_ps: float = 10.0

    @property
    def prod_steps(self) -> int:
        return round(self.total_ns * 1000 / self.dt_ps)

    @property
    def traj_steps(self) -> int:
        return round(self.traj_ps / self.dt_ps)

    @property
    def checkpoint_steps(self) -> int:
        return max(self.traj_steps * 10, 5000)
