"""MDConfig — plain production MD. Adds production fields to the shared PrepConfig."""

from __future__ import annotations

from pydantic import model_validator

from ..config import PrepConfig


class MDConfig(PrepConfig):
    """Plain production MD."""

    total_ns: float = 100.0
    traj_ps: float = 10.0
    # Resume-checkpoint cadence, counted in trajectory FRAMES (not raw steps or ps): the resume
    # checkpoint must land on a frame boundary so a resumed run appends cleanly, so we count
    # frames and derive the step interval as a multiple of traj_steps by construction. Default 10
    # frames is a fraction of a percent of wall time; lower loses less work on a crash at a little
    # more overhead, higher the reverse.
    checkpoint_every_n_frames: int = 10

    @model_validator(mode="after")
    def _check_md(self) -> "MDConfig":
        assert self.checkpoint_every_n_frames >= 1, "checkpoint_every_n_frames must be >= 1"
        return self

    @property
    def prod_steps(self) -> int:
        return round(self.total_ns * 1000 / self.dt_ps)

    @property
    def traj_steps(self) -> int:
        return round(self.traj_ps / self.dt_ps)

    @property
    def checkpoint_steps(self) -> int:
        # An integer multiple of traj_steps BY CONSTRUCTION, so the resume checkpoint always lands
        # on a trajectory-frame boundary (required for clean append on resume).
        return self.checkpoint_every_n_frames * self.traj_steps
