"""md — plain production MD backend (produce side of the seam).

`run(eq, cfg)` consumes an `EquilibratedSystem` from `prepare()` and runs one
trajectory. One-shot: ``md.run(prepare(cfg), cfg)``.
"""

from .config import MDConfig
from .production import run

__all__ = ["MDConfig", "run"]
