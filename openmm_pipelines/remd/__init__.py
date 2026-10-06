"""remd — temperature replica-exchange MD (T-REMD) backend (produce side of the seam).

`run(eq, cfg)` consumes an `EquilibratedSystem` from `prepare()` and runs replica
exchange over a geometric temperature ladder via `openmmtools.multistate`. One-shot:
``remd.run(prepare(cfg), cfg)`` with a `REMDConfig`.

Production is **NVT** (the equilibrated box is already at the right density; replicas
hold volume fixed), so the barostat is stripped from the production system and the
replicas differ only in temperature.

Requires the optional ``[remd]`` extra (``openmmtools``).
"""

from .config import REMDConfig
from .ladder import geometric_ladder
from .production import run

__all__ = ["REMDConfig", "geometric_ladder", "run"]
