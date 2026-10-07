"""rest2 — Replica Exchange with Solute Tempering (REST2) backend (produce side).

`run(eq, cfg)` consumes an `EquilibratedSystem` from `prepare()` and runs Hamiltonian
replica exchange where every replica is at the same physical temperature but replica m's
solute is scaled by λ_m (only solute DOF enter the exchange acceptance). One-shot:
``rest2.run(prepare(cfg), cfg)`` with a `REST2Config`. Production is **NVT** (barostat
stripped; equilibrated box already at density).

The one genuinely new piece is `scale_solute` (`rest2/scaling.py`) — the solute Hamiltonian
scaling (√λ charges, λ·ε LJ, λ bonded + CMAP). Requires the optional ``[remd]`` extra
(``openmmtools``).
"""

from .config import REST2Config
from .ladder import lambda_ladder
from .production import run
from .scaling import scale_solute

__all__ = ["REST2Config", "lambda_ladder", "scale_solute", "run"]
