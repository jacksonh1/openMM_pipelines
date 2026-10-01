"""Energy minimization — the static first stage on the built, restrained system.

Not equilibration: no dynamics, no time evolution. Relaxes steric clashes in the
solvated / added-atom geometry before any dynamics runs.
"""

from __future__ import annotations

import math


def minimize(sim, cfg) -> None:
    sim.minimizeEnergy()
    # Fail loud: a non-finite potential energy after minimization means a bad
    # topology/params (bad template, clashing added atoms) — do not let it reach
    # dynamics.
    pe = sim.context.getState(getEnergy=True).getPotentialEnergy()
    assert math.isfinite(pe.value_in_unit(pe.unit)), (
        f"non-finite energy after minimization: {pe}"
    )
