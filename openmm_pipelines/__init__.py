"""openmm_pipelines — explicit-solvent protein MD on OpenMM.

Plain MD is wired end to end: ``prepare(cfg)`` builds + equilibrates, ``md.run(eq,
cfg)`` produces. One-shot: ``md.run(prepare(cfg), cfg)``. T-REMD / REST2 are future
backends behind the same prepare/produce seam.
"""

from . import md
from .config import BoxShape, PrepConfig
from .forcefield import ForceFieldError, NonbondedSpec, ProteinFF, WaterModel, resolve
from .md.config import MDConfig
from .preparation import EquilibratedSystem, prepare

__all__ = [
    "BoxShape",
    "PrepConfig",
    "MDConfig",
    "ProteinFF",
    "WaterModel",
    "NonbondedSpec",
    "ForceFieldError",
    "resolve",
    "prepare",
    "EquilibratedSystem",
    "md",
]
