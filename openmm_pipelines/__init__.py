"""openmm_pipelines — explicit-solvent protein MD on OpenMM.

One prepare/produce seam serves every backend: ``prepare(cfg)`` builds + equilibrates,
then the engine's ``run(eq, cfg)`` produces. One-shot: ``md.run(prepare(cfg), cfg)``.

Just ``import openmm_pipelines`` and navigate the tree — every public object hangs off
the package or an engine subpackage, so an editor's completion is the index:
``openmm_pipelines.rest2.REST2Config``, ``openmm_pipelines.md.run``,
``openmm_pipelines.analysis.rest2.analyze_rest2``. The short form still works too
(``from openmm_pipelines import rest2`` -> ``rest2.REST2Config``).
"""

from . import analysis, md, remd, rest2
from .config import BoxShape, PrepConfig
from .forcefield import ForceFieldError, NonbondedSpec, ProteinFF, WaterModel, resolve
from .md.config import MDConfig
from .lifecycle import is_complete, run_or_resume
from .preparation import EquilibratedSystem, prepare
from .resume import resume

__all__ = [
    # enums / prep physics
    "BoxShape",
    "PrepConfig",
    "ProteinFF",
    "WaterModel",
    "NonbondedSpec",
    "ForceFieldError",
    "resolve",
    # the prepare/produce seam
    "prepare",
    "EquilibratedSystem",
    "resume",
    "run_or_resume",  # idempotent fresh-vs-continue-vs-done entrypoint for requeued SLURM jobs
    "is_complete",
    # engine subpackages (each exports its own Config + run/resume)
    "md",
    "remd",
    "rest2",
    "MDConfig",  # kept at top level for the common case; REMDConfig/REST2Config live on their subpackage
    # analysis subpackage
    "analysis",
]
