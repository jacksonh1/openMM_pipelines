"""core — the shared equilibration pipeline (prepare side of the seam).

Everything through equilibration lives here; it is identical across MD, REMD and
REST2. Public entry points (prepare, EquilibratedSystem, PrepConfig) are re-exported
as their modules land. Import submodules directly meanwhile, e.g.
``from openmm_pipelines.core.density import density_converged``.
"""
