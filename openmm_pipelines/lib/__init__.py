"""lib — the package toolkit: reusable building blocks the pipeline composes.

These are not pipeline stages; they are helpers the stages and composer call
(`restraints` is an OpenMM force helper; `density_convergence` is pure math). Kept at
the package root so any engine — preparation, production, future REMD/REST2 — can use
them without reaching into another engine's folder.
"""
