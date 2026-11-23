"""run_or_resume(cfg) — the one idempotent entrypoint a (re)submitted SLURM job calls.

Fresh outdir -> ``prepare`` + the engine's ``run``. Interrupted (manifest present, not yet at
target) -> ``resume`` to the planned length. Already complete -> no-op. This is the single call a
requeue sbatch invokes every time: with ``#SBATCH --requeue`` SLURM reruns the SAME job on
preemption / node-failure and this entrypoint continues it from the last checkpoint. The library
never talks to SLURM — eligibility for requeue is SLURM's job, and a plain non-zero crash is NOT
requeued (so a run that cannot make progress fails once, loudly, instead of looping).

``is_complete`` is read-only and context-free: it reads the real run artifacts (MD ping-pong State
XML time; multistate ``.nc`` last iteration) against the manifest's target WITHOUT building a
System/Context or touching the global context cache (which `context_cache.py` exists to protect).
A plain resume of an already-complete run fails loud in the engine backend, so the no-op
short-circuit here is what makes re-running a finished job safe — e.g. a preemption that lands
after the final step but before the process exits.

This layer sits one above ``resume.py`` (which deliberately disclaims the fresh-run and
re-submission driver); it is the fresh-vs-continue-vs-done decision, nothing more.
"""

from __future__ import annotations

from pathlib import Path

from .lib.run_layout import RunLayout, read_manifest
from .resume import _reconstruct_config, resume

_MANIFEST = "resume_manifest.json"


def run_or_resume(cfg) -> Path:
    """Start `cfg`'s run fresh, continue it, or no-op if it is already at its target.

    `cfg` is a frozen production config (``MDConfig`` / ``REMDConfig`` / ``REST2Config``). Safe to
    call on every (re)submission of the same job: the decision is driven entirely by what is on
    disk under ``cfg.outdir``.
    """
    outdir = Path(cfg.outdir)
    if (outdir / _MANIFEST).exists():
        if is_complete(outdir):
            return outdir  # already at target: idempotent re-run has nothing to do
        return resume(outdir)
    return _fresh_run(cfg)


def _fresh_run(cfg) -> Path:
    """prepare + the engine's run, dispatched on the concrete config type (fails loud otherwise)."""
    from . import md, remd, rest2
    from .md.config import MDConfig
    from .preparation import prepare
    from .remd.config import REMDConfig
    from .rest2.config import REST2Config

    eq = prepare(cfg)
    if isinstance(cfg, MDConfig):
        return md.run(eq, cfg)
    if isinstance(cfg, REMDConfig):
        return remd.run(eq, cfg)
    if isinstance(cfg, REST2Config):
        return rest2.run(eq, cfg)
    raise TypeError(f"run_or_resume: unsupported config type {type(cfg).__name__}")


def is_complete(outdir) -> bool:
    """True iff the run in `outdir` has reached its manifest target. Read-only + context-free.

    Fails loud if the manifest is absent (a run without one predates resume support). Reads the
    target from the manifest's frozen production config and the progress from the run's real
    artifacts — never from a side progress file (a second source of truth that can desync).
    """
    outdir = Path(outdir)
    manifest = read_manifest(outdir)  # fails loud if absent
    engine = manifest["engine"]
    cfg = _reconstruct_config(engine, manifest["production_config"])
    if engine == "md":
        return _md_complete(cfg)
    if engine in ("remd", "rest2"):
        return _multistate_complete(cfg, engine)
    raise ValueError(f"unknown engine {engine!r} in {outdir / _MANIFEST}")


def _md_complete(cfg) -> bool:
    """MD: newest loadable ping-pong State XML's sim time vs `cfg.prod_steps`. No Context built.

    No checkpoint yet -> not complete (a fresh run killed before its first checkpoint). Mirrors
    `md.production._load_latest_state`'s torn-file fallback (newest may be a SIGKILL-torn write)
    but deserializes a bare State instead of loading into a Simulation.
    """
    from openmm import XmlSerializer, unit

    from .md.production import _CKPT_A, _CKPT_B

    layout = RunLayout(cfg.outdir, cfg.storage_directory)
    present = [p for rel in (_CKPT_A, _CKPT_B)
               if (p := _existing_or_none(layout, rel)) is not None]
    if not present:
        return False
    present.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    errors: dict[str, str] = {}
    for p in present:
        try:
            state = XmlSerializer.deserialize(p.read_text())
        except Exception as exc:  # torn/invalid checkpoint: try the older ping-pong file
            errors[str(p)] = repr(exc)
            continue
        elapsed_ps = state.getTime().value_in_unit(unit.picosecond)
        return round(elapsed_ps / cfg.dt_ps) >= cfg.prod_steps
    raise RuntimeError(f"no loadable ping-pong checkpoint (all torn/invalid): {errors}")


def _multistate_complete(cfg, engine: str) -> bool:
    """REMD/REST2: last stored `.nc` iteration vs `cfg.n_iterations`. Opens read-only, no sampler.

    No `.nc` yet -> not complete (a fresh run killed before any iteration landed).
    """
    from openmmtools.multistate import MultiStateReporter

    from .remd.production import _NC as _REMD_NC
    from .rest2.production import _NC as _REST2_NC

    nc_name = _REMD_NC if engine == "remd" else _REST2_NC
    layout = RunLayout(cfg.outdir, cfg.storage_directory)
    nc = _existing_or_none(layout, f"{engine}/{nc_name}")
    if nc is None:
        return False
    reporter = MultiStateReporter(str(nc), open_mode="r")
    try:
        last_iteration = reporter.read_last_iteration()
    finally:
        reporter.close()
    return last_iteration >= cfg.n_iterations


def _existing_or_none(layout: RunLayout, relpath: str) -> Path | None:
    """`layout.existing(relpath)` or None if the artifact has not been written yet.

    `existing` fails loud for a resume (a missing file there IS an error); here a missing artifact
    is the legitimate "fresh run, not started/checkpointed yet" case, so translate it to None.
    """
    try:
        return layout.existing(relpath)
    except FileNotFoundError:
        return None
