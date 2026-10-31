"""resume(outdir, ...) -> Path — continue or extend an interrupted production run.

Reads the authoritative ``resume_manifest.json`` at the outdir root (written by ``run`` at the
start of production), reconstructs the exact frozen production config of the right engine (no
silent defaults — the manifest carries every field), applies optional ``platform`` /
``extend_to_ns`` overrides, and dispatches to that engine's ``resume``. The engine backend does the
continuation from its persisted state: MD reloads the ping-pong State + appends the trajectory;
REMD/REST2 reopen the ``.nc`` and continue the sampler via openmmtools ``from_storage`` + ``extend``.

Library-only: this continues a run from disk; re-submitting the SLURM job (so it calls ``resume``
instead of a fresh ``run``) is the driver's job. Completion is a derived query — a plain resume of
an already-complete run fails loud in the engine backend (pass ``extend_to_ns`` to extend).
"""

from __future__ import annotations

from pathlib import Path

from .lib.run_layout import read_manifest


def _reconstruct_config(engine: str, cfg_dict: dict):
    """Rebuild the exact frozen production config of `engine` from the manifest dict, validated like
    construction (a bad field or FF pairing fails loud here, as at first launch)."""
    if engine == "md":
        from .md.config import MDConfig

        return MDConfig.model_validate(cfg_dict)
    if engine == "remd":
        from .remd.config import REMDConfig

        return REMDConfig.model_validate(cfg_dict)
    if engine == "rest2":
        from .rest2.config import REST2Config

        return REST2Config.model_validate(cfg_dict)
    raise ValueError(f"unknown engine {engine!r} in resume_manifest.json")


def resume(outdir, *, extend_to_ns: float | None = None, platform: str | None = None) -> Path:
    """Continue the run in `outdir` to its planned length, or (with `extend_to_ns`) to a larger one.

    `platform` optionally overrides the compute platform for this segment (e.g. resume a CUDA run on
    CPU). Engine-specific backends (lazily imported so REMD/REST2's openmmtools dependency is only
    required when actually resuming one) carry out the continuation and fail loud if the run is
    already at/past the target.
    """
    outdir = Path(outdir)
    manifest = read_manifest(outdir)
    engine = manifest["engine"]
    cfg = _reconstruct_config(engine, manifest["production_config"])

    updates: dict = {}
    if platform is not None:
        updates["platform"] = platform
    if extend_to_ns is not None:
        updates["total_ns"] = extend_to_ns
    if updates:
        cfg = cfg.model_copy(update=updates)

    if engine == "md":
        from . import md

        return md.resume(cfg)
    if engine == "remd":
        from . import remd

        return remd.resume(cfg)
    if engine == "rest2":
        from . import rest2

        return rest2.resume(cfg)
    raise ValueError(f"unknown engine {engine!r} in {outdir / 'resume_manifest.json'}")
