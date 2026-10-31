"""Shared resume/extend logic for the openmmtools multistate backends (REMD, REST2).

Both backends reopen the existing ``.nc`` and continue the *same* sampler via openmmtools' native
``from_storage`` + ``run``/``extend`` — the documented resume idiom. openmmtools rebuilds the full
sampler (ladder, move, thermodynamic states) from the ``.nc``, so resume needs no `eq`/ladder, only
the target iteration count. The run-vs-extend decision and the "already complete" guard live here so
REMD and REST2 share one tested path; each backend keeps only its own ladder + QC report.
"""

from __future__ import annotations

import datetime as _dt
import time

from openmmtools.multistate import MultiStateReporter, ReplicaExchangeSampler

from .lib.run_layout import RunLayout


def utc_stamp() -> str:
    """A filesystem-safe UTC stamp for per-segment resume output filenames."""
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def open_existing_reporter(
    layout: RunLayout, subdir: str, nc_name: str, checkpoint_interval: int
) -> MultiStateReporter:
    """Reopen a run's existing ``.nc`` (+ its ``_checkpoint.nc``) in place for continuation.

    NOT via ``offload`` (which refuses an existing file); ``RunLayout.existing`` resolves the storage
    symlink and the post-SIGKILL case. Constructing a ``MultiStateReporter`` with the default
    ``open_mode=None`` does not open/truncate the file — ``from_storage`` opens it for continuation.
    The checkpoint name is derived as ``"<stem>_checkpoint.nc"`` to match what ``run`` offloaded.
    """
    stem = nc_name[:-3] if nc_name.endswith(".nc") else nc_name
    return MultiStateReporter(
        str(layout.existing(f"{subdir}/{nc_name}")),
        checkpoint_interval=checkpoint_interval,
        checkpoint_storage=str(layout.existing(f"{subdir}/{stem}_checkpoint.nc")),
    )


def continue_to(
    reporter: MultiStateReporter, target_iterations: int, total_ns: float
) -> tuple[ReplicaExchangeSampler, float]:
    """``from_storage`` then ``run``/``extend`` the sampler to ``target_iterations``.

    Returns ``(sampler, wall_seconds)``. The caller MUST keep the returned sampler referenced until
    it has finished reading `reporter` for the QC report: openmmtools closes the reporter when the
    sampler is garbage-collected, and a closed reporter raises ``NetCDF: Not a valid ID`` on read.

    Fails loud if the run is already at/past the target (pass a larger ``extend_to_ns`` to extend).
    When the target exceeds the stored target we ``extend`` past it; otherwise ``run`` finishes to the
    stored target (plain resume).
    """
    sampler = ReplicaExchangeSampler.from_storage(reporter)
    current = sampler.iteration
    if target_iterations <= current:
        raise RuntimeError(
            f"run at iteration {current} already >= target {target_iterations} "
            f"({total_ns} ns/replica); pass a larger extend_to_ns to extend."
        )
    t0 = time.perf_counter()
    if target_iterations > sampler.number_of_iterations:
        sampler.extend(n_iterations=target_iterations - current)  # extend past the stored target
    else:
        sampler.run()  # finish to the stored target
    return sampler, time.perf_counter() - t0
