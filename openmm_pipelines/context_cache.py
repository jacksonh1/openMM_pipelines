"""openmmtools global-context-cache platform setup, shared by the multistate backends.

openmmtools picks the compute platform off its *global* context cache, not per call, and
`set_platform` refuses a cache that is already populated (`RuntimeError: Cannot change
platform of a Context cache already in use`). So `run()` would not be re-entrant — a driver
that runs two backends (or two replica-exchange runs) in one process, and the test suite
running both slow multistate tests, would crash on the second. Emptying the cache first makes
the setup idempotent (and frees the prior run's cached contexts). Needs openmmtools (the
`[remd]` extra); no mdtraj/matplotlib, so the produce side imports it freely.
"""

from __future__ import annotations

from openmm import Platform
from openmmtools import cache

from .lib.platforms import gpu_props, is_gpu


def configure_global_platform(platform_name: str) -> None:
    """Point openmmtools' global context cache at `platform_name` (mixed precision on GPU).

    Empties the cache first so repeated `run()` calls in one process don't trip
    `set_platform`'s populated-cache guard.
    """
    platform = Platform.getPlatformByName(platform_name)
    # None (not {}) on CPU: set_platform's validator short-circuits on None and skips
    # building a throwaway probe Context.
    props = gpu_props(platform_name) if is_gpu(platform_name) else None
    cache.global_context_cache.empty()
    cache.global_context_cache.set_platform(platform, props)
