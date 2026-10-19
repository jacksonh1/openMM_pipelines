"""Compute-platform helpers shared across the prepare/produce seam.

OpenMM (and openmmtools) want mixed precision on a GPU platform and the default
(double) on CPU. This is the single definition of "which platforms are GPUs" and
the precision decision, rather than a copy in every module that builds a Context.
Pure Python — no OpenMM import, so it loads on any env.
"""

from __future__ import annotations

_GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")


def is_gpu(platform_name: str) -> bool:
    return platform_name in _GPU_PLATFORMS


def gpu_props(platform_name: str) -> dict[str, str]:
    """OpenMM platform properties: mixed precision on GPU, empty (double) on CPU."""
    return {"Precision": "mixed"} if is_gpu(platform_name) else {}


def precision_label(platform_name: str) -> str:
    """'mixed' on GPU, 'double' on CPU — for provenance records."""
    return "mixed" if is_gpu(platform_name) else "double"
