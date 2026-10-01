"""Analysis figures (matplotlib, Agg backend — writes PNGs, no display)."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def _save(fig, path) -> Path:
    path = Path(path)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_timeseries(time_ps, y, ylabel, path, title=None):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(time_ps, y, lw=1.0)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    return _save(fig, path)


def plot_rmsf(resseq, rmsf_ang, path):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(resseq, rmsf_ang, lw=1.0)
    ax.set_xlabel("residue")
    ax.set_ylabel("Cα RMSF (Å)")
    ax.set_title("per-residue flexibility")
    return _save(fig, path)


def plot_ss_fractions(time_ps, fractions, path):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for name in ("helix", "sheet", "coil"):
        ax.plot(time_ps, fractions[name], lw=1.0, label=name)
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("fraction of residues")
    ax.set_ylim(0, 1)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("secondary structure")
    return _save(fig, path)
