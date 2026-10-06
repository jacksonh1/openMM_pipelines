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


def plot_timeseries(x, y, ylabel, path, title=None, xlabel="time (ps)"):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(x, y, lw=1.0)
    ax.set_xlabel(xlabel)
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


def plot_neighbor_acceptance(neighbor_acceptance, path):
    """Exchange-acceptance per adjacent ladder rung (i <-> i+1). A near-zero bar is a
    broken ladder link where replicas stop diffusing through temperature space."""
    acc = np.asarray(neighbor_acceptance, dtype=float)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.bar(np.arange(len(acc)), acc)
    ax.axhline(0.2, color="k", lw=0.8, ls="--", label="0.2 rule-of-thumb floor")
    ax.set_xlabel("rung pair (i, i+1)")
    ax.set_ylabel("exchange acceptance")
    ax.set_ylim(0, 1)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("neighbour exchange acceptance")
    return _save(fig, path)


def plot_replica_state_walk(replica_states, path, max_replicas=4):
    """State index vs iteration for a few replicas — the visual proof that replicas are
    continuous coordinate sets walking through temperature space (not fixed-T ensembles)."""
    replica_states = np.asarray(replica_states)
    n_replicas = replica_states.shape[1]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for r in range(min(max_replicas, n_replicas)):
        ax.plot(replica_states[:, r], lw=0.7, alpha=0.8, label=f"replica {r}")
    ax.set_xlabel("iteration")
    ax.set_ylabel("ladder state index")
    ax.legend(loc="best", fontsize=8)
    ax.set_title("replicas walking through temperature space")
    return _save(fig, path)


def plot_ss_fractions(x, fractions, path, xlabel="time (ps)"):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for name in ("helix", "sheet", "coil"):
        ax.plot(x, fractions[name], lw=1.0, label=name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("fraction of residues")
    ax.set_ylim(0, 1)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("secondary structure")
    return _save(fig, path)
