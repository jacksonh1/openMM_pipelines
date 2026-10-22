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


def plot_series(x, series, ylabel, path, title, xlabel="time (ps)"):
    """One line per labelled series on a shared x-axis (e.g. a metric per chain pair)."""
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for label, y in series.items():
        ax.plot(x, y, lw=1.0, label=str(label))
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if len(series) > 1:
        ax.legend(loc="best", fontsize=8)
    return _save(fig, path)


def plot_cluster_populations(populations, path, design_like_cluster=-1):
    """Cluster sizes (rank-ordered, descending) with percent labels; highlight the
    design-like cluster. One dominant bar near rank 0 = a stable, homogeneous ensemble."""
    pops = np.asarray(populations, dtype=float)
    total = pops.sum()
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ranks = np.arange(len(pops))
    colors = ["tab:green" if r == design_like_cluster else "tab:blue" for r in ranks]
    bars = ax.bar(ranks, pops, color=colors)
    for r, bar in enumerate(bars):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                f"{100 * pops[r] / total:.0f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("cluster (rank by population)")
    ax.set_ylabel("frames")
    ax.set_xticks(ranks)
    ax.set_title("conformational cluster populations"
                 + (" (green = design-like)" if design_like_cluster >= 0 else ""))
    return _save(fig, path)


def plot_cluster_timeseries(x, labels, path, xlabel="time (ps)"):
    """Cluster assignment vs the x-axis of the clustered frames. For a demuxed ensemble the
    x-axis is a frame index, not time — reads as cluster occupancy, not kinetics."""
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.scatter(x, labels, s=8, alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("cluster (rank)")
    ax.set_title("cluster assignment across the ensemble")
    return _save(fig, path)


def plot_dssp_timeline(codes, path, xlabel="frame"):
    """Per-residue (y) × per-frame (x) secondary-structure map, coloured by DSSP category.
    `codes` is the (n_frames, n_residues) int category array from `dssp.dssp_timeline`.

    The figure size is BOUNDED on both axes (cf. the sibling GROMACS pipeline's dssp plot):
    `imshow(aspect="auto")` stretches the data to fill the axes, so the figure need not grow with
    frame count or system size — without the caps the width scaled unbounded (25k frames -> a
    500" sliver)."""
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    from .dssp import DSSP_CATEGORIES

    cmap = ListedColormap([c for _, _, c in DSSP_CATEGORIES])
    n = len(DSSP_CATEGORIES)
    n_frames, n_res = codes.shape
    fig_w = min(max(8.0, n_frames * 0.02), 16.0)
    fig_h = min(max(4.0, n_res * 0.05), 16.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(codes.T, aspect="auto", origin="lower", cmap=cmap,
              vmin=-0.5, vmax=n - 0.5, interpolation="nearest",
              extent=[0, n_frames, 1, n_res + 1])  # residue INDEX on y (robust for multi-chain)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("residue")
    ax.set_title(f"secondary structure ({n_res} residues, {n_frames} frames)")
    handles = [Patch(facecolor=c, label=f"{code.strip() or '~'}  {label}")
               for code, label, c in DSSP_CATEGORIES]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
              fontsize=8, title="DSSP", frameon=False)
    return _save(fig, path)


def plot_ss_fractions(x, fractions, path, xlabel="time (ps)"):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    for name in ("helix", "sheet", "coil", "na"):
        ax.plot(x, fractions[name], lw=1.0, label=name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("fraction of residues")
    ax.set_ylim(0, 1)
    ax.legend(loc="best", fontsize=8)
    ax.set_title("secondary structure")
    return _save(fig, path)
