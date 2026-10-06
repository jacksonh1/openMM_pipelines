"""Drift-from-design report — the shared core of every analysis.

RMSD / Rg / per-residue RMSF / secondary-structure fractions, all measured against the
input design pose, plus their four figures. This is the common block between plain MD
(`analyze`) and the de-multiplexed REMD ensemble (`analyze_remd`); each composer wraps it
with its own extras (clustering, exchange diagnostics) and its own x-axis (production time
vs demuxed-frame index).

No OpenMM / openmmtools import — pure mdtraj, so it is shared freely across engines.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import plots
from .dssp import ss_fractions
from .metrics import radius_of_gyration, rmsd_to_reference, rmsf_per_residue


def drift_report(traj, reference, adir, x, xlabel, rmsd_title) -> tuple[dict, dict]:
    """Compute the four drift metrics vs `reference`, write their PNGs into `adir`, and
    return ``(summary, arrays)``.

    `x` is the per-frame x-axis (production time in ps, or demuxed-frame index) and
    `xlabel` its label. `summary` is the JSON-ready stats block; `arrays` is the raw
    per-frame data the caller folds into its own ``.npz`` alongside engine-specific fields.
    """
    adir = Path(adir)
    rmsd = rmsd_to_reference(traj, reference)
    rg = radius_of_gyration(traj)
    resseq, rmsf = rmsf_per_residue(traj)
    fractions = ss_fractions(traj)

    plots.plot_timeseries(x, rmsd, "Cα RMSD to design (Å)", adir / "rmsd.png",
                          rmsd_title, xlabel=xlabel)
    plots.plot_timeseries(x, rg, "Rg (Å)", adir / "rg.png",
                          "radius of gyration", xlabel=xlabel)
    plots.plot_rmsf(resseq, rmsf, adir / "rmsf.png")
    plots.plot_ss_fractions(x, fractions, adir / "dssp.png", xlabel=xlabel)

    summary = {
        "rmsd_to_design_ang": {
            "mean": float(rmsd.mean()), "max": float(rmsd.max()), "final": float(rmsd[-1]),
        },
        "rg_ang": {"mean": float(rg.mean()), "std": float(rg.std())},
        "rmsf_ang": {
            "mean": float(rmsf.mean()), "max": float(rmsf.max()),
            "max_residue": int(resseq[int(rmsf.argmax())]),
        },
        "ss_fraction_mean": {name: float(frac.mean()) for name, frac in fractions.items()},
    }
    arrays = {
        "rmsd_ang": rmsd, "rg_ang": rg, "resseq": resseq, "rmsf_ang": rmsf,
        **{f"ss_{name}": frac for name, frac in fractions.items()},
    }
    return summary, arrays
