"""analysis — mdtraj-based trajectory analysis, all measured vs the input design pose.

Entry points: the single-purpose functions (`load_protein_trajectory`,
`rmsd_to_reference`, `radius_of_gyration`, `rmsf_per_residue`, `ss_fractions`,
`cluster_conformations`) and `analyze(outdir)`, which runs them all and writes
`outdir/analysis/` (report JSON, data npz, PNGs).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..reports import write_json
from . import plots
from .clustering import cluster_conformations, pairwise_rmsd
from .dssp import secondary_structure, ss_fractions
from .metrics import radius_of_gyration, rmsd_to_reference, rmsf_per_residue
from .trajectory import (
    load_input_reference,
    load_protein_trajectory,
    load_system_trajectory,
    superpose_to_reference,
    to_protein,
)

__all__ = [
    "analyze",
    "load_protein_trajectory",
    "load_system_trajectory",
    "to_protein",
    "load_input_reference",
    "superpose_to_reference",
    "rmsd_to_reference",
    "radius_of_gyration",
    "rmsf_per_residue",
    "secondary_structure",
    "ss_fractions",
    "cluster_conformations",
    "pairwise_rmsd",
    "plots",
]


def analyze(outdir, cluster_cutoff_ang: float = 2.0) -> Path:
    """Full analysis of a completed run: drift (RMSD to design), Rg, per-residue RMSF,
    secondary structure, conformational clusters. Writes `outdir/analysis/`."""
    outdir = Path(outdir)
    adir = outdir / "analysis"
    adir.mkdir(parents=True, exist_ok=True)

    traj, reference = load_protein_trajectory(outdir)
    time_ps = traj.time

    rmsd = rmsd_to_reference(traj, reference)
    rg = radius_of_gyration(traj)
    resseq, rmsf = rmsf_per_residue(traj)
    fractions = ss_fractions(traj)
    labels = (
        cluster_conformations(traj, cluster_cutoff_ang)
        if traj.n_frames >= 2
        else np.array([1])
    )

    plots.plot_timeseries(time_ps, rmsd, "Cα RMSD to design (Å)",
                          adir / "rmsd.png", "drift from input pose")
    plots.plot_timeseries(time_ps, rg, "Rg (Å)", adir / "rg.png", "radius of gyration")
    plots.plot_rmsf(resseq, rmsf, adir / "rmsf.png")
    plots.plot_ss_fractions(time_ps, fractions, adir / "dssp.png")

    np.savez(
        adir / "analysis_data.npz",
        time_ps=time_ps, rmsd_ang=rmsd, rg_ang=rg,
        resseq=resseq, rmsf_ang=rmsf, cluster_labels=labels,
        **{f"ss_{name}": frac for name, frac in fractions.items()},
    )

    summary = {
        "n_frames": int(traj.n_frames),
        "n_residues": int(traj.topology.n_residues),
        "rmsd_to_design_ang": {
            "mean": float(rmsd.mean()), "max": float(rmsd.max()), "final": float(rmsd[-1]),
        },
        "rg_ang": {"mean": float(rg.mean()), "std": float(rg.std())},
        "rmsf_ang": {
            "mean": float(rmsf.mean()), "max": float(rmsf.max()),
            "max_residue": int(resseq[int(rmsf.argmax())]),
        },
        "ss_fraction_mean": {name: float(frac.mean()) for name, frac in fractions.items()},
        "n_clusters": int(len(set(labels.tolist()))),
    }
    write_json(summary, adir / "analysis_report.json")
    return adir
