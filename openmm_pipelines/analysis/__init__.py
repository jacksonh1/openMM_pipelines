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
from .drift import drift_report
from .dssp import secondary_structure, ss_fractions
from .metrics import radius_of_gyration, rmsd_to_reference, rmsf_per_residue
from .trajectory import (
    export_processed_trajectory,
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
    "export_processed_trajectory",
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
    "drift_report",
    "plots",
]


def analyze(outdir, cluster_cutoff_ang: float = 2.0) -> Path:
    """Full analysis of a completed run: drift (RMSD to design), Rg, per-residue RMSF,
    secondary structure, conformational clusters. Writes `outdir/analysis/`."""
    outdir = Path(outdir)
    adir = outdir / "analysis"
    adir.mkdir(parents=True, exist_ok=True)

    export_processed_trajectory(outdir)  # viewable processed.xtc + processed_topology.pdb

    traj, reference = load_protein_trajectory(outdir)
    time_ps = traj.time

    drift, arrays = drift_report(traj, reference, adir, time_ps, "time (ps)",
                                 "drift from input pose")
    labels = (
        cluster_conformations(traj, cluster_cutoff_ang)
        if traj.n_frames >= 2
        else np.array([1])
    )

    np.savez(adir / "analysis_data.npz", time_ps=time_ps, cluster_labels=labels, **arrays)

    summary = {
        "n_frames": int(traj.n_frames),
        "n_residues": int(traj.topology.n_residues),
        **drift,
        "n_clusters": int(len(set(labels.tolist()))),
    }
    write_json(summary, adir / "analysis_report.json")
    return adir
