"""T-REMD analysis — thin wrappers that pin the REMD run subdir + `.nc` name.

The substance (exchange diagnostics, de-multiplex into a fixed-temperature ensemble, drift
metrics) is shared with REST2 in `analysis/multistate.py`; here we only fix `remd/` +
`remd.nc` and the default reference state 0 (= T_min, the design temperature). See that
module's docstring for the de-multiplex concept.

Reads only; never mutates the run. Requires the `[remd]` extra (openmmtools).
"""

from __future__ import annotations

import mdtraj as md

from . import multistate

_SUBDIR = "remd"
_NC = "remd.nc"


def exchange_diagnostics(remd_dir) -> dict:
    return multistate.exchange_diagnostics(remd_dir, _NC)


def demux_state(remd_dir, state_index: int) -> md.Trajectory:
    return multistate.demux_state(remd_dir, _NC, state_index)


def analyze_remd(outdir, reference_state: int = 0, cluster_cutoff_ang: float = 2.0):
    """Full T-REMD analysis on the fixed-T ensemble of `reference_state` (0 = T_min).

    Writes `outdir/remd/analysis/`: exchange_report.json, mixing plots, the de-multiplexed
    ensemble, the drift report (RMSD/Rg/RMSF/SS vs the input pose), and its conformational
    clusters (Cα-RMSD within `cluster_cutoff_ang` Å).
    """
    return multistate.analyze_multistate(
        outdir, _SUBDIR, _NC, prefix="remd",
        reference_state=reference_state, cluster_cutoff_ang=cluster_cutoff_ang,
    )
