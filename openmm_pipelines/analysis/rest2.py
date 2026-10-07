"""REST2 analysis — thin wrappers that pin the REST2 run subdir + `.nc` name.

The substance (exchange diagnostics, de-multiplex, drift metrics) is shared with T-REMD in
`analysis/multistate.py`; here we only fix `rest2/` + `rest2.nc` and the default reference
state 0 (= λ=1, the unscaled physical ensemble — the one REST2 ensemble that is a real
Boltzmann sample at the design temperature). See that module's docstring for the demux.

Reads only; never mutates the run. Requires the `[remd]` extra (openmmtools).
"""

from __future__ import annotations

import mdtraj as md

from . import multistate

_SUBDIR = "rest2"
_NC = "rest2.nc"


def exchange_diagnostics(rest2_dir) -> dict:
    return multistate.exchange_diagnostics(rest2_dir, _NC)


def demux_state(rest2_dir, state_index: int) -> md.Trajectory:
    return multistate.demux_state(rest2_dir, _NC, state_index)


def analyze_rest2(outdir, reference_state: int = 0, cluster_cutoff_ang: float = 2.0):
    """Full REST2 analysis on the fixed-λ ensemble of `reference_state` (0 = λ=1).

    Writes `outdir/rest2/analysis/`: exchange_report.json, mixing plots, the de-multiplexed
    λ=1 ensemble, the drift report (RMSD/Rg/RMSF/SS vs the input pose), and its conformational
    clusters (Cα-RMSD within `cluster_cutoff_ang` Å).
    """
    return multistate.analyze_multistate(
        outdir, _SUBDIR, _NC, prefix="rest2",
        reference_state=reference_state, cluster_cutoff_ang=cluster_cutoff_ang,
    )
