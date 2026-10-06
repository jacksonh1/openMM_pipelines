"""T-REMD analysis — de-multiplex replica trajectories into fixed-temperature ensembles.

The one concept that makes or breaks REMD analysis: a NetCDF reporter stores **replica**
trajectories — continuous coordinate sets that *walk through temperature space* as
exchanges swap their thermodynamic state. A single replica's frames are NOT an ensemble
at one temperature; they are a mix of every temperature the replica visited. To measure
structural drift of the design at its reference temperature you must **de-multiplex**:
at each iteration pick the replica currently occupying the target state (e.g. state 0 =
T_min) and collect *those* frames. That collection IS the canonical fixed-temperature
(NVT) ensemble, and only it should feed RMSD/Rg/RMSF.

This module provides both views:
- `exchange_diagnostics` — per-neighbour acceptance, round-trip count, and the
  replica->state permutation over time (proof the replicas mix), from the full-resolution
  state record.
- `demux_state` — the fixed-T ensemble for one ladder rung, as an mdtraj trajectory built
  from the checkpointed frames, imaged/stripped/superposed exactly like plain MD.
- `analyze_remd(outdir)` — runs both, then the standard drift metrics on the T_min
  ensemble, and writes `outdir/remd/analysis/`.

Reads only; never mutates the run. Requires the `[remd]` extra (openmmtools).
"""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
from openmm import unit
from openmm.app import PDBxFile
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer

from ..reports import write_json
from . import plots
from .drift import drift_report
from .trajectory import load_input_reference, superpose_to_reference, to_protein

_NC = "remd.nc"


def _open_reporter(remd_dir: Path) -> MultiStateReporter:
    nc = remd_dir / _NC
    assert nc.exists(), f"no REMD storage at {nc} — run remd.run first"
    return MultiStateReporter(str(nc), open_mode="r")


def exchange_diagnostics(remd_dir) -> dict:
    """Mixing health from the full-resolution state record (every iteration).

    Mixing quality comes from openmmtools' `ReplicaExchangeAnalyzer.generate_mixing_statistics`
    — the empirical state-transition matrix, its eigenvalues, and the state-index
    statistical inefficiency. The **subdominant eigenvalue** (second largest) is the
    canonical scalar: the spectral gap `1 - eigenvalue[1]` sets how fast replicas
    decorrelate across temperature space (near 1 = poor mixing). We add `neighbor_acceptance`
    (accepted/proposed per adjacent rung) as the simple, directly ladder-tunable signal, and
    carry the replica->state permutation matrix for the walk plot and the de-multiplex.
    """
    remd_dir = Path(remd_dir)
    reporter = _open_reporter(remd_dir)
    try:
        analyzer = ReplicaExchangeAnalyzer(reporter)
        mixing = analyzer.generate_mixing_statistics(number_equilibrated=0)
        n_accepted, n_proposed = reporter.read_mixing_statistics()
        # replica_states[i, r] = state index occupied by replica r at iteration i
        replica_states = np.asarray(reporter.read_replica_thermodynamic_states())
    finally:
        reporter.close()

    # Per-iteration matrices (n_iter, n_states, n_states); sum to cumulative counts.
    n_accepted = np.asarray(n_accepted, dtype=float).sum(axis=0)
    n_proposed = np.asarray(n_proposed, dtype=float).sum(axis=0)
    n_states = n_accepted.shape[0]
    neighbor_acceptance = [
        float(n_accepted[i, i + 1] / n_proposed[i, i + 1])
        if n_proposed[i, i + 1] > 0 else 0.0
        for i in range(n_states - 1)
    ]
    eigenvalues = np.real(np.asarray(mixing.eigenvalues))
    return {
        "n_states": n_states,
        "n_iterations": int(replica_states.shape[0]),
        "neighbor_acceptance": neighbor_acceptance,
        "mean_neighbor_acceptance": float(np.mean(neighbor_acceptance)),
        "min_neighbor_acceptance": float(np.min(neighbor_acceptance)),
        "subdominant_eigenvalue": float(eigenvalues[1]),  # mixing: spectral gap = 1 - this
        "state_statistical_inefficiency": float(mixing.statistical_inefficiency),
        "replica_states": replica_states,  # not JSON-written (large); for plots/demux
    }


def demux_state(remd_dir, state_index: int) -> md.Trajectory:
    """Fixed-temperature ensemble for one ladder rung, as a full-system trajectory.

    For each checkpointed iteration, takes the positions of the replica occupying
    `state_index`. Frame resolution = the run's checkpoint interval (positions live in
    the checkpoint file, not the per-iteration analysis file).
    """
    remd_dir = Path(remd_dir)
    reporter = _open_reporter(remd_dir)
    try:
        replica_states = np.asarray(reporter.read_replica_thermodynamic_states())
        checkpoints = list(reporter.read_checkpoint_iterations())
        assert checkpoints, "no checkpointed frames to de-multiplex"

        frames = []
        boxes = []
        for it in checkpoints:
            sampler_states = reporter.read_sampler_states(it)
            assert sampler_states is not None, f"no positions stored at iteration {it}"
            # Exactly one replica occupies this state per iteration.
            matches = np.where(replica_states[it] == state_index)[0]
            assert matches.size == 1, (
                f"state {state_index} not uniquely occupied at iteration {it}: {matches}"
            )
            ss = sampler_states[int(matches[0])]
            frames.append(ss.positions.value_in_unit(unit.nanometer))
            boxes.append(ss.box_vectors.value_in_unit(unit.nanometer))
    finally:
        reporter.close()

    # prepare() writes the topology under <outdir>/prepared; remd_dir is <outdir>/remd.
    cif = remd_dir.parent / "prepared" / "topology.cif"
    assert cif.exists(), f"no topology at {cif} — expected prepare() output alongside remd/"
    topology = md.Topology.from_openmm(PDBxFile(str(cif)).topology)

    # mdtraj stores coordinates/box as float32.
    traj = md.Trajectory(xyz=np.asarray(frames, dtype=np.float32), topology=topology)
    traj.unitcell_vectors = np.asarray(boxes, dtype=np.float32)  # triclinic box per frame
    protein = topology.select("protein")
    assert len(protein) > 0, "no protein atoms in topology — cannot anchor imaging"
    anchor = {topology.atom(int(i)) for i in protein}
    traj.image_molecules(inplace=True, anchor_molecules=[anchor])
    return traj


def analyze_remd(outdir, reference_state: int = 0, cluster_cutoff_ang: float = 2.0) -> Path:
    """Full REMD analysis: exchange/mixing diagnostics + drift metrics on the fixed-T
    ensemble of `reference_state` (default 0 = T_min, the design-reference temperature).

    Writes `outdir/remd/analysis/`: exchange_report.json, mixing plots, the de-multiplexed
    T_min trajectory, and the standard drift report (RMSD/Rg/RMSF/SS vs the input pose).
    """
    outdir = Path(outdir)
    remd_dir = outdir / "remd"
    adir = remd_dir / "analysis"
    adir.mkdir(parents=True, exist_ok=True)

    diag = exchange_diagnostics(remd_dir)
    replica_states = diag.pop("replica_states")
    write_json(diag, adir / "exchange_report.json")

    plots.plot_neighbor_acceptance(diag["neighbor_acceptance"], adir / "acceptance.png")
    plots.plot_replica_state_walk(replica_states, adir / "replica_walk.png")

    # Fixed-T ensemble at the reference rung, measured against the design pose.
    traj = to_protein(demux_state(remd_dir, reference_state))
    reference = load_input_reference(outdir)
    superpose_to_reference(traj, reference)

    traj.save_xtc(str(adir / f"state{reference_state:02d}_ensemble.xtc"))
    traj[0].save_pdb(str(adir / f"state{reference_state:02d}_topology.pdb"))

    # Demuxed frames are an ensemble ordered by iteration, not a continuous time series.
    frame_index = np.arange(traj.n_frames, dtype=float)
    drift, arrays = drift_report(
        traj, reference, adir, frame_index, "demuxed frame",
        f"T_min ensemble (state {reference_state})",
    )

    np.savez(
        adir / "remd_analysis_data.npz",
        frame_index=frame_index, replica_states=replica_states,
        neighbor_acceptance=np.asarray(diag["neighbor_acceptance"]), **arrays,
    )

    summary = {
        "reference_state": reference_state,
        "mean_neighbor_acceptance": diag["mean_neighbor_acceptance"],
        "min_neighbor_acceptance": diag["min_neighbor_acceptance"],
        "subdominant_eigenvalue": diag["subdominant_eigenvalue"],
        "state_statistical_inefficiency": diag["state_statistical_inefficiency"],
        "n_frames_demuxed": int(traj.n_frames),
        **drift,
    }
    write_json(summary, adir / "remd_analysis_report.json")
    return adir
