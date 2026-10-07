"""Shared analysis for the openmmtools multistate backends (T-REMD and REST2).

Both backends write the same NetCDF layout (a `MultiStateReporter` storing per-replica
trajectories plus the replica<->state permutation at every iteration), so the exchange
diagnostics and the de-multiplex are identical — only the run subdirectory and the `.nc`
filename differ. `analysis/remd.py` and `analysis/rest2.py` are thin wrappers that pin
those two names; everything substantive lives here.

The one concept that makes or breaks this analysis: a reporter stores **replica**
trajectories — continuous coordinate sets that *walk through state space* as exchanges
swap their thermodynamic state. A single replica's frames are NOT an ensemble at one state;
they are a mix of every state the replica visited. To measure the design's structural drift
at the reference state you must **de-multiplex**: at each iteration pick the replica
currently occupying the target state and collect *those* frames. For T-REMD the reference
is state 0 = T_min; for REST2 it is state 0 = λ=1 (the unscaled, physical ensemble). That
collection is the canonical fixed-state (NVT) ensemble, and only it should feed RMSD/Rg/RMSF.

Reads only; never mutates the run. Requires the `[remd]` extra (openmmtools).
"""

from __future__ import annotations

from pathlib import Path

import mdtraj as md
import numpy as np
from openmm import unit
from openmm.app import PDBxFile
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeAnalyzer

from ..exchange import neighbor_acceptance
from ..reports import write_json
from . import plots
from .clustering import cluster_conformations
from .drift import drift_report
from .trajectory import load_input_reference, superpose_to_reference, to_protein


def _open_reporter(run_dir: Path, nc_name: str) -> MultiStateReporter:
    nc = run_dir / nc_name
    assert nc.exists(), f"no multistate storage at {nc} — run the backend first"
    return MultiStateReporter(str(nc), open_mode="r")


def exchange_diagnostics(run_dir, nc_name: str) -> dict:
    """Mixing health from the full-resolution state record (every iteration).

    `ReplicaExchangeAnalyzer.generate_mixing_statistics` gives the empirical state-transition
    matrix, its eigenvalues and the state-index statistical inefficiency. The **subdominant
    eigenvalue** (second largest) is the canonical scalar: the spectral gap `1 - eigenvalue[1]`
    sets how fast replicas decorrelate across state space (near 1 = poor mixing). We add
    `neighbor_acceptance` (accepted/proposed per adjacent rung) as the directly ladder-tunable
    signal, and carry the replica->state permutation for the walk plot and the de-multiplex.
    """
    run_dir = Path(run_dir)
    reporter = _open_reporter(run_dir, nc_name)
    try:
        analyzer = ReplicaExchangeAnalyzer(reporter)
        mixing = analyzer.generate_mixing_statistics(number_equilibrated=0)
        acceptance = neighbor_acceptance(reporter)  # per adjacent rung (shared primitive)
        # replica_states[i, r] = state index occupied by replica r at iteration i
        replica_states = np.asarray(reporter.read_replica_thermodynamic_states())
    finally:
        reporter.close()

    eigenvalues = np.real(np.asarray(mixing.eigenvalues))
    return {
        "n_states": len(acceptance) + 1,
        "n_iterations": int(replica_states.shape[0]),
        "neighbor_acceptance": acceptance,
        "mean_neighbor_acceptance": float(np.mean(acceptance)),
        "min_neighbor_acceptance": float(np.min(acceptance)),
        "subdominant_eigenvalue": float(eigenvalues[1]),  # mixing: spectral gap = 1 - this
        "state_statistical_inefficiency": float(mixing.statistical_inefficiency),
        "replica_states": replica_states,  # not JSON-written (large); for plots/demux
    }


def demux_state(run_dir, nc_name: str, state_index: int) -> md.Trajectory:
    """Fixed-state ensemble for one ladder rung, as a full-system trajectory.

    For each checkpointed iteration, takes the positions of the replica occupying
    `state_index`. Frame resolution = the run's checkpoint interval (positions live in the
    checkpoint file, not the per-iteration analysis file).
    """
    run_dir = Path(run_dir)
    reporter = _open_reporter(run_dir, nc_name)
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

    # prepare() writes the topology under <outdir>/prepared; run_dir is <outdir>/<subdir>.
    cif = run_dir.parent / "prepared" / "topology.cif"
    assert cif.exists(), f"no topology at {cif} — expected prepare() output alongside {run_dir.name}/"
    topology = md.Topology.from_openmm(PDBxFile(str(cif)).topology)

    # mdtraj stores coordinates/box as float32.
    traj = md.Trajectory(xyz=np.asarray(frames, dtype=np.float32), topology=topology)
    traj.unitcell_vectors = np.asarray(boxes, dtype=np.float32)  # triclinic box per frame
    protein = topology.select("protein")
    assert len(protein) > 0, "no protein atoms in topology — cannot anchor imaging"
    anchor = {topology.atom(int(i)) for i in protein}
    traj.image_molecules(inplace=True, anchor_molecules=[anchor])
    return traj


def analyze_multistate(
    outdir,
    subdir: str,
    nc_name: str,
    prefix: str,
    reference_state: int = 0,
    cluster_cutoff_ang: float = 2.0,
) -> Path:
    """Exchange/mixing diagnostics + drift metrics + conformational clustering on the
    fixed-state ensemble of `reference_state` (default 0 = the design-reference rung: T_min
    for REMD, λ=1 for REST2). That demuxed ensemble is a Boltzmann sample at the design
    temperature, so clustering it (Cα-RMSD, `cluster_cutoff_ang` Å) partitions the sampled
    conformations. Writes `outdir/<subdir>/analysis/`.
    """
    outdir = Path(outdir)
    run_dir = outdir / subdir
    adir = run_dir / "analysis"
    adir.mkdir(parents=True, exist_ok=True)

    diag = exchange_diagnostics(run_dir, nc_name)
    replica_states = diag.pop("replica_states")
    write_json(diag, adir / "exchange_report.json")

    plots.plot_neighbor_acceptance(diag["neighbor_acceptance"], adir / "acceptance.png")
    plots.plot_replica_state_walk(replica_states, adir / "replica_walk.png")

    # Fixed-state ensemble at the reference rung, measured against the design pose. Fully
    # processed for both metrics and viewing: PBC-imaged (in demux_state), desolvated
    # (to_protein), aligned to the design pose, then origin-centered on the frame-0 centroid.
    # This writes new analysis artifacts only; the raw sim output (the .nc) is read-only.
    traj = to_protein(demux_state(run_dir, nc_name, reference_state))
    reference = load_input_reference(outdir)
    superpose_to_reference(traj, reference)
    traj.xyz -= traj.xyz[0].mean(axis=0)  # center frame-0 geometric center at the origin

    traj.save_xtc(str(adir / f"state{reference_state:02d}_ensemble.xtc"))
    traj[0].save_pdb(str(adir / f"state{reference_state:02d}_topology.pdb"))

    # Demuxed frames are an ensemble ordered by iteration, not a continuous time series.
    frame_index = np.arange(traj.n_frames, dtype=float)
    drift, arrays = drift_report(
        traj, reference, adir, frame_index, "demuxed frame",
        f"reference ensemble (state {reference_state})",
    )

    # Conformational clusters of the demuxed ensemble (Cα-RMSD, average-linkage).
    cluster_labels = (
        cluster_conformations(traj, cluster_cutoff_ang)
        if traj.n_frames >= 2
        else np.array([1])
    )

    np.savez(
        adir / f"{prefix}_analysis_data.npz",
        frame_index=frame_index, replica_states=replica_states,
        neighbor_acceptance=np.asarray(diag["neighbor_acceptance"]),
        cluster_labels=cluster_labels, **arrays,
    )

    summary = {
        "reference_state": reference_state,
        "mean_neighbor_acceptance": diag["mean_neighbor_acceptance"],
        "min_neighbor_acceptance": diag["min_neighbor_acceptance"],
        "subdominant_eigenvalue": diag["subdominant_eigenvalue"],
        "state_statistical_inefficiency": diag["state_statistical_inefficiency"],
        "n_frames_demuxed": int(traj.n_frames),
        "cluster_cutoff_ang": cluster_cutoff_ang,
        "n_clusters": int(len(set(cluster_labels.tolist()))),
        **drift,
    }
    write_json(summary, adir / f"{prefix}_analysis_report.json")
    return adir
