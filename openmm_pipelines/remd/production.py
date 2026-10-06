"""remd.run(eq, cfg) -> Path — the T-REMD production backend (openmmtools multistate).

Consumes an `EquilibratedSystem` (restraints already released) and runs replica exchange
over a geometric temperature ladder. Production is **NVT**: the equilibrated box is
already at the target density, so the barostat is stripped and replicas differ only in
temperature (one `ThermodynamicState` per rung). Each replica equilibrates at its own
temperature first (discarded), then production iterations are recorded to a NetCDF
reporter (`remd.nc`) that stores both the per-replica trajectories AND the
replica<->state permutation at every iteration — the raw material analysis de-multiplexes
into fixed-temperature ensembles.

Writes the produce side of the output contract under `outdir/remd/`:
remd.nc (+ checkpoint), ladder.json, run_report.json (QC: exchange acceptance, mixing).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
from openmm import MonteCarloBarostat, Platform, System, XmlSerializer, unit
from openmmtools import cache, mcmc, states
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeSampler

from ..reports import REMDRunReport, write_json
from .ladder import geometric_ladder

_GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")
_NC = "remd.nc"


def _nvt_system(system: System) -> System:
    """An independent copy of `system` with every MonteCarloBarostat removed.

    Round-trips through the serializer to avoid mutating the shared `eq.system`, then
    drops barostats so production samples NVT at the equilibrated box volume.
    """
    copy = XmlSerializer.deserialize(XmlSerializer.serialize(system))
    # Remove from the back so indices stay valid while deleting.
    for i in reversed(range(copy.getNumForces())):
        if isinstance(copy.getForce(i), MonteCarloBarostat):
            copy.removeForce(i)
    assert not any(
        isinstance(copy.getForce(i), MonteCarloBarostat) for i in range(copy.getNumForces())
    ), "barostat still present after strip — NVT production needs it gone"
    return copy


def _ladder(cfg) -> list[float]:
    if cfg.temperatures_k is not None:
        return list(cfg.temperatures_k)
    return geometric_ladder(cfg.temperature_k, cfg.max_temperature_k, cfg.n_replicas)


def run(eq, cfg) -> Path:
    outdir = cfg.outdir / "remd"
    outdir.mkdir(parents=True, exist_ok=True)

    # openmmtools picks the compute platform off its global context cache, not per-call.
    platform = Platform.getPlatformByName(cfg.platform)
    props = {"Precision": "mixed"} if cfg.platform in _GPU_PLATFORMS else None
    cache.global_context_cache.set_platform(platform, props)

    temperatures = _ladder(cfg)
    write_json({"temperatures_k": temperatures}, outdir / "ladder.json")

    system = _nvt_system(eq.system)
    thermo_states = [
        states.ThermodynamicState(system=system, temperature=t * unit.kelvin)
        for t in temperatures
    ]
    # One sampler state (the equilibrated pose) replicated across replicas; velocities
    # are (re)assigned per replica at its own temperature by the move.
    sampler_state = states.SamplerState(
        positions=eq.positions, box_vectors=eq.box_vectors
    )

    move = mcmc.LangevinDynamicsMove(
        timestep=cfg.dt_ps * unit.picosecond,
        collision_rate=cfg.friction_ps / unit.picosecond,
        n_steps=cfg.steps_per_iteration,
        reassign_velocities=True,
    )
    sampler = ReplicaExchangeSampler(
        mcmc_moves=move,
        number_of_iterations=cfg.n_iterations,
        replica_mixing_scheme=cfg.replica_mixing_scheme,
        online_analysis_interval=None,
    )

    reporter = MultiStateReporter(
        str(outdir / _NC),
        checkpoint_interval=cfg.checkpoint_interval_iterations,
    )
    sampler.create(thermo_states, [sampler_state] * len(temperatures), storage=reporter)

    t0 = time.perf_counter()
    if cfg.equilibration_iterations > 0:
        sampler.equilibrate(cfg.equilibration_iterations)
    sampler.run()
    wall_seconds = time.perf_counter() - t0

    report = _run_report(reporter, cfg, temperatures, wall_seconds)
    write_json(report, outdir / "run_report.json")
    reporter.close()
    return outdir


def _mean_neighbor_acceptance(reporter: MultiStateReporter) -> float:
    """Mean exchange-acceptance over adjacent ladder rungs (i <-> i+1).

    A healthy ladder keeps every neighbour pair swapping; a near-zero rung is a broken
    link that stops replicas diffusing through temperature space.
    """
    # read_mixing_statistics returns per-iteration matrices (n_iter, n_states, n_states);
    # sum over iterations for the cumulative accepted/proposed counts.
    n_accepted, n_proposed = reporter.read_mixing_statistics()
    n_accepted = np.asarray(n_accepted, dtype=float).sum(axis=0)
    n_proposed = np.asarray(n_proposed, dtype=float).sum(axis=0)
    neighbor = [
        n_accepted[i, i + 1] / n_proposed[i, i + 1]
        for i in range(n_accepted.shape[0] - 1)
        if n_proposed[i, i + 1] > 0
    ]
    assert neighbor, "no neighbour exchanges were proposed — mixing statistics empty"
    return float(np.mean(neighbor))


def _run_report(reporter, cfg, temperatures, wall_seconds) -> REMDRunReport:
    n_replicas = len(temperatures)
    aggregate_ns = cfg.total_ns * n_replicas
    return REMDRunReport(
        n_replicas=n_replicas,
        temperatures_k=list(temperatures),
        n_iterations=cfg.n_iterations,
        steps_per_iteration=cfg.steps_per_iteration,
        exchange_attempt_ps=cfg.exchange_attempt_ps,
        total_ns_per_replica=cfg.total_ns,
        replica_mixing_scheme=cfg.replica_mixing_scheme,
        mean_neighbor_acceptance=_mean_neighbor_acceptance(reporter),
        wall_seconds=wall_seconds,
        aggregate_ns_per_day=(
            aggregate_ns / (wall_seconds / 86_400.0) if wall_seconds > 0 else 0.0
        ),
        platform=cfg.platform,
    )
