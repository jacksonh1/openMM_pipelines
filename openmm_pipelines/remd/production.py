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
from openmm import unit
from openmmtools import mcmc, states
from openmmtools.multistate import MultiStateReporter, ReplicaExchangeSampler

from ..context_cache import configure_global_platform
from ..exchange import neighbor_acceptance
from ..lib.system_edits import strip_barostats
from ..reports import REMDRunReport, write_json
from .ladder import geometric_ladder

_NC = "remd.nc"


def _ladder(cfg) -> list[float]:
    if cfg.temperatures_k is not None:
        return list(cfg.temperatures_k)
    return geometric_ladder(cfg.temperature_k, cfg.max_temperature_k, cfg.n_replicas)


def run(eq, cfg) -> Path:
    outdir = cfg.outdir / "remd"
    outdir.mkdir(parents=True, exist_ok=True)

    configure_global_platform(cfg.platform)

    temperatures = _ladder(cfg)
    write_json({"temperatures_k": temperatures}, outdir / "ladder.json")

    system = strip_barostats(eq.system)
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
        mean_neighbor_acceptance=float(np.mean(neighbor_acceptance(reporter))),
        wall_seconds=wall_seconds,
        aggregate_ns_per_day=(
            aggregate_ns / (wall_seconds / 86_400.0) if wall_seconds > 0 else 0.0
        ),
        platform=cfg.platform,
    )
