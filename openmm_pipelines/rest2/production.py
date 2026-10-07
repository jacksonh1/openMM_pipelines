"""rest2.run(eq, cfg) -> Path — the REST2 production backend (openmmtools multistate).

Replica Exchange with Solute Tempering (Wang-Friesner-Berne 2011). Consumes an
`EquilibratedSystem` from `prepare()` and runs Hamiltonian replica exchange where every
replica is at the SAME physical temperature (`cfg.temperature_k`) but replica m's *solute*
is pre-scaled by λ_m (path A of the kickoff note: N independent pre-scaled Systems, one
plain `ThermodynamicState` each, at T_0). Because the replicas share T_0 but carry
different (scaled) Hamiltonians, openmmtools' reduced-potential exchange criterion becomes
exactly the REST2 criterion — only solute DOF enter the acceptance.

Production is **NVT** (the equilibrated box is already at density; barostat stripped).
State 0 is λ=1 — the unscaled, physical ensemble that analysis de-multiplexes.

Writes under `outdir/rest2/`: rest2.nc (+ checkpoint), ladder.json (lambdas +
effective temperatures), run_report.json (QC: exchange acceptance, mixing).

Requires the optional ``[remd]`` extra (openmmtools).
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
from ..lib.selections import protein_atoms
from ..lib.system_edits import strip_barostats
from ..reports import REST2RunReport, write_json
from .ladder import lambda_ladder, lambdas_from_effective_temperatures
from .scaling import scale_solute

_NC = "rest2.nc"


def _ladder(cfg) -> tuple[list[float], list[float]]:
    """`(lambdas, effective_temperatures_k)` from config (explicit λ ladder or geometric)."""
    if cfg.lambdas is not None:
        lambdas = list(cfg.lambdas)
        effective = [cfg.temperature_k / lam for lam in lambdas]
        # Keep the two representations consistent (λ = T_0 / T_eff).
        assert lambdas == lambdas_from_effective_temperatures(cfg.temperature_k, effective)
        return lambdas, effective
    return lambda_ladder(cfg.temperature_k, cfg.max_effective_temperature_k, cfg.n_replicas)


def run(eq, cfg) -> Path:
    outdir = cfg.outdir / "rest2"
    outdir.mkdir(parents=True, exist_ok=True)

    configure_global_platform(cfg.platform)

    lambdas, effective_temperatures = _ladder(cfg)
    solute = protein_atoms(eq.topology, cfg.solute_chain_index)
    write_json(
        {
            "lambdas": lambdas,
            "effective_temperatures_k": effective_temperatures,
            "temperature_k": cfg.temperature_k,
            "n_solute_atoms": len(solute),
        },
        outdir / "ladder.json",
    )

    base = strip_barostats(eq.system)
    # One pre-scaled System per replica; λ=1 reproduces `base` exactly (see scaling tests).
    thermo_states = [
        states.ThermodynamicState(
            system=scale_solute(base, solute, lam),
            temperature=cfg.temperature_k * unit.kelvin,
        )
        for lam in lambdas
    ]
    # One sampler state (the equilibrated pose) replicated; velocities are (re)assigned per
    # replica at the shared T_0 by the move.
    sampler_state = states.SamplerState(positions=eq.positions, box_vectors=eq.box_vectors)

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
        str(outdir / _NC), checkpoint_interval=cfg.checkpoint_interval_iterations
    )
    sampler.create(thermo_states, [sampler_state] * len(lambdas), storage=reporter)

    t0 = time.perf_counter()
    if cfg.equilibration_iterations > 0:
        sampler.equilibrate(cfg.equilibration_iterations)
    sampler.run()
    wall_seconds = time.perf_counter() - t0

    report = _run_report(reporter, cfg, lambdas, effective_temperatures, len(solute), wall_seconds)
    write_json(report, outdir / "run_report.json")
    reporter.close()
    return outdir


def _run_report(reporter, cfg, lambdas, effective_temperatures, n_solute, wall_seconds):
    n_replicas = len(lambdas)
    aggregate_ns = cfg.total_ns * n_replicas
    return REST2RunReport(
        n_replicas=n_replicas,
        temperature_k=cfg.temperature_k,
        lambdas=list(lambdas),
        effective_temperatures_k=list(effective_temperatures),
        n_solute_atoms=n_solute,
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
