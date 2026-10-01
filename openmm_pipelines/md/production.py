"""md.run(eq, cfg) -> Path — the plain-MD production backend.

Consumes an `EquilibratedSystem` (restraints already released) and runs one trajectory
via `app.Simulation.step()`, writing the production side of the output contract:
trajectory.xtc, production.log, final_state.xml (portable), production.chk (exact,
same-hardware restart), and run_report.json (QC). Production starts exactly at
restraint release, so the trajectory records all drift from the input pose.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import CheckpointReporter, Simulation, StateDataReporter, XTCReporter

from ..reports import RunReport, write_json

_GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")


def run(eq, cfg) -> Path:
    prod = cfg.outdir / "production"
    prod.mkdir(parents=True, exist_ok=True)

    integ = LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        cfg.friction_ps / unit.picosecond,
        cfg.dt_ps * unit.picosecond,
    )
    if cfg.seed is not None:
        integ.setRandomNumberSeed(cfg.seed)
    platform = Platform.getPlatformByName(cfg.platform)
    props = {"Precision": "mixed"} if cfg.platform in _GPU_PLATFORMS else {}

    sim = Simulation(eq.topology, eq.system, integ, platform, props)
    sim.context.setPeriodicBoxVectors(*eq.box_vectors)
    sim.context.setPositions(eq.positions)
    sim.context.setVelocitiesToTemperature(
        cfg.temperature_k * unit.kelvin, cfg.seed if cfg.seed is not None else 0
    )

    log_path = prod / "production.log"
    sim.reporters.append(XTCReporter(str(prod / "trajectory.xtc"), cfg.traj_steps))
    sim.reporters.append(
        StateDataReporter(
            str(log_path), cfg.traj_steps, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            speed=True, totalSteps=cfg.prod_steps,
        )
    )
    sim.reporters.append(
        CheckpointReporter(str(prod / "production.chk"), cfg.checkpoint_steps)
    )

    t0 = time.perf_counter()
    sim.step(cfg.prod_steps)
    wall_seconds = time.perf_counter() - t0

    sim.saveState(str(prod / "final_state.xml"))  # portable restart
    sim.saveCheckpoint(str(prod / "production.chk"))  # exact, same-hardware restart

    write_json(_run_report(log_path, cfg, wall_seconds), prod / "run_report.json")
    return prod


def _col(df: pd.DataFrame, needle: str) -> pd.Series:
    """Pick the StateDataReporter column whose verbose header contains `needle`."""
    matches = [c for c in df.columns if needle.lower() in c.lower()]
    assert matches, f"no column matching {needle!r} in {list(df.columns)}"
    return df[matches[0]]


def _run_report(log_path: Path, cfg, wall_seconds: float) -> RunReport:
    df = pd.read_csv(log_path)
    df.columns = [c.lstrip('#').strip().strip('"') for c in df.columns]

    pe = _col(df, "Potential Energy").to_numpy()
    tot = _col(df, "Total Energy").to_numpy()
    temp = _col(df, "Temperature").to_numpy()
    dens = _col(df, "Density").to_numpy()

    target_t = cfg.temperature_k
    stable = bool(
        np.all(np.isfinite(pe)) and np.all(np.isfinite(tot))
        and abs(float(temp.mean()) - target_t) < 25.0
    )

    return RunReport(
        n_steps=cfg.prod_steps,
        total_ns=cfg.total_ns,
        n_frames=len(df),
        wall_seconds=wall_seconds,
        ns_per_day=cfg.total_ns / (wall_seconds / 86_400.0) if wall_seconds > 0 else 0.0,
        temperature_mean_k=float(temp.mean()),
        temperature_std_k=float(temp.std()),
        density_mean=float(dens.mean()),
        density_std=float(dens.std()),
        potential_energy_mean_kj_mol=float(pe.mean()),
        potential_energy_min_kj_mol=float(pe.min()),
        potential_energy_max_kj_mol=float(pe.max()),
        total_energy_mean_kj_mol=float(tot.mean()),
        total_energy_std_kj_mol=float(tot.std()),
        stable=stable,
        platform=cfg.platform,
    )
