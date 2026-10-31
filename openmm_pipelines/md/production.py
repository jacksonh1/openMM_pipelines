"""md.run(eq, cfg) -> Path and md.resume(cfg) -> Path — the plain-MD production backend.

`run` consumes an `EquilibratedSystem` (restraints already released) and runs one trajectory
via `app.Simulation.step()`, writing the production side of the output contract:
trajectory.xtc, production.log, final_state.xml (portable), run_report.json (QC), plus a
pair of ping-pong resume checkpoints. Production starts exactly at restraint release, so the
trajectory records all drift from the input pose.

`resume` continues an interrupted (walltime/crash) run — or extends it to a larger target —
from the latest ping-pong checkpoint, WITHOUT re-solvating or re-equilibrating. It reloads
`prepared/` via `EquilibratedSystem.from_disk`, restores the portable State, truncates the
trajectory + log back to the checkpoint step (so the append is contiguous — no duplicate or
divergent frames), and runs the remaining steps with `append=True` reporters.

Resume correctness rests on the portable State XML (`simulation.loadState`): it works on any
node but drops the RNG stream, so the continued trajectory differs from an uninterrupted run
(physically valid — GPU runs are not bit-reproducible anyway). The resume point is written
ping-pong to two alternating files so a SIGKILL mid-write never destroys the prior good one.
(An exact same-hardware `.chk` fast path, gated by an environment fingerprint, is a documented
future enhancement — NOT implemented here; resume is portable-State-only.)
"""

from __future__ import annotations

import datetime as _dt
import time
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import Simulation, StateDataReporter, XTCReporter

from ..lib.platforms import gpu_props
from ..lib.run_layout import RunLayout, write_manifest
from ..reports import RunInfo, RunReport, write_json

if TYPE_CHECKING:
    from ..preparation.handoff import EquilibratedSystem
    from .config import MDConfig

_TRAJ = "production/trajectory.xtc"
_LOG = "production/production.log"
_FINAL = "production/final_state.xml"
# Two alternating resume checkpoints (ping-pong): a torn write never destroys the prior good file.
_CKPT_A = "production/checkpoint_state.A.xml"
_CKPT_B = "production/checkpoint_state.B.xml"


class _PingPongStateReporter:
    """Write a portable State XML every `reportInterval` steps, ALTERNATING two files.

    OpenMM's `CheckpointReporter` overwrites a single file in place, so a SIGKILL mid-write would
    leave no valid resume point. Alternating two files guarantees the previous good checkpoint
    survives any single interrupted write; `resume` loads the newest of the two that deserializes.
    Mirrors the OpenMM reporter protocol (`describeNextReport`/`report`); `saveState` pulls its own
    state, so the report requests no state data.
    """

    def __init__(self, path_a, path_b, reportInterval: int):
        assert reportInterval >= 1, "reportInterval must be >= 1"
        self._paths = (str(path_a), str(path_b))
        self._interval = int(reportInterval)
        self._next = 0

    def describeNextReport(self, simulation):
        steps = self._interval - simulation.currentStep % self._interval
        return {"steps": steps, "periodic": None, "include": []}

    def report(self, simulation, state):
        simulation.saveState(self._paths[self._next])
        self._next ^= 1


def _manifest(cfg) -> dict:
    return {"engine": "md", "production_config": cfg.model_dump(mode="json")}


def _make_simulation(eq: EquilibratedSystem, cfg: MDConfig) -> Simulation:
    """Build the production `Simulation` (box set, positions/velocities left to the caller: fresh
    seeds them, resume `loadState`s them)."""
    integ = LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        cfg.friction_ps / unit.picosecond,
        cfg.dt_ps * unit.picosecond,
    )
    if cfg.seed is not None:
        integ.setRandomNumberSeed(cfg.seed)
    platform = Platform.getPlatformByName(cfg.platform)
    sim = Simulation(eq.topology, eq.system, integ, platform, gpu_props(cfg.platform))
    sim.context.setPeriodicBoxVectors(*eq.box_vectors)
    return sim


def _add_reporters(sim: Simulation, layout: RunLayout, cfg: MDConfig, *, append: bool) -> Path:
    """Register production reporters and return the real log path.

    The ping-pong checkpoint reporter is registered BEFORE the trajectory/log reporters, so at a
    shared frame-boundary step the checkpoint is written first: a crash between them then loses at
    most one not-yet-written frame (a gap), never a duplicated/divergent frame. For a fresh run the
    trajectory is offloaded (symlink into storage) and the log recorded; on resume both reopen the
    existing real files (`existing`) in append mode.
    """
    sim.reporters.append(
        _PingPongStateReporter(layout.record(_CKPT_A), layout.record(_CKPT_B), cfg.checkpoint_steps)
    )
    traj_path = layout.existing(_TRAJ) if append else layout.offload(_TRAJ)
    log_path = layout.existing(_LOG) if append else layout.record(_LOG)
    sim.reporters.append(XTCReporter(str(traj_path), cfg.traj_steps, append=append))
    sim.reporters.append(
        StateDataReporter(
            str(log_path), cfg.traj_steps, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            speed=True, totalSteps=cfg.prod_steps, append=append,
        )
    )
    return log_path


def run(eq: EquilibratedSystem, cfg: MDConfig) -> Path:
    layout = RunLayout(cfg.outdir, cfg.storage_directory)
    prod = cfg.outdir / "production"
    prod.mkdir(parents=True, exist_ok=True)
    write_manifest(cfg.outdir, _manifest(cfg))  # authoritative engine/target record, for resume

    sim = _make_simulation(eq, cfg)
    sim.context.setPositions(eq.positions)
    sim.context.setVelocitiesToTemperature(
        cfg.temperature_k * unit.kelvin, cfg.seed if cfg.seed is not None else 0
    )
    log_path = _add_reporters(sim, layout, cfg, append=False)

    try:
        t0 = time.perf_counter()
        sim.step(cfg.prod_steps)
        wall_seconds = time.perf_counter() - t0

        sim.saveState(str(layout.record(_FINAL)))  # portable end state
        write_json(_run_report(log_path, cfg, wall_seconds), layout.record("production/run_report.json"))
    finally:
        layout.mirror("production")  # copy the small files into storage (storage = full record)
    return prod


def resume(cfg: MDConfig) -> Path:
    """Continue (or extend to a larger `total_ns`) an interrupted MD run from `cfg.outdir`."""
    from ..preparation.handoff import EquilibratedSystem

    layout = RunLayout(cfg.outdir, cfg.storage_directory)
    eq = EquilibratedSystem.from_disk(cfg.outdir / "prepared")
    write_manifest(cfg.outdir, _manifest(cfg))  # restamp the (possibly extended) target

    sim = _make_simulation(eq, cfg)
    elapsed_ps = _load_latest_state(sim, [cfg.outdir / _CKPT_A, cfg.outdir / _CKPT_B])
    start_step = round(elapsed_ps / cfg.dt_ps)
    sim.currentStep = start_step  # restore the global step (fixes the log Step column + reporter phase)

    if cfg.prod_steps <= start_step:
        raise RuntimeError(
            f"run at step {start_step} ({elapsed_ps / 1000:.4f} ns) already >= target "
            f"{cfg.prod_steps} ({cfg.total_ns} ns); pass a larger extend_to_ns to extend."
        )

    # Heal the trajectory + log to the checkpoint step so the append is contiguous (the checkpoint
    # step is authoritative; frames/rows written after it, pre-crash, would otherwise overlap).
    _truncate_xtc(layout.existing(_TRAJ), start_step // cfg.traj_steps)
    _truncate_log(layout.existing(_LOG), start_step)

    log_path = _add_reporters(sim, layout, cfg, append=True)

    try:
        t0 = time.perf_counter()
        sim.step(cfg.prod_steps - start_step)
        wall_seconds = time.perf_counter() - t0

        sim.saveState(str(layout.record(_FINAL)))  # portable end state
        _write_resume_reports(layout, cfg, start_step, wall_seconds, log_path)
    finally:
        layout.mirror("production")
    return cfg.outdir / "production"


def _load_latest_state(sim: Simulation, candidates) -> float:
    """Load the newest ping-pong State XML that deserializes into `sim`; return its sim time (ps).

    Falls back to the older file if the newest is torn (a documented, intentional recovery for the
    SIGKILL-mid-write case). Fails loud if none exist or none are loadable.
    """
    present = [Path(p) for p in candidates if Path(p).exists()]
    if not present:
        raise FileNotFoundError(f"no resume checkpoint among {list(candidates)}; nothing to resume.")
    present.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    errors: dict[str, str] = {}
    for p in present:
        try:
            sim.loadState(str(p))
        except Exception as exc:  # torn/invalid checkpoint: try the older ping-pong file
            errors[str(p)] = repr(exc)
            continue
        return sim.context.getState().getTime().value_in_unit(unit.picosecond)
    raise RuntimeError(f"no loadable resume checkpoint (all torn/invalid): {errors}")


def _truncate_xtc(path, keep_frames: int) -> None:
    """Rewrite `path` keeping only its first `keep_frames` frames (also heals a torn final frame)."""
    from mdtraj.formats import XTCTrajectoryFile

    with XTCTrajectoryFile(str(path), "r") as fh:
        xyz, times, steps, box = fh.read(n_frames=keep_frames)
    with XTCTrajectoryFile(str(path), "w") as fh:
        fh.write(xyz, time=times, step=steps, box=box)


def _truncate_log(path, max_step: int) -> None:
    """Keep the StateDataReporter header + the rows whose Step column is <= `max_step`."""
    path = Path(path)
    lines = path.read_text().splitlines(keepends=True)
    if not lines:
        return
    kept = [lines[0]]  # the single '#'-prefixed header row
    for line in lines[1:]:
        if not line.strip():
            continue
        step = int(float(line.split(",")[0].strip().strip('"')))
        if step > max_step:
            break  # rows are sequential
        kept.append(line)
    path.write_text("".join(kept))


def _write_resume_reports(layout, cfg, start_step: int, wall_seconds: float, log_path: Path) -> None:
    ts = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    # Per-segment provenance: this segment may run on a different node/GPU/openmm version than the
    # prepare-time run_info.json records.
    write_json(RunInfo.collect(cfg), layout.record(f"run_info.resume.{ts}.json"))

    seg_steps = cfg.prod_steps - start_step
    seg_ns = seg_steps * cfg.dt_ps / 1000.0
    report = {
        "resumed": True,
        "resumed_from_step": start_step,
        "resumed_from_ns": start_step * cfg.dt_ps / 1000.0,
        "target_step": cfg.prod_steps,
        "target_ns": cfg.total_ns,
        "segment_steps": seg_steps,
        "segment_ns": seg_ns,
        "segment_wall_seconds": wall_seconds,
        "segment_ns_per_day": seg_ns / (wall_seconds / 86_400.0) if wall_seconds > 0 else 0.0,
        "n_frames_total": len(Path(log_path).read_text().splitlines()) - 1,  # rows minus header
        "platform": cfg.platform,
    }
    write_json(report, layout.record(f"production/run_report.resume.{ts}.json"))


def _col(df: pd.DataFrame, needle: str) -> pd.Series:
    """Pick the one StateDataReporter column whose verbose header contains `needle`.

    Fails loud on zero or >1 matches — an ambiguous match (e.g. two "Energy" columns)
    would otherwise silently take the first and report the wrong statistic.
    """
    matches = [c for c in df.columns if needle.lower() in c.lower()]
    assert len(matches) == 1, f"expected exactly one column matching {needle!r}, got {matches}"
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
