"""The equilibration steps: restrained NVT, then NPT density equilibration.

These are the dynamics stages `prepare()` runs after minimization to bring the built,
restrained system to a production-ready equilibrated box. Engine-agnostic (every
backend shares them). `add_barostat` is the NVT->NPT transition that starts density
equilibration — glue for that step, so it lives here rather than in `lib/`.

One `LangevinMiddleIntegrator` (created by the caller) is both integrator and
thermostat for the whole system — no per-group thermostats. The `MonteCarloBarostat`
uses its default volume-move frequency, so it samples the correct NPT distribution
without a coupling-time knob. Each step acts on an `app.Simulation` in place.
"""

from __future__ import annotations

from openmm import MonteCarloBarostat, unit

from ..lib.density_convergence import density_converged


class DensityNotConverged(RuntimeError):
    pass


def equilibrate_nvt(sim, cfg) -> None:
    seed = cfg.seed if cfg.seed is not None else 0
    sim.context.setVelocitiesToTemperature(cfg.temperature_k * unit.kelvin, seed)
    sim.step(cfg.nvt_equil_steps)


def add_barostat(sim, cfg) -> None:
    sim.system.addForce(
        MonteCarloBarostat(cfg.ref_p_bar * unit.bar, cfg.temperature_k * unit.kelvin)
    )
    sim.context.reinitialize(preserveState=True)  # required after addForce


def equilibrate_density(sim, cfg) -> list[float]:
    """Run NPT in segments until the box-volume plateaus; return the per-segment mean
    volumes (nm^3, chronological). The segment count is ``len(volumes)``, and the
    volumes feed the prep QC report.

    Fail loud (`DensityNotConverged`) if no plateau is reached in `density_max_seg`
    segments — never silently ship an unconverged box.
    """
    volumes: list[float] = []
    for seg in range(1, cfg.density_max_seg + 1):
        sim.step(cfg.density_seg_steps)
        v = sim.context.getState().getPeriodicBoxVolume().value_in_unit(
            unit.nanometer**3
        )
        volumes.append(v)
        if seg >= cfg.density_min_seg and density_converged(
            volumes, cfg.density_tol_rel, cfg.density_min_seg
        ):
            return volumes
    raise DensityNotConverged(
        f"no volume plateau in {cfg.density_max_seg} segments: {volumes}"
    )
