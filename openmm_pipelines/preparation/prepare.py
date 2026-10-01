"""prepare(cfg) -> EquilibratedSystem — the shared pipeline composer.

Runs build -> restrain -> minimize -> reanchor -> NVT -> barostat -> NPT density ->
release, on a throwaway equilibration `Simulation`, then hands back plain OpenMM
objects and writes the output contract (config / run_info / input copy / prepared/ +
prep_report). Takes a `PrepConfig` (or any subclass), so it is the same call for MD,
REMD and REST2. The order of the pipeline lives here and nowhere else.
"""

from __future__ import annotations

import shutil
import time
from pathlib import Path

from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import Simulation

from ..lib.density_convergence import assess_plateau
from ..lib.restraints import add_restraint, reanchor
from ..lib.selections import protein_heavy_atoms
from ..reports import PrepReport, RunInfo, write_json
from . import equilibration_steps as eq
from .build import _count_disulfides, build_system
from .handoff import EquilibratedSystem
from .minimize import minimize

_GPU_PLATFORMS = ("CUDA", "OpenCL", "HIP")
_ION_RESIDUES = frozenset(
    {"NA", "CL", "K", "MG", "CA", "ZN", "BR", "CS", "RB", "LI", "F", "IOD", "SOD", "CLA"}
)


def _make_simulation(topology, system, cfg):
    integ = LangevinMiddleIntegrator(
        cfg.temperature_k * unit.kelvin,
        cfg.friction_ps / unit.picosecond,
        cfg.dt_ps * unit.picosecond,
    )
    if cfg.seed is not None:
        integ.setRandomNumberSeed(cfg.seed)
    platform = Platform.getPlatformByName(cfg.platform)
    props = {"Precision": "mixed"} if cfg.platform in _GPU_PLATFORMS else {}
    return Simulation(topology, system, integ, platform, props)


def _potential_kj(sim) -> float:
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )


def _instantaneous_temperature_k(sim) -> float:
    """Temperature from kinetic energy and the system's degrees of freedom."""
    system = sim.system
    dof = sum(
        3 for i in range(system.getNumParticles()) if system.getParticleMass(i) > 0 * unit.amu
    )
    dof -= system.getNumConstraints()
    if any(system.getForce(i).__class__.__name__ == "CMMotionRemover"
           for i in range(system.getNumForces())):
        dof -= 3
    ke = sim.context.getState(getEnergy=True).getKineticEnergy()
    return (2 * ke / (dof * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)


def prepare(cfg) -> EquilibratedSystem:
    t0 = time.perf_counter()
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    built = build_system(cfg)
    heavy = protein_heavy_atoms(built.topology)
    restraint = add_restraint(built.system, built.positions, heavy, cfg.restraint_k)

    sim = _make_simulation(built.topology, built.system, cfg)
    sim.context.setPositions(built.positions)

    pe_initial = _potential_kj(sim)
    minimize(sim, cfg)
    pe_minimized = _potential_kj(sim)

    minimized_positions = sim.context.getState(getPositions=True).getPositions()
    reanchor(restraint, sim.context, minimized_positions, heavy)  # restrain to minimized pose

    eq.equilibrate_nvt(sim, cfg)
    nvt_temperature = _instantaneous_temperature_k(sim)

    eq.add_barostat(sim, cfg)
    segment_volumes = eq.equilibrate_density(sim, cfg)

    # Release restraints. Bake k=0 into the FORCE DEFAULT, not just this context: a new
    # Context (production, or an eq reloaded from disk) initializes globals to the
    # force's default, which would otherwise re-enable restraints at full strength.
    restraint.setGlobalParameterDefaultValue(0, 0.0)
    sim.context.setParameter("k", 0.0)

    state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    equilibrated = EquilibratedSystem(
        system=built.system,
        topology=built.topology,
        positions=state.getPositions(),
        box_vectors=state.getPeriodicBoxVectors(),
    )

    _write_outputs(cfg, built, equilibrated, heavy, state,
                   pe_initial, pe_minimized, nvt_temperature, segment_volumes,
                   time.perf_counter() - t0)
    return equilibrated


def _write_outputs(cfg, built, equilibrated, heavy, state, pe_initial, pe_minimized,
                   nvt_temperature, segment_volumes, wall_seconds) -> None:
    write_json(cfg.model_dump(mode="json"), cfg.outdir / "config.json")
    write_json(RunInfo.collect(cfg), cfg.outdir / "run_info.json")
    shutil.copyfile(cfg.pdb_in, cfg.outdir / "input.pdb")

    equilibrated.to_disk(cfg.outdir / "prepared")

    topology = built.topology
    box = state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
    volume = state.getPeriodicBoxVolume().value_in_unit(unit.nanometer**3)
    n_waters = sum(1 for r in topology.residues() if r.name == "HOH")
    n_ions = sum(1 for r in topology.residues() if r.name in _ION_RESIDUES)
    plateau = assess_plateau(segment_volumes, cfg.density_tol_rel, cfg.density_min_seg)

    report = PrepReport(
        n_atoms=built.system.getNumParticles(),
        n_protein_heavy=len(heavy),
        n_waters=n_waters,
        n_ions=n_ions,
        box_vectors_nm=[[v[0], v[1], v[2]] for v in box],
        box_volume_nm3=volume,
        disulfides=_count_disulfides(topology),
        pe_initial_kj_mol=pe_initial,
        pe_minimized_kj_mol=pe_minimized,
        nvt_final_temperature_k=nvt_temperature,
        density_segment_volumes_nm3=segment_volumes,
        density_segments=len(segment_volumes),
        density_converged=plateau.converged,
        density_drift_rel=plateau.drift_rel,
        prep_wall_seconds=wall_seconds,
    )
    write_json(report, cfg.outdir / "prepared" / "prep_report.json")
