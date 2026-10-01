"""Equilibration primitives against a small real solvated system (CPU platform).

Marked slow (runs real dynamics). The density plateau *math* is unit-tested in
test_density; here we check the stage wiring: minimize lowers energy and stays
finite, NVT advances steps, the barostat attaches, and equilibrate_density returns a
segment count when the tolerance is met and raises when it is not.
"""

import os
from pathlib import Path

import openmm
import pytest
from openmm import LangevinMiddleIntegrator, XmlSerializer, unit
from openmm.app import Simulation

from openmm_pipelines.preparation import equilibration_steps as equilibration
from openmm_pipelines.preparation.build import build_system
from openmm_pipelines.preparation.minimize import minimize
from openmm_pipelines.md.config import MDConfig

pytestmark = pytest.mark.slow

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"

# CPU by default (works anywhere); set OPENMM_TEST_PLATFORM=CUDA on a GPU node.
PLATFORM = os.environ.get("OPENMM_TEST_PLATFORM", "CPU")


@pytest.fixture(scope="module")
def built():
    cfg = MDConfig(pdb_in=STRUCT, outdir="/tmp/_stages_unused", padding_nm=0.8,
                   platform=PLATFORM)
    b = build_system(cfg)
    # serialize once so each test can deserialize a FRESH system (add_barostat mutates it)
    return XmlSerializer.serialize(b.system), b.topology, b.positions, cfg


def _fresh_sim(built):
    system_xml, topology, positions, cfg = built
    system = XmlSerializer.deserialize(system_xml)
    integ = LangevinMiddleIntegrator(cfg.temperature_k * unit.kelvin,
                                     cfg.friction_ps / unit.picosecond,
                                     cfg.dt_ps * unit.picosecond)
    sim = Simulation(topology, system, integ,
                     openmm.Platform.getPlatformByName(PLATFORM))
    sim.context.setPositions(positions)
    return sim, cfg


def _potential(sim):
    return sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    )


def test_minimize_lowers_finite_energy(built):
    sim, cfg = _fresh_sim(built)
    before = _potential(sim)
    minimize(sim, cfg)
    after = _potential(sim)
    assert after <= before
    assert after == after  # finite (NaN != NaN)


def test_nvt_advances_steps(built):
    sim, cfg = _fresh_sim(built)
    minimize(sim, cfg)
    cfg = cfg.model_copy(update={"nvt_equil_ns": 0.002})  # 500 steps @ 4 fs
    start = sim.context.getStepCount()
    equilibration.equilibrate_nvt(sim, cfg)
    assert sim.context.getStepCount() - start == cfg.nvt_equil_steps == 500


def test_add_barostat_attaches_force(built):
    sim, cfg = _fresh_sim(built)
    minimize(sim, cfg)
    n_before = sim.system.getNumForces()
    equilibration.add_barostat(sim, cfg)
    assert sim.system.getNumForces() == n_before + 1
    assert any(isinstance(sim.system.getForce(i), openmm.MonteCarloBarostat)
               for i in range(sim.system.getNumForces()))


def test_density_converges_with_loose_tol(built):
    sim, cfg = _fresh_sim(built)
    minimize(sim, cfg)
    cfg = cfg.model_copy(update={
        "density_seg_steps": 50, "density_min_seg": 2, "density_max_seg": 6,
        "density_tol_rel": 1.0,  # anything plateaus under a 100%-drift tolerance
    })
    equilibration.add_barostat(sim, cfg)
    volumes = equilibration.equilibrate_density(sim, cfg)
    assert len(volumes) == cfg.density_min_seg == 2  # declared at the earliest allowed segment


def test_density_raises_when_never_converges(built):
    sim, cfg = _fresh_sim(built)
    minimize(sim, cfg)
    cfg = cfg.model_copy(update={
        "density_seg_steps": 50, "density_min_seg": 2, "density_max_seg": 3,
        "density_tol_rel": 1e-12,  # impossible tolerance -> never plateaus
    })
    equilibration.add_barostat(sim, cfg)
    with pytest.raises(equilibration.DensityNotConverged):
        equilibration.equilibrate_density(sim, cfg)
