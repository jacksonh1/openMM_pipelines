"""End-to-end REST2 on a small real system (CPU by default).

Two things proved here:
1. **The solute is actually scaled** on the REAL prepared System — protein atoms get
   √λ·q / λ·ε and their intramolecular bonded energy scales by λ, while water is untouched.
   This is the correctness crux of REST2 (scaling the wrong region silently gives a wrong
   result), so it is asserted directly on the force-field parameters, not just inferred.
2. The prepare/produce seam reused, Hamiltonian replica exchange run at one physical
   temperature, and the de-multiplex of the λ=1 reference ensemble.

Marked slow (real dynamics). Tiny ladder + lengths — wiring + scaling, not sampling quality.
Set OPENMM_TEST_PLATFORM=CUDA on a GPU node.
"""

import json
import os
from pathlib import Path

import pytest
from openmm import NonbondedForce, unit

from openmm_pipelines import prepare
from openmm_pipelines import rest2
from openmm_pipelines.analysis.rest2 import analyze_rest2
from openmm_pipelines.lib.selections import protein_atoms
from openmm_pipelines.lib.system_edits import strip_barostats
from openmm_pipelines.rest2.config import REST2Config
from openmm_pipelines.rest2.scaling import scale_solute

pytestmark = pytest.mark.slow

STRUCT = Path(__file__).parents[1] / "examples" / "input_structures" / "helix_fusion.pdb"
PLATFORM = os.environ.get("OPENMM_TEST_PLATFORM", "CPU")
LAM = 0.49  # sqrt = 0.7


@pytest.fixture(scope="module")
def prepared():
    import tempfile
    outdir = Path(tempfile.mkdtemp(prefix="rest2_prep_"))
    cfg = REST2Config(
        pdb_in=STRUCT, outdir=outdir, platform=PLATFORM, padding_nm=0.6,
        nvt_equil_ns=0.002, density_seg_steps=100, density_min_seg=2,
        density_max_seg=4, density_tol_rel=0.5,
        n_replicas=3, max_effective_temperature_k=450.0,
        total_ns=0.001, exchange_attempt_ps=0.1, equilibration_ns=0.0002,
        checkpoint_ps=0.1,
    )
    return prepare(cfg), cfg


@pytest.fixture(scope="module")
def rest2_run(prepared):
    eq, cfg = prepared
    rest2.run(eq, cfg)
    return cfg


def _nonbonded(system) -> NonbondedForce:
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    nb = [f for f in forces if isinstance(f, NonbondedForce)]
    assert len(nb) == 1, f"expected exactly one NonbondedForce, got {len(nb)}"
    return nb[0]


def test_solute_is_actually_scaled_on_real_system(prepared):
    eq, _ = prepared
    base = strip_barostats(eq.system)
    solute = set(protein_atoms(eq.topology))
    assert solute, "no solute atoms selected"

    scaled = scale_solute(base, solute, LAM)
    base_nb, scaled_nb = _nonbonded(base), _nonbonded(scaled)
    sqrt_lam = LAM**0.5

    # A protein atom: charge -> √λ·q, ε -> λ·ε, σ unchanged.
    p = next(iter(solute))
    bq, bs, be = base_nb.getParticleParameters(p)
    sq, ss, se = scaled_nb.getParticleParameters(p)
    assert sq.value_in_unit(unit.elementary_charge) == pytest.approx(
        bq.value_in_unit(unit.elementary_charge) * sqrt_lam, abs=1e-9)
    assert se.value_in_unit(unit.kilojoule_per_mole) == pytest.approx(
        be.value_in_unit(unit.kilojoule_per_mole) * LAM, rel=1e-6, abs=1e-9)
    assert ss.value_in_unit(unit.nanometer) == pytest.approx(bs.value_in_unit(unit.nanometer))

    # A water/solvent atom (not in the solute): completely untouched.
    w = next(i for i in range(base_nb.getNumParticles()) if i not in solute)
    for getter in ("charge", "sigma", "epsilon"):
        idx = {"charge": 0, "sigma": 1, "epsilon": 2}[getter]
        bval = base_nb.getParticleParameters(w)[idx]
        sval = scaled_nb.getParticleParameters(w)[idx]
        assert str(sval) == str(bval), f"solvent {getter} changed: {bval} -> {sval}"


def test_run_writes_contract(rest2_run):
    cfg = rest2_run
    rdir = cfg.outdir / "rest2"
    assert (rdir / "rest2.nc").exists()

    ladder = json.loads((rdir / "ladder.json").read_text())
    assert len(ladder["lambdas"]) == cfg.n_replicas
    assert ladder["lambdas"][0] == pytest.approx(1.0)
    assert ladder["effective_temperatures_k"][0] == pytest.approx(cfg.temperature_k)
    assert ladder["effective_temperatures_k"][-1] == pytest.approx(
        cfg.max_effective_temperature_k)
    assert ladder["n_solute_atoms"] > 0

    report = json.loads((rdir / "run_report.json").read_text())
    assert report["n_replicas"] == 3
    assert report["temperature_k"] == pytest.approx(cfg.temperature_k)
    assert report["n_iterations"] == cfg.n_iterations
    assert 0.0 <= report["mean_neighbor_acceptance"] <= 1.0


def test_analyze_demuxes_reference_lambda(rest2_run):
    cfg = rest2_run
    adir = analyze_rest2(cfg.outdir, reference_state=0)

    exchange = json.loads((adir / "exchange_report.json").read_text())
    assert exchange["n_states"] == 3
    assert len(exchange["neighbor_acceptance"]) == 2

    summary = json.loads((adir / "rest2_analysis_report.json").read_text())
    assert summary["reference_state"] == 0
    assert summary["n_frames_demuxed"] >= 1
    assert summary["rmsd_to_design_ang"]["mean"] >= 0.0
    assert summary["n_clusters"] >= 1  # λ=1 ensemble clustered (Cα-RMSD)

    assert (adir / "state00_ensemble.xtc").exists()     # PBC-imaged, desolvated, aligned, centered
    assert (adir / "state00_topology.pdb").exists()
    assert (adir / "acceptance.png").exists()
    assert (adir / "replica_walk.png").exists()
