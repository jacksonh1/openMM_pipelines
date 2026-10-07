"""Typed metadata + QC outputs written to an `outdir` (the output contract, draft).

Three records, all dataclasses (not raw dicts, per the fail-loud style):
- `RunInfo`    — software/env/platform provenance (reproducibility).
- `PrepReport` — QC from the prepare pipeline (build -> minimize -> equilibration).
- `RunReport`  — QC from a production run (timing, stability, energy/T/density stats).

See `knowledgebase/DECISIONS.md` (output contract, first draft).
"""

from __future__ import annotations

import dataclasses
import datetime as _dt
import json
import os
import platform as _platform
import socket
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def write_json(obj, path: Path) -> None:
    """Serialize a dataclass (or dict) to pretty JSON."""
    path = Path(path)
    data = dataclasses.asdict(obj) if dataclasses.is_dataclass(obj) else obj
    path.write_text(json.dumps(data, indent=2, default=str))


def _pkg_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


@dataclass(frozen=True)
class RunInfo:
    """What the config alone can't reconstruct — the environment that produced a run.

    GPU runs are NOT bit-reproducible even with a fixed seed (nondeterministic
    reduction order); `seed` pins setup + RNG streams, not the trajectory. A
    state-resumed trajectory also differs (saveState omits RNG state).
    """

    timestamp_utc: str
    openmm_pipelines_version: str | None
    openmm_version: str
    python_version: str
    mdtraj_version: str | None
    numpy_version: str | None
    pdbfixer_version: str | None
    platform: str  # requested compute platform (CUDA/CPU/...)
    precision: str
    hostname: str
    slurm_job_id: str | None
    forcefield_xmls: list[str]
    water_model: str
    packing_box: str

    @classmethod
    def collect(cls, cfg) -> "RunInfo":
        import openmm

        from .forcefield import resolve

        xmls, packing = resolve(cfg.protein_ff, cfg.water)
        precision = "mixed" if cfg.platform in ("CUDA", "OpenCL", "HIP") else "double"
        return cls(
            timestamp_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(),
            openmm_pipelines_version=_pkg_version("openmm_pipelines"),
            openmm_version=openmm.version.version,
            python_version=_platform.python_version(),
            mdtraj_version=_pkg_version("mdtraj"),
            numpy_version=_pkg_version("numpy"),
            pdbfixer_version=_pkg_version("pdbfixer"),
            platform=cfg.platform,
            precision=precision,
            hostname=socket.gethostname(),
            slurm_job_id=os.environ.get("SLURM_JOB_ID"),
            forcefield_xmls=xmls,
            water_model=cfg.water.name,
            packing_box=packing,
        )


@dataclass(frozen=True)
class PrepReport:
    """QC from the prepare pipeline."""

    n_atoms: int
    n_protein_heavy: int
    n_waters: int
    n_ions: int
    box_vectors_nm: list[list[float]]
    box_volume_nm3: float
    disulfides: int
    pe_initial_kj_mol: float
    pe_minimized_kj_mol: float
    nvt_final_temperature_k: float
    density_segment_volumes_nm3: list[float]
    density_segments: int
    density_converged: bool
    density_drift_rel: float
    prep_wall_seconds: float


@dataclass(frozen=True)
class RunReport:
    """QC from a production run."""

    n_steps: int
    total_ns: float
    n_frames: int
    wall_seconds: float
    ns_per_day: float
    temperature_mean_k: float
    temperature_std_k: float
    density_mean: float
    density_std: float
    potential_energy_mean_kj_mol: float
    potential_energy_min_kj_mol: float
    potential_energy_max_kj_mol: float
    total_energy_mean_kj_mol: float
    total_energy_std_kj_mol: float
    stable: bool
    platform: str
    notes: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class REMDRunReport:
    """QC from a T-REMD production run.

    The health check for replica exchange is NOT single-trajectory stability but
    *mixing*: `mean_neighbor_acceptance` is the average exchange-acceptance over
    adjacent ladder rungs; a near-zero value means a broken ladder link where replicas
    stop diffusing through temperature space (add rungs there). Round-trip diffusion is
    reported separately by `analysis/remd.py`.
    """

    n_replicas: int
    temperatures_k: list[float]
    n_iterations: int
    steps_per_iteration: int
    exchange_attempt_ps: float
    total_ns_per_replica: float
    replica_mixing_scheme: str
    mean_neighbor_acceptance: float
    wall_seconds: float
    aggregate_ns_per_day: float
    platform: str
    notes: list[str] = field(default_factory=list)


@dataclass(frozen=True)
class REST2RunReport:
    """QC from a REST2 production run.

    Like REMD, the health check is *mixing* (`mean_neighbor_acceptance` over adjacent λ
    rungs), not single-trajectory stability. REST2-specific: every replica is at the same
    physical `temperature_k`; the ladder is reported as both the solute-scaling factors
    `lambdas` (λ[0]=1 = the reference, unscaled ensemble) and the `effective_temperatures_k`
    they correspond to. `n_solute_atoms` records how much of the system was tempered.
    """

    n_replicas: int
    temperature_k: float  # physical temperature all replicas run at
    lambdas: list[float]
    effective_temperatures_k: list[float]
    n_solute_atoms: int
    n_iterations: int
    steps_per_iteration: int
    exchange_attempt_ps: float
    total_ns_per_replica: float
    replica_mixing_scheme: str
    mean_neighbor_acceptance: float
    wall_seconds: float
    aggregate_ns_per_day: float
    platform: str
    notes: list[str] = field(default_factory=list)
