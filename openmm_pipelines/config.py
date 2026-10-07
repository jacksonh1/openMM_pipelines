"""PrepConfig — everything through equilibration. Shared by MD, REMD and REST2.

Pure Python — no OpenMM import, so config construction and validation run anywhere
(no GPU, no openmm env). Engine configs (`MDConfig`, future `REMDConfig`/
`REST2Config`) subclass `PrepConfig`; inheritance (not nesting) keeps attribute
access flat while `prepare()` depends only on `PrepConfig` fields.
"""

from __future__ import annotations

import warnings
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, model_validator

from .forcefield import ProteinFF, WaterModel, _RECOMMENDED_WATER, resolve


class BoxShape(str, Enum):
    CUBE = "cube"
    DODECAHEDRON = "dodecahedron"
    OCTAHEDRON = "octahedron"


class PrepConfig(BaseModel, frozen=True):
    """Everything through equilibration — shared by MD, REMD and REST2."""

    # --- system ---
    pdb_in: Path
    outdir: Path
    protein_ff: ProteinFF = ProteinFF.AMBER14SB
    water: WaterModel = WaterModel.TIP3P
    ph: float = 7.4  # physiological; PDBFixer protonation heuristic
    box_shape: BoxShape = BoxShape.DODECAHEDRON
    padding_nm: float = 1.0  # minimum solute-to-box-edge distance
    ionic_strength_molar: float = 0.15  # NaCl added beyond charge neutralization
    neutralize: bool = True

    # --- dynamics ---
    temperature_k: float = 300.0  # reference/equilibration temperature (= T_min for REMD)
    dt_ps: float = 0.004  # 4 fs -> requires HMR (default on; see below)
    hydrogen_mass_amu: float = 4.0  # HMR; set 1.008 + dt_ps=0.002 to disable
    friction_ps: float = 1.0  # Langevin friction coefficient
    ref_p_bar: float = 1.0
    restraint_k: float = 1000.0  # kJ/mol/nm^2 on protein heavy atoms

    # --- equilibration lengths ---
    nvt_equil_ns: float = 0.2  # NVT equilibration length

    # --- NPT density equilibration (adaptive plateau protocol) ---
    density_seg_steps: int = 10_000
    density_min_seg: int = 8
    density_max_seg: int = 20
    density_tol_rel: float = 0.005

    # --- runtime ---
    platform: str = "CUDA"
    seed: int | None = None  # None = fresh RNG per run; int pins setup + RNG streams

    @model_validator(mode="after")
    def _check(self) -> PrepConfig:
        assert self.pdb_in.exists(), f"pdb_in not found: {self.pdb_in}"
        assert self.density_min_seg <= self.density_max_seg, (
            f"density_min_seg ({self.density_min_seg}) > density_max_seg "
            f"({self.density_max_seg})"
        )
        assert self.density_min_seg >= 2, "density plateau fit needs density_min_seg >= 2"
        resolve(self.protein_ff, self.water)  # raises ForceFieldError on illegal pairing
        if self.water is not _RECOMMENDED_WATER[self.protein_ff]:
            warnings.warn(
                f"{self.protein_ff.name} is recommended with "
                f"{_RECOMMENDED_WATER[self.protein_ff].name}; you chose {self.water.name}. "
                "Proceeding — verify this is intentional (e.g. ff19SB+TIP3P is a "
                "known-degraded pairing).",
                stacklevel=2,
            )
        if self.dt_ps > 0.0025:  # HMR sanity: 4 fs without repartitioned H is unstable
            assert self.hydrogen_mass_amu >= 3.0, (
                f"dt_ps={self.dt_ps} needs HMR (hydrogen_mass_amu >= 3); "
                f"got {self.hydrogen_mass_amu}"
            )
        return self

    # Derived equilibration step counts.
    @property
    def nvt_equil_steps(self) -> int:
        return round(self.nvt_equil_ns * 1000 / self.dt_ps)


class MultiStateProductionConfig(PrepConfig):
    """Shared production fields for the openmmtools multistate backends (REMD / REST2).

    openmmtools runs in **iterations**: one iteration = `steps_per_iteration` MD steps per
    replica then one exchange-attempt sweep. Production length is a per-replica MD-time wall
    (`total_ns`) plus the MD time between exchange attempts (`exchange_attempt_ps`); the
    iteration counts derive from those. Ladder specifics live in the engine subclasses.
    """

    # --- production (per replica) ---
    total_ns: float = 20.0  # MD time collected per replica
    exchange_attempt_ps: float = 1.0  # MD time per iteration (between exchange attempts)
    equilibration_ns: float = 0.5  # per-replica equilibration, discarded

    # --- replica exchange ---
    replica_mixing_scheme: str = "swap-all"  # openmmtools: 'swap-all' | 'swap-neighbors'
    checkpoint_ps: float = 100.0  # full-state checkpoint interval (restart granularity)

    @model_validator(mode="after")
    def _check_multistate(self) -> "MultiStateProductionConfig":
        assert self.exchange_attempt_ps > 0, "exchange_attempt_ps must be > 0"
        return self

    @property
    def steps_per_iteration(self) -> int:
        return round(self.exchange_attempt_ps / self.dt_ps)

    @property
    def n_iterations(self) -> int:
        return round(self.total_ns * 1000 / self.exchange_attempt_ps)

    @property
    def equilibration_iterations(self) -> int:
        return round(self.equilibration_ns * 1000 / self.exchange_attempt_ps)

    @property
    def checkpoint_interval_iterations(self) -> int:
        return max(round(self.checkpoint_ps / self.exchange_attempt_ps), 1)
