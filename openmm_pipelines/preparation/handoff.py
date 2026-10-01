"""The handoff dataclasses across the prepare/produce seam.

`BuiltSystem` is the output of `build.py` (pre-equilibration). `EquilibratedSystem`
is the output of `prepare()` — plain OpenMM objects that BOTH `app.Simulation` and
`openmmtools.multistate` can consume, so equilibration runs once and any production
backend (MD / REMD / REST2, same or separate job) starts from it.

`to_disk` / `from_disk` use the PORTABLE serialization (System XML + a saveState-style
State XML + an mmCIF topology), so `from_disk(prepared_dir)` reloads on any node with
no extra arguments. mmCIF (not PDB) because solvated boxes routinely exceed PDB's
99,999-atom serial limit. See the output contract in knowledgebase/DECISIONS.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import openmm
from openmm import System, XmlSerializer, unit
from openmm.app import PDBxFile, Topology

_SYSTEM_XML = "system.xml"
_STATE_XML = "state.xml"
_TOPOLOGY_CIF = "topology.cif"


@dataclass(frozen=True)
class BuiltSystem:
    """Solvated, parameterized system before equilibration."""

    system: System
    topology: Topology
    positions: object  # list[Vec3] with units


@dataclass(frozen=True)
class EquilibratedSystem:
    """Output of prepare(): plain OpenMM objects for any production backend."""

    system: System
    topology: Topology
    positions: object  # list[Vec3] with units
    box_vectors: object  # 3x3 with units

    def to_disk(self, prepared_dir) -> Path:
        prepared_dir = Path(prepared_dir)
        prepared_dir.mkdir(parents=True, exist_ok=True)

        (prepared_dir / _SYSTEM_XML).write_text(XmlSerializer.serialize(self.system))

        # Serialize a portable State (positions + box) via a throwaway Reference
        # context — no dynamics, no GPU allocation.
        context = openmm.Context(
            self.system,
            openmm.VerletIntegrator(1.0 * unit.femtosecond),
            openmm.Platform.getPlatformByName("Reference"),
        )
        context.setPeriodicBoxVectors(*self.box_vectors)
        context.setPositions(self.positions)
        state = context.getState(getPositions=True)
        (prepared_dir / _STATE_XML).write_text(XmlSerializer.serialize(state))
        del context

        with open(prepared_dir / _TOPOLOGY_CIF, "w") as fh:
            PDBxFile.writeFile(self.topology, self.positions, fh)
        return prepared_dir

    @classmethod
    def from_disk(cls, prepared_dir) -> "EquilibratedSystem":
        prepared_dir = Path(prepared_dir)
        system = XmlSerializer.deserialize((prepared_dir / _SYSTEM_XML).read_text())
        state = XmlSerializer.deserialize((prepared_dir / _STATE_XML).read_text())
        topology = PDBxFile(str(prepared_dir / _TOPOLOGY_CIF)).topology
        return cls(
            system=system,
            topology=topology,
            positions=state.getPositions(),
            box_vectors=state.getPeriodicBoxVectors(),
        )
