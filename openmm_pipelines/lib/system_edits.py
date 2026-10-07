"""System-level edits shared by the production backends.

Pure OpenMM (no openmmtools). `clone_system` is the serializer round-trip used to get an
independent copy before mutating — so a backend never scribbles on the shared
`eq.system`. `strip_barostats` drops every MonteCarloBarostat for NVT production (the
equilibrated box is already at density; replicas hold volume fixed).
"""

from __future__ import annotations

from openmm import (
    MonteCarloAnisotropicBarostat,
    MonteCarloBarostat,
    MonteCarloMembraneBarostat,
    System,
    XmlSerializer,
)

# Every barostat flavour OpenMM ships — NVT production must carry none of them.
_BAROSTATS = (MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat)


def clone_system(system: System) -> System:
    """An independent deep copy of `system` (via the serializer), safe to mutate."""
    return XmlSerializer.deserialize(XmlSerializer.serialize(system))


def strip_barostats(system: System) -> System:
    """A copy of `system` with every barostat removed (NVT production system).

    Round-trips through the serializer so `system` itself is untouched, then deletes
    barostats from the back so force indices stay valid during removal.
    """
    copy = clone_system(system)
    for i in reversed(range(copy.getNumForces())):
        if isinstance(copy.getForce(i), _BAROSTATS):
            copy.removeForce(i)
    assert not any(
        isinstance(copy.getForce(i), _BAROSTATS) for i in range(copy.getNumForces())
    ), "barostat still present after strip — NVT production needs it gone"
    return copy
