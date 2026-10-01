"""Secondary structure over the trajectory via mdtraj's built-in DSSP.

`compute_dssp` is a native mdtraj implementation — no external `mkdssp` needed.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np


def secondary_structure(traj: md.Trajectory, simplified: bool = True) -> np.ndarray:
    """(n_frames, n_residues) array of DSSP codes. Simplified: 'H'/'E'/'C' (+'NA')."""
    return md.compute_dssp(traj, simplified=simplified)


def ss_fractions(traj: md.Trajectory) -> dict[str, np.ndarray]:
    """Per-frame fraction of residues in helix / sheet / coil (simplified DSSP)."""
    ss = secondary_structure(traj, simplified=True)
    return {
        "helix": (ss == "H").mean(axis=1),
        "sheet": (ss == "E").mean(axis=1),
        "coil": (ss == "C").mean(axis=1),
    }
