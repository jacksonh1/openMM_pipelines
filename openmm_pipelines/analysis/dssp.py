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
    """Per-frame fraction of residues in helix / sheet / coil / NA (simplified DSSP).

    `NA` is mdtraj's code for residues it cannot assign — chain termini and caps
    (ACE/NME, common in these designs). It is reported as its own fraction rather
    than folded into coil, so the four fractions sum to 1 and coil is not inflated.
    """
    ss = secondary_structure(traj, simplified=True)
    codes = set(np.unique(ss).tolist())
    assert codes <= {"H", "E", "C", "NA"}, f"unexpected DSSP codes: {codes - {'H', 'E', 'C', 'NA'}}"
    return {
        "helix": (ss == "H").mean(axis=1),
        "sheet": (ss == "E").mean(axis=1),
        "coil": (ss == "C").mean(axis=1),
        "na": (ss == "NA").mean(axis=1),
    }
