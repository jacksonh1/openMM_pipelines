"""Secondary structure over the trajectory via mdtraj's built-in DSSP.

`compute_dssp` is a native mdtraj implementation — no external `mkdssp` needed.
"""

from __future__ import annotations

import mdtraj as md
import numpy as np

# Full (non-simplified) DSSP categories for the per-residue timeline, as
# (mdtraj code, label, colour). mdtraj's native DSSP emits H/G/I/E/B/T/S/' '/NA — it does
# NOT produce DSSP-4's PP-II ('P') or '=' states (those need mkdssp 4). ' ' (space) is coil.
DSSP_CATEGORIES: list[tuple[str, str, str]] = [
    ("H", "α-helix", "#1f4ea1"),     # dark blue
    ("G", "3₁₀-helix", "#5b9bd5"),   # light blue
    ("I", "π-helix", "#7b5fb0"),     # purple
    ("E", "β-strand", "#c0392b"),    # red
    ("B", "β-bridge", "#e59b3a"),    # orange
    ("T", "turn", "#2e8b57"),        # dark green
    ("S", "bend", "#9cd3a0"),        # light green
    (" ", "coil", "#d9d9d9"),        # light grey
    ("NA", "other", "#d4b106"),      # gold (termini / caps mdtraj can't assign)
]
_CODE_TO_INDEX = {code: i for i, (code, _, _) in enumerate(DSSP_CATEGORIES)}
_OTHER_INDEX = _CODE_TO_INDEX["NA"]


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


def dssp_timeline(traj: md.Trajectory) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue × per-frame secondary structure as category INDICES into DSSP_CATEGORIES.

    Uses the full (non-simplified) DSSP codes for the richer helix/strand/loop palette. Returns
    ``(codes, resseq)`` where ``codes`` is ``(n_frames, n_residues)`` int8 and ``resseq`` the
    per-residue numbers. An unrecognised code maps to the 'other' category.
    """
    raw = md.compute_dssp(traj, simplified=False)  # (n_frames, n_residues) of DSSP strings
    codes = np.vectorize(lambda c: _CODE_TO_INDEX.get(c, _OTHER_INDEX))(raw).astype(np.int8)
    resseq = np.array([res.resSeq for res in traj.topology.residues], dtype=int)
    assert codes.shape[1] == len(resseq), "DSSP residue count != topology residue count"
    return codes, resseq
