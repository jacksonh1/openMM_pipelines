"""Drift-from-design report — the shared core of every analysis.

Measured against the design pose: RMSD / Rg / per-residue RMSF / secondary-structure fractions,
plus three alignment-free fold-retention metrics — native-contact fraction Q, distance-RMSD, and
lDDT (one curve per configured atom set, global + per-residue) — and their figures. This is the
common block between plain MD (`analyze`) and the de-multiplexed REMD ensemble (`analyze_remd`);
each composer wraps it with its own extras (clustering, exchange diagnostics) and its own x-axis
(production time vs demuxed-frame index).

No OpenMM / openmmtools import — mdtraj + biotite (lDDT), so it is shared freely across engines.
"""

from __future__ import annotations

from pathlib import Path

from . import plots
from .config import AnalysisConfig
from .dssp import dssp_timeline, ss_fractions
from .metrics import (
    LDDT_SPEC,
    ca_drmsd,
    lddt_to_reference,
    native_contact_fraction,
    radius_of_gyration,
    rmsd_to_reference,
    rmsf_per_residue,
)


def drift_report(
    traj, reference, native_reference, adir, x, xlabel, rmsd_title, config: AnalysisConfig
) -> tuple[dict, dict]:
    """Compute the drift / fold-retention metrics vs the design, write their PNGs into `adir`, and
    return ``(summary, arrays)``.

    Two references: `reference` is the raw `input.pdb` pose (Cα-mapped RMSD); `native_reference` is
    the topology-matched prepared design (`load_prepared_reference`) that the alignment-free metrics
    (native-contact Q, dRMSD, lDDT) need for 1:1 atom correspondence. `x` is the per-frame x-axis
    (production time in ps, or demuxed-frame index) and `xlabel` its label. `summary` is the
    JSON-ready stats block; `arrays` is the raw per-frame data the caller folds into its own ``.npz``.
    """
    adir = Path(adir)
    rmsd = rmsd_to_reference(traj, reference)
    rg = radius_of_gyration(traj)
    resseq, rmsf = rmsf_per_residue(traj)
    fractions = ss_fractions(traj)

    plots.plot_timeseries(x, rmsd, "Cα RMSD to design (Å)", adir / "rmsd.png",
                          rmsd_title, xlabel=xlabel)
    plots.plot_timeseries(x, rg, "Rg (Å)", adir / "rg.png",
                          "radius of gyration", xlabel=xlabel)
    plots.plot_rmsf(resseq, rmsf, adir / "rmsf.png")
    plots.plot_ss_fractions(x, fractions, adir / "dssp.png", xlabel=xlabel)
    dssp_codes, _ = dssp_timeline(traj)  # per-residue × per-frame SS map
    plots.plot_dssp_timeline(dssp_codes, adir / "dssp_timeline.png")

    summary = {
        "rmsd_to_design_ang": {
            "mean": float(rmsd.mean()), "max": float(rmsd.max()), "final": float(rmsd[-1]),
        },
        "rg_ang": {"mean": float(rg.mean()), "std": float(rg.std())},
        "rmsf_ang": {
            "mean": float(rmsf.mean()), "max": float(rmsf.max()),
            "max_residue": int(resseq[int(rmsf.argmax())]),
        },
        "ss_fraction_mean": {name: float(frac.mean()) for name, frac in fractions.items()},
    }
    arrays = {
        "rmsd_ang": rmsd, "rg_ang": rg, "resseq": resseq, "rmsf_ang": rmsf,
        "dssp_codes": dssp_codes,
        **{f"ss_{name}": frac for name, frac in fractions.items()},
    }

    # Alignment-free fold-retention metrics vs the (topology-matched) design pose.
    q = native_contact_fraction(
        traj, native_reference, cutoff_ang=config.q_native_cutoff_ang,
        min_seq_sep=config.q_min_seq_sep, beta_per_nm=config.q_beta_per_nm, lam=config.q_lambda,
    )
    drmsd = ca_drmsd(traj, native_reference, selection=config.drmsd_selection)
    plots.plot_timeseries(x, q, "Q (fraction native contacts)", adir / "q.png",
                          "native-contact fraction", xlabel=xlabel, ylim=(0.0, 1.02))
    plots.plot_timeseries(x, drmsd, "dRMSD to design (Å)", adir / "drmsd.png",
                          "distance-RMSD to design", xlabel=xlabel)
    summary["q_native_contacts"] = {
        "mean": float(q.mean()), "min": float(q.min()), "final": float(q[-1]),
    }
    summary["drmsd_ang"] = {
        "mean": float(drmsd.mean()), "max": float(drmsd.max()), "final": float(drmsd[-1]),
    }
    arrays["q_native_contacts"] = q
    arrays["drmsd_ang"] = drmsd

    lddt_summary = {}
    for atom_set in config.lddt_atom_sets:
        selection, label = LDDT_SPEC[atom_set]
        glob, per_res, lddt_resseq = lddt_to_reference(
            traj, native_reference, selection,
            inclusion_radius_ang=config.lddt_inclusion_radius_ang,
        )
        per_res_mean = per_res.mean(axis=0)
        plots.plot_timeseries(x, glob, f"lDDT ({label})", adir / f"lddt_{atom_set.value}.png",
                              f"lDDT to design ({label})", xlabel=xlabel, ylim=(0.0, 1.02))
        plots.plot_per_residue(lddt_resseq, per_res_mean, f"mean lDDT ({label})",
                               f"per-residue lDDT ({label})",
                               adir / f"lddt_{atom_set.value}_per_residue.png", ylim=(0.0, 1.02))
        arrays[f"lddt_{atom_set.value}"] = glob
        arrays[f"lddt_{atom_set.value}_per_residue"] = per_res_mean
        arrays["lddt_resseq"] = lddt_resseq  # shared across atom sets (all protein residues)
        lddt_summary[atom_set.value] = {
            "label": label, "mean": float(glob.mean()),
            "min": float(glob.min()), "final": float(glob[-1]),
        }
    summary["lddt"] = lddt_summary

    return summary, arrays
