"""Interactive HTML overview report (Plotly) — the human-facing view of an analysis run.

A pure function of what the composers already produce: the `summary` dict (what goes to the
report JSON) and the `arrays` dict (what goes to the npz). So the report can be rebuilt from an
analysis in hand and never re-derives anything. Figures are interactive Plotly; the output is a
single **self-contained** HTML (plotly.js embedded once), so it opens offline and archives cleanly
into `storage_directory`.

Plotly is an optional dependency (`[report]` extra) — lazily imported, fail-loud with the install
hint if absent.
"""

from __future__ import annotations

import datetime as _dt
import html
import re
from pathlib import Path

import numpy as np


def _require_plotly():
    try:
        import plotly  # noqa: F401
        import plotly.graph_objects as go
        import plotly.io as pio
    except ImportError as exc:  # optional dependency — fail loud with the fix
        raise ImportError(
            "the HTML analysis report needs plotly (optional). "
            "Install it: pip install 'openmm_pipelines[report]'."
        ) from exc
    return go, pio


def _line(go, x, y, xlabel, ylabel):
    fig = go.Figure(go.Scatter(x=np.asarray(x), y=np.asarray(y), mode="lines"))
    fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)
    return fig


def _multiline(go, x, series: dict, xlabel, ylabel, yrange=None):
    fig = go.Figure()
    for name, y in series.items():
        fig.add_trace(go.Scatter(x=np.asarray(x), y=np.asarray(y), mode="lines", name=str(name)))
    fig.update_layout(xaxis_title=xlabel, yaxis_title=ylabel)
    if yrange is not None:
        fig.update_yaxes(range=yrange)
    return fig


def _dssp_heatmap(go, codes):
    """Interactive per-residue × per-frame DSSP map with a discrete category colourbar."""
    from .dssp import DSSP_CATEGORIES

    n = len(DSSP_CATEGORIES)
    labels = [f"{code.strip() or '~'} {label}" for code, label, _ in DSSP_CATEGORIES]
    colorscale = []
    for i, (_, _, c) in enumerate(DSSP_CATEGORIES):  # discrete bands, one per category
        colorscale += [[i / n, c], [(i + 1) / n, c]]
    z = np.asarray(codes).T  # (n_residues, n_frames)
    fig = go.Figure(go.Heatmap(
        z=z, zmin=-0.5, zmax=n - 0.5, colorscale=colorscale,
        customdata=np.array(labels)[z],
        hovertemplate="residue %{y}<br>frame %{x}<br>%{customdata}<extra></extra>",
        colorbar=dict(title="DSSP", tickmode="array",
                      tickvals=list(range(n)), ticktext=labels),
    ))
    fig.update_layout(xaxis_title="frame", yaxis_title="residue (index)")
    return fig


def _collect(arrays: dict, pattern: str) -> dict:
    """{label: array} for every arrays key matching `pattern` (one capture group = the label)."""
    out = {}
    rx = re.compile(pattern)
    for key, val in arrays.items():
        m = rx.fullmatch(key)
        if m:
            out[m.group(1)] = val
    return out


def _sections(go, summary: dict, arrays: dict, x_key: str, x_label: str) -> list[tuple[str, object]]:
    x = np.asarray(arrays[x_key])
    sections: list[tuple[str, object]] = []

    if "rmsd_ang" in arrays:
        sections.append(("Cα RMSD to design", _line(go, x, arrays["rmsd_ang"], x_label, "RMSD (Å)")))
    if "rg_ang" in arrays:
        sections.append(("Radius of gyration", _line(go, x, arrays["rg_ang"], x_label, "Rg (Å)")))
    if "resseq" in arrays and "rmsf_ang" in arrays:
        sections.append(("Per-residue flexibility (RMSF)",
                         _line(go, arrays["resseq"], arrays["rmsf_ang"], "residue", "Cα RMSF (Å)")))

    ss = {k[3:]: arrays[k] for k in ("ss_helix", "ss_sheet", "ss_coil", "ss_na") if k in arrays}
    if ss:
        sections.append(("Secondary structure",
                         _multiline(go, x, ss, x_label, "fraction of residues", yrange=[0, 1])))
    if "dssp_codes" in arrays:
        sections.append(("Secondary structure per residue", _dssp_heatmap(go, arrays["dssp_codes"])))

    clustering = summary.get("clustering")
    if clustering and clustering.get("populations"):
        pops = clustering["populations"]
        design_like = clustering.get("design_like_cluster", -1)
        colors = ["#2ca02c" if i == design_like else "#4c78a8" for i in range(len(pops))]
        figp = go.Figure(go.Bar(x=list(range(len(pops))), y=pops, marker_color=colors))
        figp.update_layout(xaxis_title="cluster (rank by population)", yaxis_title="frames")
        label = "Cluster populations" + (" (green = design-like basin)" if design_like >= 0 else "")
        sections.append((label, figp))
    if "cluster_labels" in arrays and "cluster_frame_indices" in arrays:
        xi = x[np.asarray(arrays["cluster_frame_indices"])]
        figs = go.Figure(go.Scatter(x=xi, y=np.asarray(arrays["cluster_labels"]), mode="markers",
                                    marker=dict(size=5)))
        figs.update_layout(xaxis_title=x_label, yaxis_title="cluster (rank)")
        sections.append(("Cluster assignment across the ensemble", figs))

    if summary.get("interchain"):
        mindist = _collect(arrays, r"interchain_mindist_ang_(.+)")
        contacts = _collect(arrays, r"interchain_contacts_(.+)")
        com = _collect(arrays, r"interchain_com_ang_(.+)")
        per_chain = _collect(arrays, r"chain(\d+)_rmsd_ang")
        if mindist:
            sections.append(("Interchain minimum distance",
                             _multiline(go, x, _relabel(mindist), x_label, "min distance (Å)")))
        if contacts:
            sections.append(("Interchain contacts",
                             _multiline(go, x, _relabel(contacts), x_label, "contact count")))
        if com:
            sections.append(("Interchain centre-of-mass distance",
                             _multiline(go, x, _relabel(com), x_label, "COM distance (Å)")))
        if per_chain:
            sections.append(("Per-chain fold retention (RMSD to design)",
                             _multiline(go, x, {f"chain {c}": v for c, v in per_chain.items()},
                                        x_label, "Cα RMSD (Å)")))

    if "neighbor_acceptance" in arrays:
        acc = np.asarray(arrays["neighbor_acceptance"])
        figa = go.Figure(go.Bar(x=[f"{i}-{i+1}" for i in range(len(acc))], y=acc))
        figa.update_layout(xaxis_title="adjacent rung pair", yaxis_title="exchange acceptance")
        figa.update_yaxes(range=[0, 1])
        sections.append(("Neighbour exchange acceptance", figa))
    if "replica_states" in arrays:
        states = np.asarray(arrays["replica_states"])
        series = {f"replica {r}": states[:, r] for r in range(min(4, states.shape[1]))}
        sections.append(("Replicas walking through state space",
                         _multiline(go, np.arange(states.shape[0]), series,
                                    "iteration", "ladder state index")))
    return sections


def _relabel(pair_series: dict) -> dict:
    """'0_1' -> 'chains 0-1' for interchain pair labels."""
    return {f"chains {k.replace('_', '-')}": v for k, v in pair_series.items()}


def _fmt(v) -> str:
    return f"{v:.3g}" if isinstance(v, float) else str(v)


def _metric_rows(summary: dict) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    frames = summary.get("n_frames", summary.get("n_frames_demuxed"))
    if frames is not None:
        rows.append(("frames", str(frames)))
    rmsd = summary.get("rmsd_to_design_ang")
    if rmsd:
        rows.append(("RMSD to design (Å)",
                     f"mean {_fmt(rmsd['mean'])} · max {_fmt(rmsd['max'])} · final {_fmt(rmsd['final'])}"))
    rg = summary.get("rg_ang")
    if rg:
        rows.append(("Rg (Å)", f"mean {_fmt(rg['mean'])} ± {_fmt(rg['std'])}"))
    ssf = summary.get("ss_fraction_mean")
    if ssf:
        rows.append(("secondary structure (mean)",
                     " · ".join(f"{k} {_fmt(v)}" for k, v in ssf.items() if k != "na")))
    clustering = summary.get("clustering")
    if clustering:
        rows.append(("clusters",
                     f"{clustering['n_clusters']} · design-like basin "
                     f"{_fmt(clustering['design_like_fraction'])} of ensemble "
                     f"({clustering['method']}, {_fmt(clustering['cutoff_ang'])} Å)"))
    inter = summary.get("interchain")
    if inter:
        for pair, m in inter.get("chain_pairs", {}).items():
            rows.append((f"interchain {pair} min distance (Å)",
                         f"mean {_fmt(m['min_distance_ang']['mean'])} · "
                         f"min {_fmt(m['min_distance_ang']['min'])}"))
    if "mean_neighbor_acceptance" in summary:
        rows.append(("mean neighbour acceptance", _fmt(summary["mean_neighbor_acceptance"])))
    return rows


_CSS = """
:root { color-scheme: light; }
body { font-family: -apple-system, system-ui, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
       margin: 0 auto; max-width: 980px; padding: 24px 16px 64px; color: #1a1a1a; }
h1 { font-size: 1.5rem; margin: 0 0 4px; }
h2 { font-size: 1.05rem; margin: 28px 0 6px; color: #333; }
.sub { color: #666; font-size: .85rem; margin-bottom: 20px; }
table.metrics { border-collapse: collapse; width: 100%; margin: 12px 0 8px; font-size: .9rem; }
table.metrics td { padding: 6px 10px; border-bottom: 1px solid #eee; vertical-align: top; }
table.metrics td.k { color: #555; width: 34%; }
section { margin-bottom: 8px; }
details { margin-top: 24px; } summary { cursor: pointer; color: #555; }
pre { background: #f6f8fa; padding: 12px; border-radius: 6px; overflow-x: auto; font-size: .8rem; }
footer { margin-top: 40px; color: #999; font-size: .78rem; }
"""


def build_report_html(summary: dict, arrays: dict, x_key: str, x_label: str, title: str) -> str:
    """Build a self-contained interactive HTML report from a run's analysis `summary` + `arrays`."""
    import json

    import plotly

    go, pio = _require_plotly()
    sections = _sections(go, summary, arrays, x_key, x_label)

    divs = []
    for i, (heading, fig) in enumerate(sections):
        fig.update_layout(template="plotly_white", height=360,
                          margin=dict(l=60, r=20, t=10, b=50))
        body = pio.to_html(fig, full_html=False, include_plotlyjs=(i == 0),
                           config={"displaylogo": False, "responsive": True})
        divs.append(f"<section><h2>{html.escape(heading)}</h2>{body}</section>")

    metric_rows = "".join(
        f"<tr><td class='k'>{html.escape(k)}</td><td>{html.escape(v)}</td></tr>"
        for k, v in _metric_rows(summary)
    )
    generated = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    full_json = html.escape(json.dumps(summary, indent=2, default=str))

    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>{html.escape(title)}</title><style>{_CSS}</style></head><body>"
        f"<h1>{html.escape(title)}</h1>"
        f"<div class='sub'>openmm_pipelines analysis overview · generated {generated}</div>"
        f"<table class='metrics'>{metric_rows}</table>"
        f"{''.join(divs)}"
        f"<details><summary>Full summary (JSON)</summary><pre>{full_json}</pre></details>"
        f"<footer>Interactive figures: Plotly {plotly.__version__}. "
        "This report is a quick default overview, not the final analysis.</footer>"
        "</body></html>"
    )


def write_html_report(adir, summary: dict, arrays: dict, x_key: str, x_label: str,
                      title: str, name: str = "report.html") -> Path:
    """Write the interactive HTML report into `adir`. Returns the path."""
    path = Path(adir) / name
    path.write_text(build_report_html(summary, arrays, x_key, x_label, title))
    return path
