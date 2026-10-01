# openmm_pipelines Knowledgebase

A living record of the `openmm_pipelines` library: what we've done, why we decided
what we decided, the shape of the code, and the pitfalls we've hit. Scoped to
**this package only** (not the sibling repos), git-tracked, maintained by Claude
across sessions.

## Files

| File | What it holds |
|------|---------------|
| [`CODEBASE_SHAPE.md`](CODEBASE_SHAPE.md) | Architecture map + mental model — the prepare/produce seam, module layout, data flow. The *where-things-are* and *why-architectural*. Updated when structure changes. |
| [`DECISIONS.md`](DECISIONS.md) | Design decision records. One entry per significant decision, dated, with the *why* and the alternatives rejected. Highest-value part of the KB. |
| [`HISTORY.md`](HISTORY.md) | Dated session log — one entry per working session/day. Reverse chronological (newest on top). Terse and factual. |
| [`GOTCHAS.md`](GOTCHAS.md) | Full explanation + fix for each discovered pitfall. The one-line index in the repo-root `CLAUDE.md` points here. |
| [`TODO.md`](TODO.md) | Concise backlog — known missing features, gaps, improvements. One terse bullet per item. |

The repo root keeps only `README.md`, `CLAUDE.md`, and `SKELETON_v2.md` (the
design sketch the package is being built from).

## Maintenance protocol (for Claude)

**Skim `HISTORY.md` + `DECISIONS.md` at the start of substantive work.**

At the **end of a working session** where something non-trivial happened:

1. **`HISTORY.md`** — add a dated entry: what we worked on, what changed, what's
   still open. A log, not a report.
2. **`DECISIONS.md`** — if a real design decision was made (an architectural
   choice, a tradeoff, a reversal), add a record. Capture *why* and what was
   rejected, not just what was chosen.
3. **`CODEBASE_SHAPE.md`** — update only if the code's shape or public API
   actually changed. Don't churn it otherwise.
4. **`GOTCHAS.md`** — when a new pitfall is discovered (a surprising OpenMM
   behavior, a cluster quirk, a wrong assumption that caused a failure), add its
   detail here immediately AND add a one-line entry to the `CLAUDE.md` index.
   Both, every time.
5. **`TODO.md`** — jot known-but-not-yet-done work as a terse one-liner.

Rules: **be terse** (low cognitive load is a project value). Prefer appending
over rewriting history. Never fabricate rationale — if the *why* isn't known, say
so. For "where is X", use `codegraph_explore` / `rg -n` against the live tree
rather than a hand-maintained index.
