# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

---

## Knowledgebase — read at session start, update at session end

`knowledgebase/` (git-tracked) is the project's persistent memory for openmm_pipelines:
- `HISTORY.md` — dated session log
- `DECISIONS.md` — design decision records (the *why*)
- `CODEBASE_SHAPE.md` — architecture map
- `GOTCHAS.md` — discovered pitfalls (full detail; the index below points here)
- `TODO.md` — known missing features, gaps, and improvements

**Skim `HISTORY.md` + `DECISIONS.md` at the start of substantive work**, and follow
the maintenance protocol in `knowledgebase/README.md` at the end of a session where
something non-trivial happened. Keep entries terse.

---

## Environment

This is an **isolated tool** with its own conda env — **`openmm_env`**. Use it for
everything here; do **not** use the sibling repos' envs (`fragforge`,
`xbl_pep_design`, `basic`). It has openmm 8.6.1, pdbfixer, mdtraj, numpy, scipy,
**openmmtools 0.26.0** (so the future REMD/REST2 backends can be tested in-place),
plus pydantic + pytest. The package is installed editable (`pip install -e .`).

```
conda activate openmm_env
python -m pytest -q
```

Note: `openmm` is *not* importable in `fragforge`/`basic`. The pure units
(`forcefield.py`, `config.py`, `density.py`) import no OpenMM and run on any env, but
run the suite in `openmm_env` so the solvation/dynamics modules are covered too.

---

## Project overview

**T-REMD**, **REST2** and **plain production MD** — for characterizing
**designed protein structures** on single-node GPU clusters (SLURM). Developed for the
Keating lab at MIT; configurable for any cluster via `site_config.sh`.

**Every run starts from a folded/designed input structure — never from an unfolded or
extended state.** The tools serve four objectives:

1. Stability/rigidity characterization of a given starting (designed) structure
2. Identifying flexible regions (per-residue)
3. Variant comparison — run the same protocol on several design variants and see which best retains its designed conformation
4. Bound-state sampling — simulate a complex in its bound pose and sample the bound ensemble

**Tool mapping:** plain MD is mainly for **#4** (and optionally #1, #2); T-REMD/REST2 is
primarily for **#1–3**. Neither is a folding-from-unfolded tool — the input pose is the
reference the analysis is measured against (RMSD = drift from the design, RMSF = local
flexibility).

### Critical REMD/REST2 concept

figure out whether trajectories are single coordinate sets moving through temperature space or whether they are thermodynamic ensembles at fixed temperatures

---

## Working style

### Push back when something doesn't make sense

If a request is scientifically or technically wrong — the wrong GROMACS flag, a misunderstanding of T-REMD, a workflow that would silently produce bad results — say so clearly before implementing it. A wrong simulation or analysis run wastes real compute time and can produce results that look plausible but are wrong. See the **Known gotchas** section below for specific discovered pitfalls.

### Suggest better alternatives

If there's a cleaner, faster, or more correct approach than what was asked for, say so and explain the tradeoff. Don't just implement what was asked if something clearly better exists. Include enough context for an informed decision.

### Prefer simple and scalable solutions

Solutions should work for 8, 48, or 128 replicas without special-casing. Prefer shell/Python idioms that stay readable as the codebase grows. If a task has a five-line solution and a fifty-line solution, understand why the complexity is or isn't justified before recommending it.

Don't hand-roll a standard algorithm — add the well-known library (scipy, etc.) and use it. I am fine with adding a dependency if it is widely used and well-maintained. Tell me when you add a dependency and why it is justified, then update the pyproject and/or yaml environment files.

If you see an existing style or convention in the codebase, don't just blindly follow it, particularly if it is a bad one. Bad patterns in quick exploratory/pilot analysis should not be propogated to production code. Ask for guidance or propose a better approach.

---

## Tests


---

## Code style: fail loudly

This project follows a "fail loudly" philosophy. Bugs that crash immediately are strongly preferred over bugs that silently produce wrong results.

### Core philosophy

- **Crashes are cheap; silent bugs are expensive.** Prefer code that crashes obviously when assumptions are violated over code that "handles" the violation by producing degraded output.
- **Don't paper over uncertainty.** If you're unsure whether a value can be None, empty, or wrong-typed, either ask, add an assertion, or leave a clearly-marked comment — never add a default to make the question go away.
- **Make illegal states unrepresentable.** Prefer types and structures where the invalid case can't be expressed, over runtime checks for the invalid case.

### Error handling

- **No bare `except:` or `except Exception:` clauses** unless the exception is logged AND re-raised, or the recovery is documented and intentional.
- **Don't catch exceptions just to log and continue.** If the operation failed, the caller needs to know.
- **No `.get(key, default)` patterns** unless the default is semantically meaningful, not just a way to avoid a `KeyError`.
- **No `value or fallback` shortcuts** (`x or []`, `x or {}`, `x or 0`) unless `None`/empty/zero is genuinely interchangeable with the fallback. These hide bugs where `x` was unexpectedly empty.
- **Don't add defensive `if x is not None:` checks** unless `None` is a real expected case. If `None` would be a bug, let it crash.

### Indexing and iteration

- **Prefer iteration over indexing.** Use `for item in items`, not `for i in range(len(items))`. When you need the index too, use `enumerate`.
- **Use `zip(strict=True)`** (Python 3.10+) so mismatched-length iterables crash instead of silently truncating.
- **Assert invariants before code that relies on them.** E.g., `assert len(a) == len(b)` before zipping when the lengths must match.

### Types and structure

- **Use dataclasses or TypedDicts, not raw dicts**, when the shape matters and is fixed.
- **Parse, don't validate, at boundaries.** Convert untrusted input into typed structures at the edge; the rest of the code should be able to assume the data is valid.

### When in doubt

- **Ask before adding error handling.** If tempted to wrap something in try/except, ask what the intended behavior is when it fails.
- **Flag assumptions explicitly.** If making an assumption about input shape, range, or type that isn't enforced by the types, leave a comment like `# ASSUMES: items is non-empty`.

---

## Known gotchas

Specific pitfalls discovered in this environment. **The full explanation + fix for each lives in
`knowledgebase/GOTCHAS.md`** — the list below is a one-line index (phrased as the trigger) so the
tripwire stays visible in every session; open the doc when one applies. **When a new pitfall is
discovered during work — a surprising GROMACS behavior, a cluster quirk, a wrong assumption that
caused a failure — add its detail to `GOTCHAS.md` immediately AND add a one-line entry here**
(both, every time). The goal is that the same mistake is never made twice.

- **SLURM: sbatch scripts are copied to a temp path** — `BASH_SOURCE[0]` inside an sbatch does *not* point to the repo; use `$SLURM_SUBMIT_DIR` or pass the path via `--export`.
- **Density convergence must be a plateau (slope) test, not a consecutive-segment difference** — segment-to-segment noise is well under `DENSITY_TOL_REL`, so a sustained drift just under tolerance passes every comparison while the box contracts several percent. `openmm_pipelines/core/density.py` (`density_converged`) fits a least-squares line through the trailing `density_min_seg` window.
- **Bash: a `[[ … ]] && echo` as the LAST line of a `set -e` script exits 1** — a false test makes SLURM report a successful job as FAILED; use a real `if` (also applies to `((n++))` / `grep -q` as the final statement).
- **Modeller.addSolvent rejects most water-model names** — `model='opc'` raises `ValueError: Unknown water model: opc`; pack by SITE COUNT (`_PACKING_BOX`), parameters come from the FF xml (`model='tip4pew'` + `opc.xml` yields real OPC). `_assert_water_topology` backstops it.
- **CHARMM36 != CHARMM36M** — different force fields (`charmm36.xml` vs `charmm36_2024.xml`); 36m is the Huang-2017 CMAP refinement, usually right for folding/stability. Pick deliberately.
- **CHARMM nonbonded is a potential-switch, not CHARMM's force-switch** — our 1.2/1.0 PME settings match the canonical ParmEd example but OpenMM's built-in `switchDistance` only does potential-switching (≈, not exact CHARMM vfswitch). Energy-match vs a reference (choderalab/OpenMMEnergyComparisons) before a big CHARMM campaign.
- **OpenMM checkpoint (.chk) is NOT portable** — loads only on identical System + Platform + OpenMM version + hardware; use a `saveState` XML (`state.xml`/`final_state.xml`) to move an equilibrated system or restart across nodes.
- **A new Context resets globals to the force DEFAULT** — `context.setParameter("k",0)` only affects that context; a fresh production context re-enables restraints at full strength. Bake release into the System with `force.setGlobalParameterDefaultValue(idx, 0.0)` (done in `prepare()`).
- **openmmtools `read_mixing_statistics()` returns PER-ITERATION matrices** — default shape is `(n_iter, n_states, n_states)`, not the documented `(n_states, n_states)`; `.sum(axis=0)` for cumulative acceptance before indexing, or it crashes on the ambiguous-truth-value comparison.
- **mdtraj Trajectory needs float32 xyz/box and takes no `unitcell_vectors` kwarg** — build from de-multiplexed REMD coords with `dtype=np.float32` and set the box via the `traj.unitcell_vectors` attribute; float64 makes `image_molecules` raise `Buffer dtype mismatch`.
