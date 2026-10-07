# REST2 kickoff note

Handoff for building the REST2 backend in a fresh session. Written 2026-10-06, right
after the T-REMD backend landed. REST2 reuses almost all of the T-REMD machinery; this
note captures what carries over, the genuinely new (solute λ-scaling) part, the
openmmtools API template verified live this session, and the open decisions.

Read first: `DECISIONS.md` (T-REMD: NVT + geometric ladder + demux), `CODEBASE_SHAPE.md`
(the prepare/produce seam), and the built `remd/` as the template.

---

## Why REST2 (the motivation you just saw)

Plain T-REMD hit the explicit-solvent wall in the demo: adjacent-replica acceptance
~ exp(-(Δβ·σ_E)²/2) and σ_E grows with the box heat capacity (∝ *all* degrees of
freedom), so a 5.4k-atom box needed ~14 rungs over just 300–360 K. The replica count
scales with √(total DOF) — impractical for big/hot ranges.

**REST2 (Replica Exchange with Solute Tempering, Wang–Friesner–Berne 2011)** fixes this
by tempering only the **solute**: instead of raising temperature, it scales down the
solute's interactions so only *solute* DOF enter the exchange acceptance. A handful of
replicas then span a wide effective-temperature range. Same objectives as T-REMD
(stability / flexibility / variant comparison) — this is the tool for explicit solvent.

## What carries over from T-REMD (do NOT rebuild)

- **prepare/produce seam** — `prepare(cfg)` is engine-agnostic; `REST2Config(PrepConfig)`
  subclasses the same base. Reuse `prepare()` unchanged.
- **Backend shape** — `rest2/production.run(eq, cfg) -> Path`, same signature as
  `md.run`/`remd.run`. `remd/production.py` is the direct template:
  - strip the barostat for NVT (`_nvt_system` pattern),
  - `openmmtools.multistate.ReplicaExchangeSampler` + `MultiStateReporter`,
  - `LangevinDynamicsMove(reassign_velocities=True)`,
  - set platform on `cache.global_context_cache.set_platform(platform, {"Precision":"mixed"})`,
  - `sampler.equilibrate()` (discarded) for per-replica equilibration,
  - `REMDRunReport`-style QC.
- **Production is NVT** (same reasoning as T-REMD; box already at density from prepare).
- **Analysis** — `analysis/remd.py` is almost entirely reusable: exchange diagnostics
  (`ReplicaExchangeAnalyzer.generate_mixing_statistics` → transition matrix, subdominant
  eigenvalue, statistical inefficiency), `demux_state(k)` to pull the fixed-λ reference
  ensemble (state 0 = full-strength solute = the design temperature), and
  `analysis/drift.drift_report` for RMSD/Rg/RMSF/SS. Factor the shared REMD/REST2 parts
  out of `analysis/remd.py` rather than copy them.
- **`[remd]` extra + gating** — same openmmtools dependency; keep `rest2/` and its
  analysis out of the eager `__init__` imports; keep `config.py`/ladder pure.
- **Gotchas** — `read_mixing_statistics` returns per-iteration matrices (sum axis 0);
  mdtraj wants float32 xyz/box. Both already in GOTCHAS.md.

## The new part: solute λ-scaling (the one genuinely involved piece)

Replica m has effective temperature T_m (ladder from T_0 = reference up). Define
**λ_m = T_0 / T_m ∈ (0, 1]** (λ=1 at the reference replica = unscaled = the physical
ensemble you analyze). Scale ONLY the solute (protein) region so that:

- solute–solute interactions  ∝ λ_m
- solute–solvent interactions ∝ √λ_m
- solvent–solvent interactions unscaled

Standard implementation (scale the force-field terms, run every replica at the *same*
physical T_0):

- **Charges:** q_solute → √λ_m · q_solute. (Coulomb solute–solute ∝ λ; solute–solvent
  ∝ √λ automatically.)
- **LJ:** ε_solute → λ_m · ε_solute. (Lorentz–Berthelot ε_ij = √(ε_i ε_j): solute–solute
  LJ ∝ λ; solute–solvent ∝ √λ automatically.) σ unchanged.
- **Bonded solute-only terms** (bonds/angles/propers/impropers fully inside the solute):
  energy × λ_m.
- **CMAP (CHARMM36m):** the backbone CMAP correction is a solute intramolecular term —
  it MUST also be scaled by λ_m. Forgetting it silently gives a wrong REST2. This is the
  **CMAP safety gate**: assert the system's CMAP (if present) is included in the scaled
  set, or refuse. (See the CHARMM36m gotcha — pick the FF deliberately.)

Only solute–solvent cross-terms need √λ; getting the cross-term scaling right is the
crux. With the charge-√λ + ε-λ trick above it falls out of the combining rules, so the
hard case is really the solute–solute bonded + CMAP terms and the exceptions/1-4 scaling.

### Two implementation paths (recommend A first)

- **A — pre-scaled Systems (simplest, matches `_nvt_system`).** Build N independent
  Systems, each with the solute pre-scaled by its λ_m, wrap each in a plain
  `ThermodynamicState(system_m, T_0)`. `ReplicaExchangeSampler` over the list. Decompose
  eagerly: a pure `rest2/scaling.py::scale_solute(system, solute_indices, lam) -> System`
  primitive (single purpose, unit-testable on energies), composed by production. Downside:
  N systems in memory; fine for our replica counts.
- **B — CompoundThermodynamicState (scalable refactor, later).** One System whose solute
  forces carry a global λ parameter; each replica is a
  `states.CompoundThermodynamicState(base_state, [rest_state_m])` where the REST composable
  state is built on `states.GlobalParameterState` / `IComposableState` and sets λ. The
  "openmmtools way," less memory, more upfront force-authoring. Note: openmmtools ships NO
  turnkey REST2 — `openmmtools.alchemy` is for alchemical free-energy, **not** REST; use
  `GlobalParameterState` + `CompoundThermodynamicState` to roll it.

Validate either path by energy: at λ=1 the scaled System must match the unscaled System
(per-force-group), and solute–solvent should scale as √λ while solvent–solvent is
invariant. Tie this into the deferred force-field energy-match tool (see TODO).

## Config / ladder

- `REST2Config(PrepConfig)`: `max_effective_temperature_k`, `n_replicas` (or explicit
  λ ladder), same iteration/exchange fields as `REMDConfig`. `temperature_k` stays the
  physical T_0 (unscaled reference). Ladder in λ: geometric in effective-T then λ=T_0/T,
  so reuse `geometric_ladder` on effective temperatures.
- Need a **solute selection** (which atoms are "hot"): default = protein heavy + H
  (whole protein), reuse/extend `lib/selections.py`. For binder work you may want only
  one chain tempered — make the selection a config knob.

## Open decisions to make

- Solute region definition: whole protein vs a chain/interface subset (binder campaigns
  may want only the peptide hot).
- 1-4 (scaled exception) handling under charge/ε scaling — exceptions must track the
  scaled charges/ε consistently; verify against energy decomposition.
- Water model with CHARMM vs AMBER for REST2 campaigns (CMAP only matters for CHARMM36m).
- Whether to expose multi-λ demux (all rungs) or just the reference (λ=1), as in T-REMD.

## First concrete steps

1. `rest2/scaling.py::scale_solute(system, solute_indices, lam)` + energy tests
   (λ=1 identity; √λ / λ / invariant scaling; CMAP included-or-refuse).
2. `rest2/config.py::REST2Config` + pure tests (λ ladder, derived counts).
3. `rest2/production.py::run` from the `remd/production.py` template (path A).
4. Factor shared exchange/demux out of `analysis/remd.py`; add REST2 analysis.
5. `examples/REST2_demo/` (should mix with far FEWER replicas than the 14-rung T-REMD demo
   — that contrast IS the showcase).
