# Open question: which force-field energy terms should REST2 scale, and by what power of λ?

## Purpose

This document poses an open methodological question for literature review. It is
intentionally neutral: it does not assume an answer. The goal is to determine, from the
primary literature and the documented behaviour of established implementations, which
terms of a molecular-mechanics force field should be scaled in REST2 (Replica Exchange
with Solute Tempering, "version 2" / solute-scaling form), and with what dependence on the
scaling factor.

## Setting

- Method: REST2 / replica-exchange solute tempering, in the solute-scaling formulation
  (Wang, Friesner, Berne, *J. Phys. Chem. B* 2011), as opposed to the original 2005 REST.
- System: a folded/designed protein (optionally a protein–peptide complex) in explicit
  solvent, NVT, single physical temperature T₀ across all replicas; replica *m* carries a
  scaling factor λ_m = T₀/T_m ∈ (0, 1], with λ = 1 the unscaled physical Hamiltonian.
- Engine: OpenMM. Force fields of interest: AMBER (ff14SB, ff19SB), CHARMM36m, and
  a99SB-disp. Analysis is performed only on the λ = 1 (reference) state.

In REST2 the potential energy is commonly partitioned by which subsystem the interaction
involves:

```
E(X) = λ · E_solute–solute  +  √λ · E_solute–solvent  +  E_solvent–solvent
```

The practical implementation question is: **which specific force-field terms are taken to
belong to E_solute–solute (and E_solute–solvent), and therefore get scaled — and is each
scaled by λ, by √λ, or left unscaled?**

## The core question

For each force-field energy term below, what does the original REST2 formulation, and what
do established, widely-used implementations, actually do — scale it, and by what power of
λ, or leave it unscaled? Where implementations disagree, what is the stated rationale, and
is there a recognized "standard" or "recommended" choice?

| Term | Scaled? | Power of λ? | Notes to resolve |
|---|---|---|---|
| Electrostatics / partial charges | ? | ? | commonly discussed as √λ per solute atom (→ λ for solute–solute, √λ for solute–solvent) — confirm |
| Lennard-Jones (van der Waals) | ? | ? | which parameter is scaled (ε? σ? C6/C12?) and with what power |
| Bonds (harmonic stretch) | ? | ? | **principal open point** — some tools scale, some do not |
| Angles (harmonic bend) | ? | ? | **principal open point** — some tools scale, some do not |
| Urey–Bradley (CHARMM 1-3) | ? | ? | treated as bond/angle or separately? |
| Proper dihedrals | ? | ? | |
| Improper dihedrals | ? | ? | scaled like proper torsions, grouped with bonds/angles, or unscaled? |
| CMAP (CHARMM backbone cross-term) | ? | ? | scaled, unscaled, or implementation-dependent? |
| 1-4 scaled nonbonded (exceptions/pairs) | ? | ? | electrostatic and LJ 1-4 components |
| Bonded terms straddling the solute boundary | ? | ? | e.g. a torsion with some atoms in solute, some not — unscaled, √λ, or excluded by construction? |

## Implementation-landscape question

What is the documented behaviour of the established REST2 implementations, and do they
agree? Of particular interest:

- The original REST2 paper (Wang, Friesner, Berne 2011) — what did it scale?
- GROMACS + PLUMED `partial_tempering` (the Bussi hamiltonian-replica-exchange workflow).
- NAMD solute scaling (including the `soluteScalingAll` / separate charge and vdW factor
  options).
- OpenMM-ecosystem implementations: the OpenMM Cookbook REST tutorial, `perses`, and
  `femto` (Psivant).
- AMBER and CHARMM native REST/solute-scaling facilities, if any.

Where these differ (e.g. on bonds, angles, impropers, CMAP, or boundary-straddling terms),
is the difference a deliberate methodological choice with a cited justification, an
accident of what each tool's topology parser can reach, or simply unimplemented
functionality?

## Consequences question

- Does the choice of which terms to scale affect only **sampling efficiency** (ladder
  acceptance, rate of barrier crossing, number of replicas needed), or can it affect the
  **correctness of the reference (λ = 1) ensemble** that is actually analyzed?
- What are the documented physical arguments for or against scaling the stiff internal
  coordinates (bonds, angles, impropers)? For example: do they contribute to the
  conformational barriers REST2 aims to lower; does scaling them perturb equilibrium
  geometry, vibrational density of states, or integration stability; and is there evidence
  of artifacts from scaling or from not scaling them?
- Is there a difference in recommended practice between tempering for conformational
  sampling of a folded structure versus other goals (binding, disordered states)?

## Force-field-specific questions

- **CHARMM36m in OpenMM** represents the nonbonded LJ as a `CustomNonbondedForce` with an
  NBFIX (off-diagonal, per-atom-type) table, the 1-4 LJ as a `CustomBondForce`, impropers
  as a `CustomTorsionForce`, and the backbone cross-term as a `CMAPTorsionForce`. Given
  that LJ parameters are keyed by shared atom types rather than by individual atoms, what
  is the correct and standard way to scale only the solute's contribution (especially when
  the solute is a subset, e.g. one chain of a complex) without perturbing non-solute atoms
  that share the same atom types?
- **a99SB-disp** includes a modified backbone hydrogen-bond / dispersion correction (an
  off-diagonal, pair-specific Lennard-Jones / NBFIX-type term) in addition to its modified
  TIP4P-D-based water. How should such a term be treated under REST2 solute scaling, and
  how have modified-dispersion / IDP force fields been combined with REST2 in the
  literature?

## What a useful answer would contain

- The term-by-term scaling actually used in the foundational REST2 references.
- A clear statement of where the major implementations agree and where they diverge, with
  citations.
- The stated rationale for scaling vs. not scaling bonds, angles, impropers, and CMAP.
- Any evidence on whether these choices are consequential for results versus only for
  efficiency.
- Guidance specific to CHARMM-family (NBFIX/CMAP) and modified-dispersion (a99SB-disp)
  force fields under REST2.
