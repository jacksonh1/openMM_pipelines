# What to scale in REST2: CHARMM36m and a99SB-disp

λ = β_m/β_0 = T_0/T_m, with λ = 1 the unscaled reference replica that is analysed.
"Solute" means the atom set you have chosen to heat, which may be a subset of one
chain.

## The rule all of it reduces to

Two per-atom factors and one counting rule generate every entry in the tables below.

1. Every solute atom gets a charge factor **√λ** and a Lennard-Jones ε factor **λ**. σ is never touched. Standard combination rules then deliver λ for solute–solute pairs, √λ for solute–solvent pairs and 1 for solvent–solvent, for both electrostatics and dispersion, with no separate code path.
2. A bonded term is scaled by **λ** if all of its atoms are solute, by **√λ** if only some are, and by 1 if none are.
3. **Bonds and angles are not scaled.** This is the original REST2 prescription (Wang, Friesner & Berne 2011, who report that scaling them does not help sampling), the PLUMED default, the NAMD default and the femto default. Scaling them is exact but costs acceptance ratio: it adds of order k_BT of fluctuating energy per stiff coordinate to the quantity that sets the exchange probability, for coordinates that do not gate conformational transitions.

Scaled terms therefore are: electrostatics, Lennard-Jones, 1-4 electrostatics,
1-4 Lennard-Jones, proper dihedrals, impropers, and CMAP where present.
Unscaled: bonds, angles, Urey–Bradley, everything involving no solute atom, and
any restraint you have added.

## CHARMM36m in OpenMM

| OpenMM force | What it carries | Action |
|---|---|---|
| `HarmonicBondForce` | bonds, and Urey–Bradley 1-3 terms (OpenMM's CHARMM readers emit UB as additional harmonic bonds) | **leave unscaled** |
| `HarmonicAngleForce` | angles | **leave unscaled** |
| `PeriodicTorsionForce` | proper dihedrals | **λ** intra-solute, **√λ** straddling |
| `CustomTorsionForce` | harmonic impropers | **λ** (both the central atom and its substituents are in the same residue, so straddling terms arise only if your solute boundary cuts a residue) |
| `CMAPTorsionForce` | backbone φ/ψ cross-term | **λ**. Every CMAP term lies wholly inside one residue, so there is no √λ case. `CMAPTorsionForce` has no per-torsion coefficient, so implement it by adding λ-multiplied copies of the energy grids with `addMap` and pointing the solute's torsions at them — which gives one fixed-λ `System` per replica rather than an openmmtools global-parameter ladder. If you need the global-parameter ladder, leave CMAP unscaled, as PLUMED's `partial_tempering` does. |
| `CustomNonbondedForce` | LJ with the NBFIX off-diagonal table | **do not edit the table.** Add an `is_solute` per-particle parameter and prepend a multiplicative factor to the existing energy expression: `scale1*scale2*(<original expression>)` with `scale_i = 1 - is_solute_i + sqrt_lambda*is_solute_i`. The product gives λ, √λ and 1 for the three pair classes automatically, the tabulated per-type ε values are untouched, and because the flag is keyed to atom index rather than atom type, a solute that is one chain of a complex needs no type duplication. |
| `CustomBondForce` | 1-4 LJ | **λ** intra-solute, **√λ** straddling, applied as the same multiplicative prefactor per pair |
| `NonbondedForce` | electrostatics, plus 1-4 electrostatics as exceptions | charge offset of **√λ** per solute particle; exception offsets of **λ** when both atoms are solute, **√λ** when one is |

Two cautions specific to this force field. First, which forces actually appear
depends on how the `System` was built: without an NBFIX table you may get a plain
`NonbondedForce` carrying both ε and charge, in which case the per-particle offset
route covers LJ as well and no `CustomNonbondedForce` work is needed. Check
`system.getForces()` before writing any of this. Second, `femto` will refuse a
CHARMM36m system outright — it raises `ValueError` on any force type outside its
supported list, which includes `CMAPTorsionForce`, `CustomTorsionForce` and
`CustomBondForce` — so borrow its `CustomNonbondedForce` pattern but do not expect
`apply_rest` to run.

## a99SB-disp

| Term | Action |
|---|---|
| Modified backbone carbonyl-O / amide-H LJ pair | **λ**. Both partners are solute backbone atoms, so this is a pure solute–solute contribution and takes the same factor as any other intra-solute dispersion term. If your parameter files carry it as an off-diagonal pair entry, use the multiplicative `CustomNonbondedForce` treatment above; if they carry it as dedicated atom types, the per-atom ε factor handles it with no extra work. |
| a99SB-disp water (TIP4P-D-derived), including its virtual site charges | **leave unscaled** — it is solvent. Solute–water dispersion and electrostatics pick up √λ automatically through the per-atom factors. |
| Proper dihedrals **and impropers** | **λ** intra-solute, **√λ** straddling. In Amber-family force fields impropers are periodic torsions living in the same `PeriodicTorsionForce` as propers, so scaling that force scales both; separating them requires inspecting the topology, and there is no reason to. |
| Bonds, angles | **leave unscaled** |
| CMAP | not applicable — a99SB-disp is Amber-family and has no CMAP term |
| Electrostatics, LJ, 1-4 terms | as in the general rule: √λ on charge, λ on ε, exceptions λ/√λ |

The one published REST study with this force field is the REST3 work of Zhang, Liu
& Chen (2023), which ran a99SB-disp under GROMACS patched with PLUMED. That route
handles the pair-specific LJ table by duplicating atom types in the
`[nonbond_params]`/`[pairtypes]` sections and is a reasonable alternative to doing
it in OpenMM.

## Before production

Validate the reconstruction rather than trusting it. Build the scaled `System` at
λ = 1 and confirm that the per-force energies match the unmodified system force by
force, not just in total. Then build it at λ = 0.5 and confirm that exactly the
groups you intended have halved — electrostatics, LJ, 1-4 and torsions down by a
factor of two for a wholly-solute system, bonds and angles bit-identical.

Make sure the exchange criterion evaluates the actual energy of each configuration
under each replica's Hamiltonian, as the openmmtools multistate samplers and
GROMACS `-hrex` both do. The λ = 1 ensemble is exact for any choice of scaled
terms under that condition. It is **not** exact if you use the closed-form REST2
acceptance expression with an E_pp that does not match the set of terms you
actually scaled — that is the only route by which these choices can corrupt the
reference ensemble, and it fails silently.

## Deliberately omitted

Scaling bonds and angles. The OpenMM Cookbook REST tutorial and perses do scale
them (`CustomBondForce` and `CustomAngleForce` with prefactor
`is_rest*λ_rest² + is_inter*λ_rest + is_nonrest`, λ_rest = √λ), so you will see it
in reference code. It is not wrong, and with `constraints=HBonds` it affects only
heavy-atom bonds, but it buys nothing for conformational sampling and costs
replicas. Note also that scaling a harmonic term leaves r₀ and θ₀ where they
were — only the width of the distribution grows, as 1/√λ — so the objection is
about acceptance ratio, not about distorted geometry.
