# REST2 scaling for CHARMM36m and a99SB-disp

λ = β_m/β_0 = T_0/T_m, with λ = 1 the unscaled reference replica that is analysed.
"Solute" means the atom set you have chosen to heat, which may be a subset of one chain.

**The construction: one `System` per replica, every factor baked in.** Build the
unmodified `System` once, then for each λ in the ladder make a deep copy and multiply the
relevant parameters in place. No global parameters, no `addParticleParameterOffset`, no
runtime λ anywhere. Every number in the finished `System` is the number the integrator
will use, which is what makes it checkable.

## The factors

Two per-atom factors and one counting rule generate everything below.

1. Every solute atom gets a charge factor **√λ** and a Lennard-Jones ε factor **λ**. σ is never touched. Standard combination rules then give λ for solute–solute pairs, √λ for solute–solvent pairs and 1 for solvent–solvent, for both electrostatics and dispersion.
2. A bonded term is scaled by **λ** if all of its atoms are solute, by **√λ** if only some are, and by 1 if none are. Classify on the term's full atom set, never by residue.
3. **Bonds and angles are not scaled.** This is the original REST2 prescription (Wang, Friesner & Berne 2011, who report that scaling them does not help sampling), and the default in PLUMED, NAMD and femto. Scaling them is exact but costs acceptance ratio: it adds of order k_BT of fluctuating energy per stiff coordinate to the quantity that sets the exchange probability, for coordinates that do not gate conformational transitions.

### CHARMM36m

| OpenMM force | What it carries | Factor |
|---|---|---|
| `HarmonicBondForce` | bonds, and Urey–Bradley 1-3 terms (OpenMM's CHARMM readers emit UB as additional harmonic bonds) | 1 — leave alone |
| `HarmonicAngleForce` | angles | 1 — leave alone |
| `PeriodicTorsionForce` | proper dihedrals | λ intra-solute, √λ straddling |
| `CustomTorsionForce` | harmonic impropers | λ intra-solute, √λ straddling (straddling arises only if your solute boundary cuts a residue) |
| `CMAPTorsionForce` | backbone φ/ψ cross-term | λ intra-solute, √λ straddling. Each term spans **three** residues — φ = (C′ᵢ₋₁, Nᵢ, Cαᵢ, Cᵢ) and ψ = (Nᵢ, Cαᵢ, Cᵢ, Nᵢ₊₁) — so a boundary at a chain, domain or loop edge produces straddling terms at both ends of the selection. |
| `CustomNonbondedForce` | LJ with the NBFIX off-diagonal table | √λ per solute particle, applied multiplicatively — never edit the table |
| `CustomBondForce` | 1-4 LJ | λ intra-solute, √λ straddling |
| `NonbondedForce` | electrostatics, plus 1-4 electrostatics as exceptions | √λ on each solute charge; λ / √λ on exceptions |

Which of these forces exists depends on how the `System` was built: with no NBFIX table
you may get a plain `NonbondedForce` carrying both ε and charge, in which case it covers
LJ as well and there is no `CustomNonbondedForce` work to do. Call `system.getForces()`
and enumerate before writing anything.

### a99SB-disp

| Term | Factor |
|---|---|
| Modified backbone carbonyl-O / amide-H LJ pair | λ. Both partners are solute backbone atoms, so this is a pure solute–solute contribution. If the parameter files carry it as an off-diagonal pair entry, treat it like the NBFIX table above; if as dedicated atom types, the per-atom ε factor handles it with no extra work. |
| a99SB-disp water (TIP4P-D-derived), including virtual-site charges | 1 — it is solvent. Solute–water terms pick up √λ automatically through the per-atom factors. |
| Proper dihedrals **and impropers** | λ intra-solute, √λ straddling. In Amber-family force fields impropers are periodic torsions in the same `PeriodicTorsionForce` as propers, so scaling that force scales both; separating them would require inspecting the topology and there is no reason to. |
| Bonds, angles | 1 — leave alone |
| CMAP | not applicable — a99SB-disp is Amber-family and has no CMAP term |
| Electrostatics, LJ, 1-4 terms | general rule: √λ on charge, λ on ε, exceptions λ / √λ |

## Applying the factors

Compute each term's factor once, from the atom sets, and write it into the parameters.

- **`NonbondedForce` particles.** `setParticleParameters(i, sqrt_lam*q, sigma, lam*eps)` for solute atoms; leave the rest.
- **`NonbondedForce` exceptions.** `setExceptionParameters(i, p1, p2, f*chargeProd, sigma, f*eps)` with f = λ when both atoms are solute, √λ when one is. Both `chargeProd` and `eps` take the same f, which is what the per-atom factors imply.
- **`PeriodicTorsionForce`.** `setTorsionParameters(i, p1..p4, periodicity, phase, f*k)`.
- **`CustomTorsionForce` (impropers).** Multiply the force-constant per-torsion parameter by f, leaving the equilibrium angle alone. Read the parameter names off the force rather than assuming their order.
- **`CustomBondForce` (1-4 LJ).** Multiply the pair's energy coefficients by f. For an `acoef/r^12 - bcoef/r^6` form, multiply **both** coefficients, which scales ε by f and leaves σ unchanged.
- **`CustomNonbondedForce` (NBFIX).** Add one per-particle parameter, `s`, set to √λ for solute atoms and 1.0 for everything else, and prepend `s1*s2*` to the existing energy expression. The product gives λ, √λ and 1 for the three pair classes automatically; the tabulated per-type values are never touched; and because `s` is keyed to atom index, a solute that is one chain of a complex needs no type duplication. This is a per-particle number, not a global parameter.
- **`CMAPTorsionForce`.** Terms sharing a map may need different factors, so duplicate the maps rather than scaling in place. For each original map, add the copies you need with `addMap(size, [f*v for v in values])`, then point each torsion at the right copy with `setTorsionParameters(i, new_map_index, a1..a4, b1..b4)`. You only ever multiply the flat energy list element-wise, so the grid's index convention and periodicity are never your problem.

## Validation

Three checks, in this order. Each one catches a class of error the others cannot.

**At λ = 1, compare force by force.** Put every force in its own force group in both the
unmodified and the rebuilt `System` and confirm the energies match group by group, not
merely in total. A total-energy match hides compensating errors.

**At λ = 0.5 with the solute set to the whole system, check what halved.** Electrostatics,
LJ, 1-4 and torsions should come out at exactly half; bonds and angles bit-identical. This
confirms the factors reach the terms you intended and no others.

**At λ ≠ 1 with the real, partial solute, check the classification.** The first two checks
cannot see a misclassified term: at λ = 1 every factor is 1 regardless, and with the whole
system as solute every term is intra-solute. So derive the classification a second time,
independently, from the topology — set intersection of each term's atoms with the solute
selection — and assert the per-term factors match. Assert the counts too: the number of
√λ terms should equal the number of atom sets that genuinely cross the boundary. This
matters most for CMAP, which straddles at both ends of a contiguous selection.

Finally, make sure the exchange criterion evaluates the actual energy of each
configuration under each replica's Hamiltonian, as the openmmtools multistate samplers
and GROMACS `-hrex` both do. The λ = 1 ensemble is exact for any choice of scaled terms
under that condition. It is **not** exact if you use the closed-form REST2 acceptance
expression with an E_pp that does not match the set of terms you actually scaled — the
only route by which these choices can corrupt the reference ensemble, and it fails
silently.

A ladder of separately built `System`s is supported by openmmtools without special
measures: `MultiStateSampler._compute_replica_energies` groups states with
`states.group_by_compatibility` and takes one `Context` per group.

## Not scaled, deliberately

Bonds and angles. The OpenMM Cookbook REST tutorial and perses do scale them, so you will
see it in reference code. It is not wrong, and with `constraints=HBonds` it affects only
heavy-atom bonds, but it buys nothing for conformational sampling and costs replicas.
Scaling a harmonic term also leaves r₀ and θ₀ where they were — only the width of the
distribution grows, as 1/√λ — so the objection is about acceptance ratio, not distorted
geometry.
