# Validating a proposed OpenMM XML port of a99SB-disp

Takes as input a candidate ffxml claiming to implement a99SB-disp, and specifies what must
be established about it. Nothing here depends on how the file was produced.

## The reference

Truth is the authors' published parameter set. Both of their distributions may be used, and
agreeing with both is stronger than agreeing with either, since they were produced
separately:

- the Amber distribution — `leaprc.protein.ff99SBdisp` with its main-chain, N-terminal and C-terminal libraries and `frcmod.ff99SBdisp`
- the GROMACS distribution — `a99SBdisp.ff` with `ffnonbonded.itp`, `ffbonded.itp` and `aminoacids.rtp`

Build the reference `System` from whichever distribution with an established reader —
`AmberPrmtopFile` on a `tleap` prmtop, or `GromacsTopFile` on a `grompp -pp` topology —
and compare the candidate against it.

The two distributions agree exactly on the defining off-diagonal term, which fixes the
single most important reference value. The Amber `LJEDIT` record
`OB HB 0.841846 0.30 0.841846 0.30` combines as Rmin_ij = 1.683692 Å, giving
σ = Rmin_ij/2^(1/6) = 1.500000 Å = **0.150000 nm** and ε = 0.30 kcal/mol =
**1.2552 kJ/mol**; the GROMACS `[ nonbond_params ]` entry gives σ = 0.150 nm,
ε = 1.2552 kJ/mol directly.

### Water reference values

The water model is small enough to specify completely, and both distributions agree on
every value:

| Quantity | Value |
|---|---|
| O σ, ε | 0.316500 nm, 0.998989 kJ/mol (Amber Rmin/2 = 1.776296 Å, ε = 0.238764 kcal/mol) |
| H and M σ, ε | zero |
| charges O / H / H / M | 0, +0.59, +0.59, −1.18 e |
| geometry | OH 0.9572 Å, HOH 104.52°, M on the HOH bisector 0.1546 Å from O |
| virtual site | three-particle average, weights 0.736124464 / 0.131937768 / 0.131937768 |

The virtual-site weights and the M-site distance are mutually consistent and should both be
asserted: 2 × 0.131937768 × 0.9572 × cos(52.26°) = 0.1547 Å against the Amber `OW-EP` bond
length of 0.1546 Å.

## Completeness

What the candidate must contain, from the reference parameter set:

| Element | Count |
|---|---|
| atom types | 79 |
| off-diagonal LJ pairs | 1 (OB–HB) |
| explicit 1-4 pair types | 0 — 1-4 terms are generated |
| bond types | 130 |
| angle types | 330 |
| dihedral types | 333 |
| residue templates | 95, including 23 N-terminal and 24 C-terminal variants, 4 caps and ions |
| water | 4-site, rigid, with a three-particle-average virtual site |

A candidate missing whole categories will usually fail loudly at template matching. A
candidate missing individual entries will not, which is what the per-term comparisons below
are for.

## Tier 0 — self-consistency

No reference needed; catches construction errors before any comparison.

1. Every class referenced in `<Residues>` exists in `<AtomTypes>`, and every class appearing in any system — protein, water and ions alike — has an `<Atom>` entry in a `<LennardJonesForce>` section. OpenMM's generator looks up every present type in that table, so a missing entry fails at `createSystem` rather than silently; it merges multiple `<LennardJonesForce>` tags but raises if their `lj14scale` or `useDispersionCorrection` disagree, so a separately supplied water file must declare the same `lj14scale`.
2. Each residue template's net charge equals its formal charge to 10⁻⁶ e, and total system charge is integral to 10⁻⁵ e. This is the only cheap guard on roughly 1400 per-template partial charges.
3. `ForceField.createSystem` succeeds on every fixture without a template-match error.
4. The off-diagonal term is present as `<NBFixPair>` under `<LennardJonesForce>` with class1/class2 equal to OB and HB and the reference σ and ε above, to 10⁻⁶ relative.
5. Classes OB and HB are assigned to backbone carbonyl oxygens and backbone amide hydrogens in every protein residue template, including all terminal variants. Count per residue — a reference DHFR topology carries 183 OB and 171 HB atoms across 184 residues, so the counts are a sharp check.

## Tier 1 — parameter equivalence

The core tier, and the only one that can establish the candidate is faithful.

**Compare `System` objects, not energies.** Build each fixture from the candidate ffxml and
from the reference topology, both inside OpenMM with identical `createSystem` arguments,
then compare parameters directly. Energy comparison across engines requires matching
cutoffs, Ewald coefficients, dispersion-correction conventions and energy grouping;
parameter comparison requires matching nothing, so any difference is unambiguously a fault
in the candidate.

Compare, per fixture:

1. Per-atom `charge`, `sigma`, `epsilon`, matched by atom index.
2. The complete set of `NonbondedForce` exceptions as `(i, j, chargeProd, sigma, epsilon)`. This one comparison covers the 1-4 scale factors, the exclusion scheme, and whether the off-diagonal term reaches 1-4 interactions.
3. `HarmonicBondForce` and `HarmonicAngleForce` as sets of `(atoms, equilibrium, k)`.
4. `PeriodicTorsionForce` as a multiset of `(atoms, periodicity, phase, k)`. Dihedral wildcard-precedence errors surface only here.
5. Which atom pairs the off-diagonal term applies to, and its parameters.
6. Constraints: count and `(i, j, length)`, under whichever `rigidWater` and `constraints` settings you will actually run.
7. Virtual sites: for each, the site type, its parent particles and its weights. For the 4-site water these must be a three-particle average with weights summing to 1.

Item 5 is separate from item 1 for a reason: the off-diagonal term is worth −103.6 kJ/mol
on DHFR, consistently in vacuum and in water, so its total absence would show up in an
energy comparison — but its application to the wrong set of atom pairs would not, since the
magnitude is unchanged. Only the pair list detects that.

Items 6 and 7 matter because they are invisible to an energy comparison at fixed
coordinates — a wrong constraint length or virtual-site weight changes dynamics without
changing a single-point energy. Comparing them statically is strictly better than probing
them with a simulation.

**Fixtures.** Capped tripeptides ACE-X-NME for all 20 residues; the separate-template
variants HID/HIE/HIP, CYS/CYM/CYX/CME, ASH, GLH, LYN, HYP, ORN, DAB; GLY and PRO
explicitly for their distinct backbone templates; one disulfide-bonded CYX pair; and every
N- and C-terminal variant, since those are 47 of the 95 templates and an error there
affects every real protein. Tripeptides rather than one long construct, so a failure
localises to a residue. Use the reference topology's own coordinate output as the input to
the candidate, so atom names and ordering agree and comparison can go index by index.

## Tier 2 — conventions to pin down

Four places where two descriptions of the same force field can differ without any term
count changing.

**Whether the off-diagonal pair applies to 1-4 interactions.** A real fork: applying the
off-diagonal values scaled by the 1-4 factor gives σ = 0.150 nm, ε = 0.6276 kJ/mol, while
generating 1-4 terms from the diagonal parameters gives σ = 0.2015, ε = 0.1201.
**Energy comparison cannot resolve this** — summed over all 171 OB–HB 1-4 pairs in DHFR the
two conventions differ by 0.263 kJ/mol, against a residual cross-engine gap of 0.733 kJ/mol
on a 22,227 kJ/mol total, because at real 1-4 separations (0.245–0.316 nm, median 0.313)
both potentials lie in a shallow attractive tail. Determine the reference convention by
reading the generated 1-4 parameters out of the reference topology directly, then assert
the candidate matches. Tier 1 comparison 2 enforces it thereafter.

**1-4 scale factors.** `lj14scale` and `coulomb14scale` must correspond to the reference
SCNB and SCEE of 2.0 and 1.2 — equivalently `fudgeLJ` 0.5 and `fudgeQQ` 0.8333333.

**Improper atom ordering.** Amber-family impropers place the central atom third in the
GROMACS convention and first in OpenMM's `<Improper>`. Getting this wrong relocates the
planarity restraint without changing the term count. Check on a fixture containing an
aromatic ring and a backbone carbonyl.

**Exclusions around rings and disulfides.** Bond-graph exclusion generation agrees with
`nrexcl = 3` for linear chains but can diverge where a bond closes a ring or crosses a
disulfide. Tier 1 comparison 2 on the CYX pair and the aromatic fixtures covers it.

## Tier 3 — energy equivalence in OpenMM

Same fixtures, same `createSystem` settings, compare per-force-group potential energies at
a fixed conformation. If Tier 1 passed this should agree to round-off, so its value is as a
regression guard rather than as evidence. Disagreement here with Tier 1 passing points at
something outside the compared parameters — constraints, virtual sites, or force grouping.

## Tier 4 — assembly sanity

Not an audit of the published force field, which is validated against its own benchmark set
and is not the candidate's responsibility. The purpose here is narrow: confirm the
assembled `System` integrates, which a single-point parameter comparison cannot establish.

1. **Short NVE energy conservation**, on a small solvated box. A virtual site whose weights are correct but whose force redistribution is misconfigured, a failing constraint, or a unit error in one place will reproduce every parameter in Tier 1 and still drift. This is the sharpest cheap check for that class of fault. Note that it also tests integrator, timestep and PME settings, so a failure needs isolating before it is attributed to the candidate.
2. **Water-box density**, against the published a99SB-disp value. Tests the water model, its virtual site and the dispersion treatment together against a sharp reference number, independently of the protein.

Deliberately excluded: folded-protein stability runs. Ubiquitin and GB3 remain folded under
force fields with substantial backbone-balance errors, so a short run has little power to
discriminate a faulty candidate from a correct one, and a run long enough to discriminate
costs more than the information is worth. Peptide helicity is discriminating but requires
converged sampling, and tests the force field rather than the port.

## Tier 5 — negative controls

A suite that has never failed is unproven. Introduce each of the following into a copy of
the candidate and confirm a named test fails:

1. Delete the `<NBFixPair>` entry.
2. Perturb one partial charge by 0.01 e in one residue template.
3. Change `lj14scale` from 0.5 to 1.0.
4. Retarget one residue's backbone O from class OB to class O — the likeliest real fault, since it is what a single missed residue produces.
5. Swap the improper atom ordering in one residue.
6. Change one dihedral's `k` by 1%.

## Acceptance criteria

| Quantity | Tolerance |
|---|---|
| per-atom charge | exact (10⁻¹² e) |
| per-atom σ, ε | 10⁻⁶ relative |
| exception parameters | 10⁻⁶ relative |
| NBFixPair σ, ε | 10⁻⁶ relative against the reference values |
| bond, angle, torsion parameters | 10⁻⁶ relative |
| term counts, per force | exact |
| residue net charge | 10⁻⁶ e of formal charge |
| OB/HB atom counts per residue | exact |
| constraint lengths and count | exact |
| virtual-site weights and parents | exact |
| per-force-group energy, candidate vs reference | 10⁻⁸ relative |
| NVE energy drift | within the drift of the same protocol run on a stock OpenMM force field |
| water density | within the published a99SB-disp value and its uncertainty |

## Sufficiency

Tiers 0, 1 and 2 establish that the candidate is a faithful translation of the reference
parameters, and are sufficient for that claim. Tier 5 establishes that those tiers test
what they claim, and is therefore not optional. Tier 3 is a cheap regression net. Tier 4
confirms the result integrates, which the static tiers cannot show.

Using two independent references — the Amber and GROMACS distributions, which were
produced separately and agree exactly on the defining off-diagonal term — covers the
residual risk that the reference topology itself was built wrong, which is the only failure
mode a complete parameter comparison would otherwise miss.
