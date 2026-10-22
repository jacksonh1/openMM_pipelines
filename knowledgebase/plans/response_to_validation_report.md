# Response to the a99SB-disp OpenMM port validation report

Reviewing `validation_report.md` dated 2026-10-21, candidate `ffxml/a99SBdisp.xml` +
`a99SBdisp_water.xml`.

## Assessment

The validation found a real conversion fault — the generic `X-CT-CT-X` torsion over-applied
to sidechain χ on seven residues — that a single-point energy check would have passed and
shipped. The fault sat in the 1–2.5 kJ/mol range: below what a total-energy comparison
resolves, well above what biases a rotamer population. The two design choices that separated
it from the benign zero-force constant layered on top of it were a second independent
reference and a multi-geometry spread test. Both earned their cost and should stay.

The report's diagnosis is plausible and consistent with its own data. It is not yet
demonstrated: the repair depends on a specific claim about OpenMM's torsion-matching
semantics that the report infers from a post-fix energy match rather than testing directly.
That claim must be established on its own terms before the fix can be considered sound.

Three items block closure of the report itself; two are defects in the record rather than in
the port.

## Corrections required to the report

**§5.2 contradicts §6a and must be rewritten.** It classifies `cap_ARG/LYS/MET/PHE/TRP/TYR`
and `term_PHE/TRP/TYR` as benign, asserting the torsion energy matches to 1e-5 — these are
precisely the residues §6a then shows carried a real geometry-dependent error of 1–2.5
kJ/mol. Both cannot be true. §5.2 is not a finding; it is the record of the misjudgement that
nearly shipped the bug, and leaving it standing as a live conclusion is more dangerous than
an ordinary inconsistency. Rewrite it to state that those flags were initially misclassified
as representation differences and that §6a corrected the classification. Determine and state
which of the two measurements was wrong and why, rather than deleting one.

**Tier 4 is not done, and its acceptance criterion was widened to close it.** The stated
criterion was drift within that of the same protocol on a stock OpenMM force field. The
result is 2.46× stock (+0.0154 vs +0.0062 kJ/mol/ns/dof), and a "5× acceptance bar" appears
with no provenance, accommodating the result after the fact. The drift may be short-run CPU
noise, but a mis-wired virtual site and a marginal constraint tolerance produce exactly this
signature, and detecting that class of fault is the only reason the check exists. Revert the
bar and close the check properly.

**§12 is internally inconsistent** — it marks Tier 4 done while recommending the longer GPU
re-check. Pick one.

## Open questions on the harness

**Does `tier1_compare.py` treat an extra or missing nonzero term as a hard failure by
construction, or does it route such differences through the energy-spread classification?**
The bug escaped because an added nonzero term was judged benign. The spread test is the right
instrument for characterising a difference — conformation-independent means zero-force,
geometry-varying means a force error — but it must not be able to excuse a term-set mismatch.
If the classifier can still return "benign" for an added nonzero term, the same failure can
recur. Answer this by inspecting the comparator, and if the answer is the wrong one, change
it and re-run.

**Is the zero-force exclusion symmetric?** The comparator "excludes zero-force terms from both
sides", and after the fix the candidate deliberately carries k = 0 shadow entries. Establish
that this exclusion cannot mask a case where the candidate holds a shadow and the reference
holds a real term, or the reverse.

## Reference values to assert against

Derive each of these from the authors' source files yourself and assert the candidate matches.
They are listed here as targets, not as results — none of them should be taken as already
checked.

| Quantity | Source | Target |
|---|---|---|
| OB–HB off-diagonal | Amber `frcmod` `LJEDIT OB HB 0.841846 0.30 …`, combining as Rmin_ij = Rmin_i + Rmin_j, σ = Rmin_ij/2^(1/6) | σ and ε in nm and kJ/mol, cross-checked against the GROMACS `[ nonbond_params ]` entry |
| Generic `X-CT-CT-X` | `frcmod.ff99SBdisp` | the k value, in kJ/mol, that the suppressed generic would contribute per affected χ |
| Water O σ, ε | Amber `frcmod.tip4pd_disp` NONBON `OW`, Rmin/2 and ε in kcal/mol | converted to nm and kJ/mol, cross-checked against GROMACS `OW_tip4pd` |
| Water charges | Amber OFF library / GROMACS water itp | O, H, H, M |
| Water virtual site | GROMACS `[ virtual_sites3 ]` a and b | weights, cross-checked for consistency against the Amber `OW-EP` bond length given the OH distance and HOH angle |
| 1-4 scale factors | Amber SCNB and SCEE | `lj14scale` and `coulomb14scale` |

Where a quantity appears in both distributions, assert against both. Agreement between two
independently produced distributions is the only available check on the reference itself, and
it is the reason the dual-reference design exists.

## Remaining work

### 1. Demonstrate the torsion-suppression semantics directly

The repair rests on the claim that a specific k = 0 override suppresses a matching generic
wildcard under OpenMM's most-specific-wins rule, contributing nothing. Establish this in
isolation rather than inferring it from a post-fix energy match: construct a minimal force
field containing the generic alone, and the generic plus a specific k = 0 override, on the
same topology, and report the resulting torsion terms and force constants in each case.
Confirm the override both wins the match and contributes no term. If it does not behave this
way, the fix is wrong regardless of the energy agreement.

### 2. Close Tier 4 on the original criterion

Re-run NVE on GPU at a length where the stock comparison is not noise-limited, same protocol
and water count for candidate and stock, reporting drift per ns per degree of freedom with an
uncertainty. If the candidate still sits several-fold above stock, isolate before accepting
it: re-run with the M-site removed in favour of a 3-site water to test whether the virtual
site is implicated, tighten the constraint tolerance, and halve the timestep. Those three
discriminate virtual-site force redistribution from constraint tolerance from integrator
error.

### 3. Publish the post-fix Tier 1 record per fixture

The headline reports 42 pass / 7 flags / 0 fail, but §5's per-fixture detail is pre-fix.
Produce the post-fix table at the same granularity, including the seven C-terminal flags with
their term-count evidence, so the record matches the verdict.

### 4. Validate the water by parameter, not only by energy

The water was checked as `Nonbonded(total)` to 3.1e-7, which is defensible given the M-site
representation difference, but it is not the per-atom comparison the plan specifies — and the
water is the hand-written component, so transcription is its dominant risk. Compare per-atom
charge, σ and ε, virtual-site type, parents and weights, and constraint lengths directly
against both distributions.

### 5. Validate the ion set

Not addressed in the report. §8 notes ions are converted separately and are absent from the
protein xml, but any solvated production system will contain neutralising Na⁺/Cl⁻, so those
parameters come from a file that has not been validated. Check them on the same footing:
per-atom charge and σ/ε against the Amber ion set, and confirm the ion xml declares its types
in a `<LennardJonesForce>` section with a matching `lj14scale`. Because the protein xml carries
LJ in `<LennardJonesForce>` for the NBFixPair, a type missing from that table fails at
`createSystem` — loud, but only once a system containing ions is actually built, which no
fixture so far does. Add an ion-containing fixture.

### 6. Add regression guards for the fault that was found

Two standing assertions, both cheap:

- For every residue, the candidate's torsion term set equals the reference's as an atom-keyed multiset, with zero-force terms canonicalised identically on both sides rather than dropped from either.
- Every specific k = 0 override present in `frcmod.ff99SBdisp` has a surviving counterpart in the generated ffxml. This guards directly against `fix_zero_periodicity.py` regressing, and tests the cause rather than the symptom.

### 7. Add a negative control for this failure mode

The six existing controls do not include it. Add: delete one specific k = 0 override from the
ffxml and confirm a named test fails. This is now the known failure mode and the suite should
be proven against it, not merely patched for it.

### 8. Unit-test the periodicity-0 post-processor

`fix_zero_periodicity.py` changed behaviour and is load-bearing. Give it direct tests on three
inputs: a torsion whose only term is periodicity-0 with k = 0; one whose only term is
periodicity-0 with k ≠ 0; and a multi-term torsion with a trailing periodicity-0 term. Assert
the output in each case.

### 9. Make the coverage claim programmatic

§8 states 76 templates and lists four absent residues. Derive the inventory from the leaprc
source and compare against the ffxml programmatically rather than by count, and confirm that a
structure containing CME, HYP, ORN or DAB raises rather than matching another template by name
collision.

### 10. Add a boundary-spanning disulfide fixture

`disulfide_CYX` passes and `cap_CYX` was correctly dropped as unphysical, but a disulfide whose
S–S crosses a chain boundary exercises the exclusion scheme differently. One fixture.

## What is not required

No audit of the published a99SB-disp parameters. It is validated against its own benchmark of
21 systems and 9000+ experimental data points, and auditing it is not this port's job.

No folded-protein stability run. Ubiquitin and GB3 remain folded under force fields with
substantial backbone-balance errors, so a short run cannot discriminate a faulty candidate
from a correct one, and a run long enough to discriminate costs more than the information is
worth. The constraint and virtual-site comparisons in Tier 1 test the same failure modes
statically and more sharply.

No cross-engine absolute energy matching beyond what is already done. Matching cutoffs, Ewald
coefficients and dispersion-correction conventions between engines tests the settings mapping
rather than the port, and the parameter-level comparison already settles the question it would
be asked to answer.

## Sequencing

Items 1 through 5 gate the claim that the port is validated. Items 6 through 8 gate the claim
that the validation suite is trustworthy going forward, and matter more than their size
suggests, because they are what prevents a silent regression the next time the conversion is
regenerated. Items 9 and 10 are completeness.

None of this should begin on the REST2 layer. The λ = 1 identity and λ = 0.5 halving checks
compare a scaled system against the candidate's own unscaled system, so they pass regardless
of whether the candidate is a faithful a99SB-disp — a conversion error and a scaling error
arriving together would be indistinguishable. Close the force field first.
