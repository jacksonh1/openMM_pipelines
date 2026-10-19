# REST2 protocol for CHARMM36m and a99SB-disp

Scientific specification, engine-independent where possible. λ = β_m/β_0 = T_0/T_m;
λ = 1 is the unmodified reference replica and the only one analysed.

## What REST2 scales

This is settled convention, not a choice you need to make. Wang, Friesner & Berne (2011)
define the modified potential as

E_m = λ·E_pp + √λ·E_pw + E_ww

for solute–solute, solute–solvent and solvent–solvent contributions, and realise it by
giving each solute atom a **charge factor √λ** and a **Lennard-Jones ε factor λ** (σ
untouched). Combination rules then produce λ for solute–solute pairs and √λ for
solute–solvent pairs automatically. Bonded terms follow by counting: **λ** if all of a
term's atoms are solute, **√λ** if only some are. 1-4 terms inherit the per-atom factors.

**Bonds and angles are not scaled.** The 2011 paper reports that scaling them does not
help sampling, and PLUMED, NAMD and Amber's REAF all default to leaving them alone (Amber
never scales bonds at any setting). The reason is the premise of solute tempering: exchange
acceptance is governed by the fluctuation of the scaled energy, so adding thousands of
stiff coordinates that do not gate conformational transitions costs replicas for nothing.

So the scaled set is: electrostatics, Lennard-Jones, 1-4 electrostatics, 1-4 LJ, proper
dihedrals, impropers, and CMAP where the engine supports it.

## Engine

Use **GROMACS patched with PLUMED**, with `plumed partial_tempering` writing the scaled
topologies and `gmx mdrun -replex N -hrex` running the exchanges. For both force fields
this is the directly precedented route:

- a99SB-disp's first author publishes exactly this workflow as a REST2 tutorial (`paulrobustelli/IDP_REST_tutorial`: GROMACS 2022.5 + PLUMED 2.9.1, a GROMACS-format `a99SBdisp.ff`, `plumed partial_tempering` on a `grompp -pp` topology, `gmx mdrun -multidir -replex 800 -plumed plumed.dat -hrex`, with Zenodo-archived reference data).
- The published REST study of a99SB-disp IDP ensembles (Zhang, Liu & Chen 2023) used the same combination.
- Bussi's implementation paper states that, applied to the whole solute, it reproduces REST2 exactly.

`partial_tempering` handles the force-field-specific terms you would otherwise have to
reason about: it rewrites `[nonbond_params]` and `[pairtypes]`, which covers CHARMM36m's
NBFIX off-diagonal table and a99SB-disp's modified backbone carbonyl-O / amide-H LJ pair,
by duplicating atom types for hot atoms. It scales proper and improper dihedrals and the
1-4 `[pairs]`. It does not touch bonds, angles or CMAP.

Running both force fields through one toolchain also means the sampling scheme is held
fixed across them, so differences you observe are attributable to the force field.

## The CMAP question, and why it is not a correctness question

`partial_tempering` leaves CHARMM36m's backbone φ/ψ cross-term unscaled. NAMD's native
`soluteScaling` does scale it (its default covers "dihedral, improper, and crossterm"),
as does GENESIS's gREST. The divergence looks alarming and is not.

Hamiltonian replica exchange samples the correct distribution at λ = 1 for **any** set of
modified Hamiltonians, provided the λ = 1 replica carries the unmodified potential and the
acceptance test uses the actual energies of each configuration under each Hamiltonian.
GROMACS `-hrex` does exactly that — it performs an unconditional coordinate swap,
recomputes the energy with the local force field, swaps back, and uses the stored
cross-energies in a corrected Metropolis test. So leaving CMAP unscaled changes **how fast
you converge, not what you converge to**.

The practical consequence is worth planning around. In CHARMM36m a substantial part of the
backbone conformational preference sits in CMAP; leaving it at full strength in every
replica means the effective backbone tempering at a given λ is weaker than the nominal
effective temperature suggests. Expect to need a wider λ range or more replicas for
CHARMM36m than for a99SB-disp, and set the ladder empirically from measured acceptance
rather than from a target temperature.

If you want CMAP scaled, the standard way is NAMD, not a custom OpenMM build. I would not
introduce a second toolchain for it.

## Ladder

Geometric spacing in λ, targeting roughly 20–40% acceptance between neighbours. Published
anchors, all with exchange attempted every 2 ps:

- Original REST2: 10 replicas, effective solute temperatures 300–572 K, ≈30% acceptance (trpcage, β-hairpin).
- Bussi: 5 replicas, λ = 1 → 0.3 geometric, 35–50% acceptance (alanine dipeptide); 16 replicas over the same λ range for an RNA tetraloop.
- Zhang, Liu & Chen: 16 replicas spanning 298–500 K, or 8 spanning 298–450 K, for IDPs in a99SB-disp; 1–2 µs per replica.

Set the hot region to the protein (all of it, for a folded or designed protein). For a
protein–peptide complex, consider restricting it to the peptide plus interface side
chains — FEP/REST and gREST+FEP both narrow the region this way for binding problems, and
a smaller hot region needs fewer replicas.

## Compaction at high effective temperature: how solid is it?

Weaker than it is usually relayed, including by me. The claim that REST2's √λ
solute–solvent factor drives artificial collapse in the hot replicas comes from a single
paper — [Zhang, Liu & Chen 2023](https://doi.org/10.1021/acs.jctc.2c01139) — on two IDPs
(p53-NTD and KID) in one force field (a99SB-disp), from one group, with no independent
replication among its forward citations.

Critically, that paper ran no temperature-REMD control at matched temperature range; it
invokes T-RE only to argue such a run would need over 100 replicas. The scheme-artifact
hypothesis predicts the collapse should be absent under T-REMD, which applies no
solute–solvent weakening. Where that control has been run at matched range, the collapse
has been reproduced under T-REMD — which is evidence against the REST2-specific mechanism
and for a force-field or physical origin.

The framing that the weakening was an intentional design choice to aid refolding is also
not supported by the primary source. The [2011 REST2 paper](https://doi.org/10.1021/jp204407d)
contains no discussion of compactness, collapse, hydrophobicity or radius of gyration; it
justifies the √λ factor on acceptance-ratio grounds, via an approximate cancellation of
E_pp and the scaled E_pw in the exchange criterion. Its complaint about REST1 was the
reverse problem — hot replicas trapped in *extended* structures.

The force-field explanation carries much more weight. a99SB-ILDN with TIP3P over-compacts
disordered states at 300 K before any heating, which is the stated motivation for TIP4P-D
([Piana et al. 2015](https://doi.org/10.1021/jp508971m)), for rebalanced protein–water
interactions ([Best, Zheng & Mittal 2014](https://doi.org/10.1021/ct500569b)) and for
a99SB-disp ([Robustelli, Piana & Shaw 2018](https://doi.org/10.1073/pnas.1800690115)).
Heating an under-solvated force field worsens this across roughly 300–400 K, under any
tempering scheme. Expect CHARMM36m to compact more than a99SB-disp for the same reason:
[Gong, Zhang & Chen 2024](https://doi.org/10.1021/acs.jpcb.4c04777) find c36m
systematically over-stabilises residue pairs relative to a99SB-disp — 0.2 kcal/mol
nonpolar, 0.6 polar, 0.8 salt-bridge, up to 2.9 for Arg–Glu — and attribute an
over-compaction bias to it directly. Both of these claims originate from the same group,
so treat the method and force-field explanations as entangled at the source rather than
independently corroborated.

**What follows for practice.** Do not pick a sampling scheme to dodge compaction; pick the
force field deliberately and treat the scheme as fixed. The residual concern that is worth
monitoring is not biased dimensions but **replica segregation**: if the hot replicas become
compact enough that they stop exchanging downward, the ladder splits into two populations
and the λ = 1 replica stops being refreshed. That does not bias the ensemble — it prevents
it from converging, which is more dangerous because the output still looks like data. It
is also free to detect from the exchange record, so make it routine per-system QC rather
than a reason to modify the method.

## Validation, in two tiers

**Tier 1 — the implementation is right. Cheap, run once, before any production.** These are
the tests `partial_tempering`'s own documentation prescribes; its help text says plainly
"It's not very robust and there might be force-field dependent issues! A few tests are
strongly suggested."

1. Scale with λ = 1.0 and `gmx mdrun -rerun` the same trajectory under the scaled and unscaled topologies. Every energy group must match, not just the total.
2. Mark **all** atoms hot, scale with λ = 0.5, and rerun. Electrostatic, LJ and dihedral energies must come out at exactly half; bonds and angles must be bit-identical.

Two reruns per force field. If these pass, the Hamiltonian is what you think it is.

**Tier 2 — the sampling is converged. This is where your expensive compute is at risk.**

1. Demultiplex and plot each replica's walk through λ space. You want every replica visiting both ends repeatedly, not two separated populations. A bimodal occupancy with few crossings is the segregation failure above, and it invalidates the run no matter how long it is.
2. Measure the acceptance ratio between every neighbouring pair, not the average. A single bottleneck pair breaks the ladder.
3. Run two independent simulations from deliberately different starting conformations and require that the λ = 1 observables agree. Zhang, Liu & Chen ran one from a helical and one from an extended structure for exactly this purpose. This is the single most informative check you can make and the hardest to fake.
4. Block-average your reported observables over the second half and confirm the blocks agree.

## Methods paragraph

Something of this shape reads as entirely standard:

> Enhanced sampling was performed with replica exchange with solute tempering (REST2) in
> the solute-scaling formulation of Wang et al. (J. Phys. Chem. B 2011, 115, 9431),
> implemented as Hamiltonian replica exchange in GROMACS <version> patched with PLUMED
> <version> (Bussi, Mol. Phys. 2014, 112, 379). Scaled topologies were generated with
> `plumed partial_tempering`, which scales the partial charges of the solute by √λ, the
> Lennard-Jones ε by λ, and the proper and improper dihedral and 1-4 terms accordingly;
> bond, angle and CMAP terms were left unscaled. The solute region comprised <atoms>. N
> replicas were distributed geometrically over λ = 1 to <λ_min>, corresponding to effective
> solute temperatures of <T0> to <Tmax> K; exchanges were attempted every <n> ps and the
> mean acceptance ratio was <x>%. All replicas were simulated at <T0> K and only the λ = 1
> replica was analysed. Scaling was validated by confirming that per-term energies at λ = 1
> reproduced the unscaled force field and that the scaled terms halved at λ = 0.5.

State the acceptance ratio and the λ = 1-only analysis explicitly. Those are the two things
a referee familiar with REST2 looks for.
