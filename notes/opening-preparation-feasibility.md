# Preparation resources and the creation boundary

S079 early feasibility inventory, 2026-10-13. The selected S072 route is a
conditional opening problem on a supplied two-ended bridge. Preparation of
its physical joint state and creation from disconnected initial geometry
remain separate unsolved obligations.

## Inputs the opening calculation currently assumes

| Input | Existing scope | Preparation obligation |
| --- | --- | --- |
| Geometry | Near-extremal magnetic pair with an existing bridge; a product throat is used in intermediate calculations | Specify regular initial data and how the pair/bridge is obtained; include the exterior and its boundary conditions |
| Quantum state | Thermofield-double-like initial pair in the opening benchmark; Gaussian/Friedrichs states in separate response estimates | Match a single physical interacting joint state and its correlations; none of those reference states automatically equals another |
| Matter and charges | The selected family has one four-dimensional fermion species, magnetic flux q and nontrivial charge/spin branching | Account for charge conservation, state accessibility and the joint sector weights |
| Boundary access | An explicit controllable interaction joins two AdS boundary systems | Specify the apparatus/Hamiltonian, causal control resources and switching work |
| Clocks | Both endpoints use a common matched boundary time with positive redshift | Implement the waveform on the actual evolving clocks |
| Preparation tolerances | Free-JT redshift and velocity intervals; separate quantum marginal bounds | Obtain the physical joint norm and tolerances from S075 and S077 |
| Thermodynamics | Temperature, extremal entropy and the supporting charge ensemble are supplied model parameters | Include cooling/entropy removal, energy exchange and duration in the preparation protocol |
| Physical completion | Finite local couplings, boundary terms and state-dependent stress are unmatched | Carry these through S074 to S076 before asserting an operational initial state |

The parameter family remains eta=k^-1, g=k^-4, L4/l_p=k^5,
q=k^4, N_f=1 for integer k>=2000, as recorded by the planning story.
This inventory neither supplies those resources nor demonstrates that they
are inconsistent.

## Primary preparation comparisons

[Maldacena and Qi, section 4.3](https://arxiv.org/html/1804.00491v3#S4.SS3)
starts its quench from the coupled ground state. It relates the subsequent
low-energy solution to a thermofield double. Using that argument to prepare
the coupled ground state in the first place would leave an assumed resource
unresolved. Its low-energy identification also does not determine our
four-dimensional charge/spin state and boundary apparatus.

[Cottrell et al., section 2](https://arxiv.org/html/1811.11528v2#S2)
construct positive parent Hamiltonians that annihilate the exact TFD and
prove uniqueness in finite dimension with sufficient operator generators.
The construction depends on the original Hamiltonian and imaginary-time
dressing of operators. Their discussion of nonlocality and section 7.3
make locality and the gravitational interpretation substantive issues.
For this route we would still have to realize those interactions, control
the physical gap and preparation time, respect gauge sectors, and show
that the prepared state lies in the physical opening basin.
The source is a candidate method, not a four-dimensional preparation result.

The existing [S071 transfer](formation-boundary-flow.md) is parked because
the direct free-channel argument fails. Its missing interacting boundary
kernel and entropy normalization cannot be supplied by importing an SYK
ground-state claim. The [S070 saddle record](ads-wall-ensemble.md)
must likewise retain its unresolved physical gravitational/matter response.
No previously closed or parked branch is reopened by this inventory.

## Topology and the precise scope of creation

The existing differential equations evolve fields on a fixed manifold.
Within that formulation, smoothly changing field values does not change the
manifold's connected components: the same points and topology are the domain
at every step. Hence its opening solution cannot also be a demonstration
that initially disconnected components become connected.
This is a restriction of the specified calculation, not a theorem excluding
all quantum creation processes.

[Friedman, Schleich and Witt, Theorem 1](https://arxiv.org/pdf/gr-qc/9305017v2)
assumes asymptotic flatness, global hyperbolicity and averaged null energy.
The current AdS boundary-coupled negative-energy setup has not been shown
to meet those premises. That theorem cannot be substituted for either a
creation mechanism or a universal impossibility result here.

Any proposal to extend this route to disconnected initial geometry must
identify new dynamics or a sum over topology sectors, its physical boundary
and state meaning, and a controlled transition calculation. Merely preparing
a TFD in a finite quantum simulator does not establish those gravitational
steps.

## Feasibility decision and next work

Continue the conditional opening chain S074 to S078, preserving the initial
bridge and joint-state requirements in every result. S079.1's early resource
inventory is complete. S079's preparation and creation gates remain open.

The first useful preparation calculation after S075 is to compare a specific
realizable preparation Hamiltonian with that story's full physical state,
derive the gap/energy/time and joint-error bounds, and test them against
S077's accepted tolerance. If such a Hamiltonian cannot be matched, record
the failed route and assess a distinct preparation mechanism under a dated
story. No formation rate, physical realization or shortcut follows from
the present comparison.

Remaining to-do list: match the target joint state and tolerances; derive a
controlled preparation protocol and its resources; resolve the separate
topology/creation mechanism before claiming the broad objective complete.
