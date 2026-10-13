# Compact transfer on the fermion Hilbert space

S072, 2026-10-13. **The compact thermal transfer bound survives actual
fermionic charge/spin branching in a quadratic-charge Hamiltonian.**
It does not require replacing those branches by the Hilbert space of
an independent group rotor. An exact anticommutation sum rule also
identifies the additional correlation needed when other charged degrees
of freedom or a radial Hamiltonian are present. Those additions remain
essential for the physical initial state and driven continuum stress.

## Why the rotor weights cannot simply be imported

[Iliesiu, section 2.2](https://arxiv.org/html/1909.05253v2) derives the
group-particle partition sum from its action and Hilbert space `L^2(G)`.
The two representation indices give its dimension-squared multiplicity.
General invariant Hamiltonians instead give general functions of the
Casimirs. This fixes the pure gauge comparison, not the determinant of
light charged matter in our fixed-source exterior.

[Iliesiu and Turiaci, sections 3.1, 3.4 and 5.1](https://arxiv.org/html/2003.02860v3)
relate the compact action to near-extremal thermodynamics and keep track
of the representation multiplicity. Their discussion leaves charged
matter contributions as a further problem. We therefore do not infer
the required state merely by inserting the screened inertias into their
pure gauge partition function. The earlier angular source-ensemble
distinction also remains in force.

The result below takes a different, explicit route: specify a charge
Hamiltonian on canonical fermions and keep all its existing branches.
It is not a claim that this is already the full black-hole Hamiltonian.

## Exact CAR sum rule

Let there be `D` canonical complex fermion modes. Their one-particle
electric charges are `e_A=+1` or `-1`; their spin matrices are copies of
spin `j` or its dual, with `sum_a t_a^2=C_j 1`, `C_j=j(j+1)`.
Use centered electric charge and traceless spin generators,

```
Q = sum_A e_A (c_A^dagger c_A-1/2),
J_a = sum_AB c_A^dagger (t_a)_AB c_B,
K = x Q^2 + a sum_b J_b^2,       x,a>0,
lambda = x+a C_j.                                         (1)
```

Equal particle and antiparticle multiplets make the centering constant
zero and preserve an integer electric lattice with neutral spin-zero
vacuum. This is the implementation below: `D=2 M q`, `q=2j+1`, with
`M` orbital copies. The necessary center relation
`2ell=(q-1)Q mod 2` follows directly. It is not used as a substitute for
the actual Fock-space representation multiplicities.

Introduce `2D` Majoranas `Gamma_(2A)=c_A+c_A^dagger` and
`Gamma_(2A+1)=i(c_A^dagger-c_A)`. Each is Hermitian and squares to one.
The unital, trace-preserving, self-adjoint channel

```
E(O) = (1/(2D)) sum_r Gamma_r O Gamma_r
     = (1/D) sum_A (c_A^dagger O c_A+c_A O c_A^dagger)       (2)
```

acts on an even Majorana monomial of degree `p` by `1-p/D`: exactly
`p` of the conjugations anticommute relative to the others. A square
of a centered Hermitian quadratic fermion generator contains only a
scalar and degree-four terms. Its scalar is its infinite-temperature
trace. Wick contraction at occupation `1/2` gives

```
tr_0 Q^2=D/4,       tr_0 sum_b J_b^2=D C_j/4,
tr_0 K=D lambda/4.
```

Consequently the following is an operator identity, before taking a
thermal expectation or choosing a representation basis:

```
E(K)-K = lambda 1-4K/D.                                 (3)
```

This incorporates the exclusion principle and all charge/spin branches.
For example, at `q=4`, the maximal charge `Q=4` in one particle and one
antiparticle multiplet is a unique spin-zero state. The center condition
alone also permits `ell=1` there, but that branch is absent. The checker
extracts spin multiplicities from `Q,J_z` weight counts and checks them
against the Casimir spectrum, without using a group-particle measure.

## Thermal transfer and free evolution

Set `rho=exp(-K)/Z` on this Fock space. The normalized, direction-averaged
midpoint transfer is

```
C(tau)=(1/(2D)) sum_r Tr[sqrt(rho) Gamma_r(tau)
                                      sqrt(rho) Gamma_r],
Gamma_r(tau)=exp(i K tau) Gamma_r exp(-i K tau).
```

At zero time this is `2/D sum_A Tr[sqrt(rho)c_A^dagger sqrt(rho)c_A]`.
In the conjugate-side fermionic TFD, with its phases fixed as below, it
equals the expectation of the charge-conserving interboundary hopping
divided by `D`. It is one for the maximally mixed single-side state.
Its normalization is for these canonical modes, not an unsmeared
continuum boundary field.

In an energy eigenbasis put
`W_mn=(1/(2D))sum_r |(Gamma_r)_mn|^2`. It is symmetric and doubly
stochastic. The same scalar thermal-affinity argument used in
[quantum-junction-state.md](quantum-junction-state.md) now gives

```
S = sum_mn (p_n-p_m)(E_m-E_n) W_mn
  = 2 Tr[rho(E(K)-K)] = 2lambda-8<K>/D >=0,
0 <= 1-C(0) <= S/8 = lambda/4-<K>/D <= lambda/4,
M2 = sum_mn sqrt(p_m p_n)(E_m-E_n)^2 W_mn <= S,
Re C(tau) >= 1-S/8-S tau^2/2
           >= 1-lambda/4-lambda tau^2.                   (4)
```

For the first inequality use
`(sqrt(p)-sqrt(r))^2 <= (p-r)(log p-log r)/4`, and the row/column
sums of `W`. For the moment bound use `2sinh(v/2)>=v` for `v>=0`;
the last line follows from `cos(v)>=1-v^2/2`. Every sum is finite.
Equation (3) also gives `<K><=D lambda/4`; it does not give a
dimension-independent thermal energy bound.

The bound (4) is uniform in `D` and spin when `lambda` is fixed. It
does not require a small-D limit or extrapolation of numerical data.
For the candidate static coefficients on `1/4<=z0<=7/25`, their
previous electric/angular bounds give `x<1/(48 k z0)` and
`a C_j<1/(16 k z0)`. The latter follows after the flux substitution from
`(1+6X)(X+5/6)>2X(1+3X)`. Hence `lambda<1/(12 k z0)<=1/(3 k)`.
Thus (4) gives the initial loss bound `1/(12 k)`, below
`4.167e-5` at `k=2000`, **if this charge Hamiltonian and canonical-mode
transfer are the prepared sector**. This does not update the physical
vertex amplitude or assert a post-switch bound. A sharper family bound
can be used after its coefficients and operator are matched.

## Exact remaining terms for a physical extension

Suppose the same canonical modes have additional charges in a commuting
environment, `Q_tot=Q+Q_e`, `J_tot=J+J_e`. Let
`K_tot=x Q_tot^2+a J_tot^2`. Applying (2) only to the fermions gives

```
E(K_tot)-K_tot
 = lambda 1-(4/D)[x Q_tot Q+a J_tot dot J].               (5)
```

The charge/spin products on the right are not automatically positive.
For an additional Hamiltonian `R`, the Gibbs state of `H=K_tot+R`, at
a finite regulator or where the displayed thermal traces converge,
still has the exact affinity identity

```
S_H=2lambda-(8/D)<x Q_tot Q+a J_tot dot J>
                   +2<E(R)-R>.                         (6)
```

The inequalities `1-C_H(0)<=S_H/8` and `M2_H<=S_H` remain true, but
`S_H<=2lambda` no longer follows. This is a concrete correlation and
commutator target for the radial, coset, horizon and gauge dressing
calculation. Nonuniformly weighted boundary wave packets likewise use
their own channel; averaging all canonical modes is not an identity
for an arbitrary boundary insertion.

An explicit spin-coupling check illustrates the sign issue: a fermion
spin `1/2` coupled to an environment spin `1` in total spin `1/2`
has `J_tot dot J=-1/4`. Center parity alone cannot remove that possibility.
No physical black-hole occupation of this example is asserted.

## A driven check that keeps the joint state

For a finite diagnostic put the same modes on both sides, use

```
H(t)=K_L+K_R-A(t) X,
X=sum_A i(c_L,A^dagger c_R,A-c_R,A^dagger c_L,A).
```

Start from the exact charge-conjugate TFD
`|TFD> proportional exp(-K_L/2) product_A(c_L,A^dagger-i c_R,A^dagger)|0>`.
The product order is fixed. The reduced left density matrix is `rho`;
the total electric charge and total spin annihilate the state. The
Hamiltonian preserves both constraints, but `K_L+K_R` does not commute
with `X`. Evolution therefore includes the joint quantum backreaction.

The checker evolves a smooth monotone drive, checks the initial transfer
against the thermal spectral formula, and independently integrates
`d<H>/dt=-A'(t)<X>`. It retains the full allowed number sector and all
spin branches; no projection onto a selected spin or a reset thermal
state is used during evolution. Negative final energies in this finite
model are not gravitational capture energies. It contains no spatial
stress, radial bath, continuum boundary subtraction or JT geometry.

In the `q=3, x=0.2, a=0.05, A_f=1` diagnostic, normalized transfer
starts at `0.983577` and reaches a sampled minimum `0.728695` during
the drive. These are diagnostic parameters, not the large-flux candidate.
They exhibit why an initial thermal affinity must not be held fixed in
the evolving state. The negative final finite-model energy does not
eliminate this change or supply a continuum negative-stress theorem.

## Decision

The independent-rotor multiplicity assumption is unnecessary for (4)
in the stated fermionic charge model. Retain its proof as a possible
state-matching route, and use (5)--(6) to test the actual physical
extension. Do not infer the missing mixed terms from the small electric
midpoint normalization found in the previous turn. That calculation
controls a different, explicitly stated Gaussian comparison.

Evidence: [checker](../checks/fermionic_compact_transfer.py),
[receipt](../evidence/fermionic-compact-transfer.json), and
[source/dependency record](../evidence/fermionic-compact-source-review.json).
The initial driven solver missed the independent work-balance tolerance;
tightening its tolerance alone did not resolve every case. Splitting the
integration at the known smooth-piece boundaries resolved the discrepancy
without changing the Hamiltonian, drive or acceptance threshold. Both
failed source/log snapshots are retained. A smaller-step, tighter-tolerance
repeat checks the final accepted calculation.
No novelty acceptance or Lean certification is claimed.

Remaining to-do list: match the physical radial boundary transfer and
bound the mixed terms in (6); control the interacting continuum state,
stress and geometry below the capture and area margins; establish
preparation, quantum fluctuations, sufficient stability and finite-energy
transport.
