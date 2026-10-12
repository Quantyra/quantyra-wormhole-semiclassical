# Finite-time fermion state under the prescribed moving load

2026-10-11, E003/S053. The bounded one-particle source in
[the moving-reflection calculation](moving-reflection-load.md) did not
control the filled sea. This note supplies that missing estimate in the
defined free LLL return model, for arbitrary noncommuting classical matrix
sources on the common inner interval. It constructs the causal quasifree
state and bounds its total free excitation energy without a fermion cutoff
or a gap uniform in flux. Unequal returns are retained throughout.

For an inner rotational source of amplitude `epsilon=Omega ell` and
bounded total time variation, the energy is `O(q^3 epsilon^2)` when
`q epsilon log L_max` is small. In the retained family this is
`O(log^6 q)=o(q delta)`. The estimate is independent of how long a static
load is held. It does not claim adiabatic following of a ground state.
It also does not bound local stress, construct the fluctuating gauge-field
state, or match the complete physical exterior. The proof uses established
CAR, energy and external-field methods; no consequential novelty is accepted.

## Specified model and source norms

Measure lengths in units of ell and frequencies in units of ell^-1.
Take q free complex massless Dirac channels on circles of lengths
`L_i=T+D_i`, with a shared coordinate interval `I=[-T/2,T/2]` and D_i>0.
Each circle has the neutral antiperiodic reference vacuum. The right-mover
one-particle energies are `2 pi (n+1/2)/L_i`, n integer; left movers are
independent copies with reversed spatial momentum. Let P_- be the filled
negative-energy projection, and H_0 the normal-ordered free many-body
Hamiltonian, with vacuum energy set to zero. Thus H_0>=0 and measures
particle plus hole energy, not the background Casimir energy.

A prescribed Hermitian q by q matrix potential A_sigma(t,s) acts on I
for each chirality sigma. Outside I it is zero. Assume bounded spatial
multiplication, continuously differentiable time dependence in both
L-infinity operator norm and L2 Hilbert-Schmidt norm, and finite total
variation in those norms on the time interval considered. C2 time
profiles suffice for all the finite-regulator manipulations below.
No spatial derivative bound is needed. In particular the estimate does
not grow with the gradient of a narrow spatial cutoff. The source is
classical; it is not the dynamical electromagnetic or rotational field.

For fixed T put

```text
C_R=1+T/(2pi)[2+log(1+L_max/pi)],
m_T=max(2/pi^2,(2/pi)(1/(2pi)+1/T)),
C_pair=4pi m_T(1+T^2/4),
b(A)=C_R max_sigma ||A_sigma||_(L-infinity,op),
a(A)=C_pair sum_sigma integral_I Tr(A_sigma^2) ds.             (1)
```

The choice of energy one in the first constant is the stated reference
unit. All source norms are evaluated at the displayed time. Take any
constant `a_*>=sup_t a(A(t))`, and assume

```text
sup_t b(A(t)) <= 1/4,
B(t)=integral_0^t b(dot A(r)) dr,
W(t)=integral_0^t sqrt(a(dot A(r))) dr.                        (2)
```

Starting from P_- at time zero, the exact one-particle evolution U(t,0)
of h_0+A gives the pure quasifree covariance
`Gamma(t)=U(t,0) P_- U(t,0)^dagger`. The excitation energy satisfies

```text
0 <= E(t) <= 2 exp(2B(t)) [2sqrt(a_*)+sqrt(2)W(t)]^2.          (3)
```

This is an all-orders statement in the prescribed external source under
(2), not a quadratic current-influence truncation. No commutation of A
at different positions or times, or with the exterior return, is assumed.
The source may be left on after a finite switch: subsequent constant
evolution adds nothing to B or W. An additive source-dependent c-number
changes only the state phase and leaves E(t) unchanged. Its variation
can still matter for physical renormalized stress.

## A fermionic energy inequality

First use any finite symmetric mode cutoff. Write positive energies as
e_p and hole energies as e_h. Every fermionic state, including a state
that is not quasifree, has one-body covariance `0<=Gamma<=1`. Its energy is

```text
E=sum_p e_p Gamma_pp + sum_h e_h(1-Gamma_hh).
```

For a Hermitian one-particle matrix A define its exact finite quantities

```text
a_exact=sum_ph |A_ph|^2/(e_p+e_h),
b_exact=max(||h_+^(-1/2) A_++ h_+^(-1/2)||,
            ||h_-^(-1/2) A_-- h_-^(-1/2)||).
```

Here a_exact includes both chiralities. Diagonal expectations are bounded
by b_exact E. The contraction inequalities Gamma^2<=Gamma and
(1-Gamma)^2<=1-Gamma give

```text
sum_ph (e_p+e_h)|Gamma_ph|^2 <= E.
```

Indeed the p row sum is bounded by Gamma_pp and the h column sum by
1-Gamma_hh. Weighted Cauchy-Schwarz then bounds the two cross terms.
For `V_A=:psi^dagger A psi:` one obtains

```text
|<V_A>| <= b_exact E+2sqrt(a_exact E).                        (4)
```

This proof uses covariance positivity, not Wick factorization or a bound
on the norm of V_A acting on the vacuum. The latter norm can diverge for
a spatially sharp source even when a_exact is finite.

Both finite quantities have the cutoff-independent bounds in (1).
For the diagonal one let R restrict a positive-energy circle to I.
In `R h_+^-1 R^dagger`, energies above one contribute operator norm at
most one. The remaining positive finite-rank operator has trace at most

```text
T/(2pi) sum_(n>=0, n+1/2<=L_i/(2pi)) 1/(n+1/2)
 <= T/(2pi)[2+log(1+L_i/pi)].
```

The direct sum over channels takes a maximum, not a sum. The same
calculation holds for holes. Since `-||A||op <= A <= ||A||op` pointwise,
this proves b_exact<=b(A), including angular channel mixing.

For the pair term, the p,h matrix element is the spatial Fourier
transform of A_ij at energy `E_ph=e_p+e_h`, divided by sqrt(L_i L_j).
The pair measure

```text
mu_ij=(1/(L_i L_j)) sum_(n,m>=0) delta_(2pi[(n+1/2)/L_i+(m+1/2)/L_j])
```

has `mu_ij((0,R])<=R^2/pi^2`. The measure nu_ij=dmu_ij/E has mass at
most m_T in every nonnegative unit bin, as shown by counting the two
half-odd grids in [the free-current bound](propagating-photon-polarization.md).
For F(E)=int_I exp(-iEs) f(s)ds, the unit-bin Sobolev estimate and
Plancherel yield

```text
integral |F(E)|^2 dnu_ij(E)
 <= 2m_T integral_R (|F|^2+|F'|^2) dE
 <= C_pair ||f||_2^2.                                      (5)
```

In the last step F' is the Fourier transform of -isf and |s|<=T/2.
Summing over matrix entries and chiralities proves a_exact<=a(A).
The denominators in (5) are essential: a bound on the bare current
vacuum column would be a different and stronger assertion.

## Work estimate and removal of the cutoff

At a finite regulator set `F(t)=<H_0+V_A(t)>` and `Y(t)=F(t)+4a_*`.
Normal ordering is always relative to the original P_-. Exact unitary
evolution gives `F'(t)=<V_(dot A)(t)>`. Equations (2),(4) imply

```text
Y >= (3/4)E-2sqrt(a_* E)+4a_* >= E/2,
Y' <= 2 b(dot A) Y+2sqrt(2a(dot A)) sqrt(Y).
```

The last step in the first line is completion of the square. The initial
vacuum has F(0)=0, whether the initial classical source is zero or not.
For Y>0, dividing the second line by 2sqrt(Y), then integrating, gives
`sqrt(Y(t))<=exp(B(t))[2sqrt(a_*)+sqrt(2)W(t)]`. Adding an arbitrarily
small positive constant handles zero Y. This proves (3) uniformly in
the regulator. No spectral gap enters this estimate.

To construct the unregulated state, use bounded one-particle sources
`A_N=P_N A P_N` on the full one-particle space, where P_N are increasing
free spectral cutoffs, and keep the exterior modes in their vacuum.
The norm bound on A_N is uniform. Strong convergence of A_N to A,
the bounded-perturbation Duhamel formula and dominated convergence give
strong convergence of U_N and U_N^dagger to U, uniformly on compact
time intervals. Hence Gamma_N converges strongly to Gamma.

Excitation energy is the sum of two positive traces:

```text
E=Tr(h_+^(1/2) Gamma_++ h_+^(1/2))
 +Tr(h_-^(1/2) (1-Gamma)_-- h_-^(1/2)).                     (6)
```

Testing finite spectral subspaces and taking their supremum proves lower
semicontinuity, so (3) passes to (6). This constructs a pure quasifree
CAR state with finite excitation energy. For each fixed finite graph the
NS gap pi/L_max is positive; it then also has finite particle-plus-hole
number and lies in the reference Fock representation. This last
fixed-system statement is not a gap uniform in q. No four-dimensional
Hadamard or interacting-state conclusion follows from finite energy alone.

## Consequence for the moving inner load

For rotational matrices J_a in spin j=(q-1)/2, write a classical source
`A_sigma=epsilon sum_a f_sigma,a(t,s) J_a`. Then exactly

```text
||sum_a f_a J_a||op=j |f|,
Tr[(sum_a f_a J_a)^2]=T_q |f|^2,   T_q=q(q^2-1)/12.
```

Assume the profile sup norms, L2 norms and total time variations are
bounded independently of q. This includes a fixed smooth preparation
followed by a held source and bounded rotations of its axis. Equation
(3) becomes

```text
b_* , B = O(q epsilon [1+log L_max]),
E = O(q^3 epsilon^2) when b_*<=1/4 and B stays bounded.       (7)
```

An electric identity source can be retained in the same proof: its
Hilbert-Schmidt contribution is q times its amplitude squared and its
operator norm is its amplitude. No source-current commutant is removed.

The preceding moving-load note supplies the inner amplitude scale and
the retained family

```text
delta=(log q)^-2,  epsilon=Omega ell=(log q)^3/(2sqrt(pi)q^(3/2)),
L_max=O(1+q delta).
```

For prescribed inner profiles satisfying the stated variation hypothesis,

```text
b_*,B=O((log q)^4/sqrt(q))->0,
E=O((log q)^6),       E/(q delta)=O((log q)^8/q)->0.           (8)
```

The bound holds in particular through the proposed time ell/delta if
the total source variation remains bounded. It improves the previous
one-particle Duhamel comparison by controlling energy in the filled sea.
It does not make the evolving state close to the vacuum in trace norm:
many soft excitations can have small total energy. It also does not
control a local or sharply differentiated stress observable. The physical
exterior may source additional fields outside I, and its time-dependent
matching is not supplied by this conditional inner calculation.

## Verification, predecessors and decision

[The checker](../checks/driven_fermion_sea.py) tests (4) with general finite
Fock vectors, Slater states and mixed covariances; compares the continuum
form bounds to finite unequal-return matrices; and checks the spatial
Fourier entries by independent quadrature. Noncommuting J_z and J_x
pulses are evolved at three time and three spectral resolutions. Work
integration is compared with excitation energy, and a separate adaptive
ODE solver checks selected evolutions. Both chiralities have equal energy
in this test family by an explicitly checked fixed channel symmetry.
[The receipt](../evidence/driven-fermion-sea.json) pins all inputs.

The numerical pulse uses `h1(u)=sin^4(pi u)` and
`h2(u)=h1(u) sin(2pi u)`, with compact time support. Their variations are
2 and `4 max|h2|`, independent of the pulse duration. Spatial profiles
are cos^2(pi s/T) and its product with sin(2pi s/T). Initial development
checks exposed cancellation error from subtracting filled-sea energies;
the final calculation uses positive particle/hole transition energies.
It does not infer convergence from preserving the one-particle norm.

The explicit bound is conservative. For this test profile with the
retained asymptotic amplitude it exceeds q delta at q=10^12, even though
the smallness condition passes; it is about 0.00774 q delta at q=10^16.
These are checks of (3), not useful physical flux thresholds or numerical
solutions of an actual wormhole at those values.

[Marcelli, arXiv:1904.11300v3](https://arxiv.org/html/1904.11300v3),
Assumptions 2.1-2.2, Theorem 2.5 and its section 3.1 proof, is an
established comparator for energy estimates uniform in switching time.
Its operator/commutator hypotheses are not silently assigned to the
current form here. The finite-regulator proof above is given separately.
[Brunetti, Fredenhagen and Pinamonti, CMP 406, 292 (2025)](https://link.springer.com/article/10.1007/s00220-025-05466-4),
section 2.2, reviews external-field CAR evolution and proves an
implementability estimate for smooth compact massive four-dimensional
sources. Its mass-dependent bound is not used for these massless channels.
The [review record](../evidence/driven-fermion-sea-source-review.json)
states exact reading and retrieval scope. These methods and the general
existence of external-field states are prior work.

Continue the physical long-return target and its matched-pulse fallback.
The prescribed inner classical drive no longer requires a uniform free
gap to obtain this energy bound. Quantizing the rotational and electric
fields in the same state, retaining their exterior determinant and
quadratic vertices, and controlling renormalized stress remain the hard
physical step. No new physical endpoint or novelty is accepted. No
Lean/Lake, GCP, outreach, publication or DOI action occurred.

Remaining to-do list: derive the coupled fluctuating gauge/fermion state
and full source-dependent exterior response; bound renormalized stress
and omitted modes; complete smooth matching, the physical endpoint and
novelty review; independently verify, GCP-certify suitable mathematics
and prepare the reproducible manuscript.
