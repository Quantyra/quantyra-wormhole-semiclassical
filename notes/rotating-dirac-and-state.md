# Rotating Dirac reduction and the equilibrium-state obstruction

2026-10-10. S051/S052/S054. **Decision:** retain local rotating scattering
as an ingredient, but retire the proposed identification of the global
corotating equilibrium determinant with the physical support energy as
a model-closure option. That identification was not previously proved.
For the prescribed free Dirac field, a global Hadamard corotating ground
or finite-temperature equilibrium state is incompatible with the spacelike
rotating time direction beyond the light cylinder. This is an application
of established state results, not a new wormhole obstruction. It does not
rule out the metastable, approximately adiabatic MMP construction.

An exact local spinor congruence does control the previously dangerous
linear shift term. The resulting **bulk** estimate is useful, but neither
it nor probability-current contraction chooses the quantum state or
matches the curved throat. The next scientific question is finite-time
support in a specified admissible state, within the existing program.

## Gauge and tetrad conventions

Work first on the flat exterior, in exact rigidly rotating coordinates:

```text
ds^2=-dt^2+|dx+b(x)dt|^2, b=Omega cross x,
A=A_t dt+A_j dx^j, D_mu=partial_mu-i A_mu, Pi_j=-i D_j.
```

Use the orthonormal coframe `e^0=dt`, `e^j=dx^j+b^j dt` and its rotating
spatial spin frame. A right-handed Weyl equation has Hamiltonian

```text
H_b=sum_j (sigma_j-b_j I) Pi_j + Z,
Z=-Omega.sigma/2-A_t I.                                      (1)
```

The other chirality is treated separately with its corresponding signs.
Equation (1) follows by transforming an inertial spinor with the spatial
rotation and its spin lift: time differentiation supplies `-Omega.J`.
The gauge one-form transforms with `A_t^rot=A_t^inertial+b.A`; hence
`-b.Pi-A_t^rot=-b.p-A_t^inertial`. This checks both the spin term and the
use of the **coordinate** electric potential. Dropping `b.A` while using
`A_t^rot` would create an erroneous flux-enhanced rotation term.
[Ambrus and Winstanley v2](https://arxiv.org/html/1401.6388v2), (9)-(12),
gives the same ungauged rotating-tetrad convention. The screened `A_t`
in [the previous note](rotating-electric-screening.md) is a coordinate
component, not an independently chosen inertial scalar potential.

## Exact local congruence

On any region with `|b|<=v0<1`, set

```text
gamma=(1-|b|^2)^(-1/2),
S=sqrt((gamma+1)/2) I + gamma/sqrt(2(gamma+1)) sigma.b,
T=I_3-gamma/(gamma+1) b b^T, M_j=sigma_j-b_j I.
```

Pauli multiplication, without a large-flux expansion, gives

```text
det S=1, S=S*, S>0, S^2=gamma(I+sigma.b),
S M_j S=sum_i sigma_i T_ij,
S H_b S=sum_ij sigma_i T_ij Pi_j+C,
C=-i sum_j S M_j (partial_j S)+S Z S.                        (2)
```

This is a congruence of differential expressions, not a unitary similarity
or a claim that the quantum vacua coincide. In particular the stationary
frequency term becomes `s S^2`, not `s I`. For rigid rotation,

```text
||S^2-I||<=C_v0 |b|, ||T-I||<=C_v0 |b|^2,
||C||<=C_v0 (|Omega|+|A_t|).                                 (3)
```

The derivative of `S` is bounded by `C_v0 |Omega|`, because the boost
depends smoothly on `b` on this closed subluminal ball. No derivative of
the magnetic potential enters `C`. The apparently large `O(|b| Pi)`
term is absent from (2); its remaining spatial correction is quadratic.
For comparison, expanding `S=I+sigma.b/2+O(b^2)` cancels the linear
principal term by `{sigma_j,sigma.b}=2b_j I`. Differentiating the boost
gives `sigma.curl(b)/2=Omega.sigma`, leaving `Omega.sigma/2` after the
spin term in (1), as the linear expansion of (2) requires.

The transformed current through a surface with unit normal `n` is

```text
(S Phi)* (sigma.n-b.n I) (S Phi)=Phi* sigma.(T n) Phi.         (4)
```

The eigenvalues are `+/-|T n|`, with `gamma^-1<=|T n|<=1`.
Also `C-C*=-i sum_ij sigma_i partial_j T_ij`, so the transformed
expression has the correct Green identity. The positive density is
`Phi* S^2 Phi`. These formulas specify the actual local flux metric;
retaining flat port traces without this transformation is not justified.
Nor does (4) conjugate an arbitrary physical throat into the old spherical
load. Its angular operator, boundary state and fields still need matching.

## Bulk consequence for the joint trial

Use `d=1`, fixed trim and frequency window, and the static trial `Phi`
of [joint mouth matching](joint-mouth-matching.md), supported in a bounded
trimmed tube/collar region. In addition to its proved residual bound,
its explicit profiles give

```text
||Phi||<=C, ||Pi Phi||<=C sqrt(q/p), p=q^(-3/4).              (5)
```

Here `Pi Phi` denotes all three Cartesian covariant derivatives. To see
the stronger norm in (5), a normalized spherical LLL profile has angular
covariant derivative squared at most `C q`; after dividing by mouth
radius `r`, this costs `C q/r^2` per radial section. The deformed profile
has relative field correction `O(r^2)` and the explicit radial corrector
is polynomially smaller in the joint regime. Its differentiated terms
obey the same upper bound. In the tube the Gaussian width is
`O(r_*/sqrt(q))`; normal and covariant azimuthal derivatives have squared
norm at most `C q/r_*^2`, while longitudinal, spin-frame and cutoff
terms are smaller or exponentially suppressed. Integrating gives
`integral_p^C q dr/r^2<=C q/p`. The overlap cutoff adds at most `C/a`,
which is absorbed since `p<a`. Axial symmetry makes this squared
covariant norm diagonal in total azimuthal angular momentum, so the
bound holds for normalized superpositions, without a channel-count factor.

Let `f0=(Xi+i sigma.Pi)Phi`, and put `Psi=S Phi`. Then exactly

```text
S(Xi+i H_b)Psi
 =f0+Xi(S^2-I)Phi+i sum_ij sigma_i(T_ij-delta_ij)Pi_j Phi+i C Phi.
```

Since `S^-1` is uniformly bounded, (2)-(5) imply

```text
|| (Xi+i H_b)Psi ||
 <= C [q^(-3/8)+v+v^2 sqrt(q/p)+|Omega|+||A_t||_infinity],  (6)
```

where `v=sup_support |b|`. Restoring units and using the prescribed
screening bound, the dimensionless terms are
`q^(-3/8), |Omega|d/g^2, (|Omega|d)^2 sqrt(q/p)` at fixed small `g`.
On the established orbital scale `d/r_e=q^nu`, `Omega d=2q^(-nu/2)`,
the last term is `O(q^(7/8-nu))`. Thus **`7/8<nu<1` suffices for this
bulk residual to vanish** while retaining the previous geometric and
orbital hierarchies. For `nu=15/16` this sufficient bulk rate is only
`O(q^(-1/16))`; it is not the old `q^(-3/8)` rate or a finite threshold.

Keeping only a middle-region gradient estimate would miss the shrinking
collar cost `sqrt(q/p)` and incorrectly weaken the condition. Conversely,
acting with the unmodified linear shift on the static trial gives the
insufficient bound `v sqrt(q/p)`. The cancellation in (2) is essential.

Equation (6) is **not a rotating return theorem**. It has no physical
electrostatic boundary data, transformed high-mode load, exterior
radiative state, gravitational matching or full energy remainder. A
bulk residual alone must not be promoted through the old boundary estimate.

## Why global equilibrium does not close the state

Let `chi=partial_t+Omega partial_phi` in inertial coordinates. In flat
space `g(chi,chi)=-1+Omega^2 r_perp^2`; it is spacelike beyond the light
cylinder for every nonzero `Omega`. The same issue occurs in an
asymptotically flat region with a putative complete rigid helical flow.

[Pinamonti, Sanders and Verch v2](https://arxiv.org/html/1806.02124v2),
section 2, proves that a finite-temperature KMS state for a complete
smooth lifted flow is not microlocally admissible where its generator is
spacelike. Its assumptions include a distributional positive two-point
function and the induced field-algebra automorphisms. For a free Dirac
field, Hadamard states meet the admissibility requirement. The proof
uses the KMS strip to restrict singular covectors to `k(chi)>=0`, while
admissibility requires all past null covectors; a spacelike `chi` violates
that inclusion. Section 3 explicitly distinguishes existing CAR states
from physically admissible ones. Mere algebraic existence is insufficient.

The same local argument excludes an admissible ground state for this
flow: its GNS generator is nonnegative, so `exp(-eta H)` is bounded for
every `eta>0` and supplies the required half-strip boundary distribution.
The one-sided wavefront argument is unchanged. The analytic boundary-value
estimate is the established result of
[Strohmaier, Verch and Wollenberg v1](https://arxiv.org/html/math-ph/0202003v1),
Theorem 2.8, with the distributional extension explained in the 2019 proof.
This ground-state specialization is not a new theorem. Smooth magnetic
connections change lower-order coefficients, not the Dirac null cone.
This application assumes a complete symmetry of the prescribed field
problem; an actually evolving radiating binary does not acquire that
symmetry by freezing its orbit in a calculation.

An independent local diagnostic is the massless translation model
`chi=partial_t+v partial_y`, `v>1`. Future modes have
`k_chi=p(1-v mu)`, where `mu=p_y/|p|`. On `mu>1/v`, their corotating
energy is negative at arbitrarily large `p`. A thermal Fermi occupation
there tends to one rather than zero, as does the ground-state limit.
For one particle polarization the leading inertial excitation energy
with momentum cutoff `Lambda` is

```text
integral_{p<Lambda} d^3p/(2 pi)^3 p 1_{mu>1/v}
                  =Lambda^4 (1-1/v)/(16 pi^2).               (7)
```

This diagnostic is not the renormalized wormhole stress tensor. It shows
why filling the negative corotating modes is not a harmless choice of
subtraction relative to an inertial Hadamard state. The general local
state conclusion above comes from the cited theorem, not this cutoff.

The local boost (2) cannot evade the state obstruction: it becomes
unbounded at the light cylinder. A small `Omega d` controls a compact
mouth region, not a uniform expansion on the entire exterior.
Probability-current contraction is a distinct property and can hold
without a Hadamard corotating equilibrium state. Therefore the existing
conditional scattering and static Gaussian energy arguments are preserved,
but their physical identification cannot be supplied by global rigid
equilibrium. No obstruction to finite-time traversability or a theorem
about the actual travel-time advantage follows from this observation.

## Verification, attribution and next question

Run `python checks/rotating_dirac_and_state.py --output
evidence/rotating-dirac-and-state.json`. Exact Pauli/gauge/current and
scale checks are supplemented by direct finite-difference evaluation
of the original differential operator against the congruence formula,
and independent quadrature of the Fermi occupation diagnostic. They do
not test a self-consistent state, a rotating loaded PDE or semiclassical
Einstein evolution. The analytic estimates and imported state theorem
carry their own stated hypotheses.

All 11 exact identities and ten rational scale checks pass. The 24
differential-operator cases compare two finite-difference steps, including
two magnetic strengths and two imaginary frequencies; the largest final
relative discrepancy is `1.77e-8`. Green/current identities agree to
floating precision. Nine independent occupation-integral comparisons
agree within `8.68e-19`, and approach the proved cutoff coefficient.
Those numbers measure these algebra/diagnostic checks, not a physical
support-error bound. The verification was rerun after adding the three
finite-time buffer scale checks; the initial run had already passed.

The local boost is ordinary Lorentz/Pauli algebra; the equilibrium
obstruction and inequivalent rotating fermion vacua are established.
Ambrus-Winstanley sections 3-4 also show why states regular inside the
light cylinder must not be described as globally regular rotating vacua.
The main proof in PSV section 2 and the complete proof of the analytic
boundary-value dependency were read. Textbook analytic/GNS dependencies
and publisher-edition comparisons remain distinct pending review items.
No first rotating-vacuum or wormhole no-go claim is available here.

The next target is the original physical-support question over a finite
time interval: specify an initial Hadamard state, exterior incoming data
and throat electrostatic/holonomy sector, propagate the two-point function,
and compare its renormalized stress and variations with the static support
budget over at least one throat transit. An arbitrary equilibrium mirror
inside the light cylinder would alter the problem and is not selected.

There is a useful finite-time scale window already compatible with MMP.
Write `ell0` for its leading throat optical scale and suppose
`Omega ell0<<1`. Choosing a candidate exterior buffer
`R=sqrt(ell0/|Omega|)` gives both
`ell0/R=sqrt(|Omega| ell0)<<1` and
`|Omega| R=sqrt(|Omega| ell0)<<1`. Thus a time interval of order `ell0`
can fit well inside a timelike buffer; the light-cylinder objection does
not eliminate this finite-time approach. At `nu=15/16` both ratios have
exponent `-13/64`, and `d/R` has exponent `-17/64`, at fixed coupling.
These are necessary scales, not a construction of a physical buffer or
an accepted state-approximation error.

The relevant locality statement is exact for a prescribed globally
hyperbolic Dirac problem: if two admissible initial covariances agree on
`U cross U` on a Cauchy surface, and their coefficients agree on its
domain of dependence, their two-point functions agree in that domain.
Indeed backward evolution of each compactly supported field test in the
domain has Cauchy support in `U`; both covariance pairings are then equal.
This does not assert that arbitrary desired interior covariance data have
a positive globally Hadamard extension. Constructing such a state, and
checking that the actual throat/exterior causal domain has the proposed
buffer, are required parts of the next model-closure test. A geometric
buffer alone does not set the initial occupation numbers or quantum stress.

Primary/fallback novelty selection and the complete energy/observer
comparison remain open. No Lean campaign for this known obstruction is
selected before a distinct complete physical consequence survives review.

Remaining to-do list: close that finite-time quantum state and physical
matching, establish the full energy/transport and observer consequence,
complete primary/fallback novelty selection, then independently verify,
GCP-certify suitable complete endpoints and prepare the manuscript.
