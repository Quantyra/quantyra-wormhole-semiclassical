# Angular gluing, Casimir support and the odd electric load

2026-10-11, E003/S053 continuation. The antipodal identification chosen
in the classical binary calculation cannot be carried into an untwisted
product throat while also assuming the common antiperiodic fermion return.
In the leading free diagonal-channel model it introduces a half-turn,
removes the order-q negative opening source at zero exterior length and
gives the wrong sign of the averaged opening source in the retained
point-dipole scaling. This is a restriction on that specified reference,
not an obstruction to all wormholes or to antipodal coordinate charts.

Use the exterior reflection as the zero-twist reference instead. The
magnetic drive remains even, but the electric/axial drive becomes odd.
The corresponding prescribed classical junction is solved below,
including the rotation datum required by equal fixed outer frames.
Neither the full physical quantum state nor a smooth junction is proved.
The phase cancellation mechanism is established prior work; consequential
novelty is not accepted.

## The invariant distinction between the two reference gluings

Freeze the leading point-dipole exterior at an orbital phase and put the
separation direction along n. The field-line endpoint map in
[the explicit exterior trial](global-mouth-trial.md) is

```text
F = I-2 n n^T;  for n=e_z: (theta,phi) -> (pi-theta,phi).
```

It reverses the sphere orientation, as required to identify the opposite
magnetic charges. The continuous spinor frame used there makes the open
transfer positive in the reflected receiving basis. That convention
alone does not remove the physical phase of a closed loop.

The [classical binary reference](binary-mouth-drive.md) instead chose
O_A=-I. With a product throat carrying no extra angular connection and
the same fixed seam lifts, the relative angular map is

```text
O_A^(-1) F = rotation by pi around n,
U_A = exp(i pi J_n), up to a common phase.                       (1)
```

Reversing the sign in the exponent only changes the common phase at
this angle. The q-dimensional lowest Landau level has j=(q-1)/2 and
m=-j,...,j. Its explicit e^(im phi) modes give (1) directly. The spin
lift and electromagnetic gauge choice can change an overall phase;
they cannot erase the alternating relative phases. In the leading
axisymmetric return, the delay and J_n commute.

Equation (1) specifies the *total* seam/transport combination. A throat
rotational Wilson line can compensate the half-turn. Likewise a change
of angular coordinates also transforms the connection and seam lift.
Thus O_A by itself is not gauge-invariant evidence for an obstruction.
The excluded reference is O_A with zero compensating connection. The
alternative O_0=F with zero additional twist has identity relative angular
return at this order. This is a reference to be matched, not a proof that
an evolving binary selects or preserves it.

## Exact support bound for the half-turn model

Take q free complex massless Dirac channels, optical throat length T,
diagonal exterior delays D_j>=0 and L_j=T+D_j. Specify the common vector
boundary condition psi_j(x+L_j)=-exp(i(phi+theta_j))psi_j(x). This excludes
reflection, escape, mode coupling and interactions. Each theta_j is a
phase of (1). There are n_+=ceil(q/2) phases in one class and
n_-=floor(q/2) in the other, separated by pi.

Writing a_j=wrap(phi+theta_j) in [-pi,pi], the circle energy and the
reference-subtracted averaged throat stress are

```text
E_loop(phi) = sum_j [-pi/(6 L_j)+a_j^2/(2pi L_j)],
rho(phi) = sum_j [-pi/(6 L_j^2)+a_j^2/(2pi L_j^2)] + q pi/(24T^2),
chi(phi) = -8T^2 rho(phi)/(pi q).                                (2)
```

Here rho denotes both Delta T_tt and Delta T_xx in conformal coordinates,
summed over angular channels. The reflecting AdS strip is the reference;
its anomaly-corrected null source is zero. Positive chi is the desired
negative opening source. These are angular averages, not a pointwise
four-dimensional source. The local anomaly is unchanged between states.
The spectral regularization and determinant derivation are retained in
[the return-phase note](lll-return-holonomy.md). At fixed phase,
-partial_T E_loop is the circle part of rho.

For positive weights w_j let W_+,W_- be their sums in the two classes.
On either of the two smooth phase branches the minimum is

```text
min_phi sum_j w_j a_j^2 = pi^2 W_+ W_-/(W_++W_-).                (3)
```

For example, class phases 0,pi are minimized at
phi=-pi W_-/(W_++W_-); the other branch has the same minimum.
This also checks the branch endpoints. Set W1_+/-=sum_+/- 1/L_j,
W1=W1_++W1_-, Delta W1=W1_+-W1_-. Then

```text
min_phi E_loop = -pi W1/24 - pi (Delta W1)^2/(8 W1).             (4)
```

Stress uses different weights: b_j=(T/L_j)^2,
B_+/-=sum_+/- b_j, B=B_++B_-, Delta B=B_+-B_-. It follows that

```text
max_phi chi = B/(3q)-1/3 + (Delta B)^2/(q B).                    (5)
```

The phase that minimizes energy generally does not maximize support.
Equation (5) bounds *every* common phase, including an energy-minimizing
Wilson phase; it does not identify the actual gauge-field ground state.

If all D_j=0, (5) is exactly 0 for even q and 1/q^2 for odd q. The
remaining odd-q negative stress is order 1/q, rather than order q.
It cannot supply the retained leading MMP opening. This conclusion
does not assume an independently adjustable phase on each two-pass loop.

## Monotone point-dipole delays and a uniform finite-q bound

Take midpoint flux labels u_j=(j+1/2)/q, j=0,...,q-1 and D_j=d f(u_j),
with the center-to-center point-dipole arc length f. This is the specified
semiclassical channel discretization, not an exact finite-flux scattering
spectrum. The prolate parametrization used in
[the collective-clock calculation](collective-clock-response.md) gives

```text
f(u)=g(u)/u,       g(u)=integral_0^1 sqrt(F(u,eta)) d eta,
a=1-eta^2,        Q=sqrt(a^2+4u^2 eta^2),    R=(a+Q)/2,
F=R+4eta^2 R^2/Q^2-4eta^2 R u^2/Q^2
 =R[a^2+4eta^2 R]/Q^2.                                        (6)
```

The following sign identity makes monotonicity global, extending the
earlier endpoint derivative argument. Put z=a/Q. In the interior,

```text
u partial_u log F
 = -(1-z)[2z+a^2/(a^2+2eta^2(a+Q))] <= 0.                      (7)
```

It follows by differentiating Q, with u partial_u Q=Q(1-z^2).
Endpoint values follow by continuity for u>0. Since F(1,eta)=1 and
F(0+,eta)=1+3eta^2 almost everywhere, monotonicity and bounded convergence
give

```text
1 <= g(u) <= kappa=1+asinh(sqrt(3))/(2sqrt(3)),
f is strictly decreasing;  1/u <= f(u) <= kappa/u.               (8)
```

One can bound the integrand by sqrt(1+3eta^2) using (7), so no singular
differentiation under the endpoint integral is needed.

Consequently the b_j are increasing, lie in [0,1], and their alternating
sum satisfies |Delta B|<=1. For even q, pair consecutive terms and bound
the sum of their increments by the total increment; for odd q, isolate
one endpoint and pair the remainder. This proves the bound without a
continuum replacement. There are at least floor(q/2)>=q/3 labels at or
above 1/2 when q>=2. With r=d/T, (8) gives

```text
B/q <= (1+r)^(-2),
B >= (q/3)(1+2kappa r)^(-2),
max_phi chi <= -[1-(1+r)^(-2)]/3 + 3(1+2kappa r)^2/q^2.         (9)
```

Thus, whenever the right side is negative, even the most favorable
common phase has positive reference-subtracted averaged null source.
For r=(pi log^2 q)^(-1), the negative term is of order 1/log^2 q and the
positive remainder is order q^-2. The old untwisted-antipodal free
reference therefore fails asymptotically in the retained family.

This is not yet an error bound for finite ports, finite-flux projections,
curved propagation, non-Abelian mixing or an interacting quantum state.
Such corrections must be smaller than the margin in (9) to transfer
the conclusion. The exact zero-delay cancellation alone has no margin
against omitted terms. The reflection reference instead has no relative
phase penalty: chi_0=4B/(3q)-1/3, which tends to 1 as r->0 in this model.
The proof is bounded convergence plus the monotone midpoint sum; it
does not assert a uniform operator approximation to the physical return.

## Reflection gluing changes the electric boundary problem

The leading binary magnetic coefficients are epsilon_+=-chi_b^2 n and
epsilon_-=+chi_b^2 n, where chi_b=P/d. Both O_A and F send n to -n;
their common-chart magnetic coefficients are equal in either reference.
But the electric fields are E_+=+E l and E_-=-E l with orbital normal
l perpendicular to n. Since F l=l, the scalar electric potentials
have opposite common-chart data under F. The axial angular velocity
transforms with det(F)F^T; the physical leading electric-induced rotations
also become opposite. The stationary reference therefore has V and omega
odd through the core. The earlier even zero-charge electric solution is
valid only for its stated antipodal classical boundary problem.

Use the classical product/RN equations of
[the stationary electric note](stationary-electric-dipole.md), with
P=1 in each RN chart, lambda=P/ell, throat ends rho=+/-R,
s=lambda R, r_s=1+s. Upper-mouth RN quantities obey

```text
u=r^2 V_RN'=e u0+b u1-2J/r,
u0=(r-1)^2(r+2)/r,
u1=(2r-1)(5r^2-5r+2)/[30r(r-1)^3],
omega_RN'=(-6J-2u)/r^4.
```

The odd throat solution, with beta=0 by parity, is

```text
H_o=rho(1+rho^2),   g_o=1+3rho^2,
V=k H_o-2J rho,     omega=-2k H_o-2J rho,
u_th=k g_o-2J,      omega'+2V'=-6J.                             (10)
```

The conserved outward charges are J and -J at the two ends. Match the
electric flux and invariant Phi=V-omega as in the earlier auxiliary
problem. At r_s and rho=R this gives

```text
b u1-k g_o = -e u0+2J(1/r_s-1),
b f u1'-lambda(1+R^2)k g_o' = -e f u0'-2J f/r_s^2,
f=(r-1)^2/r^2.                                                 (11)
```

For every lambda,R>0, the coefficient determinant has nonzero magnitude
D_o=lambda(1+R^2)g_o' u1-g_o f u1'>0. At J=0, k=e/D_o;
in the overlap limit s->0,R->infinity this is k~e lambda^2.

Parity removes the common rotating constant, so matching the outer
angular frame supplies an additional equation. For a fixed outer annulus
radius r_out, the upper coefficient is

```text
omega_out=lambda[-2k H_o(R)-2JR]
             +integral_(r_s)^r_out [-6J-2u(r)]/r^4 dr.           (12)
```

The lower common-chart coefficient is its negative. Equal fixed outer
frames require omega_out=0. Solve (11) separately for unit e and unit J
to obtain k_E,b_E and k_J,b_J, then (12) gives omega_out=e w_E+J w_J.
For s->0 and R->infinity at fixed r_out>1,

```text
w_E -> -[1-w0(r_out)],      w0=2/r-(2r-1)/r^4,
w_J -> -[1-wKN(r_out)],     wKN=(2r-1)/r^4,
J/e -> -[1-w0(r_out)]/[1-wKN(r_out)].                           (13)
```

For the electric response, (11) gives k_E~lambda^2 and the leading
b_E cancels; the throat contribution in (12) vanishes like s^3.
For the spin response, the previously checked expansion gives
k_J~(8/15)lambda/R and b_J~-6s^4. Its throat term is order s^2.
The b integrals are controlled by integral u1/r^4=O(s^-2).
The regular RN antiderivatives then give (13). In this asymptotic regime
w_J is nonzero, so the three boundary conditions determine J uniquely.
The limit tends to -1 as r_out grows. Taking r_out large still requires
the uniform-field perturbation to remain small; it is not a physical
uniform field imposed to spatial infinity.

This replaces the zero-conserved-charge even reference by a specific
odd classical load. It does not identify a physical binary ADM spin or
the quantum rotation sector. In particular, after including the LLL
angular current, the conserved displacement is
Q=(P^2/G)J_geom+sqrt(T_q/pi) chi_current, as derived in
[the rotational-current note](rotational-current-response.md).
J_geom alone is then not constant. Equations (10)-(13) supply a classical
comparison, not a license to ignore the matter displacement or the
time-dependent angular connection of a rotating binary.

## Verification, literature and next decision

[The checker](../checks/gluing_support_compatibility.py) verifies the
algebraic monotonicity identity, both odd throat equations, the RN
Wronskian and inhomogeneous equation. It independently integrates the
Euclidean arclength, enumerates all smooth common-phase branches, checks
the winding determinant numerically, and compares the weighted energy
and stress optima. Its separate boundary-value solve uses the original
inhomogeneous RN flux equation and angular equation, treating J as an
unknown boundary parameter, at two tolerances on five annuli. It starts
from the regular zero-spin field and compares with (11)-(12).
[The receipt](../evidence/gluing-support-compatibility.json) pins exact
inputs; this numerical corroboration supplements the ordinary proofs.

[MMP v3](https://arxiv.org/html/1807.04726v3) already establishes the LLL
representation, anomaly reference and long-cycle weakening of Casimir
support (sections 2, 5 and 6.1.1).
The explicit twist cancellation methods were also retained in the prior
return-holonomy review. The present result corrects an internal reference
compatibility choice and supplies its finite-delay bound. It is not
claimed as a novel general obstruction, a new support mechanism or a
stability theorem. [The review record](../evidence/gluing-support-source-review.json)
states the exact edition, scope and remaining proof-review gaps.

Continue long-return damping as primary and finite-duration band transfer
as fallback using the reflected zero-twist reference, or an explicitly
equivalent connection/seam choice. Derive the physical angular connection
and the coupled gauge/current state in that same construction. No Lean,
cloud, publication or outreach action occurs in this milestone.

Remaining to-do list: determine transported angular connection and the
coupled quantum state for the corrected reference; control omitted
response and smooth semiclassical matching; establish the physical
endpoint and consequential novelty, independently verify, GCP-certify
suitable mathematics and prepare the reproducible manuscript.
