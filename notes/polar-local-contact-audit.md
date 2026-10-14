# Polar local contact audit and missing constrained terms

S074, 2026-10-13. This audit recovers the specified aligned dipole block
from the preserved contact draft and identifies why its local terms do
not complete the physical one-loop action. The full S074 gate remains open.

## Provenance and conventions

The original untracked `checks/polar_dipole_contacts.py` has SHA-256
`14ca4c38524b91cf231ef50d50df29f59277cca99614241dd2289196c278427a`.
It lacked its referenced note and receipt. It is preserved unchanged.
The new `checks/polar_local_contact_audit.py` copies its algebra with that
provenance and adds independent determinant, constraint and angular tests.

Use the action and magnetic conventions in
[the canonical polar note](polar-canonical-vacuum.md) and the dipole
coefficients in [the source note](dipole-quadrupole-source.md).
Let c=cos(theta), Y=P2(c), V=(1-c^2)Y', and write the prescribed sphere
radius as b exp(chi Y). The base metric is g_ab+Y h_ab with tracefree h.
At a Lorentzian orthonormal base point set
h=[[p,q],[q,p]], J=p^2-q^2 and d(c)=sqrt(1-J Y^2).
Set alpha=a2/P=Acal/(b^2 E0). These are off-shell prescribed backgrounds;
the second-order gravitational constraints have not been eliminated.

## Restricted dipole block

The aligned gravitational/Maxwell curvature matrix has background
K0=diag(kg,ke), where kg=b^4/(6G), ke=4pi b^2/(3g^2).
The scalar dualization uses v=(-P,1) and mu=8pi/(3g^2).
The original horizontal Maxwell curvature includes
c Gcal + P alpha V F_A. Direct metric jets give the gravitational
connection-curvature coefficient, and direct contraction of F gives the
Maxwell block. Angular integration yields, through second polar order,

```
K11=kg[1-4chi/5+8chi^2/7+J/14]+ke P^2 (72alpha^2/35),
K22=ke[1+4chi/5+22chi^2/35+11J/70],
K12=ke P[6alpha/5+12chi alpha/35].
```

The unreduced mixed-vector coefficient follows from
(3/4) integral_(−1)^1 (1-c^2) d(c) (g+Yh)^(-1) dc.
Dualizing that two-dimensional tensor and expanding gives

```
Z^ab=g^ab+h^ab/5-11J g^ab/350,
f=sqrt(det(Z g))=1-9J/175.
```

Define r=mu P^2/kg and d0=mu/ke, so the original dipole mass is r+d0.
Inverting the entire K matrix before expanding gives M=mu v^T K^-1 v:

```
M0=r+d0,
M1=4(r-d0)chi/5+12r alpha/5,
M2=(-88r+2d0)chi^2/175+24r chi alpha/35
   +r[36/25-(108/175)r/d0]alpha^2-(5r+11d0)J/70.
```

The K1 K0^-1 K1 term in the inverse is necessary. The contact coefficients
cannot be obtained by retaining K2 alone. Integration by parts on free
dipole data recovers the previous cubic source coefficient
-(8pi/15) W R in the original (chi,Acal) variables. This checks the
on-shell cubic overlap only; it does not determine all off-shell contacts.

For a constant positive Euclidean principal matrix Z_E the scalar heat
density is exp(-t M)/(4pi t sqrt(det Z_E)). Its mass contribution to the
quadratic heat coefficient is -(M2+9(r+d0)J/175)/(4pi), with the
appropriate continuation of the Lorentzian invariant J.
The checker verifies this Gaussian integral by independent Fourier sums.
This constant-coefficient scalar calculation omits variable-background
curvature/connection terms, boundary heat coefficients, the constrained
measure and the other physical fields. It does not fix finite counterterms.
The general coefficient and boundary distinction follow
[Vassilevich, equations 2.14–2.17, 2.31–2.32 and 4.27–4.28](https://arxiv.org/pdf/hep-th/0306138v3).
No full four-dimensional divergence is extracted from this scalar block.

## An omitted angular channel occurs at the same perturbative order

Exactly,

```
P1 P2 = (2/5)P1+(3/5)P3,
integral_(−1)^1 P1 P2 P3 dc=6/35.
```

For unit-normalized m=0 harmonics, multiplication by P2 has
l=1 to l=3 matrix element 3sqrt(21)/35, whose square is 27/175.
Its total squared norm on Y10 is 11/35, of which only 4/25 remains
in l=1. This is an angular overlap, not a loop-size estimate. In particular,
it identifies an allowed l=1/l=3 mixing in the unreduced Maxwell curvature
block. Constraint reduction must determine the complete physical vertex.

For a regulated positive Hessian split into retained and omitted fields,

```
H(t)=[[A,t B],[t B*,D]],
(1/2)log det H(t)
 = (1/2)log det A+(1/2)log det D
   -(t^2/2)Tr(A^-1 B D^-1 B*)+O(t^4).
```

Thus a retained-block determinant misses mixed loops at the same quadratic
background order as local two-polar/two-dipole contacts.
This identity holds on the stated invertible domain; it is not a claim
that an unconstrained gravitational Euclidean Hessian is positive.

Nor is every omitted angular field heavy. The accepted product spectrum
for l=3 has

```
mu_3,-=120/[C(C+12+sqrt(C^2+24C+24))] ~60/C^2,
mu_dipole=2+4/C.
```

There is no mass hierarchy that makes this l=3 light branch heavier than
the retained dipole as C grows. Its physical mixing may still be suppressed;
that requires a vertex and norm bound. No such suppression or instability
is inferred from this kinematic observation.

## Nonlinear constraint and renormalization obligations

After linear reduction, suppose a second-order constraint variable n enters
as (1/2)n A n+n J2, with A invertible on a specified constrained domain.
Solving n=-A^-1 J2 adds -(1/2)J2 A^-1 J2 at fourth field order.
A multiplier with a noninvertible Hessian instead imposes its constraint
and cannot be eliminated by this formula. Boundary data and zero modes
must be resolved before using either operation.

Consequently the next action calculation must retain the off-diagonal
angular blocks and the second-order lapse/shift/Gauss sources as well as
the displayed K2 contacts. Reusing free equations inside a cubic vertex
does not supply these terms:
[Criado and Perez-Victoria, section 3, equation 10](https://arxiv.org/html/1811.09413v2#S3).

Use a common four-dimensional covariant regulator or a demonstrated
equivalent reduction prescription. Fix the cosmological, gravitational,
Maxwell and required higher-derivative bulk/boundary coefficients by explicit
renormalization conditions. Their finite parts are inputs until matched.
Normal-ordering a single dipole tadpole to zero does not determine them.

## Decision

Retain the exact aligned contacts and the cubic consistency check as a
component of S074. Reject the proposal that this block plus the previous
pair bubble exhausts the physical quadratic polar effective action.
This is an approximation-completeness decision, not a physical instability,
a no-go theorem or acceptance of a new wormhole result.

Remaining to-do list: derive the constrained mixed angular vertices and
second-order constraint sources; match the full local action and common
subtraction; quantify the physical error before S075 state/drive matching.
