# Coercive interior jets of the fixed-seed middle source

2026-10-21, E003/S167. Ordinary argument before supporting verification.

The complete fixed-seed middle trace has the improved bound

    |nabla^j Theta_seed|/(j! 100^j) < 1e22 q^2,
    0<=j<=3, s>=.001, Y>=.05.                         (1)

This replaces the 1e26 upper bound on precisely this source and region.
It does not replace the existing cold or outer estimates. Separately,
the actual gluing contribution has explicit fixed normal-chart coordinate
jets through nine, given by (11) below. Higher local UV/restoration jets,
the finite-comparison quantum source H6 norm and the positive transfer
constants are NOT evaluated here. No cancellation certificate passes.

The improvement comes from retaining the positive mass in an interior
Dirac energy identity. A derivative does not then require a factor M.
The previously proved three-derivative row bound is not extended by
assumption, and no spatial analyticity is assigned to the finite metric.

## 1. Unchanged seed, row conventions and observation neighborhood

Use the actual seed and middle normal charts of
[the canonical response proof](stationary-centered-canonical-response.md),
especially its compact non-bolt coverage and supporting-source sections.
Keep the actual HH state, paired domain, full angular tower, common
subtraction and reference scale. In A units,

    R0=1/32, |sigma_N|<.002 on the complex normal chart,
    M_n^2=n(n+q)/C, d_n=2(q+2n), q=2000^4,
    M_n>=M_1>18, (C/q)^2<31/2500000, C/q<.00352.

The fixed gluing cutoffs have slope at most 256 and source separation
D>=1/65. Both commutators remain. The plain zeroth evaluation row is
bounded by 3e4 exp(-M D/2); the rank-two spin trace is separate. Only
cap rows are used. The middle factor is the actual inverse, whose
M^-1 bound cancels the trace's M. No actual off-diagonal kernel value
is presumed. Fix these cutoffs while moving the observation point.

Set R=1/256. For |x|<=R the remaining proper separation satisfies

    D_x >= 1/65-1.003/256 > 1/88 > .19/18.          (2)

Thus the already established zeroth-row evaluation neighborhood still
fits, including at the smallest mass. After a conformal endpoint
multiplier the point bound is at most 31000 exp(-M/176) throughout
this observation ball. In particular its flat L2 norm F0 is at most
sqrt(pi) R times that bound. A row may be treated as a Hilbert-valued
spinor. If using its Hilbert--Schmidt norm, its initial sqrt(2) factor
produces exactly the rank-two factor in the eventual quadratic trace.

The chart is analytic in the two independent complex coordinates z,zbar.
On B_R, complexifying each Cartesian coordinate by less than
(R0-R)/2 stays within that chart. Cauchy's estimate, including the
Hilbert norm of the ordered Cartesian derivative tensor, therefore has
base at most 2sqrt(2)/(R0-R)<128. For a=exp(sigma_N), and also for
the endpoint multipliers exp(+-sigma_N/2), this gives

    a>=exp(-.002),
    |partial^b(a-1)| <= A_b=(21/10000) b!128^b, b>=1. (3)

The same A_b bounds the nonconstant endpoint jets. The zeroth endpoint
value is one at the normal-chart center. These estimates use the seed
inside the ball; they require no higher derivatives of the distant cap.

## 2. Interior energy with its mass retained

In the flat conformal spin frame let D*= -D and D*D=-Delta. Conformal
covariance turns a homogeneous massive row into

    (D+M a)v=0.

This is the two-dimensional instance of the conformal Dirac formula;
see [Walpuski, Spin Geometry, Proposition 15.5, section 15.2](https://walpu.ski/Teaching/SpinGeometry.pdf),
printed page 79, PDF index 78. The factor i converting the self-adjoint
geometric convention to the skew-adjoint Euclidean convention does not
change its conformal weights. Physical and half-density row conventions
differ by the endpoint factors already included in (3).

For any compactly supported spinor w in the ball, direct integration
by parts gives

    ||(D+Ma)w||^2 = ||grad w||^2 + M^2||a w||^2
                                  -M<w,c(da)w>.

The sign of the last term is not discarded. Its norm is bounded by
M A_1||w||^2. Since exp(-.004)>.996 and M>=18,

    ||(D+Ma)w||^2 >= .98 (||grad w||^2+M^2||w||^2). (4)

This coercivity is local to compact tests and does not change the actual
global boundary condition. It also holds for adjoint rows, by reversing
the harmless sign in the cross term.

Fix a target j in {0,...,9}. Put m=j+3, w0=R/m, L=2/w0=512m,
B=600m, r_k=R-k w0. Use nested cutoffs with slope at most L between
B_(r_k) and B_(r_(k+1)), for k=0,...,j+1. Let F_k denote the flat L2
norm of the full kth Cartesian derivative tensor on B_(r_k), and E_k
the combined (k+1)st derivative and M times kth derivative norm on
B_(r_(k+1)). Differentiating the local equation, then applying (4), gives

    E_0 <= (11/10)L F_0,
    E_k <= (11/10)[L E_(k-1)
                  +sum_(b=1)^k binom(k,b) A_b E_(k-b)].             (5)

For the last inequality, F_k<=E_(k-1) on the required ball, and
M F_(k-b)<=E_(k-b) on its larger ball. The b=1 term is already known;
there is no circular use of E_k. Tensor Leibniz estimates retain all
multiindices through their Hilbert norms. Exact rational arithmetic gives

    (11/10)[L/B+(21/10000)
       sum_(b=1)^k k!/(k-b)! (128/B)^b] < 1,
    0<=k<=j+1, 0<=j<=9.                            (6)

Induction proves E_k<=B^(k+1)F0, hence F_k<=B^k F0 through j+2.
The outer point-row exponential was weakened to exp(-M/176) in (2),
so there is no contradiction with the growth of massive homogeneous
solutions. An estimate at an unchanged separation is not asserted.

## 3. Point evaluation and physical endpoint factors

The last ball has radius w0. Take rho=w0/2 and a radial cutoff which
is one through rho and zero at 2rho. The quintic transition
1-10t^3+15t^4-6t^5 has |chi'|<=2/rho and |Delta chi|<=8/rho^2.
The flat Fourier identity in two dimensions gives, for H2 functions,

    |f(0)| <= ||(1-rho^2 Delta)f||_2/(2sqrt(pi)rho).

Apply it to chi partial^j v. Since the norm of a two-dimensional
Laplacian is at most twice the full second derivative norm,

    |partial^j v(0)|
      <=[9F_j+4rho F_(j+1)+2rho^2 F_(j+2)]/(2sqrt(pi)rho)
      <=31000 m [9+4(75/64)+2(75/64)^2] B^j exp(-M/176). (7)

The cutoffs may be approximated in H2 if smooth compact tests are
desired. No pointwise estimate is inferred just from an H1 norm.

Leibniz with the physical endpoint factor changes the coefficient in
(7) by at most

    1+(21/10000) sum_(b=1)^j j!/(j-b)! (128/B)^b < 1.01.

Thus the plain physical row has the following conservative constants:

    ||row_j|| <= c_j exp(-M/176),
    c_0=30000,
    c_j=520000(j+3)[600(j+3)]^j, 1<=j<=9.           (8)

For j=0 use the sharper original zeroth-row bound, weakened to the
new separation. For higher jets (7) and the endpoint conversion prove
(8). If Hilbert--Schmidt rows are used, retain sqrt(2)c_j instead.
These are coordinate jets in a FIXED chart, not derivatives of a moving
normalizer or a family of quantum states.

## 4. Full angular sum and the trace ledger

No angular truncation is needed for the proof. Write c=1/88 and
h(n)=2(q+2n)exp(-c sqrt(n(n+q)/C)) for real n>=0. In the mass variable
its logarithmic derivative is

    4CM/(q^2+4CM^2)-c <= sqrt(C)/q-c < 0.

The last inequality follows from C/q<.00352 and q=2000^4. Therefore
h is decreasing and the right-endpoint sum is bounded by its integral:

    sum_(n>=1) d_n exp(-M_n/88)
       <= integral_0^infinity h(n)dn
       =4C integral_0^infinity M exp(-M/88)dM
       =4C 88^2.                                  (9)

Here d_n=2C d(M_n^2)/dn exactly. This retains the paired multiplicity.
For a finite diagnostic sum through N the uncomputed tail is at most
4C(88^2+88M_N)exp(-M_N/88); diagnostics do not replace (9).

Let G be the actual finite local/global gluing part of Theta_seed in
the SAME decomposition used by the earlier source proof. Two rows,
the rank-two trace, both commutators and the original middle inverse give

    |partial^j G| <= 2C 256^2
       sum_(k=0)^j binom(j,k)c_k c_(j-k)
       sum_(n>=1)d_n exp(-M_n/88).                 (10)

All observation derivatives fall on the two evaluation rows. The
middle inverse and the cutoff commutators are fixed in this chart.
Equation (9) proves the explicit dimensionless coordinate bound

    |partial^j G|/(j!100^j q^2) <= G_j,
    G_j=8*256^2*88^2*(31/2500000)
            *sum_(k=0)^j binom(j,k)c_k c_(j-k)/(j!100^j),
    0<=j<=9.                                      (11)

The supporting checker retains every rational G_j, including j=4,...,9.
These higher numbers describe only G. They are not a nine-derivative
bound for the complete renormalized trace or the finite-metric residual.

## 5. Complete third source jet and quantitative consequence

At the center sigma_N=0 and grad sigma_N=0. Thus coordinate and physical
covariant scalar derivatives agree through order two. At order three
their difference is -(partial Gamma) partial G. The conformal
Christoffel formula bounds its tensor coefficient by four times the
Hessian norm of sigma_N. The same Cauchy guard gives

    |partial Gamma| <=4(.002)2!128^2 <300.

Consequently the third normalized covariant gluing bound is
G_3+G_1/200. Use the PREVIOUSLY PROVED complete local UV and common
finite-restoration bounds, 2e13 and 1e12 respectively, through order
three. Exact arithmetic gives

    max(G_0,G_1,G_2,G_3+G_1/200)+2e13+1e12 < 1e22. (12)

This proves (1), an improvement of four orders over the earlier complete
1e26 allowance. The comparison concerns the same seed, region, derivative
normalization and source decomposition. It is not a change of units.

The established [metric transfer](stationary-completed-metric-response-norm.md)
uses this middle source-three bound as an input. Its source-dependent
term .06 gamma K_m, gamma<6.8e-28, can now be bounded by 4.08e-7 instead
of 4.08e-3. Its direct high-curvature source term 6.6 gamma K_m is below
4.488e-5 instead of .4488. Other low-jet terms must be added separately;
the full force budget 10 gamma K_m from that older ledger is NOT simply
rescaled. No improved endpoint C5 contract is adopted in this note.

The next quantitative work should propagate these bounds through the
actual invariant metric equations and quantify the remaining higher
local/cold/outer terms. The target is still the actual finite-family
R0 R2 R4 N6(p) or a sharper actual q3 bound from
[the graph admissibility proof](finite-graph-admissibility.md). None of
those four quantities is supplied by (1) or (11).

## Verification contract and limits

Freeze this ordinary proof, checker and review before execution. Check
the separation, coefficient coercivity, every finite recurrence and
endpoint conversion, Fourier cutoff constants, exact angular monotonicity,
all rational G_j, the physical third-derivative conversion, and (12).
Independently integrate a variable-mass two-dimensional Dirac energy
identity at two precisions, retaining its gradient-of-mass cross term.
Check two guard-point angular sums with an explicit uncomputed tail.
They are diagnostics, not the actual charged geometry or quantum source.

The continuum PDE estimates, actual cap-row application and complete
source decomposition remain ordinary proofs. No Lean run or physical
source evaluation is claimed. Retain every scientific run and any failed
source snapshot; do not change a successful scientific snapshot afterward.

Remaining to-do list: propagate the improved middle estimate and quantify
higher local/cold/outer and finite-metric source/transfer bounds; actual
cancellation and signed compatibility; preserved nonlinear existence;
physical boundary/state and validity, stability, matching, charge/return,
opening, transport, preparation, separate creation and full audit.
