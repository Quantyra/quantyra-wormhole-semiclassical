# Outer canonical response at finite backgrounds

2026-10-19, E003/S146 supporting S139/S074. Ordinary argument and
verification protocol before execution.

The outer canonical scalar trace has a second-response estimate in
the weighted real-C5 coefficient space below. This is the proper
collar counterpart of the interior conformal result S145. It does
not identify the two coefficient spaces or solve the exterior equations.

## 1. Fixed coordinates and the input contract

Retain the centered seed, angular masses, paired domain, scale and
physical period of [the centered canonical response](stationary-centered-canonical-response.md).
In orbit units,

    M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18,
    A(t)=sum d_n M_n^2 exp(-t M_n^2)<=24C/t^2,
    q/C<400, eta_G=GC/B^2<8.5e-33.

All coordinates are fixed throughout this argument. Near infinity the
half-density operator is

    D=gamma_Y partial_Y+gamma_t f(Y) partial_t.

The covariant scalar heat diagonal includes the outside factor f(Y).
The constant complementary projectors, geodesic boundary and identical
pair intertwiner are retained. The even-double spin trace equals the
half-domain trace; no extra factor two or one half is inserted.

Write f_b=f_0 exp(b), b=p psi+v. The fixed analytic mass direction psi
is the function v_m/p in the centered response note, with

    psi(Y)=gamma Y^3+O(Y^5), |gamma|<.07.

It has the retained analytic majorant ten on the outer .01 strip,
and is capped after .4, to zero by .5. Its cap value/gradient bound
is 100. The remainder v satisfies, on 0<Y<=.1,

    |v^(j)(Y)|<=V j!100^j Y^(4-j)L(Y), j<=5,
    L(Y)=1+|log Y|.

The fixed even cap extension has |v|+|v'|<=100V. Define
N(b)=|p|+V, and define U=N(u), W=N(w) for directions. The actual
global extensions have relative first-order coefficient/connection
bound 1e4 N, on the same fixed paired domain. They agree with their
caps near the outer observation; the actual-to-cap comparison uses
the distance-24 construction and graph bound ten in the centered
response note. These are hypotheses, not an inferred transfer from
an arbitrary C0, C2 or product-inverse ball.

The real background satisfies N(b)<=1e-7. The global cap has
|f_0'/f_0|<=2 and ratio R_even<=1e10, as established for the centered
seed. Complexification below concerns parameters only. The spatial
function v remains real-C5 on the open collar; its fourth and fifth
derivatives need not extend continuously to Y=0.

Let Theta[b]=B^2 tau_hat, with the Lorentzian canonical trace
convention of the boundary-layer note. If a_b=1/f_b and a_b(0)=1,
the cubic coefficient of a_b is exactly

    c[b]=c_0-gamma p.                                  (1)

The weight-four logarithmic remainder has zero cubic jet. Products
of b with itself start at order six, so (1) is exact at finite p,v,
not a first-order approximation. Define

    R[b](Y)=Theta[b](Y)-c[b] C^2/(5pi Y),
    ||F||_out=sup_(0<Y<=.05) |F(Y)|/[1+log^2 Y].

The result, under the stated contract, is

    ||D^2 R[b][u,w]||_out <1e30 C^2 U W.               (2)

For the explicitly rescaled scalar Q=2G R/(B^2 C), (2) gives

    ||D^2 Q[b][u,w]||_out <.017 U W <.02 U W.          (3)

Thus DQ is .02-Lipschitz along a real segment in the ball, and
its Taylor remainder is bounded by .01 N(h)^2. This is a scalar
estimate in the displayed norm, not a full stress/radius estimate
or a contraction after composition with an exterior inverse.

## 2. Exact parameter responses of the bridge integrand

For a variance-two unit bridge beta, put

    X=Y+sqrt(t) beta, F=<f(X)^2>, G=<f'(X)>,
    J=f(Y) F^(-1/2) exp(t G^2/(4F)),
    k(t;Y,Y)=E[J]/(2pi t).

Brackets are time averages along the path. For f=f_b exp(su+rw),

    F_u=<2f^2 u>, F_w=<2f^2 w>, F_uw=<4f^2 u w>,
    G_u=<f'u+f u'>, G_w=<f'w+f w'>,
    G_uw=<f'u w+f(u'w+uw')>.

Let ell=log J. Exactly,

    ell_u=u(Y)-F_u/(2F)
          +t[2G G_u/F-G^2 F_u/F^2]/4,

    ell_uw=-(F_uw/F-F_u F_w/F^2)/2
       +t/4 [2(G_u G_w+G G_uw)/F
          -2G(G_u F_w+G_w F_u)/F^2
          -G^2 F_uw/F^2+2G^2 F_u F_w/F^3],

    J_uw=J(ell_uw+ell_u ell_w).                        (4)

The last product and the observation terms in ell_u,ell_w are
essential. No indicator for a crossing or exit event is differentiated.
These identities retain the exponential seagulls and scalar density.
They also provide independent finite-path tests of the parameter
calculus used below.

For real f the formula follows by Fourier transforming the two
scalar squares and doing the tangential Gaussian integral after
Feynman--Kac. [Simon, theorem 1.1](https://arxiv.org/html/math-ph/9907022)
supplies the semibounded real bridge formula; his generator is
Delta/2, so the bridge here has twice his variance. The complex
parameter extension below follows from domination and uniqueness
of analytic continuation, not a claim that his real theorem alone
proves a complex or renormalized estimate.

## 3. A uniform cap tube, before differentiating a remainder

For this section only allow N(b)<=rho=1e-4 with complex p,v.
The actual global resolvent is not continued to this larger tube.
The auxiliary cap value/gradient perturbation is at most .01.
Use the looser bound delta=.02 for both.

Let F_0=<f_0^2>. On every real path,

    Re F>=.95 F_0, |F|<=1.05 F_0,
    |G|<=2.1 sqrt(F_0),
    |J|<=1.1 R_even exp(2t).                           (5)

Indeed e^(-2delta) cos(2delta)>.95,
e^(2delta)<1.05 and (2+delta)e^delta<2.1.
The exponent is bounded by 2.1^2 t/(4*.95)<2t.
Choose the square root branch continuous from positive F.

For a direction u of norm U the global bounds are |u|,|u'|<=100U.
Consequently

    |F_u|<=210U F_0, |G_u|<=310U sqrt(F_0),
    |ell_u|<=1000U(1+t),
    |D_u k|<=1e4 U R_even t^-1(1+t) exp(2t).           (6)

These are absolute bounds for complex backgrounds. They establish
holomorphy and domination of each unrenormalized cap kernel away
from t=0. For t>=.01, split the angular exponential at half its
argument, using the full tower rather than a finite truncation:

    A(t)<=96C t^-2 exp(-162t).

Since integral_(.01)^infinity t^-3(1+t) dt=5100, (6) contributes
less than 1e20 U C^2. The use of exp(-160t)<=1 is deliberately
loose. The cap ratio is the centered seed's 1e10; a perturbation
factor is already included in 1.1 in (5).

Here is the small-time extension of the grouped real-C5 proof;
the details matter because real first-variation bounds do not by
themselves imply (2). Split b=p psi+v and first hold p fixed.
On a safe path at Y use r=Y/[1000(1+S)], S=sup|beta|.
The weighted v finite-jet norm is at most 15 V Y^4 L, and the
v' norm needed by the groups is at most 3750 V Y^3 L.
For V<=rho their maxima on Y<=.1 are below 5e-7 and .0013.
The analytic factor p psi has norm below .001; after halving
its .01 strip its derivative norm is at most .2. Thus the
normalized local coefficient cone has |f/f(Y)|<1.01,
Re(F/f(Y)^2)>.95 and |G/f(Y)|<3.

These stronger margins fit the earlier cone Re F>.5, |G|<6.
For a weight-four direction, retain the following enlarged
finite-jet bounds for the three parameter-differentiated groups:

    F^-1/2:                         3000 U Y^4 L,
    F^-3/2 G^2/4:                  1e6 U Y^3 L,
    F^-5/2 G^4/32:                 2e7 U Y^3 L.

The outside f(Y) and its parameter derivative are included.
The bounds follow by the product rule, |F^-1|<=2 and |G|<=6
in the finite-jet algebra, exactly as in the real calculation;
the factor ten covers exp(p psi+v) and its tagged product/inverse
factors. The maximal additional exponential factor is
exp(2*.001001)<1.003. The smallest inverse-cone margin is .95;
even ten factors of its reciprocal cost less than 1.68.
Together with the derivative-cone margin, these are below ten.

Expand the first group through degree four in sqrt(t), the
second through degree two and the third through degree zero.
Their real Taylor remainders use v5, v4 and v2 respectively.
The remaining exponential group is O(t^3). The direction is
never evaluated at a complex spatial point. Odd retained
coefficients average to zero by bridge reflection. Applying
E(1+S)^5<3000 gives the following ledger, with the common
a0,a1,a2 terms removed before angular integration:

| Contribution for a weight-four direction | Uniform complex-cap bound |
| --- | --- |
| Safe local, 0<t<Y^2 | 1e22 U t^(3/2) L/Y |
| Unsafe local paths, 0<t<Y^2 | 1e14 U Y^4 L t^-1 exp(-Y^2/(8t)) |
| Y^2<=t<=.01, unsubtracted kernel | 1e14 U t(1+abs(log t)) |
| Fixed collar exits | 1e11 U t^-1 exp(-.00125/t) |
| Local k1 variation | 1e7 U Y^2 L |
| Local k2 variation | 1e13 U L |
| Finite restorations and long polynomial subtractions | 1e17 U C^2 L |

The first row's leading coefficient before dividing by 2pi is
3000*1000^5*3000=9e21. The other two grouped remainders fit
inside 1e22 after that division. The unsafe and exit rows use
the original bridge tail P(S>a)<=2exp(-a^2), with (5)-(6);
the larger-short row uses E[(1+S)^10 exp(S)]<1e9.
Each is ten times the corresponding real coefficient. The cone
and domination calculation above justify this enlargement for
complex backgrounds; it is not obtained by differentiating an
unsigned O(log Y) statement.

For the analytic mass direction at v=0, use the fixed radius
1/[1000(1+S)] and the majorant ten of psi. The same grouped
calculation gives a smooth subtracted first-response bound
5e23 |dp| t^(3/2). After its cubic profile is removed, the
reflected remainder is bounded by 1e16 |dp| t^(3/2).
The first new reflected products have parabolic weight five:
the degree-three insertion multiplied by the degree-two seed;
two degree-three insertions start at weight six. There is no
unremoved weight-three term. The perturbed analytic factor is
below .001 on its strip, so the same product/inverse margin
above covers the first-response constants uniformly in p.

We only need that analytic mass estimate at v=0. At nonzero
v use the exact path integral

    R_cap[p,v]-R_cap[0,0]
      =R_cap[p,0]-R_cap[0,0]
        +integral_0^1 D_v R_cap[p,tau v][v] d tau.      (7)

This avoids pretending that the nonanalytic background is an
analytic mass curve, or differentiating an unproved asymptotic.

## 4. Common subtraction, cubic profile and complete angular limit

The cubic profile is the exact parabolic-weight-three model

    K3(t,Y)=c/(2pi t) integral_0^1
       D(Y,sqrt(2t z(1-z))) dz,
    D(Y,sigma)=E[|Y+X|^3-(Y+X)^3].

Using A(t)=2C/t^2+O(1), its leading full angular integral is
c C^2/(5pi Y). The elementary identity is

    integral_0^infinity t^-3 D(Y,sqrt(vt)) dt
      =3v^2/(2Y),
    integral_0^1 [2z(1-z)]^2 dz=2/15.

This is the same rank-two scalar normalization as the
[boundary-layer calculation](stationary-massive-boundary-layer.md).
Equation (1) implies D^2 c=0 in the chosen coefficient coordinates.
The singular profile therefore has zero mixed derivative. It is
still present in the value and first derivative and must be paired
with the opposite Weyl density before physical mass integration.

For the remaining angular integral, A(t)<=24C/t^2 and the ledger give
the following coefficients after extracting U C^2(1+log^2 Y):

    safe weighted term:        <=9.6e23,
    unsafe weighted term:     <=3.072e17,
    larger-short/local terms: <=1e18,
    fixed exits:              <=1.536e18,
    long term:                <=1e20,
    common finite terms:      <=2e17.

Here L<=2(1+log^2 Y) and L^2<=2(1+log^2 Y).
For example the safe integral is 48e22 L;
the unsafe integral uses integral_0^infinity
t^-3 exp(-Y^2/(8t))dt=64/Y^4.
For the analytic mass branch, the smooth short integral is
24*(5e23)*.2=2.4e24, and the reflected, long and local
terms fit the remaining reserve. Consequently (7) proves

    ||R_cap[b]-R_cap[0]||_out
       <1e25 C^2 N(b),  N(b)<=1e-4.                  (8)

The common a0 scalar restoration has zero metric variation,
because B, masses and mu0 are fixed. It has not been assigned
zero value. The common a2 finite restoration is included in
the table; it is not separately subtracted from the observable.
The finite-circle windings retain the physical period. The
inherited distance d>=1e17 in their Gaussian estimate and
M1>18 make their contribution smaller than the unused reserve.

Apply these bounds first with angular and proper-time cutoffs
and the same common local counterterms. The integrable absolute
remainders are uniform in the complex tube after subtracting
the explicit cubic. They give locally uniform convergence in
the displayed weighted norm. Thus the limit is holomorphic in
each parameter; its Cauchy derivatives commute with the common
subtraction and the full angular limit. This is the missing
logical step between a first-response estimate and (2).

## 5. Actual-domain correction and the second-response bound

The cap tube is larger than the tube used for the actual domain.
For the latter, fix a real background with N(b)<=1e-7 and vary
it on a parameter disk of norm radius 1e-8. The relative
coefficient bound is at most .0011; graph bound ten makes the
Neumann ratio below .02. Scalar weights preserve the paired
domain. The separated cap rows and both cutoff commutators
therefore remain valid with their inherited exponential rate.

The centered response's two-sided geometric-resolvent calculation
gives the directional bound 1e62 U M^5 exp(-9.6M) at the seed.
Use 1e64 on this tube. All changed inverse factors have norm
ratio at most 1/.98; the real and complex local evaluation
reserves and the coefficient-direction bounds are those already
used there. The extra factor 100 exceeds their product. This
comparison retains the actual global middle inverse, not an
unmatched half-line or complete-H2 replacement.

The complete angular inequality from that calculation is

    sum d_n M_n^5 exp(-9.6 M_n)<1.2e-67 q.

It can also be checked without a numerical angular truncation:
M_n>=18sqrt(n), d_n/q<=2+4n, and x^5 exp(-9.6x) decreases
for x>=18. The resulting summand
18^5(2+4n)n^(5/2)exp(-172.8sqrt(n)) decreases for n>=1.
Its sum is bounded by its first term plus the integral from one
to infinity. Substituting z=sqrt(n) reduces that integral to
18^5 integral_1^infinity (4z^6+8z^8)exp(-172.8z) dz.
Finite integration by parts and a rational lower Taylor bound
for exp(172.8) prove the displayed inequality.

The first-response actual-cap difference is therefore below
.48 U C^2. Cauchy on radius 1e-8/W bounds its mixed response
by 4.8e7 C^2 U W, below 1e8 C^2 U W.

For the cap, set r_s=1/(40000U), r_t=1/(40000W).
The background plus the two parameter directions has norm at
most .0000501, so (8) is uniform on the whole bidisk. Subtract
the constant R_cap[0] before Cauchy. The mixed coefficient is
bounded by

    1e25 * .0000501 * 40000^2 C^2 U W
      =8.016e29 C^2 U W.

There is no extra factorial for one derivative in each parameter.
Adding the actual-domain correction remains strictly below
1e30 C^2 U W, proving (2). The cases U=0 or W=0 follow directly.
The fundamental theorem of calculus on a real segment gives
the Lipschitz and Taylor statements in (3).

## 6. Verification and limitations

The focused checker is to verify (4) by symbolic differentiation,
including the product and observation terms; the affine cubic
allocation and absence of a mixed singular coefficient; the
complex cone and integral budgets; and the weighted Cauchy
coefficient. Independent deterministic bridge-path quadratures,
mixed finite differences at two steps, angular sums with a
controlled tail, and the Gaussian cubic-profile integral test
the implementation. Tolerances are frozen in the checker before
execution. A deliberately omitted product term must fail its
diagnostic. Numerical checks are not continuum operator proofs.

This result controls the canonical outer scalar under the
explicit fixed-coordinate contract. The physical metric-to-
coefficient map, global exterior norm and derivative requirements,
Weyl/LLL and radius sources, global charge and coupled solve are
not supplied by (2). In particular an interior C5 source bound
and this weighted scalar bound cannot simply be multiplied by
the product inverse's 1.3e8 constant.

Remaining to-do list: transfer the interior and outer estimates
to the actual coupled exterior norm; include all local, massless
and radius responses; solve the exterior/state and signed source;
establish opening, stability, transport, preparation and creation.
