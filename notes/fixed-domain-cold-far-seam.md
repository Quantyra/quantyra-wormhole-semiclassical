# Fixed-domain cold far-seam response

2026-10-22, E003/S167. Ordinary argument before frozen supporting checks.

This closes the DISTANT actual-to-H2-cap contribution left open by the
[admission audit](finite-cold-admission-and-response.md), for its explicit
finite-centered conformal coefficient family. It does not close the local
flat-to-H2-cap contribution, cold admission, or the full metric solve.
The actual paired inverse remains on the physical domain throughout.
No complete H2 transvection is applied to that domain.

The gain comes from substituting the actual global low metric bounds
before parameter differentiation. It is not a change of norm units.
All angular masses, charge, period, subtraction scale and physical density
conventions remain those of the retained construction.

## 1. The family and the two fixed domains

Fix ell in [0,1] and v in the preceding explicit D_ell space, with
d=D_ell(v)<=5e-10. Put

    f=ell u+v, epsilon=1e-6 ell+1000d<=1.5e-6<2e-6.

The actual metric is sigma0+z f. Its auxiliary cap is
sigma0,cap+z Ef, where E is the preceding fixed width-one C2 extension.
For every complex z these agree on the fixed physical patch through
chi1, the auxiliary s=.01 seam. Neither seam nor physical boundary moves.
The sphere masses M_n and independent radius direction are held fixed:
this is the massive canonical sector, not a full metric/radius theorem.

The retained actual low estimates and the explicit cap construction give

    sup_actual |f|<=.012 epsilon,
    sup_cap |Ef|<=.03 epsilon,
    sup_actual |grad f|<=8 epsilon,
    sup_cap |grad Ef|<=.6 epsilon.                    (1)

For example 7.7e-6 ell+d<=7.7 epsilon, and
3e-8 ell+d<=.03 epsilon. These use the actual global norms in D; they
do not infer a global bound from a cold local jet. In a fixed SEED normal
chart at any cold observation, on radius 1/32,

    |grad f|, |Hess f|<=100 epsilon a_x,
    a_x=s[1+log(.01/s)], exp(-80)<h<=s<=.001.           (2)

For u, integrate its retained radial curvature equation before converting
to the seed normal frame. The proper gradient and Hessian bounds are
1.2e-6 s(L+2) and 3.6e-6 s(L+2); L>=1+log 10. On this small chart the
weight changes by less than 1.1. The fixed frame factors are below 1.2,
and its connection contributes at most four times the gradient. These
terms fit below 100 epsilon0 a_x. For v, its first/second divided jets
are <=100d a_x and 20000d a_x; even a factor two for the full real
Hessian leaves <=40(1000d)a_x. No sixth derivative is used.

On each parameter disk centered at t in [0,1] choose

    |z-t|<=r_z=1e-5/epsilon.                           (3)

When epsilon=0 the asserted difference is zero. Otherwise (1) gives
| (z-t) f |<=3e-7 and actual gradient <=8e-5. On the cap
|sigma_cap,z|<.043. Equation (2) gives additions below 4e-6 to the local
first/second seed-normal jets, even including the real center t.
The seed's normal conformal factor has value <.025, gradient <=4 and
Hessian <=400. Thus in the fixed normal chart the constant-principal
Dirac equation has scalar mass potential Q=M exp(sigma_N,z) with

    |Q|<=1.1M, |grad Q|<=10M, |Hess Q|<=2000M.         (4)

The gradient/Hessian bounds include the Clifford contraction norms.
The complexification is only in z; every unknown spatial field remains
real C5, and the cap is C2 at its two joins.

## 2. Common domains and the weighted inverse

In the fixed seed Hilbert space the first-order operator is

    A_z=e^(-z f/2) A_0 e^(-z f/2),
    S_z=e^(z f/2) (A_0+M e^(z f))^-1 e^(z f/2).       (5)

Use Ef on the cap. The real paired A_t is anti-selfadjoint in its
canonical representation. The scalar multipliers in (5) are bounded
with bounded first derivatives and commute with the fixed paired trace
projector. They preserve its first-order domain. Thus the middle inverse
is a bounded-potential perturbation on ONE actual domain. There is no
assumption of a fixed second-order Robin expression under metric change.

Relative to any real t, the mass perturbation obeys
|exp((z-t)f)-1|<=exp(3e-7)-1. The actual inverse therefore has norm
<2/M uniformly on (3). This follows directly from the Neumann identity
for the mass multiplier, not a Killing gap or an unspecified graph norm.
It proves holomorphy on the physical paired domain.

For the complete cap use the fixed H2 Dirac operator, cap factor <.043,
and a distance weight with Lipschitz slope M/4. The Clifford commutator
is bounded by M/4 in the H2 representation. The accretive mass term
has real part at least M exp(-.043) cos(3e-7). Consequently

    M ||e^w S_cap,z e^-w||
       <=exp(.043)/[exp(-.043)cos(3e-7)-1/4]<2.        (6)

The same statement holds for the adjoint column and bounded truncations
of the weight. Taking the truncation limit gives the separated L2 bound.
No derivative of the distance is required beyond its Lipschitz bound.
This use of the COMPLETE cap does not transport the finite actual disk.

## 3. Separated point rows, with their evaluation argument

Let the support of the seam commutator be at H2 distance at least D from
the observation. First omit that commutator and take an arbitrary L2 input
supported there. In the fixed normal chart the transformed cap solution
psi solves (D_flat+Q)psi=0 on the evaluation ball. Use nested radii
r,2r,3r,4r with r=1/(20M); M>=18. They fit within radius 1/32.
Radial quintic cutoffs may be bounded by 80M for their gradients and
11200M^2 for their Laplacians. The normal-to-H2 distance conversion is
<1.1 on this ball. Formula (6) costs at most exp(.055)<1.2 when the
outer ball is subtracted from D.

Write U for the local L2 norm on the outer ball. The constant-principal
Dirac energy identity ||D_flat phi||_2=||grad phi||_2 and (4) give,
successively on the inner balls,

    ||grad psi||<=82M U,
    ||Hess psi||<=6730M^2 U.                         (7)

For the second bound differentiate the equation once before applying the
cutoff: 82^2 M^2+10M<=6730M^2 for M>=18. This uses the tensor norm
of the full first derivative; there is no missing component sum.

With Hermitian Euclidean gamma matrices D_flat^2=Delta. Hence

    (-Delta+M^2)psi=(M^2-Q^2+gamma.grad Q)psi.

After cutoff, the right side has norm <25000M^2 U. For the cutoff
first derivative the corresponding bound is <2.1e6 M^3 U: its terms
are (M^2-Q^2+gamma.grad Q)grad psi,
(gamma.Hess Q-2Q grad Q)psi, and the two cutoff terms. In the same
units their sum is bounded by

    (2.21+10/18)82 + 22/18+2000/18^2
                      +160(6730)+11200(82)<2.1e6.     (8)

For scalar or finite-vector-valued functions the whole-plane point norm
of (-Delta+M^2)^-1 is exactly 1/(2 sqrt(pi)M), from its Fourier integral.
Combining (6)-(8), the distance loss <1.2 and an overall factor four for
local measure/spinor/frame conversion gives the separated row bounds

    ||ev_x S_cap,z 1_seam|| <=3e5 exp(-MD/4),
    ||grad_x S_cap,z 1_seam|| <=1e10 M exp(-MD/4).      (9)

Indeed the value constant is bounded by
25000/(2sqrt(pi))*2*1.2*4<3e5. The derivative principal part
is below 2.1e6/(2sqrt(pi))*2*1.2*4; the conversion connection
costs at most 100 times the value row and fits well inside 1e10M.
The adjoint-column estimates use the same calculation with -D_flat.

These are bounds for physical observation derivatives in a chart held
fixed around the point where the bound is evaluated. They do not
differentiate the family of chart choices or an actual-domain transvection.
The local cutoff equations also lift L2 parameter holomorphy to point
and first-derivative holomorphy, with uniform bounds (9). Thus pointwise
parameter Cauchy below is justified; it is not applied to an unbounded
point-evaluation functional on bare L2.

## 4. Differentiate the fixed-domain identity

Fix a physical radial cutoff chi equal to one through chi1-.25 and zero
at chi1, with commutator C_z=[A_z,chi] of norm <9. The width-.25
quintic has slope <=7.5 in chi units; the cap factor costs less than
exp(.043). The actual and cap coefficients coincide on its support.
Both inverses, the commutator and the cutoff use the same common chart.

On the diagonal inside chi=1 the exact identity is

    S_actual,z-S_cap,z
                  =-S_cap,z C_z S_actual,z C_z S_cap,z. (10)

To justify it, apply the first-order geometric-resolvent identity on each
side. The remaining term is one half S_cap[A,chi^2]S_cap. This equals
one half (chi^2 S_cap-S_cap chi^2), whose kernel is zero at the indicated
diagonal. No divergent trace is cycled. The argument is valid for complex
z as well as real z. The actual inverse in (10) is retained, with its
actual physical boundary data.

Let H_M(z,x) be the regular per-mode canonical scalar from (10), including
the external M and physical diagonal-density factors. Rank two, both
commutators, the actual inverse <2/M, (9), and a factor four for the
diagonal-density conversion give

    |H_M(z,x)|<=2e14 exp(-MD/2),
    |grad_x H_M(z,x)|<=1e20 M exp(-MD/2).              (11)

For value the product is 4(2)(2)(9^2)(3e5)^2<2e14.
For gradient, either cap row can be differentiated; the product before
the density derivative is 4(2)(2)(2)(9^2)(3e5)(1e10)M<1e19M.
The derivative of the local density factor is below twenty times the
value product, and fits the stated 1e20M reserve. The fixed actual
middle inverse and cutoff have NO observation derivative. Parameter
variations of all five factors and of the density are instead included
by applying Cauchy to the entire analytic expression (10).

Cauchy on (3), followed by integration over t in [0,1], proves

    |H_M(1,x)-H_M(0,x)|<=2e19 epsilon exp(-MD/2),
    |grad_x[H_M(1,x)-H_M(0,x)]|
                         <=1e25 epsilon M exp(-MD/2). (12)

The subtracted seed far term is kept. No vacuum term is assigned zero
value. Since the actual and cap metrics coincide near x, their local
subtractions and finite restorations cancel in H_M at each z. This is
local equality between two operators, not metric independence of the
full product determinant or permission to omit its restoration elsewhere.

## 5. Substitute the cold distances and sum every angular mode

Use the slightly weakened uniform guard

    D(s)>=.945 log(.01/s)-.25.                        (13)

The retained .946 guard and finite metric value perturbation leave this
margin. It gives D>1.9 for s<=.001 and D>60 for chi<=1, since
s<=1.62h and h<1e-30. In particular no distance-60 estimate is used on
the entire cold patch.

On the actual finite h range, a useful global conversion is

    1/[a_x exp(D)]<120.                              (14)

To verify it set L=1+log(.01/s). The lower expression for a_x exp(D)
is .01 exp(.055(1-L)-.25)L, with
1+log 10<=L<81-log 100. Its logarithm is concave (second derivative
-1/L^2), so its minimum is at an endpoint. Both endpoint values exceed
1/120. This explicitly uses h>exp(-80); it is not an extremal-limit
claim or a unit renormalization.

Equations (12)-(14) establish the previously unknown far-seam prefactors

    T0=2.4e21, T1=1.2e27

in epsilon a_x exp[-(M/2-1)D] and epsilon a_x M times that exponential.
The preceding value threshold is 6.8463423100324883e22. Thus the proved
value prefactor fits with a factor greater than 28. Near the bolt the
gradient has much more room. In fact both rows fit the stronger allocation

    |delta H_M|<=1e20 epsilon a_x M^-3,
    |grad delta H_M|<=1e20 epsilon a_x M^-3, chi<=1.   (15)

For value maximize M^3 exp[-(M/2-1)D] at M=18,D=1.9. For gradient use
M^4 and D=60; the maximum is again at M=18. The logarithmic derivatives
3/M-D/2 and 4/M-D/2 are strictly negative in their stated ranges.
These are analytic maxima over all M, not sampled angular cutoffs.

The exact full counting estimate N(E)<=4CE^2 and M1>=18 imply by
Stieltjes integration sum d_n/M_n^3<=12C/18=2C/3. Thus the full massive
tower in (15), including every angular multiplicity, is absolutely and
locally uniformly summable together with its near-bolt first derivative.
Multiplication by the original extra C gives the physical Theta units.
Using a_x<=2hL_h and allowing the fixed proper-to-chi factor, sufficient
canonical normalized far-source bounds are

    |delta S_far|<=1e21(2 eta_G)epsilon a_x,
    |delta S_far(chi)-delta S_far(0)|
                    <=1e21(2 eta_G)epsilon hL_h chi,
                    0<=chi<=1.                       (16)

The distant state/boundary correction therefore spends at most one
percent of the proposed cold value allowance and less of its square-root
oscillation allowance. This establishes this contribution for the explicit
family without assuming the complete K_tr theorem.

## 6. Consequence, limitations and next calculation

The actual-domain far seam is no longer an unknown prefactor in the
finite-centered admission audit. It is controlled with the true inverse,
fixed physical cutoff, explicit cap, finite h range and complete tower.
The local flat-to-H2-cap term and its normalized endpoint interpolation
remain open. Their target is C_glue R^2<=1.8161835325509682e16, or
R<=1347.65853707494 when C_glue<=1e10. The complete local Born/restoration
and mixed weighted coefficient contract must still be verified for the new
family before claiming the real endpoint theorem or cold admission.

Next address that local COMPLETE-H2 problem. A transvection there acts
on a complete domain, unlike the finite physical disk. Retain the local
spin/frame connection and all row derivatives, or derive its weighted
response in a fixed chart. Do not repeat the already closed distant
boundary calculation or reopen the failed uniform source/inverse product.

The complete solver X, embeddings, actual correction/moment, nonlinear
contraction, exact compatibility, metric boundary/action, semiclassical
frequency validity and on-shell stability are not proved by (16).
All S167 acceptance tasks remain open and the full objective is unchanged.

## 7. Frozen verification

Freeze this note, checker and review before execution. Check the complex
cone, local energy/evaluation products, distance/weight extrema, per-mode
budgets and full counting identity. Independently verify the Fourier
point norm and constant-principal Dirac identities. In finite skew-adjoint
models compare the general two-sided identity, all parameter product
derivatives, complex Cauchy and real finite differences. A lattice
commutator C=[A,chi] need not commute with chi: retain its additional
term (1/2)S_cap[C,chi]S_cap. The first frozen run exposed precisely this
omission in the manufactured test; that failed run is retained. Check the
continuum Leibniz identity separately: C=gamma.grad chi is multiplication,
so [C,chi]=0 and (10) applies. Preserve both seam commutators and the
diagonal-density variation. Use 60/90 digits and normalized errors.
These fixtures are not physical states or continuum proof substitutes.
Pin inherited inputs, preserve every run and invoke no Lean.

Remaining to-do list: local cold quantum-source proof and admission;
actual correction in admitted norms and signed moment; complete nonlinear
map, projected solution and exact compatibility; physical boundary/state,
validity, matching, stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
