# Coincident massive response in fixed interior normal coordinates

2026-10-19, E003/S144 supporting S139/S074. Ordinary argument and
verification protocol before execution.

This proves a finite-background, coincident point-source estimate on
a specified normal-coordinate sector of the actual paired exterior.
Unlike the remote estimate, the perturbations may overlap. The normal
coordinate restriction is essential: their values and first derivatives
vanish at the observation point. Returning arbitrary perturbations to
this sector, controlling the outer collar and obtaining the full
exterior inverse norm remain separate tasks.

## 1. Geometry, directions and claim

Retain the prepared period, magnetic bundle, paired representative,
common scale and normalization of
[S143](curved-remote-dirac-response.md). Thus

    M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18,
    eta_G=GC/B^2<8.5e-33, Gamma_hat=(2G/B^2) sum_n d_n Gamma_Mn.

The orbit trace has rank two on one physical half. Pointwise traced
diagonals agree with the even-double convention; there is no additional
half or angular multiplicity in a point source.

Fix an interior observation x=0 in one of the cold or compact-middle
normal charts of
[the canonical response](stationary-centered-canonical-response.md)
and [the curved form](stationary-curved-quantum-tt.md), section 3.
The chart must contain |z|<1/32 and stay before the actual outer
boundary. The seed metric is exp(2 sigma_0)|dz|^2, with
sigma_0(0)=grad sigma_0(0)=0. In this patch its analytic normal mass
potential has the retained complex bidisk bound
|exp(sigma_0)-1|<.026 at radius 1/32 and the real quadratic/cubic
bounds 8r^2 and 1e4 r^3. This is an interior hypothesis, not a
regularity claim across the cubic outer seam.

Choose a fixed smooth radial cap equal to one through 1/64 and zero
by 1/32. Cap sigma_0 and the scalar directions with this same cap.
The cap metric is Euclidean outside that disk. It is only an auxiliary
operator; the actual-domain difference will be bounded, not omitted.

For a direction u define U as the maximum of its actual global
supremum and gbar-gradient norms, the corresponding two cap norms,
and all real coordinate derivatives of orders 0 through 5 in the
uncapped chart. Define W similarly for w. Zero directions are
understood by continuity. Directions preserve the identical paired
domain and even seam condition. At x impose

    u(0)=w(0)=0, grad u(0)=grad w(0)=0.                 (1)

The real background b has the same zero first jet and norm at most
1e-3 in this norm. No analytic spatial extension of b,u,w is required.
Use the complex PARAMETERS

    phi = sigma_0+b+s u+t w,
    |s|<=1/(40U), |t|<=1/(40W).                       (2)

The same parameters act on the actual exterior and its local cap.
This does not complexify a real spatial point. Globally the real and
imaginary perturbation sizes, and the gradients needed below, are
bounded in modulus by .051<.06.

Let F_hat[b](0) denote the canonical massive action source per fixed
normal-chart area, with the common local allocation. The theorem is

    |D^2 F_hat[b](0)[u,w]| < 7e-9 U W.                (3)

It holds uniformly for this background ball, with all massive angular
levels and the actual outer domain retained. It excludes the independent
four-dimensional Weyl functional and the LLL. It is not a bound on
the whole coupled metric operator.

On any real segment b+tau w remaining in the ball, integration gives

    |D F_hat[b+w](0)[u]-D F_hat[b](0)[u]| < 7e-9 U W.  (4)

This is a finite-background statement within (1); no changing normalizer
has been differentiated.

## 2. A uniform local coefficient cone

Write h=exp(phi_cap)-1. Throughout the parameter bidisk (2) use

    h(0)=Dh(0)=0, ||h||_infinity<.03,
    |h(z)|<=9r^2, |T2 h(z)|<=9r^2,
    |h-T2 h|<=2e4 r^3, |h-T4 h|<=1e12 r^5.          (5)

Tj is the Taylor polynomial of degree j at zero, extended as a
polynomial to the whole plane. The inequalities in (5) are global
for the capped function.

Here are sufficient quantitative checks for (5). On any real ray,
the seed's analytic bidisk bound implies Taylor coefficient bounds
.026*32^j and, for r<=1/64, fifth remainder coefficient at most
.052*32^5<1.75e6. The stronger retained quadratic and cubic bounds
are used for degrees two and three.

On this inner disk, let zeta=b+s u+t w. Every coordinate derivative
through order five is bounded by .051. On a unit ray the bound is
.051*(sqrt(2))^j, and zeta(0)=zeta'(0)=0. In particular
|zeta|<=.051 r^2. Its exponential has value below 1.001 there.
The product rule, or the elementary exponential Bell polynomials
through degree five, bounds every ray derivative through degree five
of (1+h_0)(exp(zeta)-1) by 1e8: the seed derivative reserve
.026*j!*64^j for j<=5 and |zeta^(j)|<1 suffice. Taylor's integral
remainder therefore costs less than 1e8/5! in degree five.
For degree three, use instead
h_0=O(8r^2), zeta=O(.051r^2), and the cubic Taylor remainder of
zeta, at most .051*(sqrt(2))^3 r^3/6. These bounds enlarge the
seed cubic constant by less than one. The quadratic coefficient
increases by less than .052.

Outside the inner disk, capping leaves |phi_cap|<.028 and the same
quadratic bound. The triangle inequality with Taylor coefficient
reserves |h_2|<=9r^2, |h_3|<=2e4r^3, |h_4|<=1e6r^4 gives

    |h-T2h|/r^3 <=18*64<2e4,
    |h-T4h|/r^5 <=.03*64^5+9*64^3+2e4*64^2+1e6*64
                   <1e12.

The fourth coefficient reserve follows from the preceding inner
derivative bounds, more sharply from the degree-four product rule
and .026*32^4<3e4. Thus (5) holds also for the nonanalytic, real-C5
background. The cap's large higher derivatives do not enter these
global Taylor-division estimates. Its value and gradient do enter U,W.

## 3. Complete local series and the local subtractions

Use Dflat=gamma.partial, S=(Dflat+M)^-1. The nonlocal action is
-Tr log(Dflat+M exp(phi)). At the normalized center the source
difference from the flat vacuum is

    F_M,cap(0)=-M tr[(Dflat+M+Mh)^-1-S](0,0)
      = sum_(N>=1) (-1)^(N+1) M^(N+1)
          tr[S h S ... h S](0,0),                    (6)

with N insertions of h. The first two local polynomials are

    F_local = -A/(12pi) -(B-3A^2)/(120pi M^2),
    A=Delta phi(0), B=Delta^2 phi(0).                 (7)

The minus signs in (7) are action-source signs. The retained canonical
TRACE convention in the earlier source note has the opposite sign.
Confusing these would reverse the restored curvature contact.

For completeness the first Born term gives

    -A/(12pi) -(B+A^2+2 tr H^2)/(120pi M^2),
    H=Hess phi(0).

Its radial integrals are
integral t^3(K0^2-K1^2)dt=-1/3 and
integral t^5(K0^2-K1^2)dt=-8/15.
The quadratic part of the second Born term gives
(2A^2+tr H^2)/(60pi M^2). The Hessian-anisotropy terms cancel.
An independent momentum representation of this second term is

    -(M^3/4) integral d^2k/(2pi)^2
        tr[(H:partial_k partial_k S) S
           (H:partial_k partial_k S)].               (8)

Two integrations by parts produce (8). The exact polynomial integral
in the checker verifies the displayed tensor coefficient, including
an off-diagonal Hessian. The first Born cubic Taylor polynomial has
zero angular integral. Thus (7) subtracts exactly the terms through
four derivatives, including the first-Born quartic and second-Born
quadratic product.

These signs also follow by varying
-integral sqrt(g) R^2/(960pi M^2): its source at a normal point is
-(B-3A^2)/(120pi M^2). Equation (7) is a subtraction for estimating
the remainder; the common local action is restored in section 6.

## 4. An absolute all-orders coincident remainder

The positive kernel bound is

    ||S(x)|| <= k_M(r)=(sqrt(2)/pi) exp(-Mr/sqrt(2))/r.

It follows from K0<=K1 and the real integral representation
[DLMF 10.32.9](https://dlmf.nist.gov/10.32#E9):
K1(s)=integral_0^infinity exp(-s sqrt(1+t^2))dt
<=sqrt(2) exp(-s/sqrt(2))/s. The primary representation was
inspected on 2026-10-19.

Put k_*=exp(Mr/10)k_M. Direct radial integration gives

    ||k_*||_1<5/M, ||k_*||_(3/2)^3<13/M.

For an (N+1)-edge loop, N>=2, Young's inequality for the first three
edges and L1 bounds for the others imply

    integral l^p product_edges k_M
      <=13 p!10^p 5^(N-2)/M^(N+p-1),                 (9)

where l is the total edge length. Use l^p<=p!(10/M)^p exp(Ml/10).
The trace costs at most two. For the two-edge loop, direct radial
integration instead gives

    2 integral d^2z k_M(r)^2 r^5 <16/M^5.             (10)

Subtract (7). The first-Born remainder costs at most
16*1e12 M^-3, using (5),(10). For the second Born term,

    |h(z)h(y)-h_2(z)h_2(y)| <=2*9*2e4 l^5.

Thus its remainder costs at most
(2*13*5!*10^5)(2*9*2e4) M^-3.
For every N>=3, bound the first three insertions by 9^3 l^6
and the others by .03^(N-3). Equation (9) sums them to

    [26*6!*10^6*9^3*5/(1-5*.03)] M^-4.

Together, for M>=18, these three terms are strictly below

    |E_M,cap| := |F_M,cap-F_local| <2e14 M^-3.        (11)

This is a complete resolvent series, not an asymptotic truncation.
The operator Neumann series converges since ||M S h||<.03;
the displayed moment bounds justify its coincident, subtracted
evaluation and locally uniform holomorphy in (2).

Cauchy's formula in the TWO parameter disks gives

    |D^2 E_M,cap[u,w]| <3.2e17 U W M^-3.             (12)

There is no extra 2!: each of the two independent parameters is
differentiated once. The subtraction in (7) is holomorphic as well.

## 5. Actual domain versus the local cap

Use a cutoff chi equal to one through radius 1/128 and zero by
1/64, with |grad_gbar chi|<256. The two operators and parameter
families agree on its support. A local geometric-resolvent identity
therefore expresses their kernel difference at (0,0) as

    R_cap [D,chi] R_actual,

up to a sign irrelevant to the norm. It is an identity between
operators on their own global spaces, restricted to the common chart.
The middle multiplication is supported at distance
at least exp(-.05)/128>1/256 from zero.

Both global resolvents obey the complex weighted bound 3/M at
weight M/2. Indeed

    Re exp(b+s u+t w) >= exp(-.06) cos(.06),
    exp(-.06)cos(.06)-1/2 >1/3.

The same holds for the cap. This uses actual anti-self-adjointness
and scalar domain-preserving multiplication. Reverse rows are
controlled by the conjugate-parameter adjoint, so no real-only
grading identity is silently used at complex parameters.

For local evaluation, use the two cutoff balls and argument in S143,
now without its unnecessary factor-100 reserve. On a ball of radius
1/(10M), the flattened mass has modulus below 1.13M and gradient
below 12M. Its squared Dirac potential is bounded by
(1.13^2+12/18)M^2<3M^2. The graph estimate is still 30000M^2,
and the Fourier point-evaluation norm and density factor give

    |psi(0)| <1e4 M ||psi||_local.

The resulting separated row bound is
3.2e4 exp(-M distance/2). The evaluation ball fits before the
commutator support even at M=18. Rank two and the cutoff gradient
therefore bound the SOURCE difference, including the exterior
boundary, by

    |F_M,actual-F_M,cap|
      <=2*(3.2e4)^2*256 M exp(-M/256)
      <6e11 M exp(-M/256).                           (13)

Their local subtractions are identical. Cauchy's formula gives

    |D^2(E_M,actual-E_M,cap)[u,w]|
       <9.6e14 U W M exp(-M/256).                    (14)

No scalar source has been truncated at the cap boundary. The actual
regular bolt, thermal period, paired domain and global resolvent
remain in (13). This comparison applies to interior charts; it is
not a boundary-chart theorem.

## 6. Full angular sum and common local restoration

The exact count N(E)<=4CE^2 and M_1>18 give

    sum_n d_n/M_n^3 <12C/18=2C/3.

For any delta>0, integration by parts with E exp(-delta E) gives

    sum_n d_n M_n exp(-delta M_n) <=24 C/delta^3.      (15)

To justify the sign in this step, discard the negative contribution
from E<1/delta, then bound the positive integrand on E>=1/delta
by 4C delta E^3 exp(-delta E). Do not replace N(E) by its bound
where the integration-by-parts weight is negative.

Equations (12),(14),(15), with delta=1/256, show that the normalized
massive remainder has second source variation below

    2 eta_G [(2/3)*3.2e17 +9.6e14*24*256^3] U W
      <6.58e-9 U W.                                 (16)

This is absolute convergence, uniform on the parameter bidisk, so
the angular sum, two variations and finite real path integral commute.

Restore the same common a2 functional as S142:

    Gamma_a2=alpha integral exp(-2phi)(-2Delta phi)^2,
    alpha=-C F_eff/(960pi), |F_eff|<100.

At the normal point its action source is 8alpha(B-3A^2).
Under (1) its mixed second variation is exactly

    -48alpha Delta u(0) Delta w(0).                  (17)

Since |Delta u|<=2U and |Delta w|<=2W, the normalized contact
cost is below 14 eta_G U W<1.2e-31 U W. The a1/Polyakov source
is linear in the conformal factor at the point; its mixed second
variation vanishes here. Vacuum source terms depend on the point
value, also unchanged by (1). These statements concern this mixed
variation only; the local terms stay in the action.

Combining (16),(17) proves (3). The fourth derivative contact was
restored at the common scale, not summed as an unsigned divergent
inverse-mass term.

## 7. Verification and precise remaining task

Before execution, require exact agreement of the one- and two-Born
local coefficients, direct curvature-action variation, loop-integral
constants, coefficient-cone reserves, angular integrals and the final
rational budget. Independently integrate (8) by radial/angular
quadrature for three real Hessians, tolerance 1e-10+1e-8*abs(value).
Use two resolutions; each must meet that tolerance.

A finite graded-matrix diagnostic compares the full resolvent source
with its convergent Born series and the mixed source variation with
central finite differences at steps 1e-4 and 5e-5. Require absolute
error below 1e-8 and a refinement ratio below .35, unless both
errors are below 1e-45. A finite real-background source-derivative
change must agree with quadrature of the mixed variation to 1e-30.
These diagnostics do not certify a continuum domain or chart.

The argument supplies overlapping, real-C5, finite-background source
control on a nonempty interior normal-coordinate sector. The local
coefficient proof replaces a formal heat asymptotic series with an
absolutely convergent resolvent construction. The boundary effect is
bounded separately using the actual operator. No independent reviewer
or formal verification is claimed.

It remains necessary to restore arbitrary value and first-jet
directions through the full metric-coordinate chain rule, and to prove
the outer weighted estimate at finite backgrounds. A bound on a point
source in these C5 norms is not the bound on the combined exterior
energy or the cold cosh-weighted inverse. The scalar signed source,
radius/shape coupling, Weyl and LLL terms and physical metric return
are not resolved by (3).

Remaining to-do list: general metric-coordinate and outer-norm transfer;
complete nonlinear exterior/state and signed-source solve; driven
geometry, opening, stability, transport, preparation and separate creation.
