# Auxiliary metric return: quantitative feasibility and a solver decision

2026-10-22, E003/S169. Ordinary argument before frozen verification.
The full stable semiclassical wormhole objective remains incomplete.

Decision: retain the specified c_R=10^22 construction and the exact
auxiliary reformulation, but require a source-specific solve in the
original radius/conformal graph and admission norms. A bound N2<=1 is
useful for the scalar row ONLY. The executed return calculation below
shows that using N2 alone to certify even the radius VALUE would need
N2 below approximately 3.4e-29. No positive finite N2 bound implies the
required raw C5 bound. Thus evaluating just N2 is not the next complete
milestone. Neither conclusion is nonexistence or a lower bound on the
actual correction. No new full metric response is claimed here.

## 1. Physical variables and an exact cancellation

Use the exact product-centered finite seed of
[the product audit](recentered-product-and-cold-residual.md), its actual
fixed-period disk, and [the anomaly-complete action](anomaly-complete-auxiliary-response.md).
Write B for its physical product sphere, r0=B t, s=t-1,
tau=s/t, h=(r_plus-B)/B, and w=[(s+h)/t]^(2/5). All centers and
weights are fixed during variations. B is not a newly adjustable flux
or matching parameter. The retained uniform bounds include

    h_min=1/(8*2000^10)<=h<2000^-10,
    L_h=1+log(.01/h)>64,
    a_s=s[1+log(.01/s)], h<=s<=.001.

In the common conformal disk write the physical metric as

    g_4=(r/B)^2 A^2 exp(2 sigma) g_disk+r^2 dOmega^2,
    r=B+(r0-B)exp(V), b=sigma-sigma0.

The auxiliary frame is g_E=exp(2u)g, u=u_b+v. Its sphere radius
is r_E=exp(u)r. Let ell_E=log(r_E/r_E,0). EXACTLY,

    zeta=ell_E-v=log(r/r0),
    V=log[1+(exp(zeta)-1)/tau].                         (1)

Moreover (B/r_E)^2 g_E,2=(B/r)^2 g_2. Thus the increment b of
the product-frame conformal field is independent of v, even for finite
changes. A different constant canonical sphere in the auxiliary frame
only changes the background constant, not this cancellation. The usual
frame transformation keeps matter coupled to the original metric;
see [De Felice and Tsujikawa, section 2.3](https://arxiv.org/html/1002.4928v2#S2.SS3).
Here the quantum functional and its paired domain remain exactly as in
the preceding note; this identity does not replace them by Gamma_D[g_E].

For a spherical diagonal tangent in this fixed conformal disk let
h_E=diag(a,a,c,c) in the orthonormal frame. Then

    delta b=(a-c)/2,
    delta V=(c/2-v)/tau.                              (2)

Both |a-c| and |c| are <=|h_E| in the Hilbert--Schmidt norm.
A radial diffeomorphism used to reach this disk must be included in h_E
before applying (2). A result in areal or TT gauge is not already this
input. In particular, a small lapse alone is not a small disk correction.

For finite increments a useful exact sufficient test is

    sup |zeta|/tau <=1-exp(-rho)  ==>  sup |V|<=rho.    (3)

To prove it, put d=1-exp(-rho). For 0<tau<1,
exp(-tau*d)>=1-tau*d and exp(tau*d)<=1+tau(exp(d)-1)
<=1+tau(exp(rho)-1). Insert these bounds into (1). The logarithm's
argument stays positive. Monotone radius still requires
1+s partial_t V>0; (3) alone does not impose that derivative condition.
Equation (3) is a geometric statement, not a nonlinear extension of
the preceding LINEAR scalar Green estimate.

## 2. Substitute the available scalar and thermal constants

The previous scalar theorem, for the actual product-centered seed and
physical fast-phase input, gives

    M=N2(h_E),
    ||v||_w <= a0+b0 M,
    a0=1.51e-28/(3.2e8)=4.71875e-37,
    b0=1e-27/(3.2e8)=3.125e-36.                       (4)

N2 uses weighted covariant derivatives through order two with lambda
units. Its actual value on the metric response is UNKNOWN. The seed
values |u_b|<5.4e-37, |grad u_b|/sqrt(lambda)<1e-35 and
|Hess u_b|/lambda<1e-33 have already been included in (4). They are
background values, not the inverse response. The scalar principal
perturbation is <1e-32 and its positive zeroth order part >3.3e8 lambda;
the weaker 3.2e8 denominator includes the barrier losses.

The return has two distinct weight factors:

    W_c=sup_cold w/a_s <=2^(2/5)/(64 h_min^(3/5))<4.7e18,
    W_r=sup_(t>=1+h) w/tau
        <=2^(2/5)(1+1/h_min)^(3/5)<3e20.               (5)

The W_c monotonicity is proved in the preceding scalar cold-response
note. For W_r, the logarithmic derivative of
(1+s)^(3/5)(s+h)^(2/5)/s has numerator
-(3s+2hs+5h)/5 over the positive denominator s(1+s)(s+h).
Its maximum is therefore at the bolt. Its bolt value decreases with h,
which proves the displayed exact uniform expression. No extremal limit
or change of thermal ensemble has been taken.

Equations (2)-(5) give only the following VALUE certificates:

    sup_cold |delta b|/a_s <=(W_c/2) M,
    sup |delta V| <=W_r[(1/2+b0)M+a0].                 (6)

For a component target y=5e-9 the exact sufficient thresholds using
these envelopes are

    M <= 2y/W_c,
    M <= [y/W_r-a0]/(1/2+b0).                         (7)

The second threshold is approximately 3.4e-29; the rounded W_r=3e20
threshold is 3.333333238958333...e-29. M<=3.3e-29 suffices for both
value tests. These are thresholds for THIS upper certificate, not
necessary bounds on the true response. The finite-geometric analogue
replaces y in the numerator by 1-exp(-y), if a finite zeta estimate is
available. No such nonlinear estimate is claimed from (4).

For M<=1, the scalar contribution alone to |delta V| is <1.1e-15;
with zero metric input it is <1.42e-16. Its cold logarithmic-radius
value remains <2e-17. Thus the scalar coefficient is not the present
value bottleneck. The independent metric sphere component c/2 and the
correlation c/2-v remain unbounded in the required radius norm.

These facts explain what a useful improvement must measure: bound
c/2-v with the factor tau, and a-c with a_s, before taking absolute
values. The linear fixed-temperature constraint formula in
[the earlier reconstruction](stationary-centered-first-order-geometry.md)
already exhibits correlated moving-bolt cancellations. Its source was
the fixed-seed Gamma_rest in a different local comparison theory.
Its lapse <3.3e-10, pulled-back metric <7.5e-10, gap <3.4e-10,
disk |V|<1.2e-8 and cold b bound <1e-6 times its barrier cannot be
assigned to the present coupled auxiliary inverse. Even there those
are upper allowances, not actual response lower bounds.

## 3. Execute the remaining N2-to-admission test

The decisive derivative test is negative, with no missing constant.
Choose a compact radial annulus strictly between the bolt and infinity,
a smooth cutoff chi equal to one on a smaller annulus, and a smooth
normalized radial coordinate y there. Vary only the disk conformal field
by b_n=n^-2 chi(y) sin(n y), leaving the sphere field fixed.

These are smooth spherical physical metric directions, zero at both
ends, with the original clock and regular bolt. The frame, weight and
connection coefficients on this compact set are fixed and bounded.
The corresponding metric N2 is O(1): its j-th derivatives for j<=2
are O(n^(j-2)), plus bounded cutoff/connection terms. On the smaller
annulus the fifth raw derivative has supremum n^3 (for large n that
annulus contains a sine/cosine maximum). Fixed j!100^j chart units
multiply this by a fixed positive factor. Therefore

    sup_(physical smooth directions) raw_C5/N2 = infinity. (8)

Multiplying every b_n by any fixed small positive amplitude preserves
this conclusion and can meet any positive N2/curvature smallness bound.
Thus reducing the unknown N2 constant, even below (7), cannot by itself
close the raw-C5 admission requirement. These fields are TEST directions,
not solutions of the metric equation or a sequence of physical instabilities.
Equation (8) rejects a proposed embedding; elliptic estimates using the
actual source and a stronger graph can still control the actual solution.

Likewise, the retained reference energy estimate E<=1.44e-10 would need
an embedding <=625/18=34.722... to imply the 5e-9 target. The existing
cold test proves its embedding is >2.6e15, and the energy-to-C5 map is
also unbounded. A transfer defect <=1/4 would give energy <=1.92e-10
and require embedding <=625/24=26.0416..., not cure the mismatch.
These are certificate failures, not bounds from below on the solution.

## 4. Complete solution inequalities and missing-estimate targets

Fix a COMPLETE physical correction space X that includes b,V, their
cold value/gradient and raw C5/cap contracts, outer conormal remainder,
correlated mass column and the graph of the full quantum metric operator.
Completeness, the common paired domain on the ball and bounded return
maps have not yet been established with numerical constants in this space.
For a valid fixed projection Pi and preconditioner B define

    Y0=||B Pi F(0)||_X,
    Z0=||I-B Pi F'(0)||_(X->X),
    ||B Pi(F'(x)-F'(y))||<=Z2 ||x-y||_X.

Preservation and contraction require

    Y0+Z0 rho+(Z2/2)rho^2<=rho,
    Z0+Z2 rho<1.                                      (9)

At rho=1e-8, Y0<=5e-9, Z0<=.25, Z2<=5e7 give respectively
1e-8 and .75. For fixed Y0,Z0 the missing nonlinear constant must satisfy
Z2<=2[(1-Z0)rho-Y0]/rho^2 and Z2<(1-Z0)/rho. A negative numerator
means these upper allowances do not certify this radius.

B must be injective on the projected source range. A projected solution
is a full solution only if (I-Pi)F=0 is proved exactly. S169 does NOT
inherit the old zero-c signed projection. The current auxiliary algebraic
constraint bracket >1-2e-8 and Ward/current recovery are useful for
equivalence, but do not construct a new inverse or permit deleting a
source component. Use Pi=I unless an actual physical range condition
has been derived. If a projection is needed, compute its signed residual
and parameter transversality, not just its absolute upper bound.

For every admission component i, require b_i+E_i rho<=a_i with positive
slack; the missing embedding must meet E_i<=(a_i-b_i)/rho. Retained
contracts give the following explicit conditional targets. Their old
background estimates are not silently transferred to the new seed.

| Requirement | Target at rho=1e-8 | Status |
| --- | --- | --- |
| Interior P45 | b_45+1e-8 E_45<=.001; if b_45<=1.3e-4, E_45<=87000 | Old path allowance; new complete return not certified |
| Cold intersected N_c | b_c+1e-8 E_c<=5e-7; if E_c<=1, b_c<=4.9e-7 | Actual b_c and E_c unknown; raw C5 not supplied by N2 |
| Outer response N_o | b_o+1e-8 E_o<=1e-7; if E_o<=1, b_o<=9e-8 | Old bound 1.01e-7 fails to show positive slack, not actual nonadmission |
| Radius derivative | |s partial_t V|<1; budget its X embedding times rho plus background | Value return alone insufficient |
| Curvature constraint branch | |phi+12lambda|<=lambda, |R-phi|<=lambda | Exact branch margins must be preserved by X ball |
| Domain/inverse transfer | ||L_ref^-1(L_b-L_ref)||<=.25 in the SAME full graph/domain | Unknown; product displacement divided by h is not a transfer bound |
| Full source response | Z2<=5e7 with all quantum/radius/domain terms | Unknown; component source-value constants 1e-9, 2e-7, .02 are not Z2 |
| State/bolt/outer | common paired graph, regular even bolt, fixed clock/period, physical fast branch and mass relation | Exact conditions plus open-ball coercivity; not scalar inequalities alone |

An auxiliary block solve also needs derivative-bearing scalar space U,
not just weighted scalar values. If in these spaces a scalar inverse gives
||v||_U<=A+C||h_E||_X and a metric inverse gives
||h_E||_X<=M0+P||v||_U, then

    PC<1,  (M0+P A)/(1-PC)<=5e-9                     (10)

are explicit metric-response targets. Equivalently
P<=[5e-9-M0]/[A+5e-9 C], with positive numerator, together with PC<1.
One sufficient allocation is PC<=.25, M0<=3e-9, P A<=7.5e-10.
All four constants for these complete spaces are unknown. Substituting
a0,b0 for A,C is invalid because (4) is only a VALUE estimate. The full
original-metric return must then fit Y0, including the scalar component.

## 5. Latest large allowances and where they do and do not apply

The original zero-additional-c audit is preserved, including its latest
interior bound. It is a diagnostic against restarting that certificate,
NOT a source theorem for the present c_R=10^22 theory:

| Retained quantity | Upper allowance |
| --- | ---: |
| C_I, complete interior | 9.1591500000000000015e103 |
| C_O, complete outer | 2.365283519093731664e186 |
| C_tau | 6.87674442436215072e19 |
| York C_sigma; C_b4 | 1421.40132786776365; 4.52219350981899516e14 |
| C_F,I after York | 4.141944868545865e118 |
| C_F global after York | 1.069626977892751e201 |
| Cstart after scalar recombination; mhi | 3.123143830205975e28; 2.577600186401461e12 |

With g+(11/(360pi))g^2=2000^-10, g approximately 9.765625e-34,
the old sufficient budget is

    R[Cstart+mhi C_D4 C_F/g]<=1e57 g epsilon_c/800.

At epsilon_c=1 its low term alone requires R<=3.90857e-8 while the
proved norm lower bound is R>.97. If that low allowance were eliminated,
the global high term requires R C_D4<=4.32377e-226, or C_F<=2.38393e-25
at R=.97,C_D4=2. The interior alone misses by >1.7374e143 and the global
allowance by >4.4868e225. At epsilon_c=1e-4 the thresholds tighten by
1e-4. The direct signed target is |beta|<=1e-49 epsilon_c; the exact
moment budget relaxes the conservative one by about 18, not 143 orders.
See [the full old audit](construction-feasibility-audit.md) for exact
rational endpoints, residual allocation and the operator lower proofs.

C_D4 and the transfer UPPER constants remain unknown; the C_I,C_O and
resulting failed allowances are known upper bounds. The two categories
must not be merged. Removing the outer estimate does not fix the interior
loss. A change of units scales Y0 and rho together and leaves M^2 L2 e
invariant. Any improvement must use actual source/operator structure.

## 6. Continuum, matching and physical preflight

The quantum metric source is the full traceless response with the LLL,
massive angular tower, common subtractions, local Weyl/contact terms and
both adjoint partner domains. The local trace and its scalar Green
function do not determine it. The reference TT form lower bound >.994
and mixed relative bound <.336 are retained in their reference energy
scope. The order-four logarithmic metric response must remain in the
operator graph; a local Einstein constraint iteration does not become
a same-regularity contraction because the auxiliary scalar is second order.

The hardest unclosed continuum step is a quantitative product-centered
resolvent/source estimate in that full graph, through the actual paired
endpoint realization, followed by a differentiated fixed-disk return.
Reference form closedness, a normal Mellin multiplier, and smoothness of
a kernel away from its diagonal do not give the required global weighted
constants. Separated gluing rows need both endpoints, angular summability
after common subtraction, differentiated contacts and an integrable
majorant uniform on every real metric in the neighborhood. Path-only
or value-only results cannot supply that majorant.

The current product shifts |delta sphere^2|<5.5e-23 and orbit shift
<1.71e-29, and the finite seed's horizon residual <2e-12 in thermal
units, are useful source-value facts. They do not bound an operator
difference in the full graph or a derivative of the metric-dependent
domain. The exact R2 action/boundary repartition and phase-specific mass
momentum cancellation do not establish that open-ball transfer either.
Do not divide an absolute product shift or mass allowance by the tiny
thermal gap: keep the matched product, fixed temperature and correlated
outer mass/clock data in one calculation.

No newly established physical obstruction forces a different construction
in this audit. An earlier mathematical obstacle DOES change the order:
define and estimate the physical graph/return before claiming an N2 metric
bound closes admission. Keep the curvature-constraint recovery, exact
current/Ward compatibility and actual outer phase in that first solve.

Small curvature (G lambda<1e-33), the positive scalar mass
m^2/lambda approximately 3.39530545e8 and F_*>1-1.179e-8 do not by
themselves prove semiclassical validity or retarded stability. Omitted
terms, induced fluctuations, physical smearing/scales and all perturbation
sectors need independent bounds. The distinction between intrinsic and
induced metric fluctuations is explicit in
[Hu, Roura and Verdaguer](https://arxiv.org/abs/gr-qc/0508010).
No numerical omission/noise allowance or retarded spectral gap is available
here from which to compute a physical stability margin. Those are unknown
physical inputs, not known failing upper bounds. A static Euclidean form
cannot stand in for this Lorentzian test. Opening, transport, preparation
and separate creation remain obligations after the stationary construction.

## 7. Research decision and the next acceptance result

Pursue the sharper method: solve the full coupled product-centered
nonlocal boundary problem directly for the physical pair (b,V), or an
equivalent auxiliary graph with the exact joint return (1). Keep the
scalar block as a controlled component. The next acceptance result must
bound the actual correction in cold values/gradients AND raw C5/cap,
outer/mass norms, with an explicit residual and continuum tail enclosure.
Use Y0<=5e-9, transfer/preconditioner defect <=.25 and Z2<=5e7 as the
numerical design targets, and prove all positive admission margins above.

The decisive calculation for the proposed N2-only milestone has been
executed: (5)-(7) quantify its radius loss, (1)-(2) retain its genuine
conformal cancellation, and (8) disproves its required C5 embedding.
Do not commission another chain of generic N2 embeddings or re-evaluate
the old gigantic source allowances. A compact-domain numerical solve
may guide the source-specific estimate only with actual kernel/contact,
endpoint and tail errors retained; a checker or toy solve is not that
calculation. Overall construction feasibility remains undecided because
the actual full metric correction has not been enclosed. The full goal,
S169.2-S169.5 and all independent physical requirements stay open.

Remaining to-do list: actual source-specific physical graph/metric solve,
complete admission and domain transfer, nonlinear preserved solution and
exact matching; semiclassical validity, all-sector stability, charge/return,
opening, transport, preparation, separate creation and full-result audit.
