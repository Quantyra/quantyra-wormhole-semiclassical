# Whole heat remainder: boundary cancellation and the bolt transfer gate

2026-10-23, E003/S169. Ordinary argument, frozen before verification.
The preceding turn was progress. The full stable-wormhole objective is
active and incomplete. The closed feasibility audit is not reopened.
This is a method/domain test, NOT execution of the complete signed screen.

## 1. The missing estimate and what the older constants actually bound

Retain the selected-pair identity and local coefficient of
[the smearing argument](current-smeared-heat-projection-domain.md):

    delta_h delta_v K(t)=-2t partial_t delta_v Tr[h exp(-tL)],
    r_hv^prod(t)=a_hv t^2+E_hv(t), v=-z^4,
    |a_hv|<=2.3454197e-13.

The common product-gradient contact, both square partners, original HH
period, complete angular tower and actual mass/log h=h1 are retained.
The required numerical estimate is still

    |E_hv(t)|<=B_heat t^3, 0<t<=.001,
    B_heat<=2.1332863e-5.                         (1)

The preceding qualitative half-domain argument does not supply B_heat.
In particular, the 1e15 t^2 short-time bound in
stationary-reference-nonbolt-stress.md is for a different reference
quantity. It is NOT a bound on a metric derivative of that quantity.
An inequality cannot be differentiated to obtain one. Its auxiliary
caps and finite-cut comparisons must also be differentiated and matched.
The 6e24 resolvent-gluing coefficient is an excessive KNOWN allowance
in its stated norm; it is not a numerical B_heat. No conversion between
these objects has been proved. None is a lower bound on the source.

We first test whether the boundary and bolt prevent a useful ordinary
estimate, before choosing a numerical solver for the remaining trace.

## 2. Exact paired boundary cancellation on the actual half-domain

Use inward proper distance Y from the outer seam and normalize the
time warping to one there. Write a(Y)=sqrt(Q(z)/beta), f=1/a,
z'=sqrt(Q), z(0)=0. The seed has Q'(0)=0. On the actual paired domain,
the half-density square has differential expression

    -partial_Y^2-f^2 partial_T^2
                       -gamma_Y gamma_T f' partial_T.          (2)

The two complementary Dirichlet/Neumann projectors have equal rank.
Reflection interchanges the two scalar potentials k^2 f^2 +/- k f'.
The image term is off diagonal in this spin decomposition and its trace
vanishes. Thus the paired scalar-smeared trace equals the restriction of
the diagonal trace for the EVEN reflected coefficient. This is the exact
operator identity established in stationary-massive-boundary-layer.md
and stationary-reference-nonbolt-stress.md, not a smooth-double theorem.
It holds for every thermal frequency before summation. The reflected
metric remains C2,1. A smooth analytic continuation is used below ONLY
to identify its difference from a local interior parametrix; this
difference is retained. No physical domain is replaced by that continuation.

At the seam let Q=b^2+u z^2+c3 z^3+c4 z^4+c5 z^5+c6 z^6, b=sqrt(beta).
On the genuine smooth family g_s=exp(2s v)g0, set
V(Y)=int_0^Y v(z(x))dx. At fixed proper coordinate Y the warping variation is

    dot f=-v f-f' V.                              (3)

This follows from Y_s=Y+s V+O(s^2) and
f_s(Y_s)=exp(-s v(Y))f(Y). It is essential to include the coordinate
variation; varying f by -v f alone would give the wrong odd jets.
The ODE z''=Q'(z)/2 gives, as identities of actual seam jets,

    dot f=b^4 Y^4-(u b^4/30)Y^6
                    -(3 c3 b^5/10)Y^7+O(Y^8).     (4)

In particular dot f has NO odd jets of weights 1,3,5. The first possible
odd variation is weight seven. This is an exact jet statement, not an
approximation of the metric by a polynomial in Y. The actual coefficients
c3=-4beta-2u-tau and c4,c5,c6 are retained in the execution.

Here is the consequence for the ENTIRE boundary remainder's order.
For an even germ the paired trace has no boundary correction to its
interior continuation. Every boundary-difference coefficient must
therefore contain an odd metric jet. The seed's first such jet has
weight three. The variation's first even jet has weight four and its
first odd jet has weight seven, by (4). A differentiated boundary
difference cannot occur below weight seven. In a diagonal heat kernel
of dimension two that is t^(7/2-1)=t^(5/2).

At fixed z one must ALSO differentiate the observation coordinate and
volume density. The shift is V(Y)=O(Y^5), multiplying a spatial derivative
of the seed's weight-three profile; the density variation is O(Y^4),
multiplying that profile. Both have weight seven. Neither term may be
dropped. Applying -2t partial_t preserves this order. Smooth half-domain
parametrix remainder estimates, with the original boundary operator,
give for some FINITE constants, uniformly for 0<t<=.001,

    |R_boundary(t,Y)| <= C t^(5/2) (1+Y/sqrt(t))^N
                                     exp[-Y^2/(c t)].           (5)

The constants in (5) are not quantified here. Nonlocal thermal images
and cut-off/gluing errors are exponentially small at t=0 and are INCLUDED
in the full remainder constant, not set to zero. Positive-time continuity
extends the small-time estimate to the stated compact time interval.
The regular bolt is an interior point in this heat construction.

Integrating (5) over the seam collar gives O(t^3) WITHOUT using a numerical
coefficient in |h|<=H z^3. The known global |h|<1.2e-8 is enough. This
strengthens the selected smearing argument; it does not prove a C2 action
on a neighborhood of rough base metrics. The common heat regulator can
still be removed for the selected response using Theta(t)<=4C/t.

The smooth half-domain heat expansion and local jet weights used in this
ordinary argument are described in Vassilevich, sections 2 and 5:
https://arxiv.org/html/hep-th/0306138v3. They provide qualitative orders,
not the missing numerical constants. The checks below verify (3)-(4)
and the budget arithmetic; they do not prove (5).

## 3. A numerical whole-remainder target with no unknown mass-tail constant

Write the complete time-reduced response after the retained local
subtractions and leading coefficient as

    E_hv(t)=int_0^zh h(z) R_v(t,z) dz.

Here R_v includes the product contact. The local integrations by parts
are the ones already justified for the actual mass/log phase; their
outer and bolt fluxes vanish. R_v includes the ENTIRE actual heat
response, not just its boundary expansion or finitely many modes.
Define an intrinsic, absolute-diagonal norm in the fixed z measure:

    D_heat=sup_(0<t<=.001) t^-3 int_0^zh |R_v(t,z)| dz.          (6)

The interior expansion after its t^2 term is O(t^3). Equation (5) makes
its boundary part integrable with the same order. Thus D_heat is FINITE
for the selected smooth test and actual seed. No numerical value follows.
The previously established global correction bound gives the exact
sufficient implication

    B_heat <= (1.2e-8) D_heat,
    D_heat <= 1777.7385 => B_heat < 2.1332863e-5.                (7)

The guard 1777.7385 is rounded DOWN. It is not a rescaling improvement:
the new content is cancellation of the lower boundary profiles, which
makes the absolute-diagonal norm finite and removes the unknown H.
In the older weight-four estimate the resulting upper allowance for
that norm diverges as t^-3/2, so multiplying that estimate by the global
correction bound would not have been valid. The full numerical value in
(6), including interior, boundary, thermal and finite-geometry errors,
remains UNKNOWN. Finiteness is not evidence that it meets (7).

## 4. Execute the bolt transfer test before importing cylinder estimates

Let rho be proper distance FROM the bolt, z=zh-s(rho), and
Q(zh)=0, Q'(zh)=-2gamma, gamma>0. Put Qj=Q^(j)(zh).
The exact ODE s''=-Q'(zh-s)/2 gives

    s=gamma rho^2/2+gamma Q2 rho^4/48
                      +(gamma Q2^2/2880-gamma^2 Q3/480)rho^6
                      +O(rho^8),
    A(rho)=sqrt(Q)/gamma=s'/gamma
           =rho+(Q2/12)rho^3
                    +(Q2^2/480-gamma Q3/80)rho^5+O(rho^7).    (8)

In theta=gamma tau, with the ORIGINAL period 2pi, the metric is
drho^2+A(rho)^2 dtheta^2. A/rho is smooth and even, with value one.
Its Cartesian metric extends smoothly across the bolt. There is no
conical defect. This statement does not change the period or replace
the finite geometry by its product approximation.

The reciprocal warping used by the radial cylinder method nevertheless
satisfies

    f=1/(gamma A),
    partial_rho log f=-1/rho-(Q2/6)rho+O(rho^3).                (9)

Thus sup |f'/f| is INFINITE on the actual full domain, for every positive
gamma, including the retained gamma approximately 5.0578e-34. No finite
L in the cylinder Feynman--Kac guard |f'/f|<=L reaches the bolt. Taking
gamma small does not cure (9); its singular leading coefficient is one.
This is an ACTUAL obstruction to that direct transfer of estimates,
not an obstruction to a heat kernel or a wormhole.

The complete thermal half-integer family p_j=gamma(j+1/2), j>=0,
gives the two radial indicial potentials

    j^2-1/4, (j+1)^2-1/4, divided by rho^2.                    (10)

The regular half-density powers are rho^(j+1/2) and rho^(j+3/2).
For j=0 the first channel has repeated exponent 1/2; its second local
solution is sqrt(rho) log(rho). Smooth disk spinors select the log-free
branch. Square integrability alone would not select it. This calculation
uses symbolic j for the full family; the critical j=0 case identifies a
domain issue and is not a finite-mode replacement for the determinant.
The nonnegative paired square and its regular disk extension are retained.

Consequently a calculation on a complete auxiliary cylinder must still
prove the metric-differentiated comparison to this regular disk domain.
The scalar stress comparison in the older reference note does not supply
that comparison for the current mixed response. This is a missing
CONTINUUM transfer estimate, not a rounding error or a failed numerical
quadrature. The earlier smooth half-domain argument already treats the
bolt correctly and is not invalidated by this coordinate test.

## 5. Decision and one next calculation

Executed in this turn: the exact boundary variation jets, their first
possible odd weight, the bolt regularity/indicial transfer test and the
numerical conversion (7). These tests distinguish a real obstruction to
the direct cylinder proof from an unknown intrinsic heat estimate.
They do NOT evaluate (6), produce B_heat, or execute the signed projection.
There is no evidence here for a physical nonexistence claim or for
changing the specified construction. There is also no certified numerical
basis yet for saying the full remainder target is realistically attained.

Research decision: retain the construction and geometric smearing
definition; sharpen the method to an actual-domain heat-defect estimate
in coordinates regular at the bolt, with the exact paired seam response.
Park attempts to import an undifferentiated reference/cylinder constant
or to infer the residual from additional local coefficients.

The SINGLE next calculation is a certified evaluation/bound of (6)
against D_heat<=1777.7385. Use a boundary-compatible parametrix P_s on
the actual smooth half-family and its full heat-equation defect
F_s=(partial_t+L_s)P_s. It must have the exact initial identity and
paired boundary condition, with a regular Cartesian bolt patch.
The ordinary error identity is

    H_s-P_s=-int_0^t H_s(t-u)F_s(u)du.

Differentiate this identity in s, including BOTH H'_s F_s and H_s F'_s,
and apply -2t partial_t before bounding its time-reduced diagonal L1
norm. Include the known parametrix remainder and all gluing/boundary
defects. A trace-norm bound per original time period is a sufficient
upper bound for this diagonal norm. A numerical defect alone, without
its continuum propagation bound, is insufficient. This calculation is
PROPOSED, NOT EXECUTED; its implementable parametrix and numerical norm
certificate are still missing.

Only after that gate passes should the remaining positive-time signed
projection be evaluated with angular/thermal tails and fixed contacts.
The user screen |P_hat_rest|+E_rest<1.287e-40 remains UNEXECUTED and
sufficient only. The leading heat allowance <5.4970775e-48, magnetic
<1.284506e-40, complete cold <1.693010e-43 and scalar <2.331e-44 remain
unchanged. No cold/scalar reoptimization or total-source gain is claimed.

No actual controlled solution has been established. At rho_solution=1e-8,

    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1

remain unproved, including the numerical physical inverse/domain
estimates. The full C2 action on the rough nonlinear phase, exact
compatibility, physical matching and causal boundary realization,
semiclassical validity/noise and all-sector stability are unresolved.

Remaining to-do list: certify the entire actual-domain heat remainder
and complete signed quantum/contact projection; forward/adjoint control
and nonlinear admission; rough-base action/domain and exact compatibility;
physical matching, causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
