# Short-travel comparison and the finite-power domain test

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: NO PROGRESS on the actual coupled estimate. Its exact
tail arithmetic is retained, but no physical operator was applied there.

OUTCOME: the QUANTITATIVE coupled return remains UNRESOLVED. This changes
the auxiliary comparison, with a complete numerical form-error bound, and
establishes a fixed-seed route to finite strong powers. It does not evaluate
their constants, coupled signed centers or the actual correction. No
nonlinear admission or stability requirement is declared complete.

## 1. The first unsupported step and why the old citations do not fill it

The proposed T(H0^-1 R)^4 needs high regularity of the GLOBAL geodesic
comparison, not only of the physical Dirac Hessian. Its C2,1 even metric
has seam crossings and trapped trajectories. The established strong
interior remainder is H3->L2 with coefficient 1e48 C; it does not supply
its differentiated output bounds. Repeating that bound would stop at
the same output order. It cannot justify a strong-power bootstrap.

The actual fixed-seed Hessian H has a separate weighted elliptic and
original-row phase argument. The all-order source/graph argument uses
scaled smooth HALF charts, the actual paired boundary kernel and completed
contacts; it does not make the even double smooth. Its constants at higher
orders are finite and unevaluated. This argument cannot simply be applied
to arbitrary long geodesic returns, which are not Dirac resolvent kernels.

[Kirkkopelto--Manu--Salo, Definition 1.1 and Theorem 1.2](https://arxiv.org/pdf/2606.14221v1)
require simplicity and a C^(4+l,alpha) metric for their stated W^(l,p)
regularity result. The present trapped, C2,1 doubled geometry does not
satisfy those hypotheses. That theorem is not used. [Grubb, Theorem
2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf) bounds identified C^tau
symbols, with |s|<tau and N>dimension/2+1. It does not identify this global
flow compression or give its strong inverse. Its legitimate use below
is only output order zero for differentiated, normalized LOCAL symbols.

Instead remove the long geodesic term from the COMPARISON and put its
ENTIRE value in the remainder. This is a change of preconditioner, not a
change of quantum theory, state, boundary condition or subtraction scheme.

## 2. Exact short-travel operator and its complete transfer cost

Keep P=iX on the actual unit tangent bundle, the base lift J, and

    w(r)=r sum_(n>=1) d_n[-B_Mn(r)] >=0,
    B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2).

The full angular sum obeys w(r)<=72C/r^5 at EVERY r>0. Fix a=1/128.
Let theta(t)=exp(-1/t) for t>0 and zero otherwise, and define

    zeta(r)=theta(2-r/a)/[theta(2-r/a)+theta(r/a-1)].

Thus zeta is smooth, one through a, zero from 2a=1/64. Set

    D_a=J* integral_0^infinity (1-zeta(r))w(r)[1-cos(rP)] dr J,
    H_s=H0-D_a, R_s=H-H_s=R+D_a.                    (1)

The same fixed canonical/York congruence and action factors are understood
in each term. The scalar and mixed rows are unchanged. On the form core,
the short comparison replaces the scalar bracket

    1-cos(rw)-(rw)^2/2+(rw)^4/24

by zeta(r)[1-cos(rw)]-(rw)^2/2+(rw)^4/24. Its polynomial terms and the
original fourth restoration remain EXACT. The finite contact is still
Ha2-Q4geo when the flow multiplier includes Q4geo. The new bounded operator
annihilates constants. Its diagonal row sum and nonlocal part must be
transferred together; there is no new choice of physical vacuum contact.

Spectral calculus gives 0<=1-cos(rP)<=2. Tonelli on nonnegative quadratic
forms, then polarization, includes all flow frequencies, trajectories and
angular masses. Since ||JT||SM^2=2pi||T||base^2,

    0<=D_a[T,T]<=72pi C a^-4 J0[T].                 (2)

The integral to infinity is 72C integral_a^infinity r^-5 dr; the two
additional factors are 2 from the spectral bound and 2pi from angles.
No small trapped measure or finite mode list is used. Using the RETAINED
physical conversion and actual positive diagonal form,

    C J0/E0<=j0=1000 C^2/(lambda B^4), l>=.99 Fl E0,
    delta_a=72pi a^-4 j0/(.99 Fl)<3.2e-18.           (3)

The actual parameter enclosure evaluates this upper-bound expression near
2.156202783669969564785e-18. This is not the actual norm of D_a; 3.2e-18
also covers the coarser inherited guards without depending on its last digits.

This is in the original action-form units. In particular

    ||R_s||_(V->V*)<=eta_s=.0022000000000000032,
    ||H_s^-1||<=1/(1-eta_s),
    ||K_s||<=q_s=eta_s/(1-eta_s), K_s=H_s^-1 R_s.   (4)

The slightly larger eta_s is not an improvement of the physical bound.
It is the full price for avoiding long-flow regularity in the comparison.
All actual thermal, angular, boundary and finite-domain returns remain in
H and R_s. In proving regularity of R_s, combine R+D_a FIRST: the artificial
long-flow term cancels exactly. Do not differentiate the two terms separately.

The module operators/short_travel_comparison.py defines these exact brackets
and the complete transfer bound. It is not a numerical Hessian application.

## 3. Local regularity after the exact cancellation

Here all regularity statements are for this FIXED analytic one-sided seed
and its actual paired domain. They assert finite constants, not numerical
ones or a theorem on an arbitrary nonlinear C5 neighborhood.

After (1), R_s consists of the actual local Dirac-minus-short-geodesic
kernel, its common contacts and the actual separated Dirac return. The
latter is smooth on separated physical-half charts, including boundary
charts with the original complementary projectors and regular bolt charts.
The original massive gap controls positive time and the full thermal trace;
the all-angular heat weight is summed before the common finite part. Its
finite-order endpoint derivatives follow by the same localized resolvent
and heat-parametrix argument as finite-graph-admissibility.md section 2,
applied to the current polynomial coefficients. No numerical constants
from that older finite metric path are transferred to the current seed.

The remaining geodesics have length <1/64. In the retained normal caps they
are local and have no conjugate point or long return. On a physical side
their coefficients are analytic. At a crossing use proper distance Y and
gbar=dY^2+A(Y)^2 dt^2, with A even, A'(0)=0, analytic for Y>=0.
Write c(Y)=A'(Y)/A(Y). A crossing with initial normal speed v>=0 satisfies

    Y''=c(Y)(1-Y'^2), Y(0)=0, Y'(0)=v.

On its positive leg set Y=v Z. Then

    Z''=[c(vZ)/v](1-v^2 Z'^2), Z(0)=0, Z'(0)=1.     (5)

Because c(0)=0 and the positive-side c is analytic, the divided coefficient
in (5) extends analytically to v=0. The negative leg is obtained by odd
reflection of Y, not analytic continuation of the metric through the seam.
On the short interval the coefficient bounds and Jacobi equation give
Z=t[1+O(t^2)], with the corresponding finite-jet bounds. Thus a crossing
between Y and reflected -Y' has

    v comparable to (Y+Y')/r,
    crossing fraction comparable to Y/(Y+Y').        (6)

In the variables r, v and that fraction, the short endpoint map is a small
perturbation of the flat one after dividing out these displayed factors.
Its inverse is nonsingular; the two ordinary angular charts handle v=1.
This is the local shooting/Jacobi argument, using r<1/64 and the retained
curvature bound <1.5. It is not a claim about the full geodesic flow.
At v=0 formula (5) supplies the limiting jets instead of division by an
uncontrolled transverse speed. Euler derivatives Y partial_Y and
Y' partial_Y' become bounded vector fields in the resolved variables.
Consequently every specified FINITE number of those derivatives preserves
the radial degree of the short direct/reflected kernel.

For example, the potentially nonsmooth cubic contribution on a crossing
segment has the elementary integral

    integral_0^1 |(1-t)Y-tY'|^3 dt=(Y^4+Y'^4)/[4(Y+Y')].

It is homogeneous of degree three and has bounded Euler derivatives on
the resolved positive quadrant. Replacing it by a smooth even polynomial
would be incorrect. Equation (5), rather than this example alone, supplies
the higher finite-jet argument for the actual polynomial metric.

The common flat fourth-log singularity and fourth contact cancel in R_s.
The retained Born comparison therefore has radial degree -3, instead of
-5, and local distribution order at most two. Differentiating its finite
parametrix on the just described scaled charts preserves the two-power
gain. The contacts must be differentiated through their actual compact
probe/heat expressions; differentiation of their unsigned envelopes is
not used. Fixed-seed finite-jet existence and the common regulator give
finite contact seminorms in this construction. Their values are UNKNOWN.

The scaling of these contacts is essential. On an annulus of size h, the
normalized metric differs from its flat boundary jet by O(h^2), with all
the finite scaled derivatives in use. In the combined degree-four heat
restoration this gives, for the coefficient a_l of a derivative of order
l<=2, bounds of the form |D^j a_l|<=C_jl h^(l-2)(1+|log h|)^N.
This follows by keeping the common scaled parametrix and probe together
before extraction. The older coarse 1/Y^2 first-contact envelope alone
would lose an outer power and does not establish this scaling assertion.
No numerical value for these sharper contact seminorms is claimed here.

For clarity, this is the estimate needed for each finite j and the finite
number of difference derivatives used to normalize the local symbol:

    |(Y d_Y+Y' d_Y')^j d_d^k K_s,rad|
        <=C_jk C |d|^(-3-k) (1+|log|d||)^N_jk,       (7)

with the analogous reflected Y+Y' expression, after including the local
finite-part extension and contact coefficients. The constants exist by
the finite local constructions above; (7) is not an explicit bound for
them. Normalize the order-(2+epsilon) symbol, then use Grubb's theorem at
output order zero on each scaled chart. For stationary radial charts N=2
frequency derivatives suffice; use N=3 in disk charts. Positive Holder
regularity comes from their smooth one-sided/resolved coefficients.
Commutation gives higher output orders. This respects |s|<tau; it does
not invoke that theorem at arbitrary output order on the C2,1 double.
Its direct-symbol use stays on annuli strictly inside the physical half;
the reflected term is handled by the separate weighted kernel/Schur argument.
No smooth extension of the physical coefficients through the seam is used.

Separated end terms have the original z^4 moment. When that moment is
integrable, keep it explicitly; its next Taylor remainder gains another
power. The moment-subtracted, differentiated shape row gains two outer
powers. The original undifferentiated moment stays in the matching row.
The short-travel cutoff does not alter the normal family or its exact roots.

## 4. Finite powers using the ACTUAL inverse; the bolt derivative cost

Put J_s=H^-1 R_s. The original weak inverse gives ||J_s||V<=eta_s<1,
and algebra on V gives

    K_s=J_s(I-J_s)^-1,
    K_s^p=J_s^p(I-J_s)^-p.                           (8)

This avoids assuming a strong inverse theorem for H_s. The strong bootstrap
uses the ACTUAL H inverse and its already identified physical normal family.
The positive fourth reserve and the fixed York map control critical-weight
H_b^2 from V. The latter is only a finite fixed-seed comparison here.

Take epsilon=1/4 to absorb the lower logarithmic symbol. From (7), R_s has
order at most 2+epsilon. H^-1 gains four derivatives in the scalar/shape
graph, with the shape logarithm retained. Thus each J_s gains at least
7/4 derivatives. The first application improves the critical outer weight
3/2 to 7/2, retaining the mass root at three. The second passes to sigma=
4.005, retaining the scalar quartic and EXACT slow root 4+e_*.

There is no unremoved affine z^4 log forcing in this error range: R_s has
zero canonical trace row and its separated original shape moment is z^4
times a constant. The combination (D-4) kills that moment only for the
regularity argument. The undifferentiated equation fixes its scalar quartic
return, whose coefficient is 20 beta Fstar. All phase coefficients are
extracted from the original H solution, never assigned zero. Subsequent
applications keep weight sigma and these finite phases; no further root
crossing is needed. The strip argument and exact implicit root are those
of current-global-boundary-moment-return.md and current-implicit-phase-
graph-domain.md, with the current seed's higher finite jets.

The successive regularity orders can therefore be taken as

    2, 15/4, 11/2, 29/4, 9, 43/4, 25/2, 57/4, 16.  (9)

The estimates are uniform over a dense energy core, including the actual
mass/log fields, because (7), the weighted inverse and phase extraction
are bounded linear maps with fixed finite constants. Completion extends
the composites to V. This is stronger than regularity of one chosen source.

Four powers reach H9, enough for the ORIGINAL H3 fast output (whose normal
row has order five with a logarithm), the finite phase coefficients and
bolt VALUE. Hence T_graph/phase/value K_s^4 is bounded at this fixed seed,
with a numerical constant still UNKNOWN.

That argument alone does not cover five AREAL metric derivatives at the
bolt. In regular disk radius rho, r-R is a smooth function with leading
term proportional to rho^2. The exact areal return contains one derivative
of the canonical lapse. Five areal output derivatives consequently require
six areal derivatives of the canonical fields, or twelve disk derivatives
and their point traces. H9 is not that embedding. This is a deficiency of
the four-power regularity certificate, not a proof that its actual map is
unbounded. Use EIGHT powers: H16 in two dimensions embeds into C13, and
radial evenness converts the required twelve disk jets to six derivatives
in rho^2. At the fixed positive bolt gap the subsequent areal coefficients
and the original weights are smooth with finite, possibly huge constants.
Thus the original physical C5 output T_X K_s^8 is bounded, with its actual
normalization retained. No numerical R_X or smallness follows.

Equation (8) costs at most (1-eta_s)^-p in this transfer. It supplies neither
the unknown constants in (7) nor the numerical actual-H graph inverse.
This establishes fixed-seed domain inclusion/finiteness, not quantitative
feasibility or uniform nonlinear-family admission.

## 5. Quantitative decision still required

The complete signed comparison responses have NOT been evaluated. For
each ORIGINAL output T_i, if M_i,p=||T_i K_s^p||_(E->Z_i) is numerically
bounded and U_j bounds the energy of the actual response to the full input
operator F_j, the retained exact identity gives

    error_ij <= M_i,p q_s^(N-p) U_j.                 (10)

All U_j except the affine E<=1.1e-10 remain unknown here. Both original
equations, all phase/interior columns and the common regulator are part
of F_j. Numerical application errors and the signed trial centers must
also be certified. No unit-size independent phase input is assumed.

The original reference mass lift and cutoff matching remain as in
current-coupled-global-return-attempt.md. Its local cancellations do not
remove the transition or make a global fixed-temperature zero mode.
The >.00906808164 mass-basis conversion remains in the ORIGINAL output;
both 1e-11 graph ceilings still require |m|<2.20553814952e-9. The actual m
is unknown. The sufficient kappa_m<=5e-5 target is also uncomputed.
Its conditional budgets are unchanged: 5.34235e-12 affine and .0012069
graph feedback remain, with an additional .0195239 energy coefficient
requiring its own comparison. No other entry is shown to fit.

The fast screen remains .135973<.1599164 only for its two-fast-input
admitted-family scope. Still unproved are the physical norm comparisons,
E_center+1e-8 alpha<=1e-11, R_X alpha*.028913+Z0_other<=.25,
R_X alpha^2*1.0706e10+Z2_other<=5e7, and

    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

RESEARCH DECISION: retain the construction and use the short-travel
comparison for the strong return. Its full extra form cost is negligible
on the existing scale and removes an unnecessary global-flow regularity
obligation. Use four powers for graph/phase/value and eight for the original
physical C5 return. This changes the METHOD, not the admission norm.

ONE NEXT CALCULATION: quantify the original-norm composite constants in
(10), using the differentiated, moment-subtracted R_s kernel and the actual
coupled H inverse, then apply the same bound to all columns. The numerical
strong constants and signed solve remain PROPOSED, NOT EXECUTED. The new
form transfer and fixed-seed domain argument are executed results; neither
is an executed controlled-solution feasibility test.

Remaining to-do list: numerical complete strong return and coupled signed
response; physical-norm comparisons and nonlinear admission/contraction;
exact compatibility; physical matching/causal realization; semiclassical
validity/noise; all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
