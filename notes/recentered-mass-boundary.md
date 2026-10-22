# Recentered stationary mass boundary cancellation

2026-10-22, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress, satellite a30f27a and planning 3ec526e.
The full stable semiclassical wormhole objective remains incomplete.

The complete stationary outer mass momentum is zero on the selected
fast phase of the SAME c_R=10^22 theory. Its variation within that
phase is also zero. No extra mass boundary equation is required there.
This identifies the stationary metric boundary part of the preceding
bulk reference realization; the full current/constraint problem,
nonlinear admission and Lorentzian stability remain open.

The mass cancellation and restricted C1 determinant theorem already
exist for the zero-additional-c action in
[stationary-fast-phase-variation.md](stationary-fast-phase-variation.md).
They were not absent from the retained research. The work here checks
their hypotheses on the new reference, includes the added R^2 bulk and
boundary action exactly, and fixes its same-source stationary charge.
No old scalar phase projection or numerical nonlinear constant transfers.

## 1. Phase, physical data and the actual paired domain

Keep G, q, g_R, mu0=1/b_class, c_R=10^22, fixed canonical sphere B,
paired theta=0 domain, and the prescribed positive period beta. The
physical lambda satisfies lambda0-lambda=11G lambda^2/(360pi).
Use the same fast stationary phase, in areal coordinates,

    ds^2=exp(2alpha)F dU_E^2+dr^2/F+r^2 dOmega^2,
    F=lambda r^2+1-2m/r+R_F,
    alpha=alpha_infinity+R_alpha,
    partial_r^j R_F=O(r^(-2-j)L^p),
    partial_r^j R_alpha=O(r^(-4-j)L^p), j<=5,
    L=1+log(r/B), p fixed and finite.                       (1)

Work locally in this weighted topology, with a regular nonextremal bolt
and smooth disk charts. Allow the clock source alpha_infinity in an open
interval, but set its tangent to zero for the physical fixed-clock solve.
The magnetic bundle and leading scalar curvature source are fixed.
The new RN reference satisfies (1) exactly with R_F=Z/r^2 and alpha=0.
The forced r^-4 log r conformal tail is allowed. Growing scalar roots
are excluded by the physical leading data and energy; the decaying
fast scalar belongs to the remainder. No arbitrary finite wall is used.

The new massive-domain proof already gives M_n>18 sqrt(n) and the
actual complementary-adjoint resolvent bound 1/M_n. The remaining
LLL kernel hypothesis can be checked directly, rather than importing a
gap from the old seed. In the canonical disk the orbit metric is
exp(2V)(d xi^2+xi^2 dvarphi^2), with V finite and positive conformal
factor at the bolt and outer circle. The spin structure is the smooth
disk one. The orbit has no gauge field in this purely magnetic sector.
Two-dimensional conformal covariance maps a zero mode to a flat disk
zero mode, preserving the first-order projector and its adjoint.

In a Cartesian Pauli frame a flat zero spinor has components f(z)
holomorphic and g(bar z) antiholomorphic, regular at the origin. On
the unit circle the actual tau3 gamma_normal projector gives, for
each flavor, f=+/-exp(-i varphi)g. The first trace has Fourier powers
n>=0; the second has powers -(n+1)<=-1. They cannot agree unless
both vanish. The complementary adjoint has the other sign and the
same conclusion. This proves absence of LLL zero modes, not a
numerical lower bound for the gap. The sphere LLL degeneracy does not
alter the argument. The positive paired singular-value prescription
and its phase remain those of the retained global intertwiner.

The conformal domain is compact and elliptic, so at fixed positive
temperature its kernel-free square has a positive discrete gap. On a
sufficiently small relative C1 graph neighborhood this remains true.
Neither the size of that neighborhood nor a temperature-uniform gap
is assigned a numerical value. The massive direct gap and this LLL
argument check the large-time hypothesis of the old C1 action theorem.

## 2. Exact local action and boundary decomposition

Use the retained Lorentzian action sign when displaying boundary work;
the overall Euclidean sign reversal does not change vanishing. Define
Z_R=R+12lambda and F_lambda=1-384piG c_R lambda. Then identically

    (R+6lambda0)/(16piG)+c_R R^2
      =F_lambda(R+6lambda)/(16piG)
         +c_R Z_R^2+3(lambda0-lambda)/(8piG).           (2)

The local scalar-tensor boundary action, with auxiliary field R on
its algebraic equation, likewise splits as

    [1/(8piG)+4c_R R]K
       =F_lambda K/(8piG)+4c_R Z_R K.                (3)

This is the actual F(R) boundary completion; it retains the curvature
datum, as in [Dyer--Hinterbichler (1.4), section 5.3](https://arxiv.org/pdf/0809.4033v2).
It is not metric-only Dirichlet at a finite cut. The first term in
(2)-(3) has the standard Einstein AdS counterterms multiplied by
the fixed F_lambda. The residual square has no divergent action or
surviving boundary term on (1). No new finite coefficient is selected.

In fact Z_R=O(r^-4 L^p), with its differentiated bounds from (1).
Its radial bulk square density is O(r^-6 L^(2p)), integrable. Its
boundary value 4c_R sqrt(h) Z_R K is O(r^-1 L^p), and its first
variation has the same vanishing power for all allowed tangents.
The square bulk endpoint coefficients have sizes

    Q_F=O(r^-2 L^p), P_F=O(r^-3 L^p),
    Q_alpha=O(L^p), P_alpha=O(r^-1 L^p),              (4)

where Theta=P_F delta F+Q_F delta F'+P_alpha delta alpha+
Q_alpha delta alpha'. Products with mass/fast tangents vanish. Even
an open constant clock variation leaves only O(r^-1 L^p). The finite
c_R multiplies these limits; it is not sent to zero. Locally uniform
families obey the same bounds for the variation of the boundary
one-form. Thus both its endpoint value and its linearization vanish.

## 3. All counterterm finite pieces and the exact mass sum

Keep the existing minimal intrinsic loop subtraction, per exterior,

    B_ct,loop=4pi exp(alpha)sqrt(F)(c0 r^2+2c1),
    c0=-221lambda^(3/2)/(15120pi^2),
    c1=29sqrt(lambda)/(20160pi^2).                    (5)

These are the retained divergent counterterms, not adjustable finite
parameters. To see why they do not change with c_R, the residual
vacuum volume in (2) contributes
11lambda^2 r_cut^3/(720pi)=231lambda^2 r_cut^3/(15120pi)
per unit U. The complete Weyl bulk-plus-boundary action contributes
653lambda^2 r_cut^3/(15120pi), as before. Their sum is
221lambda^2 r_cut^3/(3780pi). The linear-cut coefficient is
67lambda/(3780pi). Equation (5) cancels both powers. The added square
is integrable; the effective Einstein part has its own F_lambda-scaled
subtraction. The compact canonical determinant adds no physical-cut
power divergence. The exact quantum-vacuum recentering is essential
to this allocation; the scale mu0 is unchanged.

For clarity, the mass endpoint calculation must vary the full radial
action BEFORE restricting to RN. For the complete Weyl bulk density

    L_W=exp(alpha)/(4pi r^2){
      sigma(Rbar^2/60-Rbar/18+1/15-q^2/3)
      +Rbar T/30-T^2/10+11TF/90-11F^2/180-F/18},
    sigma=log(r/B), T=r(F'+F alpha'),
    Rbar=-r^2F''-3r^2F'alpha'-2r^2F(alpha''+alpha'^2)
                                                   +2T-2F,

put J=sigma(Rbar/30-1/18)+T/30 and
K_Rbar=Rbar/30-T/5+11F/90. At alpha=0 its derived endpoint is

    Q_F=-J/(4pi), P_F=[(2J+K_Rbar)/r+J']/(4pi).       (6)

For delta F=-2delta m/r, delta alpha=0, (6), the full literal-boundary
representative B_W=-B^2 N Gcal/(360pi), and (5) give

| Contribution to fixed-clock mass momentum | Coefficient in lambda delta m/(3780pi) |
| --- | ---: |
| Weyl bulk integration endpoint | -105 |
| Complete Weyl boundary variation, with conversion fluxes | -116 |
| Covariant loop counterterm, INCLUDING its finite mass piece | +221 |
| Canonical determinant endpoint | 0 |
| Effective Einstein bulk+GHY+counterterms | 0 |
| Residual curvature square and its boundary term | 0 |
| Total | 0 |

The boundary polynomial and all its radial metric variations are those
of [the full Weyl calculation](stationary-exterior-boundary-weyl.md),
not an action value used as a substitute for a first variation.
The coefficient 221 includes vacuum centering. The square decomposition
checks the new theory rather than merely carrying over a zero-c number.
No fixed-cut normal derivative is set to zero by hand.

The allowed fast remainders give O(r^-1 L^N) endpoint work in the
Weyl and intrinsic pieces, by the explicit bounds of the retained fast
phase theorem. Their weights are unchanged. The scalar and mixed
limits from the preceding S169 notes agree with this result.

## 4. Canonical determinant: retained distribution and limit hypotheses

The canonical boundary contribution in the table is not obtained by
dropping the quantum source. The actual orbit trace is a distribution
whose cubic contribution is

    A0 FP_ell(1/y)
      +A0[log(mu0 ell)+gamma_E+31/30]delta(y),
    A0=c3 B^2/(5pi), c3=-m sqrt(lambda)/B^3,          (7)

plus a locally integrable remainder. Equation (7) includes the full
massive tower, both partner squares and the boundary anomaly; the LLL
has zero mass weight and its local anomaly is included. The retained
Mellin proof fixes the delta coefficient. At the canonical geodesic
boundary all normal-smearing derivative rows vanish in the balanced
mixed domain; tr chi=tr(chi E)=0. In particular there is no delta'''
term that could couple to a mass cubic. This uses the complete mixed
coefficient, [Vassilevich (5.33)](https://arxiv.org/html/hep-th/0306138v3#S5.E33),
with the actual complementary first-order domains, not a guessed wall.

For a general allowed fixed-clock metric tangent, its canonical half
relative components h,j have order y^3. In a fixed proper chart,

    xi=a int_0^y (h-j)/a, omega=j-(a'/a)xi,
    xi=O(y^4), omega=O(y^3).                           (8)

The delta in (7) has zero pairing with omega; FP(1/y)omega is an
ordinary integrable function. The remaining distribution is L1 and
has no additional surface support. Canonical pressure is finite and
its xi boundary work vanishes. The nonzero bulk determinant response
is retained in the equations. Vectorlike covariance and the paired
phase convention justify the boundary-preserving diffeomorphism in (8).

The differentiability/finite-cut argument already retained in
stationary-fast-phase-variation.md applies here after the checks in
section 1. It uses the FIXED compactified areal coordinate, whose
radial and temporal coefficients have opposite mass cubics. Both
coefficients and their weighted five jets are present. The exact
cubic/a4 profiles and every lower finite-cut mixed heat profile are
removed with their common meromorphic finite parts before taking the
cut to zero. The remaining parameter-differentiated orbit heat trace
has L1 bound

    C t^(3/2)(1+|log t|)^(N+1).

Multiplying by the angular ACTION heat trace O(t^-1) and dt/t leaves
the integrable majorant C t^-1/2(1+|log t|)^(N+1)dt. The cubic
source remainder uses the corresponding integrated boundary profile;
it does not acquire another delta distribution. Low finite-cut
linear jets are O(d), and their full meromorphic profiles are handled
before the limit. No unsubtracted ultraviolet sum is dominated just
because a coefficient tends to zero. Local constants depend on the
new reference and neighborhood and are not assigned old numerical
values. The actual gap proved above handles the large-time part.

This proves C1 action transfer on (1), including the limiting first
variation. The boundary one-form is identically zero on every
fixed-source fiber, so its linearization within that fiber is zero.
We do NOT infer a quantitative C2 source map, a prescribed-power
nonlinear closure, or an admission ball from this C1 result. The
reference bulk Hessian and graph closure remain the separate actual
two-vertex result of the previous note; the present calculation removes
a possible extra mass boundary row on their smooth mass/fast class.

In particular W_(7/2) membership alone does not imply the pointwise
O(z^4 log^p z) phase (1). A reference inverse solution for the actual
source still needs that endpoint regularity/expansion and the complete
admission norm. The Mellin pole at weight 4 is not crossed by fiat.
Boundary cancellation on (1) does not silently upgrade the preceding
weighted inverse into a bijection of this nonlinear physical phase.

## 5. Stationary charge and physical matching consequence

The same calculation with an open clock source gives

    Theta_infinity=-exp(alpha_infinity)
           [Z_c m/G+E_sources]delta alpha_infinity,
    Z_c=F_lambda+(G lambda)/(5pi)
       [log(mu0/sqrt(lambda))+gamma_E+151/180].       (9)

It follows from the retained canonical/Weyl clock distributions and
the F_lambda-scaled Einstein term. The residual square gives zero
endpoint charge by (4). A finite intrinsic functional of the fixed
conformal sources can change E_sources, not the mass slope. The
normalizable mass cancellation holds for nearby clock values too,
which excludes an invisible clock-times-mass finite adjustment.

For the actual benchmark, F_lambda>1-1.179e-8 and

    |Z_c-F_lambda|<1e-33, Z_c>0.9999999882.            (10)

The leading charge change from the added theory coefficient is about
1.1781e-8 relative to m/G, not the old quantum-only 1e-33 allowance.
Use same-source differences about the NEW full reference:
Delta E_stat=Z_c Delta(m/G). Do not discard an absolute correction
proportional to the large extremal mass when comparing an extremely
small preparation gap. This is a stationary boundary conjugate; it
is not a constructed exact state, a full dynamical Hamiltonian or a
Lorentzian stability result. No measured physical parameter is changed.

## 6. Decision and verification contract

The previously decisive mass boundary target is met EXACTLY on the
declared stationary fast phase. The preceding unit-mass trace
discontinuity remains true; its dangerous coefficient is zero here.
Continue with the joint current/Ward and radial constraint realization,
the actual inverse solution's physical fast-tail regularity, then the
finite geometry/state, full residual and admission-norm inverse. The
reference energy bound 1 and relative mixed bound .336
retain their scopes. No seed-to-finite defect <=1/4 is yet proved.

This phase is a sufficient specified stationary variational class, not
a theorem for arbitrary finite-wall, nonspherical or driven boundary
data. The work does not select a different ensemble, remove the
curvature scalar, tune temperature/flux, or complete S169.2. The full
stable-wormhole goal remains unchanged. The preserved-ball requirements
rho=1e-8, Y0<=5e-9, Z0<=.25, Z2<=5e7 and EVERY positive admission
margin still require actual estimates. Compatibility, semiclassical
omission/noise control and all-sector causal stability remain open.

Before running, freeze this note, checker and source review. Independently
differentiate the unexpanded radial bulk action before RN substitution,
include the complete boundary polynomial and intrinsic counterterms,
and verify all limiting mass coefficients exactly. Check the local
bulk/boundary square identities and their fast weights, vacuum-divergence
allocation, common cubic distribution constants and LLL Fourier-domain
argument. At 70/100 digits compare direct finite-cut endpoint variations
at increasing radii with their limits, and evaluate the actual charge
factor without replacing the tiny quantum-vacuum shift by zero. These
are supporting diagnostics; they do not numerically certify the heat
limit, current equations, nonlinear solution or stability. Retain every
run, failure, source snapshot and raw log. No Lean is selected.

Remaining to-do list: complete stationary current/Ward and radial
constraint realization and physical fast-tail inverse regularity;
finite state, residual, admission and nonlinear
solution; exact compatibility, same-source physical matching, validity,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and full-result audit.
