# Anomaly-complete auxiliary reduction and scalar/metric coupling

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The auxiliary Einstein-frame method is viable only if its curvature
constraint is retained. With the original split, replacing that constraint
by a scalar wave equation leaves a differential constraint-defect equation
with no coercivity argument. An EXACT repartition of the existing local
R^2 term removes that defect derivative without changing the total action.
The curvature constraint is then algebraically recoverable in an explicit
neighborhood. No finite theory parameter or subtraction prescription changes.

In the resulting variables, the large local charged scalar/metric mixing
vanishes at the matched constant scalar. On the actual product-centered
seed, the complete scalar-row coupling obeys

    ||B_E h_E||_w/lambda <1e-27 N2(h_E).                (1)

The independent auxiliary scalar row has a positive massive inverse.
For its actual forcing and metric input it gives

    ||delta u||_w <= [1.51e-28+1e-27 N2(h_E)]/(3.2e8).
                                                               (2)

For N2(h_E)<=1, its cold VALUE norm is <2e-17. This is a conditional
coupled-row estimate, not a bound on the unknown Einstein-frame metric
response, its higher derivatives, or the full nonlinear solution.

## 1. Exact repartition and the constraint that cannot be omitted

Use the retained Lorentzian stress convention and full massless Dirac
functional Gamma_D. Put

    epsilon=G/(60pi), H=7 C4^2/12+11 J^2/6-20I,
    8piG tau_D=epsilon[H-11R^2/72+Box R],
    delta_c=1/(5760pi^2), c_*=c_R-delta_c,
    Gamma_*=Gamma_D+delta_c integral sqrt(-g) R^2,
    D=96piG c_*=96piG c_R-epsilon.                    (3)

The local R^2 stress has trace -12 delta_c Box R. Thus

    c_R R^2+Gamma_D = c_* R^2+Gamma_*,
    8piG tau_*=epsilon[H-11R^2/72].                   (4)

These are exact action and bulk-trace identities. The common R^2
boundary completion is repartitioned with the same coefficient. There
is no independently renormalized mode, discarded anomaly, or new theory.
c_*>0, delta_c/c_R<1.8e-27. The Dirac operator, state and paired domain
are unchanged by adding this explicit local functional to Gamma_*.

Introduce phi for c_* R^2, with Phi=1+32piG c_* phi and C=R-phi.
Let E_ab be the complete auxiliary metric residual:

    E_ab=Phi G_ab+Lambda g_ab+8piG c_* phi^2 g_ab
         +32piG c_*(g_ab Box-Hess_ab)phi
         -8piG(T_EM,ab+T_*,ab), Lambda=-3lambda0.

Define the scalar row

    S=phi+12lambda0-D Box phi+epsilon[H-11phi^2/72].

Tracing the actual metric equation before imposing C=0 gives exactly

    -tr E-S=C[Phi-11epsilon phi/36-11epsilon C/72].    (5)

For |phi+12lambda|<=lambda and |C|<=lambda, the bracket exceeds
1-2e-8; its reciprocal is below 1.00000003. Therefore tr E=S=0
forces C=0 in this neighborhood. Conversely C=0 reduces S to the
original full trace equation. The numerical neighborhood is an explicit
constraint-recovery condition, not a certified nonlinear solution ball.

If one uses c_R and Gamma_D without (3), the right side of (5) also
contains epsilon Box C. That term cannot be silently dropped. The
equation epsilon Delta C+approximately C=0 is not the positive-mass
maximum-principle equation. No additional solution or instability is
claimed from that diagnostic; it shows why the naive row substitution
does not establish equivalence.

The row equivalence also keeps Newton contacts. With fixed reference
density and C=0 at the seed, differentiation gives

    delta[sqrt(g)(-tr E)]-delta[sqrt(g)S]
       =sqrt(g)[Phi-11epsilon phi/36] delta C.          (6)

No derivative of delta C occurs, and the varying trace/density on the
left must be retained. One can use the three stationary rows -tr E,
E_x+E_phi and S for the two disk metric fields and phi. Their vanishing
first forces C=0 by (5). The auxiliary Ward identity then becomes the
original conserved metric identity, and the regular-bolt recovery in
[the stationary-constraint note](recentered-stationary-constraints.md)
recovers every metric component. The unchanged magnetic current still
vanishes. This is a local row-equivalence argument, not an inverse bound
for the whole three-row system.

## 2. Exact vacuum-matched Einstein-frame action

Define

    F_*=1-384piG c_* lambda, Phi=F_* exp(2u),
    g_E=exp(2u)g, phi(u)=[F_* exp(2u)-1]/(32piG c_*),
    p_v=-3(lambda0-lambda)/(8piG),
    Gamma_*ren=Gamma_*-p_v Vol_g.

The frame is normalized at infinity: u=0 there, so the prescribed outer
clock and period do not change. F_*-F_lambda=Glambda/(15pi)<2.1e-35,
where F_lambda is the original c_R reference coefficient. This very
small difference is kept algebraically, not used to retune the geometry.

With the transformed local boundary action, the exact local bulk part is

    F_*/(16piG) integral sqrt(-g_E)
      [R_E+6lambda-6(grad_E u)^2
                      -(1-exp(-2u))^2/(64piG c_*)]
       -integral sqrt(-g_E) I_E/(4g_R^2)
       +Gamma_*ren[exp(-2u)g_E,A].                    (7)

In particular the potential has a stationary point at u=0. A direct
check starts with K=1/(1024pi^2 G^2 c_*) and
F_0=1-384piG c_* lambda0. The uncentered potential is
K[F_*^2-2F_* exp(-2u)+F_0 exp(-4u)]. Subtracting
p_v exp(-4u) replaces F_0 by F_*, because p_v/K=F_0-F_*.
This gives (7), including its constant cosmological term.

For any Einstein-frame metric, the local scalar/metric mixed Hessian
at u=0 is exactly zero: grad u=0 and the potential derivative vanishes.
On a nonconstant background u_b its compact mixed form is

    F_*/(16piG) integral sqrt(g_E)
     {-12[(tr h_E/2)g_E-h_E](du_b,dv)
                              -(tr h_E/2) W'(u_b)v},
    W(u)=(1-exp(-2u))^2/(64piG c_*).                  (8)

This displays the residual/potential contact instead of treating the
seed as on shell. The full quantum mixed term still exists and is
calculated from the retained trace, not discarded by conformal covariance.
The conformal map in [De Felice--Tsujikawa, section 2.3](https://arxiv.org/html/1002.4928v2#S2.SS3)
is the standard comparator; (3)-(8) supply the present anomaly, vacuum
and constraint bookkeeping.

## 3. Complete scalar row in frame variables

The conformal identities give

    Box_g phi=F_* exp(4u)/(16piG c_*) Box_E u,
    J_g=J_E+2(Hess_E u)_TF+2(du tensor du)_TF=:Z
             as covariant two-tensors,
    exp(-4u) H_g=7 C_E^2/12+11|Z|_E^2/6-20I_E.

Consequently the rescaled row S_E=exp(-4u)S/(6F_*) is exactly

    S_E=-Delta_E u+P(u)
       +epsilon/(6F_*)[7 C_E^2/12+11|Z|_E^2/6-20I_E],
    P(u)=exp(-4u)[phi(u)+12lambda0-11epsilon phi(u)^2/72]/(6F_*),
    P(0)=0, P'(0)=C_v/D,
    C_v=1+11epsilon lambda/3.                         (9)

This is second order, although nonlinear in Hess u. The full Dirac
trace is local; the remaining metric equation still has the complete
nonlocal determinant. Formula (9) is used with the metric rows and
constraint recovery (5), not as a replacement for the whole system.

The metric remains the exact product-centered seed already constructed.
Set phi_b=R[g_b] and

    u_b=(1/2)log[1+(32piG c_*/F_*)(R[g_b]+12lambda)].

The inherited polynomial geometry, |d|<8.6e-19, x>8e9 and
Glambda<2000^-10 give the following explicit guards:

    |u_b|<5.4e-37,
    |grad_E u_b|/sqrt(lambda)<1e-35,
    |Hess_E u_b|/lambda<1e-33,
    |Riem_E|<100lambda, |Z_b|<20lambda,
    I_E<q^2/b^4, Dlambda<3e-9.                        (10)

For the first three, use |r(t)|<=10, |t r'|<=56,
|t^2 r''|<=372 in R+12lambda=d r(t)/b^2. The logarithm's
denominator differs from one by less than 1.1e-36. On t>=1+h,
f<=7x(t-1)^2, 0<f'_t<=14x(t-1+h), so f/r^2<=7lambda and
f'_r/r<=14lambda. The normalized first two u derivatives then give
(10), including the conformal correction to the Hessian.

For curvature one may use the actual orthonormal sectional bounds
15lambda, 6.1lambda and 8.1lambda for the radial, mixed and sphere
sections. Their full Riemann norm is below 50lambda before the tiny
frame correction. The traceless Ricci and Hessian corrections fit well
inside |Z_b|<20lambda. These bounds concern this explicit seed.

## 4. Scalar inverse and the quantitative mixed operator

At fixed g_E, the complete linearized scalar row on the seed is

    A_E v=-[g_E^(ab)-11epsilon Z_b^(ab)/(9F_*)]Hess_ab v
        +22epsilon Z_b(du_b,dv)/(9F_*)+P'(u_b)v.       (11)

Its principal perturbation is <1e-32, its drift divided by sqrt(lambda)
is <1e-66, and P'(u_b)>3.3e8 lambda. For the SAME physical weight
w=[(t-1+h)/t]^(2/5), the preceding radial barrier proof and (10) give

    Delta_E w<3.38lambda w,
    |grad_E w|<2sqrt(lambda)w, |Hess_E w|<20lambda w,
    A_E w>3.2e8 lambda w.                             (12)

For the Hessian bound, w'/w=p(1-h)/[t(t-1+h)], p=2/5, and
w''<0. The f and f' guards above bound the three independent radial
Hessian components by 8.4lambda w, 2.8lambda w and 2.8lambda w.
The frame correction is negligible within the stated allowance.
The positive maximum principle with a regular bolt and zero auxiliary
source at infinity therefore proves

    ||A_E^-1 F||_w <=||F||_w/(3.2e8 lambda).           (13)

Existence is the regular radial positive-elliptic Green realization;
the massive fast branch and decaying forced tail are selected at the
end. This is the scalar block of the three-row problem. It does not
give the metric inverse or the full transformed boundary domain.

To quantify the remaining coupling, define in the Einstein-frame seed

    N2(h_E)=sup w^-1 max(|h_E|,
                lambda^-1/2|nabla h_E|,lambda^-1|nabla^2 h_E|).

Here h_E is a regular fixed-period metric direction in the pulled-back
physical fast phase, with the leading outer metric, clock and curvature
data fixed. N2 alone does not impose those endpoint conditions. The
source and mixed forcing then decay at the end, as required by (13).
All tensor norms are orthonormal Hilbert--Schmidt norms. This is an
explicit INPUT norm, not one inferred from the old energy estimate.
The metric derivative B_E of (9) satisfies (1). Here is its ledger,
including an allowance for the fixed-density contact:

| Contribution to ||B_E h_E||_w/lambda | Upper coefficient times N2 |
| --- | --- |
| Metric Laplacian acting on u_b | 2e-33 |
| Curvature anomaly terms | 2e-30 |
| Fixed magnetic field variation | [40/(360pi F_*)] Gq^2/(lambda b^4) <1.5e-28 |
| Fixed-density reference-residual contact | 1.51e-28 |
| Sum | <1e-27 |

For the curvature row, varying the connection and curvature in an
orthonormal tensor norm gives |delta C_E|<=5000lambda N2 w and
|delta Z_b|<=6000lambda N2 w from (10). The bounds include metric
contractions and connection variations. Hence the variations of the
squares are at most 1.04e6lambda^2 N2 w and 241600lambda^2 N2 w.
Multiplying by 7/12, 11/6 and epsilon/(6F_*) gives the displayed
allowance. The magnetic term uses |delta I_E|<=2I_E|h_E|; despite
small curvature it dominates the curvature anomaly contribution.
The density term uses |tr h_E|/2<=|h_E| and w<=1.

On the seed phi_b=R, so S equals the original complete trace residual.
The preceding full-source theorem gives

    ||S_E,b||_w/lambda<1.51e-28.                       (14)

Equations (1), (13) and (14) prove (2), without solving for h_E.
The same cold weight conversion as before has

    sup_cold w/a_s <=2^(2/5)/[64 h_min^(3/5)],
    h_min=1/(8*2000^10).

Thus N2(h_E)<=1 implies sup_cold |delta u|/a_s<2e-17;
with h_E=0 the bound is <2.2e-18. These are VALUE bounds for the
auxiliary conformal contribution. They are not its full C5 estimate
or a claim that N2 of the actual metric response is known.

## 5. Returning to physical geometry and boundary data

The original metric variation is exactly

    delta g=exp(-2u_b)[h_E-2(delta u)g_E].              (15)

The auxiliary contribution to each logarithmic radius is -delta u;
its tensor Hilbert--Schmidt contribution is 4|delta u|. The remaining
h_E term must itself meet the actual cold, outer, raw-jet and mass
admission conditions. Equation (2) does not prove the older
||B xi||_w/lambda<=2.8e-27 target in a different decomposition. It
supplies a sharper joint organization with its return map displayed.

Do not compare u_b to zero as though the frame change were an admitted
cold tangent: its product-center value is finite relative to the much
smaller h. Use the transformed exact product and seed as the center,
and estimate delta u about u_b. The smooth positive conformal factor
preserves the regular cone and fixed period. u_b tends to zero, and
the curvature datum remains fixed at infinity.

The boundary repartition in (3)-(4) includes 4 delta_c R K and its
counterterm allocation. After auxiliary elimination it restores exactly
the original [1/(8piG)+4c_R R]K completion. The local c_* scalar-tensor
GHY transformation cancels the normal derivative from the bulk frame
change; Gamma_* retains the complementary local piece. It is not an
independent boundary condition or an omitted surface term. The full
functional is always evaluated at exp(-2u)g_E with the pulled-back
actual paired domain and state. It is not replaced by Gamma_D[g_E].
Full nonlinear domain and boundary estimates still need proof in the
chosen solver spaces; the algebraic equivalence does not supply them.

## 6. Research decision

Use this anomaly-complete auxiliary organization for the next METRIC
response calculation. The scalar/metric coupling is now quantified
from an explicit weighted two-derivative input, and the local charged
mixing has been eliminated without losing the curvature constraint.
The dominant unresolved work is the actual Einstein-frame metric
response to the full traceless quantum source, its cold derivatives
and the outer/mass return. Do not substitute the local trace for that
source or claim that the small auxiliary block solves the metric rows.

At rho=1e-8 the original complete-space targets remain Y0<=5e-9,
Z0<=.25 and Z2<=5e7, with every admission margin and exact compatibility.
The block factorization, source C2 bounds, raw C5 admission and all
finite-state/domain transfers remain required. Semiclassical omissions,
fluctuations, all-sector retarded stability and the separate opening,
transport and creation obligations are unchanged. No positivity result
here is a Lorentzian stability or wormhole existence theorem.

Remaining to-do list: actual Einstein-frame metric response and full
physical admission/return, nonlinear source and domain estimates,
global mass/kernel moments; nonlinear solution and exact matching;
semiclassical validity, all-sector stability, charge/return, opening,
transport, preparation, separate creation and full-result audit.
