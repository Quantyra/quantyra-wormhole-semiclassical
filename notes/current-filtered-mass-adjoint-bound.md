# Filtered mass functional: coupled adjoint bound and residual-estimator test

2026-10-29 HST, E003/S169. Ordinary argument and execution contract frozen
before interval execution. Previous goal turn: NO PROGRESS, a plan review.
Preserve e7abfb5, the closed audit, P32, all accepted results and original
physical norms. The full stable semiclassical wormhole objective is open.

This attempts the complete signed mass term selected in e7abfb5. It does
not claim that an adjoint norm alone evaluates that term. The execution
will either bound the required filtered mass functional or leave that
implication unsupported, then test a specified additive estimator for the
two residual pairings. No signed coefficient is supplied by a ceiling.

## 1. Exact spaces, coupled source and normalization

Use the FULL normalized Hessian h=[-a,b;b*,ell]/N, N=|W| beta, and
E_total^2=Es^2+E^2=(a+ell)/N. Below a,b,ell in pairings denote these
normalized forms. T=ell^-1 b* a^-1 b, F=T(I+T)^-1, C_op=(I+T)^-1,
M=C_op F^2. The same original paired mass/log matter realization,
complementary adjoint/squares, common regulator and contacts are retained.
There is no mode truncation, new wall, clock or temperature variation.
All angular and thermal sums, the ultraviolet logarithm, finite-domain
returns and the regular disk bolt remain in the FULL forms.

For canonical (nu,t), the frozen York inverse is

    L_y y=-y''+[6/z^2+2(Q''-4Q'/z)/(3Q)]y
         =2t'/(3z^2)-8t/(3z^3),
    xi=-t-zy-z^2 y'/2, s=nu+t/2+z^2 y'/4-zy/2.

The regular bolt condition is y(L)=0. On physical mass/log inputs the
outer condition is the existing admissible branch. On the AUXILIARY
root-zero adjoint below it instead has the forced asymptotic
y=-2/(3z)+O(z log^k z), with no z^-2 homogeneous branch. This is the
same affine adjoint extension as the accepted mass Green identity, not a
new physical metric direction. The physical forward York conditions
remain J(0)=0 and (Q^2 J'/z^4)(L)=0. The already proved whole-interval
potential inequality is L_y's potential >=5/z^2.

Let Psi_ref=(n_ref,t_ref) be the EXACT correlated reference adjoint in
current-reference-clock-adjoint-trial.md. Its two reference rows vanish,
including the canceled bolt logarithm and regularity. Its root-zero
canonical trace is (-1,1), ||t_ref||infinity<=B0=1.000001, and its
complete FULL residual r=h Psi_ref satisfies ||r||E_total*<2.1e-5.
The latter is the complete commonly restored result of
current-complete-weak-mass-residual.md, not a new local approximation.

The missing implication in an attempted reuse is numerical control of
the adjoint's scalar York coordinate s_ref. The earlier single-source
linear mass bound supplies no bound on this different filtered map.

## 2. Bound that scalar coordinate without an invalid root-zero energy norm

First prove ||xi_ref||L2(dz)<=8 B0. Write
A_y=||z^2 y'||2, B_y=||zy||2. The actual endpoint expansions give zero
boundary terms and finite A_y,B_y. Weighted Hardy, obtained by expanding
||z^2 y'+3zy/2||2^2, gives A_y>=3B_y/2. Multiply L_y y by z^4 y and
integrate the t' term BEFORE estimating:

    A_y^2+int(z^4 potential-6z^2)y^2
       =-(2/3)int t(z^2 y'+6zy),
    left >= A_y^2-B_y^2 >=5 A_y^2/9.

It follows that A_y<=6||t||2, B_y<=4||t||2, and
||xi||2<=||t||2+B_y+A_y/2<=8||t||2. L<1 gives the claimed bound.
Truncate first and use the displayed endpoint rates; this does not
require a smooth reflected cubic or an open-clock nonlinear family.

Use BOTH reference equations and the exact frozen York congruence.
Their entire scalar row is

    a_ref s_ref = z^4 B_ref(z) xi_ref,
    B_ref(z)=2 F_l[6 beta+18 eps z^2-15 eps z+2 eps+2/C_geom]/beta,
    a_ref s=F_l[-6 z^2 L_s/beta+V_ref s],
    L_s=Q s''+(Q'-2Q/z)s',
    V_ref=24+12 eps z^5/beta-24 eps z^6/beta.             (1)

Here a_ref is written as its strong operator in Hmu=L2(z^-4 dz).
The checker derives (1) from the COMPLETE reference density in
operators/signed_local_core.py; in particular the reference mixed row
has no remaining xi derivative, radial-return or J coefficient. This
uses reference stationarity, not stationarity of the full physical seed.
The physical quantum action has not been replaced by this reference:
its entire difference is still in the accepted full r.

The root-zero leading scalar vanishes under York; s_ref=O(z^2 log^k z)
and is regular at the bolt. Thus it is in the scalar energy domain.
The source in (1) is in Hmu by the preceding L2 bound. Its scalar
Laplacian is also Hmu by (1); no energy-to-point or mass-trace theorem is
being assumed. These facts justify testing (1) by s_ref.

Use the deliberately conservative guards

    F_min=1-1.179e-8, |B_ref|<12.001,
    23.99<V_ref<24.001, D_R lambda<2.946e-9,
    J_ref=12.001*8*B0,
    S0=J_ref/(F_min*23.99), A0=J_ref/sqrt(F_min*23.99),
    D0=(24.001*S0+J_ref/F_min)/6.

Then ||s_ref||Hmu<=S0, ||s_ref||a_ref<=A0, and
||Delta4 s_ref/lambda||Hmu<=D0. The ACTUAL scalar upper form from
current-canonical-interior-schur-bound.md is

    a(s)<=6(2.946e-9)||Delta4 s/lambda||Hmu^2
                       +6.0001 G(s)^2+24.001||s||Hmu^2,
    G(s)^2=int Q z^-2 s'^2/beta dz.

Consequently, in the ORIGINAL Es normalization,

    ||s_ref||Es <= S_ref,
    S_ref^2=max(6.0001/(6F_min),24.001/(23.99F_min))*A0^2
                                             +6(2.946e-9)D0^2. (2)

The numerical contract is S_ref<20. This is not a full energy norm for
Psi_ref: that root-zero field is outside E_total.

## 3. Transfer the complete adjoint to the actual filtered mass functional

Set Psi=Psi_ref-h^-1 r. The accepted coupled inverse has norm <=1, so
its scalar coordinate obeys ||s_Psi-s_ref||Es<=R, R=2.1e-5. The exact
physical mass Green flux remains K_m=2 Z_c>1.9999999764.

For w in the original E core choose the COMPLETE coupled source
G_w=(-b w,0). Both rows of h u=G_w give

    s_u=a^-1 b(xi_u+w), (ell+b*a^-1b)xi_u=-b*a^-1b w,
    xi_u=-F w.

The original York mass multiplier is -3/2. Thus the canonical mass of
u is +(3/2)m_xi(Fw), and the SAME mass Green identity gives

    m_xi(Fw)=-(2/(3K_m)) b(s_Psi,w).                     (3)

For core w, the full finite-domain Green argument keeps the boundary
commutator and common-regulator contacts. The b pairing in (3) is
continuous because s_Psi is in Es and |b(s,w)|<=sqrt(q)||s||Es||w||E,
q=.338. Extend from the original core. This agrees with the already
accepted full graph's mass coefficient; no retained Ft trace theorem
or trace on arbitrary unfiltered energy vectors is inferred.

Let d_F be the E Riesz vector of (3). Let d_F,ref use s_ref in (3).
Then

    ||d_F||E <= K_F=2sqrt(q)(S_ref+R)/(3K_lower),
    ||d_F-d_F,ref||E <= e_F=2sqrt(q)R/(3K_lower),
    g=C_op F d_F,
    ||g||E <= K_g=q K_F/(1+q)^2,
    ||g-C_op F d_F,ref||E <= e_g=q e_F/(1+q)^2.         (4)

These are numerical bounds for the ACTUAL filtered functional, including
the full quantum action through r and the coupled inverse. They are not
numerical applications of those operators or a signed center for g.
Test K_g<1 and e_g<1e-6.

For the original cut reference mass pair chi_m(n_m,b_m), let h_m be
-3/2 times its TT coordinate. Its mass coefficient is exactly one and
V h_m=-3 chi_m b_m/2. Keep chi_m=1 through .1, zero by .2, its entire
transition, original regular bolt and physical conversion 2133<cost<2174.
No numerical derivatives are assigned to an unspecified cutoff profile.
The entire original reference return f_cut=2xz chi_m-2F n_m Dchi_m,
alpha_cut=n_m Dchi_m remains; the cut lift is not a global zero mode.

Use the SAME averaged structured v0,v1 and FULL residuals of e7abfb5:
a0=r0(h_m)=j_r(h_m)-ell(v0,h_m),
a1=r1(h_m)=B(v0,h_m)-(ell+B)(v1,h_m). The mass contribution is EXACTLY

    L_mass=(a0 T+a1 I)g.                                (5)

Keep the signed polynomial before taking norms. Equation (4) gives
||L_mass||<=K_g(q|a0|+|a1|), or the sharper multiplier enclosure
K_F sup_(0<=t<=q) t|a0 t+a1|/(1+t)^2. A numerical approximation using
d_F,ref must add e_F times that SAME supremum. The two coefficients,
the signed vector and their application errors remain unknown.

## 4. Execute a failure test for the available global-gluing estimator

Attempt to enclose (5) with the available complete kernel decomposition.
The actual ROW-SUBTRACTED gluing part, including its diagonal row, obeys
|Q_glue(v,h)|<=alpha_q (6e24)||Vv||L2(dz)||Vh||L2(dz),
alpha_q=2G C_geom/(B^2 beta). This is the SAME normalized coefficient
used in the complete weak mass residual. It includes all modes and actual
returns, not just a short-distance contact. Do not confuse alpha_q with
the original radius/source factor alpha in j_r.

For the prescribed averaged trials the already established whole-spectrum
bounds are ||Vv_k||2<=S_rec L_k, k=0,1. No r=1e-110 change is made.
Suppose H_m is ANY certified upper bound for ||Vh_m||2. The additive
estimator of (5), using (4) and separating these gluing pairings, contains

    E_glue,est=K_g alpha_q (6e24) H_m S_rec(q L0+L1).     (6)

Other local, principal, contact, mixed, source and application errors
are NONNEGATIVE entries in this specified additive estimator. They are
not claimed evaluated. A lower bound on (6) can reject every completion
of that estimator; it does NOT lower-bound the signed residual or error.

No arbitrary cutoff choice is needed for this test. On 0<z<.1, where
the ORIGINAL chi_m=1, use .16<beta<.17, .15<Q<.17, |Q'|<.3 and
J_m<=z^4/(4*.15^2) to obtain b_m>.88 z^3. Therefore necessarily

    H_m >= H_min=1.32 sqrt(.1^7/7).

Substitute H_min into (6), the certified full-spectrum L0,L1,S_rec and
the original alpha_q. Test that this positive estimator expression
exceeds the COMPLETE conditional Dj target 1.49e34 by >1e5. Report
the exact ratio, original bolt contribution and the corresponding
optimistic replacement ceiling for 6e24. That ceiling grants ALL other
errors zero and uses the SMALLEST possible H_m: it is not a margin.
This is an executed test of the remaining estimator's complete budget,
not a certificate for (5) or a declaration that the actual gluing is large.

## 5. What remains unsupported and the method decision

The scalar-coordinate and filtered-mass bounds address the missing
adjoint transfer. The first remaining unsupported numerical implication
is the pair of complete signed residual evaluations a0,a1 on the
original cut lift. Continuity makes them finite at the fixed averaging
radius; it supplies neither their signs nor useful numerical constants.

If (6) fails, do not tighten K_g, cold bounds, degree or transport errors
to rescue this estimate. Compare one specified alternative on the SAME
obstacle: keep the spatial kernel in the two joint residual pairings.
Vh_m is supported in z<=.2; the very large trial L2 norms are concentrated
near the bolt. A useful calculation must retain both local row subtraction
and far correlations, rather than multiply their global L2 norms.
The existing resolvent proof has only fixed-buffer factors e^-M/270 and
e^-M/530. It does NOT yet provide a numerical kernel weighted by the
actual separation of these two profiles. Replacing those factors by
e^-M times that separation without another argument would be invalid.
Common local contacts must remain in a0,a1 even if a far kernel improves.

ONE next calculation: enclose the two COMPLETE signed pairings a0,a1,
using the original lift and a spatially weighted actual-domain resolvent
bound for their gluing term, then evaluate the polynomial in (5) and
all application errors. The executed ratio from (6) quantifies the
minimum gain this particular obstacle requires, before any other debit.
This next calculation remains PROPOSED, not executed by (2)-(6).

The unchanged whole comparison is

    Ebridge+sqrt(1+epsilon)[||L_mass||+Efast+Erest+Eapplication]
                                                        <=1.49e34.

No bracketed unknown is zero. The existing bridge <6.669e13 is charged
once. No signed Dj/Gt/Zr or bolt allowance is closed. Original P32 defect,
all source/graph/physical/boundary comparisons, actual kappa_m, nonlinear
admission/contraction/existence, exact compatibility, matching/opening,
semiclassical validity/noise and all-sector stability remain unresolved.
The limited quantum-fast screen is unchanged. No construction obstruction
or need for new physics follows. Author arithmetic review supports this
ordinary argument; it is not independent or Lean verification.

Remaining to-do list: complete signed pairings and response with errors;
all original coupled/admission and physical/stability requirements; the
full stable semiclassical wormhole objective.
