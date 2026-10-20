# Classical disk-chart positivity at the actual finite comparison

2026-10-20, E003/S167 supporting S139/S074. Ordinary argument before
frozen verification. **Component progress; the full goal remains open.**

The off-shell contact left by the
[exponential TT theorem](finite-comparison-classical-tt-form.md) costs
less than .06 of the candidate TT energy. Thus, throughout the actual
finite path, for a genuine TT tensor h and its outer-fixed disk return k,

    I0,disk''[k,k] >= .93 E_ell[c],
    (I0+Gamma_hat,LLL),disk''[k,k] >= .92 E_ell[c].       (1)

I0=S_classical+Gamma_ref is the fixed reference-potential theory. The
second line includes only the retained canonical LLL allocation. The
remaining massive/Weyl/reference terms, mixed/density blocks and
nonlinear solve are still required for the full physical theory.

## 1. Geometry and the radius equation

Keep the original physical temperature, clock, flux, state/domain,
paired tower, common subtraction and fixed potential. Use exactly the
finite path and notation of the exponential TT proof:

    g_ell=e^(2w)(d rho^2+a0^2 dtheta^2)+r0^2 e^(2S)dOmega^2,
    y=sqrt(lambda)rho, t=r0/B, x=lambda B^2, lambda=lambda_e,
    S=log T, T=1+tau(e^(ell V)-1), tau=(t-1)/t,
    nu=tau V, w=ell u+S, |w|,|S|<eps=3e-8.

A0=(log a0)_y and R0=(log r0)_y are dimensionless seed quantities.
The preceding proof supplies

    |u|,|nu|,|V|<1.2e-8, |nu_y|<7e-7,
    |S_y|<7.1e-7, |A0 S_y|<3.45e-6, |w_y|<1.03e-6,
    |u_yy+A0 u_y|<3.05e-7/t^2,
    |Delta Rtilde^2|<3.2e-6, |Delta ARtilde|<6e-6,
    Rtilde=R0+S_y, ARtilde=(A0+w_y)(R0+S_y).

The invariant radius equation in the
[completed metric proof](stationary-completed-metric-response-norm.md)
is, with phi0=log t and A_phys=B/sqrt(C),

    Delta_hat nu+2 grad(phi0).grad(nu)
      =(delta Rhat_common-t^2 delta R4_common)/6
        +2u[Delta_hat phi0+|grad phi0|^2]-(nu/3)t^2 R40.

Divide by t^2 lambda, retaining both common-coordinate curvature
variations. Since the canonical orbit Laplacian pulls back to
t^2 lambda[partial_y^2+(A0-R0)partial_y], this becomes

    nu_yy+(A0+R0)nu_y
      =delta Rhat_common/(6t^2 lambda)-delta R4_common/(6lambda)
        +2u E20-(nu/3)R40/lambda,
    E20=R0,y+A0 R0=2A0 R0-R0^2, |E20|<=22.           (2)

The already retained local trace and curvature data are

    |A_phys^2 delta Rhat_common|<5.2e-8, C/x<6.1,
    delta R4_common=-8piG Arest+nu r partial_r R40,
    |8piG r0^2 Arest/C|<1e-20,
    R40/lambda=-12+(d/x)(4t^-5-6t^-6), x>=1000, |d|<=1.

The trace bound is the complete vacuum-cancelled polynomial from the
[constraint residual proof](stationary-centered-constraint-residual.md),
not a differentiated unknown quantum remainder. In particular

    |R40|/lambda<12.01, |r partial_r R40|/lambda<=.056,
    |delta R4_common|/lambda<6.1e-20/t^2+6.72e-10.

The right side of (2) has modulus less than 6.3e-7. With
|R0 nu_y|<2.24(7e-7), this gives

    |nu_yy+A0 nu_y|<2.3e-6.                          (3)

Put k_e=e^(ell V), B_e=k_e-1-ell V k_e. The exact finite-radius
identity, not a linear approximation, is

    T_yy+A0 T_y=ell k_e(nu_yy+A0 nu_y)
          +B_e(tau_yy+A0 tau_y)+ell^2 k_e tau V_y^2,
    S_yy+A0 S_y=(T_yy+A0 T_y)/T-S_y^2.              (4)

Here k_e<1.001, T>.999, |B_e|<2e-16. The seed bounds give
|tau_yy+A0 tau_y|<=22 and tau V_y^2<=7 b_gap^2,
b_gap=(t-1)V_t, |b_gap|<3e-7. To see the latter, use
F0<=F_ext and
F_ext/[x t(t-1)]=(1-t^-1)[1+2/t+(3+1/x)/t^2+d/(x t^4)]<7.
Thus

    |S_yy+A0 S_y|<2.31e-6.                         (5)

This is the precise second derivative needed below. It uses source
values and the invariant equation, with no new source derivative.

## 2. A uniform bound for the classical gradient

In this and later sections primes mean actual physical proper-radius
derivatives and all unmarked norms use the actual dmu_ell. Set

    A=a_ell'/a_ell, R=r_ell'/r_ell,
    K=a_ell''/a_ell, E2=R'+AR, E1=K-E2.

The exact pullbacks and seed curvature bounds give

    K/lambda=e^(-2w)[a0,yy/a0+ell(u_yy+A0 u_y)+(S_yy+A0 S_y)],
    E2/lambda=e^(-2w)[E20+(S_yy+A0 S_y)],
    |K|<10.501lambda, |E2|<22.001lambda.            (6)

Also A>0, R>0 and R<2.25sqrt(lambda). Indeed A0>1/sqrt(2)
and |w_y|<1.03e-6. The exact radius derivative is
r_ell,y=e^(ell V)r0,y(1+ell b_gap)>0. These statements extend
continuously at the bolt without asserting that A itself is bounded.

Define the physical action gradient by

    I0'[v]=1/(8pi G) integral e:v dmu_ell.

In the orthonormal (proper radial,time,sphere,sphere) frame, direct
first variation of the retained radial action gives e=diag(e_n,e_t,e_s,e_s):

    e_n=AR+R^2/2-V0/(4r_ell^2),
    e_t=E2-AR+3R^2/2-V0/(4r_ell^2),
    e_s=[E1+2E2+R^2-V0,X/2]/2.                    (7)

The potential V0(X) is held fixed, X=r_ell^2. These components
vanish at the seed. They follow from
I0=-(pi/G) integral [a r'^2/N+2a'rr'/N+aN V0/2]d rho
with its retained endpoint terms. For regular fixed-period variations,
the bolt entropy cancels the surface variation, as in the preceding
proof. No finite-comparison field equation is imposed in (7).

For U(r)=V0(r^2)/(4lambda r^2), H(r)=V0,X(r^2)/(2lambda), the
explicit potential of the preceding proof gives, for every r>=B,

    |U|<4, |d_logr U|<7, |H|<7, |d_logr H|<13,
    |e^(2w)U(r_ell)-U(r0)|<4.6e-7,
    |e^(2w)H(r_ell)-H(r0)|<8.2e-7.                 (8)

For example the sharper coefficient sums are
|U|<=3+5/x, |d_logr U|<=6+24/x,
|H|<=6+14/x, |d_logr H|<=12+74/x. The finite logarithmic
radius change is |S|<eps and e^(2eps)-1<6.1e-8.

Subtract the exactly zero seed components in (7) after multiplying
by e^(2w)/lambda. The preceding coefficient differences, (5) and (8)
bound the three components by

    |e_n|<8.2e-6 lambda,
    |e_t|<1.36e-5 lambda,
    |e_s|<4.5e-6 lambda,
    ||e||_infinity<1.8e-5 lambda.                  (9)

The last norm is the four-component Frobenius norm, including both
sphere components. It bounds the complete gradient of I0, not the
gradient of the full physical action.

## 3. Cancellation of the combined gauge and chart contact

Use the genuine TT input and its original outer-fixed return:

    h=diag(b,-b-2c,c,c), b'=-2(A+R)b-2(A-R)c,
    Y=eta partial_rho, eta'-A eta=b+c,
    k=h-L_Yg, k1=2h-L_Yg.

The previous proof identified the exact contact tensor

    Z=L_Y(2h-L_Yg)+k^2-h^2,
    I0,disk''[k,k]-Q_exp[h]=I0'[Z].                 (10)

There is an additional exact cancellation in the combined Z. In this
diagonal relative frame set s=(eta',A eta,R eta,R eta), so
(L_Yg)_i=2s_i. For any diagonal covariant tensor v with relative
components v_i, (L_Yv)_i=eta v_i'+2s_i v_i. Consequently

    Z_i=eta(2h_i-2s_i)'+2s_i(2h_i-2s_i)
                          +(h_i-2s_i)^2-h_i^2
       =2eta(h_i'-s_i')=eta k1,i'.                (11)

Thus the two exponential accelerations cancel the nondifferentiated
products from the gauge contact. The estimate below needs no spatial
derivative of e. In this radial parallel orthonormal frame h_i' and
k1,i' are their radial covariant derivatives.

## 4. Bound the remaining radial derivative

Put Cc=||c||, D=||c'+2Rc|| and
e_c=1.5D^2+lambda Cc^2=8pi G E_ell[c]. The preceding TT lift and
measure comparison give

    ||b||<2.002 Cc, ||h||<4.8 Cc,
    ||h'||<75 sqrt(e_c), ||eta||<4.83 Cc/sqrt(lambda). (12)

For the first bound multiply the seed-measure bound 2.001 by
e^(4eps). For the third use the retained E_h<5500 e_c; for the
fourth use ||eta||<1.006||h||/sqrt(lambda). No derivative of
the unknown nonlinear correction is inferred by these input bounds.

Let d1=b+c and X_A=||A d1||. Multiply the exact TT equation
d1'+2(A+R)d1=c'+4Rc by A d1 and integrate in the physical measure:

    integral [2A^2+AR-K/2]d1^2=integral A d1(c'+4Rc),
    X_A^2<=||c'+4Rc||^2/4+(10.501/2)lambda||d1||^2.

This uses AR>=0 and Young's inequality. The bolt term vanishes
because d1=O(rho^2), while the permitted outer tail has zero
surface term. Since ||d1||<3.002 Cc,

    X_A <= D/2+[2.25+sqrt(10.501/2)(3.002)]sqrt(lambda)Cc
         <9.6 sqrt(e_c),
    ||c'||<5.32 sqrt(e_c).                         (13)

The exact regular derivative identities are

    eta''=K eta-A d1-2R(b-c)+c',
    (A eta)'=K eta+A d1,
    (R eta)'=E2 eta+R d1.                          (14)

For the four-component s' split these into the eta terms, the
A d1 terms, the R terms and c'. The first coefficient is
sqrt(2K^2+2E2^2)<34.5lambda. The exact TT point norm
|h|^2=2(b+c)^2+4c^2 gives
|b-c|^2<=1.5|h|^2 and |d1|^2<=.5|h|^2; hence the R terms
have norm at most sqrt(7)||Rh||. It follows that

    ||s'|| <[34.5(4.83)+sqrt(2)(9.6)
                    +sqrt(7)(2.25)(4.8)+5.32]sqrt(e_c)
           <215 sqrt(e_c).                       (15)

Only the polar product A d1 is bounded; no estimate treats A as
uniform at the bolt. From (11)-(12),(15),

    integral |Z|dmu_ell
       <=2||eta||(||h'||+||s'||)
       <2(4.83)(75+215) Cc sqrt(e_c)/sqrt(lambda)
       <3000 e_c/lambda.                         (16)

Regular bolt parity makes Z=O(rho^2), including the sphere
component, so this first variation has no omitted entropy term.
All calculations start on smooth regular TT data with the selected
outer behavior and extend through the retained H1 form closure.
The inequalities are strict for nonzero inputs; zero maps to zero.

## 5. Consequence and verification scope

Combining (9),(16) with the physical normalization,

    |I0'[Z]|<(1.8e-5)(3000) E_ell=.054 E_ell<.06 E_ell.

The established Q_exp>=.99 E_ell therefore proves the first line
of (1). The already established canonical LLL disk bound
Gamma_hat,LLL,disk''>=-8e-24 E_ell gives the second line. Both
forms now use the same actual disk chart and input energy.

The frozen checker pins the preceding exponential TT receipt and its
input union, including both invariant curvature source notes. It
verifies the radius-equation coordinate conversion and finite-radius
identity, first variation of the action and zero seed gradient,
potential bounds, combined contact cancellation, polar integration
and derivative identities, and every numerical margin. The action
gradient and contact algebra are checked directly, not just asserted
from a precomputed norm. Continuum domain and closure statements
remain ordinary proofs; this is not a formal proof or a numerical
solution of the full quantum equations. Every scientific run and any
failure/source snapshot is retained.

Next transfer the remaining massive/Weyl/reference TT form and
mixed/density terms to the candidate. The fixed-potential split must
be respected when recombining the full action. Complete logarithmic
domains and the coupled candidate range before claiming a full inverse.
Pi_cut stays fixed at the seed. All S167 tasks stay unchecked and
no physical gate is promoted.

Remaining to-do list: remaining massive/Weyl/reference TT and mixed/
density forms; logarithmic domains and coupled range; fixed projection
direction and cold/higher response; nonlinear preservation and projected
solve; signed compatibility; exterior/state, charge and return; opening,
stability, transport, preparation, separate creation and full-result audit.
