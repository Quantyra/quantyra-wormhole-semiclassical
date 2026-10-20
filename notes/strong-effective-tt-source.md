# Strong fast forcing for the remaining TT inverse

2026-10-20 HST, E003/S166. Ordinary proof and contract before execution.
**A quantitative reduction of the coupled inverse; the global fast TT
inverse and mass coefficient remain unproved numerically.**

Retain the analytic seed, full physical Hessian, original state, clock,
paired domain and mass column. Use [S165](tt-form-and-scalar-fourth-graph.md)
and the [weighted York conversion](weighted-york-source-conversion.md).
Put g=G lambda, W=(r/B)^2, Tbar=(-Box_4)/lambda,
gradbar=lambda^-1/2 grad, Dbar=D/lambda^2 and Bbar=Bmix/lambda^2.
All norms use dmu=4pi beta r^2 dr. The dimensionless York source is

    (fbar_tau,fbar_c)=lambda^-2 Rcal* f,
    Fsrc=||W fbar_tau||+||W fbar_c||.

For the actual physical compatible solution, the effective TT forcing
f_TT=fbar_c-Bbar* tau obeys

    ||W Bbar* tau|| <6.5e7 Fsrc,
    ||W f_TT||     <6.6e7 Fsrc
                       <=8.184e9 lambda^-2 Fdisk,  (1)
    Fdisk=||W(f_S,f_T)||_L2(dmu;R^2).

This supplies a strong fast source for
c=Dbar_m^-1 f_TT, with the original correlated mass datum.
The constants do not bound Dbar_m^-1 in the required fast fourth/log
norm, its mass coefficient or the actual residual/inverse composition.
Compatibility is retained, not imposed on the preparation source.

## 1. Weighted scalar second derivatives

For the physical compatible tau, S164 gives its weighted value,
gradient and second scalar graph. These also imply

    ||W gradbar^2 tau||
      <=||W Tbar tau||+22||W gradbar tau||+167||W tau||
      <3.3e9 sqrt(g) Fsrc.                        (2)

Here is the domain and constant calculation. Let v=W tau. The regular
scalar Bochner identity with |Ric|<=14lambda gives
||gradbar^2 v||<=||Box_4 v/lambda||+4||gradbar v||.
This is the scalar specialization of the Ricci identity already used
in S165; see [equation (7)](https://arxiv.org/pdf/1708.04775).
The seed has R^2<=5lambda and 0<=AR<=11lambda. Therefore

    |gradbar W|/W<4.5,
    |Box_4 W|/(lambda W)<=74,
    |Hess W|/(lambda W)<=75.

The Hessian's distinct frame entries divided by W are
2AR+2R^2, 2AR, 2R^2, 2R^2, bounded by
32lambda,22lambda,10lambda,10lambda.
Expand Box(W tau), grad(W tau) and Hess(W tau).
The gradient coefficient is 4(4.5)+4=22 and the value coefficient
74+4(4.5)+75=167, proving the first inequality in (2).

The already proved weighted gradients make v and grad v square
integrable, while the product formula makes Box v square integrable.
The regular scalar graph realization and cutoff approximation give
the Bochner inequality on this domain. In particular Hess v is
obtained as a weak limit; it is not assumed to justify its own bound.
The bolt is an ordinary smooth interior point and W is smooth there.
For the second inequality in (2), use S164 and g^(1/4)<6e-9:

    3.2e9+22(4e8)(6e-9)+167(4.4e7)sqrt(g)<3.3e9.

## 2. Curvature derivatives without a bolt singularity

Let P be the orthogonal projector onto the orbit two-plane,
S=g_metric-P, and R=sqrt(F)/r, A=F_r/(2sqrt(F)).
In an orthonormal frame, the only nonzero first derivatives are

    (nabla_alpha P)_(rho beta)
      =(nabla_alpha P)_(beta rho)=R delta_alpha_beta,

where alpha,beta are sphere directions. A second differentiation gives

    |nabla P|^2=4R^2,
    |nabla^2 P|^2=4(R')^2+4(AR)^2+24R^4.          (3)

The radial, time and two sphere derivative contributions are
4(R')^2, 4(AR)^2 and 24R^4. This is a full covariant derivative
calculation, not a derivative of frame components alone.
Since R'=AR-R^2, |R'|<=16lambda. Equation (3) gives, conservatively,
|nabla^2 P|<300lambda, including the regular bolt.

Use the exact outgoing-proof curvature decomposition

    J=j(P-S), C=w_C W0,
    W0=-P wedge P-S wedge S+P wedge S,
    nabla W0=-3(P-S) wedge nabla P.

Here |P-S|=2, |W0|=sqrt(48)<7,
|nabla W0|<=48R<108sqrt(lambda), and
|A wedge B|<=4|A||B|. Differentiation gives
|nabla^2 W0|<=24|nabla P|^2+24|nabla^2 P|<7680lambda.

The seed polynomials for j/lambda and w_C/lambda have coefficients
in u=B/r with degrees 4..6 and 3..6, respectively. S165 proved
coefficient sums 3.6,20 for jbar,u jbar_u and 5.1,40 for
wbar,u wbar_u. Hence u^2 jbar_uu and u^2 wbar_uu have sums
at most 108 and 153. The scalar Hessian entries are

    F p_rr+(F_r/2)p_r, (F_r/2)p_r, F p_r/r, F p_r/r.

Using F/r^2<=5lambda and F_r/(2r)<=11lambda gives
|Hess j|<1500lambda^2u^4 and |Hess w_C|<3000lambda^2u^3.
The existing sharper first derivative bounds are
|grad j|<=30lambda^(3/2)u^4 and
|grad w_C|<=43lambda^(3/2)u^3. Therefore

    |nabla^2 J|<1e4lambda^2u^4,
    |Box J|<2e4lambda^2u^4,
    |nabla^2 C|<1e5lambda^2u^3.                   (4)

For example the J budget is 2(1500)+2(30)(9)+3.5(600)<1e4;
the C budget is 7(3000)+2(43)(108)+5.1(7680)<1e5.
All connection terms are included by (3). No norm of A alone occurs.

## 3. Formal transpose of the complete local trace

The [exact mixed block](stationary-weighted-tt-inverse.md), section 5,
is expressed through the genuine TT lift H:c->h_c and return V:c->eta:

    Bmix* tau=-(8piG)^-1[H* Ttrace* tau+V*(E_trace grad tau)]. (5)

Use the given TT trace expression as a differential operator on
arbitrary test tensors to define Ttrace*. Applying H* then recovers
the original TT restriction. This extension adds no physical sector.
The fixed-density and off-shell return terms in (5) are essential.

There is a useful derivative cancellation: r_dot=-J:h and the
epsilon Box r_dot term transpose directly to -epsilon J Box tau.
No derivatives of J arise from that term. The J:Delta_L h term
transpose is Delta_L(tau J); its product rule uses (4).
The Weyl derivative transpose is a sum of second covariant derivatives
of tau C, with norm at most

    6[|C||Hess tau|+2|grad C||grad tau|+|nabla^2 C||tau|].

The algebraic Weyl cost is the existing
12320lambda^2u^3|tau|. These estimates follow by moving derivatives
in the exact connection-variation formula, with its two tensor slots.
Using |C|<=40lambda u^3, |grad C|<=1000lambda^(3/2)u^3,
|J|<=7lambda u^4 and |grad J|<=100lambda^(3/2)u^4, all terms give

    lambda^-1 |Ttrace* tau|
      <=(7u^4+20v u^4)|tau|
         +e[200u^3|gradbar^2 tau|
               +8000u^3|gradbar tau|+5e5u^3|tau|], (6)

where e=epsilon lambda=g/(60pi) and
v=epsilon q^2/(lambda B^4)<4e-27/180.
The Hessian coefficient before rounding is at most
14+(7/12)6(40)+(11/6)14<200.
The gradient coefficient is
(7/12)12(1000)+(11/6)2(100)<8000.
For the value coefficient include

    7(11)(13)/36+4+(7/12)[6(1e5)+12320]
      +(11/6)[20000+7(2*54+2*14)+(13/2)7+2(7^2)]<5e5.

This includes the trace prefactor, Hess R4, Weyl, Lichnerowicz,
scalar-curvature and quadratic-J contacts. The magnetic term is
retained separately. The nonlocal massive determinant remains in D.

## 4. Weighted transpose bounds and actual TT forcing

The reciprocal-weight TT integral from S166.1 gives
||W^-1 b||<=2||W^-1 c||. Thus

    ||W^-1 Hc||<=sqrt(22)||W^-1 c||<5||W^-1 c||.

The return bound gives
||W^-1 sqrt(lambda)Vc||<=9||W^-1(b+c)||<=27||W^-1 c||.
Duality in the fixed physical pairing yields

    ||W H* q||<=5||W q||,
    ||W Vbar* q||<=27||W q||, Vbar=sqrt(lambda)V. (7)

By (5), after the exact dimensionless normalization,

    Bbar* tau=-(8pi g)^-1[
       H*(lambda^-1 Ttrace* tau)+Vbar*(Ebar gradbar tau)],
    Ebar=E_trace/lambda, |Ebar|<6e-22u^4.

Multiplication by W and (6),(7) therefore give

    ||W Bbar* tau|| <=(8pi g)^-1 {
      5[(7+20v+5e5e)||W tau||
           +200e||W gradbar^2 tau||+8000e||W gradbar tau||]
                     +27(6e-22)||W gradbar tau|| }. (8)

No derivative of the trace defect is required in its transpose
contact. Insert (2) and S164's value and gradient bounds. With
alpha=sqrt(g), 3<pi<22/7,
g<2000^-10 and g>.999*2000^-10, one has g^(-1/4)<2e8.
The coefficient of Fsrc in (8) is bounded by

    (1/24) {
      5[(7+20v+5e5e)(4.4e7)
                  +(200/180)(3.3e9)alpha
                  +(8000/180)(4e8)g^(3/4)]
                  +27(6e-22)(4e8)(2e8) }<6.5e7.

The lower g guard is used only in the final return term; omitting it
would not justify that numerical reserve. Adding ||W fbar_c||<=Fsrc
gives the second inequality in (1). S166.1's conversion
Fsrc<=124lambda^-2 Fdisk gives the last.

Initially all transposes are paired on smooth compact test fields.
The weighted bounds extend them by density, and (8) identifies the
result as a strong fast L2 source on the original realization. The
actual compatible tau and correlated TT mass remain the same fields.
The operator Dbar_m^-1 still needs its numerical global fourth/log
and mass estimate before (1) becomes a complete metric inverse bound.

## Frozen verification contract

Freeze this note and checker before execution. Revalidate the exact
receipts/current hashes for the weighted TT inverse, outgoing scalar
inverse, S164, S165 and S166.1.

Compute the first and second covariant derivatives of the orbit
projector independently in the coordinate metric
diag(a(rho)^2,1,r(rho)^2,r(rho)^2 sin(theta)^2).
Contract with the inverse metric and verify both identities (3).
Verify scalar Hessian and Laplacian product rules, the W entries,
and the curvature coefficient envelopes by exact arithmetic.
Check the transpose of Box acting on a coefficient times a test
field by symbolic integration-by-parts identities, and retain the
physical reciprocal-weight transpose convention.

Check every displayed rational budget, including the lower-g use,
the Hessian, trace-transpose and final effective-source constants.
No numerical seed scan or actual preparation-source fixture is used.
The geometric/domain limits, product-rule norm bounds and signed
compatibility remain ordinary proofs. Finite checks do not certify
the complete TT inverse or nonlinear existence. Preserve every
scientific run and exact failed sources if any.

Remaining to-do list: S166's quantitative global fast fourth/log TT
inverse and mass coefficient, higher reconstruction and full residual
composition; nonlinear family and signed compatibility; exterior/state,
charge and return; opening, stability, transport, preparation, creation
and full-result audit.
