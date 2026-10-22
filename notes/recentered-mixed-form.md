# Charge-weighted TT energy and the recentered mixed form

2026-10-22, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress, technical c93a9e3 and planning 0815480.
The full stable semiclassical wormhole objective remains incomplete.

The actual scalar/returned-TT mixed form at the SAME c_R=10^22
recentered RN reference obeys the coordinate-independent relative bound

    |b(s,xi)|^2 < .336 a(s,s) l(xi,xi).                 (1)

Here a is the positive conformal TRACE comparison form (the negative
of the action's conformal Hessian), l is the full returned-TT action
form, and both include the retained density/chart contacts. Initially
the fields are smooth at the bolt and compact at the outer end; the
bound extends in the stated energy/form completions below. This is
strictly inside the selected 3/8=.375 mixed-form target.

Equation (1) is not yet the claimed physical inverse product: if l has
the required coercive physical realization, it implies that product is
<=.336 in the explicit comparison spaces. The complete metric/current/
constraint domain and its admission-norm inverse are still unproved.
The present result changes the next priority from estimating an unknown
mixed constant to establishing that physical realization.

The improvement comes from retaining the radial charge potential in
the TT energy, not from rescaling fields or norms. The coarser TT gap
throws away precisely the weight that controls delta R=-4(Z/r^4)xi.
The direct form statement (1) makes this distinction independent of
the normalization of either variable.

## 1. Fixed reference, scalar form and exact mixed row

Keep the geometry, physical period beta_T, fixed magnetic two-form,
paired Dirac domain, c_R=10^22 and mu0 of
[the quantum TT audit](recentered-quantum-tt-form.md). Use

    g=G lambda, epsilon=G/(60pi), D=96piG c_R-epsilon,
    C_v=1+11epsilon lambda/3, F_lambda=1-384pi c_R g,
    f=1-2M/r+Z/r^2+lambda r^2,
    t(r)=Z/r^4, J=diag(-t,-t,t,t), R0=-12lambda.

The scalar density operator and its lower energy are

    A=6(-Delta+4lambda)(-D Delta+C_v)
                             -(22epsilon/3)div(J grad),
    a(s,s)=(8piG)^-1 int s A s dvol >= S[s]^2,
    S[s]^2=(8piG)^-1 int [6D(Delta s)^2+6|grad s|^2
                                             +24lambda s^2]dvol. (2)

The factor 8piG is the physical action normalization, displayed on
BOTH sides of every estimate. It does not change a Schur product.
The complete-end scalar comparison realization from the recentered
trace proof supplies its energy inverse; its identification as part
of the full physical metric domain remains separate.

For a genuine TT tensor k=diag(-b-2xi,b,xi,xi),

    (r^2 f b)'=-(r^2 f'-2rf)xi,
    delta R[k]=-J:k=-4t xi.                            (3)

Return to the fixed disk slice with vector eta partial_rho,
eta=-a_orb int_rho^infinity (b+xi)/a_orb d rho, a_orb=sqrt(f)/kappa.
The exact mixed density row, specializing equations (14)-(17) of
[the weighted TT analysis](stationary-weighted-tt-inverse.md), is

    B_tr xi=(C_v-D Delta)delta R[k]
       +epsilon[(7/12)delta C4^2[k]+(11/6)delta J^2[k]
                                                    -20delta I[k]]
       -div(T0 eta),                                  (4)
    b(s,xi)=-(8piG)^-1 int s B_tr xi dvol,
    T0=epsilon[-10q^2/r^4+28M^2/r^6-56MZ/r^7
                                             +(106/3)Z^2/r^8].

I=F_em^2=q^2/(2r^4), with the physical two-form held fixed.
The canonical determinant has zero conformal/mixed derivative in the
fixed quotient variables; the anomaly variations in (4) remain.
The final divergence is essential. Scalar covariance gives
delta_Ttrace[L_eta g]=eta.grad T0, while the density contact on
k-L_eta g is -T0 div eta. Together they give -div(T0 eta), even
though the reference is off shell. No source derivative or full-action
stationarity is assumed when bounding it.

## 2. The stronger TT comparison that controls the charge weight

The exact classical TT identity and Hardy cross-term bound give

    l_classical >= F_lambda/(8piG) int {
         (3/2)f(xi'+2xi/r)^2+V(r)xi^2}dvol,
    V=3lambda-12M/r^3+13Z/r^4-7/(3r^2).                (5)

The prime is areal-r differentiation; dvol=4pi beta_T r^2 dr.
This is the retained stationary effective Einstein-Maxwell comparison,
with the nonlocal TT b tail kept. Write x=lambda B^2>8e9,
u=B/r in (0,1), M=B(1+2x+delta/2), Z=B^2(1+3x), 0<delta<1.
We now prove the stronger uniform inequality

    V >=1.99 t^2/lambda+.08lambda.                     (6)

Indeed V/lambda-1.99(t/lambda)^2 is bounded below by
P(u)-33/x-2/x^2, where

    P(u)=3-24u^3+39u^4-(1791/100)u^8.

For transparency, the degree-eight Bernstein coefficients of P on the
two affine intervals are the following exact rational certificate:

| i | [0,1/2] | [1/2,1] |
| --- | --- | --- |
| 0 | 3 | 60609/25600 |
| 1 | 3 | 30609/12800 |
| 2 | 3 | 114063/44800 |
| 3 | 165/56 | 64563/22400 |
| 4 | 3159/1120 | 38253/11200 |
| 5 | 591/224 | 22713/5600 |
| 6 | 549/224 | 12513/2800 |
| 7 | 75/32 | 759/200 |
| 8 | 60609/25600 | 9/100 |

Bernstein basis functions are nonnegative and sum to one, so P>=.09
on the whole interval. The omitted negative coefficients have sum
12+6delta+7/3+6(1.99)<33, and the x^-2 term is <2. Since x>8e9,
.09-33/x-2/x^2>.08 proves (6). This is not a sampled radial minimum.

The preceding quantum/contact loss is <.005 E_eff and the additional
local square loss <7.07e-8 E_eff, where

    E_eff=F_lambda/(8piG) int [(3/2)f(xi'+2xi/r)^2+lambda xi^2].

Both disk and affine TT versions fit this conservative allowance: the
previous ledger already includes the full off-shell return/chart cost.
The local square has zero gradient and delta R[L_eta g]=0. Subtracting
those losses from (5)-(6) proves l(xi,xi)>=X[xi]^2 for

    X[xi]^2=F_lambda/(8piG) int {
       .99(3/2)f(xi'+2xi/r)^2
                +1.99 t^2 xi^2/lambda+.07lambda xi^2}dvol. (7)

This lower norm is fixed by a proved physical potential. In particular
X^2>=.07 E_eff. The stronger-form quantum completion embeds in X;
this inequality does not assert that the actual form is closed in X.

## 3. Leading local mixed term: a direct bound below the target

Integrating the Laplacian in (4) onto s gives the leading pairing
of -4t xi with C_v s-D Delta s. Two-component Cauchy--Schwarz in
the scalar energy (2) gives

    ||C_v s-D Delta s||_2
       <=sqrt(8piG) sqrt(C_v^2/(24lambda)+D/6) S[s],
    ||t xi||_2<=sqrt(8piG lambda/(1.99 F_lambda)) X[xi].

Thus its mixed constant satisfies

    K_main^2 <= [2/(3*1.99 F_lambda)](C_v^2+4Dlambda)
                   <.33501 <.579^2.                   (8)

No cancellation of a signed source was assumed. Substitutions use
g<9.765625e-34, F_lambda>1-1.179e-8, Dlambda<2.946e-9 and
C_v=1+11g/(180pi). Keeping the TT potential rather than merely its
constant gap is the quantitative gain. A rescaling of (2) or (7)
would cancel between the two mixed factors and cannot create (6).

## 4. Remaining anomaly terms and the off-shell return

The remaining estimates use full tensor norms in physical four-volume.
Let R_p=sqrt(f)/r. The reference has R_p^2<=5lambda,
M/r^3<=(2+beta0)lambda, t<=(3+beta0)lambda,
beta0=1.223e-10. Set W=M/r^3-Z/r^4. The sectional curvatures are

    K_01=-lambda+2W-t,
    K_0a=K_1a=-lambda-W,
    K_23=-lambda+2W+t.

Consequently safe Frobenius bounds are

    |Rm|<60lambda, |Ric|<13lambda, |J|<7lambda,
    |C4|<35lambda, |P_Schouten|<9lambda,
    |grad J|<100lambda^(3/2),
    |grad C4|<1e4lambda^(3/2).                         (9)

The derivative bounds do not estimate a singular bolt connection.
If P_s is the sphere projector and Q_o=1-P_s, then
J=t(P_s-Q_o), |grad P_s|=2|R_p| and t'=-4R_p t.
Also C4=W(P_s wedge P_s+Q_o wedge Q_o-P_s wedge Q_o), where
wedge is the Kulkarni--Nomizu product. Its unit coefficient tensor
has norm sqrt(48), and its derivative is <=48sqrt(2)|R_p|.
With |W'|<41lambda^(3/2), (9) follows, with ample room. Rotations
within the orbit plane leave these projectors unchanged at the bolt.

For clarity the local variation identities used are

    delta J^2=J:Delta_L k-(R0/2)J:k-2k:J^2,
    Delta_L k=-Delta k-2Riem(k)+Ric o k+k o Ric,
    K^a_bc=(nabla_b k_c^a+nabla_c k_b^a-nabla^a k_bc)/2,
    Rdot_abcd=k_ae R^e_bcd+g_ae(nabla_c K^e_db-nabla_d K^e_cb),
    delta C4^2=2C4^abcd Rdot_abcd-8C4^abcd k_ac P_bd
                               -4k^ae C4_abcd C4_e^bcd. (10)

Integrate the one divergence in each curvature variation onto s and
the background curvature. Since |K|<=(3/2)|nabla k|, (9)-(10) imply
for EACH of delta C4^2 and delta J^2 the safe bound

    |int s delta(curvature square)| <=1e5[
       lambda ||grad s|| ||grad k||
       +lambda^(3/2)||s|| ||grad k||+lambda^2||s|| ||k||].

For example the Weyl derivative coefficients are <=6|C4| and
6|grad C4|; its algebraic coefficient is at most
2(35)(60)+8(35)(9)+4(35)^2=11620, all within 1e5.
The corresponding Ricci-tracefree coefficients are 7,100 and <2000.
Multiplying by 7/12 and 11/6 leaves the combined reserve 3e5.

The retained TT reconstruction gives ||k||^2<=22||xi||^2 and
||grad k||^2<=5200 E, E=int[ f(xi'+2xi/r)^2+lambda xi^2].
Using (7), therefore,

    ||k||^2 <=(8piG/F_lambda)(315/lambda) X^2,
    ||grad k||^2 <=(8piG/F_lambda)75000 X^2,
    ||eta||^2 <=(8piG/F_lambda)(258/lambda^2) X^2.       (11)

The curvature-anomaly mixed constant is consequently at most

    K_curv <=(3e5 epsilon lambda/sqrt(F_lambda))
          [sqrt(75000/6)+sqrt(75000/24)+sqrt(315/24)]
        <3e-28.                                      (12)

For the fixed magnetic two-form, delta I[k]=-q^2 xi/r^4 exactly.
Using its SAME t weight in (7),

    K_EM <=20epsilon q^2/[Z sqrt(24*1.99 F_lambda)]
           <3e-29,                                   (13)

since epsilon q^2/Z=g_R^2 F_lambda/(60pi^2). The geometric charge
Z has not been substituted for the physical flux in the anomaly.

Finally the exact trace source in (4), rather than a generic full-stress
value bound, gives

    ||T0||_infinity/lambda <2.1e-28.

One sufficient substitution is epsilon lambda times
10U+28(2+beta0)^2+56(2+beta0)(3+beta0)+(106/3)(3+beta0)^2,
where U<=g_R^2(3+beta0)/(pi g). Integration of the divergence in
(4) uses no derivative of this source and gives

    K_return <=(||T0||_infinity/lambda)
                            sqrt(258/(6F_lambda)) <1.5e-27. (14)

Equations (12)-(14) sum to <2e-27. Thus

    |b(s,xi)| <=(sqrt(.33501)+2e-27) S[s] X[xi],
    |b(s,xi)|^2 <.336 S[s]^2 X[xi]^2
                          <=.336 a(s,s)l(xi,xi),       (15)

which proves (1) while retaining every term in (4).

## 5. Boundaries, inverse interpretation and next decision

The integrations above start on a smooth bolt-regular compact-outer
core. The bolt is an interior point. For the selected physical mass/
fast TT tails and scalar r^-4 log^p r tails, all new integration-by-parts
fluxes vanish: curvature-mixed fluxes have at least r^-7 log^p r
physical boundary weight, and the leading t-weighted scalar flux
decays faster. The complete local Weyl mixed boundary functional was
already evaluated, including radial metric variation, in
[the boundary audit](recentered-canonical-tt-audit.md). Its limit is
zero on this class. It is not zero on arbitrary finite-cut Navier data.

The bounded bilinear form in (15) extends uniquely to the scalar S
completion and the TT X completion; its restriction to the stronger
quantum form completion agrees with the original mixed action by the
retained TT-preserving approximation. This does not complete arbitrary
outer shape-shape momenta or select the full physical boundary action.

In the natural action pairing the two diagonal signs are important:
the Hessian has blocks [-a,b;b*,l]. Multiplying the scalar equation
by minus one gives the positive trace row. Signs and 8piG factors
must accompany that operation; no Euclidean minimum or physical
stability conclusion is inferred. In a justified symmetric saddle
realization the scalar Schur sign is favorable, but such a realization
is not established for all physical constraints here.

If the required physical TT realization has inverse with bound one
from X* to X (as l>=X^2 would give for its coercive form realization),
then the explicit comparison norms give

    ||L^-1 Cmix|| ||Bmix|| <=.336 <3/8.                (16)

With the still-unproved transfer allowance ||E||<=1/4, the absolute
Schur loop is <=(4/3).336=.448<1/2. A conditional scalar inverse
bound is (4/3)/(1-.448)<2.42. Equation (16) is conditional on the
actual realization; the verified statement independent of that inverse
is the relative form domination (15). These are not bounds in the
complete metric-admission norm or on the actual nonlinear correction.

DECISION: retain the c_R construction and the charge-weighted form.
The mixed numerical target is no longer an unexplained constant in
the reference form problem. Next establish the physical renormalized
coupled domain, actual closed operator and constraint-compatible inverse.
Do not extend the old uniform source-constant chain or declare a new
Friedrichs comparison to be the physical boundary problem by fiat.

The nonlinear certificate still requires a complete admission norm,
source mapping and all positive margins: at rho=1e-8 the sufficient
targets remain Y0<=5e-9, Z0<=.25, Z2<=5e7. None of that full ball
is certified here. The scalar particular <1.82e-29 and the failed old
interior allowance 9.1591500000000000015e103 retain their previous
scopes. Exact compatibility, fixed-flux/temperature matching, the actual
finite state and geometry, omitted-operator/noise control and all-sector
causal stability remain required. The full objective is not reduced to
positivity or invertibility of a reference comparison.

## 6. Frozen supporting verification contract

Freeze this note, checker and review before running. Check the uniform
potential inequality by exact Bernstein reconstruction, all physical
budget substitutions by rational arithmetic, the tensor/projector
contractions and density contact algebra symbolically. At 70/100 digits
differentiate the unexpanded static spherical metric and full trace on
the actual RN reference, with TT fields reconstructed from their radial
integral; compare the mixed row at three radii and keep all anomaly
terms. Verify the linearized row by central differences at 1e-4 and
5e-5 with quadratic refinement, retaining the prescribed temperature.
These are diagnostics of the ordinary argument, not curved spectral
data, physical boundary realization or a coupled solution. Preserve
every run, failure, source snapshot and log row. No Lean is selected.

Remaining to-do list: full physical coupled domain and constraint-
compatible inverse; finite state, residual, admission and nonlinear
solution; exact compatibility, matching, validity, all-sector stability,
charge/return, opening, transport, preparation, separate creation and
full-result audit.
