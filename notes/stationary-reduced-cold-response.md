# Coupled cold response and its quantum product exponent

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary proof before
focused verification. **PROGRESS; full exterior inverse and nonlinear
physical control remain open.**

The actual classical/reference seed has an explicit conserved-source
inverse for radius, lapse, geometric mass and cold shape, at its fixed
temperature and laboratory clock. The complete stationary quantum
product has a coupled reduced normal matrix with a Ward factor. Its
physical cold root satisfies |nu-1|<1.3e-7, and the corresponding
finite-neck exponent factor is below 1.00002.

The last bound includes the full massive tower, LLL, common finite
subtraction and the fixed covariant reduction. It follows from a
uniform complex low-frequency value estimate, not from an isolated
scalar mass or the TT energy bound. It does not give a whole-frequency
product resolvent, the actual exterior seam inverse, a nonlinear branch
or its physical EFT error. The original signed source condition remains
undecided.

Continue the [covariant reduced equations and seed inverse](stationary-covariant-reduced-inverse.md).
The companion [parameter form bounds](stationary-reduced-parameter-form-bounds.md)
controls an actual common finite contact and supplies auxiliary outer
logarithm estimates; the full nonlinear response remains to be proved.

## 1. Fixed inputs and normalizations

The satellite notes used below are [stationary-centered-covariant-seed](stationary-centered-covariant-seed.md),
[stationary-centered-first-order-geometry](stationary-centered-first-order-geometry.md), [stationary-centered-exterior-source](stationary-centered-exterior-source.md),
[stationary-centered-complete-residual](stationary-centered-complete-residual.md), [stationary-reference-cold-bolt](stationary-reference-cold-bolt.md),
[stationary-seed-quantum-form](stationary-seed-quantum-form.md) and [stationary-centered-canonical-response](stationary-centered-canonical-response.md).
All are in the frozen quantyra-wormhole-semiclassical repository. Their
actual centered radius, temperature, state, common angular heat subtraction,
scale and complementary paired domain are retained.

Write lambda=lambda_e, x=lambda B^2, C=B^2/A^2=1+6x+d,
R=B(1+h), gamma=kappa A^2/B, epsilon=G/(60pi), c_R=1/(5760pi^2).
The accepted benchmark gives

    x>7.9e9, |d|<1e-12, h<1e-30, log(1/h)<84,
    |h/gamma-1|<=12gamma, G C/B^2<8.5e-33,
    M_n^2=n(n+q)/C, d_n=2(q+2n), q=2000^4, M_1>18,
    N(E)=sum_(M_n<=E)d_n<=4CE^2.

The weaker log bound 84 is sufficient throughout. It follows from the
actual preparation transfer in [stationary-scalar-outgoing-inverse](stationary-scalar-outgoing-inverse.md); no
extremal limit or temperature selection is taken. One physical exterior
uses a rank-two orbit Dirac kernel. The d_n already count the two nonzero
sphere-Dirac signs. A folded rank-four description must be divided by two.

The reduced action is I_red=I_full-c_R integral(R4+U(r))^2 dvol4, where

    U(r)=12lambda-(d/B^2)(4u^5-6u^6), u=B/r.

The homogeneous product H2_A x S2_B is stationary for I0 and I_full, and
R4+U=0 there. It is consequently stationary for I_red too. The finite
thermal analytic seed is stationary for I0, but not generally for I_full
or I_red. This distinction is essential in the Ward and contact arguments.

## 2. The actual classical/reference observable inverse

Let rho,p be any smooth spherical source in the same sphere-integrated
convention as [stationary-centered-first-order-geometry](stationary-centered-first-order-geometry.md)(2). Include its
angular partner so the joint Ward identity holds, require HH regularity
rho(R)+p(R)=0, and assume the integrals and norms below are finite. This
is a statement about the I0 constraint inverse; it is not a declaration
that a raw off-shell quantum variation is separately conserved.

Define the source numbers

    sigma_h=|G rho(R)/(kappa R)|,
    L_s=G integral_R^infinity |rho+p|/(r F0) dr,
    Q_rho=(G B^2/x) sup_(r>=R) |rho'(r)|/r,
    K_h=1.001(sigma_h+L_s).

The exact fixed-clock, fixed-temperature solution is

    a(r)=-G integral_r^infinity (rho+p)/(t F0(t)) dt,
    D_h=F0''(R)/2+kappa/R,
    s_R=[G rho(R)/R-kappa a(R)]/D_h,
    f(r)=-[2kappa R s_R+2G integral_R^r rho(t)dt]/r.       (1)

These enforce a(infinity)=0, f(R)+2kappa s_R=0 and
f'(R)+F0''(R)s_R+2kappa a(R)=0. In particular they fix the physical
surface gravity on the laboratory clock, not a subsequently rescaled
coordinate period. The accepted seed guards |D_h A^2-1|<7h and
|h/gamma-1|<=12gamma give kappa/(D_h B h)<1.001. Therefore

    |s_R|/(R-B)<=K_h,             ||a||_infinity<=L_s.    (2)

This is the required thermal-gap cancellation: a horizon source measured
by sigma_h is not divided by an additional h.

Put fbar=f+s_R F0'. The exact twice-subtracted numerator in the accepted
geometry proof gives, for this same arbitrary source,

    ||fbar/F0||<=L_s+3Q_rho+80h K_h.                     (3)

For clarity, that proof uses only |2F0''+rF0'''|<10C/B^2,
F0>=2kappa(r-R)+.95x(r-R)^2/B^2 and the definition of Q_rho;
it does not use a special sign or spectral property of the Dirac source.

In the fixed common disk, let xi=delta r and u_c=delta sigma_hat. Then

    H=a+fbar/F0, J(r)=integral_r^infinity H/F0,
    xi=s_R-F0 J,
    u_c=a+fbar/(2F0)-s_R/r-(F0'/2-F0/r)J.                (4)

Fixed temperature gives H(R)=0, so J is regular at the bolt and
xi(R)=s_R. The actual seed inequality

    [F0/(r-B)] integral_r^infinity dt/F0(t)<10

now proves the explicit global radius bound

    ||xi/(r-B)||<=K_V,
    K_V=1.002 sigma_h+21.002 L_s+30Q_rho.                (5)

Indeed the unrounded expression is
(1+800h)K_h+20L_s+30Q_rho, which is smaller than (5). This is
an actual-seed bound, not a product approximation or a scalar-only
estimate. Radius and orbit have already been solved together.

At the physical cold seam s=r/B-1=.01, the accepted direct integrals
give F0 integral_r^infinity dt/F0<=1.75(r-B) and
|(F0'/2F0-1/r)(r-B)|<1.01. Substitution in (4) yields

    |u_c|_seam<=U_c,
    U_c=5.1L_s+6.9Q_rho+183h K_h.                      (6)

The unrounded coefficients are 5.035, 6.8025 and 182.4. Thus (6)
does not introduce an independent freely chosen seam amplitude.

The correlated geometric mass is also an output of (1):

    delta mu_infinity=kappa R s_R+G integral_R^infinity rho dr,
    delta mu_h,moving=[mu_seed'(R)+kappa R]s_R.           (7)

The integral must be controlled if an absolute mass allowance is wanted.
The second expression retains the derivative of the reference geometric
mass at the moving horizon. No independent cubic lapse or mass datum is
added to the prepared-temperature solve.

## 3. The actual cold-shape barrier

Let A_s denote the four-dimensional trace source in the convention
delta R4=-8piG A_s. On h<=s<=s_c=.01 define

    P_tr=sup |8piG r^2 A_s|/(C s),
    P_st=sup G(|rho|+|p|)/(C s).

The exact cold curvature identity, including its radius partner, is

    delta Rbar=-8piG r^2 A_s-6f+6G(p-rho),
    delta Rbar|disk=delta Rbar|areal+xi Rbar0'.           (8)

Here Rbar=B^2 Rhat and |B Rbar0'|/C<=8.1. Formula (1) also gives

    sup_cold |f|/(C s)<=s_c P_st+3h K_h.

To check this last inequality, integrate |G rho|<=C P_st s from h
to s. Its contribution to |f|/C is at most P_st s^2. The horizon
term is below 3h^2 K_h because kappa B/(Ch)<7/6 and R/r<=1.
Consequently put

    R_c=P_tr+6.06P_st+18h K_h+8.1K_V.

Then A^2 |delta Rhat|disk<=R_c s. This includes the Lie displacement
term; bounding the areal curvature alone would miss the largest term.

In the accepted old cold chart, ghat=A^2 e^(2sigma_cold) g_H2,
|sigma_cold|<.05, .992h cosh chi<=s<=1.008h cosh chi and

    (-Delta_H2+V_c)u_c=(A^2/2)e^(2sigma_cold)delta Rhat,
    V_c=-A^2 e^(2sigma_cold)Rhat0>=2-20s>0.             (9)

The fixed dilation between this local chart and the full common disk
is kept. In particular u_c is the same scalar under that transition;
the auxiliary H2 cap is not substituted for the full exterior.

Let L_c be the actual chi coordinate of the physical seam s=.01.
The preceding bounds give 60<L_c<84. Define

    W=h[(L_c+1)cosh chi-chi sinh chi].

Then W>0, W(L_c)>=h cosh L_c>=.01/1.008, and
(-Delta+2)W=h[3cosh chi+chi/sinh chi]. The same elementary bound
used in the accepted cold proof gives 20sW/(h cosh chi)<.3.
One direct verification writes t=L_c-chi and
W/h=(t+1)cosh chi+chi exp(-chi). Since (t+1)e^-t<=1,
h cosh L_c<=.01/.992 and h(L_c+1)<85e-30, one obtains
sW/(h cosh chi)<=1.008[.01/.992+2h(L_c+1)]<.0102.
This is stronger than the .015 used for the .3 reserve. Hence
(-Delta+V_c)W>=2.7h cosh chi.

The right side of (9) is bounded by .56 R_c h cosh chi. The maximum
principle, including the actual seam datum (6), therefore proves

    |u_c|<=K_c W,
    K_c=max(102 U_c,.22 R_c),
    |u_c(0)|/h<=85K_c.                                 (10)

The deliberately rounded 102 covers the distinction between the actual
seam and arccosh(.01/h). No source estimate is extended across a false
identification of these two surfaces. Bolt smoothness gives u_c'(0)=0;
the radial equation gives its second bolt jet in the same h ideal.

For the retained first source, the accepted inputs are
sigma_h<5.1e-12, L_s<3.3e-10, Q_rho<1.3e-10,
P_tr<5.1e-11 and P_st<5.72e-11. Equations (2),(5),(6),(10) give

    K_h<3.36e-10, K_V<1.09e-8, U_c<2.59e-9,
    R_c<8.9e-8, |u_c(0)|/h<2.7e-5.                     (11)

These are ordinary consequences of the retained source bounds; no new
source or simulation was evaluated. The previous larger published bounds
remain valid. Equations (1)-(10), rather than only these example numbers,
are the useful classical observable inverse for a future defect source.

## 4. The exact classical product coupling

At the homogeneous product write X=r^2=B^2 e^(2S) and
g2=A^2 e^(2omega)g_H2. The disk shape variable is T=omega-S. Define

    P2=-Delta_H2+2,
    beta=A^2 B^2 V0''(B^2)=-(4+12x+10d)/C.

In particular |2-beta|<5. Expanding the covariant spherical I0 integrand
X R2+(grad X)^2/(2X)+V0(X) through quadratic order gives, in unit H2 area,

    B^2[2|grad S|^2+4grad S.grad omega+2beta S^2+8S omega].

Here V0(B^2)=0 and V0'(B^2)=2/A^2 cancel the linear terms. After the
Euclidean prefactor and physical four-volume are included, the Hessian
in (S,T) is exactly

    -(4piG A^2)^-1 [[3P2+beta-2, P2], [P2, 0]].         (12)

For a source normalized by this common prefactor, (12) first solves
P2 S=f_r and then

    P2 T=f_s-3f_r+(2-beta)S.                            (13)

Thus the product radius and shape use P2 successively. The isolated
conformal operator with mass 4lambda is not the radius inverse.

There is an exact finite-neck estimate. Let p0(chi)=cosh chi, take regular
bolt data and u(L)=0, and assume |f|<=a h p0 in P2 u=f. Writing u=p0 v
gives

    (p0^2 sinh chi v')'=-p0 sinh chi f,
    integral_0^chi p0^2 sinh t dt=(p0(chi)^3-1)/3.

The positive inverse and direct integration imply

    |u|/p0<=a h B_L/3,
    B_L=1-sech L+log[(cosh L+1)/2]<=L+1.                 (14)

A nonzero endpoint b adds exactly b p0/p0(L). The two solves (13)
therefore cost at most two factors (L+1)/3 in the h cosh chi norm.
For example, if normalized forcing bounds are a_r,a_s and endpoint
ratios b_r=|S(L)|/[h p0(L)], b_s=|T(L)|/[h p0(L)], then

    N_S<=b_r+q_L a_r,
    N_T<=b_s+q_L[a_s+3a_r+5N_S], q_L=(L+1)/3<85/3.     (15)

No power h^(Delta-2) occurs in this coupled classical estimate. For zero
endpoint data its derivative bounds follow from the same divergence
identity: |u'|<=a h(2L+5)sinh chi/6 and
|u''|<=a h(L+5/2)cosh chi. This is still a classical product statement;
the next sections retain the actual reduced quantum product response.

## 5. The exact reduced product matrix

Use the physical scalar pairing integral dvol4; after sphere integration
dvol4=4pi B^2 A^2 dvol_H2. Set p=P2 in invariant spectral notation,
gamma0=(8piG A^2)^-1, and define the real coefficients

    C_orb=6+24epsilon lambda-22epsilon/(3B^2)-28epsilon d/B^2,
    b=1-8epsilon lambda-epsilon(11/3+4d)/B^2,
    V_U=4lambda-2d/B^2,
    P_S=24lambda-4d/B^2+88epsilon lambda^2+6epsilon V_U^2,
    a(p)=C_orb(p-2)+A^2 P_S.                            (16)

The full reduced product Hessian, including all quantum shape response,
has the exact form

    H_red(p)=gamma0 [[-a(p), -2b p], [-2b p, q_cold(p)]]. (17)

To verify the local row without approximating the quantum determinant,
use the exact reduced trace

    E_red=R+12lambda0+epsilon[(3/2)C_W^2-(11/12)Euler4-20M_em
                         -Box U+(2U+rU')(R+U)/6].

At the product R=-U=-12lambda-2d/B^2, rU'=-16d/B^2,
J_orbit=-(A^-2+B^-2)/2, R_sphere=2/B^2. The conformal principal
tensor in the reduced-equation theorem gives C_orb. A constant conformal
variation has potential -2U+48lambda0+epsilon(2U+rU')^2/6;
using lambda0=lambda+11epsilon lambda^2/6 gives P_S in (16).
This step uses the exact stationary PRODUCT trace equation. It must not
be replaced by the classical I0 potential or applied at the off-shell
finite thermal seed without its density contact.

For an orbit-only T, delta R=delta R2=(2/A^2)P2 T,
delta C_W^2=(2R/3)delta R2 and delta Euler4=2R_sphere delta R2.
Consequently delta E_red=b delta R2, proving the mixed entry. Symmetry
of the actual stationary action supplies the transpose with the same
pairing. To recover the fixed-reference I0 system, use C_orb=6, b=1,
P_S=H_U=24lambda-12d/B^2 and q_cold=0. Then

    a_I0(p)=6(p-2)+A^2 H_U=2(3p+beta-2),

so (17) gives (12) exactly. Merely setting epsilon=0 in the FULL
product trace formula while holding the quantum-induced reference d
fixed is a different operation: it leaves an extra 8d/B^2 in P_S.
The product stationarity relation used in deriving P_S must be varied
consistently with the theory, rather than silently kept in that limit.

## 6. The product Ward factor and its domain

On H2, Hess(cosh chi)=cosh chi g. Thus T=cosh chi, S=0 is the Lie
derivative generated by grad(cosh chi). Covariance and the zero full
PRODUCT equation imply that the invariant Hessian annihilates this
direction. For the massive nonlocal part this statement can be obtained
by cutting off the vector field, applying the compact Ward identity,
and then sending the cutoff outwards. The discarded annular variations
grow at most exponentially, while every massive Green product decays
at least exp(-M_n d), M_n>18. The angular tail is absolutely summable
at positive separation. The LLL conformal response and retained contacts
are local. Thus the cutoff limit is legitimate for the invariant kernel.

Equivalently, let nu(p)=(-1+sqrt(9-4p))/2 near p=0 and
phi_p(d)=P_nu(p)(cosh d). Then Delta phi_p=(2-p)phi_p,
phi_p(0)=1 and phi_0=cosh d. The complete product multiplier is analytic
near p=0 and

    q_cold(0)=0, q_cold(p)=p g_cold(p).                  (18)

The square root uses the branch nu(0)=1. Its analyticity here follows
from the massive kernel gap and the same finite local subtraction; the
massless LLL contribution is the explicit local Polyakov polynomial.

The growing field cosh chi is not asserted to satisfy the actual finite
exterior boundary condition, fixed seam, or fixed clock. The analytic
seed has nonconstant radius, and an orbit diffeomorphism then also moves
that radius. Accordingly (18) is not an actual global kernel vector or
a new allowed mass/temperature parameter. It is the correct product
normal-factor constraint on the local coupled calculation.

Eliminating T from (17), with scaled sources F_S,F_T, gives exactly

    D_r(p)=4b^2 p+a(p)g_cold(p),
    S=-[g_cold(p) F_S+2b F_T]/D_r(p),
    det(H_red/gamma0)=-p D_r(p).                        (19)

In the I0 system just specified, g_cold=0 and D_r=4p, and (13) is
recovered. A small
quantum g_cold(0) moves the physical radius root linearly,
p_root approximately -a(0)g_cold(0)/(4b^2), with a(0) approximately -8.
The exact factor p in the determinant is retained; a generic square-root
splitting of a double root is not the relevant model here.

If nu=1+delta_nu is the physical root, propagation over this finite neck
changes a classical cosh profile by at most a constant times
exp(|delta_nu|L_c). This observation explains why an explicit small
product response bound would be useful. It is not a bound on that
response or on the actual coupled exterior inverse.

## 7. The convergent Dirac-bubble calculation

Write the invariant rank-two Green kernel on unit H2 using its positive
coefficients A_M(d),B_M(d), normalized by B_M~1/(2pi d). They obey

    A_M'+tanh(d/2)A_M/2=-M B_M,
    B_M'+coth(d/2)B_M/2=-M A_M,
    A_M/B_M ->1 at infinity.                            (20)

These are the actual complete-H2 Green equations already justified from
Muck's invariant construction in [stationary-reference-cold-bolt](stationary-reference-cold-bolt.md)(15),
[Muck, invariant spinor Green function](https://arxiv.org/html/hep-th/9912059) . Use new names for the bubble:

    K_M(d)=2M^2[A_M(d)^2-B_M(d)^2].

The spin parallel transports cancel in G(x,y)G(y,x), so this is exactly
M^2 tr[G(x,y)G(y,x)]. Its flat limit agrees with the rank-two Bessel
bubble in [stationary-seed-quantum-form](stationary-seed-quantum-form.md). The one-vertex contact is
retained through the stationary-product p=0 subtraction below.

For one massive mode define

    Pi_M(p)=2pi integral_0^infinity sinh d K_M(d)
                          [phi_p(d)-cosh d] dd-p/(12pi),
    R_M(p)=Pi_M(p)+p^2/(120pi M^2).                     (21)

The last term of Pi_M is the specified 2D conformal anomaly; its sign
agrees with Gamma_LLL''=-q integral|grad T|^2/(12pi).
The vacuum/reference contact is fixed by product stationarity, which
is why the test difference is phi_p-phi_0 rather than phi_p-1. No
mode-dependent physical finite coefficient has been introduced.

The common a2 restoration is exact on this product. Before the
stationary vacuum contact, its orbit Hessian is
alpha[8Delta^2-32Delta+16], alpha=-C F_eff/(960pi).
Subtract its value on phi_0: the resulting polynomial is
8alpha p^2=-C F_eff p^2/(120pi), where

    F_eff=4gamma_E-2H_q+4log(mu0 B).

The per-mode a2 term is correspondingly -p^2/(120pi M^2).
The accepted local normal-cap expansion after that term has a summable
inverse-cubic (or better) massive remainder on a fixed smooth bounded
test. In (21) the test and its finite jets are uniformly bounded on
each fixed ball for p in a small complex neighborhood of zero. Outside
that ball the massive gap controls its at-most-exponential growth.
Thus sum d_n R_Mn(p) converges locally uniformly there. This argument
proves local analytic continuation; it supplies no numerical norm on
the whole weighted frequency contour required for a seam inverse.

At r=B the Weyl field Phi=log(r/B) is zero for every orbit-only T,
so W4 has zero pure-T Hessian. The reference potential supplies the
stationary zero contact already used. The exact reduction term has
delta e0=(2/A^2)P2 T, hence subtracts
32pi c_R C p^2=C p^2/(180pi) in unit-H2 orbit area. Therefore the
COMPLETE reduced product shape response is

    Qhat_red(p)=-q p/(12pi)
               -C[F_eff/(120pi)+1/(180pi)]p^2
               +sum_(n>=1)d_n R_Mn(p),
    q_cold(p)=(2G/B^2)Qhat_red(p).                      (22)

The factor 2G/B^2 follows from dvol4=4pi B^2 A^2 dvol_H2 and
gamma0=(8piG A^2)^-1. It is not G/A^2 without the corresponding
C in the orbit response. All common-scheme fourth contacts and the
exact reduction subtraction have been retained in (22).

In particular the coefficient which moves the physical cold root is
the following convergent scalar moment:

    g_cold(0)=(2G/B^2){-q/(12pi)+sum_n d_n [
      -(4pi M_n^2/3) integral_0^infinity sinh d
         (A_Mn^2-B_Mn^2)
         [cosh d log((cosh d+1)/2)+cosh d-1] dd
      -1/(12pi)]}.                                    (23)

Neither the a2 restoration nor the reduction p^2 term contributes to
(23). This is more specific than an unnamed missing cold constant:
the known positive solutions of two first-order ODEs, one radial moment,
and the declared angular spectrum determine it. Section 9 bounds its
low-frequency continuation without evaluating its signed value.

The weight in (23) needs no special-function derivative assumption.
If W(x)=x log((x+1)/2)+x-1, then
[(x^2-1)W']'-2W=3x and W(1)=0 with the regular endpoint jet.
Differentiating the Legendre equation at nu=1 therefore gives
partial_p phi_p|0=-W(cosh d)/3=:psi(d). In particular

    psi(0)=0, Delta psi(0)=-1, Delta^2 psi(0)=-4.

The a2 operator at the center acts by 8(-4)-32(-1)=0. This is the
explicit cancellation which makes the massive slope summable after
the anomaly cancellation; a bound on the unrenormalized bubble alone
would miss it.

## 8. An explicit spatial tail

The accepted H2 bound ||G_M(d)||<=2e^(-Md/2)/(pi d) implies
|K_M(d)|<=8M^2 e^(-Md)/(pi^2 d^2). Exact counting gives
sum d_n M_n^2 e^(-M_n d)<=72C/d^4. For a direct count proof,
Stieltjes integration by parts and discarding its negative part bound
the sum by (4C/d^4) integral_2^infinity (v^4-2v^3)e^-v dv
=368e^-2 C/d^4<72C/d^4. Extracting half the gap gives

    sum d_n |K_Mn(d)|<=1024 C d^-6 e^(-9d).

Also sinh d<=d e^d and
|psi(d)|<=d^2 cosh d/4<=d^2 e^d/4. For the latter use
2log cosh(d/2)<=d^2/4 and cosh d-1<=d^2 cosh d/2.
Since 2pi<7, for any a>0 the slope's bubble tail satisfies

    2pi sum d_n integral_a^infinity sinh d |K_Mn psi| dd
         <=1792 C integral_a^infinity d^-3 e^(-7d) dd
         <=896 C/a^2.

For a=1/64 this is below 4e6 C. Its dimensionless contribution to
|g_cold(0)| is below 8e6 G C/B^2<6.8e-26. The unresolved contribution
is consequently a LOCAL renormalized bubble moment, not an uncontrolled
long-distance integral or the massless LLL. The LLL contribution itself
has magnitude Gq/(6piB^2)<1.9e-31 using q/C<400.

This tail bound alone does not control the local renormalized moment.
The next section supplies that separate estimate using a single compact
C5 probe and the frozen VALUE-derivative ledger. No TT form bound is
promoted to a bound on a growing Legendre test.

## 9. A uniform low-frequency product bound

The following argument proves the loose target required for the cold
root, on the whole complex disk |p|<=.01. It is a product-only pointwise
response calculation. The additional global contour and actual seam
estimates in section 10 are not consequences of it.

Use the exact H2 conformal normal chart

    g_H2=(1-|z|^2/4)^-2 |dz|^2,
    t=(z1^2+z2^2)/4, cosh d=(1+t)/(1-t).

Here complex z1,z2 are used only to bound the explicit product probe,
not to spatially continue an unknown C5 candidate. The regular radial
Legendre solution F_lambda(t), lambda=2-p, satisfies

    (t F_lambda')'=lambda(1-t)^-2 F_lambda, F_lambda(0)=1,
    F_lambda=1+lambda K F_lambda,
    (KF)(t)=integral_0^t log(t/s)(1-s)^-2 F(s) ds.

On the analytic sup norm |t|<=.1, ||K||<=.1/.9^2 and
|lambda|<=2.01. Since F_2=(1+t)/(1-t), the regular quotient
Psi_p=(F_(2-p)-F_2)/p obeys

    Psi_p=-K F_2+(2-p)K Psi_p,
    ||Psi_p||<= [(.1/.9^2)(1.1/.9)]/[1-2.01(.1/.9^2)]<.21.

The formula extends through p=0 by analyticity. For real |z|<=1/8,
a complex polydisk of radius1/8 in each Cartesian coordinate remains
inside |t|<.023<.1. Multivariable Cauchy and the Euclidean tensor norm
therefore give |D^j Psi_p|<=.21 j!12^j for j<=5.

Choose a fixed C5 cutoff equal to one for |z|<=1/16 and zero for
|z|>=1/8. In its annulus set

    v=(|z|^2-1/256)/(3/256), chi=1-I_v(6,6).

The beta polynomial has degree11 and its first five derivatives vanish
at the joins. Its Bernstein derivative bound is
|partial_v^k chi|<=2^k(11)_k, where (11)_k is falling factorial.
Also |Dv|<=64/3, |D^2v|<=512/3 and higher derivatives vanish.
Consequently the exact multilinear Bell bound for
|D^j chi|/(j!100^j) is the coefficient of z^j in

    [1+(32/75)z+(32/1875)z^2]^11.

For j=0,...,5 these coefficients are bounded respectively by
(1,5,11,14,13,9). The product probe U_p=chi Psi_p therefore satisfies

    U_p(0)=DU_p(0)=0,
    |D^j U_p|/(j!100^j)
         <=.21 sum_(k=0)^j .12^k B_(j-k)<4, j<=5.       (24)

Its global value/first-derivative bounds also satisfy the frozen C1
coefficient contract (even the coarse bounds 4 and400 are below the
corresponding 1e4 times4 graph reserve). The cutoff is supported well
inside the unchanged normal cap. At the observation point the endpoint
normalizing map is exactly the identity for the entire family
sigma_tau=sigma_H2+tau Re(U_p), or tau Im(U_p), because both the value
and first jet of U_p vanish there. No moving-normalizer derivative or
unknown sixth metric jet is introduced.

Apply the local calculation in [stationary-centered-canonical-response](stationary-centered-canonical-response.md),
sections 'Real quadratic normal gauge and the cold local loops' and
'Cold finite gluing and its mixed response', to these TWO real
infinitesimal families. This is a direct application of that local
VALUE lemma, with effective direction coefficient4 for each family.
One can first take tau arbitrarily small so every endpoint smallness
condition holds, apply the finite endpoint inequality and divide by
|tau|. The local lemma is homogeneous in this direction coefficient.
It does not require a finite-amplitude direction of size four.

All the local metric and gluing hypotheses can be checked here without
an actual-seed transfer. On |z|<=1/4 the background normal mass potential
is |exp(sigma_H2)-1|<=1/63<.026. The actual global comparison operator is
the complete H2 massive Dirac operator, with its already proved inverse
and weighted kernel bounds. The local flat cap is the same radius1/4
normal cap as in the frozen proof; its actual and capped families agree
on the support of the gluing cutoff. The probe is supported inside1/8.
The fixed-center coordinate-transport terms are zero. There is no
distant actual-to-H2 seam term in this product calculation. The flat
cap itself is still distinct from H2, and its finite two-sided gluing
term is retained with the frozen row and graph constants.

After the SAME a0/a1/a2 extraction, that calculation gives per real
family the trace-value remainder bound 4e22 M^-3. Splitting real and
imaginary parts thus gives 8e22 M^-3 for U_p. The a0/density term is
zero at the center. The massive a1 term cancels the ordinary 2D anomaly
before summing. The common a2 restoration is retained; it is not set
to zero for p different from zero. Its action on Psi_p at the center
is exactly 8alpha p, since Delta Psi_p(0)=-1 and
Delta^2 Psi_p(0)=-4+p. The frozen restoration bound 8e13 C in orbit
units is a conservative allowance for this same explicit term.

To check the trace/Hessian conversion, write the orbit effective-action
variation as delta Gamma_hat=-integral T tau_orbit dvol_H2 and
Theta=C tau_orbit. At a center where U_p=0 the varied density term is
zero. The Hessian applied to U_p at that center is therefore
- delta Theta/C. For the massive resolvent it is exactly
M^2 tr[G U_p G] there; the one-vertex seagull also has factor U_p(0)
and vanishes. The spin rank remains two and d_n already has its two
sphere signs. There is no extra normalization factor in this conversion.

The angular estimate sum d_n/M_n^3<2C/3 gives the total local plus
finite-cap-gluing remainder below (16/3)e22 C. The common restoration
allowance just stated is much smaller. The explicit LLL response adds
at most q/(12pi)<12C; the reduction term adds at most C on |p|<=.01.
Reference and other zero contacts act by U_p(0)=0.

Finally replace the compact probe by the true Psi_p. Its omitted part
is supported at proper distance at least1/16, hence at least1/64. The
uniform complex tail has the same reserve as section 8. For completeness,
the positive radial Volterra equation gives

    |Psi_p(d)|<=partial_lambda phi_lambda(d)|_(lambda=2.01)
               <=d^2 exp(1.01d)/4.

To verify the last comparison, choose nu_+>0 with
nu_+(nu_++1)=2.01, so nu_+<1.01. The regular solution phi_2.01 is
bounded by exp(nu_+ d), and
(Delta-2.01)[d^2 exp(nu_+ d)/4]>=exp(nu_+ d).
The positive radial Volterra inverse, with zero regular initial data,
then bounds the lambda derivative. This also bounds the complex divided
difference, since the Taylor coefficients of the Volterra series in
lambda are positive. In section 8 replace exp(-7d) by exp(-6.99d);
its bound by 896C/a^2 is unchanged. Thus the omitted bubble tail is
below4e6 C. All local contacts vanish on that omitted direction.

Combining the local, common-contact, finite-gluing and direct-tail terms,
with a further generous reserve, proves

    sup_(|p|<=.01) |Qhat_red(p)/p|<1e24 C<1e25 C.         (25)

This is a finite explicit low-frequency enclosure, not an evaluation of
the signed slope (23), a whole-frequency coercivity theorem, or a
same-C5 operator inverse. In particular the use of a smooth explicit
probe in a VALUE-response lemma does not remove the full shape
operator's fourth/logarithmic order.

## 10. Consequence and remaining scope

The actual classical source-to-observable map (1)-(10) is now explicit.
It retains the coupled radius/shape equations, fixed temperature and
clock, the mass output, the moving-bolt pullback and actual seam data.
For a full reduced linear solve, it can be applied to the effective
classical forcing after all quantum, density and Ward contacts have
been placed on the right side. That forcing is not known from the
published unsigned seed source bounds alone.

The next product calculation is a
resolvent bound for D_r(p)=4b^2p+a(p)g_cold(p) on the weighted H2
frequency contour. A small value of g_cold(0) controls a nearby root;
it does not by itself control the whole inverse. The full shape block
retains its fourth/logarithmic high-frequency order and its additional
regularity requirements.

The proved low-frequency estimate (25), even weakened to its1e25 C
reserve, gives |g_cold|<1.7e-7 on that disk. The explicit coefficients
(16) have |a|<9 and |b|>.999 there. Rouche on |p|=1e-6 gives exactly
one zero of D_r inside, and its own equation gives |p_root|<3.84e-7.
Reality of the multiplier makes this unique root real. Since
p=-(nu-1)(nu+2), it satisfies

    |nu-1|<1.3e-7, exp(|nu-1| L_c)<1.00002.             (26)

This explicitly closes the proposed exponent-loss issue for the actual
quantum REDUCED PRODUCT normal matrix. It is not a complete weighted
inverse or an actual-seed seam bound. The stricter1e24 reserve in (25)
would improve (26), but that improvement is not needed here.

The actual-seed transfer then needs the variable-radius cold perturbation
of this coupled product resolvent and the complete matching trace map
through the cold seam to the actual middle/outer problem. All required
shape jets/nonlocal exterior terms must be retained; treating this as
two freely prescribed second-order Dirichlet problems would be wrong.
The outer mass and clock constraints are those in (1),(7), not new
parameters used to tune a cold homogeneous solution. The accepted
critical energy inverse constant 1000/989 does not supply this weighted
trace-map norm.

There is also no same-C5 contraction hidden in (1): the classical
constraint inverse gains two derivatives, while the actual quantum
shape response has fourth/log order. One must use the reduced mixed
linear inverse and a proved parameter-dependent nonlinear map. The
explicit source-to-observable map becomes useful once the genuine
defect is bounded in its admissible Ward-completed source class.

Thus the scalar-only throat warning is resolved at the classical coupled
level and sharpened at the exact quantum product level. This is not a
nonlinear exterior, a physical reduction-error theorem, or an evaluation
of the original full-action preparation moment. The original source
moment, exterior/state/charge matching and E003/S074-S080 physical gates
remain open.

The shared [focused checker](../checks/stationary_reduced_cold_response.py)
audits finite identities and rational reserves in this theorem and its
parameter-form companion. Its [receipt](../evidence/stationary-reduced-cold-response.json)
and [raw log](../evidence/stationary-reduced-cold-response.log) retain the
actual execution outcome and input hashes. Operator domains, Ward limits,
source-lemma application and analytic continuation remain ordinary proofs;
finite algebra checks do not replace them.

Remaining to-do list: bound the full weighted product resolvent from
(21)-(23); transfer it through the actual coupled seam;
prove the nonlinear mixed map and physical remainder/observable bounds;
separately determine the original signed source moment and complete the
retained exterior/state/charge and S074-S080 requirements.
