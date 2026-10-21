# The normal quantum shape action at a moving cutoff

2026-10-21, E003/S168. Ordinary argument before focused verification.
The complete stable-wormhole goal remains incomplete.

The complete FLAT GEODESIC angular-tower model, combined with its
leading local Weyl logarithm and boundary contact, has a finite weak
cutoff limit on the compensated logarithmic form domain. Its nonzero
normal derivative at a finite cut must be retained. The principal
contact is twice the first times the second derivative; discarding it
breaks the scale cancellation. A sharper tower remainder estimate,
of order p^2 log p instead of the previous p^3 allowance, makes the
remaining model multiplier continuous on these cut profiles.

This resolves the moving-cut problem for the specified normal model.
It does NOT identify the original charged, thermal Dirac Hessian with
that model. In particular the actual potential/second-vertex contact,
curved boundary kernels, state and remote geometry remain in the
actual-minus-model difference. No physical fast branch is excluded.

## 1. The finite-cut principal representative

Use inward normal distance y>0 and a cut y=d. Initially let v be real,
smooth at every positive y, with fixed compact support at the far end.
Set h=v_yy, a=v_y(d), b=h(d), and w_d(x)=v(d+|x|). Fourier transforms
on R are unitary with phase exp(-ipx). Define the homogeneous principal
representative by meromorphic continuation,

    R_d(s)=(1/2) integral_R |p|^(4-2s) |w_d_hat(p)|^2 dp,
    L_d=-(1/2) R_d'(0),   J_d=R_d(0).                 (1)

The starting strip is 1/2<Re(s)<5/2. Far-end support can instead have
rapid decay. This is a normal-model representative, not an assertion
that an arbitrary nonsmooth double has the original determinant.
The logarithmic multiplier with scale mu is L_d-log(mu)J_d.

Let H_d(x)=h(d+|x|). Distributionally w_d''=H_d+2a delta. Continuing
the Fourier pairing, with the scaleless delta-square integral zero,
gives

    J_d=integral_d^infinity h^2 dy+2ab,
    L_d=(1/2)<H_d,log|D| H_d>-2a T[h(d+.)],
    T[f]=integral_0^infinity [f(x)-f(0)1_(x<ell)] dx/x
                                      +(gamma_E+log ell)f(0). (2)

ell is an arbitrary fixed positive subtraction length; T is independent
of it. To verify the continuation rather than multiply distributions,
the cross term uses the Riesz kernel

    c_s=Gamma(1-2s)sin(pi s)/pi=s+2gamma_E s^2+O(s^3),
    I_s[f]=2c_s integral_0^infinity x^(2s-1)f(x)dx,
    I_0[f]=f(0), I_0'[f]=2T[f].                       (3)

Subtract its constant Taylor term before continuation. Splitting the
delta-square momentum integral at an arbitrary positive scale gives
opposite analytic pieces, whose sum is zero. Equations (2) follow
without treating the cusp as an ordinary H2 full-line function.
For example w=exp(-k|x|) gives J=-3k^3/2 and
L=-3k^3 log(k)/2-k^3/2. The negative J is a boundary finite part,
not a negative compact-core energy or a physical instability.

## 2. The boundary logarithm is fixed, not adjustable

The coefficient of the model principal symbol is c=C/(30pi), as in
[the retained full flat tower](outer-flat-tower-normal-remainder.md).
The same leading contact follows independently from the actual mixed
heat coefficient in [Vassilevich, (5.33)](https://arxiv.org/html/hep-th/0306138v3)
and the [retained boundary action](stationary-exterior-boundary-weyl.md).
Here is a sign and normalization check before using it.

In outward coordinate x, take canonical orbit metric
exp(2T)(dx^2+dtau^2), with a fixed sphere. Retain only derivative
weight four in the bulk and three at the boundary. In units
P_f times sphere area, the bulk density is -3 Box R-(3/2)R2^2.
Its integrated Box R flux has quadratic boundary term
-12 T T'''-12 T'T''. The mixed boundary coefficient has
12 T T'''+24 T'T''. Their sum is 12 T'T'', while the remaining
bulk is -6 integral (T'')^2. Switching to inward y therefore gives
-6[integral h^2+2ab]. With sphere area and rank restored, the
quadratic coefficient of a4 is -C J_d/(60pi), and its Hessian
is -c J_d. Lower derivative weights are not set to zero.

The Euclidean local form is minus the Lorentzian boundary action
used in that retained note. Its explicit sigma term consequently
has Hessian -2c sigma(d)ab, accompanying -c integral sigma h^2.
For sigma=-log y plus its fixed constant, the compensated principal
normal form, in dimensionless units with the fixed scale absorbed,
is

    P_d=L_d+integral_d^infinity log(y)h^2 dy+2log(d)ab. (4)

A constant in the logarithm simply adds the corresponding fixed J_d.
No original finite coefficient is changed. In particular under
v(y)=d^(3/2)u(y/d), L_d(v)=L_1(u)-log(d)J_1(u).
Both log(d) pieces in (4) cancel that scale term. Keeping only the
bulk subtraction would leave the spurious -2log(d)ab.

## 3. A complete weak-limit argument for (4)

Assume v(0)=v_y(0)=0, h in L2(dy), and put

    z=-log y, Z=-log d,
    g(z)=exp(-z/2)h(exp(-z)),
    A(Z)=exp(Z/2)v_y(exp(-Z))
        =integral_Z^infinity exp(-(s-Z)/2)g(s)ds.

Then ||g||_2=||h||_2, A is in H1(R), and A'=A/2-g. The full
normal theorem uses the compensated logarithmic domain

    E_D(g)=(1/4) integral_R^2 k(|z-z'|)|g(z)-g(z')|^2 dz dz'<infinity,
    k(r)=1/[2sinh(r/2)].                                (5)

This is an explicit extra form-domain hypothesis, not a consequence
of H2 alone. Its identification with the actual auxiliary energy
completion remains part of the curved finite-cut transfer.

For zero extension from y>0, log|D_y|+log y becomes the Mellin
multiplier L_D with symbol Re psi(1/2+i xi). One direct derivation
applies the Riesz kernel to y^(-1/2-i xi); its coefficient is

    [Gamma(1/2-i xi)/Gamma(1/2-i xi+2s)
       +Gamma(1/2+i xi-2s)/Gamma(1/2+i xi)]/[2cos(pi s)]. (6)

Differentiating at zero yields the stated multiplier. The digamma
integral also gives <g,L_D g>=psi(1/2)||g||_2^2+E_D(g).
The underlying integral and reflection identities are established
[DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16) and
[5.5.4](https://dlmf.nist.gov/5.5.E4); no new uncertainty constant
is claimed.

Let g_Z=1_(z<Z)g and let C_H be the Carleman operator
(C_H f)(x)=integral_0^infinity f(x')/(x+x') dx'. Equation (2) gives
the exact decomposition

    P_d=<g_Z,L_D g_Z>-(1/2)<h(d+.),C_H h(d+.)>+B(Z),
    B(Z)=2 A(Z) psi(1/2+D_Z)g(Z).                       (7)

The last expression is a distributional product defined by duality.
For smooth data it is the literal finite part in (2), after the
2log(d)ab term cancels its logarithmic scale. Indeed, substituting
x=d(exp(r)-1) changes T to the one-sided kernel k(r), with multiplier
-psi(1/2+i xi). This fixes the sign and Euler constant in (7).

All three terms have the needed limit; none requires h(d) to converge.

First, in E_D(g_Z) the part with both variables below Z tends to
E_D(g) by dominated convergence. Its extra crossing contribution is

    q(Z)=(1/2) integral_(z<Z) |g(z)|^2 L(Z-z)dz,
    L(r)=integral_r^infinity k(s)ds=-log tanh(r/4),
    integral_0^infinity L(r)dr=pi^2/2.                 (8)

Thus q is in L1(dZ). It tends to zero in translated local L1,
even if it has no pointwise limit. The zeroth-order mass converges
ordinarily. This proves the weak limit of the first term.

Second, ||C_H||_(2->2)=pi: the unitary logarithmic transform turns
it into convolution with 1/[2cosh(r/2)], whose multiplier is
pi/cosh(pi xi). Right translation h(d+.) tends to h in L2(dy).
The middle term therefore has an ordinary finite limit.

Third, for any smooth compact test phi and phi_Z0(Z)=phi(Z-Z0),
|psi(1/2+i xi)|<=2+5|xi| follows from
|psi(1/2)|<2 and sum_(n>=0)(n+1/2)^(-2)=pi^2/2<5.
Moving the multiplier onto phi_Z0 A gives

    |<B,phi_Z0>| <=2||g||_2 [2||phi_Z0 A||_2
                                      +5||(phi_Z0 A)'||_2] ->0. (9)

This uses A,A' in L2 and no pointwise control of Bg. Equations
(7)-(9) establish, for every such test,

    lim_(Z0->infinity) <P_(exp(-Z)),phi_Z0>
                    =P_0 integral phi,
    P_0=integral_R [Re psi(1/2+i xi)-pi/(2cosh(pi xi))]
                                                |g_hat(xi)|^2 dxi. (10)

The lower symbol value is psi(1/2)-pi/2, agreeing with the retained
even logarithmic uncertainty constant. In particular the normal model
cannot provide a persistent nonconstant sinusoid in Z. This weak test
is a necessary condition on an ordinary cutoff limit; it is not an
action-averaging prescription.

A fixed fourth contact J_d has the same weak limit integral h^2:
its boundary part is 2A g, which is in L1(dZ). Other constant local
third-derivative contacts are also harmless in the critical profile
v=y^(3/2)V with V in H2(dz). For example v v''' is a sum of L1
products and one derivative of V V'', using
V V'''=D(V V'')-V'V''. This statement does not remove cut-dependent
logarithmic coefficients; those required the exact cancellation above.

## 4. The entire flat tower remainder is sufficiently lower order

The previous bound E(p)=K(p)-K_infty(p)=O(p^3) is insufficient for a
cusp: its Fourier tail is p^-2 and a cubic quadratic form is borderline.
Improve that bound before extending the model. Keep its fixed q>=1,C,
mu and kappa4. Let F(t), phi_t and psi_t be the real angular sums in
the retained normal-remainder note. With u=s+q/2 and v0=t-q^2/4,

    phi_t=2u/(u^2+v0), psi_t=2u/(u^2+v0)^2.

For t>=q^2, v0>=3t/4. Direct differentiation gives
integral |phi_t''|<=5/(2v0) and
integral |psi_t''|<=3/v0^2. The composite trapezoid error is at most
one eighth of these integrals. Keeping its half endpoint correction
therefore gives

    |F(t)+2log t| <=(q+5/6)/t,
    |t F'(t)+2| <=(q+4/3)/t.                           (11)

For 0<t<=q^2, the earlier bounds 4q/t+4/sqrt(t) are at most 8q/t.
Consequently both errors in (11) are <=8q/t for ALL t>0. Inserting
them into the same exact angular integral for k''(z), z=p^2, yields

    |(k-k_infty)''(z)| <=4q/(3pi z).                   (12)

The already fixed E(1), E_z'(1) have total magnitude <=20C.
Integrating (12), and using the retained small-p estimate, proves
|E(p)|<=K_q,C (1+p^2)log(e+p^2) for a finite fixed constant.
No fourth polynomial or scale is changed. This improves the growth
order; it is not a new calibrated physical benchmark.

Even extension is continuous from H2(0,infinity) to H^(5/4)(R).
For completeness, its derivative is the odd extension of v'. The
fractional order 1/4 cross-boundary seminorm is controlled by
integral |v'(x)|^2 x^(-1/2)dx, bounded by the half-line H1 norm.
The same-side seminorm is bounded by that norm as well. Translation
v(d+.) tends to v in H2, so their even extensions converge in H^(5/4).
The remainder multiplier E is continuous in this norm. Its quadratic
form thus converges ordinarily at the moving cut, without a finite
part. Along with (10) and the fixed -31/30 fourth term, this proves
the weak-limit result for the FULL retained flat geodesic tower model.

Fixed temporal momentum can be included in this translation-invariant
model. K(sqrt(p^2+gamma^2))-K(|p|)=O_gamma(p^2 log p), using the
retained pK'(p)<=6K(p) and its large-p bound; low p is bounded.
It is another continuous remainder. The lower derivative Weyl
logarithms are integrable for the zero-trace H2 profiles; their
boundary contacts have extra powers of d. This does not establish
continuation of the ACTUAL thermal Dirac form or its constraints.

## 5. What this establishes and the actual remaining transfer

The persistent vector in [the coupled energy-class trace result](energy-shape-nonlinear-trace-forcing.md)
cannot be supplied by the compensated normal model on (5). However,
neither (5) for every actual auxiliary response nor a uniform
finite-cut actual-minus-model estimate has been proved here. Existing
geodesic-seam/core comparisons do not supply either implication at
nonzero boundary jets. Matching the local a4 contact alone also does
not fix the full nonlocal Hessian. The actual second mass vertex must
remain in that comparison; K itself omits its potential/tadpole row.

The next required theorem is a finite necessary limit of the actual
quadratic action minus this normal representative, with a justified
extension from the actual energy space. Preserve the charged seed,
paired projector, frequency, state, angular tower, fixed background
cutoff and original subtraction. A failure of that transfer is not
excluded by the normal result. No independent physical averaging,
new finite counterterm, causal mode or nonlinear solution is selected.

Verification checks the independent heat/contact calculation, analytic
continuation against direct subtracted Fourier integrals, Mellin and
one-sided kernels against digamma identities, moving-cut polynomial-
exponential profiles against their Mellin limit, and the improved
angular sum bound. Continuum convergence is an ordinary argument,
not inferred from a finite list of cutoffs. No Lean run is selected.

Remaining to-do list: actual finite-cut quantum comparison and its
energy-domain transfer; physical boundary phase and source compatibility;
continuous causal constraints, exact nonlinear existence and all-sector
stability; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
