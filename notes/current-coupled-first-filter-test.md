# Complete first coupled filter on the accepted analytic input

2026-10-29 HST, E003/S169. Ordinary argument and decision contract BEFORE
execution. The previous goal turn was NO PROGRESS (a status/plan answer).
Preserve 2d8ef11, P32, the closed audit, original norms and the full stable
semiclassical wormhole objective. This tests one COMPLETE filter on that
input, not another improvement to its already sufficient transport error.

## Domain, coupled residual and the first new implication

Use the original scalar form space Es, TT form space Et and Hmu=L2(z^-4 dz),
with es=a/N, et=ell_t/N, N=|W_action| beta. The positive scalar and TT
realizations are those of current-seed-full-hessian-realization.md. Keep
the actual paired mass/log phase, complementary first-order matter adjoint
and composed squares, fixed period/spin/scale, physical cutoff transitions
and regular disk bolt. York uses J(0)=0, (Q^2 J'/z^4)(L)=0. No independent
phase perturbation or nonlinear metric-neighborhood estimate is implied.

The complete normalized coupling has norm <=sqrt(qt), qt=.338/(1-epsilon).
For xhat of 2d8ef11 define the actual filter by

    (ell_t+B)y=B xhat, B=b* a^-1 b, F_t=(I+T_t)^-1 T_t.

Choose scalar sigmahat and TT yhat jointly. Define on the ENTIRE form duals

    rs(w)=b(w,xhat-yhat)/N-es(sigmahat,w),
    rt(v)=b(sigmahat,v)/N-et(yhat,v).

Eliminating the scalar equation ON ITS ACTUAL DOMAIN gives

    ||F_t xhat-yhat||Et <= Rt+sqrt(qt) Rs = Rfilter.       (1)

Indeed the residual of (I+T_t)yhat=T_t xhat is rt plus the coupling
applied to the scalar Riesz solution of rs. The inverse norm is <=1.
There is no reference-coercivity conversion in (1). Both equation errors
are necessary; bounding only rt does not certify this filter.

The first new domain implication is that the jointly chosen rational
profiles and cutoff scalar are in these fixed form spaces. Establish it
before interpreting intervals. Put H=-partial_x((2x+a x^2)partial_x),
x=(L-z)/gamma, with its already accepted bolt Green spectral realization.
H is only an auxiliary profile generator, not the matter operator. Let
T be the EXACT finite-Q two-row reconstruction from 1b04eae/f384d0f,
including its entire Volterra return. The endpoint bounds for T extend by
interpolation from the original core; no new boundary is introduced.

Write A,muT,muS,muY,K0,Cv,lambda_l,lambda_n as in 8dc41ff and set

 c=2 tau0/(3 L^2), delta=d L^2, kappa=4 tau0 c/(beta L^4 A),
 D(p)=(p+muT)(p+muS)+kappa(Cv+delta p),
 Hfac(p)=1+(Cv-1)lambda_n/(p+lambda_n),
 f1=At K0 Hfac(H)/[(H+muT)(H+muS)(H+lambda_l)] delta_b,
 f2=kappa(Cv+delta H) D(H)^-1 f1,
 fdiff=f1-f2=At K0 Hfac(H)/[D(H)(H+lambda_l)] delta_b,
 xi0(f)=-(2/3)(H+2a)(H+muY)^-1 f,
 sigma0=c(H+muS)^-1 xi0(fdiff).

Take yhat=T f2, xhat-yhat=T fdiff, sigmahat=chi_S sigma0, S=1e-12,
with the same quintic scalar cutoff (one through S, zero after 2S,
|chi'|<=2/S, |chi''|<=15/S^2). T means the TT field, and V T f its
canonical quantum input. The new cutoff is computational and its complete
commutator is charged; it changes no physical boundary condition.

D has positive coefficients. The calculation verifies 0<=Fref(p)<=.251
for all p>=0 by positivity of the coefficients of .251 D-kappa(Cv+delta p).
Thus f2 is in the accepted order-9/2 input space; T f2 is in Et by the
SAME fixed-seed extension. Its leading ultraviolet profile is O(p^-4),
with retained rho^6 log(rho) behavior, while fdiff remains O(p^-3).
sigma0 is O(p^-4), in the scalar second-order form domain; multiplying
by C2 chi preserves that domain and introduces no distribution beyond
the already bounded weak rows. The original outer mass/log realization
and bolt flux follow from the same T formulas, not from a smooth surrogate.

The common-regulator form theorem and kernel/contact bounds of 2d8ef11
apply on this Et trial. No new heat-rate theorem, finite matter-mode
selection, or seed-to-neighborhood transfer is used. The zero canonical
scalar/mixed quantum row uses the SAME fixed-domain Ward identity as
the accepted scalar-source argument; all local anomaly/contact remainders
remain in Kbrem and KA. It is not a freely chosen stress tensor.

## Cancellation and quantitative profile control

Keep the exact cancellation of BOTH factors (p+muT)(p+muS) in fdiff.
The two reference equations then hold algebraically before norms:

    (L^2 H+4 beta)sigma0=(2 tau0/3)xi0(fdiff),
    A(H+muT)xi0(f2)=4 tau0/(beta L^4)(Cv+delta H)sigma0.

Retain the rational classical remainder (r1 H+r0)/(H+2a)^2 separately.
No root separation of D is used. Write D(-mu)=(mu-h)^2+e,
h=d1/2, e=d0-h^2. Verify |h-4/3|<.01, |e|<.001 and group its two
poles, together with muY when present, on |mu-4/3|=.2. The denominator
is bounded below by (.2-|h-4/3|)^2-|e|. Other contours use their positive
distance from h and the same subtraction of |e|. The existing near-2
circle of radius .25 is retained, as are isolated muS and lambda_n.

The positive Legendre integral yields |Q_k|<=Q_Re(k). On the new circle,
Re(k)>.65 from Re(mu)>4/3-.2 and a<1.001. Its differentiated integral
is bounded by twice the real .65 flux: |k+1|<2(.65+1), while the
exponential weight is bounded by that for .65. The old circle retains
.9 and flux factor two. Cauchy bounds use the FULL rational numerator
and denominator; isolated residues are combined before absolute values.
These contours represent complete auxiliary profiles, not truncated
angular/thermal quantum sums. Every small high-mass coefficient is retained.

For scalar spectral norms set p=t/(1-t). If sigma0/At=(1-t)^4 C(t),
the accepted complete spectral density <=1/2 gives, for j=0,1,2,

    ||H^(j/2)sigma0||2^2
       <=At^2 gamma/2 int_0^1 t^j(1-t)^(6-j)|C(t)|^2 dt. (2)

Enclose every whole interval, including t=1 (the infinite spectral tail).
All parameters and arithmetic use 512-bit balls. No finite quantum tower
is substituted. For profile f retain the complete Nf=||f||2, Mf=||s f||2,
Hf=||Q2 f_s||2 and full tails from the established Q integral. Pointwise
cutoff bounds T_f,U_f retain |f| and |Q2 f_s| on [S,2S].

## Actual reconstruction and complete scalar row

In the notation of 2d8ef11, V T f=A0 f+B0 r+C0 h_r+E0 W,
r=(H+muY)^-1 f, h_r=q2 r_x. Retain exactly its A0-1,B0,C0,E0 and W.
The actual TT coordinate is

 xi(Tf)=-(2/3)A0 f-(B0+4a A0/9)r-(2/3)C0 h_r-E0 W.   (3)

This also follows from xi=2nu-2QJ/z; nu alone is NOT xi/2 in finite
geometry. Subtract xi0(f) before norms. To bound the actual Hmu norm,
perform the same subtraction on xi/z^2, replacing A0 by z^3 R/L^5,
factoring z^3 R/L^5-1/L^2, and dividing B0,C0,E0 by z^2 algebraically.
There is no singular division at z=0. Integrate all whole physical cells
and retain dz=gamma dx. Both equations and the entire W return remain.

Let exi, exiw bound the two differences for fdiff; then
||xi(T fdiff)||Hmu <= ||xi0(fdiff)||2/L^2+exiw. The full outside-S
source is bounded by the model xi tail plus exi. With zm=L-2S,
r=M3*2S/(6a), M3=1000, CE=2L+L^2 M3/(6a),
CB=5L^2 M3/(6a)+2L(1+r), the complete second-order scalar residual has

 rP <= (2 taumax/3) exi/zm^2
       +(2 taumax/3)(Tail_xi+exi)
       +2 M_xi/zm^2+(CE M_Hsigma+CB H_sigma)/zm^2+Cut_P. (4)

Here |tau(z)-tau0|<=3s, taumax=3.01 beta, tau/z^2<=taumax; the latter
controls the ENTIRE exterior source although sigmahat is cut off.
Cut_P is the full product-rule commutator from the accepted scalar proof:
sqrt(S)/zm^2 times
L^2 Qmax(15 T_sigma/S^2+4 U_sigma/(S Q2(S)))
 +(L^2 Qsmax+2L Qmax)2 T_sigma/S.

Use the original scalar form comparisons S0=1/sqrt(23.99),
S1=sqrt(beta/5.99), S2=sqrt(beta/(6d)), and all five Kbrem components:

 Rs <=6(Cv S0+d S2)rP/beta
       +(KA S1 G_sigmahat+2w0 S0 N_sigmahat)/beta
       +(Kbrem/beta)||xi(T fdiff)||Hmu.                (5)

No accepted initial-source replacement is charged again here.

## Complete TT row and full quantum response

Compute the cutoff scalar Es norm from (2), geometry and its full
commutator, exactly as in 8dc41ff. The complete source comparison is

 Esource =1.01 Kbrem ||sigmahat||Es/beta
   +4*1.01/beta [T0(Cv Tail_sigma+d Ndiff)
       +14|epsseed|(Cv N_sigmahat+d D_sigmahat)],
 Ndiff=CE M_Hsigma+CB H_sigma+Cut_Delta+L^2 Tail_Hsigma.

Re-evaluate the complete local two-row and York-adjoint return of
1b04eae on the NEW f2 profile, using both rows, all outer log coefficients
and the full finite interval. Reuse the linear endpoint theorem only to
bound Egeo and the classical rational remainder by .251 times their
accepted f1 bounds. This follows from commuting Fref(H) with the auxiliary
spectral weights, not from scaling an actual-domain pointwise field.

Recompute actual U0,U1 and the complete cold-complement tail for V T f2
with the whole-cell reconstruction. The prior uniform Hessian comparison
scales by .251 by that same INPUT spectral inequality. Charge all ten
entries of 2d8ef11: local gain-two, full annular kernel, finite second
contact, row-subtracted global gluing, cold-cap and entire outside contact,
constant-contact gluing, curvature polynomial, vacuum uncertainty and
retained energy comparison. The UV logarithm, common regulator, complete
angular/thermal sums and finite-domain returns are unchanged.

Do not reuse the old trial's energy for the last two entries. If their
sum is cE ||yhat||Et, the actual identity gives

 ||yhat||Et <=sqrt(qt)||sigmahat||Es+Rt_other+cE||yhat||Et.

Verify cE<1 and solve this scalar inequality. This supplies its actual
energy bound before charging these two entries. Finally

 Rt=Esource+Emodel+Elocal+Egeo+Eactual_difference,
 Rfilter=Rt+sqrt(qt) Rs.                               (6)

Every actual form functional is bounded on ALL tests. No unknown source,
boundary, graph or application constant is replaced by a ceiling.

## Original output contract and limits

There remains one exact F_t after this filter. With ft=qt/(1+qt), the
new error has bolt coefficient Kfilter=pref_t ft=Kfirst/ft. The prior
first-solve error still has coefficient Kfirst=pref_t ft^2. The initial
scalar-source/input errors remain charged once in their accepted subtotal.

Before execution choose the SAME physical error amount as the preceding
1e33 first-solve screen: Kfirst*1e33. Thus this new sufficient analytic
filter screen is Rfilter<=ft*1e33, NOT 1e33 or 9.97795e32. The latter
was a conditional inner-iteration allowance for a DIFFERENT algorithm.
None of these chosen accounts establishes free physical margin. Report
the known cumulative error subtotal and all still-unknown correlated
signed terms, remaining filter, numerical caller and other outputs.

An exact analytic trial has no sampled interpolation error. Coefficient
and bound arithmetic are enclosed here; an implemented numerical caller
still needs its own Et error. If the analytic screen passes, report the
remaining caller threshold but leave that error UNKNOWN. If excessive,
quantify the loss against the unchanged screen. A residual certificate
on this exact trial is not completion of signed G_t or geometric admission.

Decision rule: retain the coupled method if its full bound is affordable;
if a particular estimator is excessive, quantify its responsible entries
and select a test of their cancellation. An undefined domain step remains
unresolved rather than receiving a numerical value. No result is an
obstruction to the construction unless a separate argument establishes it.

## Executed first estimator and the specified sharper test

The immutable first receipt current-coupled-first-filter-test.json and its
author review enclose the COMPLETE analytic filter error by <4.681e32.
The chosen screen is 2.52615844544095665e32, so this report exceeds it by
about 1.853. The TT row costs <1.9173e32; the scalar residual costs
<4.752574e32 before its sqrt(qt) transfer. Its mixed remainder is
4.5398351e32, because the first estimator used
||xi0(fdiff)||2 <= ||f1||2 <2.837410e48. This is a valid but excessive
upper estimator, not a failure of the coupled inverse or construction.

BEFORE the next execution, replace ONLY that comparison by the complete
structured spectral integral. From the already cancelled rational profile,

 xi0(fdiff)/At = -(2 K0/3)(p+2a) Hfac(p)
                   /[(p+muY)(p+lambda_l)D(p)].

With p=t/(1-t), write this as (1-t)^3 Cxi(t). Then

 ||xi0(fdiff)||2^2 <= At^2 gamma/2 int_0^1 (1-t)^4 |Cxi(t)|^2 dt.

Use the same 512 whole spectral cells including the infinite endpoint.
No high-mass coefficient, quantum mode or boundary return is deleted.
This is the norm of the actual correlated source model, followed by the
already enclosed finite-Q difference exiw, not a smooth surrogate or a
normalization change. The two-equation identity and actual trials are
unchanged. All other full-map entries are reused from the immutable
completed calculation; their source and argument snapshots are retained.

The complete-ledger threshold for xiH is

 [screen-Rt-sqrt(qt)(Rs_core+Rs_scalar_remainder)]
                                      /(sqrt(qt) Kbrem/beta).

Compute this threshold and the actual new xiH; never insert the threshold
as a bound. This isolated change is justified by its direct contribution
to (6), which the next receipt must reassemble completely. Preserve the
first estimator's excessive result and distinguish the numerical caller
error (still unknown) from this exact analytic residual.

Remaining to-do list: complete filter/caller errors; signed G_t/Z_r and
every original coupled, source, graph, physical, boundary and mass budget;
nonlinear admission/contraction/existence, exact compatibility, matching,
semiclassical validity/noise and all-sector stability; the full objective.
