# Location and proper-radial derivatives of the fixed quantile cutoff

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary argument
before frozen verification. **The original cutoff has quantified
support and five proper-radial derivatives. The projected source's
adjoint derivatives and signed coefficient remain to be controlled.**

Keep the ORIGINAL cutoff of the
[endpoint-preserving projection](endpoint-preserving-compatible-projection.md).
No cutoff, physical parameter, adjoint normalization or projection is
changed. Write t_h=1+h, y=sqrt(lambda)rho0, G0=F0/x and N=||chi0||_D.
The exact quantiles a,b and their transitions obey

    .06<y(a)<.2,
    .01h<a-t_h<.13h,
    1.01<b<2,
    y(a_-)>.04, b_+<3.                              (1)

Here a_-=(t_h+a)/2 and b_+=2b-t_h. In particular theta=1 throughout
1+1.13h<=t<=1.01. The projection generally subtracts a source in the
named cold region; equality of source germs at the bolt is a strictly
local endpoint statement. No nonzero physical residual is inferred.

In the same dimensionless seed proper coordinate,

    0<=theta<=1,
    |partial_y^j theta|<j! 10000^j, 1<=j<=5.          (2)

These constants contain no inverse power of h. They are derivatives
of theta ONLY, not of theta chi0, p_cut or the metric correction.
Conversion to the canonical normal-chart P45 norm is not asserted.

## 1. Transfer of the actual quantiles from the reference field

The [cold reference theorem](cold-reference-overlap-and-adjoint-budget.md)
uses the same normalization and supplies

    I=integral_0^infinity t^-4 dy, 17.5<I<31.7,
    Phi0=(phi_ref,0), ||Phi0||_D<4.

The [on-shell adjoint theorem](on-shell-adjoint-error-and-trace-gradient.md),
completed by the [finite trace/projection audit](finite-trace-gradient-and-projection-transversality.md),
gives ||chi0-Phi0||_D<.007. All components of chi0, including its
returned TT component, stay in this norm. For ANY measurable set J,

    | ||1_J chi0||_D-||1_J Phi0||_D |<.007.           (3)

For an initial interval, a tail, or an interval, the reference proof's
averaging estimate extends to its indicator:

    | ||1_J Phi0||_D^2-(1/2)integral_J t^-4 dy |<1e-5. (4)

This extension requires no derivative of an indicator. On y>=1e-8,
the scalar reference differs from cos(kbar y-pi/4) by <1e-9 in its
Liouville amplitude. The omitted interval contributes <3e-8. The
weighted indicator has total variation at most two; Stieltjes
integration by parts bounds its oscillatory integral by 10/kbar,
where kbar>4e17. The same error budget as the earlier smooth-cutoff
argument therefore suffices, uniformly in the interval endpoints.

Equations (3),(4) imply

    8.7<N^2<16.2,
    .0435<(1/200)N^2<.081.                          (5)

For the lower bound, sqrt(8.75-1e-5)>2.958 and
(2.958-.007)^2>8.7. For the upper bound use (4+.007)^2<16.2.
The probability C(t) defining a,b is exactly the D norm mass divided
by N^2: the constant sqrt(lambda) and Riesz normalization cancel.
Thus both the initial mass at a and tail mass at b equal N^2/200.

At y=.06, the initial reference squared norm is below .03001 and
its norm below .174. The actual initial mass is below (.174+.007)^2,
which is below .0435. Hence y(a)>.06.

At y=.2 the cold metric comparison gives t<1+1.13h, and t^-4>.99999.
The reference squared norm is therefore above .099989 and its norm
above .316. The actual initial mass exceeds (.316-.007)^2>.081.
Hence y(a)<.2. The geometric assertion used here follows independently
from the cold comparison in section 2; it does not assume the quantile.

On t>=1.001 the retained geometric bounds give

    5.99(t-1)^2/t^2 < G0(t) < 6.002(t-1)^2.          (6)

For the tail starting at t=1.01, restrict its proper-length integral
to 1.01<=t<=1.1. It is bounded below by

    log(10)/(1.1^4 sqrt(6.002))>.6.

Its reference squared mass exceeds .3-1e-5, its reference norm
exceeds .54, and its actual mass exceeds (.54-.007)^2>.081.
Consequently b>1.01. At t=2, instead use

    integral_2^infinity t^-4 dy
      < (1/sqrt(5.99)) integral_2^infinity dt/[t^3(t-1)]
      <=1/(12sqrt(5.99)).

Half this bound plus 1e-5 is below .0172<.132^2. The actual tail
mass is below (.132+.007)^2<.0435. Thus b<2. Continuity of the
cumulative mass justifies the strict quantile conclusions even if
it has flat intervals. No microscopic oscillations are numerically
resolved and no reference quantile is substituted for an actual one.

## 2. Thermal geometry and transition coordinate derivatives

The exact seed has s=t-1, v=1+h,

    G0=s^2 fbar(t)-deltabar/t,
    fbar=1+2/t+(3+1/x)/t^2+(d/x)/t^4,
    deltabar=v h^2 fbar(v), x>=1000, |d|<=1,
    1/(8*2000^10)<=h<1/2000^10<1e-30.

On h<=s<=.001 the earlier proof gives

    5.98(s^2-h^2)<G0<6.02(s^2-h^2).

Integration and inversion give

    cosh(sqrt(5.98)y)<=s/h<=cosh(sqrt(6.02)y).        (7)

The endpoint s=.001 is reached after y>.2, as follows already
from the lower bound acosh(.001/h)/sqrt(6.02). Hence (7) is valid
at y=.2. Since cosh(.2sqrt(6.02))<1.13, it proves the assertion
used in section 1. At y(a)>.06 it gives
s(a)/h>1+5.98(.06)^2/2>1.01. Finally
cosh(.04sqrt(6.02))<1.005 proves y(a_-)>.04.

For explicit derivative bounds it is convenient to expand

    G0=t^2+e-(4+2e+deltabar)/t
          +(3+e+j)/t^2-2j/t^3+j/t^4,
    e=1/x, j=d/x, 0<e<=.001, |j|<=.001,
    0<deltabar<7h^2.

The exact monotonicity proof gives G0'>0. Also fbar'<0 and
fbar<=6.002, so G0'<12.004s+7h^2<13s for s>=h.
Termwise differentiation of the expanded expression gives, for t>=1,

    |G0''|<29, |G0'''|<100, |G0''''|<460.           (8)

The radial equation t'=sqrt(G0), prime now meaning y derivative,
gives the EXACT identities

    t''=G0'/2,
    t'''=G0'' t'/2,
    t''''=[G0'''(t')^2+G0''t'']/2,
    t'''''=[G0''''(t')^3+3G0'''t't''+G0''t''']/2.    (9)

For the rising transition write q=(a-t_h)/h, .01<q<.13, and
xi=(t-a_-)/(qh/2). Here q/2< s/h-1 <q. The cold comparison gives

    |xi'|<=2sqrt(6.02(2+q)/q)<72,
    |xi''|<13(1+q)/q<1500,
    |xi'''|<1044, |xi''''|<22000, |xi'''''|<16000.   (10)

For the last three bounds use (8),(9), t'<2h, t''<8h and h<1e-30.
For example |xi''''|<(100/2)(2h)(72)+(29/2)(1500), and
|xi'''''|<(460/2)(2h)^2(72)
 +(300/2)(2h)(1500)+(29/2)(1044).

For the falling transition put d_b=b-t_h, 0<d_b<1,
xi=(b_+-t)/d_b. Then s<=2(b-1) and d_b=b-1-h.
Since b-1>.01 and h<1e-30, (b-1)/d_b<1.000001.
Consequently the safe bounds from (6),(8),(9) are

    |xi'|<5, |xi''|<14, |xi'''|<73,
    |xi''''|<1500, |xi'''''|<41000.                (11)

For instance t'<5d_b, t''<14d_b, and d_b<1 give fourth and
fifth bounds (100*25+29*14)/2<1500 and
(460*125+300*5*14+29*73)/2<41000.

Both transitions therefore have divided jets, through order five,
bounded by the coefficients of

    Q(X)=100X+1000X^2+400X^3+1300X^4+500X^5.         (12)

## 3. The fixed smooth step and cutoff derivatives

Let H(v)=exp(-1/v) for v>0, zero otherwise, and
S(v)=H(v)/[H(v)+H(1-v)], exactly as in the original projection.
On 0<v<1 put u=1/v and define polynomials

    P0=1, P_(j+1)=u^2(P_j-P_j'),
    H^(j)(v)=exp(-u) P_j(u).

For j>=1 a rigorous rational majorant is

    A_j=(1/j!) sum_m |[u^m]P_j| (m/2.7)^m,

because sup_(u>=1) u^m exp(-u)<=m^m exp(-m) and e>2.7.
The denominator of S is at least exp(-2)>1/8: one of v,1-v
is at least one half. Its divided j-th derivative is bounded by
2A_j. Set B0=1 and recursively

    B_j=8[A_j+sum_(i=1)^j 2A_i B_(j-i)].            (13)

Differentiating (H(v)+H(1-v))S(v)=H(v) proves
|S^(j)|/j!<=B_j. Exact rational arithmetic gives B_j<50^j for
1<=j<=5. Flatness extends the bounds to v=0,1 and the constant
parts outside [0,1]. This argument does not presume a uniform
complex analytic radius for the nonanalytic smooth step.

Taylor's chain rule for finite jets now bounds the j-th divided
proper derivative of S(xi(y)) by

    [X^j] sum_(n=1)^5 50^n Q(X)^n <10000^j,
    1<=j<=5.                                      (14)

The transitions are disjoint: a<1+1.13h<1.01<b. On each transition
the other factor defining theta is identically one. Elsewhere theta
is constant. There is consequently no extra product factor in (14).
This proves (2) for the EXACT fixed quantile cutoff, despite the
very small t-coordinate width of its first transition.

## 4. Multiplier and commutator consequences, with their limits

For any positive weight w(y) define the proper-radial norm

    X_(5,M,w)(v)=max_(0<=j<=5) sup |partial_y^j v|/[w(y) j! M^j].

The pointwise Leibniz rule, with no derivative of w, yields

    X_(5,M,w)(theta v)
      <=sum_(i=0)^5 (10000/M)^i X_(5,M,w)(v).         (15)

At M=100 the factor is 10101010101<1.011e10. In particular a
cold weight can be used if v is already controlled in that same
proper-radial norm. This is not a proof that the actual adjoint or
correction lies in the cold metric class, nor a chart conversion.

For the seed radial scalar Tbar=-partial_y^2-p partial_y,
p=G0'/(2sqrt(G0))+2sqrt(G0)/t, the transition bounds also give
|p|<40. On the first transition use s/h<=1.13,
(s/h)^2-1>.005(2.005), and the cold lower bound. On the second,
use t<3, (6), G0'<13s and sqrt(G0)<sqrt(6.002)s.
Let J be the union of the transition intervals. Then exactly

    [Tbar,theta]v=-2theta' v'-(theta''+p theta')v,
    ||[Tbar,theta]v||_(L2(w dmu0))
      <2e4 ||1_J v'||_(L2(w dmu0))
          +2.1e8 ||1_J v||_(L2(w dmu0)).            (16)

For zero right side read the strict inequality as <=. No regularity
of the positive measure weight is needed. In particular

    (Tbar-E)(theta v)=theta(Tbar-E)v+[Tbar,theta]v.

This identity retains the commutator if used in an on-shell source
split. The FULL coupled chi0 is not silently replaced by a solution
of (Tbar-E)chi0=0.

The unresolved source is p_cut=theta t^-4 chi0/(alpha sigma0),
sigma0=||t^-4 chi0||_H, .99<=alpha<=1. Equations (2),(15) remove
one unquantified multiplier from its higher-regularity estimate.
They do not estimate the rapidly oscillating chi0 derivatives,
physical normalization, or signed coefficient L_ell^0(f)/sigma_ell.
The existing finite value/H1 response remains unchanged. No cold
metric bound, full higher/log inverse, nonlinear preserved ball,
projected nonlinear solution or physical compatibility value follows.

## Frozen verification contract

Freeze this note and checker before execution. Reconcile the finite
range response, full adjoint error, reference overlap and fixed
projection source hashes. Check all rational quantile thresholds,
cosh enclosures by positive series with a geometric remainder,
thermal derivative coefficient sums, the exact chain-rule identities,
step derivative/quotient recursions, finite-jet composition and
commutator constants. An independent direct differentiation of H
checks its polynomial recurrence; a nonconstant symbolic test checks
the scalar commutator and detects omission of its first-derivative
term. These are algebraic fixtures, not physical adjoint sampling.
The interval averaging and continuum norm transfers are ordinary
proofs. Preserve every frozen execution and failed source. No Lean
campaign is selected and no independent review is asserted.

The first frozen execution caught an erroneous constant -3 in the
displayed seed expansion. Expanding the original product cancels
that constant exactly. It is removed in this source; the derivative
majorants and all stated cutoff constants are unchanged. The failed
run and its exact note/checker sources are retained in the receipt.

Remaining to-do list: actual adjoint/source higher derivatives and
signed projection coefficient; cold/local and higher correction;
nonlinear mapping, preservation and solve; signed physical
compatibility; exterior/state, charge/return, opening, stability,
transport, preparation, separate creation and full-result audit.
