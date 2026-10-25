# Full transition cost of the reference mass lift in the original graph

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: NO PROGRESS. The full stable semiclassical wormhole
objective remains active and incomplete; the feasibility audit stays closed.

OUTCOME: the complete physical coupled return is UNRESOLVED. The current
attempt to quantify its differentiated core did not produce that bound.
Before assigning any improved core estimate to the original admission
budget, this calculation evaluates the FULL already specified boundary-basis
conversion. Its cutoff transition, omitted from the preceding LOWER bound,
has a large cost. This is an earlier requirement for the proposed certificate,
not a physical source, a new boundary condition, or a construction obstruction.

## 1. Same mass lift, endpoints, cutoff and original norm

Keep the exact reference profiles and actual finite-seed Q:

    J_m(z)=integral_0^z t^3/Q(t)^2 dt,
    n_m=-2 beta Q J_m/z,
    b_m=beta z^3/Q+beta Q' J_m.

Both local reference rows cancel. The independently retained mass cutoff
chi_m equals one through .1 and zero by .2. Its complete areal return is

    f_cut=2xz chi_m-2F n_m D chi_m,
    alpha_cut=n_m D chi_m, D=z partial_z.

It has zero reference bolt displacement and total reference mass -xB per
unit coefficient. It is NOT a global fixed-temperature zero mode. This
calculation changes neither chi_m nor its original matching equations.

The ORIGINAL graph cutoff chi is the already specified beta(11,11)
logarithmic cutoff, one through a=.05, zero from b=.1. Put

    t=log(z/a)/log(2),
    chi(t)=1-[21!/(10!)^2] sum_(k=0)^10 (-1)^k binom(10,k)t^(11+k)/(11+k)

on 0<t<1. It is C10 at the endpoints, sufficient for the H9 quantities
below. No derivative of a sharp indicator is used. Set sigma=4.005 and

    w= z^-sigma chi(b_m-z^3),
    d_chi=B_sigma w,
    g_old=g_lifted+m d_chi.                              (1)

The known-source subtraction and every other phase are the same in (1).
This exact equality holds before any physical operator or regulator is
estimated. Neither graph is silently substituted for the other. The norm
is the original Fourier H3(dx), x=log z:

    ||v||H3^2=sum_(j=0)^3 binom(3,j) integral |partial_x^j v|^2 dx.

## 2. Complete normal operator, analytic profile and certified tails

Use the SAME complete normal symbol, including its digamma/cotangent
finite part and fixed contact. Split it only for estimation:

    B_cl(s)=K L(s-4)s(s-3),
    K L=2 beta Fstar/3,
    Fstar=F_l+23G lambda/(180pi),
    B=B_cl+B_q.

The existing all-frequency normal-symbol proof and H9 bound for THIS
.05/.1 cutoff give

    ||B_q,sigma w||H3 <2.260e-10.                       (2)

This is the complete normal quantum part, not a finite angular or Fourier
selection. It is not the unknown actual finite-domain quantum return.

Write q=Q/beta and q^-1=sum a_n z^n. The exact recurrence is

    a_0=1, a_n=-sum_(j=1)^min(6,n) q_j a_(n-j),
    q^-2 coefficient c_n=sum_(j=0)^n a_j a_(n-j).

Then b_m=z^3 q^-1+q' integral_0^z t^3 q^-2 dt. This gives its coefficients
without numerically integrating a nearly cancelling expression. Retain
all coefficients through N=80, subtract the exact mass cubic, and call
the resulting profile w_N. No physical spectrum is truncated by this
finite representation of a known analytic reference profile.

The preceding complex-disc proof gives |(b_m-z^3)/z^5|<2 on |z|=.2.
With R=.2, this implies |b_n|<=2 R^(5-n) for every n>=5. For k<=6,
the L2(dx) norm on 0<z<b of the omitted kth x derivative is bounded by

    T_k=2 R^5 b^-sigma (N+1+sigma)^k (b/R)^(N+1)
          /[sqrt(2(N+1-sigma))(1-r_k)],
    r_k=(b/R)(1+1/(N+1+sigma))^k<1.                  (3)

This follows by integrating EACH power on the whole half-line first,
then comparing the positive series geometrically. Using n+sigma in the
numerator only enlarges |n-sigma|. For the cutoff use B_0=1 and

    B_l=21!/(21-l)! *2^(l-1)/log(2)^l, 1<=l<=6.    (4)

Indeed chi has Bernstein coefficients eleven ones followed by eleven
zeros. Its lth derivative coefficients are the lth differences of that
sequence times the falling factorial; their absolute values are at most
2^(l-1). The Bernstein basis is nonnegative and sums to one. Thus (4)
bounds every x derivative, including the full transition.

Expand P(D)=(D+sigma-4)(D+sigma)(D+sigma-3)=sum p_r D^r. For j=0,...,3
the omitted classical output derivative is at most

    A_j=|K L| sum_(r=0)^3 |p_r| sum_(l=0)^(j+r)
                              binom(j+r,l) B_l T_(j+r-l).

Consequently E_series=sqrt(sum binom(3,j) A_j^2) bounds the ENTIRE
omitted classical H3 norm. Add (2) by both triangle inequalities.
All bounds use the actual enclosed seed constants, not decimal point
approximations or changes of normalization.

## 3. Whole-interval enclosure of the retained norm

On z<a, chi=1. With k_n=n-sigma and
A_n=K L b_n(n-4)n(n-3), exact integration gives

    I_inner=sum_(n,m=5)^N A_n A_m a^(k_n+k_m)
                                   (1+k_n k_m)^3/(k_n+k_m).  (5)

All cross terms retain their signs. On the transition define
E_k(t)=z^-sigma sum_(n=5)^N b_n(n-sigma)^k z^n. The exact jets are

    g_j=K L sum_(r=0)^3 p_r sum_(l=0)^(j+r)
                                  binom(j+r,l) (D^l chi) E_(j+r-l).

Partition the ENTIRE interval t in [0,1] into 2048 closed cells. On each
cell, evaluate z=a exp(log(2)t) and all E_k by outward Arb arithmetic.
For each polynomial derivative of chi, Taylor expansion about the exact
cell midpoint is FINITE and EXACT; use the sum of absolute nonconstant
terms as a radius. This avoids treating a sampled midpoint as an integral.
It encloses chi's derivatives on every point of the cell. Square the
resulting g_j intervals, using a nonnegative lower endpoint, and integrate
the interval enclosures with dx=log(2)dt. Summation encloses I_transition.
The outer interval z>b contributes zero. Endpoint regularity creates no
additional distribution through the orders in this calculation.

Taking sqrt(I_inner+I_transition) and then subtracting/adding E_series
and (2) gives a certified lower and upper bound for ||d_chi||H3. The
execution contract is to verify a lower bound above 2000 and an upper
bound below 2400; the actual, narrower interval is retained in the receipt.
The previous lower bound .00906808164 remains true but is insufficiently
informative because it used only the interval 0<z<.02.

## 4. Consequence for the COMPLETE coupled certificate

Let D_- and D_+ be the certified conversion bounds. Retaining BOTH original
and lifted graph ceilings r=1e-11 necessarily requires, by (1),

    |m| D_- < 2r.                                    (6)

This is a necessary condition for these two sufficient graph targets,
NOT for wormhole existence. Conversely the explicit split

    ||g_lifted||H3<=r/2, |m|<=r/(2D_+)               (7)

is sufficient for the graph conversion alone. Neither input in (7) has
been proved for the actual solution. In particular a lower bound on D
was never an upper bound that could establish (7).

Test the SPECIFIED positive estimator that substitutes the existing
conditional coefficient allowance

    |m|<=m0+49.462 E_lifted, m0=4.2953e-8

into ||g_old||<=E_lifted+D_+|m|. Even its best possible D consistent with
D>D_- makes the reported affine conversion allowance exceed D_- m0,
and its reported graph slope exceed 1+D_-49.462. These are lower bounds
on THAT COMPUTED UPPER-BOUND EXPRESSION. They are not lower bounds on
the actual m, actual physical graph, or actual feedback coefficient.
They show this independent-positive-bound conversion cannot certify the
target with the available mass allowance, even if a remaining quantum
column could be made arbitrarily small. Its failure precedes another
optimization of the controlled cold, scalar or L2 angular remainder.

The original mass-column target kappa_m<=5e-5 remains UNKNOWN. If proved
in the intended post-elimination norm, its unchanged debits are
2.14765e-12 affine, .0024731 graph and .0195239 per unit energy. After
the controlled far entry, the other coupled entries must still satisfy

    Y_other<=5.34234999994958e-12,
    eta_other+(.0195239+4.583e-13)gamma_E<.0012069.

Those inequalities do not themselves establish (6) or (7). Nor is D
silently added a second time to a column already certified in the signed
original output. The correct complete output is (1), with its correlations.
The scalar quartic's ORIGINAL row coefficient 20 beta Fstar, exact mass
compatibility, correlated slow pair and regular-bolt matching all remain
in the unsolved coupled system. The fast screen .135973<.1599164 keeps
its two-fast-input admitted-family scope, including the previously proved
fixed reference-lift background extension. No independent phase bound,
nonlinear ball comparison, actual signed center or stability follows.

EXECUTED SCOPE: full original-normal-graph coordinate conversion, all
cutoff-transition derivatives, analytic-series error, complete normal
quantum allowance, and the resulting certificate thresholds. The actual
common-regulated finite-domain core and coupled solve remain UNEXECUTED.
Their first unsupported quantitative step remains the differentiated
kernel/contact bound with the coupled inverse. This calculation does not
claim to have replaced that step or completed the requested milestone.

RESEARCH DECISION: retain the physical construction, but stop using
independent positive mass and lifted-graph allowances to support original
admission. The measured transition cost makes a direct signed original-
output certificate or a much sharper ACTUAL mass enclosure necessary for
this route. A different physical embedding could avoid imposing both
graph ceilings, but would need a proof at the unchanged physical radius.
Changing cutoff, basis, units or normalization alone supplies none.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): enclose the actual mass and
the correlated original graph g_lifted+m d_chi in the SAME coupled solve,
with both original equations, scalar/slow/interior/bolt rows and the common
regulator. Resolve whether the actual mass meets (6), or whether a signed
original-output/physical embedding can certify admission without the
additional lifted ceiling. Further unsigned core constants alone cannot
settle this now-quantified conversion requirement.

Remaining to-do list: actual coupled return and signed correction; original
physical-norm/neighborhood comparisons; nonlinear admission/contraction and
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
