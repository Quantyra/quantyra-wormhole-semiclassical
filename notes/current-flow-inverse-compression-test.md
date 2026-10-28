# Test the direct flow-inverse preconditioner in the original energy norm

2026-10-27 HST, E003/S169. Ordinary argument before interval execution.
The preceding email review supplied no new goal bound. Preserve the full
T_t-(I+T_t)F_n objective, P32, original physical norms, closed audit and
all accepted results. This tests the missing complementary inverse, not
the already affordable outer-input omission or transport error.

## Candidate, spaces, and the first domain implication

Use the ACTUAL paired mass/log form domain E_t with norm ell_t/N,
N=|W_action| beta, beta=lambda B_radius^2/C. Its frozen York shape map is

    V xi=-xi-(Q J)'/2,
    -(Q^2 J'/z^4)'=Q(4xi/z-xi')/z^4,
    J(0)=0, (Q^2 J'/z^4)(z_h)=0.

Let X=V E_t have the transported norm, on the fixed-gauge range and its
form completion, as in current-schur-form-preconditioner-test.md. Its
dual uses the density pairing integral f T dz. This is not the original
geometric admission norm and does not include a new output comparison.
Both original rows remain in B=b* a^-1 b. In these transported coordinates
the complete reference A_X=(ell_Q+B)/N satisfies A_X>=a R_X,
a=1-e_Q, where R_X is the ANALYSIS Riesz map of X. R_X^-1 is not used as
a numerical subroutine. The actual-minus-reference comparison and the
separate K application retain their existing scopes.

Write V0 T=(2pi)^(-1/2)T for the normalized base lift to the ACTUAL unit
tangent bundle. P=iX_flow is the complete generator on its C2,1 paired
double. In particular it is not the scalar radial Laplacian. Keep the
complete positive multiplier

    A_flow(w)=C w^4 chi(w)/(90pi^2),
    a_flow(w)=4pi G A_flow(w)/(B_radius^2 beta),
    r(w)=1/(1+a_flow(w)),
    C_flow=V0* r(P) V0.                                  (1)

Here chi is the full-tower multiplier, with the exact common fourth
contact and UV logarithm. chi(0)=-F_eff>0, r(0)=1 and 0<r<=1. The factors
in (1) keep the original action normalization: the quantum shape form
is (2G/(B_radius^2 beta)) times the unnormalized angular compression.
There is no changed temperature, scale or norm. P includes every
crossing/trapped trajectory and regular-bolt return. Flow frequency zero
is not a thermal zero mode or an omitted matter sector.

This gives a bounded positive map on H0=L2(dz). The proposed inverse
step is to use it as a map X* -> X, optionally localized as chi_I C_flow
chi_I and combined with bounded coarse/endpoint corrections. Here chi_I
is smooth, equals one on [.49,.51], and vanishes outside (.48,.52).
Its transitions are part of the candidate. The first unsupported
implication is precisely that dual-to-energy boundedness, before any
claim about an implemented complete inverse or original physical output.

The test below rules out ANY bounded extension agreeing with this
candidate on compact sources. If an output is outside X, it fails the
domain requirement already. If outputs belong to X, their operator norm
is unbounded. Adding a fixed energy-bounded correction cannot repair
unboundedness. No global inverse for A_X is assumed in this argument.

## Actual-domain localized sources and their full dual norm

Take z0=1/2, 0<h<=1e-3, phi(x)=(1-x^2)^8 on [-1,1], zero outside, and

    f_h(z)=h^(-1/2) phi''((z-z0)/h),
    Phi2=integral_(-1)^1 phi(x)^2 dx.                    (2)

These are ordinary compact source densities in the actual domain's
dual, not prescribed quantum stresses or replacements of B's range.
Their two vanishing moments are exact. For every actual shape T,
integration by parts on this compact support gives

    |integral f_h T dz| <= h^2 sqrt(Phi2) ||T''||L2(I),
    I=[.49,.51].                                        (3)

The original joint coercivity retains its fourth reserve on I, where
the original inner partition equals one. With eps the accepted retained-
transport error, the physical period/action conversion gives

    k2=(1-eps) G C/(50pi B_radius^2 beta),
    ||Delta T||L2(I) <= k2^(-1/2)||T||X,
    integral Q |T'|^2 dz <= g1 ||T||X^2,
    g1=1e5 beta/[(1-eps) .99 F_lambda].                 (4)

Equation (4) uses C S2/(100pi) and C J1/E0<=2e5 GC/B_radius^2
from current-geodesic-two-sided-remainder.md. The conversion of either
canonical integral divided by N is 2G/(B_radius^2 beta); original period
and number of identical halves cancel explicitly, not by changing norms.
The complete actual form and its comparison, rather than a local quantum
model, give (4). Since Delta T=Q T''+Q'T', put

    D2=[k2^(-1/2)+q1max sqrt(g1/qmin)]/qmin,
    ||f_h||X*^2 <= h^4 Phi2 D2^2.                       (5)

Here qmin<=inf_I Q and q1max>=sup_I |Q'|. This argument extends from
the actual core by (4) and weak differentiation. No smooth-even H4
assumption or nonlinear-neighborhood theorem enters. The compact source
annihilates endpoint boundary terms; the test functions still range over
the ENTIRE mass/log domain with its original York and bolt returns.

## Exact constant-flow contribution, with finite returns retained

Let Pi0 be the orthogonal projection onto ker P. Spectral positivity
in the COMPLETE formula (1) gives r(P)>=Pi0. This is a lower-bound
test of the proposed inverse, not a replacement of the quantum action
by a selected set of modes. The established actual orbit decomposition
gives, on crossing orbits with outer turning point b,

    j^2=Q(b),
    T_per(b)=2 integral_0^b dz/sqrt(Q(z)-Q(b)),
    (f_h)_(j,0)=2/T_per(b) integral_0^b
                         f_h(z)/sqrt(Q(z)-Q(b)) dz,
    <f_h,C_flow f_h> >= (1/(2pi)) integral
                              2T_per(j)|(f_h)_(j,0)|^2 dj. (6)

Both physical halves/directions have the original orbit normalization.
For the lower bound restrict ONLY the positive integral in (6) to
b=z0+h y with -1/2<=y<=-2/5. Such trajectories cross the actual seam,
have j^2<beta and turn before the bolt. Other returns and all trapped
orbits remain in C_flow; dropping their nonnegative contributions is
valid for this lower bound alone.

On -1<=x<=-2/5, phi''(x)=16(1-x^2)^6(15x^2-1)>=0. On
-3/4<=x<=-3/5 it is at least Lphi=16(7/16)^6(22/5).
Also 0<y-x<=7/20 there. Thus, defining

    Aphi=(3/20)Lphi/sqrt(7/20),

the exact source integral in (6) is >= Aphi/sqrt(q1max).
No oscillatory asymptotic or unknown onset frequency is needed.
The whole-interval mean-value theorem bounds Q(z)-Q(b) from above
by q1max(b-z) on this support.

For all b in I, certify -Q'(b)>=q1min>0 and Q(b)<=qmax. To upper-bound
the FULL period, use Q(z)-Q(b)>=(b-z) Dmin for 0<=z<=b, obtained by
factoring the actual polynomial before estimates. Writing Q=sum c_n z^n,
c3=-a3<0 and c4=b4>0, one sufficient coefficient bound is

    Dmin=(a3-(4/3)b4*.51)*.49^2-2|c2|*.51
                     -sum_(n=5,6) n|c_n|*.51^(n-1).

This uses (z^3+z^2 b+zb^2+b^3)/(z^2+zb+b^2)<=4b/3
for 0<=z<=b. The small initial increase of Q above beta is retained
in the c2 debit. Execute Dmin>0 and set Tmax=4 sqrt(.51/Dmin).
There is no wall at z=0, artificial short return, or sampled period.

The Jacobian dj/db=-Q'(b)/(2sqrt(Q(b))) then gives

    <f_h,C_flow f_h> >= L0 h,
    L0=Aphi^2 q1min/(5pi q1max sqrt(qmax) Tmax)>0.        (7)

The same bound holds for chi_I C_flow chi_I since chi_I f_h=f_h.
Combining (5),(7), EVERY bounded X*->X extension would have to satisfy

    ||C_flow|| >= R(h)=L0/(Phi2 D2^2 h^3).              (8)

The right side diverges as h->0. This proves nonexistence of the
required bounded extension, on the actual paired domain. Any bounded
coarse/endpoint correction C_reg of norm M changes the lower expression
by at most M. No finite M or scalar normalization removes divergence.

For the unmodified candidate, a hypothetical finite full preconditioner
defect would satisfy

    ||I-A_X C_flow||_(X*->X*) >= (1-e_Q) R(h)-1.        (9)

This follows from A_X>= (1-e_Q) R_X and the triangle inequality; B is
the FULL scalar Schur term. It is not obtained by omitting a coupled row.
Execute (5)-(9) at h=1e-16, after freezing the displayed inequalities.
No numerical application of r(P) is needed to prove this lower bound:
its exact invariant subspace and a positive finite-domain integral suffice.

## Quantitative role and one specified alternative

This test addresses the proposed UNIFORM complementary preconditioner,
whose conditional defect target was .4. It does not give a counterexample
on the structured B range. In particular, (8) must NOT be read as a lower
bound on ||S_Q-S_hat||_(Et->Et), the full residual, or physical output.
The outer-input omission remains <6.116646e-36, its full residual debit
<8.202157e-36, and the remaining sufficient S_inner application ceiling
4.9993883354e-32 remains an UNKNOWN actual error, not a physical margin.
The K error, signed G_t, mass/phase/graph/output budgets remain unknown.

The calculated loss is at least three powers of localization scale:
the inverse-compressed flow has quadratic response >=L0 h while the
original source dual norm squared is <=Phi2 D2^2 h^4. An acceptable
uniform inverse must suppress that response by at least R(h)/M for
its intended norm M. A finite coarse correction cannot supply this
regularity gain. This is stronger evidence than merely noting that
compression and inversion do not commute.

The specified alternative is to construct an inverse AFTER angular
compression, retaining its elliptic logarithmic principal part and
the global boundary reconstruction. On compact interior fields its
transported shape symbol has the already fixed leading expression

    a_shape(z,p) = G C Q(z)^2 p^4/(15pi B_radius^2 beta)
          *[log(|p| z sqrt(C Q(z)))-log(mu0 B_radius)-6/5]. (10)

To obtain (10), angular-average the full A_flow UV expansion, THEN
add the original Weyl coefficient log z-1/6. The angular finite part
is -31/30, so the total is -6/5. The scalar/R2 fourth part cancels on
the exact York principal relation nu=-T/3; B has order zero. All finite
contacts and lower terms still have to be assembled. Equation (10)
is a principal-symbol specification, not a scalar-Laplacian substitution
or a certified global inverse. It has inverse order p^-4/log|p|;
thus it addresses the three-derivative loss in (8). Its quantitative
remainder, transitions, paired mass/log endpoints, bolt return and
numerical application errors are NOT supplied by symbol matching.

Compute the exact interior coefficient and finite part at z=1/2 as
an implementation target, retaining the current mu0. Do not compare
this intermediate coefficient directly with a physical output ceiling.
For a specified comparison against the SAME localized source, freeze
only this principal part at z0, writing its positive coefficient as kappa
and finite part as b0. On the line define the Fourier multiplier

    C_pr(p)=[1+kappa p^4(b0+log(1+p^2)/2)]^(-1).

This is an explicitly labeled FROZEN PRINCIPAL comparison, not the actual
paired inverse or a substituted physical domain. Its complete frequency
integral, without a Fourier cutoff, satisfies by Plancherel and (2)

    <f_h,C_pr f_h> <= h^4 Phi2/(kappa b0).            (11)

Indeed the Fourier transform of f_h is -h^(5/2) p^2 phi_hat(hp).
The denominator is >=kappa b0 p^4. Compute the ratio of the attained
actual-flow LOWER bound (7) to the model UPPER bound (11). It quantifies
how this alternative removes the h^-3 loss on the test source. It does
NOT bound C_pr in the actual X norm: an upper bound for ||f_h||X* cannot
be used as a lower bound, and the full reconstruction/remainder is absent.
Thus no conditional operator ceiling is promoted to a certified constant.

The decisive next calculation is the full weak defect of such an
angular-compressed parametrix, including its global reconstruction,
against the existing .4 complement screen. That calculation remains
PROPOSED. No component/degree/tail optimization is selected.

## Executed result and research decision

The initial and alternative receipts retain the frozen arguments before
these reported results. The final 512-bit receipt identifies the final
sources. All 21 inherited/closed paths remain unchanged. Whole-interval
polynomial enclosures, exact bump integration and the ordinary positive
orbit argument above give, at h=1e-16,

    T_per < 10.199916 on the entire selected family,
    <f_h,C_flow f_h> > 2.7132441279e-20,
    ||f_h||X*^2 < 1.5424981260e-29,
    ||C_flow|| > 1.7589934680e9 if a bounded extension existed,
    ||I-A_X C_flow|| > 1.7552478997e9 if defined in that norm.

In fact (8) diverges as h^-3, so NO bounded extension exists. This rejects
the specified uniform inverse-compressed-flow estimator even with exact
quantum application and exact arithmetic. At this one explicit scale its
response needs suppression by more than 1.2537485005e9 to fit the necessary
inverse-norm ceiling associated with defect .4. The suppression required
over the whole domain is unbounded. It cannot be repaired by choosing
a smaller nonzero normalization or more digits, or by an energy-bounded
coarse/boundary correction. A different principal action is required.

The compress-first comparison has kappa in (9.7949972099e-36,
9.7949972100e-36) and b0 in (9.3323891833,9.3323891834). Its FULL Fourier
integral on the same source obeys (11) <4.7374670160e-31, more than
5.7272042611e10 below the attained flow-inverse lower bound. This is an
executed comparison of the required regularity mechanism, not a global
actual-domain application error. It removes this witness's h^-3 loss;
no finite coupled defect or source budget follows without the remainder
and reconstruction estimate. In particular, the ratio of two upper bounds
is not used to claim boundedness in X.

Decision: REJECT direct flow inversion followed by angular compression
as the uniform preconditioner. Retain the original energy-dual residual
method and construction; develop an angular-compressed logarithmic
parametrix instead. This is a verified failure of the specified estimator,
not failure of the structured B correction or of the wormhole construction.
It changes the preconditioner choice, not feasibility status: no actual
S_inner/K error, complete residual, signed G_t or admission margin is
obtained. No unknown source, graph or boundary comparison is replaced
by a conditional ceiling. The limited two-fast-input quantum screen
is unchanged and supplies no independent-phase or neighborhood estimate.

One next calculation: bound the complete weak defect of the compress-first
parametrix on the actual paired range, including endpoint reconstruction
and quantum/numerical errors, against .4. It remains PROPOSED. The first
missing implication for that alternative is a quantitative bound for its
non-principal and global-return operator; the symbol (10) alone is not it.
The full goal remains ACTIVE/INCOMPLETE.

Remaining to-do list: executed complete inverse and S/K application errors;
full residual and signed G_t; original admission, mass, independent-phase
and output budgets; nonlinear existence/contraction, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
