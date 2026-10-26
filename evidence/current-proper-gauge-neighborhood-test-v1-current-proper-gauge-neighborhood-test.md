# A quantitative limit on the proper-coordinate output comparison

2026-10-26 HST, E003/S169. Ordinary argument before frozen calculation.
Previous goal turn: PROGRESS (original norm recovery, exact structured
source cancellation and quantitative generic-estimator rejection).
Preserve the closed audit, accepted seed results, P32 and full objective.

## 1. The comparison being attempted, and what is already known

The direct complete return still requires the map from affine canonical
fields to the ORIGINAL outer proper coefficient. Its uniform nonlinear
extension cannot be obtained by differentiating a moving proper gauge on
only a C5 ball. stationary-centered-canonical-response.md already warns
that such a tangent can require a sixth derivative and instead fixes each
endpoint identification before interpolating. finite-outer-admission-and-
response.md likewise keeps ONE fixed coefficient manifold and radial log
zero in its proper collar. We retain those results. This note EXECUTES
the missing quantitative test of the moving-gauge shortcut; it does not
claim to discover that general composition loses derivatives.

The complete fixed-seed polynomial output is still
 T3=V_Y F^3 u,
 nu3=-A^-1 C T3+(C*)^-1(V_Y*)^-1 Bop F^2 u.
Its exact H source is (-A(C*)^-1(V_Y*)^-1 Bop F^2 u,0), with c4=0.
The actual paired first-order projectors, mass/log phase, common regulator,
full angular/thermal sums, exact ultraviolet logarithm, coincidence and
density contacts, cutoff transitions and regular-bolt returns remain in
these operators. No new quantum coefficient is supplied here.

The proposed shortcut would first bound T3 in the original five-jet
geometric norms, then transfer first/second variations uniformly through
the proper-gauge map. The counterexample below rejects THAT implication.
It does not reject finite endpoint admission, the original fixed-endpoint
proof, or the smoother restricted polynomial maps. We check this step
before numerical interpretation of a purported nonlinear budget.

## 2. Exact map on the actual outer collar

Let Y be the actual seed proper distance from the paired boundary in A
units. In fixed seed coordinates
 gbar0=dY^2+f0(Y)^(-2) dtau^2, l0=log f0,
 gbar_T=e^(2T(Y)) gbar0.
f0 is INVERSE warping, so the sign of T below is minus. Define
 X_T(Y)=int_0^Y exp(T(t))dt, Y_T=X_T^-1,
 B[T](X)=-T(Y_T(X))+l0(Y_T(X))-l0(X).                (1)
This is exactly log(f_T/f0) in the proper metric, with the SAME thermal
coordinate and period. No rescaling of the source/admission units occurs.
The identification can be patched to the identity between Y=1.1 and 1.35;
our fields are supported before .1. The original bolts and clock are
then unchanged. Equation (1) is only used where it is the actual proper
identification. Its inverse exists for bounded T by strict monotonicity.

For J_v(Y)=int_0^Y e^T v and h_T=e^-T, implicit differentiation gives
 Y_v=-h_T J_v,
 DB[T]v=[-v-h_T(l0'-T')J_v] composed with Y_T.        (2)
The first variation in C5 contains T^(6), and the second contains T^(7).
This is a statement about the comparison, not a new physical demand for
six or seven areal derivatives. All derivatives here are REGULAR proper
collar derivatives, away from the paired seam and the bolt.

For compact directions v supported strictly before an observation patch,
where T=0 on supp(v), set J=int v and J2=int v^2 over that earlier support.
There v vanishes near every observation. Put W=l0-T, L_T=e^-T partial_Y.
Differentiating the EXACT implicit coordinate change, for the affine path
T+s v, gives on the observation patch
 partial_X^5 DB[T]v=-J L_T^6 W,
 partial_X^5 D2B[T][v,v]=J^2 L_T^7 W-J2 L_T^6 W.     (3)
The J2 contact in (3) is retained. It is not permissible to approximate
the nonlinear coordinate shift by its first-order translation alone.

## 3. Actual-domain fields and finite quantitative test

Use the actual seed Q=(1-z)^2[beta(1+2z+3z^2)+u z^2+eps z^4]-tau z^3,
z(0)=0, dz/dY=sqrt(Q). Pick the observation z*=1/40 and its EXACT proper
coordinate Y*. The original parameter enclosures imply .06<Y*<.063:
Q<=.17 and, on z<=.025, Q>=beta[(1-z)^2-(1-L)^2]>.15.
No sampled endpoint or different seed is used.

Take a C11 cutoff chi supported in [.045,.080], equal to one on
[.050,.075], with beta(12,12) steps on its .005 ramps. Put
 T_N(Y)=a N^-5 chi(Y) cos(N(Y-Y*))
for the first test, and replace cos by sin for the second test. Set
a=1e-20, N=1e200. These are real fields on both original physical halves,
identical under the paired identification. nu=0, so the physical sphere
radius is unchanged. They vanish near the seam, all mass/log coefficients,
the original .05/.1 z graph transition, .1/.2 mass transition, cold caps
and the regular bolt. Smooth approximants preserve the stated estimates.
No artificial boundary, selected physical spectral modes or modified
quantum state is introduced. Each field has its complete two original
equation sources H(0,T_N); those sources are NOT evaluated or assigned
zero. These fields are domain witnesses, not candidate solutions or
claims about semiclassical validity of high-frequency perturbations.

Let p(t)=4^12 t^12(1-t)^12 on [0,1], zero otherwise, and
 v(Y)=1e-12 p((Y-.01)/.01)-1e-12 p((Y-.085)/.01).
The second bump is after the observation and keeps int v=0. All endpoint
jets through eleven vanish; the mass/log and regular-bolt constraints
remain unchanged. At the observation
 J=1e-14 4^12 (12!)^2/25!,
 J2=1e-26 16^12 (24!)^2/49!.                         (4)
Both are exact rational moments. For affine s the total proper length
can change at second order, but the physical bolt and its period are
unchanged; the fixed exterior patch restores the original identification.

For an explicit INPUT test norm on these compact collar fields use
 U5(v)=max_(j<=5) sup |v^(j)|/[j!100^j Y^(4-j)(1+log(1/Y))],
together with global C1/100. This is the original outer derivative
convention applied to the canonical input T, NOT a completed global X.
Exact Bernstein coefficients bound all derivatives of p and both chi
ramps on their full intervals. They give finite U5(v), independent of N,
and uniformly small U5(T_N). Every additional ORIGINAL fixed-chart
geometric seminorm is also uniformly bounded: its transitions are fixed
and the support is compact away from the endpoints. We do not assign a
numerical global atlas constant absent a comparison proof.

Here is a deliberately loose but proved seed-jet bound for executing
(3) without treating the seed as constant. On the actual proper collar
Y<=2.01, Q<=(1-z)^2 implies z<1-exp(-2.01)<.87. The lower polynomial
barrier gives Q>.002. Coefficient sums give |partial_z^j Q|<1e6,
j<=6. On a complex z disk of radius 1e-10 its variation is <.000101;
the positive square-root branch remains holomorphic and has modulus <1.
On the half disk its derivative is <2e10 by the elementary Cauchy bound.
Picard iteration for dz/dY=sqrt(Q), on a Y disk of radius 1e-12,
is a contraction (<.02) into that half disk. The nonconstant part of
l0=-log(Q)/2 has modulus <10 there. Cauchy's integral formula therefore
gives |l0^(j)|<10 j!10^(12j)<1e100 for j<=7. Constant clock factors cancel
in (1). This is an actual seed bound, not a globally analytic variation
hypothesis or a smooth reflection of the mass cubic.

Evaluate L_T^6 W and L_T^7 W by formal jets through seven, retaining
exp(-T) and ALL product terms. Treat every seed derivative l0^(j), j>=1,
as an interval [-1e100,1e100]. The finite result is consequently uniform
over the actual seed enclosure. Freeze tests that the original weighted
output norms of DB and D2B exceed 1e100 and 1e250 respectively. At the
actual X*=X_T(Y*), .04<X*<.09, so the fifth-jet denominator is <2e14.
Report the enclosed ratios to U5(v) and U5(v)^2 as well as raw outputs.

The backgrounds themselves remain in the outer admission. Finite chain
and inverse rules through FIVE give |partial_X^j B[T_N]|<=3a for j<=5
on X<=2; a sufficient explicit error allowance is 1e120 a/N. To see its
size, T_N^(j)<=2a N^(j-5), j<=7 after the certified ramp factors; every
nonleading inverse-coordinate term through order five contains a jet
through T_N^(4), hence <=2a/N, or the coordinate displacement <=a N^-5.
There are fewer than 10^10 product-rule terms through seven; seed jets
cost <1e100. The coefficient and exponential factors fit 1e120. The
direct fifth T jet costs <=2a. Where the field vanishes before .045 the
map is exactly zero; after the support its only tail is the small proper
translation of l0. The original cap slope <=4 and fixed patch fit the
same allowance. Thus Nout(B[T_N])<2e-14<1.5e-7. This is not merely a
sequence of inadmissible backgrounds with large fifth derivatives.

As N increases, (3) has leading magnitudes J a N for the cosine test and
J^2 a N^2 for the sine test; the lower terms grow strictly more slowly.
There is NO uniform C5-to-C5 first or second derivative bound for (1)
on the actual admitted neighborhood. The moving-gauge shortcut cannot
supply finite uniform nonlinear comparison constants from those five
jets alone. The finite witness and this asymptotic prove estimator
failure, not a failed gravitational inverse or wormhole obstruction.

## 4. Restricted polynomial maps, budgets and a specific alternative

This does NOT invalidate R_S V F^3. The accepted fixed-seed F^2 domain
already has ten local derivatives and its actual mass/log phases. It
has sufficient qualitative reserve for (2)-(3) on compact collars.
The corresponding M6,M7 quantitative output maps are UNKNOWN, as are
their weighted end counterparts. The fixed-seed output VALUE only needs
the original five geometric derivatives. Do not elevate the moving-gauge
method's six/seven requirements into necessary geometric admission.

For orientation, a stronger LOCAL comparison space can be specified
without concealing constants. Define H^(8,log) by the multiplier
<xi>^8 log(e+|xi|^2) on a regular two-dimensional chart. Then
 sup|D^j u|_F <=1/sqrt(pi) ||u||H^(8,log), j<=7.       (5)
Indeed the squared Fourier Cauchy constant is at most
(4pi)^-1 int_1^infinity dt/[t log(e+t-1)^2]<=1/pi,
using log(e+t-1)>=(1+log t)/2 and int_1^infinity
dt/[t(1+log t)^2]=1. This supplies a local kinematic repair only. It
does not complete X, its outer phase/weight embeddings, quantum C3,
common-regulator variations or nonlinear action-coordinate contacts.
The native strong source needs its exact UV-log reserve independently.

Consequential budget statement: the proved local radius multiplier and
nonlinear radius remainder, structured c4=0, form defect and limited
two-fast screen remain unchanged. Original-source defect, graph/physical
outputs, actual mass column and all application errors remain unknown.
In particular the finite .5/(h3)=1.09652849856e32 defect-product ceiling
and graph/physical ceilings do not control the moving-gauge derivatives.
If an admission/contraction chain inserts a uniform C5-only constant for
(2) or (3), that estimator has an INFINITE coefficient, not just a large
finite bound. The executed witness exceeds even a proposed finite 1e100
first-comparison allowance after its explicit input normalization. This
does not make the actual alpha_embed or quantum C3 infinite: their full
strong-space definitions and comparisons remain unproved. Higher degree
or smaller bolt-tail errors cannot repair this neighborhood implication.

The specified alternative is the ORIGINAL fixed-endpoint/fixed-coordinate
method: use the full two diagonal metric coefficients on one fixed
manifold, or prove stronger graph-space embeddings before differentiating
a proper-gauge map. The actual first-order expression with coefficients
a,b is already available in finite-outer-admission-and-response.md (3).
Small bounded C1 perturbations preserve its paired domain; obtaining the
full strong quantum source and metric-variation bounds for both coefficient
directions is still a continuum task. One may not replace it with the
proper-collar scalar theorem, whose radial coefficient is fixed to one.
Any change of nonlinear action variables also retains the gradient times
chart-acceleration contact, because the seed is not a full stationary point.

DECISION: retain P32 and the construction; retain the original fixed-chart
method for nonlinear admission. Reject an unproved moving-proper-gauge C5
transfer. The direct complete polynomial map remains UNRESOLVED; the new
evidence specifies a false shortcut and a regularity requirement of that
shortcut, not a certificate for the full map.

ONE NEXT CALCULATION (PROPOSED): derive a uniform complete outer strong
quantum response on the fixed two-coefficient chart, retaining its paired
mass/log realization and common contacts, and bind it to the original
graph neighborhood. This is the first required bridge for transferring
the canonical polynomial estimate without the failed gauge derivative.
Its numerical constant is NOT computed by the present kinematic test.

The prepared regular-bolt sector is a stationary initial bridge; the
existing plan already assigns opening to a separate drive. That scope
check changes no result and is not a new obstruction. Nonlinear existence,
exact compatibility, matching, semiclassical validity/noise, all-sector
stability and the subsequent opening/transport requirements remain open.

Verification: freeze proof and checker, execute exact coordinate identities,
whole-interval jet and actual-parameter guards, finite interval witnesses
and the embedding integral; ordinary review then unchanged acceptance.
No numerical quantum source or actual polynomial correction is claimed.

Remaining to-do list: fixed-chart strong quantum/neighborhood bridge;
complete structured original-output/source maps and all coupled entries/
application errors; nonlinear admission/contraction; exact compatibility;
physical matching; semiclassical validity/noise; all-sector stability;
opening, transport, preparation, separate creation and the full objective.
