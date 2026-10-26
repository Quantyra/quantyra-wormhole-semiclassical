# Actual-domain ultraviolet test of the retained quantum evaluator

2026-10-26 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: PROGRESS (failure of the global radial source-norm
comparison, and a finite mixed-chart replacement). The full stable
semiclassical wormhole objective remains active; the audit stays closed.

OUTCOME TO BE TESTED: the complete gravitational estimate remains
UNRESOLVED. Test a specific failed interchange in a proposed retained
quantum evaluator, and a replacement reference inverse on the ACTUAL
domain. The new parametrix below is for the massive MATTER operator A,
not for the coupled gravitational S, and is not P32. Its defect must not
be reported as ||I-SP32||. The already closed transport-to-bolt slice is
not a reason to optimize another tail. The purpose is to decide whether
the retained quantum center can be assembled term by term as proposed.

## 1. Domain and the failed interchange

Use the same finite seed Q, paired complementary realization, exact
mass/log fields, original temperature, and C2,1 identical-field double
as current-global-dirac-transport-series.md. Its two endpoints are the
original regular bolts. There is no smoothing, artificial physical wall,
or change of the coupled outer cutoffs. All angular levels n>=1 and
thermal labels j in Z remain, with M_n^2=n(n+q)/C, d_n=2(q+2n),
q=2000^4 and p=gamma(j+1/2). M>18. The massless sector, Weyl conversion,
exact ultraviolet logarithm, local/coincidence contacts and finite heat
restoration remain in their original allocations under ONE regulator.
Independent odd-mouth variations and the nonlinear neighborhood are not
covered by this seed representation.

In the fixed half-density frame, A is orthogonally equivalent to
 L=partial_y-M sigma3-(p/sqrt(Q))sigma1,
and, with the actual rotation O, to B0+alpha J. Here
 k=sqrt(M^2+p^2/Q), alpha=(atan[p/(M sqrt(Q))])'/2,
 B0=partial_y-k sigma3, J=[[0,-1],[1,0]].
Its whole-interval inverse G0 and actual regular endpoint domains are
those already proved. Write S0=G0 alpha J on Xk=L2(k dy).
The unitarily represented S0 is Hilbert--Schmidt for each fixed mode.
Indeed its kernel near either bolt has weights proportional to
(s/y)^(1/2)(y/s)^|omega| on its ordered triangle, bounded at the
borderline |omega|=1/2; alpha is bounded for that fixed mode. Away from
the endpoints all coefficients are bounded. Thus S0^2 is trace class,
and the product-kernel formula gives

 Tr S0^2=2 int_(y<s) alpha(y)alpha(s) exp[-2 int_y^s k] dy ds.       (1)

The weights cancel. Absolute integrability follows also from
int|alpha|dy<8. This bound follows by at most five monotone pieces of Q
on each half, each contributing at most pi/4 to half the theta variation.
It justifies (1), but not summation over all modes after regulator removal.

For the entire infinite angular subsequence M>=100 take the allowed
p=gamma/2. At either actual bolt let
 w=sqrt(2(z_h-z)/gamma), x=Mw, Q=gamma^2 w^2 R(x).
Taylor's theorem for the actual polynomial, |Q''|<2.1, gives on 1<=x<=4
 |R-1|<=2.1*16/(8*100^2)=.00042,
 ||Q'|/gamma-2|<=2.1*16/(2*100^2)=.00168.
Consequently
 1/2 < k |dy/dx| = sqrt(1+1/(4x^2 R))/sqrt(R) < 2,
 |alpha| |dy/dx|=|Q'|/[8 gamma sqrt(R)(x^2 R+1/4)] >1/100.
One ordered pair of subintervals [1,2] and [3,4] contributes at least
2 exp(-12)/10000 to (1), positively.

This is a bound on the FULL finite-domain trace, not only a local word.
On both central pieces z<=1/2, Q>.04 and |Q'|<.3 imply
 mu=int_central |alpha|dy <=10p/M<=gamma/20.
All potentially negative pairs touching these pieces cost at most 16mu.
On each remaining half alpha has one sign, because Q'<0 for z>=1/2.
Opposite halves have proper separation >=1/sqrt(.17), since Q<=.17;
their entire negative contribution is <=128 exp(-200/sqrt(.17)). Thus

 Tr S0^2 >= 2 exp(-12)/10000-16 gamma/20
                         -128 exp(-200/sqrt(.17)) > 1e-9.        (2)

The polynomial bounds used here are verified from the enclosed actual
coefficients: |Q'|<=16 beta/9+u/2+4|eps|+3tau<.3;
|Q''|<=12beta+2u+82|eps|+6tau<2.1; Q<=beta+u+|eps|<.17;
Q>=beta/2-|eps|-tau>.04 on [0,1/2]. In Q' the u term and -3tau z^2
are nonpositive on [1/2,z_h]. The eps term has absolute value at most
4|eps|z^2(1-z), strictly dominated by 12beta z^2(1-z). Hence the sign
above does not come from a sampled mesh.

Since sum_(n=n0)^N d_n=2q(N-n0+1)+2[N(N+1)-(n0-1)n0], (2) proves
 sum_(n,j) d_n |Tr S0^2| = infinity.
This rejects an ordinary separately unregulated wordwise trace sum.
It does NOT prove divergence of the commonly renormalized action, its
source or Hessian, or failure of the accepted adaptive tail. A single
thermal branch is an infinite-subseries witness, not a physical mode
truncation. No action projection is inferred from its sign.

There is also an operator obstruction to adding ultraviolet decay to
THIS reference by improving a constant. Take one input component equal
to one on x in [3,4]. Its squared Xk norm is <=2; the opposite output
component on [1,2] is >=exp(-6)/100. Its squared norm there is at least
exp(-12)/20000. Smooth approximants give the same operator lower bound:

 ||S0||Xk >= exp(-6)/200 > 1.2e-5, for every M>=100, p=gamma/2.    (3)

No A M^-r bound with finite A and r>0 can hold for this S0. A smaller
upper constant alone cannot fix that proposed ultraviolet gain.

## 2. One replacement: retain the regular-bolt Dirac block exactly

Use ACTUAL proper distance rho from each bolt, so the principal derivative
is unchanged. Set a(rho)=sqrt(Q)/gamma. Then a(0)=0, a'(0)=1 and
a''=(Q''/2)a. The Volterra integral equation and |Q''|<=2.1 give
 |a-rho|<=sinh(sqrt(1.05)rho)/sqrt(1.05)-rho.
On 0<=rho<=.2 this proves .99rho<=a<=1.01rho and
 |1-rho/a|<=.22rho^2. The power series has positive coefficients, so
checking its endpoint bounds controls the entire interval.

Let d=M+|p|>=100 and r=d^-1/2. On each bolt patch rho<=2r use
 L_f=partial_rho-M sigma3-(omega/rho)sigma1, omega=p/gamma,
with the ORIGINAL regular branch at zero and the decaying branch at
infinity for this auxiliary inverse P_f. Embed compactly supported patch
sources in that half-line and cut the output back to the physical patch.
This is an auxiliary reference, not a replacement physical domain.
Its defect below restores the entire true geometry and finite returns.

L_f is orthogonally equivalent to the flat polar D+M, with the regular
disk realization. Thus ||P_f||<=1/M. Squaring gives scalar potentials
M^2+(omega^2 +/- omega)/rho^2. For |omega|>=3/2 their quadratic forms
give ||(omega/rho)P_f||<=sqrt(3). For |omega|=1/2 use the gap bound
directly: do not apply a false Hardy inequality to the critical branch.
Since L-L_f is the zero-order multiplication
 V_f=omega(1/rho-1/a)sigma1,
the localized defect is <=.88sqrt(3)/d<1.53/d for |omega|>=3/2.
For |omega|=1/2 it is <=.22/(M sqrt(d))<1.53/d as well.

Outside the two patches rho<r use the existing WHOLE-interval G0,
converted back by its actual left/right orthogonal frames. Monotonicity
of Q on the deep half and Q>.04 centrally give a>=.99r there.
The accepted Q'^2<=4gamma^2+4.2Q implies
 |alpha|<=|Q'|/(8sqrt(Q))<=sqrt(d)/(4*.99)+sqrt(4.2)/8.
Moreover ||G0||L2<=1/kmin<=sqrt(2)/d, by its scalar kernels and
Q<=.17. Thus the exterior localized defect is <.4/sqrt(d) for d>=100.

Choose a squared partition eta_L^2+eta_R^2+eta_o^2=1. At each bolt
eta_b=1 through r, then cos(pi(rho-r)/(2r)) through 2r, then zero;
eta_o is the complementary sine, or one off both patches. These fixed
Lipschitz cutoffs preserve the first-order domains. They are not declared
high-order graph cutoffs; the original metric cutoffs are untouched.
The patches are disjoint. Their derivative row norm is <=pi sqrt(d)/2.

With each inverse and its left/right frames embedded in the SAME physical
spinor Hilbert space H=L2(dy;C2), define
 P=sum_i eta_i P_i eta_i,       E=I-A P.                           (4)
Here P_o is the actual-frame G0 and the bolt P_i are the mapped P_f.
P:H->Dom(A) uses the original paired transmission and regular bolts.
Products and closure follow first on compact sources, then in the graph
norm from A P=I-E. The singular difference V_f is bounded on each patch
for fixed mode. There is no new boundary distribution at the seam or at
patch edges. The actual paired projectors remain encoded by the double.

The localized inverse bound is <=2/d for every patch. For the bolt,
if |p|<=M use 1/M<=2/d; otherwise use
 ||1_patch P_f||<=2sqrt(3)/(sqrt(d)|omega|)<=1/|p|<=2/d.
The outer bound is sqrt(2)/d. Row/column Cauchy--Schwarz with the squared
partition, using localized outputs on the patch supports, proves

 ||P||<=2/d,
 E=-sum_i eta_i V_i P_i eta_i-sum_i [A,eta_i]P_i eta_i,
 ||E||<=.4/sqrt(d)+pi/sqrt(d)<4/sqrt(d)<=.4.                       (5)

The commutators, cutoff transitions and global returns have not been
discarded. The actual inverse is P sum_(k>=0) E^k, by uniqueness of A.
The sector d<100 retains the existing exact adaptive transport; every
mode belongs to one of these two sectors. This proves no gravitational
or physical-norm defect statement.

## 3. Quantitative role in the complete error budget

The following calculation tests whether (4) is usable in a retained
quantum evaluator without consuming the already closed bolt error slice.
It is not a further optimization of that slice and does not evaluate
the retained signed center.

For each mode, sum_i ||P_i eta_i||HS^2 is bounded by
 sqrt(2)*395/d+24/(M sqrt(d)) <58/sqrt(d).                         (6)
For G0 this follows directly by squaring its exponential kernel on the
actual interval of length <395. For a bolt, P_f^*P_f=(L_f L_f^*)^-1.
If |omega|>=3/2, each scalar potential is nonnegative, so its diagonal
Green function is <=1/(2M) by positive-potential comparison on the
half-line. For |omega|=1/2 the critical scalar diagonal is
rho I0(Mrho)K0(Mrho), with its regular (no log) endpoint condition;
the other is <=1/(2M). The exact Bessel Wronskian fixes the coefficient.
For all x>0, x I0(x)K0(x)<5: for x<=1 use I0<=e and
K0<=log(2/x)+1, from its positive cosh integral; for x>=1 use
cos(theta)<=1-2theta^2/pi^2 and cosh(t)>=1+t^2/2 in the two positive
integral formulas to get xI0K0<=pi/4. Thus the spin trace is <=6/M.
Each bolt source support has length 2/sqrt(d), giving (6).
These formulas and Wronskian are [DLMF 10.32.1](https://dlmf.nist.gov/10.32.E1),
[10.32.9](https://dlmf.nist.gov/10.32.E9) and
[10.28.2](https://dlmf.nist.gov/10.28.E2). They are used for positive
arguments with the stated regular realization, not a smooth-seam theorem.

The coefficient row eta_i V_i+[A,eta_i] has norm
 <=(.44/gamma+.3+pi/2)sqrt(d)<sqrt(d)/gamma.
Therefore ||P||HS||E||HS<=58/gamma, including all patch transitions.
Let R_N=P sum_(k=0)^(N-1)E^k. Its trace-class inverse remainder obeys
 ||R-R_N||1<=(58/gamma)delta^(N-1)/(1-delta), delta=4/sqrt(d).
In the frozen mass chart A_b=D+M exp(b), P stays fixed and
||E'[v]||<=2||v||infinity at b=0. Differentiating the remainder keeps one
Hilbert--Schmidt E and P and gives
 ||D(R-R_N)[v]||1<=(116/gamma)||v||infinity
                              sum_(k=N)^infinity k delta^(k-2).
The derivative of the exponential source vertex supplies an additional
R-R_N term. It is retained in the Hessian bound below.

The ENTIRE angular and half-integer thermal sums satisfy
 F5=gamma sum_n d_n M_n sum_j (M_n+|p_j|)^-5
                              <=C/3+4C gamma/81.
Indeed the complete thermal sum is <=1/(2M^4)+2gamma/M^5, and the
accepted count N(M)<=4CM^2, M>18 gives sum d_n/M_n^3<=2C/3 and
sum d_n/M_n^4<=2C/81. On d>=100,
delta^(N-1)<=.4^(N-11)4^10 d^-5. With the original conservative
doubled-trace action conversion G gamma/(pi B^2), G=1, define
 Kq=58G 4^10(C/3+4C gamma/81)/(pi B^2 gamma).
For N>=12 the complete HIGH-sector normalized source and Hessian tails are

 e1(N)=Kq .4^(N-11)/.6,
 e2(N)=Kq .4^(N-12)[.4/.6+2(N/.6+.4/.6^2)].                     (7)

These are seed mass-field L-infinity source/Hessian bounds. They include
both inverse and coincidence-vertex differentiation. All excluded low
sector errors are bounded by the EXISTING complete adaptive certificates
e1old=9.204e-79*.3961^128 and e2old=9.204e-77*.3961^128; this harmlessly
overbounds that subset. No finite set of physical modes is substituted.

Common regulator removal is used only for these trace-class seed jets.
At a positive COMMON cutoff the inverse identities and product rules
hold on the original composed domain. The relative form bounds for
A'_b=Mv and A_b^{-1} give uniformly bounded heat-cutoff derivatives,
strongly tending to zero on fixed vectors. Their products with the
trace-class seed tails tend to zero in trace norm, by finite-rank
approximation followed by the uniform bound. The majorants (7) justify
the complete mode sums. Thus the difference between exact and retained
source/Hessian has this common-cutoff limit. The retained pieces still
have to be combined with ALL original local heat subtractions and finite
contacts at the same cutoff. Equation (2) forbids deleting that step.
No trace-class family on a whole mass/log neighborhood is claimed:
E_b need not decay with d off seed. No strong spatial graph jets follow
from the Lipschitz partition. The nonlinear/domain bridge stays unproved.
The seed source and Hessian error are variational jets: differentiate
the commonly regulated tail sum_(k=N+1)^infinity Tr(E_b^k)/k before
removing the cutoff. Cyclicity of its trace-class seed jets gives
-Tr[(R-R_N) A'_b] and its derivative, including the contact. This does
not assert existence of that unregulated tail action away from the seed.

For total e_j=e_jold+e_j(N), the accepted COMPLETE two-row comparisons give
 epsH=e2 Cnew^2/beta, ef=e1 Cnew/beta,
 U(N)=Knew[ef+epsH*(1.1e-10)]/(1-epsH),                         (8)
where Cnew<2.592386349e26 and Knew<1.960636035e70 retain original units.
The frozen computation finds the first N with epsH<1 and
U(N)<8.335e-10, records the previous failing estimator, and both source
and Hessian debits. This is a mathematical inverse-comparison error;
quadrature, signed centers and numerical polynomial applications are
NOT included. Do not confuse N in (7) with gravitational polynomial degree.

## 4. What this can and cannot change

Equations (2)-(3) establish failure of a specified evaluator interchange
and of a proposed ultraviolet norm gain, not a construction obstruction.
Equations (4)-(8), if verified, supply a concrete reference with the missing
ultraviolet decay and a complete-spectrum error compatible with one
already controlled physical slice. The small matter defect is not an
original-source or physical-output comparison for the coupled system.

The complete P32 bolt estimator remains <9.216793497e25 against 1.667e-9,
an excessive specified upper estimator by 5.5289703037e34. Its dominant
known loss is the scalar trace Cs near 9.9165e36. Degree 63 passes only
the exact bolt-tail allowance. Original-source defect, graph/physical
maps, independent phases, numerical applications and every coupled entry
remain UNKNOWN. In particular kappa_m<=5e-5 is not certified; its
conditional 2.14765e-12 affine and .0024731 feedback debits are not actual
columns. No ceiling for Cdef3 IY or the graph/physical maps is assigned
to an unknown coefficient. The two-fast quantum screen .135973<.1599164
keeps its limited scope. Neither the 10^4000 C interior estimate nor
10^74739.22 absorption is invoked to obtain (8); their unknown strong
composition is not thereby resolved.

RESEARCH DECISION: retain the construction and P32. Reject separately
unregulated retained-word evaluation. Use the actual regular-bolt block
as a candidate ultraviolet reference, keeping the common regulator and
all low-order subtractions. Do not increase gravitational degree or
optimize scalar/cold/contact/tail components on this evidence.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): assemble the signed
common-restored retained quantum source and its Hessian application to
one actual phase-conforming coupled inverse trial, then enclose the
COMPLETE residual triple (R_o,c4,R_i) from current-bolt-source-domain-test.
All low-order ultraviolet cancellations, original local contacts and
finite returns must be combined before tracing. This simultaneously
tests the missing center and its original source data; (7) alone does not.

Remaining to-do list: the retained quantum center and complete coupled
residual; original source/graph/physical comparisons and every coupled
entry; nonlinear admission/contraction; exact compatibility; physical
matching; semiclassical validity/noise; all-sector stability; opening,
transport, preparation, separate creation and the full stable
semiclassical wormhole objective.
