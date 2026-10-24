# Weighted third comparison with every separated mixed tail

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
Full stable-wormhole objective ACTIVE/INCOMPLETE; feasibility audit CLOSED.

The specified seed comparison now has a continuous bilinear return
from TWO original normal-inverse inputs to the weighted H3 source row,
with bound C3_flat<1.295e10. At the proposed fast radius 1e-11 it costs
<.1295 of the unchanged .159916395999 budget. This closes the sharp-split
weighted transfer for this comparison, including every mixed tail. It
does not close the actual curved or uniform nonlinear response.

The first development bound was 2.2837863174e10 and exceeded the budget.
The executed improvement tightens a logarithmic norm inequality, with
identical physical fields, weights, cutoffs and action normalization.
Neither bound is a measured operator norm or a nonexistence theorem.

## 1. Precisely the component bounded

Retain x=log z, sigma=4.005, a=.01, chi with H3 multiplier allowance
100, and the SAME full B_sigma inverse. Put w_i=B_sigma^-1 g_i and

 u=1_(0<z<.02) z^sigma w_1(log z),
 v=theta(z) z^sigma w_2(log z).

theta is the already specified fast extension, one through .02, between
zero and one thereafter, and zero beyond .03. Both inputs are extended
evenly to the comparison line. For remote bounds use b=.03 for both
supports. This overestimates the first support; it does not move its
sharp split. The physical paired boundary realization is unchanged.

Let T=T_ge4 be EXACTLY the complete stationary comparison third kernel
in current-mixed-frequency-third-response.md: the entire angular tower,
three pair contacts, both triangle orientations, and the prescribed
fourth contact d4=C(1-2 kappa4)/(30pi). In its notation,

 T(p1,p2,p3)=sum A(p_i)+sum B(p_i)/(p_j p_k)
                         +d4 sum_(i<j) p_i^2 p_j^2,
 (pB)'=-p^3 A',
 |T|<=7 C P^4 log(2+P),  P=max|p_i|, sum p_i=0.       (1)

The response bounded here, in the original action units, is

 R3(g1,g2)=-(2G/B^2) chi z^-sigma z^5 partial_z T[u,v]. (2)

G=1 is the retained numerical unit. This is a seed flat comparison
component, not the nonlinear map on an arbitrary rough metric. T_low,
the actual curved difference, finite boundary completions, and remaining
conversion terms are separate parts of the EXACT allocation. No low
contact is set to zero. The already bounded thermal and relative pieces
remain retained and parked. No finite angular or physical-mode model is
used here. In particular, a whole-line H8 norm of the sharp input u is
not used: it would generally be infinite.

## 2. Infinite annuli, the original inverse and local inputs

Use I_j=[h_j,2h_j], h_j=(a/2)2^-j, for EVERY j>=0. The C8 cutoff eta
is supported in [.5,3], is one on [.75,2.5], and has beta(9,9) transitions
on [.5,.75] and [2.5,3]. Set u_n=eta(z/h_j)u, v_n=eta(z/h_j)v. All local
supports are below .015, where both original input restrictions are one.
There is no physical wall or smooth replacement of the boundary domain.

The exact decomposition on an observation annulus is

 T(u,v)=T(u_n,v_n)+T(u-u_n,v)+T(u_n,v-v_n).           (3)

In each remote term one factor f vanishes on [.75h_j,2.5h_j]. All its
local jets on I_j vanish. For remote y, |z-y|>=max(z,|y|)/5 and also
|z-y|>=z/8, including the even reflected side. Derivatives are taken
on the true output first and then estimated with this fixed split;
the choice of annulus is never differentiated.

Use the existing full-normal multiplier bounds, not sampled inverses:

 K=GC beta^2/(15pi B^2), L=10pi(Fl+23g/(180pi))/g,
 Cb=[2/(KL .005 sigma(sigma-3)),2/(KL sigma(sigma-3)),
                                      2/(KL sigma),2/(KL)],
 Cb4=max(Cb0,8/(KL),1/[K sqrt(L log6)]), Cq=2.328,
 W=[Cb0,Cb1,Cb2,Cb3,Cb1,Cb2,Cb3,Cb4,Cq/(K log2)].     (4)

W_m bounds the pointwise Fourier multiplier for D^m w relative to H3
data. In particular W_7<1.206e18 and W_8<9.725e35. The existing bound
K ||w||H8,log<=Cq ||g||H3 includes the eighth LOGARITHMIC derivative.
With P0=1e18, a frequency split gives

 WL_m=max(W_m log(2+P0),(Cq/K)P0^(m-8)), m<=7,
 WL_8=Cq/K.                                         (5)

These are pointwise multiplier bounds on the two frequency sets; the
maximum is therefore valid. They are not inferred by taking a maximum
of two unrelated norm estimates. No analyticity of the mass/log fields
is assumed. The top derivative remains in L2 throughout.

## 3. Logarithmic coordinate change with an explicit finite constant

Write ell(D)=log(2+|D|). For a whole-line L2 function define

 J(f)=int_(|h|<1) ||f(.+h)-f||2^2 log(1/|h|)/|h| dh,
 Q(f)=||f||2^2+J(f).

The Fourier multiplier of J is 4 int_0^1(1-cos(xi h))log(1/h)dh/h.
Splitting at 1/|xi| proves Q(f)<=8||ell(D)f||2^2. For the converse put
l=log|xi|. Integrating the cosine integral by parts on [1,|xi|] gives
J_multiplier>=2l^2-8l. If l>=17/4, this is >=2l^2/17 and
log(2+exp l)<=l+2exp(-l). The resulting ratio is <19. For |xi|<=exp(17/4),
log(2+|xi|)^2<19 directly. Consequently

 ||ell(D)f||2^2 <=19 Q(f) <=152||ell(D)f||2^2.        (6)

The loose development converse used 65 in place of 19. Sharpening this
PROVED inequality is the only material numerical refinement here; no
input, output, inverse, normalization or physical assumption changes.

Extend phi(t)=log t on [.5,3] affinely outside, with slopes 2 on the
left and 1/3 on the right. It is a C1 bi-Lipschitz bijection, phi' in
[1/3,2]. Changing both variables in J gives an inverse-Jacobian factor
9 and kernel bound 2 log(2/|x-y|)/|x-y| on |x-y|<2. Splitting the extra
log2 term at |x-y|=1/2 and the outer part at 1 gives

 Q(f composed phi)<=36 J(f)+(3+216(log2)^2)||f||2^2
                   <=107 Q(f).

For a Lipschitz multiplier psi, with M_i=||psi^(i)||infinity,
Q(psi f)<=2M_0^2 Q(f)+M_1^2||f||2^2. Combining these estimates directly
and using (6), rather than multiplying norm equivalences twice, proves

 ||ell(D)(f composed phi)||2 <=128||ell(D)f||2,
 ||ell(D)[psi(f composed phi)]||2
             <=181 M_0||ell(D)f||2+8 M_1||f||2.       (7)

The numerical guards follow from 19*107*8<128^2,
2*19*107*8<181^2, and 19*3<8^2. Translation in x costs nothing.

For the TOP derivative use a commutator, avoiding an unnecessary M_1
times W_8. The elementary Fourier inequality

 |ell(xi)-ell(eta)|<=|xi-eta|(1+|xi-eta|)/<eta>

and Young's inequality imply

 ||[ell(D),psi]F||2 <=Cpsi ||F||H^-1,
 Cpsi<=(1/pi)[(5/6)||psi||1+(3/2)||psi''''||1]
       <=2(M_0+M_4),                                (8)

when psi is supported on an interval of length 2.5. To obtain the first
constant split its Fourier integral at |xi|=1 and integrate by parts
four times in the tail. If F=D_x^8 w composed phi, then
F=(1/phi')partial_t(D_x^7 w composed phi). Since |1/phi'|<=3 and its
weak derivative has magnitude <=1,
||F||H^-1<=4sqrt(3)W_7||g||H3. The top bound is thus

 128 M_0 WL_8+8sqrt(3)(M_0+M_4)W_7.                  (9)

Scaling z=h t adds at most sqrt(3)log(1/h)M_0 W_m, since
log(2+|xi|/h)<=log(2+|xi|)+log(1/h), h<1. These estimates extend by
completion from smooth graph data. No highest-derivative L-infinity
bound and no smooth physical metric approximation is inferred.

## 4. Evaluate every local annulus

Let M(i,d)=sup_t |eta^(i)(t)t^(sigma+d)|. The beta CDF has degree-17
Bernstein controls [0 repeated 9, 1 repeated 9]. Differentiate these
controls and bisect by de Casteljau six times. On each of the 64 CLOSED
cells take the maximum absolute control, times the larger endpoint
value of the monotone power t^(sigma+d), with the transition scale
factor. Include the plateau for i=0. This encloses entire intervals;
it is not evaluation of a finite set of field modes or point samples.

Put F_n(X)=product_(r=0)^(n-1)(sigma+X-r)=sum f_nm X^m. Exactly,

 partial_z^n(eta(z/h)z^sigma w(log z))
 =h^(sigma-n) sum_i binom(n,i) eta^(i)(t)t^(sigma-n+i)
                                  F_(n-i)(D_x)w(log h+log t). (10)

The lower input Fourier L1 factor (including the unitary convolution
normalization) is <=sqrt(||v||2 ||v'||2). For one positive local input
its constant, after removing h^sigma, is V=sqrt(V0 V1), where

 V0=sqrt(3)M(0,0)W0,
 V1=sqrt(3)[M(1,0)W0+M(0,-1)(sigma W0+W1)].

Numerically V<231124. Even extension costs sqrt(2) here. For the high
logarithmic input use the triangle bound 2 for its reflected copy. With
both allocations in the bilinear convolution, the total factor is
4sqrt(2). All comparison images are thus retained.

For r=1,...,4, k=4+r, the symbol bound (1), convexity of
F(t)=<t>^k log(2+t), and F(2t)<=2^(k+1)F(t) give

 ||partial_z^r T[u_n,v_n]||2 <=7*2^k C
       [||u_n||Hk,log ||v_n_hat||1/sqrt(2pi)+symmetric]. (11)

F is increasing and convex: its second derivative is positive for
k>=5, as is seen by retaining the positive (k log2-1) term. Expanding
the Hk,log square gives derivative weights sqrt(binom(k,n)), n<=k.

Define c_00=1, c_(l+1),j=(5+j-sigma)c_lj+c_l,j-1, absent entries zero.
They are positive for 0<=j<=l<=3, and
(D-sigma)^l[z^5 partial_z T]=sum_j c_lj z^(5+j)partial_z^(1+j)T.
The output weight on I_j costs 2^(4+r-sigma) and the dz-to-dx conversion
h^-1/2. Together with (10)-(11), the surviving scale is h^p,
p=sigma+k-n>=sigma>0. Sum EVERY annulus by

 S_p=(a/2)^p/(1-2^-p),
 SL_p=S_p[log(2/a)+log2*2^-p/(1-2^-p)].              (12)

For each i,n set M0=M(i,-n+i),
M1=M(i+1,-n+i)+|sigma-n+i|M(i,-n+i-1).
For m<=7 let A_m=181M0 WL_m+8M1 W_m. The only m=8 term is n=8,i=0;
use (9) with
M4=sum_(j=0)^4 binom(4,j)|(sigma-8)_(4-j)| M(j,-12+j).
The complete near bound is the following finite expression:

 C3_near=(2G/B^2)100*4sqrt(2) C V
  sum_(l=0)^3 sqrt(binom(3,l)) sum_(j=0)^l c_lj
    2^(4+r-sigma) 7*2^k sum_(n=0)^k sqrt(binom(k,n))
      sum_(i=0)^n binom(n,i) sum_m |f_(n-i),m|
        [A_m S_p+sqrt(3)M0 W_m SL_p],
 r=1+j, k=4+r, p=sigma+k-n.                         (13)

This uses a conservative triangle sum across disjoint output annuli.
It requires neither summable input-localization overlap nor any claim
that a sharp input is globally H8. The entire eighth derivative is in
L2 with its logarithm. Its contribution dominates (13).

## 5. Complete remote terms without differentiating the input jump

Away from zero the kernel k_A equals the complete flat bubble kernel;
their difference is local. For the two-dimensional bubble,

 B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2),
 partial_M B_M=M^3[2K0(Mr)^2-K1(Mr)^2]/pi^2.

The representation [DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9) bounds
|K_nu(Mr)| by K_nu(M Re r) for nu=0,1 and Re r>0. With N(M)<=4CM^2
and the positive squared-Bessel moments
int s^5 K0(s)^2 ds=16/15, int s^5 K1(s)^2 ds=8/5,
Stieltjes integration by parts gives (224/(15pi^2)) C (Re r)^-6.
The moments follow by squaring the positive integral and performing
the gamma/beta integrals; all integrals are absolutely convergent.
For |zeta-s|<=s/4, Re sqrt(zeta^2+t^2)>=sqrt(s^2/2+t^2).
Integrating t over its ENTIRE line gives coefficient
(224/(15pi^2))(3pi/sqrt(2))<16. Cauchy therefore proves

 |k_A^(r)(s)|<=C A_r |s|^(-5-r), A_r=16 r!4^r.       (14)

The beta rescaling costs beta^2<=1. This analyticity concerns the
translation-invariant comparison kernel, not the rough actual metric.

The crucial mixed B part follows from (pB)'=-p^3 A'. If k_F is the
inverse Fourier kernel of (pB)''', then, away from zero,

 k_F(y)=3y^2 k_A''(y)+y^3 k_A'''(y),
 Kp(y)=k_F(y)/y^3,
 |Kp^(n)(y)|<=C H_n y^(-8-n),
 H_n=A_(n+3)+3 sum_(j=0)^n binom(n,j)j!A_(n-j+2).    (15)

H_0,...,H_3=[7680,118272,2300928,54088704]. The cube divided difference
in the preceding note has physical average
y^-3 int_(z-y)^z [int_s^(s+y)u][int_s^(s+y)v]ds. Its FIRST output
derivative consequently has the exact primitive formula

 int_(t>z) Kp(t-z) U(z,t)V(z,t)dt
             -int_(t<z) Kp(z-t) U(t,z)V(t,z)dt,
 U(s,t)=int_s^t u, V(s,t)=int_s^t v.                 (16)

Evenness of k_F and the original factor 1/2 cancel, so no factor is
missing. If one input f is remote, its primitive vanishes for |t-z|<z/8.
Further output derivatives through order four fall on Kp and LOCAL
jets of the other input only. They never differentiate the remote jump
at .02. All local polynomial contacts vanish for this remote factor;
they remain present in the local estimate (13).

For explicit constants use N0=1801/sqrt(2) and
N_j=N_(j-1)(sigma+2-j), j=1,...,4. The existing H6 inverse and Sobolev
evaluation give |v^(j)(z)|<=N_j z^(sigma-j)||g||H3 on the observation
set, and |u(y)|,|v(y)|<=N0 |y|^sigma||g||H3 globally on their supports.
For each r=1,...,4 bound one remote term by C times the following sum
of coefficient/power pairs (d,e), meaning d z^e:

* From partial^r[v A f], for j=0,...,r and l=r-j:
  l=0: d=binom(r,j)2 A0 5^5 N0 N_j b^(sigma-4)
             [1/(sigma+1)+1/(sigma-4)], e=sigma-j;
  l>=1: d=binom(r,j)2 A_l 5^(5+l)N0 N_j
             [1/(sigma+1)+1/(4+l-sigma)], e=2sigma-4-r.
* From partial^r A(fv):
  d=2 A_r 5^(5+r)N0^2 b^(2sigma-4-r)
             [1/(2sigma+1)+1/(2sigma-4-r)], e=0.
* From (16), with no derivative on a local jet:
  d=2 H_(r-1)N0^2 9^(4+r)b^(2sigma-4-r)
             [1/(2sigma-4-r)+1/(4+r)], e=0.
* From (16), j=1,...,r-1 derivatives on that jet:
  d=2 binom(r-1,j)H_(r-1-j)N0 N_(j-1)9^sigma
              8^(5+r-j-sigma)/(5+r-j-sigma), e=2sigma-4-r.

For the first two rows integrate y^sigma or y^(2sigma) separately below
and above z, using (14) and the factor-5 separation. For the last two
use |U|<=N0 s min(b^sigma,(9s)^sigma), s=|t-z|>=z/8; split the first
integral at b/9. The smallest positive denominator is 2sigma-8=.01.
Both t directions and even input images are included. Local f A(v)
vanishes. No unlisted mixed term remains in (1).

The weighted L2(dx) norm of z^q on (0,a) is a^q/sqrt(2q). Every
q=e+4+r-sigma above is positive. Thus the full remote bound is

 C3_far=(2G/B^2)100 C*2 sum_l sqrt(binom(3,l)) sum_j c_lj
                  sum_(d,e at r=1+j) d a^q/sqrt(2q). (17)

The final factor 2 is for the TWO remote bilinear terms in (3).
Equations (14)-(17) are bounds on actual derivatives of the complete
separated comparison kernel, not derivatives of a real envelope.

## 6. Executed result and what remains unproved

Interval evaluation at the unchanged correlated physical parameters gives

 C3_near=12946370060.2512027481...,
 C3_far=6.20567943831811042312...e-15,
 ||R3(g1,g2)||H3 <=C3_flat ||g1||H3 ||g2||H3,
 C3_flat<=C3_near+C3_far<1.295e10.                  (18)

The m=8 input term is 12946370060.2510767119...; m=7 is .00012603525
or less and all lower terms together are <1e-9. The earlier converse
65, composition 236 and product 334/14 gave 22837863173.8066102069... .
Both are upper bounds for the SAME map. Tightening (6) suffices to
change this comparison-only screening decision. Further improvement
of its tiny mixed tails, or the parked thermal/cold/scalar components,
cannot materially change the remaining budget.

Define eta_zero as the remaining actual seed nonlocal defect and
C3_add as a UNIFORM bound for the complete metric response minus this
seed comparison, including its background variation, low contacts,
conversion and actual boundary completion. A sufficient test is now

 eta_zero+1e-11 C3_add < .030416395999,               (19)

using the rounded allocation .1295 for (18). Even with eta_zero=0,
C3_add must be <3.0416395999e9. An equal split of (19) requires
eta_zero<.0152081979995 and C3_add<1.52081979995e9.
These are thresholds for UNKNOWN estimates, not achieved constants.
If one requires the seed comparison alone to consume at most half of
the original .159916395999 budget, the present guard fails that optional
allocation. It DOES satisfy the full remaining comparison-only budget.

The bound (18) defines a continuous comparison response on the actual
specified fast inputs using localization and separated primitives. It
agrees with the already identified finite part on smooth test inputs.
It does NOT prove convergence of the raw common-cutoff third response
in this graph. In particular no operator-norm cutoff limit is inferred
from distributional identification or density. That distinction is
essential after the earlier matched-normal cutoff obstruction.

FIRST UNSUPPORTED CONTINUUM STEP: bound the complete actual curved
zero-winding defect and T_low uniformly on the admitted background
GRAPH/mass-log family, with its original paired domain and a matched
strong common-regulator removal. The first-order domain and weak action
results are not this strong theorem. A-only/C5 background control was
already disproved; merely appealing to small low metric jets is invalid.
The original X/graph comparison and global graph admission are unknown.

The finite comparison estimate and the earlier complete mixed-frequency
identity show that a sharp analysis split and opposing frequencies do
not by themselves prevent this graph return. The remaining 0.0304
budget is positive but modest. Curvature and local contacts can change
the logarithmic principal coefficient, so treating them as an arbitrary
lower-order small error is unjustified. Their combined leading symbol
and exact contacts must be evaluated before declaring (19) attainable.
There is no certified actual upper bound to label excessive, and no
actual solution obstruction is proved. This supports a conditional
continuation of the full-graph method, not a claim of feasibility.

Retain M_A<=3.18e27, M_V<=2.66e30 and remaining global A<=.00094319 as
UNKNOWN contact/admission targets. At the original rho=1e-8 the exact
nonlinear conditions are
 Y0+1e-8 Z0+5e-17 Z2<=1e-8,  Z0+1e-8 Z2<1,
with unproved allocations Y0<=5e-9, Z0<=.25 and Z2<=5e7. No full
|P_hat_rest|+E_rest<1.287e-40 screen has been executed; it is sufficient,
not necessary for existence. The earlier selected 2.123e-42 enclosure
retains its limited scope and zero was not a measured full projection.
No actual controlled solution, exact compatibility, physical matching,
semiclassical validity/noise or all-sector stability is established.

RESEARCH DECISION: retain the construction and the combined full-graph
method. Park this completed comparison transfer; its bound fits the
unchanged sufficient budget. Give priority to the actual uniform domain
and regulator bridge before evaluating the remaining signed projection.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): compute the complete
actual-minus-comparison zero-winding third response in the same graph,
combining its variable logarithmic principal coefficient with T_low and
the finite paired-boundary contacts BEFORE taking norms. Prove its
matched strong regulator limit on the admitted mass/log graph family
and bound its seed/variation pair against (19). This is one actual
defect estimate, not a new collection of unconnected unsigned vertices.
The weighted comparison test (18) is executed; that uniform actual
defect calculation and the full signed projection are not.

Remaining to-do list: actual uniform curved/low-contact graph response
and common regulator removal, global admission/correction and original
X comparison, nonlinear closure, exact compatibility, physical matching/
causal realization, semiclassical validity/noise, all-sector stability,
charge/return, opening, transport, preparation, separate creation and
the full stable-wormhole objective.
