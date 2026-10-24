# Uniform near-domain return on an admitted mass/log family

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
Full stable-wormhole objective ACTIVE and INCOMPLETE; audit CLOSED.

This executes the localized calculation selected after the bolt frequency
test. No global clock gap is used. The result is for the actual massive
canonical block minus its identical-collar cylinder comparison, including
the seagull and diagonal jumps. The physical operator remains the finite
geometry with its prescribed pair and regular bolt. Admission of the
unknown nonlinear solution to the family below is a separate requirement.

## 1. Domain, norm and scope of uniformity

Keep a=.01, bcut=.02, c=.04, sigma=4.005 and all original units.
Q_c=Q through c, followed by the retained C8 taper to beta outside .05.
For a real physical mass field b define b_c by the SAME b through c,
then a fixed scalar taper to zero at .05. Retain the cubic paired seam;
this extension does not smooth or replace the physical domain.

Use the following coefficient norm, with derivatives on the physical
side of the common collar:

 ||b||_A=max{||b||infinity on the actual geometry,
             max_(1<=j<=3) sup_(0<z<c)|partial_z^j b|,
             sup_(0<z<c) z|partial_z^4 b|}.          (1)

The cylinder extension has no larger supremum norm. No derivative bound
on its remote taper enters the numerical kernel estimate. Fields and
directions must also belong to the existing admitted mass/log class,
be identical on the paired copies, and be regular even disk fields at
the bolt. In particular b(0)=b'(0)=0. For the qualitative regulator
argument retain bounded admitted interior/bolt charts and a fixed finite
log degree, as in the existing weak-action theorem. An unbounded union
of coefficient classes is not a uniform family. The numerical result
below is uniform over real backgrounds with ||b||_A<=.001.

Both exp(b) and its inverse preserve the fixed paired first-order
domain and its complementary adjoint domain. Use Dsharp and
A_M=Dsharp+M exp(b); the positive square is its composed square.
This is the canonical mass representation, not a guessed geometric
square domain. The local conversion between representations has the
same coefficients and boundary completion on the two common collars,
so its relative localized variations cancel before regulator removal.
It remains present in the full physical action.

Let Delta H_b be the difference of the two canonical Hessians, with
input restricted to 0<z'<bcut and output 0<z<a. Input restriction splits
an operator; it does not impose a new boundary on metric fields. Define

 R_b g=chi_a z^(5-sigma) partial_z
         [Delta H_b (1_(0<z'<bcut) z'^sigma B_sigma^-1 g)]. (2)

The original action factor -2G/B^2 is INCLUDED. This is the shape block
in the existing coupled H3 graph; the scalar row of this relative block
is zero. It is not the complete quantum Hessian or its complete metric
variation. The LLL stays in its existing separate allocation; here the
massive tower is n>=1. We will bound R_b and D R_b[v] with ||v||_A explicit.

## 2. Actual bolt data are holomorphic and bounded, not guessed

At fixed M and p=gamma(j+1/2)>0 use the retained Riccati equation in
each stable propagation direction,

 Udot=W(1-U^2)-2m U, W=p/sqrt(Q), m=M exp(b).         (3)

The physical regular bolt has U tending to 1. Its regular-singular
Volterra construction, with mass-independent leading normalization,
gives holomorphic parameter dependence. This retains the endpoint
sensitivities in current-endpoint-thermal-restoration.md; it does not
start them at zero at a finite cut. The cylinder has an exactly constant
end with its decaying root; its inward solution is equally holomorphic.
Reflection gives the other branch V(0)=U(0).

For a scalar complex parameter lambda with |lambda|<1/2 and a direction
||v||infinity<=1, Re(M exp(b+lambda v))>0. On |U|=1,

 d|U|^2/dY = -4 Re m <0.                            (4)

The W term has zero real radial component. Thus the unit disk is
invariant, for the actual regular branch and the cylinder branch, as
well as for the reflected branch V. Regular-tip limits are taken from
their Volterra solutions; compact-interval uniqueness and the invariant
disk give the same limit. There is no extra bolt extension parameter.
Bounded holomorphic dependence then gives, on |lambda|<=r=.001,

 |U_lambda-U_0|, |V_lambda-V_0|
              <=r/(.5-r)=:delta_U<.002005.          (5)

This is a Cauchy coefficient estimate in a METRIC PARAMETER. It does not
assume spatial analyticity of a mass/log field or shift a thermal mode.

On the common collar W>=sqrt(5)p, W<=4p. Put
k=sqrt(M^2+5p^2), m_+=exp(.001)M and m_-=exp(-.001)M.
For a real background the constant subsolution starting at zero gives

 integral (m+WU) dY >=[k-(m_+-m_-)]ell-log 2.        (6)

Indeed use U_* tanh(k_+ t), k_+=sqrt(m_+^2+5p^2),
U_*=sqrt(5)p/(m_++k_+). Integrating 1-tanh costs at most log 2;
k_+>=k gives (6). No derivative of b is used in this comparison.

For two complex solutions with the same collar coefficients, their
difference obeys EXACTLY

 Delta Udot=-[2m+W(U_act+U_c)] Delta U.              (7)

The loss from (5) and |exp(lambda v)-1|<=exp(r)-1 is controlled by

 2[1-(exp(.001)-exp(-.001))]
       -2 exp(.001)(exp(r)-1)-8 delta_U/sqrt(5)>1.98. (8)

This is the rate after dividing by k. Thus, writing Y(z)=int_0^z Q^-1/2,

 |Delta U(Y)|<=8 exp[-1.98k(Y(c)-Y)],
 |Delta V(Y)|<=32 exp[-1.98k(Y(c)+Y)].               (9)

The factors are the endpoint difference <=2, the exp(2 log 2)=4
propagation cost, and the second reflection propagation. These are
uniform actual-boundary estimates. The individual complex propagators
also obey |exp(-2I_AB)|<=4 exp[-1.98k|Y_A-Y_B|] and the two-leg image
bound with prefactor 16. Inequality (8) supplies these bounds too.

## 3. Spatial derivatives without a false analytic rough-domain extension

The old seed argument used spatial analytic germs. Here a physical
mass/log field need not be spatially analytic. At each positive output
or input point replace b LOCALLY by its Taylor polynomial through degree
three, only to calculate its finite jet. The analytic ODE with that
polynomial has exactly the actual U,V jets through order four, since
(3) is first order. Its bubble and seagull jets agree through order
three. The only fourth-order discrepancy is the explicit vertex
exp(b) b''''; it is accounted for below. This use of a finite Taylor
jet is not a surrogate domain, an analytic approximation to the action,
or permission to differentiate a real envelope.

Use the bidisk radius R=1/[1000(M+p)]. The retained Q polynomial has
|Q/beta-1|<.001 there, Q^-1/2<4 and Q^-1<16. Since the complex family
has ||b||_A<=.002, each Taylor polynomial has modulus <.003 in its
disk. In the local Riccati ODE, the value and Lipschitz guards remain
36(M+p). The total displacement from the real branch is less than
delta_U+.036<.1. Consequently |1+UV|>=.79, the rational bubble
prefactor is <2, each of its three partial derivatives is <6, and
each partial derivative of (1-UV)/(1+UV) is <4.

Set E_A=exp[-1.98k(Y(c)-Y_A)]. Propagation along the short germ arms
inflates (9) by less than 1.1, giving Delta U guards 9E_A,9E_B and
Delta V guard 36E_B. The rational-prefactor difference is <400E_A.
Integrating W Delta U or W Delta V on the REAL collar gives constants
32/(1.98 sqrt(5))<8 and 160/(1.98 sqrt(5))<37 for the direct and image
exponent differences. Germ arms cost <.3 and <.72 respectively. Use
the guards 10E_A and 40E_A. Individual two-arm exponential inflation
is <exp(.072); use direct/image prefactors 5 and 20.

The full bubble difference constant, including both density factors
and the two exponential mass vertices, is therefore bounded by

 exp(.006) 16 * 4[(400*5+2*2*10*5)
                         +(400*20+2*2*40*20)] <2^20. (10)

The seagull difference is bounded by
exp(.003) *4*4*(9+36)<2^10. No relative seagull is deleted.
The physical-half and both thermal signs retain their old normalization.

Let d=sqrt(5)(c-a), c_M=.99 sqrt(2)d, c_p=.99 sqrt(10)d and

 F=8G C(1+gamma c_p)/(pi B^2),
 Z_(r,j)=sum_(l=0)^j binom(j,l)(r+j-l+2)! l!
                       /[c_M^(r+j-l+2)c_p^(l+1)],
 K_j=F 2^20 j!1000^j Z_(2,j),
 S_j=F 2^10 j!1000^j Z_(1,j), 0<=j<=4.             (11)

These use the COMPLETE angular counting and discrete thermal sums of
the seed transfer, now with the rate 1.98 instead of 2. No modes are
discarded. The formal-germ Cauchy bounds imply ordinary derivative
guards K_j,S_j through order three. At order four they are

 K_4(z)<=K_4+.002 K_0/z,
 S_4(z)<=S_4+.002 S_0/z.                            (12)

The added terms are precisely the missing b'''' vertex contributions.
They admit the cubic and logarithmic boundary fields in (1). Mixed
derivatives of total order <=3, used for the diagonal jumps, agree
exactly with the germs and need no such addition.

## 4. Common regulator identification in the RELATIVE problem

This step must not be inferred merely from an envelope on the limiting
Green function. Here there is a separate locality argument. On real
admitted backgrounds use the self-adjoint doubled first-order operator
K_M=[[0,A_M^*],[A_M,0]]. Its square consists of the two prescribed
composed squares and has gap >=M^2 exp(-.002). Its principal symbol
has finite propagation speed one in the canonical proper metric.
The elementary local energy proof applies to Lipschitz coefficients
and the paired transmission traces: boundary Green terms cancel at
the common seam, and the actual tip is the regular disk. Multiplication
by the mass changes no propagation speed.

For a precise abstract version, use the bounded-commutator criterion in
[McIntosh--Morris, Theorem 3.1](https://arxiv.org/pdf/1201.1818).
Here a bounded Lipschitz scalar preserves the paired first-order domain,
the first commutator is Clifford multiplication by its differential with
norm <=Lip(eta), the second commutator is zero, and the doubled real
operator generates a unitary group. These verify that criterion on the
rough fold. The cited theorem supplies finite propagation only, not a
mass/log strong quantum-response theorem or any constant in (16)-(19).

If observations and the near input are in z<=bcut, two localized wave
groups on the actual and cylinder geometries agree for |s|<d0, where
d0=Y(c)-Y(bcut)>0. The same is true for their first three metric
variations when the parameter fields agree through c. Duhamel's wave
formula inserts multiplication vertices M exp(b) times products of
directions. Every path has total propagation length at most |s|; it
cannot encounter the differing coefficients or endpoint before d0.
All multiple and coincident insertions are retained by differentiating
the complete exponential, not an isolated triangle.

Now use the exact spectral identity

 exp(-tK_M^2)=(4pi t)^-1/2 int_R exp(-s^2/(4t)) exp(isK_M) ds. (13)

Wave Duhamel norms add only finite powers of |s| and M through three
variations. The relative integral begins at |s|>=d0 and hence has a
Gaussian small-time factor exp(-d0^2/(8t)), times a polynomial in M,t.
Independently, a small complex mass change is a bounded perturbation
V of K_M with ||V||<=c M r. For its SQUARE use the closed form on Dom K_M:
the two cross terms are bounded by 2||K_M psi||||V psi||, and the last
term by ||V||^2||psi||^2. The original gap absorbs these terms for a
fixed small parameter radius. This gives a uniform sectorial square
semigroup bound C exp(-c M^2 t), with c>0. Cauchy in that parameter
then gives the same mass Gaussian for its first three variations.
The complex extension uses the same holomorphic exp(b) in both off-diagonal
blocks; it does not conjugate the parameter in a purported analytic family.
This argument uses bounded mass insertions; it does not require a
derivative of an input indicator. On smooth directions it agrees with
the complete square vertices
P_I=2^|I| M^2 exp(2b) product u_i+M[exp(b) product u_i,Dsharp].
Interpolating the wave-locality and massive-semigroup bounds retains
BOTH a positive Gaussian mass factor and exp(-d0^2/(16t)).

For localized traces and the strong output, apply local elliptic
estimates to this doubled first-order system before taking the trace.
Two local Sobolev gains supply trace ideals in the two orbit variables;
four normal differentiations use at most three coefficient derivatives.
This follows by repeatedly solving the first-order normal equation for
the next normal derivative, with tangential derivatives controlled by
powers of K_M. The resulting constants grow only polynomially in M
and t^-1. Coefficients are W3,infinity on the paired collar; no bounded
fourth derivative across the cubic seam is required. Differentiating
the external mass vertex a fourth time adds (12), whose weighted
output L2 powers are positive, as verified in section 5. Fixed admitted
interior/bolt charts bound the remaining localization constants.
Thus, for some finite N and uniform C on EACH stated bounded admitted
family, the differentiated relative localized heat traces, including
their weighted strong-output norms, are bounded by

 C t^-N (1+M)^N exp(-c1 M^2 t-d0^2/(c2 t)), c1,c2>0. (14)

No numerical value of C is asserted or needed to identify this limit.
The angular count N_ang(M)<=4C_radius M^2 makes the entire angular sum
of (14) integrable at zero even after dt/t; the thermal trace is already
included in the localized two-dimensional trace. The bound is uniform
in the admitted physical-side mass/log approximations. At positive time
the fixed-domain form limit applies. At large time the same mass gap
and localized heat trace bound give domination. This proves removal of
the COMMON regulator for the relative source, Hessian and its first
metric variation in the stated strong output norm.

All local heat/contact coefficients cancel because the physical fields,
projectors and observation/input support agree before reaching c.
The nonlocal seagull and all jumps still survive and are (11)-(12).
Equations (13)-(14) identify the limit; the EXPLICIT Green/Riccati
bounds (11) give its numerical size. An unevaluated heat constant is
not being used as a numerical error certificate for a finite cutoff.
This relative argument does not prove strong regulator removal for the
singular local cylinder response, whose contacts do not cancel.

## 5. Strong return and a derivative bound with the norm comparison explicit

Use the unchanged normal inverse ||B_sigma^-1||_(H3->H6)<1801.
For v=z^sigma B_sigma^-1 g, its derivatives through four obey the
retained bound (1801/sqrt(2))L_m z^(sigma-m)||g||H3, with
L_m=product_(i=0)^(m-1)(sigma+1-i). The regular integral, every diagonal
jump and the seagull are combined exactly as in equations (9)-(13) of
current-coincident-domain-transfer.md. Replace its K,S by (11).
For the fourth response derivative add TWO terms to its positive ledger:

 .002 bcut^(sigma+1) K_0 z^-1,
 .002 S_0 z^(sigma-1).                              (15)

After chi_a z^(5-sigma)partial_z and three Euler derivatives, these
terms have L2(dz/z) powers 2.995 and 7. Both are positive. The old
regular, jump and seagull powers remain positive. This is the explicit
mass/log regularity check; fourth derivatives are not silently bounded.

The complete entire-mode evaluation gives, uniformly on ||b||A<=.001
and its complex parameter tube of radius .001,

 ||R_b||_(H3->H3)<2.205e-13.                        (16)

The exact computed enclosure is 2.204633065426858e-13 (rounded upward).
Normal convergence of (11)-(15), analytic actual endpoint data and the
fixed first-order realization give a Banach-valued holomorphic relative
operator. Cauchy in the metric parameter, NOW with the strong norm
already established, gives

 ||D R_b[v]g||H3 <2.205e-10 ||v||A ||g||H3.         (17)

No five spatial dilation derivatives are taken. This closes a uniform
near-domain return requirement on the stated admitted family.

For an executed graph comparison take the actual fast direction extended
by eta(z)=1 for z<=.02, eta=1-I_t(9,9) on .02<z<.03 with t=100z-2,
and eta=0 beyond .03. This is an auxiliary split of a fast direction,
not a replacement of its physical far part. It preserves the paired
phase and vanishes at the actual bolt. H6 point evaluation gives

 ||partial_x^m w||infinity<=1801 sqrt(a_m/512)||g||H3,
 a_m=(63,7,3,3,7,63)_m, m<=5, w=B_sigma^-1 g.       (18)

Apply the exact falling polynomial product(sigma+partial_x-l) for
ordinary z derivatives and Leibniz for eta. Bernstein coefficients on
64 subintervals bound each taper derivative on its ENTIRE transition.
This gives

 ||eta z^sigma B_sigma^-1 g||A <5.681e6 ||g||H3,
 ||D R_b[eta z^sigma B_sigma^-1 g1]g2||H3
                                   <.001253 ||g1||H3||g2||H3. (19)

At fast radius 1e-11 the A cost is <5.681e-5. The rest of the background
must therefore have A norm <=.00094319 to fit (1)'s .001 hypothesis.
This is an explicit admission target for the remaining mass/slow/global
field, NOT a measured norm of the unknown correction.

For an arbitrary global graph direction, if ||v||A<=K_A,X||v||Xgraph,
(17) would meet the proposed C_return<1e9 target when
K_A,X<4.535e18. That comparison is still UNKNOWN. The explicit fast
extension (19) satisfies its own target by many orders of magnitude;
this is evidence for that component only. It does not set the global
comparison to one or bound a far/slow direction by its outer trace.

## 6. Feasibility consequence and research decision

EXECUTED: a uniform strong near-domain relative Hessian, its first metric
variation per coefficient norm, and the explicit fast-graph comparison,
with the actual bolt realization and common regulator. The complete
angular/thermal tower, seagull, images and diagonal jumps are retained.
The global clock-gap shortcut remains rejected. Cold/scalar and the
controlled flat angular component are not reoptimized.

The admitted-family uniform bound (16) exceeds the old seed bound by
less than 2e-14. Conservatively reserve an ADDITIONAL 1e-12 from the
remaining nonlocal budget, preserving all old allocations. The remaining
target is then

 eta_local_curved_seed+1e-11 C3_other<.159916395999. (20)

This bookkeeping requires the stated A admission and covers the uniform
near return instead of inferring it from its seed value. C3_other still
includes every unbounded local/far/global component. Equation (19)
does not certify the entire C3. Both quantities in (20), seed contact
M_A<=3.18e27 and M_V<=2.66e30, and original X/graph comparison remain
unknown. No excessive complete bound or actual obstruction for them is
established. The new admission target is quantitatively compatible with
the fast ball, but its remaining global field has not been computed.

The original rho=1e-8 conditions remain
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1.
Their Y0<=5e-9, Z0<=.25, Z2<=5e7 allocations are unproved. The earlier
selected E_rest<2.123e-42<1.287e-40 screen keeps its existing scope;
zero is an enclosure center, not a measured full signed projection.
No actual global correction or controlled nonlinear solution is obtained.

RESEARCH DECISION: retain the construction and the direct strong method.
Park this now controlled NEAR relative-domain component on the explicit
admitted family. Do not call its conditional admission a preserved X ball.
The remaining hard local response must be handled with its contacts
combined before norms and its highest derivatives kept in L2.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete
local cylinder massive Hessian minus the matched normal operator and
its first metric variation in the coupled H3 graph on the admitted
mass/log family. Retain the thermal tower, finite local contacts and
the already certified actual-domain return; test (20) together with
the contact targets and explicit original norm comparisons. The local
common-regulator strong limit, unlike the relative one above, still
requires its full ultraviolet cancellation argument.

Remaining to-do list: local curved strong response and contacts, far/global
metric response and admission, actual correction, original X/graph
comparison and nonlinear closure, exact compatibility, physical matching/
causal realization, semiclassical validity/noise, all-sector stability,
charge/return, opening, transport, preparation, separate creation and the
full stable-wormhole objective.
