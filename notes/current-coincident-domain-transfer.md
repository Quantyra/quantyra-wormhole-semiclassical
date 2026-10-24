# Actual-domain transfer of the near quantum Hessian, with its diagonal terms

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
Full stable-wormhole objective ACTIVE and INCOMPLETE; feasibility audit CLOSED.

The separated-input result does not bound the change of the coincident
Green function caused by the distant geometry. That change affects both
the bubble and the local exponential-mass seagull. We now bound this
COMPLETE domain-transfer difference through the required output order.
The singular local quantum comparison itself remains open.

## 1. Exact comparison, preserving the actual paired collar

Keep the actual seed, q,C,B,gamma, mu0, paired seam, mass/log phase and
regular-bolt realization. Put a=.01, b=.02, c=.04 and sigma=4.005.
Define an AUXILIARY complete even cylinder by

 Q_c(z)=Q(z), 0<=z<=.04,
 Q_c(z)=chi(z)Q(z)+(1-chi(z))beta, .04<z<.05,
 Q_c(z)=beta, z>=.05,

with chi=1-I_t(9,9), t=(z-.04)/.01. Reflect using Q_c(|z|) and retain
the original thermal period. The positive beta taper is C8; 0<Q_c<.2.
It leaves the actual cubic seam and its C2,1 regularity UNCHANGED.
The complete cylinder Dirac operator has its usual closed H1 realization;
the mass M>0 makes its inverse unique. Its two decaying endpoint branches
give the same Riccati invariant interval 0<=U,V<=1 as the actual bolts.

More explicitly, conjugation by the radial density gives J partial_rho
plus a bounded skew-adjoint potential on the complete cylinder, with
J constant, self-adjoint and invertible. Its adjoint domain is H1: a seam jump would
create a delta term outside L2. Integration by parts cancels the two seam
boundary forms precisely for this continuous transmission domain; there
is no free seam extension parameter. The constant-coefficient operator
at each complete end has the same H1 domain, and the bounded interior
perturbation preserves it. In the skew-adjoint Dirac convention,
||(D+M)f||^2=||Df||^2+M^2||f||^2; the same estimate for the adjoint
gives surjectivity and the inverse bound 1/M. This supplies the cylinder
realization used here, rather than relying on formal symmetry alone.

This cylinder is an explicitly controlled comparison for the canonical
massive block, NOT a physical replacement geometry or a smooth surrogate
domain. All other physical local operators are retained as actual operators.
Write H_act=H_c+DeltaH for the canonical massive Hessians on observations
z<a and inputs z'<b. The comparison is useful only with the bound on
DeltaH proved below. No full determinant value on an infinite cylinder
is needed: the localized first and second variations define the difference.

The backgrounds, variation directions and paired projector coincide on
the entire support of those variations. Their local ultraviolet heat
coefficients, seam representatives and prescribed finite restorations
therefore cancel in the difference, in the SAME common regulator scheme.
This does not authorize deleting them from H_c or the actual full action.
The remaining seagull depends on the global Green function and does NOT
cancel merely because the metrics agree locally.

## 2. Propagate the difference of actual endpoint data

Use Y(z)=int_0^z Q^-1/2 and Y_c=Y(c) on the common collar. For each
M_n^2=n(n+q)/C and p_j=gamma(j+1/2), let k=sqrt(M2+5p2).
The stable equation is Udot=W(1-U2)-2MU, W=p/sqrt(Q)>=sqrt(5)p.
Two solutions obey exactly

    DeltaUdot=-[2M+W(U_1+U_2)]DeltaU.

The constant-coefficient subsolution U_* tanh(kt),
U_*=sqrt(5)p/(M+k), gives integral(M+WU)>=k ell-log2 on every
propagation segment. The unknown actual and cylinder data at Y_c both
belong to [0,1]. At 0<=Y<=Y_c this proves

 |DeltaU(Y)|<=4 exp[-2k(Y_c-Y)],
 |DeltaV(Y)|<=16 exp[-2k(Y_c+Y)].                    (1)

The second formula propagates the exact reflected datum V(0)=U(0).
It retains the original paired realization, rather than guessing the
bolt datum or identifying the physical Green function with the cylinder.

For inner A and outer B on the common collar, the direct and image
bubbles are the previously derived expressions

 Bdir=-4 U_A V_B exp(-2I_AB)/(1+U_B V_B)^2,
 Bimg=-4 U_A U_B exp(-2I_img)/(1+U_B V_B)^2.          (2)

Their normalization includes both thermal signs and the physical half.
The diagonal first-order spin trace is

    Tdiag=M(G_++G_-)=(1-UV)/(1+UV).                 (3)

For Gamma=-Tr log(D+M exp(v)), its second variation is bubble minus
M Tdiag times the product of the two directions. Thus DeltaH is the
difference of (2) PLUS the local seagull from the difference of (3).
These are the full domain-dependent terms in the same canonical chart.

## 3. Analytic germs give derivative bounds; no real envelope is differentiated

For real output z<=a and input z'<=b use the bidisk radii

    R=1/[1000(M+p)].                                (4)

M>=18 keeps these disks within |zeta|<.05 and away from the cap change
at .04. The actual polynomial gives |Q/beta-1|<.001 on |zeta|<=.05,
hence |Q^-1/2|<4 and |Q^-1|<16, with a fixed square-root branch.
As before, in |U-U_real|,|V-V_real|<=.1 the complex z-ODE has value
and Lipschitz bounds 36(M+p). Picard gives displacement <=.036 and
|1+UV|>=.79. Difference propagation on each short complex arm inflates
the REAL differences (1) by less than 1.1.

Here is an explicit complex comparison ledger, avoiding an inference
from a real-axis bound. Write E_A=exp[-2k(Y_c-Y_A)] with Y_A>=Y_B.
On the bidisks,

 |DeltaU_A|<=5E_A, |DeltaU_B|<=5E_B,
 |DeltaV_B|<=20exp[-2k(Y_c+Y_B)]<=20E_B.

For P=U_A V_B/(1+U_B V_B)^2, or its image version U_A U_B/(...)^2,
|P|<2 and each of its three first partial derivatives is <6 on the
convex interpolation of the two solution germs. Thus |DeltaP|<=180E_A.

The real direct exponent difference is at most 4E_A: integrate
W|DeltaU| and use W<=4p and p/k<=1/sqrt(5). Its two complex arms cost
less than .16E_A, so |DeltaI_AB|<5E_A. The two real image legs cost
at most 4E_A+15E_0; their complex arms cost <E_A, so take
|DeltaI_img|<24E_A. The separate complex exponents satisfy

 |exp(-2I_AB)|<=5exp[-2k(Y_A-Y_B)],
 |exp(-2I_img)|<=18exp[-2k(Y_A+Y_B)].

These follow from the real path bounds and exp(.072)<1.125 for the
two endpoint arms. The exact integral identity for the difference of
two exponentials is valid also at complex arguments. In (2) it gives
direct and image difference constants at most

 4[180*5+2*2*5*5]=4000,
 4[180*18+2*2*24*18]=19872.

Multiplying by both Q^-1/2 density factors costs 16. The sum is <2^20.
Both exponents are bounded by exp[-2k(Y_c-Y_output)]. For (3), each
partial derivative is <4, so its difference is <=100E_output; the
single output density costs four, below the reserve 2^10.

Consequently the piecewise analytic bubble density difference, per dz',
and the diagonal seagull germ have the bounds

 |DeltaB_density|<=2^20 exp[-2k(Y_c-Y_output)],
 |DeltaTdiag/sqrt(Q)|<=2^10 exp[-2k(Y_c-Y_output)].   (5)

Cauchy on (4) bounds any mixed derivative of the bubble by i!j!R^(-i-j)
times (5), and each seagull derivative by j!R^-j times (5). The germs
are those of the ACTUAL one-sided coefficients and endpoint solutions.
They do not change the physical seam or assume smooth rough-background
fields. We only claim a fixed-seed operator result, not uniform nonlinear
X-ball admission or a Z2 estimate.

## 4. Entire-mode derivative constants, including the seagull

Let d=sqrt(5)(c-a), c_M=sqrt(2)d and c_p=sqrt(10)d. Then
2k(Y_c-Y_output)>=c_M M+c_p p. Reuse the proved complete count and
discrete thermal bounds, for integers r>=1 and k>=0,

 sum_n d_n M_n^r exp(-c_M M_n)<=4C(r+2)!/c_M^(r+2),
 gamma sum_j p_j^k exp(-c_p p_j)
                      <=k! c_p^(-k-1)(1+gamma c_p). (6)

No frequency or angular cutoff is used. Put

 F=8G C(1+gamma c_p)/(pi B2),
 Z_(r,j)=sum_(k=0)^j binom(j,k)(r+j-k+2)! k!
                         /[c_M^(r+j-k+2)c_p^(k+1)],
 K_j=F 2^20 j!1000^j Z_(2,j),
 S_j=F 2^10 j!1000^j Z_(1,j), j=0,...,4.            (7)

The factor -2G/B2 is the original canonical-to-normalized action factor.
K_j bounds each ordinary j-th output derivative of the bubble difference
on either side of z=z'; mixed derivatives with total order j replace j!
by the corresponding factorial product. S_j bounds the complete local
seagull coefficient's j-th derivative. These definitions retain the
original physical units and common finite part.

The sums of all these bounds converge absolutely. Local heat contacts
are identical, so their difference is zero before regulator removal.
Dominated convergence therefore removes the COMMON regulator in the
difference, including its first two metric variations. This resolves
the transfer step even at coincidence. It does not supply the singular
local operator H_c, whose contacts still require estimation.

This cancellation can be read directly as a relative prescription: use
the same regulator for both resolvents, subtract their explicit kernels
and the seagull, and then remove it under the absolutely convergent
majorants (5)--(7). No smooth-seam heat expansion is needed for this
relative difference. Renormalization and regularity of either full
functional separately remain the earlier admitted-family obligations.

## 5. The diagonal jumps must be retained in the strong response

The bubble kernel is continuous across z=z', but its first derivative
need not be. For example in the proper propagation coordinate (2) gives

    [partial_x Bdir]_diag=4W(U+V)/(1+UV).            (8)

The image is smooth across this diagonal: using
(log(1+UV))_Y=W(U-V), its expression factors symmetrically into a
function of Y_A times the same function of Y_B, divided by (1+U(0)^2)^2.
For the difference kernel let J_r(z) be the jump of its r-th output
derivative, including its density factors and entire-mode sum. J_0=0.
For Fv(z)=int_0^b DeltaK(z,z')v(z')dz', distributional differentiation gives

 (Fv)^(n)=int_0^b (partial_z^n DeltaK)_reg v dz'
                  +sum_(r=1)^(n-1) partial_z^(n-r-1)(J_r v),
                                                   n=1,...,4.    (9)

Thus one cannot replace this response by an integral of the fourth
pointwise derivative alone. The local terms in (9) are determined by
the actual kernel; they are not new counterterms or adjustable contacts.
The bidisk Cauchy bound and the binomial identity give

 |J_r^(l)|<=2(r+l+1)/(r+1) K_(r+l), r+l<=3.          (10)

Indeed sum_(t=0)^l binom(l,t)(r+t)!(l-t)!
=r!l! binom(r+l+1,l). This counts both sides of the diagonal and every
mixed derivative when differentiating a diagonal value.

## 6. Complete H3 return in the original coupled graph

Use v=z^sigma B_sigma^-1 g restricted to 0<z'<b. The input indicator is
only a split of this bounded difference; it is not an inadmissible metric
used in the full action. Its discontinuity at b is separated from output.
The actual normal return and H1 point evaluation imply, for m<=4,

 |v^(m)(z)|<=C3/sqrt(2) L_m z^(sigma-m)||g||_H3,
 L_0=1, L_m=product_(i=0)^(m-1)(sigma+1-i), C3<1801. (11)

All coefficients of the shifted falling polynomial are positive in this
range; H_(m+1)<=H6 suffices. In particular the regular bubble integral
uses ||v||_infinity<=b^sigma C3||g||_H3/sqrt(2).

For n=1,...,4, (9)--(11) bound the n-th response derivative, after
dividing by C3||g||_H3/sqrt(2), by the positive expression

 R_n(z)=b^(sigma+1)K_n
  +sum_(r=1)^(n-1) sum_(l=0)^(n-r-1) binom(n-r-1,l)
       [2(r+l+1)/(r+1)]K_(r+l)L_m z^(sigma-m), m=n-r-1-l
  +sum_(l=0)^n binom(n,l)S_l L_(n-l) z^(sigma-(n-l)).               (12)

The three lines retain the regular kernel, every diagonal jump and the
seagull separately. Signs are allowed to cancel physically; no cancellation
is needed for this upper bound.

The outer differentiated row is chi_a z^(5-sigma)partial_z DeltaH v.
Apply the same output cutoff guards (1,3,10,60) and Euler product
recurrence C_(m+1),j=(5-sigma+j)C_mj+C_m,j-1, C_00=1.
For each k=0,...,3 sum

 binom(k,l)h_l C_(k-l),j z^(5-sigma+j)R_(j+1)(z)

over 0<=l<=k, 0<=j<=k-l. Replace each positive term A z^p by
A a^p/sqrt(2p), its exact L2(dz/z) norm on (0,a). Every p is positive:
the regular powers begin at .995, seagull powers at 4, jump powers at 6.
Call the resulting upper bounds N_k. Then

 eta_domain <=C3/sqrt(2)[sum_(k=0)^3 binom(3,k)N_k^2]^(1/2).        (13)

This is the maximum row sum of this shape-only coupled feedback block.
The scalar row is zero; all other actual local terms remain unchanged.
Freeze the target eta_domain<1e-10. Record the sharper bound, its three
parts, and the remaining allocation .18699-eta_domain.

## 7. Research consequence and the next calculation

This executes the actual-to-collar transfer INCLUDING the coincident
state dependence, seagull and differentiation jumps. It is distinct from
the already controlled far-input operator. No common ultraviolet contact
is inferred from a separated kernel or omitted from the local problem.
The two differences are allocated once: the previous result covers
physical input z'>=.02; (13) covers actual-minus-cylinder on input z'<.02.

RESEARCH DECISION: retain the construction and use this certified transfer
to study the near operator on the explicit cylinder with the SAME cubic
paired seam. The physical operator remains the actual finite geometry plus
the rigorously retained difference. No smooth-domain surrogate, artificial
physical wall or freely chosen quantum source is introduced.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete local
cylinder massive Hessian plus the actual Weyl/LLL/finite-contact terms
against the full normal operator, through the required H3 output order,
including the original cubic seam. Match the logarithmic principal and
finite fourth contacts before estimating the gain-two remainder and its
lower contacts. The conservative combined remaining target is <=.18698
if (13) leaves it available; other still-unknown original-domain feedback
must fit the same remainder allocation.

The coincident transfer is EXECUTED only by its retained receipt. The full
local quantum comparison, joint affine/known-source bounds, original-row
mass/slow moment and independent Mellin Efar<=1e-10 remain UNKNOWN.
The joint targets S_coupled<=7.49e-12 and k_coupled<=1e-14 are not measured.
No actual controlled correction or nonlinear solution is certified.

The sufficient projection screen E_rest<2.123e-42<1.287e-40 retains its
scope and unmeasured sign. At rho=1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and
Z0+1e-8 Z2<1 remain unclosed; Y0<=5e-9, Z0<=.25, Z2<=5e7 are allocations.
The C2 action on bounded admitted stationary mass/log families does not
prove nonlinear X-ball admission. Exact compatibility, physical matching,
causal realization, semiclassical validity/noise and all-sector stability
remain open. The full stable-wormhole objective and closed audit are intact.

Remaining to-do list: local paired quantum comparison and remaining global
assembly; joint sources, matching coefficients and independent far-tail
return; actual correction and nonlinear admission; compatibility, physical
matching, validity/noise and all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole goal.
