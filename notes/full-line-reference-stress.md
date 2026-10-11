# Full-line current reference with spherical photon propagation

2026-10-11, S053. Keeping the unwrapped fermion current on the entire
optical line removes the artificial mouth term from the first photon
vacuum diagram. For a specified spherical test-field reduction, the
diagram is a positive scalar resolvent trace. A fixed bulk subtraction
then defines its energy and two static derivatives, including positive
inverse-square photon tails. The resulting bound is small relative to
the working free-clock stiffness even with the stated cutoff dependence
of profile derivatives. The actual two-mouth interaction, complete
finite-return diagram, higher orders and physical state remain open.

## 1. Model and the complete order-g^2 reference

Prescribe a static globally spherical metric
`ds^2=N(s)^2(-dt^2+ds^2)+R(s)^2 dOmega^2`, s in R, and magnetic flux q.
After the usual spinor rescaling, the q angular zero modes are massless
Dirac fields on the optical line. Take their free neutral vacuum and
the constant Hermitian LLL angular vertices T_A, with

    kappa_l=l(l+1), tau_l=q beta_l/(4pi),
    beta_l=product_(j=1)^l (q-j)/(q+j), 1<=l<q.

Here tau_l=Tr(T_A^2), repeated 2l+1 times. Noncommuting vertices are
not replaced by independent bosons; only the exact free two-point
function enters this diagram. On a common unwrapped line it is diagonal
in A. The electric Maxwell test-field scalar has

    H_l=-partial_s^2+V_l(s),  V_l=kappa_l w(s), w=N^2/R^2.

This uses a prescribed metric. Coupled metric/isometry gauge modes,
higher Landau levels and changes in the supporting state are absent.
A globally spherical line is a reference geometry, not the actual
nonspherical exterior joining the two mouths in one asymptotic region.

At k>0 let A0=-partial_s^2+k^2, L u=(u',-ik u),
R_k=(A0+V)^-1 and U=L A0^-1/2. Then U*U=1. The full-line free current
projection is P=UU*, and the weighted photon propagator is
D=1-L R_k L*. Thus

    U* D U=1-A0^(1/2) R_k A0^(1/2)
          =B_k(1+B_k)^-1,
    B_k=A0^-1/2 V A0^-1/2 >=0,
    F_V(k)=Tr[B_k(1+B_k)^-1]
          =Tr[V^(1/2) R_k V^(1/2)].                       (1)

For V>=0 in L1 intersect L-infinity, B_k is positive trace class, with
Tr B_k=I/(2k), I=int_R V. These are traces of positive finite operators.
The isometry proves (1) without cyclic manipulation of divergent
identity traces. The full order-g^2 unwrapped reference is

    E_ref^(2)=g^2/(2pi^2) sum_A (tau_A/kappa_A) J_mu[V_A],
    J_mu[V]=int_0^infinity [F_V(k)-I/(2sqrt(k^2+mu^2))] dk, (2)

with mu>0 fixed before variation. A specified finite local term is
discussed below. All current vertices, including exterior and cross-mouth
ones in this reference, are included in (1). Splitting the full line
into regions does not change it. In particular there is no separate
interface subtraction. Restricting P to the throat before forming the
trace gives the extra term in the [previous counterexample](reference-mouth-completion.md).

The use of spherical separation and conformally rescaled zero modes
is established in [MMP, Appendix C](https://arxiv.org/html/1807.04726v3#A3)
and the Maxwell reduction cited in [the open-load note](open-photon-load.md).
MMP Appendix A and [Georgi-Warner, section 1](https://arxiv.org/html/1907.12705v2#S1)
already supply the commuting massive-vector baseline. Equation (1)
is the full-line form of its quadratic Ward identity, not a new
bosonization of all interacting angular currents.

## 2. Infrared estimate with inverse-square tails

Assume V<=M^2 on a central interval of length L, and at distance x>=0
from either end assume

    V(s)<=U_tail^2/(a+x)^2,  a>0.                         (3)

Increase M if necessary so V<=M^2 globally. All quantities in (3)
may depend on the finite model; the constants below do not. Define

    Q0=1+LM/2
       +2(1+U_tail/2)[2+(log_+(U_tail/a)+1)/log(2)].       (4)

The zero-tail case is interpreted with its tail contribution omitted.
Then

    int_0^1 F_V(k) dk <= Q0.                             (5)

Proof: impose independent Neumann cuts at the core endpoints and at
tail distances a(2^n-1), n>=0. The decoupled form is smaller than
the full form, so its resolvent bounds (1) above. On a cell of length
b with 0<=V<=v^2, monotonicity of Tr[B/(1+B)] in V and the Neumann
eigenvalues give

    Tr[V(-partial_N^2+k^2+V)^-1]
       <= v^2/(k^2+v^2)+b v^2/(2sqrt(k^2+v^2)).           (6)

The nonzero mode sum is bounded by the integral of its decreasing
summand. The central cell contributes at most 1+LM/2. On tail cell
n, take b=a 2^n and v=U_tail/(a 2^n). Its contribution is at most
(1+U_tail/2) min(1,U_tail/(ak 2^n)). The geometric sum is bounded by
2+log_+(U_tail/(ak))/log(2). Integrating over k in (0,1), and using
int_0^1 log_+(A/k) dk<=log_+(A)+1, proves (5).
Every finite Neumann partition has the claimed form ordering. The
infinite positive trace sum is justified by monotone limits. No gap,
finite exterior volume or finite first spatial moment of V is needed.

## 3. Fixed subtraction and two general profile derivatives

Let phi range in a compact interval, with V_phi twice continuously
differentiable in L1 and L-infinity, and suppose pointwise almost
everywhere

    |D V|<=c1 V, |D^2 V|<=c2 V, D=partial_phi.             (7)

The same central interval and tail bounds may be used locally in phi.
The free optical operator A0 is fixed; any changing coordinates must
first be expressed in this common optical line. Define
C0=1, C1=c1 and C2=c2+2c1^2, allowing an absolute numerical factor
in the following estimate. Then

    |D^j J_mu[V]| <= C Cj {Q0+I[1+log(1+M)+|log mu|]},
                                              j=0,1,2.   (8)

For j with Cj=0 the corresponding derivative vanishes. To prove this,
write Ck=(1+B_k)^-1. Since A0 is fixed, -c1 B<=B'<=c1 B and
-c2 B<=B''<=c2 B in form order. The exact trace derivatives are

    F'=Tr(Ck^2 B'),
    F''=Tr(Ck^2 B'')-2 Tr(Ck B' Ck^2 B').                 (9)

They imply |F'|<=c1 F and |F''|<=(c2+2c1^2)F. For example write
B'=B^(1/2) K B^(1/2), ||K||<=c1. The quadratic trace is the squared
Hilbert-Schmidt norm of Ck B' Ck^(1/2), bounded by c1^2 F using
||Ck B^(1/2)||<=1 and Tr(B Ck)=F. Together with (5) this controls
the differentiated infrared integrals.

For k>=K=max(1,sqrt(2)M), ||B||<=M^2/k^2<=1/2. The convergent
trace series for F-Tr B and its first two derivatives yields

    |D^j(F-I/(2k))| <= C Cj I M^2/k^3.                  (10)

Indeed the n-th trace and its differentiated versions are bounded by
C Cj n^2 ||B||^(n-1) Tr B; sum from n=2. Schatten Holder with
the relative factors in (7) justifies the bound when derivatives do
not commute with B. On 1<=k<=K, use F<=Tr B and (9), at a cost
C Cj I log K. Finally the integrals of the two elementary subtraction
differences on (0,1) and (1,infinity) are bounded by
C I(1+|log mu|). This proves (8) and a locally integrable majorant
for each derivative. The expansion never subtracts the full Born term
at k=0, where that would introduce an artificial infrared divergence.

A fixed finite local term int V f_phi changes (8) by at most
C I times the corresponding bounds on f, Df and D^2f, combined with
(7). This covers a supplied finite coefficient and a supplied local
optical-to-proper subtraction term. Their physical four-dimensional
coefficients must still be matched; they are not chosen here to alter
the curvature. No independent bare-photon vacuum energy is included
in the interaction diagram (2).

For an independent normalization check, global dilation
V_phi(s)=exp(2phi)V0(exp(phi)s) gives exactly

    J_mu[V_phi]=exp(phi)[J_mu[V0]-I0 phi/2].              (11)

Thus the optical subtraction is not silently treated as dilation
invariant. At phi=0 its first two derivatives are J_mu[V0]-I0/2
and J_mu[V0]-I0. This check alone is not the general clock theorem;
the theorem uses (7), which must be verified for the intended profile.

## 4. Angular sum and the declared weak-coupling family

For the LLL weights, exactly

    sum_A tau_A/kappa_A=q(H_q-1)/(4pi),
    sum_A tau_A=q(q-1)/(4pi),
    sum_A tau_A/sqrt(kappa_A)<=C q^(3/2).                 (12)

For the last inequality, beta_l<=exp(-l(l+1)/q) follows by taking
logs of its product. Also (2l+1)/sqrt(l(l+1))<=3; bound the resulting
Gaussian sum by its integral. This avoids an unnecessary extra angular
factor or an assumption that the noncommuting currents decouple.

For a family of reference profiles assume, uniformly in phi,

    L=O(1), ||sqrt(w)||_infinity=O(epsilon^-1),
    int w=O(epsilon^-1),
    w outside the core <= C/(epsilon+x)^2,
    |D^j w|<=C epsilon^-j w, j=1,2.                     (13)

These are explicit profile and derivative assumptions. The class is
nonempty: for |s|<=pi/2 take w_epsilon=1/(cos^2(s)+epsilon^2), and
outside take w_epsilon=1/[epsilon^2+(|s|-pi/2)^2]. The values and first
two spatial derivatives match at the joins; the integral is
pi/epsilon [1+1/sqrt(1+epsilon^2)]. The tails obey (3) with
U_tail^2=2, a=epsilon. For w_phi(s)=exp(2phi)w_epsilon(exp(phi)s),
the logarithmic derivative factors and the second relative derivatives
are O(epsilon^-1) and O(epsilon^-2), respectively, on compact phi
intervals. A fixed larger central interval covers the moving joins.
This example establishes consistency of (13), not the actual corrected
two-mouth geometry. Deriving those bounds for its physical clock
variation remains a matching task; no global dilation is identified
with a physical throat-only perturbation.

Apply (8) with V_l=kappa_l w, so I_l=kappa_l int w,
M_l=O(sqrt(kappa_l)/epsilon), U_tail,l=O(sqrt(kappa_l)), a=epsilon.
For polynomial reference scales and matched local terms whose relevant
coefficient bounds grow at most logarithmically with q, (12) gives

    |D^j E_ref^(2)|
       <= C_j epsilon^-j g^2 {
          q log q + q^(3/2)[epsilon^-1+log(q/epsilon)]
          +q^2 epsilon^-1[1+log(q/epsilon)+|log mu|]},
                                                   j=0,1,2. (14)

At the working g^2=q^-2, epsilon=q^-1/4, delta=q^-1/16,

    |E_ref^(2)|=O(q^(1/4) log q),
    |D E_ref^(2)|=O(q^(1/2) log q),
    |D^2 E_ref^(2)|=O(q^(3/4) log q).                    (15)

These absolute quantities need not tend to zero. The curvature is
O(q^-3/16 log q) relative to the free damping-scale stiffness q delta.
Against a free static curvature >=c q, this diagram alone shifts a
local equilibrium by O(q^-1/2 log q). The previously derived free
length-sensitivity bound would make that recentering error
O(q^-5/16 log q) relative to q delta at imaginary resolution eta~delta.
These are comparisons of this component, not a completed interacting
decay theorem or a statement at fixed four-dimensional coupling.

## 5. Verification, prior work and decision

The [script](../checks/full_line_reference_stress.py) and
[receipt](../evidence/full-line-reference-stress.json) verify the full
projection identity and noncommuting trace derivatives independently
in compatible finite matrices. Continuum checks compare rectangular
barrier Green traces with a Jost determinant derivative, and include
inverse-square exteriors without truncating their integration domain.
The renormalized barrier integral is also checked against the exact
dilation identity and its clock derivatives. Angular sum checks
support the quantitative bound; the ordinary proof gives
uniformity beyond the finite sample.

The current anomaly, massive-vector reduction and imaginary-frequency
Born subtraction are established tools. [Graham et al.,
hep-th/0207120v2](https://arxiv.org/html/hep-th/0207120v2), especially
sections 2.3-2.4 and Appendix A, supply the comparison with fixed
renormalization conditions and Jost determinants. The present positive
trace and tail bounds are derived above. The source supplement records
which original derivations were read or reused; the complete literature
and forward-citation review still has not passed its exit conditions.

**Decision:** continue the physical damping primary and matched
finite-duration transfer fallback. The complete unwrapped current
reference in this spherical model has a controlled first static
interaction diagram. The former artificial mouth boundary is absent.
Do not add this result to the earlier throat-only winding diagram and
call the sum complete: the full finite-return interaction across the
actual nonspherical exterior is still missing. Also retain the physical
renormalization, metric/isometry sectors, higher-current diagrams and
state requirements. No novelty, physical damping, finite-disturbance
survival or travel-time advantage is established. No Lean/Lake was run.

Remaining to-do list: match the full finite-return mouth interaction and
physical renormalization, bound higher-order/dynamical stress, close the
coupled geometry and state, complete novelty and the consequential result,
then independently verify, GCP-certify suitable endpoints and prepare the manuscript.
