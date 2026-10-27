# Signed weak source: interpolation test and an actual-norm transposition

2026-10-27 HST, E003/S169, second turn of the new research goal. Previous
turn: PROGRESS in method selection, rejecting direct relative-trace transfer.
The full stable semiclassical wormhole objective remains ACTIVE/incomplete.
Closed audit, accepted results, P32 and all original physical norms remain.

OUTCOME: the proposed tau^(5/4) complete remainder is still UNRESOLVED.
The existing strong/positive-time estimates cannot establish it by taking
their best bound, even frequency by frequency: the resulting energy-dual
upper estimator is nonintegrable. A transposed source representation below
has a proved numerical comparison to the actual E_t dual, avoiding an
unavailable C5 embedding. Its complete signed L2 decay is still missing.
No regulator error, signed Zr, admission margin or controlled solution is
certified by this calculation.

## 1. Actual input, target and first unsupported implication

Use the exact fixed phi=chi_N,r-ell_b z^4, r=1e-110, from the accepted
observation replacement, with u=Vphi and w=Vv, ||v||E_t<=1. V is the
unchanged York return; J(0)=0 and (Q^2 J'/z^4)(L)=0. The original paired
mass/log field, complementary matter adjoint, clocks/period/spin/scale,
cutoff transitions and regular-disk bolt are retained. E_t is the original
closed form energy, not a C5 coefficient space or nonlinear solution ball.

The complete positive-time geometric Hessian on this domain is

    h_tau[u,w]=gamma/(2pi) sum_jk q_tau(a_j,a_k) u_jk w_kj,
    q_tau(a,b)=-tau(a+b)^2
       [g_tau(a^2)-g_tau(b^2)]/(a^2-b^2),
    g_tau(lambda)=lambda exp(-tau lambda).             (1)

Use the diagonal limit where needed. A=iJD0 is the actual self-adjoint
paired first-order operator, A^2=L0; the complementary/composed domains
are unchanged. Formula (1) includes the coincidence vertex and complete
thermal spectrum and finite returns. It is not a finite-mode evaluator.

Write R=-Q'', Delta w=(Qw')'. The full completed reduced difference is

    sigma_tau[phi,v]=h_tau[u,w]-2/(pi tau) int uw dz
      +tau/(480pi) int[4R^2uw+8R(u Delta w+w Delta u)
                                                +8 Delta u Delta w]dz. (2)

This is the accepted common heat subtraction, including its boundary
representative. The retained exact UV logarithm and finite restoration
remain in the action. The complete massive remainder uses

    -(G/(B^2 beta)) int Theta(tau) sigma_tau d tau/tau,
    Theta(tau)=sum_(n>=1)2(q+2n) exp[-tau n(n+q)/C], q=2000^4. (3)

The immediate proposed implication is (2) <= A_phi tau^(5/4) in E_t*,
uniformly at small tau. No value of A_phi or such decay has been proved.
The prior C2 admitted-family theorem instead controls both inputs in a
weighted C5 coefficient norm A, with tau^(3/2). It does not supply the
required energy test norm. The fixed phi has a smooth mollified source,
not an evaluated A or fractional Sobolev norm. We do not assume one.

For this test grant, as an optimistic hypothesis, any finite admitted norm
of the fixed input and a finite positive constant for its strong estimate.
If these have not been evaluated or justified, that is an earlier unknown;
granting them cannot repair the interpolation failure proved next.

## 2. Executed test of the best of the available upper estimates

At positive time the accepted full spectral proof gives

    |h_tau[u,w]|<=64 ||u||infinity ||w||infinity k(tau/2),
    k(tau)=gamma Tr(exp[-tau L0])/(2pi).              (4)

The local terms in (2) have their explicit L2/Delta bounds. Form estimator
E_low by adding their absolute bounds to (4). Form E_high from the complete
admitted-family bound C_phi tau^(3/2)||w||A. The putative improvement is
the pointwise minimum min(E_low,E_high), followed by the energy-dual
supremum and common heat integration. The MINIMUM is a valid upper bound
when both estimates apply; it is this specified estimator that we test.

Use the actual-domain stationary TT sequence already proved in
current-continuum-energy-dual-estimator-test.md. On a compact interval I,
let F be the physical lapse and W=x^2 F'-2xF. With fixed smooth compact
chi,rho and int W rho!=0, set

    c_n=int W chi cos(nx) dx / int W rho dx,
    xi_n=chi cos(nx)-c_n rho,
    b_n=-(x^2F)^-1 int_left^x W(s)xi_n(s)ds,
    k_n=diag(-b_n-2xi_n,b_n,xi_n,xi_n).              (5)

Both original TT constraint terms are present. The exact zero moment
makes the field zero near both original ends and every lift transition;
the regular bolt and paired boundary conditions hold. This is a test
sequence IN the full mass/log domain, not a smooth replacement domain.
Its n is a metric-test frequency, not a quantum mode cutoff. The actual
fixed input phi is unchanged.

The retained complete-action symbol and York inverse give

    ||xi_n||E_t asymp n^2 sqrt(log n),
    ||Vxi_n||infinity asymp 1,
    ||Vxi_n||A asymp n^5,
    ||Delta Vxi_n||2=O(n^2).                         (6)

Here asymp means two fixed strictly positive bounds on this one interval.
The moment correction is super-polynomially small; the York principal
coefficient is -3/2. Coordinate/weight constants are fixed and nonzero.
No numerical value of these constants is used to pay a budget.

For v_n=xi_n/||xi_n||E_t, the leading heat asymptotic k(tau/2)~L/(pi tau)
and the positivity of every term in the reported E_low imply, for all
sufficiently small tau and large n,

    E_low(phi,v_n;tau)>=a tau^-1 n^-2/sqrt(log n),
    E_high(phi,v_n;tau)>=b tau^(3/2) n^3/sqrt(log n), (7)

where a,b>0 for every nonzero finite choice of the fixed-input bounds.
The fixed u is nonzero by the accepted actual-input argument. Additional
positive contact allowances cannot reduce the first line. Choose
n=floor(tau^-1/2). Equations (6)-(7) give, for ALL sufficiently small tau,

    sup_(||v||E_t<=1) min(E_low,E_high)
                      >=c/sqrt(log(1/tau)), c>0.     (8)

Consequently its full angular error expression has lower growth

    int_eta^tau0 Theta(tau) E_best(tau) d tau/tau
                  >= c'/(eta sqrt(log(1/eta)))        (9)

for small eta, after reducing c'>0 if needed. Use Theta(tau)~2C/tau
and integrate only tau in [eta,2eta] to obtain (9). This is a verified
nonintegrable UPPER ESTIMATOR, not a lower bound on the physical response.
It improves on the earlier quadratic divergent upper expression but
still cannot be used for regulator removal. No finite constant reduction
changes (9). The original log energy weight is included throughout.

There is a useful stronger diagnostic: even if a COMPLETE C3 test estimate
of order tau^2 were granted, instead of the available C5/tau^(3/2) bound,
the same sequence with n=floor(tau^-1) gives the best-of-two expression
>=c tau/sqrt(log(1/tau)). Its angular integral diverges as
sqrt(log(1/eta)). The proper-coordinate C3 estimate is not being assumed
valid; this optimistic grant shows why it would still be insufficient
when used only through these separate upper norms.

More generally, a strong test bound tau^s C^k, combined this way with
tau^-1 Linfinity and the actual n^2 sqrt(log n) energy, has envelope power

    beta_heat=2(s+1)/k-1.                            (10)

Integrability needs beta_heat>1; equality still fails with this logarithm.
For the target beta_heat=5/4, k=3 would require s>=19/8 (rather than 2).
At s=2 it would require k<=8/3. This quantifies the gain missing from the
specified estimator. It is not a necessary condition for the actual
response, which can keep the correlation between fixed and test inputs.

## 3. A defined complete source map with an explicit energy comparison

There is an alternative that does not take the strong test norm first.
Introduce ONLY an auxiliary scalar resolvent on H0=L2((0,L),dz):

    D=Friedrichs realization of -partial_z Q partial_z,
    S=(1+D)^-1.

Its form domain imposes value zero at z=0; at the regular bolt its finite
form energy excludes the logarithmic branch. Its operator domain has
QT'(L)=0 and Delta T=-D T in L2. Original energy images w=Vv belong to
this operator domain: the accepted Delta comparison, zero outer trace,
and regular-flux limit prove this. Their additional mass/log constraints
are NOT removed. S is an auxiliary transposition, not a new matter
boundary condition, a new physical metric family or a substitute for V.

From self-adjoint positivity,

    ||S f||2<=||f||2, ||Delta S f||2<=||f||2.

The already proved global radial estimate, with .99<L<1 and
Q>=gamma(L-z)/2, therefore supplies

    ||S f||infinity <= C_S ||f||2,
    C_S=1/sqrt(.99)+4/gamma.                         (11)

For each positive tau, f -> h_tau[u,S f] is a bounded functional on H0
by (4),(11). Thus there is a UNIQUE H_tau[u] in H0 with

    <H_tau[u],f>=h_tau[u,S f],
    ||H_tau[u]||2<=64 ||u||infinity C_S k(tau/2).     (12)

The spectral series (1) extends to bounded scalar smearings by its absolute
positive-time estimate, so using S f requires no new metric-domain
variation theorem. It is a functional extension for the transposition.
It keeps every spectral mode, contact vertex and finite return.

Define the local L2 coefficients

    a_tau=-2u/(pi tau)+tau(4R^2u+8R Delta u)/(480pi),
    b_tau=tau(Ru+Delta u)/(60pi),
    g_tau=H_tau[u]+S a_tau+(S-I)b_tau in H0.         (13)

They are defined because u,Delta u are L2 and |R|<20. There is no separate
fourth derivative of a reflected mass cubic. No boundary integration by
parts is needed for (13): Delta S=S-I follows from the actual scalar
spectral realization. For every original energy test w=Vv,

    sigma_tau[phi,v]=<g_tau,(I-Delta)w>.             (14)

Indeed apply (12) with f=(I-Delta)w, then (13), self-adjoint S and
(S-I)(I-Delta)w=Delta w. This recovers EXACTLY both local weak rows in (2).
The originally completed boundary term is retained in those rows and in
the physical h_tau. The bolt is never a wall or a puncture.

With the existing ORIGINAL comparisons ||w||2<=C_L||v||E and
||Delta w||2<=C_D||v||E, (14) gives the new usable comparison

    ||sigma_tau[phi,.]||E_t* <= C_E ||g_tau||2,
    C_E=(C_L+C_D)/sqrt(1-epsilon).                  (15)

C_E is evaluated from accepted intervals, not assigned one or replaced
by a ceiling. This closes the norm arrow from a COMPLETE signed L2
source enclosure to the original energy-dual error. It does not prove
decay of g_tau, a nonlinear neighborhood map, a strong source comparison,
or physical admission. The triangle bound from (12)-(13) is still
nonintegrable; the three terms in g_tau must be enclosed together.

## 4. Full budget consequence and the remaining calculation

The precise proposed replacement is now

    ||g_tau||2 <= B_phi tau^(5/4), 0<tau<=.001.      (16)

Neither B_phi nor this rate is proved. Equation (15) would give
A_phi=C_E B_phi in the previous test. The complete angular and original
physical conversion yields

    e_bolt,UV <= D_(1/4) C_E B_phi,
    D_(1/4)=(8.1e-45) [G/(B^2 beta)] 16C(.001)^(1/4). (17)

Every factor outside B_phi is known and evaluated by the frozen receipt.
The hypothetical ceiling

    B_phi <= [1.667e-9-6.447981645e-11]/(D_(1/4) C_E) (18)

allocates all remaining bolt allowance to this UV error and assumes zero
for every unknown debit. It is not a proved margin or bound on B_phi.
The signed center, positive-time evaluation, inverse/Riesz/application
errors and all other coupled columns remain. E32 and Ereg stay paid once.

The failure (8)-(10) and the finite norm comparison (15) are the executed
results. The alpha=1/4 actual-input remainder test ends UNRESOLVED: its
first unsupported implication is the L2 norm of the COMPLETE three-term
expression (13) decaying as (16). No calculation here evaluates that norm.
Qualitative regularity of the mollified input does not provide B_phi.
No arbitrary-source counterexample or finite quantum spectrum is used.

DECISION: reject interpolation of the existing separated strong and
positive-time bounds. Retain the construction and use direct transposition
of the complete signed source (13). Its advantage over that failed
estimator is an established energy comparison without a C5 test norm;
its cancellation/decay still requires proof and numerical constants.

Single next calculation: enclose (13) in L2 for the ACTUAL mixed-Green
input, with its common regulator, full thermal/angular tails and paired
boundary/bolt returns, against (16)-(18). This is PROPOSED, not completed.
An independently proved bound on the same combination would also suffice.
Further tail, degree, cold or isolated-contact optimization is not justified.

The original-space defect target, physical and graph outputs, mass column
kappa_m<=5e-5 and every application error remain unclosed where previously
unproved. The seed quantum-fast screen retains its limited scope. No exact
compatibility, physical matching, validity/noise or stability result follows.

Verification: the ordinary domain and frequency arguments above support
(8)-(15). Arithmetic only evaluates the exponents and known factors in
(15),(17),(18), with source/runtime hashes. It does not verify the missing
continuum bound, evaluate the full action or certify a solution.

Remaining to-do list: complete signed regulator estimate and Zr; every
original coupled/source/graph/physical budget; nonlinear admission and
contraction; exact compatibility; physical matching/opening; semiclassical
validity/noise; all-sector stability; transport, preparation, separate
creation and the full stable semiclassical wormhole objective.
