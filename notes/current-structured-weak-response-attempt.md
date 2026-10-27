# Structured weak response: positive-time extension and failed absolute majorant

2026-10-27 HST, E003/S169, execution of R1-R3 of the approved research
plan. The preceding turn delivered the plan; this is the first goal turn
of its execution. Full stable-wormhole objective ACTIVE and incomplete.
Closed audit, accepted results, action, state, physical norms and P32 stay
unchanged. This is NOT completion of the requested signed milestone.

OUTCOME: the complete positive-time geometric heat Hessian has an explicit
energy-space bound on the original domain. The absolute spectral estimator
defined below nevertheless fails regulator removal: its integrated upper
expression diverges quadratically. Its finite positive-time contribution
is also very large. The actual commonly renormalized signed response is
still UNRESOLVED, not proved unbounded. No original admission gate closes.

## 1. Composition, spaces and the first missing implication

The fixed input remains phi=chi_N,r-alpha*z^4, r=1e-110,
chi_N,r=N V^-1 C^-1(ell_b mu_r), alpha=ell_b=(1+h)/h,
N=|W_action| beta. The original smoothed observation and its already paid
error are unchanged. E_t is the completion of the accepted paired form
core, with norm squared ell_t/N and

    (1-eps)||v||E^2 <= ||v||E_t^2 <= (1+eps)||v||E^2,
    eps < 1.227e-74.                                           (1)

The core has the one-sided mass cubic and fixed finite quartic log phase,
the original lift cutoffs, fixed clocks/period/spin/projector and smooth
regular-disk matching at the bolt. Its energy completion is not declared
to be a weighted C5 nonlinear metric family. York retains J(0)=0 and
(Q^2 J'/z^4)(L)=0. Both original equations enter the accepted joint output
and the definition of phi; no new scalar source or boundary is introduced.

The immediate composition is

    (phi,v) in E_t x E_t --(V,V)--> (u,w)
      -- complete commonly restored massive Hessian --> R,
    v -> L_t[phi,v] in E_t*,
    Riesz and correlated F_t,F_t^2,F_t^3 transfers -> Zr,
    original prefactor, E32, Ereg and other debits -> bolt certificate.

Only a fixed-positive-time version of the second arrow is bounded below.
The full common-regulator limit, its signed evaluation, Riesz/application
errors and original strong source/graph comparisons are still missing.
The original global geometric Banach X is not completed by this result;
its cold/outer/bolt/radius/boundary requirements remain mandatory.

For u=Vphi and w=Vv the accepted energy comparisons give

    ||u||2 <= C_L ||phi||E_t/sqrt(1-eps),
    ||Delta u||2 <= C_D ||phi||E_t/sqrt(1-eps),
    ||u||infinity <= C_inf ||phi||E_t/sqrt(1-eps),
    Delta u=(Q u')',  C_inf <7.477e50.                       (2)

The same bounds hold for w. These maps include the original regular-bolt
flux and radius/temperature conversion costs; C_inf is not set to one.
The input norm ||phi||E_t is not numerically enclosed in this attempt.
No conditional ceiling is inserted for it.

The first missing implication tested here is extension of the differentiated
positive-time trace to this full test ball without W3,infinity or C5 test
norms. Common-regulator removal is a further implication, tested in section 5.
The boundary realization itself is reused, not inferred from formal symmetry.

## 2. Actual-domain spectral identity, including the coincidence vertex

Let A=i J D0 on the already proved paired first-order realization. The
accepted identities D0*= -J D0 J and (J D0)*= -J D0 make A self-adjoint
with A^2=L0=D0*D0. This does not declare the unpaired D0 self-adjoint or
choose another extension. On an admissible real Weyl core,

    A_b=exp(-b/2) A exp(-b/2),  L_b=A_b^2.

Scalar multiplication preserves the original first-order domain. The
square always has its composed domain. At positive heat time tau the
compact finite-seed realization has trace-class heat kernel; its spectrum
contains the full thermal tower and all finite-boundary/bolt returns.
Both partners give the same physical trace by the accepted J relation.
Use the physical-half convention throughout, with no extra doubling.

On the core, for a fixed scalar u,

    delta_u L = -A u A - (u L+L u)/2,
    delta_u Tr exp(-tau L)
      = -2 tau partial_tau Tr[u exp(-tau L)].                (3)

This follows from differentiating A_b on the common first-order domain
and cyclicity at positive time, interpreted through the closed forms.
Differentiate (3) in w. Let A e_j=a_j e_j, lambda_j=a_j^2 and
g_tau(lambda)=lambda exp(-tau lambda). The result is the complete Hessian

    D2 K(tau)[u,w] = sum_(j,k) q_tau(a_j,a_k) u_jk w_kj,
    q_tau(a,b) = -tau (a+b)^2
                  [g_tau(a^2)-g_tau(b^2)]/(a^2-b^2).         (4)

At equal squared eigenvalues use g_tau'; if a+b=0 the entry is zero.
At a=b this gives 4 tau a^2(tau a^2-1)exp(-tau a^2), agreeing with
the exact constant-Weyl second variation. In particular (4) retains the
second vertex, not just the bubble. Deriving it from (3) avoids differentiating
the arbitrary test field in space. No radial integration by parts discards
a boundary flux. All cutoff and bolt effects remain in the actual spectrum.

To justify the infinite sum and its core extension, use the bound proved
next. Finite spectral projections can justify intermediate algebra, but
their limits are taken with an absolute summable majorant; no finite
selection of physical modes defines (4). This proof requires only the
accepted first-order domain and positive-time trace properties. Smooth
boundary heat-coefficient theorems are not being imported into the mass/log
family. The literature's smooth-smearing conformal identities, for example
[Esposito-Gilkey-Kirsten](https://arxiv.org/abs/math/0510156), are background,
not a theorem supplying this energy-dual or regulator estimate.

## 3. An explicit whole-spectrum bound without test-field derivatives

For all real a,b and tau>0,

    |q_tau(a,b)| <=32[exp(-tau a^2/2)+exp(-tau b^2/2)].       (5)

Here is a direct scalar proof, covering all frequencies. Put lambda=a^2,
mu=b^2 and order mu>=lambda. If mu>=2lambda, use
(a+b)^2<=2(lambda+mu) and (lambda+mu)/(mu-lambda)<=3. Then the left side
is <=6 tau(lambda exp(-tau lambda)+mu exp(-tau mu));
x exp(-x)<=2 exp(-x/2)/e proves (5), with constant 12/e<32.
If lambda<=mu<=2lambda, the mean value theorem gives
|g_tau(lambda)-g_tau(mu)|/|lambda-mu|
<= (1+2tau lambda)exp(-tau lambda). Thus the left side is at most
6x(1+2x)exp(-x), x=tau lambda. Its ratio to exp(-x/2) is bounded by
12/e+192/e^2<32. The zero/equal cases follow by continuity.

For bounded real scalar multipliers, row Cauchy-Schwarz gives
sum_k |u_jk w_kj| <=||u||infinity ||w||infinity.
Apply this to both terms in (5), summing every j. Consequently

    |D2 K(tau)[u,w]|
      <=64 ||u||infinity ||w||infinity K(tau/2).              (6)

The same bound controls differences of core approximants, by bilinearity
and (2). It proves a unique continuous positive-time extension on E_t x E_t.
It does not assert a nonlinear metric family for every energy vector.
The needed derivative loss is removed at FIXED positive tau; no uniform
small-time regularity gain has yet followed.

Write k(tau)=gamma K(tau)/(2pi), the original time-reduced physical trace.
Then (6) and (2) give

    |D2 k(tau)[Vphi,Vv]|
      <= C_heat k(tau/2)||phi||E_t ||v||E_t,
    C_heat=64 C_inf^2/(1-eps) <3.578e103.                    (7)

The huge factor is the actual previously established global supremum
comparison, including 1/gamma; it is not an unknown constant.

## 4. Complete positive-time subtracted remainder and its budget contribution

Retain the exact common subtraction from
current-canonical-heat-restoration-assembly.md. On the actual phase the
reduced heat Hessian has area subtraction
2 int u w dz/(pi tau), zero a1 Hessian, and the completed a2 subtraction

    -tau/(480pi) int[4R^2 u w+8R(u Delta w+w Delta u)
                                      +8 Delta u Delta w]dz. (8)

This includes the paired divergence/contact completion. Its continuation
to E_t follows from (2), not from an unproved pointwise fourth derivative.
The seed polynomial supplies |R|=|Q''|<20 on 0<=z<=L<1; the arithmetic
checks the coefficient-sum bound, retaining every u, eps and tau parameter.

Define

    C0=2 C_L^2/[pi(1-eps)],
    C2=[4(20)^2 C_L^2+16(20)C_L C_D+8 C_D^2]/[480pi(1-eps)].

If r is the time-reduced remainder after these exact subtractions, then

    |D2 r(tau)[Vphi,Vv]|
      <=[C_heat k(tau/2)+C0/tau+C2 tau]||phi||E_t||v||E_t.   (9)

The full massive Hessian remainder in normalized energy units is

    -(G/(B^2 beta)) int Theta(tau) D2 r(tau) d tau/tau,
    Theta(tau)=sum_(n>=1) 2(q+2n) exp[-n(n+q)tau/C].          (10)

The beta divisor is necessary. The same common regulator applies to the
whole expression. The extracted UV logarithm and finite a0/a2 restorations
are kept in the original action; (10) does not delete or reassign them.
LLL, Weyl, classical/scalar and other restored terms remain separate exact
parts of L_t and are not declared evaluated by (10).

For tau>=tau0=.001 use the accepted complete seed trace
k(.00025)<1521, hence k(tau/2)<1521 by positivity. Keep the entire angular
sum: Theta(tau0)<=4C/tau0 and
Theta(tau)<=Theta(tau0)exp[-324(tau-tau0)], since M1^2>324.
These yield a finite bound for ALL times tau>=tau0:

    ||L_mass,positive[phi,.]||E_t* <= B_positive ||phi||E_t,
    B_positive=(G/(B^2 beta))*(4C/tau0)/324
                  *[1521 C_heat/tau0+C0/tau0^2+C2]
                <2.362e79.                                (11)

This includes the thermal tower, the whole angular tower, the local
subtractions and global returns. It is an upper bound for this specified
positive-time remainder, not its measured value or the complete action.

Its explicit role if used as an unevaluated third-functional error is

    e_bolt,positive <= K3 B_positive ||phi||E_t,
    K3=pref_t f_t^3, K3<8.1e-45,
    K3 B_positive <1.910e35.                               (12)

Even assigning the ENTIRE bolt allowance to (12), this estimator's
arithmetic would require ||phi||E_t<8.731e-45. That is a CONDITIONAL
input threshold, not an upper bound on phi or an available budget. The
exact sufficient threshold is evaluated as 1.667e-9/(K3 B_positive);
the displayed rounded value is descriptive, not a sufficient rounded guard.
For a unit energy input the estimator allowance is about 1.1454e44 times
the entire bolt allowance. That comparison is to the fully composed error,
not to an intermediate constant. It is not a lower bound on the action,
nor a test of the unknown actual phi norm. No debit is paid or margin gained.

## 5. Executed ultraviolet test of this estimator

The scalar majorant in (9) is finite for each tau>0. It is NOT integrable
in (10) at zero. On the actual seed, the accepted full heat/area asymptotics
give k(tau/2)~L/(pi tau), .99<L<1, and Theta(tau)~2C/tau. Therefore the
positive upper expression already from the first term is asymptotic to

    D_L/tau^2 in the d tau/tau integrand,
    D_L=(G/(B^2 beta))*C_heat*2CL/pi,
    D_1 approximately 8.0058370615e71, D_L=L D_1.            (13)

Its integral from a common heat cutoff eta to tau0 diverges as
D_L/(2 eta^2). This is a VERIFIED FAILURE of the specified absolute
majorant, not an infinite value of the subtracted functional. The other
positive terms in (9) cannot repair that estimator.

The cancellation needed is now explicit: retain the signed difference
of (4) and (8) before norms and obtain, uniformly on the actual E_t ball,

    ||D2 r(tau)[Vphi,V(.)]||E_t* <= A_phi tau^(1+alpha),
    alpha>0,                                               (14)

or another majorant integrable after the exact Theta(tau)/tau weight.
Relative to the tau^-2 whole-action majorant, a power-law certificate
needs more than TWO powers of heat-time improvement. A finite reduction
of C_heat alone cannot make its ultraviolet integral finite. The exact
UV logarithm must remain in the restored principal action; (14) is for
the matched remainder, not the raw cutoff action.

## 6. Method decision, next calculation and preserved gates

R2 has advanced: (4)-(7) give a full-spectrum, derivative-free positive-time
weak representation and its actual energy comparison. The completed
contacts also extend by (8). R2-R3 have NOT supplied the regulator-uniform
matched estimate (14), its numerical constant, or the signed Zr. The
overall requested signed feasibility milestone remains incomplete.

Retain the physical construction and change the estimation step from
absolute spectral/heat terms to the MATCHED signed spectral kernel minus
its completed local terms. The next calculation is (14) for the actual
phi, with the boundary-return remainder kept in the same expression.
The smoothing of the observation supplies no rate for (14) by itself.
Do not optimize C_inf, the cold/scalar parts, P32 degree or transport tail:
none fixes the nonintegrable exponent of this estimator.

The specified alternative boundary log-determinant parametrix faces the
same obstacle: its retained principal logarithm must be separated from a
quantitatively controlled differentiated boundary remainder in E_t*.
A bare trace-class whole-action argument was already rejected. The new
positive-time extension supports keeping the weak spectral route for the
next attempt; it does not establish superiority of its UV remainder or
justify changing construction.

Unknown: ||phi||E_t as a usable numerical input, A_phi/alpha, signed Zr,
the other residual errors, all unproved original source/graph/boundary
comparisons and remaining coupled entries. The existing mass-column
target, limited quantum-fast screen and all conditional ceilings retain
their scopes. Neither seed estimates nor this extension admit a nonlinear
X ball. No earlier physical matching or compatibility theorem is supplied.

Verification: the ordinary infinite-spectrum and domain argument above
is the evidence for (4)-(14). The small accompanying Arb calculation
checks its scalar guards, loads accepted interval inputs and evaluates
the complete error composition (11)-(13). It does not numerically solve
the orbit spectrum, verify the missing continuum theorem or enclose the
signed response. No Lean, simulation or finite-mode physical surrogate.

Remaining to-do list: matched regulator-uniform remainder and signed Zr;
complete original coupled/source/graph/physical budgets; nonlinear admission
and contraction; exact compatibility; physical matching and opening;
semiclassical validity/noise; all-sector stability; transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
