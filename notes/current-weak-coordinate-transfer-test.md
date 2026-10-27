# Actual-input test of a separated proper-coordinate weak estimator

2026-10-27 HST, E003/S169, second goal turn executing the R1-R8 plan.
Previous goal turn: PROGRESS, in the positive-time domain step only.
The full stable-wormhole objective is active and incomplete. Preserve
the closed audit, all accepted results, P32 and original physical norms.

OUTCOME: a specific proposed transfer of the collar weak estimate has an
INFINITE energy-ball upper estimator. This is proved on the actual fixed
input and admissible physical test fields. It is not a counterexample to
the complete quantum functional, the structured inverse or admission.
The matched regulator bound and signed Zr are still unresolved. This
attempt does not complete the requested signed feasibility milestone.

## 1. The attempted full estimate, not a change of physical problem

Continue phi=chi_N,r-alpha*z^4 and the exact E_t certificate of
current-structured-weak-response-attempt.md. Let u=Vphi, w=Vv. The
positive-time spectral identity and bound in that note remain accepted.
The goal is to bound the complete signed heat remainder after its area
and completed curvature subtractions, before the common angular integral.
For fixed phi a sufficient target remains

    ||D2 r(tau)[u,V(.)]||E_t* <= A_phi tau^(1+alpha), alpha>0. (1)

The operator-uniform version of (1) with only ||phi||E_t on its right is
already ruled out if the nonlocal principal logarithm is omitted: see
current-heat-cutoff-form-obstruction.md. A fixed-input version can use
additional regularity of phi, but requires a proved numerical comparison.
Neither extra regularity nor a free constant is assumed here.

The proposed transfer uses the exact proper-coordinate collar kernel in
current-curved-paired-heat-bridge.md. A possible asymmetric Taylor argument
would distribute five spatial derivatives as three on the fixed direction
and two on the arbitrary energy test. At the coefficient level this
suggests an H3 x H2 estimate. This is only a candidate; it is NOT a proved
bound for the complete metric response. It must first pass the actual
metric-to-proper-coordinate chain rule. This note tests that necessary
step, rather than computing a tangent kernel or a chosen mode response.

Keep the actual paired first-order/adjoint realization, both partner
squares, original period/spin/clock/projector, mass/log core, common
regulator, exact UV logarithm and all contact restorations. The full
thermal and angular sums and global returns remain in the target. Nothing
in the test replaces the regular bolt by a wall or changes a cutoff.

## 2. Exact coordinate variation at the physical seed

Let x be seed proper distance from the outer end and f(x)>0 the seed
inverse time radius on any compact interval away from the bolt. This is
the same coordinate used by the accepted collar argument. For a Weyl
factor b set

    Y_b(x)=integral_0^x exp(b(s))ds,
    x_b=Y_b^-1, f_b(y)=exp[-b(x_b(y))] f(x_b(y)).             (2)

No physical domain is changed by this identification. At b=0 write
I_u(x)=integral_0^x u, I_w=integral_0^x w, I_uw=integral_0^x uw.
Then the accepted exact chain rule gives

    xi_u=-I_u, xi_w=-I_w,
    xi_uw=-I_uw-u xi_w-w xi_u,
    eta_u=D_u f_b=-u f+f' xi_u,
    chi_uw=D_uw f_b=uw f-(u'f+uf')xi_w-(w'f+wf')xi_u
                                  +f''xi_u xi_w+f'xi_uw.   (3)

Differentiate chi_uw twice in x, keeping I_w'=w and I_uw'=uw.
The unique third-test-derivative term is

    partial_x^2 chi_uw = f I_u w''' + lower_test_orders.     (4)

Every other local test derivative has order at most two; the remaining
terms contain I_w or I_uw. All their fixed-input coefficients are smooth
and bounded on a compact interval. Equation (4) is an exact consequence
of (3), not a guess based on differential order. In particular the term
cannot be dropped as a coordinate normalization. The symbolic check
retains all five terms of chi_uw and verifies its coefficient f I_u.

## 3. Why that coefficient is nonzero for the actual fixed input

The accepted trial argument proves phi is not zero. If it were zero,
C V(alpha*z^4)=N ell_b mu_r. The left side is analytic on the open physical
interval, whereas the right is a nonzero compact smooth density that
vanishes on an open subinterval. The analytic identity theorem rules out
that equality. V is the original invertible York return, so u is not zero.

The mixed Green equation has smooth coefficients and a smooth source on
compact interior intervals; its local regularity puts this fixed u in
C-infinity there. Its mass/log outer phase makes I_u well-defined from
the original endpoint. If I_u vanished everywhere in the interior then
u=I_u' would vanish there, contradicting the preceding argument. Therefore
there is a compact interval I away from both physical ends on which

    |f I_u| >= m_I >0.                                    (5)

No numerical value of m_I is needed to decide boundedness. Its positivity
will force divergence for any value. It is not assigned a value for an
admission estimate. This argument uses the actual phi, not a convenient
replacement input. If an earlier global H3 comparison for phi is unknown,
it remains unknown; the local test below does not assume that comparison.

## 4. Admissible full-equation test and the unbounded comparison

Use the physical areal coordinate r on I and write F for the actual seed
lapse. Set W=r^2 F'-2rF. It is analytic and not identically zero; restrict
I if necessary so a smooth compact rho can satisfy integral W rho !=0.
Choose a nonzero smooth compact chi in I and let x(r) be the same proper
coordinate as above. For integer n define

    c_n=integral W chi cos(n x(r))dr / integral W rho dr,
    xi_n=chi cos(n x(r))-c_n rho,
    b_n=-(r^2 F)^-1 integral_(left endpoint)^r W(s)xi_n(s)ds,
    k_n=diag(-b_n-2xi_n,b_n,xi_n,xi_n).                    (6)

Both original TT constraint terms are kept:
(r^2 F b_n)'=-(r^2 F'-2rF)xi_n. The exact zero moment makes b_n vanish
on both sides of I. Thus the physical perturbation vanishes near the
original paired boundary, all lift transitions and the regular bolt.
It belongs to the actual smooth core with zero mass/log coefficients;
it does not replace the full mass/log domain by a smooth surrogate.
Its canonical w_n=Vxi_n uses the original global York realization.

Since x'(r) never vanishes on I, repeated integration by parts gives
c_n=O(n^-m) for every fixed m. The accepted complete common-regulator
local symbol and York calculation give, on this fixed interval,

    ||xi_n||E_t <= C_I n^2 sqrt(log n),
    w_n=-(3/2)chi cos(nx)+one-order-lower terms.             (7)

Here the second statement is the local symbol statement, including
derivatives through order three. The lower-order global York return is
retained, not declared zero. These are the same actual finite-seed facts
used by current-continuum-energy-dual-estimator-test.md; no continuum
conclusion is drawn from samples of n. The index n is a test-field
frequency, not a truncation of quantum angular or thermal labels.

Insert (7) into (4). By (5), the leading oscillatory term has L2 size
bounded below by a positive constant times n^3. Integration of sin^2
against the fixed nonzero amplitude gives a nonzero limit after division
by n^6. All other terms are O(n^2), including the primitive returns.
Hence, for sufficiently large n and suitable positive c_I,C_I,

    ||chi_(u,w_n)||H2(I) >= c_I n^3,
    ||D2 f_b[u,V(xi_n/||xi_n||E_t)]||H2(I)
                     >=(c_I/C_I) n/sqrt(log n) -> infinity. (8)

Thus there is NO finite constant for

    ||D2 f_b[u,Vv]||H2(I) <= C(u) ||v||E_t.                (9)

The infinite norm is an executed analytic counterexample to (9), not an
unknown or excessively large finite constant. Changing normalization or
choosing a smaller fixed heat cutoff cannot repair this derivative loss.

## 5. Consequence for the complete regulator and output budget

Let R_tau(f) denote the coefficient-coordinate subtracted trace. Even
granting finite candidate coefficient estimates, the separated chain-rule
upper bound would contain

    |D2 R_tau[eta_u,eta_w]| + |D R_tau[chi_uw]|
      <= C_32(tau)||eta_u||H3||eta_w||H2
                            + C_2(tau)||chi_uw||H2.        (10)

For any positive C_2(tau), its second term has infinite supremum on the
actual E_t unit ball by (8). If C_2 or the first term is unestablished,
that is a separate preceding unknown, not a cure for (9). No vanishing
of D R_tau or C_2 is proved. Terms from the moving coordinate cutoffs,
measure and local restorations also remain in the complete target.

This estimator therefore supplies NO finite A_phi in (1). Its contribution
to a certified regulator error cannot be charged to the original bolt
ledger. The accepted coefficient K3=8.08582851697e-45 multiplying that
third-functional error does not make an infinite upper estimator finite.
No physical/source/graph margin, mass-column bound or other coupled entry
improves. The former positive-time bound, E32 and Ereg are unchanged.

The exact full response need not inherit (8): the coordinate terms in
D2 R, D R, the moving cutoffs and the trace measure can cancel BEFORE
their norms are taken. Indeed the previous fixed-Hilbert positive-time
identity already gives a finite answer on this same test ball. This
confirms the scope of the failure: separated coordinate factors, not
the actual functional. A repair must recover at least this ONE derivative
of the test through an explicit combined identity and estimate, in
addition to proving the ultraviolet heat-time gain.

## 6. Method decision and one next calculation

Reject the proposed separated proper-coordinate H3 x H2 transfer. Retain
the already defined FIXED-HILBERT spectral response, where the test is a
bounded multiplication operator, and combine its exact local subtractions
before estimation. This avoids introducing (9). It does not by itself
bound the remaining diagonal/contact difference or evaluate the UV log.

Compare the specified alternative: a fixed-coordinate boundary
log-determinant/resolvent representation also avoids the transport term
w' I_u. Its common first-order domain and bounded mass vertices are
available. However their operator-norm bounds do not give a summable
differentiated trace over the complete angular tower; the matched local
and boundary remainder still needs an energy-dual estimate. No available
theorem or executed calculation supplies that estimate. This comparison
does not justify changing the physical construction.

The single next calculation is the full fixed-Hilbert spectral Hessian
minus its completed local terms, with a quantitative E_t* bound on the
actual phi and common regulator, including boundary/bolt returns. It is
PROPOSED, not completed here. The current result rules out one attempted
regularity transfer but does not supply the requested A_phi, positive
alpha, signed Zr, or consequential admission bound. There is still no
certified controlled solution or construction obstruction.

Remaining to-do list: matched regulator estimate and signed response;
all original coupled/source/graph/physical budgets; nonlinear admission
and contraction; exact compatibility; physical matching and opening;
semiclassical validity/noise; all-sector stability; transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
