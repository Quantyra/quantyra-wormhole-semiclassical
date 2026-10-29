# First inverse residual and the uniform source estimator

2026-10-28 HST, E003/S169. First turn of the NEW execution goal. The
preceding plan supplied no bound; earlier blocked runs are historical.
Preserve the full stable semiclassical wormhole objective, P32, closed
audit, original physical norms and all accepted results.

This tests one complete first-inverse residual estimator. It does not
evaluate the full signed Z_r. A finite excessive upper report and a
lower limit on this CLASS OF UPPER REPORTS reject taking the source norm
separately. They do not lower-bound the actual residual or correction.

## Actual source, trial, and full-space residual

Use H_mu=L2((0,L),z^-4 dz), e_t=ell_t/N on the original paired E_t,
N=|W_action| beta, and e_s=a/N on the scalar energy E_s. The fixed
mass/log phase, complementary first-order adjoint and composed squares,
clock/period/spin/scale, original cutoff transitions and regular disk bolt
are unchanged. York retains J(0)=0 and (Q^2 J'/z^4)(L)=0.

Let t_hat be the analytic Green trial accepted at 0459e57. The complete
adjoint estimates of current-complete-bolt-energy-bound.md imply

    g_s(t)(w)=<t,(V*)^-1 bbar* w>H_mu/beta,
    ||g_s(t)||E_s* <= K_s ||t||H_mu, K_s=K_V K_b/beta.

These are the complete local mixed rows, including the Ward/anomaly and
return terms; the canonical nonlocal determinant has zero scalar/mixed
row in the fixed chart. No independent quantum stress is chosen. Let y
be the scalar Riesz solution e_s(y,w)=g_s(t)(w). Its induced first source is

    g(t)(v)=b(y,v)/N
          =<t,(V*)^-1 Bbar v>H_mu/beta,
    ||g(t)||E_t* <= K_g ||t||H_mu,
    K_g=K_s sqrt(q)/sqrt(1-epsilon), q=.338.                 (1)

The identities follow first on the original form core and extend by the
quoted bounds. In particular they do not require V^-1 t to be a C5
metric path, or an unproved unfiltered H_mu-to-E_t embedding. The source
exists, but neither y nor its full numerical source functional has been
computed on this t_hat. That is an earlier implementation gap than the
independent quantum inverse application.

The first inverse is e_t(x,v)=g(t_hat)(v), namely the filtered extension
of T_t V^-1 t_hat. Choose the explicit nonzero trial

    x_hat=c psi, c=alpha q, psi=z^4, alpha=(1+h)/h.

For psi, J=0, Upsi=psi/2 and Vpsi=-psi. This is the actual admissible
quartic with the original outer and bolt data. The retained full-energy
bound vcap>=||psi||E_t from current-signed-adjoint-trial-attempt.md includes
the entire commonly renormalized quantum form, exact UV logarithm,
complete angular/thermal towers, contacts and finite-domain returns.
It is not a local-only energy. Do not transfer its separate quartic
quantum-action certificate to the Green input.

The COMPLETE signed residual is r(v)=g(t_hat)(v)-e_t(x_hat,v), for all
v in E_t. The fallback being tested takes its two terms separately:

    ||r||E_t* <= R_upper=K_g ||t_hat||H_mu+|c| vcap.          (2)

This upper bound needs no unevaluated quantum matrix element set to zero:
the second functional's exact dual norm is ||x_hat||E_t. It supplies no
signed center. Parameter/rounding errors enter the interval expressions;
the trial and input are analytic functions, not sampled interpolants.
Unknown application errors for any later implemented solver remain unknown.

For F_t=(I+T_t)^-1 T_t and f_t=q_t/(1+q_t), q_t=q/(1-epsilon), coercivity
and the two remaining filters put (2) in ORIGINAL bolt-output units:

    E_first <= K_first R_upper, K_first=pref_t f_t^2.        (3)

The previously accepted input replacement debit <7.432e-12 is charged
once in the complete certificate; it is not included again in (2).
The correlated j_r term, quartic subtraction, later inverse errors and
all other outputs are still required. No available margin is assumed.

## Whole-interval norm of the actual computational input

Write s=L-z, SP=1e-16 and

    t_hat(s)=chi_SP(s) At [Q_l(1+a s/gamma)-Q_n(1+a s/gamma)],

with exactly the accepted correlated parameters l>=1, n>l and At>0.
The existing cutoff is one through SP and zero beyond 2SP, with values
in [0,1]. Standard real Q is used, not Olver's differently normalized Q.
The physical interval and its cutoff remain unchanged.

The [Legendre equation](https://dlmf.nist.gov/14.2.E1) gives

    [(w^2-1)(Q_n' Q_l-Q_l' Q_n)]'
                        =[n(n+1)-l(l+1)] Q_n Q_l.

At w=1 the standard Q expansion has common leading -log(w-1)/2 and
constant log(2)/2-EulerGamma-digamma(k+1); see
[DLMF 14.8.9](https://dlmf.nist.gov/14.8.E9) with the normalization in
[14.3.10](https://dlmf.nist.gov/14.3.E10). At infinity the boundary term
vanishes for these positive degrees. Integrating and taking the coincident
degree limit therefore proves

    I(n,l)=int_1^infinity Q_n Q_l dw
          =[digamma(n+1)-digamma(l+1)]/[n(n+1)-l(l+1)],
    I(k,k)=trigamma(k+1)/(2k+1),
    J=int_1^infinity (Q_l-Q_n)^2 dw=I(l,l)+I(n,n)-2I(n,l). (4)

The coincident limit is justified by the positive integral representation
and its logarithmic endpoint/integrable infinity bounds on a compact
degree interval. No matter mode is being selected in (4).

For x>=1, Q_1(1+x)<=1/(3x^2), by its positive inverse-power series.
Monotonicity in degree of the same positive integral gives
0<=Q_l-Q_n<=Q_1. With x0=a SP/gamma>1,

    tail=int_(x0)^infinity (Q_l(1+x)-Q_n(1+x))^2 dx
                          <=1/(27 x0^3).

Since the cutoff is one before SP and at most one afterwards, and the
ORIGINAL measure has L^-4<=z^-4<=(L-2SP)^-4 on its support,

    At L^-2 sqrt((gamma/a)(J-1/(27x0^3)))
      <= ||t_hat||H_mu
      <= At (L-2SP)^-2 sqrt((gamma/a)J).                     (5)

The infinite comparison integral is only a bound for the computational
profile; it does not extend the physical geometry. Equation (5) covers
the whole input, including its cutoff transition and regular bolt limit.
Its purpose is to test (2), not to optimize the closed input contract.

## Why sharpening the uniform constant cannot rescue this estimator

Let C_g be ANY valid norm bound for the map g:H_mu->E_t*. Test it on
t=Vpsi=-z^4, an actual admissible field. The existing two-scalar-test
variational lower matrix in current-quartic-ritz-rank-test.md gives
B[psi,psi]/N>=B_low,22>0. Its finite matrix is a LOWER variational bound
for the full scalar inverse, not a replacement for that inverse. Then

    C_g >= ||g(Vpsi)||E_t*/||Vpsi||H_mu
        >= B_low,22/(vcap sqrt(L^5/5)) = C_floor.           (6)

Both the numerator and complete quantum energy in the denominator retain
their actual domains. This tests a UNIFORM source-map constant. It does
not lower-bound g(t_hat) on the special Green input.

Every separated report of the form C_g||t_hat||+||x_hat||, even with an
exact best uniform constant and arbitrary trial coefficient, is therefore
at least C_floor times the LOWER endpoint of (5). Its original-output
report is at least K_first C_floor H_lower. If this exceeds 1.667e-9,
no better uniform constant or quartic coefficient can certify that report.
The zero trial is the no-data baseline for (2); nonzero c only adds to
this particular upper expression. A signed or input-specific calculation
can escape this estimator limitation.

## Quantitative alternative and execution decision

Compute the scalar source y before estimating its norm. For an admissible
y_hat define its complete normalized scalar residual

    rho_s(w)=g_s(t_hat)(w)-e_s(y_hat,w),
    ||y-y_hat||E_s<=||rho_s||E_s*.

The original coupling bound gives

    ||g(t_hat)-b(y_hat,.)/N||E_t*
                       <=sqrt(q)/sqrt(1-epsilon) ||rho_s||E_s*,
    E_source <= K_scalar ||rho_s||E_s*,
    K_scalar=K_first sqrt(q)/sqrt(1-epsilon).               (7)

This step contains the full known scalar and mixed operators and actual
boundary returns; it does not need evaluation of the independent quantum
shape determinant. A scalar approximation and its residual are NOT yet
available. The subsequent quantum residual must retain
b(y_hat,v)/N-e_t(x_hat,v) before taking norms. Formula (7) is its explicit
budget role, not a claim that completing the scalar solve closes Z_r.

Execute (4)-(7) with 512-bit outward intervals. Test the current report
and the best-uniform-constant report against the WHOLE bolt allowance,
granting all other debits zero only to reject this estimator. Record the
conditional scalar-residual ceiling 1.667e-9/K_scalar, not an available
allocation. Preserve every uncertainty and no numerical application error
is silently zeroed.

Decision at this test: if those report floors fail, abandon the separated
uniform-source estimator, retain the construction, and compute the actual
coupled scalar source followed by its signed quantum residual. Do not
optimize the quartic column, Green input, tails, polynomial degree or the
uniform coefficient. The next calculation is the complete scalar-source
solve with a certified rho_s on the actual domain, charged through (7).
The actual signed Z_r and full coupled certificate remain unresolved.

Remaining to-do list: actual source and signed quantum residual applications;
all coupled admission/source/graph/physical/boundary/mass requirements;
nonlinear admission/contraction; exact compatibility; physical matching/
opening; semiclassical validity/noise; all-sector stability; the full
stable semiclassical wormhole objective.
