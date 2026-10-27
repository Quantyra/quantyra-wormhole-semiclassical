# Finite-rank inverse test on the actual energy domain

2026-10-27 HST, E003/S169. The resumed goal is ACTIVE and incomplete.
Preserve the full stable semiclassical wormhole objective, closed audit,
accepted results, P32, original norms and the sufficient transport depth.
This is a test of a specified numerical estimator, NOT an evaluation of G_t.

The previous response supplied a requested goal statement and no research
evidence. This attempt identifies a genuine failure of one possible
implementation: approximating a bounded inverse uniformly on the whole
energy space by a finite-rank matrix before composing with its source.
The preceding spatial/Neumann argument did not require that implementation;
its accepted conditional bounds are not contradicted or replaced.

## Actual spaces and complete operators

Use the accepted paired first-order realization, its complementary adjoint
and composed squares. Keep the mass cubic and finite logarithmic phase,
fixed period, clock, spin and renormalization scale, original mass-lift
cutoff transitions, and regular-disk bolt. York has J(0)=0 and
(Q^2 J'/z^4)(L)=0. Let E_t be the actual retained positive TT form
completion, with e_t=ell_t/N and N=|W_action| beta. Here and below all
operators written as E_t -> E_t use THIS Riesz identification.

    T_t=ell_t^-1 B,  0<=T_t<=q_t=.338/(1-epsilon),
    L_t=I+T_t,      F_t=L_t^-1 T_t,
    G_t=J0 F_t^3+J1 F_t^2 : E_t -> R.

Both original equations give J0,J1 as in current-spatial-inverse-output-
test.md (1). Their output is the original bolt radius observation, including
ell_b=(1+h)/h. The full common-renormalized ell_t contains the exact fourth
ultraviolet logarithm, all angular and thermal sums, coincidence/density
contacts, cutoff transitions and finite/regular-bolt returns. They are not
replaced by a matrix or discarded in the arguments below.

E_t is infinite-dimensional. For example the accepted genuine TT core
contains infinitely many independent interior-supported directions. In the
York construction choose compactly supported J with its one return moment
zero; xi=(QJ)'-z^4 integral_0^z (Q'/z^4)'J then has compact support as well.
The scalar/shape reconstruction has zero exterior phases and is regular at
the actual bolt. An infinite-dimensional test subspace inside the actual
domain does not replace the domain or its nonzero mass/log sectors.

The seed form realization and 0<=T_t<=q_t justify boundedness and inversion
of L_t. This proof invokes no smooth-double theorem, quantitative graph
embedding, or unproved extension to a nonlinear neighborhood.

## A uniform inverse approximation is the wrong intermediate target

Let R_n:E_t->E_t be ANY bounded finite-rank linear map. This includes a
finite Galerkin inverse followed by finite-basis reconstruction, when its
input functionals are bounded. If they are not bounded, an operator-norm
certificate is not defined in the first place. There is a unit vector v
in ker R_n. For any bounded invertible L on this space,

    ||(L^-1-R_n)v||=||L^-1 v|| >= 1/||L||.

Consequently the ACTUAL retained operator obeys, for every n,

    ||L_t^-1-R_n|| >= 1/(1+q_t).                       (1)

The complete geodesic preconditioner in the accepted spatial method obeys
||L0||<=1+q_t+e0, e0=(.0022+epsilon)/(1-epsilon), so also

    ||L0^-1-R_n|| >= 1/(1+q_t+e0).                    (2)

These bounds are approximately .747384 and .746157, respectively. No
increase of matrix rank or floating-point precision makes either error
small in this topology. This is an ordinary infinite-dimensional proof;
arithmetic only encloses its constants. It is not a measured condition
number or a lower bound on the error of F_t, G_t, or an affine correction.
In particular a restricted composed map may be approximable even though
the inverse on the whole energy space is not. Qualitative smoothing alone
still supplies no effective approximation constant for that composed map.

For a quantitative comparison test this SPECIFIED separated estimator:

    F_hat=R_n T_t,
    a_sep=q_t alpha_n, alpha_n>=||L_t^-1-R_n||.

It estimates ||F_t-F_hat|| by multiplying separate norm bounds. Equation
(1) forces a_sep>=f_t=q_t/(1+q_t), about .252616. Thus this expression
cannot certify the previously proposed 1e-30 scale for a uniform F
application: it misses by at least 2.52616e29. That scale remains a
sufficient target, not an attained error or a necessary existence test.

Use the accepted physical weights W0,W1, including the original radius
conversion and pref_t, from current-spatial-inverse-output-test.md. Its
telescoping error expression is

    Gamma(z)=W0(3z+3z^2+z^3)+W1(2z+z^2).

For the separated estimator z=a_sep/f_t>=1, so its reported physical
upper expression is at least Gamma(1)=7W0+3W1. Compare that expression
with the ENTIRE 1.667e-9 bolt allowance, granting every other cost zero.
This is a floor on an EXCESSIVE UPPER ESTIMATOR, not on the physical error.
No claim is made that the retained Neumann implementation must use this
separation. Equations (1)-(2) do not reject analytic infinite-dimensional
preconditioners, finite-rank corrections to them, or range-restricted solves.

## Test a specified alternative against the same obstacle

Use conforming Ritz approximation of the COMPOSED map, not a uniform
approximation of L_t^-1. Choose a finite-dimensional V_n contained in the
actual paired form domain, with the original phase lifts retained when
needed. For each u in E_t define F_n u in V_n by

    <L_t F_n u,v>_t=<T_t u,v>_t for every v in V_n.     (3)

L_t>=I proves (3) is uniquely defined. It is a finite-rank map, but its
complete residual is

    R_full,n=T_t-L_t F_n,
    F_t-F_n=L_t^-1 R_full,n,
    ||F_t-F_n|| <= ||R_full,n||.                       (4)

It is the FULL residual on E_t, not the finite residual in V_n. Galerkin
orthogonality sets the latter to zero without controlling the former.
Writing P_n for the orthogonal E_t projection gives
R_full,n=(I-P_n)(T_t-L_t F_n). Thus a norm enclosure of this omitted
component, with quadrature/application errors, is the first unsupported
quantitative implication in this alternative.

There is also an exact range approximation interpretation. If Q_n is the
L_t-orthogonal Ritz projection, F_n=Q_n F_t. For any E_t-orthogonal P_n
onto V_n, coercivity, orthogonality and boundedness give

    ||F_t-F_n|| <= sqrt(1+q_t) ||(I-P_n)F_t||.         (5)

Proof: for each u, Ritz minimizes the L_t energy error over V_n; compare
it with P_n F_t u, then use ||w||_t<=||w||_L and
||w||_L<=sqrt(1+q_t)||w||_t. No graph constant is silently assumed.

The lower bound (1) does NOT apply to (4)-(5): their input has already
passed through T_t, and they do not approximate the identity component of
the inverse on arbitrary inputs. This addresses the specific topology
obstacle. To reach 1e-30 by (5) would require the actual range tail below
1e-30/sqrt(1+q_t), approximately 8.64514e-31. That is UNKNOWN. The
coefficient sqrt(1+q_t) is proved; the tail is not a conditional ceiling
renamed as a known constant. No theorem inspected here provides its
numerical value on the actual mass/log domain.

Implementation stops BEFORE assembly of a claimed certified matrix:
the full ell_t action and the complement in (4) have not been enclosed.
The existing evaluator still covers compact directions and only two
angular levels. Finite trial residuals, the known .0022 comparison, and
the accepted small spatial/Neumann errors do not fill this missing map.
The alternative was tested analytically against the finite-rank obstruction;
its complete signed numerical application was NOT executed.

## Outcome and decision

Executed: a quantified failure of the separated finite-rank inverse
estimator, valid on the actual infinite-dimensional energy space, and a
proof that range-restricted Ritz certification avoids that particular
obstacle. This rules out a solver design; it closes NO physical allowance
and supplies NO new signed value of G_t. The full requested milestone
remains unresolved. No new claim about construction feasibility follows.

Decision: retain the construction and original P32, reject separate uniform
finite-rank inverse certification, and require complete range-residual or
direct joint-adjoint certification before treating numerical output as
evidence. The next single calculation is an enclosure of R_full,n in (4)
for an explicit actual-domain Ritz trial, including the complete signed
quantum action and its full complement. It is proposed, not executed.
If that enclosure is unavailable, publishing more matrix precision or
conditional error budgets would not advance the signed certificate.

Every original defect, graph, physical, mass-column and independent-phase
entry remains unknown where previously unknown. No remaining margin is
declared available. Preserve the limited two-fast quantum screen, all
nonlinear admission/contraction, exact compatibility, matching/opening,
semiclassical validity/noise, all-sector stability, and the full objective.

Remaining to-do list: complete signed action and range residual; actual G_t
and all application errors; every original coupled/admission budget;
nonlinear existence and compatibility; matching, validity/noise and stability;
transport, preparation and separate creation.
