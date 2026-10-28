# Test a finite-output implementation in the complete strong residual norm

2026-10-27 HST, E003/S169. Ordinary argument frozen before execution;
the initial receipt retains its exact pre-result text and arithmetic.
The previous turn bounded the complete source and associated operator
graph, but did not implement S_Q. Test the proposed implementation class
before assembling a large matrix: any linear map S_hat with output in a
fixed space of dimension at most one million, certified through

    beta_H=||B_H-A_Q S_hat||_(Et->H),
    H=L2(z^-4 dz), A_Q=(ell_Q+B)/N in its H realization.       (1)

The safe sufficient target is beta_H<=4.964e-32. A coefficient or action
rounding error cannot improve a lower bound valid for EVERY such output
map, even with exact coefficients. If its outputs are not in Dom(A_Q),
(1) is unavailable already. An adaptive method with input-dependent
output spaces, or an analytic infinite-rank solver, is outside this test.
This is not a test of the full physical inverse's existence.

## 1. Actual-domain fields, not arbitrary prescribed quantum sources

Keep the actual paired mass/log closure E_t, original period/scale,
complementary Dirac domains, common regulator, full angular/thermal sums,
contacts, cutoff transitions and regular bolt. Work on I=[.4,.6], where
all relevant coefficients are smooth. No operator or boundary is replaced.
Set M=1000001 and divide I into M disjoint cells of width h=.2/M.
On cell j let phi_j(z)=(1-x_j^2)^8 for |x_j|<1 and zero otherwise,
where x_j=2(z-center_j)/h. Let Q be the actual full seed polynomial and

    m=(Q'/z^4)', J=m phi'+2m' phi,
    xi=(QJ)'-z^4 m^2 phi,
    nu=xi/2+QJ/z, T=-xi-(QJ)'/2.                         (2)

Since mJ=(m^2 phi)', the full York moment is zero. Direct differentiation
gives the original York equation. J,xi,nu,T vanish outside their cell,
including every endpoint trace needed here. In particular the original
outer and bolt conditions hold, with zero mass/log coefficients on these
test directions. These compact directions lie INSIDE the actual domain;
they do not replace its nonzero phases or truncate physical matter modes.
The fields are C5 and have more than the required weak derivatives. Smooth
core approximation preserves the strict inequalities below. All cell
transitions are part of the fields and their derivative bounds.

Write xi=sum_(i=0)^2 a_i phi^(i). Exactly,

    a2=Qm, a1=Q'm+3Qm', a0=2(Qm')'-z^4m^2.              (3)

The corresponding nu and T coefficient arrays follow from (2), without
solving a new boundary problem. Qm has a positive minimum on I. Thus
the nonzero disjoint fields xi_j span an M-dimensional real space U_M.

## 2. Uniform bounds on this entire space

Let f_k be the EXACT L2(-1,1) norms of the derivatives of (1-x^2)^8,
k=0,...,5, computed by polynomial integration. Put w=2/h and

    D0=(inf_I a2) w^2 f2 - ||a1||inf w f1 - ||a0||inf f0.

For each cell, ||xi||_L2(dz)>=sqrt(h/2) D0; execute D0>0.
For any of the three coefficient arrays c_i for xi,nu,T and k<=3,

    N_k(c)=sum_i sum_(j=0)^k binom(k,j)
                  ||c_i^(j)||inf w^(i+k-j) f_(i+k-j),
    ||field^(k)||_L2(dz)/||xi||_H <= .6^2 N_k(c)/D0.       (4)

The same inequalities hold for every linear combination in U_M, since
supports are disjoint. For H derivative norms multiply (4) by .4^-2.
For xi itself use its exact H norm, not the loose k=0 bound.
All coefficient bounds are whole-interval enclosures, not samples.

The physical scalar Laplacian and gradient are exactly

    Deltabar xi=(z^2/beta)[Q xi''+(Q'-2Q/z)xi'],
    ||gradbar xi||_H^2=(1/beta) int Q z^-2 |xi'|^2 dz.

Using (4) gives upper constants L2,L1 for these two quantities per
unit H norm. Cell Dirichlet Poincare gives the LOWER constant

    ||gradbar xi||_H^2 >= g_low ||xi||_H^2,
    g_low=(inf_I Q)*.4^4*pi^2/(beta*.6^2*h^2).            (5)

No estimate on an arbitrary-source correction is used. The scalar test
in the original first equation is s=xi. The complete scalar energy obeys

    a(xi,xi)/N <= A_up ||xi||_H^2,
    A_up=6d L2^2+6.01 L1^2+24.01, d=D_R lambda.          (6)

The coefficients follow from the current exact scalar row and its fixed
curvature/coupling guards: ||K||<6.01 and
|2R4+48lambda0|/lambda<24.01. This is a global physical scalar test
with compact support, not a replacement of a^-1 by a scalar matrix.

Let tau_bar=tau/lambda, tau_min=inf_I tau_bar>0, and
t2=||Deltabar tau_bar||inf. Integrate the main mixed row ONCE:

    b_main(xi,xi)/N
      =4 int [C_v tau_bar xi^2+d tau_bar |gradbar xi|^2
                                    -(d/2)(Deltabar tau_bar)xi^2] dmu_H.

The COMPLETE nonleading mixed bound from the previous turn is
K_rest<1.683278e-27 in A x H. It includes all curvature, magnetic,
off-shell and nonlocal York-return terms, including their cross-cell
couplings. Consequently

    b(xi,xi)/N >= b_low ||xi||_H^2,
    b_low=4(C_v tau_min+d tau_min g_low-d t2/2)
                                                    -K_rest sqrt(A_up).

Execute b_low>0. Maximizing 2b(s,xi)/N-a(s,s)/N over the SINGLE scalar
direction s=c xi then proves the full coupled variational lower bound

    (B_H xi,xi)_H >= b_low^2/A_up ||xi||_H^2=:B_low||xi||_H^2. (7)

The exact full scalar inverse remains in B_H. Restricting its variational
maximization supplies a lower bound, not an approximate scalar solve.

## 3. Upper-bound the COMPLETE quantum energy before using dimension

A local classical energy alone cannot upper-bound E_t. Use the accepted
complete geodesic comparison ell_0 with |ell-ell_0|<=.0022 ell.
Its physical quantum logarithm, full common local restoration and all
state/thermal/boundary differences stay in the comparison and its error.

Evaluate ALL FOUR exact affine local densities (2)-(3) of
current-finite-local-action-assembly.md: Einstein/cosmological, magnetic,
additional R2 and Weyl. Use the six ORDINARY independent jets
(nu,nu',nu'',T,T',T''). This retains the off-shell density terms directly;
no unquantified chart-conversion reserve or stationarity assumption enters.
For each allocation alpha, interval differentiation of its displayed
density gives its symmetric Hessian matrix H_alpha(z). If C_alpha,ij
bounds |H_alpha,ij| on I, (4) gives

    E_local <= (1/beta) sum_alpha sum_ij
                            C_alpha,ij Ljet_i Ljet_j.    (8)

Here the off-diagonal terms are counted twice, exactly as in the quadratic
form v^t H v. The interval operation computes the second derivative, not
the Taylor coefficient divided by two. The negative physical action
factor is harmless for this absolute upper bound; beta is the ORIGINAL
conversion from I/W to ell/N. This deliberately conservative local bound
may discard beneficial cancellations only after retaining the full
two-row returned fields. Its looseness weakens, never strengthens, the
residual lower bound (11).

Add canonical a2, LLL and the shared bounded vacuum coefficient exactly
once, using (5)-(6) of current-canonical-heat-restoration-assembly.md.
With R=-Q'', Delta T=QT''+Q'T', their absolute quadratic bounds are
|beta2|[4Rmax^2 LT0^2+16Rmax LT0 LDelta+8LDelta^2],
|Gq/(6pi B_radius^2)| Qmax LT1^2, and U0 LT0^2,
respectively, all divided by beta. Set LDelta=Qmax LT2+Q1 LT1.
The vacuum value remains one bounded shared coefficient; it is not zero.

The fixed log z and F_eff retain their actual values; no scale is reset.
The original boundary Hessians vanish on (2), rather than being omitted
from the operator. This assembly is only an upper energy bound for the
test space, not a signed evaluation of the full residual.

The remaining positive geodesic form has the already proved complete
all-angular bound P_M[T]<=||X^3 T||^2/(525pi^2 M^4),
sum d_n/M_n^4<4C^2/q. The actual full flow retains trapped trajectories
and regular-bolt returns; the generator argument applies because these
compact C5 fields lie in Dom(X^3). With Qmax=||Q||inf, Q1,Q2 its
first two derivative bounds, and the T-jet constants from (4), set

    A=Qmax LT3, D=1.5 Q1 LT2+.5 Q2 LT1,
    E_geo=2(4C^2/q)/(525pi B_radius^2 beta)
                           *Qmax*(5A^2/8+1.5 A D+D^2).   (9)

This follows by angular integration of the exact third-flow derivative;
it bounds the ENTIRE positive logarithmic form, not selected frequencies.
Combining all terms and the actual retained-transport comparison gives

    ||xi||_Et^2 <= E_up ||xi||_H^2,
    E_up=(1+epsilon)/(1-.0022)
                       *(E_local+E_a2+E_LLL+E_vacuum+E_geo). (10)

Disjoint-support derivative norms bound every cross-cell quantum pairing
through (9), so (10) is uniform on U_M, not only its individual basis fields.

## 4. The actual finite-output residual lower bound

For a fixed output space of dimension at most M-1, choose nonzero
u in U_M orthogonal in H to its A_Q image. Then

    (u,(B_H-A_Q S_hat)u)_H=(u,B_H u)_H>=B_low||u||_H^2.

Cauchy--Schwarz and (10) prove

    beta_H >= B_low/sqrt(E_up).                          (11)

This uses the actual structured source B_H and both original equations.
It is not an arbitrary-source counterexample and does not assert failure
of the actual coupled inverse. The analogous full energy-dual test yields
only rho_S>=B_low/E_up. Compute BOTH bounds, their ratios to the actual
sufficient screens, and the resulting floor of the proposed strong-norm
upper estimator's complete source debit. Do not call that estimator floor
a lower bound on the physical output or on the full weak residual.

If (11) exceeds 4.964e-32, this rejects the SPECIFIED strong-residual
finite-output implementation at this rank, even with perfect arithmetic.
A weak lower bound below its screen is NOT a weak certificate or proof
of attainability. No degree, transport-depth or tail optimization is
authorized by that inconclusive comparison.

## 5. Execution and decision boundary

Execute the complete inequalities once with exact polynomial norm
integrals and 512-bit outward intervals, preserving every outcome.
The supporting arithmetic verifies this ordinary full-domain argument;
it neither implements a million-column solver nor evaluates the signed
quantum correction. The next method choice must use the resulting bounds.

The source regularity bound, P32, all accepted results and closed audit
remain unchanged. No mass-column, original geometric admission, nonlinear
existence/contraction, compatibility, matching, semiclassical validity/noise
or all-sector stability requirement is closed by a finite-output estimator
failure. An infinite-rank structured inverse and the original energy-dual
certificate remain distinct alternatives with their own unproved errors.

## 6. Executed result and method decision

The initial 512-bit execution is retained in
[the initial receipt](../evidence/current-strong-residual-rank-test.json).
Its whole-interval coefficient bounds and exact bump integrals give

    A_up < 6.063886e24,
    b_low > 8267.08677664,
    B_low > 1.1270781644e-17,
    E_up < 1.200446e24,
    beta_H > 1.0286858190e-29.

The last bound is over 207.22 times the sufficient 4.964e-32 screen for
EVERY linear map into a fixed at-most-million-dimensional output space
in Dom(A_Q). In particular this failure persists with exact coefficients
and exact quantum applications. It cannot be repaired by reducing their
rounding errors. If a proposed output is outside Dom(A_Q), the stronger
residual is undefined rather than small.

Through the existing strong-residual conversion and complete 12-term
accounting, this estimator's source debit alone is >1.3893307489e-29:
over 13.89 times the ENTIRE sufficient 1e-30 residual screen. Thus giving
it all of the conditional residual room, instead of its prior source
allocation 6.704319e-32, still cannot certify that screen. This is a lower
bound on the report required by this upper estimator; it is NOT a lower
bound on the actual weak application error, signed G_t or physical output.

The mathematical loss is explicit: the complete B_H source has a positive
response on every direction of U_M, whereas a fixed smaller output image
leaves an H-orthogonal direction uncancelled. No arbitrary external source
or finite selection of matter modes was substituted. The actual Q and
York identities, scalar maximization, all local Hessians, common
restoration, full spectral geodesic bound and two-sided complete quantum
comparison were used. The original boundary and regular-bolt conditions
hold on the test directions, and the operator retains the full domain.

The first required continuum implication for this test is membership of
these returned fields in the actual form domain. Compact C5 fields from
(2) have the required scalar H2 and flow H3 control. Smooth approximation
of phi preserves the exact York construction, zero outer moment and bolt
traces; (8)-(10) control the complete energy convergence. This proves
membership without a smooth surrogate domain. It does not assert that a
particular numerical output has the stronger Dom(A_Q) regularity.

The specified alternative is to certify the original E_t* residual
rho_S of the complete reference solve, while retaining an analytic
infinite-dimensional principal action in the approximation. The same
test yields only rho_S >= 9.3888298768e-42, about 1.882e-10 of its
4.989e-32 sufficient screen. This small LOWER bound neither certifies the
alternative nor proves it attainable. The already proved coercivity gives
a_S <= rho_S/(1-e_Q); meeting that screen would spend <6.705e-32 in the
complete residual. This is a conditional allowance, not an attained error.
Returning to the original dual norm is not a normalization improvement:
the original conversion must still be used, and the error must actually
be bounded in that norm.

Decision: reject a fixed at-most-million-dimensional strong-residual
certificate as the next implementation method. Retain the full geodesic
principal operator and original energy-dual residual. Merely enlarging
a finite matrix beyond one million dimensions has not been certified or
recommended. The dominant missing bound remains the complete reference
application and its infinite-dimensional complement, not the controlled
P32 transport error or quartic column.

One decisive next calculation, still PROPOSED: construct and bound the
complete weak residual B-(ell_Q+B) S_hat of an approximation that retains
the analytic principal action on its complementary space, including
the scalar solve, common quantum action, returns and arithmetic errors,
against rho_S<=4.989e-32. No such implemented map or upper certificate was
obtained here. K still needs its separate application certificate before
the full T_t-(I+T_t)F_n claim can close.

The final receipt freezes this outcome and the added comparison with the
ENTIRE residual screen. Neither execution supplies a new physical margin.

Remaining to-do list: actual full-domain implementation and numerical
quantum/source/solve/complement errors; complete residual and signed G_t;
all original coupled/admission and physical/stability requirements.
