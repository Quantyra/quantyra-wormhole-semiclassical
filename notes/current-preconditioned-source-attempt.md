# Complete preconditioned source attempt: global application still unresolved

2026-10-27 HST, E003/S169. The previous goal turn was PROGRESS: its
full-operator lower bound rejected a one-dimensional correction range.
This follow-up does NOT obtain a new complete application or residual bound.
It audits the first attempted implementation of S and rejects unnecessary
spectral-quadrature refinement. The active goal, closed audit, original
physical norms, P32 and all accepted results remain intact.

## 1. Actual composition and its domain

Use the same E_t form closure, norm squared ell_t/N, N=|W_action| beta,
paired first-order realization, complementary adjoint and both squares.
Retain the outer mass/log phase, original cutoff transitions, fixed
clock/temperature/spin/scale and regular disk bolt. York still obeys
J(0)=0 and (Q^2 J'/z^4)(z_h)=0. Both original equations give
B=b* a^-1 b, with the complete scalar and mixed forms.

Distinguish the Riesz maps from the weak forms. If R_t:E_t->E_t* is
ell_t/N, write

    B_N=B/N, L0_form=(ell_0+B)/N,
    T_t=R_t^-1 B_N, L0_Riesz=R_t^-1 L0_form.

Then exactly

    S=L0_Riesz^-1 T_t=L0_form^-1 B_N=(ell_0+B)^-1 B.       (1)

Thus implementing S does not require a separate numerical application of
the physical ell_t inverse. This is an operator identity, not a computed
solution. B_N:E_t->E_t* is bounded by the accepted mixed estimate;
L0_form is coercive by the complete geodesic comparison. Lax--Milgram
defines the inverse on this actual form space. The common regulator,
complete angular/thermal state and coincidence/boundary contacts remain
in the physical forms and their comparison. No new domain theorem is
inferred from the formal identity.

For any proposed numerical S_hat the COMPLETE source residual is

    r_S(u)[v]=B[u,v]/N-(ell_0+B)[S_hat u,v]/N,
    rho_S=sup_(||u||Et=||v||Et=1) |r_S(u)[v]|.

If ||ell_t^-1(ell_0-ell_t)||<=e0<1, then

    ||S_hat-S||_(Et->Et)<=rho_S/(1-e0).                   (2)

The supremum in rho_S includes all inputs and tests, including the entire
complement of a numerical trial space. No finite Ritz matrix supplies it.
The first unsupported implication in this attempt is a computable bound
on rho_S for an actual implementation. The underlying form operator is
defined; its numerical application constant is unknown.

## 2. Test the apparent multiplier-accuracy bottleneck

The existing current-geodesic-flow-multiplier.md constructs a positive
contact-corrected quadrature of the complete flow multiplier A(P), where
P=iX generates the actual C2,1 doubled geodesic flow and the base lift is
retained. It is not the flat angular-averaged multiplier K(P). The weights
and common node shift use exact formulas, preserving both the zero-frequency
fourth contact and the ultraviolet logarithm's finite constant. The full
angular tower is resummed before this quadrature. Actual thermal/state
and global return differences remain in the comparison remainder.

The retained receipt gives

    eta_quad <2.535721e-5,
    eta_actual_to_exact_reference <=.0021040238,
    eta_actual_to_quadrature_reference <.002129381007<.0022.

This permits defining ell_Q as the quadrature comparison form and using
it CONSISTENTLY in both parts of the algorithm:

    S_Q=(ell_Q+B)^-1 B,
    L_Q=ell_t^-1(ell_Q+B),
    K_Q=L_Q^-1(I+T_t-Q_delta-L_Q),
    F_delta,12,Q=sum_(k=0)^11 (-K_Q)^k S_Q.                (3)

The exact inverse identity and ordered Neumann remainder hold for any
such comparison. The target ell_t, source B, physical output and admission
allowance do not change. This is a consistent change of computational
reference, not an improved physical constant by normalization.

If instead one insists on approximating the FIXED exact-reference S by
S_Q while leaving the rest of the algorithm fixed, the generic bound
(with e_quad=eta_quad/(1-epsilon) in E_t)

    ||S_Q-S|| <= e_quad f_t/[(1-e0)(1-e0-e_quad)]

reports about 6.434074e-6, over 1.28681e26 times the proposed 5e-32
application target. This is an excessive upper estimator, not a measured
error or an obstruction. Refining the multiplier to remove that number
is unnecessary when (3) is used. The same already accepted .0022 guard
still proves the <8.953511e-31 EXACT-APPLICATION residual bound. Using
the tighter recorded comparison would only change it to <8.939032e-31;
this does not resolve the numerical application problem.

Accordingly no higher-order quadrature campaign is launched. Arithmetic
rounding of weights, spectral evaluation and actual geodesic propagation
still belong in the numerical errors of S_hat_Q and K_hat_Q. Exact
quadrature formulas are not declared to be a full numerical inverse.

## 3. Attempted global realization and the precise missing estimate

Inspection of current_global_geodesic_return.py and its ordinary argument
shows orbit geometry, trapped trajectories and a cubic-seam Fourier-tail
formula, but no implemented global inverse or complete source residual.
Its result explicitly leaves full_profile_J4_certified and
complete_global_operator_evaluated false. The multiplier implementation
likewise leaves geodesic_flow_matrix_evaluated false. These statements
agree with the inspected function bodies; they are not just status flags.

The orbit-tail estimator requires a quantitative bound on the regular
fourth flow derivative after subtracting the explicit mass-cubic jump,
together with the mass coefficient, for the actual correction range.
No uniform bound for V S_Q applied to the whole unit E_t ball is supplied.
The earlier failure of an unmodified periodic-H4 test on a mass lift does
not prove failure on this structured range. Conversely, qualitative graph
regularity does not supply the missing constant or the uniform dual source
residual. This auxiliary fourth-flow-derivative estimate is not silently
made into an original geometric admission requirement.

Small trapped-orbit measure cannot bound arbitrary localized unit inputs.
A finite set of periods or matrix entries also does not bound the omitted
test space. The common contacts and exact UV logarithm prevent replacing
the global flow by a smooth periodic surrogate or a finite matter spectrum.
These are precisely the hypotheses a complete rho_S calculation must cover.
No applicable quantified approximation theorem closing that step has been
established in this attempt. The first unsupported step is an UNKNOWN
uniform approximation/residual estimate, not an undefined seed operator.

The executed numerical contract for (2), with the actual quadrature
reference and retained transport comparison, is

    e_Q <.002129381007,
    rho_S <=(1-e_Q)*5e-32,
    (1-e_Q)*5e-32 =4.98935309496752...e-32.

This is a sufficient required residual accuracy, not an attained rho_S.
No actual S_hat_Q, K_hat_Q, coefficient error or full numerical residual
was obtained. There is NO new consequential bound toward the requested
complete certificate in this follow-up. The arithmetic prevents an
unnecessary component optimization; it is not a completed research milestone.

## 4. Decision and one next calculation

Retain the complete geodesic preconditioner, using its existing quadrature
consistently in S and K; do not refine that already affordable component.
The specific next calculation remains a conforming implementation of S_Q
with a full input-to-dual residual bound rho_S, including the complete
geodesic return and a proved complementary-space error. Test
rho_S<=4.989e-32, which safely implies a_S<=5e-32. A numerical matrix
without its complement is insufficient. This calculation is PROPOSED,
not executed. The separate K application and signed G_t remain open.

This is the first no-new-bound continuation after the rank-one progress
turn. The goal stays ACTIVE/INCOMPLETE; neither completion nor an actual
construction obstruction is established. No physical budget changes.

Remaining to-do list: complete source and quantum applications and full
residual; signed G_t and every original coupled/admission, mass and boundary
phase entry; nonlinear existence, compatibility, physical matching,
semiclassical validity/noise and all-sector stability; the full objective.
