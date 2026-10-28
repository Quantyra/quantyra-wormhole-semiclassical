# Original-dual complement attempt: no application certificate obtained

2026-10-27 HST, E003/S169. The preceding strong-residual rank test was
progress. This is the first follow-up without a new consequential bound.
The full T_t-(I+T_t)F_n goal remains ACTIVE/INCOMPLETE. The calculation
below is conditional accounting for an attempted block construction, not
an implementation or a completed feasibility milestone.

## Actual spaces and the attempted construction

Keep the paired mass/log E_t form domain, ||u||_Et^2=ell_t[u,u]/N,
N=|W_action| beta, with the original York return, both equations, boundary
phases and regular bolt. Keep the exact UV logarithm, common regulator,
all angular/thermal modes, contacts, cutoff transitions and global returns.
The complete scalar elimination is B=b* a^-1 b. Define the ANALYSIS
Riesz map R_t=ell_t/N and the bounded operator

    A=R_t^-1 (ell_Q+B)/N=I+T_t+Q,
    ||Q||<=e_Q, ||T_t||<=q_t, A>=a I, a=1-e_Q.

This A is the energy-normalized form operator, not the H-associated
unbounded A_Q of the previous test. They cannot be interchanged. The
accepted comparison and source inequalities justify these analysis maps;
R_t^-1 is not thereby a numerical routine.

Take an E_t-orthogonal projection Pi onto a conforming finite space
(the existing quartic direction z^4 is one admissible choice). Write
A in blocks A00,A01,A10,D on Pi E_t plus its full complement. All cross
blocks stay. The coercivity gives D>=a I and ||A01||,||A10||<=b=q_t+e_Q.
This uses no finite replacement of the quantum spectrum.

For an ACTUAL complement approximation C, suppose its full defect were

    d_C=||I-D C||<=d, C:E_c->E_c, E_c=(I-Pi)E_t.       (1)

In weak coordinates the required map is from the complement energy dual
to the complement energy space; its Riesz identification must be explicit.
The exact inverse D^-1 exists, but assigning C=D^-1 simply restates the
missing solve. No such implemented C or bound (1) is available here.

The attempted joint block inverse would be

    H_C=A00-A01 C A10,
    P_C = [ H_C^-1,             -H_C^-1 A01 C;
           -C A10 H_C^-1, C+C A10 H_C^-1 A01 C ].       (2)

If h=a-b^2*d/a>0, the exact Schur complement is >=a I and
||H_C^-1||<=1/h, ||C||<=c=(1+d)/a. Direct multiplication, without
commuting any blocks, gives

    I-A P_C = [0,0;
              -(I-D C) A10 H_C^-1,
               (I-D C)(I+A10 H_C^-1 A01 C)],

so its complete energy-norm defect obeys

    delta<=d sqrt[(b/h)^2+(1+b^2*c/h)^2].               (3)

The quartic block alone does not bound (1). Equation (3) retains the
source/complement couplings; it does not make an unknown defect small.
The coefficient, projection and finite Schur inversions would also need
their numerical errors, not exactness inferred from their definitions.

## First unsupported implication and evidence from the attempt

The first missing constructive implication is a NONCIRCULAR complement
map C on the actual domain with a certified full defect (1). The current
geodesic multiplier code evaluates scalar spectral functions and their
contact-corrected quadrature. Inspection of its functions and of the
global-return implementation finds orbit geometry and seam-tail bounds,
but no inverse for the compressed, York-returned coupled operator.
Its existing flags agree with the implementation bodies: the global
operator/forward-adjoint solve and full-profile regularity are unevaluated.

In particular, if V lifts a radial shape field to the unit tangent bundle,
the known form V* A_flow(P) V is a COMPRESSION. Inverting the scalar
function A_flow does not invert that compression or the sum with the
other local and B terms. Preservation of the physical boundary domain,
the complement output map, and a bound on the resulting joint defect
are still required. No commutation or domain theorem supplies them here.

A tempting shortcut would set C=I in energy coordinates, since
||D-I||<=q_t+e_Q. But the identity inverse in those coordinates represents
the full complement R_t^-1 in weak coordinates. Applying it would require
the unresolved quantum Riesz solve. The small analysis norm is valid;
calling it an implemented preconditioner would be circular.

The previous strong-source graph constant does not supply C: it bounds
A_Q S_Q, not an approximation to S_Q or to D^-1. Neither a fourth-flow
jet nor an actual signed application error follows from it. The prior
rank rejection is also retained; it is not used as an arbitrary-source
counterexample to this proposed infinite-dimensional approximation.

Thus the missing quantity is an UNKNOWN constructive approximation/defect
estimate. The actual seed operator is defined. No excessive actual error,
failed weak complement map, or construction obstruction is proved by
this follow-up. No new consequential evidence toward the complete
numerical bound was obtained.

## Conditional quantitative test, explicitly not an attained constant

Execute (3) at the trial requirement d=.4 using the accepted e_Q and q_t.
This checks whether an attainable modest complement defect WOULD suffice
as a starting point. It is not a measurement of d_C. If delta<1, the
ordered correction

    P_N=P_C sum_(j=0)^(N-1) (I-A P_C)^j,
    S_N=P_N T_t,

would satisfy rho_S<=q_t delta^N with EXACT applications. The factors are
ordered as written; no source or block commutation is assumed. For an
implemented S_hat_N with complete map error epsilon_app,

    rho_S,implemented<=q_t delta^N+(1+q_t+e_Q)epsilon_app. (4)

Select the first integer N making the first term <= half the sufficient
4.989e-32 source-residual screen, then compute the remaining epsilon_app
ceiling in (4). All scalar, quantum, projection, solve, caller, return and
rounding errors belong in epsilon_app. The separate outer K application
and the full F_n/G_t estimates still remain. No conditional ceiling is an
available geometric allowance, nor an original physical-norm comparison.

The accompanying receipt executes only this requirement arithmetic and
checks source identities. It does not constitute progress on the missing
complement estimate. No new checker campaign or component optimization
is selected, and the full goal is not marked complete or blocked here.

Decision: retain the original-dual method conditionally, with the complete
principal operator and all couplings. One decisive next calculation is
to construct a domain-preserving complement preconditioner in weak
coordinates and bound its full defect against (1), including the actual
compressed quantum action. This remains PROPOSED, not executed. The
method cannot be called numerically viable until that implication is proved.

Remaining to-do list: actual complement and S/K applications with errors;
full residual and signed G_t; every original coupled/admission, mass and
phase budget; nonlinear existence/contraction, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
