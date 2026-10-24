# Remove the reference-prefactor defect in the joint source test

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous turn: progress. The full stable-wormhole objective remains active
and incomplete; the feasibility audit stays closed. This is a test of a
sharper physical first approximation. It does not evaluate the complete
quantum source or certify a controlled solution.

## 1. Exact action identity and a changed trial correction

Use the original-G reference I0 and affine canonical fields of
[the cold action defect](current-seed-cold-action-defect.md). Use the
corrected PLUS trace row of [the local assembly](current-finite-local-action-assembly.md).
Let F_l=1-e_R, e_R=384 pi G c_R lambda. The reference with prefactor
F_l/G in [the full Hessian realization](current-seed-full-hessian-realization.md)
has precisely the same covariant seed potential V0. Thus I_ref=F_l I0,
including its consistently repartitioned reference boundary action.
In the unchanged physical action and original renormalization scheme,

    I=F_l I0+J,
    J=-I_pot+I_sq+Gamma_q,
    I_sq=-c_R integral (R+12lambda)^2 dvol,
    I0'(g0)=0, f=I'(g0)=J'(g0), H=I''(g0).

The old actual first correction satisfies I0'' h1=-f. Consequently

    r_old=f+H h1=e_R f+J'' h1,
    h_tilde=h1/F_l,
    r_new=f+H h_tilde=(r_old-e_R f)/F_l=J''h1/F_l.       (1)

This removes an entire known multiple of the seed source, including BOTH
canonical rows and their boundary representatives. Subtracting only a
trace row would not give (1). It is not a second reference iteration,
a change of c_R, a rescaling of a norm or a change of the physical H.
It changes the first approximate metric. The exact linear full-Hessian
solution is unchanged: h_tilde-H^-1 r_new=h1-H^-1 r_old.

The actual h1 remains the evaluated full-source correction, not a local
surrogate. No uncomputed stress is inserted in its place. Fixed leading
clock, period, scalar source, mass/log phase and regular bolt conditions
are linear homogeneous tangent conditions, hence preserved by this
constant multiplication. It supplies no nonlinear domain-neighborhood
theorem. The physical metric boundary identities and common quantum
subtractions are still required on the original realization.

The old bolt allowance is |L_b h1|<3.33e-10. The executed test will bound

    |L_b(h_tilde-h1)| < [e_R/F_l] 3.33e-10 <4e-18,
    |L_b h_tilde| <3.33000004e-10.

The remaining 2e-9 total bolt budget is therefore >1.666999996e-9.
Use the conservative remainder budget 1.66699e-9. The selected working
reduced-source ceiling is 2.60e-40, slightly stricter than the prior
2.6035780e-40 entire-budget necessary ceiling scaled by
1.66699e-9/1.667e-9 for the selected remainder budget. Passing this
screening allocation is insufficient for a physical error certificate.
Missing the deliberately stricter allocation alone does not establish a
violation of that necessary ceiling; compare the actual signed value with
the recalculated necessary ceiling before rejecting that certificate.
The previous half-budget number is not transferred without recalculating
the scalar source return for r_new.

## 2. Execute the known local part of the signed combination

Keep xi=z^4,J_Y=0 and the controlled GLOBAL scalar inverse
u=A^-1 B xi from [the scalar projection](current-schur-scalar-projection.md).
The corresponding affine test is (v,-z^4), v=z^4/2+u. Its scalar row
under the full H vanishes. No scalar solver refinement is selected.
With W=-N pi B^2/(G gamma), the quantity now required is

    P_new=<r_new,(u,xi)>/W
         =[R_nu(v)-R_b(z^4)-e_R(F_nu(v)-F_b(z^4))]/F_l. (2)

F_nu and F_b are the FULL seed gradient rows; R_nu and R_b are the old
defect rows. In particular F_b is S_rest, not a separately chosen pressure.
The complete value (2) is still unknown.

Put x=lambda B^2, t=r/B, z=1/t, d_u=d/x, chi=16 pi G c_R lambda,
tau=T0/lambda, r=4z^5-6z^6, r_n=-12+d_u r and
t_n=(3+1/x+d_u)z^4-d_u(5z^5-9z^6/2). Then O/lambda=r_n/2-2t_n.
All tau coefficients, including the exact product relation and actual
positive-temperature mass term, are those of the preceding scalar test.
Write D_t=-z^2 partial_z, F_n=F/x, and

    K=d_l z^4 partial_z(F_n partial_z tau)+eta K_t tau,
    d_l=(96 pi G c_R-G/(60pi))lambda,
    eta=G lambda/(60pi), K_t=11r_n/36-7c_n/3+11t_n/3.

Equation (5) of the cold action defect gives the explicit T0-only shape
projection, after multiplication by z^4 cancels the r^4 density:

    B_T=chi (x/C) integral [2(O/lambda)tau+2F_n D_t^2 tau
                                  +2(D_t F_n+4zF_n)D_t tau] dz.  (3)

The covariant R2 variation behind (3) retains both Hessian terms; see
[De Felice--Tsujikawa, section 2.1](https://arxiv.org/html/1002.4928v2#S2.SS1).
That reference is used only for the local field equation. It does not
prove our quantum boundary/domain realization or a source bound.

The known part of (2) is therefore

    P_T={ (x/C) integral z^-4 v(K-e_R tau) dz-B_T }/F_l. (4)

Use the already certified fractional trial for v, exact power primitives,
and all actual parameter intervals. Bound the numerical inverse error by

    E_T <=[2.3300059e-44+e_R C_seed(3.6891662e-12)]/F_l,
    C_seed=(x/C) sqrt(integral z^-4 tau_abs^2 dz),

where tau_abs is the absolute-coefficient polynomial. Its outer order is
at least four, so the norm is finite. The previous full trace majorant
also majorizes K, which justifies using its larger pairing-error bound.
Freeze the acceptance target |P_T|+E_T<1e-42. This tests an explicit
signed component on the actual geometry, not the entire source.

## 3. Cancel the unknown principal shape stress algebraically

The old shape row contains -2chi r_n S_rest. Its seed subtraction in
(1) gives EXACTLY

    -2chi r_n S_rest-e_R S_rest=-2chi d_u r S_rest.       (5)

This holds for the complete stress without evaluating it. The known
global bound |G(rho-p)/(lambda r^2)|<2.01e-7 implies
|z^4 S_rest|<(2x/C)2.01e-7. Since
integral_0^1 (4z^5+6z^6)dz=32/21, the remaining projection satisfies

    E_principal <=[4chi |d_u| (x/C)(2.01e-7)(32/21)]/F_l. (6)

Freeze E_principal<1e-43. This closes this specific unknown-stress
coefficient at the required scale. It does not bound the other unknown
quantum, density, radius-transport or f-dependent local terms.

With (4) and (6), a sufficient allocation for the screening projection is

    |P_remaining| <2.589e-40,
    |P_T|+E_T+E_principal <1.1e-42,
    hence |P_new|<2.60e-40.                               (7)

P_remaining is an UNKNOWN SIGNED joint quantity. The inherited scalar
contact allowance, after division by F_l, remains about 5.81e-39, over
22 times the entire screen. It is a known upper allowance that does not
certify (7); it is not an observed error. No corresponding complete
shape value has been evaluated. A factor above 22 reduction in that
one unsigned allowance alone would still not close the other terms.
Correlated action evaluation can outperform those separate allowances,
but no numerical evidence yet establishes that it does.

## 4. Decision, continuum limits and verification scope

Adopt h_tilde as the sharper first approximation for subsequent source
tests, retaining the original h1 and all its evidence. The executed local
screen removes the reference-prefactor loss by an exact physical trial
change. Do not claim that the total source is smaller until P_remaining
is bounded. The finite-seed H, its scalar Schur elimination, domain and
model-comparison norm remain the same. A new nonlinear center still
needs its actual residual, admission and full-ball estimates.

The SINGLE next calculation is the complete residual projection
J''(h1,(v,-z^4))/(F_l W), with the remaining common quantum/contact terms
evaluated together on the actual paired boundary domain, and a certified
error small enough to decide (7). This is proposed, not executed here.
Before calling a finite-mode calculation that value, it must include the
full angular tail, common a0/a2 restoration, domain derivatives and the
one-sided mass/log regularity. Trace anomaly data cannot supply the
independent shape response. Fixed-seed kernel bounds are not full-ball
quantum-source differentiability bounds. Endpoint or compatibility failure
in this evaluation takes priority over further source optimization.

The same full nonlinear requirements remain: rho=1e-8,
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with sufficient
Y0<=5e-9, Z0<=.25, Z2<=5e7. The X/admission inverse, exact complement
and constraints, physical matching, causal energy balance, semiclassical
validity/noise and all-sector stability remain open. No physical solution
or nonexistence result follows from this component.

Freeze this argument, implementation and source review before execution;
retain every run and source snapshot. Recheck the action identities and
the covariant shape variation symbolically, then the actual interval
projection and bolt-budget cost. Independently integrate the retained
fixed-trial polynomial/fractional expressions. Preserve the exploratory
attempt missing the +22 term in Ps and its corrected replay as diagnostics;
neither is acceptance evidence. No Lean is selected.

The first frozen run passed the numerical component targets but failed an
exact decimal-budget identity tested with Arb interval equality. Test that
rational identity with exact rationals instead; this changes no bound or
scientific allocation. Preserve the failed run and replay the full package.

Remaining to-do list: complete the joint quantum/contact source projection;
actual Schur forward/adjoint solve and physical bolt/exterior certificate;
X/admission, nonlinear closure and exact compatibility; matching and causal
realization, semiclassical validity/noise, all-sector stability, charge/return,
opening, transport, preparation, separate creation and the full stable-wormhole objective.
