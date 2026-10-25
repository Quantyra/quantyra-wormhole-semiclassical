# A regular reference trial for the current mass adjoint

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: PROGRESS (current coupled mass-adjoint certificate).
The full stable-wormhole objective is incomplete; the audit remains closed.

OUTCOME: the actual coupled estimate remains UNRESOLVED. An explicit
regular-bolt adjoint trial is constructed, with BOTH reference rows zero.
Its complete seed-source pairing has an available absolute enclosure, but
that specified enclosure is much too wide for the mass certificate. The
actual signed pairing and the full quantum adjoint residual are unknown.
This is failure of the specified unsigned estimator, not a mass lower
bound, instability, or obstruction to the physical construction.

## 1. Actual reference trial, with its clock and temperature distinguished

Keep the current finite Q, z_h=1/(1+h), Q(z_h)=0, Q'(z_h)=-2gamma,
beta=lambda B^2/C and I_ref=F_l I0. Define

    I(z)=integral_0^z dt/Q(t), J(z)=integral_0^z t^3/Q(t)^2 dt,
    n_c=-Q I/z, b_c=1+Q'I/2,
    n_m=-2beta QJ/z, b_m=beta z^3/Q+beta Q'J,
    D_A=z_h^2 Q''(z_h)/2+3gamma z_h,
    m_T=gamma^2/(beta z_h D_A),
    Psi_ref=(n_c,b_c)+m_T(n_m,b_m).                    (1)

These are the actual rational-seed reference profiles, not a product
geometry or smoothed seam. The exact areal return is

    f=2F(-D n-b)+(D F)n, alpha=n+2b+D n,
    F=CQ/z^2, D=z partial_z.

It gives f_c=0, alpha_c=1, f_m=2(lambda B^2)z, alpha_m=0.
Thus (r f)'=0 and alpha'=0 for Psi_ref: both reference constraint rows
vanish, and the original joint Ward identity supplies their angular row.

The unchanged reference temperature derivative per unit mass is
-D_A beta z_h/gamma^2. With alpha_infinity=1, (1) gives

    delta kappa_coordinate/kappa=1-m_T D_A beta z_h/gamma^2=0.

The coordinate period and regular bolt therefore match. This is an
AUXILIARY open-clock adjoint direction. It is not a physical fixed-clock
zero mode. It introduces no clock or mass degree of freedom in the solve.
In particular the original cut mass lift, its .1/.2 transition source,
and the .05/.1 original graph cutoff remain unchanged. An uncut mass
profile is used here only in the correlated auxiliary trial (1).

## 2. Exact cancellation at the bolt before numerical evaluation

Let delta=z_h-z, q2=Q''(z_h), U(z)=Q(z)/(z_h-z), and

    c=z_h^2/(2D_A), A=I+2beta m_T J.

Separately I and J have logarithmic terms at the bolt. Their combined
logarithmic coefficient in A is exactly

    -1/(2gamma)+2beta m_T[
             3z_h^2/(4gamma^2)+q2 z_h^3/(8gamma^3)]=0.     (2)

There is a stable analytic representation on the whole actual interval:

    V(z)=[(z_h-z)U(z)+2beta m_T z^3-c U(z)^2]/(z_h-z)^2,
    A_reg(z)=-c/z_h+integral_0^z V(t)/U(t)^2 dt,
    A(z)=c/(z_h-z)+A_reg(z).                             (3)

The numerator defining V has a DOUBLE zero at z_h by (1), so V is a
polynomial after exact division. U is a polynomial, positive on the
closed real interval, with U(z_h)=2gamma>0. No separately rounded Q root
is divided by at the bolt. Equation (3) is an exact antiderivative
identity, not an endpoint extrapolation from sampled values.

The trial components can be evaluated without subtracting poles:

    n_ref=-U[c+(z_h-z)A_reg]/z,
    W(z)=[c Q' U/2+beta m_T z^3]/(z_h-z),
    b_ref=1+Q'A_reg/2+W/U.                              (4)

The numerator of W has a simple zero at z_h. All functions in (4) are
analytic there. In particular n_ref(z_h)=-gamma z_h/D_A, agreeing with
s_R/R for s_R=-kappa/D_h, D_h A_radius^2=D_A. Analyticity in z_h-z is
regular rotational-disk smoothness, since z_h-z is a smooth function of
rho^2 with positive leading coefficient at this fixed nonzero temperature.

At zero, write Q=beta+u z^2-a z^3+..., a=4beta+2u+tau. Then

    n_ref=-1-2u z^2/(3beta)+[3a/(4beta)-m_T/2]z^3+O(z^4),
    b_ref= 1+u z^2/beta+[-3a/(2beta)+m_T]z^3+O(z^4).       (5)

Thus with P=(-1,1), v_ref=P-Psi_ref is in the ACTUAL critical energy
domain: it starts at weight two and is regular at the original bolt.
The fixed-seed graph/form argument of current-coupled-mass-adjoint.md
applies. It is an auxiliary adjoint correction, not a proposed nonlinear
metric in the selected mass/fast phase.

The tiny areal mass parameter m_T must not be called the trial's canonical
cubic coefficient: that coefficient is m_T-3a/(2beta), near -6. Neither
quantity is the unknown physical solution's mass. The calculation checks
this explicitly to avoid a misleading small-mass claim.

## 3. Apply the full operator: what vanishes and what remains

Use exactly the physical energy normalization h,F,q_P of the current
coupled mass-adjoint certificate. With h=h_ref+j, j=J_rest''/(|W|beta),
and v_ref from (5), its exact residual is

    r_ref=h v_ref-q_P=-h Psi_ref=-j Psi_ref.             (6)

The reference contribution is zero GLOBALLY, including regular-bolt
matching, not merely on a collar. The remaining j contains -I_pot,
I_sq, the entire renormalized Dirac action, density/York contacts and
every angular/thermal and finite-domain return. None is assigned zero.

Equation (6) is well-defined in E*: q_P is in E* by the current auxiliary
domain argument, v_ref is in E, and h is a bounded form on E. The first
unsupported quantitative step is a norm bound for this COMPLETE r_ref.
The existing strong interior estimate concerns a different restricted
output norm; the fast screen does not cover this independent root-zero
input. The source's radial derivative bound Theta'(t) is not a metric
or temperature derivative and does not bound (6). Likewise the old
interior first-Cauchy estimate requires its direction norms, fixed-area
density term and the boundary part before it could bound this residual.
No unknown constant, old-theory constant or family analyticity is assumed.

For zero forward trial the signed center is F(Psi_ref). It can be
identified with an ACTUAL complete source moment without evaluating a
quantum Hessian. Let h1 be the existing original-G first correction,
I0''h1=-I'(g0). The current reference response is h1/F_l. Its reference
mass Green constant is 2F_l in the SAME energy units, so

    F(Psi_ref)=2 m(h1)
      =-2[ kappa R s_R+G integral_R^infinity rho_rest dr ]/(lambda B^3).
                                                               (7)

Here rho_rest is the full evaluated seed source. This uses both original
reference equations and fixed-clock horizon matching for h1. Its known
|s_R|/(Bh)<3.33e-10 is retained. The rescaling F_l cancels exactly in
(7); it improves neither the physical mass target nor the source error.
F is the compact-test source extended by the convergent moment of the
preceding adjoint argument. Replacing (7) by an unrestricted open-clock
derivative of the whole action would require its additional boundary
one-form. It must not be omitted or counted again in this bulk moment.

The actual coupled mass certificate consequently reads

    |m_actual| <= [|F(Psi_ref)|+
                    (1.1e-10)||j Psi_ref||E*]/K_m,
    K_m=2Z_c>1.9999999764.                            (8)

Numerical quadrature error for (7) must also be added when it is evaluated.
This is a specific full coupled trial, with its remaining error explicit.

## 4. Execute the available COMPLETE source-moment enclosure

Do not integrate the massive and Weyl tails separately. The current
rational-family guards and fixed-scale restoration verified in
fixed-clock-seed-bolt-feasibility.md and current-global-return-energy-attempt.md
give the same complete Dirac, vacuum-subtracted estimate proved by the
explicit combination in stationary-centered-exterior-source.md:

    B integral |rho_D-rho_vac| dr <8.001e15 q_flux^2.   (9)

Its required K=1e15, K0=1e11, outer remainder 4e15, |Fang|<100,
M1>18, rational coefficient bounds and actual paired/state realization
hold at the current seed. The leading L/t terms cancel BEFORE absolute
integration. The bound includes all angular and thermal modes, the actual
finite-domain gluing, LLL and Weyl terms. The old reference-subtraction
mass allowance is not transferred to the new c_R theory.

Instead add the current local gradients directly. With delta I=int E:delta g,
the stress has magnitude at most 2|E|. The actual bounds from
current-global-return-energy-attempt.md, section 2, are

    (G/lambda)|E_pot| < e_p t^-4, e_p=5.5e-23,
    (G/lambda)|E_sq| <= e_s t^-5,
    e_s=c_R g[200w0+w0^2+10w2],
    w0=1.1e-27, w2=4.3e-25, c_R=1e22, g=G lambda.

Sphere integration and integration in r give the dimensionless complete
absolute source-mass bound

    M_abs=8.001e15 L_q+8pi e_p+4pi e_s,
    L_q=G q_flux^2/(lambda B^4).                       (10)

The factors 8pi and 4pi include the stress/action factor two and the
integrals int_1^infinity t^-2=1, int_1^infinity t^-3=1/2. The actual
lower limit v=1+h only decreases these upper bounds. The same vacuum
and original finite coefficients are retained throughout.

The unchanged bolt allowance costs at most

    M_bolt=(gamma/beta)(1+h)h(3.33e-10),
    |F(Psi_ref)| < I_abs=2(M_abs+M_bolt).              (11)

Evaluate (10)-(11) with the existing Arb enclosure of the current seed.
The execution contract is I_abs<7e-11 and I_abs/(3.3e-15)>1e4.
The actual tighter reported interval is retained. This certifies an
enclosure about zero using the already available absolute source bound;
it does not locate the signed value inside that interval.

The SPECIFIED unsigned estimator (11) therefore cannot certify the
selected 3.3e-15 center allowance, even if the adjoint residual were zero.
This is a known excessive upper bound. The actual signed center is UNKNOWN,
and the residual norm in (6) is a separate UNKNOWN estimate. Neither is
proved large. As a quantified alternative, an absolute bound replacing
8.001e15 in (9) by

    C_mass <= [3.3e-15-2(M_bolt+8pi e_p+4pi e_s)]/(2L_q)       (12)

would fit the selected center budget. A direct SIGNED evaluation of (7)
can do so without improving the entire absolute-source norm. This
calculation selects that signed evaluation, not generic source optimization.
The already controlled cold value contributes negligibly to (11); no cold
or scalar refinement is proposed. There is still no estimate for (6).

## 5. Decision and execution boundary

EXECUTED: the full reference adjoint trial, exact cancellation of the
bolt logarithm, regular auxiliary-domain realization, complete source-
moment identity, original-unit absolute envelope and its required
improvement. The checks verify those identities and arithmetic only.
NOT EXECUTED: a narrow signed source-moment enclosure or a certified
full quantum adjoint residual. Thus (8) has not certified actual mass,
and the requested full coupled feasibility milestone remains incomplete.

The original full conversion 2133<||d_chi||H3<2174, conditional mass
target 2.29e-15, kappa_m<=5e-5 target, residual global budgets and the
.135973<.1599164 two-fast-input screen retain their exact prior scopes.
No original graph, physical-norm embedding, nonlinear neighborhood,
compatibility, matching, validity/noise or all-sector stability gate closes.

In particular kappa_m is the post-elimination column of the FULL coupled
return, not m_T or the scalar mass observable in (8). If that column were
<=5e-5, its affine cost would be 2.14765e-12 and its direct graph feedback
.0024731. After the controlled far contribution, the remaining entries
would still have to satisfy

    Y_other <=5.34234999994958e-12,
    eta_other+(.0195239+4.583e-13) gamma_E <.0012069.    (13)

Here gamma_E is the required source-family energy/graph comparison, not
Euler's constant. Even with eta_other=0 it must be <.06181654. The
quartic, correlated slow pair and fast/interior couplings in (13) have
not been bounded by this mass trial. No remaining entry is set to zero.
The scalar original row 20beta Fstar c4 and the exact slow difference
z^4-z^(4+e_star), e_star about 1.865097e-34, remain in the full system.

The first unsupported FULL GRAPH step remains numerical control of the
differentiated common-regulated core composed with the coupled inverse,
including its finite-domain and transition returns. The fixed-seed
weak/graph realization is qualitative. Its energy estimate is not an H3
estimate. The physical input/output comparisons and uniform nonlinear
neighborhood estimate are also unknown. In particular a two-fast-slot
bound does not bound independent mass/quartic/slow phase variations.
The actual-original-graph criterion and preserved-ball/contraction
requirements in the closed audit are unchanged.

Equation (12) quantifies the required improvement of this one estimator;
it does not make that improvement plausible by itself. There is no
retained signed measurement or complete residual estimate from which to
infer attainability. Classical cancellation, cold smallness and finite
energy alone do not answer it. Conversely, the positive envelope gives
no lower bound on the source moment, so it supplies no reason to reject
the construction. This is the precise feasibility decision still open.

RESEARCH DECISION: retain the construction and this regular coupled
adjoint trial; sharpen the SIGNED complete mass integral rather than
substitute the available unsigned source envelope. Its excess is an
estimation failure, not evidence against the physical mass. The full
residual in (6) remains mandatory before accepting the mass certificate.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): enclose the full common-
restored quantum source moment in (7), using this explicit trial and
retaining the massive/Weyl cancellation inside the integral; in the same
operator assembly bound (6) against 1e-5. Test their combined error in
(8) against 2.29e-15. Actual original graph control remains independent.

Remaining to-do list: executed signed mass pairing and coupled adjoint
residual; actual original graph return; physical-norm/neighborhood
comparisons; nonlinear admission/contraction; exact compatibility;
physical matching; semiclassical validity/noise; all-sector stability;
charge/return, opening, transport, preparation, separate creation and
the full stable semiclassical wormhole objective.
