# Physical weak source conversion on the finite seed

2026-10-23, E003/S169. Ordinary argument before frozen execution.
The full stable semiclassical wormhole objective remains incomplete.

The canonical source rows and quantum action cannot be inserted directly
as scalar/TT dual rows. This supplies their exact fixed-seed transpose,
including the bolt cancellation, outer mass phase and action measure.
It also bounds the kinematic return of a source supported in z<=1/2.
That return is NOT the response of the gravitational Hessian inverse.
The full curved quantum remainder and physical correction remain open.

## 1. Freeze the coordinates and action normalization

Keep the current c_R=1e22 action, q=2000^4, original mu0, paired state,
fixed clock and finite-seed parameters. Set z=B/r, z_h=1/(1+h),
gamma=kappa B/C. The canonical disk has

    Q=(1-z)^2[beta(1+2z+3z^2)+u z^2+epsilon z^4]-tau z^3,
    F=C Q/z^2, d rho=-dz/sqrt(Q), a=sqrt(Q)/gamma,
    Q(z_h)=0, Q'(z_h)=-2gamma, Q'(0)=0.                (1)

The last two endpoint equations are exact identities with the original
temperature root. Correlated endpoint expressions must use them before
interval estimation; cancellation of independent polynomial intervals
does not establish a positive bolt denominator.

For N identical exteriors, the canonical area measure is 2pi N dz/gamma.
The prior normalized affine rows (f,g) multiply -B^2/(2G). Thus the
physical action pairing with canonical directions (nu,T) is

    W integral_0^z_h (f nu+g T)dz,
    W=-N pi B^2/(G gamma).                              (2)

There is no extra sphere degeneracy or thermal factor. The same W
multiplies the normalized Hessian and forward defect. It cancels from
the forward equation, but NOT from the physical adjoint right-hand side.
This is a unit conversion, not an improved error estimate.

Use the frozen linear York return of
[the seed form](stationary-seed-quantum-form.md), not a field-dependent
projection differentiated along a nonlinear path. A TT tangent is
k=diag(-b-2xi,b,xi,xi), with

    (r^2 F b)'_r=-(r^2 F'_r-2rF)xi,
    I(r)=integral_r^infinity (b+xi)/F dt,
    S=xi/2+F I/r, T=-(b+3xi)/2+(F'_r/2-F/r)I.          (3)

Here nu=s+S for the physical scalar s. These formulas follow from the
outer-fixed radial gauge return and keep the induced r^-4 radial tail.
The finite-seed full-gradient/chart contacts used by the coercivity proof
remain in the canonical Hessian. If another nonlinear tensor chart or a
moving York map is used, its gradient-times-second-variation contact must
be added explicitly. The seed is not stationary for the full action.

## 2. One positive differential auxiliary problem

Put J=C I/B and P=Q^2/z^4. Equation (3) is equivalent to

    L J=Q/z^4 (4xi/z-xi'), L=-(P partial_z)',
    J(0)=0, (P J')(z_h)=0,
    b=Q J'-xi,
    S=xi/2+Q J/z, T=-xi-(QJ)'/2.                       (4)

The two auxiliary conditions follow from the original outer-fixed return
and regular bolt; they add no physical mass or wall equation. The
explicit inverse for a source q is

    (L^-1 q)(z)=integral_0^z t^4/Q(t)^2
                                  integral_t^z_h q(v)dv dt. (5)

This auxiliary L is a kinematic operator, not either gravitational
diagonal. For smooth one-sided directions xi=O(z^3) and regular at the
bolt, (3) or (5) gives J=O(z^4) at zero and finite J(z_h). In particular
nu(z_h)=s(z_h)+xi(z_h)/2, while
T(z_h)=-xi(z_h)+gamma J(z_h). The nonzero J endpoint term is retained.
Weight-four logarithmic remainders can be included one-sided. Nothing
here upgrades the reflected mass cubic from C2,1 to smoothness.

## 3. Exact weak transpose and cancellation before division

For ordinary canonical source densities f,g with respect to dz, define

    q=Q(f/z+g'/2), w=L^-1 q,
    J_Y^*(f,g)=(f, f/2-g+(Qw)'/z^4).                  (6)

Indeed integrate -g(QJ)'/2 once and use the symmetric Green pairing
for L. Integrating Qw(4xi/z-xi')/z^4 once more gives (Qw)'xi/z^4.
All resulting endpoint terms vanish on the declared source/tangent phase.
For combined sources f,g=O(log^k z), g'=O(z^-1 log^k z), one has
w=O(z^5 log^(k+1) z) at zero. These bounds suffice for those endpoint
terms; no numerical norm bound on the full actual source is assumed.

A derivative-free implementation is preferable when only an oscillation
bound is available. Write

    R_g(z)=integral_z^z_h [Q(v)f(v)/v
                             -(Q'(v)/2)(g(v)-g(z))]dv,
    w'(z)=z^4 R_g(z)/Q(z)^2, w(0)=0,
    [J_Y^*(f,g)]_xi=f/2-g+R_g/Q+Q'w/z^4.              (7)

The equivalent expression
R_g=int_z^z_h [Qf/v-Q'g/2]dv-Q(z)g(z)/2 has two separately
first-order terms that cancel at the bolt. Dividing their independent
unsigned bounds by Q^2 loses this cancellation. Formula (7) keeps it.
For g constant and f=0, R_g=w=0 EXACTLY and the returned shape row
is -g. This is also the integrated identity int T=-int xi.

At the bolt let delta=z_h-z. Since Q=2gamma delta+O(delta^2),
bounded f and a modulus |g(v)-g(z)|<=omega(|v-z|) give
R_g=O(gamma delta^2)+O(gamma delta omega(delta)). Hence w' is
integrable provided int_0 omega(delta) ddelta/delta is finite.
Positive Holder regularity suffices. The retained cold oscillation
sqrt(chi) yields exponent 1/4 in delta at fixed positive h, since
delta=h(cosh chi-1)/[(1+h)(1+h cosh chi)]. The dependence on h is
retained; this does not prove a uniform new inverse norm. With this
condition, R_g/Q tends to zero and the actual bolt row is

    f(z_h)/2-g(z_h)-2gamma w(z_h)/z_h^4.               (8)

The global moment w(z_h) cannot be discarded from a local source bound.
For general V* data, (6)-(7) are interpreted by the transpose pairing
itself, not by requiring a pointwise derivative of a distribution.

## 4. Physical bolt point source and a bounded exterior conversion

The physical functional is L_b(s,xi)=ell_b[s(z_h)+xi(z_h)/2],
ell_b=(1+h)/h. Since Q delta_(z_h)=0, the transpose of the canonical
bolt evaluation is exactly the TWO sources ell_b(delta,delta/2).
For the normalized equations in (2) the coefficient instead is

    ell_b/W=-G(1+h)gamma/(N pi B^2 h)
           =-G(1+h)K(1+h)/(N pi B^2), gamma=h K(1+h). (9)

Evaluate both expressions on the actual parameter enclosure for N=1,2.
The 1/h has not disappeared from the observable: W carries its compensating
factor. Adjoint pairing with the physical action is unchanged. The
previous scalar and TT singular profiles remain required, and no smooth
homogeneous adjoint is substituted.

For a source supported in 0<=z<=z_c=1/2, set
A(z)=int_0^z t^4/Q(t)^2 dt. Fubini in (7), first on smooth cutoffs
and then by integrable limits, gives

    w(z_h)=int_0^z_c [Q A f/z-(Q A)'g/2]dz,
    |w(z_h)|<=M_f ||f||_infinity+M_g ||g||_infinity,
    M_f=int_0^z_c Q A/z dz,
    M_g<=int_0^z_c [|Q'|A+z^4/Q]dz/2.                (10)

Thus the returned TT source value divided by h is at most
2gamma/(h z_h^4)(M_f||f||+M_g||g||). This is the exact h-sensitive
kinematic return, not the independent gravitational exterior-return target
<6.67e-10. For the actual seed, verify Q>.11, Q<.17 and |Q'|<.3 on
this whole interval. Even the elementary global bounds give

    M_f<.17 z_c^5/(25*.11^2),
    M_g<.3 z_c^6/(60*.11^2)+z_c^5/(10*.11).

The intended conversion coefficients are <.04 for f and <.08 for g.
Execute stronger whole-cell integrals on 128 and 256 cells using the
same positive kernel. These are constants for the canonical-to-York
source conversion ONLY; they are not K_m, K_c, or a solved metric error.

## 5. Global physical trial directions and independent transpose fixtures

For a polynomial Q=q0+sum_(k=2)^6 q_k z^k, the auxiliary direction
J_m=z^(m+4) produces a polynomial TT tangent

    xi_m=(m+4) z^(m+3) Q
           +(m+4)z^4 sum_(k=2)^6 k q_k z^(k+m-1)/(k+m-1),
    b_m=Q J_m'-xi_m.                                  (11)

This follows by integrating
(xi/z^4)'=[QJ''+(2Q'-4Q/z)J']/z^4. The missing q1 is essential for
the displayed polynomial primitive at m=0. A separate direction
xi=z^4 has J=0, b=-xi, S=xi/2, T=-xi. The m=0 direction has the
nonzero physical mass cubic 4Q(0)z^3; m>=1 starts at z^4 or faster.
All are bolt regular, with b(z_h)=-xi(z_h), and return to the fixed
outer clock. They need no artificial inner cutoff. This supplies actual
finite-seed weak trial functions; no density/completeness theorem for a
whole numerical approximation space is claimed, and logarithmic tails
must still be included where the actual equation requires them.

An independent exact transpose family is w=z^(m+5), f=0,

    g=-4(m+5)z^m Q
          +2m(m+5)int_0^z t^(m-1)Q(t)dt,              (12)

where the integral term is zero for m=0. It satisfies Lw=Qg'/2 and
the endpoint conditions. In particular m=0 gives g=-20Q and its
returned shape source is 25Q+zQ'. Verify the weak equality by integrating
both original canonical and returned polynomial pairings across the
entire actual interval, with the root boundary term retained exactly.
These are manufactured source fixtures, NOT an evaluation of r1 or H psi.

## 6. Continuum scope, decision and execution contract

[Dabrowski--Dossena](https://arxiv.org/pdf/1209.2021) establishes spin-Dirac
equivariance under diffeomorphisms on smooth manifolds. It does not alone
justify the paired outer domain or reflected C2,1 seam. Here the explicit
fixed-seed return, original boundary representative and vanishing endpoint
pairings supply the restricted weak conversion. The full quantum operator
still needs the retained finite-kernel/contact construction and actual
common curved/boundary remainder. The local lower-form estimate is not
an evaluated remainder or a numerical dual norm.

Decision: keep the complete angular logarithmic principal model and use
this frozen weak transpose in the physical primal-dual assembly. Do not
feed raw canonical rows into the TT equation, drop the global moment in
(8), or differentiate a moving chart without its off-shell contacts.
The next decisive calculation is unchanged: assemble the actual common-
subtracted curved/boundary remainder on the forward/bolt-adjoint trial pair
and evaluate the two residual norms and signed bolt pairing. The targets
<1.667e-9, exterior return <6.67e-10, Y0<=5e-9 and all preserved-ball
inequalities remain. That physical solve is still unexecuted.

Freeze argument, checker and source review before the run; retain all
attempts and snapshots. Use python-flint 0.9.0 at 256 bits, symbolic
identities, and direct polynomial antiderivatives on the whole finite
domain. Verify action measure, two-source transpose, endpoint identities,
mass-mode retention, both grid conversion targets and polynomial duality.
No scientific edit follows a successful run without full replay. No Lean
or proof of a controlled solution is supplied by this package.

Remaining to-do list: full common curved/boundary operator and physical
forward/adjoint residual certificate; X/admission, nonlinear closure and
exact compatibility; physical matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and stable wormholes.
