# Complete structured residual: retain the curvature cancellation before norms

2026-10-27 HST, E003/S169. Ordinary argument before arithmetic. The user
resumed the unchanged full-domain residual goal. The preceding turn was
NO PROGRESS and recorded an implementation impasse. This calculation tests
the proposed finite-output implementation on the ACTUAL structured range,
in the ORIGINAL energy residual norm, not the stronger H residual norm.
The full stable-wormhole objective, closed audit and accepted results stay.

## Certificate under test and its domain

Keep E_t with norm ell_t/N, N=|W_action| beta, R_t=ell_t/N and
T_t=R_t^-1 B_N, B_N=b* a^-1 b/N. The positive self-adjoint T_t has
norm <=q_t=.338/(1-epsilon). The target remains a full implemented
R_n=T_t-(I+T_t)F_n with all quantum/application errors; 1e-30 is the
recorded sufficient residual screen, not a necessary existence condition.

Test fixed output ranks n=10^6,10^7,10^8. These are implementation classes,
not a finite selection of physical modes. All coefficients and quantum
actions are granted exactness in the LOWER-bound test. A rejection would
therefore also apply to numerical implementations in the same range.
This does not test P32, which is a polynomial in a full operator and need
not have finite-dimensional range.

Use exactly the paired mass/log completion and core of the accepted
strong-residual rank test. Pi_- psi=0 for D, Pi_+ phi=0 for its adjoint,
the complementary Robin square, original thermal spin structure, common
regulator, exact ultraviolet logarithm, all angular/thermal sums and
coincidence contacts remain. York has J(0)=0 and (Q^2 J'/z^4)(z_h)=0;
the bolt is regular, not a finite wall. No nonlinear neighborhood or
independent phase control is asserted.

For M=n+1 cells in I=[.4,.6], h=.2/M, let phi_j=(1-x_j^2)^8 inside
its cell and zero outside, x_j=2(z-center_j)/h. Set

    m=(Q'/z^4)', K=Q(m phi'+2m' phi), w=z^4 m^2 phi,
    J=K/Q, xi=K'-w, nu=xi/2+K/z, T=-xi-K'/2.           (1)

The accepted exact York and zero-moment identities apply. These fields
are compact admissible directions INSIDE the full paired domain. They
retain every cell transition and vanish at both physical ends, so their
original boundary Hessians vanish by evaluation, not operator deletion.
There is no arbitrary-source surrogate. The full scalar inverse remains
in B_N. The disjoint xi_j span an M-dimensional subspace U_M of E_t.

## First loss: bound the actual source without a global coefficient ratio

Use f_k=||partial_x^k (1-x^2)^8||L2(-1,1), k=0,...,5, from exact
polynomial integration. Write xi=sum_(i=0)^2 a_i phi^(i), with a2=Qm>0,
and put omega=2/h. On a coefficient interval [l,r] containing an entire
cell, define

    N_k(c)=sum_i sum_(j=0)^k binom(k,j)
                 ||partial^j c_i||inf omega^(i+k-j) f_(i+k-j),
    D0=(inf a2) omega^2 f2-||a1||inf omega f1-||a0||inf f0,
    D1=(inf a2) omega^3 f3
                 -sum_(i=0)^2 ||(xi')_i||inf omega^i f_i.

Execute D0,D1>0. The leading coefficient of xi' is again a2. Exact
triangle inequalities, with no sampling, give

    ||field^(k)||L2(dz)/||xi||H <= r^2 N_k(c)/D0=:L_field,k,
    ||gradbar xi||H^2/||xi||H^2
      >= (inf Q) l^4/(beta r^2) (D1/N0(xi))^2=:g_j,
    ||gradbar xi||H/||xi||H <= sqrt(sup Q)/(sqrt(beta) l) L_xi,1,
    ||Deltabar xi||H/||xi||H
      <= [(sup Q)L_xi,2+||Q'-2Q/z||inf L_xi,1]/beta.    (2)

Here H=L2(z^-4 dz). In the last inequality the z^2 in Deltabar cancels
the H weight BEFORE taking a norm. This uses both original normalizations.
Cover cell centers by 128 macro intervals, expand each by h/2 and clip
to I. Every cell is contained in at least one resulting coefficient
interval. Whole-interval Laurent evaluation supplies all extrema; it does
not enumerate a million cells or substitute representative sample fields.

Let A_j=6d L2_j^2+6.01 L1_j^2+24.01, d=D_R lambda. The complete scalar
energy obeys a(xi_j,xi_j)/N<=A_j||xi_j||H^2. Keep the exact main mixed
identity and its full known remainder K_rest:

    b_main(xi_j,xi_j)/N >= p_j||xi_j||H^2,
    p_j=4[C_v inf(tau/lambda)+d inf(tau/lambda) g_j
                          -(d/2)||Deltabar(tau/lambda)||inf],
    |b_rest(s,u)|/N<=K_rest ||s||A ||u||H.              (3)

The complete scalar row is local; different compact cells have zero
scalar cross energy. For u=sum u_j, take the scalar trial
s=sum (p_j/A_j)u_j. Its a energy is at most
sum (p_j^2/A_j)||u_j||H^2 and its main mixed pairing is at least the
same sum. Thus the full dual main-source norm is >=sqrt(b0)||u||H,
where b0=min_j p_j^2/A_j. Triangle inequality for the COMPLETE remainder
gives, provided sqrt(b0)>K_rest,

    B_N[u,u] >= b_* ||u||H^2,
    b_*=(sqrt(b0)-K_rest)^2.                           (4)

This retains the full a^-1 and all global mixed/York returns. It is a
variational LOWER bound on the actual coupled source, not a scalar solve
approximation. Use the minimum over the whole macro cover in (4).

## Second loss: an exact two-row curvature cancellation

The previous independent-jet estimator bounded the full R2 Hessian using
nu'',T'' separately. On (1) their combination has lower order. Put

    d1=-2 Delta T-6 Delta nu+12Q nu'/z,
    Dseed=-Q''+6Q'/z-12Q/z^2, Delta=(Q partial_z)'.

Since T+3nu=-w/2+3K/z and xi'=K''-w', direct differentiation gives

    d1=Q w''+(Q'-6Q/z)w'
                +(24Q/z^2-6Q'/z)(K'-K/z).             (5)

The K'' terms cancel exactly. All terms in (5) are at most second
derivatives of phi. This is an identity on the actual fields, including
their cutoff transitions. No integration by parts, coordinate change,
new chart, altered domain or discarded off-shell acceleration is used.

The COMPLETE affine R2 diagonal from the accepted local generator is

    alpha_R {16 T^2+4C d2+C^2[
        4T^2 Dseed^2-8Dseed T d1+2d1^2+2Dseed d2]},
    alpha_R=8pi c_R/(C B_radius^2), d2=-12Q(nu')^2.     (6)

Use (5) to enclose ||d1||L2/||xi||H from its three coefficient functions
before bounding (6). In particular, with Z=||Dseed||inf, Ld the new
d1 norm ratio, and L_T,0,L_nu,1 from (2), the safe full R2 upper bound
per H norm squared is

    alpha_R/beta * [16 L_T,0^2
      +(4C+2C^2 Z)12(sup Q)L_nu,1^2
      +C^2(4Z^2 L_T,0^2+8Z L_T,0 Ld+2Ld^2)].          (7)

The original action sign is handled by this absolute upper bound.
All terms in (6), including the density/second-vertex term d2, stay.
The disappearance of the artificial fourth-derivative loss is algebraic,
but its quantitative consequence is tested in the COMPLETE energy below.

## Complete quantum energy and the original residual consequence

Reuse the accepted whole-I coefficient matrices for Einstein, magnetic
and Weyl terms; combine them with the sharper cell norm ratios (2).
Use (7) for R2. Add the SAME a2, LLL and shared vacuum allowances and the
complete positive geodesic logarithmic bound of the strong-rank proof:

    Egeo=2(4C^2/q)/(525pi B_radius^2 beta) (sup Q)
            [5A^2/8+3AD/2+D^2],
    A=(sup Q)L_T,3, D=(3/2)||Q'||inf L_T,2
                                      +(1/2)||Q''||inf L_T,1.

This bounds ALL angular levels and full geodesic/thermal returns, not a
selected spectral set. The actual-to-geodesic error and retained transport
give

    ||u||Et^2 <= E_* ||u||H^2,
    E_*=(1+epsilon)/(1-.0022)
          max_j(E_local,j+E_a2,j+E_LLL,j+E_vac,j+Egeo,j). (8)

All derivative-square bounds sum over disjoint supports. The nonlocal
energy is bounded by its full local third-flow-derivative majorant, so
cross-cell quantum pairings are retained in (8). Source/core approximation
is exactly that of the accepted strong-rank proof; no new smooth-seam
theorem is needed. The common regulator and contacts stay in its comparison.

For any F_n with fixed output dimension <=M-1, choose a nonzero
u in U_M E_t-orthogonal to (I+T_t) range(F_n). Equations (4),(8) give

    ||T_t-(I+T_t)F_n||Et->Et >= mu_M=b_*/E_*.          (9)

For S_hat with fixed output dimension <=M-1, impose instead the at most
M-1 conditions <A_Q range(S_hat),u>=0. All inputs have chi_outer u=0,
so B_i u=B_N u. The ORIGINAL weak structured residual satisfies

    ||B_i-A_Q S_hat||Et->Et* >= mu_M.                 (10)

This is not the previous stronger H-residual rejection or the arbitrary-
source uniform-inverse obstruction. Execute (9) against 1e-30 and (10)
against 4.9887427328e-32 at all three preselected ranks. A lower bound
below a screen proves neither success nor attainability.

If a routine has additional output L outside its fixed range, (9) weakens
by at most (1+q_t)||L||Et->Et. Thus an analytic complementary action must
have norm at least (mu_M-1e-30)/(1+q_t) when that expression is positive.
This is a necessary size, not an upper error certificate. Increasing
precision or correcting entries inside the same range cannot repair (9).

Decision rule: if the test rejects practical finite ranges, retain the
full logarithmic principal action and require a constructed analytic
complement with its complete weak defect. Do not infer that the inverse,
P32 or the physical construction fails. If inconclusive, report that
without treating a conditional cutoff or margin as an attained bound.

The proposed alternative remains a solver with analytic principal range
and finite correction. It must pass the same original weak residual test;
no implementation or quantum-application error is supplied by this rank
argument. The complete signed residual goal is not reduced to this test.

Remaining to-do list: implemented full-range S/K, their complete residual
and all application errors; signed G_t and original physical/source/graph/
mass/phase/admission budgets; nonlinear existence, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
