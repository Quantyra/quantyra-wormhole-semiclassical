# Physical curvature reconstruction and the charged optical reduction

2026-10-20, E003/S168. The actual seed admits a bounded reconstruction
of its conformal metric amplitude from the physical curvature variable
at sufficiently large real causal frequency. This is a weighted Green
map, not the earlier volume-norm inverse. Substitution gives an exact
curvature/shape equation with bounded conformal feedback and a compact
derivative feedback term. The actual retarded shape equation and its
physical domain still have to be solved. The full goal is incomplete.

Retain the action, state, charge and temperature of
[the explicit metric-source calculation](explicit-curvature-metric-source.md).
All background functions below are those of the same off-shell seed;
none is substituted for an exact nonlinear solution.

## 1. Curvature in optical variables

Use inward tortoise coordinate y, r_y=-F, and the optical perturbations
h^t_t=h^y_y=2(p+v), h^theta_theta=h^phi_phi=2p, h_ty=0.
These are the conformal and shape variables of the full retarded two-row
problem. Here v denotes shape, not a variation parameter. Put

    H_R=2R+rR', q_p=r p, u=r(delta R-r R' p),
    D=F'/r-H_R/6, K=-partial_y^2+F D,
    u_v=r[-2F^-1(-v_tt+v_yy)
                         +(2F''+4F/r^2+8F'/r)v].          (1)

Primes in radial coefficients mean r derivatives. Direct conformal
variation delta_p R=-6Box p-2Rp and the orbit-only curvature variation
give the exact kinematic identity

    u=6F^-1[(q_p)_tt+K q_p]+u_v.                         (2)

This uses no metric equation. In particular u_v contains the shape's
second time derivative, which must be kept when taking a Laplace transform.

For the rational seed write z=B/r, x=lambda B^2, mu=m/B,
qhat=1+3x+d and f=z^2F(B/z). The previous seed theorem gives
mu>=1+2x, x>1, and |d|<x/2. The last bound follows already from its
stronger 3x-1-34|d|>0 temperature-monotonicity guard. Then

    f=x+z^2-2mu z^3+qhat z^4-2d z^5+d z^6,
    D=B^-2 P(z),
    P=6x+2mu z^3-2qhat z^4+8d z^5-8d z^6 >=3x.        (3)

Indeed its classical x term is bounded below by 4x because
6+4z^3-6z^4-4=(1-z)(6z^3+2z^2+2z+2)>=0 on [0,1].
The remaining classical term is 2z^3(1-z)>=0, the thermal increment
is nonnegative, and the d term is -2d z^4(1-2z)^2>=-2|d|.
Thus D>=3lambda everywhere outside the nonextremal horizon.

K is the unique positive scalar realization: its potential is
6/y^2+O(1) at zero and decays at the horizon. Both endpoints are
limit point; the integrable outer homogeneous solution is O(y^3).
For real gamma>0, (K+gamma^2)^-1 has a positive Green kernel G_gamma.
It can be constructed by exhaustion with finite Dirichlet intervals
and the maximum principle. This is an auxiliary curvature-to-metric
operator, not the full gravitational response.

## 2. A global bound with the singular lapse retained

The required map on physical curvature L2(dy) is

    R_gamma u(y)=(1/6) integral G_gamma(y,s) F(s) u(s)ds. (4)

Multiplication by F is unbounded near y=0. Therefore (4) cannot be
justified by multiplying the ordinary unweighted resolvent estimate
by a nonexistent bound for F. Define it first on compact sources,
then extend by the kernel estimate below.

There is an explicit finite threshold. With z derivatives in this
paragraph define

    N=(2f-zf')^2-f[2f-2zf'+z^2f''+P-3x],
    Hpol=(N+x^2)/z^2,
    Gamma^2=(1+|Hpol|_1)/B^2.                           (5)

Hpol is a polynomial: N has constant term -x^2 and no linear term.
The coefficient norm |.|_1 is the sum of absolute coefficients, with
the actual fixed seed parameters inserted. Direct differentiation gives

    F'^2-F(F''+D-3lambda)
          =B^-2[-x^2/z^2+Hpol(z)] <=Gamma^2.            (6)

No smallness or sharpness of this threshold is claimed. For real
gamma>=Gamma the two positive supersolutions are

    (K+gamma^2)1=F D+gamma^2 >=3lambda F,
    (K+gamma^2)(1/F)=F''-F'^2/F+D+gamma^2/F>=3lambda.   (7)

On any finite Dirichlet interval, comparison with 1/(3lambda) and
1/(3lambda F) bounds the row and column integrals of the kernel in
(4) by 1/(18lambda). Passing through exhaustion preserves the bounds.
The Schur criterion therefore proves

    ||R_gamma||_(L2(dy)->L2(dy)) <=1/(18lambda),
                     gamma>=Gamma, gamma real.         (8)

[Teschl, Lemma 0.32 and its proof](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf)
is the established Schur criterion used here. The seed-specific
supersolutions and threshold are (3)-(7), derived in this note.
For compact input (4) agrees with the positive resolvent applied to
F u/6. Density extends it to all L2 inputs, and local elliptic estimates
show (K+gamma^2)R_gamma u=F u/6 in distributions. The resulting q_p
is the unique L2 solution: an L2 difference is a homogeneous solution
of the positive limit-point operator at eigenvalue -gamma^2.

This extended q_p need not lie in K's usual energy domain. For example
in the leading model F=1/(lambda y^2), q_p=y^(1/4)exp(-y) has both
q_p and u=6F^-1(K+gamma^2)q_p in L2, while q_p' and F u are not L2
near zero. This distinction is necessary for the physical norm; it
does not exclude such curvature inputs by an auxiliary energy choice.

For later derivative control choose ell=lambda^-1/2,
rho(y)=ell y/(ell+y), and M0=sup rho^2 F<infinity. Multiplying the
equation for q=R_gamma u by rho^2 conjugate(q) and integrating gives

    ||rho q_y||_2^2
      <=(M0/6)||u||_2||q||_2+||q||_2^2,
    ||rho partial_y R_gamma||
      <=[M0/(108lambda)+1/(324lambda^2)]^(1/2).          (9)

Here (rho^2)''/2<=1 and FD+gamma^2>=0. First use compactly truncated
weights, then pass to the limit; their second derivatives are uniformly
bounded and the endpoint L2 tails vanish. This proves the weighted
energy step without assuming q is in the unweighted energy domain.
Equations (8)-(9) are uniform on the stated real-frequency ray. No
bound near the imaginary frequency axis or causal metric inverse follows.

## 3. Exact feedback after conformal reconstruction

Use the previously derived c_R,d_R,e_R,A1,A0 and Atilde, without
changing any trace coefficient. A pure optical shape has the areal
representative a=2v, b=-2Fv, beta=-2v. Hence its remaining source is
exactly -2A0 v: the spatial-derivative terms cancel.

For a pure conformal p, the radial gauge displacement is xi^r=r p.
Removing the mixed component locally requires partial_r xi^t=r p_t/F^2.
The required local quantities, independent of the time integration
constant, are

    b=2rF p'-rF'p,
    a'=-r p''-r p_tt/F^2,
    b'+F a'=-r Psi_p/6-L p,
    Psi_p=-6Box p-H_R p,
    L=F'+rF''+rH_R/6.                                 (10)

The gauge need not be globally regular; these equations calculate the
gauge-invariant trace row and keep its off-shell contact. Set

    B1=2rF e_R,
    B0=2e_R F^2+rF e_R F'+F d_R L.

Adding conformal and shape variations gives exactly

    g_mix=-A1(u-u_v)/6-B1(q_p)_y-B0 q_p-2A0 v.         (11)

On the actual seed, B1=O(y^3), B0=O(y^2) at zero, and both decay
exponentially at the horizon. Explicit leading radial terms are
B1=-(88/3)lambda Qhat/r^3+O(r^-4) and
B0=-(176/3)lambda^2 Qhat/r^2+O(r^-3). In particular B1/rho and B0
are bounded and vanish at both endpoints. Define at gamma>=Gamma

    T_gamma=(B1 partial_y+B0)R_gamma,
    C_gamma=A1/6+T_gamma.                              (12)

Then C_gamma is uniformly bounded on curvature L2 by (8)-(9):

    ||C_gamma|| <= ||A1||_infinity/6
      +||B1/rho||_infinity [M0/(108lambda)+1/(324lambda^2)]^(1/2)
      +||B0||_infinity/(18lambda).                     (13)

For each fixed gamma, T_gamma is compact. On a compact interior
interval R_gamma maps a bounded L2 set into bounded H2, so both its
value and first derivative are compact in local L2. Outside such
intervals the norm is arbitrarily small, by (8)-(9) and the vanishing
coefficients B0 and B1/rho. Taking the norm limit proves compactness.
Multiplication by A1 is bounded, not asserted compact on L2.

For an exponential mode, or a transform with zero initial jets, define

    U_v(gamma)V=r[2F^-1(gamma^2 V-V_yy)
                             +(2F''+4F/r^2+8F'/r)V].

The exact curvature row after eliminating q_p is

    (Atilde+gamma^2+C_gamma)U
                  =C_gamma U_v(gamma)V-2A0 V-rF f_A.  (14)

It is meaningful in the stated spaces when U,V,U_v(gamma)V and rF f_A
are L2 and U has the required scalar graph-domain regularity. A natural
shape source class here is V in L2 with U_v(gamma)V in L2. This does
not impose the physical domain of the retarded shape operator. The
latter still contains the full nonlocal Dirac response; it is not
replaced by U_v or by the positive auxiliary K.

## 4. Initial contacts, compatibility and the remaining decision

For a general Laplace transform set I_q=gamma q_p(0)+q_p,t(0),
I_v=gamma v(0)+v_t(0), I_u=gamma u(0)+u_t(0), and

    u_v_hat=U_v(gamma)V-2r I_v/F,
    Q_p=R_gamma(U-u_v_hat)+(K+gamma^2)^-1 I_q.          (15)

All initial sources must belong to the relevant spaces; no free
choice of independent curvature and metric jets is inferred. With
J_q=(K+gamma^2)^-1 I_q the full version of (14) is

    (Atilde+gamma^2+C_gamma)U
      =I_u+C_gamma u_v_hat-2A0 V-rF f_A_hat
                                      -(B1 partial_y+B0)J_q. (16)

If both fast curvature traces are removed, the preceding scalar
range theorem imposes exactly

    Mtilde_gamma[I_u-C_gamma(U-u_v_hat)-2A0 V
                         -rF f_A_hat-(B1 partial_y+B0)J_q]=0. (17)

This keeps the original horizon functional for Atilde; it is applied
to the full reconstructed source. No perturbed Green kernel is silently
substituted, and (17) remains a projected trace equation. It is not an
independent physical condition after (16) and its boundary data hold.

The conformal feedback cannot cancel the inverse-square negative
quadratic directions merely by being included: for any fixed admitted
real gamma, (13) changes a normalized compact curvature test's real
quadratic value by a bounded amount, whereas the Atilde value tends
to minus infinity along the preceding outer tests. This is a statement
about the reconstructed curvature row with shape fixed. C_gamma depends
on frequency and need not be symmetric. It neither constructs a
full-system pole nor proves physical instability. The other metric row,
Ward constraints, admissible boundary data and exact background are still
required. No numerical physical collar or survival tolerance is supplied.

Next substitute the bounded reconstruction into the actual retarded
shape row and solve the remaining coupled problem on a justified
physical domain. This is now a concrete reduction in the physical
curvature norm, with all conformal feedback and initial contacts kept.

## Verification scope

Exact checks cover the kinematic identities, positive seed polynomial,
finite threshold, comparison functions, source cancellation and endpoint
coefficients. Independent four-dimensional curvature contractions test
general optical perturbations with nonzero time derivatives and the
off-shell radial gauge subtraction. Half-integer Bessel Green quadrature
tests the weighted reconstruction, derivative and Schur bounds in the
declared leading model, including an input outside the auxiliary energy
domain. Those fixtures are not an actual charged spectrum. The global
bound and compactness proof remain ordinary mathematics; no Lean is used.

Remaining to-do list: actual retarded shape reduction, physical initial
and boundary data, total moment and constrained stability transfer or
justified exclusion; S167 signed source, cold/higher response, nonlinear
solve and compatibility; exterior/state, charge/return, opening, full-sector
stability, transport, preparation, separate creation and final audit.
