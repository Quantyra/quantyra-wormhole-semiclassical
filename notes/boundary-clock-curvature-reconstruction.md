# Curvature reconstruction with the boundary clock retained

2026-10-21, E003/S168. The boundary contribution exposed by the
[global gauge calculation](global-optical-gauge-and-boundary-clock.md)
can be included in the curvature reconstruction explicitly. A weighted
Green estimate admits its non-L2 conformal amplitude, an exact boundary
subtraction reduces smooth boundary data to the previous L2 inverse,
and the resulting curvature forcing is bounded. This completes a
kinematic reconstruction and bulk-curvature coupling calculation. It
does not supply the boundary shape/state response or a physical mode.

Keep the original action, flux, clock, temperature and quantum state.
All geometry is the same off-shell rational seed. The physical curvature
amplitude U=r Psi remains in its previously derived L2(dy) space. The
new weight below is an AUXILIARY reconstruction space, not a replacement
for that physical norm or a selected gravitational boundary condition.

## 1. A global weighted extension of the positive reconstruction

Use r_y=-F, z=B/r, x=lambda B^2 and the previously established

    K=-partial_y^2+F D, D>=3lambda,
    D=F_r/r-(2R+rR_r)/6,
    R_gamma u=(1/6) integral G_gamma(y,s) F(s)u(s)ds.

G_gamma is the positive Green kernel of K+gamma^2 for real gamma>0.
Let Z=L2(z^2 dy), with norm ||u||_Z=||z u||_2. Since z~B lambda y at
the outer end and tends to a positive constant at the horizon, Z admits
1/y but excludes the homogeneous reconstruction branch 1/y^2.

There is an explicit finite frequency threshold. Write

    f=z^2 F(B/z), P=B^2 D,
    N1=f[(15/4)f-(3/2)z f_z+x-P],
    N2=(2f-z f_z)^2-f(2f-2z f_z+z^2 f_zz)
                       +(1/2)f(2f-z f_z)+(3/4)f^2+xf-fP,
    H_i=(N_i+5x^2/4)/z^2,
    Gamma_Z^2=[1+max(|H_1|_1,|H_2|_1)]/B^2.           (1)

For the actual polynomial seed both H_i are polynomials; |.|_1 is the
sum of absolute coefficients after the fixed parameters are inserted.
Define h=z^(-3/2), b=z^(1/2)/F. Direct differentiation gives

    h_yy/h=(3/2)F F_r/r+(3/4)F^2/r^2,
    b_yy/b=F_r^2-F F_rr+(1/2)F F_r/r+(3/4)F^2/r^2.    (2)

The two quantities h_yy/h+lambda F-FD and
b_yy/b+lambda F-FD are respectively B^-2 N_i/z^2. Equation (1)
bounds them above by Gamma_Z^2. Thus for gamma>=Gamma_Z,

    (K+gamma^2)h>=lambda F h,
    (K+gamma^2)b>=lambda F b=lambda z^2 h.             (3)

Apply the positive finite-interval comparison principle and then
exhaust the half-line. In the measure z^2dy the kernel is
G_gamma(y,s)F(s)/(6z(s)^2). The first inequality in (3) bounds its
h-weighted row integrals by h/(6lambda); the second bounds its
h-weighted column integrals by the same quantity. Weighted Schur gives

    ||R_gamma||_(Z->Z)<=1/(6lambda).                  (4)

For clarity, weighted Schur here follows directly by applying Cauchy-
Schwarz to the factors sqrt(kernel*h(s)/h(y)) and
sqrt(kernel*h(y)/h(s))*u(s), then integrating in y. This is the
established [Schur criterion, Teschl Lemma 0.32](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf),
with the seed-specific supersolutions (2)-(3). No boundedness of the
singular multiplication F on unweighted L2 is presumed.

Density and local elliptic regularity extend the integral to Z and
give (K+gamma^2)R_gamma u=Fu/6 in distributions. It is the unique Z
solution. A homogeneous difference has outer exponents 3 and -2;
the second is excluded by Z. The first and the horizon L2 branch lie
in the old positive K realization, which has no eigenvalue -gamma^2.
On old L2 inputs this extension agrees with the preceding R_gamma.
Use gamma>=max(Gamma_Z,Gamma_R) whenever also invoking that note's
unweighted norm estimate; Gamma_R is its explicit threshold.

## 2. Derivatives and the full conformal feedback

Set M_F=sup z^2 F and
M_w=max(0,sup (z^4)_yy/(2z^2)). Both are finite on the actual seed:
z_y tends to B lambda at zero, z_yy is bounded, and both derivatives
decay at the horizon. Multiplication of the reconstruction equation by
z^4 conjugate(q), with q=R_gamma u, yields

    ||z^2 q_y||_2^2
       <=(M_F/6)||u||_Z||q||_Z+M_w||q||_Z^2,
    ||z^2 partial_y R_gamma||_(Z->L2)
       <=C_D=(1/6)sqrt(M_F/lambda+M_w/lambda^2).        (5)

Prove this first for compact sources on finite Dirichlet intervals.
The nonnegative potential term is retained, and the weight derivative
is bounded as displayed. Exhaustion, local convergence and weak lower
semicontinuity give (5). No unweighted energy regularity of q is assumed.

Use the exact coefficients A1,A0,B1,B0 from
[the physical reconstruction](physical-curvature-reconstruction.md).
Their outer orders are A1=O(y), B1=O(y^3), B0=O(y^2), A0=O(1), and
all decay at the horizon. Consequently

    C_gamma^Z=A1/6+(B1 partial_y+B0)R_gamma : Z -> L2,
    ||C_gamma^Z||<=||A1/z||_infinity/6
          +||B1/z^2||_infinity C_D+||B0/z||_infinity/(6lambda). (6)

This bound is uniform on the stated real frequency ray. The derivative
feedback is compact for each fixed gamma: on compact interior intervals
the output is locally H2, and at both ends the coefficients B1/z^2
and B0/z tend to zero. This proves compactness by local compact embedding
and (4)-(5). The multiplication A1/6 is not claimed compact.
On physical curvature L2, (6) restricts to the original C_gamma.

## 3. An exact boundary subtraction using only the old L2 inverse

Let V be twice continuously differentiable near zero, with trace v_0,
and let chi be a smooth cutoff equal to one there and zero sufficiently
far into the exterior. Assume the remaining fields have the stated
source-space membership at the horizon; no physical boundary domain
is selected by this smoothness example. Put

    V_reg=V-v_0 chi,
    H_R=2R+rR_r,
    C_v=2F_rr+4F/r^2+8F_r/r,
    U_v(gamma)V=r[2F^-1(gamma^2 V-V_yy)+C_v V].

The simple pair (p,v)=(-chi,chi) has exact curvature amplitude

    u_chi=4r F^-1(chi_yy-gamma^2 chi)-12chi_y
                                      +r(C_v+H_R)chi. (7)

For F=lambda r^2+1-2m/r+Qhat/r^2+c3/r^3+c4/r^4,

    C_v+H_R=4/r^2+10c3/r^5+36c4/r^6.                 (8)

Hence u_chi=4(lambda-gamma^2)y+O(y^3) near zero; it is L2 and
compact at the horizon. Also U_v V_reg is locally L2 because V_reg=O(y).
For U in physical L2, with U_v V_reg in L2 globally, the full exact
reconstruction is

    q=r p=-v_0 r chi
                  +R_gamma[U-U_v V_reg-v_0 u_chi].     (9)

Equation (9) is a boundary lift plus the previously controlled L2
operator. It is not a new homogeneous fast-scalar solution. Indeed
(K+gamma^2)(-r chi)=F(u_chi-U_v chi)/6, so direct substitution proves
the curvature identity. Uniqueness in Z shows that (9) is the same
as the extended formula R_gamma[U-U_v V].

It is independent of the chosen cutoff. If chi changes by a compact
interior function, the induced change in the bracket in (9) is exactly
6F^-1(K+gamma^2) applied to the compensating lift; the old inverse
cancels it. No artificial wall or additional boundary physics is added.

The old Green estimates imply R_gamma g=O(y^-1/2)||g||_2 locally for
g in L2, which follows by Cauchy-Schwarz against the regular y^3 and
singular y^-2 fundamental solutions. Therefore y R_gamma g->0 and

    p->-v_0,
    q=-v_0/(lambda y)+o(1/y).                          (10)

For the zero-initial clock displacement c, v_0=-gamma c and (10)
reproduces the boundary contribution found in the gauge theorem.
It is retained rather than set to zero. On general initial data the
time-domain relation is v(t,0)=-partial_t c(t), so the transform is
v_0=-gamma c_hat+c(0); the initial boundary clock is not discarded.

## 4. The boundary contribution in the physical curvature equation

For zero initial jets the exact charged trace row now extends to
V in L2 with U_v V in Z:

    (Atilde+gamma^2+C_gamma)U
       =C_gamma^Z U_v V-2A0 V-rF f_A.                 (11)

All terms on the right are in physical L2 when rF f_A is in L2.
The new weight has not changed U's kinetic norm, Atilde's boundary
realization, the fast branch, or the total signed curvature condition.
It only makes the previously omitted shape-boundary forcing meaningful.

Equivalently write V=V_reg+v_0 chi. The curvature-free boundary column is

    q_boundary=-r chi-R_gamma u_chi,
    g_boundary=C_gamma^Z U_v chi-2A0 chi
      =A1 U_v chi/6-(B1 partial_y+B0)q_boundary-2A0 chi. (12)

Then (11) has the extra source v_0 g_boundary. This is a derived bulk
curvature source, not a freely adjustable physical source. Its value
depends on the actual boundary/metric problem; the total expression is
cutoff-independent together with V_reg. The remaining shape equation,
its boundary action and quantum-state contacts still have to determine it.

Nonzero initial jets can be retained without ambiguous products. With
I_q=gamma q(0)+q_t(0), I_v=gamma v(0)+v_t(0), I_u=gamma u(0)+u_t(0),
use the single extended source

    q_hat=R_gamma[U-U_v V+(6I_q+2rI_v)/F],
    u_v_hat=U_v V-2rI_v/F,
    (Atilde+gamma^2+A1/6)U
       =I_u+A1 u_v_hat/6-(B1 partial_y+B0)q_hat
                                      -2A0 V-rF f_A_hat. (13)

The bracket in (13) must belong to Z and the displayed trace forcing
to physical L2. Horizon membership involving 1/F is a real requirement;
it is not automatic for arbitrary jets. Where the earlier separated
resolvent expressions are defined, (13) agrees with them exactly.

## 5. A uniform bound on the clock-induced curvature source

There is also a useful large-real-frequency consequence. For a fixed
smooth profile chi supported in [0,b], equal to one near zero, choose
chi_gamma(y)=chi(gamma y) as a decomposition profile. On its support,
F^-1<=C y^2 and C_v is bounded for gamma sufficiently large. Exactly,

    z U_v chi_gamma
      =B[2 gamma^2 F^-1(chi-chi'')+C_v chi](gamma y),

where the geometric coefficients are evaluated at y. Thus

    ||chi_gamma||_2=gamma^-1/2 ||chi||_2,
    ||U_v chi_gamma||_Z<=C_chi gamma^-1/2,
    ||g_boundary,gamma||_2<=C_b gamma^-1/2.            (14)

C_b is independent of gamma on a sufficiently large real ray, by (6).
For a gauge forcing f_+ the previous clock estimate is
|c|<=||f_+||_2/sqrt(2gamma). Since v_0=-gamma c in the zero-initial
case, (14) gives the full rank-one composition bound

    ||v_0 g_boundary,gamma||_2<=C_b ||f_+||_2/sqrt(2). (15)

The boundary contribution is controlled; it is not assumed negligible
or compatible with the fast curvature moment. These frequency-dependent
cutoffs are only used to estimate an exact decomposition. They do not
define a causal time-domain projection, and (4)-(15) are not a full
right-half-plane metric resolvent or a physical stability theorem.

## Verification and next decision

The ordinary proof above supplies the global estimates. Focused checks
verify the seed polynomials, both weighted supersolutions, derivative
weights, the exact boundary lift and cutoff cancellation, initial
contacts, and norm scaling. Independent Bessel Green quadratures test
the non-L2 boundary reconstruction in the normal model. Actual charged
seed boundary-value fixtures compare two numerical refinements and two
subtraction cutoffs for the SAME shape source, retaining the resulting
curvature feedback. Finite-domain numerical boundaries are convergence
fixtures, not physical walls or a numerical replacement for the proof.
No Dirac kernel, physical boundary action or full metric spectrum is
evaluated. Freeze scientific sources before acceptance; retain every
failed source snapshot and raw outcome. No Lean is selected.

The next physical step is now explicit: include the column (12) and its
transported boundary/state response in the global shape problem, or
derive why the admitted physical data have v_0=0. Preserve the actual
curvature compatibility moment, Ward/initial constraints, all off-shell
contacts and the eventual exact-background transfer. No physical gate
is advanced merely by extending this auxiliary inverse.

Remaining to-do list: physical boundary shape/state response and clock
closure, global coupling and charged stability; S167 signed source,
higher response, nonlinear existence and compatibility; exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and final audit.
