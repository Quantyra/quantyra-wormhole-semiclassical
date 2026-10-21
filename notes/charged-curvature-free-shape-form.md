# The actual charged curvature-free shape form, including its tails

2026-10-20, E003/S168. The remaining shape row has an exact congruence
formula after physical curvature reconstruction. For shapes supported
in a fixed rescaled outer annulus, its actual charged curvature-free
form converges in the full logarithmic form norm to the positive normal
form. The accompanying conformal perturbation is global and its tails
are included. This proves positivity and a restricted variational inverse
for sufficiently small annuli, not a global physical metric inverse,
constraint-satisfying mode, or stable nonlinear wormhole.

Use the unchanged model, quantum state, common contacts and positive
temperature of [the physical reconstruction](physical-curvature-reconstruction.md).
The seed remains off shell. The positive normal comparison and actual
quantum convergence come from [the all-frequency estimate](curvature-free-all-frequency-bound.md)
and [the annular retarded theorem](annular-retarded-operator-limit.md).

## 1. The remaining row is a congruence, with transformed sources

At real positive Laplace frequency gamma, let H_gamma be the full
fixed-coordinate-density two-row response in optical variables (p,v).
Write its blocks H_pp,H_pv,H_vp,H_vv. Use the radial action pairing,
normalized per unit time-mode norm, including the sphere integration.
The local columns are the retained action Hessian; the canonical quantum
response occurs only in H_vv. Reciprocity makes this real-frequency form
symmetric on its core. This does not assert time-domain self-adjointness.

With the preceding notation set

    P_gamma=r^-1 R_gamma,
    Q_gamma=-P_gamma U_v(gamma),
    J_gamma=[[P_gamma,Q_gamma],[0,1]],
    p=P_gamma U+Q_gamma V.                              (1)

Here U is the physical curvature amplitude r Psi. Q_gamma V is the
unique L2(r^2dy) conformal amplitude making Psi zero, with the auxiliary
limit-point reconstruction. Define the transformed form by J* H J;
the stars denote adjoints in these form/dual pairings. Its shape block is

    S_gamma=H_vv+H_vp Q_gamma+Q_gamma* H_pv
                              +Q_gamma* H_pp Q_gamma,   (2)
    B_gamma=Q_gamma* H_pp P_gamma+H_vp P_gamma.

The transformed shape equation is

    S_gamma V+B_gamma U=b_v+Q_gamma* b_p.                (3)

Simply substituting p into the original v row omits the Q* trace-row
contribution. Equation (3) is equivalent to that substitution together
with the transformed trace equation on the admitted reconstruction range.
It is not a Schur complement using an inverse of H_pp or the fast factor.

All chart/density contacts remain in H. For nonzero initial jets the
kinematic reconstruction has an additional fixed term

    p_init=r^-1(K+gamma^2)^-1[I_q+(r/3)I_v],
    I_q=gamma q_p(0)+q_p,t(0), I_v=gamma v(0)+v_t(0).

Replace b by b-H(p_init,0) before applying J*. The source b itself must
already include every initial contact of the full retarded equations and
independent state source. This formula does not choose independent
curvature and metric jets. The map J is a fixed linear tangent map at the
seed; no nonlinear curvature-coordinate Hessian with omitted second
coordinate derivatives is being used.

## 2. An exact contact identity for evaluating the form

Let T be the retained equation trace, epsilon=G/(60pi), and
lambda0=lambda+11epsilon lambda^2/6. In optical variables its conformal
chart row per physical volume is

    H_Phi[p,v]=-[delta T+(4p+2v)T]/(8piG).

For Psi=0, delta T=epsilon g_mix/(rF)+r T' p. Therefore, with q=r p,

    Q_pp[p]+Q_pv[p,v]
       =-epsilon/(2G) integral q g_mix dy
        -1/(2G) integral F r^2[(rT'+4T)p^2+2Tpv]dy.     (4)

The same identity polarizes. Consequently S_gamma[v]=Q[p,v] with
p=Q_gamma v can be evaluated as Q_vv[v]-Q_pp[p] plus twice the
right side of (4). All residual contacts in (4) are kept at finite seed.
There is no assumption that the seed solves the physical equations.

For explicit tail control the complete conformal form is

    Q_pp[p]=(1/(2G)) integral {
       6epsilon (r^2/F)(p_yy-2F p_y/r-gamma^2 p)^2
       -r^2 k_o(p_y^2+gamma^2 p^2)
       -F r^2(2R+48lambda0)p^2 }dy,
    k_o=6+epsilon R/6+(22epsilon/3)j.                    (5)

For complex profile arguments use polarization of the real quadratic
formula. Equation (5) follows from the complete conformal Hessian,
including its density contact, integrated twice by parts. The endpoint
decay established below makes those integrations legitimate.

The trace itself is a known local function of the seed:

    T=R+12lambda0+epsilon[7Wc^2/36+22j^2/3
                                  -11R^2/72-20M+Box R],
    M=q_magnetic^2/(2r^4),
    Box R=F R''+(F'+2F/r)R'.                            (6)

Here q_magnetic is the fixed flux parameter, not the conformal amplitude
q. The vacuum relation cancels the constant term. Thus T=O(r^-4),
rT'+4T=O(r^-5) at infinity; both are smooth and bounded at the positive-
temperature horizon. These contacts vanish in the scaled limit for a
proved reason; they have not been removed from (2) or (4).

## 3. Global reconstruction under outer dilation

Use the physical inward tortoise coordinate and time, and set

    y=d eta, gamma=kappa/d,
    v_d(y)=d^-1/2 v(eta), supp v in I=(ell,b),
    p_d=Q_(kappa/d) v_d,
    z_d(eta)=lambda d^(3/2) r(d eta) p_d(d eta),
    R_d(eta)=lambda d r(d eta).

Keep 0<kappa0<=kappa<=kappa1 fixed. The reconstruction equation is

    (K_d+kappa^2)z_d=f_d,
    K_d=-partial_eta^2+d^2 V_K(d eta), V_K=F D,
    f_d=R_d(v''-kappa^2 v)/3
           -R_d d^2 F(d eta) C_v(d eta) v/6,
    C_v=2F''+4F/r^2+8F'/r.                              (7)

The forcing is compact in I, although z_d is not. The seed asymptotics
give R_d=eta^-1+O(d^2) in every required coefficient norm on I and
d^2 F C_v=24/eta^2+O(d^2). Thus

    f_d -> f_0=B_kappa v/(3eta^3) in L2,
    B_kappa=eta^2(partial_eta^2-kappa^2)-12,
    ||f_d-f_0||_2<=C_I d^2||v||_H2.                     (8)

There is a GLOBAL potential comparison. The preceding seed formulas give
V_K=6/y^2+2lambda+O(y) at zero and V_K->0 at the horizon. Hence

    W_K=V_K-6/y^2 is bounded on (0,infinity),
    K_d=K_0+d^2 W_K(d eta), K_0=-partial_eta^2+6/eta^2,
    ||(K_d+kappa^2)^-1-(K_0+kappa^2)^-1||
                              <=d^2||W_K||_infinity/kappa0^4. (9)

Both K_d and K_0 are nonnegative. Their limit-point domains agree after
the bounded perturbation. The estimate is the second resolvent identity,
whose established framework is [Teschl, Lemma 6.5](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf).
It controls output over the whole half-line, not just annular pairings.
Equations (8)-(9) give z_d-z_0=O(d^2)||v||_H2 in L2 and in the graph
norm of K_0. The graph upgrade follows by applying K_0 to the difference
equation. On the compact core,

    ||K_0 z||_2^2=||z''||_2^2+12||z'/eta||_2^2,

and the weighted Hardy inequality controls z/eta^2 as well. Closure
extends these estimates to the resolvent outputs.

The limiting conformal profile is p_0=eta z_0=-T_kappa v, where
T_kappa=(1/3)A_kappa^-1 B_kappa is precisely the previously separated
normal auxiliary inverse. It retains its regular outer and decaying
interior tails. No truncation of that inverse is made.

## 4. Why the entire conformal action converges

The following bounds are uniform for v in a bounded H2(I) set and kappa
in the stated compact interval. Near zero the regular solution of (7)
is O(eta^3), with its first two conormal derivatives; hence the scaled
conformal amplitude p_d(d eta)d^(1/2)=z_d/R_d is O(eta^4). Its coefficient
converges to the normal coefficient. To see this uniformly, use the
regular Frobenius solution below ell/2, the bounded remainder in (9),
and the graph-norm convergence of value and derivative at ell/2. The
regular branch is the only integrable one there.

For eta>b, positivity and comparison with the free massive Green kernel
give

    |z_d(eta)| <= ||f_d||_1 exp[-kappa(eta-b)]/(2kappa).

The homogeneous ODE bounds its first two derivatives by a fixed multiple
of the same exponential for eta>=b+1. For any fixed small physical y_c,
on b+1<=eta<=y_c/d the differentiated seed estimates give
R_d comparable to eta^-1. The rescaled integrands of (5), and the
contacts in (4), are bounded by a fixed polynomial in eta times
exp(-2kappa0 eta). Near zero they are bounded by integrable positive
powers of eta. These estimates include both derivatives in (5).

On the rest of the exterior, eta>=y_c/d, the actual rational seed and
simple horizon imply that all coefficients in (4)-(5) grow at most
as a fixed power of d^-1 times exp(C d eta). The conformal amplitudes
and their differentiated coefficients have the same kind of bound
times exp(-kappa0 eta). Choosing d small makes the integral of this
region at most C' d^-N exp(-kappa0 y_c/d), which tends to zero.
This keeps the finite-temperature horizon tail; no zero-temperature
limit or uniform smallness of the geometry on the entire exterior is used.

On compact eta intervals, (8)-(9) and the differentiated seed expansions
give convergence in every local norm needed in (4)-(5). The uniform
endpoint bounds just proved then give convergence of their full bilinear
forms, uniformly on the H2 unit ball. The contacts in (4) tend to zero
after multiplication by lambda d^4, using (6) and the exact source
g_mix=A1 u_v/6-B1 q_y-B0 q-2A0 v on Psi=0. All local conformal and mixed
terms therefore converge with their global tails, not just pointwise.

## 5. The actual shape limit and its coercivity

Let X_I be the closure of compact smooth shape profiles in

    ||v||_X^2=(1/2) integral (1+xi^2)^2[1+log(2+xi^2)]
                                           |v_even_hat(xi)|^2 dxi.

Only the canonical quantum V,V term is nonlocal in the original variables.
Its two shape arguments remain compact in I, so the preceding annular
retarded theorem supplies convergence in this norm. Combining it with
section 4 yields the actual result

    sup_(kappa in [kappa0,kappa1])
      ||lambda d^4 S_(kappa/d)[v_d,w_d]-R_kappa[v,w]||
                                            _(X_I -> X_I*) ->0. (10)

The notation in (10) denotes the polarized profile form. The physical-
coordinate normalization is lambda d^4. In the preceding canonical
coordinates y=L Y, t=L tau, L=ell_o/sqrt(x), the profile amplitude
and time-norm conversion give C L^2=1/lambda, proving this normalization.
It also follows directly from (5)'s leading quartic coefficient.

R_kappa is the already established positive curvature-free normal form.
Positivity in H2 alone would not absorb an error measured in X. Here is
the needed stronger bound. Write h(v)=||v||_H2^2 in the same even-half
normalization and set

    m=min(1,kappa0^2), M=max(1,kappa1^2), c0=min(1,kappa0^2/2),
    a0=19m^2/(150pi), A=m^2/(60pi),
    B=A+[M^2/(30pi)] max(0,-Lbar+6/5-log ell-(log c0)/2),
    c_I=a0 A/(a0+B)>0.                                (11)

The previous theorem gives R_kappa>=a0 h. Its explicit logarithmic
formula also gives R_kappa>=A||v||_X^2-B h, since
log(sqrt(xi^2+kappa^2))>=(log(2+xi^2)+log c0)/2 and
log eta>=log ell on the support of H_kappa v. Combining these two
bounds with weight a0/(a0+B) proves

    R_kappa[v]>=c_I||v||_X^2.

Thus for sufficiently small d, uniformly on this compact frequency set,

    lambda d^4 S_(kappa/d)[v_d]>=c_I||v||_X^2/2.         (12)

The bounded coercive form on X_I has its unique variational inverse to
X_I*, with norm at most 2/c_I in the normalized profile variables.
This is a restricted shape problem with a globally reconstructed
conformal tail. Completing X_I introduces no physical wall at ell or b.
The quantum comparison constants remain uncalibrated, so (12) supplies
no numerical physical collar radius. It does not prove positivity or
invertibility for arbitrary shapes occupying the whole exterior.

## 6. What remains physically unresolved

Equation (2) specifies the actual shape operator and (10)-(12) control
its escaping-annulus restriction. The U,V cross block B_gamma, the
global shape response and the physical metric/state domain still have
to be controlled to solve (3). A restricted inverse cannot be substituted
for the global inverse in a Schur complement. The trace retains the fast
branch and its total horizon compatibility condition. The Ward initial
and incoming constraints and the off-shell-to-on-shell transfer remain.

No actual charged growing mode, complete stability theorem, nonlinear
exterior, cutoff or order reduction is obtained. S168.3-S168.4 and all
S167 tasks remain open. The next task is the global reduced shape and
cross-block estimate, or a physical domain argument that excludes the
fast branch under the original objective.

## Verification scope

The checker verifies the congruence/source algebra, density-contact
identity, potential and trace asymptotics, coordinate normalization,
normal reconstruction identity and logarithmic coercivity constants.
Actual-form seed ODE fixtures compare the globally reconstructed scaled
conformal field with an independently integrated Bessel Green solution,
refine the boundary-value solve and dilation, and retain explicit tail
estimates. These fixtures do not evaluate the full Dirac shape spectrum.
The continuum transfer uses the pinned earlier quantum theorem; it is
ordinary mathematics, not a Lean certificate. No Lean is selected.

Remaining to-do list: global reduced shape and cross-block control,
physical initial/boundary data, total compatibility and charged stability;
S167 signed source, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, full-sector stability, transport,
preparation, separate creation and final audit.
