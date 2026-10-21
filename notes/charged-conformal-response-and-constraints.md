# Charged conformal response and Lorentzian constraints

2026-10-20, S168. Ordinary argument before focused verification.

The exact AdS Klein--Gordon lift is not an identity on a non-Einstein
charged exterior. This note derives its missing local terms, identifies
the gauge-invariant spherical curvature perturbation, and derives the
dynamic constraint propagation equations. It supplies part of S168.3;
it does not close the coupled retarded problem or decide stability.

The [vacuum preflight](ads-fast-scalar-stability-preflight.md) remains
valid in its stated scope. The earlier [static conformal response](stationary-coupled-asymptotic-obstruction.md)
already establishes conformal locality. The additional result here is
the causal column with its off-shell contacts and the precise failure
of direct Klein--Gordon transfer and static constraint recovery.

## 1. Model and causal state prescription

Retain the complete selected massless Dirac action, its full angular
tower, common subtraction, finite coefficients, magnetic flux, paired
domain and prepared temperature. Put epsilon=G/(60pi),
lambda0=lambda+11epsilon lambda^2/6, E=1/epsilon+11lambda/3. Use
Lorentzian signature (-+++) and Box=g^ab nabla_a nabla_b. Define

    calE_ab=G_ab-3lambda0 g_ab-8piG(T_Maxwell,ab+T_Dirac,ab),
    T=-g^ab calE_ab
     =R+12lambda0+epsilon[7C^2/12+11J^2/6-11R^2/72
                          -20M+Box R],                       (1)
    J_ab=Ric_ab-Rg_ab/4, M=F_ab F^ab=q^2/(2r^4).

T denotes the equation trace, not a stress tensor. The formula is the
full-action trace, without the stationary bookkeeping parameter. At
the finite comparison g1 the residual is not zero. Statements below
that require calE=0 are conditional on an exact stationary solution;
this note does not substitute g1 for that solution.

First take a smooth spherical conformal change g_s=exp(2s)g, with s
zero on an open neighborhood of the initial Cauchy surface and near
the two exterior endpoints. Start in the selected HH state. Transport
the Dirac field by Psi_s=exp(-3s/2)Psi, keeping A fixed. Conformal
covariance transports its causal Green functions and preserves the
initial density matrix and boundary projector. This is the evolved
initial HH state, not an instantaneous thermal-state assumption.
The massless-field identity on the two in-in contour branches gives
the local Weyl cocycle, with no change of the past state. In particular

    delta_s T_Dirac,ab=-2s T_Dirac,ab+A_ab[s;g],              (2)

where A_ab is the metric variation of the linear anomaly cocycle.
It is local and of differential order at most four in s. The
Weyl-invariant functional contributes the base-stress rescaling in
(2); it has not been discarded. Maxwell stress has the same rescaling
without A_ab. The gauge-anomaly cocycle varies the current by a constant
multiple of nabla_a(s F^ab). For the spherical magnetic background,
F has only angular components, nabla_a F^ab=0, and s is independent
of angle. This current variation vanishes. Thus this conformal column
does not require an electric or magnetic perturbation to repair Maxwell.

Independent initial-state changes add their own causal stress and
current response. At fixed metric their bulk stress trace is zero by
the same anomaly identity. They cannot cancel a nonzero variation of
(1). They remain relevant in the other tensor equations. A general
shape perturbation requires the full retarded stress and current
response; (2) gives only its conformal column.

[Anderson, Molina-Paris and Mottola](https://arxiv.org/html/0907.0823v2)
separate causal metric response, local contacts and independent state
variations in sections I-III, particularly (3), (27)-(33). Their FRW
calculation is not a charged HH transfer theorem. Their removal of
Planck-scale branches is not adopted here. Source scope and document
identity are retained in the accompanying source-review receipt.

## 2. Exact conformal column, including residual contacts

Let Lh=DT_g[h], with h a covariant metric variation and fixed A.
Let L* be its LOCAL formal transpose with respect to the spacetime
volume and h^ab pairing, using compact variations. Then the complete
causal conformal column in the state prescription above is

    D_ret calE_ab[2s g]
       =-2(L*s)_ab-s T g_ab-2s calE_ab.                      (3)

To derive the contact factors, suppress a common action normalization
and write delta I=int sqrt(-g) calE^ab h_ab. Its conformal derivative
is -2 int sqrt(-g) s T. Differentiate this in h. The Weyl vector 2s g
also varies, giving delta I[2s h]. On the other side the volume varies
by 4s and raising the two equation indices contributes -4s. The result
is (3). For the causal quantum part this argument applies to the LOCAL
cocycle in (2), retaining the base-stress rescaling. It is not an
assertion that a general retarded kernel is self-adjoint in time.
There are no boundary integrations in this compact-support statement.

An explicit local representation useful for implementing (3) is

    Lh=(1-11epsilon R/36+epsilon Box) DR[h]
        +epsilon[7 D(C^2)[h]/12+11 D(J^2)[h]/6
                  -20 DM[h]+K_R[h]],                       (4)
    DR[h]=nabla_a nabla_b h^ab-Box H-Ric_ab h^ab,
    H=g^ab h_ab,
    K_R[h]=-h^ab nabla_a nabla_b R
           -(nabla_a h^ab-nabla^b H/2)nabla_b R.

Here K_R is the variation of the Box operator acting on the background
R, not Box DR. Each D in (4) is an ordinary metric derivative at fixed
A. In particular DM[h]=-2F_ac F_b{}^c h^ab. Its transpose form is

    L*s=DR*[(1-11epsilon R/36)s+epsilon Box s]
          +epsilon[7 D(C^2)*s/12+11 D(J^2)*s/6
                    -20 DM*s+K_R*s],                       (5)
    (DR*f)_ab=nabla_a nabla_b f-g_ab Box f-Ric_ab f,
    (K_R*s)_ab=nabla_(a s nabla_b)R
                 -g_ab[(nabla s).(nabla R)+s Box R]/2.

No unknown nonlocal kernel occurs in (5). The remaining two transposes
can equivalently be evaluated as the metric derivatives of the weighted
local functionals int sqrt(-g)s C^2 and int sqrt(-g)s J^2, subtracting
their volume terms s C^2 g_ab/2 and s J^2 g_ab/2. This specifies the
tensor without a curvature-sign convention hidden in a named tensor.

At calE=0, T=0, (3) reduces to -2L*s. At exact AdS, C=J=M=0,
R=-12lambda, so L=(1+11epsilon lambda/3+epsilon Box)DR and

    D_ret calE_ab[2s g]=-2epsilon Q_ab(Box+E)s,
    Q_ab f=nabla_a nabla_b f-g_ab Box f+3lambda g_ab f.       (6)

For a charged non-Einstein background, (4)-(5) retain additional
curvature tensors and their derivatives. Off shell the two contacts
in (3) cannot be dropped. Already the Einstein limit tests their signs:
with T=R+12lambda0 and calE=G-3lambda0 g, (3) gives exactly
-2 Hess s+2g Box s-6lambda0 s g.

## 3. Gauge-invariant trace and failure of the direct lift

In a monotone exterior use the background areal radius r. Under a
spherical infinitesimal diffeomorphism, delta R changes by xi^r R'
and delta r changes by xi^r. Therefore

    Psi=delta R-R' delta r                                 (7)

is gauge invariant. The bare delta R is invariant only if R'=0.
Equivalently choose h_A=h-Lie_xi g with xi^r=delta r, so delta r_A=0.
The residual time gauge does not affect any static background scalar.
On an exact background, DT[Lie_xi g]=Lie_xi T=0, hence (4) becomes

    DT[h]/epsilon=(Box+E-11Z/36)Psi+K_R[h_A]
                     +7 delta_A(C^2)/12+11 delta_A(J^2)/6,
    Z=R+12lambda.                                         (8)

The magnetic term vanishes in (8) because fixed q and delta r_A=0
give delta_A M=0. The other terms depend on the remaining metric
variables. Equation (8) is not a closed Klein--Gordon equation for Psi.
It does not require the numerical quantum stress to be known to expose
the missing coupling. For h=2s g one has delta r=rs and

    Psi=-6 Box s-(2R+rR')s.                                (9)

The exact conformal trace, using DT from (1), is

    DT[2s g]=-6epsilon Box^2 s-(6+epsilon R/6)Box s
       -(22epsilon/3)J^ab nabla_a nabla_b s
       -2epsilon(nabla R).(nabla s)
       +[-2R-4epsilon B_an]s,                             (10)
    B_an=7C^2/12+11J^2/6-11R^2/72-20M+Box R.

On shell the final coefficient is 2R+48lambda0. For completeness,
the curvature identities behind (10) are

    delta_s R=-6Box s-2Rs, delta_s C^2=-4s C^2,
    delta_s M=-4s M,
    delta_s Euler=-4s Euler+8J:Hess s-2R Box s,
    delta_s(Box R)=-6Box^2 s-2R Box s
                     -2 grad R.grad s-4s Box R.

There is a local obstruction to importing the SAME equation
(Box+E)s=0 as a conformal metric lift. At a point where s=0 and
nabla s=0, a solution of this equation has Box s=Box^2 s=0.
Equation (10) then gives

    DT[2s g]=-(22epsilon/3)J:Hess s.                        (11)

For spherical fields on a static spherical background, the angular
Hessian at this point is zero. In an orthonormal frame choose
H_tt=H_rr=H0, with the other Hessian entries zero. This is tracefree
and is compatible with local spherical KG Cauchy jets. If J^t_t=j_t
and J^r_r=j_r, then J:H=(j_r-j_t)H0. When j_r != j_t, (11) is
nonzero. Even when j_r=j_t, a different jet with s or first derivative
nonzero can contribute in (10); (11) alone does not test that case.

For unrestricted perturbations, any nonzero J permits a tracefree
Hessian with nonzero contraction. This is a pointwise algebraic
obstruction, not a constructed mode or an instability theorem. In
particular, the classical RN tensor has j_r=j_t: it would be wrong
to claim that the special spherical zero-first-jet test proves failure
on that metric. An independent explicit RN test below supplies that
case without calling classical RN an on-shell semiclassical solution.

For classical RN-AdS, write F=1+lambda r^2-2m/r+Q2/r^2,
R=-12lambda, J^a_b=diag(-Q2/r^4,-Q2/r^4,Q2/r^4,Q2/r^4).
For any spherical KG solution,

    J:Hess s=(Q2/r^4)[E s+4F s'/r].                       (12)

At a point with s=0, s' !=0, (10) on this classical reference is
-88epsilon Q2 F s'/(3r^5), nonzero for Q2,F,s' !=0.
All other terms there vanish by R'=0 and the KG equation. Thus the
vacuum conformal KG lift fails even the trace test on the charged
classical reference. The value does not depend on dropping an
off-shell multiple of s. It is a reference diagnostic, not the
charged semiclassical background's spectrum. Independent state
variations cannot repair its nonzero trace.

## 4. Constraints propagate in time; static recovery is insufficient

On a stationary exact solution write

    ds^2=-N(r)^2dt^2+A(r)^2dr^2+r^2dOmega^2,
    e_t=delta calE^t_t, e_r=delta calE^r_r,
    e_o=delta calE^theta_theta, j=delta calE^r_t,
    delta calE^t_r=-A^2 j/N^2.

When the joint Maxwell/current equations hold, the linearized Ward
identity is nabla_a delta calE^a_b=0. Direct differentiation gives

    partial_t e_t+[partial_r+(log(N A r^2))']j=0,
    -(A^2/N^2)partial_t j+e_r'
                 +(N'/N)(e_r-e_t)+2(e_r-e_o)/r=0.          (13)

These are on-shell identities; variations of the connection multiply
the background residual off shell. They also explain why the static
two-row recovery cannot simply be continued to time-dependent fields.
If the trace and orbit trace vanish, then e_o=0 and e_t=-e_r. The
remaining equations are

    -partial_t e_r+j'+(N'/N+A'/A+2/r)j=0,
    -(A^2/N^2)partial_t j+e_r'+2(N'/N+1/r)e_r=0.            (14)

They are a homogeneous hyperbolic constraint propagation system with
characteristic speeds +/-N/A on a regular exterior patch. They do not
set e_r and j to zero. Their initial data and incoming boundary data
must vanish to recover all equations. For example in a flat exterior
patch N=A=1,

    e_r=j=r^-2 f(t+r), e_t=-e_r, e_o=0                   (15)

solves both Ward equations for any smooth nonzero f. This is a tensor
counterexample to recovery from the two rows alone, not a physical
metric solution. In the stationary case j=0, (14) instead yields
(N^2 r^2 e_r)'=0, and regularity at a nonextremal bolt kills its single
constant, as in the previous static argument. Lorentzian initial
constraint data cannot be replaced by that static constant argument.

## 5. Consequence and verification contract

The next transfer must use the full retarded shape response with (3),
(8) and the initial constraints (13). A static Euclidean shape gap
does not estimate that causal response. No controlled remainder in
the normalized scalar tail of the vacuum preflight has been proved.
Conformal locality alone supplies neither a constraint-satisfying
scalar evolution nor a physically justified exclusion of its branch.
No order reduction, finite-counterterm change or independent state
adjustment is used to remove the branch. S168.3-S168.4 stay open.

Before execution freeze this note, checker and source review. Reconcile
the previous vacuum receipt and its transitive input hashes. Check
the contact factors, local transpose, gauge cancellation, trace
decomposition, RN contraction and KG counterjet exactly. Independently
form the Lorentzian Christoffel symbols and covariant divergence of a
general symmetric spherical equation tensor to test (13)-(15).

For the conformal column, independently vary the Einstein equation
of a nonconstant static warped metric at two rational radii, retaining
the lapse variation, and compare (3) in the epsilon=0 limit. For the
anomaly transpose contact, independently differentiate the one-dimensional
radial density action int sqrt(-g) s Box R with arbitrary lapse,
time and radius variations, and compare its Euler derivatives with
the transpose of Box DR+K_R plus the volume contact. These exact
local fixtures are neither exterior solutions nor spectral data.
Detect omission of the residual contacts and of K_R. Retain all runs
and exact failed sources. Continuum covariance, the causal state
identification and local KG solvability remain ordinary arguments;
no Lean campaign is selected and no simulation is claimed.

Remaining to-do list: full constrained retarded shape/scalar transfer
or justified physical exclusion; S167 signed coefficient, cold/higher
response, nonlinear solve and compatibility; actual exterior/state,
charge/return, opening, stability, transport, preparation, creation
and full-result audit.
