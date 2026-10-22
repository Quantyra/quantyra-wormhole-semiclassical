# Charged R-squared reference and its anomaly-complete bulk trace block

2026-10-22, E003/S169. Ordinary charged-reference theorem, not a
semiclassical exterior or a coupled physical inverse. S169.2 stays open.

Keep the [selected additional Lorentzian coefficient](finite-r2-scalar-window.md)
c=10^22, all other model data and the complete Dirac determinant. The
previous goal turn made progress: commits 7cb6831 and planning 3f3539e
selected and quantitatively tested this candidate. This calculation
addresses its charged equations before extending the old estimate chain.

Two useful results hold. The LOCAL Einstein-Maxwell-R^2 theory has an
exact charged reference and an on-shell scalar tensor lift. On that
reference the FULL BULK Dirac trace anomaly, including the volume-density
contact, gives a coercive conformal trace block. Its comparison inverse
has energy-norm bound 1, independent of cutoff and positive temperature.
The reference is off shell for the complete Dirac theory. Its quantum
boundary terms and remaining shape equation are not solved by this result.

## 1. Complete bulk equations and the actual local charged reference

Write the electromagnetic two-form as F_em, I=F_em^2 and
T_EM,ab=g_R^-2(F_em,ac F_em,b^c-g_ab I/4). Introduce varphi without
changing the selected action:

    S_local=int sqrt(-g)[(R-2Lambda)/(16piG)
                      +c(2varphi R-varphi^2)-I/(4g_R^2)],
    F_s=1+32piG c varphi,
    F_s G_ab+Lambda g_ab+8piG c varphi^2 g_ab
       +32piG c(g_ab Box-nabla_a nabla_b)varphi
                         =8piG(T_EM,ab+T_D,ab),
    R=varphi.                                             (1)

Here T_D is the variation of the UNCHANGED renormalized full Dirac
functional, with its paired domain, scale and state; it is not supplied
by its trace. Define j_D=-delta Gamma_D/(sqrt(-g) delta A), so the other
equation is nabla_a(F_em^ab/g_R^2)=j_D^b. The metric/gauge Ward identity
must use the same functional and boundary realization. The divergence
of the local gravitational left side is
16piG c(varphi-R)nabla_b varphi. Thus its conservation follows with the
auxiliary equation and matter/current Ward identities, not an independent
choice of quantum pressures.

Put lambda0=-Lambda/3=1/(G k^10), k=2000, q=k^4, g_R=k^-4.
The LOCAL reference (T_D=j_D=0 only for defining this reference) is

    R0=-12lambda0, F0=1-384piG c lambda0,
    Q_c^2=piG q^2/(g_R^2 F0),
    f(r)=1-2M/r+Q_c^2/r^2+lambda0 r^2,
    ds^2=-f dt^2+dr^2/f+r^2 dOmega^2,
    F_em=(q/2)sin(theta)dtheta wedge dphi.                  (2)

For constant R0, (1) becomes F0(Ric-Lambda g)=8piG T_EM;
direct substitution proves (2). Q_c is a geometric charge coefficient,
not a redefinition of the fixed physical flux or Maxwell coupling.
The established f(R)-Maxwell solution is the comparator:
[Moon--Myung--Son, section 2, (2.2)-(2.13)](https://arxiv.org/pdf/1101.1153v2).
Its normalization is converted explicitly above; no novel black-hole
solution or physical quantum state is claimed.

Let C_c=sqrt(1+12lambda0 Q_c^2), B_c^2=(C_c-1)/(6lambda0).
The extremal radius is B_c. Every exterior with outer radius r_+>=B_c
obeys Q_c^2/B_c^4=3lambda0+B_c^-2. For r_+>B_c the smooth Euclidean
period is set by

    T(r_+)=[1/r_+-Q_c^2/r_+^3+3lambda0 r_+]/(4pi),
    T'(r_+)>= [3lambda0-1/(12Q_c^2)]/(4pi)>0.             (3)

The derivative bound minimizes 3Q_c^2 x^2-x over x=1/r_+^2.
Since T(B_c)=0 and T tends to infinity, each prescribed positive
temperature has exactly one local reference radius. This avoids using
temperature as a free compatibility parameter. It does not transfer
the old prepared gap to the new quantum geometry.

On every such exterior, define t=Q_c^2/r^4, W=M/r^3-Q_c^2/r^4.
The Euclidean traceless Ricci eigenvalues are (-t,-t,t,t), and

    J^2=4t^2, C^2=48W^2, I=q^2/(2r^4),
    0<=t<=3lambda0+B_c^-2,
    |W|<=5lambda0+2B_c^-2.                              (4)

For the last inequality use M/r_+^3=(r_+^-2+Q_c^2/r_+^4+lambda0)/2
and r>=r_+. These bounds hold for the entire exterior, not sampled radii.

## 2. Local charged scalar lift, with constraints retained

For the local theory the trace is R-4Lambda-96piG c Box R=0.
Let psi solve (Box-m0^2)psi=0, m0^2=1/(96piG c), on (2). Set

    h_ab=2s g_ab, s=-16piG c psi/F0, delta A=0.           (5)

The conformal curvature identity gives delta R=-6(Box-4lambda0)s=psi.
Moreover delta F_s=-2F0 s. In the metric variation of (1), the
F0 delta G terms cancel 32piG c(g Box-Hess)psi. Maxwell conformal
covariance gives delta T_EM=-2s T_EM. The remaining algebraic tensor
coefficient reduces to 4Lambda-R0=0 after using the background equation.
Thus ALL local tensor components and the Maxwell equation are solved.
This is an on-shell charged lift, not the failed off-shell seed lift.

Equivalently the Einstein-frame metric variation
delta(F_s g)=F0 h+delta F_s g vanishes. The local scalar and metric
blocks therefore decouple to first order at this reference in those
variables. The original metric observables still require the inverse
conformal map. In the full theory the Dirac anomaly couples the scalar
back to the shape; conformal invariance of the massless classical field
does not remove this anomaly or the shape determinant.

## 3. Exact density contact and the full bulk conformal trace block

Set epsilon=G/(60pi), D=96piG c-epsilon>0 and
H=7C^2/12+11J^2/6-20I. The exact bulk trace on a metric is

    T=R+12lambda0+epsilon(H-11R^2/72)-D Box R.           (6)

On (2), T0=epsilon(H0-22lambda0^2) generally is NONZERO.
For a genuine conformal coordinate g_s=exp(2s)g0 the fixed-background
density row is exp(4s)T[g_s]. Its derivative is
delta T+4s T0, not delta T alone. Constant R0 and div J=0 give

    A s=6D Delta^2 s-div(K grad s)+24lambda0 s,
    K=(6+24D lambda0+22epsilon lambda0)g
                                      +(22epsilon/3)J. (7)

Delta is the nonpositive Euclidean Laplacian. This identity includes
the complete bulk anomaly of the full angular tower; it is not a
lowest-level truncation or an estimate on a spatially analytic family.
For verification without rearranging a preexisting Hessian:

    sqrt(g_s)R_s=sqrt(g0)exp(2s)(R0-6Delta s-6|grad s|^2),
    delta[sqrt(g_s) C_s^2]=delta[sqrt(g_s) I_s]=0,
    delta[sqrt(g_s) J_s^2]=-4sqrt(g0) J:Hess s,
    delta[sqrt(g_s) R_s^2]=-12sqrt(g0) R0 Delta s,
    delta[sqrt(g_s) Box_s R_s]
                  =sqrt(g0)(-2R0 Delta s-6Delta^2 s).

These identities derive (7) directly. The apparent gauge and Weyl
zeroth-order potentials cancel with the density contact. No physical
source has been set to zero: T0 is still the forcing. The cancellation
cannot be used for a nonconformal variation or for the shape-shape
Hessian. With the repository Euclidean action convention the conformal
action Hessian is -A/(8piG); positivity of A is an elliptic inverse
statement, not positivity of the full gravitational energy.

## 4. Realized bulk comparison inverse and numerical margins

First work on the smooth Euclidean bolt exterior truncated at r=L.
Use the stated comparison domain s=0, Delta s=0 at L and regularity
at the bolt. There is no physical wall at the bolt. For the LOCAL
metric/varphi Dirichlet data these conditions follow from s=0 and
delta R=0. The quantum boundary completion is discussed below.

The form domain is H^2 intersect H_0^1 on this compact smooth manifold,
with no independent condition on the normal first derivative. Its form is

    a(s,v)=int[6D Delta s Delta v+K(grad s,grad v)
                                            +24lambda0 s v],
    ||s||_E^2=int[6D(Delta s)^2+6|grad s|^2+24lambda0 s^2]. (8)

The natural second boundary condition is Delta s=0. Imposing both
s=0 and its normal derivative zero would be a different realization.
Using (4),

    K >= [6+24D lambda0-(22epsilon/3)B_c^-2]g.           (9)

Rational benchmark bounds give

    1-1.179e-8 < F0 < 1,
    B_c^-2/lambda0 < 1.223e-10,
    7.06e-8 < 24D lambda0-(22epsilon/3)B_c^-2 < 7.08e-8.
                                                                    (10)

Consequently a(s,s)>=||s||_E^2. On each finite cut, standard elliptic
norm equivalence makes the form domain complete; boundedness follows
from the bounded coefficients in (4). The Riesz/Lax--Milgram argument
gives one weak solution for each E-dual source and

    ||A^-1 f||_E <= ||f||_(E*).                         (11)

This is a bound on an actual infinite-dimensional comparison operator,
not on a sampled matrix. It is uniform in L and every positive T.
For the complete end take the E-completion of smooth compactly supported
functions regular at the bolt. The same bounded coercive form gives its
inverse for E-dual sources. No cut-dependent Poincare constant is used.
Elliptic regularity gives the interior strong equation. At finite cuts
the variational natural condition gives the displayed Navier condition.

Equation (11) does NOT apply to the raw constant asymptotic forcing T0:
T0 tends to -22epsilon lambda0^2, which is not in this global dual.
For example unit-amplitude tests on large annuli have volume-growing
pairing while E grows only as the square root of volume. The outer
geometry must first be centered at the exact quantum vacuum lambda,
where g=G lambda solves g+11g^2/(360pi)=k^-10. Its curvature shift is
12(lambda0-lambda), not zero. The original reference also needs recentering
at the new full-Dirac product. A local black-hole reference does not
supply the old finite geometry, prescribed state or admissible residual.

## 5. What is still required for the coupled physical inverse

The local boundary action is
epsilon_n int sqrt(|h|)[1/(8piG)+4c varphi]K with fixed h,varphi;
see [Dyer--Hinterbichler, (1.4)](https://arxiv.org/pdf/0809.4033v2).
This is only the EH plus added R^2 sector. The retained
[spherical anomaly boundary calculation](spherical-anomaly-boundary.md)
exhibits independent local boundary momenta, and the
[paired Dirac heat calculation](curved-fermion-boundary.md) contains
normal-smearing terms. Fixing s and delta R alone does not establish
that all these terms vanish. Navier is therefore a PROVED BULK COMPARISON
REALIZATION, not an asserted physical determinant boundary condition.
No finite-cut boundary functional or cutoff limit is silently changed.

The complete mixed conformal row at this reference is explicit. For an
arbitrary covariant metric variation k, with gauge field fixed, it is

    Bk=[1+11epsilon lambda0/3-D Delta]delta R[k]
         +epsilon[7delta C^2[k]/12+11delta J^2[k]/6-20delta I[k]]
         +(tr k)T0/2.                                  (12)

For TT k, delta R[k]=-J:k, but neither the remaining mixed terms nor
the shape response may be discarded. In a fixed conformal frame the
Weyl-invariant determinant has no conformal derivative; its full nonlocal
shape-shape derivative remains. This locality is an actual cancellation
that can be used in computing mixed estimates, not a norm rescaling.

After a justified physical domain and recentering, write the Jacobian
blocks as [A+E, B; C, L]. With (8) FIXED and compatible spaces for the
shape/current/constraint variables, sufficient numerical targets are

    ||E||_(E->E*) <= 1/4,
    ||L^-1 C||_(E->X) ||B||_(X->E*) <= 3/8.             (13)

Then ||(A+E)^-1||<=4/3 and the Schur loop
L^-1 C(A+E)^-1 B has norm <=1/2. L must be invertible on its actual
gauge-fixed domain, including constraints and physical boundary data.
The scalar response obeys ||s||_E<=(4/3)(||f||_E*+||B|| ||k||_X);
the shape response obeys
||k||_X<=2[||L^-1 g||_X+(4/3)||L^-1 C|| ||f||_E*].
These are requirements, not numerical claims about the uncomputed
quantum shape operator. Boundary perturbations must be bounded in the
same forms or separately realized. Unbounded trace terms cannot simply
be included in E. Source and mixed norms in (13) are UNKNOWN, not known
upper bounds already meeting the targets. Rescaling E or X cancels in
the Schur loop and does not improve feasibility.

## 6. Decision and verification

Continue the c=10^22 candidate. The charged bulk trace block passes its
quantitative inverse test with the full local anomaly, and the local
charged scalar has a constraint-respecting tensor lift. Next compute
the renormalized paired-domain boundary quadratic form and the remaining
coupled shape response, preferably in the auxiliary Einstein-frame
variables where the local scalar/metric mixing vanishes at (2). Retain
the anomaly mixed terms and original metric observables. Recenter the
outer vacuum BEFORE applying (11) to a residual. Do not extend the old
uniform source-times-inverse chain or import the old signed projection.

The checker independently derives RN curvature from sectional curvatures,
the local tensor cancellation, density variation, global curvature and
coercivity bounds, and temperature monotonicity. Rational arithmetic
certifies (10). At 60/90 digits it differentiates the unexpanded
conformal bulk density on charged annuli and compares with (7).
It also checks radial form/integration-by-parts consistency with the
actual charged coefficients. These are supporting checks, not quantum
state samples, formal continuum verification, or a coupled inverse.
Frozen source snapshots, input hashes, logs and failures are retained.

Remaining to-do list: physical quantum boundary form and coupled shape
inverse; recentered finite geometry, full residual and nonlinear map;
controlled solution and exact compatibility; matching, validity, full
stability, charge/return, opening, transport, preparation, separate
creation and the full-result audit.
