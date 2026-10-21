# Coupled fast tails and the quadratic boundary obstruction

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The independent second-order optical shape cannot repair the preceding
trace/action mismatch in the stated expansion. The leading trace and
boundary action depend on the same combination of the two metric tails.
The canonical determinant contributes a finite second-order coefficient,
but no persistent cutoff oscillation. Every trace-compatible pair therefore
retains the previously derived nonzero oscillation.

This is an obstruction to a convergent quadratic coefficient of the
unchanged finite-cut action on this family. It is not an exclusion of
every possible fast metric or boundary theory, a Lorentzian instability,
or a wormhole no-go theorem. Full nonlinear existence remains unresolved.

## 1. Scope and the necessary trace equation

Use the unchanged action, paired domain and fixed background cut of the
[finite-cut calculation](fast-trace-finite-cut-action.md). Write
z=log(r/B), epsilon=G/(60pi), t=epsilon lambda>0,
k^2=1/t+17/12 and S0=lambda B^3/G. In the optical chart put

    w=Phi-Phi_seed=eta exp(-3z/2) f(z)H(tau)
                                      +eta^2 exp(-3z)h(z,tau),
    v=V-V_seed=eta^2 exp(-3z)j(z,tau),
    f=a cos(kz+phi), a>0.                                (1)

H is a fixed nonzero thermal cosine. The first shape variation is zero.
Second-order h,j have bounded differentiated radial coefficients and
fixed finite thermal frequencies. Smaller differentiated remainders do
not change the displayed leading coefficients. No theorem that every
exact solution, or the full first-order auxiliary shape response, has
this expansion is asserted. In particular, arbitrary first-order shape
tails are not silently included.

Time integration only sees the second-order mean. Define hbar=<h>/<H^2>
and jbar=<j>/<H^2>; the former special choice h(z,tau)=hbar(z)H^2 is
not required. Fixed-frequency derivatives and charged corrections gain
two radial powers in the equation at this order. The
[retained nonlinear trace calculation](fast-nonlinear-trace-boundary-compatibility.md)
therefore gives the necessary mean equation

    L_h hbar+L_j jbar+N_t=0,
    C(D)=3tD^2-9tD+11t+3, D=d/dz,
    L_h=-(D-4)(D+1)C(D),
    L_j=-(D-6)(D+1)C(D)/3,
    N_t=(471t/4+289/4+6/t)f^2+(35t/3-5)ff'
                                         -(21t+3)(f')^2. (2)

Both zero and doubled temporal frequencies remain necessary in the full
equations. Taking their mean here extracts a necessary condition; it
does not discard the other equation or invoke a positive-frequency
inverse to solve the static compatibility condition.

## 2. Complete local shape contribution at the fixed cut

At first order in an independent shape variation v, the background
areal chart stays fixed. The physical areal variables are exactly

    F_new=exp(-2v)F_seed, alpha_new=2v.

Differentiate the original Einstein and converted Weyl bulk densities,
the GHY term, both Einstein counterterms, the actual Weyl boundary
polynomial and the existing loop counterterms. In units S0 exp(3z),
their leading first-order static bulk and boundary coefficients are

    L_v=(1/2+2t)(6v-v'-v''),
    B_v=t v''+(1/2+9t/7)v'-(1+76t/7)v.                  (3)

The classical magnetic density -exp(alpha)Q2/(2G r^2), with
Q2=pi G q^2/g_c^2, is retained in the finite-radius fixtures. Its shape
variation is subleading here. This physical magnetic coefficient is
distinct from the comparison metric's charge coefficient.

The bulk primitive is exp(3z)(1/2+2t)(2v-v'). Combining it with
the boundary value leaves exp(3z)(t/7)(7v''-5v'-48v). The intrinsic
curvature terms are retained until their lower-power limit; none is
deleted from the finite-radius calculation. Varying the proper radial
derivative in the Weyl boundary polynomial is essential for (3).

For v=exp(-3z)jbar, its persistent contribution is

    M_j jbar=(t/7)(7D^2-47D+30)jbar
                          =(t/7)(7D-5)(D-6)jbar.       (4)

The conformal contribution already derived is

    M_h hbar=(3t/7)(7D^2-33D+20)hbar
                         =(3t/7)(7D-5)(D-4)hbar.

Thus the complete LOCAL quadratic cutoff term is

    S0 <H^2>[Q_t+M_h hbar+M_j jbar],
    Q_t=-(3001t+546)f^2/504-(71t+21)ff'/14
                                      -(167t/42)(f')^2. (5)

This expression must still be checked against the canonical determinant.
Treating that determinant as unchanged when v varies would be incorrect.

## 3. Canonical first variation on the critical cubic directions

Only the FIRST canonical variation enters [eta^2]Gamma_hat[V_seed+v],
because v starts at eta^2 and Gamma_hat is independent of w. Its value
is not set to zero. We prove that it has a finite cutoff limit, even
for the oscillatory cubic j in (1). This is a fixed-seed first-variation
statement, not C2 dependence on a whole family of rough baselines.

Use background canonical proper distance Y at the limiting boundary.
The smooth seed coefficient has

    A(Y)=1+Y^2/(2B^2)+c3 Y^3+O(Y^4 log^p Y).

The shape direction u(Y)=exp(-3z)j is O(Y^3), with
u'=O(Y^2), u''=O(Y), u'''=O(1). It belongs to C^(2,1) on a closed
collar, and all of its first three boundary jets u,u',u'' vanish.
Its third derivative need not have a boundary limit. Hence the earlier
[stationary C1 phase theorem](stationary-fast-phase-variation.md), which
fixes a convergent cubic coefficient, cannot simply be applied to the
perturbed metric family.

Instead use the actual fixed-seed orbit-Weyl variation identity from the
[canonical boundary charge proof](stationary-same-source-boundary-charge.md):

    delta zeta_L(s)=2s[Tr(u L^-s)-Tr(u m_S^2 L^(-s-1))]. (6)

It retains both partner squares, their induced domains, the complete
angular sum and the common four-dimensional subtraction. The massless
angular sector has zero mass weight, with its anomaly still included.

At a cut Y=d, pull back the inward collar by Y=d+x. The normalized
background coefficients and their jets converge smoothly. The
[balanced reflection identity](stationary-massive-boundary-layer.md)
annihilates the flat reflected spin trace. Consequently there is no
weight-zero boundary layer. For a boundary heat profile of weight n,
the orbit diagonal scales as heat_time^(n/2-1). Multiplying by the
angular mass weight O(heat_time^-2) and integrating heat_time gives a
meromorphic x^(n-4) profile. The remaining low weights n=1,2,3 therefore
produce only finite parts x^-3,x^-2,x^-1, possible logarithmic partners,
and contacts of order at most two. Their coefficients and fixed-scale
finite parts remain bounded as d tends to zero. Odd-linear-jet terms
tend to zero; the cubic profile has its retained finite limit.

The full mixed heat coefficient [Vassilevich, (5.33)](https://arxiv.org/html/hep-th/0306138v3)
has a possible third smearing derivative proportional to tr chi. It
vanishes identically in the balanced paired domain. The remaining local
boundary anomaly uses u,u',u'' with bounded coefficient jets, and has
zero limit for the present direction. At the geodesic boundary only
the previously derived u(0) contact remains, and that also vanishes.

All higher-weight mass-Green remainders are uniformly L1 after the
common subtraction, by the fixed-seed boundary-layer estimate. The
smooth bulk restoration, thermal/state and remote-attachment pieces
are in the same integrable class; the retained positive gap controls
large heat time. This estimate concerns the baseline heat kernel and
its smearing, so it does not require five bounded derivatives of u.

For completeness, each FP x^-p, p<=3, with fixed finite-part scale is
continuous on C^(2,alpha), 0<alpha<1: subtract the quadratic Taylor
polynomial. The remainder is O(x^(2+alpha)), integrable against x^-3.
Its contacts require at most the second jet. Logarithmic partners have
the same continuity. The translated tests u_d(x)=u(d+x) converge to u
in C^(2,alpha), since they converge in C2 and have a common Lipschitz
bound on the second derivative. Smooth extension at the far collar end
is fixed throughout. Pairing the finite profiles, local anomaly and
uniformly L1 remainder therefore gives

    lim_(d->0) delta Gamma_hat,d[u_d]
                         = delta Gamma_hat,0[u], finite. (7)

This extends the finite-cut first variation to these fixed nonsmooth
tests; it is not an inference from pointwise trace asymptotics alone.
No collar width is tied to d and no log(d) renormalization scale is
inserted. The limit can depend on the interior and state. It contributes
a finite constant to (5), not a persistent radial cutoff oscillation.
Nonzero temporal components of u have zero linear pairing with the
stationary seed, so only the mean used in (2) is needed here.

An explicit check of the continuity step is available on [0,1]. For
u(Y)=Re Y^nu, nu=3+i omega, and finite-part length one, let T_p=FP Y^-p.
Then T_p[u]=(4-p)/[(4-p)^2+omega^2], p=1,2,3. For 0<d<=1/2,

    |T_p[u(d+.)]-T_p[u]| <= C_omega sqrt(d),
    C_omega=(16/15)|nu(nu-1)(nu-2)|
                              +1+|nu|+|nu(nu-1)|.       (7a)

Indeed the difference of second derivatives has supremum at most
|nu(nu-1)(nu-2)|d and Lipschitz constant at most twice that coefficient.
Interpolation gives a half-Holder seminorm at most twice the coefficient
times sqrt(d). The quadratic Taylor remainder integrates to at most
(16/15)|nu(nu-1)(nu-2)|sqrt(d) against the worst Y^-3 kernel. The
remaining finite-part polynomial terms obey the other terms in (7a).
This proof does not assume that u''' has a limit.

## 4. Both equations see the same combination

Define the single radial combination

    X=3(D-4)hbar+(D-6)jbar.

Equations (2),(4),(5) become

    (D+1)C(D) X=3N_t,
    quadratic cutoff term=S0 <H^2>[Q_t+(t/7)(7D-5)X]. (8)

Equivalently the determinant of the two commuting operator rows is
L_h M_j-L_j M_h=0 identically. Treating hbar and jbar as two freely
adjustable cancellation variables would miss this exact degeneracy.
The shape contribution must be derived from the unchanged boundary
action to obtain the common combination.

For bounded hbar,jbar and their stated jets, X is bounded. The
homogeneous roots of (D+1)C(D) are -1 and 3/2 +/- i k. Thus its bounded
normal solution is the unique periodic forced part plus an exp(-z)
term. If the leading equation has an o(1) error, the bounded solution
and its first derivative differ from that periodic part by o(1):
variation of constants uses an integrable exponentially decaying kernel
after excluding the two growing homogeneous components. Splitting the
forcing at a large fixed z proves the o(1) statement. This supplies the
needed endpoint control without assuming convergence of hbar or jbar
separately.

Choose the particular solution with jbar=0. The preceding trace result
already evaluates the forced X and gives, for unit radial amplitude,

    Q_t+(t/7)(7D-5)X
       =constant+R_c cos(2kz+2phi)+R_s sin(2kz+2phi)+o(1),
    R_c=-(78815t^3+74391t^2+16137t+81)
                              /[168(5t+3)(11t+3)]<-3/56,
    R_s=-k(6655t^3+111t^2-3015t-567)
                                  /[84(5t+3)(11t+3)].   (9)

Every bounded trace-compatible coupled pair has the same persistent
oscillation. Phase rotation cannot annihilate a nonzero coefficient
vector, and its amplitude exceeds (3/56) S0 <H^2> a^2. The finite
canonical first variation in (7) cannot cancel it. Changes in the
individual tails preserving X can change a finite action constant and
other equations, but they do not change (9).

Consequently the prescribed finite-cut quadratic coefficient has no
limit on a family of the form (1) with nonzero first-order fast amplitude
and trace-compatible second-order coefficients. In particular this
family does not supply a C2 variational completion obtained with a
convergent quadratic finite-cut coefficient. We do not infer convergence
or differentiability of a finite-amplitude rough determinant family from
(7), exchange unproved limits, or declare all possible physical fast
directions excluded.

## 5. Verification and physical consequence

The checker differentiates the original finite-radius shape action,
including the varied radial metric and literal boundary terms; checks
the exact row degeneracy and the shared forced oscillation; and tests
the shape coefficient on refined charged annuli at two precisions.
Independent complex-amplitude extraction checks its linear Taylor
coefficient against the unexpanded local action.

Finite-part translation controls use u(Y)=Y^3 cos(omega log Y), which
has no limiting third derivative. They compare the continuous FP
Y^-p actions, p=1,2,3, against the exact d=0 values and the proved
C^(2,alpha) continuity bound. These are tests of the distributional
step; they are not numerical evaluation of the full determinant.
The finite-cut canonical statement remains an ordinary proof using
the retained actual-domain heat profiles. No Lean run is selected.

The outcome closes this coupled second-order repair within (1).
The remaining physical decision is to determine whether the complete
first-order shape response and a justified gravitational boundary phase
admit the fast branch under the unchanged action. Stronger exclusion
must be matched to actual source and initial-data compatibility. An
additional boundary functional, boundary degrees of freedom or changed
cut prescription would require its own explicit physical derivation;
none is selected here. The full shape equation, causal constraints,
nonlinear solution and physical stability remain unresolved.

Remaining to-do list: full first-order shape and gravitational boundary
phase, static source compatibility, continuous causal response and
constrained stability; nonlinear existence and signed compatibility;
exterior/state, charge/return, opening, all-sector stability, transport,
preparation, separate creation and full audit.
