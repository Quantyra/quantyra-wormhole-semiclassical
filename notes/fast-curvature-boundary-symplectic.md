# The variational pairing of the fast curvature traces

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The two outer fast curvature amplitudes carry a nonzero antisymmetric
boundary pairing of the retained gravitational Hessian. Its coefficient
is fixed by the original action. The complete auxiliary thermal shape
elimination leaves this pairing unchanged. Thus a fixed-source
variational domain cannot allow both amplitudes independently merely by
adding a boundary functional. A maximal linear reflecting trace space
has one real phase relation; the action does not select its phase.

This is a necessary boundary-domain result at a fixed admitted thermal
frequency. It does not construct a renormalized action on every fast
metric, identify a Lorentzian physical phase space, or prove stability.

## 1. The full conformal Hessian and its boundary momenta

Keep the charged seed, physical H=L2(dy), state/domain, action coefficients
and contacts of the [metric trace result](coupled-fast-metric-boundary-traces.md).
The inward optical coordinate has r_y=-F. Fix an admitted positive
Matsubara frequency gamma, and write

    epsilon=G/(60pi), E=1/epsilon+11lambda/3,
    d0=E+4lambda, k^2=E/lambda-9/4,
    H_R=2R+rR_r, t=F/r, V=gamma^2+F F_r/r,
    a=6epsilon/F, b=6+epsilon R/6+(22epsilon/3)j,
    c=F(2R+48lambda0), lambda0=lambda+11epsilon lambda^2/6.

For q=r p with the independent shape variable v zero, the complete
conformal Hessian from the [Schur calculation](thermal-curvature-shape-schur.md)
is the polarization of

    c_gamma[q,q]=(1/(2G)) integral {
          a(q_yy-Vq)^2-b[(q_y+tq)^2+gamma^2 q^2]-c q^2 }dy. (1)

The quadratic linearized action is one half of (1). The canonical Dirac
functional is independent of this four-dimensional conformal variable;
its contribution is not omitted from the shape block. The local Weyl
action, classical action and seed density contacts all enter (1).

Set Jq=q_yy-Vq and Dq=q_y+tq. Integration by parts gives the increasing-y
endpoint one-form for that quadratic action,

    Theta_y(q;delta q)=(1/(2G))[Pi1(q) delta q+Pi2(q) delta q_y],
    Pi2(q)=a Jq,
    Pi1(q)=-b Dq-(a Jq)_y.                              (2)

The outer boundary has the opposite orientation, so Theta_out=-Theta_y.
For two real variations define its antisymmetric derivative by

    Omega_out(q1,q2)=(1/(2G))[
       q1 Pi1(q2)+q1_y Pi2(q2)-q2 Pi1(q1)-q2_y Pi2(q1)]. (3)

For complex amplitudes conjugate every occurrence belonging to q1.
This fixes the sign and the normalization; a Hessian and a quadratic
action must not be confused by an extra factor of two.

## 2. Evaluation on the actual charged maximal domain

Use U in the maximal scalar graph domain and q=R_gamma U. The exact
reconstruction equation, not a vacuum substitution, gives

    K_gamma q=FU/6,
    Pi2(q)=-epsilon(U+H_R q),
    Pi1(q)=epsilon U_y+(epsilon H_R-b)q_y
                              +(epsilon (H_R)_y-bt)q.   (4)

This also defines the momenta for graph-domain data without demanding
a pointwise third derivative of an arbitrary L2 source. The preceding
trace theorem gives

    U=sqrt(y)[A_c cos(k log(y/y0))+A_s sin(k log(y/y0))]
                                                     +o(y^(3/2)),
    q=U/(6d0)+z, z=o(y^(3/2)), z_y=o(y^(1/2)),          (5)

with the corresponding differentiated remainder for U. The exact charged
coefficients obey

    H_R=-24lambda+O(y^5), b=6-2epsilon lambda+O(y^3),
    t=1/y+O(y), (H_R)_y=O(y^4).                         (6)

All subleading terms in (3)-(6) tend to zero. For example, a momentum
remainder o(y^(1/2)) pairs with q=O(sqrt(y)), while a Pi2 remainder
o(y^(3/2)) pairs with q_y=O(y^(-1/2)). The potentially singular t term
cancels antisymmetrically even before taking the limit. The leading
momenta can equivalently be written

    Pi2=-(6+22epsilon lambda)q+o(y^(3/2)),
    Pi1=24epsilon lambda q_y-(6-2epsilon lambda)q/y
                                                     +o(y^(1/2)).

Since the cosine/sine basis has Wronskian k, the limit is

    Omega_infinity(U1,U2)=C_fast
          [conj(A1_c) A2_s-conj(A1_s) A2_c],
    C_fast=epsilon k/(12G d0)>0.                       (7)

The two amplitudes have nonzero rank-two pairing. Its coefficient does
not depend on gamma, the seed mass or the subleading charge coefficients.
Those terms were retained until their controlled boundary limit.
The formula is per exterior and per the same time-mode norm as (1).

For precision about operator domains, let C_gamma,max be the strong
conformal row paired back with R_gamma^*: on the maximal graph domain,

    C_gamma,max U=epsilon R_gamma^*
            [Atilde+gamma^2+C_gamma,strong-epsilon^-1 P R_gamma]U/(2G),
    P=F(r calT_r+4calT).

Here C_gamma,strong is the preceding bounded local trace feedback,
not the variational operator on the left. The residual contacts are
unchanged. Its source is in H, so this expression is well defined.
The horizon limit-point data have zero endpoint pairing. Integration
on finite intervals followed by (4)-(6) gives the Green identity

    <U1,C_gamma,max U2>-<C_gamma,max U1,U2>
                                             =Omega_infinity(U1,U2). (8)

The individually divergent terms of (1) are not assigned values for
all maximal-domain inputs. Equation (8) uses strong H rows and the
finite limit of their boundary pairing instead.

## 3. Why the complete auxiliary shape cannot remove the pairing

The retained mixed map K_cross:E_gamma->H is bounded, with adjoint
K_cross^dagger:H->E_gamma. In the shape energy inner product the
auxiliary block rows are

    (C_gamma,max U+K_cross w, K_cross^dagger U+w).

Their off-diagonal pairings cancel in the antisymmetric Green identity,
and the identity shape block contributes zero. Thus (8) holds for the
full pair on D_max x E_gamma. This is a form-dual calculation; it does
not assume that every energy-class shape has pointwise higher normal
derivatives at the boundary.

Equivalently, shape elimination subtracts

    M_gamma=K_cross K_cross^dagger=M_gamma^*,
    C_gamma,reduced=C_gamma,max-M_gamma.

A bounded self-adjoint term on H contributes zero to (8), regardless
of its norm or detailed kernel. Compactness, although already proved,
is not needed for this cancellation. The complete canonical response
therefore remains included in the auxiliary variational statement.
This is not a replacement of the strong nonsymmetric curvature feedback
by a positive operator: it uses the actual congruent Hessian pairing.

## 4. What boundary functionals can and cannot determine

At each finite cutoff, adding any twice differentiable boundary
functional B changes the integrated boundary one-form by delta B.
Its antisymmetric second derivative is zero. This includes the existing
GHY, Weyl-transgression and counterterm functionals. Tangential exact
ambiguities integrate to zero over the thermal circle and sphere; for a
finite time strip the corner terms must instead be retained. These are
the standard distinctions in [Iyer--Wald, equations (40)-(47)](https://arxiv.org/pdf/gr-qc/9403028).

Consequently a C2 variational completion whose fixed-source boundary
one-form vanishes on all admitted variations must have zero pullback
of (7). It cannot admit independent A_c and A_s on one exterior with
no additional boundary equations or boundary degrees of freedom. A
counterterm alone cannot cancel their rank-two pairing. This necessity
does not claim that the full action is already finite on those metrics,
or forbid a different boundary theory with explicitly matched dynamics.

For a complex linear trace subspace preserving real perturbations,
maximal isotropy of (7) is exactly one real homogeneous relation,

    cos(theta) A_c+sin(theta) A_s=0, theta modulo pi.    (9)

To see this, a complex line spanned by (a,b) is isotropic iff
conj(a)b-conj(b)a=0: its ratio is real, including the vertical line.
No two-dimensional subspace is isotropic. The zero trace subspace is
isotropic but not maximal. Combining with the preceding index result,
the maximal reflecting trace choices have auxiliary index zero, while
the zero trace choice retains its index -1 compatibility condition.
Neither statement proves existence, uniqueness, positivity or stability.

Equation (9) is necessary for this maximal, linear, separated reflecting
choice at a fixed thermal frequency; the physical value of theta is
undetermined. A relation linking the two exteriors, or boundary dynamics
with frequency dependence, requires its own matched variational problem.
No such choice is inferred from the matter reflector or leading metric.

## 5. A cutoff does not silently select a phase

The established inverse-square normal model makes the remaining choice
explicit. For U=sqrt(y)cos(k log(y/y0)+phi), define its projective
logarithmic derivative wherever U is nonzero:

    rho(y)=y U_y/U-1/2=-k tan(k log(y/y0)+phi),
    d rho/d log y=-(rho^2+k^2).                         (10)

It has logarithmic period pi/k and poles that are coordinate singularities
of the trace line. The normalized pair (U/sqrt(y),
(yU_y-U/2)/(k sqrt(y))) passes smoothly through these poles as a rotation.
A fixed rho at a shrinking cutoff instead changes phi modulo pi.
Two geometric cutoff sequences can yield two distinct limiting trace
lines. Thus a fixed cutoff Robin number does not specify a convergent
boundary phase. This is a normal-model statement, not a computed charged
resolvent or an imposed boundary counterterm.

This limit-cycle behavior is established for the inverse-square problem;
see [Bawin--Coon, equations (3)-(5)](https://arxiv.org/pdf/quant-ph/0302199).
Here its relevance is tied to the gravitational pairing (7). At the
retained benchmark G lambda<(4e-17)^2, k^2>1.125e35 and pi/k<1e-17.
This bounds the phase period in log cutoff, not a physical length,
growth time or validity scale of the effective theory.

## 6. Verification and physical scope

The focused checker differentiates the original nonlinear optical
classical/Weyl densities to recover the conformal boundary pairing,
independently of the integrated form (1). It compares their antisymmetric
momenta for arbitrary charged radial jets, checks the reconstruction
identities and exact asymptotic coefficient, and tests actual rational
charged exterior coefficients at successively larger radii and two
precisions. Symmetric Schur controls, phase-line algebra and independent
integration of the normal-model phase rotation test the domain claims.
The fixture parameters are not a calibrated physical benchmark and do
not evaluate the full quantum kernel. The continuum endpoint and
variational statements remain ordinary proofs; no Lean run is selected.

The new result supplies the full auxiliary variational meaning of the
previous scalar trace relation. It still does not derive a physical
choice of boundary ensemble, a causal flux law on a nonlinear solution,
the remaining tensor/horizon conditions or physical stability. The seed
remains off shell. S168.3-S168.4 and every S167 task stay unchecked.

Remaining to-do list: matched gravitational boundary rule and source
compatibility; continuous causal response and constrained stability;
nonlinear existence and signed compatibility; exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and full audit.
