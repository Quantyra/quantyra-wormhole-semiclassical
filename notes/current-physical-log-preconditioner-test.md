# A physical logarithmic preconditioner fails the coarse uniform defect test

2026-10-28 HST, E003/S169. Ordinary argument followed by outward arithmetic.
The newly resumed goal is ACTIVE and incomplete. Preserve the full stable
semiclassical wormhole objective, P32, closed audit and original norms.
This tests one explicit infinite-range reference before implementing its
inverse to extreme precision. It is a uniform-preconditioner rejection,
not a lower bound on the actual structured correction or physical error.

## 1. Complete operator and the proposed reference

Use H=L2((0,z_h),z^-4 dz), E_t with squared norm ell_t/N,
N=|W_action| beta, and the complete weak operator

    A_Q=(ell_Q+B)/N : E_t -> E_t*, B=b* a^-1 b.

The scalar inverse and mixed rows are those of BOTH original equations.
The accepted quadrature comparison has |ell_Q-ell_t|/N <= e_Q R_t,
e_Q<.0022. Nothing here discards its quantum action. Keep the paired
first-order domains Pi_- psi=0 for D and Pi_+ phi=0 for Ddagger,
Pi_+=(1+tau_3 gamma_n)/2. The square has Pi_- psi=0 and
Pi_+(nabla_n+K/2)psi=0, with exchanged partner projectors. The regular
thermal bolt has no wall. The common regulator, exact ultraviolet log,
all angular/thermal modes, coincidence contacts, cutoff transitions and
finite returns are inherited from the accepted COMPLETE form comparisons.

Define the comparison scalar Laplacian on the ACTUAL physical exterior:

    L=-Delta_4/lambda
     =-z^2/beta [Q d_z^2+(Q'-2Q/z)d_z].

It acts in H, with its complete-end Friedrichs realization and smooth
rotational disk behavior at the bolt. It is not the canonical-double
Laplacian and not a replacement quantum field or state. Put g=G lambda,
Fstar=F_lambda+23g/(180pi), kappa=-6/5-log(mu0/sqrt(lambda)), and

    H_M(s)=psi(s-1)+(pi/2)cot(pi s/2)+kappa,
    s(R)=3/2+sqrt(9/4-R),
    f(R)=(3/2)Fstar R+(3g/(20pi))R(R-2)H_M(s(R)),
    P_phys=f(L)^-1.                                      (1)

For R>9/4 choose the positive imaginary square root. Reflection makes
the value real and branch-independent; remove the pole at R=2 in the
PRODUCT. The accepted scalar-coordinate TT normal family is
9 lambda^2 N_T(s). Its H normalization multiplies it by 8piG/lambda,
giving exactly (1) when s(s-3)=-R. Thus (1) retains the full outer
normal logarithm and its finite constant, not just its leading growth.
The actual finite-geometry difference A_Q-f(L) remains in the defect.
Spectral notation specifies an exact comparison, NOT an implemented solve.

## 2. Domain before numerical interpretation

The physical seed satisfies F''>=lambda, F'>0 and F(r_h)=0. Integrating
2F'F''>=2lambda F' gives F'^2>=2lambda F+F'(r_h)^2. In proper radius,
the radial volume density m is proportional to sqrt(F) r^2, and

    (log m)'=F'/(2sqrt(F))+2sqrt(F)/r >= sqrt(lambda/2).

Integration of (m u^2)' for a compact-outer, bolt-regular scalar gives
sqrt(lambda/2)||u||^2 <= 2||u||||grad u||. The bolt term vanishes
because m=O(rho); no Dirichlet condition is imposed there. Closure gives

    L>=1/8.                                               (2)

For 1/8<=R<=9/4 the accepted whole-strip bound is
|(q+2)H_M|<400, q=-R. For R>=9/4 set t=sqrt(R-9/4). Then

    H_M=Re psi(1/2+it)-(pi/2)sech(pi t)+kappa>0.

The digamma series [DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6), valid
since 1/2+it avoids its poles, gives
Re psi(1/2+it)-psi(1/2)
=sum_(n>=0) t^2/[(n+1/2)((n+1/2)^2+t^2)] <=7 zeta(3)t^2<9t^2.
Use 39/5<kappa<54/5, psi(1/2)=-EulerGamma-2log2, and
zeta(3)<1+1/8+integral_2^infinity x^-3 dx=5/4. Positivity follows
from monotonicity of both displayed t-dependent terms and their value
at zero. Consequently, for every R>=1/8,

    a R <= f(R)
       <=(3/2)Fstar R+(3g/(20pi))(400R+14R^3),
    a=(3/2)(1-1.179e-8)-60(1e-33)/pi>0.                  (3)

In particular P_phys:H->H is defined and has norm <=8/a<5.34.
This is not an original physical-output comparison.

Here is the additional space identification needed for a WEAK defect.
The actual physical exterior is smooth in disk charts; its compactified
mass/log boundary is at infinity, not replaced by a smooth canonical
double. On rescaled end charts L is uniformly elliptic in log z after
the critical H conjugation. Its coefficients and all conormal derivatives
are bounded; in the middle and regular disk charts the scalar elliptic
estimates apply. Integer graph norms of L are the corresponding patched
physical Sobolev norms. Applying the spectral dyadic decomposition to
these graph norms supplies the logarithmically weighted norms: f(R)
has positive growth comparable to R+g R^2 log(2+R), so

    Dom(f(L))=W_(3/2)=z^(3/2) H_b^(4,log)

with the original regular disk charts. Constants are finite, not numerical.
The fixed-seed full-Hessian argument, section 5, identifies the same
operator domain for the actual TT block. Its genuine-TT core includes
the outer tail compensator and the mass cubic; it is not the punctured-
bolt core. Closed graph gives equivalence of the two graph norms.
Interpolation between H and these positive graph domains identifies
their square-root domains. The retained transport/quadrature form
comparisons then identify this common form space with E_t. Thus f(L)
is a bounded positive coercive map E_t->E_t*, and its inverse maps
E_t*->E_t. This supplies qualitative domain membership only. No
numerical E_t inverse norm, physical radius return, boundary coefficient
norm or nonlinear-neighborhood comparison is inferred.

No theorem about smooth metric variations of the canonical double is
used. The finite seed realization and its accepted complete comparisons
are inputs, not a new regulator-removal theorem on a neighborhood.

## 3. An admissible whole-operator test

Use psi=z^4, J=0, Upsi=z^4/2, Vpsi=-z^4. It solves the exact York
equation -(Q^2 J'/z^4)'=Q(4psi/z-psi')/z^4 with J(0)=0 and zero
bolt flux. It is smooth in the actual regular disk and has zero mass
coefficient, an allowed member of the paired mass/log domain. The
accepted TT-preserving cutoff/compensator argument includes its outer
transitions and returns. No cutoff error is dropped by treating a
finite interval as the original end.

The accepted full-form quartic receipt provides a known diagonal k,
a positive complete geodesic remainder, and a scalar variational lower
bound b. Restricting the scalar maximization to z^4,z^5 LOWER-bounds
the full B, rather than replacing its inverse by a matrix. Therefore

    A_Q[psi,psi] >= (1-.0022)(1-eps)/(1+.0022) k+b,
    eps<1.227e-74,
    k=2.26666663457616857965...,
    b=0.259726029048268612....                            (4)

The first .0022 allows the actual ell_Q comparison, the denominator
allows the complete geodesic-to-actual form comparison. This deliberately
conservative combination includes all retained quantum error. Positivity
of the whole geodesic remainder is used only in this lower bound.

From (3), compute WITHOUT spectral truncation

    m1=<psi,Lpsi>H=beta^-1 integral Q z^-2 (psi')^2 dz,
    m3=<psi,L^3psi>H=beta^-1 integral Q z^-2 ((Lpsi)')^2 dz,
    <psi,f(L)psi>H
      <=(3/2)(1+1e-30)m1+(3e-33/(20pi))(400m1+14m3).   (5)

Both moments are finite polynomial integrals on the exact correlated
geometry. The endpoint identity Q(z_h)=0 is retained exactly; subtracting
uncorrelated tiny endpoint intervals is not the bolt proof. The massless
trial and its L image are in the scalar form domain, so (5) is legitimate.

For the defect D=I-A_Q P_phys on E_t*, let C and M be the f(L) and A_Q
Riesz representatives in E_t. They are bounded positive and invertible
by section 2 and accepted coercivity. The defect is similar to
I-C^-1/2 M C^-1/2. Its spectral radius, hence its norm in the ORIGINAL
dual, is at least any generalized Rayleigh quotient minus one. Thus

    ||D||_(E_t*->E_t*)
      >= A_Q[psi,psi]/<psi,f(L)psi>H-1.                 (6)

Unknown numerical norm-equivalence constants do not weaken this spectral
lower bound. They still matter for any upper bound or output certificate.

## 4. Executed result and quantified repair requirement

The accompanying 512-bit outward calculation, with frozen accepted
inputs and runtime, gives

    m1 = 0.5333333333410889287...,
    m3 = 9.42492801360908014...,
    A_Q[psi,psi] >2.516441223600354,
    <psi,f(L)psi>H <0.800000000011634,
    ||I-A_Q P_phys|| >2.14555152945470.                 (7)

The exact reference therefore fails the sufficient .5 test. Application
errors are set to ZERO for this rejection; no numerical solve is claimed.
This is a whole-operator witness using both equations and the complete
quantum comparison, not a selected-mode approximation to that operator.

For a positive corrected reference f(L)+Delta to have defect <=.5, its
quadratic value on this same trial must be at least (2/3)A_Q[psi,psi].
Consequently its increase must exceed .87762748238860 in these unchanged
energy units, a factor >2.09703435296980 over the old reference value.
This is NECESSARY on this one trial, not sufficient on the whole space.
Merely increasing correction count cannot validate a Neumann expansion
whose exact defect already has spectral radius >2.1455. Rescaling is
not a physical-norm improvement and is not certified by this test.

No weak-source allowance, complete F residual, signed G_t or physical
admission margin closes. The sufficient source screen remains
4.9887427328e-32; the actual implemented S/K errors are UNKNOWN.
The failed uniform reference does not refute the restricted map on B_i,
nor the coupled inverse, P32, or existence of the construction.

## 5. Method decision and single next calculation

Retain the construction; change this preconditioner. The test says the
finite-geometry local/York energy must be retained before high accuracy
is pursued. Compare the already accepted COMPLETE compensated scalar-log
form A_sc: it keeps those terms and B and has complete relative difference
e_SQ<.000654 from A_Q and coercivity >=(1-.002784)R_t. Its EXACT inverse
would have defect <.000656, but that inverse is still unimplemented.
These established form comparisons are evidence for this choice, not a
claim of an implemented solver or attained application error.

The single next calculation is a coarse certified full-range solve for
A_sc: construct an explicit analytic logarithmic complement plus a
finite correction retaining the complete local/York rows, and bound
||I-A_sc P0|| in the original dual INCLUDING complement, interface,
boundary and numerical errors. Target .49 leaves room for transfer to
A_Q: e_SQ<.000654 and coercivity >.997216 imply
||I-A_Q P0|| <=.49+.000654*(1+.49)/.997216 <.491 <.5.
This calculation is PROPOSED. Neither f(L)^-1 nor A_sc^-1 is a ready
numerical implementation. A fixed finite output alone is already rejected;
the analytic complement is an essential part of the proposed solver.

Only after that gate may ordered corrections be budgeted toward the
source screen, with separately certified applications and output maps.
The closed audit and original physical norms remain unchanged.

Remaining to-do list: full-range implemented inverse and complete weak
source certificate; S/K application errors and signed G_t; original graph,
physical, radius, mass, phase and coupled budgets; nonlinear admission and
contraction, exact compatibility, physical matching, semiclassical
validity/noise and all-sector stability; the full stable-wormhole objective.
