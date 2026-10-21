# First thermal frequency: full Dirac horizon flux

2026-10-21, S168. This is an ordinary argument and a reproducible algebraic
and numerical audit, not a formal proof certificate. It extends the
[zero-mode cancellation](causal-horizon-and-zero-mode-response.md) to one
horizon component of the full retained free Dirac response.

## Statement and fixed hypotheses

Use the actual smooth, simple-horizon
[centered seed](stationary-centered-covariant-seed.md), its prepared HH
state with beta=2pi/kappa_h, identical paired spherical metric sources,
the unchanged zero-source paired boundary realization, and the common
four-dimensional subtraction and fixed local actions. The seed is not
an exact nonlinear semiclassical solution. A source is smooth and has
radial support in a compact subset of the open exterior, separated from
both the horizon and outer endpoint. Its frequency is gamma=kappa_h>0.
There is no independent initial-state perturbation or boundary forcing.

The frequency response means the causal matter resolvent at omega=i gamma,
or the response to exp(gamma t) prepared from t=-infinity. It is not the
entire response to a finite-time switched source; the latter has transients.
The nonzero Matsubara identity identifies this response with the Euclidean
first harmonic, including the prescribed contacts. Then, in the source-free
horizon collar and with the unchanged background used to raise indices,

    lim_(r->rh+) exp(-kappa_h v) delta T^r_v(r,v)=0.       (1)

Here T is the full renormalized Dirac stress, including all massive angular
channels and the declared zero-mode component. The fixed local restoration
has zero variation on this collar. Other components and the bulk response
need not vanish. Equation (1) fixes an exact value; it does not evaluate
the general massive kernel, supply a gravitational boundary condition,
exclude a charged mode, or establish physical stability.

## 1. Thermal response and separated readout

The [thermal spectral derivation](thermal-shape-response-and-subtractions.md)
and [ADM vertex construction](adm-dirac-boundary-response.md) fix the causal
and Euclidean conventions. At finite spectral regulator let H be Hermitian,
C=(1+exp(beta H))^-1, B the input vertex, O a stress readout and W its
mixed second vertex. For gamma=2pi l/beta, l a positive integer,

    (gamma+i ad_H) delta C=-i[B,C],
    K_R=Tr(O delta C)+Tr(C W),
    K_E=Tr(C W)-integral_0^beta exp(i gamma tau)
         Tr[C exp(tau H) O exp(-tau H)(1-C) B] d tau,
    K_R=K_E.                                           (2)

To verify (2), diagonalize H. In each ordered energy pair the tau integral
has denominator E_a-E_b+i gamma. The identity
f_a(1-f_b)(exp(beta(E_a-E_b))-1)=f_b-f_a gives the covariance resolvent
after interchanging the pair labels. Equal-energy blocks integrate to zero
at nonzero l. The static equilibrium derivative has a separate equal-energy
term and cannot be substituted at l=0. The same local contact and common
subtraction accompany both sides when the regulator is removed.

For compact separated input and readout there is no coincident ultraviolet
singularity between the two vertices. Below we also control the angular
sum in a whole bolt collar, so this use of (2) does not assume an optical
mass gap or interchange an uncontrolled angular sum with a horizon limit.
The Euclidean source has harmonic exp(-i kappa_h tau); because its support
avoids the bolt it is a smooth complex metric variation at fixed period
and spin structure. Real cosine and sine variations give its real and
imaginary parts. The background operator, not a complex deformed operator,
is used for the linear response.

## 2. Smoothness of the massive sum at the bolt

In the retained constant-sphere conformal frame, the proper angular masses
and multiplicities are

    M_n=sqrt(n(n+q)/C), d_n=2(q+2n), n>=1, q>0, C>0.    (3)

These are not gaps of the Lorentzian optical Hamiltonians. The latter
retain their gapless horizon spectrum. The paired Euclidean realization
and its local half-density/domain transport are fixed in the
[canonical comparison](stationary-centered-canonical-response.md).
Write S_M=(Dsharp+M)^-1. At the real background Dsharp is anti-self-adjoint,
||S_M||<=1/M, and scalar Lipschitz multiplication preserves its domain.
For a weight phi with |d phi|<=M/2, the conjugated operator differs by
the bounded Clifford multiplication -c(d phi), so its inverse has norm
at most 2/M by the resolvent Neumann series. A bounded weight which changes
by M d/2 between two separated sets gives an exponential off-diagonal
bound. Approximating the distance weight by smooth weights, or decreasing
the exponent slightly, removes any cut-locus issue.

Take a compact bolt neighborhood K and a larger one K1, both separated
from the source support J by positive proper distance. Nested cutoffs in
these neighborhoods and interior elliptic estimates bound all prescribed
derivatives of the separated kernel by a polynomial in M times exp(-c M d),
with c>0. To spell out the parameter dependence: the equation
Dsharp u=-M u off J, differentiated finitely many times, has fixed smooth
geometric coefficients; each use of the equation or a derivative estimate
costs only finitely many factors of (1+M). Sobolev evaluation on the fixed
compact neighborhoods converts these local L2 estimates to kernel
derivatives. The same argument on the other kernel variable gives the
two-sided estimates. No endpoint regularity estimate is needed at the
observation point, since the Euclidean bolt is an interior point.

In the fixed representation delta S=-S delta D S, together with the local
endpoint/density terms of the representation. The latter vanish at an
off-source readout. The differential insertion delta D is supported in J;
its finite order and any angular mass factors cost another polynomial.
Applying the local stress differential operator to delta S on the diagonal
therefore gives, for every finite derivative order j on the smaller K,

    ||delta T_n||_(C^j(K)) <= A_j (1+M_n)^p_j exp(-t M_n), (4)

for finite A_j, an integer p_j>=0 and some t>0. Integration over the compact
source does not change the exponent type. The background conformal factor
and its derivatives are smooth on K; restoring the four-dimensional tensor
preserves (4). The common local counterterms, Weyl restoration, density and
index-raising contacts all vary by zero there because the metric variation
vanishes on a neighborhood. The paired phase convention stays fixed; this
does not introduce an independently adjustable determinant phase.

Here is an explicit summable majorant. Set alpha=sqrt(q/C), b=t alpha/2 and

    K_(0,t)=1, K_(p,t)=(1+2p/t)^p for p>0.

For M>=0, (1+M)^p exp(-tM)<=K_(p,t) exp(-tM/2). This follows by splitting
at M=2p/t: before it the polynomial is bounded by K, and after it the
logarithmic derivative of (1+M)^p exp(-tM/2) is negative. Since
M_n>=alpha sqrt(n), and 2(q+2z) exp(-b sqrt(z)) decreases for z>=4/b^2,
the integral test gives, for N>=max(1,4/b^2),

    sum_(n>N) d_n(1+M_n)^p exp(-t M_n)
      <=4 K_(p,t) exp(-b sqrt(N)) [
          q(sqrt(N)/b+1/b^2)
          +2(N^(3/2)/b+3N/b^2+6sqrt(N)/b^3+6/b^4)].    (5)

The integral follows by z=sqrt(n). Thus the massive sum and every finite
derivative converge uniformly on K. This is smoothness of the entire
massive response, not just of finitely many angular channels. The finite
zero-mode multiplicity is covered by the previous evaluated response;
no positive-mass estimate is applied to M=0. This proves smoothness of the
full Euclidean stress variation at the bolt in the stated source class.
The A_j in (4) are not numerically calibrated here, so (5) alone is not
a numerical error bar for the full physical kernel.

## 3. Smooth harmonic tensor and the mixed flux

Suppress the subscript on kappa. In the unchanged collar use

    ds^2=-F dv^2+2dv dr+r^2 dOmega^2, F(rh)=0, F_r(rh)=2kappa,
    u=t+y, v=t-y, r_y=-F,
    U=-exp(-kappa u), V=exp(kappa v), x=-UV=exp(-2kappa y).

Then F x_r=2kappa x and x=a(r-rh)+O((r-rh)^2), a>0. The Killing field
is xi=kappa(V partial_V-U partial_U). A frequency n kappa response has
Lie_xi delta T=n kappa delta T. In Euclidean signature, U and V become
oppositely rotating complex coordinates at the bolt. A smooth scalar
harmonic of integer weight k is the appropriate complex coordinate to
power |k| times a smooth function of squared radius. This is a smooth
equivariance fact: its Taylor jets contain only the allowed monomials;
division of the flat remainder by the required power stays smooth.
It does not assume a convergent Taylor series or a real-analytic source.

Including the weights of the two covariant indices, the first harmonic has

    T_UU=V^3 A(x), T_UV=V B(x), T_VV=U C(x),             (6)

where A,B,C are smooth at x=0. Angular scalar components have the form
V D(x). Continuation of this single harmonic and its radial amplitudes
gives the exterior Lorentzian frequency response by (2). This step is not
a claimed analytic continuation theorem for arbitrary time-dependent
smooth functions or for a different Dirac state.

Direct covariant transformation uses

    partial(U,V)/partial(v,r)=[[-kappa U,-x_r/V],[kappa V,0]],
    V^-1 T_vv=kappa^2 x(x A+2B-C),
    V^-1 T_vr=-kappa x_r(x A+B),
    V^-1 T_rr=x_r^2 A.                                 (7)

The physical mass flux is the mixed component, not an arbitrarily chosen
covariant component. Raising its first index before taking the limit gives

    V^-1 T^r_v=V^-1(T_vv+F T_rv)
               =-kappa^2 x(C+x A)=O(r-rh).              (8)

This proves (1). In contrast, T_rr/V and T^v_v/V may tend to nonzero values.
The spherical energy Ward identity on this source-free collar is
partial_r(r^2 T^r_v)+kappa r^2 T^v_v=0 after the common exp(kappa v)
factor is removed. Equations (7)-(8) then require C(0)=-B(0); they do
not set A(0) or B(0) to zero. Background off-shell metric contacts vanish
here because the source and all its derivatives vanish here.

## 4. Scope tests and quantitative consequence

Three controls prevent overinterpreting this result.

* Future-horizon regularity alone does not suffice. The conserved ingoing
  spherical null tensor exp(kappa v) dv^2/(4pi r^2) has nonzero mixed flux.
  It is regular at the future horizon for V>0, but T_VV is proportional
  to 1/V and is singular at the bifurcation surface. This is a tensor
  counterexample, not a construction of a different Dirac state.
* Smoothness does not erase the full stress. For the flat product metric
  -2dU dV plus a constant sphere, T_UU=0, T_UV=V, T_VV=-U is smooth,
  divergence-free and of first boost frequency. Its trace is nonzero.
  It is a geometric counterexample, not a value of the charged kernel.
* The next positive thermal harmonic permits T_VV=C(x). More generally
  n>=2 has T_UU=V^(n+2)A, T_UV=V^n B, T_VV=V^(n-2)C, giving
  T^r_v/V^n=kappa^2(C-x^2 A). It can be nonzero at the horizon.

The usable consequence is one exact horizon response value, including
the massive sum, and the regular approach (8). This supplies a boundary
consistency check for future numerical response work. It does not give
a calibrated norm of that response, a factor (gamma-kappa) in a neighborhood
of this frequency, a horizon/infinite-time Laplace interchange, an averaged
null-energy sum rule, a time-integrated full mass-return result, or a full
gravitational mode. Sources reaching the horizon require their local
contacts and possible moving-horizon terms and are outside this theorem.
Independent initial-state forcing and boundary sources are also outside it.

## 5. Source comparison and frozen verification protocol

[Wang and Wang](https://arxiv.org/pdf/2208.01047), section II equations
(6)-(9) and (17)-(18), use established Kruskal boost weights to obtain
horizon restrictions and pole-skipping conditions in local holographic
bosonic theories. No pole-skipping, chaos exponent or gravitational mode
follows here. Our contribution is applying smooth thermal harmonics and
the retained separated Dirac sum to this specific matter horizon flux.

[Sanders](https://arxiv.org/pdf/1310.5537), section 2.3, Proposition 2.11
and Corollary 2.12, explains why regular coordinates matter at the
bifurcation surface. His field existence theorem is for a scalar field;
it is not imported as a theorem about our paired Dirac boundary problem.
[Romatschke and Son](https://arxiv.org/pdf/0903.3946v1), equations (13)-(17),
distinguish metric-response contacts from a bare stress product. No QGP
sum rule is transferred. Our retained source/domain and spectral identities
supply the model-specific assumptions. No novelty is claimed for generic
boost-weight counting, KMS continuation or elliptic off-diagonal decay.

Before acceptance freeze this note, the checker and the source review.
The checker independently transforms the covariant tensor and raises its
index, verifies weights, the horizon Ward coefficient and the controls,
and differentiates the explicit tail integral. It compares the covariance
resolvent with direct Euclidean matrix-exponential quadrature in four
noncommuting three-level fixtures, including degeneracy and the second
thermal harmonic, with 48 and 96 nodes and absolute error below 2e-11.
A nonzero-contact control must detect dropping W. Three tail-envelope
fixtures compare direct finite sums plus an analytic remainder at two
cutoffs with (5). These matrices and envelopes check conventions and
majorants; they are not numerical evaluations of the charged Dirac kernel.
Retain every run and exact failed or superseded source. No Lean campaign
is selected and no scientific source changes after final success.

Remaining to-do list: evaluate the remaining massive response and full
physical domain; decide charged stability; establish nonlinear existence
and compatibility; complete exterior/state, charge/return, opening,
all-sector stability, transport, preparation, separate creation and audit.
