# Uniform outer logarithmic frequency comparison

2026-10-20, E003/S168. Ordinary continuum argument before focused checks.
The full research goal remains incomplete.

At the retained centered seed, there are finite constants A,B, independent
of thermal frequency and distance of the support from the outer seam,
such that every smooth real spherical Euclidean tangent T supported in
0<Y<.02 satisfies

    Qcan[T]-C/(30pi) integral(omega+19/60)(Delta T)^2
       >= C[J2/(100pi)-A J1-B J0].                         (1)

The integrals are in the normalized canonical orbit metric; the physical
half convention and Jj are those of
[the canonical frequency core bound](canonical-frequency-core-bound.md).
Qcan includes the massive tower, fixed finite contacts and LLL. It does
not include the classical action or the full four-dimensional Weyl action.
The displayed weighted subtraction is the specific logarithmic fourth
term, not an assertion about the rest of that coupled action.

This resolves the missing comparison at the level of uniform finite
constants and asymptotic annuli. A and B are not numerically calibrated;
no numerical collar radius or physical disturbance threshold follows.
The result concerns this seed and its fixed state/domain. It does not
solve the constrained retarded equations or transfer to a nonlinear
semiclassical exterior. S168.3-S168.4 remain open.

## 1. The two-dimensional scalar comparison

Use the same complete even auxiliary cylinder as in
[the curved quantum form](stationary-curved-quantum-tt.md), section 7.
Rescale the tangential coordinate by the constant a(0), so that

    g_aux=dY^2+exp(2h(Y)) d tau^2,
    h=h_ext, |h|<=.04, |h'|<=1.2, |h''|<=42,
    |h'''|<=3270 almost everywhere,
    L=-Delta_aux=-partial_Y^2-h' partial_Y-exp(-2h)partial_tau^2.

The unchanged thermal circle has length Ltau; it is positive and fixed
(the retained benchmark has circumference far above one). This constant
rescaling does not change the clock or thermal period. The metric is
C2,1 across the cubic seam. It agrees with the actual double for
|Y|<=.1. Work on the full double with even extensions of T and divide
forms by two at the end. No boundary image contribution is silently
discarded. T vanishes near the seam for each core element; estimates
below are uniform as that neighborhood shrinks.

Define the scalar comparison form by the established multiplier,

    Qsc[T]=<T,Kflat(sqrt(L))T>,
    Kflat(p)=-C F_eff p^4/(120pi)
       +sum_n d_n M_n^2/(2pi) integral_0^1 R(alpha(1-alpha)p^2/M_n^2)d alpha,
    R(z)=log(1+z)-z+z^2/2.

The following is the key comparison. For each epsilon>0 there are finite
K_epsilon,L_epsilon depending only on the fixed geometry, angular data,
period, cutoff and common subtraction, such that

    |Qgeo,actual[T]-Qsc[T]|
       <= C[epsilon J2+K_epsilon J1+L_epsilon J0].          (2)

It holds for all thermal Fourier sums, without a frequency-dependent
constant. We prove it before taking Fourier coefficients. This is why
the growing potential in an individual Fourier sector is not estimated
as a bounded perturbation of the one-dimensional flat Laplacian.

## 2. Heat kernel and geodesic kernel have the same principal part

For one mass put a_alpha=alpha(1-alpha) and

    beta_M(t)=M^2/(2pi t) integral_0^1 exp(-t M^2/a_alpha)d alpha.

The exact summed weight obeys beta(t)<=.8 C pi^-1 t^-3 exp(-648t).
Away from the diagonal, Qsc has kernel -integral beta(t)K_L(t;x,y)dt.
In flat two-space, G(t,r)=(4pi t)^-1 exp(-r^2/(4t)) gives

    integral beta_M(t)G(t,r)dt
       =M^3/(2pi^2 r) integral_0^1 a_alpha^-1/2 K1(Mr/sqrt(a_alpha))d alpha
       =M^4/(2pi^2)[K1(Mr)^2-K0(Mr)^2] = -B_M(r).          (3)

For the second equality put alpha=(1+tanh u)/2. The integral becomes
2 integral_0^infinity sech(u)K1(2Mr cosh(u))du. The product formula
[DLMF 10.32.17](https://dlmf.nist.gov/10.32.E17), followed by
K2(z)-K0(z)=2K1(z)/z, identifies it with
Mr[K1(Mr)^2-K0(Mr)^2]. Thus the normalization matches the full rank-two
Dirac bubble exactly, not merely its leading high-momentum power.

Let j(x,r,e) be the exponential-map volume ratio, so
dvol=r j dr dtheta. For flow lengths below a fixed normal radius rho,
the local geodesic kernel is B_M(r)/j. Choose rho small enough that
all balls centered at |Y|<=.02 lie in |Y|<.1 and are normal, uniformly
also across the seam. The established C2,1 normal charts provide this.
Curvature bounds and the Jacobi equation imply, with fixed finite k,

    |Delta r-1/r|<=k r,  |j-1|<=k r^2,  k rho^2<1.          (4)

Here and below constants may be enlarged by fixed factors. None depends
on a Fourier frequency or the center's distance from Y=0.

We give the local heat comparison rather than assume a smooth heat
expansion across the seam. Write delta=Delta r-1/r. For
F=exp(b r^2+c t)G, direct differentiation yields

    (partial_t-Delta)F/F
       =c-4b-4b^2 r^2+2b r^2/t
                           +delta[r/(2t)-2b r].           (5)

Using |delta|<=k r, the choices

    F_plus=exp(2kt+k r^2/4)G,
    F_minus=exp(-2kt-k r^2/4)G

are respectively super- and subsolutions in the normal ball. They have
the same unit delta initial mass as G. The maximum principle therefore
compares them with the heat kernel killed at the ball boundary; a
spatially constant correction of size t^-1 exp(-c_rho/t), increasing
for sufficiently small t, accommodates the boundary values. For the
full heat kernel the same correction bounds returns from that boundary.

For completeness, those full boundary values have a uniform Gaussian
bound without an assumed spectral gap. Conjugating by exp(h/2) gives
-partial_Y^2-exp(-2h)partial_tau^2+V, with |V|<=22. On the universal
cover its heat kernel is exactly

    [exp(h(Y)+h(Y'))]^-1/2 (4pi t)^-1/2 exp(-(Y-Y')^2/(4t))
       times E[exp(-integral_0^t V(Y_s)ds)
               (4pi A_t)^-1/2 exp(-(tau-tau')^2/(4A_t))],
    A_t=integral_0^t exp(-2h(Y_s))ds.                       (6)

Y_s is the variance-two Brownian bridge. This follows by Fourier
integration of the bounded-potential
[Feynman-Kac formula](https://arxiv.org/pdf/math-ph/9907022v1),
Theorem 1.1, with its generator rescaled. Since
exp(-.08)t<=A_t<=exp(.08)t, (6) is bounded by a fixed multiple of
exp(22t)t^-1 exp(-c distance^2/t). Sum all thermal windings to obtain
the cylinder kernel; they remain in this estimate. At small t their
nonlocal contribution is exponentially small, and at large t the
bound grows at most exp(22t)(t^-1+t^-1/2). This is integrable against
beta(t). The Brownian formula also establishes the conservative scalar
semigroup used below; no actual Dirac state is replaced by this scalar
model.

Combining the barriers, boundary correction and j estimate, for r<rho/4,

    |K_L(t;x,y)-G(t,r)/j|
       <= c1(t+r^2)t^-1 exp(-r^2/(c2 t))
                         +c3 t^-1 exp(-c4 rho^2/t)         (7)

at sufficiently small t. All ci are finite and uniform. The calculation
uses only bounded curvature and the fixed normal radius, so smooth
even approximations with uniformly bounded C2,1 norms pass to the seam
limit. It does not assume fourth or higher seam jets.

Integrating (7) against beta gives a separated-kernel error bounded by
K C/r^4 plus a bounded kernel. Indeed the integrals of t^-3 and
r^2 t^-4 against exp(-r^2/(c t)) both scale as r^-4. Large times and
the boundary correction contribute finite constants. This is a gain of
two powers relative to the principal r^-6 kernel.

## 3. Fixing the finite second contacts

A separated-kernel bound alone does not prove (2). The fourth contacts
are the same common restoration. Explicitly,

    Qgeo,4-Qsc,4 = C F_eff/(360pi) integral Rbar |grad T|^2, (8)

on the common local geometry. This is already a bounded second-order
form. No new fourth coefficient is selected. Both complete comparison
models annihilate the constant BASE scalar, which fixes the zero contact.
Their restricted-collar escape terms still have to be retained.

Fix the remaining second contact using a single compact coordinate
quadratic at each center, q_x(z)=chi(z)(e.z)^2, with a uniformly bounded
smooth cutoff. Use the global warped coordinates locally, rescaled by
the center's bounded positive warp. Its gradient vanishes at the center.
The bounds on h through its third weak derivative imply

    |L^2 q_x(z)-L^2 q_x(x)|<=H distance(x,z)                (9)

locally, with H finite uniformly in x, e and the seam distance. The
potentially discontinuous h''' occurs multiplied by partial_Y q_x,
which vanishes at x and is O(distance). All other coefficients in this
expression are Lipschitz. Thus (9) does not require h'''' or Lipschitz
h'''. Away from the cutoff center, boundedness extends the estimate
after enlarging H. This also defines the center value through the seam.

The scalar diffusion has bounded drift h' and bounded uniformly elliptic
diffusion coefficients (1,exp(-h)). Its first displacement moment is
O(sqrt(t)) for t<=1. For larger times, boundedness of L^2 q_x suffices.
There is therefore a uniform finite H1 such that

    |(exp(-sL)L^2 q_x)(x)-L^2 q_x(x)|<=H1 sqrt(s).

The twice-integrated semigroup identity now gives

    |(exp(-tL)q_x)(x)-q_x(x)+t Lq_x(x)-t^2 L^2q_x(x)/2|
       <=(4H1/15)t^(5/2).                                (10)

This identity is valid for the weak bounded L^2q_x as well, by the
conservative kernel formula and approximation. Integrating the exact
single-mass heat weight bounds its subtracted remainder at q_x by

    H1 sqrt(pi)/(2048 M^3).                               (11)

The coefficient follows from Gamma(5/2) and
integral_0^1 [alpha(1-alpha)]^(5/2)d alpha=5pi/1024.
Since sum d_n/M_n^3<2C/3, the scalar contact is uniformly finite.

For the geodesic remainder, the same cutoff has uniform jets. The
geodesic equation for a C2,1 metric gives bounded fourth weak curve
derivative. Since grad q_x(x)=0, Taylor expansion of q_x along the
curve through degree four has remainder O(r^5), uniformly. Odd jets
cancel in the angular cosine average. The exact absolute Bessel moment
integral r^6|B_M(r)|dr=225/(4096 M^3) therefore controls the local
remainder. The long-flow part has exponential mass decay and bounded
cutoff/polynomial jets; for M>=18 it too is bounded by a fixed multiple
of M^-3. This is the finite-contact mechanism already used in section 6
of the curved seed proof, now stated for every coordinate quadratic.
It proves that the second-contact tensor of the difference is bounded
uniformly. It does not infer a diagonal extension from (7).

The geodesic long-flow row mass beyond rho/4 is finite by Liouville
invariance and the Bessel tail bound, including returns and conjugate
points. The scalar long-distance row mass is finite by (6) and the
heat weight. When the two different global completions are identified
only in the collar, keep both escape terms and their row-sum contacts.
They cost a finite multiple of C J0. Constant annihilation does not
erase them.

We now have a symmetric gain-two kernel, fixed bounded second contacts,
the known curvature contact (8), and bounded zero-order tails. Apply
the two-dimensional geodesic difference estimate from section 6 of
the curved seed proof. Splitting at a sufficiently small fixed radius
b absorbs its fourth-order cost into any chosen epsilon J2, leaving
finite K_epsilon J1+L_epsilon J0. The angular Hessian identity controls
the second directional derivatives. This proves (2) with constants
uniform in every thermal Fourier sector and every outer core support.

## 4. Actual logarithmic bound and escaping annuli

The preceding canonical result supplies

    |Qcan-Qgeo,actual|<=C[.001 J2+3e24 J1+3e26 J0].

The scalar logarithmic theorem from that same note, applied to each
thermal mode and summed, gives

    Qsc-C/(30pi) integral(omega+19/60)(Delta T)^2
       >= C[(33/(600pi))J2-(12/(30pi))J1].

All operators and weights are stationary, so the Fourier sectors are
orthogonal, with their actual time norms. Choose epsilon=.001 in (2).
Since 33/(600pi)-.002>1/(100pi), (1) follows with

    A=3e24+K_.001+12/(30pi), B=3e26+L_.001.                 (12)

These are finite definitions, not evaluated numerical constants. The
result is stronger than a fixed-support estimate: the same A and B
apply on all annuli tending to the outer end and at all Matsubara
frequencies. No continuation to arbitrary complex frequencies follows
merely from this Euclidean statement.

For a concrete asymptotic consequence take nonzero
phi in C_c^infinity((1,2)) and T_d=d^-1/2 phi(Y/d)cos(nu_d tau).
Choose a nonzero Matsubara sequence with nu_d d->kappa finite, kappa>=0.
For kappa>0 nearest-integer rounding gives error at most pi d/Ltau;
for kappa=0 a fixed nonzero mode suffices. Since h(0)=h'(0)=0,

    J0 = time_norm [integral phi^2+o(1)],
    d^2 J1 = time_norm [integral(phi'^2+kappa^2 phi^2)+o(1)],
    d^4 J2 = time_norm [integral(-phi''+kappa^2 phi)^2+o(1)].

The last limit is strictly positive. Hence J1/J2=O(d^2), J0/J2=O(d^4),
and the liminf of the left side of (1) divided by C J2 is at least
1/(100pi). The actual canonical form supplies the logarithmic outer
compensation on this Euclidean seed core. This is not a normalized
physical energy or a permitted arbitrarily short-distance experiment.

The next stability step is to combine this result with the full local
frequency polynomial, including mixed/conformal terms, and impose the
initial/current constraints. The constraint-reduced retarded problem,
its boundary realization, transfer to a nonlinear exterior and numerical
control at an admitted physical scale remain required. S167's signed
coefficient, nonlinear solve and compatibility remain unresolved.

## Verification scope

The checker pins the preceding evidence and this ordinary argument.
Exact algebra checks the radial barriers, anchored quadratic jets,
heat moments, common fourth contact, final margin and annular scaling.
Independent positive-integral/Bessel fixtures check the kernel identity;
an exactly solvable nonconstant warped metric checks the half-density
heat normalization and generator. These fixtures do not
compute the actual exterior spectrum or calibrate A,B. No Lean run or
formal continuum certification is claimed.

Remaining to-do list: coupled retarded reduction and initial/current
constraints; transfer or justified exclusion at an admitted physical
scale; signed coefficient, cold/higher response, nonlinear solve and
compatibility; exterior/state, charge/return, opening, stability, transport,
preparation, separate creation and full audit.
