# Global long returns in a logarithmic scale family

2026-10-11, S053. The fixed-trim mouth and photon-source construction
extends to shrinking resonant labels by a polynomial bound on its geometry
constants. Selecting delta=(log q)^-2 makes every finite such loss harmless.
This is a conditional return theorem for the flat dipole with the prescribed
spherical higher-mode load. It does not construct that loaded realization,
match the physical curved metric, or evaluate a vacuum determinant or stress.

## 1. Statement and scope

Keep the operator, gauges, ports, higher-mode load, and positive-resolution
current/resolvent assumptions of [joint mouth matching](joint-mouth-matching.md)
and [source matching](dipole-source-matching.md). All lengths here are in
units of separation d. In particular z=d omega_phys+i nu, nu=d eta_phys.
Assume the loaded identities uniformly over the sectors being retained;
the all-angular-sector collar loss bound supplies the same boundary estimate.
The uncommitted exterior-boundary draft is not an input to this theorem.

Fix u0<1, Omega<infinity, and finitely many axisymmetric Hermitian Pauli
source profiles with fixed compact support away from both centers and bounded
C2 coefficients. Their support, norm bounds and mouth exclusion radii are
independent of q and u. The compact support is retained from the source proof;
no general time-dependent source or angular-mixing theorem is asserted.
For small sigma>0 retain physical labels

    sigma <= u <= u0,       q(1-u) in Z+1/2.

There exist finite positive C, M, N and c, depending on these fixed data,
such that for q>=C sigma^-N, |Re z|<=Omega and 0<nu<=c sigma, the matched
return R and Wilson comparison R_W of the linked source note obey

    ||partial_t^alpha (R-R_W)|_(t=0)||
      <= C_alpha sigma^-M nu^(-|alpha|-1/2)
            [q^-3/8 + q^1/8 exp(-c q sigma^N)],     |alpha|<=2.       (1)

M,N may be increased to cover all these inequalities simultaneously.
They are existence constants, not computed numerical exponents. No explicit
finite-flux threshold follows. Sources vanish near the shrinking overlap
annuli. The norm includes all retained ground inputs with higher modes
still attached; axisymmetry keeps the different azimuthal sectors separate.

Select the same weak-coupling theory family g^2=q^-2 but now take

    delta=d/ell=(log q)^-2, p=q^-3/4, a=c_a q^-1/4,
    eta_phys=delta/ell, nu=d eta_phys=delta^2=(log q)^-4.              (2)

For any fixed finite set of harmonics n=2,...,K, take sigma=c_K delta
below their limiting resonant labels

    u_n = kappa delta/[pi(2n-1)] + o(delta),
    kappa=1+asinh(sqrt(3))/(2 sqrt(3)).                              (3)

Rounding to the nearest physical half-integer sector changes u by at most
1/(2q)=o(delta). Equations (1)-(2) give, for every fixed k and |alpha|<=2,

    ||partial_t^alpha (R-R_W)|| = o(delta^k).                       (4)

The order of limits is fixed K first, then q to infinity. This does not
bound the actual physical contribution of all discarded harmonics. The
previous scalar series tail remains a tail estimate for that comparison
series only. Equation (4) also does not control a frequency/source-dependent
vacuum trace, or source profiles whose dimensionless norms grow with q.

## 2. Geometry with controlled dependence on the trim

Use prolate coordinates xi>=1, -1<eta<1 in the meridional half-plane:

    r=1/2 sqrt[(xi^2-1)(1-eta^2)], z=xi eta/2,
    rho_-=(xi+eta)/2, rho_+=(xi-eta)/2,
    U=xi(1-eta^2)/(xi^2-eta^2).

On U=u,

    eta^2=xi(1-u xi)/(xi-u),
    1-eta^2=u(xi^2-1)/(xi-u),
    r=sqrt(u)(xi^2-1)/(2 sqrt(xi-u)),     1<=xi<=1/u.                (5)

The two signs of eta join smoothly at the midplane. They are the single
embedded source-to-receiver curve in r>0. Their only limiting axis points
are the two centers. Direct differentiation in the prolate metric gives

    b^2=|grad U|^2/r^2
       =16[xi^2+eta^2(3xi^2-1)+eta^4]/(xi^2-eta^2)^4,
    b>=4/xi^3>=4u^3.                                              (6)

These identities are checked symbolically in the accompanying script.
Let rho_*=(1+rho_-^-2+rho_+^-2)^-1/2. For sigma<=u<=u0 the following
properties hold, with constants allowed a finite power of sigma^-1:

* b rho_*^2 is bounded above and bounded away from zero; the lower bound
  can be taken c sigma^3. Near either center the monopole term dominates
  uniformly, while away from fixed small balls (6) applies.
* r/rho_* >= c sqrt(sigma). Near a center, rescale by its distance:
  the flux ray has sin(theta)=2 sqrt(u(1-u)); on the rest (5) applies.
* Every fixed finite collection of scaled derivatives of b, the frame,
  curvature, m=(U-u)/r, and rho_* is bounded by C sigma^-A. For example
  the quantities needed in the old proof are c rho_*^3, c_s rho_*^4,
  omega rho_*, omega_s rho_*^2, kappa rho_*, kappa_s rho_*^2, and
  rho_*^4 times the third normal derivative of m. Here c=m_yy and
  omega=b_s/(4b); subscripts s are unit-tangent derivatives.

For the last assertion, Cartesian derivatives are rational expressions in
r,z,rho_+,rho_-,r^-1 and b^-1, with square roots for the positive norms.
In a center chart x=center+rho v, multiplying by the displayed powers
of rho_* removes the monopole singularity. For a fixed trim the resulting
expressions extend to rho=0, with v separated from the axis. The rest of
the curve is compact (coordinates are bounded by C/sigma). Their suprema
are thus finite semialgebraic functions of sigma. The elementary growth
argument below bounds each by a power. This establishes finite powers,
without assuming that unspecified smooth-family constants are polynomial.

## 3. Global tube radius and integration

The normal n=(-t_z,t_r), t=B/|B|, is algebraic on the curve. Consider the
normal map x+y n(x) for |y|<epsilon rho_*(x), including the curve inside
the mouth balls down to, but excluding, its limiting centers. Define an
admissible epsilon by the following requirements for every u in [sigma,u0]:
the map is injective, its Jacobian 1-kappa y is at least 1/2, its image
avoids the axis and centers, and distances to each center stay within a
fixed factor of their values on the curve. Replacing epsilon by a smaller
value preserves these requirements.

For a fixed sigma a positive epsilon exists. On a fixed compact middle
portion this is the ordinary embedded-curve tubular argument. In each
center chart the curve is a smooth perturbation of a ray, the scaled frame
converges to its ray frame, and rho times curvature tends to zero. The
ray is separated from the axis uniformly for u in [sigma,u0]. Thus the
normal maps in the two small center cones are injective for one positive
relative width. The cones are disjoint; shrinking their width and that of
the compact middle portion excludes nonlocal intersections. This argument
also gives a continuous positive choice over the compact interval of u.

All the requirements just stated can be written with quantifiers over real
coordinates and polynomial equations/inequalities: describe the center
distances and b by positive square-root auxiliary variables, impose U=u,
and quantify two points and their normal offsets for injectivity. Neither
arclength nor an exponential enters that formula. Restrict epsilon<=1/10
and take half its admissible supremum. Quantifier elimination makes this a
positive semialgebraic function of sigma. Consequently

    epsilon(sigma)>=c sigma^A                                  (7)

for some finite A. Using a fixed smooth transverse bump with width
epsilon(sigma) rho_* costs only additional powers of sigma^-1. Its first
longitudinal derivative uses rho_*', which is bounded. Finite additional
derivatives used in the residual have the same polynomial property.

For completeness the elementary semialgebraic fact is as follows. A
positive semialgebraic function f(t) on (0,t0] has a polynomial lower bound
near zero, and its reciprocal has a polynomial upper bound. A quantifier-free
description of its graph uses finitely many nonzero polynomials. At every
graph point at least one must vanish: otherwise their strict signs persist
in a two-dimensional open neighborhood, contradicting that it is a graph.
Their product P(t,f(t)) therefore vanishes. Write P=sum a_j(t)y^j, divide
out its lowest nonzero y power, and shrink t0 so that its nonzero lowest
coefficient has |a_j(t)|>=c t^m. For 0<f<=1, the remaining bounded
coefficients give c t^m<=C f. For f>1 the lower bound is automatic.
Apply this to 1/f for the upper bound. Projection/quantifier elimination
is the established Tarski-Seidenberg input, not a new result.

The full curve length L_p<=C/sigma: coordinates have magnitude <=C/sigma
by (5), and each coordinate has a uniformly bounded number of monotonic
pieces, because its critical points satisfy fixed-degree algebraic
equations. Total variation bounds arclength by variation(r)+variation(z).
Near a center ds/drho is bounded for each trim, with polynomial dependence
by the same algebraic argument. On the middle portion rho_* is bounded
away from zero and the above length estimate applies. Hence, for small a,

    integral_(rho_*,mouth >= a) ds/rho_*^2 <= C sigma^-A/a.         (8)

Changing the precise collar exclusion by fixed factors changes only C.
This bounds an integral by algebraic quantities; it does not assert that
arclength or the integral itself is semialgebraic.

## 4. Gaussian, overlap and boundary estimates

On the tube, (6)-(7) and the scaled derivatives above give the same Taylor
remainder and Gaussian-moment proof as in joint mouth matching. The
uncut section residual is bounded by C sigma^-A/(sqrt(q) rho_*).
The cutoff annulus has q b y^2>=c q sigma^A. Its polynomial prefactors
can be absorbed into sigma^-A exp(-c q sigma^A), after changing A,c.
Equation (8) consequently gives the previous bulk estimate with
polynomial constants and a weakened exponential. A finite number of
such losses can be bounded using common powers M,N.

The remaining issue is the collar/tube overlap, not addressed by the
earlier local long-return tube. At a fixed small center radius rho, the
peak theta0 of the exact deformed angular ground state solves U=u. Its
distance from either pole is bounded below by a power of sigma (near
rho=0 it is comparable to sqrt(u(1-u))). The derivatives of its logarithmic
phase, divided by q, are algebraic in sin(theta), cos(theta), rho and u
after differentiation. Its peak is nondegenerate; at rho=0 the curvature
is 1/2 in these conventions. Peak curvature, derivatives through the
finite order needed for the Gaussian comparison, and the normal-projection
map have polynomial upper/lower bounds on a neighborhood of polynomial
width in sigma. Choose a<=c sigma^A, which follows by increasing N in
the hypothesis q>=C sigma^-N.

Taylor expansion and normalized Gaussian moments on this neighborhood
give overlap errors

    ||tube-collar|| <= C sigma^-A (q^-1/2+rho^2)
                        + C sigma^-A exp(-c q sigma^A).            (9)

Here the tube contains the cubic and width correctors, and the collar
contains its integral lower-component corrector, as in the original
overlap calculation. There is no extra sqrt(q) in the normalized estimate.
For the tails, the angular phase has its unique maximum at U=u; its
derivative has the sign pointing toward that root. The Taylor decrement
at the boundary of the chosen neighborhood is at least c sigma^A.
Monotonicity then bounds the rest by the same decrement. Powers from the
spinor/normalization factors are harmless after Gaussian normalization.
This argument uses bounds on phase derivatives and its monotonicity;
the exponential ground state and its normalization integral are not
being declared semialgebraic. These facts explain why compactness alone
would not justify uniform Laplace constants and why polynomial bounds do.

The radial collar residual and the loaded higher-mode boundary estimates
are already uniform over all angular sectors in the prior note. Combining
(8)-(9) with those estimates gives

    ||f|| <= C sigma^-M [(qa)^-1/2+a^3/2+sqrt(q)a^7/2
                         +a^-1/2 exp(-c q sigma^N)],
    boundary_loss_cost <= C nu^-1 q^-3/4,                         (10)

for p=q^-3/4. The latter is the squared, loss-weighted boundary norm;
ground basis errors are smaller. A possible polynomial enlargement of
its constant is absorbed in (1). At a=c_a q^-1/4 all three algebraic
bulk terms are O(q^-3/8). Current and trace estimates convert ||f|| to
a return error with factor nu^-1/2, and the square root of the boundary
cost has the same order. The exact powers of nu are not hidden in sigma.

## 5. Source derivatives on a long path

Use the previous full source trial, including P=P0+i w(s)y, not only
the Wilson phase. Its exact Pauli residual identity and cancellation
i(partial_y+qby)(iwy g)=-w g are unchanged. C2 norms in the field frame,
Gaussian moments and path integrals now introduce only polynomial trim
losses. A crude bound using L_p<=C/sigma already suffices for its bulk norm.
The source is zero near both overlap regions, so (9) is unchanged apart
from its longitudinal phase.

For j source directions of bounded operator norms M_j, choose the complex
polydisc |t_j|<=nu/(8 j M_j), omitting zero sources. The perturbation of
the loaded resolvent has norm less than 1/4 in its Neumann-series test;
the dissipative current estimate retains a fixed fraction of nu. The
trial phase grows by at most exp(C nu L_p)<=exp(C nu/sigma), uniformly
bounded under nu<=c sigma. Thus (10) and the return estimate hold on
that polydisc with the same nu power. Cauchy's formula costs nu^-1 for
each source derivative, with the usual repeated-derivative factorial.
This proves (1). A fixed-radius source polydisc would instead allow
exponential growth in 1/sigma and would not prove the claim.

## 6. Compatibility of the chosen physical scaling

Put x=log q. Then sigma is comparable to x^-2, nu=x^-4, and every term
in (1) is q^-3/8 times a finite power of x, plus an exponentially small
term with q/x^(2N) in its exponent. This proves (4), including after
any fixed power of delta^-1 from a later finite response calculation.
It also shows q>=C sigma^-N and nu<=c sigma eventually.

The tempting choice delta=x^-1 is rejected: the prior soft-pair bound
delta log(2+q delta^2) would tend to 1. For delta=x^-2 it tends to zero
like x^-1. The previously required scale tests also remain satisfied:

| Existing component condition | Value/order in this family |
| --- | --- |
| q delta^3 -> infinity | q/x^6 |
| q delta^5 -> infinity | q/x^10 |
| log(2+q delta)/(q delta^2) -> 0 | O(q^-1 x^5) |
| 1/[delta sqrt(q delta)] -> 0 | q^-1/2 x^3 |
| sqrt(r_e ell)/d -> 0 | O(q^-3/2 x^2) |
| Omega_rotation t_d -> 0 | O(q^-3/2 x^5) |
| q Omega_rotation t_d/2 -> 0 | O(q^-1/2 x^5) |
| eta_ang/delta -> 0 | O(q^-1 x^3) |
| (M ell)^2/delta -> 0 | O(q^-1 x^2) |
| open photon weighted error -> 0 | O(q^-1 x^6) |
| spherical reference relative error -> 0 | O(q^-1/4 x^3) |
| singlet path correction delta log q -> 0 | x^-1 |

These are checks against the established component estimates in the
weak-coupling, long-return and spherical-response notes, not an assertion
that their conditional inputs now constitute a physical stress theorem.
In particular a general photon trace need not have only polynomial
feedback losses and requires its own proof.

## 7. Verification, sources and remaining gap

Run `checks/logarithmic_long_return.py --output
evidence/logarithmic-long-return.json` with the recorded Python environment.
It checks (5)-(6) symbolically, compares a cancellation-free exact dipole
flux formula with 70-digit arithmetic, and tests both global bulk residual
integrals and two overlap annuli at q=2^20,2^24,2^28 and n=2,3,5 using
physical half-integer sectors. Two longitudinal/normal quadrature orders
and geometry steps are compared. Bulk attenuation is omitted to give a
larger norm. The numerical frequency is purely imaginary; the bounded
real-frequency and source-derivative assertions use the analytic proof
above and the prior exact residual checks. No finite sample proves the
polynomial-uniform theorem or determines its constants.

The sampled bulk residuals are large (squared norms about 1.2e4 to 2.3e5),
and do not show decreasing return error. The far-field Gaussian is wider
than the narrow, unit-capped tube used by this verification trial. Direct
integration over the supported interval, split at the taper junctions,
replaces the old fixed Hermite rule, which could have no nodes inside this
support. Relative bulk refinement is below 0.00125. The largest overlap
norm discrepancy is below 0.000311, with refinement below 1e-5. These are
checks of formula evaluation and integration, not evidence that the
unspecified asymptotic q threshold has been reached. An explicit useful
finite-flux window remains unestablished.

The standard real-algebraic input is Michel Coste, *An Introduction to
Semialgebraic Geometry* (October 2002), Theorems 2.3 and 2.6 and Proposition
2.11. The authored text is available through this
[academic mirror](https://blogs.mat.ucm.es/jesusr/wp-content/uploads/sites/52/2020/03/SAG.pdf).
The source receipt records provenance, exact pages read and dependency
limits. The root-coefficient argument above independently explains the
one-variable growth step. This is established machinery; applying it to
this matching construction does not establish publication novelty.

This removes the fixed-trim restriction for a fixed finite set of resonant
long paths in an admissible scale family, subject to the loaded assumptions.
It does not control the discarded physical tail, general angular photon
mixing, the exterior determinant, a renormalized closed-system stress or
finite-disturbance stability. Those are the next physical questions.
No Lean or Lake process was invoked; suitable endpoint certification
remains exclusively on GCP after a physical endpoint survives review.

Remaining to-do list: construct/justify the physical loaded boundary and
curved mouth interface; obtain the renormalized state-dependent response,
including angular mixing, exterior determinant and discarded tail; complete
the literature/proof/novelty review and consequential endpoint; independently
verify it, certify suitable mathematics on GCP, and prepare the manuscript.
