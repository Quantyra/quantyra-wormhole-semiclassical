# Leading MMP monopole throat constraints

2026-10-12, S053 primary-target gate. The prescribed proper-mouth metric
kernel is not itself a gravitational normal mode. This note derives the
leading spherically symmetric dilaton constraint on the MMP
`AdS_2 x S^2` throat and reproduces the static radial equation used to
match the throat to the mouths. It identifies the homogeneous throat
solution space, but does not solve the sourced dynamical problem or impose
global mouth boundary conditions.

## Leading throat equation and components

In MMP's long-throat expansion, write the unit-radius `AdS_2` metric as

```text
ds_2^2=-f(rho) d tau^2 + d rho^2/f(rho),
f(rho)=1+rho^2.
```

Let `Phi(tau,rho)` be the first-order fractional perturbation of the
two-sphere area factor in `ds_4^2=r_e^2[ds_2^2+(1+Phi)dOmega_2^2]`.
At leading order in the near-product/LLL approximation, fixed magnetic
flux and after projection onto the traceless two-dimensional LLL source,
the spherical Einstein equation takes the JT constraint form

```text
L_ab[Phi] = S_ab,
L_ab[Phi]=g_ab Box Phi - nabla_a nabla_b Phi - g_ab Phi.       (1)
```

`S_ab` denotes the source in MMP's dimensionless throat normalization;
trace-proportional four-dimensional vacuum terms, higher Landau levels,
and corrections to the product background are outside (1). The overall
sign and normalization are fixed by MMP equation (5.36), rather than
chosen independently here. The three independent coordinate components
of the left side, with the `rho rho` and `tau tau` equations rescaled for
clarity, are

```text
f L_rhorho = -Phi_tautau/f + rho Phi_rho - Phi,
L_taurho   = -Phi_taurho + rho Phi_tau/f,
L_tautau/f = -f Phi_rhorho - rho Phi_rho + Phi.                (2)
```

For a static `Phi`, the first equation reduces to `rho Phi'-Phi`. MMP's
static LLL source gives

```text
rho Phi' - Phi = -alpha/f,
Phi(rho)=alpha(1+rho atan(rho)),                              (3)
```

which is their equations (5.36)-(5.37). Direct substitution verifies
equation (3); matching the large-`rho` coefficient to the asymptotically
flat near-extremal Reissner-Nordstrom mouth then fixes the throat length
in equations (5.39)-(5.40). This reproduces the relevant static result;
it does not add a new static solution.

## Homogeneous solutions are constrained, not an arbitrary clock history

For a source-free perturbation, tracing (1) gives `Box Phi=2 Phi`, so
`nabla_a nabla_b Phi=g_ab Phi`. On connected `AdS_2`, such a solution is
fixed by its value and two first derivatives at one point, hence the
solution space has dimension at most three. The three independent
solutions are

```text
Phi_0=rho,
Phi_c=sqrt(f) cos(tau),
Phi_s=sqrt(f) sin(tau).                                      (4)
```

Substitution into all three equations (2) verifies that they span the
homogeneous leading-throat solution space. They are the constrained
embedding/dilaton modes of the product throat, not an arbitrary function
`phi(t)` multiplying a prescribed metric profile. Whether a combination
survives as a physical boundary mode depends on mouth matching, boundary
time normalization, and the Schwarzian/boundary action. This local result
does not declare the modes pure gauge.

This interpretation is consistent with the September 2026
[Sadhukhan perturbation paper](https://arxiv.org/html/2609.33511v1), which
states that the `l=0` sector contains the JT dilaton controlling throat
length and is “not trivial” despite not propagating. In its exact
`alpha=0` throat limit, its equation is `Box K=2K`, and it identifies the
static MMP sphere-radius profile among its solutions. The paper restricts
its master equations and `O(alpha)` stability result to `l>=2`; it does not
solve the `l=0` boundary-value problem. The present three-dimensional local
solution space agrees with the equation it displays, but the asymptotic
clock and two-mouth classification remain open.

For a time-dependent matter source, (2) must be solved simultaneously,
with the stress/current response evaluated in that same perturbed state.
The source must satisfy `nabla^a S_ab=0`, as follows from the divergence
of (1) on `AdS_2`; otherwise the linearized Bianchi identity fails. The
MMP static solution only checks this system at zero frequency. The
prescribed Gaussian proper-mouth response has not yet supplied the
constraint-solved `S_ab(tau,rho)` or the matching data.

## Verification and scope

Run `checks/l0_throat_constraint.ps1`. It substitutes the three candidate
homogeneous modes into each of the three differential expressions in (2)
on a deterministic coordinate grid and checks the static identity (3) at
multiple radii. Its receipt binds this note and script. These checks
verify the stated local formulas; they do not verify the 4D reduction,
quantum source conservation, coupled state response, or global matching.

The closest-source comparison is MMP [section 5.3.2](https://arxiv.org/html/1807.04726v3#S5.SS3), equations (5.34)-(5.40): a static spherical perturbation sourced by their specified LLL stress, matched to the near-extremal exterior. Their section 5.4 discusses the Schwarzian validity range and the energy scale at which horizons form, not a time-dependent linear monopole eigenmode. The September 2026 [Sadhukhan perturbation analysis](https://arxiv.org/html/2609.33511v1) is more direct: it explicitly says the `l=0` throat contains the nontrivial, nonpropagating JT dilaton but leaves it untreated; in its exact `alpha=0` limit its equation reduces to `Box K=2K`, whose static solutions include MMP's sphere-radius profile. Its perturbation master equations and first-order stability result are restricted to `l>=2`. That paper already derives the LLL stress response and induced Hall current required by the Bianchi identities in that higher-multipole sector. Neither its `l>=2` stability theorem nor its linear-response formulas establish the constrained `l=0` two-mouth boundary problem, and the higher-multipole equations cannot simply be continued to `l=0` because the harmonic gauge fixing degenerates there.

Remaining to-do list: derive the `l=0` gauge-invariant boundary data and solve the constraints through both mouths; incorporate the matched LLL response and induced current in the state/boundary convention used by Sadhukhan; include massive-mode and full-anomaly corrections; only then interpret a clock kernel as a physical stability operator.
