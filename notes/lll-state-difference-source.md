# Conserved LLL state-difference source in the spherical throat

2026-10-12, continuation of the monopole constraint reduction. This
isolates the nonlocal, state-dependent lowest-Landau-level (LLL) source on
the fixed global `AdS2 x S2` throat, proves its conservation and zero
area-dilaton source within the fixed-boundary LLL model, and reproduces the
static MMP source profile. It does not calculate the full four-dimensional
renormalized stress or its response when the mouth boundary conditions move.

## State difference and its source dependence

Take `ds4^2=g_ab dx^a dx^b+X dOmega_2^2`, fixed monopole flux `q`, and
restrict the free charged Dirac field to its `q` sphere zero modes. The
canonical rescaling `Psi=X^(-1/2) chi` cancels the radial spin-connection
term. Because the angular eigenvalue is zero, the resulting two-dimensional
LLL action depends on `g_ab` and the two-dimensional gauge field, but has
no explicit `X` dependence. Thus for two LLL states `s_1,s_0` on the same
fixed two-dimensional region and with identical boundary conditions,

```text
Delta Gamma_LLL[g,A] = Gamma_LLL[g,A;s_1]-Gamma_LLL[g,A;s_0],
delta(Delta Gamma_LLL)/delta X = 0.                         (1)
```

The common local determinant Jacobian and four-dimensional counterterms
cancel in this state difference. Equation (1) is only about the finite
state-dependent LLL difference. It says nothing about the absolute
four-dimensional vacuum source: the full angular tower, spinor-rescaling
Jacobian and trace anomaly still contribute to the reference geometry.
Nor does (1) hold unchanged if varying `X` changes the domain, boundary
conditions, magnetic holonomy, or which state is being compared.

With vanishing two-dimensional electric field and no state-dependent
current, diffeomorphism invariance of this two-dimensional state difference
gives

```text
nabla^a Delta T_ab = 0,       g^ab Delta T_ab=0.             (2)
```

The trace vanishes because the two states have the same local anomaly; their
difference is traceless. A nonzero current instead adds the Lorentz-force
term in the Ward identity and must be accompanied by the Maxwell constraint.

## MMP global-throat value

Use `x=atan(rho)` and `f=1+rho^2`, so
`ds2^2=r_e^2 sec^2(x)(-d tau^2+dx^2)` on `-pi/2<x<pi/2`. For `q` complex
massless LLL channels, compare the antiperiodic ground state on the full
circle of circumference `pi` with the reflecting-strip reference. Their
finite Casimir stress difference in conformal coordinates is

```text
Delta T_tau tau = Delta T_xx = -q/(8 pi),
Delta T_tau x = 0.                                          (3)
```

The coefficient is the difference between circle and reflecting-strip
vacuum energies divided by the common length; the local Weyl-anomaly stress
is common and cancels. Transforming the covariant tensor to `rho` gives

```text
Delta T_tau tau = -q/(8 pi),
Delta T_rho rho = -q/(8 pi f^2),
Delta T_tau rho = 0.                                        (4)
```

Direct substitution into the two-dimensional metric verifies both
identities in (2). Substituting (4) into the spherical metric constraint
reproduces MMP's static `rho Phi'-Phi=-alpha/f` and
`Phi=alpha(1+rho atan rho)` after matching their normalization. Therefore
the familiar static throat opening is precisely sourced by a conserved
LLL state difference with zero direct `X` variation, rather than by an
independent angular pressure.

This does not determine the dynamical response needed by the stability
question. On a fixed conformal strip the CFT anomaly fixes the local Weyl
part of the stress response, but the state functions and their variation
when the two mouth clocks/boundary conditions move must be obtained from
the globally matched state. The massive angular tower and absolute local
terms also change when `X` varies. Reusing (1) while varying those data
would omit the corresponding boundary/state response.

Run `checks/lll_state_difference_source.ps1`. It checks the coordinate
transformation, trace and divergence of (3)-(4) over a deterministic
coordinate grid. Its receipt records the script and note hashes. The test
verifies the fixed-throat Ward identities only; it does not test the
four-dimensional anomaly or a varying-boundary state.

MMP [section 2.3](https://arxiv.org/html/1807.04726v3#S2.SS3) derives the
`q` massless two-dimensional modes, and [section 5.3.2](https://arxiv.org/html/1807.04726v3#S5.SS3)
uses the corresponding static quantum stress in the spherical Einstein
equations. The state-difference interpretation and conservation check here
are a reproduction of that fixed-throat source, not a new stability result.

Remaining to-do list: compute the matched LLL stress/state response when
both mouth conditions vary; include the renormalized massive tower and
absolute anomaly/dilaton terms; solve the coupled dynamical constraints and
match to the asymptotically flat exteriors.
