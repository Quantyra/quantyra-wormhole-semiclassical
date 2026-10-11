# Trace matching before physical wormhole stability

2026-10-10. S051/S052. This is a necessary background condition for a
specified one-Dirac-loop theory, not a wormhole impossibility theorem.

The traceless lowest-Landau-level (LLL) throat lift cannot by itself supply
the complete renormalized four-dimensional Dirac stress. On a homogeneous
magnetic throat the missing trace requires an unequal-radius correction.
Its leading fractional size is of order `g_c^2`, whereas the retained MMP
Casimir backreaction is of order `g_c^2/q`. The next physical calculation
must therefore specify a quantum-corrected reference background before
assigning full-theory significance to the alpha-order perturbation result.
The earlier compact-minus-line Casimir calculations remain valid within
their prescribed-background scope.

## Model, normalization and sources

Use mostly-plus signature and one free, massless, unit-charge four-component
Dirac determinant with `D_mu = nabla_mu - i A_mu`. The classical action is

```text
I_cl = integral sqrt(-g) [R4/(2 kappa) - F^2/(4 g_c^2)],
kappa = 8 pi G,       Lambda = 0,
A_phi = q cos(theta)/2,       q an integer >= 1.
```

Only the Dirac field is integrated out. There are no photon or graviton
loops in this model. Renormalized couplings and the usual local
curvature-squared counterterms must be fixed in any complete stress
calculation. The bulk trace test below is independent of their finite
coefficients on a constant-curvature product. Extra matter, a cosmological
constant or higher-dimension operators would change the specified problem.

For this stress convention the state-independent, non-derivative anomaly is

```text
A = <T^mu_mu>
  = [C_Weyl^2/20 - 11 E4/360]/(16 pi^2) - F^2/(24 pi^2)
    + a scheme-dependent coefficient times box(R4),
E4 = Riem^2 - 4 Ric^2 + R4^2.
```

The negative gauge sign is important. [Bastianelli and Broccoli,
v2, section 4.2](https://arxiv.org/html/1808.03489v2#S4.SS2), equations
(81) and (91), obtain it with two regulators; their Appendix A uses
mostly-plus signature. Set their axial gauge field to zero. [Shukla,
v3, section III](https://arxiv.org/html/1906.02334v3#S3), equations (49)-(51),
gives both curvature coefficients and the same gauge sign with the stress
definition in equation (4). Equation (49) and its coefficients were also
checked in v2 and on page 7 of the published article,
[Physical Review D 100, 096010](https://doi.org/10.1103/PhysRevD.100.096010).
These are established anomaly coefficients, not a new result of this project.

## Independent contraction of the local coefficient

As a check separate from substituting the published anomaly, take Euclidean
orthonormal sectional curvatures `-x,y`, with magnetic `F_23=B`. The
four-component Dirac Laplace operator has endomorphism and curvature

```text
E = -R4 I/4 - (i/2) gamma^i gamma^j F_ij,
Omega_ij = (1/4) R_ijkl gamma^k gamma^l - i F_ij I.
```

Explicit Clifford matrices give

```text
tr E = -R4,
tr E^2 = R4^2/4 + 4 B^2,
tr Omega_ij Omega^ij = -Riem^2/2 - 8 B^2.
```

Insert these into the universal local heat-kernel coefficient

```text
a4 = tr[60 R4 E + 180 E^2 + 30 Omega^2
        + (5 R4^2 - 2 Ric^2 + 2 Riem^2) I] / [360 (4 pi)^2],
```

with derivative terms absent on the product. The fermion determinant has
anomaly `-a4`. This yields

```text
A = (7 Riem^2 + 8 Ric^2 - 5 R4^2)/(5760 pi^2) - F^2/(24 pi^2),
```

which equals the expression above. The Laplace-type coefficient, its
derivation and its relation to the anomaly are imported from
[Vassilevich, v3](https://arxiv.org/html/hep-th/0306138v3), sections 3.3,
4.1 and 7.1, especially (4.28); the fermionic sign is also checked against
Bastianelli and Broccoli's regulator derivation. The verification script
does not prove the general heat-kernel expansion theorem.

## Product trace and required radius ratio

Consider the local product

```text
ds^2 = a^2[-(1+rho^2) dtau^2 + drho^2/(1+rho^2)] + b^2 dOmega^2,
x = 1/a^2,       y = 1/b^2.
```

Direct Christoffel and Riemann contractions give

```text
R4 = 2(y-x),          Ric^2 = 2(x^2+y^2),
Riem^2 = 4(x^2+y^2),  E4 = -8xy,
C_Weyl^2 = (4/3)(y-x)^2,       F^2 = q^2 y^2/2.
```

Consequently

```text
A(a,b,q) = [3x^2 + 5xy + (3-15q^2)y^2]/(720 pi^2).
```

At equal radii `a=b=R`, the trace is

```text
A(R,R,q) = (11-15q^2)/(720 pi^2 R^4) < 0    for integer q >= 1.
```

Classical Maxwell stress is traceless. Thus the trace of the semiclassical
Einstein equation, `-R4 = kappa A`, requires

```text
2(x-y) = kappa [3x^2 + 5xy + (3-15q^2)y^2]/(720 pi^2).
```

The equal-radius Bertotti-Robinson product fails this equation at the first
Dirac-loop order for every nonzero integral magnetic flux. Changing the
state cannot remove this local massless trace. This excludes that exact
uncorrected product as a solution of the specified model, not corrected
products, global wormholes or other theories.

For a perturbatively nearby homogeneous reference, set `delta=a^2/b^2-1`.
Expand with an explicit loop parameter multiplying `A`, taking delta to be
of the same order. The trace equation fixes the necessary first-order term

```text
delta = kappa (15q^2-11)/(1440 pi^2 b^2) + higher backreaction orders.
```

Use the classical radius relation only inside this first-order term:

```text
r_e^2 = kappa q^2/(8 g_c^2),
delta = g_c^2/(12 pi^2) [1-11/(15q^2)] + higher orders.
```

This determines one combination of the two radii. It neither supplies
their common scale nor the two independent quantum pressures. In
particular, it is not an existence proof or a prescription for inserting
these homogeneous radii into an inhomogeneous Casimir wormhole.

For the existing throat notation,

```text
alpha = g_c^2/(4 pi^2 q),       s0 = g_c^2 q/(4 pi^2),
delta_leading/alpha = q/3 - 11/(45q),
delta_leading = s0/(3q) [1-11/(15q^2)],       alpha = s0/q^2.
```

The ratio is about 0.919, 2.636, 21.330 and 85.332 at `q=3,8,64,256`.
Asymptotically, vacuum matching precedes the retained Casimir backreaction
even at fixed `s0`. Making the coupling weak reduces both corrections; it
does not reverse this relative ordering. These are comparisons of
perturbative orders, not a claim that weak-coupling corrections are large
in absolute size.

## What the W22 truncation does and does not include

For the [W22 baseline](https://arxiv.org/html/2610.09847v1), section 2.1,
write the metric in `r_e` units with

```text
S = f + alpha gamma,       C = 1 + alpha phi,       f = 1+rho^2,
phi = rho atan(rho) + 4/3,
gamma = -4rho^2/3 - rho(rho^2+3)atan(rho) + log(f).
```

The scalar curvature of `diag(-S,1/S,C,C sin(theta)^2)` is

```text
R4 = -S'' + 2/C - 2 S C''/C - 2 S' C'/C + S (C')^2/(2C^2).
```

Both its zeroth-order term and its alpha coefficient vanish. The latter
identity is `-gamma''-2phi-2f phi''-2f'phi'=0`.
Likewise W22's displayed section 6.2 quantum stress,

```text
T_Q = q/(96 pi^2 r_e^2)
      diag(rho^2-2, -(rho^2+4)/f^2, 1, sin(theta)^2),
```

has zero four-dimensional trace on the zeroth-order throat. More generally,
the section 6.1 lift defines its sphere pressure so that the total trace
cancels. A two-dimensional anomaly in the reduced tensor does not supply
the four-dimensional anomaly of all Dirac modes.

This identifies an additional matching requirement when using the LLL
model for full one-Dirac-loop accuracy. It does not invalidate W22's
expressly restricted LLL/local/gradient calculation. [MMP,
v3](https://arxiv.org/html/1807.04726v3), section 5.2, explicitly allows
other four-dimensional modes to contribute to the trace; its weak-coupling
and magnetic-scale renormalization assumptions must be preserved. Neither
paper is being recast here as claiming an exact all-mode quantum background.

## Renormalization and prior-work limits

A fixed finite Maxwell counterterm has zero metric trace in four dimensions.
Finite curvature-squared terms shift the bulk trace by a multiple of
`box(R4)`, which is zero on this homogeneous product. These counterterms
therefore cannot cancel the displayed trace while holding the model fixed.
Their role in the remaining stress equations is a separate calculation.
Changing the cosmological constant or adding supporting fields is a change
of physical inputs, not a subtraction-scheme resolution of this test.

Evaluating the running coupling at a magnetic scale does not remove the
gauge trace either. At one Dirac loop,

```text
d(1/g_c^2)/d log(mu) = -1/(6 pi^2),       mu proportional to sqrt(B).
```

Under a constant Weyl scaling, `B -> exp(-2sigma) B`; hence
`delta(1/g_c^2)=sigma/(6 pi^2)`. Varying the complete Maxwell density
`-F^2/[4g_c^2(sqrt(B))]` gives `-F^2/(24 pi^2)` again. Holding the inserted
number `g_c(B)` fixed during metric variation would miss this contribution.
This running-coupling argument checks the leading magnetic term, not the
complete curved-space determinant.

Quantum corrections to extremal reference geometries have extensive prior
work. [Matyjasek and Zaslavskii, v2, section IV.1](https://arxiv.org/html/gr-qc/0102109v2#S4.SS1)
solve a massive-field Schwinger-DeWitt problem. At zero cosmological
constant their retained perturbative branch keeps equal radii and changes
the radius-charge relation. That result is not the massless charged model
used here; sending the mass to zero leaves their approximation regime.
[Abe, Noumi and Yoshimura, v2](https://arxiv.org/html/2305.17062v2),
sections 2.2.1 and 3.1, derive magnetic extremality equations with the
Euler-Heisenberg action. Their stated slowly varying massive-field
condition also cannot be imported by setting the mass to zero. These
comparisons preclude a broad claim to first quantum-corrected extremality.

## Verification and research decision

[The script](../checks/throat_trace_anomaly.py) independently contracts the
metric curvature and Clifford matrices, reproduces the anomaly, verifies
the W22 trace cancellations, and checks the first-order ratio and running
coupling variation. [The receipt](../evidence/throat-trace-anomaly.json)
contains 20 exact-zero identities, one positive shifted polynomial and
four exact rational order comparisons. The universal heat-kernel formula
and quantum anomaly framework remain literature inputs. There is no Lean
certificate and no physical stability or instability claim.

Complete the full renormalized vacuum stress on the magnetic curved
reference, with a state and renormalization conditions, before extending
the stability calculation. A useful next feasibility test is whether a
spectral determinant on the homogeneous product supplies both pressures
and a controlled expansion; its spectral normalization, boundary state
and closest prior derivations must first be checked. Then add the
state-dependent wormhole Casimir term and determine the total response.
The anomaly cancels in compact-minus-line differences with common local
subtraction, so the earlier finite-size kernel cannot supply the missing
reference stress by itself. Formalizing the present necessary condition
would not substitute for this unresolved physical question.

Remaining to-do list: complete the source and novelty review, compute the
quantum reference stress and its matched background, then establish and
verify a consequential total-model stability, transport or obstruction result.
