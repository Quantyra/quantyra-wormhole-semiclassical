# Dirac mode reduction and source terms for the spherical constraints

2026-10-12, continuation of the radial-constraint gate. This reduces the
classical four-dimensional Dirac operator on a spherically symmetric
fixed-magnetic-flux background to two-dimensional angular modes and writes
the stress, area-dilaton source and electric current that enter the exact
spherical Einstein-Maxwell equations. It does not evaluate those operators
in the MMP two-ended state.

## Geometry and canonical mode fields

Use the fixed magnetic flux and general spherical metric

```text
ds4^2=g_ab(x) dx^a dx^b + R(x)^2 dOmega_2^2,
F_theta phi=(q/2) sin(theta),       X=R^2.
```

Let `D_S,q` denote the self-adjoint unit-sphere Dirac operator twisted by
the monopole connection. In a warped orthonormal frame its four-dimensional
Dirac operator has a two-dimensional derivative term, an additional
`gamma^a partial_a(log R)` spin-connection term, and the angular operator
`R^-1 gamma_* D_S,q`. Rescaling the four-dimensional spinor as
`Psi=R^-1 chi` cancels the extra radial spin-connection term in the
quadratic action. Expanding `chi` in an orthonormal eigenbasis of `D_S,q`
then gives canonical two-dimensional fields.

For unit charge the angular spectrum is `q` zero modes of one sphere
chirality and, for each `n>=1`, eigenvalues
`+/-lambda_n`, `lambda_n=sqrt(n(n+q))`, with multiplicity `q+2n` for
each sign. Pairing the signs by a constant two-dimensional unitary
transformation gives the reduced action

```text
I_D,red = sum_(A=1)^q I_2[chi_(0,A); g_ab,A_a]
        + sum_(n>=1) sum_(A=1)^(2(q+2n))
              I_2[chi_(n,A); g_ab,A_a,m_n(X)],
m_n(X)=lambda_n/sqrt(X),
I_2[chi;g,A,m]=integral sqrt(-g) chi_bar
                       [i gamma^a( nabla_a-i A_a)-m] chi.     (1)
```

The displayed Lorentzian `i` convention can equivalently be shifted into
the gamma-matrix convention; the invariant statements are the mass
magnitude and multiplicities. The `n=0` fields are massless two-dimensional
Dirac fields in the vectorlike four-dimensional Dirac theory. For a single
four-dimensional Weyl field, the same zero modes are chiral and their
chirality changes with flux sign. This agrees with MMP's explicit LLL
reduction and prevents a silent factor-of-two change between a Weyl and a
Dirac model.

For the conformal two-metric `g_ab=N(x)^2(-dt^2+dx^2)`, the covariant mass
`m_n=lambda_n/R` becomes the optical mass
`m_n^opt=N lambda_n/R` after the standard two-dimensional Weyl rescaling.
That is the mass used in the mode-transfer calculations. It is not a
constant `AdS2` mass on a general spherical perturbation.

## Sources for the dilaton and Maxwell equations

Let `Gamma_D,ren[g,X,A;state]` be the renormalized in-in or Euclidean
effective action for a specified state, including the complete angular
tower, local counterterms and the Jacobian of `Psi=R^-1 chi`. Use the
source conventions of `exact-spherical-em-reduction.md`:

```text
T_ab^(2) = -2/sqrt(-g) delta Gamma_D,ren/delta g^ab,
J_X      =  4G/sqrt(-g) delta Gamma_D,ren/delta X,
j^a      =  1/sqrt(-g) delta Gamma_D,ren/delta A_a.          (2)
```

These are the quantum sources in the metric and angular Einstein equations;
`j^a` sources the two-dimensional electric Maxwell constraint. If the
spherical electric potential is retained, its reduced Maxwell action is
`I_M,2=-pi/g_c^2 integral sqrt(-g) X F_ab F^ab`. It has no local propagating
mode in two dimensions, but its Gauss constraint must be solved with the
charge density in `j^a`. With the current convention in (2), its equation
is `4 pi/g_c^2 nabla_a(X F^ab)+j^b=0`. The classical reduced Einstein
equations are therefore (2)-(3) of the exact spherical reduction note
with these sources.

At the formal mode-action level, `T_ab^(2)` is the sum of the two-dimensional
stress tensors of all modes plus local renormalization terms. The explicit
mass dependence contributes to `J_X`: since
`dm_n/dX=-m_n/(2X)`, each massive mode contributes a renormalized scalar
bilinear proportional to `m_n <chi_bar chi>/(2X)`. The LLL action (1) has
no classical explicit `X` dependence after the canonical rescaling, but
the regulated four-dimensional determinant can still generate local
dilaton terms from the rescaling Jacobian and from the massive tower.
Consequently setting the LLL nonlocal contribution to `J_X` to zero is
not a valid statement about the complete renormalized four-dimensional
theory without evaluating those local terms in the same scheme.

The mode determinants are not separately UV finite: their large-`n`
contribution must be combined and subtracted as a four-dimensional
determinant. The homogeneous product calculation already shows that the
full Dirac trace/anomaly shifts the reference radii at order `g_c^2`,
which can precede the retained `g_c^2/q` LLL flare at large flux. A
constraint calculation that retains only the LLL stress but omits this
reference correction is therefore a specified approximation, not the
complete one-loop spherical equation.

## What remains to compute for the selected stability target

Equation (1) identifies the field content and the exact area dependence
of the angular masses, but it does not determine any expectation value in
(2). For the MMP stability problem those values must be evaluated in the
two-ended state with the actual mouth conditions, magnetic holonomy,
boundary time normalization, and the same renormalization prescription
used for the reference geometry. Linear response requires the coupled
variations `delta T_ab^(2)`, `delta J_X`, and `delta j^a` with respect to
`delta g_ab`, `delta X`, `delta A_a`, and the state. They must obey the
reduced diffeomorphism and gauge Ward identities together; a prescribed
stress history need not satisfy these identities.

This is the matter-side input that the classical lapse and shift
constraints need. Until it is computed and matched through both mouths,
the current free proper-mouth Gaussian kernel remains a response to a
prescribed metric family, not a physical radial stability operator.

MMP [section 2.3](https://arxiv.org/html/1807.04726v3#S2.SS3) gives the
`q` massless angular zero modes and their two-dimensional propagation;
Appendix C gives the separated Dirac equation on a general conformal
two-metric and varying sphere radius. Its [section 5.3.2](https://arxiv.org/html/1807.04726v3#S5.SS3)
uses the corresponding static LLL stress in the spherical Einstein
equations. The full magnetic spectrum and determinant multiplicities used
here are also recorded in [the homogeneous product calculation](magnetic-product-vacuum.md).
These are reproduced source inputs, not a novel angular reduction.

Remaining to-do list: evaluate the full renormalized source variations in
the two-ended MMP state; impose current conservation and both gravitational
constraints; match the resulting radial solutions at both asymptotically
flat mouths; then reassess primary/fallback novelty and the physical
stability claim.
