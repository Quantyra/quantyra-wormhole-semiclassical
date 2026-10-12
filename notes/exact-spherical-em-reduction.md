# Exact classical spherical reduction and its JT throat limit

2026-10-12, continuation of the S051 radial-constraint gate. This derives
the exact two-dimensional classical Einstein-Maxwell equations at fixed
magnetic charge before making a near-product approximation. It identifies
the Hamiltonian and momentum constraints covariantly and reproduces the
leading JT equation used in `l0-throat-constraint-audit.md`. It does not
reduce the quantum Dirac determinant or solve the sourced dynamical problem.

## Four-dimensional model and reduced action

Use mostly-plus signature and

```text
I_EM = integral d^4x sqrt(-g4) [ R4/(2 kappa) - F^2/(4 g_c^2) ],
kappa=8 pi G,    F=(q/2) sin(theta) dtheta wedge dphi,
ds4^2=g_ab(x) dx^a dx^b + R(x)^2 dOmega_2^2.
```

The magnetic flux `q` is fixed. Set `X=R^2`; the warped-product curvature is
`R4=R2+2/R^2-4 Box(R)/R-2 (grad R)^2/R^2`. Integrating over the sphere and
integrating the `Box(R)` term by parts (with the corresponding boundary
term retained or fixed separately) gives

```text
I_EM = 1/(4G) integral d^2x sqrt(-g) [ X R2
              +(grad X)^2/(2X)+2-2 r_e^2/X ],
r_e^2 = pi G q^2/g_c^2.                                      (1)
```

The coefficient follows directly from `F^2=q^2/(2R^4)` and agrees with the
extremal radius in the MMP normalization. The result is exact within the
spherically symmetric metric and fixed magnetic-flux ansatz; it is not yet
the full quantum-effective action. A spherically symmetric electric field
is a two-dimensional Maxwell constraint sector: its flux is constant in
the absence of charged sources, and at fixed zero electric charge it adds
no local propagating degree of freedom. If electric charge or charged
radial currents are allowed, that constraint and the matter current must be
included rather than silently set to zero.

## Exact equations and constraints

Write `V(X)=2-2 r_e^2/X`. Varying (1) with respect to `g^{ab}` and `X`
gives

```text
E_ab := g_ab Box X - nabla_a nabla_b X
       + (nabla_a X nabla_b X - g_ab (nabla X)^2/2)/(2X)
       - g_ab V(X)/2 = 2G T_ab^(2),                         (2)

R2 - Box X/X + (nabla X)^2/(2X^2) + V'(X) = J_X.          (3)
```

Here `T_ab^(2)=-2/sqrt(-g) delta I_matter^(2)/delta g^ab` and
`J_X=4G/sqrt(-g) delta I_matter^(2)/delta X`; these definitions fix the
source normalizations. In an ADM decomposition
`ds2=-N^2dt^2+Lambda^2(dr+N^r dt)^2`, the projections
`u^a u^b(E_ab-2G T_ab^(2))=0` and
`u^a s^b(E_ab-2G T_ab^(2))=0`, for unit normal `u` and unit radial vector
`s`, are respectively the lapse/Hamiltonian and shift/momentum constraints.
Equation (2) supplies both constraints and the remaining evolution
component; it must not be replaced by only a radial equation for `X`.
The dilaton equation (3) is the angular Einstein equation. Diffeomorphism
invariance relates the divergence of (2) to (3) and to the matter identity;
the reduced quantum stress, dilaton source, and any electric current must
obey that identity together.

## Product throat and the leading JT constraint

The magnetic product has `X_0=r_e^2`, `R2_0=-2/r_e^2`, and `V(X_0)=0`.
Write `X=r_e^2(1+Phi)` and retain first order on the fixed background.
The derivative-squared term in (2) is second order, while
`V'(X_0) delta X/2=Phi`. Dividing (2) by `r_e^2` therefore yields

```text
g_ab Box Phi - nabla_a nabla_b Phi - g_ab Phi
       = (2G/r_e^2) delta T_ab^(2) + higher-order/background terms. (4)
```

Thus the homogeneous operator and its sign in the local constraint audit
follow from the exact classical reduction, rather than an assumed JT
equation. For the prescribed leading LLL static source, the radial
projection of (4) is MMP's `rho Phi'-Phi=-alpha/(1+rho^2)` and its solution
`Phi=alpha(1+rho atan rho)`. The previous deterministic checker verifies
the homogeneous component identities and this static radial equation;
the separate [full-tensor reference calculation](casimir-reference-matching.md)
checks all components for its static first-order source family.

The deterministic check `checks/spherical_reduction_consistency.ps1`
evaluates the exact metric equation (2) on arbitrary first/second-derivative
jets of `X` at 15 throat points and compares all three coordinate
components with (4). It checks 45 components at `epsilon=10^-5`; the maximum
residual divided by `epsilon^2` is `10.56`, consistent with the discarded
quadratic terms. This is an algebra check of the written reduced equation
and its linearization, not an independent derivation of the four-dimensional
quantum source.

The limit in (4) freezes the background metric, neglects quadratic
perturbations, trace-proportional reference shifts, and higher angular
fermion modes, and treats `delta T_ab^(2)` as a source. For the actual
Einstein-Maxwell-Dirac problem this is not sufficient: the quantum source
also varies with `g_ab`, `X`, the electric potential and the state. Its
variation contributes to the linear operator and must be calculated from
the two-ended state, with local counterterms and the anomaly fixed. The
classical reduction determines the constraints but does not decide whether
the proposed boundary clock survives matching as a physical degree of
freedom.

MMP's [section 5.3.2](https://arxiv.org/html/1807.04726v3#S5.SS3) derives
the static spherical Einstein-Maxwell equations and matches the throat to
the near-extremal exterior. Equation (1) is their underlying exact
spherically reduced classical sector in an area-radius variable, written
before their throat expansion; equation (4) recovers the already-known
leading limit. No novelty is claimed for this reduction.

Remaining to-do list: derive the spherical reduction of the Dirac action and
its renormalized state-dependent stress, dilaton source and current; solve
the coupled linearized lapse/shift/area constraints with both mouth
matching conditions; then determine whether a physical clock mode remains.
