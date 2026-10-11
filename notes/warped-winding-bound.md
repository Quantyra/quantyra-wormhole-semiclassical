# A winding energy and null-stress bound on a warped magnetic circle

2026-10-10. S051/S052 feasibility result. For a prescribed static spherical
circle, the antiperiodic winding part of a free charged Dirac determinant
has non-positive radial null stress. The entire higher-angular-mode sector
is bounded by the integrated throat barrier, even when the local angular
mass becomes small elsewhere on the circle. This removes the constant-radius
assumption of the earlier cylinder test for this observable.

The reference subtraction, global geometry and absence of angular mixing
are essential. This result does not yet bound the complete MMP support
stress, establish its stability, or constitute an accepted novelty verdict.

## Defined model and observable

Let `x` have coordinate period `X`, and prescribe

```text
ds^2 = -N(x)^2 dt^2 + B(x)^2 dx^2 + R(x)^2 dOmega^2,
F = (q/2) sin(theta) dtheta wedge dphi,    integer q >= 1.
```

Take positive, smooth, periodic `N,B,R`, a fixed gauge field, no electric
potential, and one free massless four-dimensional Dirac field in its static
ground state. Piecewise smooth finite coefficients can also be used in the
reduced transfer problem with continuous rescaled-spinor matching; these
are not asserted to be smooth Einstein geometries. No photon/graviton loops or interacting
multi-flavor dynamics are included. Spherical symmetry makes angular modes
independent. The known magnetic Dirac spectrum gives `q` zero modes and,
for `n>=1`, masses `lambda_n=sqrt(n(n+q))` of both signs with degeneracy
`q+2n` per sign. The two signs give identical determinants.

For spatial phase `psi(x+X)=exp(i theta) psi(x)`, define

```text
E_w(theta) = E(theta) - (1/(2 pi)) integral_0^(2 pi) E(vartheta) d vartheta.  (1)
```

Use the identical metric, field and local renormalization prescription in
each term. All local ultraviolet subtractions and anomaly terms cancel.
Equivalently, the average is the Bloch average per cell of the periodic
universal cover. This isolates the winding contribution. It is not the
reflecting global-AdS2 vacuum or the homogeneous product reference used in
the previous notes. We evaluate (1) at `theta=pi`.

Define the optical circumference and integrated angular barrier

```text
L = integral_0^X B/N dx,       J = integral_0^X B/R dx.                    (2)
```

After the usual spinor rescaling, one angular mode is a two-dimensional
Dirac field of optical mass `lambda N/R`. At imaginary frequency `xi>=0`,
a constant spinor basis change gives the real transfer equation

```text
M' = A M,  M(0)=1,
A = a sigma_a + b sigma_b,
a=lambda B/R >= 0, b=xi B/N >= 0,
sigma_a=diag(-1,1), sigma_b=[[0,1],[1,0]].                                (3)
```

Here `lambda` denotes the positive angular mass. Negative masses are
unitarily equivalent. Write `I=lambda J`.

## Transfer estimate and finite determinant

The generator is traceless, so `det M(X)=1`. Expand its evolution around
the diagonal part `a sigma_a`. Every remaining insertion of `b sigma_b`
and every diagonal propagator has non-negative entries. The zeroth term
therefore yields `tr M(X)>=2 cosh(I)`. A constant orthogonal basis change
interchanges the roles of `a` and `b` while retaining non-negative
off-diagonal entries. The same argument gives `tr M(X)>=2 cosh(xi L)`.
These Dyson-series arguments extend by approximation to the stated smooth
coefficients. Thus

```text
zeta(xi) = arcosh(tr M(X)/2) >= max(I, xi L).                              (4)
```

For completeness, the phase-dependent spectral characteristic function is
`det(M(X)-exp(i theta)1)`, proportional to `tr M(X)-2 cos(theta)`.
It determines the twisted Dirac eigenvalues. The argument-principle
vacuum-energy difference, rotated to imaginary frequency, is
`-(1/pi) integral log[(tr M-2cos(theta))/(tr M-2cos(theta_0))] dxi`.
The sign and factor are those of a complex two-component Dirac field:
each positive/negative pair contributes minus its positive frequency.
The spectrum is paired because the mass-only Hamiltonian anticommutes with
a constant Pauli matrix preserving the twisted domain. All operators use
the same smooth coefficients; the high-frequency phase-dependent ratio
tends to one, fixing the normalization without a boundary counterterm.

Factor the characteristic function as
`exp(zeta)(1-2cos(theta)exp(-zeta)+exp(-2zeta))`. Its phase-averaged logarithm
is `zeta`; expansion in Fourier modes proves this for `zeta>0`. Consequently
the antiperiodic winding energy of one angular mode is exactly

```text
E_w(pi) = -(2/pi) integral_0^infinity log(1+exp(-zeta(xi))) dxi.             (5)
```

Equation (4) makes (5) absolutely convergent. Using `log(1+y)<=y` and
splitting the integral at `xi=I/L` gives

```text
0 <= -E_w(pi) <= 2 (I+1) exp(-I)/(pi L).                                  (6)
```

For a zero mode, `a=0`, `zeta=xi L` and `E_w(pi)=-pi/(6L)` exactly.
For constant optical mass `m`, (5) also equals
`(2m/pi) sum_{k>=1} (-1)^k K_1(k m L)/k`. These independent normalizations
are checked numerically. Neither a pointwise exterior mass gap nor a
low-real-frequency scattering approximation was used in (4)-(6).

## Local metric variation and radial null stress

Let `K_a(x)=delta zeta/delta a(x)` and likewise for `K_b`. Duhamel's formula
inserts the corresponding Pauli matrix into the cyclic transfer product.
In the original basis that product has non-negative entries, so
`|tr(sigma_a P)|<=tr P`. In the interchanged basis the same is true for
`sigma_b`. Positivity of the off-diagonal insertion in each respective
basis also fixes the lower signs. Since
`delta zeta=delta(tr M)/(2 sinh zeta)`, this proves, for `I>0`,

```text
0 <= K_a(x), K_b(x) <= coth(zeta) <= coth(I).                              (7)
```

For the angularly averaged stress of one mode, variation of the conserved
coordinate-time energy defines

```text
rho_w = [delta E_w/delta N]/(4 pi B R^2),
p_w   = -[delta E_w/delta B]/(4 pi N R^2).
```

The angular degeneracy restores exactly uniform stress when summed over
each multiplet. Differentiating (5) using (3) gives the radial orthonormal
null combination

```text
rho_w+p_w = -1/(2 pi^2 N B R^2)
              integral_0^infinity (a K_a+2b K_b)/(exp(zeta)+1) dxi <= 0.    (8)
```

This is a statement about the finite difference (1), not the full
renormalized stress. Metric variation is legitimate on compact sets of
positive coefficients: (4) and (7) dominate the differentiated integrals
by an integrable exponential times a linear frequency polynomial.

Suppose `J>=T>0` is known from just a throat segment, and set `u=lambda T`.
Both moment bounds decrease as their argument increases:

```text
integral exp(-max(u,xi L)) dxi = (u+1)exp(-u)/L,
integral xi exp(-max(u,xi L)) dxi = (u^2/2+u+1)exp(-u)/L^2.
```

Together with (7), these imply the pointwise estimate

```text
0 <= -(rho_w+p_w)
 <= coth(u) exp(-u)/(2 pi^2 N R^2)
       [ lambda(u+1)/(R L) + (u^2+2u+2)/(N L^2) ].                        (9)
```

The `q` zero modes have the exact sum
`(rho_w+p_w)_LLL = -q/(12 N^2 R^2 L^2)`.

## Sum over every higher angular mode

Put `s=sqrt(q+1)`, `u_1=T s`, and assume `u_1>=3`. Let
`gamma=N(x)L/(R(x)T)` at the observation point, and define

```text
P(u)=(gamma+1)u^2+(gamma+2)u+2,
Q(u)=(gamma+1)(u^3+3u^2+6u+6)
       +(gamma+2)(u^2+2u+2)+2(u+1).
```

Including both signs of each nonzero angular eigenvalue, (6) and (9) yield

```text
|E_w,HLL|/|E_w,LLL|
 <= 24 exp(-u_1)/(pi^2 q)
       [(q+2)(1+u_1)+2(s^2+3s/T+3/T^2)],                                (10)

|(rho_w+p_w)_HLL|/|(rho_w+p_w)_LLL|
 <= 12 coth(u_1) exp(-u_1)/(pi^2 q)
       [(q+2)P(u_1)+2 Q(u_1)/T^2].                                      (11)
```

Here is the tail argument, rather than a finite truncation. For `v>=1`,
`lambda(v)=sqrt(v(v+q))` obeys
`(q+2v)dv=2lambda d lambda`. The summands have the form
`(q+2v) p(Tlambda(v)) exp(-Tlambda(v))`, where `p` is a polynomial of degree
at most two with positive coefficients. Its logarithmic derivative is at
most `lambda'(3/lambda-T)`, using
`2/(q+2v)<=lambda'/lambda` and `p'/p<=2/u`. Thus it decreases when
`Tlambda>=3`; the sum is bounded by its first term plus the integral from
one to infinity. Elementary integration gives (10)-(11). Bound `coth(u_n)`
by `coth(u_1)` before applying this comparison. This also proves absolute
convergence of the mode-summed stress difference and justifies the retained
variations. The local factor `gamma` is not assumed universally bounded.

For an exact product throat segment `-rho_c<=rho<=rho_c`,
`B/R=1/sqrt(1+rho^2)`, so `T=2 asinh(rho_c)` is available regardless of the
rest of the positive spherical circle. At a point with `NL/R=pi`, examples
of the analytic upper bounds are:

| q | rho_c | Higher-mode energy / LLL | Higher-mode null stress / LLL |
| --- | --- | --- | --- |
| 3 | 2 | 0.117682 | 0.820956 |
| 3 | 8 | 0.000864 | 0.008088 |
| 8 | 8 | 0.000003506 | 0.000047896 |
| 16 | 8 | 0.000000008118 | 0.000000150294 |

These are evaluations of (10)-(11), rounded upward, for the specified model
parameters. They are not certified measurements or bounds for the actual
MMP geometry. A long exterior optical path changes `L` and hence `gamma`.

## Source comparison and unresolved physical bridge

[Freivogel, Fumagalli and Tomasevic, v2](https://arxiv.org/html/2606.12528v2),
section 4.3 and appendix C, give the warped spherical Dirac reduction and
distinguish propagating probes from support energy. Their cited
[Bolognesi, Bucciotti and Luzio, v2](https://arxiv.org/html/2406.15552v2),
section 2.1, estimates the summed higher-angular-wave suppression at a
monopole core. Neither inspected derivation supplies the observable (1)
or a bound on the local winding null stress (8). Their larger proof chains,
editions and other physical sectors have not all been audited.

[Sundberg and Jaffe](https://arxiv.org/html/hep-th/0308010v1), sections 2-4,
establish Dirac transfer/scattering methods, mass-barrier Casimir energies
and the need to distinguish total energy from a finite interaction
subtraction. These methods are precedents, not new claims of this work.
Comparison with other inhomogeneous determinant and Casimir-stress results
remains open; search non-results do not establish novelty.

The decisive physical gaps are specific:

1. The real exterior connecting two MMP mouths is not a globally spherical
   `S1 x S2` product. It can mix angular modes. The scalar positive transfer
   argument (3)-(7) has not been proved for that boundary operator.
2. The Bloch/holonomy average in (1) differs from the reflecting-AdS reference.
   Local vacuum polarization and the remaining reference-state difference
   must be matched before (9) bounds the complete source opening the throat.
3. A bound on the static source is not a bound on its dynamical response,
   nor on the coupled gravitational spectrum or finite-energy transport.
   The prescribed gauge-field assumption also leaves photon response out.

Decision: retain this as a bounded feasibility argument for a quantitative
LLL-support target. The integrated barrier does control energy and local
null stress in a nonconstant spherical model. The actual exterior operator
and reference-state comparison, rather than another local mass-gap estimate,
are now the hardest steps. Do not accept it as an MMP result or start Lean
formalization before the contribution comparison resolves these boundaries.

## Verification

[The check](../checks/warped_winding_bound.py) verifies constant-mass
normalizations against Bessel sums, transfer and local-insertion bounds,
independent metric finite differences, and a smooth AdS transfer against an
ODE solver with spatial refinement. Profiles include an exterior angular
mass reduced by a factor of 100. Exact polynomial primitives and all-mode
bound samples are checked separately. An independent Fourier Hamiltonian
calculation for `m(x)=0.7+0.3 cos(x)` uses an energy-based heat regulator
and a spatial-phase average. Its regulator extrapolation agrees with the
transfer winding energy within `3.7e-10`; the spectral truncation and
extrapolation error are not certified. The ordinary positivity and summation
arguments above supply the inequalities; a finite sample does not prove them.
The [receipt](../evidence/warped-winding-bound.json) retains results and scope.
Floating-point/ODE errors are not all certified. No Lean invocation occurred.

Remaining to-do list: finish the precise contribution/source comparison,
resolve angular mixing and the reference-state bridge, then establish the
physical consequence, verify the complete endpoint and prepare the manuscript.
