# Exterior Dirac reduction and a finite window error target

2026-10-10. S052/S054. The static, axisymmetric, flat point-dipole Weyl
equation reduces exactly to independent two-dimensional Dirac differential
operators. Their mass-zero contours are the magnetic field lines, with
the midpoint flux labels previously used as a comparison discretization.
This identifies an actual exterior operator to analyze. It does not yet
identify its finite-mouth scattering map or control its projection error.

A second ordinary result makes the energy comparison tolerant of an
approximate return operator on a closed imaginary-frequency window away
from zero. This replaces the need for an exact unequal-delay inequality
on every low frequency. The [preceding gapless audit](dipole-transverse-modes.md)
and its assumptions remain valid; the present argument does not restore
a global transverse spectral gap.

## Exact axial reduction in the leading flat exterior

Use charge convention `e=1` and one massless Weyl sector
`H = sigma . (-i grad-A)` in cylindrical coordinates `(r,phi,z)`.
The monopoles lie at `z=+-a`, with `d=2a` and positive integer `q`.
Write

```text
R_minus = sqrt[r^2+(z+a)^2],   R_plus = sqrt[r^2+(z-a)^2],
u(r,z) = [(z+a)/R_minus-(z-a)/R_plus]/2,
A = A_phi dphi,              A_phi = q(1-u),
B_r = (q r/2)(R_minus^(-3)-R_plus^(-3)),
B_z = (q/2)[(z+a)R_minus^(-3)-(z-a)R_plus^(-3)].       (1)
```

This gauge is regular on the axis between the monopoles. Outside that
segment, the axis has `A_phi=q`; the usual integer gauge change removes
the string in another patch. All differential formulas below are on
`r>0`, away from the monopoles. Regularity on the axis and matching at
physical mouth surfaces are domain conditions, not freely selectable
reflecting walls.

The eigenvalues `ell` of `J_z=-i partial_phi+sigma_z/2` are half-integers.
The unitary angular and radial-measure transformation is

```text
Psi_ell(r,phi,z) = (2 pi r)^(-1/2) exp(i ell phi)
                         exp(-i sigma_z phi/2) chi_ell(r,z).
```

It preserves the norm after integrating over `phi`. Direct substitution
gives, without a Landau-level approximation,

```text
H_ell = -i sigma_x partial_r - i sigma_z partial_z
                                  + sigma_y M_ell(r,z),
M_ell = [ell-A_phi(r,z)]/r.                            (2)
```

The `1/(2r)` terms from the spin rotation and radial measure cancel.
The reduced radial and axial currents are `chi^dagger sigma_x chi` and
`chi^dagger sigma_z chi`; therefore the transformation also preserves the
corresponding integrated flux. The other four-dimensional Weyl sector
has the opposite Hamiltonian sign. Mixing different `ell` sectors requires
an axisymmetry-breaking term or nonaxisymmetric boundary matching.

For `0<ell<q`, the equation `M_ell=0` is precisely the field-line contour

```text
u_ell = 1-ell/q,
{u_ell: ell=1/2,3/2,...,q-1/2}
        = {(k+1/2)/q: k=0,...,q-1}.                  (3)
```

This derives the midpoint labels as mass-zero contours. It does not prove
that each contour supports an isolated propagating quantum channel.
In particular, `q` contours do not contradict the `q-1` normalizable
zero modes of the full-plane transverse problem. The latter is a
different spectral problem. At large `r` for fixed `z`,
`M_ell ~ (ell-q)/r`, so the mass tends to zero rather than opening a
uniform bulk gap.

The square of (2) is also exact:

```text
H_ell^2 = -partial_r^2-partial_z^2+M_ell^2
                   +sigma_z partial_r M_ell-sigma_x partial_z M_ell. (4)
```

On the mass-zero contour,
`grad M_ell=(-B_z,B_r)`, and hence `|grad M_ell|=|B|`.
This is the local magnetic confinement scale, not a positive spectral
gap for the complete exterior.

## Field-aligned coordinates and the coupling terms

Orient a contour by the magnetic field and parametrize it by arclength
`s`. Let `t=(sin theta,cos theta)` and `n=(-cos theta,sin theta)` be its
unit tangent and normal in the `(r,z)` plane. Write `t'=kappa n`, so
`theta'=-kappa`. In a tubular patch,

```text
x(s,n) = gamma(s)+n n(s),             h(s,n)=1-kappa(s)n>0.
```

Rotate by `S(s)=exp[-i theta(s) sigma_y/2]` and absorb the Jacobian by
`tilde_chi=sqrt(h) S^dagger chi`. The resulting operator on the flat
`ds dn` measure is exactly

```text
H_F = i sigma_x partial_n
          -i sigma_z h^(-1/2) partial_s h^(-1/2)
          +sigma_y M_ell(s,n).                        (5)
```

The normal spin-connection term cancels the normal Jacobian derivative.
The symmetric tangential derivative retains the curvature dependence;
it must not be replaced by `partial_s` before estimating the error.
Also, `M_ell(s,n)` is not generally linear in `n`.

The local straight linear model has `M=beta n`, `beta=|B|` on the contour.
Its normal operator annihilates
`exp(-beta n^2/2)(1,0)`, and its square is
`-partial_n^2+beta^2 n^2-beta sigma_z`. The associated positive-spin
branch travels at unit speed along `s`. This reproduces the leading
field-line picture, while (5) explicitly exposes the curvature, changing
width and nonlinear mass terms needed for an error estimate.

In the prolate coordinates used in the [dipole-delay note](dipole-delay-bound.md),

```text
r/a=sqrt[(mu^2-1)(1-eta^2)], z/a=mu eta,
u=mu(1-eta^2)/(mu^2-eta^2),
|B|^2=(q^2/a^4)
 [mu^2+3mu^2 eta^2+eta^4-eta^2]/(mu^2-eta^2)^4.      (6)
```

For `mu>=1` and `|eta|<1`, the numerator is at least `mu^2`, and
`mu^2-eta^2>=mu^2(1-eta^2)`. Substitution of the expression for `u`
therefore gives

```text
|B| >= q u^3/[a^2(1-eta^2)^2] >= q u^3/a^2.          (7)
```

Equality in the second comparison occurs at the midplane. The local
linear normal width there is `a/(sqrt(q) u^(3/2))`. It is not uniformly
small on the outermost contours `u=O(1/q)`. Near the axial endpoint
`u=1`, the distance to the cylindrical axis is another required scale.
Uniform tube estimates should initially use `u` in a compact subinterval
of `(0,1)` and trim the neighborhoods of the monopoles. Equation (7)
alone does not justify a global oscillator approximation.

## Comparison with established curved interface results

[Aoki and Fukaya, section 2](https://arxiv.org/html/2203.03782v3)
derive induced spin connections and edge localization on curved mass
walls. Their large step-mass model differs from the smooth mass in (2).
[Bal et al.](https://arxiv.org/html/2106.00729v1) already prove approximate
unit-speed propagation of well-prepared wavepackets along smooth curved
mass-zero interfaces. [Bal's later analysis](https://arxiv.org/html/2206.09439v2)
uses tubular coordinates, transverse mode expansions and energy estimates
for relativistic and dispersive packets. These are direct methodological
predecessors, not new mechanisms supplied by our reduction.

The inspected results assume suitable globally smooth coefficients,
nondegenerate interfaces and prepared initial data. Our mass has axis and
mouth singularities, and vanishes at infinity. A smooth extension of a
trimmed tube can support a local comparison, but does not automatically
preserve the physical port scattering problem. Their real-time packet
bounds also do not automatically give an operator norm on the entire
incoming channel space. The complete relevant propagation proofs and the
physical matching specialization remain under review.

## An approximation-tolerant Gaussian energy theorem

Retain the leading fixed-exterior energy model

```text
E(T)=A/T^2+q pi/(24T)+E_loop(T), A>0,
E_loop=-(1/pi) sum_sigma integral_0^infinity
               Re log det[I+exp(-xi T) R_sigma(xi)] dxi.
```

For each chirality let the physical-model return `R_sigma` be causal and
passive with common front `d`; this is a hypothesis about flux-normalized
ports. On a window `lambda<=xi<=Lambda`, suppose a comparison `R0_sigma`
obeys the unequal-delay exterior-power condition for `d_j=d/u_j>=d`,
and the approximation satisfies

```text
||R_sigma(xi)-R0_sigma(xi)|| <= epsilon exp(-xi d).    (8)
```

The comparison norm is also at most `exp(-xi d)`. Its delay-label CDF is
within `delta` of uniform. `R0` need not share eigenvectors with `R`, and
both may be nonnormal. All exterior data are independent of `T` in this
derivative, as in the preceding energy notes.

Put `X=exp(-xi T)R`, `X0=exp(-xi T)R0`, and `L=T+d`. The resolvent
identity and `||X||,||X0||<=exp(-xi L)` give

```text
|Tr[X(I+X)^(-1)-X0(I+X0)^(-1)]|
       <= q epsilon exp(-xi L)/(1-exp(-xi L))^2.       (9)
```

This follows by writing `X(I+X)^(-1)=I-(I+X)^(-1)` and using
`|Tr Y|<=q||Y||`. No eigenvector conditioning assumption is needed.
For `xi<lambda`, common-front passivity gives the one-sided bound
`Re Tr[X(I+X)^(-1)]<=q/2`. For `xi>Lambda`, use the preceding exponential
tail bound. On the window use (9) and the positive unequal-delay
comparison, extending the latter to all frequencies to obtain an upper
bound. Since `L>=d`, the integrated perturbation is at most

```text
(2q epsilon)/(pi d^2) F(y),     y=lambda d,
F(y)=integral_y^infinity t exp(-t)/(1-exp(-t))^2 dt
    =y/(exp(y)-1)-log(1-exp(-y)).                     (10)
```

Consequently, for every `0<T<=d`,

```text
E'(T) <= -2A/d^3 -(q/d^2) C(delta,y,x,epsilon),
x=Lambda d,
C = pi[8 log(2)-5-delta]/24
       -(2/pi) exp(-x)(x+1)
       -y^2/(2pi)-(2 epsilon/pi)F(y).                 (11)
```

The infrared cost `q lambda^2/(2pi)` comes directly from integrating the
`q/2` trace bound for the two chiralities. It needs no approximation
at `xi=0`. Negative low-frequency trace values only improve this
one-sided estimate. Thus resonant phase choices do not invalidate (11)
when the stated passive-return assumptions hold.

One explicit sufficient target is

```text
delta <= 1/20,     lambda d=1/20,     Lambda d=4,
epsilon <= 1/1000,
C >= 0.003576370769... >0.                            (12)
```

The decimal is a numerical evaluation, not the proof of positivity.
An exact rational certificate follows from `F(1/20)<=1+log(21)<81/20`,
positive lower series for `log(2)` and `exp(4)`, and `pi>3.14159`.
For `log(21)<61/20`, the first twenty terms of `exp(61/20)` already
exceed 21. The check retains the rational certificate after multiplying
the conservative lower bound for `C` by `24pi`.

An additional energy term preserves the exclusion if its positive
`T` derivative is uniformly less than `2A/d^3+q C/d^2`. No physical
omitted contribution has yet been bounded that way.

Replacing `b` midpoint labels in (3) by `u=1` has CDF error at most
`(b+1/2)/q`. In a fully axisymmetric port problem, sectors remain
independent, so uncontrolled sectors can be placed in the comparison
using their own passive returns and delay `d`; no approximation to them
is required. With coupling between sectors, (8) must control the full
return operator. A fraction of modes can be reserved to avoid both weak
outer confinement and the near-axis tube degeneracy, provided the CDF
budget and the full return hypotheses are verified.

## What the closed frequency window changes

The next estimate need only cover `0.05<=xi d<=4`, rather than approach
the zero-energy threshold. To see the mathematical benefit, let `H` and
`h` be self-adjoint bulk and effective Hamiltonians, `J` a bounded
identification, and suppose an operator estimate on the chosen input
space is proved for all `t>=0`:

```text
||exp(-itH)J-J exp(-ith)|| <= C0 alpha (1+t).
```

The Laplace representation of the resolvent immediately gives

```text
||(H-i xi)^(-1)J-J(h-i xi)^(-1)||
       <= C0 alpha (1/xi+1/xi^2), xi>0.              (13)
```

Time and energy here are in consistent fixed dimensionless units.
Equation (13) follows by integrating the evolution estimate against
`exp(-xi t)`, so no global spectral gap is used. The input estimate is
substantive: a bound for one packet with an uncontrolled Sobolev constant
does not establish it uniformly. Nor is a bulk resolvent matrix element
automatically a flux-normalized return between mouth ports.

For a trimmed contour family, rescaling (2) by `a/q` gives the usual
semiclassical Dirac form with small parameter `1/q` and mass
`[ell/q-A_phi/q]/(r/a)`. Establish a uniform prepared-input estimate,
include the actual mouth matching, and translate it to (8). This is a
specific route to test, not a claimed application of the cited theorems.
For example a valid bound with a fixed constant times `q^(-1/2)` would
eventually meet (12), but neither that constant nor that port estimate
is established here.

## Physical problem and contribution decision

Use matching surfaces of radius `r_p` in the overlap `r_e<<r_p<<d` around
the magnetic mouths, where `r_e` is the throat/mouth gravitational scale.
Derive the incoming/outgoing low angular modes and flux normalization
from the throat and exterior fields, including the evanescent higher
modes; do not impose artificial reflecting boundary conditions. The
flat differential expression (2) applies in the leading overlap/exterior
approximation. Gravitational lapse, physical front distance, finite-port
geometry, gauge interactions and the vacuum subtraction remain separate
terms in the final energy/observer comparison.

There is an exact finite-port geometric simplification in this flat model.
Let the two spherical port radii be `r_p<d/2`, and put `D=d-2r_p`.
A classical field line of label `u` crosses the midplane a distance
`d/(2u)` from each monopole center. By the triangle inequality, its
length between the two port surfaces obeys

```text
L_port(u) >= d/u-2r_p >= D/u,
(d/u-2r_p)-D/u = 2r_p(1/u-1) >=0.                   (14)
```

Thus the same comparison distribution applies with the actual flat
surface-to-surface front `D`. Apply (8)-(12) with every front-time `d`
replaced by `D`; the physical matching target is `0.05<=xi D<=4`.
This geometric statement does not identify a quantum mode with its
center contour. It also does not settle the common-emission and
receiver-clock comparison: the minimum distance between two surfaces
is not the exterior travel time for every fixed pair of observer events.
Lapse and event/port matching still require an explicit final calculation.

The exact axial reduction and ordinary window bound are verified. The
primary candidate remains controlled robustness of the existing MMP
no-shortcut comparison; the final primary/fallback and novelty verdicts
remain pending. A leading curved-edge reproduction alone would be known
physics. The needed additional result is the actual port return/error
bound and its total-energy consequence. Keep the remaining work on that
question rather than expand a generic domain-wall theory campaign.

## Verification

[The check](../checks/exterior_dirac_window.py) and
[receipt](../evidence/exterior-dirac-window.json) verify 12 exact identities,
18 independent Cartesian finite-difference reductions, 30 field-line
geometry cases with three derivative resolutions, 24 noncommuting matrix
perturbations, five integral comparisons and the positive rational margin.
The initial coarse gradient stencil slightly exceeded its tolerance;
refinement reduces the maximum residual from about `2.2e-8` to `3.9e-10`.
These supplementary checks do not certify a complete numerical error or
simulate physical mouth scattering. No Lean invocation occurred.

Remaining to-do list: prove the prepared-input and physical-port estimates
on the closed window, control total energy and observer corrections,
complete contribution selection, then certify and prepare the manuscript.
