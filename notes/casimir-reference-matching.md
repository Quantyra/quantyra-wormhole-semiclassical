# Invariant vacuum shifts and leading Casimir opening

2026-10-10. S051/S052 baseline and contribution decision. The combined
first-order throat is determined for arbitrary invariant reference pressures
and a specified leading lowest-Landau-level (LLL) state difference. The
invariant part shifts the product radii; the Casimir part supplies the flare.
This reproduces established ingredients and their linear superposition. It
does not establish a new wormhole, corrected transit time or stability result.

## Model and order

Use mostly-plus signature, one massless unit-charge Dirac field, fixed magnetic
flux `F=(q/2) du wedge dphi`, and Einstein-Maxwell normalization
`kappa=8 pi G`, `r_e^2=kappa q^2/(8 g_c^2)`. The sign of the flux depends on
the choice `u=cos(theta)` orientation and does not affect this calculation.
Introduce a common loop bookkeeping parameter `epsilon` and write

```text
ds^2 = r_e^2 [-S(rho) d tau^2 + d rho^2/S(rho) + C(rho) d Omega^2],
f = 1+rho^2,   S=f+epsilon s,   C=1+epsilon c.
```

Reflection symmetry fixes `c'(0)=s'(0)=0`; choose `s(0)=0` as the central time
normalization. Evaluate the quantum source on the zeroth-order product:

```text
kappa T_Q = epsilon diag(-U f-chi, U/f-chi/f^2,
                        V/(1-u^2), V(1-u^2)).                         (1)
```

Here `U,V` are constants. If `p_A,p_S` are the renormalized invariant
pressures from [the magnetic product reference](magnetic-product-vacuum.md),
then `epsilon U=kappa r_e^2 p_A` and `epsilon V=kappa r_e^2 p_S`.
Using `epsilon=alpha=kappa q/(32 pi^2 r_e^2)=g_c^2/(4 pi^2 q)` gives
`U=32 pi^2 r_e^4 p_A/q` and `V=32 pi^2 r_e^4 p_S/q`.

The Casimir amplitude is `chi=1` for the specified long-throat, antiperiodic
LLL state difference with negligible ambient optical length. In conformal
coordinates `x=atan(rho)`, the circle-minus-reflecting-strip difference is
`Delta T_2,tau tau=Delta T_2,xx=-q/(8 pi)`. Dividing by `4 pi r_e^2` and
transforming back gives the two entries proportional to `chi` in (1).
The common local anomaly cancels in this state difference; it remains in
the invariant reference pressures. This is a prescribed leading LLL source,
not a calculation of all modes in the actual exterior-connected state.

At fixed charge, the quantum source's response to an order-epsilon metric
change enters at the next order of this bookkeeping expansion. That fact
does not control a simultaneous large-charge limit or justify omitting a
state-dependent source of the same retained order. On a finite throat patch
require `|epsilon U|, |epsilon V|, |epsilon chi rho| << 1`. No uniform error
estimate for the full semiclassical problem is supplied here.

## Direct Einstein-Maxwell derivation

With primes denoting rho derivatives and the common radius scaled out, direct
metric contraction gives

```text
G_tau tau/S = 1/C - S C''/C - S'C'/(2C) + S(C')^2/(4C^2),
G_rho rho S = -1/C + S'C'/(2C) + S(C')^2/(4C^2),
G_theta theta/C = S''/2 + S C''/(2C) + S'C'/(2C) - S(C')^2/(4C^2).
```

The corresponding magnetic stress combinations are `1/C^2,-1/C^2,1/C^2`.
Consequently the three independent first-order equations are

```text
c - f c'' - rho c' = -U - chi/f,
-c + rho c'        =  U - chi/f,
s''/2 + f c''/2 + rho c' + 2c = V.                         (2)
```

Adding the first two equations gives `c''=2 chi/f^2`. Integrate and impose
reflection symmetry; the radial equation fixes the remaining constant.
Substitution into the angular equation then gives

```text
c = chi (1+rho atan(rho)) - U,
s = (V+2U) rho^2 + chi s_C,
s_C = -rho^2-rho(rho^2+3) atan(rho)+log(1+rho^2),
s_C'' = -6(1+rho atan(rho)).                                (3)
```

The two possible homogeneous lapse terms are removed by `s(0)=s'(0)=0`.
All four diagonal equations, six off-diagonal equations, magnetic Maxwell
equations and background covariant conservation of (1) are checked by
[direct tensor contraction](../checks/casimir_reference_matching.py).
The contraction starts from Christoffel symbols and the magnetic two-form,
rather than assuming the reduced equations (2).

The scalar and radial null checks are

```text
r_e^2 R_4 = -2 epsilon (U+V) + O(epsilon^2),
k=(1/f,1,0,0),
kappa T_Q,mu nu k^mu k^nu = -2 epsilon chi/f^2,
R_areal = r_e [1+epsilon c/2] + O(epsilon^2),
R_areal'' = epsilon r_e chi/f^2 + O(epsilon^2).              (4)
```

Thus the leading opening depends on `chi` alone. For positive `chi`, the
reflection point is a local area minimum within the perturbative solution.
For `chi=0`, the sphere has constant radius and there is no first-order
flare. Neither statement proves global traversability.

## Reference radii and known special cases

For `chi=0`, put `A=1+epsilon(V+2U)` and rescale
`rho_hat=sqrt(A) rho`, `tau_hat=sqrt(A) tau`. The reference product has

```text
a_ref^2/r_e^2 = 1-epsilon(V+2U) + O(epsilon^2),
b_ref^2/r_e^2 = 1-epsilon U + O(epsilon^2),
a_ref^2/b_ref^2 = 1-epsilon(U+V) + O(epsilon^2).              (5)
```

These agree with the independently varied pressure equations of the
previous reference calculation. Keeping the full Dirac trace therefore
changes the reference geometry without removing the leading Casimir flare.
This resolves the first-order baseline question raised by the
[trace audit](throat-trace-anomaly.md); it does not complete quantum closure.

Setting `U=-1/3,V=1/3,chi=1` recovers the W22 throat:
`c=rho atan(rho)+4/3` and
`s=-4rho^2/3-rho(rho^2+3)atan(rho)+log(f)`.
That traceless invariant lift is a special case of (3). Replacing its
constants by full reference pressures changes the constant and quadratic
terms, while the leading opening term is unchanged.

The relevant precedents are explicit:

- [Maldacena, Milekhin and Popov, v3](https://arxiv.org/html/1807.04726v3#S5.SS2),
  section 5.2 and appendix F: LLL Casimir support and anomaly accounting;
  massive invariant contributions shift the two radii. Their discussion
  also anticipates suppression of massive boundary effects deep in a throat.
- [Kanai, Maeda and Yoshida, v1](https://arxiv.org/html/2511.21017v1),
  equations (2.84)-(2.92): the Casimir area and lapse solution. Sections
  3.2-3.3: invariant sources have equal mixed AdS pressures, no leading
  flare, and a product geometry after coordinate rescaling. Their all-order
  local-EFT argument is not an all-order theorem for a nonlocal quantum
  state coupled through the exterior.
- [Sadhukhan, v1](https://arxiv.org/html/2610.09847v1), section 2:
  the special case above, independently reproduced in the earlier
  [throat audit](axial-baseline-and-throat-energy.md).

W25 sections 2.2-2.3 and 3.1-3.3 have been read at derivation level. Its full
dependency/edition audit and five-dimensional branch remain outside this
comparison. The source register does not mark its full proof review complete.

## What matching does and does not determine

The large-rho coefficient of (3) is

```text
d R_areal/d rho = epsilon r_e chi pi/4 + lower overlap terms.
```

This is an overlap coefficient in `1 << rho << 1/|epsilon chi|`, not a limit
of a uniformly valid solution at fixed nonzero epsilon. Classical MMP
matching with `tau=t/ell` gives `R-r_e approximately r_e^2 rho/ell`, hence

```text
ell_0 = 4 r_e/(pi epsilon chi),
ell_0(chi=1,epsilon=alpha) = 16 r_e^3/(G q).                 (6)
```

The global transit time beyond (6) is not obtained by simply replacing
`r_e` with one of (5)'s corrected radii. To see the missing information,
suppose an exterior near an extremal horizon has the form

```text
ds_ext^2 approximately -Z^2 (R-b)^2/a^2 dt^2
                      +a^2 dR^2/(R-b)^2 +R^2 dOmega^2.
```

Use canonical product coordinates for the throat, with `tau_hat=t/ell`.
Matching then requires `R-b approximately a^2 rho_hat/(Z ell)`. If `B`
denotes its opening slope in those same coordinates, `ell=a^2/(Z B)`.
The exterior fixes `Z`; the corrected slope requires the mixed vacuum/state
response and consistent coordinate normalization. Ambient optical lengths
also affect the LLL state. These data are absent from a homogeneous
determinant. W25's equations (2.75)-(2.77) already separate corrected
exterior matching from exact Reissner-Nordstrom matching.

## Decision and verification

Classify homogeneous pressure shifts plus the leading opening solution as
**already known ingredients / reproduced baseline**, and remove this route
as a standalone candidate for a new-physics paper. The calculation remains
a prerequisite for a total response or global-matching result. Do not start
a separate Lean campaign for this elementary linear superposition.

The next feasibility question is quantitative: can the omitted nonlocal
state response be controlled in the actual throat/exterior geometry tightly
enough to change or establish a stability or transport conclusion? The
existing constant-radius cylinder bounds alone do not answer it; the sphere
and its local angular mass gap vary outside the throat. The conditional
drift energy also has a gapless channel, so small omitted terms cannot be
declared harmless by a spectral-gap argument. These are requirements for
the next candidate comparison, not a newly accepted primary or fallback.

The [receipt](../evidence/casimir-reference-matching.json) records 44 exact
symbolic identities, including conservation, the scalar trace, radius
rescaling, the W22 specialization and the MMP leading length. This is
reproducible computer algebra supporting the ordinary derivation. There is
no global state calculation, controlled higher-order error, physical
stability result, Lean certificate or novelty claim.

Remaining to-do list: complete decisive source comparisons and primary/fallback
selection; bound the relevant global state/response corrections; establish
the consequential physical result, then verify and prepare the manuscript.
