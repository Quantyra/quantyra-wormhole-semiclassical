# Static formation support: a radial pass defeated by warping

2026-10-12, E003/S068. Units `G=c=1`. This calculation concerns the
background that would hold prospective mouths apart. It constructs no
wormhole and supplies no quantum stress or formation rate.

**Decision:** reject the infinitesimally thin perfect-fluid support shell.
There is a nonempty overlap of endpoint confinement, dominant energy, and
radial stability. It nevertheless has a known nonspherical instability.
An arbitrarily thin regularization with only bounded subleading restoring
terms does not remove the short-wavelength obstruction. This is an
application of existing methods and instability results; novelty is not
accepted.

## 1. Formation question and predecessor

The orbiting MMP construction remains metastable; S066/S067 did not prove
formation or full stability. The new avenue uses the test-string action of
[Horowitz, Marolf, Santos and Wang (HMSW), 1904.02187v2](https://arxiv.org/html/1904.02187v2),
equations (3.4)–(3.10). With `ds²=-A(r)²dt²+dr²/g(r)+r²dΩ²`, define

```
I(r) = integral_0^r A(u)/sqrt(g(u)) du,
U(r) = m A(r) - mu I(r),             k = mu/m > 0.
```

`2U` is the static pair energy relative to the unbroken test string. A
smooth zero-energy resting endpoint requires `U=U'=0` and a positive
quadratic restoring term. The full relativistic endpoint energy, rather
than a Newtonian force balance, fixes the subtraction and redshift here.

HMSW section 2 distinguishes its nearly static instantons from the exactly
degenerate limit, in which the decay prefactor vanishes. Its black-hole
replacement and quantum-opening argument in section 7 must be established
on any different background. A probe equilibrium is not that argument.

## 2. A shell that passes the limited mechanical tests

Take a material shell at radius `R`, with Minkowski interior and
Schwarzschild exterior, normalized to the exterior clock. Put

```
C = M/R,   s = sqrt(1-2C),   0<C<1/2,
A = s, g = 1                         (r<R),
A = sqrt(1-2M/r), g = 1-2M/r         (r>R).
```

The center is regular; the exterior starts outside `2M`. This is an ordinary
interior/exterior junction, not two exterior regions glued into a wormhole.
The string is a probe. Its infinite extent would introduce a conical
asymptotic defect once backreaction is included.

Direct integration gives `I=sr` inside and `I=sR+r-R` outside. Thus
`U(R)=ms(1-kR)`, and zero energy selects `kR=1`. At that value,

```
(R/m) U'(R-) = -s,
(R/m) U'(R+) = C/s - 1.
```

The shell location is a strict local **cusp** minimum precisely when
`C>sqrt(2)-1`. Equality is not a minimum: the outer second derivative is
negative. No smooth-orbit formula is evaluated at the discontinuity. A cusp
does not itself supply HMSW's limit of arbitrarily small acceleration.

The Israel junction conditions give the surface energy density and pressure

```
Sigma = (1-s)/(4 pi R),
p = (1-s)^2/(16 pi R s),
p/Sigma = (1-s)/(4s).
```

They obey the dominant energy condition when `s>=1/5`, equivalently
`C<=12/25`. Notice that `p` is a compressive pressure, not a tension.

For an isolated barotropic shell let `eta=dp/dSigma` at equilibrium, keep
the ADM mass `M` fixed, and write its proper-time equation as

```
Rdot² + V(R) = 0,
V = 1 - M/R - m_s²/(4R²) - M²/m_s²,
m_s = 4 pi R² Sigma,
Sigma' = -2(Sigma+p)/R.
```

Differentiating these equations before substituting the static junction
values gives

```
R² V'' = [12 eta s³ + 4 eta s² + 3s³ - s² - s - 1]/(2s²),
eta_crit = (1-s)(3s²+2s+1)/[4s²(3s+1)].
```

This is a radial test only. `V''>0` requires `eta>eta_crit`; imposing the
local causal range `0<=eta<=1` requires `s>s_*`, where `s_*` is the unique
positive root of `15s³+3s²-s-1=0`. Descartes' rule gives uniqueness, and
the root is between `0.39529` and `0.39530`. The joint interval is

```
sqrt(2)-1 < C < (1-s_*²)/2 = approximately 0.421871,
eta_crit(s) < eta <= 1.
```

The dominant energy condition holds throughout. A strictly subluminal
example is `s=2/5`, `C=21/50`, `eta=99/100`. It has
`p/Sigma=3/8`, `eta_crit=171/176`, `R²V''=81/1000`, and endpoint slopes
`-2/5` and `1/20`. These exact values are reproducible diagnostics, not a
claim of complete material stability.

The checker differentiates the Israel potential symbolically and also
evaluates it using the independently integrated conservation equation for
the local affine equation of state. It compares `eta_crit` with the
independent equilibrium-sequence criterion (33) of
[Pitre, Schneider and Poisson, 2604.05980v1](https://arxiv.org/html/2604.05980v1#S3.SS4):
their adiabatic index is `Gamma=eta(Sigma+p)/p`, not `eta Sigma/p`.

## 3. Why the apparent window fails

[Yang, Bonga and Pan, 2207.13754v2](https://arxiv.org/html/2207.13754v2),
equations (18)–(19), include the displaced matter and gravitational response.
For local proper frequency and proper tangential wave number their leading
relation is

```
omega_proper² = -(p/Sigma) k_proper² + lower-order terms.
```

Positive `p/Sigma` therefore produces high-wave-number growth. The radial
sound derivative does not change this leading bending term. In the entire
endpoint-confining interval,

```
p/Sigma > 1/(2 sqrt(2)).
```

The dangerous coefficient cannot be tuned to zero while retaining this
particular endpoint support. The explicit example gives `omega²/k² -> -3/8`.
This sign has an order-one margin, so the radial pass cannot be promoted to
full stability. The conclusion is an application of the cited instability,
not a new generic shell-instability theorem.

Pitre, Schneider and Poisson study precisely the neutral
Minkowski/interior–Schwarzschild/exterior perfect-fluid system. Their
sections I and VII report growing even-parity modes at the sampled
multipoles and parameters. We do not turn that sampling into an all-finite-
multipole proof, or claim to have reproduced their full mode solver. The
analytic short-wavelength obstruction already suffices for this rejection.

## 4. Next variation tested: making the layer slightly thick

After the thin-fluid rejection, the next comparison is with the proposed
finite-thickness/rigidity cure. Yang et al. already discuss both after (23).
The following conditional inequality makes the relevant scale explicit; it
does not derive a new material law.

Let `epsilon=h/R -> 0`, `K=k_proper R`, `beta=p/Sigma>=beta0>0`. Suppose a
conservative normal-mode branch exists and its *real* squared frequency
obeys, uniformly on `1 << K << 1/epsilon`,

```
W² <= -beta0 K² + A K + B + D epsilon² K⁴,   W=omega_proper R,
```

with nonnegative finite `A,B,D` independent of `epsilon`. This assumption
includes only bounded lower-order gravity/curvature terms and a bounded
leading bending correction. It does not cover an arbitrary interacting
material, dissipative frequency-dependent kernel, or a rigidity coefficient
that diverges as the layer is thinned.

Choose integer angular wave numbers with `K=epsilon^(-1/2)(1+o(1))`.
They lie within the assumed band. Then

```
epsilon W² <= -beta0 + A sqrt(epsilon) + (B+D)epsilon + o(1) < 0.
```

Thus some branch still grows for sufficiently small thickness, under these
explicit hypotheses. A continuum sufficient threshold, when `A+B+D>0`, is

```
sqrt(epsilon) < 2 beta0 / [A + sqrt(A²+4 beta0(B+D))].
```

If all three coefficients vanish, the sign is already negative. The checker
also uses `K=sqrt(l(l+1))` for integer `l`, so the numerical examples do not
rely on an exactly selectable continuum wave number. Its illustrative
coefficients are declared test inputs, not measured or derived shell
coefficients. Only the conditional implication is accepted. A finite
thickness capable of changing the result requires a full constitutive and
gravitational calculation, not interpolation of this shell metric.

## 5. Source discrepancy and scope

The inspected [1708.04588v1](https://arxiv.org/html/1708.04588v1) equations
(38)–(39) agree with the static junction data. Its printed (40), specialized
to `psi_,R=0`, does not agree with direct differentiation of the potential
and the conservation equation above: it gives `2s(2eta+1)` for `R²V''`.
We preserve this discrepancy and do not use that expression as validation.
The independent sequence criterion from 2604.05980 and finite-difference
checks agree with the expression in section 2. No published erratum or
comprehensive source correction is claimed.

The original research target remains unachieved. This support model fails
before black-hole replacement, renormalized quantum support, and formation
probability are addressed. The finite-thickness variation above supplies a
necessary modeling condition, not a positive construction. A stable support
object would still not prove a stable wormhole.

Remaining to-do list: select a full finite-thickness or other explicit
matter action, test endpoint creation and all relevant support modes,
construct the backreacted mouths and quantum state, and control formation
and finite-disturbance traversability. No publication or Lean run is selected.
