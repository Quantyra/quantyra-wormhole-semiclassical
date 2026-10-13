# Higher Landau modes at the magnetic AdS boundary

S072, 2026-10-13. **The exterior radial barrier can suppress the boundary
response of higher Landau modes even when the collective strip frequency
exceeds their AdS2 mass.** An exact zero-frequency matching coefficient
and a conditional finite-frequency bound quantify this statement. The
condition is a bounded inner retarded load at a specified radial cut; it
has not been proved for the evolving wormhole. Bulk gauge polarization is
a different contribution and is not bounded by this boundary result.

This continues [collective electric response](collective-electric-response.md).
Its mass comparison remains correct, but was not an excitation-weight
calculation. The present result neither establishes a physical switching
energy bound nor rules out the candidate.

## Radial problem and boundary convention

Use the fixed, purely magnetic, extremal RN-AdS4 exterior, with extremal
radius `b`, `X=b^2/L4^2>0`, and `C=1+6X`. Set `u=b/r` and

```
f(r) = (1-u)^2 H_X(u)/u^2,
H_X(u) = X+2Xu+(1+3X)u^2,
dx/dr = 1/f,                         x(infinity)=0.
```

For one free massless four-dimensional Dirac field with unit electric
charge and monopole flux `q>=1`, the sphere eigenvalues and multiplicities
are

```
lambda_n = sqrt[n(n+q)],
d_0=q,              d_n=2(q+2n) for n>=1,
M_n=lambda_n/sqrt(C).
```

The nonzero-level multiplicity includes both signs of the sphere Dirac
eigenvalue. Spinor and radial-measure rescaling, followed by a constant
unitary rotation, gives for a positive eigenvalue

```
H_n = -i sigma_2 partial_x + W_n(x) sigma_1,
W_n = lambda_n sqrt(f)/r,
partial_x f_s = -W_n f_s + omega g_s,
partial_x g_s =  W_n g_s - omega f_s.                 (1)
```

Here `f_s,g_s` are spinor components, distinct from the metric function.
There is no electric `A_t` background. This is a free-field radial
calculation; it does not integrate fluctuating Maxwell or metric fields.

The radial gamma projector is proportional to `sigma_1` in this basis.
Choose source `J=(f_s+g_s)/sqrt(2)` and response
`A=(f_s-g_s)/sqrt(2)`, with the retarded sign convention

```
G=A/J=(r_s-1)/(r_s+1),       r_s=f_s/g_s,
rho=Im G=2 Im(r_s)/|1+r_s|^2 >=0.                     (2)
```

The flux is `2 Im(f_s^* g_s)=-2|J|^2 Im G`. A massless ingoing channel
has `(f_s,g_s)=(1,-i) exp(-i omega x)` and `G=i`, so `rho_0=1`.
This fixes the per-channel spectral normalization used below; it is not
a new identification of the physical inter-boundary coupling. Reversing
the sphere eigenvalue or the complementary source convention interchanges
the projectors. The other spectral denominator is `|1-r_s|^2`; the leading
weight and the upper bound below are the same. Real boundary contact
terms do not change `rho`.

The separation and gamma projectors are checked against
[Bintanja, Freivogel and Rolph, sections 2.1.1–2.1.2 and appendix A](https://arxiv.org/html/2308.00871v2).
The sphere spectrum is the one already used in the preceding checkpoint.

## Exact zero-frequency exterior coefficient

Define the positive radial action

```
I_X(u) = integral_0^u dv/[(1-v)sqrt(H_X(v))]
       = 1/sqrt(C) log{ [sqrt(C H_X(u))+C-(1+4X)(1-u)]
                       /[(1-u)(sqrt(CX)+2X)] },
Lambda_X = 2C/[sqrt(CX)+2X].                          (3)
```

For example, putting `z=1/(1-u)` changes the integral into
`integral dz/sqrt[C z^2-(2+8X)z+(1+3X)]`, which gives (3), with
`I_X(0)=0`. Direct differentiation provides an independent check.

The zero-frequency solutions normalized at the AdS4 boundary are
`f_s=exp(lambda_n I_X)` and `g_s=exp(-lambda_n I_X)`. Near the horizon,

```
sqrt(C) I_X = -log(1-u)+log Lambda_X+O(1-u),
zeta = b/[C(1-u)] [1+o(1)],
f_s ~ (C Lambda_X/b)^M zeta^M,
g_s ~ (b/[C Lambda_X])^M zeta^(-M).                 (4)
```

If the inner solution has coefficients `A_IR zeta^(-M)` in `g_s` and
`B_IR zeta^M` in `f_s`, the leading outer coefficient ratio is therefore

```
r_s(boundary) = [b/(C Lambda_X)]^(2M) B_IR/A_IR.      (5)
```

The regular zero-frequency branch is `g_s`. Its AdS4 source is nonzero
(`J=g_s/sqrt(2)`), so this exterior has no zero-frequency source zero
in the stated free, purely magnetic problem. The factors of `Lambda_X`
are fixed by the full exterior rather than by a freely chosen throat cut.
At large `X`, `Lambda_X -> 12/(sqrt(6)+2)=2.696938...`.

Equation (5) fixes the leading matching coefficient. It is not an exact
finite-frequency formula: real analytic outer response and inner matching
corrections must also be retained at nonzero frequency.

## Leading absorptive weight, including its normalization

For a Poincare-AdS2 throat with `W=M/zeta`, an exact ingoing solution is

```
g_s = sqrt(zeta) H^(1)_(M+1/2)(omega zeta),
f_s = sqrt(zeta) H^(1)_(M-1/2)(omega zeta).
```

Its small-argument expansion gives, for positive `omega`,

```
Im(B_IR/A_IR) = pi/Gamma(M+1/2)^2 (omega/2)^(2M),
rho_n(omega) ~ 2pi/Gamma(M+1/2)^2
                    [omega b/(2 C Lambda_X)]^(2M).  (6)
```

The equality for the inner absorptive coefficient is independent of real
analytic subtractions. At half-integer `M`, logarithmic terms replace a
singular real coefficient; the absorptive expression is continuous.
The asymptotic expression for `rho_n` is for `n>=1`, fixed geometry and
fixed `M>0`, as frequency tends to zero. It must not be extrapolated to
`M=0`: that channel has `r_s=i`, not a small coefficient, and (2) gives
`rho_0=1` exactly.

For an initial near-extremal thermal throat, put `y=omega/(2pi T)`.
With the same eigenvector normalization the finite-temperature result is

```
Im(B_IR/A_IR) = (pi T)^(2M) cosh(pi y)
                 |Gamma(M+1/2+i y)|^2/Gamma(M+1/2)^2,
rho_n(omega,T) ~ 2[pi T b/(C Lambda_X)]^(2M) cosh(pi y)
                 |Gamma(M+1/2+i y)|^2/Gamma(M+1/2)^2. (7)
```

This follows either from the inner spinor Green function or from
reflection and duplication identities applied to its gamma factors.
Specifically the ratio before taking its imaginary part is
`i (pi T)^(2M) Gamma(1/2-M)/Gamma(1/2+M)` times
`Gamma(1/2+M-i y)/Gamma(1/2-M-i y)`. The positive form in (7) avoids
subtracting large real contact terms. The reference is
[Faulkner et al., sections III–IV and appendices A.4, D.2, D.4](https://arxiv.org/html/0907.2694v2).
Only its matching method and neutral-electric AdS2 spinor result are
imported, not its electrically charged planar background or Fermi surfaces.
The independent Hankel/flux calculation fixes the convention and the
factor of two in (6); relying on the quoted eigenvector ratio without
checking its normalization would not suffice.

In the frozen clock comparison `omega=2pi T E`, so `y=E` is the global
strip frequency. For `k=2000`, `M_1=18.0481`, `Tb=7.46039e-24` and
`E=Omega=460.659`, (7) gives `log10(rho_1)=-1151.861`. At
`E=2 Gamma=2666.667` it gives `-1124.338`. These are **leading matched
estimates**, not error-certified physical numbers. Taking `k` large also
changes `M` and `E`; fixed-parameter low-frequency asymptotics alone do
not justify that simultaneous limit. In particular the tiny absorptive
part cannot be inferred by numerically subtracting two approximate real
Green functions.

## Exact conditional propagation at finite real frequency

A bound can avoid subtracting real parts. Fix an outer cut
`u_c=1-delta`, with `0<delta<1`. Assume for a particular inner retarded
solution that

```
Im r_s(u_c)>=0,    |r_s(u_c)|<=R<1,
|omega|(1+R^2) <= 2 R W_min,
W_min = (lambda_n/b) min[sqrt(X), delta sqrt(H_X(u_c))]. (8)
```

The value of `W_min` is exact: the only interior stationary point of
`(1-u)^2 H_X(u)` is a maximum, since its derivative is
`2(1-u)u[1-2(1+3X)u]`. Thus the minimum on the retained exterior is
at an endpoint.

From (1), along increasing `x`,

```
r_s' = omega(1+r_s^2)-2W r_s,
(Im r_s)' = 2[omega Re r_s-W] Im r_s.
```

At the boundary of the disk `|r_s|<=R`, its outward radial derivative
is at most `|omega|(1+R^2)-2W R<=0`. The disk is invariant. Consequently
there is no exterior denominator zero, and positive flux propagates as

```
Im r_s(0) = Im r_s(x_c) exp[-2 S_n
                                  +2 omega integral_(x_c)^0 Re r_s dx],
S_n=lambda_n I_X(u_c),
tau_out=-x_c <= (b/X)(1/delta-1).
```

Combining this identity with (2) proves the finite-frequency bound

```
rho_n <= 2R/(1-R)^2 exp[-2 S_n+2|omega|R tau_out].    (9)
```

This is an exact ordinary bound on the frozen extremal exterior,
conditional on (8). It is not a WKB estimate, nor does it infer (8)
from causality alone. A general passive inner load has positive imaginary
part but can have arbitrarily large modulus near a resonance. A reflecting
evolving throat cannot simply be replaced by an absorbing RN horizon.
The small finite-temperature deformation of the exterior also requires
its own comparison before (9) is used on the actual initial RN solution.

There is also a bound on the complete higher-level boundary spectral
trace, if the same inner-cut disk condition holds for **every** level in
the specified frequency band. Since `lambda_n>=sqrt(n) lambda_1`, put
`a=2 S_1>=4` and `epsilon=2|omega|R tau_out`. Then

```
(1/q) sum_(n>=1) d_n rho_n
 <= 4R/(1-R)^2 exp(-a+epsilon) B(a,q),
B(a,q) = 1+2/q+2(1/a+1/a^2)
                +(4/q)(1/a+3/a^2+6/a^3+6/a^4).      (10)
```

To prove it, the decreasing function `(1+2x/q) exp(-a sqrt(x))`
bounds the sum by its first term plus its integral from one to infinity.
The substitution `v=sqrt(x)` gives the elementary integrals in (10).
Degeneracy therefore does not defeat the radial barrier under the stated
common load condition. This is a bound on the boundary spectral trace,
not a stress tensor or a bilinear noise convolution.

For the existing candidate choose `delta=k^(-4)`, `R=0.1`, and the band
`|E|<=2 Gamma`. At `k=2000` the exterior frequency condition has
`|omega|/W_min=2.257e-18` at worst. The logarithmic frequency correction
in (9) is at most `4.887e-17`. Conditional on all the inner load bounds,
the full tower trace per LLL channel obeys
`log10[(sum d_n rho_n)/q] <= -492.4778`.
The corresponding bounds at `k=4000,8000` are in the receipt.
The numerical margin makes checking the inner condition worthwhile;
it does not establish that condition. The leading inner small-argument
solution suggests `Re r_s=omega zeta/(2M-1)+...`, comfortably within the
disk here, but an error bound uniform in the retained modes and evolving
state remains to be supplied.

## What this does and does not resolve

The preceding mass comparison did not by itself establish a large
boundary load. The exterior matching calculation now identifies a
strong suppression mechanism, and (8)–(10) reduce a finite-frequency
boundary estimate to an explicit inner response condition. That is
progress toward testing the candidate, without accepting its formation.

Two qualifications matter for the next calculation. First, a boundary
bilinear response involves spectral convolutions, state occupations and
the actual vertex. A smooth finite-time switch is not band limited;
thermal and switching tails, local counterterms, and the matched value
of the coupling still need bounds. The trace in (10) is not that work.
Second, a collective electric or angular field can excite massive
fermions **inside** the throat. Its bulk polarization loop has no need
for each fermion to propagate to the AdS boundary. Neither (6) nor (10)
bounds that contribution. Renormalized massive-mode polarization at the
collective scale is therefore a separate immediate calculation.

[Checker](../checks/higher_landau_boundary.py),
[receipt](../evidence/higher-landau-boundary.json), and
[source/dependency record](../evidence/higher-landau-boundary-source-review.json)
check the radial integral, independent Bessel flux normalization, thermal
identity, exact finite-frequency propagation identity, and tower bound.
Numerical checks support the stated ordinary arguments; they do not
replace the open inner-load, state or gravitational matching proofs.

Remaining to-do list: bound bulk massive-mode polarization and verify the
inner-load condition through the evolving geometry; control spectral tails,
quantum noise, mixed energy terms and higher-order interactions; match
physical clocks, vertices, states and complete stress; establish preparation,
sufficient stability and finite-energy transport.
