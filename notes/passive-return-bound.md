# Passive return bounds and throat transit time

2026-10-10. S051/S052/S054 conditional obstruction. In the leading MMP length
energy, a causal passive Gaussian return cannot give a stationary throat
whose transit time is at or below the earliest exterior arrival time.
The bound permits frequency-dependent mixing, unequal additional delays
and contraction from loss. It has an explicit margin for omitted energy
derivatives. The physical reduction of the complete four-dimensional
quantum energy remains unproved. The argument extends MMP's existing
shortest-field-line comparison using established passive-scattering methods;
its novelty and suitability as a standalone contribution remain unresolved.

## Defined model and comparison times

Use units `hbar=c=1`. Let `T>0` be the throat optical transit time and `D>0`
the earliest exterior arrival between the same port surfaces, measured by
the same time coordinate. `D` is a signal front time, not a group delay or
a scattering phase derivative. It is held fixed when varying `T`.

Take `q` complex massless Dirac channels. Each of their two chiral sectors
has a `q` by `q` exterior return `R_sigma(s)`, analytic for `Re s>0`, with
the untwisted antiperiodic minus sign separated out. Impose

```text
R_sigma(s) = exp(-s D) G_sigma(s),
||G_sigma(s)||_op <= 1,                    Re s > 0.       (1)
```

No normality, reciprocity, simultaneous diagonalization or frequency
independence is assumed. Loss is admissible in this defined Gaussian
ground-state problem only with its passive environmental completion;
replacing an interacting quantum state by a lossy amplitude is not enough.
Any constant Wilson/spin phase can be absorbed in the return without
changing (1). The returns and all their internal parameters are independent
of `T` in the variation below.

The length-dependent energy functional is

```text
E(T) = A/T^2 + q pi/(24 T) + E_loop(T),       A>0,

E_loop(T) = -(1/pi) sum_{sigma=+,-} integral_0^infty
                Re log det[I + exp(-xi T) R_sigma(xi)] dxi. (2)
```

The loop energy subtracts the disconnected exterior self energies and the
same-regulator infinite-line contribution. Equation (2) specifies the
Gaussian return model being studied; a causal response matrix alone does
not establish that a given four-dimensional quantum state has this energy.
For two identical sectors it reduces to the determinant checked in the
[equal-length phase calculation](lll-return-holonomy.md).

For `R_sigma=exp(-s D) I`, the loop term is `-q pi/[6(T+D)]`. With
`T=pi ell` and `A=pi^2 r_e^3/G_N`, (2) is exactly MMP's shortest-return
comparison in equations (5.44)-(5.45). The positive `q pi/(24T)` term is
their throat conformal-anomaly/reference contribution. Changing its
coefficient, adding other matter or changing the classical energy changes
the conclusion. This is a specified leading-order model.

## Why a physical front and passivity give (1)

For completeness, consider a time-translation-invariant linear return on
flux-normalized amplitudes with convolution kernel zero before time `D`.
Assume this map is a contraction on the full time-domain `L2` space and its
kernel admits the Laplace transform for `Re s>0`. Advance the output by `D`.
The resulting map `G` remains causal, and the full-line time shift preserves
its `L2` norm. This step would fail without a genuine front bound.

For `Re s>0`, input `f(t)=exp(s t) v` for `t<=0` and zero for `t>0`.
Causality gives `Gf(t)=exp(s t) G(s)v` for `t<=0`. Contractivity implies

```text
||G(s)v||^2/(2 Re s)
   = integral_{-infinity}^0 ||Gf(t)||^2 dt
  <= ||Gf||_L2^2 <= ||f||_L2^2 = ||v||^2/(2 Re s).
```

This proves (1). In a scattering application, both the flux normalization
and the causal port reference planes must be justified. Unitarity on real
frequencies alone, without analytic growth and front control, is insufficient.
For example `exp(+s D)` is unit modulus on the imaginary axis but is an
advance. Conventional far-field channel conventions can contain reference
advances and cannot be substituted for a physical port-to-port map.

## Sharp matrix derivative bound

For an eigenvalue `z=t exp(i theta)` with `0<=t<=r<1`,

```text
-r/(1-r) <= Re[z/(1+z)] <= r/(1+r).                      (3)
```

The lower bound follows from `|z/(1+z)|<=r/(1-r)`. For the upper bound,
multiply the difference by `(1+r)|1+z|^2`. Its numerator is

```text
r + (r-1)t cos(theta) - t^2
  = (r-t)(1+t) + (1-r)t[1-cos(theta)] >= 0.
```

At imaginary frequency define `X_sigma=exp(-xi T) R_sigma(xi)` and
`r=exp[-xi(T+D)]`. The norm bound implies that every eigenvalue of `X_sigma`
lies in the disk of radius `r`. Analytic functional calculus gives
`Tr[X_sigma(I+X_sigma)^(-1)]` as the sum of `z/(1+z)`, including algebraic
multiplicities. This remains true for nondiagonalizable matrices.

Differentiating (2) at fixed return and applying (3) gives

```text
-q pi/[3(T+D)^2] <= E_loop'(T) <= q pi/[6(T+D)^2].         (4)
```

Here `integral_0^infty u/(exp(u)+1) du = pi^2/12` and the corresponding
minus-sign integral is `pi^2/6`. They follow by expanding the geometric
series and integrating, using the convergent sum of reciprocal squares.
The logarithm is well-defined for `xi>0`; its real part is `log|det|`.
Near zero its possible logarithmic divergence is integrable. A dominating
derivative is proportional to `xi/[exp(xi(T+D))-1]`, justifying the
derivative on any compact positive `T` interval.

Both constants are sharp. The upper bound is attained by the ideal
antiperiodic return, the lower by its periodic counterpart
`R_sigma=-exp(-sD)I`. This is a one-sided upper limit on the support-driving
derivative, not a claim that every fermionic Casimir force has one sign.

## No stationary shortcut and its margin

Combining (2) and (4),

```text
E'(T) <= -2A/T^3 - q pi/(24T^2) + q pi/[6(T+D)^2].        (5)
```

For `0<T<=D`, the right side equals

```text
-2A/T^3 - q pi (D-T)(D+3T)/[24 T^2 (T+D)^2] < 0.        (6)
```

Consequently no stationary point, including a local minimum, exists with
`T<=D`. The energy derivative is uniformly at most `-2A/D^3` throughout
that interval. Thus the same exclusion survives an omitted contribution
`Delta E(T)` whenever `Delta E'(T) <= epsilon < 2A/D^3` there. This is an
explicit error budget; no actual four-dimensional remainder is bounded here.

A stronger necessary condition follows by defining `d0=16A/(q pi)`.
Every stationary point satisfies

```text
4 T^3/(T+D)^2 - T >= 3 d0.                              (7)
```

On `T>=D`, the left side starts at zero and increases strictly to infinity:
its derivative multiplied by `(T+D)^3`, with `T=D+x`, is
`8D^3+24D^2 x+18D x^2+3x^3>0`. There is a unique `T_min>D` where equality
holds. Thus `T>=T_min`, and also `T>d0` for `D>0`. The ideal return attains
this threshold, so mixing or passive dispersion cannot beat it within (2).
The zero-exterior-length limit gives the familiar MMP value `T=d0`.

The same no-shortcut conclusion holds if a common Wilson phase is minimized:
the bound is uniform in that phase. Integrating the uniform negative
derivative between two `T` values and then minimizing preserves strict
decrease, even where the minimizing phase changes branches.

## Physical applicability and prior comparison

[MMP section 3.1 and appendices D-E](https://arxiv.org/html/1807.04726v3#S3.SS1)
already use adiabatic magnetic field-line confinement, identify exceptional
escaping modes and give unequal exterior lengths. Section 5.5 already
excludes an external shortcut using the shortest-return comparison. Its
equations (5.44)-(5.47) and conclusion agree with printed pages 18-20 of the
[2023 accepted manuscript hosted by OSTI](https://www.osti.gov/servlets/purl/2419919)
(PDF pages 20-22). This is an accepted-manuscript comparison, not a complete
comparison with the publisher's typeset version or all proof dependencies.

[Lambrecht, Jaekel and Reynaud](https://arxiv.org/html/quant-ph/9801055v1)
sections 2-4 derive passive scattering contraction and a perfect-mirror
Casimir-force bound, including frequency dependence. That method is an
explicit predecessor of (4). [Beenakker](https://arxiv.org/html/2402.02477v2#S2)
section II gives the fermionic scattering-density and imaginary-frequency
determinant prescription. The circle normalization in (2) is independently
checked against the known Dirac spectrum; it is not copied from a two-wall
geometry with a different round-trip length.

[Han et al.](https://arxiv.org/abs/2608.10427v1) discuss channel reference
advances and conditional Schur representations. The selected main-text
discussion is read; its cited supplementary proof is not audited. The
time-domain argument above supplies the needed restricted result directly
and does not import their sum rules. [Aoki and Fukaya](https://arxiv.org/html/2203.03782v3#S3.SS1)
derive the spin-connection phase for circular domain-wall modes. That
precedent helps assess an adiabatic phase argument but does not determine
the open two-mouth exterior return or its gluing.

[Fernandez, Fosco and Hansen](https://arxiv.org/html/2505.13682v2#S3)
give a Dirac force-sign argument under reflection and positivity conditions.
Its general sign claim is not imported here. The inference after their
(51), from positivity of a difference of positive Laplace transforms to
its decreasing with separation, requires further justification: take
`I_+(a)=exp(-a)` and `I_-(a)=(3/4)exp(-2a)`. At `a=log(4/3)` their difference
is `21/64>0` and its derivative is `3/32>0`. This refutes that generic
inference, not the physical theorem: the example does not realize their
specific operator/eigenvectors. Edition comparison and any additional
constraints that repair the step remain to review. Bound (4) needs neither
reflection symmetry nor a universal force sign.

To apply (5) to the physical geometry one must still show that the projected
return has the stated front, norm and Gaussian energy, and that varying
the throat holds that exterior fixed to the required accuracy. Same-mouth
reflection, channel projections with nonlocal time dependence, extra light
degrees of freedom, singlet/gauge interactions, background radius shifts
and exterior backreaction can require a different determinant or additional
`T` derivatives. If a shorter return path exists, `D` must be reduced to
its true front. An assumption failure alone does not establish a shortcut.
Formation, nonlinear stability and the existence of any stationary solution
are separate questions.

Decision: retain the sharp conditional bound as a mathematical endpoint
for contribution review. In this model, abandon passive mixing or dispersion
as a route to beating the exterior travel time. The next physical test is
whether the free-sector determinant and its positive error budget survive
the actual exterior and omitted quantum sectors. MMP already supplies the
qualitative no-shortcut result and passive-force bounds are established;
the full novelty exit and primary/fallback selection remain open. Do not
start a large Lean campaign or a manuscript claiming a new physical no-go
until that comparison and physical consequence are resolved.

## Independent verification

[The executable check](../checks/passive_return_bound.py) and
[receipt](../evidence/passive-return-bound.json) check five exact algebraic
identities, 90 pointwise matrix cases, 18 integrated cases and 18 independent
energy variations. The matrices include nonnormal, frequency-noncommuting,
lossy and resonant returns. Both sharp endpoints reproduce the circle result.

A separate real-spectrum calculation uses
`R(s)=exp(-sD)(s-a)/(s+a)`. Its positive roots satisfy
`omega(T+D)+2 atan(omega/a)=2 pi n`. Differentiating these roots gives
`omega'=-omega/[T+D+2a/(a^2+omega^2)]`. Summing the regulated vacuum
derivative and subtracting the identical-regulator continuum density gives
an independent estimate of (4)'s integrand integral. The latter density is
`(T+D)/(2pi)+a/[pi(a^2+omega^2)]`; its isolated-resonance part is independent
of `T`. Successive regulator extrapolations agree in three cases, with
refined residual at most `2.2e-9`. The zero mode contributes no energy.
Numerical cutoff and extrapolation errors are not certified. The ordinary
proof establishes the bound under (1)-(2); these checks test its formulas
and normalizations. No Lean invocation or physical novelty acceptance.

Remaining to-do list: audit the physical return and energy reduction, bound
omitted derivatives against the margin, complete decisive source comparisons,
and select the primary result and fallback before certification and manuscript preparation.
