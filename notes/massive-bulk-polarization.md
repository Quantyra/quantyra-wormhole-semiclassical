# Massive bulk polarization: subtraction, thresholds and finite pulse work

S072, 2026-10-13. **The complete higher-Landau electric bubble admits a
finite, positive, all-frequency pulse-work bound after the static charge
coefficient is matched.** Individual real-frequency thresholds are
unbounded, so a uniform pointwise smallness argument fails. An integrated
bound avoids that problem. This result is for a free fermion vacuum on
`R^(1,1) x S2(b)` with a prescribed classical electric source. The sphere
spectrum is exact; replacing the longitudinal flat space by the actual
AdS2/RN throat, and the source by the interacting quantum phase, remains
open. The result is not a physical capture or stability theorem.

The [preceding boundary calculation](higher-landau-boundary.md) does not
control this bulk loop: both fermions can stay inside the throat. This
note computes that distinct contribution without using the exterior
barrier. The [collective electric comparison](collective-electric-response.md)
also needs care: temporal frequency alone is not invariant pair energy.

## Exact product reduction and Ward normalization

Take one massless, unit-charge four-dimensional Dirac field with integer
magnetic flux `q>=1`. A spherically constant potential `A_a(t,x)` couples
diagonally to the sphere harmonics. The resulting two-dimensional Dirac
channels have

```
m_n^2=n(n+q)/b^2,     d_0=q,     d_n=2(q+2n), n>=1.
```

Both nonzero sphere eigenvalue signs contribute the same quadratic vector
response. Their opposite pseudoscalar masses can be rotated to ordinary
masses; possible linear topological terms do not alter this bubble.
Nonspherical gauge fields mix harmonics and are outside this calculation.

In Euclidean longitudinal momentum `z=p_E^2>=0`, the Ward-preserving
one-channel determinant gives the electric-field coefficient

```
c_m(z) = (1/pi) integral_0^1 dx v/[m^2+z v],
v=x(1-x),
Pi_ab=(p_E^2 delta_ab-p_a p_b)c_m(z),
Gamma_E^(2)=(1/2) integral F(-p)c_m(z)F(p).             (1)
```

Here the potential has unit-charge normalization; the coupling resides
in the Maxwell coefficient `kappa=4pi b^2/g^2`. The massive limit is
`c_m(0)=1/(6pi m^2)`. The massless channels instead have
`z c_0(z)=1/pi`, supplying the already retained `q/pi` Schwinger term.
They must not be included again in the massive sum.

Equation (1) follows by combining the two fermion denominators and taking
the two-component trace with a vector Ward regulator. The finite contact
from the regulated divergent tensor integral is essential: dropping it
would violate transversality. The independent reference is
[Hattori and Itakura, appendix A.1, equation 91 and section II.2](https://arxiv.org/html/2205.04312v2).
Only its vacuum two-dimensional loop is imported. The planar Landau
degeneracy and its finite-temperature medium are replaced here by the
stated sphere spectrum and vacuum; there is no import of their medium
dispersion relation.

## Static renormalization and the finite frequency difference

Every massive channel has a finite coefficient, but their static sum does
not. With a common sphere cutoff `1<=n<=N`,

```
c_H,N(0) = b^2/(3pi) [H_N+H_(N+q)-H_q].             (2)
```

This follows from
`(q+2n)/[n(n+q)]=1/n+1/(n+q)`. With physical cutoff
`Lambda_N=sqrt[N(N+q)]/b`, subtract `2b^2 log(Lambda_N/mu)/(3pi)`.
The finite static coefficient in this explicitly specified spectral
scheme is

```
c_H,ren(0;mu)=b^2/(3pi)[2 log(b mu)+2 gamma_E-H_q],
kappa_* = 4pi b^2/g_R(mu)^2+c_H,ren(0;mu),
g_*^2=4pi b^2/kappa_*.                               (3)
```

Assume `kappa_*>0`. The subtraction uses one common local Maxwell
counterterm, not a separate subtraction chosen for each Landau level.
Its scale cancellation requires
`d(1/g_R^2)/d log mu=-1/(6pi^2)`. The finite magnetic part in (3)
is retained. Defining `g_*` as the measured electric coefficient at this
fixed product background does not identify it automatically with the
Maxwell coupling in a different geometry. A finite scheme change must
be accompanied by the corresponding change in `g_R`.

The logarithmic all-level divergence and the distinction between it and
finite magnetic response are consistent with
[Hattori and Itakura, appendix B](https://arxiv.org/html/1209.2663v2).
The earlier [curved product vacuum](magnetic-product-vacuum.md) uses its
own stated gravitational/Maxwell subtraction; (3) is not silently equated
with that prescription or with the RN response.

The frequency-dependent difference is independent of the local constant:

```
Delta c_H(z) = sum_(n>=1) d_n [c_mn(z)-c_mn(0)]
 = -(z/pi) sum_(n>=1) d_n integral_0^1
                v^2/[m_n^2(m_n^2+z v)] dx.           (4)
```

It converges absolutely. A telescoping sphere identity gives

```
sum_(n>=1) d_n/m_n^4 = (2b^4/q) H_q^(2),
0 <= -Delta c_H(z)/kappa_*
   <= g_*^2 z b^2 H_q^(2)/(60pi^2 q).                (5)
```

Indeed `(q+2n)/[n^2(n+q)^2]=[1/n^2-1/(n+q)^2]/q`, and
`integral_0^1 v^2 dx=1/30`. This bound covers the full tower without
assigning it a hard ultraviolet validity cutoff. It is a statement about
the specified renormalized one-loop model, not about its UV completion.

There is also an exact convergent digamma expression, useful when the
bound in (5) is loose. Let `r_++r_-=q`, `r_+r_-=b^2 z v`. Then

```
Delta c_H(z)=(2b^2/pi) integral_0^1 dx v [
  psi(1)+psi(q+1)-psi(1+r_-)-psi(1+r_+)].             (6)
```

Above the quadratic discriminant crossing the two roots are conjugate,
and their digamma sum is real. The expression is symmetric in the roots.
The checker compares (6) with a direct finite sum plus its analytic tail,
rather than using the digamma expression as its own verification.

For numerical comparison only, set `a=b/sqrt(C)`, `C=1+6X`, using the
existing RN parameters as a length normalization, and write `Z=z a^2`.
Set `g_R(1/b)=k^-4` in the scheme (3); this scale convention is part of
the comparison. At `k=2000` the relative static correction is of order
`10^-27`, while (6) gives

```
-Delta c_H/kappa_* = 1.77918e-28 at Z=Omega^2,
                   = 2.93600e-28 at Z=(2 Gamma)^2.    (7)
```

These are **Euclidean invariant momentum** values. They do not identify
an inhomogeneous phase drive at temporal frequency `Omega` with zero
radial momentum, or imply a uniform bound on the real axis.

## Retarded thresholds: a pointwise bound fails

For positive temporal frequency continue `z=-(s+i0)`, where
`s=omega^2-p_x^2`. The absorptive coefficient of one massive channel is

```
Im c_m,R(s) = 2m^2/[s^2 sqrt(1-4m^2/s)] theta(s-4m^2). (8)
```

The static subtraction does not change it. One derivation applies the
delta function to (1). An independent pair-state derivation uses
`h(p)=p sigma_3+m sigma_1`: the squared interband matrix element of
`j_x=sigma_3` is `m^2/E_p^2`. The two roots of `omega=2E_p`, with their
Jacobian, give `Im chi_jj=2m^2/[s sqrt(1-4m^2/s)]` in the pair rest frame;
division by `s` reproduces (8). These separate normalizations agree in
the checker.

There is no massive-pair absorption below invariant threshold. Above it,
the square-root singularity means that for every nonzero coupling the
pointwise ratio `Im c_H,R/kappa_*` is unbounded. Small `g^2` does not
justify an operator bound uniform over all real invariant frequencies.
This is the specific failed approximation, not a wormhole instability.
It agrees with the threshold analysis in
[Hattori and Itakura, sections IV.2–IV.3](https://arxiv.org/html/1209.2663v2).

## A positive cumulative measure and an all-frequency work bound

The singularities are integrable. Equation (8) has the exact primitive

```
R(S) := integral_0^S Im c_H,R(s) ds
      = sum_(4m_n^2<=S) d_n sqrt(1-4m_n^2/S)
      <= (b^2 S/2)(1+1/q).                           (9)
```

To prove the last inequality let `N` be the largest allowed level.
Dropping the square roots gives `2N(q+N+1)`. The threshold condition
`4N(N+q)<=b^2 S` and `N<=b^2 S/(4q)` give (9). This exact counting
bound includes arbitrarily high levels and needs no separation between
the frequency of the source and the masses. In particular,
`R(S)/(kappa_* S)<=g_*^2(1+1/q)/(8pi)`.

Use Fourier convention `E(t,x)=integral d omega d p/(2pi)^2
E_tilde(omega,p) exp(-i omega t+i p x)`. Let the real classical field be
off at early and late times. At quadratic order its work deposited in
the massive vacuum channels is

```
W_H = integral dp/(2pi) integral_0^infinity d omega/pi
          omega Im c_H,R(omega^2-p^2) |E_tilde(omega,p)|^2 >=0. (10)
```

This is dissipative pair energy after the pulse. The real static
coefficient returns its reversible energy when the pulse is off. The
massless LLL, other fields, higher source orders and a dynamical source's
own quantum energy are not included in `W_H`.

For each spatial Fourier component define
`A0(p)=||E(t,p)||_L1(dt)` and `A2(p)=||partial_t^2 E(t,p)||_L1(dt)`.
Assume the derivatives have no endpoint boundary term and
`integral dp A0(p)A2(p)<infinity`. Integration by parts gives
`|E_tilde|<=min(A0,A2/omega^2)`. On the absorptive support,
`omega^2=s+p^2>=s`, so its square is bounded by
`w(s)=min(A0^2,A2^2/s^2)`. Changing variables in (10) gives

```
W_H <= (1/2pi) integral dp/(2pi) integral_0^infinity w(s) dR(s)
    <= b^2/(2pi) (1+1/q) integral dp/(2pi) A0(p) A2(p). (11)
```

For the second step integrate the positive measure by parts, use
`R(s)<=b^2(1+1/q)s/2` and
`integral_0^infinity w(s) ds=2 A0 A2`. This proves a bound on the whole
frequency axis, including all thresholds and the ultraviolet tail. It
does not replace the retarded kernel by its zero-frequency value.

An explicit smooth compact pulse is

```
E(t,x)=E0(x) [1-(t/tau)^2]^3 for |t|<=tau, otherwise zero,
U_peak=(kappa_*/2) integral dx E0(x)^2.
```

The shape has `||f||_1=32/35` and `||f''||_1=384/(25sqrt(5))`.
Parseval in (11) gives, for any square-integrable spatial profile,

```
W_H/U_peak <= 3072/[875sqrt(5)pi^2] g_*^2(1+1/q)
            = 0.159... g_*^2(1+1/q).                 (12)
```

The bound is uniform in positive duration `tau`. At `k=2000` the
specified flat-product comparison gives `6.21425e-28`. This denominator
is the pulse's Maxwell electric energy, **not** the negative wormhole
opening energy. The actual S072 source is an interacting quantum field;
neither its profile nor its energy is supplied by this illustrative pulse.

## The light pole also uses invariant mass

Retain the already included LLL Schwinger term and resum the massive
quadratic bubble. In this flat vacuum, its positive light mass squared
solves

```
s_* [kappa_*+Delta c_H,R(s_*)]=q/pi,
s_0=q/(pi kappa_*),
r_0=s_0/(4m_1^2)=g_*^2 q/[16pi^2(q+1)].              (13)
```

If `r_0<1`, positivity and monotonicity of the subthreshold integral
give one root `0<s_*<s_0`. On `0<=s<=s_0`, the analogue of (5) has an
extra factor `1/(1-r_0)` and positive sign. Consequently

```
0 <= 1-s_*/s_0
   <= g_*^4 H_q^(2)/[240pi^4(1-r_0)].                (14)
```

There is no higher-mode absorptive width at this subthreshold pole in
the specified one-loop vacuum model. At `k=2000`, `r_0=2.47366e-29`
and the fractional shift bound is `1.07364e-57`. These extremely small
numbers do not bound higher-loop or geometry errors. A propagating mode
can have large temporal frequency and almost equal spatial momentum
while retaining this small invariant mass. Conversely, an externally
forced field can have timelike off-shell components that reach (8).
Neither statement identifies the collective boundary phase with a free
bulk photon. Its response must be obtained with the radial constraints.

## Scientific decision and verification

The feared conclusion from a large temporal charge rate is not established
by the massive thresholds. Their actual invariant response is now explicit,
with a common renormalization, a controlled light pole and an integrated
work estimate. A uniform pointwise real-axis loop bound is rejected;
the positive cumulative spectral estimate is the useful replacement.

The necessary next transfer must retain AdS2 curvature, the initial
thermal/entangled state and the evolving radial geometry. Flat-space
translation invariance made (8)–(11) possible; it is not present globally
in the RN problem. A large `m_1 a` is a comparison parameter, not a proved
curvature remainder. The physical source is quantum and self-consistent,
so replacing it by the prescribed pulse would omit precisely the mixed
current/noise work already identified in [the phase action](phase-influence-action.md).
Angular and nonspherical response are additional sectors.

[Checker](../checks/massive_bulk_polarization.py),
[receipt](../evidence/massive-bulk-polarization.json), and
[source/dependency record](../evidence/massive-bulk-polarization-source-review.json)
verify the subtraction, sum identities, independent pair cut, cumulative
measure, explicit pulse and subthreshold pole. The pulse integration uses
the exact transform `F(y)=96[(15-6y^2)sin(y)+(y^3-15y)cos(y)]/y^7`.
For `y>=Y`, `|F(y)|<=D(Y)/y^4`, where
`D(Y)=96[1+6/Y+15/Y^2+15/Y^3]`. In the numerical units `b=E0=1`,
the uncomputed work beyond level `N` is bounded by
`D(Y)^2/[1890pi tau^6 N^6]`, with `Y=2m_(N+1)tau`;
this follows from `sum_(n>N) d_n/m_n^8<=2/(3N^6)` and
`integral_0^1 (1-t^2)^4 dt=128/315`. The tail is separate from
quadrature refinement. No curved-throat error,
all-order quantum switching stress, physical preparation or stable
creation claim follows from the checks.

Remaining to-do list: transfer the massive response and integrated work
bound to the curved, evolving throat and actual quantum source; verify
inner matching and include noise, mixed work and complete stress; establish
physical state/clock/vertex matching, preparation, sufficient stability
and finite-energy transport.
