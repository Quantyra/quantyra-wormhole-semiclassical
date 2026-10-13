# Cancelling electric returns and the minimum switch duration

S072, E003, 2026-10-13. **A nonzero switch shorter than one strip return
period cannot leave the quadratic restoring response stationary. A
monotone switch lasting just over that period can cancel every restoring
return in the fixed-strip model.** After that longer switch the positive
quadratic Hamiltonian is autonomous, so its regulated energy is conserved.
This removes one source of parametric driving without assuming that the
switched quantum state becomes the ground state.

The longer switch lies outside the existing pre-return JT capture proof.
These results therefore give a control option and a duration constraint,
not a completed physical opening or a full stability theorem. The rapid
candidate remains selected; its nonstationary stability has not been
replaced by the stationary theorem for a different protocol.

## Model and the reviewed stationary result

Use the exact finite-channel electric-neutral vacuum spectra from
[neutral-transfer-current.md](neutral-transfer-current.md), with angular
dynamics fixed. In the common global strip clock the restoring spectrum
has frequencies `E=2,4,...` and strictly positive weights `v_E(q)` for
`q>2`. The [quadratic feedback proof](passive-electric-feedback.md) gives

```
mu=q/(2pi),  a_E=2v_E/E,  0<a_E<=2q/pi^2,
k(t)=mu h(t)+h(t)^2 delta_q-h(t) R_q(t),
R_q(t)=sum_even a_E integral_0^t cos[E(t-u)] h'(u) du,
delta_q=chi_XX(0)-chi_YY(0)>0.                         (1)
```

The static difference uses a common Abel regulator before taking the
limit. Indeed its regulated value is `2 integral_epsilon^infinity
(C_XX-C_YY) d tau`; positivity and integrability of the difference for
`q>2` justify the limit. The digamma expression in the preceding note
agrees with an independent integral. A paired hard truncation is not the
same subtraction. At finite common regulator the exact algebra uses
`delta_reg`; the finite-band demonstrations match their local restoring
coefficient to `delta_q`, rather than asserting that their raw truncated
static difference already equals it.

The positive-real proof excludes stationary growing linear poles for
any of the accepted positive radial measures in the same clock and
normalization. The all-return cubic stiffness bound also survives
review. Neither proof supplies the physical clock, vertex, continuum
state, stress subtraction or nonlinear remainder. The existing checker
was rerun unchanged and reproduced its receipt byte for byte.

## Why the fast switch keeps ringing

Let real `f=h'` be supported in `[0,T]`, with integral `h_f!=0`, and define

```
F(E)=integral_0^T f(u) exp(-i E u) du.
R_q(t)=sum_(n>=1) a_(2n) Re[exp(i2nt) F(2n)], t>=T.    (2)
```

For example `f` may be continuous, vanish at its support endpoints, and
have derivative of bounded variation. Then `F(E)=O(E^-2)`, so (2)
converges absolutely and uniformly. It is periodic with period `pi`,
has zero mean, and can tend to zero at late times only if it is identically
zero. Positivity of every `a_(2n)` makes that equivalent to

```
F(2n)=0 for every n>=1.                                (3)
```

Periodize `f` with period `pi`. Its Fourier coefficients are `F(2n)/pi`.
If (3) holds, reality and Fourier uniqueness in L1 imply that the
periodization equals `h_f/pi` almost everywhere. If `T<pi`, it instead
vanishes on a nonempty interval. This is a contradiction. The argument
also covers signed switches; monotonicity is not needed. For `T=pi`,
the only possible derivative is `h_f/pi` almost everywhere, which is
incompatible with continuous zero endpoint derivatives. Smooth switches
with complete cancellation therefore require `T>pi`.

There is a useful quantitative version for monotone switches. If
`f>=0` and `T<pi/2`, projecting `exp(iT)F(2)` onto its real axis gives

```
|F(2)| >= h_f cos T,
v_2=2(q-1)/pi^2,
first harmonic amplitude of [k-k_stat]/(mu h_f)
  = a_2 |F(2)|/mu >= (4/pi)(1-1/q)h_f cos T.           (4)
```

The equality for `v_2` follows directly from the exact generating
coefficient `f_2(a)=8a^2`. The first harmonic cannot be canceled by
higher harmonics. In particular the root mean square of the normalized
stiffness oscillation is at least the right side of (4) divided by
`sqrt(2)`. This is a coefficient bound, not a growth exponent.

For the existing cubic switch,

```
|F(2)|/h_f = 3(sin T-T cos T)/T^3.
```

At `q=2000^4`, `h_f=4pi/2000`, `T=0.25`, the first normalized stiffness
harmonic is about `0.0079501`; every monotone switch of that duration
has the lower bound `0.0077513`. It persists in this reference even
though the earlier uniform instantaneous stiffness ratio exceeds
`0.9246`. This does not establish a parametric instability.

## A switch that cancels the entire even spectrum

Let `p_tau` be any nonnegative normalized smooth pulse supported in
`[0,tau]`. Average it over one full period:

```
f(t)=h_f (B*p_tau)(t),  B(t)=1_[0,pi](t)/pi,
h(t)=integral_0^t f(u) du,  T_star=pi+tau,
F(E)=h_f (1-exp(-i pi E))/(i pi E) p_tilde(E).          (5)
```

This is monotone and has final value `h_f`. The box transform vanishes
at every positive even frequency, independently of `q` and of the
weights. Therefore `R_q(t)=0` exactly for `t>=T_star`, including the
complete infinite restoring spectrum. Cancellation concerns the
coherent `X` restoring response; it does not remove the `Y` bath noise,
radial fluctuations, quantum correlations or excitations.

The explicit checked pulse is

```
p_tau(t)=30 x^2(1-x)^2/tau,  x=t/tau in [0,1],
S(x)=10x^3-15x^4+6x^5, extended by 0 and 1,
f(t)=h_f/pi [S(t/tau)-S((t-pi)/tau)].                   (6)
```

Here `p_tau` is C1 after extension, `f` is C2, and `h` is C3; a C-infinity
bump may be used in (5) instead. No instantaneous quench is needed.

For `0<tau<pi`, positivity during the whole longer ramp follows from
the same all-frequency comparison used in the previous note. With
`Z_E=integral_0^t cos[E(t-u)]f(u)du`, two integrations by parts give

```
|Z_E(t)| <= A/E^2,  A=75 h_f/(8pi tau),
|R_q-R_infinity| <= (16 zeta(3)/pi^2) A.                (7)
```

To see the first constant, `||f'||_infinity=15h_f/(8pi tau)`,
`f'(0)=0`, and `||f''||_1=15h_f/(2pi tau)`. The comparison spectrum
has the exact Dirac-comb expression

```
R_infinity/mu = f(t)+2 sum_(m>=1)f(t-m pi)-2h(t)/pi.
```

The periodized nonnegative `f` in (5) equals `h_f/pi`, so the sum with
coefficients at most two is bounded by `2h_f/pi`. Using `delta_q>0`
in (1) proves

```
k(t)/(mu h(t)) >= 1-2h_f/pi-300 zeta(3)h_f/(pi^2 q tau)
                                                     (8)
```

for every `t>0`. For the same `q=k^4`, `h_f=4pi/k`, with `tau=0.25`,
this lower bound is greater than `0.9959999999999` for every integer
`k>=2000`. The final Hamiltonian has `k_stat=mu h_f+h_f^2 delta_q>0`.

## What final energy conservation proves

At every finite oscillator regulator with the matched static term,
the Hamiltonian after `T_star` is the positive completed-square
Hamiltonian of the preceding note with constant coefficients. Thus
for any evolved state with finite regulated energy,

```
<H_quad(t)>=<H_quad(T_star)>,
<alpha(t)^2> <= 2<H_quad(T_star)>/k_stat, t>=T_star.     (9)
```

The classical finite-dimensional statement is identical. Every normal
frequency is real and all solutions are bounded. Correlations created
during the switch are retained; no reset to a factorized or stationary
state is used. The constants in (9) do not explicitly count modes, but
the unsubtracted vacuum energy on its right can diverge with the
regulator. Consequently (9) alone does not provide a continuum quantum
bound or a finite renormalized capture error. Those require the ordered
work, dressing and stress calculation that remains open.

Exact cancellation also assumes the unchanged integer strip spectrum.
For a displaced frequency `2n+d`, (5) gives the elementary bound

```
|F(2n+d)| <= h_f |d|/(2n-|d|), |d|<2n.                 (10)
```

The smooth pulse supplies additional high-frequency decay. Physical
mode shifts, varying clocks and geometry, and their aggregate weighted
error have not been bounded. Fourier cancellation is an established
control method: [Guery-Odelin and Muga, sections II/III and appendix B](https://arxiv.org/html/1410.4957v2)
derive its role in oscillator transport. Equations (2)--(10) apply that
method to this complete restoring spectrum; no novelty is accepted.

## Capture compatibility and verification

The [JT switching proof](switching-radiation-capture.md) bounds its
global duration by `48/127<pi`, before the first returned signal.
The same source convention gives the return period used here. Its
pre-return incident covariance cannot be carried unchanged through
the longer switch (5). Nor is the reference `h` already the physical
four-dimensional coupling: vertex normalization remains to be matched.
The longer protocol therefore needs a returned-state stress calculation
and a new capture/area bound before it can help establish stable opening.

[Checker](../checks/electric_return_cancellation.py),
[receipt](../evidence/electric-return-cancellation.json), and
[review record](../evidence/electric-return-source-review.json).
Independent piecewise Fourier quadrature checks 128 even modes at each
of three smoothing widths, off-lattice factorization and the detuning
bound. Twelve monotone fast pulses check (4). Coupled oscillator ODEs
propagate the switch with both noise sectors, followed by an independent
matrix exponential of the autonomous Hamiltonian; endpoint returns,
energy conservation, uncertainty and the phase bound are checked.
A tighter ODE tolerance checks numerical refinement. The ordinary
Fourier and positivity proofs, not the finite sample, cover all modes
and all future times within their stated scopes.

Remaining to-do list: control the rapid candidate's driven state or
derive capture and returned stress for a longer shaped switch; match
physical clocks, vertices and initial dressing; include nonlinear
angular/massive interactions and evolving geometry; establish preparation,
sufficient stability and finite-energy transport.
