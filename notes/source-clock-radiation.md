# Clock weight of Lhat and the radiation law on a realized angle

S074 course, E003, 2026-10-14. This note identifies `Lhat^{uu}` on the
varied clock boundary retained by the JT reduction, and it evaluates
the returned-switch radiation law on an angle that is not a prescribed
c-number. It does not close S074.2–S074.6, and it does not mark a
physical opening.

The comparison, the hardest step, the observable, and the test below
are fixed before the fit.

## Comparison

The retained rule-D clock is

```
z_tt + (1 + Lambda r0(t)^2) z = g/2,
r0(t) = 4 eta [S(t/tau) - S((t-pi)/tau)],
```

with `g` taken from the prescribed returned switch. Two inputs left
open by the contact-channel calculation are now evaluated on that
same benchmark.

The varied clock boundary of the retained near-horizon reduction is
the cutoff curve in the product throat

```
ds^2 = r_e^2 (-dt^2 + dz^2)/z^2 + r_e^2 dOmega^2,
```

at fixed induced `gamma_tt = -r_e^2/epsilon^2`. Its outward extrinsic
curvature is the retained expansion

```
K = r_e^{-1} [1 - epsilon^2 Sch(F,t) + O(epsilon^4)].
```

The sphere radius is constant, so the angular extrinsic curvature
vanishes. With the inward heat-kernel convention this is not an
umbilic surface: `Lhat_tt = (2/3) K gamma_tt`. The sphere integral
of `Lhat^{ab} u_a u_b` is

```
-(8 pi/3) epsilon (1 - epsilon^2 Sch) (d theta/dt_cutoff)^2 + O(epsilon^5).
```

In the conformal-boundary limit that defines the JT clock, this
weight vanishes. A finite `c_theta` multiplies that weight. It does
not enter the retained clock equation. This is the geometric factor
going to zero with `epsilon`. It is not an assignment of `c_theta`.

The realized angle rate remains `rho = r0(t) + kappa z`. The
returned-switch law is evaluated on that rate. With `D = 0` before
`t = 0` and `D(t) = rho(t) + D(t-pi)`,

```
theta(t) = integral_0^t rho,   Q(t) = integral_(t-pi)^t D^2,
E/N = [-theta + Q]/(4 pi),     g = (theta - Q)/(4 pi eta).
```

At `kappa = 0` this `g` is the retained force. The contact charge
keeps the varied-clock reduction already checked,

```
z_tt + (1 + Lambda r0^2) z + 4 Lambda kappa r0 z^2
    + 3 Lambda kappa^2 z^3 = g/2,
```

now with this realized `g`.

## Hardest step

`Lhat_ab` is fixed by the cutoff embedding, not left as an arbitrary
tensor. The product throat is not umbilic, and the lapse still makes
the integrated weight `O(epsilon)`. The JT benchmark is the limit
`epsilon -> 0` in which that weight drops out. A sphere radius that
varies along the normal is outside this product metric.

The radiation law makes `g` a delayed functional of `z`. The method
of steps has to reproduce `make_forces` at `kappa = 0` before any
`kappa` window is accepted. The contact terms still vanish at
`Lambda = 0`; the radiation terms do not. `Lambda = 0` is therefore
a real sample of this drive, not a hidden controller.

No replacement area target is written before the test. The area
budget stays unevaluated.

## Observable

The sample is `eta = 1/2000` and the six corners
`z(0) in {1/4, 7/25}`, `z_t(0) in {-1/2000, 0, 1/2000}`. Deviations
are taken against the retained trajectory at `Lambda = 0`,
`kappa = 0`. The frozen budgets are

```
sup |delta z| on the switch < 1/200,
sup |delta z| after the switch < 1/50,
|delta e_f| < 1/50,
```

with the retained diagnostic
`e_f = z_t(T_star)^2 + z(T_star)^2 - g_f z(T_star)`,
`g_f = 1 - 4 eta`. A trajectory that reaches `z <= 0` leaves the
retained domain and is outside the admissible set. The contact
strengths tested with the realized drive are `Lambda = 0` and the
fractions `1/2` and `1` of each prescribed-angle anchor
`-4757.61079788208` and `4875.4563331604`.

## Test

1. Direct geometry. For `L_tt = K gamma_tt`, vanishing angular
   extrinsic curvature, and the retained `K`, the sphere-integrated
   density equals the `O(epsilon)` expression above and its limit
   at `epsilon = 0` is zero.
2. Retained drive. At `kappa = 0` the delay law reproduces
   `make_forces` on `[0, T_end]`, and the clock trajectory at the
   positive anchor reproduces the retained rule-D integrator.
3. Independent integration. DOP853 and RK45 agree on the realized
   drive at the `Lambda = 0` edge.
4. Controller window. At each predeclared `Lambda`, both signs, the
   checker reports the largest `kappa` of each sign for which all
   three budgets hold on all six corners. Interior fractions `1/4`,
   `1/2`, `3/4`, and `1` of that endpoint are probed on the most
   restrictive corner of the `Lambda = 0` sample, on both signs.
   The other samples are bracketed, and their interiors are not
   probed one by one.
5. Overshoot probe. At `Lambda` equal to `1.1` times the positive
   anchor, on the corner `z(0) = 7/25`, `z_t(0) = 1/2000`, the
   magnitudes `10^{-4}`, `10^{-3}`, `10^{-2}`, `10^{-1}`, and `1`
   are tested on both signs. This probe does not map a joint boundary.

## Derivation

The inward tensor with `L_tt = K gamma_tt` and `L_angular = 0` has
trace `K` and `Lhat_tt = (2/3) K gamma_tt`. Raising the index and
integrating the round sphere produces one factor of `epsilon` from
`sqrt(-gamma) gamma^{tt}`. Substituting the retained `K` leaves the
displayed density. Its Taylor series in `epsilon` begins at order
`epsilon`. The Schwarzian correction begins at order `epsilon^3`.

The delay identity `D(t) - D(t-pi) = rho(t)` is the returned-switch
relation between the outgoing profile and the junction rate. For the
prescribed rate it rebuilds `D(t) = 4 eta S(t/tau)` and the force
`g = a - r` used by the benchmark. The contact reduction is the same
varied charge as the preceding controller calculation, evaluated on
the realized `g`.

## Benchmark result

The acceptance receipt is `evidence/source-clock-radiation.json`.
Both acceptance runs agree on the geometry, the retained drive, the
anchor trajectory, the `kappa` endpoints, and the budget comparison.

The checker confirms the product-throat density. The integrated
weight is `-(8 pi/3) epsilon (1 - epsilon^2 Sch) rate^2`. Its limit
at `epsilon = 0` is `0`. The surface is not umbilic. The Schwarzian
correction begins at order `epsilon^3`. `c_theta` has no separate
clock interval on this benchmark. The identified input is this
`O(epsilon)` weight. In the JT limit it drops out of the clock
equation. That limit is the geometry of the retained reduction.

At `kappa = 0` the delay law reproduces `make_forces`. The maximum
absolute drive error on the 401-point grid is
`5.875251396503245e-10`. The drive at `T_star` is
`0.9980000004377896`, against `1 - 4 eta = 0.998`. The clock
trajectory at `Lambda = 4875.4563331604` and `kappa = 0` reproduces
the retained rule-D integrator with gap
`1.0325795773979962e-10`.

The endpoints below are the admissible bisection values of `kappa`
that hold on all six corners at once. They are the admissible `lo`
values next to a rejected neighbor. The binding budget is the one
closest to its threshold on the corner that sets the shared endpoint.

```
Lambda                  kappa interval                                      binding (+ / -)
0.0                     [-5.16183614730835e-05, 5.072014927864075e-05]      post / post
2437.7281665802         [-4.602195620536803e-05, 5.619207024574278e-05]     post / post
4875.4563331604         [-0.0, 6.0794997215271e-05]                         post / switch
-2378.80539894104       [-5.705349445343019e-05, 4.523511528968811e-05]     post / post
-4757.61079788208       [-6.149901151657105e-05, 1.788139343261719e-11]     switch / post
```

At `Lambda = 0` both endpoints come from the corner `z(0) = 7/25`,
`z_t(0) = -1/2000`. The positive post-switch deviation is
`0.019999999835735978`, and the negative one is
`0.019999999367255057`. The same corner,
on both signs, was probed at one quarter, one half, three quarters,
and the endpoint. All eight samples are admissible. Switch,
post-switch, and energy deviations increase through those fractions
on both signs. DOP853 and RK45 agree on the positive endpoint with
gap `1.9633215053715958e-09`. The other four `Lambda` samples were
bracketed, and their interiors were not probed one by one.

At `Lambda = 4875.4563331604` the negative endpoint is `-0.0`. The
`kappa = 0` edge on the corner `z(0) = 7/25`, `z_t(0) = 1/2000` has
switch deviation `0.004999999995464445`. At
`Lambda = -4757.61079788208` the positive endpoint is
`1.788139343261719e-11`. That edge, on the same corner and velocity,
has switch deviation `0.0049999999383626564`. Both edges remain
inside the switch budget. The opposite signs at those anchors are
the post-switch endpoints in the table.

The overshoot `Lambda = 5363.00196647644` is `1.1` times the positive
anchor, on the single corner `z(0) = 7/25`, `z_t(0) = 1/2000`.
`kappa = 0` lies outside the budgets. Every probed magnitude
`10^{-4}`, `10^{-3}`, `10^{-2}`, `10^{-1}`, and `1`, on both signs,
also lies outside the budgets. That sample does not map a joint
boundary.

`Lambda = 0` is a real radiation sample. The contact terms vanish
there, and the realized drive still bounds `kappa`. The area budget
is not evaluated. S074.2–S074.6 stay open. S075 through S079 stay
blocked. This is not an impossibility theorem and not a physical
opening.

[Checker](../checks/source_clock_radiation.py),
[receipt](../evidence/source-clock-radiation.json), and
[source/run review](../evidence/source-clock-radiation-review.json).
No Lean, paid compute, publication, or outreach was used.

Remaining to-do list: on the realized angle,
`sigma = 2 c_theta Lhat^{uu}` enters through
`Lambda = (pi ell^2 / N)(zeta_s + sigma)`. At
`kappa = 2.5360074639320374e-05` the six-corner interval is
`-3957.533203125 <= Lambda <= 5497.42431640625`, switch on both signs.
At `kappa = -2.580918073654175e-05` the six-corner interval is
`-5509.62939453125 <= Lambda <= 4146.0244140625`, switch on both signs.
Both edges are the corner `z(0) = 7/25`, `z_t(0) = 1/2000`. `zeta_s`
stays unmatched inside that sum, so `sigma` has no separate number. The
product-throat weight remains the recorded `O(epsilon)` geometric factor.
The half-anchor `kappa` slices, on all 6 corners, are
`2437.7281665802: [-4.602195620536803e-05, 5.619207024574278e-05]`
and `-2378.80539894104: [-5.705349445343019e-05, 4.523511528968811e-05]`.
The one-corner full-edge probe at the `Lambda = 0` endpoints keeps a
same-sign anchor admissible and is not an empty interval. Unequal
collars, the nonlocal interacting determinant, and `cJ` stay unmatched.
S074.2–S074.6 stay open. S075 through S079 stay blocked. The area budget
stays unevaluated. None further for the S080 obstruction audit, the S081
benchmark interval, the S082 contact-channel slices, or the S083
radiation-channel slices.
