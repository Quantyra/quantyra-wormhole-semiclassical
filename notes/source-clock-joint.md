# Finite Weyl product on the realized angle

S074 course, E003, 2026-10-14. This note carries the finite product
`sigma = 2 c_theta Lhat^{uu}` and the realized angle
`rho = r0 + kappa z` into the same rule-D clock. It does not close
S074.2–S074.6, and it does not mark a physical opening.

The comparison, the hardest step, the observable, and the test below
are fixed before the fit. The benchmark numbers that follow them are
the output of that test.

## Comparison

The retained realized-drive clock is

```
z_tt + (1 + Lambda r0^2) z + 4 Lambda kappa r0 z^2
    + 3 Lambda kappa^2 z^3 = g/2,
rho = r0 + kappa z,
D(t) = rho(t) + D(t-pi),   D = 0 for t <= 0,
theta = integral_0^t rho,  Q = integral_(t-pi)^t D^2,
g = (theta - Q)/(4 pi eta).
```

`r0(t) = 4 eta [S(t/tau) - S((t-pi)/tau)]`. The diagnostic
`g_f = 1 - 4 eta` stays the frozen returned-switch value.

Two earlier calculations leave this joint entry open. On a prescribed
drive, `c_theta` enters the clock only through

```
sigma = 2 c_theta Lhat^{uu},
lambda_theta = (pi ell^2 / N) sigma,
lambda_zeta = pi zeta_s ell^2 / N,
Lambda = lambda_zeta + lambda_theta
       = (pi ell^2 / N) (zeta_s + sigma).
```

That is the u-boundary charge on `gamma = -du^2 + ell^2 dOmega^2`
with a `y`-independent `Lhat^{uu}`. Setting `Lambda = 0` removes every
contact correction of that channel. That reduction is a check. It is
not an admissible range for `sigma`.

On the product-throat cutoff the integrated weight of `Lhat^{uu}` is
`O(epsilon)` and drops out of the JT limit. A finite `c_theta`
multiplies that geometric factor. The factor is the cutoff. It is not
an assignment of `c_theta`, and it is not an assignment of `sigma`.

This test keeps both inputs in the one equation above. `Lambda` carries
`zeta_s + sigma`. The delay law carries `kappa`. Neither coefficient
is set to zero as the result.

## Hardest step

The cutoff surface and the u-boundary charge are different. The product
throat at fixed `gamma_tt = -r_e^2/epsilon^2` supplies the recorded
`O(epsilon)` weight. The clock charge used here still contains `sigma`
when `Lhat^{uu}` is the finite tensor in `I_theta`. Substituting the
cutoff limit for that tensor would replace `sigma` by zero. This test
does not do that.

At nonzero `kappa` the linear term `Lambda r0^2 z` still bounds
`Lambda`. The quadratic and cubic contact terms, and the delayed `g`,
move that bound. The shared `kappa` endpoints at `Lambda = 0` already
sit on a budget. The predeclared nontrivial `Lambda` slice is half of
each shared endpoint. A probe at the full endpoint records
admissibility at the half anchors and the full anchors. A failed sample
in that probe is not an empty-interval claim.

`lambda_zeta` and `lambda_theta` add before the variation. An interval
for `Lambda` at fixed nonzero `kappa` is an interval for
`zeta_s + sigma`. No retained identity selects a value of `zeta_s`
inside the recorded anchor interval, so the interval does not split
into a separate number for `sigma`. The equivalent product statement is
`sigma = Lambda N / (pi ell^2) - zeta_s`, with `zeta_s` left as that
unmatched input.

No replacement area target is written before the test. The area budget
stays unevaluated.

## Observable

The sample is `eta = 1/2000` and the corners returned by the shipped
radiation checker. That list is `z(0) in {1/4, 7/25}` and
`z_t(0) in {-1/2000, 0, 1/2000}`. The receipt records
`corner_count` as the length of that list. Deviations are taken
against the retained trajectory at `Lambda = 0`, `kappa = 0`. The
frozen budgets are

```
sup |delta z| on the switch < 1/200,
sup |delta z| after the switch < 1/50,
|delta e_f| < 1/50,
```

with `e_f = z_t(T_star)^2 + z(T_star)^2 - g_f z(T_star)` and
`g_f = 1 - 4 eta`. A trajectory that reaches `z <= 0` leaves the
retained domain and is outside the admissible set. The area budget is
not part of the comparison.

The declared `kappa` values for the `Lambda` bisection are half the
shared `Lambda = 0` endpoints of each sign, recomputed by the shipped
`bracket_kappa`. The declared `Lambda` values for the `kappa`
bisection are half of each recorded anchor,
`0.5 * 4875.4563331604` and `0.5 * -4757.61079788208`.

## Test

1. Direct sum. Sympy checks `sigma = 2 c_theta Lhat^{uu}`,
   `lambda_theta = (pi ell^2 / N) sigma`, the equality of that sum
   with `(pi ell^2 / N)(zeta_s + sigma)`, the controller form
   `lambda_theta = (2 pi ell^2 / N) c_theta Lhat^{uu}`, and the charge
   factor `(r0 + kappa z)(r0 + 3 kappa z) = r0^2 + 4 kappa r0 z + 3 kappa^2 z^2`.
2. Shipped regressions. The checker calls `drive_regression` and
   `anchor_regression` from the radiation checker. Those functions
   assert the retained drive and the positive-anchor trajectory.
3. Prescribed reduction. At `kappa = 0`, both recorded anchors are
   admissible on every shipped corner. At `1.1` times the positive
   anchor, on the corner `z(0) = 7/25`, `z_t(0) = 1/2000`, the sample
   `kappa = 0` is outside the three budgets. The `kappa = 0` lambda
   window is not re-bisected.
4. Radiation regression. At `Lambda = 0`, the shipped `bracket_kappa`
   on every shipped corner, both signs, reproduces the recorded
   radiation endpoints
   `[-5.16183614730835e-05, 5.072014927864075e-05]`.
5. Finite product at a nonzero angle. At `kappa` equal to half of each
   shared `Lambda = 0` endpoint, both signs, the checker bisects
   `Lambda` on every shipped corner, both signs. The reported endpoint
   is the intersection. The search starts at `1`, doubles until a
   budget fails, stops past `10^6`, and then takes 24 bisections. The
   returned value is the admissible `lo`.
6. Finite angle at a nonzero contact. At `Lambda` equal to half of
   each recorded anchor, both signs, the shipped `bracket_kappa`
   reports the `kappa` interval on every shipped corner. Those samples
   are the recorded radiation endpoints
   `2437.7281665802: [-4.602195620536803e-05, 5.619207024574278e-05]`
   and
   `-2378.80539894104: [-5.705349445343019e-05, 4.523511528968811e-05]`.
   The joint receipt contains them so both axes are one comparison.
7. Full-edge probe. At each shared `Lambda = 0` endpoint, on the
   corner that set it, the checker evaluates `Lambda` at half and full
   of each recorded anchor. It records admissibility. It does not turn
   a failed sample into an empty interval.
8. Independent integration. DOP853 and RK45 agree at the positive
   `Lambda` endpoint of the positive half-`kappa` slice, on the corner
   that sets that shared endpoint.
9. The receipt records `corner_count` from the shipped corner list.
   The area budget stays unevaluated. `opening_achieved` stays false.
   `s074_gate_closed` stays false. The result is not an impossibility
   theorem. `c_theta` and `sigma` stay explicit.

`kappa = 0` and `Lambda = 0` are reduction checks. The joint range is
the `Lambda` interval at the half-`kappa` slice together with the
`kappa` interval at the half anchors.

## Derivation

The u-boundary reduction already checked for a prescribed angle gives
`lambda_theta = (2 pi ell^2 / N) c_theta Lhat^{uu}`. The definition
`sigma = 2 c_theta Lhat^{uu}` rewrites that coupling as
`(pi ell^2 / N) sigma`. Adding `lambda_zeta = pi zeta_s ell^2 / N`
produces the displayed `Lambda`. The same reduction expands the
varied charge by the factor
`(r0 + kappa z)(r0 + 3 kappa z)`, which is the quadratic and cubic
polynomial in the clock equation. The force `g` is the returned-switch
delay evaluated on `rho = r0 + kappa z`, the law already checked at
`kappa = 0` against `make_forces`.

## Benchmark result

The acceptance receipt is `evidence/source-clock-joint.json`. Both
acceptance runs agree on the charge identities, the shipped
regressions, the `Lambda` intervals, the `kappa` intervals, and the
budget comparison. The shipped corner list has length 6.

The charge identities hold. `sigma = 2 c_theta Lhat^{uu}` enters as
`lambda_theta = (pi ell^2 / N) sigma`, and
`Lambda = lambda_zeta + lambda_theta` equals
`(pi ell^2 / N)(zeta_s + sigma)`. The varied-charge factor expands to
`r0^2 + 4 kappa r0 z + 3 kappa^2 z^2`. At `kappa = 0` the delay law
reproduces `make_forces` with maximum absolute error
`5.875251396503245e-10`, and the final drive is
`0.9980000004377896` against `0.998`. The clock trajectory at
`Lambda = 4875.4563331604` and `kappa = 0` reproduces the retained
rule-D integrator with gap `1.0325795773979962e-10`.

At `kappa = 0` both recorded anchors are admissible on all 6 corners.
At `Lambda = 5363.00196647644`, on the corner `z(0) = 7/25`,
`z_t(0) = 1/2000`, the sample `kappa = 0` is outside the budgets. The
switch deviation is `0.005493315463730761`. The switch budget is
`1/200`.

At `Lambda = 0` the shipped bracket reproduces the recorded radiation
endpoints `[-5.16183614730835e-05, 5.072014927864075e-05]`. Both edges
are the corner `z(0) = 7/25`, `z_t(0) = -1/2000`, and both bind on the
post-switch budget. The positive post-switch deviation is
`0.019999999835735978`. The negative one is `0.019999999367255057`.

The half-endpoint slices below are the admissible bisection values of
`Lambda` that hold on all 6 corners at once. They are the admissible
`lo` values next to a rejected neighbor. The binding budget is the one
closest to its threshold on the corner that sets the shared endpoint.

```
kappa                      Lambda interval                         binding (+ / -)
2.5360074639320374e-05     [-3957.533203125, 5497.42431640625]     switch / switch
-2.580918073654175e-05     [-5509.62939453125, 4146.0244140625]    switch / switch
```

Every edge in that table is the corner `z(0) = 7/25`,
`z_t(0) = 1/2000`. The switch deviations are
`0.004999999757793729` and `0.004999999802967798` on the positive
`kappa` slice, and `0.004999999656788608` and
`0.004999999523526456` on the negative `kappa` slice. DOP853 and RK45
agree on the positive endpoint of the positive slice,
`Lambda = 5497.42431640625`, with gap
`3.897244749140327e-10`.

These intervals constrain `zeta_s + sigma`. No retained identity
selects `zeta_s` inside them, so `sigma` has no separate number on
this benchmark. The product-throat weight remains the recorded
`O(epsilon)` geometric factor. It is not used as a value of `sigma`.

The shipped bracket at the half anchors reproduces the recorded
radiation endpoints. On all 6 corners,

```
Lambda                  kappa interval                                      binding (+ / -)
2437.7281665802         [-4.602195620536803e-05, 5.619207024574278e-05]     post / post
-2378.80539894104       [-5.705349445343019e-05, 4.523511528968811e-05]     post / post
```

The full-edge probe uses the `Lambda = 0` endpoints on the single
corner `z(0) = 7/25`, `z_t(0) = -1/2000`. At
`kappa = 5.072014927864075e-05`, the samples
`Lambda = 2437.7281665802` and `Lambda = 4875.4563331604` are
admissible, and the samples `Lambda = -2378.80539894104` and
`Lambda = -4757.61079788208` are outside the budgets. At
`kappa = -5.16183614730835e-05` the opposite pair is admissible and
the positive anchors are outside. That probe is not a six-corner range
and not an empty-interval claim.

The area budget is not evaluated. S074.2–S074.6 stay open. S075
through S079 stay blocked. This is not an impossibility theorem and
not a physical opening.

[Checker](../checks/source_clock_joint.py),
[receipt](../evidence/source-clock-joint.json), and
[source/run review](../evidence/source-clock-joint-review.json).
The radiation checker and the controller checker are reused and are
not edited. No Lean, paid compute, publication, or outreach was used.

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
