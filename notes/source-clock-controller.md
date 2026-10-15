# Weyl contact and a non-prescribed controller on the varied clock

S074 course, E003, 2026-10-14. This note carries the finite Weyl contact
`c_theta` and a controller whose angle is not a prescribed c-number onto
the same rule-D clock as the benchmark window. It does not close
S074.2–S074.6, and it does not mark a physical opening.

The comparison, the hardest step, the observable, and the test below
are fixed before the fit. The benchmark numbers that follow them are
the output of that test.

## Comparison

The retained rule-D clock, on the returned-switch benchmark, is

```
z_tt + (1 + lambda (Theta_t)^2) z = g/2,
Theta_t = r0(t) = 4 eta [S(t/tau) - S((t-pi)/tau)],
lambda = pi zeta_s ell^2 / N,
```

with `g` the retained returned-switch drive. Two further finite inputs
are now admitted on that same equation.

The curved-collar finite coupling `c_theta` multiplies

```
I_theta = integral sqrt(gamma) Lhat^{ab} u_a u_b,   u_a = D_a theta.
```

For a homogeneous angle and a `y`-independent `Lhat^{uu}`, the
Lorentzian sphere integral of that density is the u-Lagrangian
`A_theta y^2 (dTheta/dt)^2` with
`A_theta = 4 pi ell^2 c_theta Lhat^{uu}`. Its clock coupling is
`lambda_theta = A_theta/(2N) = (2 pi ell^2 / N) c_theta Lhat^{uu}`.
The prescribed-angle clock sees only the sum
`Lambda = lambda + lambda_theta`. Equivalently

```
Lambda = (pi ell^2 / N) (zeta_s + 2 c_theta Lhat^{uu}).
```

The `1/2` in the retained source action `I_s` is why the curvature
contact carries the factor 2 relative to `zeta_s`.

The controller replaces the prescribed rate by the affine response

```
rho(t) = r0(t) + kappa z,
```

with `kappa` a finite coupling and `z = y/eta`. The response depends
on the clock and not on `z_t`. The matter drive `g` stays the retained
force built from the prescribed pulse.

## Hardest step

`I_theta` needs `Lhat^{uu}` on the varied clock boundary. The heat
theorem treats `L_ab` as an arbitrary symmetric tensor, so the boundary
metric `gamma = -du^2 + ell^2 dOmega^2` does not determine it. The
relative RN-AdS weight tends to zero at the conformal boundary and is
the anomaly transgression, not this finite counterterm. An umbilic
boundary would make `I_theta` vanish; that has not been established
for this clock. Until `Lhat^{uu}` is supplied, `c_theta` and `zeta_s`
are not separate coefficients of this equation.

A `z`-dependent rate is not homogeneous of degree 2, so the u-energy
is no longer equal to the Lagrangian. The charge has to be varied
before any window is computed. The radiation law
`F_t/N = [theta_dot^2 - theta_dot]/(4 pi)` was derived for a
prescribed angle before the first return. Evaluating that law on a
state-dependent angle would replace the retained drive. This
calculation does not do that. The drive-channel coupling is an
unresolved input. It is not assigned the value zero.

No replacement area target is written before the test. The area budget
stays unevaluated. `Phi` still depends on `y` and `E`.

## Observable

The sample is the retained one: `eta = 1/2000` and the six corners
`z(0) in {1/4, 7/25}`, `z_t(0) in {-1/2000, 0, 1/2000}`. Deviations
are taken against the retained trajectory at `Lambda = 0`, `kappa = 0`.
The frozen budgets are

```
sup |delta z| on the switch < 1/200,
sup |delta z| after the switch < 1/50,
|delta e_f| < 1/50,
```

with `e_f = z_t(T_star)^2 + z(T_star)^2 - g_f z(T_star)` and
`g_f = 1 - 4 eta`, evaluated at `T_star`. A trajectory that reaches
`z <= 0` leaves the retained clock domain and is outside the
admissible set. The prescribed-angle anchors
`-4757.61079788208 <= lambda <= 4875.4563331604` are the comparison
interval already recorded for `kappa = 0`. They are not refit here.

## Test

1. Direct variation. For `L = A y^2 rho(y)^2` and
   `rho = r0 + kappa y/eta`, the global-time charge reduces to the
   clock equation

```
z_tt + (1 + Lambda r0^2) z + 4 Lambda kappa r0 z^2
    + 3 Lambda kappa^2 z^3 = g/2.
```

   At `kappa = 0` this is rule D. At `Lambda = 0` every `kappa` term
   drops out. A central difference of `L` at the predeclared point
   `A = 1`, `eta = 1/2000`, `y = 0.28 eta`, `r0 = 0.001`, `kappa = 0.2`
   reproduces that charge.
2. Additivity. A split of the positive anchor into a `zeta_s` piece
   and a `c_theta Lhat^{uu}` piece with the same sum produces the
   retained rule-D trajectory at that anchor.
3. Independent integration. DOP853 and RK45 agree on the controller
   trajectory at the probed edge.
4. Controller window. At each fraction `1/4`, `1/2`, `3/4`, and `1`
   of each recorded anchor, both signs, the checker reports the
   largest `kappa` of each sign for which all three budgets hold on
   all six corners. The bisection returns an admissible endpoint.
   Interior fractions `1/4`, `1/2`, `3/4`, and `1` of that endpoint
   are probed on the most restrictive corner of the half positive
   anchor, on both signs. The other samples are bracketed, and their
   interiors are not probed one by one.
5. Overshoot probe. At `Lambda` equal to `1.1` times the positive
   anchor, on the corner `z(0) = 7/25`, `z_t(0) = 1/2000`, the
   magnitudes `10^{-4}`, `10^{-3}`, `10^{-2}`, `10^{-1}`, and `1`
   are tested on both signs. This probe does not map a joint boundary.

`kappa = 0` is the prescribed-angle reduction. It is a check, not the
controller result. The contact-channel range is the family of `kappa`
endpoints in step 4.

## Derivation

The u-Lagrangian `L = A y^2 rho^2` with `rho = r0 + (kappa/eta) y`
has momentum `y dL/dy` and u-energy `H_u = y dL/dy - L`. Direct
variation gives

```
H_u = A y^2 rho^2 + 2 A (kappa/eta) y^3 rho,
E_t = H_u/y = A y rho^2 + 2 A (kappa/eta) y^2 rho.
```

Insert `E_t` beside the retained matter charge, divide by `eta`, and
set `Lambda = A/(2N)`. The cubic equation in the test is that
reduction. The retained source and the curvature contact add in `A`
before this variation whenever `Lhat^{uu}` is independent of `y` and
`rho` is still the prescribed rate. The heat polynomial `a4` is built
from `L_ab` and `u_a` and does not contain `c_theta`.

## Benchmark result

The acceptance receipt is `evidence/source-clock-controller.json`.
Both acceptance runs agree on the identities, the joint constraint,
the controller endpoints, and the overshoot probe.

`c_theta` has no separate interval on this benchmark. The identified
input that separates it from `zeta_s` is `Lhat^{uu}` of the varied
clock boundary. With a prescribed angle, the clock constrains only
the sum `zeta_s + 2 c_theta Lhat^{uu}` through `Lambda`. The recorded
anchor interval remains that constraint at `kappa = 0`. Splitting the
positive anchor into equal `lambda` and `lambda_theta` pieces
reproduces the retained rule-D trajectory with gap `0`. The central
difference of the varied charge matches the closed form with relative
gap `6.702802905747133e-09`.

The controller endpoints below are the admissible bisection values of
`kappa` that hold on all six corners at once. These are the admissible
`lo` values next to a rejected neighbor, not a proved supremum past
that neighbor. The binding budget is the one closest to its threshold
on the corner that sets the shared endpoint. That corner is
`z(0) = 7/25`. Its initial velocity is `-1/2000` when the post-switch
budget binds and `+1/2000` when the switch budget binds.

```
Lambda                  kappa interval                         binding (+ / -)
1218.8640832901         [-0.0021144577026367193, 0.0021056297302246094]   post / post
2437.7281665802         [-0.0015083452224731444, 0.001492590808868408]    post / post
3656.5922498703003      [-0.0012409776687622074, 0.0007495726108551027]   switch / post
4875.4563331604         [-0.0010815817832946781, 1.490116119384766e-10]   switch / post
-1189.40269947052       [-0.002084294509887696, 0.002086107635498047]     post / post
-2378.80539894104       [-0.0014881957054138186, 0.0014863876342773438]   post / post
-3568.20809841156       [-0.0012264337539672852, 0.0007515277385711672]   switch / post
-4757.61079788208       [-0.0010714412689208984, 1.1324882507324219e-10]  switch / post
```

At the positive anchor the switch budget is already met by the contact
term, so the positive `kappa` endpoint sits at `1.490116119384766e-10`
with switch deviation `0.004999999999556939`. The opposite sign still
has room, `-0.0010815817832946781`, and there the post-switch budget
binds at `0.019999998733850444`. Where the post-switch budget binds,
the endpoint scales approximately as `1/sqrt(|Lambda|)`, which is the
cubic term `3 Lambda kappa^2 z^3`.

On the half positive anchor, `Lambda = 2437.7281665802`, the corner
`z(0) = 7/25`, `z_t(0) = -1/2000` was probed at one quarter, one half,
three quarters, and the endpoint. All eight samples are admissible.
The post-switch deviation increases through those fractions on both
signs. The switch deviation increases on the positive sign and
decreases on the negative sign, so that one observable is not the
certificate. DOP853 and RK45 agree on the positive endpoint with gap
`1.1297329738368944e-10`. The other seven contact strengths were
bracketed, and their interiors were not probed one by one.

The overshoot `Lambda = 5363.00196647644` is `1.1` times the positive
anchor, on the single corner `z(0) = 7/25`, `z_t(0) = 1/2000`.
`kappa = 0` is not admissible there. Of the listed magnitudes, only
`kappa = -0.001` is admissible, with switch deviation
`0.003496518304417648`, post-switch deviation
`0.018637263343309618`, and `|delta e_f| = 0.002493530795375204`.
That sample is not a six-corner range, and it does not map the joint
boundary. The prescribed-angle anchors are not the full boundary once
`kappa` is free.

Setting `Lambda = 0` removes every controller correction in this
contact channel. That identity is not an admissible range for `kappa`.
The radiation-law response of `g` to a non-prescribed angle remains
unresolved. No retained input selects a point inside the anchor
interval or a preferred `kappa`.

The area budget is not evaluated. S074.2–S074.6 stay open. S075
through S079 stay blocked. This is not an impossibility theorem and
not a physical opening.

[Checker](../checks/source_clock_controller.py),
[receipt](../evidence/source-clock-controller.json), and
[source/run review](../evidence/source-clock-controller-review.json).
No Lean, paid compute, publication, or outreach was used.

Remaining to-do list: identify `Lhat^{uu}` on the varied clock
boundary, and evaluate the retained radiation law on a
non-prescribed angle if that drive channel is kept. The
contact-channel `kappa` slices are recorded. The joint boundary
outside the prescribed-angle anchors is not mapped. Unequal collars,
the nonlocal interacting determinant, and `cJ` stay unmatched.
S074.2–S074.6 stay open. S075 through S079 stay blocked. The area
budget stays unevaluated.
