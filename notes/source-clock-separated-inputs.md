# Separated inputs after the retained clock sum

S086, E003, 2026-10-14. S085 recorded that the retained rule-D clock
does not separate `zeta_s` inside `Lambda`, does not contain unequal
collars or the nonlocal interacting determinant, and does not measure
`cJ`. This note takes the four inputs named after that record. It does
not close S074.2–S074.6, and it does not mark a physical opening.

The assumption, the error comparison, the comparison, the hardest
step, the observable, and the test for each input are fixed below,
before any checker run. Measured values, if the test produces them,
are written only in the benchmark result and only from the receipt.

The recorded S084 intervals stay the latest ranges of `Lambda`. This
note does not replace them. The S085 march is not repeated. No
coefficient is set to zero. No endpoint contact energy is cancelled
as a result. No area target is written. The area budget stays
unevaluated. No successor model is opened.

## Shared benchmark

The retained clock is the equation checked in
[the joint note](source-clock-joint.md):

```
z_tt + (1 + Lambda r0^2) z + 4 Lambda kappa r0 z^2
    + 3 Lambda kappa^2 z^3 = g/2,
Lambda = (pi ell^2 / N) (zeta_s + sigma),
sigma = 2 c_theta Lhat^{uu}.
```

`r0(t) = 4 eta [S(t/tau) - S((t-pi)/tau)]`. Deviations on that clock
stay measured against `Lambda = 0`, `kappa = 0`, inside the frozen
budgets of the joint note. This test does not march that equation and
does not re-bracket `Lambda`.

## `zeta_s` and `sigma`

Assumption. The coordinate component `Lhat^{uu}` is held fixed while
the lapse `N` and the sphere radius `R` vary. The boundary metric is
`gamma = -N^2 dt^2 + R^2 dOmega^2`, and the angle is homogeneous. On
the round metric `N = 1`, `R = ell`, this is the same `Lhat^{uu}` that
enters `sigma`.

Error comparison. The collar tensor in the curved-fermion note is
`L_ab = -<nabla_a n, e_b>`. A lapse variation moves that tensor. The
total lapse derivative of `c_theta I_theta` therefore contains the
extra piece `(partial lag / partial Lhat^{uu}) (partial Lhat^{uu} / partial N)`.
This test does not evaluate the second factor. The product-throat
`O(epsilon)` weight is not used as `sigma`.

Comparison. The clock charge depends on `zeta_s` and `sigma` only
through the sum. The one alternative compares that charge with the
sphere-integrated actions. The source integrand is the published
contact `2 pi zeta_s R^2 (d theta/dt)^2 / N`. At fixed `Lhat^{uu}` the
curvature contact is `4 pi c_theta N R^2 Lhat^{uu} (d theta/dt)^2`.
The clock-sum integrand gives both coefficients the source's `1/N`
dependence. The S085 opposite shift is not marched again.

Hardest step. At `N = 1` the two actions add to the clock-sum
integrand, so the value on the round metric does not separate them.
The radius derivatives also agree there. The lapse derivatives agree
only when `sigma` vanishes, and this test does not make it vanish.
Equating the fixed-`Lhat` lapse response to the collar stress would
drop the extra piece above.

Observable. At `N = 1` and `R = ell`, the lapse derivatives of the
fixed-`Lhat` sum and of the clock-sum integrand, and the response
ratio `R (partial lag / partial R) / (N (partial lag / partial N))`
of each. The rule-D trajectory is not an output.

Test.

1. Calibrate the source. The metric variation of its integrand must
   reproduce the published `T_s` under
   `delta I = (1/2) integral sqrt(-gamma) T^{ij} delta gamma_ij`,
   and the stored energy at `N = 1`, `R = ell` must be
   `E_s = 2 pi zeta_s ell^2 (d theta/dt)^2`. Call
   `contact_variations()` in the shipped curved-fermion checker and
   require that call to pass. A finite difference of the integrand
   must agree with its analytic derivatives.
2. Reduce the fixed-`Lhat` contact at `N = 1`, `R = ell` to the
   controller Lagrangian `4 pi ell^2 c_theta Lhat^{uu} (d theta/dt)^2`.
   Require the round-metric sum to equal the clock-sum integrand, and
   the radius derivatives to agree. Differentiate in `N`.
3. Record the lapse gap and both response ratios. An opposite shift
   `zeta_s -> zeta_s + a`, `sigma -> sigma - a` must move the
   fixed-`Lhat` ratio. The same shift leaves `Lambda` fixed; that
   identity is the shipped `charge_identities()` result and is not a
   new march.
4. Retain `partial Lhat^{uu} / partial N` as a free symbol and record
   the extra lapse piece. Do not set that symbol, `zeta_s`, or
   `sigma` to zero.

A difference of lapse derivatives, with the extra piece unevaluated,
is a scoped obstruction. The class is an uncontrolled approximation.
No separate admissible range is recorded.

## Unequal collars

Assumption. Each collar, taken by itself, has the published
one-boundary energy `E_s = 2 pi zeta_s ell^2 (d theta/dt)^2`. The
candidate for unequal radii is the average of those two energies.
The average is the unique multiple that returns one `E_s` when the
radii agree.

Error comparison. The matched-collar theorem folds identical jets
onto one boundary and states that unequal collar perturbations lie
outside it. The coupled-domain note states that unequal collars
require additional invariants not covered by its product formula.
The average has no jet-cross term. This test does not estimate that
term, and it does not supply a two-collar equation of motion.

Comparison. Write `ell_L = ell_bar + delta/2` and
`ell_R = ell_bar - delta/2`. Compare the average with the retained
one-boundary energy, and compare the average with the same average
plus an unevaluated jet-cross function of `delta`.

Hardest step. The sum of the two published energies is twice `E_s`
when `delta = 0`, so that sum is not the retained contact. The
average matches `E_s` at `delta = 0` and still depends on `delta^2`.
Setting `delta = 0` would hide that dependence. This test does not
do that. A positive second derivative in `delta` is not a stability
mass and is not a range.

Observable. The average as a function of `delta`, its value at
`delta = 0`, and its second derivative. The jet-cross function stays
in the second derivative.

Test.

1. Expand the average. Require agreement with `E_s(ell_bar)` at
   `delta = 0`, and require the sum of the two energies to be twice
   that value.
2. Differentiate twice in `delta`. The second derivative of the
   average must not vanish as a symbol.
3. Add `jet_cross(delta)` and require the second derivative to gain
   `d^2 jet_cross / d delta^2`. Do not replace that derivative by
   zero. Read the matched-collar exclusion from the curved-fermion
   note and the missing-invariant sentence from the coupled-domain
   note.

A pass is a scoped obstruction. The class is an uncontrolled
approximation. No admissible range for `delta` is recorded.

## Nonlocal interacting determinant

Assumption. The free resolvent factor and the interaction share the
clock parameter. In the finite-matrix stand-in,
`G(z) = (1+z) I_2`, the free potential is `v I_2`, and the constant
interaction is `W = [[0, m], [m, 0]]`. The predeclared point is
`z = 0`, `v = 1/100`, `m = 1/2`.

Error comparison. The free static bound controls derivatives of the
free Gaussian energy on the product throat. The substitution under
test treats that bound as a bound on
`d/dz log det(1 + G(z)(v I + W))`. No norm of `W` is supplied by the
retained clock. This matrix is not the wormhole determinant.

Comparison. The free derivative is `d/dz log det(1 + G(z) v)`.
The interacting derivative is the same operation with `W` retained.
The finite-mouth note states that its Gaussian result does not
control the corrected background's nonlocal determinant. The
spherical-anomaly note states that its local sources do not supply
the independent nonlocal determinant. Neither determinant is set to
one, and the free Gaussian energy is not recomputed.

Hardest step. At `m = 0` the two derivatives agree, so the
disagreement is the interaction. A bound on the free derivative at
the predeclared point does not bound the interacting derivative.
The clock polynomial is not given a determinant factor.

Observable. The two derivatives at the predeclared point, and the
second `m` derivative of the interacting derivative at `m = 0`.

Test.

1. Differentiate both logarithms in `z`. Evaluate both at the
   predeclared point. Require the interacting magnitude to exceed
   the free magnitude.
2. Require the free derivative to be independent of `m`, and require
   the second `m` derivative of the interacting derivative at `m = 0`
   to be a nonzero explicit expression.
3. Check the interacting derivative by a finite difference of
   `log det` at the predeclared point.
4. Read the two scope sentences from the finite-mouth note and the
   spherical-anomaly note.

A pass is a scoped obstruction. The class is an uncontrolled
approximation. The physical determinant remains uncomputed.

## `cJ`

Assumption. The matching condition is the shipped parallel electric
susceptibility, solved for `cJ`. The magnetic magnitude `B0` and that
susceptibility are the data of the condition.

Error comparison. The susceptibility is an electromagnetic response
at a stated magnetic background. The rule-D march does not take those
data. Solving the shipped formula does not measure them. The
conditional axial window is a chosen mass tolerance from an earlier
note. It is not a value of this match, and it is not read here.

Comparison. Call `invariants()` in the shipped finite-local checker.
Solve its returned parallel susceptibility for `cJ`. Compare the
symbols in that solution with the symbols in the shipped
`charge_identities()` expression for `Lambda`.

Hardest step. A formula that still contains the susceptibility and
`B0` does not by itself assign `cJ`. The clock sum does not contain
`cJ`. Setting `cJ` to zero would be an extra condition. This test
does not add it.

Observable. The solved matching condition, and the parameter list of
the shipped `march`.

Test.

1. Parse the susceptibility string returned by `invariants()`. Solve
   it for `cJ`. Substitute the solution back and require the
   susceptibility to be recovered.
2. Require `Lambda` from `charge_identities()` to be independent of
   `cJ`, and require the solved `cJ` to be independent of `Lambda`.
3. Require the shipped march parameters to be
   `(eta, lam, kappa, z0, v0, method)`. Record that the
   susceptibility and `B0` are absent from that list.
4. Do not adopt the conditional axial window and do not set `cJ` to
   zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. The matching condition names the susceptibility and `B0`. The
benchmark does not supply them, so no admissible value of `cJ` is
recorded.

## What a pass does not do

These four results do not give every declared-order coefficient an
admissible range. They do not check S074.2–S074.6. They do not start
S075. They do not evaluate the area budget. They do not accept the
epic objective. A coefficient window, including the recorded `Lambda`
intervals, is not that acceptance. The S085 classifications stay the
record of the clock march. This note does not replace that march.

## Benchmark result

The acceptance receipt is
`evidence/source-clock-separated-inputs.json`, SHA-256
`38867d29ca83a15e2911af21d52d4d8dc3b76932914bfee413a348c3f3711b70`.
Both acceptance runs agree. The numbers below are copied from that
receipt. No new `Lambda` interval is recorded. The S085 march was not
repeated.

On the round metric the fixed-`Lhat^{uu}` action and the clock-sum
integrand agree, and their radius derivatives agree. Their lapse
derivatives differ by `4*pi*ell**2*theta_dot**2*sigma`. The source
response ratio is `-2`. The fixed-`Lhat` ratio is `2`. The clock-sum
ratio is `-2`. The combined fixed-`Lhat` ratio is
`2*(sigma + zeta_s)/(sigma - zeta_s)`. The opposite shift that leaves
`Lambda` fixed changes that combined ratio. The unevaluated lapse
piece is `4*pi*N*R**2*c_theta*dLhat_dN*theta_dot**2`. The finite-
difference relative gaps are `2.1690019832717698e-11` for the source
and `3.472477328216686e-10` for the contact. The shipped
`contact_variations()` call passes, with maximum discrepancy
`2.4741694804042424e-11`. `zeta_s` and `sigma` are not set to zero.
The class is an uncontrolled approximation. No separate range is
recorded.

The average of the two published collar energies equals `E_s(ell_bar)`
at `delta = 0`. The sum of those energies is `2*E_s` there. The second
derivative of the average is `pi*theta_dot**2*zeta_s`. The second
derivative of `jet_cross(delta)` is retained. `delta` is not set to
zero, and no two-collar equation is introduced. The class is an
uncontrolled approximation.

At `z = 0`, `v = 1/100`, `m = 1/2`, the free derivative is `2/101` and
the interacting derivative is `-4798/7701`. The second `m` derivative
of the interacting derivative at `m = 0` is
`-4*(z + 1)/(v*z + v + 1)**3`. The finite-difference relative gap of
the interacting derivative is `7.176209872275215e-10`. The determinant
is not set to one, and the wormhole determinant is not computed. The
class is an uncontrolled approximation.

The shipped parallel susceptibility is `1+2*cJ*kappa^4*B0^2`. Solved
for `cJ`, the matching condition is `(chi - 1)/(2*B0**2*kappa**4)`.
The shipped march parameters are `(eta, lam, kappa, z0, v0, method)`.
That list does not supply `chi` or `B0`. `cJ` is not set to zero, and
the conditional axial window is not adopted. The class is an
unavailable physical input.

The area budget is not evaluated. S074.2–S074.6 stay open. S075 through
S079 stay blocked. This is not an impossibility theorem and not a
physical opening.

[Checker](../checks/source_clock_separated_inputs.py),
[receipt](../evidence/source-clock-separated-inputs.json), and
[source/run review](../evidence/source-clock-separated-inputs-review.json).
The joint, radiation, curved-fermion, and finite-local checkers are
called and are not edited. No Lean, paid compute, publication, or
outreach is authorized. The recorded S084 intervals remain the latest
ranges of `Lambda`.
