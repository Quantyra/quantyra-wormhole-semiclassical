# Unmatched inputs on the retained rule-D clock

S085, E003, 2026-10-14. This note takes the four inputs named where
the E003 loop stands and asks what the retained rule-D clock can say
about each one. It does not close S074.2–S074.6, and it does not mark
a physical opening.

The comparison, the hardest step, the observable, and the test for
each input are fixed below, before any checker run. Measured gaps, if
the test produces them, are written only in the benchmark result and
only from the receipt.

The recorded S084 intervals stay the latest ranges of `Lambda`. This
note does not replace them. No coefficient is set to zero. No area
target is written. The area budget stays unevaluated.

## Shared benchmark

The retained realized-drive clock is the equation already checked in
[the joint note](source-clock-joint.md):

```
z_tt + (1 + Lambda r0^2) z + 4 Lambda kappa r0 z^2
    + 3 Lambda kappa^2 z^3 = g/2,
rho = r0 + kappa z,
Lambda = (pi ell^2 / N) (zeta_s + sigma),
sigma = 2 c_theta Lhat^{uu}.
```

`r0(t) = 4 eta [S(t/tau) - S((t-pi)/tau)]`. The sample used here is a
recorded admissible point of that note, read from its receipt: the
positive `Lambda` endpoint on the positive half-`kappa` slice, at the
corner that receipt names. The march is the shipped
`source_clock_radiation.march`. This note does not re-bracket `Lambda`.

## `zeta_s` inside `Lambda`

Comparison. The S084 intervals constrain `zeta_s + sigma`. The
comparator is an opposite shift of the two coefficients. The strongest
statement the clock can make is the identity of that sum.

Hardest step. `lambda_zeta` and `lambda_theta` add before the
variation. A shift that raises one and lowers the other by the same
amount leaves `Lambda` unchanged, and therefore leaves the clock
equation unchanged. A shift of the sum does change the equation. No
retained identity selects a point for `zeta_s` inside the recorded
interval.

Observable. The clock state `(z, z_t)` on the shipped march, at fixed
`kappa` and fixed `Lambda`.

Test. Direct variation: the shift
`zeta_s -> zeta_s + a`, `sigma -> sigma - a` leaves `Lambda` fixed,
and the same statement in the split
`Lambda = lambda_zeta + lambda_theta`. Neither piece of either test
split is zero. Numerical check: two splits with the same sum, both
pieces nonzero, are marched with the shipped integrator and must agree.
The two recorded endpoints of the same `kappa` slice, whose sums
differ, must not agree. That second march is a comparison of recorded
endpoints. It is not a new interval.

A pass of this test is a scoped obstruction. The class is an
unavailable physical input: the clock identifies the sum, which S084
already ranged, and no retained input separates `zeta_s`.

## Unequal collars

Comparison. The curved-collar theorem used for this contact assumes
matched collars: both collars have identical pulled-back metric and
gauge jets, and unequal collar perturbations lie outside that theorem.
The clock coupling uses one sphere radius, `ell`, in
`lambda = pi zeta_s ell^2 / N`.

Hardest step. A collar difference `ell_L - ell_R` would be a second
boundary scale. The retained equation has no such scale. Extending the
theorem is a different derivation. This test does not supply that
derivation, and it does not set the difference to zero.

Observable. The same shipped clock state. The question is whether that
state depends on a collar difference.

Test. Direct variation: differentiate the retained `Lambda` and the
retained acceleration with respect to a collar-difference symbol. Both
derivatives are zero because the symbol is absent. Numerical check: the
shipped march is a function of `(eta, Lambda, kappa, z(0), z_t(0))`.
Its trajectory at the recorded point is compared with the acceleration
reconstructed from that state. The reconstruction has no collar
difference. Agreement means the observable that was actually integrated
does not carry unequal collars.

A pass is a scoped obstruction. The class is an unresolved proof: the
matched-collar theorem does not cover unequal driven collars, and the
retained clock does not contain their entry.

## Nonlocal interacting determinant

Comparison. The local relative Weyl transfer is not a computation of
the finite nonlocal determinant. The free Gaussian determinant on the
static product throat does not control the nonlocal determinant of the
corrected interacting background. The retained clock is a local
ordinary differential equation.

Hardest step. A nonlocal determinant is a functional of the background,
not a constant already sitting in `Lambda`. No formula in the retained
equation adds that functional. This test does not evaluate the free
Gaussian determinant again, and it does not set the interacting
determinant to one.

Observable. The same shipped clock state. The question is whether that
state depends on the interacting determinant.

Test. Direct variation: the derivative of the retained acceleration
with respect to a determinant symbol is zero. Numerical check: the
acceleration reconstructed from the shipped march agrees with the local
polynomial and has no determinant factor.

A pass is a scoped obstruction. The class is an unresolved proof: the
finite nonlocal interacting determinant has not been entered on this
benchmark.

## `cJ`

Comparison. The finite operator
`Delta L = cJ kappa_em^4 J^2 / 16` is a different direction from the
clock contact `I_s = -(zeta_s / 2) integral (D theta)^2`. On a purely
magnetic background, `J = 0`, so that operator and its first field
variation vanish for every finite `cJ`. The background therefore does
not measure `cJ`. The earlier conditional axial window is a chosen mass
tolerance. It is not a range of this clock, and it is not reused here
as an admissible value.

Hardest step. The clock acceleration depends on `Lambda` and `kappa`.
It does not depend on `cJ`. Setting `cJ = 0` would be an extra matching
condition. This test does not add that condition.

Observable. The same shipped clock state, together with the algebraic
vanishing of `Delta L` at `J = 0` for a free symbol `cJ`.

Test. Direct variation: `Delta L` and `d(Delta L)/dJ` are zero at
`J = 0` while `cJ` remains free, and the clock acceleration is
independent of `cJ`. The two densities `J^2` and `(D theta)^2` are not
the same monomial. Numerical check: the reconstructed acceleration
matches the shipped march and does not move with `cJ`.

A pass is a scoped obstruction. The class is an unavailable physical
input: neither the magnetic background nor this clock identifies `cJ`.

## What a pass does not do

These four results do not give every declared-order coefficient an
admissible range. They do not check S074.2–S074.6. They do not start
S075. They do not evaluate the area budget. They do not accept the epic
objective. A coefficient window, including the recorded `Lambda`
intervals, is not that acceptance.

## Benchmark result

The acceptance receipt is `evidence/source-clock-unmatched-inputs.json`.
Both acceptance runs agree. The numbers below are copied from that
receipt. No new `Lambda` interval is recorded.

The shift identity holds. An opposite change of `zeta_s` and `sigma`,
and the same change of `lambda_zeta` and `lambda_theta`, leaves
`Lambda` fixed. At the recorded point read from the S084 receipt,
`kappa = 2.5360074639320374e-05`,
`Lambda = 5497.42431640625`, corner `z(0) = 7/25`,
`z_t(0) = 1/2000`, the two nonzero splits
`(0.25 Lambda, 0.75 Lambda)` and `(0.75 Lambda, 0.25 Lambda)` produce
the same shipped trajectory. Their maximum absolute state gap is `0.0`.
The other recorded endpoint of that slice,
`Lambda = -3957.533203125`, is a different sum. Its gap from the
positive endpoint is `0.009895015671797413`. The clock distinguishes
recorded sums and does not distinguish splits. `zeta_s` is an
unavailable physical input inside the recorded sum. It is not set to
zero, and `sigma` is not set to zero.

The acceleration reconstructed from the shipped state agrees with the
local polynomial. The maximum absolute gap is
`2.1894991625165394e-06`. The derivatives of that polynomial with
respect to a collar difference, a nonlocal-determinant symbol, and
`cJ` are zero. The derivative with respect to `Lambda` is not
identically zero. The shipped march takes only
`(eta, Lambda, kappa, z(0), z_t(0), method)`. Unequal collars and the
nonlocal interacting determinant are unresolved proofs on this
benchmark. They are not assigned a zero cost.

At `J = 0`, `Delta L` and its first derivative in `J` vanish for a free
symbol `cJ`. The clock contact is not that operator. `cJ` is an
unavailable physical input. It is not set to zero, and the earlier
conditional axial window is not adopted as a clock range.

The area budget is not evaluated. S074.2–S074.6 stay open. S075 through
S079 stay blocked. This is not an impossibility theorem and not a
physical opening.

[Checker](../checks/source_clock_unmatched_inputs.py),
[receipt](../evidence/source-clock-unmatched-inputs.json), and
[source/run review](../evidence/source-clock-unmatched-inputs-review.json).
The joint, radiation, controller, and allocation checkers are reused
and are not edited. No Lean, paid compute, publication, or outreach
was used.
