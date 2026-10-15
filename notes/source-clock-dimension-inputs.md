# Dimension-one coefficient, second bracket, strength ceiling, and the ratio `(2*B0**2)/B0**2`

S103, E003, 2026-10-15. This note takes the four inputs named by S102.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S102 curvature derivative is not recomputed as
a new result. The S102 radial derivative is not recomputed as a new
result. The S102 integral strength is not recomputed as a new result.
The S102 factor derivative is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written before
the test. The area budget stays unevaluated. `zeta_s`, `sigma`, `cJ`,
`c_theta`, `d_theta`, and `delta` are not set to zero. An endpoint
contact energy is not cancelled. A positive quadratic term is not used
as a stability mass.

## Dimension-one coefficient

Assumption. A number for `c_theta(mu0)` would be the statement that the
subtraction renormalizes a dimension-one coefficient.

Error comparison. That statement is not a value of `c_theta(mu0)`. The
clause contains no numeral. A scale derivative `c_one` does not assign
`c_theta(mu0)`. `c_theta` is not set to zero and is not set equal to
`b_theta`. The curvature derivative of the subtraction is not recomputed.

Comparison. Require the displayed clause. Form the dimension-one monomial
`c_one*mu`. Differentiate it in the scale `mu` and check the derivative
by a finite difference. Compare the result with the endpoint that still
contains `c_theta(mu0)`.

Hardest step. The dimension can be carried by the power of the scale
while the coefficient stays symbolic. A derivative that returns the
coefficient does not assign `c_theta(mu0)`.

Observable. Whether the dimension-one clause is a number, and the scale
derivative of `c_one*mu`.

Test.

1. Require the sentence `it renormalizes a dimension-one coefficient`.
2. Require that clause to contain no numeral and no equality. Record the
   derivative of `c_one*mu` in `mu`. Require a finite difference to agree.
   Require the endpoint to keep `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not differentiate the subtraction in a curvature
   symbol.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Second bracket in (4)

Assumption. The unequal-collar invariant would be the displayed
second-bracket term `-d_theta (R^2/f)f'/2`.

Error comparison. That bracket term is not a two-collar theorem. The
local relative limit is not a proof for unequal driven collars. The
derivative in `f'` is `-R**2/(120*pi**2*f)`. That derivative is not the
recorded collar mismatch and is not the radial derivative
`-1/(60*pi**2)`. `delta` is not set to zero. `d_theta` is not set to zero.

Comparison. Require the displayed bracket, the shipped coefficient
`d_theta=4/[15(4pi)^2]`, and the exclusion of unequal collars.
Differentiate `-d_theta*(R**2/f)*f_prime/2` in `f_prime` and check the
derivative by a finite difference. Compare the derivative with the
recorded mismatch `1`.

Hardest step. The bracket can stay dependent on one radius and one lapse
while the two collars remain uncovered. A one-boundary quotient is not
the missing invariant.

Observable. The derivative of the second bracket in `f'`, and whether it
equals the mismatch.

Test.

1. Require `-d_theta (R^2/f)f'/2`, the shipped `d_theta` formula, and the
   statement that this limit is not a proof for unequal driven collars.
2. Record the derivative in `f_prime`. Require a finite difference to
   agree. Require the derivative to differ from `1` and from
   `-1/(60*pi**2)`.
3. Do not adopt the derivative as a two-collar theorem. Do not
   differentiate the RN relative series. Do not store the radial
   derivative of `-d_theta R` as this result. Do not set `d_theta` to
   zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## Strength ceiling

Assumption. A number for one `D_j` would be the ceiling `2*a/h`.

Error comparison. That ceiling is not an assignment of `D_j`. It keeps
`a` and `h`. It is not the integral strength. `a` and `h` are not set to
zero. `delta` is not set to zero. The determinant is not set to one.

Comparison. Require the displayed ceiling. Form `2*a/h`. Differentiate
that expression in `D_j` and check the derivative by a finite difference.
Differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The derivative can vanish because the ceiling does not
contain `D_j`. A quotient that still contains `a` and `h` does not assign
`D_j`.

Observable. Whether `2*a/h` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the ceiling to keep `a` and `h`. Require its derivative in
   `D_j` to vanish, and require a finite difference to agree. Require the
   clock derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read the ceiling as `D_j`. Do not store the integral
   strength as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Ratio `(2*B0**2)/B0**2`

Assumption. A measured `cJ` would be the ratio `(2*B0**2)/B0**2`.

Error comparison. The integer `2` does not measure `cJ`. Cancelling
`B0**2` does not supply a coupling. The conditional axial window is not
adopted. `cJ` is not set to zero. `B0` is not set to one. A positive
quadratic term is not used as a stability mass.

Comparison. Form `(2*B0**2)/B0**2`. Simplify the ratio and check it at
one nonzero `B0` other than one. Differentiate the unsimplified ratio in
`B0` and check that derivative by a finite difference.

Hardest step. The ratio can collapse to an integer while `cJ` never
enters. The integer is not a measured coupling, and setting `B0` to one
is not this test.

Observable. The value of the ratio, and whether that value is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the simplified ratio. Require a numerical evaluation at one
   nonzero `B0` other than one to agree. Require the `B0` derivative of
   the unsimplified ratio to vanish, and require a finite difference to
   agree. Do not set `B0` to one.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store the `cJ` derivative of `2*B0**2` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
dimension-one clause is not a number and does not assign `c_theta(mu0)`.
The scale derivative of `c_one*mu` is `c_one`. The endpoint keeps
`c_theta(mu0)`. `sigma` stays `2*c_theta*Lhat^{uu}`. The class is an
unavailable physical input. `c_theta` is not set to zero and is not set
equal to `b_theta`. The curvature derivative is not recomputed as a new
result.

The derivative of the second bracket in `f'` is
`-R**2/(120*pi**2*f)`. It is not the recorded mismatch `1` and it is not
the radial derivative `-1/(60*pi**2)`. The class is an unresolved proof.
The derivative is not adopted as a two-collar theorem. `delta` is not
set to zero. `d_theta` is not set to zero. The radial derivative is not
recomputed as a new result.

The ceiling `2*a/h` is not a number. It keeps `a` and `h`. Its derivative
in `D_j` is `0`. The clock derivative in `D_j` is `0`. `numeric_D_j` is
null. The class is an unavailable physical input. The ceiling is not read
as `D_j`. The integral strength is not recomputed as a new result.

The ratio `(2*B0**2)/B0**2` simplifies to `2`. The integer `2` does not
measure `cJ`. The `B0` derivative of the unsimplified ratio is `0`. The
class is an unavailable physical input. `cJ` is not set to zero. `B0` is
not set to one. The factor derivative is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
