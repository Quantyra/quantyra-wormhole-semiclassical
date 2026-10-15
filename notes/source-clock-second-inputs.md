# Second mu0 derivative, lapse derivative of the ratio, width derivative in h, and the constant 8

S106, E003, 2026-10-15. This note takes the four inputs named by S105.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S105 mu0 derivative is not recomputed as a new
result. The S105 ratio is not recomputed as a new result. The S105
amplitude derivative is not recomputed as a new result. The S105
background derivative is not recomputed as a new result.
No new `Lambda` interval is sought. No replacement area target is
written before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Second `mu0` derivative

Assumption. A number for `c_theta(mu0)` would be the second `mu0`
derivative `-1/(20*pi**2*mu0**2)`.

Error comparison. That derivative is not a value of `c_theta(mu0)`. It
keeps `mu0`. It is not the first derivative `1/(20*pi**2*mu0)`.
`c_theta` is not set to zero and is not set equal to `b_theta`. The
first `mu0` derivative is not recomputed.

Comparison. Require the shipped `b_theta` formula. Differentiate
`1/(20*pi**2*mu0)` in `mu0` and check the derivative by the limit of a
central quotient and by a finite difference. Compare the result with the
endpoint that still contains `c_theta(mu0)`.

Hardest step. The second derivative can be exact while `mu0` remains in
the denominator. A quotient that still contains `mu0` does not assign
`c_theta(mu0)`.

Observable. Whether `-1/(20*pi**2*mu0**2)` is a number, and whether the
endpoint still contains `c_theta(mu0)`.

Test.

1. Require `b_theta=4/[5(4pi)^2]`.
2. Record the derivative of `1/(20*pi**2*mu0)` in `mu0`. Require the
   limit of a central quotient and a finite difference to agree. Require
   the derivative to equal `-1/(20*pi**2*mu0**2)`. Require it to keep
   `mu0` and to differ from `c_theta(mu0)`. Require the endpoint to keep
   `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not store the first derivative
   `1/(20*pi**2*mu0)` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## `f` derivative of the bracket ratio

Assumption. The unequal-collar invariant would be the `f` derivative
`-2/f` of the bracket ratio.

Error comparison. That derivative is not a two-collar theorem. It keeps
`f`. It is not the recorded collar mismatch. The local relative limit is
not a proof for unequal driven collars. `delta` is not set to zero.
`d_theta` is not set to zero. The ratio `-2*log(f)` is not recomputed.

Comparison. Require both displayed brackets, the shipped `b_theta` and
`d_theta` formulas, and the exclusion of unequal collars. Differentiate
`-2*log(f)` in `f` and check the derivative by the limit of a central
quotient and by a finite difference.

Hardest step. The derivative of one logarithm can be exact on one lapse
while the second collar is absent. A quotient that still contains `f`
is not the missing invariant.

Observable. The `f` derivative of the bracket ratio, and whether it
equals the mismatch.

Test.

1. Require `b_theta (log f/2)(2/3)(R^2/f)w`, `-d_theta (R^2/f)f'/2`,
   both shipped coefficients, and the statement that this limit is not a
   proof for unequal driven collars.
2. Record the derivative of `-2*log(f)` in `f`. Require it to equal
   `-2/f`. Require the limit of a central quotient and a finite
   difference to agree. Require the derivative to differ from `1`.
3. Do not adopt the derivative as a two-collar theorem. Do not
   differentiate the RN relative series. Do not store the ratio
   `-2*log(f)` as this result. Do not set `d_theta` to zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## `h` derivative of the amplitude derivative

Assumption. A number for one `D_j` would be the `h` derivative
`4/h**3`.

Error comparison. That derivative is not an assignment of `D_j`. It
keeps `h`. It is not the amplitude derivative `-2/h**2`. `a` and `h`
are not set to zero. `delta` is not set to zero. The determinant is not
set to one. The amplitude derivative is not recomputed.

Comparison. Require the displayed ceiling. Differentiate `-2/h**2` in
`h` and check the derivative by the limit of a central quotient and by
a finite difference. Differentiate the result in `D_j`, and
differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The derivative in `h` can be exact while `D_j` is absent.
A quotient that still contains `h` does not assign `D_j`.

Observable. Whether `4/h**3` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the `h` derivative of `-2/h**2` to equal `4/h**3`. Require
   the limit of a central quotient and a finite difference to agree.
   Require the derivative in `D_j` to vanish. Require the clock
   derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read `4/h**3` as `D_j`. Do not store the amplitude
   derivative as this result. Do not set `a` or `h` to zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Second `B0` derivative

Assumption. A measured `cJ` would be the second `B0` derivative `8`.

Error comparison. The integer `8` does not measure `cJ`. It does not
contain `cJ`. The conditional axial window is not adopted. `cJ` is not
set to zero. `B0` is not set to one. A positive quadratic term is not
used as a stability mass. The slope `8*B0` is not stored as this result.

Comparison. Read the shipped parallel susceptibility. Reduce it at
gravitational `kappa = 1` and substitute the integer `2` for `cJ`.
Differentiate that expression twice in `B0` and check the second
derivative by a second finite difference.

Hardest step. The second derivative can be the integer `8` while `cJ`
is absent. That integer is not a measured coupling, and a positive
quadratic is not a stability mass.

Observable. The second `B0` derivative of the substituted
susceptibility, and whether that integer is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the second `B0` derivative. Require it to equal `8`. Require a
   second finite difference to agree. Require a nonzero `B0` other than
   one. Require the `cJ` derivative of `8` to vanish. Do not set `B0` to
   one.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store `8*B0` as this result. Do not use the positive
   quadratic as a stability mass.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
second derivative of `1/(20*pi**2*mu0)` in `mu0` is
`-1/(20*pi**2*mu0**2)`. It is not a number and does not assign
`c_theta(mu0)`. The endpoint keeps `c_theta(mu0)`. `sigma` stays
`2*c_theta*Lhat^{uu}`. The class is an unavailable physical input.
`c_theta` is not set to zero and is not set equal to `b_theta`. The
first `mu0` derivative is not recomputed as a new result.

The `f` derivative of the bracket ratio is `-2/f`. It is not the
recorded mismatch `1`. The class is an unresolved proof. The derivative
is not adopted as a two-collar theorem. `delta` is not set to zero.
`d_theta` is not set to zero. The ratio is not recomputed as a new
result.

The `h` derivative of `-2/h**2` is `4/h**3`. It is not a number. It
keeps `h`. Its derivative in `D_j` is `0`. The clock derivative in
`D_j` is `0`. `numeric_D_j` is null. The class is an unavailable
physical input. The derivative is not read as `D_j`. The amplitude
derivative is not recomputed as a new result.

The second `B0` derivative of the substituted susceptibility is `8`.
The integer `8` does not measure `cJ`. Its derivative in `cJ` is `0`.
The class is an unavailable physical input. `cJ` is not set to zero.
`B0` is not set to one. A positive quadratic term is not used as a
stability mass. The slope `8*B0` is not stored as this result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
