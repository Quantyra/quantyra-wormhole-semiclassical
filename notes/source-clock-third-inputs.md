# Third mu0 derivative, second lapse derivative, width derivative in h, and a vanishing slope

S107, E003, 2026-10-15. This note takes the four inputs named by S106.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S106 second mu0 derivative is not recomputed
as a new result. The S106 f derivative is not recomputed as a new
result. The S106 h derivative is not recomputed as a new result. The
S106 second B0 derivative is not recomputed as a new result.
No new `Lambda` interval is sought. No replacement area target is
written before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Third `mu0` derivative

Assumption. A number for `c_theta(mu0)` would be the third `mu0`
derivative `1/(10*pi**2*mu0**3)`.

Error comparison. That derivative is not a value of `c_theta(mu0)`. It
keeps `mu0`. It is not the second derivative
`-1/(20*pi**2*mu0**2)`. `c_theta` is not set to zero and is not set
equal to `b_theta`. The second `mu0` derivative is not recomputed.

Comparison. Require the shipped `b_theta` formula. Differentiate
`-1/(20*pi**2*mu0**2)` in `mu0` and check the derivative by the limit
of a central quotient and by a finite difference. Compare the result
with the endpoint that still contains `c_theta(mu0)`.

Hardest step. The third derivative can be exact while `mu0` remains in
the denominator. A quotient that still contains `mu0` does not assign
`c_theta(mu0)`.

Observable. Whether `1/(10*pi**2*mu0**3)` is a number, and whether the
endpoint still contains `c_theta(mu0)`.

Test.

1. Require `b_theta=4/[5(4pi)^2]`.
2. Record the derivative of `-1/(20*pi**2*mu0**2)` in `mu0`. Require
   the limit of a central quotient and a finite difference to agree.
   Require the derivative to equal `1/(10*pi**2*mu0**3)`. Require it to
   keep `mu0` and to differ from `c_theta(mu0)`. Require the endpoint to
   keep `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not store the second derivative
   `-1/(20*pi**2*mu0**2)` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Second `f` derivative of the bracket ratio

Assumption. The unequal-collar invariant would be the `f` derivative
`2/f**2` of `-2/f`.

Error comparison. That derivative is not a two-collar theorem. It keeps
`f`. It is not the recorded collar mismatch. The local relative limit is
not a proof for unequal driven collars. `delta` is not set to zero.
`d_theta` is not set to zero. The derivative `-2/f` is not recomputed.

Comparison. Require both displayed brackets, the shipped `b_theta` and
`d_theta` formulas, and the exclusion of unequal collars. Differentiate
`-2/f` in `f` and check the derivative by the limit of a central
quotient and by a finite difference.

Hardest step. Another derivative of one lapse can be exact while the
second collar is absent. A quotient that still contains `f` is not the
missing invariant.

Observable. The `f` derivative of `-2/f`, and whether it equals the
mismatch.

Test.

1. Require `b_theta (log f/2)(2/3)(R^2/f)w`, `-d_theta (R^2/f)f'/2`,
   both shipped coefficients, and the statement that this limit is not a
   proof for unequal driven collars.
2. Record the derivative of `-2/f` in `f`. Require it to equal
   `2/f**2`. Require the limit of a central quotient and a finite
   difference to agree. Require the derivative to differ from `1`.
3. Do not adopt the derivative as a two-collar theorem. Do not
   differentiate the RN relative series. Do not store `-2/f` as this
   result. Do not set `d_theta` to zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## `h` derivative of `4/h**3`

Assumption. A number for one `D_j` would be the `h` derivative
`-12/h**4`.

Error comparison. That derivative is not an assignment of `D_j`. It
keeps `h`. It is not `4/h**3`. `a` and `h` are not set to zero. `delta`
is not set to zero. The determinant is not set to one. The `h`
derivative `4/h**3` is not recomputed.

Comparison. Require the displayed ceiling. Differentiate `4/h**3` in
`h` and check the derivative by the limit of a central quotient and by
a finite difference. Differentiate the result in `D_j`, and
differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The derivative in `h` can be exact while `D_j` is absent.
A quotient that still contains `h` does not assign `D_j`.

Observable. Whether `-12/h**4` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the `h` derivative of `4/h**3` to equal `-12/h**4`. Require
   the limit of a central quotient and a finite difference to agree.
   Require the derivative in `D_j` to vanish. Require the clock
   derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read `-12/h**4` as `D_j`. Do not store `4/h**3`
   as this result. Do not set `a` or `h` to zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## `B0` derivative of the constant `8`

Assumption. A measured `cJ` would be the `B0` derivative `0` of the
constant `8`.

Error comparison. The derivative `0` does not measure `cJ`. It does not
contain `cJ`. `cJ` is not set to zero. `B0` is not set to one. The
conditional axial window is not adopted. A positive quadratic term is
not used as a stability mass. The integer `8` is not stored as this
result.

Comparison. Read the shipped parallel susceptibility. Reduce it at
gravitational `kappa = 1` and substitute the integer `2` for `cJ`.
Differentiate that expression three times in `B0` and check the third
derivative by a difference of second differences.

Hardest step. The third derivative can vanish while `cJ` is absent. A
vanishing derivative is not a measured coupling, and it does not set
`cJ` to zero.

Observable. The `B0` derivative of the constant second derivative, and
whether that value is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the third `B0` derivative. Require it to equal `0`. Require a
   difference of second differences to agree. Require a nonzero `B0`
   other than one. Require the `cJ` derivative of that value to vanish.
   Do not set `B0` to one. Do not set `cJ` to zero.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store `8` as this result. Do not use the positive
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
third derivative of `-1/(20*pi**2*mu0**2)` in `mu0` is
`1/(10*pi**2*mu0**3)`. It is not a number and does not assign
`c_theta(mu0)`. The endpoint keeps `c_theta(mu0)`. `sigma` stays
`2*c_theta*Lhat^{uu}`. The class is an unavailable physical input.
`c_theta` is not set to zero and is not set equal to `b_theta`. The
second `mu0` derivative is not recomputed as a new result.

The `f` derivative of `-2/f` is `2/f**2`. It is not the recorded
mismatch `1`. The class is an unresolved proof. The derivative is not
adopted as a two-collar theorem. `delta` is not set to zero. `d_theta`
is not set to zero. The derivative `-2/f` is not recomputed as a new
result.

The `h` derivative of `4/h**3` is `-12/h**4`. It is not a number. It
keeps `h`. Its derivative in `D_j` is `0`. The clock derivative in
`D_j` is `0`. `numeric_D_j` is null. The class is an unavailable
physical input. The derivative is not read as `D_j`. The derivative
`4/h**3` is not recomputed as a new result.

The third `B0` derivative of the substituted susceptibility is `0`.
The derivative `0` does not measure `cJ`. Its derivative in `cJ` is
`0`. The class is an unavailable physical input. `cJ` is not set to
zero. `B0` is not set to one. A positive quadratic term is not used as
a stability mass. The integer `8` is not stored as this result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
