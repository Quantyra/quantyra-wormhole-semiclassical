# Tenth mu0 derivative, ninth lapse derivative, width derivative in h, and a repeated vanishing slope

S114, E003, 2026-10-15. This note takes the four inputs named by S113.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S113 ninth mu0 derivative is not recomputed
as a new result. The S113 f derivative is not recomputed as a new
result. The S113 h derivative is not recomputed as a new result. The
S113 ninth B0 derivative is not recomputed as a new result.
No new `Lambda` interval is sought. No replacement area target is
written before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Tenth `mu0` derivative

Assumption. A number for `c_theta(mu0)` would be the tenth `mu0`
derivative `-18144/(pi**2*mu0**10)`.

Error comparison. That derivative is not a value of `c_theta(mu0)`. It
keeps `mu0`. It is not the ninth derivative `2016/(pi**2*mu0**9)`.
`c_theta` is not set to zero and is not set equal to `b_theta`. The
ninth `mu0` derivative is not recomputed.

Comparison. Require the shipped `b_theta` formula. Differentiate
`2016/(pi**2*mu0**9)` in `mu0` and check the derivative by the limit of
a central quotient and by a finite difference. Compare the result with
the endpoint that still contains `c_theta(mu0)`.

Hardest step. The tenth derivative can be exact while `mu0` remains in
the denominator. A quotient that still contains `mu0` does not assign
`c_theta(mu0)`.

Observable. Whether `-18144/(pi**2*mu0**10)` is a number, and whether the
endpoint still contains `c_theta(mu0)`.

Test.

1. Require `b_theta=4/[5(4pi)^2]`.
2. Record the derivative of `2016/(pi**2*mu0**9)` in `mu0`. Require the
   limit of a central quotient and a finite difference to agree. Require
   the derivative to equal `-18144/(pi**2*mu0**10)`. Require it to keep
   `mu0` and to differ from `c_theta(mu0)`. Require the endpoint to keep
   `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not store the ninth derivative
   `2016/(pi**2*mu0**9)` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Ninth `f` derivative of the bracket ratio

Assumption. The unequal-collar invariant would be the `f` derivative
`-80640/f**9` of `10080/f**8`.

Error comparison. That derivative is not a two-collar theorem. It keeps
`f`. It is not the recorded collar mismatch. The local relative limit is
not a proof for unequal driven collars. `delta` is not set to zero.
`d_theta` is not set to zero. The derivative `10080/f**8` is not recomputed.

Comparison. Require both displayed brackets, the shipped `b_theta` and
`d_theta` formulas, and the exclusion of unequal collars. Differentiate
`10080/f**8` in `f` and check the derivative by the limit of a central
quotient and by a finite difference.

Hardest step. Another derivative of one lapse can be exact while the
second collar is absent. A quotient that still contains `f` is not the
missing invariant.

Observable. The `f` derivative of `10080/f**8`, and whether it equals the
mismatch.

Test.

1. Require `b_theta (log f/2)(2/3)(R^2/f)w`, `-d_theta (R^2/f)f'/2`,
   both shipped coefficients, and the statement that this limit is not a
   proof for unequal driven collars.
2. Record the derivative of `10080/f**8` in `f`. Require it to equal
   `-80640/f**9`. Require the limit of a central quotient and a finite
   difference to agree. Require the derivative to differ from `1`.
3. Do not adopt the derivative as a two-collar theorem. Do not
   differentiate the RN relative series. Do not store `10080/f**8` as this
   result. Do not set `d_theta` to zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## `h` derivative of `-725760/h**10`

Assumption. A number for one `D_j` would be the `h` derivative
`7257600/h**11`.

Error comparison. That derivative is not an assignment of `D_j`. It
keeps `h`. It is not `-725760/h**10`. `a` and `h` are not set to zero.
`delta` is not set to zero. The determinant is not set to one. The `h`
derivative `-725760/h**10` is not recomputed.

Comparison. Require the displayed ceiling. Differentiate `-725760/h**10`
in `h` and check the derivative by the limit of a central quotient and by
a finite difference. Differentiate the result in `D_j`, and
differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The derivative in `h` can be exact while `D_j` is absent.
A quotient that still contains `h` does not assign `D_j`.

Observable. Whether `7257600/h**11` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the `h` derivative of `-725760/h**10` to equal `7257600/h**11`.
   Require the limit of a central quotient and a finite difference to
   agree. Require the derivative in `D_j` to vanish. Require the clock
   derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read `7257600/h**11` as `D_j`. Do not store
   `-725760/h**10` as this result. Do not set `a` or `h` to zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## `B0` derivative of the ninth derivative `0`

Assumption. A measured `cJ` would be the `B0` derivative `0` of the
ninth derivative `0`.

Error comparison. The derivative `0` does not measure `cJ`. It does not
contain `cJ`. `cJ` is not set to zero. `B0` is not set to one. The
conditional axial window is not adopted. A positive quadratic term is
not used as a stability mass. The ninth derivative `0` is not stored as
this result.

Comparison. Read the shipped parallel susceptibility. Reduce it at
gravitational `kappa = 1` and substitute the integer `2` for `cJ`.
Differentiate that expression ten times in `B0` and check the tenth
derivative by a difference of ninth differences.

Hardest step. The tenth derivative can vanish while `cJ` is absent. A
repeated vanishing derivative is not a measured coupling, and it does
not set `cJ` to zero.

Observable. The `B0` derivative of the ninth derivative, and whether
that value is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the tenth `B0` derivative. Require it to equal `0`. Require a
   difference of ninth differences to agree. Require a nonzero `B0`
   other than one. Require the `cJ` derivative of that value to vanish.
   Do not set `B0` to one. Do not set `cJ` to zero.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store the ninth derivative `0` as this result. Do not
   use the positive quadratic as a stability mass.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
tenth derivative of `2016/(pi**2*mu0**9)` in `mu0` is
`-18144/(pi**2*mu0**10)`. It is not a number and does not assign
`c_theta(mu0)`. The endpoint keeps `c_theta(mu0)`. `sigma` stays
`2*c_theta*Lhat^{uu}`. The class is an unavailable physical input.
`c_theta` is not set to zero and is not set equal to `b_theta`. The
ninth `mu0` derivative is not recomputed as a new result.

The `f` derivative of `10080/f**8` is `-80640/f**9`. It is not the recorded
mismatch `1`. The class is an unresolved proof. The derivative is not
adopted as a two-collar theorem. `delta` is not set to zero. `d_theta`
is not set to zero. The derivative `10080/f**8` is not recomputed as a new
result.

The `h` derivative of `-725760/h**10` is `7257600/h**11`. It is not a number.
It keeps `h`. Its derivative in `D_j` is `0`. The clock derivative in
`D_j` is `0`. `numeric_D_j` is null. The class is an unavailable
physical input. The derivative is not read as `D_j`. The derivative
`-725760/h**10` is not recomputed as a new result.

The tenth `B0` derivative of the substituted susceptibility is `0`.
The derivative `0` does not measure `cJ`. Its derivative in `cJ` is
`0`. The class is an unavailable physical input. `cJ` is not set to
zero. `B0` is not set to one. A positive quadratic term is not used as
a stability mass. The ninth derivative `0` is not stored as this result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
