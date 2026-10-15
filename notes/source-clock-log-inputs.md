# Logarithmic endpoint, bracket ratio, amplitude derivative, and `8*B0`

S105, E003, 2026-10-15. This note takes the four inputs named by S104.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S104 coefficient derivative is not recomputed
as a new result. The S104 weight derivative is not recomputed as a new
result. The S104 width derivative is not recomputed as a new result.
The S104 substituted susceptibility is not recomputed as a new result.
No new `Lambda` interval is sought. No replacement area target is
written before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Logarithmic endpoint term

Assumption. A number for `c_theta(mu0)` would be the logarithmic
endpoint term `log(mu0)/(20*pi**2)`.

Error comparison. That term is not a value of `c_theta(mu0)`. It keeps
`mu0`. Its derivative in `mu0` is `1/(20*pi**2*mu0)`. That derivative
is not a number and does not assign `c_theta(mu0)`. `c_theta` is not
set to zero and is not set equal to `b_theta`. The coefficient
derivative that returns `mu` is not recomputed.

Comparison. Require the shipped `b_theta` formula. Form
`log(mu0)/(20*pi**2)`. Differentiate it in `mu0` and check the
derivative by the limit of a central quotient and by a finite
difference. Compare the result with the endpoint that still contains
`c_theta(mu0)`.

Hardest step. The logarithm can be an exact term in the endpoint while
the integration constant stays independent. A derivative that still
contains `mu0` does not assign `c_theta(mu0)`.

Observable. Whether `log(mu0)/(20*pi**2)` is a number, and its
derivative in `mu0`.

Test.

1. Require `b_theta=4/[5(4pi)^2]`.
2. Record the derivative of `log(mu0)/(20*pi**2)` in `mu0`. Require the
   limit of a central quotient and a finite difference to agree. Require
   the derivative to keep `mu0` and to differ from `c_theta(mu0)`.
   Require the endpoint to keep `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not store the coefficient derivative `mu` as
   this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Ratio of the bracket derivatives

Assumption. The unequal-collar invariant would be the ratio of the
first-bracket `w` derivative to the second-bracket `f'` derivative.

Error comparison. That ratio is not a two-collar theorem. It simplifies
to `-2*log(f)`. That value is not the recorded collar mismatch. The
local relative limit is not a proof for unequal driven collars. `delta`
is not set to zero. `d_theta` is not set to zero.

Comparison. Require both displayed brackets, the shipped `b_theta` and
`d_theta` formulas, and the exclusion of unequal collars. Form the ratio
of the two one-boundary derivatives and check it by a numerical
evaluation at one lapse other than one. Differentiate the ratio in `R`
and check that derivative by a finite difference.

Hardest step. The ratio can cancel the radius and the lapse while one
logarithm of one lapse remains. A one-boundary logarithm is not the
missing invariant.

Observable. The value of the ratio, and whether it equals the mismatch.

Test.

1. Require `b_theta (log f/2)(2/3)(R^2/f)w`, `-d_theta (R^2/f)f'/2`,
   both shipped coefficients, and the statement that this limit is not a
   proof for unequal driven collars.
2. Record the simplified ratio. Require a numerical evaluation at one
   lapse other than one to agree. Require the derivative in `R` to
   vanish, and require a finite difference to agree. Require the ratio
   to differ from `1`.
3. Do not adopt the ratio as a two-collar theorem. Do not differentiate
   the RN relative series. Do not store either bracket derivative as
   this result. Do not set `d_theta` to zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## Amplitude derivative of the width derivative

Assumption. A number for one `D_j` would be the `a` derivative
`-2/h**2`.

Error comparison. That derivative is not an assignment of `D_j`. It
keeps `h`. It is not the width derivative `-2*a/h**2`. `a` and `h` are
not set to zero. `delta` is not set to zero. The determinant is not set
to one.

Comparison. Require the displayed ceiling. Differentiate `-2*a/h**2` in
`a` and check the derivative by a finite difference. Differentiate the
result in `D_j`, and differentiate the clock expression for `Lambda` in
`D_j`.

Hardest step. The amplitude derivative can be exact while `D_j` is
absent. A quotient that still contains `h` does not assign `D_j`.

Observable. Whether `-2/h**2` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the `a` derivative of `-2*a/h**2` to equal `-2/h**2`.
   Require a finite difference to agree. Require the derivative in
   `D_j` to vanish. Require the clock derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read `-2/h**2` as `D_j`. Do not store the width
   derivative as this result. Do not set `a` or `h` to zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## `B0` derivative of the substituted susceptibility

Assumption. A measured `cJ` would be the `B0` derivative `8*B0`.

Error comparison. The expression `8*B0` does not measure `cJ`. It keeps
`B0`. The conditional axial window is not adopted. `cJ` is not set to
zero. `B0` is not set to one. A positive quadratic term is not used as
a stability mass.

Comparison. Read the shipped parallel susceptibility. Reduce it at
gravitational `kappa = 1` and substitute the integer `2` for `cJ`.
Differentiate that expression in `B0` and check the derivative by a
finite difference.

Hardest step. The derivative can be linear in `B0` while `cJ` is
absent. That slope is not a measured coupling, and setting `B0` to one
is not this test.

Observable. The `B0` derivative of the substituted susceptibility, and
whether that derivative is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the `B0` derivative. Require it to equal `8*B0`. Require a
   finite difference to agree. Require a nonzero `B0` other than one.
   Require the `cJ` derivative of `8*B0` to vanish. Do not set `B0` to
   one.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store `4*B0**2 + 1` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
logarithmic term `log(mu0)/(20*pi**2)` is not a number and does not
assign `c_theta(mu0)`. Its derivative in `mu0` is `1/(20*pi**2*mu0)`.
The endpoint keeps `c_theta(mu0)`. `sigma` stays `2*c_theta*Lhat^{uu}`.
The class is an unavailable physical input. `c_theta` is not set to zero
and is not set equal to `b_theta`. The coefficient derivative is not
recomputed as a new result.

The ratio of the two bracket derivatives is `-2*log(f)`. It is not the
recorded mismatch `1`. The class is an unresolved proof. The ratio is
not adopted as a two-collar theorem. `delta` is not set to zero.
`d_theta` is not set to zero. Neither bracket derivative is recomputed
as a new result.

The `a` derivative of the width derivative is `-2/h**2`. It is not a
number. It keeps `h`. Its derivative in `D_j` is `0`. The clock
derivative in `D_j` is `0`. `numeric_D_j` is null. The class is an
unavailable physical input. The derivative is not read as `D_j`. The
width derivative is not recomputed as a new result.

The `B0` derivative of the substituted susceptibility is `8*B0`. It
does not measure `cJ`. Its derivative in `cJ` is `0`. The class is an
unavailable physical input. `cJ` is not set to zero. `B0` is not set to
one. The substituted susceptibility is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
