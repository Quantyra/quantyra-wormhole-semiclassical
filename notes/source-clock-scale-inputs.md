# Scale `mu`, first bracket, width derivative, and the susceptibility `1+4*B0**2`

S104, E003, 2026-10-15. This note takes the four inputs named by S103.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S103 scale derivative is not recomputed as a
new result. The S103 bracket derivative is not recomputed as a new
result. The S103 ceiling is not recomputed as a new result. The S103
ratio is not recomputed as a new result. No new `Lambda` interval is
sought. No replacement area target is written before the test. The area
budget stays unevaluated. `zeta_s`, `sigma`, `cJ`, `c_theta`, `d_theta`,
and `delta` are not set to zero. An endpoint contact energy is not
cancelled. A positive quadratic term is not used as a stability mass.

## Scale in the dimension-one monomial

Assumption. A number for `c_theta(mu0)` would be the scale `mu` in
`c_one*mu`.

Error comparison. That scale is not a value of `c_theta(mu0)`. The
derivative of the monomial in `c_one` is the symbol `mu`. That symbol is
not a number. `c_theta` is not set to zero and is not set equal to
`b_theta`. The scale derivative that returns `c_one` is not recomputed.

Comparison. Require the dimension-one clause. Form `c_one*mu`.
Differentiate it in `c_one` and check the derivative by a finite
difference. Compare the result with the endpoint that still contains
`c_theta(mu0)`.

Hardest step. The scale can be the object returned by a coefficient
derivative while the endpoint coefficient stays independent. Returning
`mu` does not assign `c_theta(mu0)`.

Observable. Whether the scale `mu` is a number, and the derivative of
`c_one*mu` in `c_one`.

Test.

1. Require the sentence `it renormalizes a dimension-one coefficient`.
2. Record the derivative of `c_one*mu` in `c_one`. Require a finite
   difference to agree. Require the derivative to be the symbol `mu` and
   not a number. Require the endpoint to keep `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set `c_theta`
   equal to `b_theta`. Do not store the derivative in `mu` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## First bracket in (4)

Assumption. The unequal-collar invariant would be the first-bracket term
`b_theta (log f/2)(2/3)(R^2/f)w`.

Error comparison. That bracket term is not a two-collar theorem. The
local relative limit is not a proof for unequal driven collars. The
derivative in `w` is `R**2*log(f)/(60*pi**2*f)`. That derivative is not
the recorded collar mismatch and is not the second-bracket derivative
`-R**2/(120*pi**2*f)`. `delta` is not set to zero. `d_theta` is not set
to zero.

Comparison. Require the displayed first bracket, the shipped coefficient
`b_theta=4/[5(4pi)^2]`, and the exclusion of unequal collars.
Differentiate `b_theta*(log(f)/2)*(2/3)*(R**2/f)*w` in `w` and check the
derivative by a finite difference. Compare the derivative with the
recorded mismatch `1`.

Hardest step. The first bracket can stay dependent on one lapse and one
weight while the two collars remain uncovered. A one-boundary quotient
is not the missing invariant.

Observable. The derivative of the first bracket in `w`, and whether it
equals the mismatch.

Test.

1. Require `b_theta (log f/2)(2/3)(R^2/f)w`, the shipped `b_theta`
   formula, and the statement that this limit is not a proof for unequal
   driven collars.
2. Record the derivative in `w`. Require a finite difference to agree.
   Require the derivative to differ from `1` and from
   `-R**2/(120*pi**2*f)`.
3. Do not adopt the derivative as a two-collar theorem. Do not
   differentiate the RN relative series. Do not store the `f'` derivative
   of the second bracket as this result. Do not set `d_theta` to zero.

A pass is a scoped obstruction. The class is an unresolved proof. No
admissible two-collar coefficient is recorded.

## Width derivative of the ceiling

Assumption. A number for one `D_j` would be the width derivative
`-2*a/h**2`.

Error comparison. That derivative is not an assignment of `D_j`. It keeps
`a` and `h`. It is not the ceiling `2*a/h`. `a` and `h` are not set to
zero. `delta` is not set to zero. The determinant is not set to one.

Comparison. Require the displayed ceiling. Differentiate `2*a/h` in `h`
and check that derivative by the limit of a central quotient and by a
finite difference. Differentiate the resulting expression in `D_j`, and
differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The width derivative can be an exact rational expression
while `D_j` is absent. A quotient that still contains `a` and `h` does
not assign `D_j`.

Observable. Whether `-2*a/h**2` is a number for `D_j`, and the clock
derivative in `D_j`.

Test.

1. Require `eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h.`
2. Require the `h` derivative of the ceiling to equal `-2*a/h**2`.
   Require the limit of a central quotient and a finite difference to
   agree. Require the derivative in `D_j` to vanish. Require the clock
   derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `h` to be absent from the march
   parameters. Do not read the width derivative as `D_j`. Do not store
   the ceiling `2*a/h` as this result. Do not set `a` or `h` to zero.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Substituted susceptibility

Assumption. A measured `cJ` would be the susceptibility `1+4*B0**2`.

Error comparison. The expression `1+4*B0**2` does not measure `cJ`.
Writing the integer `2` in place of `cJ` removes `cJ` from the reduced
susceptibility. The conditional axial window is not adopted. `cJ` is not
set to zero. `B0` is not set to one. A positive quadratic term is not
used as a stability mass.

Comparison. Read the shipped parallel susceptibility. Reduce it at
gravitational `kappa = 1`. Substitute the integer `2` for `cJ` and
compare the result with `1+4*B0**2`. Differentiate the substituted
expression in `cJ` and check that derivative by a finite difference.

Hardest step. The substitution can produce a quadratic in `B0` while
`cJ` never remains. That quadratic is not a measured coupling, and
setting `B0` to one is not this test.

Observable. The substituted susceptibility, and whether that expression
is `cJ`.

Test.

1. Require the sentence that magnetic background matching alone cannot
   measure `cJ`.
2. Record the susceptibility after the integer `2` replaces `cJ`.
   Require it to equal `1+4*B0**2`. Require a numerical evaluation at
   one nonzero `B0` other than one to agree. Require the `cJ` derivative
   of the substituted expression to vanish, and require a finite
   difference to agree. Do not set `B0` to one.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`. Do not store the simplified ratio `2` as this result.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
coefficient derivative of `c_one*mu` is `mu`. That symbol is not a number
and does not assign `c_theta(mu0)`. The endpoint keeps `c_theta(mu0)`.
`sigma` stays `2*c_theta*Lhat^{uu}`. The class is an unavailable physical
input. `c_theta` is not set to zero and is not set equal to `b_theta`.
The scale derivative is not recomputed as a new result.

The derivative of the first bracket in `w` is
`R**2*log(f)/(60*pi**2*f)`. It is not the recorded mismatch `1` and it is
not the second-bracket derivative `-R**2/(120*pi**2*f)`. The class is an
unresolved proof. The derivative is not adopted as a two-collar theorem.
`delta` is not set to zero. `d_theta` is not set to zero. The second-bracket
derivative is not recomputed as a new result.

The width derivative of the ceiling is `-2*a/h**2`. It is not a number.
It keeps `a` and `h`. Its derivative in `D_j` is `0`. The clock derivative
in `D_j` is `0`. `numeric_D_j` is null. The class is an unavailable
physical input. The derivative is not read as `D_j`. The ceiling is not
recomputed as a new result.

The substituted susceptibility is `4*B0**2 + 1`. It does not measure
`cJ`. Its derivative in `cJ` is `0`. The class is an unavailable physical
input. `cJ` is not set to zero. `B0` is not set to one. The ratio is not
recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
