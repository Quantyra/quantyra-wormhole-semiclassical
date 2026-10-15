# J_theta choice, mixed curvature trace, local bound, and the cJ slope

S095, E003, 2026-10-14. This note takes the four inputs named by S094.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S094 source coefficient is not recomputed as
a new result. The S094 product trace is not recomputed as a new result.
The S094 bracket derivative is not recomputed as a new result. The S094
second `B0` derivative is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written
before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Finite `J_theta` choice

Assumption. A number for `c_theta(mu0)` would be the finite `J_theta`
choice.

Error comparison. The curved-collar note defines `J_theta` as a
functional and says the coefficient of `-d_theta J_theta/3` is scheme
dependent. `J_theta` is not that source coefficient. `c_theta` is not
set to zero and is not set equal to `d_theta`.

Comparison. Require the definition of `J_theta` and the sentence that
its coefficient is scheme dependent. Search that note for an equality
that sets `J_theta` to a number. Differentiate `-d_theta J_theta/3` in
`J_theta` and check the derivative by a finite difference.

Hardest step. The variation can name `J_theta` while the finite choice
stays unassigned. Naming the functional is not a value of
`c_theta(mu0)`.

Observable. Whether `J_theta` is a number, and whether the derivative
in `J_theta` assigns `c_theta(mu0)`.

Test.

1. Require the definition of `J_theta` and the sentence that its
   coefficient is scheme dependent.
2. Require the note to contain no equality of the form `J_theta =`
   followed by a number. Require a finite difference in `J_theta` to
   agree with the derivative of `-d_theta J_theta/3`. Require the
   endpoint to keep `c_theta(mu0)`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `J_theta` and `c_theta`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Product trace `tr((chi')^2 S)`

Assumption. The unequal-collar invariant would be the product trace
`tr((chi')^2 S)`.

Error comparison. The coupled note says that trace vanishes on the
product, together with the other product traces, and unequal collars
are not covered by equation (8). That trace is not `tr(SE)`. `delta`
is not set to zero.

Comparison. Build `S` and `chi'` from the shipped collar matrices.
Evaluate `tr((chi')^2 S)` and check it at one numerical angle. Compare
that trace with the recorded mismatch `1`.

Hardest step. This product trace can vanish while the unequal-collar
mismatch stays outside the theorem. A vanishing product trace is not
the missing invariant.

Observable. The product trace `tr((chi')^2 S)`, and whether it equals
the mismatch.

Test.

1. Require the sentence that `tr((chi')^2 S)` vanishes on the product
   and the exclusion of unequal collars.
2. Record `tr((chi')^2 S)`. Require a numerical evaluation at one
   angle to agree. Require the trace to differ from `1`.
3. Do not adopt the vanishing trace as a two-collar theorem.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Local J derivative bound

Assumption. A number for one `D_j` would be the local J derivative
bound `O(a+a epsilon[1+|log(lambda_0/epsilon)|])`.

Error comparison. The mouth note displays that bound as a big-O
statement. A big-O statement is not an assignment. The expression
inside the bound is not the bracket derivative. `delta` is not set to
zero. The determinant is not set to one.

Comparison. Require the displayed bound. Differentiate the expression
inside it in `lambda_0`, and check that derivative by the limit of a
central quotient and by a high-precision numerical difference.
Differentiate the clock expression for `Lambda` in `D_j`.

Hardest step. The expression inside the bound can be differentiated
while every symbol remains and the clock derivative in `D_j` stays
zero. Neither fact assigns `D_j`.

Observable. Whether the bound is a number, and the clock derivative in
`D_j`.

Test.

1. Require the displayed big-O bound.
2. Require the expression inside the bound not to be a number. Require
   the limit of the central quotient in `lambda_0` to agree with the
   direct derivative, and require a numerical difference to agree.
   Require the clock derivative in `D_j` to vanish.
3. Require `D_j`, `a`, `epsilon`, and `lambda_0` to be absent from the
   march parameters. Do not evaluate the big-O statement as an
   equality. Do not assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## `cJ` derivative of the second `B0` derivative

Assumption. A measured `cJ` would be the `cJ` derivative of the second
`B0` derivative of the parallel susceptibility at gravitational
`kappa = 1`.

Error comparison. Removing `cJ` by differentiation does not measure
`cJ`. The second `B0` derivative is a different quantity. The
conditional axial window is not adopted. `cJ` is not set to zero. A
positive quadratic term is not used as a stability mass.

Comparison. Call `invariants()`. Differentiate the shipped parallel
susceptibility twice in `B0` and then once in `cJ`, at gravitational
`kappa = 1`. Check that `cJ` derivative by a finite difference.

Hardest step. The derivative in `cJ` can be independent of `cJ` while
no measurement fixes `cJ`. The resulting number is not a mass and is
not a value of `cJ`.

Observable. The `cJ` derivative at gravitational `kappa = 1`, and
whether it is a value of `cJ`.

Test.

1. Require the returned susceptibility and the sentence that magnetic
   background matching alone cannot measure `cJ`.
2. Record the `cJ` derivative at gravitational `kappa = 1`. Require a
   finite difference in `cJ` to agree. Require that derivative not to
   be recorded as a value of `cJ`.
3. Require the march parameters to omit `cJ`, `chi`, and `B0`. Do not
   use the march's `kappa` as the gravitational `kappa`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The derivative
of `-d_theta J_theta/3` in `J_theta` is `-d_theta/3`. It is not a number
and does not assign `c_theta(mu0)`. The endpoint keeps `c_theta(mu0)`.
`sigma` stays `2*c_theta*Lhat^{uu}`. The class is an unavailable physical
input. `c_theta` is not set to zero. The source coefficient is not
recomputed as a new result.

The product trace `tr((chi')^2 S)` is `0`. It is not the recorded collar
mismatch `1`. The class is an unresolved proof. That trace is not adopted
as a two-collar theorem. `delta` is not set to zero. The product trace
`tr(SE)` is not recomputed as a new result.

The expression inside the local J derivative bound is not a number. Its
derivative in `lambda_0` is
`a*epsilon*sign(-log(epsilon) + log(lambda_0))/lambda_0`. The clock derivative in
`D_j` is `0`. `numeric_D_j` is null. The class is an unavailable physical
input. The big-O statement is not evaluated as an equality. The bracket
derivative is not recomputed as a new result.

The `cJ` derivative of the second `B0` derivative at gravitational
`kappa = 1` is `4`. That number is not a value of `cJ`. The class is an
unavailable physical input. `cJ` is not set to zero. The second `B0`
derivative is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
