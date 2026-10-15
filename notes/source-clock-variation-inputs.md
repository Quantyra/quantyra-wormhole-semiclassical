# Weyl coefficient, traced term, little-o bound, and the susceptibility ratio

S096, E003, 2026-10-14. This note takes the four inputs named by S095.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S095 `J_theta` derivative is not recomputed as
a new result. The S095 curvature trace is not recomputed as a new result.
The S095 local-bound derivative is not recomputed as a new result. The
S095 `cJ` derivative is not recomputed as a new result. No new `Lambda`
interval is sought. No replacement area target is written before the
test. The area budget stays unevaluated. `zeta_s`, `sigma`, `cJ`,
`c_theta`, `d_theta`, and `delta` are not set to zero. An endpoint
contact energy is not cancelled. A positive quadratic term is not used
as a stability mass.

## Coefficient `-3` in the Weyl variation of `J_theta`

Assumption. A number for `c_theta(mu0)` would be the coefficient `-3`
in the Weyl variation of `J_theta`.

Error comparison. That coefficient multiplies the Weyl density. It is
the change of the functional, not a value of the finite coupling. That
coefficient is not a value of `c_theta(mu0)`. `c_theta` is not set to
zero and is not set equal to `-3`.

Comparison. Require the sentence that `J_theta` changes by `-3` times
the Weyl density. Differentiate that variation in the density and check
the derivative by a finite difference. Compare the coefficient with the
endpoint that still contains `c_theta(mu0)`.

Hardest step. The variation can display the number `-3` while the
endpoint stays unassigned. A counting factor is not a value of
`c_theta(mu0)`.

Observable. The coefficient `-3`, and whether the endpoint equals it.

Test.

1. Require the sentence that `J_theta` changes by `-3` times the Weyl
   density.
2. Record the derivative of that variation in the density. Require a
   finite difference to agree. Require the endpoint to keep
   `c_theta(mu0)` and to differ from `-3`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`. Do not set
   `c_theta` equal to `-3`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Traced `S_:aa` term

Assumption. The unequal-collar invariant would be the traced `S_:aa`
term.

Error comparison. The curved-collar note says `tr S=-2L` is independent
of `theta`, so the traced term adds none. Independence of `theta` is
not a two-collar theorem. Unequal collars are not covered by equation
(8). This term is not `tr((chi')^2 S)`. `delta` is not set to zero.

Comparison. Differentiate the displayed trace `-2L` in `theta`. Build
`S` from the shipped collar matrices, differentiate its trace in
`theta`, and check that derivative by a finite difference. Compare both
derivatives with the recorded mismatch `1`.

Hardest step. Both derivatives can vanish on the product while the
unequal-collar mismatch stays outside the theorem. A vanishing
derivative is not the missing invariant.

Observable. The `theta` derivative of the traced term, and whether it
equals the mismatch.

Test.

1. Require the sentence that the traced `S_:aa` term adds none, and the
   exclusion of unequal collars.
2. Record both `theta` derivatives. Require a numerical evaluation of
   the shipped derivative at one angle to agree. Require each derivative
   to differ from `1`.
3. Do not adopt either derivative as a two-collar theorem. Do not
   identify the shipped trace with `-2L`.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Little-o bound on the derivatives of `E_cov`

Assumption. A number for one `D_j` would be the statement that both
derivatives of `E_cov` are `O(1)=o(q delta)`.

Error comparison. The mouth note writes a little-o comparison. A
little-o comparison is not an assignment. The product `q delta` is not
the local J derivative bound. `delta` is not set to zero. The
determinant is not set to one.

Comparison. Require the displayed statement. Differentiate the product
`q delta` in `q`, with `delta` kept as a symbol, and check that
derivative by a finite difference. Differentiate the clock expression
for `Lambda` in `D_j`.

Hardest step. The product can be differentiated while `delta` remains
and the clock derivative in `D_j` stays zero. Neither fact assigns
`D_j`.

Observable. Whether `o(q delta)` is a number, and the clock derivative
in `D_j`.

Test.

1. Require the sentence that both derivatives of `E_cov` are
   `O(1)=o(q delta)`.
2. Require the derivative of `q delta` in `q` to keep `delta`, and
   require a finite difference to agree. Require the clock derivative
   in `D_j` to vanish.
3. Require `D_j`, `q`, and `delta` to be absent from the march
   parameters. Do not evaluate the little-o statement as an equality.
   Do not assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Ratio of the second `B0` derivative to its `cJ` derivative

Assumption. A measured `cJ` would be the ratio of the second `B0`
derivative to its `cJ` derivative at gravitational `kappa = 1`.

Error comparison. Returning the symbol `cJ` does not measure `cJ`. The
`cJ` derivative is a different quantity. The conditional axial window
is not adopted. `cJ` is not set to zero. A positive quadratic term is
not used as a stability mass.

Comparison. Call `invariants()`. Form the second `B0` derivative of the
shipped parallel susceptibility at gravitational `kappa = 1`, divide by
its derivative in `cJ`, and check the quotient at one nonzero value of
`cJ`.

Hardest step. The quotient can be the symbol `cJ` while no measurement
fixes that symbol. The quotient is not a mass.

Observable. The ratio at gravitational `kappa = 1`, and whether it is
a number.

Test.

1. Require the returned susceptibility and the sentence that magnetic
   background matching alone cannot measure `cJ`.
2. Record the ratio at gravitational `kappa = 1`. Require a numerical
   evaluation at one nonzero `cJ` to agree. Require `cJ` to remain.
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
of the Weyl variation in the density is `-3`. It is not a value of
`c_theta(mu0)`. The endpoint keeps `c_theta(mu0)` and does not equal `-3`.
`sigma` stays `2*c_theta*Lhat^{uu}`. The class is an unavailable physical
input. `c_theta` is not set to zero and is not set equal to `-3`. The
`J_theta` derivative is not recomputed as a new result.

The `theta` derivative of the displayed trace `-2L` is `0`. The `theta`
derivative of the shipped trace is `0`. Neither equals the recorded
mismatch `1`. The shipped trace is not identified with `-2L`. The class
is an unresolved proof. Neither derivative is adopted as a two-collar
theorem. `delta` is not set to zero. The curvature trace `tr((chi')^2 S)`
is not recomputed as a new result.

The statement `O(1)=o(q delta)` is not a number. The derivative of
`q delta` in `q` is the symbol `delta`. The clock derivative in `D_j`
is `0`. `numeric_D_j` is null. The class is an unavailable physical
input. The little-o statement is not evaluated as an equality. The
local-bound derivative is not recomputed as a new result.

The ratio of the second `B0` derivative to its `cJ` derivative at
gravitational `kappa = 1` is the symbol `cJ`. That symbol is not a
number and does not measure `cJ`. The class is an unavailable physical
input. `cJ` is not set to zero. The `cJ` derivative is not recomputed
as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
