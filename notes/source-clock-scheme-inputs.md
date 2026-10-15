# Scheme coefficient, product trace, mouth length, and the susceptibility slope

S093, E003, 2026-10-14. This note takes the four inputs named by S092.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S092 increment is not recomputed as a new
result. The S092 collar mismatch is not recomputed as a new result.
The S092 trace gap is not recomputed as a new result. The S092 second
parallel derivative is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written
before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Scheme coefficient `d_theta`

Assumption. A number for the endpoint value of `c_theta` would be the
shipped scheme coefficient `d_theta`.

Error comparison. `d_theta` multiplies the normal term, not `I_theta`.
The curved-collar note says the finite `J_theta` choice and the finite
`c_theta` must both remain explicit. The endpoint value is not replaced
by `d_theta`. `c_theta` is not set to zero.

Comparison. Read `d_theta=4/[15(4pi)^2]`. Form that number and subtract
it from the endpoint `c_theta(mu0) + log(mu0)/(20*pi**2)`. Differentiate
the difference in the integration constant and check that derivative by
a finite difference. Compare `d_theta` with the march parameters and
with the shipped sigma string.

Hardest step. `d_theta` can be a pure number while the endpoint still
contains `c_theta(mu0)` and `mu0`. The pure number is not a
renormalization condition for `c_theta`.

Observable. Whether `d_theta` equals the endpoint value.

Test.

1. Require the shipped formulas for `b_theta` and `d_theta`, and the
   sentence that both finite choices remain explicit.
2. Record `d_theta` and the difference from the endpoint. Require the
   difference to keep the integration constant, and require a finite
   difference in that constant to agree with its derivative.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta`, `mu`, and `d_theta`.
   Do not set `c_theta` equal to `d_theta`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Product trace `tr S^3`

Assumption. The unequal-collar invariant would be the product trace
`tr S^3`.

Error comparison. The coupled note says `tr S=tr S^3=0` on the product,
and unequal collars are not covered by equation (8). The product trace
is not the collar mismatch. `delta` is not set to zero.

Comparison. Build `S` from the shipped collar matrices. Evaluate
`tr(S^3)` and check it at one numerical angle. Compare that trace with
the recorded mismatch `1`. Do not introduce a two-collar equation.

Hardest step. The product trace can vanish while the unequal-collar
mismatch stays outside the theorem. A vanishing product trace is not
the missing invariant.

Observable. The product trace, and whether it equals the mismatch.

Test.

1. Require `tr S=tr S^3=0` and the exclusion of unequal collars.
2. Record `tr(S^3)`. Require a numerical evaluation at one angle to
   agree. Require the trace to differ from `1`.
3. Do not adopt the vanishing trace as a two-collar theorem.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Identity `Lbar = T + average_j D_j`

Assumption. An equality setting one `D_j` would solve
`Lbar = T + average_j D_j`.

Error comparison. `Lbar` does not assign that length. The mouth note
states the identity as a definition, and `average D_j` remains
`O(delta log q)`. That order statement is not an equality. The
determinant is not set to one. `delta` is not set to zero.

Comparison. Require the definition `Let Lbar=T+average_j D_j`.
Differentiate `Lbar` in the average length and differentiate the clock
expression for `Lambda` in `D_j`. Evaluate the definition at the
predeclared point `T = 1`, `average D_j = 1`.

Hardest step. The definition can hold at a chosen point while the clock
derivative in `D_j` stays zero. Neither fact assigns `D_j`.

Observable. The derivative of `Lbar` in the average length, and the
clock derivative in `D_j`.

Test.

1. Require the definition and the order statement `average D_j=O(delta log q)`.
2. Require the derivative of `Lbar` in the average length to be `1`,
   and require a finite difference to agree. Require the clock
   derivative in `D_j` to vanish.
3. Require `D_j` and `Lbar` to be absent from the march parameters.
   Do not assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## `B0` derivative of the parallel susceptibility

Assumption. A measured `cJ` would be the `B0` derivative of the parallel
susceptibility at `B0 = 1` and gravitational `kappa = 1`.

Error comparison. That slope is not a number. It stays linear in `cJ`.
The second parallel derivative is a different quantity. The conditional
axial window is not adopted. `cJ` is not set to zero. A positive
quadratic term is not used as a stability mass.

Comparison. Call `invariants()`. Differentiate the shipped parallel
susceptibility in `B0` at the unit background. Check that derivative by
a central finite difference. Compare it with a pure number.

Hardest step. The slope in `B0` can be `4*cJ` while no measurement
fixes `cJ`. The coefficient is not a mass.

Observable. The `B0` derivative at the unit background, and whether it
is a number.

Test.

1. Require the returned susceptibility and the sentence that magnetic
   background matching alone cannot measure `cJ`.
2. Record the derivative at `B0 = 1`, `kappa = 1`. Require a central
   finite difference to agree. Require `cJ` to remain.
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

Filled from the receipt after the two agreeing checker runs. The shipped
scheme coefficient is `d_theta = 1/(60*pi**2)`. The endpoint minus that
number is `(60*pi**2*c_theta(mu0) + 3*log(mu0) - 1)/(60*pi**2)`. That
difference is not a number. `sigma` stays `2*c_theta*Lhat^{uu}`. The
class is an unavailable physical input. `c_theta` is not set equal to
`d_theta`. The scale increment is not recomputed as a new result.

The product trace `tr(S^3)` is `0`. It is not the recorded collar
mismatch `1`. The class is an unresolved proof. That trace is not
adopted as a two-collar theorem. `delta` is not set to zero. The collar
mismatch is not recomputed as a new result.

The derivative of `Lbar` in the average length is `1`. The clock
derivative in `D_j` is `0`. `numeric_D_j` is null. The class is an
unavailable physical input. The determinant is not set to one. The
trace gap is not recomputed as a new result.

At `B0 = 1` and gravitational `kappa = 1` the `B0` derivative of the
parallel susceptibility is `4*cJ`. It is not a number. The march
parameter `kappa` is not used as the gravitational `kappa`. The class
is an unavailable physical input. `cJ` is not set to zero. The second
parallel derivative is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S088 relative series is not differentiated
again. No new `Lambda` interval is recorded. The area budget is not
evaluated. No successor model is opened.
