# Source coefficient, mixed trace, bound constant, and the second slope

S094, E003, 2026-10-14. This note takes the four inputs named by S093.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S093 scheme coefficient is not recomputed as
a new result. The S093 product trace is not recomputed as a new result.
The S093 mouth derivative is not recomputed as a new result. The S093
susceptibility slope is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written
before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, `d_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Finite coefficient of `integral sqrt(gamma) u^2`

Assumption. A number for `c_theta(mu0)` would be the finite coefficient
of `integral sqrt(gamma) u^2`.

Error comparison. The curved-collar note names that dimension-one
source and does not assign it a formula. It says the finite source
coefficient, the finite `J_theta` choice, and the finite `c_theta`
must all remain explicit. That source coefficient is not `d_theta`.
`c_theta` is not set to zero.

Comparison. Require the sentence that the a3 subtraction renormalizes a
dimension-one coefficient of `integral sqrt(gamma) u^2`. Search that
note for an equality that sets the coefficient to a number. Compare the
endpoint `c_theta(mu0) + log(mu0)/(20*pi**2)` with the absence of that
number. Differentiate the endpoint in the integration constant and check
the derivative by a finite difference.

Hardest step. The source term can be named while its finite coefficient
stays unassigned. Naming the term is not a value of `c_theta(mu0)`.

Observable. Whether the curved-collar note assigns a number to that
coefficient, and whether the endpoint is a number.

Test.

1. Require the dimension-one sentence and the sentence that the three
   finite choices remain explicit.
2. Require the note to contain no equality of the form
   `source coefficient =` followed by a number. Require the endpoint to
   keep `c_theta(mu0)`, and require a finite difference in that constant
   to agree with its derivative.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta` and `mu`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Product trace `tr(SE)`

Assumption. The unequal-collar invariant would be the product trace
`tr(SE)`.

Error comparison. The coupled note says `tr(SE)` vanishes on the
product, together with the other product traces, and unequal collars
are not covered by equation (8). The product trace is not the collar
mismatch. `delta` is not set to zero.

Comparison. Build `S` from the shipped collar matrices and let `E` be a
general angular endomorphism of the same kind the coupled checker uses.
Evaluate `tr(S E)` and check it at one numerical angle. Compare that
trace with the recorded mismatch `1`.

Hardest step. The mixed product trace can vanish for a general angular
`E` while the unequal-collar mismatch stays outside the theorem. A
vanishing product trace is not the missing invariant.

Observable. The product trace `tr(S E)`, and whether it equals the
mismatch.

Test.

1. Require the sentence that `tr(SE)` vanishes on the product and the
   exclusion of unequal collars.
2. Record `tr(S E)`. Require a numerical evaluation at one angle to
   agree. Require the trace to differ from `1`.
3. Do not adopt the vanishing trace as a two-collar theorem.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Bound constant `C`

Assumption. A clock number for one `D_j` would be the bound constant
`C` in the mouth estimate that multiplies `Lbar`.

Error comparison. The mouth note says that constant is independent of
the individual `D_j`. Independence is not an assignment. `Lbar` remains
a definition. The determinant is not set to one. `delta` is not set to
zero.

Comparison. Require estimate (2) and the independence sentence.
Differentiate the displayed bracket in the average length inside
`Lbar`, and differentiate the clock expression for `Lambda` in `D_j`.
Require `C` to remain a symbol.

Hardest step. The bracket can depend on `Lbar` while `C` stays an
unnamed factor and the clock derivative in `D_j` stays zero. Neither
fact assigns `D_j`.

Observable. Whether `C` is a number, and the clock derivative in `D_j`.

Test.

1. Require the displayed bound and the sentence that the constant is
   independent of individual `D_j`.
2. Require the derivative of the bracket in the average length to keep
   the factor `a`, and require a finite difference to agree. Require
   the clock derivative in `D_j` to vanish. Require `C` not to be a
   number.
3. Require `D_j`, `Lbar`, and `C` to be absent from the march
   parameters. Do not assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Second `B0` derivative of the parallel susceptibility

Assumption. A measured `cJ` would be the second `B0` derivative of the
parallel susceptibility at gravitational `kappa = 1`.

Error comparison. That second derivative can stay linear in `cJ`. The
first `B0` slope is a different quantity. The conditional axial window
is not adopted. `cJ` is not set to zero. A positive quadratic term is
not used as a stability mass.

Comparison. Call `invariants()`. Differentiate the shipped parallel
susceptibility twice in `B0` at gravitational `kappa = 1`. Check that
second derivative by a second finite difference.

Hardest step. The second slope can be a constant multiple of `cJ`
while no measurement fixes `cJ`. The coefficient is not a mass.

Observable. The second `B0` derivative at gravitational `kappa = 1`,
and whether it is a number.

Test.

1. Require the returned susceptibility and the sentence that magnetic
   background matching alone cannot measure `cJ`.
2. Record the second derivative at gravitational `kappa = 1`. Require a
   second finite difference to agree. Require `cJ` to remain.
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

Filled from the receipt after the two agreeing checker runs. The curved-collar
note assigns no number to the finite coefficient of `integral sqrt(gamma) u^2`.
The endpoint is `c_theta(mu0) + log(mu0)/(20*pi**2)`. Its derivative in the
integration constant is `1`. The endpoint is not a number. `sigma` stays
`2*c_theta*Lhat^{uu}`. The class is an unavailable physical input. `c_theta`
is not set to zero and is not set equal to `d_theta`. The scheme coefficient
is not recomputed as a new result.

The product trace `tr(SE)` is `0`. It is not the recorded collar mismatch `1`.
The class is an unresolved proof. That trace is not adopted as a two-collar
theorem. `delta` is not set to zero. The product trace `tr(S^3)` is not
recomputed as a new result.

The derivative of the displayed bracket in the average length is `a`. The
bound constant `C` is not a number. The clock derivative in `D_j` is `0`.
`numeric_D_j` is null. The class is an unavailable physical input. The
determinant is not set to one. The `Lbar` derivative is not recomputed as a
new result.

The second `B0` derivative of the parallel susceptibility is
`4*cJ*kappa**4`. At gravitational `kappa = 1` it is `4*cJ`. It is not a
number. The march parameter `kappa` is not used as the gravitational
`kappa`. The class is an unavailable physical input. `cJ` is not set to
zero. The first `B0` slope is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
