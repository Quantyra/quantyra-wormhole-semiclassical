# Curved boundary, angle cross, mouth norm, and a nonzero background

S088, E003, 2026-10-14. This note takes the four inputs named by S087.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S087 product-throat derivative is
not recomputed. No new `Lambda` interval is sought. No replacement
area target is written before the test. The area budget stays
unevaluated. `zeta_s`, `sigma`, `cJ`, and `delta` are not set to
zero. An endpoint contact energy is not cancelled. A positive
quadratic term is not used as a stability mass.

## Extrinsic curvature inside the clock

Assumption. The candidate is the shipped RN relative type-B series,
not the product-throat weight. The boundary coordinate is `z = 1/R`.
A nonzero lapse derivative at the clock boundary would be a nonzero
first `z` derivative of that series at `z = 0`.

Error comparison. The product-throat weight already has limit 0, and
that calculation is not repeated. The type-B series is the curved
boundary retained by `weyl_and_asymptotic_checks()`. Its conformal
boundary is `R -> infinity`. A series that starts at `z^2` has a
vanishing value and a vanishing first derivative there. `sigma` is
not replaced by that zero.

Comparison. Call `weyl_and_asymptotic_checks()`. Differentiate the
returned relative series in `z`. Compare the analytic derivative with
a finite difference. Compare the boundary values with the shipped
`Lambda` string.

Hardest step. The first derivative at `z = 0` can vanish while the
second derivative stays nonzero. The second derivative is not a lapse
derivative on the clock boundary, and it is not a value of `sigma`.

Observable. The series, its value at `z = 0`, its first derivative
at `z = 0`, and its second derivative at `z = 0`.

Test.

1. Parse `RN_typeB_relative_series`. Require the value and the first
   derivative at `z = 0` to be 0. Require the second derivative at
   `z = 0` to be a nonzero explicit expression.
2. Check the first derivative by a finite difference at one positive
   `z`.
3. Require `Lambda` from `charge_identities()` to contain `zeta_s`
   and `sigma`. Do not substitute the boundary zero for either one.
4. Do not march `rad.march`.

If the first derivative vanishes and `sigma` remains in `Lambda`, a
nonzero clock-boundary lapse derivative is an unavailable physical
input.

## Folded paired heat coefficient

Assumption. The missing invariant is the cross of two angle gradients
`u_L` and `u_R`. Equation (8) has one angle. The candidate cross is
`u_L · u_R`, exposed by
`|u_L - u_R|^2 = |u_L|^2 + |u_R|^2 - 2 u_L · u_R`.

Error comparison. Unequal jets and unequal clocks lie outside the
matched-collar theorem and outside equation (8). The cross is not a
two-collar equation of motion. `delta` is not set to zero. A positive
square is not a stability mass.

Comparison. Expand the identity. A finite difference in one component
of `u_R` must agree with the analytic derivative of the cross. Read
the exclusion sentences.

Hardest step. The separate squares `|u_L|^2` and `|u_R|^2` can be
written down without the cross. The difference of the jets requires
the cross. Equation (8) does not contain it.

Observable. The cross `u_L · u_R` and the expanded difference identity.

Test.

1. Require the difference identity exactly.
2. Differentiate the cross in one component and check it by a finite
   difference.
3. Read the two exclusion sentences. Do not introduce a two-collar
   clock equation.

A pass is a scoped obstruction. The class is an unresolved proof.

## Norm on the wormhole `z` derivative

Assumption. The norm inputs are the symbols named by the finite-mouth
note for its nonlocal determinant: a positive `a`, `epsilon`, `D_j`,
and the exterior time normalization `lambda_0`. The S087 Schatten
value bound is not this norm and is not recomputed.

Error comparison. Those symbols are data of the mouth determinant.
The rule-D march does not take them. A stand-in radius is not
`lambda_0`. The determinant is not set to one, and the wormhole
determinant is not computed.

Comparison. Read the scope sentence from the finite-mouth note.
Compare its symbols with the parameter list of `march` and with
`Lambda`.

Hardest step. Naming `a`, `epsilon`, `D_j`, and `lambda_0` does not
assign their values. The clock polynomial is not given those inputs.

Observable. The march parameter list, and the absence of those four
symbols from it.

Test.

1. Require the march parameters to be
   `(eta, lam, kappa, z0, v0, method)`.
2. Require `a`, `epsilon`, `D_j`, and `lambda_0` to be absent from
   that list and from `Lambda`.
3. Read the mouth-note sentence that it must control the corrected
   background's nonlocal determinant. Do not call `march`.

A pass is a scoped obstruction. The class is an unavailable physical
input.

## Nonzero `B0`

Assumption. A nonzero background is the symbol `B0` kept free in the
shipped susceptibility. The four-null-photon coefficient returned by
`invariants()` is the sector that does not use that background. The
conditional axial window is not read.

Error comparison. A symbolic nonzero `B0` is not a measured
susceptibility. The photon coefficient is a pure number. The march
still does not take `chi` or `B0`. Setting `cJ` to zero is an extra
condition.

Comparison. Call `invariants()`. Record the photon coefficient.
Substitute a nonzero symbol for `B0` in the susceptibility and require
`cJ` to remain. Compare with the march parameters.

Hardest step. The solved formula still contains `chi` and `B0` when
`B0` is nonzero. The photon number does not replace either one.

Observable. The photon coefficient, the susceptibility at nonzero
`B0`, and the march parameter list.

Test.

1. Require the shipped four-null-photon coefficient to be `256` and
   to contain neither `cJ` nor `B0`.
2. Keep `B0` symbolic and nonzero. Require `cJ` to remain in the
   susceptibility.
3. Require the march parameters to omit `chi` and `B0`. Do not set
   `cJ` to zero. Do not adopt the conditional axial window.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The shipped
RN relative series is `-4*Q2*ell**2*z**3/3 + 2*ell**2*m*z**2`. Its value
and its first derivative at `z = 0` are 0. Its second derivative at
`z = 0` is `4*ell**2*m`. `sigma` is not replaced by that zero. The class
is an unavailable physical input.

The angle cross `u_L · u_R` satisfies the difference identity. The
folded coefficient is not proved. The class is an unresolved proof.
`delta` is not set to zero.

The mouth-note inputs `a`, `epsilon`, `D_j`, and `lambda_0` are absent
from the march parameters `eta, lam, kappa, z0, v0, method`. The class
is an unavailable physical input. The determinant is not set to one.

The four-null-photon coefficient is `256` and contains neither `cJ` nor
`B0`. `cJ` remains in the susceptibility at nonzero symbolic `B0`. The
class is an unavailable physical input. `cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S087 product-throat derivative is not
recomputed. No new `Lambda` interval is recorded. The area budget is
not evaluated. No successor model is opened.
