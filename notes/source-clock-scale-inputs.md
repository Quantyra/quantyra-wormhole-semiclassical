# Renormalized coupling, substitution error, channel data, and a measured pair

S090, E003, 2026-10-14. This note takes the four inputs named by S089.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S087 product-throat derivative is
not recomputed. The S088 relative series is not differentiated again.
The S089 `L00` derivative is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written
before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, and `delta` are not set to zero. An endpoint contact
energy is not cancelled. A positive quadratic term is not used as a
stability mass.

## Renormalized value of `c_theta`

Assumption. A renormalized value is a number for `c_theta` at one
stated scale. The candidate coefficient is the shipped factor in the
mixed a4 curvature term, divided by `(4pi)^2`, which the curved-collar
note calls `b_theta`. The finite coupling in `sigma = 2 c_theta Lhat^{uu}`
is a different object from that beta-function coefficient.

Error comparison. The `L00` derivative recorded by S089 is not a value
of `c_theta`, and it is not differentiated again. A beta-function
coefficient does not assign the finite coupling at a scale. `sigma` is
not replaced by `b_theta`.

Comparison. Call `clifford_checks()` and `charge_identities()`. Divide
the shipped a4 by the traceless quadratic form and by `(4pi)^2`.
Differentiate `2*c_theta*Lhat` in `c_theta` and check that derivative
by a finite difference. Compare the symbols with the march parameters.

Hardest step. The curvature prefactor can be a pure number while
`c_theta` remains a symbol in `sigma`. The pure number is not a
renormalized coupling at a stated scale.

Observable. The a4 prefactor, `b_theta`, and the `c_theta` derivative
of `sigma`.

Test.

1. Require the returned a4 of both domains to be a constant times
   `(u·L·u - tr(L)|u|^2/3)`. Record that constant.
2. Set `b_theta` to that constant over `(4pi)^2`. Require it to be a
   number that does not contain `c_theta`.
3. Require the shipped sigma string to be `2*c_theta*Lhat^{uu}`.
   Differentiate `2*c_theta*Lhat` in `c_theta` and require a finite
   difference to agree. Require the derivative to contain `Lhat`.
4. Require the march parameters to omit `c_theta` and `mu`. Do not
   march. Do not set `c_theta` to zero.

If `b_theta` is a number and `c_theta` remains in `sigma`, the
renormalized value is an unavailable physical input.

## Error of the two-rate substitution

Assumption. An error bound would compare the one-angle normal-smearing
coefficient with a proved two-collar coefficient. The candidate is the
difference between the shipped one-angle square and the same square
after the two rates are substituted. Equation (8) has one angle.

Error comparison. The product coefficient recorded by S089 is the
algebra of that substitution, not a bound on a two-collar coefficient.
The shipped return does not contain a two-collar error. `delta` is not
set to zero. A positive square is not a stability mass.

Comparison. Call `analytic()`. Record its keys and the one-angle
string. Form the difference between that square and the substituted
square. Check the difference by a finite difference in one rate. Read
the exclusion sentence.

Hardest step. The difference can be an explicit expression while the
shipped return still has no two-collar coefficient to compare it with.
An expression for a substitution is not an error bound.

Observable. The shipped return keys, and the difference created by the
substitution.

Test.

1. Require `product_Delta_a4_unsmeared` to be 0 and the normal-smearing
   string to be one angle's square.
2. Form the difference between that square and the substituted square.
   Require a finite difference in one rate to agree with its derivative.
3. Require the return keys to contain no two-collar error bound.
4. Read the sentence that unequal collars are not covered by (8). Do
   not introduce a two-collar clock equation.

A pass is a scoped obstruction. The class is an unresolved proof.

## Declared channel count and one exterior length

Assumption. The clock datum would be one declared positive `q` together
with one numeric `D_j`. The scaling formulas are the three closed
expressions in the finite-mouth note. The mean return remains the
big-O sentence.

Error comparison. The `q` derivatives recorded by S089 are not values
of `D_j`. Evaluating the three closed formulas at a positive rational
`q` is an illustration. That rational is not assigned to the rule-D
march. The big-O is not changed by the illustration.

Comparison. Evaluate `a`, `epsilon`, and `lambda_0` at `q = 4`.
Differentiate each formula once more only to check the values by a
finite difference of the formulas, not to assign them. Read the big-O
sentence. Compare `q` and `D_j` with the march parameters and with
`Lambda`.

Hardest step. The three formulas become numbers at `q = 4`, and
`D_j` does not. A triple with a missing length is not the norm's data.

Observable. The three illustrated numbers, the big-O sentence, and the
march parameter list.

Test.

1. Require each illustrated value to agree with a finite difference of
   its formula between two neighbors of `q = 4`, in the sense that the
   centered quotient's limit is the derivative and the value at `q = 4`
   is the formula.
2. Require the note to state `average D_j=O(delta log q)`.
3. Require `q` and `D_j` to be absent from the march parameters and
   from `Lambda`. Do not call `march`. Do not assign `q = 4` to the
   clock.

A pass is a scoped obstruction. The class is an unavailable physical
input.

## Measured `chi` and gravitational `kappa`

Assumption. A measurement is a number for `chi` and a number for the
gravitational `kappa` at one nonzero numeric `B0`. The candidate is
the shipped parallel susceptibility with `B0 = 1` and `kappa = 1`
both substituted. Those substitutions are not measurements.

Error comparison. Leaving `kappa` symbolic was the S089 result, and
that expression is not the measurement. Substituting a number for
`kappa` by hand leaves `cJ`. The march parameter spelled `kappa` is
the contact coefficient, and it is not substituted for the
gravitational `kappa`. `cJ` is not set to zero.

Comparison. Call `invariants()`. Substitute `B0 = 1` and `kappa = 1`.
Differentiate in `cJ` and check the derivative by a finite difference.
Compare the remaining symbols with the march parameters.

Hardest step. Both backgrounds can be numbers while `chi` remains
linear in `cJ`. A linear expression is not a measured susceptibility.

Observable. The susceptibility after both substitutions, and the march
parameter list.

Test.

1. Require the twice-substituted susceptibility to contain `cJ` and no
   `B0` or `kappa`.
2. Require its `cJ` derivative to agree with a finite difference.
3. Require that derivative to be a nonzero number.
4. Require the march parameters to omit `chi`, `B0`, and `cJ`. Do not
   use the march's `kappa` as the gravitational `kappa`. Do not set
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
a4 prefactor is `4/5`. Then `b_theta = 1/(20*pi**2)`. The shipped sigma
string is `2*c_theta*Lhat^{uu}`, and its `c_theta` derivative is
`2*Lhat`. `sigma` is not replaced by `b_theta`. The class is an
unavailable physical input. `c_theta` is not set to zero.

The substitution difference is `-4*theta_R*(2*theta_L - theta_R)/15`.
The shipped return keys contain no two-collar error bound. The class is
an unresolved proof. `delta` is not set to zero.

At the illustrated `q = 4`, which is not assigned to the clock, the
scaling values are `a = 1/(16*pi**2)`, `epsilon = 1/16`, and
`lambda_0 = 1/(1024*pi)`. `average D_j` remains `O(delta log q)`.
`numeric_D_j` is null. The class is an unavailable physical input.
The determinant is not set to one.

At `B0 = 1` and `kappa = 1` the susceptibility is `2*cJ + 1`. Its
`cJ` derivative is `2`. `chi` is not a number. The march parameter
`kappa` is not used as the gravitational `kappa`. The class is an
unavailable physical input. `cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S089 `L00` derivative is not recomputed as
a new result. No new `Lambda` interval is recorded. The area budget
is not evaluated. No successor model is opened.
