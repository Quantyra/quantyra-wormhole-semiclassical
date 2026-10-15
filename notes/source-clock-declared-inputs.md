# Curvature derivative, folded rate product, mouth scaling, and a numeric background

S089, E003, 2026-10-14. This note takes the four inputs named by S088.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S087 product-throat derivative is
not recomputed. The S088 relative series is not differentiated again.
No new `Lambda` interval is sought. No replacement area target is
written before the test. The area budget stays unevaluated.
`zeta_s`, `sigma`, `cJ`, and `delta` are not set to zero. An endpoint
contact energy is not cancelled. A positive quadratic term is not used
as a stability mass.

## Curvature contact and one boundary derivative

Assumption. The candidate is the shipped mixed a4 curvature coefficient
returned by `clifford_checks()`, for one matched collar. A first
boundary derivative of that contact is its derivative with respect to
one component of the second fundamental form.

Error comparison. The RN relative series starts at order `z^2`, and
that calculation is not repeated. On an umbilic boundary the
coefficient itself vanishes. A nonzero derivative in one curvature
component is the response to leaving that boundary. `sigma` is not
replaced by that derivative.

Comparison. Call `clifford_checks()`. Parse the returned `a4`.
Differentiate it in `L00`. Compare that derivative with a finite
difference. Evaluate the coefficient on an umbilic matrix and on one
matrix that is not umbilic. Compare the result with the shipped
`Lambda` string.

Hardest step. The derivative can be nonzero while the umbilic value
stays zero. Neither evaluation assigns `sigma`, and neither is the
lapse derivative of `Lhat^{uu}` on the clock.

Observable. The `a4` expression, its `L00` derivative, the umbilic
value, and one non-umbilic value.

Test.

1. Require the returned `a4` of both domains to equal
   `(4/5)` times `(u·L·u - tr(L)|u|^2/3)`.
2. Differentiate in `L00`. Require a finite difference to agree, and
   require the derivative to be nonzero at one rational gradient that
   is concentrated on the first axis.
3. Require the coefficient to vanish when `L` is a multiple of the
   identity. Require one other evaluation to be nonzero.
4. Require `Lambda` from `charge_identities()` to contain `zeta_s`
   and `sigma`. Do not substitute the derivative for either one.
   Do not march `rad.march`.

If the derivative is nonzero and `sigma` remains in `Lambda`, this
contact does not supply the clock coefficient. The class is an
unavailable physical input.

## Folded heat coefficient and the rate product

Assumption. The folded coefficient is the shipped one-angle a4 from
`analytic()` and from `clifford_checks()`. Equation (8) has one angle.
The two collar rates enter only by substituting their difference into
that one-angle square.

Error comparison. The algebraic identity for `|u_L - u_R|^2` is not
recomputed as a heat coefficient. The product collar has a vanishing
unsmeared `a4`. The normal-smearing square contains a product of two
rates only after that substitution. The matched-collar theorem does
not license the substitution. `delta` is not set to zero. A positive
square is not a stability mass.

Comparison. Call `analytic()`. Read the unsmeared value and the
normal-smearing string. Substitute the difference of two rates.
Isolate the coefficient of their product. Check that coefficient's
parent expression by a finite difference. On the curved `a4`, form
the polarization in two gradients and evaluate it on an umbilic
matrix. Read the exclusion sentences.

Hardest step. The product of the two rates appears only after the
substitution. The shipped return value has one angle. The curvature
polarization vanishes on the umbilic product collar, which is the
setting of equation (8).

Observable. The shipped one-angle normal-smearing coefficient, and
the product coefficient produced by the substitution.

Test.

1. Require `product_Delta_a4_unsmeared` to be 0. Require the
   normal-smearing string to be a constant times `theta_prime**2`.
2. Substitute the difference of two rates. Require the coefficient
   of their product to equal minus twice that constant. Check the
   substituted expression by a finite difference in one rate.
3. Require the curved polarization to vanish on an umbilic matrix.
4. Read the exclusion sentences. Do not introduce a two-collar clock
   equation.

A pass is a scoped obstruction. The class is an unresolved proof.

## Mouth scaling in the channel count

Assumption. The values named by the finite-mouth note, in its
particular scaling, are `a=1/(4pi^2 q)`, `epsilon=q^-2`,
`lambda_0=1/(16pi q^3)`, and `average D_j=O(delta log q)` with
`delta=log^-2 q`. The channel count `q` is not a rule-D march
parameter.

Error comparison. Absence of the four names from the march signature
was already recorded and is not this result. A formula in `q` becomes
a number only after `q` is given. The big-O statement is not a value
of `D_j`. Choosing a numerical `q` does not insert the four-tuple
into `Lambda`. The determinant is not set to one.

Comparison. Differentiate the three closed formulas in `q`. Compare
each derivative with a finite difference. Read the big-O sentence.
Compare `q` with the march parameters and with `Lambda`.

Hardest step. The three closed formulas have nonzero derivatives, so
they are not constants. The fourth member stays a big-O. A partial
tuple is not the data of the norm.

Observable. The three `q` derivatives, the big-O sentence, and the
march parameter list.

Test.

1. Require each `q` derivative to be nonzero and to agree with a
   finite difference.
2. Require the note to state `average D_j=O(delta log q)`.
3. Require `q`, `a`, `epsilon`, `D_j`, and `lambda_0` to be absent
   from the march parameters. Require `q` to be absent from `Lambda`.
4. Do not call `march`. Do not assign a numerical channel count to
   the clock.

A pass is a scoped obstruction. The class is an unavailable physical
input.

## Susceptibility at a numeric background

Assumption. A measured susceptibility is the shipped parallel
susceptibility evaluated at a nonzero numeric `B0`. The measurement
would be a number. The conditional axial window is not read.

Error comparison. A symbolic nonzero `B0` was already recorded and is
not this result. Substituting a number for `B0` leaves `cJ` and
`kappa` in the susceptibility. The four-null-photon coefficient is
not this measurement. `cJ` is not set to zero.

Comparison. Call `invariants()`. Substitute `B0 = 1` and
`B0 = 1/2`. Differentiate the first result in `cJ` and check it by a
finite difference. Compare the ratio of `(chi - 1)` at the two
backgrounds with the squared ratio of those backgrounds. Compare the
symbols with the march parameters.

Hardest step. The ratio test can pass while `chi` remains an
expression in `cJ` and `kappa`. A ratio of two expressions is not a
measured number.

Observable. The susceptibility at the two numeric backgrounds, and
the march parameter list.

Test.

1. Require both substituted susceptibilities to contain `cJ` and
   `kappa`, and to differ.
2. Require the `cJ` derivative at `B0 = 1` to agree with a finite
   difference.
3. Require `(chi - 1)` at the two backgrounds to stand in the ratio
   of the squared backgrounds.
4. Require the march parameters to omit `chi`, `B0`, and `cJ`. Do
   not set `cJ` to zero. Do not adopt the conditional axial window.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The shipped
`a4` curvature coefficient is
`4*(2*L00*u0**2 - L00*u1**2 - L00*u2**2 + 6*L01*u0*u1 + 6*L02*u0*u2 - L11*u0**2 + 2*L11*u1**2 - L11*u2**2 + 6*L12*u1*u2 - L22*u0**2 - L22*u1**2 + 2*L22*u2**2)/15`.
Its derivative in `L00` is `4*(2*u0**2 - u1**2 - u2**2)/15`. At the
rational gradient `(u0, u1, u2) = (1, 0, 0)` that derivative is `8/15`.
The umbilic value is 0. One non-umbilic value is `8/15`. `sigma` is not
replaced by that derivative. The class is an unavailable physical input.

The shipped normal-smearing coefficient is `4*theta_prime**2/15`. The
unsmeared product `a4` is 0. Substituting the difference of two rates
produces the product coefficient `-8/15`. The umbilic polarization of
the curved coefficient is 0. The folded coefficient is not proved. The
class is an unresolved proof. `delta` is not set to zero.

The scaling derivatives are `a`: `-1/(4*pi**2*q**2)`, `epsilon`:
`-2/q**3`, and `lambda_0`: `-3/(16*pi*q**4)`. Their finite-difference
relative gaps are the receipt rationals
`-1/63999980000001`,
`-21333328888889/455110826666725333328888889`, and
`-9557329493333834666646666667/87381251413363029328106667130666646666667`.
`average D_j` remains `O(delta log q)`. No channel count is assigned.
The class is an unavailable physical input. The determinant is not set
to one.

At `B0 = 1` the susceptibility is `2*cJ*kappa**4 + 1`. At `B0 = 1/2`
it is `(cJ*kappa**4 + 2)/2`. The `cJ` derivative at `B0 = 1` is
`2*kappa**4`. `chi` is not a number. The march parameter `kappa` is
not used as the gravitational `kappa`. The class is an unavailable
physical input. `cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S087 product-throat derivative is not
recomputed. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated.
No successor model is opened.
