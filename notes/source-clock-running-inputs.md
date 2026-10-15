# Running coupling, collar ratio, exterior length, and the magnetic variation

S091, E003, 2026-10-14. This note takes the four inputs named by S090.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S089 `L00` derivative is not recomputed as a
new result. The S090 substitution difference is not recomputed as a new
result. No new `Lambda` interval is sought. No replacement area target
is written before the test. The area budget stays unevaluated.
`zeta_s`, `sigma`, `cJ`, `c_theta`, and `delta` are not set to zero.
An endpoint contact energy is not cancelled. A positive quadratic term
is not used as a stability mass.

## Renormalization condition for `c_theta(mu)`

Assumption. A renormalization condition is an equality that assigns one
number to `c_theta(mu)` at one named scale `mu`. The curved-collar note
states that `c_theta(mu)` multiplying `I_theta` has beta function
`-b_theta`, with `b_theta=4/[5(4pi)^2]`. That beta function fixes the
scale derivative. It does not assign the integration constant.

Error comparison. S090 recorded `b_theta = 1/(20*pi**2)` and did not
replace `sigma`. Integrating the beta function adds `c_theta(mu_0)`.
The integration constant is not set to zero. `sigma` is not replaced
by the running coupling.

Comparison. Read the shipped beta-function statement. Integrate
`mu partial_mu c_theta = -b_theta` from a named `mu_0`. Differentiate
the integral in `mu` and check that derivative by a finite difference.
Compare `c_theta` and `mu` with the march parameters and with the
shipped sigma string.

Hardest step. The scale derivative can be a pure number while
`c_theta(mu_0)` remains an arbitrary constant. The pure number is not
a renormalization condition.

Observable. The integrated running coupling and the symbol that the
beta function does not fix.

Test.

1. Require `b_theta=4/[5(4pi)^2]` and the beta-function sentence.
   Set the scale derivative to `-b_theta`.
2. Integrate from a named `mu_0`. Require `mu` times the `mu`
   derivative of that integral to equal `-b_theta`, and require a
   finite difference to agree.
3. Require the integration constant to remain a symbol. Do not set it
   to zero.
4. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta` and `mu`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu)` is recorded.

## Two-collar ratio

Assumption. The compared coefficient is the cross term between the sum
of two equation-(8) integrals and the substituted square. Equation (8)
covers one angle. The coupled note says unequal collars require their
additional invariants and are not covered by (8).

Error comparison. The sum of two one-angle terms and the single
substituted square differ by the cross term. That ratio is not adopted
as a two-collar theorem. `delta` is not set to zero. No two-collar
clock equation is introduced.

Comparison. Call `analytic()`. Form the diagonal sum and the
substituted square from the shipped normal-smearing coefficient.
Divide the cross term by the diagonal sum. Evaluate that ratio at the
predeclared point `theta_L = 1`, `theta_R = 1/2`. Check the `theta_L`
derivative of the cross term by a finite difference.

Hardest step. The ratio can be a number at one point while the
additional invariants named by the note stay unstated. A number at one
point is not a proved two-collar coefficient.

Observable. The ratio of the cross term to the diagonal sum, and
whether `analytic()` returns a two-collar bound.

Test.

1. Require the shipped normal smearing to be a constant times
   `theta_prime**2`. Record the cross term
   `substituted square minus diagonal sum`.
2. Record the ratio at `theta_L = 1`, `theta_R = 1/2`. Require a
   finite difference of the cross term to agree with its derivative.
3. Require the returned keys to contain no two-collar error bound.
   Require the exclusion sentence. Do not introduce a two-collar
   equation.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Equality for one numeric `D_j`

Assumption. An equality assigning `D_j` would be an equation `D_j`
equals a number, written as clock data before a test. The mouth note
keeps `arbitrary fixed D_j>0` and writes `average D_j=O(delta log q)`.
An order symbol is not that equality.

Error comparison. The singlet trace depends on `T+D_j` through
`coth(k(T+D_j)/2)`. The trace derivative is not that equality. Moving
`D_j` moves the trace, so the trace does not select one number.
`delta` is not set to zero. The determinant is not set to one.

Comparison. Differentiate one mouth's contribution
`eta*coth(k(T+D_j)/2)/(2k)` in `D_j`. Check that derivative by a
central quotient and by a numerical evaluation at `eta = 1`, `k = 1`,
`T = 1`, `D_j = 1`. Compare `D_j` with the march parameters and with
`Lambda`.

Hardest step. A nonzero derivative shows that the trace does not pin
`D_j`. It does not supply the missing equality.

Observable. The `D_j` derivative of one mouth term, and whether it is
zero.

Test.

1. Require the order statement and the arbitrary-positive statement.
2. Require the derivative to equal
   `-eta*csch(k(T+D_j)/2)**2/4`, and require the central quotient to
   agree. Record the value at the predeclared point. Require it to be
   nonzero.
3. Require `D_j` to be absent from the march parameters and from
   `Lambda`. Do not assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Magnetic first variation

Assumption. A measurement of `cJ` would be a number returned at
`B0 = 1` and gravitational `kappa = 1`. The electromagnetic note says
that on a purely magnetic background both `Delta L` and its first
functional variation vanish for every `cJ`, and that this does not
establish the physical value `cJ=0`.

Error comparison. The vanishing first variation is not that number.
The quadratic correction remains proportional to `cJ`. The conditional
axial window is not adopted. `cJ` is not set to zero.

Comparison. Call `invariants()`. Rebuild `Delta L = cJ*kappa**4*J**2/16`
and differentiate it in `J` at `J = 0`. Check that derivative by a
finite difference. Substitute the unit background into the shipped
quadratic correction and record whether `cJ` remains.

Hardest step. The first variation can vanish for every `cJ` while the
quadratic response stays linear in `cJ`. A zero variation is not a
measured coupling.

Observable. The first variation at `J = 0`, and the quadratic
correction at `B0 = 1`, `kappa = 1`.

Test.

1. Require the returned magnetic background action and first variation
   to be zero. Require the note sentences that the background does not
   measure `cJ` and does not establish `cJ=0`.
2. Require `Delta L` and its `J` derivative to vanish at `J = 0`, and
   require a finite difference to agree.
3. Substitute `B0 = 1` and `kappa = 1` into the shipped quadratic
   correction. Require `cJ` to remain. Do not adopt the conditional
   axial window.
4. Require the march parameters to omit `cJ`, `chi`, and `B0`. Do not
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
beta function is `-1/(20*pi**2)`. The integrated running coupling is
`c_theta(mu0) - b_theta*log(mu/mu0)`. The integration constant
`c_theta(mu0)` remains a symbol. `sigma` is not replaced by the running
coupling. The class is an unavailable physical input. `c_theta` is not
set to zero.

The ratio of the cross term to the diagonal sum is `-4/5` at
`theta_L = 1`, `theta_R = 1/2`. The shipped return keys contain no
two-collar error bound. The class is an unresolved proof. That ratio
is not adopted as a two-collar theorem. `delta` is not set to zero.

The mouth-trace derivative is `-eta*csch(k*(T+D_j)/2)**2/4`. At the
unit point it is `-1/(4*sinh(1)**2)`, which is nonzero. `numeric_D_j`
is null. The class is an unavailable physical input. The determinant
is not set to one.

The magnetic first variation is zero for every `cJ`. At `B0 = 1` and
`kappa = 1` the quadratic correction is `cJ*e_parallel**2`. `cJ` is
not a number. The march parameter `kappa` is not used as the
gravitational `kappa`. The class is an unavailable physical input.
`cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S090 substitution difference is not
recomputed as a new result. No new `Lambda` interval is recorded. The
area budget is not evaluated. No successor model is opened.
