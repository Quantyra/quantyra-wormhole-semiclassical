# Scale increment, collar mismatch, trace gap, and the parallel coefficient

S092, E003, 2026-10-14. This note takes the four inputs named by S091.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S091 running coupling is not recomputed as a
new result. The S091 collar ratio is not recomputed as a new result.
The S091 trace derivative is not recomputed as a new result. No new
`Lambda` interval is sought. No replacement area target is written
before the test. The area budget stays unevaluated. `zeta_s`,
`sigma`, `cJ`, `c_theta`, and `delta` are not set to zero. An endpoint
contact energy is not cancelled. A positive quadratic term is not used
as a stability mass.

## Number for `c_theta(mu0)` at a named scale

Assumption. A number for `c_theta(mu0)` would be the value of the
running coupling at one numerical scale named by the curved-collar
note. That note writes `c_theta(mu)` and does not name a numerical
scale.

Error comparison. The comparison endpoints `mu = 1` and `mu = E` are
not a physical assignment of `mu0`. The definite increment between
those endpoints is not that value. The integration constant is not set
to zero. `sigma` is not replaced by the increment.

Comparison. Integrate the shipped beta function from `mu = 1` to
`mu = E`. Compare the increment with the value of the running coupling
at `mu = 1`. Check the increment against a numerical evaluation.
Compare `c_theta` and `mu` with the march parameters and with the
shipped sigma string.

Hardest step. The increment between two scales can be a pure number
while the value at either endpoint still contains `c_theta(mu0)`. The
pure number is not a renormalization condition.

Observable. The definite increment, and whether the value at `mu = 1`
is a number.

Test.

1. Require the shipped formula `b_theta=4/[5(4pi)^2]` and the absence
   of a numerical scale assigned to `c_theta`.
2. Form the increment from `mu = 1` to `mu = E`. Require it to be
   independent of the integration constant, and require a numerical
   evaluation to agree.
3. Require the value at `mu = 1` to contain the integration constant.
   Do not set that constant to zero.
4. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta` and `mu`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Unequal-collar invariants

Assumption. The additional invariants are the quantities the coupled
note says vanish on the product and are not covered by equation (8).
The curved-collar theorem excludes unequal collar perturbations. Its
traceless tensor is `Lhat_ab=L_ab-L gamma_ab/3`.

Error comparison. The predeclared collars are not the same matrix.
The left diagonal is `(1, 0, 0)` and the right diagonal is `(0, 1, 0)`,
with gradient `(1, 0, 0)`. The collar mismatch is not adopted as a
two-collar theorem. `delta` is not set to zero.

Comparison. Call `analytic()`. Contract `u·Lhat·u` on each predeclared
collar and subtract. Differentiate that mismatch in the left `L00`
entry and check the derivative by a finite difference. Compare the
mismatch with the shipped product unsmeared coefficient.

Hardest step. The mismatch can be nonzero while the product return
stays zero. A nonzero mismatch is not the missing two-collar
coefficient.

Observable. The mismatch of the two contractions, and whether the
product return contains it.

Test.

1. Require the exclusion sentences and the definition of `Lhat_ab`.
2. Record the mismatch at the predeclared collars. Require a finite
   difference in the left `L00` entry to agree with its derivative.
3. Require the shipped product unsmeared coefficient to be zero and
   the returned keys to omit the vanishing product invariants. Do not
   introduce a two-collar equation.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Equality setting one `D_j`

Assumption. An equality setting `D_j` would be a finite positive root
of the free-diagonal match `coth(k(T+D_j)/2) = 1`. The mouth note gives
the trace and says its difference from `eta/(2k)` decays at large `k`.
That difference is not the missing clock equality.

Error comparison. The trace gap is not that root. The predeclared point
is `eta = 1`, `k = 1`, `T = 1`, `D_j = 1`. `delta` is not set to zero.
The determinant is not set to one.

Comparison. Subtract `eta/(2k)` from one mouth term. Test whether
`coth(u) - 1` has a finite positive root. Evaluate the gap at the
predeclared point and check it by a numerical evaluation. Differentiate
the clock expression for `Lambda` in `D_j`.

Hardest step. The gap can be a nonzero number while no finite `D_j`
removes it, and the clock derivative in `D_j` can be zero. Neither fact
assigns `D_j`.

Observable. The trace gap at the predeclared point, and the clock
derivative in `D_j`.

Test.

1. Require the trace formula and the free-diagonal sentence.
2. Require `coth(u) - 1` to have no finite positive root. Record the
   gap at the predeclared point and require a numerical evaluation to
   agree.
3. Require the derivative of the clock `Lambda` expression in `D_j` to
   vanish. Require `D_j` to be absent from the march parameters. Do not
   assign a channel count.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Measured parallel quadratic coefficient

Assumption. A measured parallel coefficient would be a number for the
second electric derivative at `B0 = 1` and gravitational `kappa = 1`.
The electromagnetic note's `k = 2000` axial illustration is a different
quantity.

Error comparison. The axial illustration is not that coefficient.
Substituting the note's example `cJ = 1` is not a measurement. The
conditional axial window is not adopted. `cJ` is not set to zero. A
positive quadratic term is not used as a stability mass.

Comparison. Call `invariants()`. Differentiate the shipped quadratic
correction twice in the parallel electric field at the unit background.
Check that second derivative by a second finite difference. Evaluate
the note's axial ratio `16*pi*gamma*(1+1/C)` at `k = 2000` and compare
it with the parallel coefficient. The displayed `theta/cJ` digits in
that note are the comparison target.

Hardest step. The second electric derivative can stay linear in `cJ`
while an axial ratio at one `k` is a tiny number. The tiny number is
not a measurement of the parallel coefficient.

Observable. The second parallel derivative at the unit background, and
whether it equals the axial illustration.

Test.

1. Require the returned quadratic correction and the sentences that
   the axial shifts are not the polar bound and that magnetic matching
   alone cannot measure `cJ`.
2. Record the second derivative at `B0 = 1`, `kappa = 1`. Require a
   second finite difference to agree. Require `cJ` to remain.
3. Require the axial ratio at `k = 2000` to agree with the displayed
   digits and to differ from the parallel coefficient. Do not adopt
   the conditional axial window.
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

Filled from the receipt after the two agreeing checker runs. The
definite increment from `mu = 1` to `mu = E` is `-1/(20*pi**2)`. The
value at `mu = 1` is `c_theta(mu0) + log(mu0)/(20*pi**2)`. That value
is not a number. The curved-collar note names no numerical scale.
`sigma` stays `2*c_theta*Lhat^{uu}`. The class is an unavailable
physical input. `c_theta` is not set to zero. The running coupling is
not recomputed as a new result.

The predeclared collar mismatch is `1`. The shipped product unsmeared
coefficient is `0`, and the returned keys omit the vanishing product
invariants. The class is an unresolved proof. That mismatch is not
adopted as a two-collar theorem. `delta` is not set to zero. The
collar ratio is not recomputed as a new result.

The trace gap at `eta = k = T = D_j = 1` is `1/(exp(2)-1)`, recorded
as `1/((-1 + E)*(1 + E))`. `coth(u) - 1` has no finite positive root.
The clock derivative in `D_j` is `0`. `numeric_D_j` is null. The class
is an unavailable physical input. The determinant is not set to one.
The trace derivative is not recomputed as a new result.

At `B0 = 1` and gravitational `kappa = 1` the second parallel
derivative is `2*cJ`. It is not a number and it is not the displayed
axial ratio. The example `cJ = 1` is not a measurement. The march
parameter `kappa` is not used as the gravitational `kappa`. The class
is an unavailable physical input. `cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. The S088 relative series is not differentiated
again. No new `Lambda` interval is recorded. The area budget is not
evaluated. No successor model is opened.
