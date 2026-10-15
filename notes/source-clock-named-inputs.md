# Lapse derivative, jet-cross, determinant norm, and magnetic background

S087, E003, 2026-10-14. This note takes the four inputs named by S086.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. No new `Lambda` interval is sought.
No replacement area target is written before the test. The area budget
stays unevaluated. `zeta_s`, `sigma`, `cJ`, and `delta` are not set
to zero. An endpoint contact energy is not cancelled. A positive
quadratic term is not used as a stability mass.

## Lapse derivative of `Lhat^{uu}`

Assumption. The retained cutoff is the product throat in the radiation
note, `ds^2 = r_e^2 (-dt^2 + dz^2)/z^2 + r_e^2 dOmega^2`, at fixed
induced `gamma_tt = -r_e^2/epsilon^2`. The positive lapse is
`N = r_e/epsilon`. The shipped sphere-integrated weight is the object
whose lapse derivative S086 left unevaluated. `sigma` stays a symbol
in `Lambda`.

Error comparison. S083 records that this weight is `O(epsilon)` and
drops out of the JT clock. That factor is geometric. It is not an
assignment of `c_theta` or of `sigma`. A sphere radius that varies
along the normal is outside this product metric. The chain rule can
vanish at `N -> infinity` while the `epsilon` derivative stays
nonzero. `sigma` is not replaced by that limit.

Comparison. Call `geometry_identities()` and differentiate its
integrated density with respect to `N` after `epsilon = r_e/N`.
Compare that derivative with the chain rule through `epsilon`.
Compare its limit at the conformal boundary with the shipped `Lambda`
string from `charge_identities()`. The S086 fixed-`Lhat` gap is not
recomputed, and the clock is not marched.

Hardest step. The boundary that defines the JT clock is
`epsilon -> 0`. A derivative that exists only at finite `epsilon`
is the cutoff weight, not a coefficient inside the retained clock
sum. Substituting the boundary value for `sigma` would erase the
symbol the sum still carries.

Observable. The lapse derivative of the shipped integrated weight,
its limit at `N -> infinity`, the `epsilon` derivative at
`epsilon = 0`, and the shipped `Lambda` string.

Test.

1. Parse the integrated density returned by `geometry_identities()`.
   Differentiate in `epsilon` and in `N`. Require the chain rule.
   A finite difference in `N` must agree with the analytic `N`
   derivative at one finite cutoff.
2. Require the `epsilon` derivative at `epsilon = 0` to be a nonzero
   explicit expression. Require the weight and its `N` derivative to
   have limit 0 as `N -> infinity`.
3. Require `Lambda` from `charge_identities()` to contain `zeta_s`
   and `sigma`. Do not substitute the boundary limit for either
   coefficient.
4. Do not march `rad.march`.

If the boundary limit of the derivative is 0 and `sigma` remains in
`Lambda`, the separation is an unavailable physical input. The limit
is not an admissible range and is not a value of `sigma`.

## Jet-cross term of unequal collars

Assumption. The jet-cross term is the quadratic scalar in the
published square `(D_a chi)(D^a chi)`, whose trace contains
`L_ab L^{ab}`. Split one orthonormal jet as
`L_bar + (delta/2) J` and the other as `L_bar - (delta/2) J`. The
candidate is the average of the two quadratic scalars. `delta` stays
free.

Error comparison. The matched-collar theorem folds identical jets
onto one boundary. Unequal collar perturbations lie outside it.
Equation (8) of the coupled-domain note has no collar difference.
The average of two single-collar squares is not that folded
coefficient, and it is not a two-collar equation of motion.

Comparison. Expand the average. Compare its second `delta`
derivative with a finite difference. Call `clifford_checks()` so the
shipped trace that contains `L_ab L^{ab}` is the one used. Read the
exclusion sentences from the curved-fermion note and the
coupled-domain note.

Hardest step. The linear cross cancels in the average, so a first
derivative at `delta = 0` does not display the jet. Setting
`delta = 0` would hide the quadratic piece. A positive second
derivative is not a stability mass.

Observable. The averaged quadratic scalar, its value at `delta = 0`,
and its second `delta` derivative.

Test.

1. Require the average at `delta = 0` to equal the mean jet's square.
2. Differentiate twice. A finite difference of the average must
   agree. The second derivative must be the square of `J` over 2,
   as a symbol, and must not be replaced by zero.
3. Require `clifford_checks()` to pass. Read the two exclusion
   sentences. Do not introduce a two-collar clock equation.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible range for `delta` is recorded.

## Norm on the interacting determinant

Assumption. The norm is the published Schatten bound
`|log det_3(I+X)| <= r/[3(1-r^(1/3))]` for `||X||_S3^3 <= r < 1`.
The stand-in remains `G(z) = (1+z) I_2`, free potential `v I_2`, and
`W = [[0, m], [m, 0]]`. The predeclared point remains
`z = 0`, `v = 1/100`, `m = 1/2`. `X` is the interaction relative to
the free scalar factor.

Error comparison. The bound controls the value of `log det_3`. It
does not control a `z` derivative, and a bound on a value does not
become a bound on a derivative by differentiating the right-hand
side. The finite-mouth note does not let a Gaussian result control
the corrected background's nonlocal determinant. This matrix is not
the wormhole determinant. The determinant is not set to one.

Comparison. Compute the singular values of `X` at the predeclared
point and the Schatten radius `r = 2 |s|^3`. Compare `|log det_3|`
with the published right-hand side when `r < 1`. Compare the
interacting `z` derivative with the free `z` derivative. A finite
difference checks the interacting derivative.

Hardest step. A finite `r` inside the hypothesis can hold while the
interaction still moves the `z` derivative. The clock polynomial is
not given this norm as a factor.

Observable. `r`, the bound, `|log det_3|`, and the two `z`
derivatives at the predeclared point.

Test.

1. Build `X` from the predeclared matrix. Require an explicit `r`.
   If `r < 1`, require `|log det_3|` to lie at or below the published
   expression. If `r >= 1`, record that the hypothesis fails.
2. Require the interacting `z` derivative to differ from the free
   `z` derivative. Check the interacting derivative by a finite
   difference.
3. Do not treat the bound as a bound on that derivative. Do not
   recompute the free Gaussian energy. Read the scope sentence from
   the finite-mouth note and the bound sentence from the
   nonlinear-current note.

A value bound on `log det_3` of this stand-in is not a norm on the
wormhole determinant's `z` derivative. The class is an uncontrolled
approximation.

## Parallel susceptibility and `B0`

Assumption. The zero magnetic background is the point `B0 = 0` in
the shipped susceptibility `1+2*cJ*kappa^4*B0^2` and in the shipped
quadratic correction `cJ*kappa^4*B0^2*e_parallel^2`. The conditional
axial window is not read.

Error comparison. The shipped magnetic background has vanishing
action and first variation. At `B0 = 0` the susceptibility is 1 for
every `cJ`, and the quadratic correction is 0 for every `cJ`. That
identity does not determine `cJ`. Setting `cJ` to zero would be an
extra condition. The rule-D march does not take `chi` or `B0`.

Comparison. Call `invariants()`. Substitute `B0 = 0` into the
returned strings. Compare the result with the parameter list of
`march` and with `Lambda` from `charge_identities()`.

Hardest step. The solved formula `(chi - 1)/(2*B0**2*kappa**4)` is
undefined at `B0 = 0`. The reference background is that point's
first-variation statement. A formula that becomes an identity does
not assign `cJ`.

Observable. The susceptibility and the quadratic correction at
`B0 = 0`, and the march parameter list.

Test.

1. Parse both shipped strings. Require the susceptibility at
   `B0 = 0` to equal 1 with `cJ` absent, and the quadratic
   correction at `B0 = 0` to equal 0 with `cJ` absent.
2. Require `Lambda` to be independent of `cJ` and of `B0`. Require
   the march parameters to be `(eta, lam, kappa, z0, v0, method)`.
3. Do not call `march`. Do not set `cJ` to zero. Do not adopt the
   conditional axial window.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The
shipped integrated density is
`-(8*pi/3)*epsilon*(1-epsilon**2*Sch)*rate**2`. Its lapse derivative
is `8*pi*r_e*rate**2*(N**2 - 3*Sch*r_e**2)/(3*N**4)`. The `epsilon`
derivative at zero is `-8*pi*rate**2/3`. The weight and the lapse
derivative both have limit 0 as `N -> infinity`. The finite-difference
relative gap on the lapse derivative is `8.691635025882847e-11`.
`Lambda` remains
`lambda_zeta+lambda_theta=(pi*ell**2/N)*(zeta_s+sigma)`. `sigma` is
not replaced by the boundary limit. The class is an unavailable
physical input.

The averaged jet square has second derivative
`(j0**2 + j1**2 + j2**2)/2`. The finite-difference relative gap is
`6.18487904569359e-11`. The folded coefficient is not proved. The
class is an unresolved proof. `delta` is not set to zero.

At the predeclared point the Schatten radius is `250000/1030301`,
inside `r < 1`. `|log det_3|` is
`0.03606155179153857880017716417480772858782` and the published bound
is `0.2149548739358667360411194995898728439925`. The free derivative
remains `2/101` and the interacting derivative remains `-4798/7701`.
The finite-difference relative gap is `7.176209872275215e-10`. The
bound does not control the `z` derivative. The class is an
uncontrolled approximation.

At `B0 = 0` the susceptibility is 1 and the quadratic correction is
0. `cJ` is not constrained. The march parameters remain
`eta, lam, kappa, z0, v0, method`. The class is an unavailable
physical input. `cJ` is not set to zero.

`opening_achieved` is false. `s074_gate_closed` is false. The S085
march is not repeated. No new `Lambda` interval is recorded. The area
budget is not evaluated. No successor model is opened.
