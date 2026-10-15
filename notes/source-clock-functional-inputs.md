# Functional, curvature combination, self energy, and the quadratic quotient

S097, E003, 2026-10-14. This note takes the four inputs named by S096.
The comparison, the hardest step, the observable, and the test below
are fixed before any checker run. No result from a fit is stated in
this section.

The S085 march is not repeated. The S088 relative series is not
differentiated again. The S096 Weyl coefficient is not recomputed as a
new result. The S096 traced term is not recomputed as a new result. The
S096 little-o derivative is not recomputed as a new result. The S096
susceptibility ratio is not recomputed as a new result. No new `Lambda`
interval is sought. No replacement area target is written before the
test. The area budget stays unevaluated. `zeta_s`, `sigma`, `cJ`,
`c_theta`, `d_theta`, and `delta` are not set to zero. An endpoint
contact energy is not cancelled. A positive quadratic term is not used
as a stability mass.

## Weyl-invariant functional `I_theta`

Assumption. A number for `c_theta(mu0)` would be the Weyl-invariant
functional `I_theta=integral sqrt(gamma)Lhat^ab u_a u_b`.

Error comparison. That functional is the boundary integral multiplied
by the finite coupling. A functional is not a value of the coupling.
That functional is not a value of `c_theta(mu0)`. `c_theta` is not set
to zero and is not set equal to `I_theta`.

Comparison. Require the sentence that names `I_theta`. Differentiate
the product of `c_theta(mu0)` with `I_theta` in the functional and
check the derivative by a finite difference. Compare the functional
with the endpoint that still contains `c_theta(mu0)`.

Hardest step. The functional can be displayed while the coupling stays
unassigned. Naming `I_theta` is not a value of `c_theta(mu0)`.

Observable. The derivative of the product in `I_theta`, and whether the
endpoint equals the functional.

Test.

1. Require the sentence that names
   `I_theta=integral sqrt(gamma)Lhat^ab u_a u_b`.
2. Record the derivative of `c_theta(mu0) I_theta` in `I_theta`. Require
   a finite difference to agree. Require the endpoint to keep
   `c_theta(mu0)` and to differ from `I_theta`.
3. Require the shipped sigma string to stay `2*c_theta*Lhat^{uu}`.
   Require the march parameters to omit `c_theta` and `I_theta`. Do not
   set `c_theta` equal to `I_theta`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `c_theta(mu0)` is recorded.

## Curvature combination `Delta tr[(D chi)^2 S]`

Assumption. The unequal-collar invariant would be the curvature
combination `Delta tr[(D chi)^2 S]=-2L u^2-4L_ab u^a u^b`.

Error comparison. The curved-collar note displays that combination on
the product. A product heat identity is not a two-collar theorem.
Unequal collars are not covered by equation (8). This combination is
not the `theta` derivative of the traced `S_:aa` term. `delta` is not
set to zero.

Comparison. Differentiate the displayed combination in `L`. Check that
derivative by a finite difference. Compare the combination and the
derivative with the recorded mismatch `1`.

Hardest step. The combination can be differentiated while `L_ab`
remains. A curvature derivative is not the missing invariant.

Observable. The derivative in `L`, and whether the combination equals
the mismatch.

Test.

1. Require the displayed combination and the exclusion of unequal
   collars.
2. Record the derivative in `L`. Require a finite difference to agree.
   Require the derivative and the combination to differ from `1`.
3. Do not adopt the combination as a two-collar theorem. Do not
   identify it with the shipped trace.

A pass is a scoped obstruction. The class is an unresolved proof.
No admissible two-collar coefficient is recorded.

## Clock-independent mouth self energy

Assumption. A number for one `D_j` would be the clock-independent mouth
self energy inside the undifferentiated energy.

Error comparison. The mouth note says no uniform bound on that energy
is asserted. Clock independence is not an assignment. The self energy
is not the little-o statement `O(1)=o(q delta)`. `delta` is not set to
zero. The determinant is not set to one.

Comparison. Require the displayed statement. Differentiate a
clock-independent self energy in `D_j` and check that derivative by a
finite difference. Differentiate the clock expression for `Lambda` in
`D_j`.

Hardest step. The self-energy derivative can vanish because the symbol
does not depend on `D_j`. Vanishing does not assign `D_j`.

Observable. Whether the self energy is a number, and the clock
derivative in `D_j`.

Test.

1. Require the sentence that the undifferentiated energy contains large
   clock-independent mouth self energies.
2. Require the derivative of that self energy in `D_j` to vanish, and
   require a finite difference to agree. Require the clock derivative
   in `D_j` to vanish.
3. Require `D_j` to be absent from the march parameters. Do not assign
   a channel count. Do not evaluate the missing bound as an equality.

A pass is a scoped obstruction. The class is an unavailable physical
input. No numeric `D_j` is recorded.

## Quadratic correction divided by `e_parallel^2`

Assumption. A measured `cJ` would be the quadratic correction divided
by `e_parallel^2` at gravitational `kappa = 1`.

Error comparison. The quotient can keep `B0`. Returning `cJ B0^2` does
not measure `cJ`. The conditional axial window is not adopted. `cJ` is
not set to zero. A positive quadratic term is not used as a stability
mass.

Comparison. Call `invariants()`. Divide the shipped quadratic correction
by `e_parallel^2`, set gravitational `kappa = 1`, and check the
quotient at one nonzero `cJ` and one nonzero `B0`.

Hardest step. The quotient can stay proportional to `cJ` while `B0`
remains. The quotient is not a mass.

Observable. The quotient at gravitational `kappa = 1`, and whether it
is a number.

Test.

1. Require the returned quadratic correction and the sentence that
   magnetic background matching alone cannot measure `cJ`.
2. Record the quotient at gravitational `kappa = 1`. Require a
   numerical evaluation at one nonzero `cJ` and one nonzero `B0` to
   agree. Require `cJ` and `B0` to remain.
3. Require the march parameters to omit `cJ`, `chi`, `B0`, and
   `e_parallel`. Do not use the march's `kappa` as the gravitational
   `kappa`.

A pass is a scoped obstruction. The class is an unavailable physical
input. No admissible value of `cJ` is recorded.

## What a pass does not do

These four results do not match every declared-order term. They do
not check S074.2–S074.6. They do not start S075. They do not evaluate
the area budget. They do not accept the scientific objective of E003.
They are not an impossibility theorem.

## Benchmark result

Filled from the receipt after the two agreeing checker runs. The derivative
of `c_theta(mu0) I_theta` in `I_theta` is `c_theta_mu0`. It is not a number
and does not assign `c_theta(mu0)`. The endpoint keeps `c_theta(mu0)` and
does not equal `I_theta`. `sigma` stays `2*c_theta*Lhat^{uu}`. The class is
an unavailable physical input. `c_theta` is not set to zero and is not set
equal to `I_theta`. The Weyl coefficient is not recomputed as a new result.

The derivative of the curvature combination in `L` is `-2*u2`. It is not
the recorded mismatch `1`. The combination itself is not that mismatch.
The class is an unresolved proof. The combination is not adopted as a
two-collar theorem. `delta` is not set to zero. The traced term is not
recomputed as a new result.

The clock-independent mouth self energy is not a number. Its derivative
in `D_j` is `0`. The clock derivative in `D_j` is `0`. `numeric_D_j` is
null. The class is an unavailable physical input. The missing bound is
not evaluated as an equality. The little-o derivative is not recomputed
as a new result.

The quadratic correction divided by `e_parallel^2` at gravitational
`kappa = 1` is `B0**2*cJ`. That expression is not a number and does not
measure `cJ`. The class is an unavailable physical input. `cJ` is not set
to zero. The susceptibility ratio is not recomputed as a new result.

`opening_achieved` is false. `s074_gate_closed` is false. The S085 march
is not repeated. The S088 relative series is not differentiated again. No
new `Lambda` interval is recorded. The area budget is not evaluated. No
successor model is opened.
