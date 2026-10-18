# Variable-radius kernel comparison and the coupled exterior seam

2026-10-18 HST, E003/S139 supporting S074/S075. Ordinary argument before
the focused checker. **The coupled quantum seam is unsolved. The
nonlocal kernel difference is an unresolved proof. No physical gate
closes.**

This attempt uses the reconstructed product inverse in
[stationary-product-weighted-response](stationary-product-weighted-response.md)
and the prescribed-source mass and temperature map in
[stationary-reduced-cold-response](stationary-reduced-cold-response.md).
The constant-background bubble factor is the one already bounded in
[stationary-reduced-parameter-form-bounds](stationary-reduced-parameter-form-bounds.md).
The signed moment is the fast outgoing condition in
[stationary-scalar-phase-moment](stationary-scalar-phase-moment.md).
None of `zeta_s`, `sigma`, `cJ`, `c_theta`, or `delta` is set to zero.

## Frozen test

Comparison: the invariant kernel of the complete-product inverse
`M(P2)^-1` against the actual variable-radius reduced Hessian on the
cold chart.

Hardest step: bound the nonlocal part of that Hessian difference in
the weighted jet norm `X_L^2`, including the radius derivative of the
paired common-subtracted bubble.

Observable: the horizon lapse constraint and the geometric mass return,
with the kernel correction included, at the retained fixed temperature.

Test: accept a solved coupled seam only when that difference has a
proved bound in `X_L^2` and the mass and temperature constraints are
solved with the correction included. A product of the local conformal
gap with the retained weighted norm is not that bound. No area target
is written before this attempt, so the area budget is not evaluated.

## 1. Retained product inverse

On a cold ball of radius `L<=84`, the normalized product inverse obeys

    ||f(P2) F||_(X_L^2) <= 1.3e8 ||F||_(X_L^2)

for a compact radial source, with the cosh weight retained. The same
note gives the invariant kernel majorant

    ||K(d)|| <= [4e5 exp(-2.01 d)+1e3 exp(-1.999 d)]
                                        /sqrt(1-exp(-2d)).

The map is restriction to the ball, after the complete-H2 inverse,
after zero extension. It selects the decaying complete-H2 tail. It
does not supply a Dirichlet or Robin datum at `chi=L`, and it does
not invert a variable-radius operator.

The arithmetic reserve used for `1.3e8` is

    (22/7)(4e7+2000*84*(6/5)) < 1.3e8.

That reserve is rechecked here. It is not a new spatial theorem.

## 2. Local conformal coefficient

In two dimensions, if `ghat=exp(2 sigma) g` and `g` is flat in the
coordinate chart, the Laplace-Beltrami operator is exactly

    Delta_ghat = exp(-2 sigma) Delta_g.

The flux identity is `sqrt(ghat) ghat^{ij}=delta^{ij}`. Division by
the area density `exp(2 sigma)` produces the factor. This is a local
differential identity. It is not a transformation law for the nonlocal
quantum Hessian.

On the cold chart the retained guard is `|sigma_cold|<1/20`. The
coefficient `exp(-2 sigma)` therefore lies between `exp(-1/10)` and
`exp(1/10)`. Its largest absolute departure from 1 on that interval
is `exp(1/10)-1`, because

    exp(a)-1 > 1-exp(-a)

for `a=1/10>0`, which is the identity
`exp(a)+exp(-a)-2=(exp(a/2)-exp(-a/2))^2`.

Two independent enclosures give

    21/200 < exp(1/10)-1 < 11/100.

The direct enclosure truncates the exponential series after the
degree-five term and bounds the positive tail by a geometric majorant.
The reciprocal enclosure uses the alternating series for `exp(-1/10)`
and inverts the resulting rational bounds. Both are finite rational
comparisons. Neither identifies this gap with a quantum kernel.

## 3. The conformal substitution is uncontrolled

Let `A` be invertible and `A+E` be invertible. The resolvent identity

    (A+E)^-1 - A^-1 = -A^-1 E (A+E)^-1

is exact. It becomes a norm bound only after `E` is estimated in the
same operator norm as `A^-1`.

Write the variable-radius Hessian as `H_prod+E_local+E_bubble` on the
common domain where that sum is defined. `E_local` may contain the
conformal Laplacian factor above. `E_bubble` is the change in the
paired common-subtracted bubble when the radius function varies.
The retained notes construct `H_prod^-1`. They do not estimate
`E_bubble`.

Substituting `E_local+E_bubble` by `(exp(-2 sigma)-1) H_prod`, and
then multiplying the weighted norm `1.3e8` by the gap
`exp(1/10)-1`, is a formal size. The rational endpoints are exact:

    (21/200)*1.3e8 = 13650000,
    (11/100)*1.3e8 = 14300000.

The series enclosures of the same gap, multiplied by `1.3e8`, lie
strictly between those endpoints. That size belongs to the
substitution in the normalized weighted norm. It is not a bound on
`||K_variable-K_product||`. Two further reasons keep it from being
that bound. The guard `|sigma_cold|<1/20` is stated on the cold chart,
not on the whole ball `chi<=84`. And the quantum shape response is a
nonlocal function of the Dirac operator, not the local Laplacian.

Classification of this substitution: uncontrolled approximation.
The failed argument is preserved above and is not used as an input
to any later gate.

## 4. Nonlocal kernel difference

The missing derivative is the one named by the parameter-form note:
three resolvents of the actual paired, common-subtracted two-vertex
bubble, with a mass-potential insertion, the differentiated seagulls,
and the same common finite parts. Its first-jet cancellation and the
matching to the Weyl logarithm at the same scale are not proved here.

For a constant complex background shift `v`, the parameter-form factor

    R_v=exp(-4v)(1+t x)/(1+exp(-2v) t x)

has derivative `-4` at `v=t=x=0`. The same value is recovered by
differentiating in the positive coordinate `u=exp(-2v)` and
multiplying by `du/dv`. That derivative is the constant-background
factor. Replacing a constant `v` by a function `v(Y)` inside the
integrand does not estimate the variable-coefficient Dirac operator.
The constant-background derivative is therefore not the variable
three-vertex or the common contact.

The flat positive remainder at constant `v` likewise stays the
statement already proved: a pointwise replacement is not a
variable-coefficient theorem. No term in `E_bubble` is set to zero.

Classification: unresolved proof.

## 5. Prescribed-source mass and temperature map

For a prescribed spherical source, the retained solution defines

    s_R=[G rho(R)/R - kappa a(R)]/D_h,
    f(R)=-[2 kappa R s_R]/R,
    delta mu_infinity=kappa R s_R+G integral_R^infinity rho dr.

The integral that would correct `f(R)` runs from `R` to `R`, so

    f(R)+2 kappa s_R=0.

This is the horizon lapse constraint for that prescribed source.
Differentiate `(s_R, delta mu_infinity)` in the prescribed data
`(a(R), integral rho)`, holding `rho(R)`, `kappa`, `R`, `G`, and
`D_h` fixed. The Jacobian determinant is

    -kappa G/D_h.

It is nonzero for positive `kappa`, `G`, and `D_h`. The
prescribed-source map is therefore locally one-to-one. This is the
finite-dimensional mass and temperature system on a source that does
not depend on the unknown state.

Fixed temperature in the retained geometry proof also gives `H(R)=0`
for that same prescribed source, and therefore `xi(R)=s_R`. That
identity is not re-derived from the product kernel.

## 6. Coupled seam

The quantum source `rho, p` is the stress of the variable-radius
state. That state depends on `E_bubble` and on the exterior values
sampled by the nonlocal kernel. Section 5 solves the system only
after `rho` and `p` are prescribed. Substituting the product response
for that stress changes the source that the map was given, so the
output of section 5 is not a fixed point.

No boundary datum at `chi=L` repairs this. The product operator
samples both sides of the seam, and the exterior continuation changes
the interior operator. The classical Jacobian does not invert that
coupling.

Classification of the coupled seam: unresolved proof. The
prescribed-source solution keeps its original scope. It is not a
solution of the quantum seam.

No replacement area target was written before this test. The area
budget is not evaluated. No numerical value is assigned to it.

## 7. Constant-background bubble derivative

The constant-background integrand factor of the parameter-form
remainder is `exp(2v) f(exp(-2v) x)`, with

    f(x)=log(1+x)-x+x^2/2.

Differentiation at `v=0` gives

    2 log(1+x)-2x+x^2(1-x)/(1+x).

The same expression is the chain rule
`2 f(x)+f'(x)(-2x)` with `f'(x)=x^2/(1+x)`. At `x=1` the value is
`2 log(2)-2`. The integral `log(2)=integral_0^1 dt/(1+t)` and the
decreasing Riemann upper sum `5/6` give

    2 log(2)-2 < 2*(5/6)-2=-1/3<0.

The constant-background derivative is therefore not the zero function.
It is still not the variable Dirac three-vertex. The higher weighted
nonlinear map is not this derivative.

## 8. Nonlinear map, reduced branch, and observable errors

A reduced branch in the sense of the covariant reduced-inverse note
needs an inverse constant `M_s`, a nonlinear derivative bound `L_s`,
and a source defect `epsilon_s` in the matched weighted spaces, with
the certificate `2 M_s^2 L_s epsilon_s<1` inside the admitted phase.
The symbol `2 M_s^2 L_s epsilon_s` still depends on all three names.
None of the three is set equal to `1.3e8` or to zero. The product
bound is a linear weighted estimate on complete H2. It is not a `C1`
bound for the renormalized determinant in those spaces.

The reduced branch is not constructed. Physical observable errors are
not bounded. An action-value remainder is not that error theorem.

Classification of each of these three items: unresolved proof.
None is set to zero, and none is replaced by the constant-background
derivative.

## 9. Signed source condition, decided separately

The fast outgoing condition is the scalar statement

    m_num(80 pi^2 f_eff)=0.

Necessity and sufficiency for that fixed linear problem are the
retained moment theorem. The product operator is the two-by-two
multiplier `M(p)` of `P2`, with

    det M(p)=-p(4 b^2 p+a g).

Its `p` derivative at the unit positive point `(p,b,a,g)=(1,1,1,1)`
is `-9`, so this determinant is not a constant scalar moment. The
fast scattering ratio `(1+i a_mix)/(1-i a_mix)` has modulus one and
equals `i` at `a_mix=1`. That phase law is not an evaluation of
`m_num`.

The product norm `1.3e8` therefore supplies no sign for `m_num`.
The retained signed data still do not enclose the moment away from
zero and do not prove that it vanishes. The actual state contribution
and both resolvent pieces in that moment remain unenclosed.

Decision: the original signed source condition is undecided.
Classification: unavailable physical input. `m_num` is not set to zero.
This decision does not use the kernel substitution in section 3.

## 10. Transfer and the physical gates

The finite-pulse work theorem remains a statement about its static
reference. Transfer onto a matched variable-radius exterior requires
the seam solution in section 6. That solution is not available.
Driven stress, radius, controller, capture, tail, and spatial area
are the S074 through S078 obligations. They are not calculated on a
stand-in background.

S139.1 and S139.2 ask for the common exterior prescription and the
controlled exterior stress. This attempt does not provide either one.
S074.2 through S074.6 ask for every term at the declared order, a
remainder in the observables that feed the geometry gate, and a range
or an identified input for every unknown coefficient. The unresolved
bubble derivative leaves that pass unsatisfied. S075 through S079
stay blocked on that pass. The earlier obstruction audit remains the
S080 outcome. No satisfying row of the completion plan is met.

## 11. Checker

The focused checker is
[stationary_variable_radius_seam.py](../checks/stationary_variable_radius_seam.py).
It audits the conformal Laplacian, both gap enclosures, the arithmetic
of the rejected substitution, the resolvent identity, the prescribed-source
Jacobian, the constant-background derivatives, and the separation
between `M(p)` and the fast phase. Analytic continuation, spectral
domains, the variable-coefficient bubble, and the coupled fixed point
remain ordinary arguments. The checker does not assign a bound to
`E_bubble` and does not evaluate the area budget.

Remaining to-do list: estimate `E_bubble` in `X_L^2` and solve the
coupled mass and temperature fixed point with that estimate included;
prove the variable Dirac three-vertex, the higher weighted nonlinear
map, the reduced branch, and the physical observable errors without
setting them to zero; enclose the signed moment, or prove it vanishes,
from its own data; then match the exterior, transfer the finite-pulse
theorem, and complete the driven stress, radius, controller, capture,
tail, spatial area, and S074-S080 gates.
