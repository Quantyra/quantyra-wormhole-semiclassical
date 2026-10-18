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

Classification of the full nonlocal derivative, before the estimate
below: unresolved as a three-resolvent theorem. The next section
estimates the part that the retained a2 bound actually controls, and
records why that part does not close the weighted-jet norm.

## 5. Low-jet embedding and the retained a2 derivative

Let `D` be the radial hyperbolic Laplacian
`u''+coth(chi) u'`. The divergence identity

    (sinh(chi) u')'/sinh(chi)=D u

integrates, for a regular bolt value `u'(0)=0`, to

    u'(chi)=(1/sinh chi) integral_0^chi sinh(t) D u(t) dt.

Hence `|u'(chi)|<=||D u||_infinity tanh(chi/2)<||D u||_infinity`.
The weight in the transverse jet collapses exactly:

    coth(chi) tanh(chi/2)=1-1/(2 cosh(chi/2)^2)<1,

so `|coth(chi) u'|<||D u||_infinity`. Therefore
`|u''|<2||D u||_infinity`. If

    ||u||_infinity+||D u||_infinity<=eta,

every jet in the weighted norm `X_L^2` is at most `2 eta` before the
factor `1/cosh(chi)<=1`. The radial low-jet ball embeds into `X_L^2`
with constant 2. This embedding does not run in the other direction.

On `0<=chi<=84`, `cosh(84)>exp(84)/2`. Since `exp(1)>27/10`,

    cosh(84)>(27/10)^84/2>10^35.

A field of weighted norm 1 may therefore have supremum above `10^35`.
The a2 hypothesis in the parameter-form note requires
`||w||_infinity+||Delta w||_infinity<=1/100`. The unit ball of `X_84^2`
is not contained in that hypothesis ball. The ratio of the two sup-norm
allowances exceeds `10^37`.

Inside the hypothesis ball, the retained a2 estimate is

    ||Dr^(-1/2)(partial_w Q_a2) Dr^(-1/2)||<1600 eta.

The ledger `9|c2|/C<=15/(2 pi)` and `[15/(2 pi)](200 pi+1)<1503`,
rounded up to 1600, is used unchanged and covers only this contact.
At `eta=1/100` the relative factor is 16, which is not a contraction.
The factor falls below 1 only for `eta<1/1600`. On that smaller ball
the contraction is a statement about the a2 quadratic form in the `Dr`
topology. It is not a bound on `||E_bubble||` from `X_L^2` to `X_L^2`.

The a2 third-variation matrix, evaluated at vanishing curvature and
vanishing `Delta w`, is diagonal with entries proportional to
`(-2,0)`. Its determinant vanishes, so it has a kernel direction.
A scalar three-resolvent trace does not have that matrix shape. The
constant 1600 is not substituted for `E_bubble`, and `E_bubble` is not
set to zero.

## 6. Failed fixed-point composition

The classical seam map of the next section takes the two scalars
`(a(R), integral rho)`. The a2 bound returns a relative quadratic form,
not those two scalars. Composing the relative factor 1600 with the
product inverse norm `1.3e8` is a formal product:

    1600*1.3e8/100=2080000000,
    1600*1.3e8/1600=130000000.

At the contraction threshold the formal product equals the full product
norm, so it does not produce a smaller seam correction. On the
hypothesis ball it is larger still. This composition is an uncontrolled
approximation to `||E_bubble||_{X_L^2}`. It is preserved as a failed
argument and is not used as a Lipschitz constant for the seam map.

Consequently the coupled mass and temperature fixed point is not
solved. The prescribed-source Jacobian below remains nonzero, and the
missing input is a seam-source Lipschitz bound for the three-resolvent
insertion.

Classification of `E_bubble` in `X_L^2`: unresolved proof, with the
embedding, the factor 16, and the formal products above recorded.
Classification of the coupled fixed point: unresolved proof.

## 7. Prescribed-source mass and temperature map

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

## 8. Coupled seam

The quantum source `rho, p` is the stress of the variable-radius
state. That state depends on `E_bubble` and on the exterior values
sampled by the nonlocal kernel. Section 7 solves the system only
after `rho` and `p` are prescribed. Substituting the product response
for that stress changes the source that the map was given, so the
output of section 7 is not a fixed point. Section 6 shows that the
retained a2 factor does not supply the missing source Lipschitz bound.

No boundary datum at `chi=L` repairs this. The product operator
samples both sides of the seam, and the exterior continuation changes
the interior operator. The classical Jacobian does not invert that
coupling.

Classification of the coupled seam: unresolved proof. The
prescribed-source solution keeps its original scope. It is not a
solution of the quantum seam.

No replacement area target was written before this test. The area
budget is not evaluated. No numerical value is assigned to it.

## 9. Constant-background bubble derivative

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

## 10. Nonlinear map, reduced branch, and observable errors

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

## 11. Signed source condition, decided separately

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

## 12. Transfer and the physical gates

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

## 13. Three-resolvent candidate inside the classical Jacobian

Let `A(s)=D+s V` for an invertible matrix `D` and a matrix insertion
`V`. Differentiating the squared resolvent gives

    d/ds (A^-1 A^-1) = -A^-1 V A^-1 A^-1 - A^-1 A^-1 V A^-1.

Each term contains three resolvents. This is the constant-matrix
insertion. A radius function `v(Y)` replaces `V` by multiplication by
a noncommuting function, and the same algebra does not estimate that
operator.

The only scalar mass derivative available from the retained heat
coefficient is

    d/dM [1/(120 pi M^2)] = -1/(60 pi M^3).

The guards `q=2000^4`, `x>7.9e9`, `|d|<1e-12` and `C=1+6x+d` give
`C>6x>6*7.9e9`, hence

    33755/100 < (1+2000^4)/(6*7.9e9) < 33756/100.

Thus `M_1^2<33756/100` and `M_1<19`. They do not give `M_1>18`.
At `x=10^{10}` and `d=0`, `C=1+6*10^{10}` and

    (1632/100)^2 < M_1^2 < (1633/100)^2,

so `16.32<M_1<16.33`. With `333/106<pi<22/7`, the formal product of
this derivative's absolute value with `1.3e8` satisfies

    158 < product < 159.

That value is not inside `(100,124)`. It is not a Lipschitz constant
below 1.

The candidate is also not the seam source. The classical Jacobian
takes `(a(R), integral rho)` to `(s_R, delta mu)`. A scalar
coefficient derivative does not supply those two values. Using the
formal product as `||E_bubble||_{X_L^2}` is an uncontrolled
approximation. The three-resolvent insertion in `X_L^2` remains an
unresolved proof. It is not set equal to the heat candidate and it is
not set to zero. The coupled fixed point remains unsolved.

## 14. Variable Dirac three-vertex

The constant-background remainder bound `|P_{M,v}-P_M|<=.05 P_M` for
`|v|<=.01` stays a constant-shift statement. The variable three-vertex
needs the commutator of the Dirac operator with `v(Y)`, a gain-two
remainder, and the common finite parts at one scale.

The sequence `f_N(Y)=N^-5 sin(N Y)` has fifth derivative `cos(N Y)` and
sixth derivative `-N sin(N Y)`. Its `C5` amplitude is 1, while its next
derivative has amplitude `N`. A bound on five derivatives does not
bound the sixth. That sequence is a counterexample to reading a `C5`
label as the higher weighted control the three-vertex remainder asks
for. It is not a statement that the physical vertex vanishes.

Classification of the variable Dirac three-vertex and its common
contact: unresolved proof. Neither is set to zero, and neither is
replaced by the constant-background factor or by `f_N`.

## 15. Higher weighted nonlinear map

The reduced-branch certificate is `2 M_s^2 L_s epsilon_s<1`.
Substituting the only numerical inverse constant on hand,
`M_s=1.3e8`, gives

    2 M_s^2=33800000000000000.

The product `L_s epsilon_s` remains in the certificate. No retained
estimate supplies values that make the product smaller than the
reciprocal of that square. Section 14 also shows that a `C5` bound
does not control the next derivative, so the weighted product bound
does not yield the higher map.

The reduced branch is not constructed. Physical observable errors are
not bounded. Classification of each: unresolved proof. None of
`M_s`, `L_s`, or `epsilon_s` is set to zero.

The signed moment remains the unavailable input recorded in section 11.
The heat candidate does not evaluate it. The exterior match, the
finite-pulse transfer, and the driven stress, radius, controller,
capture, tail, and spatial area stay blocked. No pass condition of
S139.1, S139.2, or S074.2 through S079.6 holds. The S080 outcome
remains the earlier obstruction audit. No satisfying row is met.
The area budget is not evaluated.

## 16. Noncommuting kernel in the weighted norm

The three-resolvent source with a radial multiplier `phi` has the
integral kernel

    Lambda_phi(x,x')
      = integral K(d(x,y)) phi(y) K(d(y,z)) K(d(z,x')) dvol(y,z)
      + integral K(d(x,y)) K(d(y,z)) phi(z) K(d(z,x')) dvol(y,z).

The two insertions sit at different points. Split `phi(y)=phi(x)+(phi(y)-phi(x))`.
The first piece commutes out of the integral. The second is the
noncommuting gap.

On the ball `chi<=84`, the geodesic distance is at most `2*84=168`.
The weighted norm controls `|phi'|/cosh(chi)`. Hence

    |phi(y)-phi(x)| <= 168 ||phi||_X cosh(84).

The retained comparison `cosh(84)>10^35` therefore gives a gap larger
than `10^37 ||phi||_X`. That constant is not below 1.

The same weight shows that `X_L^2` is not an algebra. The radial
profile `cosh(chi)` has all four weighted jets bounded by 1, so its
norm is 1. Its square is `cosh(chi)^2`, and dividing by `cosh(chi)`
returns `cosh(chi)`. The weighted norm of the product is therefore
greater than `10^35`. Pulling `||phi||_X` out of `Lambda_phi` without
this factor is false. The corrected factor is not a Lipschitz constant
below 1.

The horizon observable is the bolt value. There `cosh(0)=1`, so
`|u(0)|<=||u||_X`. The evaluation constant is 1, which is not strictly
below 1. One composition with the retained inverse norm replaces 1 by
`1.3e8`, which is not below 1. The contraction threshold for a source
derivative is `1/1.3e8`. The multiplication factor `10^35` exceeds that
threshold by more than `10^43`.

The rho-integral observable in the classical Jacobian is `integral rho dr`.
Its change-of-variables factor `dr/dchi` is not a retained finite
constant on the exterior ray. The hyperbolic comparison functional
`integral_0^84 u(chi) sinh(chi) dchi`, evaluated on the same unit-norm
profile `cosh(chi)`, equals `sinh(84)^2/2` and exceeds `10^69`. A unit
ball of `X_L^2` therefore does not force this integral below 1. The
physical `dr` integral is not identified with that comparison. Its
missing Jacobian is an unresolved proof.

No ball inside the classical Jacobian receives a Lipschitz constant
below 1 for the pair `(horizon value, rho integral)`. The Jacobian
determinant `-kappa G/D_h` remains nonzero, and the source increments
are not placed in its contracting neighborhood. The tail discarded by
zero extension is also unresolved: the ball restriction of `K` does not
reproduce `Lambda_phi`.

Classification: uncontrolled approximation for any bound that omits the
`cosh(84)` weight; unresolved proof for the physical `dr` Jacobian and
the exterior tail. Neither obstruction is removed by setting a
coefficient to zero.

## 17. Exterior tail and the radial factor

The coordinate definition `v=log tanh(chi/2)` has derivative
`dv/dchi=1/sinh(chi)`. The retained radial formula has
`dv/dr=kappa/F0`. Therefore

    dr/dchi=F0/(kappa sinh(chi)).

On the region where the retained upper bound is `F0<=C s^2`,

    dr/dchi<=(B/gamma) s^2/sinh(chi),

because `kappa=gamma C/B`. The chart inequality `s<=1.008 h cosh(chi)`
rearranges to `cosh(chi)>=s/(1.008 h)`. For `cosh>2`, `sinh>cosh/2`, so

    dr/dchi<(B/gamma) s h (2*1.008)

when that lower bound on the hyperbolic sine is inserted. At `h=gamma`
this is `2.016 B s`.

The guards `x>7.9e9`, `h<1e-30` and `|h/gamma-1|<=12 gamma` admit
`h=gamma=10^{-31}`, `B=1000` and `lambda=10^4`. Then `x=10^{10}`, and
at the cold seam `s=1/100` the cap equals `20.16`, which is not below 1.
The guards do not cap `B`. This attempt therefore does not produce a
Lipschitz constant below 1 for the rho integral. The exact factor is
kept. The cap is not sharpened by setting `B` or `gamma` to zero.

The physical disk ends at `D0=log(coth(h/2))`. For `h<1e-30`,
`coth(h/2)>2/h>2*10^{30}`, and `D0>69`. A source supported at distance
at least 69 from the bolt, with unit `L-infinity` density, meets the
kernel majorant integrated against hyperbolic area. For distance at
least 1 the square-root denominator is at least `9/10`, and `sinh`
is at most half an exponential. The resulting tail is below `10^{-20}`.
That is below 1. It is a bound on one kernel leg acting on a unit
density, not on the three-resolvent composition.

The composition applies that leg to the intermediate field on
`chi>D0`. The retained theorem bounds that field in the weighted norm,
not in `L-infinity`. No `L-infinity` bound is selected here, and the
weight comparison recorded earlier is not reused as this Lipschitz
constant. The boundary layer, where the distance to the fictitious
exterior is short, is likewise not bounded by the long-distance tail.
Classification of the composition tail: unresolved proof. The one-leg
unit-density tail below `10^{-20}` remains available and is not set
equal to the composition.

No Lipschitz constant below 1 is obtained for the pair of seam
observables. The rho-integral cap exceeds 1 on an admissible guard
point, and the composition tail has no `L-infinity` input.

## 18. Exterior L-infinity bound and the guard set for B

On `chi<=84` the retained weighted estimate gives
`|U|/cosh(chi)<=1.3e8 ||F/omega||_infinity`. The fictitious exterior
begins at `D0>69`, which is inside that ball. Therefore every value
there satisfies

    |U(chi)|<=1.3e8 cosh(chi) ||F/omega||_infinity.

The series through degree 7 is a lower bound for `exp(1)`. Its 69th
power exceeds `10^29`, and `cosh(69)>exp(69)/2`, so `cosh(69)>10^29`.
The resulting `L-infinity` allowance on `chi>=69` exceeds `10^37` times
the weighted source. It is not below 1. This is the bound supplied by
the retained weighted theorem. It does not become a composition
Lipschitz constant below 1.

The factor `2.016 B s` is below 1 only if `B<1/(2.016 s)`. At the cold
seam `s=1/100` this would require `B<1000/20.16`. The guards do not
impose that restriction. For every `B>=1` the choice `h=gamma=10^{-31}`
and `lambda=8*10^9` has `x=lambda B^2>=8*10^9>7.9e9` and
`|h/gamma-1|=0<=12 gamma`. At `B=50` and `s=1/100`,
`2.016 B s=1.008>1`. The admissible set is unbounded in `B`, so it
contains no uniform cap below 1. This obstruction is read off the
guards. It does not reuse a single previously evaluated cap, and it
does not set `B` to zero.

## 19. Variable Dirac three-vertex, continued

The constant-background theorem assumes a shift with
`|v|<=1/100` and no spatial derivative. The family
`v_N=sin(N Y)/100` obeys `|v_N|<=1/100` for every `N`, while

    dv_N/dY=N cos(N Y)/100.

At `N=100` the derivative amplitude is 1. The constant-shift derivative
contains no factor of `N`. Replacing `v_N` by a constant inside the
three-vertex drops this term. That replacement is an uncontrolled
approximation. The variable three-vertex remains an unresolved proof.
The derivative term is not set to zero, and `N` is not set to zero.

The higher weighted map needs a derivative bound `L_s`. The same family
shows that the constant-shift hypothesis does not provide one: the
derivative amplitude grows with `N` while the supremum stays inside the
hypothesis. `L_s` and `epsilon_s` stay unfixed. The reduced branch is
not constructed, and the physical error bound is not obtained.

## 20. Higher weighted map, reduced branch, and observable errors

The product estimate controls the weighted jet norm `X_L^2`, whose
jets stop at order two. Define the third-jet seminorm by the weighted
supremum of `|u'''|`. For `N>=1` set

    f_N(chi)=cos(N chi)/N^2.

Then `f_N''(0)=-1`, so the weighted second jet at the bolt is 1 and
`||f_N||_X>=1`. The third derivative is `N sin(N chi)`. Its amplitude
is `N`. The identity `|sin u|<=|u|` follows from `u-sin u` vanishing
at 0 with nonnegative derivative `1-cos u`, and
`|sin(N chi)/(N sinh chi)|<=chi/sinh(chi)<=1` keeps the transverse
jet inside the same norm. Thus these fields lie in the unit ball of a
fixed multiple of `X_L^2`, while the third-jet seminorm is at least a
positive multiple of `N`.

No finite `L_s` satisfies
`third seminorm <= L_s ||f||_X` for every such field. Taking
`N` larger than that finite candidate violates it. `epsilon_s` is not
set to zero: each `f_N` has weighted second jet 1. The Newton
certificate `2 M_s^2 L_s epsilon_s<1` therefore has no finite `L_s`
on this ball. The higher weighted map is not obtained. This comparison
uses the weighted jet norm directly. It does not pass through a `C5`
bound.

The reduced branch is not constructed. Physical observable errors are
not bounded. Classification of each: unresolved proof. Neither `L_s`
nor `epsilon_s` is set to zero.

## 21. Signed-moment enclosure

The retained fast expansion has `rho_mix=8*10^{-12}` and

    error in m_num <= 2 ||g|| ||h_R|| rho_mix^(m+1)/(1-rho_mix).

At `m=0` the numerical factor is `16/(10^{12}-8)`, which is below
`2*10^{-11}`. The norms `||g||` and `||h_R||` are not retained
numbers. Their product with this factor is not an interval about
`m_num`. The enclosure does not exclude zero and does not prove zero.
`m_num` is not set to zero. Classification: unavailable physical
input for those two norms.

The higher-map obstruction does not evaluate the fast pairing. The
exterior match, the finite-pulse transfer, and the driven stress,
radius, controller, capture, tail, and spatial area stay blocked.
No pass condition of S139.1, S139.2, or S074.2 through S079.6 holds.
The area budget is not evaluated.

## 22. Moment norms

Equation (5) of the scalar phase note gives

    ||g||^2/c0=||Im L||<=||L||<4*10^{-12},

with `c0=k(E+4 lambda)`. Hence `||g||<2*10^{-6} sqrt(c0)`. The factor
`4*10^{-12}` is exactly `(2*10^{-6})^2`. The quantity `c0` is not a
retained number, so this is not a numerical interval for `||g||`.

The second factor is `h_R=V*R4_PV p` with `p=80 pi^2 f_eff`. The same
note states that the retained source data do not evaluate (27), and
that the series remainder does not certify errors in `h_R` or the
source. No retained operator norm converts the unknown `f_eff` into a
numerical `||h_R||`.

The truncation product therefore stays

    [16/(10^{12}-8)] ||g|| ||h_R||
      < [16/(10^{12}-8)] (2*10^{-6}) sqrt(c0) ||h_R||.

Both `sqrt(c0)` and `||h_R||` remain. The product is not an interval
about `m_num`. It does not exclude zero and it does not prove zero.
Classification: unavailable physical input. `m_num` is not set to zero.

## 23. Exterior matching, with the moment left open

The physical outer endpoint is defined by `tanh(D0/2)=exp(-h)`, so
`D0` is finite for `h>0`. The product inverse is the complete-H2
resolvent, whose comparison coordinate runs past every finite `D0`.
Its value at `chi=D0` is an interior sample of that resolvent, not a
matched outer datum.

The only retained pointwise allowance on `chi>=69` is
`|U|<=1.3e8 cosh(chi)` times the weighted source. That allowance
exceeds `10^37`. It does not show that the interior sample equals the
required outer datum within a tolerance below 1. The open factors
`sqrt(c0)` and `||h_R||` are not used to cancel this difference, and
`m_num` is not set to zero in their place.

Exterior and state matching therefore remains an unresolved proof.
S139.1 stays open. The finite-pulse transfer still requires this match.
Driven stress, radius, controller, capture, tail, and spatial area stay
blocked. The area budget is not evaluated.

## 24. Search of the retained notes

The scalar phase note defines

    E=60 pi/G+11 lambda/3,  k^2=E-9 lambda/4,  c0=k(E+4 lambda),

and `h_R=V*R4_PV p` with `p=80 pi^2 f_eff`. It gives
`||g||^2/c0<4*10^{-12}` and `||L||<4*10^{-12}`. A line search of that
note finds the formula `m_num(p)=m0+<g_R,J h_R>`, no assignment
`m_num=0`, no assignment of `c0` to a numeral, and no numeric bound on
`||h_R||`. It contains no `chi=D0` trace. The outgoing-inverse note
states `k_rad^2=E-9 lambda/4` and

    .999*2000^{-10}<G lambda<2000^{-10},
    k_rad/sqrt(lambda)>4*10^{17},  E>180/G.

That note contains no `c0` and no `chi=D0` trace.
The phase note states the fast condition `m_num(80 pi^2 f_eff)=0` and
the formula `m_num(p)=m0+<g_R,J h_R>`. It does not assign `m_num=0`.

The product `G lambda=(1999/2000)*2000^{-10}` lies in the open window
above. At `G=1` and at `G=4` the identities fix `E` and `k^2`. Because
`pi>3`, both points meet `k/sqrt(lambda)>4*10^{17}` and `E>180/G`.
The squared values of `c0` stand in the ratio 64, so the positive
values stand in the ratio 8. `c0` is not one number. `m_num` is not
set to zero.

The explicit enclosures are in section 27. `||h_R||` and the
`chi=D0` trace remain unavailable physical inputs.

## 25. Why S138 and the later gates stay blocked

S138 transfers the finite-pulse work theorem from its static reference
onto a matched quantum exterior. That exterior is the unmatched trace
at `chi=D0` together with the unenclosed moment. The transfer is not
performed.

S076 stress and radius need an accepted geometry. S075's controller
and drive need the S074 pass. Capture needs the matched drive. The
exterior tail and the `dr` Jacobian remain unresolved. The spatial
area budget has no replacement target, so it is not evaluated.
S074.2 through S074.6 and S075.1 through S079.6 stay unchecked. The
S080 outcome remains the earlier obstruction audit. No satisfying row
is met.

## 26. Checker

The focused checker is
[stationary_variable_radius_seam.py](../checks/stationary_variable_radius_seam.py).
It audits the spectral identities, the Machin enclosures of `c0`, the
line search of the phase and outgoing notes, and the scaling of the
moment error with `||h_R||`. It does not replace `m_num` by zero, and
it does not evaluate the area budget.

## 27. Explicit c0 values and the S074.2 gate

The comparison is the phase identity `c0=k(E+4 lambda)` at the two
window points. The hardest step is the square root of `k^2`. The
observable is `c0`. The test is Machin's enclosure of `pi` together
with a line search for a numeric `||h_R||` and a `chi=D0` trace below 1.

On `G lambda=(1999/2000)*2000^{-10}`,

    pi=16 atan(1/5)-4 atan(1/239),

and the identities `E=60 pi/G+11 lambda/3` and `k^2=E-9 lambda/4` give

    258792/100 < c0 < 258793/100    at G=1,
    32349/100 < c0 < 32350/100      at G=4.

The positive ratio is 8. A line search of the phase note finds
`h_R=V*R4_PV p` and no comparison of `||h_R||` with a numeral. The
outgoing note does not contain `h_R`. Neither note contains `chi=D0`,
so neither states a trace below 1 at that endpoint.

S074.2 asks for a frozen action, boundary ensemble, gauge, canonical
fields, perturbative order, and parameter domain, and for a map of the
closest source equations. The phase identities are mapped. They do not
freeze the domain: `c0` meets both enclosures, and `||h_R||` is not an
identified number. The exterior trace is absent. S074.2 stays
unchecked. `m_num` is not set to zero.

## 28. Per-gate outcomes after the unmet S074.2 dependency

Each gate below is blocked by section 27. The parameter domain is not
one point, `||h_R||` is not a number, and no `chi=D0` trace below 1 is
on record. No named derivation is substituted for that missing freeze.
No pass box is checked. `m_num` is not set to zero.

S074.3 asks for the cubic and quartic vertices, including lapse,
shifts, Gauss constraints, induced contacts, and boundary terms, with
angular sectors and either fermions and higher modes or a bound on
their omission. Those vertices need the frozen domain from S074.2.
Classification: unavailable physical input. The vertices are not
derived and are not set to zero. S074.3 stays unchecked.

S074.4 asks for one regulator, subtraction, and renormalization
condition, and for the finite couplings that need physical input. It
needs the S074.3 vertices. Classification: unavailable physical input
inherited from S074.2. No regulator is fixed. S074.4 stays unchecked.

S074.5 asks for a local/nonlocal match without double counting, a
regulator-independent remainder, and checked current and stress
identities. It needs the S074.4 regulator. Classification: unresolved
proof. No remainder is bounded. S074.5 stays unchecked.

S074.6 asks for a direct variation of that derivation and an
independent check. The derivation is not available. The shipped checker
audits the `c0` enclosures and the note search. It does not verify
cubic vertices. Classification: unresolved proof. S074.6 stays
unchecked.

S075.1 asks for both boundary clocks, frames, ensembles, and vertex
normalization, and for the relation of the free angle to the
interacting coupling. That relation needs the S074 action. Classification:
unavailable physical input. S075.1 stays unchecked.

S075.2 asks for the joint dressed initial state, its charge and spin
branching, and a preparation-error norm. It needs S075.1 and a numeric
`||h_R||`. Classification: unavailable physical input. The state is not
constructed. S075.2 stays unchecked.

S075.3 asks for the boundary Hamiltonian of the returned switch and for
gauge, charge, and stress identities that include the drive. It needs
S075.1. Classification: unresolved proof. S075.3 stays unchecked.

S075.4 asks for mixed current and noise correlations over the whole
switch, or a uniform remainder. It needs S075.3. Classification:
unresolved proof. S075.4 stays unchecked.

S075.5 asks for the ordered massive-mode work estimate on this state,
with endpoint dressing removed. It needs S075.2. The S138 theorem stays
on its static reference. Classification: unresolved proof. S075.5 stays
unchecked.

S075.6 asks for continuum work and charge balance and for the inputs
S076 requires. It needs S075.3 through S075.5. Classification:
unresolved proof. S075.6 stays unchecked.

S076.1 asks for the observable ledger that carries action, state, and
source errors into energy, clock, and area, frozen before any
acceptance run. Those errors are not numbers. Classification:
unavailable physical input. S076.1 stays unchecked.

S076.2 asks for the renormalized stress, currents, and gravitational
sources from that same action and state, with Ward identities. The
action and state are not available. Classification: unavailable
physical input. S076.2 stays unchecked.

S076.3 asks for the coupled response through the switch and the
returns, with backreaction. It needs S076.2. Classification:
unresolved proof. S076.3 stays unchecked.

S076.4 asks for reduction, tail, and quantum-source bounds on one
common parameter domain, including the extremal area coefficient. The
S074.2 domain is two disjoint `c0` intervals, not one domain.
Classification: unavailable physical input. S076.4 stays unchecked.

S076.5 asks for a nonempty budget window, or a replacement target
written before the test, with positive clock and area. No replacement
target was written. The area budget is not evaluated. Classification:
unresolved proof. S076.5 stays unchecked.

S076.6 asks for an independent check of whether physical opening is
achieved. No opening is claimed. Classification: unresolved proof.
S076.6 stays unchecked.

S077.1 asks for perturbation sectors, a gauge-invariant norm, and a
survival criterion before any stability calculation. No accepted S076
geometry exists. Classification: unavailable physical input. S077.1
stays unchecked.

S077.2 asks for the constrained perturbation operator and controlled
spectral tails. It needs S077.1. Classification: unresolved proof.
S077.2 stays unchecked.

S077.3 asks for a nonlinear disturbance radius inside the opening
margins. Those margins are not established. Classification: unresolved
proof. S077.3 stays unchecked.

S077.4 asks for fluctuation bounds of smeared observables in the
interacting state. That state is the missing S075.2 state.
Classification: unavailable physical input. S077.4 stays unchecked.

S077.5 asks for a weakest-mode test and an exact time scope. It needs
S077.2. Classification: unresolved proof. S077.5 stays unchecked.

S078.1 asks for sender, receiver, clocks, packet profile, energy,
charge and spin, the packet time interval, and an arrival criterion.
The clocks are the missing S075.1 clocks. Classification: unavailable
physical input. S078.1 stays unchecked.

S078.2 asks for propagation on the matched geometry and a nonzero
reception bound. The `chi=D0` trace is absent, so the geometry is not
matched. Classification: unresolved proof. S078.2 stays unchecked.

S078.3 asks for packet stress in the coupled equations and a nonzero
energy range that keeps clock and area positive until reception. It
needs S076.3. Classification: unresolved proof. S078.3 stays unchecked.

S078.4 asks for causal attribution of reception to throat transport.
It needs S078.2. Classification: unresolved proof. S078.4 stays
unchecked.

S078.5 asks for the achieved energy and travel-time scope. It needs
S078.3. No shortcut claim is made. Classification: unresolved proof.
S078.5 stays unchecked.

S079.1 stays checked. It is the earlier inventory.

S079.2 asks for one joint target state and one preparation-error norm
from the physical tolerances. Those tolerances are the open S075, S077,
and S078 gates. Classification: unavailable physical input. S079.2
stays unchecked.

S079.3 asks for a finite-resource preparation protocol that reaches the
opening basin. It needs S079.2. Classification: unresolved proof.
S079.3 stays unchecked.

S079.4 asks for that protocol to be evolved through switching and
reception, with errors carried into S076 through S078. Those stories
are open. Classification: unresolved proof. S079.4 stays unchecked.

S079.5 records that the initial bridge remains supplied. Creation from
disconnected geometry is not constructed. The topology mechanism
remains an unresolved proof. S079.5 stays unchecked.

S079.6 asks for separate preparation and creation verdicts. Preparation
is unachieved. Creation is unachieved. The unresolved requirement is
the unfrozen S074.2 domain, the missing numeric `||h_R||`, and the
missing `chi=D0` trace below 1. S079.6 stays unchecked.

S080 remains the earlier obstruction audit. No satisfying row is met.
The E003 objective is unmet.

## 29. One frozen point, attempted for S074.2

The comparison is `||h_R||` and the weighted sample at `chi=D0` on one
frozen point. The hardest step is passing from `||L||<4*10^{-12}` to
`||V* R4_PV p||`. The observables are that norm and the sample. The
test is a two-by-two source that meets the `L` scale while the image
exceeds 1, together with the disk endpoint at `h=10^{-31}`.

Freeze

    G=1,  lambda=(1999/2000)*2000^{-10},  h=10^{-31}.

Section 27 gives `258792/100<c0<258793/100` at this `G`. The phase
definition remains `h_R=V*R4_PV p` with `p=80 pi^2 f_eff`, and
`||L||<4*10^{-12}` for `L=V*R+ V/(E+4 lambda)`. The denominator
`E+4 lambda` is fixed by `c0` and `k` at this point. The source `p` is
not. The diagonal example

    V=(1,0),  R=diag(10^{-12},1),  p=(10^{20},0)

has sandwich `V* R V=10^{-12}<4*10^{-12}` and image `V* R p=10^8`. A
bound of the `L` type does not determine `||h_R||`. Substituting the
upper value of `c0` into the moment-error factor leaves the positive
symbol `||h_R||`. `m_num` is not set to zero.

The same point has `D0=log(coth(h/2))=2 atanh(exp(-h))`. The degree-3
exponential polynomial at `h=10^{-31}` lies above `1-h`, the omitted
term is positive, and the resulting ratio exceeds `10^{31}`. Since
`exp(69)<2*10^{30}`, this gives `D0>69`. The retained weighted
allowance on `chi>=69` exceeds `10^{37}` times the weighted source. It
is not a trace below 1 at `chi=D0`.

Classification: unavailable physical input, from this frozen-point
calculation. S074.2 stays unchecked. The later gates in section 28 stay
unchecked. They are not closed by this attempt.

## 30. Signed source, charge, and the finite-pulse transfer

The comparison is the fast condition `m_num(80 pi^2 f_eff)=0` on the
frozen point of section 29, kept separate from the product inverse.
The hardest step is the sign of `m0+<g_R,J h_R>`. The observable is
that moment. The test is the outgoing identity and two pairings of
equal size and opposite sign.

The phase relation is

    A_out=i exp(i delta_full)<psi,p>/c0,
    <psi,p>=m_num/sqrt(1+a_mix^2).

At the frozen point `c0>2587>0`, so `A_out=0` if and only if
`m_num=0`. The truncation bound is a positive multiple of `||h_R||`.
Section 29 leaves that norm unset. With reference value `m0=1` and
pairing size `1`, the moment `m0-<g_R,J h_R>` is `0` and
`m0+<g_R,J h_R>` is `2`. The same size admits both values. The signed
condition is not decided. `m_num` is not set to zero.

The product determinant at the unit point `(p,b,a,g)=(1,1,1,1)` equals
`-5`. That value is not the fast condition `0`, and it is not the full
charge. Exterior matching still requires the sample at `chi=D0` to lie
within `1` of the outer datum. Section 29 places the weighted
allowance above `10^{37}`. The sample is not matched, and the full
charge is not matched.

S138 transfers the finite-pulse theorem onto that matched exterior.
The match tolerance is not met, so the transfer is not performed.
Classification: unavailable physical input for the moment sign and the
trace. No replacement area target is written in this section. The area
budget is not evaluated. No pass box is checked. No satisfying row
is met.

## 31. Reference projection and the three admissible signs

The comparison is the fast condition `m_num(80 pi^2 f_eff)=0`, kept
separate from the product inverse. The hardest step is the cancellation
of `(a_tr-epsilon T)Z` against the reference field. The observable is
the sign of the remaining pairing. The test is the coefficient `1/6`,
one self-adjoint eigenvector, and three remainder vectors of equal norm.

The phase note sets `epsilon=G/(60 pi)` and `E=a_tr/epsilon`. Then

    80 pi^2 epsilon/(8 pi G)=1/6.

On a self-adjoint matrix, an eigenvector of `T` at eigenvalue `E` pairs
to zero against `(a_tr-epsilon T)Z` for every height `Z`. The curvature
remainder in equation (23a) is not fixed by that cancellation. The three
vectors `(1,0)`, `(-1,0)`, and `(0,1)` have squared norm `1` and
pairings `1`, `-1`, and `0` against the first coordinate. The same
absolute norm admits both signs and exact zero.

Equation (23a) does not remove the second term of (27). The phase note
states that the retained data do not evaluate (27) or bound it away
from zero. Section 30's pairings `0` and `2` are both nonnegative, so
they do not exhibit the negative sign. `m_num` is not set to zero. The
fast condition stays undecided.

The unit product determinant remains `-5`. That value is not one of
these three pairings. Section 29 still places the weighted sample at
`chi=D0` above `10^{37}`, so exterior, state, and full charge are not
matched. The S138 transfer is not performed. No replacement area target
is written. The area budget is not evaluated. No pass box is checked.
No satisfying row is met.

Classification: unavailable physical input, from this projection. The
missing input is the unevaluated second term of (27), together with the
unmatched trace at `chi=D0`.

Remaining to-do list: the reference projection coefficient is `1/6`, and
the trace term cancels on an eigenvector, but the curvature remainder
still admits moments `1`, `-1`, and `0` at equal norm. Equation (23a)
does not remove the second term of (27). The fast condition stays
undecided, and `m_num` is not zero. The unit product determinant is
`-5`, not that condition. The weighted sample at `chi=D0` stays above
`10^{37}`, so exterior, state, and full charge are not matched and the
S138 transfer is not performed. No replacement area target was written.
The area budget is not evaluated. S074.2 through S079.6 stay unchecked.
No satisfying row is met.
