# Nonlinear determinant current and susceptibility on a fixed domain

2026-10-15, E003/S132; finite response milestone supporting S074. This
supplies the derivative estimate missing from S087's value-bound test.
It controls the nonlinear current and its susceptibility for prescribed
sources in the existing unequal-return Dirac model. It does not compute
the corrected physical background's interacting determinant or close S074.

## Scope and predeclared acceptance

The fixed-domain source map and coefficient are those of
[the nonlinear-current note](nonlinear-current-determinant.md), equations
(1)-(2). The unresolved observable is identified in
[the S085 note](source-clock-unmatched-inputs.md#nonlocal-interacting-determinant)
and [S087](source-clock-named-inputs.md#norm-on-the-interacting-determinant).
S087 correctly says that a bound on a value cannot be differentiated into
a response bound. Here a derivative operator is constructed and bounded.
No derivative of its scalar envelope is used.

Before running the new checker, fix these numerical tolerances: value
series and derivative/Hessian series tails must satisfy their analytic
bounds plus 5e-11 absolute rounding allowance; central differences with
h=2e-4 and h/2 must agree with the first derivative within 2e-6, and mixed
Hessians and second path derivatives within 2e-5. Symmetry and algebraic
identities have tolerance 5e-11. Inequality tests allow a relative 5e-10
rounding margin. Noncommuting examples must have commutator norm >1e-6.
The equal-return limit is not substituted for the unequal-return source
example. Cutoff convergence in the theorem is an analytic result;
finite refinements are supplementary evidence with no fitted rate.

## Operator response theorem

Let K belong to the complex Schatten ideal S3. Set s=||K||3 and
a=||K||op<1. The logarithm based at zero is

```text
F(K)=log det_3(I+K)=sum_(n>=3) (-1)^(n+1) Tr(K^n)/n.       (1)
```

The condition a<1 is sufficient; s need not be less than one. For
H,J in S3, put R=(I+K)^(-1). Then

```text
DF(K)[H] = Tr(K^2 R H),
|DF(K)[H]| <= s^2 ||H||3/(1-a),                           (2)
D2F(K)[J,H] = Tr(((KJ+JK)R-K^2 R J R)H),
|D2F(K)[J,H]| <= s (2-a) ||J||3 ||H||3/(1-a)^2.          (3)
```

No matrices are assumed to commute or be normal. Every term within the
displayed traces is trace class by Schatten Holder. The two negative
and positive terms in (3) retain their ordering.

Proof. Termwise first differentiation of Tr(K^n)/n gives
Tr(K^(n-1) H) by cyclicity. Schatten Holder, with two factors K in S3
and the remaining factors in operator norm, bounds its absolute value
by s^2 ||H||3 a^(n-3). The resulting geometric series converges
locally uniformly in the S3 norm and sums to (2). Differentiating that
series again gives n-1 ordered terms per degree n. Each contains H,J,
and n-2 copies of K. One K and H,J are placed in S3; the remaining
n-3 copies have operator norm at most a. Summing
sum_(n>=3) (n-1)a^(n-3)=(2-a)/(1-a)^2 proves (3), and local uniform
convergence proves Frechet differentiability. Differentiating K^2 R
using DR[J]=-R J R gives the stated identity. The trace series proves
symmetry in H,J, even though the displayed ordered expression is not
manifestly symmetric.

For a C2 path K(z), with h=||K'(z)||3 and b=||K''(z)||3,

```text
|dF/dz| <= s^2 h/(1-a),
|d2F/dz2| <= s(2-a)h^2/(1-a)^2 + s^2 b/(1-a).            (4)
```

If ||K||3^3<=r<1, replace a,s by r^(1/3) in these sufficient bounds.
In particular |dF/dz|<=r^(2/3)||K'||3/(1-r^(1/3)). A value bound
without the independent derivative norm supplies no such conclusion.

There is a stronger, local alternative outside the operator-norm ball:
if I+K is invertible, any local branch of log det_3 has derivative (2)
and bound s^2 ||R||op ||H||3. To see the extension without continuing
an unsafe logarithm, first use finite-rank approximations and the
finite-dimensional determinant derivative. Invertibility is open;
S3 convergence gives operator-norm convergence and uniformly bounded
inverse operators near K. The ordered trace-class expression converges
by Holder, so the local derivatives converge. The small-norm regime
above already suffices for the source corollary; this alternative
requires an independently established resolvent bound.

## Explicit series and Taylor error certificates

Let F_N retain degrees 3 through N, N>=3. Direct geometric summation
gives

```text
|F-F_N| <= s^3 a^(N-2)/((N+1)(1-a)),
|DF[H]-DF_N[H]| <= s^2 ||H||3 a^(N-2)/(1-a),
|D2F[J,H]-D2F_N[J,H]|
 <= s ||J||3 ||H||3 a^(N-2) [N/(1-a)+a/(1-a)^2].         (5)
```

At a=0 these expressions extend continuously, including zero tails.
On a segment K+tH, 0<=t<=1, suppose the uniform bounds are
||K+tH||op<=a_*<1 and ||K+tH||3<=s_*. Integrating the second
derivative gives a useful first-order Taylor certificate,

```text
|F(K+H)-F(K)-DF(K)[H]|
 <= (1/2) s_* (2-a_*) ||H||3^2/(1-a_*)^2.                 (6)
```

For increasing finite-rank orthogonal projections P_m strongly tending
to identity, P_m K P_m -> K in S3. The same holds for H,J and path
derivatives. On a compact C2 path these convergences are uniform, by
finite nets in the compact images of K,K',K''. Compression contracts
the operator norm. Every finite series term converges in trace norm,
and the uniform tails (5) prove convergence of F,DF,D2F. Thus removing
the fermion regulator preserves these response bounds. No photon
regulator, geometry limit or Lorentzian continuation is removed here.

## Realization as nonlinear current and susceptibility

Fix q, all circle lengths, the common support interval, beta (possibly
infinite), the Dirac operator D0 and its antiperiodic domains. For real
eight-component sources phi, define the existing linear map

```text
A_q(phi)=D0^(-1) i sum_mu gamma_mu sum_A phi_(mu,A) E_A,
c_q=4 sqrt(2) Lbar_beta/(pi sqrt(q)),
||A_q(phi)||3 <= c_q^(1/3) ||phi||_L3.                    (7)
```

The earlier interpolation proof establishes (7) on bounded sources.
Linearity and density extend this map continuously to L3 sources on
the fixed source region. This is an operator extension, not an assertion
that every unbounded multiplier has an independently defined classical
Dirac domain. For bounded sources it agrees with the original product.
For beta=infinity the region has compact temporal support as before.
No maximum return length replaces Lbar_beta.

For a=c_q^(1/3)||phi||_L3<1, compose (1) with (7). Its first variation
is a bounded complex-valued functional on real L3. Apply real Lp duality
separately to its real and imaginary parts to define j_NL. Its two real
current components and the bilinear susceptibility satisfy

```text
max(||Re j_NL(phi)||_(3/2), ||Im j_NL(phi)||_(3/2))
 <= c_q ||phi||_3^2/(1-a),
|D2(F o A_q)(phi)[psi,xi]|
 <= c_q ||phi||_3 (2-a) ||psi||_3 ||xi||_3/(1-a)^2.        (8)
```

The full complex-vector L^(3/2) norm is at most twice the first bound,
by the triangle inequality. No equality between the norm of a complex-
valued functional on real sources and that complex-vector norm is used.

The effective action has the opposite sign to F, so its nonlinear
current has the opposite sign and the same norm bound. These are
response bounds with the fermion cutoff removed, rather than just
small nonlinear values. They bound all orders starting at cubic;
the linear and quadratic renormalized responses must still be added.
The det_3 component and its current are not separately asserted to be
gauge invariant; the full matched current includes those other terms.
For fixed normalized profile and source amplitude, c_q tends to zero
in the retained long-return scaling, so (8) supplies a quantitative
nonlinear source-response consequence, with no finite-flux physical
threshold inferred. A prescribed path phi(z) uses
K'=A_q(phi') and K''=A_q(phi'') in (4).

## Why the S087 stand-in is still insufficient

For its K(z)=s(z) sigma_x, s(z)=m(1+z)/(1+v(1+z)),

```text
F=log(1-s^2)+s^2,
F'=-2s^3 s'/(1-s^2),
(log det(I+K))'=-2s s'/(1-s^2).                           (9)
```

The difference is the derivative of the quadratic subtraction. Thus
even an exact det_3 response cannot stand in for the entire interacting
response. The new checker verifies this at the original v=1/100,
m=1/2,z=0 and tests noncommuting, nonnormal cases separately.

An exact counterexample to a derivative inferred from a uniform value
envelope is K_M(z)=(1/4+sin(Mz)/16) I_1. All z,M obey
||K_M||3<=5/16<1 and the same value bound, while at z=0
F'(K_M)=M/320. This divergence is precisely the missing K' control.

## Physical interface and remaining decision

The S085 clock contains no finite determinant functional. To enter
(4) into a clock or stress equation one must first derive a fixed-space
realization of its actual family D_z,V_z, including the fermion domain,
boundary work, induced vertices and the same subtraction. If that
realization has K_z=D_z^(-1)V_z, its formal derivative is
K'_z=-D_z^(-1)D'_z K_z+D_z^(-1)V'_z. The derivative expression is
usable only after proving the relevant compositions are bounded/in S3
and the family is C1 there. A moving boundary cannot silently be treated
as a fixed prescribed source. For a dynamical response, nonlocal time
dependence and retarded continuation must also be retained.

Integrating over a fluctuating gauge measure additionally requires
differentiability and domination for the weight and measure, including
large-field configurations. The existing global determinant value and
total-variation estimates in [the tail note](determinant-tail-control.md)
do not provide these unbounded current/stress estimates. Neither fixed
Euclidean source response nor an arbitrary small-norm matrix computes
the physical gauge covariance, nonlinear constraints, finite couplings,
unequal curved collar matching or the complete S074 remainder.

The strongest justified next mathematical input is a controlled
fixed-domain operator family for the actual geometry or clock insertion,
with its K' (and K'' for susceptibility) bounds. Repeating scalar
derivatives of unrelated inputs does not supply that map.

## Prior work and verification

[Britz et al., exact v3](https://arxiv.org/html/2007.12834v3), introduction,
equation (1.3), defines the general regularized determinant for complex
Schatten operators. Taking A=-K and k=3 gives the signs in (1). Only
that definition and its immediate surrounding scope were read for this
milestone; their product-anomaly theorem is not imported. The ordered
derivative proof above is supplied directly. No novelty claim is made.

The [checker](../checks/determinant_response_control.py) compares the
ordered identities with independent log determinants and central
differences, tests the mixed Hessian without commuting directions,
checks degree tails and Taylor error, exercises the existing unequal-
return Fourier Dirac map with distinct source directions and an
independent spacetime multiplier quadrature, and verifies (9) and the
value-envelope counterexample. Its [receipt](../evidence/determinant-response-control.json)
and [raw log](../evidence/determinant-response-control.log) record exact
inputs, source hashes, tolerances and results. Numerical examples
supplement the ordinary theorem; no Lean/Lake or physical opening is
claimed.

The retained run passes 12 matrix cases and 108 series-tail checks. The
largest first-difference error is 1.43e-10; the largest Hessian/path-second
difference error is 6.14e-7. For the specified unequal-return source,
the conservative source radius is 0.232240; (8) gives each real/imaginary
current component norm bound 1.402846 (combined bound 2.805692) and
susceptibility bilinear norm bound 11.739097.
These finite-source numbers are not assigned to the physical clock.
The first invocation failed because this embedded Python did not place
the checker directory on its import path. The explicit local import path
fix and the successful run are retained with that failure in the raw log.

Remaining to-do list: derive the physical clock/geometry operator map
and its response norms; match lower-order, local and boundary terms;
control gauge-weighted response and its tails; complete the constrained
state, retarded continuation, exterior matching and S074 error budget.
