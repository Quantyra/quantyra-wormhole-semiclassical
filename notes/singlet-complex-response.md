# The singlet preserves the full subtracted response at damping resolution

2026-10-11, S053. The spherical-photon Gaussian model now admits a
mass-uniform comparison of its **complete complex subtracted stress
response**, including the real part. For frequency imaginary part `eta`,
the error is `O(1+eta^-4)`, independent of channel count, return lengths
and bounded singlet mass profile. Consequently it is `o(q delta)` at the
damping resolution `eta~delta` if `q delta^5->infinity`.

This strengthens the earlier
[positive-band absorption result](singlet-stress-continuum.md).
It does not determine the renormalized static vacuum curvature. That
remaining quantity, nonspherical photons and physical matching must be
controlled before using the [finite-channel decay theorem](finite-channel-decay.md)
for the physical wormhole. No static counterterm is adjusted here.

## 1. Fixed Gaussian model and precise result

Use exactly the common throat, unequal return graph, and zero/charge
sector convention of [the singlet model](schwinger-return-feedback.md).
There are `q` scalar current channels, a throat interval `[0,T]`, positive
finite exterior lengths `D_j`, and a real bounded mass profile `m(s)`
coupling only the unit vector `v=(1,...,1)/sqrt(q)`. The complete throat
energy `V`, including its mass term, is the clock vertex in
`H(phi)=exp(phi)V+H_ext`. No nonspherical photon is included.

Normalize the reference frequency to one. For its positive scalar
frequencies `omega_a`, let `U_ab` be the ordered-pair vacuum stress
weights. Define the convergent subtracted response by

`Delta chi_m(w)=sum_ab 2w^2 U_ab`
`              /[(omega_a+omega_b)(w^2-(omega_a+omega_b)^2)]`. (1)

For fixed constants `0<c0<c1<infinity`, there is `C(c0,c1)` such that
for `c0<=Re w<=c1`, `0<eta=Im w<=1`,

`|Delta chi_m(w)-Delta chi_0(w)| <= C(c0,c1)(1+eta^-4)`.      (2)

The constant is independent of `q`, all `D_j`, `T`, and the magnitude
or shape of `m`, as long as every individual mass profile is bounded.
The theorem allows the mass bound to grow along a sequence of finite
models. An untruncated singular AdS endpoint is not thereby included.
The estimate is sufficient, with no claim that the power four is sharp.

If dimensional frequency is restored using reference scale `omega_ref`,
multiply the right side by `omega_ref` and replace `eta` by
`Im(w)/omega_ref`. This is an estimate of the retarded analytic function
in the upper half-plane, not a pointwise real-axis bound between discrete
spectral measures. The first harmonic and all zero modes are handled by
the actual pair weights, not by assigning them a small positive frequency.

## 2. The two-condition comparison and the stress observable

On the scalar graph Hilbert space let
`d_m f=(f',m v*f_throat)` and introduce the self-adjoint first-order
operator `D_m=[[0,d_m*],[d_m,0]]`. The second component includes the
auxiliary throat mass space. Let `X` be the throat indicator on the
scalar and derivative spaces, and the identity on that auxiliary space.
It is a projection with norm one, independent of the mass.

For positive-frequency scalar modes `f_a`, the normalized eigenvectors
are `e_a^+/-=(f_a,+/-d_m f_a/omega_a)/sqrt(2)`. Thus

`U_ab=(omega_a omega_b/2)|<e_a^-,X e_b^+>|^2`.              (3)

Restrict `d_m` by the two endpoint conditions
`v*f_throat(0)=v*f_throat(T)=0`, and call the resulting first-order
operator `D_m,K`. The common symmetric restriction has defect two.
Its resolvent difference has rank at most two, and for nonreal `z`,

`||(D_m-z)^(-1)-(D_m,K-z)^(-1)||_1 <=4/|Im z|`.             (4)

The original domain and adjoint argument is in the previous absorption
note. It uses continuity and surjectivity of the two endpoint traces in
the graph norm; no uniform mass or length bound is required. After the
cut, the singlet throat is isolated and `X=I` on it. Its negative-to-positive
matrix elements vanish by spectral orthogonality. The remaining reflecting
network and its observable are independent of `m`. Consequently the cut
stress pair statistic is identical for `m` and zero, for every bounded
spectral test supported in the negative-positive quadrant.

For `F` smooth and compactly supported away from both axes, put

`I_D(F)=sum_(lambda<0<mu) F(lambda,mu)|<e_lambda,X e_mu>|^2`.

The two-variable Cauchy-Green formula and (4) give

`|I_Dm(F)-I_D0(F)| <= (64/pi^2)`
`                         *sum_(i,j=0)^2 ||partial_1^i partial_2^j F||_1`.
                                                                  (5)

Here the cut cancels between the two comparisons. To check the constant,
an order-one almost-analytic extension has antiholomorphic derivative
`(i chi'/2)f-(y chi'/2)f'+(iy chi/2)f''`. An even cutoff supported in
`[-2,2]`, equal to one on `[-1,1]`, gives integrals against `1/|y|`
bounded by `1,1,2`. Expanding the difference of
`R_D(z) X R_D(t) X` costs `8/(|Im z||Im t|)` in trace norm per cut.
The two cuts, two scalar Cauchy formulas and those three bounds give (5).
Only the trace-class **difference** is integrated; the individual
resolvent products can have infinite trace from the auxiliary zero kernel.

Crucially (5) also holds for `D_m/r,D_0/r`, for any `r>0`, with the same
constant. It is a resolvent-rank statement, invariant under a simultaneous
change of frequency unit. This permits frequency-dependent imaginary
extension scales without imposing a global high-frequency cutoff.

## 3. Removing the axes with an integral derivative estimate

Write `a=-lambda>0`, `b=mu>0`, `s=a+b`. From (1) and (3), the test
function to be controlled is

`F_w(-a,b)=w^2 ab/[s(w^2-s^2)]`.                           (6)

Choose once a nonnegative smooth partition `sum_(k in Z) theta(s/2^k)=1`
for `s>0`, with `theta` supported in `[1/2,2]`. Choose a smooth cutoff
`beta(t)=0` for `t<=1/2`, `beta(t)=1` for `t>=1`. Initially multiply
by `beta(a/epsilon)beta(b/epsilon)` to avoid the axes.

The elementary estimate needed on the unit annulus is

`sum_(i,j=0)^2 ||partial_a^i partial_b^j`
` [ab g(a+b) beta(a/e) beta(b/e)]||_1`
` <= C_beta sum_(k=0)^4 ||g^(k)||_L1([1/2,2])`,             (7)

when `g` is supported in that annulus; the constant is independent of
`e>0`. The functions are extended by zero outside the positive quadrant.

To prove (7), expand the nine derivatives by the product rule and change
variables to `(s,a)`, `0<a<s`. Terms without differentiated cutoffs have
bounded polynomial factors on this annulus. On a cutoff transition,
`a beta'(a/e)/e` is bounded and its integral is `O(e)`, while
`[2 beta'(a/e)+(a/e)beta''(a/e)]/e` has bounded `L1(da)` norm.
These are precisely the first and second derivatives of `a beta(a/e)`.
The corresponding statements hold for `b`. If `e<1/4`, both transitions
cannot occur on this annulus simultaneously. The factor from the other
variable is then a bounded polynomial; integrating along the transition
leaves the displayed one-dimensional norms of `g` and its derivatives.
For `e>=1/4`, every cutoff derivative is uniformly bounded on the
annulus. This covers all terms and proves (7).

The second derivative's supremum does grow like `1/e`; its integral
does not. Replacing that integral by its supremum would recreate an
unnecessary soft-mode error. For each finite graph there are only finitely
many nonzero frequencies in a compact band, and a positive gap above its
zero subspace. Therefore taking `epsilon->0` in a fixed annulus is
legitimate. The uniform bound (7) survives that limit. This is not an
assumption of a gap uniform in the sequence of graphs.

## 4. Summing both frequency tails

For `r=2^k`, rescale the annular contribution by `a=r A`, `b=r B`.
Its test function for `D/r` is `AB g_(r,w)(A+B)`, where

`g_(r,w)(s)=r w^2 theta(s)/[s(w^2-r^2 s^2)]`.               (8)

Apply (5) in those frequency units, followed by (7). Uniformly on the
stated frequency rectangle,

`sum_(j=0)^4 ||g_(r,w)^(j)||_1 <= C r`       for small `r`,
`sum_(j=0)^4 ||g_(r,w)^(j)||_1 <= C/r`       for large `r`. (9)

For example, take `r<c0/4` in the first bound and
`r>4sqrt(c1^2+1)` in the second. The denominators are then separated
from zero on `[1/2,2]`, and ordinary differentiation proves both bounds.
These geometric series control the infrared and ultraviolet tails of
the **difference**, independently of mass and graph size. Truncating the
annular sum to pair frequencies between `R_low` and `R_high` leaves an
interaction-comparison error `O(R_low+1/R_high)` once those scales lie
outside the resonance band. This is not a bound on either individual
response tail uniform in `q`.

Only finitely many annuli remain. On each, partial fractions separate
the simple denominator `w-rs` from factors bounded away from zero.
For `j>=1`, the integral of its `j`th derivative is bounded by
`C eta^-j`; for `j=0` the bound is `C(1+|log eta|)`.
Thus the sum of the derivative norms through order four is
`C(1+eta^-4)`. This proves (2), subject to convergence of the individual
responses, which can be established without an extra ultraviolet premise.

At `w=i`, (6) is the positive function `ab/[s(1+s^2)]`.
The free response is finite by its known circle spectral sum. Apply
the annular comparison bounds to positive partial sums. Their differences
are absolutely summable by (9); the massive positive sum is therefore
finite as well. Axis cutoffs can be removed monotonically. For general
upper-half-plane `w`,

`|F_w(-a,b)|/F_i(-a,b) = |w|^2(1+s^2)/|w^2-s^2|`

is bounded on `s>0`, locally uniformly in `w`. Dominated convergence
therefore defines the analytic function (1) and validates all limits
used above. This proves its sharp-profile subtracted convergence for
each bounded mass profile. No two divergent unsubtracted quantities
have been subtracted after separate summation.

## 5. Consequence for the clock, with its static premise exposed

In the finite inverse-tail family let `delta->0`, `q delta^5->infinity`.
For `w=Omega_delta+delta z`, with `z` in a compact upper-half-plane set,
(2) gives

`[Delta chi_m(w)-Delta chi_0(w)]/(M delta)`
`                 =O(1/(q delta^5))->0`, `M=q/8`.           (10)

The convergence includes the real dispersive part. It is uniform in
all bounded-per-model throat singlet masses. The working
`delta=q^(-1/16)` gives the relative bound `O(q^(-11/16))`.
This also implies the earlier free-channel condition `q delta^3->infinity`.

A precise conditional corollary now isolates the remaining conservative
input. Suppose the **derived renormalized** static clock curvature obeys

`K0_m-K0_free=o(q delta)`.                                 (11)

For the same canonical clock vertex and kinetic coefficient, define
`K_m(w)=K0_m-Mw^2+Delta chi_m(w)`. Then the scaled-kernel limit,
positive spectral representation and finite-time decay conclusion in
the previous note remain valid, with its exact free-continuum center.
To justify the spectral step, (1) at `w=i` gives summability of
`2U_ab/[M(omega_a+omega_b)^3]` at high frequency; the finitely many
positive low modes cause no divergence. The same closed positive form
therefore realizes `M/K_m`, and its inverse second moment is `M/K0_m->1`.
Condition (11) makes `K0_m` positive eventually. No separate assumption
of a continuum bath is used.

This is a corollary with a stated unresolved physical premise, not a
license to choose `K0_m` to satisfy it. A frequency subtraction controls
neither the vacuum energy nor its clock derivatives. The photon reduction
fixes a mass profile and renormalization/matching terms that must supply
(11). Other current interactions need their own bounds. The result does
not prove that a massive boson with imposed graph matching is the complete
fermion/photon state of the physical mouths.

## 6. Verification, prior work and decision

Run `python checks/singlet_complex_response.py`; see its
[receipt](../evidence/singlet-complex-response.json) and
[source review](../evidence/singlet-complex-source-review.json).
The checks cover the pair kernel, rescaling, axis derivative norms,
resonant and nonresonant annuli, finite-rank trace comparisons, and
independent finite-element massive/free complex stress calculations.
The ordinary estimates, not numerical parameter scans, establish mass
and channel uniformity. Finite-element checks retain their cutoff and
mesh limits explicitly.

The functional-calculus method is established prior art. In
[Peller, 1509.02803v1](https://arxiv.org/abs/1509.02803v1), sections 1.1-1.2
provide a comparator for double operator integrals and their hypotheses.
Here `D_m-D_m,K` is not a bounded finite-rank operator: it is the
**resolvent difference** that has finite rank. Equation (5) supplies the
applicable argument directly rather than importing a bounded-perturbation
theorem. [MMP v3, Appendix A](https://arxiv.org/html/1807.04726v3#A1) and
[the later reduction, section 4.3.1](https://arxiv.org/html/2606.12528v2#S4.SS3.SSS1)
are physical predecessors for the singlet. They do not supply this
unequal-return complex stress comparison. No novelty acceptance follows
without the remaining physical result and full citation review.

Continue the existing damping primary and finite-duration energy-transfer
fallback. The Gaussian dynamical singlet gap is closed at damping
resolution; the static and nonspherical interaction gaps remain open.
No Lean, Lake, GCP or publication action occurred.

Remaining to-do list: derive the singlet static-curvature bound and control
nonspherical interactions; close physical background, complete returns and
state matching; finish novelty and the consequential result; independently
verify, GCP-certify suitable endpoints and prepare the manuscript.
