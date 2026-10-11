# An interior dipole trial field and a global return-error criterion

2026-10-10. S052/S054. The finite-collar result alone cannot control an
arbitrary attached passive load. This note makes that limitation concrete
and develops the next construction: a field-aligned trial spinor in the
actual point-dipole exterior, together with a weighted current estimate
that would control the full return from its residual and boundary errors.

The first correction cancels the order-one width and nonlinear-mass
sources. On any fixed regular, trimmed interior tube its residual has an
ordinary `O(q^(-1/2))` bound, plus an exponentially small transverse cutoff
term. The mouth-end matching and a single global trial field are not yet
constructed. The global boundary existence source is under review.
Consequently no full physical return bound or novel wormhole result is
accepted by this milestone.

## Why the collar column does not settle feedback

For `0<tau<1`, the real unitary scattering matrix

```text
S_tau = [[r,tau],[tau,-r]], r=sqrt(1-tau^2),
||S_tau (1,0)-(1,0)||^2=2(1-r) ->0
```

has a prepared-input column arbitrarily close to perfect direct
reflection. Connect its second port to the causal passive load
`K(s)=-exp(-s tau^4)`. The loaded response is

```text
S_loaded(s)=r+tau^2 K(s)/[1+r K(s)].                  (1)
```

At every fixed positive real `s`, it tends to `-1` as `tau->0`, whereas
the uncoupled approximation tends to `+1`. Indeed,
`1-r exp(-s tau^4)=tau^2/2+O(tau^4)`. A positive delay and passivity do
not prevent this amplification when no uniform lower delay or loop-loss
condition is supplied. This is an elementary network counterexample,
not evidence that the physical dipole has such a resonance.

The [collar estimate](mouth-collar-scattering.md) remains valid with its
stated incoming data. The counterexample rules out promoting that column
estimate to a loaded-return bound using passivity alone. Instead, the
following global residual estimate includes feedback in the actual
solution and requires no estimate of an interface inverse.

## Weighted current criterion for a global trial

Let `Omega` be the flat exterior of the two spherical ports. Use its
outward normal `n`, so at a hole this normal points into the removed
ball. Let `P_out=(I+sigma.n)/2`, `P_in=(I-sigma.n)/2`. Suppose the chosen
causal Dirac boundary realization supplies the `L2` solution of

```text
(H-i xi)psi=0, H=sigma.(-i grad-A), xi>0,
P_in psi=J_in a.
```

The magnetic connection is treated in smooth gauge patches of its line
bundle. Gauge strings are not extra physical boundaries. This section is
an a posteriori estimate **conditional on that boundary realization**;
its complete existence/trace/dependency review is not claimed here.

For a global trial `Phi a`, put `f=(H-i xi)Phi a`, and let `b` be the
incoming trace of the difference `psi-Phi a`. Let `h(x)` be a
nonnegative 1-Lipschitz function, zero at the emitting port and at least
`D` at the receiving port. Euclidean distance from the emitting sphere
is one choice for equal spherical ports, with `D=d-2rho_p`.
For `0<alpha<1`, set `w=exp(alpha xi h)`.

For the error spinor `e=psi-Phi a`, the current identity gives

```text
div(e^dagger sigma e)=-2xi |e|^2-2 Im(e^dagger (H-i xi)e),
xi I-sigma.grad(log w) >=(1-alpha)xi I.
```

Integration, Cauchy-Schwarz and completing the square yield

```text
||w P_out e||_boundary^2
 <= ||w b||_boundary^2 + ||w f||_bulk^2/[2(1-alpha)xi]. (2)
```

The sign of the error forcing is immaterial to the last inequality.
To justify unbounded weights, first use
`exp(alpha xi min(h,R))`, apply an additional cutoff at spatial infinity,
remove that cutoff using the unweighted `L2` solution, and then let `R`
increase. The coercive bulk term gives the uniform weighted bound needed
for monotone convergence. Lipschitz weights follow by approximation, or
by the weak product rule. No real-energy spectral gap is used.

If `J_out` is the receiving flux isometry and `R_0` the desired comparison,
suppose uniformly for unit input `a` that

```text
||w b||<=B, ||w f||<=F,
||J_out^* P_out Phi a-R_0 a||<=delta_out.
```

Then (2) implies

```text
||R-R_0|| <= delta_out
       +exp(-alpha xi D) sqrt[B^2+F^2/(2(1-alpha)xi)]. (3)
```

Every incoming mismatch, including at the receiving mouth, belongs in
`B`. The norm is on the complete outgoing state before projection, so
unobserved angular excitations are not discarded by assumption. With
axisymmetry, estimates uniform in individual azimuthal sectors extend to
all their superpositions without an additional mode-count factor.

For `alpha=1/2`, define `F_D=sqrt(D) F`. On the selected window
`0.05<=xi D<=4`, a sufficient criterion for the earlier
`0.001 exp(-xi D)` return target is

```text
exp(4) delta_out + exp(2) sqrt(B^2+20 F_D^2) <=0.001.  (4)
```

For example, `B<=5e-5`, `F_D<=2e-5`, `delta_out<=1e-6` satisfy (4).
An exact rational certificate uses
`exp(1)<65/24+1/100<11/4` and
`sqrt[(5e-5)^2+20(2e-5)^2]<103e-6`, giving upper bound
`214049/256000000 <0.001`. These are **targets for unproved physical
norms**, not measured or certified return errors.

## Corrected trial field in the actual dipole interior

Rescale lengths by `d` in this section, so the centers lie at `z=+-1/2`,
the spectral parameter is `Xi=xi d`, and derivatives are dimensionless.
Fix a flux contour `u=u_0` with `u_0` in a compact subinterval of `(0,1)`.
Trim its endpoints away from the monopoles and use a regular tube of
fixed half-width `w_0` that stays away from the axis and both centers.
The width is fixed independently of `q` for this interior assertion.

Use the exact [field-aligned reduction](exterior-dirac-window.md), with
arclength `s`, signed normal coordinate `y` and `h=1-kappa(s)y`:

```text
H_F=i sigma_x partial_y-i sigma_z h^(-1/2) partial_s h^(-1/2)
                                            +sigma_y q m(s,y),
m(r,z)=[u(r,z)-u_0]/r, ell/q=1-u_0,
m(s,0)=0, b(s)=m_y(s,0)=|B(s,0)|/q>0,
c(s)=m_yy(s,0), omega(s)=b_s/(4b).
```

Only flux/charge pairs with half-integer `ell` represent physical axial
sectors. A compact contour family supplies uniform smooth bounds on a
regular trimmed tube; the earlier field estimate gives `b>=4u_0^3` in
these length units. Tube reach and endpoint matching still need explicit
uniform choices before the global return criterion can be applied.

Let

```text
g(s,y)=(q b/pi)^(1/4) exp(-q b y^2/2),
P=1-q c y^3/6, Q=omega y,
Phi_F=exp(-Xi s) chi(y) g(s,y) (P,Q).                 (5)
```

Here `chi=1` on `|y|<=w_0/2`, vanishes for `|y|>=w_0`, and
`|chi'|<=pi/w_0`. For instance use
`cos^2(pi(|y|/w_0-1/2))` in the transition region. This `C1` cutoff is
sufficient for the first-order operator and the `L2` residual.

The changing Gaussian width creates the upper source
`-i omega(1-2qb y^2)g`. The lower corrector `Q g` cancels it under
`i sigma_x partial_y+sigma_y q b y`. The quadratic mass creates the
lower source `i q c y^2 g/2`; the upper correction `-q c y^3 g/6`
cancels it. These are the two order-one terms that the leading Gaussian
alone leaves behind.

The full residual is still evaluated using the exact dipole mass and
curvature, not the linearized operator. Put
`Delta=m-b y`, `r_3=m-b y-c y^2/2`,
`G=omega(1-2qb y^2)`, and `h_s=-kappa_s y`. Before applying `chi`, its
two components divided by `i exp(-Xi s)g` are exactly

```text
U = G[-kappa y+q c y^3/6]/h +q c_s y^3/(6h)
     +Xi(h^(-1)-1)P +h_s P/(2h^2)-q omega y Delta,
V = q r_3-q^2 c y^3 Delta/6
     +omega_s y/h+Q G/h-Xi Q(h^(-1)+1)-h_s Q/(2h^2). (6)
```

Thus all uncancelled terms are small under the Gaussian scaling. Here
is a sufficient bound with explicit dependence on the fixed tube.
Suppose `|kappa y|<=1/2` and, throughout the tube interval,

```text
0<b_0<=b<=b_1, |c|<=c_0, |c_s|<=c_1,
|omega|<=a_0, |omega_s|<=a_1,
|kappa|<=k_0, |kappa_s|<=k_1,
|r_3|<=c_3 |y|^3, 0<=Xi<=Xi_max.
```

The Taylor remainder condition follows from a bound on the third normal
derivative of the exact smooth mass on the tube. Define positive
polynomials

```text
P_U(t)=2a_0(1+2b_1 t^2)(k_0 t+c_0 t^3/6)
       +c_1 t^3/3+2(Xi_max k_0+k_1)t(1+c_0 t^3/6)
       +a_0 t(c_0 t^2/2+c_3 t^3),
P_V(t)=c_3 t^3+c_0^2 t^5/12+c_0 c_3 t^6/6
       +2a_1 t+2a_0^2 t(1+2b_1 t^2)
       +3Xi_max a_0 t+2k_1 a_0 t^2.
```

For a polynomial with nonnegative coefficients, let `M` act linearly
with `M(t^k)=sqrt[(2k-1)!!/(2^k b_0^k)]`. Then, for `q>=1`, (6) and
Gaussian moments give

```text
||uncut residual(s,.)|| over |y|<=w_0
 <= exp(-Xi s) C/sqrt(q), C=M(P_U+P_V).              (7)
```

To check (7), expand (6) using `1/h<=2`. Each resulting term has the
form constant times `q^p |y|^k g` with `k>=2p+1`; hence its Gaussian
`L2` norm is bounded by its coefficient times `M(t^k)/sqrt(q)`.
The polynomials retain exactly those coefficients with `q` set to one
after the power check. This proves an upper order, not an assertion that
the leading coefficient is nonzero at every section.

The additional cutoff residual is `i sigma_x chi' exp(-Xi s)g(P,Q)`.
On its support split the Gaussian exponent to obtain

```text
||cutoff residual(s,.)|| <= exp(-Xi s) C_cut
                              exp[-q b_0 w_0^2/16],
C_cut=(pi/w_0)2^(1/4)
       [1+c_0 sqrt(15)/(6 b_0^(3/2))+a_0/sqrt(b_0)]. (8)
```

For example, `exp(-qb y^2)<=exp(-qb w_0^2/8)exp(-qb y^2/2)` on
`|y|>=w_0/2`; the remaining Gaussian moments prove (8). Integrating
(7)-(8) in `s` gives the corresponding bulk bound on a fixed interior
tube. If `s>=0`, the longitudinal factor has squared integral at most
`1/(2Xi)` for positive `Xi`.

## The remaining global matching step

The statement just proved is an interior estimate. A longitudinal cutoff
creates another derivative term that need not decrease with `q`; it
cannot be omitted. Nor does a tube ending on an arclength cross-section
automatically have the required trace on a spherical mouth.

To connect local fields `Phi_A` and `Phi_B` smoothly in an overlap, use
`Phi=chi Phi_A+(1-chi)Phi_B`. In the original Cartesian spin frame,

```text
(H-i xi)Phi = chi f_A+(1-chi)f_B
                     -i sigma.grad(chi)(Phi_A-Phi_B). (9)
```

The overlap discrepancy, weighted by `|grad chi|`, is therefore part of
the residual budget. Both fields must be expressed in the same physical
spin frame and gauge patch. This gives a concrete next requirement:
match (5) to the spherical mouth profiles and correctors, calculate the
weighted incoming and receiving traces, and bound (9). Outer exceptional
flux labels, the axis, continuum escape and curved throat matching
remain included in the eventual full-return problem.

## Literature scope and checks

[Ginoux and Murro v6](https://arxiv.org/html/2007.02544v6) gives a general
Friedrichs-system framework on manifolds with timelike boundary. The
introduction, section 2.2 and section 5 were read. The outgoing-normal
projector has the future-admissible sign and rank in Definition 2.13.
Theorems 1.1-1.2 and the global extension depend on sections 3-4 and
further references whose full audit is pending. A smooth magnetic
connection is a lower-order bundle term, but the exact exterior
realization and trace construction are not accepted solely from the
theorem statement.

[Barbaroux et al., Lemma A.2](https://arxiv.org/html/2509.17151v1#A1.Thmtheorem2)
and its full proof were read. Its Combes-Thomas conjugation already
produces the factor `[(1-alpha)|Im z|]^(-1)` for a magnetic Dirac
resolvent. Lemma A.3's proof was also inspected. The half-plane operator
and self-adjoint boundary realization differ from this exterior port
problem. Weighted resolvent decay is established mathematics, not a
novel claim of this program. The prepared curved-interface method was
already compared with [Bal et al.](https://arxiv.org/html/2106.00729v1#S3).

Run

```powershell
python checks/exterior_trial_residual.py --output evidence/exterior-trial-residual.json
```

The [receipt](../evidence/exterior-trial-residual.json) checks the source
cancellations, Gaussian power conditions, exact conditional budget and
weighted Pauli coercivity. Thirty actual dipole normal sections/frequency
cases use three flux values each, with two geometry-derivative resolutions.
The residual is evaluated in the exact field-aligned differential operator
and independently compared with expansion (6). The initial numerical
test incorrectly required a nonzero leading coefficient everywhere;
some sections decay faster. The corrected test checks the proved upper
order, and retains the scaled ratios. Nine passive feedback examples
verify (1). These supplementary calculations are not a global scattering
simulation or an interval bound for all tube geometry constants.
No Lean invocation occurred; GCP remains the exclusive Lean route.

Remaining to-do list: complete the global boundary proof review, construct
the mouth/tube overlap and trace matching, control the weighted full-return
residual and exceptional modes, finish energy/observer and novelty
comparisons, then certify suitable endpoints on GCP and prepare the manuscript.
