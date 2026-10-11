# A joint small-mouth, large-flux trial with reflected-mode control

2026-10-10. S051/S052/S054. The fixed-port trial and the spherical-load
estimate can now be combined in a joint limit. Use exact spherical
dipole modes near each mouth and the corrected field-aligned trial in
the exterior, with an explicit overlap error. For trimmed flux labels,
ports `rho_p/d=q^(-3/4)` and overlap radius `rho_a/d=c_a q^(-1/4)` give
a conditional loaded-return error `O(q^(-3/8))` on a fixed positive
imaginary-frequency window. Every higher angular mode enters the boundary
estimate. The constants depend on the trim and window; no numerical flux
threshold is certified.

This is a result for the specified flat dipole operator coupled to a
spherical boundary load, conditional on a solution with the stated
current/trace identity. The actual matched operator, smooth physical
metric/gauge interface, full vacuum energy and its variation remain
unproved. In particular, this estimate does not restore the exterior
front bound after transmission through the throat has been attached.

## Statement and definitions

Set `d=1`; dimensional norms follow by rescaling. Fix `0<delta<1/2` and
`0<Xi_0<Xi_1<infinity`. Use the exact static point-dipole Weyl operator
outside two balls of radius `p`, with flux integer `q`, and retain its
sectors `ell=q(1-u) in Z+1/2`, `delta<=u<=1-delta`. All exterior angular
and continuum states are retained; only the prescribed ground inputs
are trimmed.

The source/receiver ground bases are the isolated spherical monopole
bases, in the common gauge and spin-frame phase convention of the
[global trace calculation](global-mouth-trial.md). Let `L_p(u)` be the
actual dipole arc length between the two radius-`p` spheres and set
`R_0(u,Xi)=exp(-Xi L_p(u))`. This is an open-path phase convention, not
a determination of a closed physical loop's holonomy.

Attach the higher-mode load of a prescribed spherical throat as in
[the load estimate](throat-load-feedback.md). Suppose its two disjoint
optical collars have length at least `c_ell p` and `N/R<=C_beta/p`,
with positive constants independent of `q`. Use

```text
beta=C_beta/p, ell_c=c_ell p,
lambda_n=sqrt(n(n+q)),
delta_n=Xi/(Xi+lambda_n beta)
                  [1-exp(-2(Xi+lambda_n beta)ell_c)].          (1)
```

Increasing the upper bound `beta` and shortening a collar to `ell_c`
preserve the proved loss estimate. Let `p<a`, `2a<1/10`. For all
sufficiently large `q`, there is one smooth global exterior trial for
each retained input, with constants uniform in `u,Xi,p,a`, such that

```text
||f|| <= C [(qa)^(-1/2)+a^(3/2)+sqrt(q)a^(7/2)]
                       + C a^(-1/2) exp(-c q),                (2)
||b_G||+epsilon_out <= C q p^4,
sum_n delta_n^(-1)||b_H,n||^2
 <= C [q p^4+q^(3/2)p^3+p^4+sqrt(q)p^3+q p^5].               (3)
```

Here `f` is the exact flat exterior residual, `b_G,b_H` are the actual
loaded boundary mismatches defined in (6) of the load note, and
`epsilon_out` is the ground receiving trace error relative to `R_0`.
The ground estimate requires `sqrt(q)p^2` sufficiently small, as it is
for the scale chosen below. Constants in this note may depend on the
fixed trim, frequency window and collar constants, but not on `q,p,a`.
They are not claimed small or optimal.

Conditional on existence of the loaded solution with `L2` traces and
the exterior current identity, the load note gives

```text
||R_loaded-R_0|| <= epsilon_out
 + sqrt[||b_G||^2+sum_n delta_n^(-1)||b_H,n||^2+||f||^2/(2Xi)]. (4)
```

Thus `p=q^(-3/4)`, `a=c_a q^(-1/4)` with fixed `c_a>0` implies

```text
||R_loaded-R_0|| <= C q^(-3/8).                               (5)
```

Equations (2)-(5) hold in operator norm on arbitrary superpositions of
the retained ground inputs. No factor `sqrt(q)` is lost in superposition:
the trial and residual sectors are azimuthally orthogonal. The load
is contractive in each angular eigenspace, and the weight (1) is scalar
there. Unitary identifications that preserve those eigenspaces leave
the boundary estimate unchanged.

## Uniform tube geometry as the mouths shrink

Let `r_*(s)` be a positive smooth scale comparable to the distance to
the nearer monopole, capped at one, along a trimmed field line. It can
be chosen with `|r_*'|<=C` and `|r_*''|<=C/r_*`. In each small mouth
neighborhood, use center distance as the curve parameter. The remote
field is smooth, whereas the monopole field is radial of size
`1/(2r_*^2)`. Consequently the field/radial tilt is `O(r_*^2)` and
`ds/dr=1+O(r^4)`.

Here are the needed uniformity details. Rescale a small annulus by its
center distance. The normalized monopole field and flux profile are
smooth away from the axis; the remote contribution and its derivatives
have uniform scaled bounds. Trimming keeps the limiting ray angle in
a compact subset of `(0,pi)`. The scaled field-line ODE is therefore
uniformly regular down to zero center distance. On the remaining compact
interior, the exact prolate field-line parametrization is embedded and
regular. The two mouth neighborhoods are disjoint, and each limiting
mouth ray stays apart from nonadjacent curve pieces. These facts give
a tubular radius `w(s)=eta r_*(s)` for one fixed sufficiently small
`eta>0`, including a uniform bound `|kappa y|<=1/2`.

For the mass `m=(u(x)-u_0)/r_cyl` and the notation of
[the interior residual](exterior-trial-residual.md), scaled differentiation
then gives

```text
b ~ r_*^(-2), |c_2|<=C r_*^(-3), |c_2,s|<=C r_*^(-4),
|omega|<=C/r_*, |omega_s|<=C/r_*^2,
|kappa|<=C/r_*, |kappa_s|<=C/r_*^2,
|m-by-c_2 y^2/2|<=C |y|^3/r_*^4.                             (6)
```

The cylindrical radius is comparable to `r_*` near each mouth by the
trim; it is bounded below on the compact interior. Thus division by
that radius causes no unaccounted singularity in (6).

Use the corrected Gaussian `(1-q c_2 y^3/6,omega y)` with the exact
spin frame and measure factors, multiplied by a smooth transverse
cutoff in `y/r_*(s)`. Every residual monomial from the displayed exact
formula in the interior note has sectional norm at most
`C/(sqrt(q)r_*)` for bounded `Xi`. For example,
`q c_2,s y^3` contributes
`q r_*^(-4)(r_*/sqrt(q))^3=1/(sqrt(q)r_*)`.
The other terms obey the same or a stronger bound using (6).

The varying cutoff contributes both a normal and a longitudinal
derivative. Their sizes are `O(1/r_*)`; both lie in the Gaussian tail
`|y|>=eta r_*/2`. Since `b r_*^2` has a uniform positive lower bound,
their sectional norm is at most `C exp(-cq)/r_*`. Polynomial correctors
are bounded by splitting the Gaussian exponent and using its moments.
On the exposed tube, center distances are at least a constant times `a`,
and `integral ds/r_*^2<=C/a`. This proves the tube contribution in (2).

This is the step unavailable from a fixed-geometry constant with no
tracked dependence on port radius. It does not assume a positive gap
for the complete noncompact exterior.

## Exact collars and overlap

Near the source use the exact deformed intrinsic zero mode `g(r)` and
its lower-spin inverse correction `h(r)` from
[the mouth calculation](mouth-mode-matching.md):

```text
Phi_c(r)=exp[-Xi(r-p)] (g(r),h(r)).
```

At the receiving mouth use the reflected convention and amplitude
`exp(-Xi L_p(u))`, with the radial exponential reversed. The exact
[collar calculation](mouth-collar-scattering.md) bounds either residual
by `C(r+sqrt(q)r^3+Xi r^2)`. Integrating from `p` to `2a` gives the
collar contribution `C(a^(3/2)+sqrt(q)a^(7/2))` in (2).

Join each collar to the tube with a radial cutoff equal to one below
radius `a` and zero above `2a`, whose gradient is at most `C/a`.
Construct the tube down to center radius `a/4`, with its longitudinal
end cutoff entirely in the region where the radial tube multiplier
vanishes. Thus no longitudinal end-cutoff residual is exposed.

On every sphere in the overlap, the tube's desired radial spin component
and `g(r)` differ in `L2(dtheta)` by
`O(q^(-1/2)+r^4)`. Its undesired component and `h(r)` differ by
`O(q^(-1/2)+r^2)`. To see the uniformity, use the scaled angular variable
`sqrt(q A_r)(theta-theta_0(r))`. The maximum is the exact flux intersection,
`A_r=r^2 B_radial/q` is bounded above and below on the trimmed family,
and the normalized logarithmic profile and its first three angular
derivatives have uniform bounds as `r->0`. Taylor's theorem with Gaussian
domination gives the `q^(-1/2)` error, as in the global trace proof but
now uniformly in small `r`. The exact spin/measure rotation gives tilt
error `O(r^2)` and Gaussian squeeze error `O(r^4)`.

The bounds extend over the whole sphere: outside a fixed angular
neighborhood of the unique maximum, monotonic flux and regular pole
powers give uniform exponential tails. The tube is zero there. Also
`||h||<=C r^2`. Its detailed cancellation with the wrong-spin trace is
not needed. The phase conventions differ by at most `O(Xi a^5)` because
`ds/dr=1+O(r^4)`, and this is included in the displayed error.

Hence the complete overlap difference is at most
`C(q^(-1/2)+a^2)` on each sphere. The commutator in the exact Cartesian
spin frame is `-i sigma.grad(chi)(Phi_c-Phi_t)`. Squaring, integrating
over radial thickness `O(a)` and taking a square root gives
`C[(qa)^(-1/2)+a^(3/2)]`. The same bound holds at the receiving mouth.
This proves (2), including the overlap term rather than omitting it.

## The boundary errors at the inner ports

Only the exact collar field reaches radius `p`. Let `P_0` be the isolated
ground projection and `M=q p^2/[4(1-p)^2]`. The load note proves

```text
||(I-P_0)g||^2<=M^2/(q+1),
||Lambda^(1/2)(I-P_0)g||^2<=M^2/sqrt(q+1).                   (7)
```

For normalized positive profiles in a fixed sector, write
`P_0 g=gamma g_0`, `gamma>0`. Then
`1-gamma <= ||(I-P_0)g||^2` whenever the small-error regime is used;
in fact this follows directly from `gamma=sqrt(1-||(I-P_0)g||^2)`.
The prescribed input and receiving projection errors are therefore
`O(q p^4)`. Normalizations and phases at the two mouths are the same
ones used to define `R_0`.

The corrector satisfies

```text
||h||<=C p^2,
||L_0 h||<=p||g_r||+M||h||
          <=C [sqrt(q)p^2+q p^4],
||Lambda^(1/2)h||^2<=||h|| ||L_0h||
          <=C [sqrt(q)p^4+q p^6].                            (8)
```

Here `L_0` is the unit-sphere operator and `Lambda=|L_0|` on its
complement. The derivative bound is the proved
`||g_r||<=(3/sqrt(2))sqrt(q)r`. The lower-spin corrector has no ground
component in the relevant mouth convention.

For (1), `delta_n^(-1)<=C[1+lambda_n/(Xi p)]` uniformly on the window.
Indeed its exponential denominator is at least
`1-exp(-2 C_beta c_ell sqrt(q+1))`. Applying this to (7)-(8) gives
respectively the first two terms and the last three terms of (3).
Finally `b_H=K out_H Phi-in_H Phi`; the weighted norm commutes with the
angular blocks of the spherical load and `K` is contractive there.
The squared triangle inequality proves (3) without truncating any
angular mode or assuming a mode-independent load loss.

With `a=c_a q^(-1/4)`, all three polynomial bulk terms in (2) are
`O(q^(-3/8))`. With `p=q^(-3/4)`, the largest boundary cost in (3) is
`O(q^(-3/4))`. Its square root has exactly the same upper order. This
proves (5); the exponential term is eventually smaller. The exponents
balance these particular sufficient bounds, not an optimal error theorem.

## Physical and novelty boundary

This removes the former incompatibility between an unweighted fixed-port
Gaussian approximation and the stronger norm required by reflected modes.
It does not prove that the mathematical spherical load is the physical
MMP throat load. The finite-sphere dipole field has an angularly deformed
normal magnetic field; a spherical throat has a uniform monopole field.
Identifying spinor port spaces does not make those gauge-invariant fields
or the metrics match smoothly. A physical interpolation, its Maxwell
sources and its effect on the estimate must be derived. No unspecified
matching shell is asserted to solve this problem.

The necessary geometric hierarchy is not algebraically inconsistent:
with fixed small `g`, choosing `d/r_e=q^nu`, `3/4<nu<1`, gives
`rho_p/r_e=q^(nu-3/4)->infinity`. The established MMP leading length
obeys `ell_0/r_e=16 pi q/g^2`, so `d/ell_0->0` and
`Omega ell_0=32 pi q^(1-3nu/2)/g^2->0`. This is only a scale comparison
on the known long-throat branch. It is not a bound on metric/gauge
perturbations of the large-flux wavefunction, gauge running or interactions.
In particular a small coefficient perturbation can act on derivatives
that grow with flux; it cannot simply be added as an unweighted small
operator norm. The asymptotic comparison is not a certified finite
parameter regime of a fixed interacting theory.

[MMP v3](https://arxiv.org/html/1807.04726v3), (2.3), (5.31), (5.44)-(5.47)
and (6.50)-(6.51), supplies those baseline scales and already finds a
long rather than shortcut wormhole. Its section 5.5 also discusses
nonuniform angular stress and the expected subleading geometric response.
The new matching estimates do not establish that expectation at the
full energy level. Established curved-interface packet methods, including
[Bal et al.](https://arxiv.org/html/2106.00729v1), remain methodological
predecessors. No novelty acceptance or primary/fallback review exit is
inferred from this construction.

The operator existence, excluded channels, all-frequency front behavior,
length dependence of the load and its separate determinant, state and
holonomy, full renormalized energy, and common observer clocks remain
required for the intended physical consequence. Do not apply the earlier
absorbing-port energy theorem to (5) without those steps. No new Lean
endpoint is selected before that consequence and its novelty survive.

## Supplementary verification

Run `python checks/joint_mouth_matching.py --output
evidence/joint-mouth-matching.json` with the pinned requirements. The
[receipt](../evidence/joint-mouth-matching.json) contains seven exact rational
exponent checks, 54 two-mouth overlap cases, 18 bulk-residual integrations
and five analytic-order budget cases. It evaluates the actual dipole
field, spherical profiles and integral corrector, and compares expanded
residuals with direct Pauli-operator evaluation. The normal cutoff includes
both normal and longitudinal derivatives. Overlap grids use 1,025 and
2,049 points; bulk integrals use 64 and 96 longitudinal nodes, 256 and
512 normal nodes, and two geometry-derivative steps. Numerical angular tails are not substitutes for
the ordinary all-mode estimates (7)-(8).

The first implementation used a root bracket that admitted both branches
on large spheres. It now selects the source-to-midplane branch. A second
attempt omitted the tube cutoff when integrating broad Gaussian tails;
the implemented cutoff and its derivatives now match the trial definition.
The lower-flux cases then showed that the cutoff residual dominates the
polynomial order. They are retained; later flux values test the limiting
upper order. Refinement was increased instead of loosening its tolerance.
These numerical failures and corrections are retained in the external
execution audit. They do not supply a finite threshold for (5).

Decision: retain this complete conditional joint-matching argument and
proceed to the matched operator and physical interface, with the energy
and novelty gates still open. Do not expand a separate toy-model paper.

Remaining to-do list: construct the loaded operator and a controlled
physical metric/gauge interpolation, bound full energy and load derivatives
with excluded modes and interactions included, finish the common-observer
and primary/fallback novelty assessment, then verify and GCP-certify suitable
complete endpoints and prepare the reproducible manuscript.
