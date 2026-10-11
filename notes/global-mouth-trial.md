# A global dipole trial and its two spherical-mouth traces

2026-10-10. S051/S052/S054. This note extends the corrected interior
trial through both spherical boundaries. Its longitudinal cutoffs are
supported entirely inside the removed balls. The resulting field on the
whole flat exterior has a vanishing weighted residual as magnetic flux
increases at fixed geometry. Both mouth traces have explicit, uniform
large-flux limits on a trimmed family of flux contours.

The leading trace is not exactly an intrinsic spherical zero mode.
Field tilt creates a wrong-spin component of order `(rho/d)^2` and a
Gaussian squeeze mismatch of order `(rho/d)^4`. These errors remain at
fixed mouth geometry as `q` increases. Keeping them in the incoming
error budget gives a conditional large-flux prepared-return estimate.
The exterior boundary realization still has an incomplete source-proof
audit; no finite numerical flux threshold, curved-throat matching,
full physical return or novelty acceptance is claimed.

## Geometry and the global trial

Use the exact point-dipole Weyl operator, its magnetic line bundle and
the conventions of [the interior trial](exterior-trial-residual.md).
Rescale `d=1`. The positive source is at `(r,z)=(0,-1/2)` and the negative
receiver at `(0,1/2)`. Remove balls of equal radius `0<rho<=1/10`.
Fix `0<delta<1/2` and contours `u in [delta,1-delta]`; retain physical
sectors `ell=q(1-u) in Z+1/2`. Their sphere intersections are unique:
the enclosed angular flux is strictly increasing with polar angle,
as established in [the mouth calculation](mouth-mode-matching.md).

Orient each contour with the magnetic field. Let `s=0` at its source
sphere intersection and `s=L(u,rho)` at the receiver intersection.
Extend it at both ends to center distance `rho/2`. The extended family
is compact, smooth, embedded and avoids the axis and monopoles. Here
compactness refers to the family parameter and the bounded extended
arclength intervals, not to the unbounded physical exterior. Smooth
dependence follows either from the field-line ODE away from its centers
or the explicit prolate parametrization. Monotonicity in prolate eta
precludes a self-intersection. Compactness, bounded curvature and
separation of nonadjacent curve segments supply a common positive
tubular radius. Shrink it to `w<rho/8`, away from the axis, with
`|kappa y|<=1/2`. No numerical value of this uniform reach is asserted.

Choose a longitudinal cutoff equal to zero before center distance
`rho/2` and equal to one from `3rho/4` onward at each end. Its derivatives
are supported where the centerline is at most `3rho/4` from the relevant
center. Every point of that normal cross-section has center distance
at most `3rho/4+w<rho`. Thus both longitudinal derivative supports lie
strictly inside removed balls. In particular, wherever the tube meets
the physical exterior the longitudinal cutoff is already one.

Multiply the interior trial by these cutoffs and extend by zero outside
the tube. After restriction to the exterior this is a single compactly
supported `H1` spinor. A smooth transverse cutoff may be used; the earlier
`C1` cutoff is also sufficient. With `h=1-kappa y`, write its Fermi-frame
components as

```text
Phi_F = N exp(-Xi s) chi(y) (q b/pi)^(1/4) exp(-q b y^2/2)
                          (1-q c2 y^3/6, omega y),                 (1)
b=|B|/q, c2=m_yy(s,0), omega=b_s/(4b), Xi=xi d.
```

The source normalization `N` is specified below. The source gauge
`A_phi=q(1-u)` is smooth on the supporting tube; the field vanishes near
the axis. It therefore defines a global bundle section without treating
a gauge string as a boundary. The azimuthal spin/phase factors obey the
half-integer `ell` condition. All frames and phases are fixed by this
common gauge and continuous field frame.

The residual formulas (6)-(8) of the interior note apply everywhere in
the exterior support: there is no additional longitudinal cutoff term.
For the distance weight `w_x=exp(alpha Xi dist(x,B_source))`, bounded
`Xi`, and fixed `rho,delta`, compactness gives finite uniform constants

```text
||w_x (H-i Xi)Phi|| <= C_1/sqrt(q)+C_2 exp(-c_0 q), c_0>0.         (2)
```

Indeed, multiply the prior sectional bound by the maximum of `w_x`
over the fixed tube and by the square root of the finite integral of
`exp(-2Xi s)` over the extended arclength interval. That integral remains
bounded uniformly on any compact frequency window, including zero.
The small part of the extended interval with `s<0` is finite. The
normalization `N` is bounded. These constants are allowed to depend on
`rho,delta`; (2) is not a uniform joint limit as `rho` or `delta` tends
to zero. Orthogonal azimuthal sectors permit arbitrary retained-sector
superpositions with no additional square-root-of-mode-count factor.

## The exact sphere trace map

At the source let the centerline intersect at polar angle `theta_0`.
Use the unit radial vector `e_r=(sin(theta),cos(theta))` and field tangent
`t=(sin(beta),cos(beta))`. Set

```text
a=beta_0-theta_0, c=cos(a)>0, tau=tan(a/2),
A=rho^2 (B.e_r)/q = rho^2 b_0 c,
N=sqrt(2c/(1+c)).                                                 (3)
```

The symbol `c` in (3) is a tilt cosine, distinct from the mass derivative
`c2` in (1). For points of the sphere inside the tube, solve
`x(theta)=gamma(s(theta))+y(theta)n(s(theta))`, with
`n=(-cos(beta),sin(beta))`. At the central intersection,

```text
y_theta=-rho c,     s_theta=rho sin(a).                           (4)
```

The unitary spherical spinor `eta` has boundary norm `integral |eta|^2
dtheta`. Comparing the cylindrical and spherical measure factors gives
the exact trace relation

```text
eta(theta)=sqrt(rho/h) R_y(beta(s)-theta) Phi_F(s,y),
R_y(v)=[[cos(v/2),-sin(v/2)],[sin(v/2),cos(v/2)]].                (5)
```

This retains the change of longitudinal coordinate along the sphere,
the spin rotation and the Fermi Jacobian. Omitting any of them changes
the endpoint error. The exact incoming/outgoing projectors use the
normal pointing out of the exterior, into the removed ball. At the
source, the radial upper component is incoming to the exterior; at the
receiver, the radial lower component is outgoing from the exterior.

Reflection across `z=0` maps paired mouth points by
`theta_receiver=pi-theta_source`, `s_receiver=L-s_source`, with
`beta_receiver=-beta_source`, equal `b,c2,kappa`, and opposite `omega`.
Thus the receiver rotation angle is `-pi-a` at its center. Choose the
receiving angular isometry to be minus the radial lower component. This
is the phase obtained from the displayed continuous field spinor; it
makes the leading transfer coefficient positive. A subsequent closed
wormhole loop must still match its physical throat/gauge holonomy.
This convention does not remove a physical loop phase.

## Uniform limiting traces and the persistent endpoint error

Let `g_q(theta)` be the exact normalized deformed spherical zero mode
from the mouth calculation in this sector. Its unique maximum solves
`A_phi(theta_0)=ell`, and

```text
(log g_q)''(theta_0)=-q A.
```

On the scale `T=sqrt(q A)(theta-theta_0)`, its normalized leading
Gaussian and that of the source desired-spin trace are respectively

```text
g_0=(q A/pi)^(1/4) exp(-T^2/2),
f_c=(q A c/pi)^(1/4) exp(-c T^2/2).                              (6)
```

The second width follows from `b_0 y_theta^2=A c`. Before applying `N`,
the desired component has norm tending to `sqrt((1+c)/(2c))` and the
wrong component to `sqrt((1-c)/(2c))`; their squared difference is one,
as required by flux. Normalization (3) therefore gives the following
`L2(dtheta)` limits, uniformly in the compact contour/frequency family:

```text
source incoming       = f_c       + O(q^(-1/2)),
source outgoing       = tau f_c   + O(q^(-1/2)),
receiver desired out  = exp(-Xi L)[f_c + O(q^(-1/2))],
receiver incoming     = exp(-Xi L)[-tau f_c + O(q^(-1/2))],
g_q                   = g_0       + O(q^(-1/2)).                  (7)
```

The receiver expressions use the reflected angular coordinate and the
desired-output phase specified above. They are statements about the
constructed trial, not the exact scattering solution.

Here is the ordinary error argument behind (7). The sphere and tube
coordinate maps are smooth uniformly near the compact family of central
intersections, and (4) is bounded away from degeneracy. Taylor expansion
of the exact normal exponent gives
`q b(s(theta)) y(theta)^2/2 = q A c (theta-theta_0)^2/2
+O(q |theta-theta_0|^3)`.
The intrinsic log-profile has the analogous expansion with `A`.
On a sufficiently small fixed angular neighborhood both exact exponents
are bounded below by a positive uniform constant times
`q(theta-theta_0)^2`. The mean-value formula bounds their exponential
difference by a constant times `q|theta-theta_0|^3` times a Gaussian.
Its normalized `L2` norm is `O(q^(-1/2))`. Smooth amplitude variations
and `exp(-Xi s(theta))` give the same order. The cubic and linear
correctors in (1) also have that order on the Gaussian scale. Outside
this neighborhood the cut-off tube trace and the intrinsic profile have
uniform exponential tails. For the latter, monotone angular flux gives
a unique maximum and a uniform logarithmic drop on the complementary
compact set; the endpoint powers suppress the polar endpoints. The
same estimates show its normalization constant differs from the
Gaussian one by `O(q^(-1/2))`. These bounds prove (7), including uniformity
when the admissible midpoint labels vary with `q`.

The normalized Gaussian overlap is

```text
gamma(c)=<g_0,f_c>=sqrt(2 sqrt(c)/(1+c)),
d_G(c)=||g_0-f_c||=sqrt(2-2 gamma(c)).                            (8)
```

Consequently the source desired-input mismatch tends to `d_G`, not
zero at fixed nonzero tilt. The receiver unwanted incoming amplitude
tends to `|tau| exp(-Xi L)` and must be included in the boundary budget.
For the *projected* desired receiving amplitude the error tends to
`(1-gamma) exp(-Xi L)`, smaller than the full profile distance `d_G`.
Small-angle expansion gives
`gamma=1-a^4/64+O(a^6)` and `d_G=a^2/(4 sqrt(2))+O(a^4)`.

## Conditional prepared-return consequence

Retain the boundary realization hypothesis of the weighted-current
criterion in the interior note. In particular, suppose the actual
solution with source input `g_q`, zero receiver incoming data, and no
incoming field from infinity has the required `L2` current/trace
realization. For `alpha=1/2`, its trial incoming mismatch satisfies

```text
B <= sqrt[d_G(c)^2+tau^2]+O(q^(-1/2)),
F_D -> 0,
delta_out <= 1-gamma(c)+O(q^(-1/2)).                             (9)
```

To justify the receiver weight in (9), its limiting value is
`exp(alpha Xi h_R)`, where `h_R` is the distance of the receiving
intersection from the source ball. The contour from the emitting
intersection has length `L>=h_R`. Hence the weighted unwanted input
factor `exp(-Xi L+alpha Xi h_R)` is at most one. Smooth weight variation
across the shrinking trace gives `O(q^(-1/2))`. At the source the distance
weight is exactly one.

A simple geometry bound makes (9) useful. The remote monopole field has
magnitude at most `k` times the local monopole field, with

```text
k=rho^2/(1-rho)^2,
|tau| <= k/[2(1-k)],
d_G <= [k/(1-k)]^2.                                             (10)
```

Indeed the radial field is at least `1-k` times the local magnitude,
the tangential part at most `k`, and
`|tau|=|B_t|/(|B|+B_r)`. Also `|tan a|<=k/(1-k)` and
`1-c<=tan(a)^2/2`. Writing `t=sqrt(c)` gives
`d_G^2=2[1-sqrt(1-(1-t)^2/(1+t^2))]<=2(1-t)^2`, which proves the
stated sufficient bound.

For any **fixed** `0<rho<=1/200`,

```text
|tau| <= 1/79200,
d_G <= 1/1568160000,
sqrt(d_G^2+tau^2) < 2e-5,
1-gamma=d_G^2/2 < 1e-18.                                       (11)
```

Fix also `delta` and the window `0.05<=xi D<=4`, `D=1-2rho`.
Equations (2), (7) and (11) imply that some finite `q_*(rho,delta)`
exists above which all retained sectors satisfy the earlier budgets
`B<=5e-5`, `F_D<=2e-5`, `delta_out<=1e-6`. Under the boundary realization
hypothesis, the projected exact exterior return therefore satisfies

```text
||R_trim(xi)-diag(exp(-xi L_j))|| <= .001 exp(-xi D).             (12)
```

This is a conditional ordinary *existence-of-threshold* consequence.
No numerical value of `q_*` or physically compatible joint scaling is
certified. The supplementary tests below use larger mouth radii to make
tilt visible; they do not certify (12) at a finite flux. The complete
untrimmed return and curved throat are not covered by (12).

The number of excluded midpoint labels is at most `2 delta q+2`.
This counts the exceptional channels; it does not establish their
transport. A later energy comparison may bound their contribution by
passivity and rank. It must actually prove that comparison and retain
the higher mouth modes and throat feedback before claiming a full result.

## Boundary literature and verification limits

The additional reading of [Ginoux and Murro v6](https://arxiv.org/html/2007.02544v6)
covers sections 2.1, 3 and 4.1-4.3, supplementing the earlier sections
2.2 and 5. The local energy argument and compact exhaustion support the
intended application. However, Proposition 4.6 imports Rauch's 1985
Theorem 8 for weak-to-strong boundary regularity, and the noncharacteristic
part of Theorem 4.8 imports Rauch--Massey (1974), Theorem 3.1. These
proof dependencies have not been read. The source's weak identity with
test fields vanishing at both time endpoints does not by itself encode
initial data, so that identity alone is not used to infer a causal
semigroup or a boundary trace domain.

The author's [course page](https://sites.lsa.umich.edu/jbrauch/course-material/)
links the 1985 article, but its migrated PDF link returned 404; the
tested AMS PDF routes returned 403. Retrieval outcomes are retained
outside Git. A complete application still needs the initial/boundary
trace argument and its positive-imaginary-frequency limit. This is a
specific remaining proof audit, not a claim that the standard existence
theory is false. No new boundary-existence theorem is claimed here.

Run `python checks/global_mouth_trial.py --output
evidence/global-mouth-trial.json`. The [receipt](../evidence/global-mouth-trial.json)
contains eight exact spin/flux/coordinate identities, nine actual dipole endpoint
geometries, eighteen two-mouth/frequency cases at three admissible fluxes,
and fifteen exceptional-channel counts. An independently integrated
field-line ODE supplies Fermi coordinates. The exact intrinsic angular
profile is compared with the transformed corrected tube trace at both
mouths; both ODE and angular quadrature resolutions are refined.
The observed errors to the analytic Gaussian limits decrease at the
predicted upper order. These are supplementary trace tests, not a full
exterior PDE solve or a proof of interval bounds on all geometry constants.

The first execution exposed an empty-mask call to the ODE interpolator,
which was corrected. The initial 240/320-point quadrature pair failed
the prescribed refinement tolerance at the transverse cutoff; raising
the pair to 480/640 passed without weakening the tolerance or the
convergence test. The final receipt pins the executed source bytes.
No Lean invocation occurred. This construction applies established
localized-mode and current-estimate methods; its eventual contribution
still requires the full model and novelty comparison.

Remaining to-do list: complete the boundary source/domain audit and make
the large-flux threshold quantitative where physically needed; control
exceptional channels in the energy comparison and match the curved
throat/state/holonomy; establish the physical consequence and novelty,
then certify suitable endpoints on GCP and prepare the manuscript.
