# AdS formation wall: matter response before a mode count

S070, 2026-10-12. **The minimal material-sheet model has negative angular
directions, but its instability cannot yet be assigned to the microscopic
formation process.** An explicit collisionless comparison has the same
spherical stress and reduced radial action with a different angular response.
This resolves a modeling ambiguity in the next calculation; it establishes
neither a wormhole creation rate nor stability of a formed wormhole.

## Source and scope

[Bintanja, Freivogel and Rolph (BFR), 2308.00871v2](https://arxiv.org/html/2308.00871v2)
assume a radiation wall. We have now inspected both ensemble calculations,
the perfect-fluid appendix D and radial-fluctuation appendix F. Their
discussion leaves its microscopic origin open. Appendix F's extra radial
negative modes near extremality are already known.

The following results specify two **probe matter models** on a fixed static
geometry. The material fluid calculation uses an independent embedding and
comoving fluid fields. The ray calculation uses independent particle phase
space. Both agree on the spherical trajectory. This agreement does not make
their nonspherical constitutive laws interchangeable.

The underlying comparisons are established physics: the fluid action of
[Dubovsky et al., sections II/III](https://arxiv.org/html/1107.0731v1),
the material-membrane instability of
[Yang, Bonga and Pan](https://arxiv.org/html/2207.13754v2), and null-orbit
Lyapunov growth as in
[Cardoso et al., section III.2](https://arxiv.org/html/0812.1806v3).
No novelty acceptance is claimed. We use the last source's local geodesic
analysis, not its asymptotically flat quasinormal-mode correspondence in AdS.

## 1. Specify a material action, rather than just an equation of state

On the 2+1 dimensional worldvolume take two comoving labels `phi^I`, induced
metric `gamma_ab`, and density `b=sqrt(det B)`, with
`B^IJ=gamma^ab partial_a phi^I partial_b phi^J`. On a sphere the equivalent
global definition pulls back a reference entropy area form. The minimal
non-dissipative radiation action is

```text
S_mat = -integral sqrt(-gamma) rho(b) d^3 y,
rho(b)=c b^(3/2),       p=b rho'(b)-rho=rho/2,
c_sound^2=dp/drho=1/2,  c>0.
```

It conserves the comoving entropy and does not freeze the surface density
when the embedding changes. That distinction is essential when varying
the angular shape.

In a locally flat patch set the normal displacement to `z(t,x,y)` and hold
the comoving labels at their equilibrium values. Directly from the induced
metric `gamma=diag(-1,1,1)+dz dz^T`,

```text
-det gamma = 1+|grad z|^2-zdot^2,
b^2/b0^2 = (1-zdot^2)/(1+|grad z|^2-zdot^2),
L2 = (rho0/2) zdot^2 + (p0/2)|grad z|^2.
```

The kinetic energy is positive. The gradient term has the instability sign,
giving `omega_proper^2=-(p0/rho0) k^2=-k^2/2`. After Wick rotation the
quadratic Euclidean form is
`[rho0 (partial_tau z)^2-p0 |grad z|^2]/2`. These are matter displacement
directions in the probe model, not gravitational conformal-factor modes.

## 2. Exact normal Hessian on the static spherical probe

Let `ds^2=-f dt^2+dr^2/f+r^2 dOmega^2`, `u(r)=sqrt(f)/r`, and let the
unperturbed radius `R` satisfy `R f'(R)=2 f(R)` with `f(R)>0`. Write
`rho=C0/r^3` on spherical configurations, `C0>0`. The Euclidean material
action for a general normal graph at fixed comoving labels has density

```text
L_E = (C0/r) X^(3/4) (X+Y)^(-1/4),
X=f+(partial_t r)^2/f,   Y=|grad_Omega r|^2/r^2.
```

This follows from the full induced determinant, including its mixed
time-angle components. For `r=R+chi`, expanding through quadratic order gives

```text
S_E2 = (C0/2) integral dt dOmega [
  (partial_t chi)^2/(R f_R^(3/2))
  + u''(R) chi^2 - |grad_Omega chi|^2/(2 R^3 sqrt(f_R)) ].
```

The normal/tangential compression mixing is proportional to `u'(R)` and
vanishes at this probe equilibrium. For example the static compression
Jacobian contributes `u(r) J^(3/2)`; its mixed `r,J` derivative is zero
there. Thus this is the normal block of the specified fluid model, rather
than a variation that imposes inconsistent frozen fluid pressure.

For `chi ~ exp(2 pi i n t/beta) Y_lm`, division by the positive kinetic
coefficient gives

```text
Lambda_nl = (2 pi n/beta)^2 - gamma_t^2 - f_R l(l+1)/(2R^2),
gamma_t^2 = -R f_R^(3/2) u''(R)
          = f_R^2/R^2 - f_R f''(R)/2.
```

When `gamma_t^2>0`, every static angular direction in this continuum model
is negative. Its high-`l` growth matches the local membrane result. A real
material cutoff and a gravitational constraint calculation are still needed
before counting physical instanton modes.

For extremal RN-AdS, define `a=bar_r`, `x=a^2/L_AdS^2`. Direct substitution
in `f=1+r^2/L_AdS^2-2a(1+2x)/r+a^2(1+3x)/r^2` gives

```text
R=2a(1+3x),
f_R=(1+2x)(1+6x)^3/[4(1+3x)],
gamma_t^2=(1+2x)(1+6x)^4/[32 a^2(1+3x)^4].
```

This reproduces BFR's general radial criterion F.5 directly. The retrieved
HTML and PDF text of F.8 instead contain this last expression divided by
`R`, which has the wrong length dimension for a squared frequency. The
checker records the discrepancy; it does not use F.8 for validation or
claim an author-confirmed correction. In the flat extremal limit we obtain
`gamma_t^2=1/(32a^2)`.

## 3. Same spherical action, different collisionless response

A massless particle with fixed angular momentum magnitude `J>0` has
coordinate-time Hamiltonian

```text
H(r,p_r,J)=sqrt(f^2 p_r^2+f J^2/r^2).
```

Legendre transforming only the radial momentum gives its radial Routhian

```text
L_J=-(J/r) sqrt(f-rdot^2/f).
```

An isotropic ensemble of orbital planes, with identical radial motion and
fixed `J_i`, therefore has exactly the same spherical action as the fluid
wall, after identifying `sum J_i=4 pi C0`. At a circular orbit its local
stress is `p_radial=0`, `p_theta=p_phi=rho/2`. It consists of rays that can
stream through one another; no material surface confines them after a
nonspherical perturbation. The massless Vlasov stress integrals in
[Andreasson, Fajman and Thaller, equations (2.8)-(2.10)](https://arxiv.org/html/1511.01290)
give this same tangential pressure in the circular-orbit limit. Their
self-gravitating existence theorem is not used as a stability theorem here.

Expand the one-particle Hamiltonian at `(R,0)`. Its radial quadratic form
is `a_p delta_p^2/2+b_r chi^2/2`, with

```text
a_p=f_R^2/H0 > 0,  b_r=J u''(R),  -a_p b_r=gamma_t^2.
```

Let `G_geo` generate geodesic motion on the unit tangent bundle of the
round sphere, with angular speed `Omega=sqrt(f_R)/R`, and define
`D_t=partial_t+Omega G_geo`. Since `u'(R)=0`, the angular speed has no
linear radial correction at fixed `J`. The radial phase-space perturbation
equations are

```text
D_t chi=a_p delta_p,   D_t delta_p=-b_r chi,
(D_t^2-gamma_t^2) chi=0.
```

On the growing branch,
`chi(t)=exp(gamma_t t) chi0 o Flow(-t)`.
The geodesic flow preserves the phase-space measure, so its pullback is
unitary on `L2` and the norm growth rate is exactly `gamma_t`, independent
of the angular scale of the initial perturbation. This is the norm of
displacement amplitudes with the invariant orbital measure, not an `L2`
norm of the delta-supported distribution itself. In a local streaming
description, `omega=k.v +/- i gamma_t`; an angular gradient adds advection,
not the fluid sheet's growth proportional to `|k|`.

This is an explicit difference in Lorentzian responses with identical
spherical data. It does **not** prove stability of the ray ensemble:
`gamma_t` is positive. Nor does it calculate its Euclidean determinant;
periodicity, ensemble and particle occupation constraints have not been
continued into a formation saddle. The ray ensemble does not remain a
barotropic perfect fluid under arbitrary nonspherical perturbations.

| Initial angular scale l, flat extremal example a=1 | Material-sheet growth in t | Collisionless phase-space norm growth in t |
| --- | ---: | ---: |
| 0 | 0.176777 | 0.176777 |
| 2 | 0.467707 | 0.176777 |
| 10 | 1.862458 | 0.176777 |
| 50 | 8.928536 | 0.176777 |

The table compares the two specified probe responses. It does not identify
collisionless phase-space norms with gravitational quasinormal modes.

## 4. Constraints and the unresolved physical scale

The local thermodynamic constraint does not by itself fix the material
sign. For `E=c S^(3/2)/sqrt(A_surface)`,
`dE=T dS-p dA_surface`. At fixed entropy, an area increase lowers energy;
at fixed energy, it increases entropy by `dS=(p/T)dA_surface`. Thus the
static Euclidean contribution `-dS` retains the negative area sign in
the local microcanonical calculation. This is not a replacement for the
full constrained gravitational Hessian.

For an actual ideal material membrane with a regular local derivative
expansion, the known high-wave-number gravitational correction is lower
order than the `k^2` pressure term. A potential controlled band would need

```text
max(R^-1, curvature_scale^-1, G rho) << k
    << min(thickness^-1, mean_free_path^-1, microscopic_cutoff).
```

Neither thickness nor a collision/transport scale follows from `rho~R^-3`.
No such band has been established for BFR's radiation state. A collisionless
gas has no fluid mean-free-path band, even if its background stress looks
isotropic. The matched Einstein-Maxwell-fermion state, rather than the
radial action alone, must determine which description applies. We have not
integrated out the bulk matter or solved the nonspherical Israel constraints.

## 5. Verification and current decision

[The checker](../checks/ads_wall_response.py) verifies the induced-metric
expansion, positive kinetic coefficient, pressure and sound derivative,
spherical normal Hessian, vanished compression mixing, fixed-angular-momentum
Routhian, independent Hamiltonian Hessian, RN-AdS substitution and local
fixed-energy identity. It also integrates the **nonlinear** material action
for `chi=epsilon P_l(cos(theta))`, at `l=0,2,5`, `a^2/L_AdS^2=0,0.1`, and
three step sizes. All 18 finite differences agree with the Hessian; the
largest final relative error is `2.85e-10` and shrinks quadratically.

Run `python checks/ads_wall_response.py` from the satellite root. The
workstation interpreter is
`C:/Users/dfred/QuantyraTools/SpaceProofs/python/python.exe`.
The [receipt](../evidence/ads-wall-response.json) pins the code and records
what was not computed. [Source provenance](../evidence/ads-wall-source-review.json)
pins the actual reading scope, hashes and limited novelty search.

**Decision:** the minimal material wall is a poor candidate for a saddle
with one negative mode. Applying that conclusion to a physical nucleation
process remains conditional. An off-equilibrium perfect-fluid constitutive
law selects the material case; specifying only its spherical moments does
not. The explicit ray comparison prevents an
unjustified rejection of all radiation-supported formation. S070 remains
active, with the next calculation narrowed by evidence to the radiation
state's phase-space response and allowed variations. The original goal
remains incomplete and no creation rate is accepted.

Remaining to-do list: construct the collisionless radiation response with
its occupation and boundary constraints, determine whether it supplies an
admissible Euclidean formation saddle, then include backreaction and test
finite-disturbance traversability if a candidate survives.
