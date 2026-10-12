# Renormalized spherical Dirac response and the clock obstruction

2026-10-12. This is the constraint-gate endpoint for the symmetric MMP
throat-clock question. One covariant prescription supplies the retarded
lowest-Landau-level stress, the massive-tower dilaton source, and the
electric current. Their Ward identities hold at that order. The sourced
Misner–Sharp mass and electric Gauss flux integrate through both mouths
on one normalized asymptotic clock. The symmetric proper-mouth clock is
not a gauge-invariant dynamical boundary mode.

No stability sign, pole, shortcut, formation, propulsion, occupancy, or
engineering claim is made.

## Prescription

The spherical Hadamard–Polyakov scheme is the renormalized effective
action of the spherically reduced Dirac field after the canonical
rescaling `Psi = R^{-1} chi` of `dirac-spherical-source-reduction.md`.
It has three pieces, subtracted together:

1. The `q` massless two-dimensional Dirac fields, central charge `c=q`.
   Their stress is the Christensen–Fulling conformal stress in the frame
   `ds^2=-e^{2 omega} du dv`, including the characteristic state maps.
   This piece does not depend on the area `X=R^2`.
2. The local Hadamard finite part of each massive angular mode, taken
   from the zeta function of the squared Dirac operator at constant mass
   and then restored as a covariant function of `X`. Nonlocal massive
   propagators are omitted and estimated below.
3. The anomaly-free vector current of the massless level, including the
   center-of-charge Schwinger response. Its mass matches Sadhukhan's
   `m_S^2=g_c^2 q/(4 pi^2)` on the sphere, so the worldsheet mass is
   `m_S^2/X`.

The same subtraction defines the static reference, the retarded
variation, and the moving-join Schwarzian contact. A prescribed
Gaussian clock kernel is not an input.

## Stress, dilaton source, and current

For total central charge `c`,

```text
T_uu = (c/(12 pi)) (partial_u^2 omega - (partial_u omega)^2) + t_+(u),
T_vv = (c/(12 pi)) (partial_v^2 omega - (partial_v omega)^2) + t_-(v),
T_uv = -(c/(12 pi)) partial_u partial_v omega.                    (1)
```

The state functions `t_pm` are the retarded pullbacks

```text
t_+(u) = R_+'(u)^2 t_+^0(R_+(u)) - (c/(24 pi)) Sch(R_+,u),
```

and likewise for `t_-`, with `Sch` the Schwarzian. On a fixed RN collar,
`omega=(1/2) log F` reproduces equations (2) of
`rn-collar-lll-anomaly-flux.md`, including `T^a_a=c R_2/(24 pi)`.
The moving-interface composition is that note's Schwarzian chain rule:
the join contact enters `t_+` and then the mass flux, and it does not
add a homogeneous metric coefficient.

A chiral vector current `j_u(u)`, `j_v(v)` is covariantly conserved on
the same frame. With `epsilon_uv=+sqrt(-g)` and the Maxwell convention
of `sourced-spherical-mass-charge-balance.md`,

```text
partial_u Q_e = -(g_c^2/(4 pi)) j_u,
partial_v Q_e = +(g_c^2/(4 pi)) j_v.                          (2)
```

The center-of-charge mode uses `q` coherent flavors,
`j^a=(q/pi) epsilon^{ab} partial_b chi` with `Box chi` equal to the
worldsheet electric field, up to the orientation that makes the photon
mass positive. That mass is `m_S^2/X` with
`m_S^2=g_c^2 q/(4 pi^2)`.

Each massive mode has `m_n=lambda_n/sqrt(X)`, `lambda_n=sqrt(n(n+q))`,
and multiplicity `N_n=2(q+2n)` as written in the mode-reduction action.
The zeta function of one two-dimensional Dirac operator gives the
Euclidean density `m^2 (log(m^2/mu^2)-1)/(4 pi)`. At vanishing area
gradients the dilaton source is

```text
J_X = - (G N_n lambda_n^2 /(pi X^2)) log(lambda_n^2/(mu^2 X)).  (3)
```

If the local effective Lagrangian also depends on `(grad X)^2`, `J_X`
is the Euler–Lagrange derivative, not the explicit partial derivative.
The metric stress is the Hilbert stress of that Lagrangian.

## Ward identities at this order

Diffeomorphism invariance of the effective action gives

```text
nabla_a j^a = 0,
nabla^a T_ab = (J_X/(4 G)) partial_b X + j^a F_ba.             (4)
```

The checker builds the Christoffel symbols from the conformal metric,
builds `T_ab` from (1) and from the Hilbert variation, builds `J_X`
from the Euler–Lagrange derivative, and simplifies the residuals.
They vanish identically. The trace identity `T^a_a=c R_2/(24 pi)` is
the same simplification. Equation (4) is not imposed by setting a
residual to zero in advance.

The fractional size of the omitted nonlocal massive condensate on the
product throat is the curvature-to-mass ratio `R_2/m_1^2=2/(1+q)`.
It changes the value of `J_X`. It is not an extra homogeneous solution
of the source-free dilaton equation.

## Both mouths, one normalized clock

On each RN collar the balance note and the collar flux note give

```text
partial_u m = -(T_uu+T_uv),
partial_v m = +(T_vv+T_uv).                                   (5)
```

For a state difference on a fixed collar, the anomaly terms cancel and
the constant-time slice from a cut `r_b*` out to that mouth's infinity is

```text
M_i(t) - m_i(t,r_{b,i}*)
  = integral_{r_{b,i}*}^{infinity}
      [t_{+,i}(t-r*)+t_{-,i}(t+r*)] dr*.                      (6)
```

The same split with the orientation (2) integrates `Q_{e,i}`. The left
mouth uses its own outward tortoise. Both charts share the Killing time
`t` fixed by `N -> 1` at each infinity and by a fixed origin. A trial
reparametrization `t+alpha(t)` multiplies the asymptotic lapse by
`1+alpha'`. Normalization forces `alpha'=0`, and the fixed origin forces
the remaining constant to vanish. That constant is pure gauge.

The checker compares (6), and its charge analogue, with an independent
Simpson integral of compact Gaussian pulses on both mouths, and compares
the time derivative of the collar energy with the cut fluxes. The
maximum mass and charge residuals on those samples are about `2.3e-10`.
The rate residual is about `9e-9`.

## The symmetric clock does not survive

On the unit throat, the source-free dilaton solutions are

```text
Phi = a rho + sqrt(1+rho^2) (b cos tau + c sin tau).          (7)
```

Around the supported slope `Phi_0=beta rho`, the linear Casimir
variation is `delta C=-2 beta a`. The RN mass shift is

```text
delta M = r_e delta C / (8 G) = - r_e beta a / (4 G).         (8)
```

With `beta=2 r_e/ell` this is the static tangent already matched to
MMP's per-mouth mass. The coefficients `b` and `c` are the Lie
derivatives of `Phi_0` along the two non-translational Killing fields
of `AdS2`. Their asymptotic time components are opposite at the two
ends and are nontrivial functions of `tau` unless `b=c=0`. Holding
either normalized clock fixed sets `b=c=0`. The mass integral sets `a`.
No continuous coefficient remains.

Thus every retarded spherical source that obeys (4) changes the physical
data only by changing `M_pm(t)` and `Q_{e,pm}(t)` through (2) and (6).
The proper-mouth history is not an independent boundary degree of
freedom. It is fixed by those fluxes and by the clock normalization,
or it is the pure-gauge Killing frame that the normalization removes.
The prescribed Gaussian kernel is a response to a metric family outside
this solution space. It is not a stability operator.

The massive Schwinger mode remains what it already is in MMP Appendix A
and in Sadhukhan's axial-metric paper: a charge mode with mass
`m_S^2`, not a gravitational clock. This closure does not compute its
scattering and does not claim a new bound on it.

## Prior-art verdict

**Already known.** Editions pinned on 2026-10-12:

- Maldacena–Milekhin–Popov, arXiv:1807.04726v3, section 5.3.2 (static
  spherical support and RN mass match), section 5.4 (Schwarzian
  validity window and the binding-energy horizon warning, not a
  retarded radial eigenmode), and Appendix A (massive center-of-charge
  mode, massless relative sector, Wilson line fixed at the Casimir
  minimum).
- Sadhukhan, arXiv:2609.33511v1. For `l=0` the mouth reduces by the
  generalized Birkhoff theorem to a change of mass plus gauge, and the
  throat Jackiw–Teitelboim dilaton obeys `Box K=2K`, does not propagate,
  and is left untreated.
- Sadhukhan, arXiv:2610.09847v1. At `l=0` the metric amplitude vanishes
  and the sector is the spherically symmetric electric perturbation.
  The paper states that the only `l=0` mode is the massive Schwinger
  mode, and quotes MMP Appendix A for the fixed Wilson line.
- Maldacena–Qi, arXiv:1804.00491v3. The coupled-Schwarzian mode is a
  pair of boundary reparametrizations of rigid nearly-`AdS2`, with
  `SL(2)` imposed as a gauge constraint. The boundaries are cutoff
  surfaces, not two asymptotically flat normalized clocks. Asymptotic
  flatness appears there only inside the topological-censorship
  citation.

Two searches on 2026-10-12, one for an MMP `l=0` retarded stability
kernel and one for an asymptotically flat spherical clock
reparametrization, return these families and unrelated thin-shell,
Euclidean-axion, and echo calculations. Absence of a further hit is
not a novelty certificate. The unread remainder is the internal
bosonization proof of MMP Appendix A, which was not re-derived in this
session; Sadhukhan v1 and the existing satellite audit quote the same
massive, Wilson-line-fixed conclusion, so that gap does not reopen a
massless clock. The bulk of Maldacena–Qi beyond the boundary action and
the `SL(2)` gauge constraint was not re-read; the construction does not
use normalized asymptotically flat mouths.

No same-emission, same-receiver exterior-only clock comparison is
proved. MMP's leading no-shortcut comparison stays the predecessor.

## Formalization boundary

No theorem is asserted as Lean-certified. The argument is the ordinary
identity and mode count checked by
`checks/retarded_dirac_spherical_response.py`. No GCP instance was
started, and no local Lake or Lean run is acceptance. A local run was
not performed.

## Check

Run `checks/retarded_dirac_spherical_response.py`. It is the evaluator
of (1)–(8), not a second implementation hidden behind the note. The
receipt is `evidence/retarded-dirac-spherical-response.json`.

Remaining to-do list: none for this endpoint. Publication, DOI minting,
journal submission, and outreach are unselected.
