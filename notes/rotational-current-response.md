# Rotational current, the dipole constraint and a positive reduced energy

2026-10-11, E003/S053. The stationary rotation constant must be coupled
to the fermion angular current before it can be used as a quantum spin
datum. This note performs that coupling for the **defined quadratic
Einstein-Maxwell plus free projected lowest-Landau-level (LLL) action**.
It eliminates the dipole gauge constraint without dividing by frequency,
obtains a positive two-field energy on a prescribed spherical geometry,
and computes the two product-throat masses, including finite-flux factors.

This is a component of the retained physical long-return/finite-transfer
investigation. It does not supply the background Casimir-state Hessian,
the omitted Landau levels, the interacting state, or a smooth quantum
wormhole. In particular it is not a physical stability theorem. Rotation
screening itself was already predicted explicitly by
[Maldacena, Milekhin and Popov (MMP), Appendix A.1](https://arxiv.org/html/1807.04726v3#A1.SS1).
No consequential novelty is accepted.

## Fields, normalization and the classical action

Use magnetic-radius units: the dimensionless areal radius is R(r),
coordinates are dimensionless, and the physical metric is P^2 times

```text
ds^2 = -S(r) dt^2 + dr^2/S(r)
       + R(r)^2 [dtheta^2 + sin^2(theta)(dphi - w_a dx^a)^2],
A_geom/P = -cos(theta) dphi + cos(theta) V_a dx^a,
A_charge = (q/2) A_geom/P,      K=P^2/G,      P^2=pi G q^2/g_c^2.
```

Here S,R>0 are prescribed static functions and a=t,r. The microscopic
fermion has unit charge under A_charge. Fix the polar Maxwell angular
gauge. The remaining small rotation changes w and V by the same exact
one-form; Phi=V-w is invariant. The square in the metric specifies the
second-order path. A matter Hessian on a different path can differ by
the background stress tadpole; see the existing
[Ward audit](casimir-cylinder-ward-completion.md).

Let E_w=dot(w_r)-w_t', E_V=dot(V_r)-V_t'. Reduction of the classical
action (16 pi G)^(-1) integral sqrt(-g_4)(R_4-F_geom^2) gives

```text
I_cl = integral sqrt(-g_2) [ -A_w F_w^2/4 - A_V F_V^2/4
                            -D Phi_a Phi^a/2 ],
A_w=K R^4/6,       A_V=K R^2/3,       D=2K/3.                  (1)
```

The gravitational kinetic term follows from the sphere Killing vector:
the curvature correction is -R^2 sin^2(theta) F_w^2/4. The integrals
of sin^2(theta) and cos^2(theta) are 8pi/3 and 4pi/3. In the horizontal
frame the Maxwell field has angular-longitudinal component proportional
to sin(theta) Phi, which supplies D. The checker independently evaluates
the full four-dimensional curvature and Maxwell action on constant-radius
test fields with nonzero curl, determining all three coefficients.

It also derives every linear Einstein-Maxwell component with arbitrary
time-dependent w_t,w_r,V_t,V_r on both exact RN (R=r,
S=(r-1)^2/r^2) and product AdS2 (R=1,S=1+r^2). Their equations agree
with (1), including the angular Maxwell constraint. At zero frequency
and with radial potentials zero they reproduce

```text
(R^2 V_t')' = 2(V_t-w_t)/S,
(R^4 w_t')' = -4(V_t-w_t)/S.                                 (2)
```

The restricted classical action makes sense for arbitrary S,R. Agreement
with a complete perturbation problem on a non-electrovacuum background
also requires the source's action and its variations. We do not replace
those variations by (1).

## The projected fermions and the correct conserved flux

There are q massless two-dimensional fermions in spin j=(q-1)/2.
In the magnetic orientation used here the normalized mode density obeys
<cos(theta)>_m=-2m/(q+1). Indeed its density in u=cos(theta) is
proportional to (1+u)^(j-m)(1-u)^(j+m); the beta-integral mean is
[(j-m)-(j+m)]/(2j+2). The charge vertex and the metric vertex,
including the spin connection, combine to give the effective potential

```text
A_eff,m = -m a,       a=(w+q V)/(q+1)=k_w w+k_V V,
k_w=1/(q+1),         k_V=q/(q+1).                            (3)
```

This follows directly from the previously checked
[full Dirac vertex](finite-flux-response-audit.md): the metric term has
expectation -m w/(q+1), and the charge term -qm V/(q+1). When V=w,
the answer is -m w, as required by a local rotation. The curl term of
the four-dimensional Dirac operator changes angular chirality and has
zero LLL-to-LLL matrix element; its virtual higher-level contribution
has not been integrated out here.

For a fixed rotation axis the vertices commute. Abelian bosonization of
the q free projected fermions leaves one normalized current field chi:

```text
T_q=sum_m m^2=q(q^2-1)/12,
gamma=sqrt(T_q/pi),
I_chi=integral dt dr [dot(chi)^2/(2S)-S chi'^2/2
                     +gamma chi(k_w E_w+k_V E_V)].           (4)
```

The sign of chi has been chosen to give the displayed coupling. The
orthogonal q-1 bosonic combinations do not couple to this one axis.
The electric singlet is orthogonal because sum_m m=0. They have not
been removed from the underlying theory or its Casimir energy. Three
axes have the same two-point normalization by rotational covariance;
they are not three independent Abelian theories at nonlinear order.
The non-Abelian current algebra, compact sectors and global spin lift
must be retained in a nonlinear/global treatment.

Equations (3)-(4) specify the nonzero-mode free-LLL current response.
They do not choose a Green function, a boundary holonomy, or a quantum
state on the actual closed exterior/throat trajectories. q=1 has T_q=0
and no rotational current field; do not retain a spurious chi in that
case. The rest of this note takes integer q>=2.

The canonical electric displacements are

```text
p_w=A_w E_w+gamma k_w chi,
p_V=A_V E_V+gamma k_V chi,
Q=p_w+p_V,             Q'=dot(Q)=0.                          (5)
```

Both the spatial and temporal constancy follow from the two common
rotation equations; no frequency division is used. For stationary
fields define the geometric spin expression

```text
J_geom=-(R^4 w_t'+2 R^2 V_t')/6.
Then Q=K J_geom+gamma chi.                                   (6)
```

In electrovacuum Q=KJ is the spin normalization calibrated to slow
Kerr-Newman in the [stationary audit](stationary-electric-dipole.md).
With fermions, the geometric expression alone is not constant. (6) is
the local angular Gauss law with its matter contribution. The value of
Q and its interpretation at the two physical ports require boundary
data; we have not identified it with the ADM angular momentum of the
complete rotating binary. Treating each mouth as an isolated quantum
rotor would omit this coupling.

## Constraint elimination and a positive energy

Set p=p_V, p_w=Q-p and pi_chi=dot(chi)/S. Up to the displayed Gauss
constraint and the boundary Legendre term, eliminating Phi_t gives

```text
H = integral dr [ (Q-p-gamma k_w chi)^2/(2 A_w)
                  +(p-gamma k_V chi)^2/(2 A_V)
                  +S(pi_chi^2+chi'^2)/2
                  +S p'^2/(2D)+DS Phi_r^2/2 ].              (7)
```

Indeed Phi_t=-S p'/D. The momentum term is p dot(Phi_r), apart from
the common gauge pair. The transformation
y=p/sqrt(D), pi_y=-sqrt(D) Phi_r is canonical up to a total time
derivative; D is constant in this normalization. Each term in (7) is
nonnegative. In a fixed-Q sector the quadratic perturbation energy is
the Q=0 Hessian even if the stationary solution has Q nonzero.
Boundary work must be included when the prescribed potentials vary.

Writing z=(chi,y), the homogeneous reduced equations are

```text
ddot(z)/S - (S z')' + M(r) z = 0,
H_2=1/2 integral dr [S(|pi|^2+|z'|^2)+z^T M z],

M11=gamma^2(k_w^2/A_w+k_V^2/A_V),
M12=gamma sqrt(D)(k_w/A_w-k_V/A_V),
M22=D(1/A_w+1/A_V),
det M=D gamma^2/(A_w A_V)>0.                                (8)
```

The determinant uses k_w+k_V=1. Thus M is pointwise positive for any
R>0. The local energy law is d_t e=d_r[dot(z)^T S z']. For a finite
interval with homogeneous Dirichlet conditions, the energy is conserved
and positive; an exponentially growing nonzero finite-energy solution
is impossible. The same conclusion holds for other well-posed domains
whose boundary flux supplies no energy, including energy-continuous
piecewise matching. On an unbounded domain positivity need not be
uniform; no spectral gap or uniform decay follows. This is a theorem
for (1),(4), not for the uncomputed full semiclassical source.

In tortoise coordinate dx=dr/S, the operator is
-d_x^2+S M, with the common radial light cone. In particular no extra
massless gravitational dipole propagates. One classical electric/axial
degree of freedom and the one current field give two local scalar
degrees of freedom. A constant Q and global holonomy are separate data.

Define h=gamma^2/K. The explicit matrix is

```text
M11=3h[2/R^4+q^2/R^2]/(q+1)^2,
M12=sqrt(6h)[2/R^4-q/R^2]/(q+1),
M22=4/R^4+2/R^2,       det M=12h/R^6.                        (9)
```

This also supplies the leading projected current model in the RN mouth.
Its light eigenvalue tends to zero as R^(-4) at large R. A mass inferred
at the throat center is not a global gap or a reason to remove boundary
spin/holonomy data.

## Product-throat spectrum and predecessor checks

For R=1 let c=(q-2)/[3(q+1)] and

```text
s0=g_c^2 q/(4pi^2),      h=(s0/3)(1-q^(-2)),
M=[[ (2+9c^2)h, -c sqrt(54h) ],[-c sqrt(54h),6]],
B=6+(2+9c^2)h,
mu_minus^2, mu_plus^2 = [B -/+ sqrt(B^2-48h)]/2.              (10)
```

These are positive for h>0. At small h,
mu_minus^2=2h+O(h^2), mu_plus^2=6+9c^2 h+O(h^2).
For h increasing without bound at fixed q, the light eigenvalue tends
to 12/(2+9c^2); it does not grow without bound. For q=2, c=0 and the
two eigenvalues are simply 2h and 6, with a crossing at h=3. These
large-h statements test the defined quadratic model, not perturbative
control of the microscopic theory at strong coupling.

As q increases at fixed s0, its characteristic polynomial becomes

```text
x^2-(6+s0)x+4s0.                                           (11)
```

This equals the nonzero polynomial obtained by algebraically setting
l=1 in [Sadhukhan, (8.4)](https://arxiv.org/html/2610.09847v1#S8).
That paper explicitly restricts its reduction to l>=2 and leaves the
dipole gauge constraint untreated (section 2.2). We use (11) as a
limiting cross-check; the third zero eigenvalue of its three-variable
formula is not a physical third dipole field. No claim that its l>=2
neutral channel disappears is made. MMP Appendix A.1 already identifies
rotational screening at order G q^3/P^2, equivalent to g_c^2 q.

The exact finite-flux trace and determinant are

```text
tr M = 6+s0 (q-1)(q^2+2)/[q^2(q+1)],
det M=4s0(1-q^(-2)).                                       (12)
```

In the retained g_c^2=q^(-2) family,
mu_minus^2=1/(6pi^2 q)+O(q^(-2)). This does not validate the earlier
operator-norm comparison across one-unit spin sectors, and it does not
prove that those sectors fluctuate or destabilize the throat. It says
that a dynamical angular-current response must replace a fixed
electrovacuum rotation parameter in this approximation.

On the ideal global product throat with the standard fast-decaying
boundary condition, the independent scalar frequencies obey
Omega_n=n+Delta, Delta=(1+sqrt(1+4mu^2))/2. Other boundary choices,
including the actual two-mouth load, need not have this spectrum.
The numerical check solves the original coupled two-field operator on
refined grids and compares with these analytic frequencies; it does
not use them as a physical wormhole boundary condition.

## What remains before a physical conclusion

The free LLL projection is exact only as a defined sector. Missing pieces
include virtual higher levels and contact terms, renormalized couplings,
the Casimir-state stress response on the specified metric path, internal
photon/graviton interactions, nonspherical/time-dependent matching, and
the non-Abelian/global spin sectors. Gauge and energy identities of the
truncated action cannot bound those omissions. Earlier finite-flux and
Ward audits remain applicable. An isolated-RN rotor partition function
is not a substitute for this coupled wormhole state.

The concrete next step is to put this angular-current sector into the
same physical exterior/throat return and state construction used for
the electric singlet, and determine its boundary charge and holonomy.
Recompute its source response, rather than add a central constant spin
shift by hand. The screening mechanism and its scaling are prior work;
the explicit constraint/finite-flux completion is supporting evidence,
not an accepted consequential endpoint or a new Lean campaign.

[The checker](../checks/rotational_current_response.py) and
[receipt](../evidence/rotational-current-response.json) retain independent
curvature, full dynamic tensor, angular projection, Hamiltonian and
spectrum checks. [The source review](../evidence/rotational-current-source-review.json)
records the searches, editions and bounded comparison. No Lean, cloud,
outreach, publication or DOI action was undertaken.

Remaining to-do list: derive the physical rotational boundary/state and
holonomy response together with the electric singlet; control omitted
quantum terms and smooth semiclassical matching; finish the consequential
physical endpoint and novelty review, independently verify, GCP-certify
suitable mathematics and prepare the reproducible manuscript.
