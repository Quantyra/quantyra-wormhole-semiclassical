# The bulk velocity contact leaves an order-epsilon joined-geometry term

2026-10-11, E003/S053. Expanding the moving-cut throat Hamiltonian to
second order identifies its explicit `phidot^2` contact exactly. In the
free conformal model, this bulk term reproduces the finite circle-Casimir
piece of the already derived trajectory coefficient. The remaining
order-epsilon piece belongs to the joined geometry and its local stress
contact; it is not supplied by the bulk Hamiltonian or the static
determinant alone.

This closes one algebraic bookkeeping question and narrows the interacting
contact problem. It does not identify the full seam term for the charged
potential, bound the interacting ultraviolet remainder, or imply a
physical damping result.

## 1. Expand the cut condition and Hamiltonian

Keep the product-throat conventions of the [moving-mouth vertex](moving-mouth-clock-vertex.md):

```text
A=acos(epsilon), b=sqrt(1-epsilon^2), r=epsilon/b, x=exp(phi), J=s/A,
cos^2(s)=epsilon^2(x^2-sdot^2).
```

For a fixed smooth history, the unique positive-cut branch expands to
quadratic amplitude as

```text
s=A-r phi-epsilon phi^2/(2b^3)+epsilon^3 phidot^2/(2b^3)+O(phi^3).
J=1+kappa phi+lambda phi^2+mu phidot^2+O(phi^3),
kappa=-r/A, lambda=-epsilon/(2 A b^3), mu=epsilon^3/(2 A b^3).
```

The `phidot^2` term in `s` follows by putting `sdot=-r phidot+O(phi^2)`
inside the cut condition; its coefficient is not an independent slow-motion
approximation. The pulled-back Hamiltonian is

```text
H=1/2 integral_th [(x/J)(pi^2+f_z^2)
                    +a x J sec^2(Jz)|v*f|^2] dz
  +1/2 integral_th (Jdot z/J){pi,f_z} dz + H_ext.                 (1)
```

Here `H_ext` has no explicit bulk clock-velocity term in the prescribed
coordinate map. The coefficients of `phidot^2` in `x/J` and
`x J sec^2(Jz)` are respectively `-mu` and
`+mu sec^2(z)(1+2z tan(z))`.
The shift term in (1) is linear in `Jdot` through quadratic amplitude.
Thus the explicit bulk seagull in the Hamiltonian is

```text
H_bulk,02=-mu integral_th T_zz dz
           +mu integral_th z V0 tan(z)|v*f|^2 dz,
T_zz=(pi^2+f_z^2-a sec^2(z)|v*f|^2)/2.                            (2)
```

With the effective-kernel convention of the preceding free-response note,
this contributes `[-2mu <integral_th T_zz>_ren
+2mu <integral_th z V0 tan(z)|v*f|^2>_ren] w^2` to the quadratic
kernel. Equation (2) is an operator identity for this bulk coefficient;
it is not the full time-ordered response. The retarded vertex covariance,
anomaly, moving-cut boundary work and exterior matching must still be
included.

## 2. Independent free-field calibration

For one nonchiral conformal field on the joined circle, central charge c,
the renormalized static stress in the throat is

```text
e0=-pi c/(6 L^2),   T_zz=e0-c tan^2(z)/(24 pi),
L=T+D.
```

Therefore

```text
integral_-A^A T_zz dz=2 A e0-c/(12 pi)(b/epsilon-A).
```

Substitution in (2) gives the bulk velocity-contact coefficient

```text
B_bulk=-epsilon^3/(A b^3) integral_-A^A T_zz dz
       =pi c epsilon^3/(3 L^2 b^3)
        +c epsilon^2(b-A epsilon)/(12 pi A b^3).                  (3)
```

The independently derived complete free proper-cut coefficient is

```text
B_free=c epsilon/(12 pi b)+pi c epsilon^3/(3 L^2 b^3).            (4)
```

The first term in (4) is the universal local part of the joined
stress/contact calculation. The difference not supplied by the bulk
seagull is exactly

```text
B_free-B_bulk=c epsilon/(12 pi b)
              *[1-epsilon(b-A epsilon)/(A b^2)].                  (5)
```

For `0<epsilon<=1/4`, the bracket is between `11/15` and `1`: use
`A>1`, `b^2>=15/16` and `0<=b-A epsilon<=1`. Hence the joined-geometry
term is positive, order c epsilon, and independent of D. Equation (5)
is a calibration of the decomposition in the free CFT, not a derivation
of its interacting continuation.

For the massive charged Gaussian theory, (2) fixes both the full
renormalized pressure term and an additional potential-profile contact,
`mu integral z V0 tan(z)|v*f|^2`. The free calibration does not determine the remaining cut/join
contact or make its potential-dependent expectation small. The existing
finite-mouth static determinant controls clock derivatives along its
static proper-cut family; it does not by itself bound this independent
velocity-squared boundary variation.

## 3. Closest prior work and the ultraviolet consequence

Franchino-Viñas and Mazzitelli, [Effective action for delta potentials](https://arxiv.org/abs/2010.11144v2), derive a heat-kernel effective action for a different model: a thin flat codimension-one delta layer, perturbative in its spacetime-dependent inhomogeneity. Their sections 3.1-3.2 separate a local divergence, a finite local choice and a nonlocal form factor; section 4 treats the dissipative pair-production part. The divergence and subtraction depend on the underlying layer theory. This agrees with the need to derive both local matching and nonlocal response from one action, but it does not give formulas for our finite smooth 1+1 throat, charged rank-one potential, or transmitting exterior join.

Their local/nonlocal separation is already known physics, so it is not a novelty claim for this project. It also shows why the finite-band estimate cannot simply be extended by formal subtraction: a boundary counterterm belongs to the chosen microscopic model. The [canonical interacting response](singlet-complex-response.md) remains a different observable and does not control the endpoint-pressure form.

## 4. Verification, evidence and decision

The [checker](../checks/proper_mouth_seagull.py) verifies the cut-condition
series, Hamiltonian coefficient extraction, stress decomposition, exact
free calibration and bracket bound at high precision. The
[receipt](../evidence/proper-mouth-seagull.json) binds the script, notes
and cached source. The [source supplement](../evidence/proper-mouth-seagull-source-review.json)
records the pinned arXiv v2, reading scope and model differences. These
checks independently verify the algebra; they do not verify an interacting
boundary contact or an ultraviolet tail estimate.

Continue long-return damping as provisional primary and matched
finite-duration transfer as fallback. The fixed-band pressure theorem
survives this audit. Next derive the joined velocity contact for the
interacting action, then bound its renormalized full retarded response.
The local term's matching coefficient and the response's physical 4D
embedding remain unresolved; physical damping and consequential novelty
are not accepted.

Remaining to-do list: derive and bound the interacting joined contact and
ultraviolet response; close constraint-compatible geometry, state and
exterior; finish primary/fallback novelty assessment; independently verify,
certify suitable mathematics on GCP and prepare the manuscript.
