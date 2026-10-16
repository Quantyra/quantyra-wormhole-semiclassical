# A centered covariant comparison exterior

2026-10-15, S139/S074. Ordinary derivation before focused verification.
There is an explicit static spherical comparison family with the exact
homogeneous quantum center, its independently determined AdS2 radius,
the corrected outer vacuum curvature, and the prescribed positive
temperature. A scalar potential in the reduced spherical action
generates this family exactly. Adding and subtracting that potential
preserves the retained theory and its stress/radius Ward identity.
This constructs a comparison geometry, not the corrected Dirac exterior
or a justified inverse for the full one-loop equations.

## Fixed data and two independently corrected radii

Keep the physical magnetic parameter Q, G, lambda0=L4^-2, the selected
finite coefficients and mu0=1/b0 fixed. In particular mu0 is not reset
to a quantum radius. Let B=b_star and A=a_star be the exact product
radii of [the homogeneous theorem](uniform-dirac-reference.md), in its
one-percent box. Select the positive outer-vacuum curvature lambda_e by

```
lambda_e+11G lambda_e^2/(360pi)=lambda0.                 (1)
```

This is the specified AdS-invariant vacuum branch of
[the outer-domain audit](stationary-exterior-conformal-domain.md),
not proof that an actual monopole HH exterior approaches that vacuum.
The reference construction assumes this outer prescription.

Put t=r/B, x=lambda_e B^2, C_star=B^2/A^2,
C_RN=1+6x and d=C_star-C_RN. No classical RN relation between A and B
is imposed. The accepted product bounds and (1) imply

```
|d|/C0<5e-26<1e-24,   |d|<C_RN,
C0=sqrt(1+12pi k^6),   k=2000.                         (2)
```

For clarity, C_star/C0=z_star/u_star, whereas
C_RN/C0=C0^-1+(1-C0^-1)(lambda_e/lambda0)z_star.
Subtracting these expressions bounds |d|/C0 by
z_star|u_star^-1-1|+|z_star-1|/C0
+(1-C0^-1)z_star|1-lambda_e/lambda0|.
Here |u_star-1|<4e-26, |z_star-1|<1.4e-22, C0>6k^3,
and |1-lambda_e/lambda0|<G lambda0=k^-10; these give (2).

## Exact centered and thermal metrics

Define

```
F_ext(t)=(t-1)^2/t^2
            [1+x(t^2+2t+3)+d/t^2],
F_delta(t)=F_ext(t)-delta/t,
ds4^2=-F_delta dU^2+dr^2/F_delta+r^2 dOmega^2.           (3)
```

The extremal metric has a double root at r=B and
F_ext=(r-B)^2/A^2+O((r-B)^3). It has outer coefficient
lambda_e r^2, constant term 1 and a finite vacuum-subtracted geometric
mass. Its 1/r^2 coefficient is a comparison coefficient, not a changed
physical magnetic flux.

To prove positivity and determine every thermal root, set P=t F_ext.
Direct expansion gives

```
P'=1+3xt^2-(1+3x)/t^2+d(-t^-2+4t^-3-3t^-4)
  =(t-1)[(t+1)(1+3x(t^2+1))/t^2-d(t-3)/t^4].           (4)
```

For t>=1, |t-3|/[t^2(t+1)]<=1. The perturbing bracket
relative to the positive first bracket is therefore at most
|d|/(1+6x)<1. Consequently P'>0 for t>1. For each v>1,
delta=P(v)>0 gives exactly one exterior root, at r=Bv;
F_delta(t)>0 for t>v. There are no additional outer roots.

Normalize the boundary clock by alpha=0. The surface gravity is

```
kappa B=[F_ext'(v)+F_ext(v)/v]/2=P'(v)/(2v).
v=1+h,  gamma=kappa A^2/B,
gamma=h K(h),
K(h)=[(1+3x)(v^-2+v^-3)+3x(1+v^-1)
                        -d(v^-4-3v^-5)]/(2C_star).    (5)
```

Here a prime on F_ext or P denotes t differentiation. K(0)=1.
The temperature is also globally one-to-one in v>1: differentiating
P'(v)/(2v) gives
[3x-v^-2+3(1+3x)v^-4+d(3v^-4-16v^-5+15v^-6)]/2,
which is at least [3x-1-34|d|]/2>0, since x>1 and (2).
It tends to zero at v=1 and to infinity with v.
On 1<=v<=1.001,
|K'|<=[5+18x+19|d|]/(2C_star)<3, using (2).
Thus (hK)' >=1-6h>0 on this interval. Every sufficiently small
positive gamma has a unique such thermal root. For gamma<1e-30,
evaluate hK at h=gamma/2 and h=2gamma to bracket the root;
then |h-gamma|<=3h^2<=12gamma^2. In particular

```
R-B=kappa A^2[1+O(gamma)],    |h/gamma-1|<=12gamma.      (6)
```

This root is measured from B, without expanding the center translation
in the much smaller preparation gap. The original k=2000,
eta=1/2000, z_i in [1/4,7/25] Killing temperature is retained.
Its accepted positive preparation interval gives gamma<1e-30;
no zero-temperature interchange is used in (5)-(6).

## An exact action-level split, with both sources

In the spherical reduction use X=r^2 and

```
S_classical=(1/(4G)) integral sqrt(-g2)
                   [X R2+(grad X)^2/(2X)+V_physical(X)],
V_physical=2+6lambda0 X-2Q/X.
V0(X)=2[F_ext(r)+r partial_r F_ext(r)],
q_e=B^2(1+3x),
V0=2+6lambda_e X-2q_e/X
           +2d[-B^2/X+4B^3/X^(3/2)-3B^4/X^2].          (7)
```

For a static alpha=0 solution, the reduced equations require
partial_r(rF)=V0(r^2)/2. Equation (7) solves this identically;
the integration constant is precisely the -delta/t term in (3).
The remaining angular equation follows by differentiating this
constraint, or directly from the X equation. Thus all thermal members
are exact solutions of this reduced potential theory, with
V0(B^2)=0 and V0'(B^2)=2/A^2.

Define fixed bookkeeping functionals on spherical fields,

```
DeltaV=V0-V_physical,
Gamma_ref=(1/(4G)) integral sqrt(-g2) DeltaV(X),
S0=S_classical+Gamma_ref,
Gamma_rest=Gamma_Dirac-Gamma_ref.                       (8)
```

The full common Dirac functional, its domain, state and subtraction
remain in Gamma_Dirac. Equation (8) is an exact addition/subtraction,
not a new finite coefficient, changed renormalization scale or flux.
Hold B, A and lambda_e fixed when functionally varying Gamma_ref.
It is covariant under two-dimensional diffeomorphisms of the spherical
orbit, with X a scalar. No arbitrary nonspherical four-dimensional
local extension of Gamma_ref is asserted. Regulate the split and
combine it in the same prescription; its separate volume terms need
not define separately finite on-shell actions.

With T2,ab=-2 delta Gamma_ref/(sqrt(-g2) delta g2^ab),
use the repository normalization
J_X=4G delta Gamma_ref/(sqrt(-g2) delta X). Its sources are

```
T2,ab=DeltaV g2,ab/(4G),   J_X=DeltaV',
nabla^a T2,ab=J_X partial_b X/(4G),
pA,ref=DeltaV/(16piG r^2),
pS,ref=DeltaV'/(16piG),   rho2,ref+p2,ref=0.             (9)
```

These signs also follow by varying the four-dimensional sphere area.
The joint Ward identity is exact, so an orbit source cannot be retained
while discarding its radius partner. At X=B^2, (7) makes (9) exactly
the two homogeneous quantum pressures selected by the product
Einstein equations. At infinity both pressures tend to
3(lambda_e-lambda0)/(8piG)=-11lambda_e^2/(960pi^2),
the chosen outer vacuum. Equality at the product is equality on the
homogeneous geometry; it does not set Gamma_rest's stress at a finite
thermal bolt to zero.

For example, with rho2,vac=A_v r^2 and
A_v=3(lambda0-lambda_e)/(2G), the seed's physical-Q mass definition
F=1+lambda0 r^2-2m/r+Q/r^2 has

```
mu=m-G A_v r^3/3,
mu_infinity=(1+2x)B+delta B/2.                          (10)
```

The d terms and the difference between Q and q_e contribute inverse
powers of r to mu. Equation (10) proves finiteness for the comparison
metric only. It is not the corrected physical mass or full-action
Hamiltonian.

## What this supplies for spatial order reduction

The covariant bookkeeping family S_lambda=S0+lambda Gamma_rest has
the seed at lambda=0 and the original retained theory at lambda=1.
The homogeneous product remains stationary for every lambda, because
both components of Gamma_rest's source vanish there. This resums the
known center translation rather than Taylor expanding it across the
old classical bolt. The asymptotic vacuum also matches.

If Y denotes reduced geometric fields in a fixed gauge, with a prescribed
covariant HH state/domain map and boundary ensemble, write the equations as
E0(Y)+lambda J_rest(Y)=0. A perturbative spatial reduction would solve

```
E0(Y0)=0,
D E0(Y0)Y1=-J_rest(Y0).                                (11)
```

Higher spatial derivatives in J_rest are evaluated on Y0 in this
equation, including the retained Box R term. They are not removed from
the theory or applied to the unknown Y1 at this order. The Ward identity
comes from the action split, rather than independent filtering of
stress components. If the differentiability and domains required for
Taylor's formula have been justified, the exact residual at
Y_app=Y0+lambda Y1 is

```
E_lambda(Y_app)=lambda^2 {
 integral_0^1 (1-t)D^2 E0(Y0+t lambda Y1)[Y1,Y1]dt
 +integral_0^1 D J_rest(Y0+t lambda Y1)Y1 dt }.           (12)
```

This is a conditional residual identity, not an established bound.
The necessary norms must control the additional spatial derivatives,
the HH family and the moving domain. Exact centering alone does not
make (12) small relative to the thermal gap. A useful quantitative
criterion is that its dimensionless horizon forcing, in the units
eA,eS of [horizon matching](stationary-horizon-product-matching.md),
is bounded by C_rem epsilon^2 L_h fT, where epsilon is a bound on Y1
in a specified strong metric norm and L_h=1+log(1/h),
and that the selected branch has a controlled inverse in that same
horizon-weighted norm. The resulting relative-gap error can then be
bounded by the inverse constant times C_rem epsilon^2 L_h.
Neither this mixed h-weighted response estimate nor that inverse is
proved here. A source vanishing at h=0 can still scale as sqrt(h),
which does not meet this criterion.

In particular the earlier fixed-RN source, mass and first-order lapse
bounds are not estimates of J_rest on (3). The companion
[horizon-source comparison](stationary-horizon-source-defect.md)
now bounds both complete Dirac/product pressure defects and the
product-minus-reference contribution at this explicit seed's bolt.
This supplies its horizon forcing value; it does not bound the full
exterior forcing or its response to an unknown corrected metric.
The corrected state, gap-weighted response/inverse, higher-order
remainder and full boundary charge remain open. Analytic radius loss or a finite
order expansion does not by itself remove the
[same-regularity solver obstruction](stationary-exterior-solver-obstruction.md).
No exact global branch or physical opening follows from this seed.

Remaining to-do list: extend the seed bolt value to the complete common
exterior source and horizon-weighted metric response; justify the spatial
approximation and its remainder or a full higher-order inverse; match
the corrected exterior/state and full boundary charge.
