# Stationary Weyl sources and the outer matching terms

S139/S074, 2026-10-15. This is a continuation of the frozen
[isolated exterior result](stationary-lll-quantum-exterior.md), not a
replacement of its action or evidence. The retained four-dimensional
Weyl functional has explicit stationary bulk sources. On a classical
RN-AdS reference they contain an outer vacuum term and a logarithmic
geometric-mass tail. The latter cannot be discarded or inferred to be
integrable from asymptotic AdS geometry. The companion massive-sector
calculation identifies the boundary contribution needed to cancel it.

## Conventions and complete compact-variation representative

Use the common one-loop subtraction of
[the anomaly action](spherical-anomaly-action.md) and
[its separate sources](spherical-anomaly-sources.md). This note gives
one exterior's normalization. Identical folded exteriors double the
action and sources; they do not double the source in one exterior's
Einstein equation. Write

```
ds2^2 = -exp(2alpha) F dU^2 + dr^2/F,
X=r^2, sigma=log(r/b0), h2=exp(-2sigma)g2,
T=r(F'+F alpha'),
Rbar=b0^2 R[h2]
    =-r^2 F''-3r^2 alpha' F'-2r^2 F(alpha''+alpha'^2)+2T-2F.
```

Primes here denote r derivatives. The fixed magnetic flux is q;
there is no two-dimensional electric field. In the product frame,
b0^2 Box_h sigma=T-F, b0^2(nabla_h sigma)^2=F and
sqrt(-h2)=b0^2 exp(alpha)/r^2. Integrating sigma Box_h R[h2]
twice by parts gives the following exact representative for compact
variations:

```
4pi W_bulk = integral dU dr L,
L=exp(alpha)/r^2 {
  sigma[Rbar^2/60-Rbar/18+1/15-q^2/3]
  +Rbar T/30-T^2/10+11TF/90-11F^2/180-F/18
}.
```

The surface terms are part of the finite-cut action. In particular,
the additional integration above contributes at each oriented endpoint

```
[sigma exp(alpha) F Rbar' - exp(alpha) F Rbar/r]/(120pi).
```

The preceding anomaly-action flux and the literal boundary Weyl
transgression must also be retained. The boundary functional is
specified in [the boundary note](stationary-exterior-boundary-weyl.md).
The expression L alone is sufficient for compact variations and is
not a boundary-charge prescription.

For z=F,alpha define E_z=sum_j(-d/dr)^j partial L/partial z^(j).
Only j<=2 is required before taking the Euler derivatives. At arbitrary
alpha, the sphere-integrated physical sources obey

```
rho2,W=-exp(-alpha) E_alpha/(4pi),
rho2,W+p2,W=-F exp(-alpha) E_F/(2pi),
alpha'_W=G(rho2,W+p2,W)/(rF)=-G exp(-alpha) E_F/(2pi r).
```

The last equality is the source's contribution to the lapse equation,
not a solution of the complete higher-derivative equations. It follows
directly by varying g_UU and g_rr with respect to alpha and F. This
normalization avoids identifying a four-dimensional stress with its
sphere-integrated counterpart: T4=T2/(4pi r^2).

## Exact RN evaluation, radius source and Ward identity

Evaluate the compact sources after varying, with

```
alpha=0, F=1+lambda r^2-2m/r+Q/r^2,
lambda=L4^-2, Q=r_m^2, m=G M.
```

The resulting expressions are

```
E_F=2[6Q sigma-8Q+3mr]/(15r^4),

E_alpha={
 (432Q^2+288Q lambda r^4-864Qmr+432Qr^2-60q^2r^4)sigma
 -593Q^2-338Q lambda r^4+1408Qmr-504Qr^2
 -33lambda^2r^8+144lambda mr^5-588m^2r^2+288mr^3-11r^4
}/(180r^6).
```

These formulae retain the gauge anomaly; its pure radius potential
contributes to E_alpha and the radius source, but not to E_F. On RN,
R4=-12lambda, Ric4^2=36lambda^2+4Q^2/r^8 and
C4^2=48(m/r^3-Q/r^4)^2. Consequently the retained anomaly is

```
16pi^2 A4=14(m/r^3-Q/r^4)^2/15
          +11Q^2/(45r^8)-11lambda^2/15-q^2/(3r^4),
trace T2,W=(E_alpha-F E_F)/(2pi),
J_X,W=(2G/r^2)[4pi r^2 A4-(E_alpha-F E_F)/(2pi)].
```

The last line is the independent radius source, obtained from the
exact joint trace identity of the Weyl functional. With
rho2,W=-E_alpha/(4pi) and p2,W=E_alpha/(4pi)-F E_F/(2pi),
direct differentiation gives the stationary Ward identity

```
p2,W' + F'(rho2,W+p2,W)/(2F) = r J_X,W/(2G).
```

Thus the radius source cannot be omitted while inserting this stress.
For a regular horizon at r_h>0 the E coefficients are finite and
rho2,W+p2,W has the requisite factor F. The local Weyl contribution
does not have a spurious lapse singularity there. This observation
does not control the separate state-dependent determinant at the bolt.

## Integrated lapse and outer mass term

At first order, on this fixed RN reference, imposing alpha_W(infinity)=0
integrates the lapse source exactly:

```
alpha_W(r)=G/(15pi){m/r^3 + Q[(3/2)log(r/b0)-13/8]/r^4}.
```

For b0=b and the classical extremal parameters
m=b(1+2xi), Q=b^2(1+3xi), xi=b^2/L4^2, this gives
alpha_W(b)=-G(5+23xi)/(120pi b^2). This is an order-reduced
component value. It neither recenters the full horizon nor determines
the prepared HH temperature in the corrected exterior.

The exact outer expansions have distinct implications:

```
E_alpha=-11lambda^2r^2/60+4lambda m/(5r)+O(log(r)/r^2),
rho2,W=11lambda^2r^2/(240pi)-lambda m/(5pi r)+O(log(r)/r^2),
rho2,W+p2,W=-lambda m/(5pi r)+O(log(r)/r^2).
```

The leading term is the AdS Dirac vacuum energy density
11lambda^2/(960pi^2) in four-dimensional normalization. It must be
included in the outer curvature equation or absorbed by an explicit
cosmological renormalization condition. Subtracting this constant
vacuum alone leaves the term -lambda m/(5pi r) in the mass equation
M'=rho2. Its radial integral is logarithmic. The Weyl component by
itself therefore does not satisfy a finite vacuum-subtracted geometric
mass condition on this reference.

The massive product determinant has a leading boundary layer
tau_hat=-lambda m r/(5pi b0^2) in its sphere-integrated product-frame
trace. Conservation reconstructs the opposite physical tails
rho2,hat=+lambda m/(5pi r) and
alpha_hat=-Gm/(15pi r^3), as derived in
[the massive boundary note](stationary-massive-boundary-layer.md).
The two leading coefficients cancel. The companion note supplies a
local bound of order (1+log r)/r^2 on the remaining mass source, proving
convergence at this fixed reference's outer endpoint. Its constants are
not quantified or uniform in the tiny HH temperature. Equality of the
leading coefficients alone would not have supplied that remainder.

## Consequence for the stationary proof

The common-scheme source must combine the Weyl term, the entire
massive tower with its boundary layer, and the LLL allocation before
asserting outer decay or bounding a geometric mass. A uniformly bounded
raw product-frame trace would exclude the required leading boundary
term. The appropriate remainder norm subtracts that term explicitly.
Finite boundary charge further requires the full finite-cut variation
and a stated boundary prescription; geometric mass and action finite
part do not by themselves determine it.

Remaining to-do list: quantify the massive boundary remainder and its
global, horizon-weighted metric response; complete the finite boundary
charge and outer vacuum matching; solve the full corrected stationary
exterior and transfer the pulse result before closing S074-S079.
