# Horizon product matching and the prepared thermal gap

2026-10-15, S139/S074. Ordinary argument before focused verification.
The exact horizon equations identify the two-radius homogeneous Dirac
solution as the possible regular extremal center. They also identify
the additional source response needed to match the very small prepared
temperature. The verified fixed-RN source and finite geometric mass,
even when both are available, do not alone establish that response or
a corrected exterior. No computation is reported here.

## Fixed theory, source normalization and exact horizon equations

Keep the one free four-dimensional Dirac field, fixed magnetic charge,
renormalized couplings and zero extra finite curvature coefficients of
[the reference theorem](uniform-dirac-reference.md). In particular
mu0=1/b0 stays fixed, where b0 is the classical extremal radius; it is
not reset to the corrected radius. The benchmark remains k=2000,
eta=1/2000 and z_i in [1/4,7/25], with the original laboratory drive.
Here kappa denotes the prescribed physical Killing temperature times
2pi in the clock normalized by alpha(infinity)=0.

Write the static spherical metric and geometric mass function as

```
ds4^2=-exp(2alpha)F dU^2+dr^2/F+r^2 dOmega^2,
F=1+r^2/L4^2-2m(r)/r+Q/r^2,
m'=G rho2, rho2=4pi r^2 rho4,
alpha'=G(rho2+p2)/(rF).                                (1)
```

Classical magnetic stress is already in Q. Quantum rho2,p2 include
the full common Dirac allocation, including Weyl sources once. These
are necessary exact Einstein constraints for a solution of the full
higher-order effective equations; writing them does not turn the
quantum stress into a prescribed second-order source.

At a regular nonextremal horizon R, F(R)=0, and regularity gives
rho2,h+p2,h=0. Define the orbit curvature radius at that horizon by

```
a_h^-2=-R2,h/2=F''_h/2+(3/2)F'_h alpha'_h,
pA,h=T^U_U,h=T^r_r,h, pS,h=T^theta_theta,h,
F'_h=2kappa exp(-alpha_h).                              (2)
```

The p's in (2) are four-dimensional quantum pressures, not the
sphere-integrated p2. The two independent horizon equations are

```
F'_h/R-1/R^2-3/L4^2=-Q/R^4+8piG pA,h,
a_h^-2+F'_h/R-3/L4^2= Q/R^4+8piG pS,h.                 (3)
```

For example G^theta_theta,h=F''_h/2+3F'_h alpha'_h/2+F'_h/R.
The mixed orbit equation also gives the exact clock law

```
kappa=exp(alpha_h)
 [R^2+3R^4/L4^2-Q-2G R^2 rho2,h]/(2R^3).                (4)
```

The value of the mass at infinity is not another term in (4).
Once R and the source are selected, it instead fixes the mass
integration constant through

```
m_h=R(1+R^2/L4^2+Q/R^2)/2,
rho2,vac=A_v r^2, rho2,ren=rho2-A_v r^2,
mu(r)=m(r)-G A_v r^3/3,
lambda_eff=L4^-2-2G A_v/3,
mu_h=m_h-G A_v R^3/3
    =R(1+lambda_eff R^2+Q/R^2)/2,
mu_infinity=mu_h+G integral_R^infinity rho2,ren dr.       (5)
```

Here A_v is the specified fixed-reference vacuum coefficient for a
first-order reconstruction, or the appropriate actual asymptotic
vacuum coefficient for an actual solution. Existence of that actual
asymptotic vacuum is a separate matching requirement. The horizon
subtraction in mu_h cannot be omitted. Equations (3)-(4) use the
original retained cosmological constant and the full quantum stress;
switching to lambda_eff requires subtracting the same vacuum there.
The full Hamiltonian
can still differ from geometric mass, as in
[the charge audit](stationary-exterior-charge-audit.md).

## What first-order reconstruction does and does not say

For a fixed classical RN/HH source, let R0 be its horizon and impose
the same temperature at first order. Linearizing (1)-(4) gives

```
delta m_h=kappa0 R0 delta R,
0=delta kappa=kappa0 delta alpha_h
 +(F0''/2+kappa0/R0)delta R-G rho2,h/R0,
delta R=[G rho2,h/R0-kappa0 delta alpha_h]/D0,
D0=F0''/2+kappa0/R0=C(R0)/R0^2-3kappa0/R0,
C(R)=1+6R^2/L4^2.                                     (6)
```

Thus D0 stays near C(b0)/b0^2 as the classical thermal gap tends to
zero. The shift formula itself is finite. It does not imply that
delta R is small relative to R0-b0. The full homogeneous reference
already has b_star<b0 with |b_star-b0|/b0>2e-28, whereas the
prepared classical gap is of order 1e-33. The complete fixed-source
linearization must therefore not be promoted to a uniformly small
metric correction across the old bolt. Nor may its exterior radial
profile simply be evaluated at a corrected root below its domain.

## Conditional extremal product necessity

Suppose an actual extremal solution has a smooth limiting horizon B,
finite alpha_B and

```
F(r)=(r-B)^2/a_star^2+O((r-B)^3),
alpha(r)=alpha_B+O(r-B).                               (7)
```

The near-horizon scaling r=B+lambda x,
U=exp(-alpha_B)t/lambda gives AdS2(a_star) times S2(B).
If the needed curvature derivatives are uniformly controlled, all
covariant derivatives of the limiting curvature vanish. In particular
Box R4 tends to zero; it is not set to zero on a finite-temperature
unknown exterior. Algebraic curvature and magnetic anomaly terms
remain and belong to the common determinant.

The state statement also has hypotheses. For a regular near-extremal
HH family require pointed C6 convergence of the Euclidean inner
metrics on every fixed proper ball, uniform local normal-chart jet
majorants, and
the actual outer projector receding to infinite proper distance.
Then each nonzero angular Dirac sector has its positive proper mass
gap. The two-sided locality identity in
[the cold-bolt proof](stationary-reference-cold-bolt.md), equation (14),
removes the distant actual boundary. On fixed balls, the local
resolvent/normal-chart construction converges. After the common
a0,a1,a2 subtraction its uniform M_n^-4 majorant is angular summable,
by that note's equation (22). Dominated convergence therefore gives
the homogeneous massive orbit trace. At the regular bolt the two
orbit pressures agree, so this also determines their limiting value.
The LLL trace/conservation formulas
in [the exact LLL note](stationary-lll-quantum-exterior.md), equations
(1)-(3), give the regular invariant limit directly. The joint
four-dimensional trace then fixes the angular stress as well.

Consequently, under these geometric and state convergence hypotheses,
the limiting pressures are the complete homogeneous common-scheme
pressures pA,prod(a_star,B), pS,prod(a_star,B). Equations (3) at zero
temperature are exactly the two homogeneous Einstein equations.
The quantitative stationary-product theorem in
[the uniform reference](uniform-dirac-reference.md), equations (9)-(13),
already proves a unique solution (a_star,b_star) in the one-percent
box, with fixed-point derivative rows below .02. Its inverse in those
dimensionless homogeneous variables has norm at most 1/(1-.02)=50/49.
Thus any such regular limiting exterior whose product lies in that
box must have B=b_star. This is an attractor necessity statement and
a zero-mode inverse bound, not construction of the global exterior.

The stress/radius subtraction and exterior clock normalization must
both be retained as in
[the centering note](quantum-reference-clock-centering.md), sections
"Exact centering with both reference sources retained" and
"The exterior normalization is still physical data". In particular
b_star^2/a_star^2 is not replaced by its classical RN formula.

## Exact finite-temperature product defects and a gap-scale target

At a finite-temperature horizon define actual defects

```
delta pA=pA,h-pA,prod(a_h,R),
delta pS=pS,h-pS,prod(a_h,R).
```

These include the nonproduct geometry, the actual state and the full
retained derivative sources. Equations (3) are equivalently

```
EA=1/R^2+3/L4^2-Q/R^4+8piG pA,prod
                 =2kappa exp(-alpha_h)/R-8piG delta pA,
ES=1/a_h^2-3/L4^2-Q/R^4-8piG pS,prod
                 =-2kappa exp(-alpha_h)/R+8piG delta pS. (8)
```

To compare with the already solved product map, put

```
z=R^2/b0^2, u=a_h^2/a0^2, a0^2=b0^2/C0,
fT=2kappa exp(-alpha_h)R^3/(C0 b0^2),
eA=8piG R^4 delta pA/(C0 b0^2), eS similarly.
```

If Tz,Tu_den denote the explicit map and inverse denominator of
uniform-reference equation (10), the exact horizon algebra is

```
z=Tz(z,u)+fT-eA,
u=[Tu_den(z,u)-2fT/z^2+(eA+eS)/z^2]^-1.                 (9)
```

This is a useful auxiliary solve when the defects are supplied; it
is not a definition that makes physical defects small. With eA=eS=0
and any supplied |alpha_h|<=1, the new thermal terms have size O(fT) and
derivative O(fT), with fT<1e-30 in the preparation range. The same
one-percent box remains invariant and its derivative row bound stays
below .021. Banach therefore gives a unique auxiliary algebraic root.
Relative to the zero-temperature product, the change has
|delta(z,u)|_infinity<2.2fT. Its first equation gives
|delta z-fT|<.045fT, hence delta z>0 and, conservatively,

```
.9 < (R-b_star)/(kappa a_star^2 exp(-alpha_h)) <1.1.     (10)
```

The last conversion also uses the accepted product-radius bounds,
which are much smaller than these margins. This auxiliary result
recenters before measuring the gap and therefore avoids delta R/h.

A sufficient physical value target is

```
|8piG delta pA|, |8piG delta pS|
 <=.01 * 2kappa exp(-alpha_h)/R,                        (11)
```

or |eA|,|eS|<=.01fT. In addition, along a chosen regular exterior
family the two-radius response rows of the complete added forcing
in (9), relative to the homogeneous map, must have norm below .1.
This includes fT(z,u,alpha_h(z,u)) and eA,eS; separately require
|alpha_h|<=1. The auxiliary
Banach statement above instead treats alpha_h as a supplied fixed
scalar; its value bound alone gives no derivative bound for an actual
family. With the product-map rows below .021, the total map then
has rows below .121. The stated full response target preserves a
contraction, a positive root and a
stable gap bracket. A global higher-order
inverse or justified spatial reduction must provide that family and
its response. The fixed-RN spatial derivative bound on tau is not
this defect-map derivative. Even a bounded finite mass in (5) cannot
supply it, since (5) is an integral constraint whereas (8) contains
local horizon jets and state response.

## A precise missing jet: the retained anomaly sees it at the bolt

This distinction can be proved without guessing a stress constant.
At any fixed positive-temperature RN horizon take alpha unchanged and

```
delta F=A(r-R)^3 chi((r-R)/Delta), chi=1 near zero,
```

with smooth compact support toward the exterior. The horizon radius,
F'_h, F''_h and a_h are unchanged. All algebraic curvature invariants
and the magnetic invariant at the horizon are therefore unchanged.
But direct differentiation of the exact scalar-curvature formula gives

```
delta R4,h=0, delta R4'_h=-6A,
delta(Box R4)_h=-6F'_h A.
```

For the retained Dirac trace coefficient b_an=1/30 this gives exactly

```
8piG delta A4,h=-G F'_h A/(10pi).                       (12)
```

The conformally invariant part of the action has zero bulk trace,
so changing its state cannot cancel this compulsory trace variation.
Meanwhile ||delta F||_C2=O(|A|Delta) and
|delta F/F|=O(|A|Delta^2/F'_h) near its support. At fixed positive
F'_h one may take A unbounded and Delta decreasing fast enough that
both norms tend to zero and F stays positive. The perturbation is
smooth in regular bolt coordinates and preserves the HH period.
Thus unchanged two-radius horizon data and small C2 coefficients do
not bound the product defects uniformly.

The actual trace is 2pA,h+2pS,h. Relative to a zero-defect base, (11)
and (12) require G R|A|<=40pi times the chosen relative error.
If both metrics merely satisfy that defect bound, their difference
instead gives the conservative 80pi factor by the triangle inequality.
Either form excludes an unbounded A. Equivalently a third areal
metric derivative / first radial curvature derivative, or the
corresponding stronger covariant solver norm,
must be controlled. One must not interchange kappa tending to zero
with this unbounded third-jet construction. Removing it by a spatial
order-reduction rule requires the controlled branch/error argument;
its absence cannot be inferred from a small fixed-reference source.

Remaining to-do list: prove the physical horizon product-defect and
metric-response estimates on a justified exterior branch; obtain the
finite mass/charge and global clock normalization; match the corrected
state and retain the original driven preparation, capture, tail and
renormalized spatial-area obligations. No corrected-exterior existence
or physical gate is closed by this horizon algebra.
