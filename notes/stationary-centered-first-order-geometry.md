# Linear reconstruction and a regular centered completion

2026-10-15, S139/S074. Ordinary proof before focused verification.
The complete source on the explicit quantum-centered seed determines
a unique prepared-temperature linear correction, with the boundary
clock fixed. In coordinates following the bolt, its relative F
correction is below 7.5e-10 and its relative thermal-gap correction
is below 3.4e-10. There is an explicit smooth, positive finite metric
completion of this tangent that preserves the bolt period, outer
clock and quantum center. This completion is a family of comparison
metrics, not a solution of the full retained one-loop equations.
Its complete metric/state response and retained-action residual
are not bounded by this construction.

## Sources, background and the linear boundary problem

Use the exact potential theory S0 and its seed in
[the comparison proof](stationary-centered-covariant-seed.md).
Write B=b_star, A=a_star, C=B^2/A^2=1+6x+d,
x=lambda_e B^2, R=B(1+h), and alpha0=0. All physical inputs,
the scale mu0=1/b0, the selected HH ensemble and paired domain
are those of [the source proof](stationary-centered-exterior-source.md).
Here rho,p are its sphere-integrated Gamma_rest energy and radial
pressure; pS is its four-dimensional angular pressure. The source
is evaluated on the seed, not on the unknown correction.

The inputs are

```
sigma_hor=|G rho(R)/(kappa R)|<5.1e-12,
L=G integral_R^infinity |(rho+p)/(rF0)|dr<3.3e-10,
Qrho=(GB^2/x)sup |rho'|/r<1.3e-10,
Qp=(GB^2/x)sup |p'|/r<1.3e-10,
B integral_R^infinity |rho|dr<1e16 q^2.                 (1)
```

HH regularity gives rho(R)+p(R)=0. The first line follows from
the orbit-pressure result in
[horizon matching](stationary-horizon-source-defect.md); its angular
partner is also retained below. The source's joint Ward identity
is the compatibility condition for the angular linear equation.
The reference potential is fixed under variation, including B,A,
lambda_e. Its mixed pressures depend only on the areal radius.

Let f=delta F, a=delta alpha in Eulerian areal coordinates, and
s_R=delta R. The complete linear constraints about S0 are

```
(r f)'=-2G rho,
a'=G(rho+p)/(rF0),
a(infinity)=0,
f(R)+2kappa s_R=0,
f'(R)+F0''(R)s_R+2kappa a(R)=0.                        (2)
```

The last line fixes the physical surface gravity, not the period
of a coordinate time subsequently rescaled. Their unique solution is

```
a(r)=-G integral_r^infinity (rho+p)/(tF0(t))dt,
D=F0''(R)/2+kappa/R,
s_R=[G rho(R)/R-kappa a(R)]/D,
f(r)=-[2kappa R s_R+2G integral_R^r rho(t)dt]/r.         (3)
```

The remaining angular equation follows from the linearized joint
Ward identity and (2). It is not an independent freely assigned
pressure. The finite mass datum is consequently fixed at this order:

```
delta mu_E(R)=kappa R s_R,
delta mu_h,moving=[mu_seed'(R)+kappa R]s_R,
delta mu_infinity=kappa R s_R+G integral_R^infinity rho dr. (4)
```

Here mu_E is the Eulerian variation and mu_seed(r) is the seed's
vacuum-subtracted geometric mass. Its derivative term cannot be omitted when
quoting a mass variation at the moving horizon. The reference
already contains its asymptotic vacuum subtraction. Equation (4)
concerns geometric mass, not the complete Hamiltonian.

## Uniform geometry bounds and a moving bolt

For the prepared seed, h<1e-30, x>.99k^3>7.9e9 and |d|<1e-12.
The x bound follows from xi0>k^3 and the accepted product-radius
and outer-vacuum bounds; the weaker x>1000 also holds.
The explicit rational formula, with t=r/B, gives

```
1.9x/B^2 <= F0'' <=12.1x/B^2,
5 < kappa B/(xh) <7,
|2F0''+rF0'''|<10C/B^2.                               (5)
```

For the sharper convexity bounds, isolate the x-dependent part of
B^2F0'' as 2x[1-4/t^3+9/t^4]. The bracket has its minimum
26/27 at t=3 and maximum 6 at t=1. The remaining constant,
thermal and d terms have absolute coefficient sum below 54,
much smaller than .01x at this benchmark. The last inequality
follows directly from
2F0''+rF0'''=4lambda_e+4m_seed/r^3-12Q_c B^2/r^4
+72dB^3/r^5-80dB^4/r^6.

Put y=r-R. In particular

```
F0(r)>=2kappa y+.95x y^2/B^2.
D A^2=1+O(h), |D A^2-1|<7h,
|s_R|/(R-B)<3.4e-10.                                  (6)
```

The second line uses D=(w+j/2)/R^2, w=R^2F0''/2,
j=2kappa R, and the explicit seed bounds
|w-C|<3Ch, j<3Ch. Thus no source estimate is divided by h
in the prepared-temperature gap bound.

Pull back by r_physical=r+s_R at first order and vary the areal
scalar X at the same time. Define fbar=f+s_R F0'. It vanishes
at the fixed coordinate bolt r=R and satisfies
fbar'(R)=-2kappa a(R). Let T=rF0'. Exact subtraction in (3) gives

```
fbar/F0={s_R[T(r)-T(R)-T'(R)y]
          -2Rkappa a(R)y-2G integral_R^r[rho(t)-rho(R)]dt}
                                                        /(rF0).
```

Here T'(R)=2RD. The first numerator has absolute value at most
5C|s_R|y^2/B^2. For the last integral use
|rho'(t)|<x Qrho t/(GB^2) and t<=r. With (5), even the
weaker .9x lower curvature bound proves

```
||fbar/F0||_infinity <= L+3Qrho+80|s_R|/B <7.5e-10,
||a||_infinity<L.                                     (7)
```

This is a norm on a moving-bolt pullback. An Eulerian expansion
across the displaced root is neither needed nor asserted. The
angular metric variation is 2s_R/r and is also uniformly small.

## An exact common conformal disk

Choose D0=2atanh(exp(-h)), so log tanh(D0/2)=-h and D0<82.
For the full seed, define

```
v(r)=log tanh(chi(r)/2)
    =-h-kappa integral_r^infinity dt/F0(t),
ghat0=(B/r0)^2 g2,0
     =A^2 exp(2sigma0)[dchi^2+sinh(chi)^2 dphi^2],
phi=kappa U_E, 0<=chi<D0,  phi identified modulo2pi.    (8)
```

This is a valid full exterior disk: v<0 everywhere and chi=D0
is the actual outer endpoint. The bolt is chi=0. The auxiliary
H2 cap used for the cold source proof is not substituted for it.

The linear displacement xi=delta r at fixed chi and the conformal
variation u=delta sigma have the exact representations

```
H=a+fbar/F0,       J(r)=integral_r^infinity H(t)/F0(t)dt,
xi=s_R-F0 J,
u=a+fbar/(2F0)-s_R/r-(F0'/2-F0/r)J.                   (9)
```

These follow by differentiating (8) with fixed kappa and fixed
D0. The identity integral_r^infinity F0'/F0^2=1/F0
isolates the moving-bolt term. Fixed temperature gives H(R)=0
exactly, so J is regular at the bolt. In particular xi(R)=s_R.

Here are quantitative real bounds, including the endpoint
cancellation. Globally ||H||<1.08e-9. From (5), elementary
integration of the quadratic lower bound gives

```
[F0/(r-B)] integral_r^infinity dt/F0(t)<10.             (10)
```

Indeed, writing z=(r-R)/(Bh) and c=kappa B/(xh), the bound
is at most log(1+4c/1.9)+12.1/1.9<10. This uses
z log(1+a/z)/(1+z)<=log(1+a) and log(1+a/z)<=a/z.
Consequently, with V=xi/(r0-B),

```
||V||<1.2e-8.                                        (11)
```

For completeness a derivative bound follows without a second
source derivative. On 0<=y<=Bh, (1), (3) and the exact
subtraction above imply
|H|<5e-10 y/(Bh). For the lapse difference use
|a(r)-a(R)|<Qsum y/(10Bh), Qsum=Qrho+Qp<2.6e-10.
For the horizon term use
0<=1-2kappa R y/(rF0)<.7y/(Bh); the density term costs
Qrho y/(10Bh), and the s_R term is negligible. Thus
|J|<[5e-10+2.8||H||]/(2kappa) in this inner interval.
Outside it, (5) gives F0'(r)y/F0<13. From
xi'=H-F0'J and (10), treating y<=Bh and y>=Bh separately,

```
||(r0-B) partial_r V||<3e-7,    ||u||<2e-7.            (12)
```

For the second region, F0' integral dt/F0<260; the first
region uses the regular J bound just stated. Substitution in (9)
gives the u bound. Smooth radiality implies V_chi(0)=u_chi(0)=0.

At the cold seam s=r/B-1=.01 the stronger bounds are

```
|xi|/(r-B)<1.9e-9,  |u|<3e-9.                        (13)
```

For .01<=s<=.1 use F0>=.8Cs^2; beyond .1 use
F0>=.98xs^2. Their reciprocal integrals give
integral_r^infinity dt/F0<175B/C at s=.01, while F0<=Cs^2
and |(F0'/2F0-1/r)(r-B)|<1.01 there. Equations (9) yield (13).

## The invariant inner curvature ideal

Raw a'(R) can be of order 1/h. It must not be bounded as if
it were a covariant curvature derivative. The linear constraints
and angular equation give instead

```
delta R2=-16piG pS+2[-f+G(p-rho)]/r^2,
delta Rbar=-16piG r^2pS-6f+4G(p-rho)
          =-8piG r^2 Arest-6f+6G(p-rho),
Rbar=B^2 Rhat,  Arest=A_Dirac-2(pA_ref+pS_ref).          (14)
```

The full four-dimensional trace at fixed areal r is particularly
simple: delta R4=-8piG Arest. The reference mixed trace is a
fixed function of r. The conformal-disk curvature variation adds
xi Rbar0' to the second line of (14).

There is an elementary local bound |Arest'|<1e5q^2/B^5.
Here are coefficient counts to make it independent of an assumed
stress derivative. The seed F is a Laurent polynomial with
exponents {2,0,-1,-2,-3,-4} and absolute coefficient sum<2C.
For dimensionless B^2R4 the monomial multiplier is
-(n+1)(n+2), with the additional 2/t^2. Its coefficient sum
is below 61C, with derivative sums 366C,2562C,20496C through
order three. Hence the coefficient sum for the derivative of
B^4 Box R4 is below 123708C^2. The Weyl-square derivative is
below 14884C^2: use the numerator
t^2F''-2tF'+2F-2, coefficient sum 61C, and divide its square
by 3t^4. The Euler derivative is below 8832C^2, from
4[(F')^2+(F-1)F'']/t^2. The retained coefficients
c=1/20,a=11/360,b_an=1/30 therefore bound
16pi^2|partial_t(B^4 A_Dirac)| by
(745+295+4124)C^2+4q^2/3<5166q^2.
Since 16pi^2>144, this is below 36q^2. The explicit reference
profiles have trace derivative below 5.6e4q^2 in these units,
which proves the asserted 1e5 bound. No derivative of a heat
remainder or unknown metric response is used.

The horizon trace value and (1) now give on h<=s<=.01

```
|8piG r^2 Arest|/C<5.1e-11 s,
G(|rho|+|p|)/C<5.72e-11 s,
|f|/C<7.6e-12 s,
A^2 |delta Rhat at fixed chi|<1e-7 s.                 (15)
```

For the last inequality (14) contributes less than 4.5e-10s;
the displacement adds at most 1.2e-8s*8.1, since
|B Rbar0'|/C<8.1 in this patch. This retains the Lie term.

Apply the curvature equation in the accepted OLD cold chart,
where the physical background is exp(2sigma_cold)gH2. Its
potential Vpot=-exp(2sigma_cold)A^2Rhat0 obeys
Vpot>=2-20s. To verify the one-sided bound, the explicit seed
curvature gives -Rbar0/C>=2-10s, and sigma_cold>=-2s;
thus Vpot>=(1-4s)(2-10s)>=2-20s.
This chart is used only on the actual inner patch. Relative to
(8) it differs by the exact fixed baseline dilation
tanh(chi_cold/2)=exp(c_h)tanh(chi/2). The seam integral
above bounds |c_h|<277h. Keep this transition fixed for the
whole correction family. Then u is a scalar under it; no
auxiliary capped metric is identified with the full exterior.

Let Dc=arccosh(.01/h) and

```
W=h[(Dc+1)cosh(chi_cold)-chi_cold sinh(chi_cold)].
(-Delta_H2+2)W=h[3cosh(chi_cold)+chi_cold/sinh(chi_cold)].
```

W(Dc)>=.01, W>0 and 20sW/(hcosh(chi_cold))<.3.
The last bound follows from s<=1.008hcosh(chi_cold) and
max s[(Dc+1)-chi_cold tanh(chi_cold)]<.015.
Thus (-Delta_H2+Vpot)W>=2.7hcosh(chi_cold).
Its positive potential also permits the maximum principle.
Equation (15) and |sigma_cold|<.05 put the absolute forcing
in the conformal curvature equation below
.56e-7hcosh(chi_cold). At the seam (13) applies. Comparison
with 1e-6W proves

```
|u|<=1e-6W,   |u(0)|<8.3e-5h.                       (16)
```

The second radial derivative at the bolt has the same h ideal
by the curvature equation and (15). This is a bound on the
geometric conformal variation, with its homogeneous dilation
mode fixed by the actual seam data. Separate bounds on a and
a' would not prove it.

## An exact regular finite completion of the tangent

In the common disk (8), for 0<=lambda<=1 define

```
sigma_lambda=sigma0+lambda u,
r_lambda-B=(r0-B)exp(lambda V),
g4,E,lambda=(r_lambda/B)^2 A^2 exp(2sigma_lambda)
                 [dchi^2+sinh(chi)^2dphi^2]
                      +r_lambda^2dOmega^2.            (17)
```

This is an explicit comparison family with precisely the tangent
(2)-(3), not a nonlinear field equation. Equations (11)-(12)
give r_lambda>B and
partial_chi r_lambda=exp(lambda V)partial_chi r0
[1+lambda(r0-B)partial_rV]>0 away from the bolt.
Its horizon gap is h_lambda=h exp(lambda V(0)), so
|h_lambda/h-1|<3.5e-10. The bolt is smooth and its U_E period
is exactly 2pi/kappa, because phi has period 2pi and both metric
coefficients use the same smooth conformal factor at chi=0.
No further lapse rescaling changes that clock.

The outer data follow exactly from (8), not a shooting assumption.
For any fixed-kappa asymptotically AdS geometry with alpha(infinity)=0,

```
sigma(D0)=log[B sqrt(lambda_e)/(A kappa sinh(D0))],
sigma'(D0)=-coth(D0),
r(chi)~c_pole/(D0-chi), c_pole=kappa sinh(D0)/lambda_e. (18)
```

The second relation is the geodesic product-frame boundary. The
linear source constraints give u(D0)=u'(D0)=0 and V(D0)=0.
More precisely their asymptotic mass term is proportional to
Y^3 and the remaining terms are O(Y^4(1+|logY|)), where
Y=D0-chi. The same first-derivative behavior follows directly
from the constraints and the quantitative trace remainder,
not by differentiating an unsigned asymptotic estimate.
Consequently (17) preserves all three relations (18) exactly.
Exponentiation introduces only order Y^6 terms at second order;
it changes neither the outer vacuum curvature nor constant
asymptotic normalization. Higher conormal metric-response bounds
are not inferred from this value and first-derivative statement.

The cubic terms have an additional necessary correlation. Put
m1=delta mu_infinity from (4). The outer constraints and (9) give

```
V=m1/(2lambda_e r^3)+O(r^-4(1+log r)),
u=-m1/(lambda_e r^3)+O(r^-4(1+log r)),
u3+2V3=0.
```

Indeed f=-2m1/r+O(r^-2log r), while a=O(r^-4log r).
The lapse reconstructed exactly from (17) is

```
alpha_lambda-alpha0=2log(r_lambda/r0)+2lambda u-lambda V
                         -log[1+lambda(r0-B)V_r].
```

Its cubic coefficient is 2lambda u3+4lambda V3=0.
Thus the finite completion preserves the cancellation of the
physical r^-3 lapse term; this is stronger than merely fixing
alpha(infinity). No independent cubic outer parameter is added.

The completion also avoids a constant second-order product error.
At the bolt, (16), its curvature equation and the radial gap
factor in (17) keep the orbit curvature and radius in the
product-vanishing h ideal. For example the exact product-frame
curvature identity is

```
Rhat_lambda=exp(-2lambda u)
                 [(1+2lambda u)Rhat0+lambda delta Rhat]. (19)
```

At chi=0 the part beyond linear order has magnitude below
1e-7h^2/A^2, using (15)-(16). The radius derivatives are
also multiplied by r0-B; (12) controls their second derivative
at the bolt. Thus no orderlambda^2 constant translation of B
or constant change of A is generated by this completion.
This is why an additive F/alpha completion with uncontrolled
raw horizon jets would be insufficient.

## Residual and physical scope

The linear reconstruction, its finite mass datum, uniform gap/
positivity bounds and the regular finite completion are established
here. They specify a candidate first-order geometry on the whole
exterior at the original prepared temperature. Equation (17) is
not certified as a controlled approximation to the full retained
equations merely because its tangent solves (2).

To obtain that conclusion one must bound E0(g_lambda)-lambda
D E0(g0)g1 in suitable invariant, horizon-weighted and outer
conormal norms, and bound the complete source difference
J_rest(g_lambda)-J_rest(g0), including the actual HH state and
domain. The conditional order-reduction identity in the seed note
still applies. The bulk trace does simplify its dangerous linear
jet: delta R4=-8piG Arest, and at the moving bolt
delta(Box R4)=2kappa[-a_h R40'-8piG Arest'+s_R R40''].
It does not bound the independent nonlocal source response or
all finite-amplitude higher-order terms. An unrestricted
same-regularity Einstein contraction remains excluded.

The allowed outer Y^4 log Y remainder does not establish membership
of this correction in a bounded analytic or C6 response ball used
by a stronger source-response argument. That hypothesis remains
unproved here. Conversely, a nonzero logarithmic coefficient has
not been established by the upper bound, so this observation is
not proof of nonmembership in such a ball. Conormal response
estimates require their own complete differentiated argument.

Remaining to-do list: bound the full retained-action residual and
horizon-weighted metric/state response in this regular family;
justify the approximation or a higher-order inverse; match the
corrected exterior/state and complete boundary Hamiltonian.
