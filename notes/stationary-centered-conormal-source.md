# Conormal source control for the centered exterior

2026-10-16, S139/S074. Ordinary supporting proof, independently reviewed
before focused verification, for the complete-response milestone.
No focused checker has been run for this note.

This note supplies the three outer source derivatives needed to test the
regular completion in [the geometry construction](stationary-centered-first-order-geometry.md).
It retains the common four-dimensional subtraction and the actual paired
domain. It does not identify three spatial derivatives with a metric
response, or certify a controlled full approximation by itself.

## Statement and normalization

Use the centered seed, parameters and fixed HH ensemble of
[the complete source proof](stationary-centered-exterior-source.md).
Write `ell=A=B/sqrt(C)`, `t=r/B`, `u=1/t`, and

```
Y=sqrt(C) integral_0^u calB(v)^(-1/2)dv,
D=Y partial_Y, L(Y)=1+log(1/Y),
Theta=B^2 tau_hat,mass,
c=-a sqrt(x)/C^(3/2), S(Y)=c C^2/(5pi Y).
```

Here Y is dimensionless proper distance to the conformal-frame outer
boundary, measured in ell units. D is the Euler derivative, not the
ordinary derivative. The result is

```
|D^j[Theta-S]| <= 1e20 j! 100^j q^2 L(Y),
0<Y<=1/10, 0<=j<=3.                              (1)
```

All physical quantities and `mu0=1/b_classical` stay fixed. In particular
the auxiliary dilation used below is not a change of the physical HH
state. The singular coefficient is differentiated together with its
profile: its auxiliary dilation is `c_lambda=exp(3lambda)c`.
It is not an independently adjustable boundary datum.

For ordinary derivatives, `Y^j partial_Y^j=D(D-1)...(D-j+1)`.
The coefficient sums at j<=3 show that (1) implies

```
|partial_Y^j[Theta-S]|
 <=2e20 j!100^j q^2 Y^(-j)L(Y).                    (2)
```

The universal `S` is removed before applying (1)-(2). It is restored
with the correlated Weyl contribution in the physical density and
pressure; removing it is not a new finite subtraction.

## Exact heat formula and common finite parts

The actual Euclidean paired projector induces complementary D/N squares
at the geodesic outer cut, with no angular-mass Robin term. The rank-two
orbit trace is the diagonal of the evenly reflected scalar problem;
there is no additional half factor. In proper coordinates put `f=1/a`.
For a variance-two bridge `X_s=Y+sqrt(2tau) beta_s`, set

```
F=<f(X)^2>, G=<f'(X)>,
J=f(Y) F^(-1/2) exp[tau G^2/(4F)].                 (3)
```

The zero-winding heat diagonal is `E[J]/(2pi tau)`, with the understood
fixed thermal-circle normalization. The outside `f(Y)` is part of the
covariant density. Every bound below includes its derivative.
The nonzero winding terms are treated separately, not discarded.

The common angular weight is

```
M_n^2=n(n+q)/C, d_n=2(q+2n),
Atilde(tau)=sum d_n M_n^2 exp(-M_n^2 tau),
Atilde<=24C/tau^2,
sum d_n/M_n^4<4C^2/q.                              (4)
```

The fixed-scheme value formula is the one in
[the common non-bolt calculation](stationary-reference-nonbolt-stress.md):
subtract `k0/tau+k1+tau k2` inside the angular proper-time integral,
and add back the common a0 and k2 finite parts. They are proportional to

```
zeta_h'(-1)+(1+2log(mu0 B))zeta_h(-1),
Fang+4log(mu0 B),  Fang=4gamma_E-2H_q.              (5)
```

The pole's `+2` cancels the corresponding local anomaly. The k1 terms
cancel. The a0 contribution is constant in position and is bounded by
167 q^2 before the negligible fixed-scale shift; use 200 q^2 here.
No separate angular renormalizations, scale reset or suppressed finite
constant is used in deriving (1).

## Auxiliary caps and finite parameter family

Move the auxiliary inner seam to proper distance 24 beyond the nearest
observation, instead of changing the physical geometry. The exact seed
bounds give caps with

```
R_even=max f/min f <=1e14,
R_smooth<=4e14,
minimum thermal circumference>=1e17.               (6)
```

For example take `zc=exp(-24)/1001` in the compact inverse radius.
The actual benchmark `h<=eta/k^9`, not merely its coarser `h<1e-30`,
gives the last inequality: the cap construction has circumference
at least `zc/h`, and `e<3` gives
`zc/h>2000^10/(1001*3^24)>1e17`.
The smooth negative continuation agrees
with the analytic centered-seed branch through Y=-1 before tapering. Only
the even cap is used at long proper time.

A C6 taper is explicit: integrate the normalized beta(7,7) polynomial
and use its decreasing complement. On width 1/2, derivatives through
order six are at most `52^j`. This follows directly from the Bernstein
derivative bound `2^(j-1)13!/(13-j)!`, followed by the width rescaling.
On width 1/10 use `260^j`. No spatial analyticity of a taper is assumed.

The needed parameter disk is `|lambda|<=1/100`. Near the actual
analytic geometry use `f(exp(lambda)Y)`. On the negative far cap use
instead the finite, globally parameter-holomorphic jet

```
P_lambda f=sum_(j=0)^3 lambda^j D^j f/j!.           (7)
```

Blend these two on `[-1/2,-2/5]`. Their lambda jets through order
three agree everywhere at zero, including all spatial derivatives
needed here. Thus this blend has the true real dilation jets without
ever evaluating a C6 taper at a complex spatial argument.
On the positive side the analytic family is used to Y=4 and is
blended to the fixed cap by Y=4.5. Its resulting distant defect is
explicitly retained below.

Here are finite guards for (7). For `w=f'/f` the capped negative
derivatives through order five are bounded respectively by

```
(2,60,4000,4e5,7e7,1.2e10).
```

The actual negative branch has the smaller real first-three bounds
`(1,3,9)`. The centered compact polynomial has degree six: its RN
quartic part is augmented by `-2e u^5+e u^6`, `|e|<1e-24`.
On the Picard disk `|u|<1.25`, `a/C<.334` and
`(Q0+d)/C<.501` bound the quartic curvature by 7.3; the two
additional curvature terms are less than 1e-20. Thus the same
complex .05 strip has `|K|<8`. For `w=f'/f`, the identity
`K=w'-w^2`, equivalently `w'=K+w^2`, gives
the remaining derivative bounds
`4 j!20^j`. The preceding cap array follows by the Leibniz rule
with the explicit taper. Exponential Bell polynomials then give

```
|f^(j)|/f <=(1,2,64,4368,5e5,8e7,1.4e10), j=0..6.
```

Direct differentiation of `D^j f` gives on `|Y|<=1`

```
|P_lambda/f-1|<.024061,
|P_lambda'|/f<2.976,
|P_lambda''|/f<149.2,
|P_lambda'''|/f<16174.                             (8)
```

For clarity, `partial_Y^a D^3 f` equals
`Y^3 f^(a+3)+(3a+3)Y^2 f^(a+2)`
`+(3a^2+3a+1)Y f^(a+1)+a^3 f^(a)`.
The analytic-minus-polynomial tails through spatial order three
are bounded by `(.0024,.012,.444,120)` on the parameter disk.
Their blend's third derivative costs at most
`42182.4+2433.6+346.32+120<45200` in the same relative units.
Consequently the entire family has the strict guards

```
|f_lambda'/f_lambda|<5,
|f_lambda''/f_lambda|<1e3,
|f_lambda'''/f_lambda|<1e5,
|arg f_lambda|<.05.                               (9)
```

The positive analytic region has the same guards. These estimates
are coefficient calculations for a finite parameter family, not an
analytic continuation claim about the artificial spatial cap.

## Complex smooth and reflected heat bounds

Let `S=sqrt(2)sup|beta|`, so `P(S>s)<=2exp(-s^2)`.
The actual analytic geometry has a uniform complex spatial strip .05.
For the anchor and the interval from it to Y<=.1, `Y*<2`; safe paths
have real span at most 2.4. Their parameter shift plus the bridge
Cauchy shift is bounded by

```
2.4(exp(.01)-1)+exp(.01)/40 <.05.                  (10)
```

On the bridge disk `rho=1/[40(1+S)]`, the finite analytic guards give
`Re F>=.6775 Q`, `|G|^2/|F|<=25/.6775`,
`|F_Y/F|<17`, and exponent derivative less than 14.
Including the outside density therefore gives

```
|J_Y|<=100 exp(.4S).
```

The sixth Cauchy remainder, after the expected odd terms vanish,
is at most `100*40^6*3000 tau^3=1.2288e15 tau^3` before the
heat prefactor. Unsafe paths have `S>.4/sqrt(tau)`. The global
guard `1e5 R(1+tau)exp(7tau)` and the exact conversion
`sup_(tau<=.01) tau^-3 exp(-.16/tau)<.113` bound their contribution
well below `1e21 tau^2`. Thus, for the subtracted smooth heat value
and gradient, use the common enlarged constant

```
Cs=1e21.                                          (11)
```

For the reflected part retain the explicit cubic profile K3 before
bounding its remainder. On a crossing path of span V<=.4 the local
complex displacement is less than .0041, hence its log ratio is
less than .01. The interpolated F,G satisfy `|F|>.1`, `H<3.5` for
`H=F^(-1/2)exp(tau G^2/(4F))`. Safe derivative guards are

```
(|H|,|H_F|,|H_G|,|H_FF|,|H_FG|,|H_GG|)
 <=(6,30,30tau,500,1000tau,50tau),
|H_F+1/2|<=2e5(V^2+tau).
```

The quartic remainder arrays are

```
|DeltaF-DeltaF3|<=5e7 V^4,
|(DeltaF-DeltaF3)_Y|<=4e8 V^3,
|DeltaF|<=3e7 V^3, |DeltaF_Y|<=2e8 V^2,
|DeltaG|<=1e7 V^2, |DeltaG_Y|<=3e7 V,
|F_Y|<=800V, |G_Y|<=100.
```

Differentiate the exact interpolated identity
`R=integral(H_F DeltaF+H_G DeltaG)dtheta+DeltaF3/2`.
Its seven product terms, including the `H_F+1/2` cancellation,
sum to `|R_Y|<=3e13(V^3+tau V)`. Multiplication by the outside
density and the correlated cubic density term gives
`2e14(V^3+tau V)`. On crossing paths `V<=2sqrt(tau)S`.
Gaussian moments and the crossing probability then give the
localized `sqrt(tau)exp[-Y^2/(4tau)]` profile after division by
the heat prefactor. No event indicator is differentiated.

For Y<=.1 the exited-path conversion is bounded by
`sup tau^(-3/2)exp(-.085/tau)<.21`. For Y>=.1 use a tighter
global gradient than the deliberately loose smooth-tail bound:
from (9), the covariant logarithmic derivative is at most
`10.1+2590tau`; hence `|J_Y|<50R` for tau<=.01. The two caps have
combined ratio at most `5.5e14`. Their crossing heat difference
is therefore less than

```
1e17 tau^(-1)exp(-Y^2/tau).
```

Its ratio to `sqrt(tau)exp[-Y^2/(4tau)]` is at most 631, attained
at the limiting Y=.1, tau=.005. Thus the enlarged reflected
constant, including both regions and density factors, is

```
Cb=1e21.                                          (12)
```

## Actual domain and the two different distant defects

Only cap rows are estimated; no hypothetical inner load or actual
off-diagonal Green kernel is substituted. The two-sided cutoff
identity for positive squares uses the actual inverse norm `M^-2`
and cap rows. In first-order paired coordinates multiplication
weights preserve the domain, and Combes--Thomas gives
`||(A+M)^(-1)||_weighted<=2/M` at weight M/2.

A local radius `1/(20M)`, nested cutoffs and the flat Fourier
evaluation identity yield a point row less than `3e4 exp(-MD/2)`.
For reference the variable-principal absorption uses .05 on the
support of the evaluation cutoff; its larger energy neighborhood
only needs .07. Tangential derivatives commute with the domain.
For k<=5, cutoffs of width `1/[20M(k+2)]` have gradient at most
280M; Lichnerowicz/Caccioppoli bounds each next derivative by
600M. The conservative bootstrap `(2000M)^k k!` suffices.
Normal derivatives through three follow from the first-order
radial equation using only f''; one never differentiates the
nonsmooth reflected f''' across the boundary.

These elementary estimates give cap-row derivatives through total
order three bounded by

```
1e12 j!1e4^j M^j C0 exp(-MD/5), C0<=1e23.          (13)
```

For the actual/cap seam `D>=24`, quadratic row products therefore
have `exp(-9.6M)`. The exact benchmark angular sum obeys

```
sum d_n M_n^5 exp(-9.6M_n)<1.2e-67 q.              (14)
```

One direct proof uses `M_n>=18sqrt(n)`, `d_n<=6qn` and a decreasing
integral comparison. The exponential guard `exp(-172.8)<1e-74`
follows from `exp(7/3)>(65/24)^2(25/18)>10`.
After division by `j!100^j`, all actual/cap terms fit the very
coarse reserve `4e11 q^2`.

The positive compact-dilation defect starts at Y=4, not at this
seam. It must use its own separation `D>=3.9` for Y<=.1.
Real finite resolvent differentiation has vertex coefficients
`q_j<=1e4^j`, graph bound ten, and hence inverse jets bounded by
`2j!1e5^j/M`, j<=3. Every difference term has a far-supported
vertex, so two point rows give the permissive per-mode bound
`1e12 j!1e5^j M^2 exp(-3.9M)`. Independently,

```
sum d_n M_n^2 exp(-3.9M_n)<2.1e-27 q.              (15)
```

The decreasing integral is bounded by
`1944q exp(-70.2)[1+2/(70.2-5)]`.
Thus the normalized defect contributes less than `8e-9 q^2`.
It is not assigned the distance 24. The negative finite-jet
construction (7) has no corresponding defect through order three.

## Thermal winding is retained

Under the exact dilation heat identity the circle is beta exp(-lambda),
not beta. The auxiliary winding bound is uniform on the parameter
disk. With `d>=1e17`, `Re F>=.99 Q`, `|G|^2<=25Q`, completing the
Fourier Gaussian gives a bound with exponential

```
exp[100tau-d^2 w^2/(20tau)], w!=0.
```

Since M^2>=324 the growing term uses less than half the mass
exponential. Extract `exp(-Md|w|/10)` by allocating
`M^2tau/4+d^2w^2/(100tau)>=Md|w|/10`.
The residual integral is bounded by one: with
`z=M^2tau/4`, `a=Md|w|/10>=1`,

```
integral_0^infinity z^-1 exp(-z-a^2/z)dz
 <=integral_0^1 z^-1 exp(-1/z)dz
   +integral_1^infinity z^-1 exp(-z)dz <2/e<1.
```

Winding signs and the inverse-tenth-power sum cost less than four.
Using `exp(-x)<=10!/x^10` and `sum d_n M_n^-8<q/1e8`, the full
normalized prefactor is less than `2e43 d^-10 q`, hence less
than `1e-127 q`. Cauchy differentiation on the same .01 parameter
disk preserves the bound after normalization by `j!100^j`.
This is not a zero-temperature replacement of the HH state.

## Anchor, common Mellin identity and finite aggregate

The anchor r=2B has `1<Y*<sqrt(8)log2<2`. Integrate the gradient
from this anchor; `log(Y*/Y)<=L(Y)` for Y<=.1. With (11)-(12),
the following normalized contributions bound the auxiliary
subtracted value, uniformly on the parameter disk:

| Contribution | Bound in units q^2 L(Y) |
|---|---:|
| Long even-cap anchor plus interval length below two | `3*1e7*1e14*(C/q)<=1.056e19` |
| Short reflected remainder | `96*1e21*(C/q)^2<1.190e18` |
| Short smooth anchor plus interval | `3*8*1e21*(C/q)^3<1.047e15` |
| Common finite terms and long polynomial restoration | `<1e8` |
| Angular cubic discrepancy and its removed long tail | `<1e8` |
| Actual/cap and positive fixed-cap defects | `<4e11` |
| Thermal windings | `<1e-100` |

For the local-restoration row no small unproved k2 guard is needed.
The actual analytic strip and `|w|<=5` give
`|w^(j)|<=5j!40^j` on the nested strip. Thus
`|R|<=450`, `|R'|<=36000`, `|R''|<=4.32e6`,
`|R'''|<=6.912e8`, `|Delta R|<=4.5e6`,
`|(Delta R)'|<7.3e8`. Consequently both
`k2=(R^2+2Delta R)/60` and its derivative are below 1e8.
Their finite restoration has normalized coefficient below 1.1e4;
long polynomial subtractions still fit the displayed 1e8 reserve.
The common a0 anchor is included explicitly.

The table is below `2e19`; call this Kaux. It remains necessary to
relate the parameter derivatives to the true Euler derivative, rather
than interpreting parameter holomorphy as spatial analyticity.
The exact heat identity is

```
k_f(exp(2lambda)tau,exp(lambda)Y;beta)
 =exp(-2lambda)k_(f_lambda)(tau,Y;beta exp(-lambda)).
```

For `Ksub=k-(k0/tau+k1+tau k2)/(4pi)`, integration by parts gives

```
D integral Atilde Ksub
 =integral Atilde delta_scale Ksub
   +2 integral tau Atilde' Ksub.                   (16)
```

Both proper-time endpoints vanish. The local k2 restoration obeys
`delta_scale k2=Dk2+4k2`, so its extra `-4k2` term is retained;
the a0 position derivative is zero. Formula (5) accounts for all
finite Mellin pieces. No repeated a0 term is generated.

For completeness define angular moments
`A_l=sum d_n M_n^2(tau M_n^2)^l exp(-tau M_n^2)`.
Their continuum/error bounds follow directly from unimodality:
for `E=x(x+q)/C`, the summand in the E variable has logarithmic
derivative `(l+1)/E+2C/(q^2+4CE)-tau`, strictly decreasing. Hence

```
|A_l-2C(l+1)!/tau^2|
 <=4q(l+1)!/tau+8sqrt(C)(l+2)!/tau^(3/2).           (17)
```

Apply (16) at most three times. The Stirling decomposition of
`(tau partial_tau)^j` and the moment factorials cost, after
normalization by `j!100^j`, at most

```
1+.12+.006+.000264 <2.                            (18)
```

The finite k2 scaling terms have the same smaller normalized
budget and are already covered by the local reserve above.
Equations (16)-(18), Kaux, and the separate distant defects give
`2Kaux+4e11+1e8+1e-100<1e20`. This proves (1).

## Consequence for the response milestone

The canonical stresses are reconstructed before adding the explicit
Weyl, LLL and reference pieces:

```
I(t)=integral_v^t calB_t Theta dt,
B^2 rho_mass=I/(2F)-Theta/t^2,
B^2 p_mass=I/(2F).
```

Their universal outer linear term is combined with the Weyl term
before any physical absolute norm is taken. The local full trace
remains the exact anomaly minus the reference trace. Thus three
outer conormal derivatives are available for the invariant
second-order reconstruction; they are not obtained by differentiating
an undifferentiated O(log Y) statement.

This still requires a quantitative transfer to the actual completed
metric's proper-gauge weighted C5 norm, the complete real-C5 source
response including the cold bolt, and combination with
[the E0 residual](stationary-centered-constraint-residual.md).
No inverse of the full higher-derivative equations, corrected physical
state, or boundary Hamiltonian is proved here.

Remaining to-do list: derive the actual metric-response norm and combine the complete
source and E0 errors; justify the retained-action approximation or
higher-order inverse; match the corrected exterior/state and charge.
