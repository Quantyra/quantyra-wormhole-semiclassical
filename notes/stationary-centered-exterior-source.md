# The complete exterior source on the covariant centered seed

2026-10-15, S139/S074. Ordinary derivation before focused verification.
The explicit covariant seed admits the full common Dirac source, not
only its horizon value. Its remainder after the action-level reference
subtraction has an absolutely integrable exterior energy source, a
small complete first-order lapse and a finite density-variation norm.
These are evaluations on the specified seed and its HH state. They
are inputs to a first-order construction, not estimates of the source
on an unknown corrected metric.

Retain the exact seed and action split of
[the centered comparison](stationary-centered-covariant-seed.md),
the complete horizon source of
[the horizon-defect theorem](stationary-horizon-source-defect.md),
and the ordinary heat/Green estimates in
[the cold proof](stationary-reference-cold-bolt.md),
[the non-bolt proof](stationary-reference-nonbolt-stress.md), and
[the finite-mass proof](stationary-reference-finite-mass.md).
No frozen note, finite coefficient, source waveform or preparation
range is changed.

## Statement and conventions

Put B=b_star, A=a_star and use one exterior's normalization:

```
t=r/B, v=R/B=1+h, x=lambda_e B^2,
C=B^2/A^2=1+6x+d, e=d/C,
Q0=1+3x, Qc=Q0+d, a=1+2x+delta/2,
F(t)=1+xt^2-2a/t+Qc/t^2-2d/t^3+d/t^4,
alpha0=0, F(v)=0, F'_r(R)/2=kappa=z_i nu,
ell=B/sqrt(C), q=k^4, k=2000, eta=1/2000,
z_i in [1/4,7/25].
```

Q0 is the extremal RN comparison coefficient. Qc is the coefficient
of the RN part of this rational function. Neither is substituted for
the fixed physical magnetic parameter Q in the retained theory.
The extra terms are part of the covariant seed's reference action.

The selected outer branch obeys
lambda_e+11G lambda_e^2/(360pi)=lambda0. The fixed common scale is
mu0=1/b0, so mu0 B is not set to one. Use the actual paired Dirac
domain, complementary-adjoint square and the regular HH state at the
seed's original prescribed kappa. Define

```
Theta(t)=B^2 tau_hat,mass(Bt),
Gamma_rest=Gamma_Dirac-Gamma_ref,
Mabs_rest=integral_R^infinity |rho2,rest(r)| dr,
Qrho=(G B^2/x) sup_(r>=R) |partial_r rho2,rest(r)|/r.
Qp=(G B^2/x) sup_(r>=R) |partial_r p2,rest(r)|/r.
```

The stress rho2 is sphere integrated: T4=T2/(4pi r^2). The
conclusions are

```
|Theta'(t)|<=Kq^2, K=1e15,
|Theta(2)|<=K0q^2, K0=1e11,
B Mabs_rest<1e16q^2,
sup |alpha_rest^(1)|<3.3e-10,
Qrho<1.3e-10.                                        (1)
Qp<1.3e-10, Qrho+Qp<2.6e-10.
```

The superscript (1) denotes the first constraint response to the
complete source evaluated on this seed, normalized by
alpha_rest^(1)(infinity)=0. The reference source has rho2+p2=0,
so the full Dirac and Gamma_rest first-order lapse responses agree.
The earlier complete horizon estimate remains available:

```
max(|8piG pA_rest,h|,|8piG pS_rest,h|)/(2kappa/R)
 <5.1e-12.                                          (2)
```

No metric-response estimate follows by replacing a spatial derivative
in (1) with a parameter derivative. The full-action boundary charge
also remains separate from the finite source integral.

## Explicit geometry guards for transferring the heat estimates

The product and temperature bounds give

```
x>=1000, |d|<1e-12, |e|<1e-24,
eta/(8k^9)<=h<=eta/k^9<1e-30, log(1/h)<80,
C/q<.0035<.00352, M_n^2=n(n+q)/C, M_1>18,
d_n=2(q+2n), |log(mu0 B)|<1e-21.
```

In particular a<2.01x, Qc<3.01x, a/C<.334 and Qc/C<.501.
The stronger .0035 guard used in the original non-bolt angular
budgets is explicit: C0<(123/20)k^3 and the accepted product
radius bounds give C/C0<1+1e-21, hence
C/q<(123/(20k))(1+1e-21)<.0035. The looser .00352 is used
only where the earlier mass and conservative ledgers already used it.
Let u=1/t, z=1-u and write the compactified orbit function as

```
calB(u)=u^2F(1/u)
 =x+u^2-2a u^3+Q0 u^4+d u^4(1-u)^2,
calB_ext=z^2[x(1+2u+3u^2)+u^2+d u^4],
calB=calB_ext-delta u^3.                              (3)
```

The bracket divided by C lies between 1/7 and 1 on [0,1]. For
the lower bound use x/C>=1/7. For the upper bound its derivative
is 2x+6xu+2u+4d u^3>0, and its value at one is C. Thus the
exact extremal factorization used by the non-bolt construction
survives, with no loss at the long-throat seam.

There is also an exact derivative control. If P=u^4(1-u)^2, then

```
calB_ext'=2u(1-u)[1-2Q0u+d u^2(2-3u)],
2C-calB_ext''
 =12(1-u)[x+Q0u]+d[2-P''(u)].                        (4)
```

The bracket in the first line strictly decreases from one because
its derivative is at most -2Q0+13|d|<0. It has one zero. In the
second line P''(1)=2 and |P'''|<=264, so
|2-P''(u)|<=264(1-u). Since 12x>264|d|, calB_ext''<=2C.
The elementary quadratic RN minimum is at least -(1+6x);
|P''|<=82 and C=1+6x+d then give calB_ext''>=-2C.

Consequently calB_ext'^2<=4C calB_ext: differentiate their
difference and use the single interior maximum from (4); the
difference is negative at zero and zero at one. The thermal term
has delta/C<=1.001h^2. With zc=exp(-16)/1001, z>=zc/2 gives
delta/(Cz^2)<5e-38. Directly, therefore,

```
Cz^2/8<=calB<=Cz^2, calB<=x+1<=C/5.                 (5)
```

The latter upper bound follows also from
-2a+Q0u+d u(1-u)^2<0. It is important that the shape term d
has not been compared with h as an unweighted number.

Here is the derivative ledger used in the frozen non-bolt proof.
All derivatives in this table are u derivatives except the last row.

| Quantity | Seed upper bound | Reason |
| --- | ---: | --- |
| abs(calB')/C | 5 | coefficient sum, a/C<.334,Qc/C<.501 |
| abs(calB'')/C | 2 | (4), thermal negative term and the lower margin |
| abs(calB''')/C | 17 | 12a/C+24Qc/C plus the degree-five/six terms |
| abs(calB'''')/C | 13 | 24Qc/C plus the degree-five/six terms |
| abs(calB^(5))/C | 960abs(e) | direct degree-six derivative |
| abs(calB^(6))/C | 720abs(e) | direct degree-six derivative |
| orbit abs(R),abs(R_y),abs(R_yy),abs(R_yyy) | 2,17,56,115 in ell units | the same polynomial formulas, including the new fifth derivative |

For the last entry the extra term is
sqrt(calB)calB calB^(5), in addition to
sqrt(calB)[3calB'calB''''/2+calB''calB'''/2]. Its normalized
size is below 960|e|. This fits inside the old 114.5-to-115 margin.

The logarithmic warping derivative w=f'/f consequently satisfies
the same real cap guards |w|<=1.01,|w'|<=2.021,|w''|<=13
before the taper. The exact same taper has
|f'/f|<=2,|f''/f|<=8,|f'''/f|<=100 and
R_even<=1e10. The smooth negative-u continuation has
R_smooth<=4e10. Its explicit bounds are unchanged: on [-1/2,0]
the degree-six addition changes the normalized polynomial and its
first six derivatives by less than 2000|e|. The strict RN margins
include

```
|calB'/(2C)|<3.523<3.54,
|calB''/(2C)|<7.202<8,
|calB'|/C<.751667<.752  on [-1/2,0],
|calB''|/C<3.507334<3.51 on [-1/2,0].
```

The 2000|e| correction is below 2e-21. Thus the actual analytic
strip, complex lower bound for the bridge Fbar, smooth Cauchy
remainder Cs=1e15, reflected remainder Cb=1e16 and cap row
constants in the non-bolt proof apply to this explicit seed. These
lemmas require the displayed guards, not a quartic field equation.
The two cap ratios remain distinct in the long-time and anchor
estimates.

The actual boundary remains geodesic because calB'(0)=0. The
balanced image trace cancellation is exact for the same paired
domain. The cubic coefficient is now

```
c=-a sqrt(x)/C^(3/2), |c|<1,                          (6)
```

since the new d terms start at u^4. The fourth-order Taylor and
two-derivative collar guards used for Cb follow from the same
analytic strip and relative f-jet bounds. No extra image sector
or physical boundary at the cap is introduced.

Actual-to-cap localization also transfers exactly. Equation (5)
gives a proper distance at least 16ell from each non-bolt
observation to the inner cap. The actual square is positive with
mass M_n^2, and the separated positive-square identity keeps that
actual inverse. Every constant in the frozen locality proof is a
cap row, cutoff or positive-inverse constant already bounded above.

## The cold patch and the common scale

Use the exact q0 perturbation and seed cap established in the
horizon-defect theorem. Relative to the RN polynomial at the same
x,h it has q_seed=q_RN+e r, with
|r|<4.2s,|r'|<4 and, on the old half s-disk, |r|<10s0.
The crucial resulting estimates are

```
|delta sigma|<=100|e|s0 L(s0)  on the quarter disk,
|delta sigma_s|<=400|e|L(s0),
|delta sigma_ss|<=3200|e|L(s0)/(s0+h).                 (7)
```

They follow by Cauchy with margin (s0+h)/4. The cold complex
normal-coordinate guards are preserved by the same explicit
degree-six derivative bound. The radial map gives |s_chi|<=2s;
from s_chichi=exp(2sigma)[sigma_s P+P'/2] and the q0
bounds, |s_chichi|<=3s for s<=.001. Its transverse Hessian is
s_chi coth(chi)<=2s. Thus the real transvection estimates
|sigma_t|<=200sL and |sigma_tt|<=4000sL retain large
margins: direct inner-patch bounds are 40sL and 500sL.
The perturbations (7) contribute less than 1e-16sL. At the distant
taper the explicit C2 perturbation is below 1e-14, with its same
sL exp(distance) weight, as in the horizon proof.

This verifies the hypotheses of the cold own-chart Cauchy remainder,
weighted actual-domain comparison and finite cap-gluing proof. The
proper gap M_1>18 and angular sums are unchanged. The only
quartic-specific local-polynomial identity acquires an explicit term:

```
partial_u(Rnorm^2+2Delta Rnorm)
 =-4 Bnorm' Bnorm''''-2 Bnorm Bnorm'''''.              (8)
```

The new term is below 2e-18 in normalized units and fits inside
the retained common-local budget 100q^2. The full cold bound
remains below 3e14q^2, and the non-bolt bound also remains below
3e14q^2. Hence K=1e15 is a single global source-gradient bound.

All common finite parts must be kept. With hat sphere B,

```
Fang_effective=4gamma_E-2H_q+4log(mu0 B),
a0 bracket=zeta_h'(-1)+[1+2log(mu0 B)]zeta_h(-1),
zeta_h(-1)=1/30-q^2/6, |zeta_h'(-1)|<1000q^2.
```

The effective Fang remains below 100. The anchor at t=2 has the
same smooth short-time, even-cap long-time, actual-domain and
angular EM bounds; its strict subtotal below 5e10q^2 has ample
room for these fixed-scale terms. This proves K0=1e11. No angular
mode is renormalized independently.

## Source reconstruction and elementary global geometry

Static HH regularity fixes the null integration constant. Let

```
I(t)=integral_v^t Bfun'(s)Theta(s) ds,
Bfun(t)=F(t)/t^2,
J(t)=integral_v^t F(s)Theta'(s)/s^2 ds.
```

The exact source formulas are

```
B^2 p2,mass=I/(2F)=Theta/(2t^2)-J/(2F),
B^2 rho2,mass=I/(2F)-Theta/t^2
             =-Theta/(2t^2)-J/(2F).                   (9)
```

F is increasing: t^3F'/2=xt^4+a t-Q0+d(-1+3/t-2/t^2)
has derivative 4xt^3+a+d(-3/t^2+4/t^3)>0, and its horizon
value is positive. Thus the anchor and gradient give exactly the
same individual bounds as before,

```
|B^2rho2,mass|,|B^2p2,mass|<=S q^2,
S=K+K0/2=1.00005e15.                                 (10)
```

For a useful lower bound compare with RN at the same x and v.
Put H(t)=(t-1)^2/t^4. The seed differs from that thermal RN
function by d[H(t)-vH(v)/t]. The bracket vanishes at v and its
derivative has magnitude at most 7(t-1). Integration gives
abs(bracket)<=7(t-v)(t-1). The RN factorization therefore gives

```
F(t)>=(x-7|d|)(t-v)(t-1)>=.99x(t-v)(t-1).             (11)
```

In particular F(t)>=xt^2/2 for t>=4. There is also a uniform
convexity bound .9x<=F''<=21x. For the lower bound maximize
4a/t^3-6Qc/t^4: its positive maximum is a^4/(8Qc^3)<.077x;
the extra inverse powers cost at most 44|d|. The upper bound is
at most 2x+6Qc+44|d|<21x. Taylor's integral formula then gives

```
(t-v)F'(t)/F(t)<=50.                                 (12)
```

Indeed F'>=F'(v)>0, F'<=F'(v)+21x(t-v), and
F>=F'(v)(t-v)+.45x(t-v)^2.

## The complete first-order lapse

The exact massive-source lapse derivative is

```
partial_t alpha_mass^(1)
 =-(G/B^2)J(t)/(t F(t)^2), alpha_mass^(1)(infinity)=0.
```

Monotonicity and (11) give

```
sup|alpha_mass^(1)|
 <=[G Kq^2/(.99B^2xv)]{log[v/(v-1)]-1/v}
 <[80*1.001/.99]K/k^8<3.2e-10.                       (13)
```

Here Gq^2/(B^2x)<1.001k^-8 follows directly from the product
radius bounds and lambda_e/lambda0; no parameter derivative is
used. The canonical LLL on this same seed has

```
Theta_LLL=q[-2+12a/t-12Qc/t^2+40d/t^3-30d/t^4]/(24pi),
|Theta_LLL|<=qC, |Theta_LLL'|<=qC/4.                  (14)
```

Its lapse is below 3e-28 by the same integral inequality.

For the Weyl piece the generic alpha=0 Euler derivative is

```
E_F={sigma r^2F''''+(4sigma+2)rF'''+2F''-2F'/r}/30.
```

Take its RN part with coefficient Qc and add DeltaF=-2d/t^3+d/t^4.
The exact integrated response is

```
alpha_W^(1)=G/(15piB^2){a/t^3
                  +Qc[(3/2)log t-13/8]/t^4}
 +Gd/(2piB^2){[22/25-(8/5)log t]/t^5
                                +[2log t-11/15]/t^6}. (15)
```

The new term is below G|d|/B^2 in absolute value. The entire local
Weyl term is below 2k^-10. The three allocations sum to less than
3.3e-10. Gamma_ref and the vacuum term both have rho2+p2=0,
so this is the complete Gamma_rest first-order lapse, not a sum
with the separately solved nonlinear LLL geometry.

## Quantitative outer trace and pressure integral

The same angular discrepancy estimate and reflected remainder now
have the explicit coefficient (6). Nothing in their derivation uses
the absence of the degree-five/six terms once the displayed guards
hold. Therefore

```
Theta_sing=c C^2/(5piY),
|partial_Y(Theta-Theta_sing)|<=2e13 q^2/Y, 0<Y<=1,
Y=sqrt(C)integral_0^(1/t) calB(u)^(-1/2)du.            (16)
```

The finite scale restoration is included in this statement. The
geometry (5) still gives sqrt(5)u<=Y<=sqrt(8)log[1/(1-u)],
so the unique Y=1 point lies at t1 in (2,4). Also

```
|calB(u)/x-1|<2.5u^2, 0<=u<=1/sqrt(5).
```

The additional normalized terms are at most
4|d|u^4/x, below the existing .099u^2 margin. The reciprocal
distance argument of the finite-mass proof hence gives

```
Theta(t)=Lt+Rtrace(t), L=-a x/(5pi),
|Rtrace(t)|<=Tq^2(1+log t), T=4e15, t>=4.             (17)
```

This is an explicit remainder estimate, not a differentiated
qualitative asymptotic expansion.

The derivative of Bfun is

```
Bfun'=-2/t^3+6a/t^4-4Qc/t^5+10d/t^6-6d/t^7.
```

On [v,4], t^5Bfun'/2=-t^2+3at-2Q0+d(-2+5/t-3/t^2)
is concave: its second derivative is below -2+28|d|<0.
It is positive at the bolt and at four. Moreover
Bfun(4)<=x+(9-13x)/256+9d/4096<x. Thus
abs(I(4))<=x(K0+2K)q^2. On t>=4 the exact coefficient bound is

```
|Bfun'|<4x/t^3;
2+(3/2)a+Q0/4+(1/2)|d|<4x.                          (18)
```

Decompose I using (17) without differentiating Rtrace. The old
RN primitive 2/t-3a/t^2+4Qc/(3t^3) receives the explicit
integrable addition -5d/(2t^4)+6d/(5t^5) for integral tBfun'.
Its extra absolute contribution is less than |L||d|. All strict
integral margins remain:

```
|I(t)|<4e15 xq^2,
|I/(2F)|<=4e15q^2/t^2,
integral_4^infinity |I/(2F)|dt<=1e15q^2.              (19)
```

The compact interval is bounded by 3Sq^2 using (10). The
canonical LLL bounds (14), the same I reconstruction and (18)
give its global energy-source norm below q^2.

## The non-RN Weyl tail, with a polynomial coefficient ledger

It is necessary to combine the massive and Weyl sources before
absolute integration. This can be proved quantitatively for the
rational seed, rather than inferred from its asymptotic form.
Let H be the braces of the compact-variation Weyl Lagrangian,
Rbar=-t^2F''+2tF'-2F and T=tF'. Put
Ebar_alpha=B^2E_alpha for the dimensionless Euler coefficient;
all primes in the following identity are t derivatives. Then exactly

```
t^2 Ebar_alpha=H-Rbar H_R
 -(t^2F'+2tF)H_R'-2t^2F H_R''
 -tF H_T'-(tF'-F)H_T.                               (20)
```

The RN part uses a,Qc,x in the frozen radial formula. Its source
is V(t)+L/t+RW(t), with

```
V(t)=11x^2t^2/(240pi),
|RW|<=2q^2(1+log t)/t^2,
integral_v^infinity |RW|dt<4q^2.                     (21)
```

The same coefficient sums 2076 and 3730 apply because a,Qc,x<q.
The extra DeltaF=-2d/t^3+d/t^4 does not change V or L.
Here is an explicit bound for its entire source, including all
quadratic cross terms. For the monomial exponents
n in {2,0,-1,-2,-3,-4}, write r_n=-(n-1)(n-2); then
|n|<=4 and |r_n|<=30. The following are coefficient norms per
absolute monomial amplitude, with a common (1+log t) factor:

| Entry | Bound |
| --- | ---: |
| quadratic H | 16 |
| linear H_R | 2 |
| t partial_t of linear H_R | 5 |
| t^2 partial_t^2 of linear H_R | 32 |
| linear H_T | 2 |
| t partial_t of linear H_T | 8 |

Equation (20) therefore has quadratic coefficient norm at most
16+60+30+64+8+10=188<200. Its linear coefficient norm is below
5; the F-independent magnetic term cancels in the difference.
The absolute RN amplitude is at most 9xt^2 and the extra
amplitude at most 3|d|t^-3. Divide (20) by 4pi t^2 to obtain

```
|B^2 Delta rho2,W|
 <=1e4 |d|x(1+log t)/t^3,
integral_v^infinity |B^2 Delta rho2,W|dt<1e4q^2.       (22)
```

Every monomial has inverse power at least three and at most ten;
differentiating its log-polynomial also gives

```
|partial_t(B^2Delta rho2,W)|/t
 <=1.1e5|d|x(1+log t)/t^5.                            (23)
```

This is a direct derivative estimate of the polynomial, not a
derivative of an unsigned remainder bound.

Combining (9),(17),(21)-(22), the renormalized full Dirac source
has its L/t terms cancel exactly inside the absolute integral:

```
B^2(rho2,Dirac-rho2,vac)
 =I/(2F)-Rtrace/t^2+RW+B^2Delta rho2,W
                                      +B^2rho2,LLL,  t>=4.
rho2,vac=11lambda_e^2 r^2/(240pi).
```

The outer trace norm is below 4e15q^2, and the massive compact
norm below 3Sq^2. On [v,4] the Weyl linear term must also be
included: its norm is at most |L|log(4/v)<q^2. Including this
term, the outer pressure 1e15, both entire
local remainders and the LLL, the full absolute norm is below
8.001e15q^2. The separate massive and Weyl mass integrals are
not assigned opposite divergent values.

## Exact reference cancellation and the finite remainder mass

Use the exact pressure profiles of the horizon theorem, with p0A,p0S
the complete center product pressures,

```
p_v=3(lambda_e-lambda0)/(8piG)
   =-11lambda_e^2/(960pi^2),
D_P=p0A+p0S-2p_v, |D_P|<3q^2/B^4.
```

Since rho2,ref=-4pi r^2pA_ref, its vacuum term is EXACTLY
rho2,vac. The remaining expression is

```
B^2(rho2,ref-rho2,vac)
 =-4pi B^4{(p0A-p_v)t^-2
                    +D_P(-t^-2+4t^-3-3t^-4)}.         (24)
```

The pressure bound |p0A|<3300q^2/B^4 and |p_v|<q^2/B^4
give its absolute global integral below
4pi(3301+12)q^2<5e4q^2. This proves the remainder mass norm
in (1) with its strict 1e16 margin.

In the first constraint response about S0=S_classical+Gamma_ref,

```
delta mu(infinity)=delta mu_h
                  +G integral_R^infinity rho2,rest(r)dr. (25)
```

Here delta mu_h is the Eulerian homogeneous datum at the fixed
background radius R. For a horizon displacement delta R it equals
kappa R delta R; the actual moving-horizon mass variation is
[mu_seed'(R)+kappa R]delta R. The companion
[first-order geometry construction](stationary-centered-first-order-geometry.md)
keeps these two variations distinct.

The reference vacuum has already been included in the seed's
lambda_e. Thus the integral is finite without an additional
cosmological retuning. The homogeneous horizon datum remains a
separate prepared-temperature condition. In particular

```
G Mabs_rest/(BC)<1e16 Gq^2/(B^2C)<7e-12.              (26)
```

This is a geometric constraint source, not the Hamiltonian of the
full renormalized action.

## Density variation and the remaining source components

The mass-source reconstruction differentiates once using Theta'
alone. From (9), exactly

```
partial_t(B^2rho2,mass)
 =Theta/t^3-Theta'/t^2+F'J/(2F^2).                    (27)
```

The anchor gives |Theta|<=4Kq^2t. Monotonicity bounds
|J|<=Kq^2F(t)(t-v)/(vt), and (12) bounds the last term. Thus

```
sup |partial_t(B^2rho2,mass)|/t <=31Kq^2.             (28)
```

No Theta'' or temperature derivative has been assumed.
For the RN Weyl remainder, direct differentiation of its finite
log-polynomial gives |RW'|/t<=20q^2 and |(L/t)'|/t<=q^2/15.
Equation (23) handles its extra seed tail. For the reference profile
(24), its derivative divided by t is below
4pi(2*3301+3*26)q^2<1e5q^2. The growing vacuum derivatives
cancel before this estimate. The LLL contribution is below 8qC,
by (12),(14),(27). These local pieces sum to less than 1e6q^2.
Therefore

```
Qrho< (31K+1e6) Gq^2/(B^2x)<1.3e-10,                (29)
```

as stated. The reference identity
partial_r rho2,ref/r=-8pi pS_ref is consistent with this
cancellation.

The analogous pressure seminorm also needs no second trace
derivative. From (9), exactly

```
partial_t(B^2p2,mass)=-Theta/t^3+F'J/(2F^2),
sup |partial_t(B^2p2,mass)|/t<=29Kq^2.                (30)
```

For the local Weyl pressure, in dimensionless notation,
B^2p2,W,ren=-B^2rho2,W,ren-F Ebar_F/(2pi), where
Ebar_F=B^2E_F. The explicit polynomial used in (15) gives

```
|Ebar_F|<=4q(1+log t)/t^3,
|t Ebar_F'|<=8q(1+log t)/t^3.
```

For example its RN terms have coefficients
2[6Qc log t-8Qc+3at]/(15t^4); the new terms are
d[(6-8log t)/t^5+(12log t-32/5)/t^6]. Their absolute
values and derivatives satisfy the displayed bounds using
a,Qc<q and |d|<1e-12. Also F<=5xt^2 and F'<=13xt,
so
abs(partial_t(F Ebar_F))/(2pi t)<16q^2. The Weyl
pressure derivative is therefore bounded by its density ledger
plus 16q^2. The vacuum terms again cancel before subtraction
of the reference; p2,ref=-rho2,ref. The LLL pressure has the
same canonical reconstruction, with derivative bound below 8qC.
All additional terms are below 1e6q^2, proving

```
Qp<(29K+1e6)Gq^2/(B^2x)<1.3e-10,
Qrho+Qp<2.6e-10.                                    (31)
```

These are fixed-seed spatial seminorms. They can be used in a
moving-horizon first-order constraint construction; they are not
renormalized metric-response norms.

The full angular pressure is fixed by the retained trace anomaly:
pS_Dirac=[A4-(p2_Dirac-rho2_Dirac)/(4pi r^2)]/2;
subtract the explicit pS_ref for Gamma_rest. In particular this
retains its radius source and the retained Box R4 term.
For reference B^2R4=-12x+4d/t^5-6d/t^6, so the extra derivative
anomaly is an explicit local integrable tail. The joint Ward identity
holds because both allocations come from the retained action split.
This supplies the full fixed-seed source needed in the first-order
constraints; it does not provide its Frechet derivative under an
unknown geometry/state change.

The companion [first-order geometry construction](stationary-centered-first-order-geometry.md)
uses these source bounds to establish the prepared-temperature linear
reconstruction and its regular finite completion. Focused verification
outcomes are recorded separately from these ordinary arguments.

Remaining to-do list: control the actual metric/state response and
retained-action residual in appropriate regularity spaces;
justify continuation or a full inverse; match the corrected exterior
and full boundary charge before the physical gates can close.
