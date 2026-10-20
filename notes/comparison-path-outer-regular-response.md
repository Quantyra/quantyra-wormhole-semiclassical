# Fixed-chart outer regular response on the retained path

2026-10-20, E003/S157 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The regular canonical outer source and its HH pressure contribution
have first and second responses on the ACTUAL comparison path in
fixed seed coordinates. Shrinking charts avoid differentiating a
moving proper-coordinate map in C5. The price is a sixth logarithmic
power, which remains integrable in the outer constraints.

The singular canonical cubic is allocated explicitly below. Its
combination with the Weyl sector, other local sources and the full
equation inverse remain required. A regular-sector bound is not a
bound for the complete outer physical source.

## 1. Fixed observations and scaled metric jets

Use [S155](comparison-path-interior-response.md)'s actual path
sigma_lambda=sigma0+lambda u and
r_lambda=B+(r0-B)exp(lambda V), 0<=lambda<=1.
Write y=Y0 for the FIXED seed inward proper coordinate and

    Y_lambda(y)=integral_0^y exp(lambda u(v))dv,
    L(y)=1+log(1/y), 0<y<=y_*=.08.                  (1)

The physical clock, actual paired domain, B,A,C and scheme are fixed.
No source is evaluated at a complex spatial point. The scalar at
fixed y is the one in the original common disk.

The finite metric proof gives the common-chart expansion

    u(y)=m y^3+v(y), |m|<1e-11,
    |v^(j)|<=2e-8 j!100^j y^(4-j)L(y), j<=5,       (2)
    |u|<1.2e-8, |grad_gbar u|<7.7e-6.

The factor two in (2) is the proper-to-common reserve in
[the metric estimate](stationary-completed-metric-response-norm.md).
The mass cubic follows from |p|<4e-11 and d^3<.075.
No sixth unknown metric derivative is used.

At an observation y use a fixed SEED conformal normal chart and
put r=y. Its cap is one through r/80 and zero by r/40. We need
the scaled high-jet norm

    max_(2<=j<=5) r^j |partial^alpha u|/(j!100^j),

with unscaled value/gradient, actual-domain and cap low norms.
Here is an explicit, loose finite-coordinate majorant for this step.
The seed a0 obeys |(log a0)'|<9 on its .01 complex strip and
.99<a0<1.01 on the real collar. In the isothermal coordinate
xi=int dy/a0, a fixed quadratic normalizer has coefficient
modulus <.1. On a complex real-coordinate bidisk of radius .001
about any point in |z|<=y/40, its inverse proper coordinate remains
within the .01 strip: the isothermal displacement is <.0042 and
|a0|<1.12 gives a proper displacement <.0048.
Only this SEED map is complexified.

Cauchy applied to its increment, of modulus <.01, bounds its
two-variable Taylor series by

    H(T)=20T/(1-2000T), T=|delta z1|+|delta z2|.

Indeed the grouped degree-j coefficient is at most
.01(j+1)1000^j <=20*2000^(j-1). Real points in the cap obey
.95y<Y0(z)<1.05y by the unit proper gradient and the seed
normal metric bound exp(|sigma0|)<1.03.
Thus the normalized scaled degree-j coefficients of (2) are
majorized, for j<=5, by the coefficient of X^j in

    4e-8 y^4 L(y) sum_(k=0)^5 [20X/(1-20yX)]^k
      +1e-11 y^3 [1.05+.2X/(1-20yX)]^3.          (3)

The first term includes the two-variable composition and the
variation of y^(4-k)L across the cap. Formula (3) has positive
coefficients increasing with y on this interval. At y=.08,
each coefficient for j=2,...,5 is below 2.6e-5.

The unscaled first jet is <8e-6. Also |u|<3e-8 y^3 on
the cap, so a cutoff of slope 160/r gives capped gradient
less than exp(.027)[8e-6+160(3e-8)y^2]<1e-5.
Consequently the scaled norm is <1.3e-4, uniformly down to
y=0. On the parameter ellipse used in
[S156](comparison-path-cold-response.md), rho=256,
|lambda|<=64.5009765625, and all these coefficient norms
remain below .011.

## 2. A complex bound that allows small M r

Reuse S154/S155's finite Born subtraction on the cap r/40.
The angular masses M themselves are NOT rescaled or truncated.
All Taylor-degree-j coefficient bounds acquire r^-j. The
complete degree-five remainder therefore acquires r^-5:

    |E_cap,M|<1e17 r^-5 M^-3.

The fifth-and-higher Born sum still has ratio 5(1.02)(.06)<1.
The complex-mass kernel envelope is the same one justified by
[DLMF 10.32.9](https://dlmf.nist.gov/10.32#E9), rechecked 2026-10-20.
The local common contact has total derivative order four and
is bounded by 1e17 C r^-4. All counterterms and fixed vacuum
allocation are the original ones.

The old point-evaluation ball 1/(10M) need not fit when M r
is small. Replace its inverse scale by

    K=M+32/r.

The evaluation ball 1/(10K) lies inside r/320, before the
gluing commutator. The flattened mass potential remains bounded
by 3M^2<=3K^2. S143/S144's same cutoff/Caccioppoli/Fourier
calculation gives point evaluation <1e4 K times the local L2
norm. The actual and cap weighted resolvents still have norm
3/M; their actual domain and scalar multipliers are unchanged.
The two separated row norms are therefore bounded by

    3.2e4(K/M) exp(-M distance/2).

The scalar-source comparison, including both rows, density
factor and cutoff slope 320/r, is bounded by

    (9e11/r)[M+64/r+1024/(M r^2)] exp(-Mr/320).   (4)

This treats ALL modes, including M r<1. It is not an
exponential-only estimate with an unfitted evaluation ball.
For delta>0, the retained count N(E)<=4CE^2 gives

    sum d_n M_n exp(-delta M_n)<=24C/delta^3,
    sum d_n exp(-delta M_n)<=8C/delta^2,
    sum (d_n/M_n) exp(-delta M_n)<=8C/delta.       (5)

The last two follow by Stieltjes integration by parts with
exp(-delta E) and exp(-delta E)/E. There are no masses below
18, so the apparent lower-end singularity in the latter has
zero counting measure.

Using delta=r/320, sum d_n/M_n^3<2C/3 and eta_G<8.5e-33,
the nonconstant normalized scalar has complex bound

    K(r)<=2eta_G {
       (2/3)1e17 r^-5+1e17 r^-4
       +9e11 r^-4[24*320^3+64*8*320^2+1024*8*320]}
       <2e-11 r^-5, 0<r<=.08.                   (6)

The bounds are locally uniform for every y>0, proving parameter
holomorphy there. They deliberately allow growth at the boundary.
They do not assert a bounded unweighted source on the closed collar.

## 3. Real regular values and actual-path response

Let S_lambda=2G Theta_lambda/(B^2 C) at fixed y and let c_lambda
be the cubic of the orbit warping in actual proper coordinates.
The exact finite cubic relation in the reconstruction gives
c_lambda=c0+lambda m. Define

    T_lambda(y)=2eta_G c_lambda/[5pi Y_lambda(y)],
    Q_lambda(y)=S_lambda(y)-T_lambda(y).           (7)

For complex lambda in the ellipse, (1) has |Y_lambda/y-1|<1e-6,
so its inverse and T_lambda are holomorphic. Here |c_lambda|<2.
The difference Q_lambda-Q_0 consequently has bound

    B(y)=1e-10 y^-5.                              (8)

The fixed vacuum scalar cancels in this difference. The cubic
is kept at its actual Y_lambda, including its parameter responses.
It is not frozen at y before differentiating.

For REAL lambda the endpoint canonical proof bounds the proper
regular difference by
1e24(1e-7+4e-11)lambda C^2[1+log^2 Y] when Y<=.05.
Pulling to y adds the frozen regular gradient contribution:
|R0,Y|<=2e13 q^2/Y and |Y_lambda/y-1|<1.3e-8.
On the remaining compact collar .05<=Y<=.081 the same endpoint
proof gives the unweighted coefficient 2e18 lambda C^2.
Its endpoint-coordinate identification is fixed once; the
frozen pullback correction obeys the same gradient bound.
The tiny explicit cubic difference on this compact collar is
included as well. Thus a single generous real bound is

    |Q_lambda(y)-Q_0(y)|<=5e-15 L(y)^2,
                  0<y<=.08, 0<=lambda<=1.        (9)

There is no omitted strip around Y=.05. On this full interval
L>=1+log(12.5)>3.5. The seed value bound
|R0|<=1e20 q^2 L and (9) also give |Q_lambda|<3e-7 L^2.

Apply S156's closed-interval lemma to (8)-(9). With t=log(1/y),

    log[B/(5e-15 L^2)]<=log(20000)+5t
                          <7(1+t)=7L,
    k<=7L/log(256).

It follows, INCLUDING lambda=0,1, that

    |partial_lambda Q|<6e-14 L^4,
    |partial_lambda^2 Q|<9e-14 L^6.               (10)

These are derivatives along the retained path, not arbitrary
mixed directions. They do not use S151's C7 coordinate contract.
They give C2 into the slightly weaker L^7 weighted supremum
space: each derivative divided by L^7 tends uniformly to zero
at y=0, while continuity away from zero follows from holomorphy.
The sharper L^6 bound in (10) remains available for integration.
No C2 claim into the L^6 supremum space is needed.

## 4. Fixed-cut pressure memory and the retained cubic sector

Use the fixed seed cut y_*=.08, which belongs to S155's interior
cover. Its HH pressure P_* has value/first/second bounds

    (P0,P1,P2)=(6e-10,5.2e-13,1.1e-14).

There is no moving-cut evaluation. With
W_lambda(y)=a0(y)^2 exp(2lambda u(y)) up to a fixed constant,
the exact Ward solution is decomposed as

    P=P_reg+P_cub,
    P_reg(y)=H_lambda(y)P_*+L_lambda[Q](y),
    P_cub(y)=L_lambda[T](y),
    H_lambda=W_lambda(y_*)/W_lambda(y),
    L_lambda[f]=-[2W_lambda(y)]^-1
                    integral_y^y_* W_lambda'(v)f(v)dv. (11)

P_cub has zero matching value at y_* by this definition;
the full, nonzero HH pressure memory stays in P_reg.
This is an exact allocation, not a new subtraction or boundary
condition. P_cub-T and 2P_cub-T must still be combined with
the correlated Weyl sector before a full physical absolute norm.
In particular T has nonzero second response from Y_lambda.

From (2), |u_y|<1e-7 y, and |a0'/a0|<2y. Write a=1.2e-8.
Then W(v)/W(y)<1.01 and
|a0'/a0+lambda u_y|<3y. The value/first/second operator bounds
on any weight L^p, p>=0, are

| Operator | Value | First lambda response | Second lambda response |
| --- | ---: | ---: | ---: |
| H_lambda | 1.01 | 5e-8 | 2.4e-15 |
| L_lambda | .01 | 8e-10 | 5.5e-17 |

For H use 4a H and 16a^2 H. For L, put d=u(v)-u(y);
the differentiated integrands are respectively
K[2d A+u_y(v)] and K[4d^2 A+4d u_y(v)], where
K=W(v)/W(y), A=a0'/a0+lambda u_y.
Using |d|<=2a and integral_0^.08 v dv=.0032 gives the table.
The weight L(v)^p is at most L(y)^p inside the integral.

The full product rule in (11) and L>=3.5 give

    |P_reg|<4e-9 L^2,
    |P_reg'|<5e-15 L^4,
    |P_reg''|<1e-15 L^6.

Thus the regular canonical normalized profiles
rho_reg=P_reg-Q and n_reg=2P_reg-Q satisfy

    value: 3.1e-7 L^2,
    first: 7e-14 L^4,
    second: 9.2e-14 L^6.                          (12)

Every changing-kernel and matching-memory cross term is retained.
The bound does not request a spatial derivative of a varying source.

## 5. Entire seed outer interval and its constraint consequence

For t=r0/B>=40 the retained seed guards give

    1/t<sqrt(5.9)/(1.001t)<=y
          <=sqrt(8)log[t/(t-1)]<.08.

Hence the whole outer interval is covered and L(y)<=L_e(t)=1+log t.
For the sixth logarithmic power, direct integration gives

    I6=int_40^infinity L_e^6/t^2 dt
       =(1/40) sum_(j=0)^6 binom(6,j) j! L_e(40)^(6-j)<1580,
    (2C/x)I6/40^3<.302,
    (2C/x)int_40^infinity L_e^6/t^5 dt<.019.       (13)

These replace only S152's OUTER column. With d_g=3e-7 define

    m2=1.003(9.2e-14)+4.02d_g(7e-14)+5.54d_g^2(3.1e-7),
    a2=1.006(9.2e-14)+10.08d_g(7e-14)+26.5d_g^2(3.1e-7).

Use L^2,L^4<=L^6 for the cross terms. The normalization cancels
one 2eta_G factor exactly as in S152/S156. For the OUTER REGULAR
CANONICAL component, the second-response bounds are

    zero-datum mass primitive: .302 m2/2 <1.4e-14,
    absolute mass density in BC units: 1580 m2/2 <7.3e-11,
    absolute lapse density L1: .019 a2/2 <9e-16.  (14)

Taylor allowances are half (14) times lambda^2. Pointwise
continuity and the integrable L^6 majorants justify twice
differentiating these density maps in their output norms.
The primitive norm follows from the same absolute integral
and its fixed positive outer denominator. No scalar derivative
of an absolute value or supremum is asserted.

## 6. Predeclared verification and limits

Revalidate S155, S156 and S152 receipts and current inputs.
Check the scaled finite-jet majorants, seed coordinate margins,
Born degree-five scaling, contact degree-four scaling, the K/M
row factor, complete angular counts and complex bounds. Check
the real/complex interval budgets, pressure product rules,
coverage and every sixth-log constraint constant.

Use 70 digits for numerical diagnostics. At y=.08,.01,exp(-20),
compare a manufactured Chebyshev family of degree
floor(log[B/a]/log(256)) and amplitude a/2 against the closed-
interval bounds. Normalize by a. Compare coefficient derivatives
using [DLMF 18.5.1](https://dlmf.nist.gov/18.5#E1)'s representation
with 128-node Cauchy contours of radius .005 at lambda=0,.4,1,
tolerance 1e-25(1+|target|); real first/second differences at
steps 1e-7,5e-8 must pass tolerance 1e-6(1+|target|) and
refinement ratio .35 unless both errors are below 1e-35.
Endpoint differences stay within [0,1].

For pressure use manufactured a0=exp(y^2/2),
u=1e-8 y^3+2e-8 y^4 L(y), a regular source
(3e-10+lambda*2e-13+lambda^2*4e-15)L(y)^2 and matching
pressure 5e-10+lambda*4e-13+lambda^2*5e-15.
At y=.001,.02,.06 and lambda=0,.4,1 compare direct 32/64-node
Gauss quadrature, after a logarithmic radial substitution,
against adaptive quadrature of the complete second-response
formula. Normalize by 1e-14, tolerance 1e-24 for quadrature;
use steps .001,.0005, tolerance 1e-7(1+|target|) and the same
refinement policy. Detect omission of matching memory and of
kernel derivatives. Check sixth-log integrals independently
after t=40exp(z). These fixtures are not quantum determinants.

Preserve every run, failed source snapshot, raw log and input
hash before/after. The infinite sums, analytic continuation and
operator estimates remain ordinary proofs.

Remaining to-do list: recombine the allocated cubic and Weyl
sector, complete other local/radius/reference sources and E0;
establish a general admissible family and usable full inverse;
solve exterior/state, charge, return and signed source; complete
opening, stability, transport, preparation and separate creation.
