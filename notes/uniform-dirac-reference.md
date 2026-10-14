# Uniform Dirac correction to the magnetic AdS reference

S074, 2026-10-13. The complete homogeneous Dirac determinant admits a
flux-uniform subtraction and pressure bound. In the specified one-loop
model, the candidate family has a unique stationary product within one
percent of its classical radii for every integer k>=2000. Quantitative
bounds put the actual fractional squared-radius shifts below 1.4e-22
for the sphere and 4e-26 for AdS2. This supplies a controlled reference
for later matching, not a wormhole solution or a stability result.

## Model and fixed scheme

Reuse the [homogeneous determinant](magnetic-product-vacuum.md): one free
massless four-component Dirac field, unit-charge potential F, magnetic
flux q, complete Euclidean H2 L2 kernel continued to the invariant AdS2
state, and zeta subtraction at fixed mu. The classical density is
(R-2Lambda)/(2kappa^2)-F^2/(4g_R^2), with kappa^2=8pi G and Lambda<0.
All additional finite local coefficients are set to zero in this
particular reference prescription, as in the retained determinant.
That is a specified input, not a determination of the full EFT's finite
coefficients. No photon/graviton loop or boundary coupling is included.

For the candidate set G=l_p^2 and

```
q=k^4, g_R(mu0)=k^-4, L4/l_p=k^5, Lambda=-3/L4^2,
C0=sqrt(1+12pi k^6), b0^2=L4^2(C0-1)/6,
a0^2=b0^2/C0, mu0=1/b0, k an integer >=2000.           (1)
```

Flux, mu0, Lambda and the renormalized coupling are fixed while varying
the two radii in each member. The family varies the theory's coupling
as k changes. This g_R is defined in the retained zeta scheme; matching
it to a separately defined electric coefficient or to the complete
boundary-driven model remains necessary.

The sphere spectrum is n(n+q)/b^2 with multiplicities q at n=0 and
2(q+2n) at n>=1. The hyperbolic density is
lambda coth(pi lambda)/(pi a^2). These are established inputs:
[Castellano, Montella and Zatti, (2.29)-(2.31), (2.84)-(2.85) at zero electric field, and (3.40)-(3.42)](https://arxiv.org/html/2603.15730v1).
Their magnetic parameter is q/2. The retained note independently checks
the hyperbolic normalization against Camporesi--Higuchi. We use the
minimal Dirac determinant, not the supersymmetric or electrically
unstable sectors of the source.

## Proper time at the magnetic scale

The old sphere remainder contained exp(q^2 s/(4b^2)), which is unsuitable
for a uniform large-flux estimate. Put

```
t=q s/b^2, r=b^2/(a^2 q), L=log(mu0^2 b^2/q),
S_q(t)=[1+2 sum_(n>=1)(1+2n/q) exp(-t n(1+n/q))]/(4pi),
h(t,r)=(r/pi) integral_0^infinity lambda coth(pi lambda)
                                  exp(-rt lambda^2) d lambda,
K4(s)=q^2 b^-4 k(t), k=h S_q.                          (2)
```

No mode is omitted. Define D=r partial_r and hP=(1-D)h. Direct
differentiation gives the positive pressure kernel
hP=r^2 t/pi integral lambda^3 coth(pi lambda) exp(-rt lambda^2).
The local product coefficients are

```
k(t)=a0/t^2+a1/t+a2+O(t),
a0=1/(4pi^2), a1=(r-1/q)/(24pi^2),
a2=(15-3/q^2-5r/q-3r^2)/(720pi^2).                    (3)
```

Here a0 in (3) is a heat coefficient, distinct from the reference radius
a0 in (1). For hP S_q the coefficients are
p0=a0, p1=-1/(24pi^2 q), p2=(15-3/q^2+3r^2)/(720pi^2).
For (D hP)S_q they are 0,0,r^2/(120pi^2).

For any of these channels define

```
I[k]=integral_0^1 [k-a0/t^2-a1/t-a2] dt/t
     +integral_1^infinity k dt/t,
F=1/2 [I[k]-a0/2-a1+a2(gamma_E+L)],
P=(1-D)F, w_E=q^2 b^-4 F,
pA=q^2 b^-4 A, A=-P,
pS=q^2 b^-4 S, S=P-a2/2, A+S=-a2/2.                 (4)
```

The pressure channel uses its own p0,p1,p2 in the same functional.
Its weighted derivative uses the third set above. The -log q in L is
required when changing the proper-time variable. It is a scale conversion
within the same subtraction, not a new finite counterterm. Equations
(3)-(4) reproduce the earlier anomaly and both radius variations.
They retain the complete LLL and higher-level determinant exactly once.

## A sphere remainder without exponential flux growth

Let f(n)=(1+2n/q) exp[-t(n+n^2/q)]. Euler--Maclaurin gives

```
S_N=1/(2pi t)+(1/(2pi)) sum_(j=1)^N B_(2j)
      sum_(v=0)^j (-1)^(2j-v) t^(2j-v-1) q^-v
                                     /[(2j-2v)! v!]. (5)
```

Use f=-(1/t)partial_n exp[-t(n+n^2/q)]. For m=2N+1, differentiating
the quadratic exponential and changing variable y=n+n^2/q in its
absolute integral yields

```
|S_q-S_N| <= zeta(2N)/[pi(2pi)^(2N)]
 sum_(j=0)^N sum_(v=0)^(N-j)
 [m! binomial(N-j,v) 4^v v!/(j!(m-2j)!)]
                     q^(-j-v) t^(2N-j-v-1).           (6)
```

Indeed 1+2n/q=sqrt(1+4y/q); its remaining even power expands into a
finite polynomial, whose Laplace moments are v!/t^(v+1). Applying the
periodic-Bernoulli bound gives (6). All coefficients are positive and
there is no exp(qt/4) penalty. The smallest t power is N-1.

The thermal expansion of h through t^M has coefficients
B_(2j+2) r^(j+1)/[(2j+2)j!pi], j=0,...,M, after 1/(2pi t).
Its remainder is at most

```
D_M r^(M+2) t^(M+1),
D_M=2(2M+3)! zeta(2M+4)/[(M+1)!(2pi)^(2M+4) pi].       (7)
```

For the pressure channel multiply (7) by M+1. For the coarse M=2
bound, the weighted pressure derivative has multiplier at most 12:
the thermal integrand is proportional to (2-x)e^-x, x=rt lambda^2,
and |(2-x)e^-x-(2-3x)|<=2x^2. This proves the required derivative
estimate without differentiating an unsigned remainder inequality.

For a numerical sphere cutoff Ncut with t Ncut(1+Ncut/q)>=80, monotonicity
past that point and the exact integral of f give
0<=S_q-S_cut<=exp(-80)/(2pi t). The retained numerical calculation uses
N=8, M=6 and t0=1/8 or 1/5. The resulting explicit UV and angular
tail bounds supplement adaptive quadrature estimates; they do not turn
the complete floating-point evaluation into an interval certificate.

## Uniform pressure and first derivative bounds

For q>=1, 0<r<=1 and any finite L, the subtracted integrals satisfy

```
|F|, |A|, |D A| <=100(1+|L|),
|partial_L A|<=1, |a2|<=1, |D a2|<=1.                 (8)
```

Here is a deliberately conservative proof. Split proper time at 1/4
and 1. On 0<t<=1/4, use S_3 and the order-two h, hP, D hP polynomials.
With pi>3 and zeta(6),zeta(8)<2, (6) bounds the sphere remainder by
(5155/69984)t^2<t^2, and the largest thermal remainder is at most
(35/4374)t^3<t^3. Both |S_3| and each exact hyperbolic factor are
bounded by 1/t there. The integrated product remainder is at most 1/2.
The absolute integral of each remaining polynomial after its three
subtractions is below 1, as verified by exact rational coefficient sums.

For 1/4<=t<=1, positivity and elementary exponential bounds give
S_q<1/t^2. Each of |h|, |hP|, |D hP| is below 1/t by its thermal
moments. Each local subtraction coefficient has absolute value below 1.
The absolute subtracted integral is therefore at most
integral_(1/4)^1 (t^-4+t^-3+t^-2+t^-1) dt<40.

For t>=1, S_q<1. Applying coth x<=1+1/x to the hyperbolic integrals
gives |h|, |hP|, |D hP|<=1/t+1/sqrt(t); for the derivative use
|(2-x)e^-x|<=(2+x)e^-x. The tail integral is at most 3.
The finite constants in (4) then give (8), with ample margin. The same
bounds justify differentiation on compact positive-radius domains.
The trace and scale bounds in (8) follow directly from (3).

## Quantitative stationary product theorem

Set z=b^2/b0^2, u=a^2/a0^2 and

```
alpha=1/C0, epsilon=kappa^2 q^2/(C0 b0^2)
                     =4 g_R^2(1+alpha),
r=(C0/q)z/u, L=-log q+log z.
```

The invariant current vanishes. The two Einstein equations, multiplied
by b^4/(C0 b0^2), become

```
0=-alpha z-(1-alpha)z^2/2+(1+alpha)/2-epsilon A,
0=z^2/u-(1-alpha)z^2/2-(1+alpha)/2-epsilon S.           (9)
```

Use A+S=-a2/2 to write the equivalent fixed-point map

```
Tz=1-epsilon A(r,L)-(1-alpha)(z-1)^2/2,
Tu=[1-alpha+alpha/z-epsilon a2(r)/(2z^2)]^-1.          (10)
```

On the closed box |z-1|,|u-1|<=rho=1/100, all k>=2000 obey

```
r<0.004, alpha<1e-10,
B=101(1+log q), |A|, |D A|<=B,
epsilon B<7e-23, alpha B<1e-6.                        (11)
```

For example C0>6k^3, C0/q<7/k and epsilon<5k^-8.
The functions k^-8(1+4log k) and k^-3(1+4log k) decrease in this range;
log(2000)<8 gives the last two bounds. On the box,
|partial_z A|,|partial_u A|<=2B. The denominator in Tu is at least
0.99. Thus

```
|Tz-1| <=epsilon B+rho^2/2<rho,
|Tu-1| <=[alpha rho/(1-rho)+epsilon/(2(1-rho)^2)]/0.99<rho.
```

The absolute derivative row sums are bounded by
rho+4epsilon B and (2alpha+3epsilon)/0.99^2, respectively; both are
below 0.02. Banach contraction therefore proves existence and uniqueness
within this box for the full stationary equations of the specified
one-loop truncation. This is quantitative and uniform in the family,
unlike an implicit-function statement at each separate flux value.

At the fixed point, (10) improves the bounds to

```
|z-1|<=2epsilon B<1.4e-22,
|u-1|<=2epsilon<4e-26.                                (12)
```

The second estimate follows from the first and alpha B<1e-6. These
are relative squared-radius bounds about the classical product. They
are not the plan's opening error budget.

The exact first-order shifts, with A0=A(C0/q,-log q) now denoting the
scaled pressure rather than a heat coefficient, are

```
delta z1=-epsilon A0,
delta u1=epsilon[a2(C0/q)/2-alpha A0],
|z-1-delta z1|, |u-1-delta u1|
                       <=100 epsilon^2 B^2<5e-43.     (13)
```

To see the remainder bound, put R=2epsilon B. The pressure Lipschitz
estimate and (10) bound the z remainder by 10epsilon^2 B^2.
Expand alpha/z about z=1, and bound
|a2(r)/z^2-a2(r0)|<=8R using (8) and the box. The u remainder is
at most alpha times the z remainder, plus
2alpha R^2+4epsilon R+2(2alpha R+epsilon)^2. Equations (11) imply (13).
This estimate is relative to the exact subtracted pressure integrals;
it does not certify numerical quadrature of those integrals.

## Numerical illustration and continuation

At k=2000, the scaled pressures from the common scheme are
A=0.0289351583549 and S=-0.0299905853617. The first-order shifts are
approximately -4.52112e-28 for z and +1.64910e-29 for u. These are
numerical illustrations; the existence and bounds (12)-(13) come from
the analytic estimates. The corrected product is still homogeneous,
so its invariant stress has zero contraction with an AdS2 null vector.
It does not supply the state-dependent opening stress.

[The checker](../checks/uniform_dirac_reference.py) verifies the
rescaled coefficients, rational remainder majorants, both contraction
rows, first-order Einstein identities and family bounds. Four existing
small-flux pressure records are recovered in the same scheme. Three
candidate values, each with two UV splits, and an independent radius
derivative check pass. Complete angular and UV tails are explicit;
remaining quadrature estimates are not interval-certified errors.

This advances the common reference/tadpole calculation in S074. The
unknown finite local inputs identified by the preceding audits remain
unknown. Adding them, changing the state, or matching the flat bubble's
electric coefficient requires the corresponding common action and finite
conversion. A global exterior, physical drive, nonlinear nonspherical
response, stability, transport and topology creation are not determined
by the stationary product theorem. Closed predecessor workstreams
remain closed.

Remaining to-do list: complete finite local and boundary matching,
independent Dirac terms, nonlinear constraints and remaining interactions;
then S075-S080 physical state, stress/geometry, stability, transport,
preparation/creation and final audit.
