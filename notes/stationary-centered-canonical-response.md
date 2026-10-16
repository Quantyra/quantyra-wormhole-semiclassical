# Canonical response on the regular centered comparison

2026-10-16, S139/S074. Ordinary argument before focused verification.

This note bounds the complete canonical determinant response, including
the massive angular tower, on a real C5 comparison. It retains the HH
ensemble, the complementary-adjoint paired domain and the common 4D
subtraction. The input is a quantitative finite-regularity metric
contract. Its application to the particular completed geometry requires
the separate metric-transfer proof; neither a trace-gradient bound on
the seed nor an undifferentiated logarithmic remainder supplies it.

The output includes a cold trace oscillation estimate. A value estimate
alone is insufficient at the bolt in the stationary Ward integral.
The Weyl/reference response and the Einstein constraint residual are
separate contributions, not included a second time in this canonical
estimate. No inverse of the full retained equations is asserted.

Here C5 is local on open orbit charts. At the outer boundary it
means precisely the WEIGHTED open-collar condition (2), not a
bounded C5 norm on the closed compactification. A Y^4 log Y
metric term is allowed: it is C3 at Y=0, while its proper-map
primitive is weighted C6 for Y>0 and ordinarily C4 at Y=0.
The fixed paired domain and parameter inverse need only the
lower boundary regularity and fixed projector stated below.

## Fixed data, directions and paths

Use B,A,C=B^2/A^2, q and the actual centered seed of
[the complete source proof](stationary-centered-exterior-source.md).
The unit of orbit length is ell=A. The angular masses and multiplicities
are

```
M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18,
Theta=B^2 tau_hat, Atilde(tau)=sum d_n M_n^2 exp(-tau M_n^2),
Atilde<=24C/tau^2, C/q<.00352, q/C<400.                 (1)
```

The orbit spin trace has rank two. The paired normalization is the one
already used in the source proofs; no additional factor two, half or
angular degeneracy is inserted here. In dimensionless orbit units the
canonical trace is C times the angular sum of the orbit traces. Thus
the response bounds below carry C^2, not q^2. The common scale mu0,
the sphere radius B in the product frame and the physical time period
are fixed throughout.

For each fixed original completion parameter ell_p in [0,1], compare
the seed and that endpoint. Make the endpoint proper-coordinate
identification once; do not differentiate it in ell_p. It is the
identity in the cold patch and is patched to the proper identification
near infinity. Interpolate the logarithms of the two endpoint metric
coefficients and of their radius fields on this fixed manifold.
The interpolation parameter below is tau_p in [0,1]. Its direction
is an endpoint difference. The proper tangent of the original moving
proper-coordinate homotopy is not substituted: its fifth derivative
could contain a sixth derivative of the original conformal tangent.
The endpoint quantum-source difference is independent of this
auxiliary path. Both endpoints and the entire path have the same
physical period and regular Euclidean bolt.

All local normalizations below are endpoint normalizations. Fix the
endpoint map once, pull the endpoint metric back, and interpolate
the endpoint logarithmic factors. Near the observation the actual
and UV-cap chords agree for every comparison parameter; away from
that patch interpolate positive relative metric matrices. Thus the
geometric-resolvent identity has the same path on its two sides.
Differentiating a moving normalizer of a merely C5 metric would
not prove that its direction is C5.

The sufficient contract is as follows. Epsilon denotes the size of
the chosen endpoint direction, and is proportional to ell_p. In the
compact middle region distinguish epsilon_hi for the C5 direction
from epsilon_lo for its C0/C1/C2 transport and global graph bounds.
There epsilon_hi may be 1e-4 while epsilon_lo is at most 1e-6;
these two quantities are not replaced by their maximum.

* In a real normalized orbit chart the common metric direction has
  C5 jets bounded by epsilon j!100^j. Its first spatial-center mixed
  jets through order four have one additional factor 100. The global
  relative C1 operator coefficients are bounded by 1e4 epsilon.
  Curvature stays below three in ell units. These are graph-norm
  coefficient bounds, including the spin connection and density.
  In the compact middle chart the first C5 bound uses epsilon_hi,
  while center value/gradient/curvature transport and the global
  relative C1 bound use epsilon_lo. In particular the large
  epsilon_hi bound is not used to perturb the global inverse.
* On the retained cold patch s=r/B-1<=.001 put
  a_x=s[1+log(.01/s)]. The auxiliary cap seam is at .01;
  it is not the endpoint of the retained cold observation region. The
  preceding local direction and mixed bounds have the additional
  factor a_x; global comparison-cap directions have the corresponding
  exp(distance) envelope. The direction size is at most 1e-6.
* Near infinity Y is inward proper orbit distance in ell units.
  After retaining the correlated mass cubic exactly, the relative
  inverse-warping direction v obeys

```
|v^(j)(Y)|<=epsilon j!100^j Y^(4-j)L(Y), j<=5,
L(Y)=1+|log Y|, 0<Y<=.1.                              (2)
```

  This is required along the fixed-endpoint proper cap chord. Its
  auxiliary cap extension has |v|+|v'|<=100 epsilon. The direction
  size is at most 1e-7. The cap-only global proper bound is not imposed
  through the moving physical bolt. The actual global inverse is
  treated in the fixed common chart instead.

Real C5 is enough. No complex spatial extension of a direction, C6
direction bound or analytic correction ball is assumed. The seed
coefficients used as comparison coefficients have their proved
analytic strips. The endpoint scalar radius direction has the
corresponding fourth-order contract for the explicit local sources.
Outer fourth and fifth derivatives in that contract are weighted
open-collar derivatives; they are not asserted bounded at Y=0.

## Actual domain and finite parameter dependence

The canonical half-density representation has a fixed paired domain.
In the proper collar its expression is

```
A=gamma_Y partial_Y+gamma_t f(Y) partial_t,
L=-partial_Y^2-f^2 partial_t^2-gamma_Y gamma_t f' partial_t.
```

The projector is constant, the complementary-adjoint square is D/N,
and f'(0)=0. Even reflection gives the exact spin-trace cancellation
of the image term. Its lack of smoothness across the cubic jet does
not change this domain statement.

For real comparison parameters the paired first-order realization is
anti-self-adjoint and ||(A+M)^-1||<=1/M. Scalar multiplication weights
preserve its domain. At weight M/2 its inverse is bounded by 2/M.
For a conformal direction u, the fixed-Hilbert-space identities are

```
A_u=exp(-u/2) A_0 exp(-u/2),
S_u=exp(u/2)(A_0+M exp(u))^-1 exp(u/2),
S_delta=(S u+u S)/2-M S u S.                          (3)
```

The endpoint and density terms in (3) must be retained. For a patched
diagonal metric the same argument uses its relative first-order
coefficient and connection bounds. Complexifying only the comparison
parameter on |mu|<=1e-8/epsilon changes these coefficients by less
than 2e-4 in graph units. The Neumann ratio is below .01 using graph
bound ten. This proves parameter holomorphy on the actual domain;
it does not evaluate a C5 spatial coefficient at a complex point.

On a real cap, local Caccioppoli and Fourier point evaluation give
the separated row 3e4 exp(-MD/2). Under the preceding complex
parameter perturbation, use the looser 3e5 exp(-MD/4). For the
actual-to-cap comparison it is therefore sufficient to retain

```
C0 exp(-MD/5), C0=1e23.                              (4)
```

To see the evaluation bound directly, use radius 1/(20M). The
evaluation cutoff has gradient 80M and Laplacian 11200M^2.
On its support |f^2-1|<.05; the larger energy neighborhood only
needs .07. The graph estimate gives gradient below 90M times
the neighboring L2 norm. The variable-principal absorption then
gives ||(-Delta+M^2)phi||<34000M^2||psi||. The whole-plane
point-evaluation norm is 1/(2sqrt(pi)M). Together with CT this
fits 3e4, including a possible square-root-two reflection factor.
The complex mass or relative-coefficient perturbation above is
far smaller than this reserve.

The two-sided geometric-resolvent identity uses only these cap
rows and the actual middle inverse. Differentiate that identity
in the fixed parameter. Parameter Cauchy costs 1e8 epsilon;
proper observation motion costs at most 100 epsilon times the
first spatial row bound 1e12 M; cutoff/connection variation costs
at most 1e4 epsilon. A permissive combined bound per mode is

```
1e62 epsilon M^5 exp(-9.6M),                         (5)
```

when the cap seam is at distance 24. The finite angular estimate
sum d_n M_n^5 exp(-9.6M_n)<1.2e-67 q gives an actual-domain
response below .005 epsilon C^2. No actual off-diagonal kernel,
Killing gap or unspecified interior load is used.

## Outer real-C5 heat response

For a variance-two unit bridge with supremum S, P(S>a)<=2exp(-a^2),
the covariant zero-winding diagonal is

```
k_f(tau;Y,Y)=E[J]/(2pi tau),
J=f(Y) F^-1/2 exp(tau G^2/(4F)),
F=<f(Y+sqrt(tau)b)^2>, G=<f'(Y+sqrt(tau)b)>.
```

The outside density factor f(Y) is essential. Differentiate this
formula first, then use crossing or exit events for bounds; no
indicator is differentiated. In particular

```
delta log J=v(Y)-delta F/(2F)
 +tau[2G delta G/F-G^2 delta F/F^2]/4,
delta F=<2f^2v>, delta G=<f'v+fv'>.                    (6)
```

Write zeta=sqrt(tau). For tau<Y^2 split off paths staying in
[Y/2,3Y/2]. On those paths (2) gives

```
|v^(j)(X)|<=12 epsilon Y^4L(Y) j!(200/Y)^j.
```

Use the finite Taylor-majorant radius

```
r=Y/[1000(1+S)].                                    (7)
```

The direction's finite jet norm is at most 15 epsilon Y^4L and
the v' jet through order three at most 3750 epsilon Y^3L.
The baseline coefficient moves by less than 1e-4 inside its
proved analytic strip. Consequently the finite coefficient cone
has Re F>.5, |G|<6 and H<2. This is a finite real-jet calculation;
only the baseline coefficients are analytically extended.

There is an important derivative-counting step. Do not take a fifth
zeta derivative of the ungrouped H(F,G): that would request v6.
Instead first write

```
H=F^-1/2+zeta^2 F^-3/2 G^2/4
             +zeta^4 F^-5/2 G^4/32+O(zeta^6).
```

For its parameter variation the three groups have coefficient
norms, including the outside density, bounded respectively by

```
300 epsilon Y^4L,
1e5 epsilon Y^3L,
2e6 epsilon Y^3L.
```

Expand the first through degree four, the second through degree
two and the third through degree zero. Their absolute remainders
use v5,v4,v2 respectively. The exponential remainder has coefficient
below 1e13 epsilon Y^3L and carries zeta6. Averaged odd polynomial
coefficients vanish under bridge reflection. Since
E(1+S)^5<3000, the leading expected fifth remainder is below
9e20 epsilon zeta5 L/Y before division by 2pi. All other grouped
remainders fit this reserve. Thus

```
|delta k_sub|<=1e21 epsilon tau^(3/2)L(Y)/Y,
0<tau<Y^2,                                        (8)
```

on the retained local paths, with k_sub denoting the common
a0,a1,a2-subtracted diagonal.

For tau>=Y^2, exact (6) and the weight-four direction give
vmax<=48 epsilon R^4 L(R)(1+S)^5 and
v'max<=2400 epsilon R^3 L(R)(1+S)^4, R=Y+sqrt(tau).
The bridge bound E[(1+S)^10 exp(S)]<1e9 then gives

```
|delta k|<=1e13 epsilon tau[1+|log tau|],
Y^2<=tau<=.01.                                    (9)
```

For an explicit moment check, integration of the bridge tail gives
E[(1+S)^10e^S]<=1+4[3072I0+2560I9+512I10], where
I0<1.3,I9=384,I10<1300; this is below 7e6.
The retained local variations satisfy
|delta k1|<=1e6 epsilon Y^2L and
|delta k2|<=1e12 epsilon L.

Take the target collar Y<=.05 and the direction collar |X|<=.1.
A fixed collar exit has S>.05/sqrt(tau). The global cap-direction
bound in (6) gives 1e10 epsilon tau^-1 exp(-.00125/tau).
Its complete angular integral is below
24e10 epsilon C^2/(.00125)^2<2e17 epsilon C^2.
The local unsafe-path estimate for tau<Y^2 is bounded by
1e13 epsilon Y^4L tau^-1 exp[-Y^2/(8tau)]; its integral
fits the same enlarged budget. No estimate is inferred by
differentiating an O(logY) remainder.

For tau>=.01 the real cap has

```
|delta k|<=1e5 epsilon R_even tau^-1(1+tau)exp(6.25tau).
```

The shifted mass-gap integral is below 1e7 epsilon R_even q.
With R_even<=1e14 and q/C<400 its contribution is below
4e23 epsilon C^2. Smooth negative-cap ratios are not used here.
The common local finite parts and long polynomial subtractions
add less than 1e16 epsilon C^2L. The common a0 restoration has
zero metric variation: its angular masses, B and mu0 are fixed.
This statement is not made about the entire homogeneous-product
determinant.

Combining (1), (8) and (9), the short safe integral is at most
48e21 epsilon C^2L. The large-short integral is below
96e13 epsilon C^2L^2. Together with the exits, long term,
restorations and (5), this proves the deliberately enlarged bound

```
|delta[Theta-c C^2/(5pi Y)]|
 <=1e24 epsilon C^2[1+log^2 Y], 0<Y<=.05.           (10)
```

The singular cubic is retained with the actual geometric mass and
is combined with its opposite Weyl density term before any physical
absolute mass integral. It is not independently varied or discarded.
Its bulk response beyond the singular profile must also be retained.
Write p=m1/(Bx), |p|<4e-11, and vary the outer analytic
coefficient calB_p=calB-2xp u^3 at fixed proper Y. Exactly,

```
v_m/p=xu^3/calB
 +[x calB'/(2sqrt(calB))] integral_0^u z^3 calB(z)^(-3/2)dz.
```

Hence v_m=p(x/C)^(3/2)Y^3+O(pY^5), cubic coefficient
less than .07p and no quartic term. On |Y|<=.4 and its
complex strip .01, |u|<.25, |calB|>x/2 and |calB'|<5C
bound the two terms by .03125 and .069. Use analytic
majorant ten. Taper only after .4, to zero by .5; the cap
C1 direction is below 100|p|. This is an auxiliary cap mode,
not a proposed physical horizon variation for p=.01.

The grouped heat argument with fixed radius 1/[1000(1+S)]
gives a smooth remainder below 5e22|p|tau^(3/2), whose
full short angular integral is below 3e23|p|C^2. After
K3 is removed, the reflected variation starts at weight five:
H_F+1/2 supplies V^2+tau, giving a bound below
1e15|p|tau^(3/2). Long and domain terms fit the preceding
4e23 reserve. Thus the mass mode's remaining bulk response
is below 1e24|p|C^2L(Y)^2. For the full direction (10)
therefore uses epsilon+|p|, rather than discarding everything
beyond its singular cubic. At the actual endpoint bounds this
is less than twice the estimate using epsilon=1e-7.

The thermal circle is fixed on the actual path. For a normalizing
coordinate rescaling its transformed period is retained in the FK
formula. The winding estimate exp[100tau-d^2w^2/(20tau)], d>=1e17,
combined with M1>18 bounds its first variation below
1e-100 epsilon q^2. Thus (10) is not a zero-temperature substitution.

## Real quadratic normal gauge and the cold local loops

The flat Born calculation needs only sigma_N(0)=grad sigma_N(0)=0.
It does not need a holomorphic direction or the exact Goursat
condition sigma_N(z,0)=0. In a baseline analytic normal chart let
nu be the perturbed real conformal factor. The holomorphic quadratic
coordinate map is

```
F(z)=a z+b z^2, a=exp[-nu(0)], b=-a^2 nu_z(0).       (11)
```

Indeed the transformed factor is nu(F)+log|F'|; its z derivative
at zero is a nu_z+b/a=0. The map uses only the center value and
first jet. Its first spatial-center mixed variation uses the second
jet, not a sixth metric jet.

Normalize the seed and endpoint separately, and keep their maps
F0,F1 fixed. The normalized difference is
delta nu(F1)+nu0(F1)-nu0(F0). Its first term uses only the
endpoint C5 difference; the second may use the analytic seed's
sixth jet. Interpolate these normalized endpoint factors. Extend
F1 globally once before interpolating the actual metric, so its
local chord agrees with the UV-cap chord. First-center mixed jets
through order four require five endpoint metric derivatives, not
the unknown sixth derivative of nu1.

The low Born cancellation is genuinely general-Hessian. With
A0=Delta sigma_N,H0=Hess sigma_N,B4=Delta^2 sigma_N, the one-
and two-vertex terms are

```
A0/(12pi)+(B4+A0^2+2trH0^2)/(120pi M^2),
                 -(2A0^2+trH0^2)/(60pi M^2).
```

Their sum at order M^-2 is (B4-3A0^2)/(120pi M^2).
Since R=-2A0 and Delta R=-2B4+4A0^2 at the normal center,
this is exactly -(R^2+2Delta R)/(480pi M^2). The traceless
Hessian cancels. The common per-mode anomaly and common a2
restoration are therefore the same as in the retained scheme.

Extend F-id with the frozen quintic cutoff on radii 5/32..7/32.
Its radial derivative bounds are 30,1500,3e5. The baseline map
has DG<1.439,DG^-1<=2,D2G<22 and D3G<5000. The latter follows
from |v_z|<.02/.018, |v_zz|<2(.02)/.018^2, giving F'''<300;
the cutoff adds less than 4000.
From the contract and (11), the map-direction bounds before
baseline composition are

```
(|G_delta|,|DG_delta|,|D2G_delta|)
               <(11,420,3e4) epsilon a_x,
(|G_delta,x|,|DG_delta,x|,|D2G_delta,x|)
               <(1100,5e4,3e6) epsilon a_x.
```

After baseline composition use (20,1e3,2e6) and
(2e3,1e5,1.2e8). For the mixed resolvent only DG_delta,x
is required in the divergence-form bilinear coefficient. There
is no hidden derivative of the baseline map of order six.

For the local potential v_N=exp(sigma_N)-1, keep |v_N|<.026.
The direction is quadratic, with coefficient below
2e4 epsilon a_x. Its global capped Taylor remainders can be
bounded by

```
|delta v_N-T4 delta v_N|<=2e10 epsilon a_x |z|^5,
|delta v_N,x-T3 delta v_N,x|<=3e10 epsilon a_x |z|^4. (12)
```

Inside radius 1/8 these follow by the real fifth/fourth Taylor
formula and the C5 contract. Outside it the bounded capped
potential and the degree-two through degree-four polynomials
give the same bounds after division by |z|^5 or |z|^4.
For products, use the enlarged coefficients 1e12 and 1e14
for the corresponding fifth and fourth remainders. All odd
retained loop polynomials integrate to zero before taking bounds.

The explicit flat Green kernel is
M[K0(Mr)+gamma dot x K1(Mr)/r]/(2pi). For an N-edge loop,
N>=3, total edge length l and p=4,5,6, the frozen Young/Holder
estimate generalizes directly to

```
integral l^p product(k_M)
 <=p!4^p 100 5^(N-3)/M^(N+p-2).                    (13)
```

For p=5 and 4 the coefficients are 12,288,000 and 614,400.
For two edges use coefficient 1000: K0<=K1 and
K1(s)=integral_0^infinity exp[-s sqrt(1+t^2)]dt
 <=sqrt(2)exp(-s/sqrt(2))/s give this bound by a gamma integral.
Equation (12), the one- and two-vertex terms and the tagged
complete Neumann sum consequently fit

```
|delta E_M,local|<=1e22 epsilon a_x M^-3,
|grad delta E_M,local|<=1e22 epsilon a_x M^-2.       (14)
```

For example Born two costs at most
2(12,288,000)(1e12) M^-3 and
2(614,400)(1e14) M^-2. Higher Born orders still have three
quadratic vertices and the degree-six loop bound; the tag costs
a polynomial in .026 and the direction's quadratic coefficient.
The Neumann denominator is 1-5(.026)>.87. Their sum is below
1e21 in either enlarged coefficient in (14).

## Cold finite gluing and its mixed response

The local flat cap and the global H2 cap are different metrics.
For this first/mixed-response argument use the original width-one
C2 H2 cap, with |sigma'|<.4 and |sigma''|<1.5. A narrower C3
cap used to obtain seed spatial derivatives is a different auxiliary
construction; its larger curvature is not assigned energy twenty.
The C2 cap is sufficient here: local C5 is needed only in the
unchanged inner metric, while the distant inverse needs C1
coefficients and the mixed conformal direction needs C2.
Use the exact two-sided resolvent identity and the transported
divergence expression B^a=exp(-nu)gamma^i(DG)^-1_i^a.
Freeze explicit conformal variation for (3), and charge the
coordinate transport once to B. The map bounds above give

```
|B_delta|<1e4 epsilon a_x,
|grad B_delta|<2e7 epsilon a_x,
|B_delta,x|<1e6 epsilon a_x.
```

For instance nu_delta<80 epsilon a_x,
nu_delta,j<2e4 epsilon a_x,
K_delta<4000 epsilon a_x and K_delta,j<9e6 epsilon a_x.
The mixed terms have nu_delta,x<8000 epsilon a_x and
K_delta,x<4e5 epsilon a_x. The baseline-center quadratic
terms are negligible on chi<=1, where a_x<1e-27.
Weighted graph energy twenty and the multiplier-only identity
(3) then give

```
M||W S_delta W||<1e6 epsilon a_x,
M||W S_delta,x W||<1e9 epsilon a_x.                 (15)
```

Here W=exp(-2distance). The uniform baseline conformal-center
source can be bounded by five. It is not replaced by the old
smaller bound three after changing a cap. Products use one
weighted direction and one uniform multiplier. No spatial
derivative of that multiplier is needed.

For this response calculation the local normal-cap row is the
commutator-weighted row R_M=||ev_x S_cap C_chi||, bounded by
400exp(-M/32). Each row includes ONE cutoff commutator, whose
slope is at most 33; the product includes the two commutators
once. Its parameter
and mixed variations can be bounded by 1e6 epsilon a_x R_M
and 1e9 epsilon a_x R_M. One direct finite check uses
||r^2G_flat||_(L2,weight M/4+2)<10/M^2; the preceding K1
bound gives this because sqrt(2)-1/2-4/M>2/3. Multiplication
by the quadratic direction, the cap inverse norm 9/M and the
cutoff slope 33 proves the first bound; the mixed quadratic
coefficient costs a factor 100 and fits the second. The
baseline-center tagged terms have the additional a_x<1e-27.

Retaining rank two and the external trace-M in the gluing identity,
(15) bounds its mixed scalar by
1e10 epsilon a_x R_M^2=1.6e15 epsilon a_x exp(-M/16).
The maxima M^2exp(-M/16)<150 and
M^3exp(-M/16)<5600 make it smaller than the reserves in

```
|delta E_M|<=K epsilon a_x M^-3,
|grad delta E_M|<=K' epsilon a_x M^-2,
K=K'=1e22.                                        (16)
```

The actual-to-H2-cap distant seam is handled by the same actual
domain identity as (5), with its original cold separation.
At chi<=1 it is bounded by a power of h of degree at least nine,
including its parameter and mixed variation. It fits (16) with
ample room. For chi>=1 its exp(distance) envelope retains a_x
for the VALUE bound. The mixed-gradient bound in (16) is needed
and asserted on chi<=1; there the explicitly indicated baseline
quadratic terms are negligible and give (18). The value proof
on the remaining cold patch uses the first endpoint variation
and uniform low coefficients, with no small-a_x^2 omission.
The common k2 finite restoration contributes below
1e13 epsilon a_x C^2 to the value and gradient; include it
inside the final scalar bounds below. Its flux/scale constants
are not set to zero.

## Cold angular sum, continuity and Ward use

Let N(E)=sum_(M_n<=E)d_n. Exact finite summation gives
N(E)<=4CE^2. Thus for R>=M1,

```
sum_(M<=R)d_n/M^2<=8C[1+log(R/M1)],
sum_(M>R)d_n/M^3<=12C/R.                            (17)
```

For two points separated by chi<=1 take R=M1/chi. The low
modes use the gradient in (16), the high modes twice its value.
After multiplication by C this gives

```
|delta Theta(chi)-delta Theta(0)|
 <=[(8/3)K+32K'] epsilon h L_c(h) C^2 sqrt(chi)
 <1e24 epsilon h L_c(h) C^2 sqrt(chi),
L_c(h)=1+log(.01/h)<81.                             (18)
```

We used chi[1+log(1/chi)]<=2sqrt(chi), an elementary maximum,
and a_x<=1.62hL_c(h)<2hL_c(h) on chi<=1.
Value alone gives |delta Theta|<=1e23 epsilon a_x C^2 after
enlargement: the complete angular inverse-cubic sum is below
2C/3, and K=1e22 in (16). The local restorations and constant product response
are included in these scalar bounds. A frozen constant Ward
reference has no oscillation; the full product determinant is
not declared metric independent.

In a fixed common conformal disk W=exp(2sigma)sinh^2chi,

```
P_hat=(2W)^-1 integral_0^chi W' Theta,
rho_hat=P_hat-Theta.
```

Separate the frozen homogeneous value before varying this formula.
Its pressure is exactly half that value for every W. The metric
kernel response of a nonconstant trace is

```
delta P=-2uP+(2W)^-1 integral W' deltaTheta
                         +W^-1 integral(Wu)'Theta.
```

In rho+p use the trace oscillation about the observation point.
Equation (18) makes the cold endpoint integrable, since the
remaining radial integral has at worst chi^(-1/2). A mere bound
|deltaTheta|<=constant epsilon h would not establish this fact.
Physical radius factors must be combined with the reference at
the same moving radius; they are not another canonical anomaly.

## Compact non-bolt coverage

The estimates do not leave the region between the cold patch and
Y=.05 uncovered. On r>=1.001B with Y>=.05 use a conformal
normal chart of radius R0=1/32, whose physical closure remains
before the boundary. For the NORMALIZED compact polynomial
calB/C the complex guards are |(calB/C)'|<7.08,
|(calB/C)''|<8 and |p|<.224. The Goursat integral equations
close |u-u0|<.02 and |sigma_N|<.002. The two boundary axes
cost 2(.224)/32=.014; their integral corrections are below
7.08e^.004/(4*1024)<.00174 and
8e^.004/(8*1024)<.00099. This chart is for the seed; the
direction is real C5 and normalized at the endpoints.

The same integral equation gives the sharper quadratic bound
|sigma_N(z)|<=2|z|^2, so a reserve c2=8 suffices for v_N
without dividing a constant bound by the shrunken radius.
Capped Taylor division gives value/mixed remainder coefficients
4e10 and 1e11; product reserves 1e12 and 1e14 suffice.
Take the gluing cutoff one through R0/2 and zero by 3R0/4,
with slope below 256. Both commutators are included. Their
Euclidean separation is 1/64, not 1/128. Since |sigma_N|<.002,
its proper separation is at least exp(-.002)/64>1/65.
The small evaluation ball plus
nested energy balls has radius below .19/M<=.01056, which
fits before that support even at M=18; the old radius 2/M
would not fit. The commutator-weighted flat row is bounded by
3200exp(-M/260), including its slope-256 commutator; its
quadratic exponent is exp(-M/130). This row convention differs
from the plain seed-evaluation rows below.
The maxima M^2exp(-M/130)<9200 and
M^3exp(-M/130)<3e6 give a two-norm bound. The local loops
use epsilon_hi. Born two costs at most
2(12,288,000)(1e12)epsilon_hi M^-3<2.46e19 epsilon_hi M^-3;
Born one costs less than 4e13 in that coefficient. For the complete
higher-Born tag take c2=8, delta c2=2e4 epsilon_hi and
delta v_infinity=1250 epsilon_hi. The geometric derivative

```
3 c2^2 delta c2/[1-5(.026)]
 +5 c2^3 delta v_infinity/[1-5(.026)]^2
 <9e6 epsilon_hi
```

multiplies the degree-six loop coefficient 2,949,120,000, hence
is below 3e16 epsilon_hi M^-4. These terms and local
restoration fit 1e22 epsilon_hi C^2 after the angular sum.

Gluing VALUE uses only the first parameter variation and low
normalizer jets. The bounds are R_delta<=1e6 epsilon_lo R_M,
M||S||<=1.2 and M||S_delta||<=1e6 epsilon_lo. With rank two
its coefficient is below 2e7 epsilon_lo R_M^2. For the new
R_M=3200exp(-M/260) this is
2.048e14 epsilon_lo exp(-M/130), and therefore less than
6.2e20 epsilon_lo M^-3. Both cutoff commutators and the
fixed-domain middle response are included in this enlarged row
bound. Since sum d_n/M_n^3<=12C/M1<2C/3, the complete
interior value response obeys

```
|delta Theta|<= [1e22 epsilon_hi+1e24 epsilon_lo] C^2. (19)
```

At epsilon_hi=1e-4 and epsilon_lo=1e-6 this is 2e18 C^2.
The actual middle inverse and its C1 direction are the same
fixed-domain ones used above. Patch transitions need no global
proper-coordinate interpretation at the moving bolt. The cold
and non-bolt regions overlap, and the outer estimate (10)
applies on the remaining boundary collar.

## Supporting fixed-seed source regularity

The metric-transfer argument needs seed source derivatives, in addition
to the endpoint response just proved. Here the directions are changes
of the observation chart on the ANALYTIC SEED, not derivatives of the
unknown completed metric or of its state. The following ledgers supply
that input without assuming it from a first radial derivative.

On the cold patch use the frozen H2 uniformization, s=h cosh(chi)
before its conformal correction, and L_c=1+log(.01/h). At the
seam s1=.01 the actual centered seed satisfies

```
|sigma|<.02, |sigma_chi|<.04,
|sigma_chichi|<.09, |sigma_chichichi|<3.             (20)
```

For the last inequality, on the complex s-disk of radius .0075
about .01 the exact divided-difference q0 has |q0-1|<.01172,
|q0'/q0-a|<.014. The uniformizing integrand J has |J'|<1.5,
J(h)=1, and its integral D has |D|<1.5log4<2.08.
The two logarithms in the conformal factor are bounded by .042
and .038, so |sigma|<.1 on that disk. Cauchy gives
|sigma_s|<13.334, |sigma_ss|<3556 and
|sigma_sss|<1.423e6. On the real seam the sharper first two
values are 3.35 and 416. The exact radial polynomial gives
s_chi<.01024, |s_chichi|<.011, |s_chichichi|<.013.
The three terms in the third chain rule are consequently below
1.423e6(.01024)^3+3(3556)(.01024)(.011)+13.334(.013)<3.
The degree-six centered-seed perturbation changes these strict
bounds by less than 1e-14.

For this SEED derivative calculation use a C3 cap, distinct from
the width-one response cap above. Multiply the cubic seam Taylor
polynomial by eta7(t)=1-35t^4+84t^5-70t^6+20t^7,
t=(chi-chi1)/.01, and set it to zero past t=1. Since
|eta7'|<=2.1875, |eta7''|<=26.25 and |eta7'''|<=210,
the cubic polynomial and its first three derivatives are bounded
by .021,.0411,.12,3. Leibniz gives the global cap guards

```
|sigma_cap|<.021, |sigma_cap'|<5,
|sigma_cap''|<6e4, |sigma_cap'''|<5e6.              (21)
```

Point-centered hyperbolic transvections turn these into weighted
observation derivatives. The uniformizing identities retain the
factor a_x=s_x[1+log(.01/s_x)] before differentiating. In the
inner cap each derivative of s costs at most its distance weight;
in the taper a_seam=.01 and the lower bounds for a_x exp(d)
give the same estimates. For first through third derivatives the
resulting multiplier bounds, with one distance weight reserved,
are respectively 500,7e6,1e9 times a_x exp(d). For example
the taper contributions are bounded by 5/.01,6e4/.01 and
5e6/.01; derivatives of the transvection vectors and the inner
uniformization fit the respective reserves 500,7e6,1e9.
This retains a_x, including at the bolt; a uniform cap derivative
bound alone would lose it.

Here the C3 cap is not assigned the old C2 curvature-based energy
bound. Its inverse is controlled by the conformal multiplier
factorization (3) on H2; the base massive weighted inverse is at
most 1.2/M. The required first-order graph bound can nevertheless
be rederived without the large cap curvature: for
T=D_H2+M exp(sigma),
||D_H2 T^-1||<=1+exp(.021)(1.2)<2.23. Distance-weight
commutators and the H2 Lichnerowicz identity give its gradient
bound below 2.5. Differentiating the endpoint factor exp(sigma/2)
adds at most (.5*5)(1.2)/18<.17. The physical gradient bound
is therefore below three, hence below the frozen reserve four.
The map bounds DG<1.439 and density connection 44(1.2)/18
then give transported energy below ten, for which twenty is
again a safe reserve. No sigma'' is used in this argument.
Derivatives
are finite ordered insertions of exp(sigma_cap). Explicitly

```
(T^-1)^(j)=sum_{l=1}^j (-1)^l
 sum_{i1+...+il=j, ia>=1} j!/(i1!...il!)
 T^-1 T^(i1) T^-1 ... T^(il) T^-1,                (22)
T=A_H2+M exp(sigma_cap).
```

Reserve one weighted insertion in every term. Use (21) for the
other insertions and the bounds 500,7e6,1e9 for the reserved
one. The endpoint multipliers in (3) are differentiated too.
For this ANALYTIC SEED calculation on s<=.001 retain the exact
Goursat normalizer, instead of differentiating a quadratic
normalizer along a C5 endpoint path. On the holomorphic axis
it obeys F'=exp[v(0)-2v(F,0)], F(0)=0. Its image has margin
.018 in the quarter disk. The uniformizing inverse is bounded
by two; differentiating the exact axis identity
h cosh(chi_y)=h cosh(chi_x)+(h sinh(chi_x)/2)z and using
Cauchy on the s quarter disk gives the deliberately enlarged
axis bounds |v_center^(j)|<=a_x(200,4000,1e9)_j.
There is also a direct global amplitude check. On the quarter
s-disk about s0<=.001 the exact logarithms obey
|log A|<=4.1s0 and
|log[cosh(hD)+(s/h)sinh(hD)]|
 <=s0[3.09 log((.01+h)/(s0+h))+6.85].
The integrand bound
|sD|<=s0[3.045 log((.01+h)/(s0+h))+6.75]<.014
gives |v_axis|<=12a_x. For
H(s)=h cosh(chi(s)), the same expression gives
|H'-1|<.05 and inverse derivative below 1.06.
A complex center displacement .01 combined with the axis
radius .232 changes H by less than .128(s0+h), within
the quarter disk. Therefore the normalizer is holomorphic
in that center disk with
|F-id|<=1.062(3)(12)(.232)a_x<9a_x.
Center Cauchy gives F_center'''<6e6(9)a_x and its
first/second z derivatives on the .018 margin; these are
strictly inside the third arrays used below. This check
retains the full nonlinear exponent, not only its linearization.
In the third differentiated normalizer equation the linear
source is v_center'''(0)-2v_center'''(F), bounded by
3e9a_x. Its coefficient of F_center''' is -2v_z,
|v_z|<.02/.018, while |F'|<1.05 and |F|<.232.
The scalar Gronwall factor is below exp(2(.232)(.02/.018))<1.68;
the integrated linear source is at most
.232(1.05)(3e9)(1.68)a_x<1.3e9a_x.
Nonlinear center terms are retained: on this region
a_x<=.001(1+log10)<.00331. Use one weighted factor and
uniform bounds for all other factors. This gives the safe
third-normalizer bounds

```
|F_center'''|<3e9a_x,
|(F')_center'''|<2e10a_x,
|(F'')_center'''|<3e12a_x.
```

The last bound follows either by differentiating its first-order
equation in z or by Cauchy on the remaining image margin.
Multiplication by the fixed quintic extension gives
(|G_center'''|,|DG_center'''|,|D2G_center'''|)
<(3e9,2e11,1e13)a_x: for the last entry the explicit sum is
3e12+60(2e10)+1500(3e9)<1e13. The first-center array
is the frozen (260,9200,6e5)a_x. For second derivatives use
the larger GLOBAL cold arrays
(|F_center''|,|(F')_center''|,|(F'')_center''|)
<(2e5,1e6,2e8)a_x and
(|G_center''|,|DG_center''|,|D2G_center''|)
<(2e5,7e6,6e8)a_x. To verify the new second bound, the
inhomogeneous second-exponent terms are at most
12000+4(200/.018)(2.6)
 +2[2(.02)/.018^2](260)(2.6)<3e5, times a_x.
Here 2.6 bounds 260a_x even on the larger range a_x<=.01.
The length and Gronwall factor give F_center''<2e5a_x;
z differentiation on the image margin gives the other entries.
The third equation uses these global bounds, in particular
|F_center''|<=2000, and its full inhomogeneous source is below
4e9a_x. Its integrated bound is below 1.7e9a_x, within
the displayed reserve 3e9a_x. These are SEED analytic
derivatives; using a quadratic map here instead would request
its fourth seed conformal jet, also available, but that is not
needed for this ledger.

For K=(DG)^-1 keep every inverse-matrix term. With
||K||<=2, |DG_center'|<=9200a_x<=92 and
|DG_center''|<=7e6a_x, the three groups in K_center'''
are bounded by 8e11a_x,
6(8)(7e6)(92)a_x and
6(16)(9200)(92)^2a_x. Their sum is below 1e12a_x.
The full chain rule similarly gives |nu_center'''|<1.5e10a_x,
using |grad v|<4, |Hess v|<400 and its next analytic
Cauchy derivative with the first/second global map arrays.
For its lower derivatives use |nu_center'|<2000a_x<=20
and |nu_center''|<1.2e6a_x. The exp(-nu) and Clifford
product rule then gives |B_center'''|<2e12a_x. Its direct
bilinear insertion costs less than 5e13a_x. The lower
transport arrays and multiplier insertions fit the safe array
M||S_center^(j)||<=a_x(2e6,1e11,1e18)_j,
j=1,2,3. To check the largest entry before coordinate transport,
the only cubic compositions in (22) have weighted factors
1e9, 3(7e6)(5), 3(500)(6e4), and 6(500)(5)^2;
their sum, endpoint factors and 1.2 powers are below 1e11.
The compact direct term is bounded above. The uniform first
inverse insertion is below 2e4: the transported B first
coefficient is below 700 and its energy factor is 24.
The uniform second inverse insertion is below 1e9 by its
ordered-insertion formula (the direct second coefficient costs
less than 2e7 and two first insertions less than 8e8).
Thus the mixed third groups cost at most
6(2e6)(1e9)a_x and 6(1e11)(2e4)a_x, and the triple
first group at most 6(2e6)(2e4)^2a_x. Their sum and
the direct term remain below 1e18a_x. No nonlinear
compact term has been discarded using near-bolt smallness.
No completed-metric sixth jet is involved.

The local UV part is analytic in (u0,p,pbar) on the frozen
polydisk of radius .01. Its complete flat Born remainder is
less than 1.4e11 M^-4, including every Born order and the
same a0,a1,a2 subtraction. Cauchy in the point-centered
parameters costs at most j!100^j, j<=3. Subtract the H2
value before taking the parameter variation: its exact
homogeneous value is invariant under the transvection. The
seed coefficients and their parameter difference have the
factor a_x from the preceding uniformization. Allowing 300
times 1.4e11 for the three parameter directions gives
4.2e13 a_x M^-4 after normalized differentiation. Thus
C sum d_n/M_n^4<4C^3/q bounds its angular sum by
1.68e14(C/q)^3 a_x q^2<7.4e6 a_x q^2. The finite local
restoration fits the enlarged local reserve 1e8 a_x q^2.

For the finite local/global gluing term, retain the two endpoint
rows, their normalizer densities and rank two. Equations
(21)-(22) and the local analytic row Cauchy bounds give its
third scalar derivative below

```
1e26 a_x exp(-M/16).                              (23)
```

One can check this by bounding the middle third derivative by
1e18 a_x/M and local center row derivatives by
j!1e4^j a_x times their commutator-weighted base row 400.
The ten ordered triples of derivative orders adding to three
have coefficients at most six. Retain one a_x in each term;
the other compact factors are bounded uniformly as above.
The largest middle term is 2(400)^2(1e18)=3.2e23, and each linear
endpoint term is below 2(400)^2(6e12). Even enlarging all
ten terms to the first bound and their coefficients to six
gives less than 2e25<1e26. The exact angular majorant
sum d_n exp(-M_n/16)<4q then makes the third derivative,
divided by 3!100^3 and multiplied by C, less than
1e26(4)(.00352)/(6e6) a_x q^2<2.4e17 a_x q^2.
The first and second distributions are smaller. The actual-domain
seam has distance growing like log(.01/h); its two cap rows
give h to a power at least nine after three derivatives and
are smaller than (23). No actual kernel evaluation estimate
is required there.

Thus the complete fixed-seed cold source, after its homogeneous
value is separated, obeys

```
|nabla^j(Theta_seed-Theta_product)|/(j!100^j)
 <=1e20 q^2 a_x, j<=3.                            (24)
```

The retained local anomaly and angular pressure have the same
finite rational-seed derivative bounds after using the exact
trace/radius Ward relation. Product-value subtraction is at the
same physical radii and fixed mu0; it is not a new scheme.

For the fixed-seed MIDDLE source take the radius-1/32 chart
already constructed above. Its cutoff is one through R0/2
and zero by 3R0/4, slope at most 256, source separation
D>=1/65 in proper distance. In this paragraph rows are PLAIN evaluation rows;
the two commutators are included explicitly below. The small
evaluation ball and nested energy neighborhoods have total
radius less than .19/M, so they fit at M>=18. Tangential
energy with nested widths 1/[20M(k+2)], k<=5, costs at
most (2000M)^k k!; the first-order radial equation supplies
the mixed jets with the same reserve. Starting from the
3e4 evaluation row gives, for j<=3,

```
||row_j||<=3e4 j!(2000M)^j exp(-MD/2).             (25)
```

For j=3 the differentiated product of the two rows is at most
4(3e4)^2 3!2000^3 M^3 exp(-M/65).
Rank two and the two slope-256 commutators give less than
2.3e25 M^3 exp(-M/65). Dividing by 3!100^3 leaves
3.9e18 M^3 exp(-M/65). The exact count N(E)<=4CE^2,
integrated against the absolute derivative of E^3exp(-cE),
gives

```
sum d_n M_n^3 exp(-cM_n)<=768C/c^5.
```

At c=1/65 this is less than 8.92e11 C. Multiplication
by C and (C/q)^2<1.24e-5 bounds the complete normalized
gluing derivative by 4.4e25 q^2. Lower derivatives are smaller.
For the local analytic UV term shrinking the frozen radius
1/4 to 1/32 costs at most 8^6 in its degree-six remainder;
the c2=8 higher-Born term costs less than its cubed ratio.
The complete coefficient 1e20 M^-4 is therefore sufficient.
Its center Cauchy factor j!100^j is the normalization already
used, and its angular sum is below 2e13 q^2. The local
polynomial restoration contributes below 1e12 q^2. Hence

```
|nabla^j Theta_seed|/(j!100^j)<=1e26 q^2,
j<=3, s>=.001, Y>=.05.                            (26)
```

This honest middle constant is larger than the cold and outer
ones; it is not reset to 1e20. The outer source counterpart is
the separately proved conormal result
[K_outer=1e20](stationary-centered-conormal-source.md), after
the correlated c/Y term is retained. Together (24),(26) and
that result are the fixed-seed source-three input for the
invariant second-order metric equations, not a claim of a
Frechet bound obtained by differentiating those equations.

## Scope and combination

The estimates prove canonical trace value control and cold
oscillation under the stated real C5 endpoint/path contract.
They give no uncontrolled proper-coordinate derivative of the
original family. Complete stress and radius forcing follow by
the stationary Ward reconstruction and the full retained anomaly,
after adding [the explicit local response](stationary-centered-local-response.md).
The outer mass cubic cancels with the Weyl term before absolute
integration, and the fixed reference vacuum is retained.

Application requires the actual metric-transfer bounds, the finite
integration ledger converting these trace bounds into mass/lapse/
horizon norms, and addition of [the E0 residual](stationary-centered-constraint-residual.md).
Until those are supplied this supporting response theorem is not
a claim of a controlled complete approximation, a stationary
solution or a Hamiltonian charge calculation.

Remaining to-do list: independently audit all finite constants and
the actual endpoint C5 transfer; combine the complete source response
and E0 residual in inverse-weighted norms; perform focused verification
only after the combined ordinary argument is accepted; justify an
approximation or higher-order inverse and complete the physical match.
