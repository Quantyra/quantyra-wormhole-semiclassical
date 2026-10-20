# LLL and reference response on the retained comparison path

2026-10-20, E003/S159 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The massless canonical Dirac sector has an exact static HH stress.
Using that formula, its response and the reference density can be
bounded globally on the actual comparison path. In the cold region
the fixed massive homogeneous density must be included before applying
the vanishing weight. In the outer region the vacuum is already in
S158 and is not subtracted again.

This supplies the LLL/reference part of the density and null source.
It does not complete the Weyl cold/middle sector, the full radius
equation, E0 response, general source map or higher-order inverse.

## 1. Exact HH stress, with the state constant fixed

Retain the action, paired domain, physical temperature, B,A and common
scale of [the centered seed](stationary-centered-covariant-seed.md).
Use normalized profiles

    rho=(2G/(B^2 C))r^2 rho2,
    N=(2G/(B^2 C))r^2(rho2+p2),
    S=(2G/C)tau_hat, beta=Gq/(12pi B^2)
      =eta_G(q/C)/(12pi)<1e-31.

For the massless sector S=beta R, where R is product curvature
in A units. This follows from the retained q massless Dirac anomaly
tau_hat=q Rhat/(24pi). The anomaly and state integration are standard
two-dimensional mechanisms; see
[Grumiller--Kummer--Vassilevich](https://arxiv.org/abs/hep-th/0204253),
inspected 2026-10-20. The finite-path bounds below are derived here.

Write the dimensionless product metric as

    d rho_proper^2-a(rho_proper)^2 d tau^2,
    a(0)=0, a'(0)=1, w=log a, R=-2a''/a.

The last bolt derivative fixes the Euclidean period 2pi in the
original angle; it does not change the laboratory Killing clock.
The Ward equation and smooth HH bolt give exactly

    P=beta(1-a'^2)/a^2,
    rho_L=beta[a^-2+2w''+w'^2],
    N_L=2beta[a^-2+w''].                         (1)

Indeed (a^2 P)'=a a' S=-2beta a'a'', whose integral is
beta(1-a'^2). Any additional constant divided by a^2 is singular
at the bolt and is excluded. Formula (1) therefore contains the
same HH matching information as the Ward integral. No independent
outer pressure value or temperature is being fitted.

In the common H2 coordinate a=exp(sigma)sinh chi these become

    rho_L=beta exp(-2sigma)[1+2sigma''-sigma'^2],
    N_L=2beta exp(-2sigma)
                  [sigma''-coth(chi)sigma'-sigma'^2]. (2)

Primes in (2) mean chi derivatives. The apparent bolt pole is the
smooth tangential Hessian of a radial scalar. Its two Hessian
eigenvalues agree at chi=0, so N_L(0)=0.

For sigma_lambda=sigma0+lambda u, use seed proper derivatives D and
the seed w0. At the same seed point, (1)-(2) imply the exact formulas

    rho_L,lambda=exp(-2lambda u)
       [rho_L,0+beta(2lambda D^2u-lambda^2(Du)^2)],
    N_L,lambda=exp(-2lambda u)
       [N_L,0+2beta lambda(D^2u-(Dw0)Du)
                         -2beta lambda^2(Du)^2].     (3)

The two terms in D^2u-(Dw0)Du are the radial and tangential seed
Hessian entries. These formulas use only two fixed spatial jets.
One spatial derivative for a bolt modulus uses a third jet.

## 2. The reference and homogeneous allocations

Let t_lambda=r_lambda/B=1+s exp(lambda V), s=t0-1.
Keep the fixed coefficients p0A,p0S,pv,DP from
[the reference pressure proof](stationary-horizon-source-defect.md):

    pA_ref=p0A t^-4+pv(1-t^-4)
                         +DP(-t^-4+4t^-5-3t^-6),
    pS_ref=-p0A t^-4+pv(1+t^-4)
                         +DP(t^-4-6t^-5+6t^-6),
    DP=p0A+p0S-2pv,
    K=8piG B^2/C, f(t)=-1+4/t-3/t^2.

All coefficients, radii and subtraction scales remain fixed.
The exact normalized minus-reference density is

    Rref=K[p0A-pv+pv t_lambda^4+DP f(t_lambda)],   (4)

and its null profile is identically zero. The vacuum part is
K pv t_lambda^4. It is precisely minus the density paired with
the Weyl vacuum in S158.

Let rho_m,*=-S_m,*/2 be S156's FROZEN massive product density,
with the same A,B and mu0. At phi=0 the Weyl orbit stress vanishes.
The homogeneous product orbit equation consequently says

    rho_m,*+beta=-K p0A.                          (5)

It is a source identity, not an assumption that the massive
determinant's curvature response vanishes.

The additional density/null pair beside the existing canonical
allocations is therefore exactly

    cold (S156 massive deviation):
      H_rho=(rho_L,lambda-beta)
                       +K[pv(t_lambda^4-1)+DP f(t_lambda)],
      H_N=N_L,lambda;
    middle (S155 full massive canonical):
      H_rho=rho_L,lambda+Rref, H_N=N_L,lambda;
    outer (S158 massive+Weyl-minus-vacuum):
      H_rho=rho_L,lambda
                       +K[p0A-pv+DP f(t_lambda)],
      H_N=N_L,lambda.                             (6)

Add each line to its named predecessor and also include the still-open
cold/middle Weyl allocation. The resulting full physical densities
agree across the regional interfaces. The cut locations are fixed seed
locations, so no differentiation of the cuts is taken. The individual
allocations in (6) need not agree across the cuts. Their agreement
after all sectors are included is an exact algebraic identity, not
a claim that the missing sector bounds have already been proved.

In particular the cold constant p0A cancels identically. Bounding
the massive homogeneous density and reference separately would
introduce a false nonvanishing cold source.

The radius partner has not been discarded: pS_ref is displayed
above and obeys pS_ref=pA_ref+(t/2)pA_ref'.
For the canonical LLL the normalized angular pressure is -S_L/2.
These facts do not yet give a cold angular cancellation with p0S.
The Weyl angular product value is nonzero and must be included in
the subsequent full radius allocation.

## 3. Fixed-jet bounds and the cold bolt

All estimates below hold on the complex parameter disk |lambda|<=33,
while spatial arguments and region weights remain real and fixed.
There is no complex spatial evaluation of the unknown metric.
The actual low bounds give

    |u|,|V|<1.2e-8, |exp(-2lambda u)|<1.000001,
    |t_lambda/t0-1|<4e-7.

On the cold region use a_x=s[1+log(.01/s)]<=.00331 and
a_x<=2h L_h for chi<=1. The reconstruction gives |u|<=1e-6 a_x.
The fixed-seed covariant jet reserve

    |grad^j u|<=10 a_x, j=1,2,3                  (7)

is deliberately generous. Here is the finite-chart justification,
rather than a derivative of an endpoint inequality. The cold
barrier bounds the value and gradient by 1e-6 a_x and 2e-6 a_x.
The endpoint normalized jets are 1e-6 a_x j!100^j.
Invert the same quadratic normalizer as S155, at the observation:
its inverse linear deviation is <=2e-6 a_x, its quadratic
coefficient <=3e-6 a_x, and cubic coefficient <=2e-11 a_x^2.
The seed analytic normal metric has coefficient majorant
.026/(1-4T)^2. Compose these finite jets only through degree three.
The unknown contribution plus analytic seed/Jacobian contributions
has third ordinary derivative <7 a_x. At the chart center the
connection is zero. The explicit seed majorant bounds its
ordinary second metric-log derivatives by 3; all third-jet
covariant conversion terms cost at most 120|grad u|.
Thus 7a_x+120(2e-6 a_x)<10a_x. This uses the analytic
chart coefficients, not a curvature-only bound on a coordinate
connection derivative.
The lower jet bounds have much more room. No fourth or sixth
unknown derivative is used in (7).

For middle points, S155's fixed normal second jets and its low
gradient bound give |D^2u|,|(Dw0)Du|<10. For outer points,
the retained expansion u=m y^3+O_(4)(2e-8 y^4 L), |m|<1e-11,
gives |D^2u|<1e-4 and |(Dw0)Du|<1e-4.
The seed |Dw0|<y there, by the analytic seed polynomial.

We also need seed LLL values; a parameter response is not inferred
from a spatial derivative bound. In the cold region the exact seed
polynomial gives

    |R0+2|<=18s, |partial_t R0|<=17,
    |rho_L,0-beta|<=27beta s, |N_L,0|<=36beta s.  (8)

The first bound follows from calB''(1)=2C+O(C h^2)
and |calB'''|<=17C. The second row follows by the positive
Ward average with norm 1/2. Its null difference at the bolt is
bounded more sharply by 34beta(s-h). For the bolt modulus use
the exact s_chi=exp(sigma0)sqrt[(s^2-h^2)q0(s)].
The cold seed guards exp(sigma0)sqrt(q0)<1.1 imply
s<=h cosh(1.1chi). Hence s-h<=2h chi^2 on chi<=1.
Thus |N_L,0|<=100beta h chi on that neighborhood.

Through the middle the seed orbit warping increases: the explicit
calB' is negative for 1/40<=z<=1/(1+h). Together with |R0|<=2,
the HH average gives |rho_L,0|<=3beta and |N_L,0|<=4beta.
On the outer interval, the direct proper formula (1) gives the
weaker bound 20beta for either component: |R0|<=2,
|Dw0|<y and a0^-2<1. The cold bolt pole is not used there.

Insert (7)-(8) in (3). Since a_x<=.00331,
the brackets for the cold deviation and null are below

    1e5 beta a_x.                                (9)

This also accounts for (exp(-2lambda u)-1) times the product
density. For the bolt null profile, compare along a radial seed
geodesic. Both Hessian eigenvalues agree at the bolt and
|grad^3 u|<=10a_x; their difference is at most
44h L_h chi. The gradient-square derivative costs at most
200 a_x^2, and the segment length is <1.1chi.
The seed estimate and (3) then give

    |N_L,lambda|<1e5 beta h L_h chi
                       <=1e5 beta h L_h sqrt(chi). (10)

At middle and outer points the same insertion is bounded by
1e-26. These very small source bounds use q/C<400 and eta_G,
not a gap assumption for a massless resolvent.

## 4. Uniform holomorphy and the constraint bounds

Use |p0A|<3300q^2/B^4, |pv|<q^2/B^4, |DP|<3q^2/B^4 and

    gamma=Gq^2/(CB^2)<6.8e-28,
    K q^2/B^4=8pi gamma<1.71e-26.

In the cold region |t_lambda-1|<=1.000001s, so
|t_lambda^4-1|<4.1s and |f(t_lambda)|<10.1s.
In the middle |t_lambda|<41 and its inverse is <1.001.
In the outer region its inverse is <.026 and |f|<2.
Equations (6),(9),(10) therefore give the following complex
value bounds. The same weights as S152 apply, including its
original square-root bolt weight.

| Region/weight | Value bound B_i | First response | Second response |
| --- | ---: | ---: | ---: |
| cold / a_x | 1e-23 | 3.2e-25 | 2e-26 |
| bolt null / (h L_h sqrt chi) | 1e-24 | 3.2e-26 | 2e-27 |
| middle / 1 | 1e-18 | 3.2e-20 | 2e-21 |
| outer / L_e^2 | 1e-21 | 3.2e-23 | 2e-24 |

The finite covariant jet formulas and rational radius functions
are holomorphic into these weighted spaces. At the bolt the
Lipschitz estimate supplies a uniform stronger modulus.
Thus Cauchy on circles of radius 32 about each real lambda in
[0,1] applies in norm and gives the last two columns. There is
no endpoint interpolation loss for this explicit local sector.

Let d_g=3e-7. For each row i S152 gives

    m1_i=1.003 N1_i+2.01 d_g N0_i,
    a1_i=1.006 N1_i+5.04 d_g N0_i,
    m2_i=1.003 N2_i+4.02 d_g N1_i+5.54 d_g^2 N0_i,
    a2_i=1.006 N2_i+10.08 d_g N1_i+26.5 d_g^2 N0_i.

Apply half of the S152 matrix to these normalized profiles.
The combined LLL/reference/homogeneous allocation in (6) has

| Output | First response | Second response |
| --- | ---: | ---: |
| zero-datum mass primitive | 4.6e-18 | 2.9e-19 |
| absolute mass density /BC | 1.7e-20 | 1.1e-21 |
| absolute lapse density L1 | 2.2e-17 | 1.4e-18 |

Taylor allowances are half the last column times lambda^2.
These are density-map bounds, not absolute-value derivatives or
Hamiltonian charges. They supplement S155/S156/S158, and do not
complete the remaining Weyl/radius/E0 allocation or nonlinear solve.

## 5. Predeclared verification

Revalidate the S155, S156, S158 and S152 receipts and current inputs.
Use exact algebra to check the Ward antiderivative, the unique
regular bolt constant, the fixed-seed conformal transformation (3),
reference density/radius variation and identity (5), regional
allocation, finite chart jets, complex guards, Cauchy budgets and
all changing-geometry constraint products. Preserve every run and
failed source snapshot.

At 70 digits use manufactured radial data on chi>=0:

    sigma0=.02 chi^2/(1+chi^2),
    u=.003 chi^2/(1+chi^2)+.002/(1+chi^2),
    t0=1+.01 cosh chi, V=.004/(1+chi^2).

These deliberately larger directions test the algebra, not the
physical comparison metric. Set beta=1, K=1, p0A=7,
pv=-.2, DP=.03 and rho_m,*=-p0A-beta.
Test chi=.001,.2,1,4 and lambda=0,.4,1.
Compare exact HH pressure with independent adaptive and 32/64-node
Gauss quadrature of the Ward integral, tolerance 1e-28. Use eight
equal panels, with the stated node count on each panel.
Compare both exact transformed components with direct geometry
evaluation, tolerance 1e-35.

For each regional reference allocation, compute second responses
by symbolic differentiation, 64/128-node Cauchy contours of radius
.01, and real differences at .001 and .0005 (second-order
one-sided stencils at endpoints). Require absolute error below
1e-6(1+|target|) and refinement ratio below .35 unless both
errors are below 1e-35. A chi series checks the smooth bolt and
N=O(chi^2). Require omission tests to detect both a free 1/a^2
state term and failure to include the frozen massive density in
the cold cancellation. The complete reference Ward identity must
also pass independently of the orbit-only fixture.

These calculations test exact local/state and parameter algebra.
Uniform metric estimates and weighted-space arguments remain
ordinary proofs. No determinant sampling, new state, full branch,
stability result or physical completion gate is claimed.

Remaining to-do list: Weyl cold/middle and full radius allocations,
E0 and complete source/residual assembly; general admissible family
and usable full inverse; coupled exterior/state, charge, return and
signed-source solve; opening, stability, transport, preparation and
separate creation.
