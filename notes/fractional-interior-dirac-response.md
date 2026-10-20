# Interior canonical response with four derivatives and a half-Hölder jet

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **A weaker sufficient interior
regularity contract; no preserved global nonlinear family or solution.**

The actual-domain canonical scalar response theorem of
[S154](divided-jet-dirac-response.md) remains valid when its divided C5
chart norm is replaced by the C^(4,1/2) norm below. Its first and mixed
response constants remain 1e-9 and 2e-7. This changes the spatial class,
not the physical action, angular tower, state, domain, clock or scheme.

The [S149 audit](spatial-response-norm-audit.md) already proves that
ordinary C5 bounds do not control a uniform spatial source gradient,
even for the exact flat comparison. Its logarithmic endpoint failure
is retained. The present result concerns source VALUES and finite
metric responses, and does not infer a spatial derivative of them.
It supplies another interior norm in which to construct the nonlinear
family required by S167; it does not prove the inverse preserves it.

## 1. Norm, actual-domain assumptions and statement

Use S154's fixed seed normal charts, observation cover, cap, original
paired domain and physical B,A,C,q, mu0. Keep its actual global and
capped value/gradient norms unchanged. On each convex real chart put

    [D^alpha u]_(1/2)
      =sup_(x!=y) |D^alpha u(x)-D^alpha u(y)|/|x-y|^(1/2).

Define P45(u) as the maximum of those global/cap low norms, the
coordinate value and first-jet norms, and

    sup |D^alpha u|/(j!100^j), 2<=j=|alpha|<=4,
    [D^alpha u]_(1/2)/(4!100^5), |alpha|=4.              (1)

All coordinates and derivatives have exactly S154's dimensionless
chart convention. Take the supremum over observations for a uniform
statement. Directions are real C^(4,1/2); only their finite parameters
are complexified. The identical actual-domain, kernel-free and low
coefficient conditions of S145/S154 remain explicit assumptions.
This is an interior conformal-shape theorem, not an arbitrary radius,
nonspherical, moving-domain or cold-endpoint theorem.

For P45(b)<=.001, U=P45(u), V=P45(w), use

    phi=sigma0+b+s u+t w,
    |s|<=1/(200U), |t|<=1/(200V).                       (2)

Zero directions are handled by continuity. The nonconstant canonical
scalar source T_hat, with S145's normalization 2G/B^2 per deformed
metric area and the same vacuum/contact allocation, obeys

    |D T_hat[b][u]|<1e-9 U,
    |D^2 T_hat[b][u,w]|<2e-7 U V.                      (3)

The source is holomorphic in the parameter bidisk and bounded there
by K<5e-12 after separating the fixed vacuum scalar. That vacuum is
not assigned zero value. The Lorentzian trace has the opposite sign
and the same absolute bounds. Along an admitted real segment, (3)
also gives derivative Lipschitz control and a Taylor remainder at
most 1e-7 P45(u)^2. These are source-value norm estimates only.

The new norm is controlled by the old divided P5 norm with constant
one. For a coordinate fourth derivative, put A=4!100^4 P5(u) and
L=sqrt(2)5!100^5 P5(u). Its difference is at most
min(2A,L|x-y|), hence its half-Hölder seminorm is at most sqrt(2AL).
After division by 4!100^5 this is less than .4P5(u), since
2sqrt(2)*5/100<.16. All other terms are already in P5. Thus

    P45(u)<=P5(u).                                      (4)

This inclusion does not manufacture regularity: the new class is
strictly larger. A compactly cut off |x1|^(9/2) in an interior chart
has a half-Hölder fourth derivative and an unbounded fifth derivative.
It can be made radial on a compact annulus and arbitrarily small in
P45, preserving the fixed domain and low coefficient requirements.
It is a test field, not a physical solution.

## 2. Fractional Taylor division on real spatial rays

Put m=M exp(phi(0)), h(z)=exp(phi_cap(z)-phi(0))-1. The unchanged
low coefficient estimates give

    Re m>.98M, |m|<1.02M, |exp(-2phi(0))|<1.03,
    |h|<=.06, |h(z)|<=2r, r=|z|.                       (5)

Every parameter field zeta=b+s u+t w has P45 norm at most .011.
On an inner real ray r<=1/64, retain the finite polynomial

    Z(X)=.016X+.011 sum_(j=2)^4 (150X)^j,
    E(X)=exp(Z(X)), through degree four,
    A(X)=1.026+.026 sum_(j=1)^5 (64X)^j.                (6)

Only the analytic SEED needs its fifth coefficient in A. The unknown
field has no fifth derivative in this proof. The bound 100sqrt(2)<150
converts coordinate derivatives to ray derivatives. The exponential
of zeta(z)-zeta(0) is below 1.001. As in S154, the normalized j-jets
of h for j=2,3,4 are bounded by coefficients of 1.001 A E, below
361, 45369, 6276299. The center Taylor coefficients obey

    |h1|<=.02r, |h2|<=300r^2, |h3|<=5e4r^3,
    |h4|<=7e6r^4.                                     (7)

The previous global integer remainder bounds through degree three
therefore remain valid: |h-T1h|<=400r^2, |h-T2h|<=5e4r^3 and
|h-T3h|<=7e6r^4. Their derivations use at most four metric jets.

Here is a direct fractional estimate for the remaining fourth
division. Write z_j=[X^j]Z and e_j=[X^j]E, with e_0=1. Bounds for
the half-Hölder seminorms of D_ray^j zeta/j!, j=0,...,4, are

    h_0=.016, h_j=(j+1)z_(j+1) for 1<=j<=3,
    h_4=4(.011)100^5=4.4e8.                            (8)

For j<4, bounded next derivatives give these bounds on a segment
of length at most one. For j=4, the sum of absolute direction
coefficients is at most (sqrt(2))^4=4, so (1) gives (8).

Let a_0=1.026, a_j=.026*64^j and ah_j=(j+1)a_(j+1).
The positive finite Bell polynomials for the exponential satisfy
partial e_j/partial z_l=e_(j-l). Thus the product/difference rules
give divided half-Hölder coefficient bounds

    eh_j=1.001[h_0 e_j+sum_(l=1)^j h_l e_(j-l)],
    hh_j=sum_(i=0)^j [ah_i(1.001e_(j-i))+a_i eh_(j-i)].

Their exact finite arithmetic gives hh_4<1e9. Taylor's integral
remainder, subtracting the fourth derivative at the center, gives

    |h-T4h|<=hh_4 r^(9/2)<1e9 r^(9/2), r<=1/64.

Indeed the coefficient multiplying the divided fourth Hölder
seminorm is 4 integral_0^1 (1-v)^3 v^(1/2)dv<=1.
For r>=1/64 use (5),(7). The coefficient is at most

    2.02*64^(7/2)+300*64^(5/2)
       +5e4*64^(3/2)+7e6*sqrt(64)<1e9.

This proves the global fractional division without a spatial analytic
extension. Interpolating the preceding integer remainder bounds gives
the full list

    |h-T_jh|<=B_(j+1) r^(j+1/2), j=0,...,4,
    B=(2,29,6000,840000,1e9).                          (9)

For j=1,2,3 the squared coefficients are bounded by
400(2+.02), 5e4(400+300), 7e6(5e4+5e4), respectively.
For j=0 use min(.06,2r). All are below the stated reserves.

## 3. Every Born order, with a convergent fractional angular tail

The local two- and four-derivative Taylor terms are unchanged from
[S145](general-interior-dirac-response.md). They use only derivatives
through four. In particular the common fourth term remains

    -Q_phi/(120pi m^2),
    Q_phi=Delta^2 phi-3(Delta phi)^2
             -4 grad phi.grad Delta phi+4|grad phi|^2 Delta phi.

It includes every first-through-fourth Born contribution and the
linear-gradient terms. No contact or common finite coefficient changes.

For the Taylor remainder of N insertions, use S145's product recursion
with p=(.02,300,5e4,7e6), A0=2 and the fractional B of (9):

    C(1,d)=B_(d+1), C(N,N-1)=2^N,
    C(N,d)=B_(d-N+2)2^(N-1)
                 +sum_(a=1)^(d-N+1) p_a C(N-1,d-a).

The exponent now is d+1/2. In the base case N=d+1, the full product
is bounded by min(2L,.06)^N<=2^N L^(N-1/2), separately for L<=1
and L>=1. In every recursive term, one fractional remainder and
the remaining Taylor or linear factors have total power d+1/2.
For total degree through four this gives

    C(1,4)=1e9, C(2,4)=4946800,
    C(3,4)=41990.4, C(4,4)=234.343432.                 (10)

For N>=5, keep the old first-five-linear-factors estimate
2^5 L^5 .06^(N-5); it needs no fifth spatial derivative.

The existing positive complex-mass kernel envelope has, for N>=2,

    J_(N,j)=integral L^j product_(N+1 edges) k
         <=14 j!10^j 5^(N-2)/M^(N+j-1), j=4,5.

Apply Cauchy--Schwarz to this positive loop measure:
J_(N,9/2)<=sqrt(J_(N,4)J_(N,5)). After the rank-two trace and
mass factors, its coefficient is bounded by

    28sqrt(4!5!)10^(9/2)<4.8e7,

and its mass power is M^-5/2. For N=1 use the sole insertion
distance r, rather than the total loop length 2r. The two-edge
fourth and fifth moments, including rank two, are below 10/M^4 and 20/M^5.
Their geometric mean is below 20/M^(9/2). These bounds follow
directly from S145's retained envelope: with r_m=1.02 the j-th
coefficient is at most (8/pi)r_m^2 Gamma(j)/(sqrt(2)*.98)^j.
The separate mass-insertion factor r_m^(N+1) is retained below.

Including the entire N>=5 geometric tail and the scalar-area factor,
the full remainder coefficient of M^-5/2 is at most

    1.03 {20r_m^2 C(1,4)
       +sum_(N=2)^4 4.8e7*5^(N-2)*r_m^(N+1)*C(N,4)
       +(1/4)*28*5!*10^5*2^5*5^3*r_m^6/(1-5r_m*.06)}
       <3e14.                                         (11)

For the last term M^-3<=M^-5/2/4, since M>=M1>18>16.
The Neumann ratio remains 1.02*.06<1, and the loop-tail ratio
5(1.02)(.06)<1. No finite Born or angular truncation is used.

The exact counting bound N(E)<=4CE^2 implies, by Stieltjes
integration and the vanishing counting function below M1,

    sum_n d_n/M_n^(5/2)
      <=(5/2) integral_M1^infinity 4CE^2 E^(-7/2)dE
      =20C/sqrt(M1)<5C.                              (12)

Thus the fractional remainder still has an absolutely convergent
whole-tower sum. The borderline exponent two would not suffice;
the half derivative is doing actual work in (12).

## 4. Actual domain, restoration and parameter derivatives

S154's actual-to-cap comparison uses only the unchanged low
coefficient norms, actual domain, separated rows and commutators.
Its per-mode scalar bound remains 7e11 M exp(-M/256). The common
fourth contact estimate also uses only the four retained jets and
remains below 1e17 C. No fifth unknown derivative is needed for
either estimate. Keep the fixed-scale vacuum and Polyakov allocation.

With eta_G=GC/B^2<8.5e-33, the complete nonconstant source bound is

    K<=2eta_G [5(3e14)+7e11*24*256^3+1e17]<5e-12.      (13)

The common subtractions, local finite contact and actual-domain
gluing are all present. The uniform majorants in (9)-(13) justify
coincident evaluation, full angular summation and parameter
holomorphy on (2). They apply directly to real C^(4,1/2) spatial
fields; no density assertion for smooth fields in the full Hölder
norm is required. Equivalently approximate in C4 with uniformly
bounded fourth Hölder seminorm and use the same summable majorants.
No complex spatial point of a nonsmooth field is evaluated.

Cauchy on (2) now gives 200KU and 40000KUV, proving (3).
Integrating the mixed bound along an admitted real segment gives
the stated Lipschitz and Taylor estimates in the observation
supremum source norm. This is a finite-background nonlinear
canonical component in the actual retained domain, not only a
flat multiplier estimate.

## 5. Consequence for S167 and limits

The existing candidate's fixed-chart P5 bounds imply the new P45
bounds by (4). A prospective interior nonlinear norm therefore need
not demand a fifth classical metric derivative just to use this
canonical source-value/response theorem. The inverse must still
control the chosen norm; P45 is not controlled by the small H1
correction established in the preceding milestone.

No spatial source gradient follows from (3); S149's obstruction
and the shape operator's fourth-order logarithm remain. A natural
mixed Hölder/logarithmic solver scale must be matched to its actual
source norm, rather than applying C5 value estimates as gradient
bounds. No norm-one property of Pi0 beyond its proved L2 forcing
space is assumed. The signed moment is still unevaluated.

This theorem does not remove S153's separate real cold endpoint,
cold direction, state/domain and clock contracts. It does not prove
their preservation for a P45 ball. The outer conormal and correlated
mass/cubic source cancellation, independent radius/local response,
full nonlinear map and candidate inverse require their own control.
All S167 acceptance tasks stay open. No physical wormhole, stability,
opening, charge, transport, preparation or creation result follows.

## Frozen verification contract

Freeze the note and checker before execution. Reconcile S145, S149,
S154 and the current S167 receipt/input identities. Check the P5
inclusion, finite exponential/difference recursions, fourth Hölder
coefficient, all integer and fractional division reserves, product
recursion, positive-measure moment interpolation, exact angular tail,
unchanged contact allocation and final parameter constants. Use an
explicit C^(4,1/2) cusp and independent finite polynomial products
to check fractional Taylor behavior; these are not physical Dirac
sources or continuum determinant simulations. Preserve failed
sources and every terminal run. Continuum resolvent, parameter and
domain arguments remain ordinary proofs; no Lean campaign is selected.

Remaining to-do list: match the correction and nonlinear source in
usable interior and cold/outer weighted norms; prove family
preservation and the projected full solve; evaluate actual signed
compatibility; exterior/state, charge and return; opening, stability,
transport, preparation, separate creation and full-result audit.
