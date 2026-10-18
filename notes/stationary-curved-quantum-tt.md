# The actual curved quantum TT form at the centered seed

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument before
focused verification. This continues the
[seed quantum-form calculation](stationary-seed-quantum-form.md).

At the analytic centered seed, in the retained common subtraction and
actual paired Dirac domain, the full action on the smooth bolt-regular,
compact-outer TT core and the stronger-norm completion specified in
section 11 satisfies

    I_full''[h,h] >= .99 E0[h]

for the static spherical TT tangents specified below. The proof includes
the complete massive tower, its local contacts and global domain, the
Weyl allocation, LLL, reference subtraction, and the off-shell return
from the disk chart. The growing outer logarithm is combined with an
even scalar spectral comparison; it is not bounded by a global constant.

This is a TT-sector lower bound at a comparison basepoint. It does not
exclude a kernel of the coupled trace/TT operator, evaluate its adjoint
phase or source moments, transfer an inverse to a nonlinear candidate,
or construct the physical exterior/state. S139.1-S139.2 remain open.

## 1. Fixed theory, basepoint and norms

Keep the actual centered seed, prepared period, laboratory clock, magnetic
flux, finite coefficients, scale mu0 and complete paired domain of the
[covariant seed](stationary-centered-covariant-seed.md),
[source theorem](stationary-centered-exterior-source.md) and
[domain audit](stationary-exterior-conformal-domain.md). In particular,

    I = I0 + Gamma_rest,   I0 = S_classical + Gamma_ref,
    Gamma_rest = W4 + Gamma_hat,LLL + Gamma_hat,mass - Gamma_ref,
    I0'(g0) = 0.

The seed is used for the FULL retained Hessian. No candidate regularity
is assumed. Let B be its sphere-center radius, lambda=lambda_e,
x=lambda B^2 and C=1+6x+d. Distinguish

    ell_o = B/sqrt(C),   ell_A = B/sqrt(x) = 1/sqrt(lambda).

The dimensionless canonical orbit is gbar=ghat/ell_o^2. Its outer inward
proper coordinate is Y, and u=B/r. With the time integral included, use
one physical half for

    J0[T] = integral_gbar T^2,
    J1[T] = integral_gbar |grad T|^2,
    J2[T] = integral_gbar (Delta T)^2.

The TT tensor has frame components (-b-2c,b,c,c), with b RADIAL. The
[global TT theorem](stationary-global-tt-reduction.md) gives, in the
physical radial norm with measure r^2 dr,

    (r^2 F b)_r = -(r^2 F_r-2rF)c,   ||b|| <= 2||c||,
    C0 = ||c||^2,   D = ||c'+2Rc||^2,   E = D+lambda C0,
    E0 = beta/(2G) [3D/2+lambda C0] <= I0''[h,h].

Here this prime is physical proper-radius differentiation, A=a'/a and
R=r'/r; the curvature Rbar below is a different quantity. The disk return is

    eta = -a integral_rho^infinity (b+c)/a d rho,
    S = c/2-R eta,
    T = -(b+3c)/2-(A-R)eta.

The preceding seed proof supplies

    ||eta||^2 <= 18 C0/lambda,
    ||S||^2 <= 100 C0,   ||T||^2 <= 500 C0,
    ||S'||^2 <= 5e4 E,   ||T'||^2 <= 1e5 E.

The weaker S bound 500 is also sufficient where indicated. It follows that

    J0 <= 500 beta C C0/B^4,   J1 <= 1e5 beta E/B^2,
    C J0/E0 < 5.2e-29,        C J1/E0 < 1.7e-27.          (1)

These use the actual benchmark bounds

    q=2000^4, M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18,
    C/q<.00352, q/C<400, C/x<6.1,
    GC/B^2<8.5e-33, Gq^2/(lambda B^4)<4e-27,
    |log(mu0 B)|<1e-21.

All inequalities below that use these small numbers concern this retained
benchmark, not every geometry satisfying only the coarse x>=1000 guard.

## 2. Actual double, normalization and positive comparison

The [paired modulus/phase allocation](stationary-paired-dirac-phase.md)
and [cold realization](stationary-reference-cold-bolt.md) supply an
anti-self-adjoint Dsharp with -Dsharp^2=L2>=0 and a preserved grading
J Dsharp J=-Dsharp. This is not an assertion that the original unpaired
Euclidean D2 is self-adjoint. Scalar cutoffs preserve the paired domain.

The [boundary-layer identity](stationary-massive-boundary-layer.md)
holds along every even static conformal path preserving the geodesic
seam and projector. The squared kernel is chirality-even and its image
insertion chirality-odd, so the image spin trace vanishes along the
whole path. Pointwise tr K_half=tr K_double, without a factor 1/2.
After integration, the physical functional and its variations are HALF
the evenly varied double functional. The common angular heat finite
part preserves this identity. The folded rank four is divided by two;
the orbit formulas here have rank two and d_n already counts both
nonzero sphere-Dirac signs.

The double has smooth bolts and a C2,1 metric at its cubic outer seam.
Its geodesic field on the unit tangent bundle is C1,1 and complete.
The Hamiltonian field of g^{ij}p_i p_j/2 has zero canonical divergence,
so its flow preserves Liouville measure dvol_gbar dtheta, with
integral dtheta=2pi. The Koopman group is strongly continuous and
unitary. Define its self-adjoint generator P=iX by that group; no
unproved smooth pseudodifferential assertion is needed.

Lift a base scalar T to this bundle. For one mass set

    B_M(r) = M^4/(2pi^2) [K0(Mr)^2-K1(Mr)^2] < 0,
    g(z) = 1-z^2/2+z^4/24-cos z.

Repeated integration of 1-cos z>=0 proves 0<=g(z)<=z^6/720.
The Bessel sign and moments follow from the positive-real integral
representations [DLMF 10.32.9-10](https://dlmf.nist.gov/10.32).
Squaring the Schwinger representation and integrating gives

    integral_0^infinity t^(s-1) Knu(t)^2 dt
      = 2^(s-3) Gamma(s/2)^2 Gamma(s/2+nu)Gamma(s/2-nu)/Gamma(s).

Thus the radial B_M moments I_j=integral r^(j+1)B_M dr are

    I2=-1/(6pi^2), I4=-4/(15pi^2 M^2), I6=-48/(35pi^2 M^4).

Define the nonnegative form

    P_M[T] = integral_0^infinity r dr [-B_M(r)] <T,g(rP)T>_SM.

It is bounded above on the H3 core by ||X^3 T||^2/(525pi^2 M^4).
The exact angular sum sum d_n/M_n^4<4C^2/q makes the tower finite
on that core. Its closed form domain is defined by spectral calculus
and compression by the bounded base lift. Closedness follows directly
from closedness of the spectral square root; no unproved identification
with a particular logarithmic Sobolev space is being used.

The same common fourth restoration as the seed calculation gives

    F_eff=4gamma_E-2H_q+4log(mu0B), kappa4=-F_eff/4>0,
    Qgeo[T] = sum d_n P_Mn[T] + Q4[T],
    Q4 = -C F_eff/(90pi^2) ||X^2 T||_SM^2
       = C kappa4/(30pi) J2
           - C kappa4/(90pi) integral Rbar |grad T|^2.       (2)

Indeed the angular Hessian identity is
integral (X^2T)^2 dtheta=(pi/4)[(Delta T)^2+2|Hess T|^2];
weak Bochner gives the second expression. It is valid at C2,1 regularity.
Physical-half bundle norms are likewise half the double norms.

There is no free fourth contact hidden in later comparisons. At each
finite angular cutoff, P_M plus its ordinary NEGATIVE fourth moment
exactly cancels the fourth Taylor subtraction, leaving cos-1+r^2P^2/2.
The scalar model below similarly reduces to log(1+aL/M^2)-aL/M^2
after its ordinary negative fourth term is added. Each finite-mode
actual bubble plus anomaly also has local order at most two. Their
common principal fourth finite parts agree. The remaining a2 difference
is the explicit curvature-gradient and zero-order polynomial

    C F_eff/(72pi) integral Rbar |grad T|^2
      + alpha2 integral (4Rbar^2+8Delta Rbar)T^2,
    alpha2=-C F_eff/(960pi).                              (3)

Consequently only second and zero contacts remain to be fixed. This
conclusion uses the common finite-mode allocation, not just an
off-diagonal kernel estimate.

## 3. Uniform normal caps and the gain of two orders

The actual degree-six seed retains the source theorem's bounds

    |Rbar|<=2, |Rbar_Y|<=17, |Rbar_YY|<=56, |Rbar_YYY|<=115,
    |w|<=1.01, |w'|<=2.021, |w''|<=13, w=f'/f, f=1/a,

on the non-bolt cover. Here gbar=dY^2+a(Y)^2 dtau^2.
These include d u^4(1-u)^2 and the thermal delta term, rather than
comparing the shape coefficient with the tiny horizon displacement.

Uniform seam-crossing normal charts can be constructed explicitly.
At Y0 put x=integral f/f0 dY, v=a(Y0)tau and sigma=-log(f/f0).
For |x|<=.07 the chart stays in |Y|<=.3 for |Y0|<=.1.
The exact derivatives are

    sigma_x=-w e^sigma,
    sigma_xx=-Rbar e^(2sigma)/2,
    sigma_xxx=-(Rbar_Y-2wRbar)e^(3sigma)/2.

Since (e^-sigma)_x=w, one may use |sigma|<.1 and derivative bounds
1.3,1.3,16. With b=sigma_x(0), the injective holomorphic map
H(z)=z-bz^2/2 removes the first jet:

    sigma_N=sigma(Re H)+log|1-bz|,
    sigma_N(0)=D sigma_N(0)=0,
    ||D^2 sigma_N||<5, ||D^3 sigma_N||<30, |z|<=1/16.

The third derivative is an essential bound, sufficient for C2,1 Taylor.
For v_N=e^sigma_N-1 this gives |v_N|<2.6r^2, sup|v_N|<.026 and
|v_N-T2v_N|<=6r^3. Capping after 1/32 and before 1/16 leaves
a global cubic remainder <=336r^3 and quadratic bound 8r^2.

The accepted cold and middle complex Goursat charts give the same
global reserve without a singular bolt frame. Restrict their analytic
v(z,w) to v(t e^(i theta),t e^(-i theta)). One-variable Cauchy at
radius R0 yields remainder <=2(.026)R0^-3 r^3 for r<=R0/2.
The cold R0=1/4 gives 3.328 and the middle R0=1/32 gives 1703.936.
Outside the inner radius, |v|,|T2v|<=c2 r^2 give 4c2/R0, at most
1024. Thus, everywhere on the actual seed, retain

    |v|<=8r^2, sup|v|<=.026, |v-T2v|<=1e4 r^3.            (4)

For the flat rank-two Dirac kernel use
k_M(r)=(sqrt(2)/pi)e^(-Mr/sqrt(2))/r. At exponential weight M/4,
its L1 norm is <=6.2/M and sup r^2 k_M<=.37/M. A quadratic-tag
Born insertion costs 8(.37)(6.2)/M, and each further insertion costs
.026(6.2)<.162. Including the endpoint factors gives

    |S_cap(0,y)-S_flat(0,y)| <=40/M e^(-M|y|/4).

These are COVARIANT kernels:
S_geom=e^(-sigma(x)/2)(Dflat+M e^sigma)^-1 e^(-sigma(y)/2).
The positive half-density factors belong to a different Hilbert
representation and are not substituted in this formula.

The bubble difference is <=1e4(M/r+1)e^(-Mr/2).
The exact count N(E)=sum_(M_n<=E)d_n<=4CE^2 gives

    sum d_n M_n e^(-M_n r/2)<=256C/r^3,
    sum d_n e^(-M_n r/2)<=(2304/324)C/r^4.

The local cap error is therefore <=1e7C/r^4. On r<=1/64 inside a
radius-1/32 controlled chart, a minimizing path stays within 1.02r;
|d-r|<=10r^3. The Jacobi equation with |sectional curvature|<=1.5
gives |Jexp-1|<=d^2/3 and |Jexp^-1-1|<=d^2/2. Since

    B_M'(r)=M^4 K1(Mr)^2/(pi^2 r)
           <=2M^2 e^(-sqrt(2)Mr)/(pi^2 r^3),

the coordinate/geodesic correction fits the enlarged 1e8C/d^4.
The usable separation is 1/64, not the full chart radius.

For radial comparison with a scalar model, freeze the warp at the
symmetric proper-Y midpoint. In coordinates with a_mid=1, endpoints
are (+/-s/2,+/-tau/2), r^2=s^2+tau^2. Write the minimizing path as
the straight path plus xi, Z=||xi'||_2. The local metric coefficient
b=a^2 has |b'|<3, |b''|<6. First-exit comparison gives Z<.75r.
The linear metric term on the straight path integrates to zero, and

    d^2 >= r^2+.95Z^2-7r^2Z-3r^4 >= r^2-16r^4,
    d^2 <= r^2+.25r^4, hence |d-r|<=17r^3.

The radial half-density factor sqrt(a(Y)a(Y'))/a_mid differs from one
by <.3s^2. The distance/Jacobi correction is <=30M^2e^-Mr per mode.
Adding it to the cap comparison and integrating tau gives

    |radial separated error| <1e9 C/|Y-Y'|^3.             (5)

All these comparisons have a finite local radius. Global returns,
thermal state and the actual resolvent are handled separately next.

## 4. Actual gluing, with two endpoint buffers

Let A(z)=(Dactual+z)^-1, B(z)=(Dcap+z)^-1 and Q=[D,chi].
The geometric resolvent identities imply

    chi A chi = B chi^2+B Q chi B-B Q A Q B,
    B Q chi B = (chi^2 B-B chi^2)/2.

The middle term has zero diagonal, and vanishes off diagonal when both
points have chi=1. Thus E=A-B=-BQAQB on that local pair. This is a
first-order identity; it does not discard the direct term in a different
positive-square IMS identity.

The local pair radius is now 1/128. Choose chi=1 through 1/64 and zero
by 3/128, with slope <260. Both endpoints have proper commutator
separation >1/132. Shrink the established evaluation/energy balls by
two: .19/(2M) fits even at M=18. Its row reserve, including the
commutator, is

    Rloc <=7.8e7 e^(-M/270).

For completeness, shrinking by k changes the cutoff bounds to
80kM and 11200k^2 M^2, the graph bound to 90kM and the second-order
evaluation bound to 34000k^2 M^2. The whole-plane Fourier evaluation
norm remains 1/(2sqrt(pi)M). Thus a k^3 row inflation for k>=2 is
conservative; the local inflation ten exceeds the required eight.

The point bound is ||E(x,y)||<=Rloc^2/M, but its full L2 row is only
260Rloc/M^2. Use these different bounds: the cap covariant L1 kernel
norm is <=10/M, so both rank-two cross terms cost <=40Rloc^2.
The E-squared term costs <=2(260)^2 Rloc^2/M^2<418Rloc^2.
It would be invalid to square the large point bound over the whole
domain. A reserve 1e3Rloc^2 gives angular local Schur bound <1e24C.

For the complementary FAR region use a different cutoff, one through
1/256 and zero by 1/128, with slope <512. Its proper separation is
>1/263. Shrink the evaluation neighborhoods by four and inflate the
plain row by 64. The commutator row is

    Rfar <=9.8304e8 e^(-M/530).

The one-sided identity ev_x A 1_far=ev_x B Q A 1_far and ||A||<=1/M
give physical rank-two bubble Schur bound 2Rfar^2. Summing
sum d_n e^(-cM_n)<=8C/c^2 gives <2e24C. The total bubble Schur
reserve is therefore 3e24C.

Reverse rows have the same bounds: adjoints and the preserved grading
return the same realization. In a paired half-chart the direct and
reflected L1 envelopes integrate over complementary halves of the
plane, so their sum has the same 10/M reserve. Rank four divided by
two gives the same physical rank-two factors.

Use the ROW-SUM-SUBTRACTED gluing form

    -.5 integral Kglue(x,y)[T(x)-T(y)]^2,

whose absolute cost is <=6e24 C J0. The actual constant probe below
already includes its row sum. Counting only the raw 3e24 bubble while
also retaining the full actual constant contact would miss this local
row correction. Gluing is kept as a bounded form, not included in the
smaller ultraviolet second-contact reserve.

## 5. Actual constant contact, including the tadpole

For S=(Dsharp+M)^-1 and P=M^2-Dsharp^2, grading cancels odd terms.
The exact constant conformal Hessian is

    M^2 Tr S^2-M Tr S = 2M^4 Tr P^-2-2M^2 Tr P^-1.

Its angular heat weight is

    B(t)=2[t sum d M^4 e^(-tM^2)-sum d M^2 e^(-tM^2)],
    |B(t)|<120C/t^2,
    |B(t)|<=800C t^-2 e^(-162t).                         (6)

The count proves sum d M^2 e^-tM^2<=12C/t^2 and
sum d M^4 e^-tM^2<=40C/t^3; extracting half the gap proves the
last inequality. The constant 2D anomaly has zero second derivative.

The finite local terms must be restored in the SAME common scheme.
With k0=2 and k2=-(Rhat^2+2Delta Rhat)/120, let

    A0=k0[zetah'(-1)+(1+2log(mu0B))zetah(-1)]/(4pi B^2),
    A2=B^2 F_eff k2/(4pi).

The constant-contact restoration is -2A0+2A2; the k1 term is zero.
This follows from the action: Gamma_a0=-A0/2 scales as e^(2T),
and Gamma_a2=A2/2 scales as e^(-2T). It does not introduce a new
finite +2 into F_eff.

Apply (6) to the complete cap decomposition of the
[non-bolt proof](stationary-reference-nonbolt-stress.md), using its
actual centered-seed transfer. Its smooth short-time remainder is
<=1e15 t^2, contributing <1.2e15C. The reflected cubic satisfies
|K3|<=2|c|sqrt(t)e^(-Y^2/(4t)), |c|<=1, and contributes <960C/Y.
The undifferentiated reflected fourth remainder is
<=1e16 t e^(-Y^2/(4t)), contributing
<240e16C(1+|log Y|).

Here the fourth-remainder estimate retains the outside density and
the six-term bound 1.3e11(V^4+tV^2) on crossing paths. With
V<=2sqrt(t) sup|sqrt(2)bridge| and its eighth moment <=48 it has
the stated reserve. Outside Y=.1 use the two-cap crossing bound
2e10 t^-1 e^(-Y^2/t) and the stronger cubic Gaussian
2sqrt(t)e^(-Y^2/(2t)); their spare Gaussian absorbs the t^-2 ratio.
This is not a smooth-double heat expansion through the cubic seam.

For t>=.01 the cap value bound is <=4e14 t^-1 e^(6.25t).
The gap part of (6) contributes <6e6(4e14)C.
Long subtracted polynomials, the fixed finite terms and the retained
thermal winding estimates are smaller. In Theta=C times the
dimensionless orbit-density convention these cap contributions fit
1e18q^2(1+1/Y).

The cold cover uses its complete normal-cap Born remainder, not a
time-Fourier estimate at a polar bolt. The real source remainder is
<=1.4e11 M^-4. On |z-M|<=M/4 the Bessel integral gives
Re z>=.75M, |z|/Re z<=5/3. The Neumann ratio increases from .13
to less than .362; including all retained mass/edge powers and the
denominator ratio gives inflation <2046. Thus 3e14M^-4 is a safe
complex source remainder. Constant mass differentiation by Cauchy
costs four, giving 1.2e15M^-4. Summing and multiplying by C costs
<2.1e8q^2. Its exact quadratic/quartic terms are the same common
a1/a2 terms; the complete loop remainder includes the homogeneous
H2 contribution and all higher curvature powers.

For actual-to-cap constant-probe gluing use the diagonal identity
of section 4, now holomorphic on |z-M|<=M/4. The actual inverse is
<=4/(3M). Repeating the REAL-metric graph proof with mass coefficient
1.25 and doubled 34000 reserve permits the conservative complex row
R(z)<=7.8e8e^(-M/260), with the original center-only separation.
For F(z)=tr(A(z)-B(z))(x,x), rank two gives
|F(z)|<=(8/3)R(z)^2/M; use eight instead. The Hessian difference is
-M^2F'(M)-MF(M), bounded by 5M sup|F|. Its sum is <3.3e24C,
or <4.1e19q^2 in Theta units. This is Cauchy applied to an actual
holomorphic localized identity, not differentiation of a value bound.

Altogether,

    |Theta_p0| <=1e20 q^2(1+1/Y), Y>0.                  (7)

The distributional Delta Rbar seam term is retained; T(0)=0 makes
its quadratic form vanish on the selected physical directions.
Since Y>=sqrt(5)B/r and ||T||^2<=500C0, the pointwise part costs

    beta 1e20 q^2 integral r^-2(1+1/Y)T^2 dr
       <=1500e20 Gq^2/(lambda B^4) E0 <6e-4 E0.          (8)

## 6. Finite second contacts and their form estimate

Test a centered normal quadratic q_v(z)=chi(z)(v.z)^2, with chi one
through 1/64 and zero by 1/32. Use 1-I_s(6,6) on the transition.
Its first five radial cutoff derivative reserves are
(11,220,3960,63360,887040). After Cartesian conversion the coefficients
are (1,22,968,42768,1907136,86359680); the quadratic's normalized
C5 norm, divided by j!100^j, is <3. This is a finite test, not a
form estimate inferred from a general C5-to-C0 source bound.

The direct local Born calculation uses (4). After the exact anomaly
and a2 terms, the one-Born fifth-degree remainder is <3e8/M^3.
The two-Born fifth loop has coefficient 12,288,000/M^6 before
mass vertices and spin rank. Two tagged positions AND rank two give
4(1e4+512)(12,288,000)/M^3<6e11/M^3.
The remaining v^2q and higher terms have degree at least six; summing
their tagged Neumann series fits 1e13/M^4. Since sum d/M^3<2C/3
and sum d/M^4<4C^2/q, the finite ultraviolet second contact fits
1e13C.

The known curvature contact is retained explicitly. At a conformal
normal center set A=Delta sigma_bg. The actual a2 action on the
quadratic is C F_eff A Delta q/(20pi); the geodesic action is
C F_eff A Delta q/(45pi). Their difference is
C F_eff A Delta q/(36pi), exactly the gradient term (3).

The local geodesic quadratic remainder must also be bounded after its
known fourth jet is removed. On flow length r<=1/128 the curve stays
within normal radius 1/64. The uniform normal-cap metric jets give
|Gamma|<=15, |D Gamma|<=6000, |D2 Gamma|<=3e6, so its unit geodesic
has speed <=1.02 and derivative bounds 16, 7000, 4e6 through order
four, the last essentially everywhere. At the center gamma''=0 and
|gamma'''|<=15. Taylor through the cubic curve jet therefore gives

    |q(gamma(r))-r^2(e.v)^2-r^4(e.v)(e.gamma'''(0))/3|
       <4e5 r^5.

The exact fifth absolute Bessel moment is
integral r^6 |B_M(r)|dr=225/(4096M^3)<.5/M^3. Including 2pi
angles and the angular sum, the local flow remainder costs <2e6C.
No fourth metric derivative is assumed across the seam.

The full geodesic model also has a finite long-flow contribution to
this test, including returns. Start it at 1/128, exactly where the
local estimate ends. For any a>=1/400 the radial tail moments obey

    |J0,tail|<=18C/a^4, |J2,tail|<=36C/a^2,
    |J4,tail|<=1152C E1(9a)<1e4C.

These follow from sum d|B_M(r)|<=72C/r^6 and the half-gap variant
1152C r^-6 e^-9r. The probe has supremum <=1/1024,
|X^2q(x)|<=2 and |X^4q(x)|<=120. Its full angular tail action is
<=2pi[(1/1024)|J0|+|J2|+5|J4|]<4e9C, inside the same reserve.
It is not identified with a compact cap.

The cap probe uses radius 1/32, whereas the localization below uses
the radius-1/400 second-moment anchor. Change anchors explicitly by
the annular second moment of the gain-two kernel. Its bound is at most
the angular constant times 2e8C log[(1/32)/(1/400)], below 1e11C.
Together with the local flow and full-flow tail, this stays inside
the 1e13C finite-contact allowance. The bounded gluing form remains
outside this singular finite-part anchor.

For a symmetric remainder bounded by KC/|s|^3 in one dimension,
fix the second moment by this compact quadratic and p0 by the constant
probe. With X=||v'||, Z=||v''||, splitting at b<=R gives

    |F[v]| <=KC[b XZ+b^2 Z^2/8+2log(R/b)X^2].

Indeed the near Taylor remainder has L1 bound
|s|^3 XZ+s^4Z^2/4; the farther difference and its subtracted square
each have bound s^2X^2. With b^2=min(R^2,epsilon/(4K)),

    |F[v]| <=epsilon C Z^2/3
       +KC[1/4+2log(R/b)]X^2.                          (9)

The two-dimensional version uses the same unitary geodesic transport
estimate and the angular Hessian identity. Angular and Jacobian factors
are included in K. This avoids a derivative assumption on a kernel
envelope. Below K=1e13 is sufficient for all such remainders.

## 7. An all-frequency bound and the auxiliary even collar

The accepted flat multiplier is

    Kflat(p)=-C F_eff p^4/(120pi)
      +sum d M^2/(2pi) integral_0^1 R(a p^2/M^2) d alpha,
    a=alpha(1-alpha), R(z)=log(1+z)-z+z^2/2.

To obtain a lower bound at ALL frequencies, put
nu_+/-=(q+/-sqrt(q^2-4t))/2 and
F(t)=-2[psi(1+nu_+)+psi(1+nu_-)].
Partial fractions identify this with the same common finite part,
F(0)=4gamma_E-2H_q, and F'(t)<0. For Re w>0,
[DLMF 5.9.13](https://dlmf.nist.gov/5.9.E13) writes psi(w) as log w
minus a Laplace integral with kernel in (0,1). Both 1+nu have real
part at least one and their product is 1+q+t. Therefore
F(t)<=-2log t+4.

For K(z)=Kflat(sqrt z), exact differentiation of the convergent remainder
gives

    K''(z)=-C/(2pi) integral_0^1 a^2
          [F(Caz)+Caz F'(Caz)+4log(mu0B)]d alpha.

Using integral a^2=1/30 and integral a^2 log a=-47/900, then
integrating twice with K(0)=K'(0)=0, proves

    Kflat(p)>=C/(30pi)p^4
                [log(p/(mu0 ell_o))-38/15].             (10)

No large-frequency asymptotic is substituted for this inequality.

Extend the actual even collar log warp h=log(a/a(0)) from |Y|<=.1
to zero at |Y|=.2 using the quintic 1-10s^3+15s^4-6s^5.
The Riccati equation and seed curvature bounds give
|h|<=Y^2, |h'|<=2|Y|, |h''|<=2 and one-sided |h'''|<=30.
The extension satisfies |h_ext'|<=1.2, |h_ext''|<=42,
|h_ext'''|<=3270. Thus

    H=-partial_Y^2+V, V=h_ext''/2+h_ext'^2/4,
    |V|<=22, Lip(V)<=2000,
    H=(partial_Y-h_ext'/2)^*(partial_Y-h_ext'/2)>=0.

This is an auxiliary COMPLETE EVEN line, conjugate to the weighted
radial Laplacian. It agrees with the actual collar on the later support.
Set v=sqrt(a_ext)T and w=Hv. Since H+24>=-partial_Y^2, the resolvent
integral for log gives

    <w,log sqrt(H)w> >= <w,log|D|w>-12<v,Hv>.

Use positive regularization and then w=Hv to remove it:
log(1+24/H)<=24/H has finite quadratic form on this range.
No spectral gap or Dirichlet substitution is assumed.
[Beckner, Corollary 2, equation (10)](https://arxiv.org/pdf/math/0612486)
then gives

    <w,log sqrt(H)w> >= integral -log|Y| |w|^2
         +D_H||w||^2-12<v,Hv>,
    D_H=log2+psi(1/4)=-gamma_E-2log2-pi/2>-4.             (11)

The full-line result applies to even functions, including their image
normalization. It does not require every actual geodesic to cross the seam.

The scalar model is matched quantitatively, not just by a symbol.
Its exact positive heat weight is

    beta(t)=sum d M^2/(2pi t)
          integral_0^1 exp[-t M^2/(alpha(1-alpha))]d alpha,
    Kflat(L)= -C F_eff L^2/(120pi)
       -integral beta(t)[e^(-tL)-1+tL-t^2L^2/2]dt.

Counting and the gap give
beta(t)<=.8C pi^-1 t^-3 e^-648t.
For H=D^2+V, the bounded-potential
[Feynman-Kac formula](https://arxiv.org/html/math-ph/9907022)
bounds the off-diagonal heat difference by
(e^(24t)-1)(4pi t)^-1/2 e^(-s^2/(4t)).
Integration against beta gives <13C/|s|^3. Together with (5)
this proves the required gain-two scalar comparison.

Its second contact is also fixed. For the scalar quadratic v(Y)=Y^2
at the observation center, Hv=-2 and H^2v=-4V0 there. After that
known fourth contact is removed, Feynman-Kac gives

    |(e^-tH-e^-tD^2)v(0)+2V0 t^2|
      <=4 Lip(V)t^(5/2)+484t^3 e^(22t).

Integrating each massive weight costs <100/M^3+3/M^4, whose tower
is <70C. A compact quadratic changes this by <1e4C: use the
Gaussian tail and beta(t)e^(22t)<=C/(4pi t^3), which follows from
integral a^2=1/30 and (1-22/1296)^-2<1.04.
Converting v to the base test u=a_ext^-1/2 v retains its explicit
cubic and quartic jets in the common contact polynomials. The remaining
cubic coefficient of a_ext^-1/2 is bounded using the displayed h jets;
after multiplication by v2 it has degree at least five. Proper-to-normal
coordinate composition has the same bounded jets. The normal-cap
1e4 remainder reserve therefore applies, and 1e13C safely covers
the finite second-contact conversion.

Both comparisons annihilate constant BASE T. For the scalar model,
H sqrt(a_ext)=0 and Kflat(0)=0; equivalently its complete weighted
heat semigroup preserves constants. Cutoffs of this bounded generalized
zero mode justify the identity. Half-density V'' and V^2 terms are
fixed by this identity, not additional independent p0 contacts.

The two complete models have different global completions. When
identifying their collars and extending the localized field by zero,
retain both escape terms. The same heat bound and the Gamma(5/2)
integral give |Kaux(s)|<=3C/(pi|s|^5), hence tail row norm
<=3C/(2pi R^4)<1.3e10C at R=1/400. Even allowing 2e11C for
density/image conventions, this fits the bounded zero-order envelope.
The geodesic tail is counted separately below; constant annihilation
alone is not used to erase either restricted-collar escape term.

## 8. The exact Weyl completion

The Euclidean W4 is minus B^2/(4pi) times the fixed bulk polynomial in
[the spherical anomaly action](spherical-anomaly-action.md).
In this section use physical canonical ghat derivatives and put

    omega=log(r/B), R=Rhat, U=Delta omega, V=|grad omega|^2,
    y=B^-2, k=y^2(1/15-q^2/3),
    d=Delta T, f=Delta(S+T/3), w=grad omega.grad S,
    j=|grad S|^2.

Its compact bulk Hessian is exactly

    W4''= B^2/(20pi) integral f^2
       -B^2/(30pi) integral(omega+1/6)d^2
       -B^2/(2pi) integral[d D+f L+J j+W T w-7w^2/45+Z],
    D=S(-R/15+y/9)+T[(6omega+1)R/45+11V/135]-11w/135,
    L=T[-R/15+2U/5+7V/45]-7w/45,
    J=R/30-7(U+V)/90-y/18,
    W=-2R/15+14(U+V)/45,
    Z=T^2[omega R^2/30+R(U+V)/15-U^2/5
             -7UV/45-7V^2/90+2k omega]+ST[-R^2/30+2k].

In particular U cancels from D. The ordinary 1/6 term remains.
With theta=S+T/3=-b/6-(A+2Rphysical)eta/3,

    theta'=-Rphysical(b/3+c)-L0 eta/3,
    L0=F_rr/2+4A Rphysical-2Rphysical^2, |L0|<=65lambda,
    ||theta||<8sqrt(C0), ||theta'||^2<=1e4lambda C0.

The seed bounds |Rhat|<=13lambda, |U|<=22lambda t^2,
V<=5lambda t^2, y<=lambda and omega<=t^2, t=r/B, give the
following integral reserves in units beta lambda E/B^2:

| Quantity | Upper bound |
| --- | ---: |
| integral D^2 | 2e4 |
| integral L^2 | 1.4e5 |
| integral abs(J) j | 1.5e5 |
| integral abs(W T w) | 1.3e5 |
| integral abs(Z without k) | 7.6e4 |

Use |dD|<=d^2/100+25D^2 and
f^2/10-|fL|>=f^2/20-5L^2. The w^2 term is favorable.
For absolute continuity in the norm of section 11 it also satisfies
integral_hat w^2<=5 beta lambda ||S'||^2/B^2, hence <=2.5e5 in
the same table units. This follows from w^2<=V|grad_hat S|^2.
The k terms cost <=[2000/(3pi)]Gq^2/(lambda B^4)E0<1e-24E0.
Consequently, now in normalized orbit units,

    W4'' >=-C/(30pi) integral(omega+19/60)(Delta T)^2
       +C/(40pi) J2[S+T/3]
       -[2e6G lambda+1e-24]E0.                         (12)

The actual finite-cut boundary representative is retained as in
[the boundary Weyl audit](stationary-exterior-boundary-weyl.md).
It contributes no additional compact-core Hessian here. The selected
mass/fast directions have both outer traces T(0)=T'(0)=0 (mass Y^3
and faster behavior included), and the decaying-end extension is taken
in the stronger finite model-plus-absolute-W4 norm. These traces and
that approximation control the logarithmic boundary flux; T(0)=0
alone would only eliminate the seam delta potential. The retained
boundary representative is not dropped by an H1 cutoff argument.

## 9. Localization, including its finite contacts

Choose a=1/100 and P(s)=I_s(6,6), explicitly
462s^6-1980s^7+3465s^8-3080s^9+1386s^10-252s^11.
Let theta=0 before Y=a, theta=pi/2 after Y=2a and
theta=(pi/2)P((Y-a)/a) between. Put chi_o=cos theta, chi_i=sin theta.
This is a C5 partition with bounded piecewise sixth derivative and
chi_o^2+chi_i^2=1. It is constant near the seam and bolt.

Bernstein differences bound P derivatives 1..6 by
(11,220,3960,63360,887040,10644480), with the sharper P'<3.
Useful theta derivative reserves are
a^-j(5,350,6300,1e5,1.4e6,1.7e7), and Bell composition gives
chi reserves a^-j(5,375,12000,650000,3.8e7,3e9).

No high seam jet is required. On the transition and its radius-1/400
neighborhood, the geodesic equation Y''=h'(Y)(1-Y'^2) and the
seed curvature derivatives give

    |h'|<1, (|h''|,|h'''|,|h4|,|h5|)<=(2,10,50,200),
    (|Y'|,...,|Y6|)<=(1,1,10,100,1000,10000).

For example the direct Y5 and Y6 upper sums are 564 and 8994.
Thus

    |X6chi| <=chi6+15chi5+245chi4+2115chi3+8500chi2+10000chi1
             <1e22,     |X4chi|<1e14.

Fix the second-contact anchor q_x(s)=s^2 xi(s), with proper radial
displacement s=Y-Y_x and R=a/4=1/400,
with xi one through R/2, zero by R, and the same beta-polynomial
transition. Its derivatives 0..6 are bounded by
(R^2,24R,970,37100/R,1.28e6/R^2,4e7/R^3,1.1e9/R^4).
Leibniz gives |X6(chi q_x)|<1e20 and |X4(chi q_x)|<1e12.
If the probe meets the transition it remains away from the seam.
If it does not, the partition is constant on its support and the IMS
operator applied to that probe is exactly zero.

The point action of P_M is <=2/(525pi M^4) sup|X6u|.
This includes the full 2pi angular integral. Its sum is
<2e-5 C sup|X6u|, and the restored Q4 point action is
<=C sup|X4u|. Hence the IMS operator Qgeo-sum chi_j Qgeo chi_j
has explicit zero and second-contact reserves EACH <1e22C.
The zero contact is -sum chi_j Qgeo(chi_j), since Qgeo(1)=0.
The fixed quadratic extracts the second contact using the same
finite-part anchor; no unspecified diagonal extension is inferred
from a kernel bound.

The separated IMS kernel gains two powers:
sum(chi_j(Y)-chi_j(Y'))^2<=25|Y-Y'|^2/a^2.
On r<=1/400 the principal envelope is <1e4C/r^6; even the rounded
2e8C/r^4 error adds only 1250C/r^6 there. The resulting IMS
kernel, radial comparison (5), and scalar heat comparison all fit
K=1e13 in (9). Beyond that radius, Liouville invariance bounds the
flow pushforward's row mass even past conjugate points. Its cosine
IMS tail is <=500C/R^4<2e13C. Polynomial contacts have already
been counted in the exact point actions.

The weighted fourth Weyl form must also be localized. With
w=omega+19/60, its exact radial IMS difference has integrand

    6w theta'^2 T'^2
     +[8w theta'theta''+4w h'theta'^2+2w'theta'^2]TT'
     +w[theta'^4+(theta''+h'theta')^2]T^2.

Here derivatives are in Y. All coefficients are supported in the
transition, where w<7, |w'|<=2/a and |h'|<=4a. Its absolute cost
is <=6e4 a^-2 J1+1e6 a^-4 J0,transition, before C/(30pi).
The unbounded outer logarithm is never differentiated at Y=0.

## 10. Positive fourth margins and the complete budget

For the outer field chi_o T the auxiliary metric agrees with the
actual one. The seed polynomial has calB>=x/2 on Y<=.02, so
u(Y)>=sqrt[x/(2C)]Y and

    omega<=-log Y-.5log(x/C)+.5log2,
    -log(mu0 ell_o)+.5log(x/C)=-log(mu0 ell_A)>9.

The last bound uses actual benchmark data. Write C_seed,0 for the
reference seed coefficient (distinct from the tangent norm C0) and
k_bench=2000. Then C/C_seed,0>1-1e-21,
C_seed,0>6k_bench^3 and |d|<1e-12 give x>k_bench^3/2=4e9.
Since e<3 and 3^20<4e9, .5log x>10; the fixed scale change is <1e-21.
Combining (10)-(12) leaves outer reserve

    [C/(30pi)] [9-38/15-19/60-4-(log2)/2] J2[chi_oT]
       >[C/(30pi)] (33/20) J2[chi_oT],                 (13)

apart from the explicit Hardy gradient cost 12C J1/(30pi).

On the interior support Y>=a, monotonicity gives
omega<=log(sqrt(2C/x)/a)<log350<6. Also
kappa4=H_q/2-gamma_E-log(mu0B)>13.9. For a rational check,
e<11/4 implies e^7<1200 and e^.5<5/3, hence log2000>7.5;
use H_q>=log(q+1), q=2000^4 and gamma_E<1.
The interior reserve from (2) and (12) exceeds
7.5 C J2[chi_iT]/(30pi), apart from the displayed curvature-gradient term.

Apply (9) with epsilon=1e-3, K=1e13, R=1/400. Its gradient
coefficient is <1e15C. The global actual/geodesic comparison, model
comparison and IMS remainders consume small portions of (13).
Use the angular Hessian identity for the global comparison and
||v''||^2<=2||Hv||^2+968||v||^2 for the auxiliary comparison.
There are two global geodesic gain-two forms (actual comparison and
model IMS) and one outer scalar comparison. The two global losses are
at most (2epsilon C/3)||X2T||_SM^2, with weak Bochner bound
||X2T||_SM^2<=(3pi/4)J2[T]+(pi/2)J1[T]. The scalar loss is at
most (2epsilon C/3)J2[chi_oT]+(968epsilon C/3)J0[chi_oT].
The exact unweighted IMS identity converts global J2 to the localized
sum with only the already bounded lower terms. The total fourth loss
per outer field is <=epsilon(pi/2+2/3)C<.00224C; its initial
margin 33C/(600pi)>.0175C therefore leaves, conservatively,
C/(100pi) sum J2[chi_jT].

There is an important physical weight in the cutoff gradients. On
Y<=2a, u<=Y/sqrt5, hence u^4<=16a^4/25. Equation (1) therefore
improves on this region to

    C J0,outer/E0 <5.2e-29(16a^4/25)<3.4e-37,
    C sum J1[chi_jT]/E0 <1.8e-27.

The extra gradient is exactly integral theta'^2T^2 and costs <8.5e-32
in the latter normalized units. Estimating it by the unweighted global
J0 bound would waste this margin.

For completely explicit aggregation, let J1* sum the gradient norms
of T, both localized fields and the auxiliary half-density v;
let J0* sum the corresponding zeroth norms with the unitary measures.
The half-density derivative costs at most twice the localized
Dirichlet norm plus .0008 J0,outer. Thus safely

    C J1*/E0 <1e-26,   C J0*/E0 <2e-28.

All finite second contacts, IMS, curvature and Hardy terms fit
1e23 C J1*. All bounded gluing/Schur/IMS zero remainders OTHER
than the actual constant contact fit 1e25 C J0*. These envelopes
include the corrected 6e24 gluing form and the finite model tails.
Their costs are <.001E0 and <.002E0. The actual p0 cost (8) is
<.0006E0. The weighted Weyl IMS cost is <1e-12E0.

The remaining action allocations are also retained. The canonical
LLL costs <1e-25E0 by the seed proof. For the reference potential,
fixed B,A,lambda_e give the exact disk Hessian

    Gamma_ref,Euclidean''=-16pi beta integral r^2 [
       (S+T)^2 pA +2S(S+T)pS+S^2pS+(r/2)S^2 partial_r pS]dr.

The minus is the Euclidean continuation of the Lorentzian reference
potential convention. The S^2pS term comes from the SECOND variation of X=e^(2S)r^2
and must not be omitted. The exact source profiles give
|pA|,|pS|<3341q^2/B^4 and |r partial_r pS|<1.4e4q^2/B^4.
With ||S||^2<=100C0, ||S+T||^2<=1100C0, their coefficient sum
is <7e6, so |Gamma_ref''|<1e-16E0. This is the spherical
bookkeeping functional, with no arbitrary nonspherical extension.
The low-order Weyl cost in (12), LLL and reference together are
<1e-12E0.

The complete displayed disk loss is <.004E0. Retain the comfortable
bound

    Gamma_rest,disk'' >= -.005 E0.                     (14)

Finally, the preceding seed proof bounds the FULL action's off-shell
gauge and exponential-chart contacts by <1e-6E0. I0 is stationary,
so its gauge/chart Hessians agree. Combining (14) with I0''>=E0
proves the stated actual full-action TT bound I_full''>=.99E0.

## 11. Domain, verification and remaining scope

Begin with smooth bolt-regular tangents compact at the outer end.
The bolt is an ordinary interior point. A logarithmic H1 cutoff that
removes a nonzero bolt value is not used in this fourth-order argument.
Use the completion of this core in the explicit stronger norm

    N[h]^2=E0[h]+sum d_n P_Mn[T]
       +C[J2[T]+J2[S+T/3]+integral omega(Delta T)^2].

Each displayed positive component is a closed spectral or weak
derivative norm, pulled back by the fixed return map; its zeroth
norm is bounded by E0. The completion therefore embeds in the
physical E0 space. Sections 4-10 make the actual Hessian continuous
in N, so its lower inequality passes to this completion. This does
not infer closability from semiboundedness of an abstract form, or
claim an unproved identification with a global logarithmic Sobolev
space. The boundary logarithms cancel at relative bound one.

The physical mass/fast directions belong to this completion by an
explicit TT-preserving approximation. Let L tend to infinity and set
c_L=chi(r/L)c+a_L psi(r/L), where chi truncates above 2L and psi is
supported in [2L,3L]. Choose a_L so
integral_(r_h)^infinity (r^2F_r-2rF)c_L dr=0, with the actual bolt
r_h=B(1+h). For c=O(r^-3) with the
selected derivative bounds, the numerator stays bounded and the
denominator is of order L^2, because r^2F_r-2rF=-2r+O(1).
Thus a_L=O(L^-2), and the exact TT equation makes b_L compact above
3L and unchanged below L. Both tensors are smooth at the bolt.

Below L the return difference is delta eta=a(r)J_L with J_L=O(L^-3).
Its scalar components satisfy delta T=-(a'-aR)J_L=O(L^-3/r),
where the leading lambda r terms cancel, and delta S=-aR J_L=O(L^-3 r).
In the outer transition delta T=O(L^-2), with Y derivatives through
order three of size O(L^(j-2)). Consequently

    ||X2 delta T||^2=O(L^-1), ||X3 delta T||^2=O(L).

To control the full model norm, use both elementary spectral bounds
g(z)<=z^4/24 and g(z)<=z^6/720. For its spectral variable p,

    P_M(p)<=min[p^4/(90pi^2 M^2),p^6/(525pi^2 M^4)].

Split the angular sum at M=|p| and use N(E)<=4CE^2. Its multiplier
is bounded by a fixed-seed constant times p^4(1+log_+|p|), and hence
by a constant times |p|^4+|p|^(4+2epsilon), 0<epsilon<1/2.
Spectral interpolation of the two displayed flow norms gives
O(L^(-1+2epsilon))->0. This establishes approximation without a
Sobolev-equivalence assertion. The largest absolute Weyl contribution
is bounded by
L^-6 integral_(1/L)^1 Y^-6(1+|logY|)dY=O(log L/L), also tending
to zero; E0 and the remaining terms converge as well.

Thus the existing mass/fast directions are included with both outer
traces T(0)=T'(0)=0. The first also kills the retained Delta Rbar
seam contact. All outer integrations and their boundary fluxes are
defined by this stronger-norm approximation, not by pointwise traces
alone or a punctured-bolt H1 argument.

All common finite parts were taken before removing the angular
regulator. Auxiliary caps serve different purposes: the source heat
caps fix the actual constant contact, normal caps estimate the
separated kernel, and the complete scalar collar proves (11).
Their actual-domain differences have been counted rather than
identifying these auxiliary functionals with the physical determinant.
No new state, temperature, finite coefficient or scale has been chosen.

The [focused checker](../checks/stationary_curved_quantum_tt.py) audits
the displayed finite algebra and rational budgets. Its immutable
[receipt](../evidence/stationary-curved-quantum-tt.json) and
[raw log](../evidence/stationary-curved-quantum-tt.log) record the actual
execution outcome and input hashes. Geodesic-flow construction,
Feynman-Kac, common finite-part passage, locality and form closure
remain ordinary proofs; a finite algebra run does not replace them.

This result supplies the missing curved TT-sector estimate at the
analytic seed. It does not establish positivity of the coupled
trace/TT system or infer a physical instability from its oscillatory
adjoint channel. The [actual weighted map](stationary-global-weighted-map.md)
and [adjoint phase theorem](stationary-adjoint-phase.md) retain their
existing scope. Physical kernel, phase and source-moment evaluation,
the compatible inverse and nonlinear transfer are the next tasks.

Remaining to-do list: determine the coupled physical kernel and adjoint
phase, evaluate source moments and quantify the compatible inverse;
transfer regularity and close nonlinear physical logarithmic phase
or justify controlled spatial reduction; match exterior/state and
charge, transfer S138, and complete driven stress, radius, controller,
capture, tail, spatial area and S074-S080.
