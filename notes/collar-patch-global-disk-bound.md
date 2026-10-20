# Joining the collar and patch graphs in the physical disk norm

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof
and verification contract before execution. **A numerical global
linear disk bound; nonlinear closure and the physical goal remain open.**

The [annulus/bolt graph](seam-safe-annulus-bolt-graph.md) and
[mass-extended outer graph](outer-mass-column-source-bound.md)
can now be joined. The missing steps are quantitative maps between
their different dimensions, coordinates, logarithmic norms and
time normalizations. This note proves both maps and absorbs their
two derivative interfaces. No derivative interface remains in the
result. Its very large constant is not a Newton contraction bound.

## 1. Norms, physical domain and theorem

Retain the action, state, clock, mass column, paired domain and
all finite/nonlocal allocations of the two predecessors. Fix

    epsilon=10^-160, r=epsilon/32, L=4r, t=10^-170,
    delta=7/2, x=log Y, k=sqrt(xi/C), xi=lambda B_phys^2,
    A0=beta sqrt(C lambda), gamma0=4pi beta B_phys^3.

B_phys denotes the seed radius B, not the interface norm below.
Keep the same maximal t/4096-separated chart cover of {Y>=r/4},
cutoffs chi_i, eta_i and partition zeta_0,zeta_i. Fix any such
cover once. All bounds below are uniform over that choice.

Write v=(theta,b), theta=S+T/3, b=2T, and w=(theta,T).
The quantified similarity is vhat=M_epsilon v. Its mass column is
ehat_m=c_L(Y)Y^3(-1/12,1), e_m=M_epsilon^-1 ehat_m. The unique
correlated coefficient is m=m_Y=-3k^3 m_c. Put

    v_f=v-m e_m, vhat_f=M_epsilon v_f,
    w_f=(theta_f,b_f/2),
    w_m=M_epsilon^-1 c_L Y^3(-1/12,1/2),
    chi=c_r, sigma=(1-c_(r/4))c_(8r).

The fast field, not the uncut mass cubic, belongs to X_(7/2).

Let P(w), B(w), N(w), F(w) be exactly the canonical patch,
interface, whole-half value and two-row source norms of the
annulus/bolt graph. In particular P uses the square sum of the
local sums H4(theta)+H4log(T). Define

    O=||chi vhat_f||_X, Z=||chi lambda^-2 M_epsilon f_dual||_Z,
    p=P(w)/sqrt(A0), b_I=B(w)/sqrt(A0),
    Q=F(w)/(C sqrt(A0)), V=V0(vhat), M_r=abs(m)r^-1/2,
    f_dual=(f_S,(f_T-f_S/3)/2), f=A(S,T).

The lowercase b_I is an interface norm, not the field b. Define
the following ORIGINAL-variable, mass-extended global input norm:

    G(S,T)^2=||chi v_f||_X^2+abs(m)^2
      +A0^-1 sum_i [||chi_i theta_f||_H4(R2)^2
                                  +||chi_i T_f||_H4log(R2)^2]. (1)

The outer X and every chart Fourier multiplier are the retained
ones, including log(2+|frequency|). The normalized time factor in
(1) is explicit. This is a fixed norm on the accepted mass-extended
disk space; it is not an arbitrary rescaling of the physical source.
Finite atlas/cutoff equivalence identifies its topology with the
accepted global space. No numerical equivalence with an unspecified
other atlas or a higher York norm is asserted.

The source/value/mass a priori result is

    G <=10^22 Z+10^5000 Q+10^386000 V+10^5100 M_r.       (2)

For the unique compatible linear solution in the accepted range
ker L0, the actual physical source norm is

    Fdisk(f)=||u^-2(f_S,f_T)||_L2(dmu4), u=B_phys/r_phys.

Then (2) and the retained value and mass estimates give

    G(A_m^-1 f) <=10^386020 lambda^-2 gamma0^-1/2 Fdisk(f). (3)

This is an upper bound for the complete actual GLOBAL disk
inverse in the specified norm, including the bolt and the fast
remainder of its mass datum. It is not the normal inverse alone.
The one signed compatibility condition is retained; (3) does
not assert surjectivity onto arbitrary sources or decide the
nonlinear preparation-source equation. The giant upper constant
does not establish a useful inverse/residual or Newton estimate.

## 2. A common collar cylinder and finite derivative reserves

Use y=a(0) tau, of period A0. On the collar the canonical metric
and measure are exactly

    gbar=dY^2+j(Y)^2 dy^2, dvol_gbar=j(Y)dY dy,
    1/2<j<2, m_0(Y)=sqrt(j(Y)), M_epsilon=m_0 for Y<=epsilon.

Only the band r/16<Y<32r is used for the coordinate transfers.
It lies inside the analytic one-sided collar, away from its seam
and the bolt. Local y coordinates may be unwrapped on each chart;
averages always use the original period A0, with no period change.

Here are sufficient derivative reserves through order five.
Use the annulus cover's native and transition bounds and its
partition constants D=10^12/t^2, H_part=2e18/t^2,
E_part=10^28/t^2. On relevant overlaps the maps between z_i
and (Y,y) have derivative bounds j!(10^4/t)^j, and Jacobian
and length conversion factors less than two. These follow either
from the explicit one-sided quadratic normal chart or its analytic
inverse with margin t/512; the real first derivatives are bounded
by the metric ratios. They do not use a polar frame at the bolt.

The factors eta_i zeta_0, chi_i c_(8r), and f(s)chi_i have
bounded component jets in the chart. In (Y,y) coordinates the
partition/transition product bound has base less than
3 E_part(10^4/t)<10^33/t^3. The beta a(s), as well as f(s),
has the same finite derivative reserve used in the predecessor.
Thus set the common enlarged base

    H=10^40/t^3=10^550.

The lift coefficients formed by multiplying these factors by
m_0^-1 Y^delta have jets bounded by 10^10 j! H^j, including
order zero. Their supports have Y<32r<1. This also follows
directly by Cauchy for m_0^-1 Y^delta on each small complex
chart, whose distance from Y=0 exceeds its radius. Choose the
positive real branch of Y^delta and its local analytic extension.

For the lift's composition x=log Y(z_i), its nonzero-order jets
are at most j!(10^6/t^2)^j<=j!H^j. One can check this without
a complex logarithm: combine the derivatives (j-1)!Y^-j of log
with the chart jets, using Y>r/16>t. Constants of x itself do
not enter derivatives of a composed input.

For the reverse map use D_x=Y partial_Y. Its repeated product
and chain rules give, at order j<=5, coefficient sums bounded by

    K_j=10^40 j!(10H)^j.                              (4)

This includes sigma, m_0, the constant two-field conversions and
the conormal derivatives of Y^-delta after removing that factor.
For clarity, a sufficient finite product budget before those
last constant reserves is

    sum_(k=0)^j binom(j,k) 10^10 (j-k)! H^(j-k) C_k,
    C_0=1,
    C_k=sum_(l=1)^k 2^l Bell_(k,l)(1!H,2!H^2,...).

It is below (4). The expansion
(Y partial_Y)^j=sum_(k=0)^j Stirling2(j,k)Y^k partial_Y^k,
the retained conormal cutoff jets <10^14, and the additional
factor two for (theta,T) to (theta,b) fit the spare 10^30.
The checker verifies the finite budgets, not just their exponents.

The number of eta supports through a point, as well as chi
supports, is <10^6. For eta this follows from packing centers
inside the ball of radius 1.1t/64 about that point, with disjoint
d/2 balls and d=t/4096. The enlarged ball stays inside a normal
chart. Use this overlap count, not the total 10^384 patch count.

## 3. Lifting radial data to the interface and mass patches

For a scalar weighted radial input U on the x-line define

    (T_I U)_i(z)=eta_i(z)zeta_0(z)m_0(Y(z))^-1
                                      Y(z)^delta U(log Y(z)),
    (T_P U)_i(z)=chi_i(z)c_(8r)(Y(z))m_0(Y(z))^-1
                                      Y(z)^delta U(log Y(z)). (5)

Only patches meeting their collar supports contribute. Extend
each output by zero in its own R2 chart. The vanishing cutoff
jets make these extensions bounded through H5. Apply the
second map only on the cover's patches; their chi supports all
have Y>r/8. The first map already has r/8<Y<r/4 support.

For integer s=3 or 5, product/chain rule and (4) bound every
output derivative by a sum of radial derivatives U^(l). Summing
their squared integrals over patches and using overlap gives

    sum_i integral |U^(l)(log Y)|^2 dvol_gbar
      <=2*10^6 A0 integral |U^(l)(x)|^2 Y dx
      <=2*10^6 A0 ||U^(l)||_2^2.

The first expression is restricted to the supports in (5).
Euclidean/canonical measure conversion costs less than two.
The coefficient sum in the two-dimensional Hs polynomial is
3^s. Thus the finite endpoint budget

    2*sqrt(10^6)*3^s*max_(j<=s) K_j <10^4000

proves both endpoint operator norms <10^4000 sqrt(A0).

The dyadic transfer in the compact-cover proof is independent of
the dimensions of the input and output Hilbert spaces. Indeed
take orthogonal frequency shells in R on input and in each R2
output, using their orthogonal direct sum. The same block bound
32A 2^-|j-k|(1+|j-k|) has Schur sum 224A. Its norm-equivalence
factors give 28672A for H4log and 1536A for H4. Hence, including
the pair and local-sum/Hilbert conversions, (5) costs at most
10^5000 sqrt(A0) in the required mixed norms.
There is no H4-to-H4log inference or unquantified trace theorem.

For the fast field, chi=1 on the interface, so T_I applied to
Y^-delta chi vhat_f gives eta_i zeta_0 w_f, after b/2 conversion.
To handle the mass in both (5), use the legitimate fast test profile

    U_m(Y)=Y^-delta (1-c_(r/32)(Y))c_L(Y)Y^3(-1/12,1).

It is supported in [r/32,8r], and its X norm is <10^40 r^-1/2.
To see this, scale Y=r z: the weighted amplitude is below
6 r^-1/2, the x-support length is log256<6, and each cutoff's
log derivatives through five are <10^14. Product rule gives
derivatives <10^33 r^-1/2, and H4log<=2H5 supplies the stated
reserve. The inner cutoff is identically one on every patch
where the mass can contribute, since those patches have Y>r/8.
The extra c_(8r) in T_P is one on its mass support. Thus

    b_I <=10^5000 O+10^5040 M_r,
    P(w_m)/sqrt(A0)<=10^5040 r^-1/2.                  (6)

This uses an auxiliary cutoff profile in X, not the uncut cubic.

## 4. Averaging patch data back to the outer transition

For arbitrary chart inputs h_i, define the radial operator

    R(h)(x)=Y^-delta sigma(Y)m_0(Y) A0^-1
      integral_0^A0 sum_i f(s)chi_i h_i(z_i(Y,y)) dy. (7)

For the second component multiply this expression by two.
Every summand is zero outside the corresponding chart; chi_i
provides the required vanishing jets. All active points satisfy
r/4<=Y<=16r. On that band zeta_0=0 and sum f(s)chi_i^2=1.
Therefore, for the actual stationary field,

    R((chi_i w)_i)=Y^-delta sigma vhat.               (8)

This is an average of scalar values, not a point trace. Its
normalization is the original stationary time period A0.

Differentiate (7) in x through order s=3 or 5. For each derivative,
(4), Cauchy-Schwarz in y, and the pointwise overlap bound give
an L2(dx) bound by

    4*sqrt(10^6) A0^-1/2 (r/4)^-4
                   3^s max_(j<=s)K_j (sum_i ||h_i||_Hs^2)^1/2.

The power four is exactly delta+1/2: Y^-delta from the weight
and Y^-1/2 from dx=dY/Y. Coordinate Jacobians and field factors
fit the four. The displayed finite endpoint constant, without
A0^-1/2, is <10^4000. Apply the same dimension-independent
dyadic transfer, now from a direct sum of R2 spaces to R. It gives

    ||sigma vhat||_X<=10^5000 P(w)/sqrt(A0)=10^5000 p. (9)

Equations (6),(9) account for both normalized time factors.
No assumption that A0 is one or bounded below is needed.

## 5. Two-sided graph absorption and the mass remainder

The annulus/bolt and outer estimates, with N(w)/sqrt(A0)<=3V,
give

    p<=10^8 Q+10^-18000 b_I+3*10^381000 V,
    O<=2e20 Z+10^-40 ||sigma vhat||_X
                         +10^2000 r^-4 V+10^2005 abs(m)L^-1/2.

Insert (6),(9), and note L=4r and r^-4<10^647:

    p<=10^8 Q+10^-13000 O+10^-12960 M_r+3*10^381000 V,
    O<=2e20 Z+10^4960 p+10^2647 V+(10^2005/2)M_r.     (10)

The feedback product is 10^4960*10^-13000=10^-8040<1/2.
Solving these two scalar inequalities removes BOTH derivative
interfaces. The norm in (1) satisfies

    G<=1.01 O+p+10^5040 M_r+abs(m),                  (11)

by the mass-patch bound in (6) and the quantified outer inverse
multiplier M_epsilon^-1. Since r<1, abs(m)<=M_r. Substitution
of (10) into (11), retaining every source, value and mass term,
gives (2). The checker verifies each coefficient after exact
division by 1-10^-8040. Neither normal-inverse smallness alone
nor support overlap alone would justify this absorption.

## 6. Exact physical forcing conversion and the linear solution

Set D=lambda^-2 gamma0^-1/2 Fdisk. Source conversion is sharper
if the original u^-2 weight is retained before estimating. The
exact canonical source relation is

    (L_theta,L_T)=(4pi B_phys^4/C)u^-4(f_S,f_T-f_S/3),
    dmu4=(4pi B_phys^4/C)u^-4 dvol_gbar.

Hence the u weights cancel exactly in the source norms, and

    F(w)<=(10sqrt(pi)/3) B_phys^2 C^-1/2 Fdisk,
    Q<=(20pi/3)k^(7/2)D<3D.                         (12)

The first bound uses ||f_S||+||f_T-f_S/3||<=5/3||f||
in the common weighted Hilbert measure. Also

    gamma0=4pi A0 B_phys^4/(kC).

On supp chi, M_epsilon^2=j, u=kY w_seed and
|w_seed-1|<.001. The dual (theta,b) matrix has squared Frobenius
norm 23/18. Direct change of variables therefore gives

    Z^2 <=(23/18)k^7 sup(w_seed^8)D^2<9D^2.          (13)

This keeps both the half-density and stationary time measure;
the physical u^-2 source weight has not been replaced by a
different source space. The loose annulus forcing conversion
remains valid but is unnecessary for (12).

The accepted compatible-solution estimates from
[the York source conversion](weighted-york-source-conversion.md)
and [the mass bound](tt-mass-numerical-bound.md) are

    ||(S,T)||_physical<=5.58e9 g lambda^-2 Fdisk,
    abs(m_c)<=5.40144e14 g D, 0<g<2000^-10<1.

Using M_epsilon<2, the two field conversions and u<1 yields
V<=24 gamma0^-1/2 ||(S,T)||_physical<2e11 g D.
Also 3k^3<1/4, so abs(m)<1.4e14 g D. Thus r^-1/2<10^82
and (2), (12),(13) prove (3), with ample explicit reserve.
This uses value and mass bounds already established independently
of the new higher graph estimate; there is no circular absorption.

The [accepted range and injectivity theorem](stationary-scalar-outgoing-inverse.md)
gives the unique inverse on ker L0. Apply the a priori estimate first on
the mass-extended smooth core and pass to its graph closure.
Boundedness of the two transfer maps just proved justifies these
limits. The mass column is included through its actual accepted
domain and full image, not through a false fast-space approximation.

## 7. Verification and remaining scope

Freeze this proof and checker before the scientific run. Pin the
annulus, mass-column, multiplier, value and mass receipts and all
their current inputs, including the accepted range/domain source.
Check the finite chain/partition/probe-profile ledgers, both
endpoint transfer budgets, dyadic Schur constant, period and
measure identities, exact two-interface coefficient absorption,
source transpose, physical-source cancellation, and final value/
mass substitution. Preserve every run and failed source snapshot.
Finite arithmetic and symbolic checks support the ordinary
continuum operator and domain arguments; they do not formalize them.

This supplies a numerical full disk inverse in one fully specified
mixed norm. Higher York extraction/reconstruction and the actual
residual/inverse composition still need quantitative estimates.
In particular no useful Newton bound follows merely because the
large number in (3) is finite. The fixed-density nonlinear map,
signed compatibility equation, admissible family and all physical
completion requirements remain unresolved.

Remaining to-do list: higher reconstruction and matched inverse/
residual composition; nonlinear family and signed compatibility;
exterior/state, charge and return; opening, stability, transport,
preparation and separate creation requirements.
