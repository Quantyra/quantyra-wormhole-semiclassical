# Common-chart transfer into outer proper coordinates

2026-10-20, E003/S151 supporting S139/S074. Ordinary argument and
verification protocol fixed before the acceptance run.

This resolves the local coordinate portion of
[S150's joint contract](averaged-ward-matching.md) on a stronger
common-chart space. It does not prove preservation by the unknown
global source/equation map. In particular it does not promote the
previously constructed C5 endpoint to C7, supply a cold estimate, or
prove a bounded inverse of the full equations.

The regularity issue is consistent with the general composition
problem studied by [Inci, Kappeler and Topalov](https://arxiv.org/abs/1107.0488v2).
Their paper concerns Sobolev diffeomorphisms; it is background, not
a citation for the weighted estimates below. Those are derived here.
No new quantum-state estimate is imported from that paper.

## 1. Exact map and its parameter jets

Fix the seed inward proper coordinate y=Y0(chi) once. On 0<y<=r,
r=.11, the common-disk conformal deformation reads

    g_b=exp(2b(y))[dy^2+a0(y)^2 d tau^2],
    F_b(y)=integral_0^y exp(b(z)) dz, Y=F_b(y),
    y=y_b(Y), a_b(Y)=exp(b(y))a0(y),
    beta_b(Y)=L0(Y)-L0(y)-b(y), L0=log a0.          (1)

The same physical clock tau is used throughout. We compare on the
fixed target interval 0<Y<=.1; the larger input collar leaves room
for inversion. There is no endpoint movement at Y=0.

Let b=b_*+s h+t k be affine in the common chart. At its current
real background define D=y partial_y, Jf=y^(-1) integral_0^y f,

    q=log(F_b/y), T=exp(q-b),
    A_h=exp(-b) J(exp(b)h),
    A_hk=exp(-b) J(exp(b)hk),
    B_hk=-A_hk+h A_k+k A_h-(Db)A_h A_k.           (2)

All right sides in the following identities are evaluated at
y_b(Y):

    y_h=-y A_h, y_hk=y B_hk, D_Y=T D,
    beta_h=-h+(Db+DL0)A_h,
    beta_hk=Dh A_k+Dk A_h
       -(D^2b-Db+D^2L0-DL0)A_h A_k
       -(Db+DL0)B_hk.                             (3)

Here A_hk means the integral of the product hk, NOT the mixed
parameter derivative of A_h. Differentiate F_b(y_b)=Y twice to
obtain y_hk: the terms are -exp(-b) integral exp(b)hk,
-h y_k-k y_h-b' y_h y_k. This checks both the sign and the
acceleration. Differentiating (1) gives the remaining identities.

## 2. A sufficient weighted space and explicit finite constants

Write ell(y)=1-log y, w3=y^3, w4=y^4 ell(y). Define the Euler jet
norm on the input collar by

    |f|_{j,w}=sum_{i=0}^j sup |D^i f|/(i! w).
    b=m y^3+v, N7(b)=21|m|+|v|_{7,w4}.            (4)

Functions have seven continuous derivatives on the open collar and
the indicated finite weighted norms. The cubic coefficient is part
of the data. This is stronger than the earlier ordinary weighted C5
contract. Since sum_{i=0}^7 3^i/i!<21 and r ell(r)<.36,
|b|_{7,w3}<=N7(b). Assume N7(b)<=.001, and write H=N7(h),
K=N7(k). Norms below are for real coefficients, without complex
spatial evaluation of b, h or k.

The norm is a product algebra with multiplied weights. D costs
at most j when mapping j jets to j-1 jets. J commutes with D;
on w3 and w6 its norms are 1/4 and 1/7, respectively. In
particular, elementary exponential and logarithmic series give

    |exp(+/-b)|_{7,1}<1.001,
    |q|_{7,w3}<.252 N7(b), |T|_{7,1}<1.001,
    |A_h|_{6,w3}<.251 H,
    |A_hk|_{5,w6}<.144 HK,
    |B_hk|_{5,w6}<.65 HK.                         (5)

For example |exp(b)-1|_{7,w3}<1.001 N7(b), so
|F/y-1|_{7,w3}<.251 N7(b). For the last bound the two
linear product terms cost .502 and the Db product costs at
most 7*.001*.251^2*r^3; their sum with .144 is below .65.
The numerical constants are uniform on the full input collar.

Let M>=1 bound both

    |DL0|_{6,y^2}, |D^2L0-DL0|_{5,y^2}.

Let P bound |psi-d^3Y^3|_{5,w4}, where psi is S146's fixed
analytic mass mode and .4<d=sqrt(x/C)<.42. These are FIXED
seed constants, not unknown geometry norms. Explicit bounds for
them are given below. Composing a j-jet with y_b costs at most

    (1+1.001)^j times the weight-ratio bound.

To see this, iterate D_Y=T D: applying D from n to n-1 jets
costs n and multiplication costs |T|. Sum the resulting
binomial majorants. Because |q|<.252*.001*r^3<4e-7,
the ratios of w3, w4, w6 and y^2 at y and Y are below 1.01.
Thus composition of five jets costs less than 33.
The image of r is greater than .1.

Subtract the cubic BEFORE using w4. Formula (3) becomes

    beta_h+m_h Y^3 =
      [-v_h+m_h(F^3-y^3)+(Db+DL0)A_h](y_b(Y)).     (6)

Here |F^3-y^3|_{5,y^6}<16 N7(b): the exponential costs
less than .757 and multiplication by y^3 costs less than 21.
Since |m_h|<=H/21, this term costs less than .762 N7(b)H.
Products y^5, y^6 and higher are bounded by w4 on this collar.
Equations (5)-(6)
therefore bound the right side in Euler w4 norm by
50(1+M)H. Equation (3) bounds the mixed response by

    33[3.514+49*.001*.251^2+7*.001*.65
                       +(.251^2+.65)M] HK
       <150(1+M)HK.                              (7)

This deliberately keeps a loose finite reserve; it is not an
optimized norm for proving contraction.

The ordinary weighted C5 seminorm used in S146 is bounded by
this Euler five-jet norm: y^j partial_y^j is
D(D-1)...(D-j+1), its absolute coefficient sum is j!, and
j!<=100^j for j<=5. Hence
|partial_y^j f|/(j!100^j y^(4-j)ell)<=|f|_{5,w4}.

The exact output mass parameter is p=-m/d^3. Its remaining
mass-basis correction is (m/d^3)(psi-d^3Y^3).
Since 21 d^3>1, define

    C1=51+50M+P, C2=150(1+M).

With N_out^loc denoting S146's mass-plus-remainder norm restricted
to Y<=.1, we have

    N_out^loc(beta)<=C1 N7(b),
    N_out^loc(beta_h)<=C1 H,
    N_out^loc(beta_hk)<=C2 HK.                     (8)

The first follows by integrating the first response on the radial
segment from zero to b. For an affine common-chart family the
mixed cubic vanishes exactly. It must not be counted as a new
mass datum.

These are Frechet C2 maps into the weighted output space, not
just formal pointwise derivatives. On [delta,r] the usual chain
and inverse differentiation is continuous in the indicated norms.
At y=0 every transported term in beta_hk has weight at least
six, and the terms containing the seventh background derivative
have weight at least eight (the b part has weight nine).
After division by w4 they tend uniformly to zero on bounded
input balls. This controls the possibly nonuniform continuity
of the highest weighted derivative near zero. The untransported
linear term has two spare derivatives. Splitting at delta and
then sending delta to zero proves continuity and the difference
quotient limits. No weighted C8 premise is used.

For clarity the fixed constants can be enclosed without fitting
sampled data. The retained seed collar has |calB|>x/2,
|calB'|<5C, C/x<6.1 on its .01 complex strip, and
L0'=calB'/(2 sqrt(C calB)); thus |L0'|<9.
On the half strip, |L0^(l)|<=9(l-1)!200^(l-1).
Use L0'(0)=0 to bound |L0'|/y<=1800. For l>=2 put

    R_l=9(l-1)!200^(l-1) r^(l-2), R_1=1800,
    E_j=sum_{l=1}^j Stirling2(j,l) R_l,
    M=max(1, sum_{j=0}^6 E_(j+1)/j!,
             sum_{j=0}^5 [E_(j+2)+E_(j+1)]/j!).   (9)

S146's |psi|<10 on the .01 strip and its zero jets through
order two with cubic d^3 imply, for f=psi-d^3Y^3,
|f^(l)|<=10*4!*100^4 Y^(4-l)/(4-l)! for l<=4;
for l=5 use 10*5!*100^5. With

    V_l=10*4!*100^4/(4-l)! (0<=l<=4),
    V_5=10*5!*100^5 r,
    P=V_0+sum_{j=1}^5 sum_{l=1}^j
                         Stirling2(j,l)V_l/j!,    (10)

we obtain an explicit rational enclosure. Using ell>=1 makes it
conservative. The checker evaluates (9)-(10); these enclosures
are not optimal and do not certify the existing endpoint.

## 3. Radius jets and the physical correlation

Let the same common-chart radius be

    u_disk(y)=u0(y)exp(n(y)), n=n3 y^3+z,
    N7(n)=21|n3|+|z|_{7,w4}.

Then u_b(Y)=u_disk(y_b(Y)). Write rho=log(u_b/u0(Y)),
R0=log u0, and allow n to be affine in the same parameters.
The logarithm R0 contains log y, but DR0 and D^2R0 are regular.

    rho_h=n_h-(Dn+DR0)A_h,
    rho_hk=-Dn_h A_k-Dn_k A_h
       +(D^2n-Dn+D^2R0-DR0)A_h A_k
       +(Dn+DR0)B_hk.                            (11)

These identities also retain the acceleration. Suppose the joint
background N7(b)+N7(n)<=.001 and use joint direction sizes
H=N7(h)+N7(n_h), K=N7(k)+N7(n_k).
Take M_R=500, bounding |DR0|_{2,1} and
|D^2R0-DR0|_{1,1}. To verify this from the retained seed,
u0'=sqrt(calB/C), u0''=calB'/(2C), and |u0'''|<=250
by Cauchy on the .01 strip. On the real collar,
calB/C<=d^2+5u0<.42^2+1.25<1.2^2 and
u0'>d/sqrt(2)>.28. Hence u0>=.28y. Put
A=y u0'/u0, B=y^2 u0''/u0, E=y^3 u0'''/u0:
|A|<5, |B|<1, |E|<11. The exact identities

    DR0=A, D^2R0=A+B-A^2,
    D^3R0=A+3B+E-3A^2-3AB+2A^3

give bounds 5,31,359. Therefore
|DR0|_{2,1}<5+31+359/2<500 and
|D^2R0-DR0|_{1,1}<26+359+31<500.

The unweighted consequences of (5), r^3<.002 and r^6<2e-6
give, after the one-jet composition,

    |rho_h|_{1,1}<=M_R H,
    |rho_hk|_{1,1}<=M_R HK.                        (12)

For example the first precomposition norm is bounded by
.002H+2M_R*.001H. The mixed one is below
.000012HK+.000002M_R HK+.000004M_R HK.
Composition for one jet costs less than three, so the reserves
in (12) are ample. The two-jet norms of the factors used here
are controlled by the seven input jets.

Under S147's retained real radius guards D_Y log u>=1/4,

    u_h/u=rho_h,
    u_h'/u'=rho_h+(D_Y rho_h)/(D_Y log u),
    u_hk/u=rho_hk+rho_h rho_k,
    u_hk'/u'=rho_hk+rho_h rho_k
       +(D_Y rho_hk+rho_h D_Y rho_k
                       +rho_k D_Y rho_h)/(D_Y log u).

Thus first relative value/derivative norms are <=5 M_R H,
and mixed norms <=14 M_R^2 HK. These are coordinate estimates;
the real low guards must also hold for the chosen family.

More importantly, the finite boundary coefficients are fixed
algebraically. If the seed has

    a0=1+Y^2/(2C)+c0 Y^3+...,
    u0=dY+dY^3/(6C)+d c0 Y^4/4+...,

then (1) gives

    F=y+m y^4/4+O(y^5 ell), y=Y-mY^4/4+...,
    c=c0+m,
    e4=d c0/4+d(n3-m/4),
    [Y^3]alpha=2m-4n3.                           (13)

Consequently the zero cubic lapse condition is precisely
n3=m/2. Then e4=dc/4 and p=delta mass/(Bx)=-m/d^3,
as in S147. The proper radius-log cubic is m/4; it is the
COMMON-CHART radius-log cubic that equals m/2. Confusing the
two leaves an unwanted cubic lapse. All these coefficients are
exact in the small family, not first-order-only correlations.

Let Q=max(C1,sqrt(C2),5M_R,1), which the enclosures above
allow us to take as C1. With U=QH and V=QK,
(8) and (12) imply the local first/mixed coordinate and radius
jet portions of S150. A background N7(b)<=1e-7/C1 meets the
local S146 size bound. This statement does not establish the
cap extension, actual global operator/domain coefficients,
interior raw norm, source differentiability in the same family,
cold weights or preservation of the input C7 ball.

## 4. Why the earlier C5 estimate cannot be substituted

At the seed with a0=1 and b=0, (3) reduces to

    D^2 beta[ h,k ]=h' integral_0^Y k
                              +k' integral_0^Y h. (14)

Choose smooth cutoffs supported strictly inside the outer collar,
with k=1 on a smaller interval and its primitive nonzero there.
Let h_N=N^-5 sin(NY) on that smaller interval, with a fixed
cutoff outside it. These are legitimate zero-cubic weighted
directions with uniformly bounded C5 norm, but on the interval
the fifth derivative of (14) has leading term
N cos(NY) integral_0^Y k (up to an immaterial phase).
The remaining derivatives are uniformly bounded. Hence the
bilinear map is unbounded C5 x C5 -> C5. For an analytic
nonconstant seed the extra L0 terms contain at most fifth
derivatives of h_N after the primitive factors are used;
they do not cancel this leading sixth derivative.

This is a coordinate regularity obstruction, not a quantum
state or physical nonexistence theorem. Seven background
derivatives are sufficient for the finite-background mixed
estimate; we do not claim that C7 is the minimal possible space
or that the existing C5 source estimates provide C7.

The next global proof must either preserve a stronger space
with its own quantitative inverse, or use a formulation that
avoids this moving-coordinate differentiation. Merely rescaling
the old C5 ball cannot make an unbounded bilinear map bounded.
The stronger local coordinate theorem does not remove S149's
separate source-gradient limitation.

## 5. Focused verification protocol

Before execution: exact symbolic chain-rule/inverse identities,
cubic/radius/lapse series and algebraic bounds (5)-(10).
The test geometry is manufactured smooth real data, not a
solution or quantum state. Use exp(b) with cubic/quartic smooth
b, and independent n and seed polynomials. Compare (3) and
(11) against independently solved coordinate inversions at
Y=.025,.05,.09, using 32/64-point Gauss integration, 55 digits,
and centered mixed steps .001,.0005. Require absolute mixed
error <=2e-6(1+|target|), refinement ratio <=.35 unless both
errors <1e-30; quadrature agreement <=1e-11(1+|target|).
Check exact omission differences for the inverse acceleration
and radius-cubic correlation. Test the unbounded sequence in
(14) analytically and at increasing integer frequencies.

Preserve all runs, failed source snapshots, input hashes before
and after and raw output. Verification supports the ordinary
proof and the limited coordinate scope. It is not an acceptance
test for a full semiclassical wormhole.

Remaining to-do list: prove a usable global common-family/source
contract and cold response; complete remaining source sectors,
global inverse and exterior/state solve; complete charge, signed
source, opening, stability, transport, preparation and creation.
