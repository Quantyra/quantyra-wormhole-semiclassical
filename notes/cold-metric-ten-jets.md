# Cold finite-comparison metric derivatives through ten

2026-10-21, E003/S167. Ordinary argument before supporting verification.
For the ORIGINAL finite comparison

    sigma_ell=sigma0+ell u, t_ell=1+(t-1)exp(ell V),
    nu=tau V, tau=(t-1)/t, w_ell=log(t_ell/t),

the full ordered physical orbit tensors satisfy

    |nabla^j(ell u)|, |nabla^j w_ell|
       <ell 1e24 j!100^j a(s), 0<=j<=10,
    a(s)=s[1+log(.01/s)], h<=s=t-1<=.001.            (1)

At ell=0 both fields vanish. Derivatives use the seed's Euclidean
orbit metric in A units. The cold weight is retained at the regular
bolt, where s=h>0. No inverse power of h appears in (1). Keep the
original smaller C5 and C0/C1/C2 contracts. This is a bound for the
explicit finite comparison, not a nonlinear solution or stability
theorem. Outer order ten, endpoint maps and finite-metric H6 remain open.

The sources and algebraic reconstruction are those of the
[middle metric proof](middle-metric-ten-jets.md). The new input is the
[complete actual cold source C8](seed-cold-actual-domain-eight-jets.md).
The [original endpoint proof](stationary-completed-metric-response-norm.md)
and [repaired dependency review](repaired-source-metric-transfer.md)
supply the smaller low radial jets. The qualitative smoothness of
the original fields at the disk center is already established by the
[finite comparison theorem](finite-comparison-conormal-regularity.md).

## 1. Product subtraction before Ward reconstruction

Use outward seed proper radius rho from the bolt and the seed metric

    g=d rho^2+N(rho)^2 dtheta^2, N(0)=0, N'(0)=1.

With z=1/t, b=calB/C, p=sqrt(b), beta=N'/N and e=d/C,

    z'=-p, p'=-b_z/2, N''=(b_zz/2)N,
    beta'=b_zz/2-beta^2, L=p/z,
    S=-b_z/z+2b/z^2,
    c_nu=z b_zzz/6+4(x/C)/z^2+2e z^3-4e z^4.        (2)

The relation N proportional to p has a nonzero fixed proportionality
constant at the prepared positive temperature. We use N's regular
normalization, never divide its initial value by a bound for p(0).

Let theta0 be the full canonical homogeneous product trace divided by
q^2, including the LLL. Set T=(Theta_can-Theta_product)/q^2 and
Q=P-theta0/2, where P=B^2 p_can/(q^2 z^2). Subtract the same product
from the explicit local Weyl-minus-reference allocations:

    E_rho=B^2 rho_loc/(q^2 z^2)-theta0/2,
    E_delta=B^2(p-rho)_loc/(q^2 z^2)+theta0,
    E_A=4pi B^4 t^2 Arest/q^2,
    E_u=E_A-3z^2 E_delta, E_nu=z^2 E_delta.

All three first-line local differences and E_A vanish at z=1,
delta/C=0. This is the exact homogeneous stationarity of Gamma_rest
in the [reference-action split](stationary-centered-covariant-seed.md),
with canonical pressure theta0/2. It is not a subtraction of a physical
nonconstant source or a change of finite renormalization.

The LLL trace difference is -C(b_zz-2)/(24pi q), an explicit local
function. Adding it to the massive cold theorem gives, for j<=8,

    |T^(j)|/(j!100^j)<=K a(s), K=2e44.              (3)

These are proper radial derivatives: radial curves are unit-speed
geodesics, so the full source tensor bound controls them directly.
The weight a(s) is increasing on h<=s<=.001 and s increases with rho.

## 2. Quantitative seed coefficients, including the bolt

The centered polynomial, with tau=1-z and delta_bar=delta/C, is

    b=tau^2-B3 tau^3+B4 tau^4-4e tau^5+e tau^6
                            -delta_bar(1-tau)^3,
    |B3|<1.334, |B4|<.501, |e|<1e-24,
    0<delta_bar<1.001h^2, h<1e-30.                 (4)

At any real cold point, tau<=s and p<=s. On a complex proper disk of
radius .005, the equations tau'=p and p'=-b_z/2 bootstrap
|tau|+|p|<4s: while |tau|<=.005, |p'|<=2|tau|+2delta_bar, and
Gronwall gives exp(.01)(2s+.01delta_bar)<2.1s. Thus |z-1|<.004.
The polynomial bounds in the [local cold proof](seed-cold-local-eight-jets.md)
give |b_z|<.011, |b_zz|<2.05, |b_zzz|<8.1, |b_zzzz|<13.
The following larger coefficient reserves therefore have kth proper
derivative at most the displayed constant times k!200^k:

    |L|<1, |pz|<1, |z^2|<2, |b_zz|<12,
    |z b_zzz/2|<10, |2S|<2, |c_nu|<10.             (5)

The local formulas used in the middle proof are holomorphic jointly
on |tau|<=.01, |delta_bar|<=1e-4. Their absolute values, in the
same q-normalizations and before product subtraction, are below 1e16;
the Laurent/log coefficient estimates there remain valid on this
smaller radius range and enlarged thermal coefficient. Subtracting
the value at (0,0), Cauchy gives difference bounds

    2e18|tau|+2e20|delta_bar|.

Use |tau|<4s and delta_bar<=1.001s^2, s<=.001, to bound these by
1e22 s on every .005 proper disk. Hence the normalized kth derivatives
of E_rho,E_delta,E_A are bounded by 1e22 2^k a(s), and those of
E_u,E_nu by 1e24 2^k a(s). This also bounds the explicit LLL difference
well inside the extra allowance in (3). The analytic continuations
are of local seed coefficients only; T is differentiated on the real
geometry with its actual-domain C8 theorem.

For rho>=.5, 0<b_zz/2<1.1 on the real cold patch. Riccati comparison
with sqrt(1.1)coth(sqrt(1.1)rho) gives beta<1.1+1/rho<4.
On the .005 complex disk beta'=b_zz/2-beta^2 bootstraps |beta|<5,
since the displacement while |beta|<=5 is at most .005(6+25)<1.
This supplies the regular away-from-center coefficient bound.

For 0<=rho<=1 a separate, uniform analytic disk calculation removes
the polar singularity. On |rho|<=2, (4) gives

    |tau''-tau|<=10|tau|^2+2delta_bar.

Starting with tau(0)=h/(1+h), tau'(0)=0, a bootstrap |tau|,|p|<=100h
and the sinh/cosh variation formula sharpen both to 5h. Consequently
|b_zz/2-1|<51h. Write N=rho n(rho^2), n(0)=1. The regular Volterra
equation for N''=(b_zz/2)N gives

    |n(rho^2)-sinh(rho)/rho|<1e4h on |rho|<=2.     (6)

For example bootstrap |N/rho|<10; the Green kernel for d^2/drho^2-1
bounds the difference by 51h*10*4<1e4h and closes the bootstrap.
Evenness removes the apparent origin quotient. Cauchy in w=rho^2,
on |w|<=4, bounds the first two w-derivatives of this difference on
|w|<=2.25 by 1e-24. The same bound applies after dividing the
corresponding squared difference by w, using its zero value at w=0
and the maximum principle before Cauchy differentiation.

The positive series for sinh(rho)/rho and
psi(w)=(sinh(sqrt(w))/sqrt(w))^2 give the following strict reserves
also for the actual n and psi=n^2:

    |rho|<=1.2: |n|<1.3, |n|>.7, |psi_w|<1;
    |rho|<=1.5: |n|<1.5, |n|>.5,
               |psi_w|<2, |K0|<1, |(K0)_w|<1,
    K0=(psi-1)/w.                                  (7)

The errors in (6) are below the margins. In particular
c=beta-1/rho=n'/n is holomorphic, |c|<10 on |rho|<=1.2.
On proper disks of radius .1 about 0<=rho<=1, its derivative bounds
are 10 k!10^k. These disks also suffice for the Ward kernel below.

## 3. Regular pressure and metric integral equations

The exact Ward formula, including its bolt integration constant, is

    Q(rho)=[2N(rho)^2]^-1 integral_0^rho (N^2)'(v)T(v)dv.

For rho<=1 it becomes integral_0^1 H(rho,t)T(t rho)dt, with

    H=[2t n(t rho)^2+rho t^2 (n^2)'(t rho)]/[2n(rho)^2].

Here n(t rho) denotes N(t rho)/(t rho), and the prime is the derivative
in that argument. Formula (7) bounds |H|<10t on the complex .1 disk:
(1.3^2+1.2^2)/(.7^2)<10. Thus the normalized jth derivative of Q is
bounded by Q_near,j a(s), where

    Q_near,j=5K sum_(i=0)^j .1^i, j<=8.             (8)

Only real derivatives of T(t rho) are used in this Leibniz formula;
their weights are at most a(s(rho)). Away from the center use
Q'+2beta Q=beta T. Its value bound is |Q|<=K a/2 because (N^2)'
is positive. Put Q_away,0=K/2 and

    Q_away,n+1=sum_(i=0)^n 5*2^i(K+2Q_away,n-i)/[100(n+1)].

Take Q_j=max(Q_near,j,Q_away,j). These estimates overlap on
.5<=rho<=1. No endpoint value is prescribed independently.

With F=f/C and gamma<g=6.8e-28, the exact equations are

    F'=-L F-2gamma pz(Q-T+E_rho),
    u''+beta u'=b_zz u-(z b_zzz/2)nu
                               +3F-3gamma z^2T+gamma E_u,
    nu''+(beta+2L)nu'=2S u+c_nu nu
                               -F+gamma z^2T+gamma E_nu. (9)

The common-coordinate anomaly cancellation and d-dependent Lie terms
are exactly those in the middle proof. Product subtraction above
leaves every term on the right in the cold ideal.

The original cold radial ledger supplies |F|<1e-10 a and normalized
u,nu derivatives through five below 1e-4 a. This generous low allowance
uses the barrier and regular radial inverse before the old quadratic
normalizer, not an inference from arbitrary chart C5 norms. The low
barrier derivatives are below 3e-6 s(L_c+2), and the old new high jets
are below 1e-8 a, both inside 1e-4 a.

For rho<=1, move c u' (respectively (c+2L)nu') to the right of the
flat radial equation v''+v'/rho=G. Regularity gives exactly

    v''=G-integral_0^1 t G(t rho)dt,
    v^(n+2)=G^(n)-integral_0^1 t^(n+1)G^(n)(t rho)dt. (10)

There is no (n+1)st source derivative. This follows by integrating
v'=rho integral tG and integrating once in t before differentiating.
The formula is valid through the center by even smoothness. Its bound
is at most twice the sup of |G^(n)| on the interior radial interval,
whose weight is at most the observation weight. For rho>=.5, ordinary
differentiation of (9) has the smaller factor one. The regular drift
near the center and the full drift away from it both have normalized
coefficient bounds 12*2^i, so one ledger covers both regions.

Normalize by j!100^j a(s). Starting F_0=1e-10, define

    F_(n+1)={sum 2^i F_(n-i)
       +2g sum 2^i[Q_(n-i)+K+1e22 2^(n-i)]}/[100(n+1)], n<=7.

Set U_j=Z_j=1e-4 through five. For n=4,...,8 define

    C_n(c,X)=sum_(i=0)^n c 2^i X_(n-i),
    D_n(c,X)=sum_(i=0)^n c 2^i 100(n-i+1)X_(n-i+1),
    U_(n+2)=2{D_n(12,U)+C_n(12,U)+C_n(10,Z)+3F_n
                      +3g C_n(2,K)+g 1e24 2^n}/[1e4(n+1)(n+2)],
    Z_(n+2)=2{D_n(12,Z)+C_n(2,U)+C_n(10,Z)+F_n
                      +g C_n(2,K)+g 1e24 2^n}/[1e4(n+1)(n+2)]. (11)

Here C_n(2,K) uses the constant sequence K. The recurrences use only
previous metric indices, so are not circular. Exact evaluation yields
U_10<1.003e15, Z_10<3.341e14 and every entry below H0=2e15.
Thus the original tangent fields have radial C10 allowance H0 a.

## 4. Original finite completion with the weight retained

The middle proof's inverse-tau argument also works at the bolt after
retaining the thermal term. Put alpha=p/tau. On a .05 complex proper
disk its equation is

    alpha'=-tau g0_z/2+delta_bar[z^3/tau^2+3z^2/(2tau)],
    g0=[x(1+2z+3z^2)+z^2+d z^4]/C, tau'=alpha tau.

Bootstrap |alpha|<3. Since real |alpha|<=1.01, tau0>=h/(1+h),
delta_bar/tau0^2<1.002 and exp(.3)<1.35, the derivative is below 2.
The .1 possible change closes the bootstrap. Thus
|tau^-1|<1.2/tau0 on this disk and its normalized coefficient sum
through ten is below 1.5/tau0<2/tau0. Also
a(s)/tau0=(1+s)L_c(s)<82 and a(s)<.00331.

In the Taylor algebra through degree ten, separate degrees <=5 and
>=6. Products of two high parts vanish. The low nu norm is at most
.0006a; hence the low V norm is below .0012*82<.1 and its high norm
is below 821H0. The exact identity

    t_ell/t=1+ell nu h(ell V), h(v)=integral_0^1 exp(tv)dt

has h,h' norms below exp(.1)<1.2 on the low part. The low part of
X=ell nu h(ell V) is bounded by ell*.00072a<ell*2.4e-6. The high part
is bounded by ell*a[.00072+(6+.00072*821)H0]<ell*7H0*a.
Using log(1+X_low+X_high)=log(1+X_low)+(1+X_low)^-1 X_high in this
algebra bounds the full radial coefficient norm by ell*10H0*a,
hence by ell*3e16*a. The orbit field ell u has the smaller reserve
ell*H0*a. The original finite family and ell factor are retained.

## 5. Full tensor conversion at and away from the bolt

For rho>=.5, use nabla(d rho)=beta(g-d rho tensor d rho). Its coefficient
derivatives are bounded by 15 k!200^k after the tensor-product norm
factor below three. The ordered connection recurrence with initial
radial reserve J=3e16 has amplification below 8.64 through ten.

Near the bolt, first convert radial derivatives to Cartesian derivatives
in x=(rho cos theta,rho sin theta). For a smooth even radial scalar f,
subtract its Taylor polynomial at zero of degree less than m. That is
a polynomial in |x|^2 of Cartesian degree less than m, so its mth
Cartesian derivative vanishes. For the remainder,

    |f_rem^(j)(rho)|<=rho^(m-j)/(m-j)!
                          sup_(0<=v<=rho)|f^(m)(v)|.

Expand ordered flat derivatives in products of n=d rho and the flat
metric. Since |n|=1, |g_flat|=sqrt(2) and
D n=(g_flat-n tensor n)/rho, positive coefficient majorants obey

    C_(1,1)=1,
    C_(m+1,j)=C_(m,j-1)+(4m-j)C_(m,j),
    A_0=A_1=1, A_m=sum_(j=1)^m C_(m,j)/(m-j)!.      (12)

Indeed differentiating rho^(j-m) costs m-j and differentiating at
most m factors n costs less than 3m. Thus full Cartesian derivatives
have normalized bounds A_m J a. The tenth A_m is below 4.868e6.
The estimate extends to rho=0 by continuity and retains a(h).

For the physical connection, the Cartesian seed metric is

    g_ij=psi delta_ij-K0 x_i x_j, w=x_1^2+x_2^2.

Its Christoffel symbols have the analytic expression

    Gamma^k_ij=A(x_i delta_kj+x_j delta_ki)+B x_k delta_ij
                                                +D x_k x_i x_j,
    A=psi_w/psi, B=-K0-psi_w,
    D=2psi_w K0/psi-(K0)_w.                        (13)

On complex coordinate polydisks of radius .1 about real |x|<=1,
|x|<1.2, so (7) gives |A|<8, |B|<3, |D|<17 and the Hilbert norm
of Gamma below 100. Cauchy with all ordered coordinate derivatives
gives |partial^k Gamma|<=100 k!20^k, since sqrt(2)/.1<20.
The positive real curvature coefficient b_zz/2 makes N>=rho, hence
g>=g_flat and physical covariant tensor norms are at most the
Cartesian component norms. No extra singular polar-frame factor occurs.

Set B_(0,k)=A_k J k!100^k and use the full ordered recurrence

    B_(r+1,k)=B_(r,k+1)
       +r sum_(i=0)^k binom(k,i)100 i!20^i B_(r,k-i). (14)

Its amplification B_(j,0)/(J j!100^j) is below 1.02e7 through ten.
The tenth entry is below 10186500. Gamma_0 is retained and no higher
tensor is symmetrized. Equations (12)-(14) treat only the region
rho<=1; the regular away estimate covers the rest with overlap.
Since 1.02e7*3e16<1e24, both give (1).

## Verification and consequence

The focused checker verifies seed-series and complex coefficient
margins, exact regular integral and product-subtraction identities,
the rational pressure/metric and tensor ledgers, the original finite
radius algebra, independent Cartesian connection formulas, and
independent regular-inverse integral and derivative evaluations.
The continuum argument is an ordinary proof, not machine verified.
No actual charged kernel is sampled. Failed runs, if any, retain their
source snapshots and raw logs.

Together with the middle result, this supplies two regional metric
inputs through ten. The outer conormal metric estimate, endpoint maps
and correlated mass cubic remain necessary before differentiated
finite-metric quantum and Schur response. The transfer target
R0 R2 R4 N6(p)<=1e57 epsilon_c is still unverified, as are exact signed
compatibility, nonlinear existence and physical stability. No S167
acceptance task or physical gate is completed by this regional result.

Remaining to-do list: outer metric derivatives through ten, endpoint
maps and correlated mass cubic; finite-metric quantum/Schur estimates;
positive or sharper actual transfer; signed compatibility, nonlinear
existence, physical state/boundary and stability; remaining matching,
charge/return, opening, transport, preparation, creation and audit.
