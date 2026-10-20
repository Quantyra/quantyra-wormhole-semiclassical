# TT form regularity and the global compatible scalar fourth graph

2026-10-20 HST, E003/S165 supporting S139/S074. Ordinary proof and
verification contract before execution. **The scalar fourth graph is
controlled; the complete metric fourth/log inverse remains open.**

Keep the analytic seed, full action, paired domain, temperature,
fixed clock and correlated mass column. Use all source conventions of
[S164](global-compatible-scalar-response.md):

    g=G lambda, W=(r/B)^2, Tbar=(-Box_4)/lambda,
    gradbar=lambda^-1/2 grad, Dbar=D/lambda^2,
    (fbar_tau,fbar_c)=lambda^-2 Rcal* f,
    Fsrc=||W fbar_tau||+||W fbar_c||,
    dmu=4pi beta r^2 dr.

For the actual physical compatible solution
u=G0 f=Rcal(tau,c), this note proves

    ||W Tbar^2 tau|| <=6.2e11 g^(-1/2) Fsrc.       (1)

The lower scalar graphs retain S164's bounds. In addition, writing
h_c for the genuine TT tensor and T_c for its disk shape return,

    ||Dbar^(1/2)c||      <=3500 sqrt(g) Fsrc,
    ||Box_2 T_c/lambda|| <=3.2e5 sqrt(g) Fsrc,
    ||gradbar^2 h_c||    <=1.8e6 sqrt(g) Fsrc.       (2)

The norms in (2) are UNWEIGHTED physical L2 norms. Here
Box_2=F partial_r^2+F_r partial_r is the physical two-dimensional
orbit scalar Laplacian; it differs from Box_4 by 2F partial_r/r.
The fast weight appears explicitly in (1). No bound for the mass
coefficient or fast fourth/log shape remainder is inferred from (2).
No numerical source transpose or disk reconstruction norm, nonlinear
solution, actual preparation moment or stability result is supplied.

## 1. Recover the positive margin on the actual form domain

Section 10 of the [curved TT proof](stationary-curved-quantum-tt.md)
retains, before its final rounding,

    <c,Dc> >= (1-.004)E0[c]+C/(100pi) sum_j J2[chi_j T_c].

These are the actual disk Hessian and its full retained allocations.
The quoted .995/.99 energy gaps discarded the positive term.
Let Q=sum |grad_gbar chi_j|^2=theta'^2. Since sum chi_j^2=1,

    Delta T=sum chi_j Delta(chi_j T)+Q T,
    sum J2[chi_j T] >= J2[T]/2-integral Q^2 T^2.

The cutoff has |theta'|<=5/a=500, with a=1/100. Its support is
in the outer transition, where C J0,outer/E0<3.4e-37. Thus the
additional loss is less than

    (500^4/(100pi))*3.4e-37 E0 <1e-28 E0.

The actual orbit metric is gbar=(C/r^2)(F dt^2+dr^2/F).
Its two-dimensional conformal transformation gives exactly

    C J2[T]=beta integral r^2(Box_2 T)^2 dr
           =||Box_2 T||_dmu^2/(4pi).

Consequently, with room also for the earlier off-shell contact reserve,

    <c,Dc> >=.99 E0[c]+||Box_2 T_c||^2/(800pi^2).  (3)

This first holds on the smooth regular-bolt compact-TT core.
Use the actual critical graph-domain density proved in
[the TT inverse](stationary-weighted-tt-inverse.md), section 3.
For core approximants, apply (3) to their differences. Continuity of
<c,Dc> in the actual graph norm makes Box_2 T_c Cauchy in physical
L2. The order-zero return converges in L2; distributional closedness
identifies the derivative limit. Thus (3) holds on Dom(D).

Next take spectral truncations of the already identified positive
selfadjoint D. For any c in Dom(D^(1/2)), (3) again makes both the
classical energy and Box_2 T_c Cauchy. Their closed limits give (3)
on the actual spectral form domain. This does not assert membership
in the older absolute Weyl/model norm N. The cancellation at the
outer end is retained, and compact-TT density at fast weight 7/2
is not assumed.

In particular, for c=Dbar^(-1/2)y,

    ||Box_2 T_c/lambda|| <=sqrt(800pi^2)||y||<89||y||. (4)

The earlier energy ledger also gives, with alpha=sqrt(g),

    ||c||<5.1alpha||y||,       ||gradbar c||<34alpha||y||,
    ||h_c||<24alpha||y||,      ||gradbar h_c||<364alpha||y||,
    ||sqrt(lambda)eta||<22alpha||y||,
    ||div eta||<150alpha||y||.                     (5)

For the c gradient use ||c'||^2<=2D_energy+40lambda C0 and
E0=beta(3D_energy/2+lambda C0)/(2G), hence its physical bound
320piG E0. For the last estimate,
div eta=2A eta+(b+c)+2R eta,
||A eta||<=3||c|| and ||eta||<=sqrt(18/lambda)||c||,
so ||div eta||<28||c||. All constants include the physical volume.

## 2. Two derivatives of the entire TT lift, including the bolt

Primes in this section mean physical proper-radius derivatives.
Set A=a'/a, R=r'/r, and

    K=F_rr/2-2AR+R^2.

The TT and return equations are

    b'=-2A(b+c)-2R(b-c), eta'=A eta+b+c,
    T_c=-(b+3c)/2-(A-R)eta,
    A'=F_rr/2-A^2, R'=AR-R^2.

Differentiate before estimating. Exact cancellation gives

    T_c'=-3c'/2+2Rb-K eta,
    Box_2 T_c=-3Box_2 c/2-(6R^2+K)b
                +(4R^2-4AR-K)c-(K'+2AK)eta.       (6)

Thus no norm of A alone, which diverges at the bolt, is needed.
The seed polynomial yields

    K=1/r^2-6m/r^3+6Q_e/r^4-20dB^3/r^5+15dB^4/r^6,
    |K|<40lambda, |K'|<300lambda^(3/2).

Here m=B(1+2x+delta/2), Q_e=B^2(1+3x+d). For x>=1000,
0<delta<=1, |d|<=1, the absolute coefficient sums for K/lambda
and r K_r/lambda are below 31 and 109. Multiply the latter by
|R|<=sqrt(5lambda) for its proper derivative.
Using ||b||<=2||c||, ||A eta||<=3||c|| and the return norm in
(5), equation (6) implies

    Qc:=||Box_4 c||
       <=(2/3)||Box_2 T_c||+1300lambda||c||
                                      +5sqrt(lambda)||c'||. (7)

The regular-disk scalar Bochner identity, with |Ric|<=14lambda,
gives ||Hess c||<=Qc+4sqrt(lambda)||c'||. In particular it
controls ||A c'||, the time-time component of Hess c.

The tensor rough Laplacian is not the scalar Laplacian applied
componentwise. Its radial and sphere frame entries are exactly

    (Box h)_rr=b''+(A+2R)b'-4A^2(b+c)-4R^2(b-c)
             =-F_rr(b+c)-2R^2(b-c)-2(A-R)c',
    (Box h)_ss=Box_4 c+2R^2(b-c),
    (Box h)_tt=-(Box h)_rr-2(Box h)_ss.             (8)

The two sphere entries coincide. Equations (7)-(8), F_rr<=21lambda,
R^2<=5lambda, and scalar Bochner imply

    ||Box h||<=7Qc+26sqrt(lambda)||c'||+276lambda||c||. (9)

For completeness, a covariant tensor estimate supplies the remaining
second derivatives. The third-derivative Ricci identity is equation
(7) of [Semmelmann--Weingart](https://arxiv.org/pdf/1708.04775),
with zero torsion. Applying it to a covariant two-tensor and integrating
by parts gives the bound

    ||grad^2 h||^2 <=||Box h||^2
       +(||Ric||inf+4||Rm||inf)||grad h||^2
       +4||grad Rm||inf ||h||||grad h||.            (10)

There is one Ricci contribution from the derivative index, two
curvature actions on each of h's two slots, and two differentiated
curvature contractions. Each latter contraction costs at most
sqrt(4)=2 in four dimensions. These counts give the displayed
constants; no unspecified elliptic constant is used.

The outgoing proof has |Rm|<=54lambda, |grad C|<=1000lambda^(3/2),
|grad J|<=100lambda^(3/2), and |grad R4|<lambda^(3/2).
Differentiating Rm=C+(g wedge J)/2+R4(g wedge g)/24, with
|A wedge B|<=4|A||B|, gives |grad Rm|<2000lambda^(3/2).
Thus (10) and Young's inequality imply

    ||grad^2 h||<=||Box h||+100sqrt(lambda)||grad h||
                                      +100lambda||h||. (11)

All integrations start on regular-bolt compact TT fields. No bolt
boundary is removed. Graph approximation, then the spectral-form
approximation used above, extend (6)-(11). In particular their
right sides bound differences, so the weak derivatives have the
claimed limits. These are global estimates on the actual exterior.

Combining (4)-(11), g<2000^-10, gives for c=Dbar^(-1/2)y

    ||Box_4 c/lambda||<61||y||,
    ||Box h_c/lambda||<430||y||,
    ||gradbar^2 h_c||<500||y||.                    (12)

## 3. The mixed TT vertex is a strong fast source

Let V_T=sqrt(80pi^2) Bbar Dbar^(-1/2) be the TT part of S164's
dimensionless vertex. We now improve its weak-source representation to

    ||W V_T y|| <=250 g^(-1/2)||y||.              (13)

Use the exact local trace and off-shell contact, not a truncated
determinant. For h=h_c, the scalar curvature variation is especially
simple: r_dot=-J:h=-4j c, where J=j(P_orb-S_sphere). The established
curvature bounds and the polynomial j give

    |j|<=3.5lambda u^4, |j'|<=30lambda^(3/2)u^4,
    |Box_4 j|<2000lambda^2 u^4, u=B/r.

Hence Box_4 r_dot costs at most
14lambda u^4 |Box_4 c|+240lambda^(3/2)u^4|c'|
+8000lambda^2u^4|c|.
The exact Weyl variation costs at most
6|C||grad^2 h|+12320lambda^2u^3|h|.
The exact J variation uses Delta_L h and its algebraic curvature
terms. Combining all terms in the trace gives, with
e=epsilon lambda=g/(60pi) and
v=epsilon q^2/(lambda B^4),

    lambda^-1 |T_trace h|
      <=7u^4|h|+e[14u^4|Box_4 c/lambda|
        +240u^4|gradbar c|+8000u^4|c|
        +140u^3|gradbar^2 h|+13u^4|Box h/lambda|
        +20000u^3|h|]+20v u^4|h|.                 (14)

The 20000 reserve includes the R4 prefactor, Hess R4, Weyl and
Lichnerowicz algebraic terms and the retained fixed-density contact
allocation of the mixed block. The separate off-shell return is
still -div(E_trace eta); it has not been set to zero.

We need its derivative, which can be computed from the exact local
trace polynomial. Define z=Z/lambda, jbar=j/lambda and wbar=w_C/lambda.
With f=F/(lambda r^2),

    f=1+u^2/x-(2+4x+delta)u^3/x
        +(1+3x+d)u^4/x-2d u^5/x+d u^6/x,
    L[p]=f u^2 p_uu+(u^2 f_u-2fu)p_u=Box_4 p/lambda,
    z=(d/x)(4u^5-6u^6),
    jbar=-(1+3x+d)u^4/x+5d u^5/x-(9d/2x)u^6,
    wbar=-(1+2x+delta/2)u^3/x+(1+3x+d)u^4/x
                          -(10d/3x)u^5+(5d/2x)u^6,
    E_trace/lambda=(1+11e/3)z+e L[z]-11e z^2/72
                       +28e wbar^2+(22e/3)jbar^2-10v u^4. (15)

For |u|<=1, the absolute coefficient sums of f,u f_u,u^2 f_uu
are below 9,50,250. Those of jbar,u jbar_u are below 3.6,20,
and those of wbar,u wbar_u below 5.1,40. For z they are
10|d|/x, 56|d|/x and 260|d|/x through two derivatives.
Thus L[z] has coefficient sum <7000|d|/x, degree at most 12,
and L[jbar] has sum <2000. Every term in E_trace has at least
four powers of u. Since |gradbar u|<=sqrt(5)u, (15) gives

    |grad E_trace|<4e-20 lambda^(3/2)u^4.         (16)

One explicit majorant is
200|d|/x+1e5 e+100v, using
|d|/x<1.3e-22, e<2000^-10/180 and v<4e-27/180.
The e|d|/x terms fit the reserves in that majorant.
Together with the earlier |E_trace|<6e-22 lambda u^4,
this bounds the full divergence contact by

    lambda^-1 |div(E_trace eta)|
      <=6e-22 u^4 |div eta|
                          +4e-20 u^4|sqrt(lambda)eta|. (17)

Finally V_T y=-(sqrt(5)/(2G lambda))
[T_trace h-div(E_trace eta)]. Multiply (14),(17) by W=u^-2.
All remaining powers of u are at most one. Insert (5),(12);
the leading coefficient is bounded by
(sqrt(5)/2)*168/sqrt(g), and all other terms fit below
the reserve to 250/sqrt(g). This proves (13) without differentiating
Dbar^(-1/2), or the square-root factors of the scalar local vertex.

## 4. Recover the full scalar fourth graph

S164 gives the actual auxiliary solution y=V* tau with

    ||y||<=3sqrt(g) A+6100sqrt(g) Bv,
    A<=800||W fbar_tau||+30000||W fbar_c||,
    Bv<=.2sqrt(g)||W fbar_c||.

Therefore ||y||<91000sqrt(g) Fsrc; its TT component obeys the
same bound. The Schur source p also has a strong representation:

    p=80pi^2 fbar_tau-sqrt(80pi^2)V_T Dbar^(-1/2)fbar_c,
    ||W p||<=36000 Fsrc.                          (18)

Here ||Dbar^(-1/2)||<=sqrt(8pi g/.99); the coefficient in (18)
uses sqrt(80pi^2)*250*sqrt(8pi/.99)<36000.

Keep the scalar local trace correction directly:

    Pbar=P0+Cbar-V_T V_T*, P0=(Tbar+4)(Tbar-E),
    Cbar=divbar(abar gradbar)-z/(3e),
    abar=z/36+11jbar/9, E=60pi/g+11/3.

Radiality is used here; J's radial component is j.
The bounds |abar|<=5u^4, |gradbar abar|<=40u^4 and
|z|<=2|d|u^4/x, combined with S164's value, gradient and second
graph estimates, give

    ||W Cbar tau||<1e-6 Fsrc.                    (19)

In detail its three constants are bounded by
5(3.2e9)sqrt(g), 40(4e8)g^(3/4), and
40pi(|d|/x)(4.4e7), respectively. This avoids differentiating
a nonsmooth sign/square-root decomposition of the local scalar vertex.

For q=P0 tau=p+V_T y_T-Cbar tau, (13),(18),(19) imply
||W q||<2.4e7 Fsrc. The exact equation now yields

    Tbar^2 tau=q+(E-4)Tbar tau+4E tau.

Since 4<E<190/g, its weighted norm is at most

    [2.4e7+6.08e11/sqrt(g)+3.344e10]Fsrc
       <6.2e11/sqrt(g) Fsrc,

proving (1). The equality is initially distributional; the actual
regular fourth-order bolt domain and the obtained strong source
identify its physical graph. No new bolt or outer condition is imposed.
The compatibility assumption remains exactly S164's physical signed
range condition. It is not replaced by the reference moment of p.

Finally,

    Dbar^(1/2)c=Dbar^(-1/2)fbar_c-y_T/sqrt(80pi^2).

Equations (4),(12) and the auxiliary bound prove (2). The mass and
fixed radial tail remain in the actual TT response throughout.

## Frozen verification contract

Freeze this proof and checker before execution. Revalidate the receipts
and exact current inputs of the curved TT, weighted TT inverse,
outgoing scalar and S164 results. Check every new normalization and
constant budget by exact rational arithmetic, using 3<pi<22/7.

Verify the two return identities (6) by symbolic differentiation
under the exact TT, return and seed identities. Independently compute
the four-dimensional Christoffel symbols and tensor rough Laplacian
for diag(a(rho)^2,1,r(rho)^2,r(rho)^2 sin(theta)^2); compare its
frame components with (8) before the TT substitution. Verify the
substituted formula, trace cancellation, the scalar radial operator,
and the polynomial K and local trace formulas.

Check all polynomial coefficient envelopes by exact coefficient
majorants on x>=1000, 0<delta<=1, |d|<=1. Check the tighter |d|/x,
g and magnetic bounds when forming the final constants.
Use the regular manufactured bolt
a=rho, r=constant,
c=c0+c2 rho^2+c4 rho^4,
b=-c0-c2 rho^2/2-c4 rho^4/3 and
eta=eta1 rho+c2 rho^3/4+c4 rho^5/6.
Check the equations, finite bolt limits and both reconstructed graph
identities symbolically. This is a geometric diagnostic, not an
actual prepared seed or a wormhole solution.

No numerical exterior or free-boundary fit is needed for these exact
identities. The continuum form/domain passage, Cauchy--Schwarz,
integration by parts and spectral approximation remain ordinary
proofs, not conclusions from a finite sampling test. Retain every
scientific run and exact failed sources if any.

Remaining to-do list: quantitative fast fourth/log shape inverse with
its mass coefficient, source/reconstruction norm conversions, nonlinear
family and actual signed compatibility; exterior/state, charge and
return; opening, stability, transport, preparation, creation and audit.
