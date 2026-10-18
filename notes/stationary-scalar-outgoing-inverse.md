# The actual scalar outgoing inverse and absence of a physical kernel

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary proof before
focused verification. **PROGRESS; physical source compatibility open.**

At the analytic centered seed, the complete scalar Schur operator has
an outgoing inverse in the source and response spaces specified below.
Its homogeneous physical fast kernel is zero. By the
[actual TT inverse and Schur theorem](stationary-weighted-tt-inverse.md),
the complete mass-extended coupled kernel is therefore zero as well.
The established index -1 then leaves exactly one actual adjoint source
condition. The phase and value of that condition are not evaluated here.

The estimate retains the complete TT operator, its global York maps,
the fixed-density contact and the off-shell return. Its contraction
constant is below 8e-12 after an exact positive-block elimination.
This is a quantitative outgoing construction at the analytic seed;
it is not a numerical fast H4/log inverse constant or a nonlinear
physical exterior/state solution.

## 1. Operators, physical measure and benchmark

Keep the positive prepared temperature, fixed laboratory clock, magnetic
bundle, common subtraction and Dirac realization of the
[centered seed](stationary-centered-covariant-seed.md). Put

    t=r/B, u=1/t, x=lambda B^2, r_h=B(1+h),
    F=lambda r^2+1-2m/r+Q_e/r^2-2dB^3/r^3+dB^4/r^4,
    m=B(1+2x+delta/2), Q_e=B^2(1+3x+d),
    epsilon=G/(60pi), E=1/epsilon+11lambda/3,
    k_rad^2=E-9lambda/4=60pi/G+17lambda/12.

Here E is a spectral energy. Write E_tr for the scalar TRACE DEFECT,
and k_bench=2000 for the preparation parameter. The physical scalar
measure is dmu=4pi beta r^2 dr. Proper radius rho starts at the bolt
and obeys dr/d rho=sqrt(F). All adjoints use this physical measure.

The accepted [geometry](stationary-centered-first-order-geometry.md),
[trace benchmark](stationary-trace-range-condition.md) and
[TT form](stationary-curved-quantum-tt.md) give

    x>.99 k_bench^3>7.9e9, |d|<1e-12, 0<delta<=1,
    .999 k_bench^-10 < G lambda < k_bench^-10,
    Gq^2/(lambda B^4)<4e-27,
    3e-17<sqrt(G lambda)<4e-17,
    k_rad/sqrt(lambda)>4e17, E>180/G.                (1)

The fixed magnetic invariant is M_em=q^2/(2r^4), from the
[exact electromagnetic reduction](exact-spherical-em-reduction.md).
It is not obtained by identifying the comparison coefficient Q_e
with a changed physical charge.

Let T=-Box>=0 have the smooth-bolt scalar realization, and set
Q=(T+4lambda)^(1/2). Let D be the actual scalar-c TT operator of the
preceding theorem, and Bmix its exact mixed block. Its scalar Schur
equation, after multiplication by 80pi^2, is

    P tau=p,
    P=(T+4lambda)(T-E)+W-F_T F_T*,
    F_T=sqrt(80pi^2) Bmix D^(-1/2),
    W=div((Z g/36+11J/9)grad)-Z/(3epsilon),
    Z=R4+12lambda, J=Ric-(R4/4)g.                  (2)

For a full coupled source, p=80pi^2 f_eff with f_eff defined in the
preceding theorem. The scalar variable is the same York coordinate
tau; no source or output transpose is silently changed.

## 2. The actual TT spectral form and global return bounds

The constructed D_crit^-1 is bounded, symmetric, positive, injective
and defined on all physical scalar L2. Thus it is selfadjoint. Its
range is dense, and its inverse D is selfadjoint on the proved
operator domain W_(3/2). This identifies the actual operator; it does
not promote an arbitrary semibounded form to an assumed realization.

The preceding critical graph-domain coercivity extends to the genuine
spectral form domain Dom(D^(1/2)). Indeed spectral truncations lie in
Dom(D), are Cauchy in the D form norm, and hence in the classical TT
energy norm by the lower bound. The latter is closed, being an
ordinary regular-disk and weighted H1 form. Consequently

    ||D^(1/2)c||^2 >= .99 E0[c].                   (3)

For c=D^(-1/2)y, h=Hc=diag(-b-2c,b,c,c), and eta=Vc the actual
outer-fixed return vector, the [seed energy ledger](stationary-seed-quantum-form.md)
gives, in physical tensor/scalar L2 norms,

    ||h||^2 <=176piG E0/lambda,
    ||nabla h||^2 <=41600piG E0,
    ||eta||^2 <=144piG E0/lambda^2,
    E0[D^(-1/2)y] <=||y||^2/.99.                  (4)

These are bounds on the entire global TT reconstruction and Volterra
return, including their fixed tails. They are not pointwise locality
claims about the York map. No weighted D inverse is needed below.

## 3. Reference scalar potential, bolt domain and thermal constant

Set m_rho=r^2 sqrt(F). The unitary Liouville map is
tau -> v=sqrt(4pi beta m_rho)tau. In this section a prime denotes a
proper-radius derivative unless it carries an r subscript. Put
a=sqrt(F)/kappa, A=a'/a and R=r'/r. Then

    H=UTU^-1=-partial_rho^2+V_L,
    V_L=F_rr/4+F_r/r-F_r^2/(16F),
    D_rho=U partial_rho U^-1=partial_rho-(A+2R)/2,
    H=D_rho* D_rho.                               (5)

At the actual bolt r-r_h=kappa rho^2/2+O(rho^4), and

    V_L=-1/(4rho^2)+F_rr(r_h)/6+2kappa/r_h+O(rho^2).

Let D0 be the graph closure of partial_rho-1/(2rho) on
C_c^infinity(0,infinity). Its two closed products are

    H_B0=D0*D0=-partial_rho^2-1/(4rho^2),
    H_B1=D0D0*=-partial_rho^2+3/(4rho^2).           (6)

The first is the critical Friedrichs Bessel realization. Its form
is integral rho |partial_rho(v/sqrt(rho))|^2 d rho. Logarithmic
cutoffs approximate the regular sqrt(rho) branch with vanishing
form cost; ordinary H1 regularity of v itself is not asserted.
The second endpoint is limit point: rho^-1/2 is not L2, whereas
the regular branch is rho^3/2. Neither product has an additional
bolt parameter or an L2 zero vector on the full half-line.

The [seed coefficient bounds](stationary-seed-quantum-form.md) are

    lambda<=F_rr<=21lambda, F/r^2<=5lambda,
    0<=F_r/(2r)<=11lambda.                         (7)

Since a''/a=F_rr/2 and a(0)=0, a'(0)=1, the Wronskian comparison
with rho and sinh(c rho)/c, c=sqrt(21lambda/2), gives

    1/rho<=A<=c coth(c rho),
    0<=A^2-rho^-2<=21lambda/2.

For example (a' sinh(c rho)/c-a cosh(c rho))' is nonpositive;
division by the positive product proves the upper comparison.
The last bound uses sinh(s)>=s. Thus

    U_geo=V_L+1/(4rho^2)-9lambda/4,
    |U_geo|<=30lambda,
    H=H_B0+9lambda/4+U_geo.                        (8)

Also |D_rho-D0|<4sqrt(lambda). The physical form domain is the same
graph domain as D0. Since U_geo is bounded, the operator domain is
the maximal Sturm--Liouville L2 domain with

    lim_(rho downarrow 0) Wronskian(sqrt(rho),v)=0. (9)

For smooth sources this gives v=sqrt(rho)[c+O(rho^2)], excluding
sqrt(rho)log(rho). It is precisely smooth even scalar disk
regularity, not an independently imposed zero scalar value at the bolt.

Here is a quantitative L1 estimate retaining the actual temperature.
Write rho_2=rho(2B). From F_r(B)=delta/B and (7),
2kappa>=lambda B h. For s=r-r_h,
F>=2kappa s+lambda s^2/2. Integration and F<=5lambda r^2 yield

    1/(5sqrt(lambda)) < rho_2
      <=2[2+log(1/h)]/sqrt(lambda).                (10)

The lower bound uses r_h<1.1B and log(20/11)>1/2. The upper
integral is 2sqrt(2/lambda) asinh sqrt(lambda(2B-r_h)/(4kappa));
asinh(1/sqrt(2h))<=log(3/sqrt(h)) gives the stated reserve.

The [prepared cold reference](stationary-reference-cold-bolt.md) has
1/(8 k_bench^10)<=h_class<=1/k_bench^10 and
.998<=gamma_class/h_class<=1.001. The centered theorem gives

    gamma_actual/gamma_class=u_star/sqrt(z_star),
    |u_star-1|<4e-26, |z_star-1|<1.4e-22,
    |h/gamma_actual-1|<=12gamma_actual,
    gamma_actual<1e-30.

Therefore h>h_class/2, log(1/h)<84 and rho_2<172/sqrt(lambda).
No uniform zero-temperature limit is being taken.

For t>=2 the exact polynomial is

    F=x(t^2-4/t+3/t^2)+(1-1/t)^2
                         +d(t-1)^2/t^4-delta/t.

Its leading parenthesis divided by t^2 is at least 11/16; x>=1000
gives F>=lambda r^2/2. If f=F-lambda r^2, absolute coefficient
sums give |f|<=3x, |r f_r|<=4x and |r^2 f_rr|<=9x. Before rounding
these are 2.75x+3.3125, 3.5x+3.5 and 8.5x+10.25. Consequently

    |F_r^2/F-4lambda|<=64x/r^2,
    |V_L-9lambda/4|<=11lambda B^2/r^2.

The first follows by expanding F_r^2-4lambda F and dividing by
F>=lambda r^2/2. The outer integral of the second bound is below
2sqrt(lambda), and integral_(rho_2)^infinity (4rho^2)^-1 d rho
is below 5sqrt(lambda)/4. Combining with (8)-(10) gives

    L_geo=||U_geo||_L1(d rho)
      <[130+60log(1/h)]sqrt(lambda)<6000sqrt(lambda),
    eta_B=L_geo/k_rad<2e-14<1e-12.                 (11)

## 4. Direct outgoing resolvent and derivative bounds

The outgoing Bessel kernel at k=k_rad>0 is

    G0(k;rho,sigma)=(i pi/2)sqrt(rho sigma)
             J0(k min(rho,sigma)) H0^(1)(k max(rho,sigma)).

The coefficient follows from the [Bessel Wronskian](https://dlmf.nist.gov/10.5.E3)
and the derivative jump -1. [Nicholson's formula](https://dlmf.nist.gov/10.9.E30)
and the positive decreasing [K0 representation](https://dlmf.nist.gov/10.32.E9)
give

    J0(s)^2+Y0(s)^2
      =(8/pi^2) integral_0^infinity K0(2s sinh v)dv
      <=(8/pi^2) integral_0^infinity K0(2s v)dv
      =2/(pi s),
    |G0(k;rho,sigma)|<=1/k.                        (12)

Here integral K0(a v)dv=pi/(2a), obtained directly by integrating
exp(-a v cosh y) first. No numerical Bessel enclosure is used.

The complex bound is |G0(zeta;rho,sigma)|<=1/|zeta| for Im zeta>0.
For fixed positive rho,sigma, zeta G0 is analytic there, tends to
zero at the indented origin and has boundary modulus at most one
on both real rays. The [large-argument expansions](https://dlmf.nist.gov/10.17)
and [continuation formulas](https://dlmf.nist.gov/10.11) give limsup
modulus at most one on expanding upper semicircles: the two leading
exponentials have coefficients of modulus 1/2 and do not grow.
Apply the expansions on closed quadrants, using continuation near
the negative ray, and then the maximum-modulus principle.

Factor U_geo=v_geo u_geo with u_geo=|U_geo|^(1/2) and
v_geo=sign(U_geo)|U_geo|^(1/2). The Hilbert--Schmidt norm of
v_geo (H_B0-zeta^2)^-1 u_geo is at most L_geo/|zeta|. Thus the
bounded-potential resolvent identity and a Neumann series give the
actual reference kernel

    G=G0-G0 U_geo G0+G0 U_geo G0 U_geo G0-...,
    |G|<=1/[|zeta|(1-L_geo/|zeta|)].                (13)

For any multiplication weights omega_1,omega_2 in L2(d rho), dominated
convergence proves convergence of the weighted kernels in
Hilbert--Schmidt norm as zeta approaches k from above. The sandwiched
inverse also converges in operator norm. This proves the actual
limiting-absorption resolvent R+=(T-E-i0)^-1 and

    ||omega_1 R+ omega_2||
       <=||omega_1||2||omega_2||2/[k(1-eta_B)].     (14)

Its regular bolt condition is (9). For smooth compact sources it is
outgoing, with Liouville asymptotic v=c exp(i k rho)+o(1) and the
corresponding derivative. A regular outgoing homogeneous solution
obeys the L-infinity contraction v=-G0 U_geo v and hence vanishes.
Equivalently its bolt flux is zero and its outgoing flux is k|c|^2.
For rough weighted sources the norm limit defines the outgoing class;
pointwise asymptotics are not asserted for every such datum.

Derivative estimates retain the bolt contact. Nicholson's formula
at order one, with p=sinh v, gives

    |H1^(1)(s)|^2<=2/(pi s)+1/(pi s^3),
    |H1^(1)(s)|^2<=2/(pi s)+4/(pi^2 s^2).

Use (1+2p^2)/sqrt(1+p^2)<=1+2p^2 or <=1+2p. The K0 moments are
integral p K0(a p)dp=1/a^2 and integral p^2 K0(a p)dp=pi/(2a^3).
The [Poisson integral](https://dlmf.nist.gov/10.9.E4) gives |J0(s)|<=1,
|J1(s)|<=s/2 for s<=1. Split k min(rho,sigma)<=k max(rho,sigma)
according as both, one or neither is below one. These bounds give

    |G1(k;rho,sigma)|<=2/k,
    |D0,rho G0(k;rho,sigma)|<=2,
    D0 (H_B0-zeta^2)^-1 D0*
       =I+zeta^2(H_B1-zeta^2)^-1.                 (15)

For example the larger-end derivative when both arguments are at
most one is bounded by sqrt(2pi+4)/2<2; when both exceed one it is
at most sqrt(6)/2. The closed products (6) prove the last identity
by polar decomposition. There is no bolt rank-one term. The same
maximum-modulus argument gives the complex kernel bounds needed for
norm limits; D0 G0 has a finite limit at zeta=0.

Differentiate the two endpoint free kernels of (13), first as
distributions. Their uniform bounds give

    |D0,rho G|<=2/(1-eta_B),
    D0 R+ D0*=I+K0,
    |K0|<=2k+4L_geo/(1-eta_B).

Since |D_rho-D0|<4sqrt(lambda),

    |D_rho G|<3,
    D_rho R+ D_rho*=I+K,
    |K|<3k.                                       (16)

Before rounding, the last bound is
2k+[4L_geo+16sqrt(lambda)+16lambda/k]/(1-eta_B).
The identity distribution remains explicit. Equations (14)-(16)
give operator-norm derivative limits and

    ||omega_1 partial_rho R+ omega_2||
       <=3||omega_1||2||omega_2||2,
    ||omega_1 partial_rho R+ partial_rho* omega_2||
       <=||omega_1 omega_2||infinity
                       +3k||omega_1||2||omega_2||2. (17)

These are PHYSICAL derivatives before the Liouville transform.
The second inequality also requires the product of weights bounded.
The adjoint one-derivative version has the same estimate.

This direct proof only needs bounded U_geo in L1. The global
Marchenko moment assumption in the critical-angular-momentum
[scattering paper](https://arxiv.org/pdf/1611.05210) is not satisfied
by its 1/(4rho^2) tail and is not invoked.

## 5. Mixed local coefficients and the small trace contact

The exact trace and mixed variation are equations (14)-(17) of the
preceding TT inverse note. For h genuinely TT, set r_dot=-J:h. Then

    T_tr h=(1-11epsilon R4/36)r_dot+epsilon Box r_dot
      -epsilon h:Hess R4+(7epsilon/12)delta C4^2
      +(11epsilon/6)[J:Delta_L h-(R4/2)J:h-2h:J^2]
      +40epsilon h^ab F_ac F_b^c,
    Bmix c=-(T_tr h-div(E_tr eta))/(8piG).          (18)

The full trace is local. The artificial reference cancels in the
full action; it must not be added again to (18).

Here are explicit differentiated curvature bounds, including the
frame variation. With P_orb the orbit-plane projector and
S_sph=g-P_orb,

    J=j(P_orb-S_sph),
    B^2 j=-(1+3x+d)u^4+5d u^5-(9d/2)u^6,
    |nabla P_orb|=2sqrt(F)/r,
    |J|<=7lambda u^4,
    |nabla J|<=100lambda^(3/2)u^4.                 (19)

Indeed |j|<=3.5lambda u^4 and |j'|<=30lambda^(3/2)u^4; differentiating
the projector as well gives (19), with no singular bound on A.
The Weyl tensor is C=w_C W0, where

    w_C=-m/r^3+Q_e/r^4-(10d/3)B^3/r^5+(5d/2)B^4/r^6,
    W0=-P_orb wedge P_orb-S_sph wedge S_sph+P_orb wedge S_sph,
    |W0|=sqrt(48), |nabla W0|<=48sqrt(F)/r.

The last bound follows from
nabla W0=-3(P_orb-S_sph) wedge nabla P_orb and
|A wedge B|<=4|A||B|. Coefficient sums give |w_C|<=5.01lambda u^3
and |w_C'|<=43lambda^(3/2)u^3, hence

    |C|<=40lambda u^3,
    |nabla C|<=1000lambda^(3/2)u^3,
    |Rm|<=54lambda, |R4|<=13lambda,
    |Ric|<=14lambda, |P_schouten|<=5lambda.         (20)

For the last Schouten bound use its orthogonal decomposition
P_schouten=J/2+R4 g/24, not the weaker triangle bound.
The scalar curvature and its radial derivatives give

    Z=(d/B^2)(4u^5-6u^6),
    |Z|<=4|d|lambda u^5/x,
    |Z|<=2|d|lambda u^4/x,
    |Hess R4|<=3652 x|d| B^-4 u^5<4lambda^2 u^5.   (21)

Both Z bounds are stated deliberately: a coefficient 2 multiplying
u^5 would be false near infinity. The Hessian entries are
F Z_rr+(F_r/2)Z_r, (F_r/2)Z_r and twice F Z_r/r; use
|r Z_r|<=56|d|u^5/B^2 and |r^2 Z_rr|<=372|d|u^5/B^2.

After one integration by parts, write

    T_tr h-div(E_tr eta)=-div g_h+f_h,
    g_h=epsilon grad(J:h)-(7epsilon/3) C^{abcd}K_adb
                  +(11epsilon/6)J^{ab}grad h_ab+E_tr eta,
    K^a_bc=(nabla_b h_c^a+nabla_c h_b^a-nabla^a h_bc)/2.

The Weyl vector above carries free index c. Its derivative term is
4C^{abcd}nabla_c K_adb before weighting; |K|<=3|nabla h|/2.
Thus the exact algebraic contacts and (19)-(21) give

    |g_h|<=epsilon[160lambda u^3|nabla h|
                   +100lambda^(3/2)u^4|h|]+|E_tr||eta|,
    |f_h|<=3700epsilon lambda^(3/2)u^3|nabla h|
       +[7lambda u^4+10000epsilon lambda^2u^3
                         +20epsilon q^2 B^-4u^4]|h|. (22)

For an audit of these reserves, the unweighted Weyl variation costs
6(|C||grad sigma|+|nabla C||sigma|)|nabla h| and at most
13000lambda^2u^3|sigma||h|. Its algebra is bounded by
2|C||Rm|+8|C||P_schouten|+4|C|^2<=12320lambda^2u^3.
The J variation costs 7lambda u^4|grad sigma||nabla h|,
100lambda^(3/2)u^4|sigma||nabla h| and
1100lambda^2u^4|sigma||h|. Apply weights 7/12 and 11/6 and add
the Box r_dot and R4-prefactor terms. The derivative totals are
7+140+77/6<160 and 3500+1100/6<3700; the algebra total is below
10000. The magnetic coefficient uses |F_ac F_b^c|<=q^2/(2r^4).

Crucially, the return is kept as div(E_tr eta), requiring only the
trace VALUE bound. Subtracting the chosen AdS vacuum before taking
absolute values gives

    E_tr=(1+11epsilon lambda/3)Z+epsilon Box Z-11epsilon Z^2/72
                       +(7epsilon/12)C^2+(11epsilon/6)J^2-20epsilon M_em,
    |E_tr|<=lambda e_* u^4,
    e_*=4.004|d|/x+1100epsilon lambda
                          +10epsilon q^2/(lambda B^4)<6e-22. (23)

Use |Box Z|<=2|Hess R4| and (19)-(21). The curvature reserves are
8, less than 1, 2800/3 and 539/6, totaling below 1100. The first
term of e_* is below 5.07e-22 by (1), and the others total below
1e-27. No derivative of E_tr or logarithmic nonlocal stress estimate
is substituted here.

## 6. Local response maps before scalar preconditioning

For c=D^(-1/2)y, equation (22) gives

    F_T y=div f_T[y]+g_T[y],
    f_T=(sqrt(5)/(2G))g_h,
    g_T=-(sqrt(5)/(2G))f_h.

The normalization sqrt(80pi^2)/(8piG)=sqrt(5)/(2G) is exact.
With w_3=u^3, (4), (22) and (23) imply

    ||w_3^-1 f_T[y]||<.001sqrt(lambda)||y||,
    ||w_3^-1 g_T[y]||<201sqrt(lambda/G)||y||.        (24)

Explicit majorants before rounding are

    A1=sqrt(5)/(2sqrt(.99)) [epsilon lambda/sqrt(G)
        (160sqrt(41600pi)+100sqrt(176pi))
                              +12sqrt(pi)e_*/sqrt(G)],
    A0=sqrt(5)/(2sqrt(.99))*sqrt(lambda/G)
       [7sqrt(176pi)+epsilon lambda
          (3700sqrt(41600pi)+10000sqrt(176pi))
                  +20epsilon q^2 sqrt(176pi)/(lambda B^4)].

Use pi<22/7, pi>3 and (1). In A1 the return term is bounded by
25e_*/sqrt(G lambda)<.0005 after division by sqrt(lambda).
The LOWER bound on G lambda is essential for this step. A0 is below
sqrt(lambda/G)[200+2e6epsilon lambda+600epsilon q^2/(lambda B^4)].

In particular F_T*Q^-1 is bounded on L2, so
K=Q^-1 F_T F_T* Q^-1 is a bounded positive selfadjoint operator.
This follows from the actual spectral form of D and one integration
by parts, not an assumed D^-1 on an unspecified negative Sobolev space.
Likewise V=Q^-1 W Q^-1 is bounded selfadjoint. Thus
H_eff=T+V-K is selfadjoint on Dom(T).

On radial scalars let a_s=Z/36+11j/9 and set

    U_s tau=(sqrt(|a_s|)grad tau,sqrt(|Z|/(3epsilon))tau),
    J_s=diag(-sign(a_s),-sign(Z)), W=U_s* J_s U_s,
    Fcal=(F_T,U_s*), Jcal=diag(I,-J_s).             (25)

Then Fcal Jcal Fcal*=F_T F_T*-W. The square roots are bounded
measurable form multipliers; their derivatives are not taken.
For w=u, |a_s|<5lambda u^4 and (21) give

    ||w^-1 f_s||<=sqrt(5lambda),
    ||w^-1 g_s||^2<=4|d|lambda/(3epsilon x).

Since w_3<=w, combine these with (24) by Cauchy on the auxiliary
direct sum. There are bounded linear maps f,g with

    Fcal y=g[y]+div f[y],
    ||w^-1 g[y]||<=1000sqrt(lambda/G)||y||,
    ||w^-1 f[y]||<=3sqrt(lambda)||y||.              (26)

All nonlocal York and TT effects remain inside the bounded maps f,g.
The decaying multipliers, rather than weighted norms of D^-1, supply
localization for the following resolvent estimates.

## 7. Exact absorption of the low-frequency contribution

The fourth-order reference outgoing inverse is the exact identity

    R4+=(E+4lambda)^-1[R+-(T+4lambda)^-1]
        =Q^-1 R+ Q^-1.                            (27)

Both resolvents are of the SAME selfadjoint T. The subtracted term
retains all low frequencies. Replacing Q^-1 globally by 1/sqrt(E)
would not establish (27) or the estimate below.

For y=Fcal* tau, equation (2) becomes

    [I+Acal Jcal-Lcal Jcal]y=Fcal* R4+ p,
    Acal=Fcal*(T+4lambda)^-1 Fcal/(E+4lambda)>=0,
    Lcal=Fcal*R+ Fcal/(E+4lambda).                 (28)

These are bounded auxiliary-space maps by (17) and (26). Only the
local scalar block of Acal needs to be small:

    ||Acal_22||<= [||a_s||infinity
                    +||Z||infinity/(12epsilon lambda)]/(E+4lambda)
                <b0=1e-21.                       (29)

Here the global bound ||Z||infinity<=2|d|lambda/x is valid, and
||a_s||infinity<5lambda. The TT block may be large.
Put T1=(I+Acal_11)^-1. Positivity gives

    0<=S22=Acal_22-Acal_21 T1 Acal_12<=Acal_22,
    ||T1 Acal_12||, ||Acal_21 T1||<=sqrt(b0)/2.

For the cross bounds factor a positive block matrix through its
diagonal square roots and use sup_(a>=0) sqrt(a)/(1+a)=1/2.
The local Schur block is I+S22 Jcal_22; its inverse norm is at most
(1-b0)^-1. Triangular factorization therefore yields

    ||(I+Acal Jcal)^-1||
       <=(1+sqrt(b0)/2)^2/(1-b0)<2.                (30)

This bound is independent of the norm of Acal_11.

The actual spatial weight has ||w||infinity<=1 and, by (10)-(11)
and F>=lambda r^2/2 for r>=2B,

    ||w||_L2(d rho)^2<200/sqrt(lambda).

Thus the zero-, one- and two-derivative bounds (14),(17) are bounded
respectively by 201/(k_rad sqrt(lambda)), 600/sqrt(lambda), and
601 k_rad/sqrt(lambda). The last includes the identity contact.
Using (26),

    ||Lcal|| <=(E+4lambda)^-1[
          201*1e6 lambda/(G k_rad sqrt(lambda))
          +1200*3000lambda/(sqrt(G)sqrt(lambda))
          +601*9lambda k_rad/sqrt(lambda)]
       <2.8e-12+8e-13+1.4e-14<4e-12.             (31)

Use E>180/G in the first two terms, and k_rad/E<1/k_rad in the
third. Equation (1) gives the displayed numbers. Consequently

    ||(I+Acal Jcal)^-1 Lcal Jcal||<8e-12<1,
    ||[I+Acal Jcal-Lcal Jcal]^-1||
                            <=2/(1-8e-12).         (32)

This is the quantitative outgoing inverse bound in the declared
auxiliary response space. The TT derivative bound and the small exact
trace contact were necessary; the TT L2 gap alone did not suffice.

## 8. Limiting domain and full-operator reconstruction

Let P_w be the physical weak-source space

    P_w={p=g+div f: w^-1g,w^-1f in L2(dmu)},

with the quotient sum norm. Let V_w be the regular-bolt weighted
H1 output space with norm ||w tau||+||w grad tau||. In the
Liouville variables a source is g_tilde-D_rho* f_tilde. Equations
(14),(17) give R+:P_w->V_w, including its operator-norm limit.
Because w<=1, P_w embeds in the ordinary scalar H^-1 form dual;
the positive elliptic resolvent (T+4lambda)^-1 maps it into the
ordinary H1 form domain and hence V_w. Thus R4+ is defined on P_w.
Equation (26) gives Fcal:H_aux->P_w and its paired map
Fcal*:V_w->H_aux, without a weighted estimate for Q or D^-1.

For p in P_w solve (28) and define

    tau=R4+(p+Fcal Jcal y).                        (33)

Applying Fcal* returns y exactly, so (33) solves the complete
scalar equation weakly. Smooth compact source approximations and
the norm limits above justify each identity. The outgoing domain
for rough data is this limiting-absorption realization, rather than
a pointwise asymptotic postulate for every rough source.

This realization is also the one of the actual compressed operator
H_eff. For nonreal z near E keep Fcal and Q fixed and replace (27)
by (z+4lambda)^-1[(T-z)^-1-Q^-2]. The corresponding Acal(z) is a
scalar multiple of Acal(E). The bounded sandwiches are norm
continuous up to E; (30)-(32) remain invertible in a neighborhood.
Away from the real axis the ordinary resolvent factorization
identifies (33) with Q^-1(H_eff-z)^-1 Q^-1 on smooth sources.
Taking the norm limit identifies the outgoing map, rather than
introducing a different boundary operator.

The fourth-order bolt condition is retained as well. For p in P_w,

    T R4+ p=[E R+ p+4lambda(T+4lambda)^-1p]/(E+4lambda).

The distributional copies of p cancel in this identity. Both tau
and T tau therefore have the regular scalar form boundary condition;
locally tau is H3 for weak H^-1 data. For smooth fast physical
forcing the actual local elliptic regularity and the mixed block's
order at most two bootstrap this solution into the previously fixed
scalar and TT domains, allowing the outgoing oscillatory component.
The reconstructed TT field is

    c=D_crit^-1(f_c-Bmix* tau).                    (34)

Its right-hand side belongs to the original global critical L2
input domain. Indeed w tau,w grad tau in L2 is outer weight delta=1/2
in the accepted physical b-Sobolev convention. Global weighted
elliptic regularity for (T+4lambda)(T-E), with its divergence-L2
source, gives tau in z^(1/2)H_b^3; ordinary disk regularity handles
the bolt. The mixed adjoint has order at most two and gains at least
two powers of outer decay. The established sharp York-transpose
map in -1<delta<4 therefore gives Bmix* tau in z^(5/2)H_b^1 plus
its retained z^4 tails, which is physical L2.
The weak reconstruction first uses the genuine spectral D form;
elliptic regularity identifies it with (34). On fast sources the
critical inverse agrees with the mass-extended inverse by uniqueness.

## 9. Physical kernel and the one remaining source condition

There is no unaccounted free homogeneous remainder for a physical
fast kernel vector. Such a tau is smooth at the bolt by the existing
elliptic theorem, belongs to Z_(7/2), and has Fcal* tau in H_aux.
Its source Fcal Jcal Fcal* tau has the controlled smooth outer tails.
Use its own y=Fcal* tau, before solving any auxiliary equation, and set
q=tau-R4+[Fcal Jcal Fcal* tau]. This difference satisfies

    (T+4lambda)(T-E)q=0,
    q_osc=(T+4lambda)q/(E+4lambda),
    q_ell=-(T-E)q/(E+4lambda), q=q_osc+q_ell.

Both components are regular at the bolt. The elliptic component
has no growing outer branch, so is its decaying L2 branch;
T+4lambda>0 forces it to vanish. The oscillatory component is
outgoing, and section 4 forces it to vanish. Thus q=0, and the
physical fast kernel really obeys (28). Equation (32) gives y=0
and hence tau=0. The TT inverse then gives c=0.

Therefore the [mass-extended map](stationary-global-weighted-map.md)
and its scalar Schur reduction satisfy

    ker A_m=ker S_eff={0},
    index A_m=index S_eff=-1,
    dim coker A_m=dim coker S_eff=1.               (35)

They have bounded inverses on their closed ranges. Equation (32)
quantifies an outgoing response norm; it is not a numerical bound
for the full fast weighted H4/log inverse.

By the [actual adjoint theorem](stationary-adjoint-phase.md), the
mass-orthogonal adjoint space is now exactly one real line, with a
nonzero scalar oscillatory leading coefficient. If that coefficient
also vanished, the field would be a physical mass/fast kernel vector,
contradicting (35). The unique adjoint line therefore has the actual
scalar phase, together with its forced TT corrections. Its TT field
is -D_crit^-1 Bmix* psi, and the necessary and sufficient fast-range
condition is still the actual convergent moment

    <chi,f>=<psi,f_tau-Bmix D_crit^-1 f_c>=0.       (36)

The phase angle and value of (36) for the retained preparation source
have not been evaluated. In particular an outgoing solution of an
arbitrary source can have a nonzero oscillatory amplitude and fail
the physical fast condition. Neither the small contraction constant
nor a small equation residual sets this moment to zero.

## 10. Verification and remaining scope

The [focused checker](../checks/stationary_scalar_outgoing_inverse.py)
audits finite curvature, Liouville, normalization, block and rational
budget identities. Its immutable
[receipt](../evidence/stationary-scalar-outgoing-inverse.json) and
[raw log](../evidence/stationary-scalar-outgoing-inverse.log) record
the exact executed inputs. Endpoint realizations, spectral form
extension, complex maximum-modulus/Born limits, weak reconstruction
and physical kernel exclusion remain ordinary analytic arguments.
The checker does not solve a physical exterior or evaluate (36).

All results concern the analytic centered seed and the retained
theory. Transfer to a finite nonlinear candidate, its regularity,
physical logarithmic phase or controlled spatial reduction, and
exterior/state and charge matching remain open. No S139.1-S139.2
or S074.2-S074.6 gate closes by this linear theorem.

Remaining to-do list: evaluate the unique actual scalar phase and
source moment, quantify the compatible fast weighted inverse;
transfer regularity and close nonlinear physical logarithmic phase
or justify controlled spatial reduction; match exterior/state and
charge, transfer S138, and complete driven stress, radius, controller,
capture, tail, spatial area and S074-S080.
