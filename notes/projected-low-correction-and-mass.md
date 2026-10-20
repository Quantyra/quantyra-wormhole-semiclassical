# Sharper projected correction in value and mass norms

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **Small value and mass upper
bounds in the declared normalization; no preserved nonlinear ball.**

Keep the original physical action, analytic seed, positive temperature,
clock, state, paired domain, full angular tower and correlated mass
column. Let f_ell be the actual fixed-volume comparison residual of
[S166](matched-residual-inverse-composition.md), 0<=ell<=1, and Pi0 its
orthogonal projection onto the actual seed range in Fdisk. For
h_ell=-A_m^-1 Pi0 f_ell, write h_ell=Rcal(tau,c)=(S,T). Then

    ||W tau||/sqrt(gamma0) <2.2e-5 ell^2,
    ||(S,T)||/sqrt(gamma0) <1.2e-6 ell^2,
    |m_c|                  <1.2e-4 ell^2.                  (1)

Here g=G lambda, W=(r/B)^2, dmu=4pi beta r^2 dr and
gamma0=4pi beta B^3. The full metric value norm is UNWEIGHTED and
includes the mass column; the scalar value norm is fast weighted.
m_c is the actual sphere-component TT mass coefficient of
[the mass functional](tt-mass-numerical-bound.md). Its current disk
column conversion is m_Y=-3 k_geom^3 m_c as specified by
[the higher York theorem](quantitative-higher-york-isomorphism.md).
These are upper bounds, not values or errors to an exact solution.
They do not certify the cold C5 neighborhood or a candidate Hessian
inverse. Every S167 acceptance checkbox remains open.

## 1. Extract the positive elliptic part of the actual inverse

Use [S164's physical source and operator conventions](global-compatible-scalar-response.md):

    Tbar=(-Box_4)/lambda, gradbar=lambda^-1/2 grad,
    (fbar_tau,fbar_c)=lambda^-2 Rcal* f,
    Fsrc=||W fbar_tau||+||W fbar_c||,
    H0=Tbar+4, E=60pi/g+11/3, e0=E+4, k^2=E-9/4,
    Pbar=H0(Tbar-E)-V J V*, V=(V_T,V_s), J=diag(I,-J_s).

All Hilbert norms retain the physical volume. The retained guards give
E>180/g, sqrt(E)<15/sqrt(g), k sqrt(g)>13, sqrt(g)<4e-17,
g^(1/4)<6e-9 and g^(-1/4)<2e8. The actual Schur forcing is

    p=80pi^2(fbar_tau-Bbar Dbar^-1 fbar_c)=a+divbar b,
    A=||W a||<=30000 Fsrc, Bv=||W b||<.2sqrt(g)Fsrc.

The exact outgoing solve has y=V*tau, q=p+VJy and
tau=(R+q-R_e q)/e0, where R_e=H0^-1. Compatibility belongs to q,
not necessarily p. The existing estimates give

    ||W tau||<=4.4e7 g Fsrc,
    ||W Tbar tau||<=3.2e9 sqrt(g)Fsrc,
    A_q<=3001A+6100000Bv<9.004e7 Fsrc,
    B_q<=9sqrt(g)A+1.001Bv<270001sqrt(g)Fsrc.

Put z=R+q=H0 tau. The partial-fraction identity proves the equality,
including weak forcing. Thus ||z||<=||Wz||<3.201e9 sqrt(g)Fsrc.
It is useful to obtain its unweighted gradient without differentiating
the source. Since (Tbar+1)z=q+(E+1)z is in H^-1, the regular scalar
resolvent realizes z in H1. The known regular-bolt realization removes
the logarithmic homogeneous extension; the L2 tail removes an infinity
boundary contribution. Equivalently approximate in the regular form
domain and use uniqueness of its shifted positive resolvent. Testing
(Tbar-E)z=q and taking absolute values gives

    ||gradbar z||^2<=E||z||^2+A_q||z||+B_q||gradbar z||,
    ||gradbar z||<=sqrt(E)||z||+sqrt(A_q||z||)+B_q
                  <4.9e10 Fsrc.                         (2)

Separate the positive TT vertex from the small signed local vertex:

    B_T=V_T/sqrt(e0), B_s=V_s/sqrt(e0),
    K=H0+B_T B_T*-B_s J_s B_s*.

[The outgoing proof](stationary-scalar-outgoing-inverse.md), (29), gives
||B_s* H0^-1/2||^2<=b0=1e-21. The vertices are bounded from H1 to
their auxiliary spaces; their adjoints are form maps into H^-1. Thus
the displayed K is a closed coercive form, with

    K>=(1-b0)H0+B_T B_T*.

No pointwise or commuting interpretation of the TT vertex is used.
Factorization through K^-1/2 and H0^-1/2 gives

    ||K^-1(a+divbar b)||<=(A/4+Bv/2)/(1-b0),
    ||K^-1 B_T||<=1/[2sqrt(1-b0)],
    ||K^-1 B_s||<=sqrt(b0)/[2(1-b0)].

The exact rearrangement, as an identity of form maps, is

    tau=z/e0-K^-1 p/e0-K^-1 V J V* z/e0^2.              (3)

Indeed the outgoing identity first gives K tau=(H0 z-p)/e0;
then K^-1 H0 z=z-K^-1 VJV*z/e0. This is valid without commuting
H0 with either vertex. Their separate coefficient bounds imply

    ||V_T*z||<=201||z||/sqrt(g)+.001||gradbar z||
                <6.5e11 Fsrc,
    ||V_s*z||<=1000||z||/sqrt(g)+3||gradbar z||
                <3.4e12 Fsrc.

Substitute in (3), using e0>180/g. The leading coefficient after
division by g Fsrc is 30000/(4*180); all z and vertex contributions
have an additional sqrt(g), while the Bv contribution is smaller.
With (1-b0)^-1<1.001, sqrt(b0)<1e-10 and 180^(3/2)>2400,

    ||tau||<43 g Fsrc,
    ||gradbar tau||<3.8e5 g^(3/4)Fsrc,
    ||Hessbar tau||<3.3e9 sqrt(g)Fsrc.                    (4)

The gradient follows from ||gradbar tau||^2=<tau,Tbar tau>.
The unweighted Bochner bound ||Hessbar tau||<=||Tbar tau||+
4||gradbar tau|| uses the retained |Ric|<=14lambda. The regular
scalar graph and cutoff approximation justify both identities.

## 2. Retain coefficient decay in the strong TT forcing

Write w=B/r (distinct from comparison u below). Equation (6) of
[the strong-source proof](strong-effective-tt-source.md) gives

    lambda^-1 |Ttrace*tau|
      <=(7+20v)w^4|tau|
        +e[200w^3|Hessbar tau|+8000w^3|gradbar tau|
                                     +5e5w^3|tau|],
    e=g/(60pi), v<4e-27/180, |Ebar|<6e-22w^4.

Use its actual transpose norms 5 for H* and 27 for the return Vbar*.
Multiplication by W=w^-2 leaves w or w^2, both at most one. Therefore
UNWEIGHTED norms from (4) suffice in the right side:

    ||W Bbar*tau||<=(8pi g)^-1 {
      5[(7+20v+5e5e)||tau||+200e||Hessbar tau||
                                      +8000e||gradbar tau||]
                   +27(6e-22)||gradbar tau||}<63 Fsrc.

The last return term costs g^(-1/4), bounded by 2e8; it is included.
The leading budget is 5*7*43/24<62.709. Consequently

    ||W f_TT||<64 Fsrc, f_TT=fbar_c-Bbar*tau,
    ||c||<1664 g Fsrc,
    ||(S,T)||<=||tau||+25||c||<42000 g Fsrc,
    |m_c|<4.224e6 g gamma0^(-1/2)Fsrc.                    (5)

Here ||Dbar_crit^-1||<26g, the actual unweighted York return is
bounded by 25, and the actual mass functional has norm
6.6e4 g gamma0^(-1/2) on the strong weighted TT source. The physical
mass column is included, with the same critical inverse and fixed
return tails. These facts are retained from the cited mass and
[source-conversion](weighted-york-source-conversion.md) proofs.

For a sharper fast scalar value, use the vertices separately once more.
The TT a,b coefficients are 201/sqrt(g), .001 with w^3 decay.
The local ones are sqrt(40pi|d|/(gx))<1.3e-10/sqrt(g), sqrt(5)<2.24
with w^2 decay, since |d|/x<1.3e-22. Equation (4) implies

    ||y_T||<8644sqrt(g)Fsrc, ||y_s||<.006sqrt(g)Fsrc,
    A_q<1.768e6 Fsrc, B_q<8.86sqrt(g)Fsrc.

Apply S164's compatible reference estimates to this SAME q:
||W tau||<=(81A_q+1100B_q)/E and
||W Tbar tau||<=348A_q/k+701B_q. Thus

    ||W tau||<8e5 g Fsrc,
    ||W Tbar tau||<4.8e7 sqrt(g)Fsrc,
    ||W gradbar tau||<6.3e6 g^(3/4)Fsrc.                 (6)

The last line follows from the weighted scalar form identity with
|gradbar log W|^2<=20. No reference moment is imposed on p alone.

## 3. Improve the actual comparison's low jets

Use the exact formulas and guards in
[the completed metric proof](stationary-completed-metric-response-norm.md)
and [the first geometry reconstruction](stationary-centered-first-order-geometry.md).
Here t=r0/B, s=t-1, V=xi/(Bs), u=delta sigma, nu=sV/t, b=s V_t.
Let a,f,s_R be the Eulerian tangent and moving-bolt displacement,
fbar=f+s_R F0_r and H=a+fbar/F0. Then

    |H|<1.08e-9, |a|<3.3e-10, |fbar/F0|<7.5e-10,
    |s_R|/B<3.4e-10 h, h<1e-30,
    J=B integral_t^infinity H/F0,
    xi=s_R-F0 J,
    u=a+fbar/(2F0)-s_R/(Bt)
      -(F0_t/2-F0/t) integral_t^infinity H/F0.

The seed is F0=s^2 f(t)-delta/t, with
f=x+2x/t+(3x+1)/t^2+d/t^4, C=1+6x+d, x>=1000, |d|<1e-12.
The root is 1+h and delta=(1+h)h^2 f(1+h). Since f decreases,
delta/C<=1.001h^2. An exact polynomial identity is

    2s F0_t-3F0
      =s^2/t^3 [x p(t)+4-3t+d(8-7t)/t^2]
                           +delta(5t-2)/t^2,
    p(t)=t^3-2t^2-5t+12.

For 1<=t<=2, p decreases to 2, hence the bracket is at least
1998-15|d|>0. For t=2+y>=2,

    x p(t)+4-3t=x y^3+4x[y-(x+3)/(8x)]^2
                          +(31x^2-38x-9)/(16x)>x.

The d term has magnitude at most 15|d|. Therefore sF0_t/F0>=3/2
everywhere outside the bolt. Integrating
-(2s/F0)_t>=1/F0, whose boundary value at infinity is zero, proves

    F0(t) integral_t^infinity dz/F0(z)<=2s.             (7)

For s>=.001, f'<0 gives sFext_t/Fext<=2. Also Fext>=Cs^2/7;
the thermal term gives F0>=Cs^2/8 and
sF0_t-2F0<=delta(3t-1)/t^2. Hence sF0_t/F0<=2.0001.
If D=F0_t/2-F0/t is nonnegative, D integral 1/F0<=2.0001;
if negative, its absolute value times that integral is at most 2s/t.
These are difference bounds, so no two positive terms are added.
Substitute (7) into the exact tangent formulas to obtain

    |V|<2.2e-9, |u|<3e-9, |nu|<2.2e-9,
    |f/F0|<7.6e-10,
    b=a+f/F0+(sF0_t/F0-1)V, |b|<3.3e-9.                (8)

The existing invariant curvature identity then gives
|Delta u|<=2|u|+5e-8<5.6e-8. Keep the existing sharper outer
t>=4 bounds U<=8e-9/t^2, |V|,|nu|<=2e-9/t^2,
|b|<=1.4e-8/t^2 and the original weighted cold bounds. Equation
(8) is only needed on the middle patch. No bolt logarithm is
estimated by a false uniform bound.

## 4. Actual finite-endpoint source and geometry residual

Use the finite-endpoint profiles of
[the complete centered residual](stationary-centered-complete-residual.md),
not the larger parameter-derivative allowance. Before final mass/lapse
Jacobians, H_rho=r^2 rho2 and H_N=r^2(rho2+p2) have endpoint
coefficients T_c=1.1e17, T_m=3e18, T_o=2e17. Their differences
are bounded by 2T_c C^2 ell a_x on the cold patch, 2T_m C^2 ell
on the middle patch, and C^2 ell[2T_m+4T_o L(Y)^2] at t>=t_o,
where 40<t_o<60, a_x=s[1+log(.01/s)] and L(Y)<=L_e=1+log t.

The separate Ward metric contribution requires care. In its own hats,
H_rho=P-Theta and H_N=2P-Theta, so H_N-2H_rho=Theta exactly.
Every pure pressure memory or metric variation cancels from this
combination. Thus apply the finite trace/local profile envelopes to
this combination without adding the separate mass/lapse metric-kernel
bound. The local Weyl and reference pieces are already included in
T_c,T_m,T_o; the Ward identity is not imposed separately on them.
No mass/lapse Jacobian is appropriate to the orbit equation.

The exact normalized orbit source is
O_rest=(2G/(B^2C))(H_N-2H_rho). With GC/B^2<8.5e-33, triangle
bounds on the allocated profiles give

| Region | abs(Delta O_rest)/ell |
|---|---:|
| cold, s<=.001 | 1.2e-14 a_x |
| 1.001<=t<=40 | 3.1e-13 |
| t>=40 | 3.1e-13+4.1e-14 L_e^2 |

The last row covers 40..t_o by the middle bound. The budgets are
12(GC/B^2)T_c, 12(GC/B^2)T_m and 24(GC/B^2)T_o. Keep S160's
full local trace bounds: |Delta scriptT|/ell<=3.2e-16 a_x,
3.2e-17 and 3.2e-14 L_e^2 on these three regions.

Now use the exact geometric remainders N_O of S166 and N_S of
[S161](comparison-path-complete-equation-residual.md). No new
equation is obtained by differentiating a constraint primitive.
For explicit reproducibility, put z=U+N, fac=2U+Vabs+1.001N,
CP=10z^2+5z Vabs+2tau_geom Vabs^2. Their positive majorants are

    N_O/ell^2 <=1.01 fac |O_*|+2|Hpot| CP
      +2.02 D1 N fac+1.002 D1 tau_geom Vabs^2+1.004 D2 N^2
      +1.01 Vabs^2 t^2|L2 tau_geom|
      +2.02 t^3(F0/C)b^2/s+2.02t^2(F0/C)(Vabs/t+|b|)^2,
    N_S/ell^2 <=1.01 fac |A_*|+|P| CP
      +1.01 D1 N fac+.501 D1 tau_geom Vabs^2+.502 D2 N^2
      +1.01t Vabs |Delta u|+1.01Vabs^2 t^2|L tau_geom|
      +2.02t^3(F0/C)b^2/s+1.01t^2(F0/C)(Vabs/t+|b|)^2.

These display absolute upper majorants. On 1.001<=t<=4 use (8),
|Hpot|<=1.1t^4, |P|<=.51(t^4+1), D1<=2.1t^4, D2<=6.3t^4,
F0/C<=s^2, tau_geom<=1, |O_*|<=5.201e-9, |A_*|<=2.6e-9,
t^2|L2 tau_geom|<=6t^3 and t^2|L tau_geom|<=3t^3.
Every resulting monomial is increasing in t,s; evaluate at t=4,s=3.
This gives |N_O|<2.5e-13 ell^2 and |N_S|<7.1e-14 ell^2.
The cold bounds remain 3e-12 a_x, 2e-12 a_x, and the outer
t>=4 bounds remain 5e-15 L_e^2, 2e-15 L_e^2, all times ell^2.

The full normalized profiles satisfy the exact source cancellation

    R_V=N_O-ell Delta O_rest,
    R_Phi=N_O+2N_S-2ell Delta scriptT.

The resulting sufficient rows, divided by ell^2, are

| Region | abs(R_Phi)/ell^2 | abs(R_V)/ell^2 |
|---|---:|---:|
| cold | 7.1e-12 a_x | 3.1e-12 a_x |
| 1.001..4 | 4e-13 | 5.7e-13 |
| 4..40 | 9.1e-15 L_e^2 | 3.1e-13+5e-15 L_e^2 |
| 40..infinity | 7.4e-14 L_e^2 | 3.1e-13+4.6e-14 L_e^2 |

## 5. Physical forcing integral and projected composition

The exact fixed-volume conversion of S166 is

    f_ell=C/(8pi G B^2) exp(2ell u)t^-4 (R_Phi,R_V).

Its density factor cancels the moved radius exactly. The weighted
squared integral is gamma0 times integral t^-2(R_Phi^2+R_V^2)dt,
with exp(2ell u)<1.001. For I_n(T)=integral_T^infinity
t^-2(1+log t)^n dt, integration by parts gives

    I_n(T)=T^-1 sum_(j=0)^n binom(n,j)j!(1+log T)^(n-j).

Use I_0(4)=1/4, I_2(4)<3.12, I_4(4)<60;
I_0(40)=1/40, I_2(40)<.835, I_4(40)<33. The log bounds
log4<1.387 and log40<3.689 follow from positive exponential
series. On the cold patch use a_x<=.00331 and length<=.001;
on the middle patch extend its t^-2 integral to 1..4, giving 3/4.
Extend the positive 4..40 integrand to infinity and add the separate
40..infinity estimate. This harmless double counting gives

    1.001^2 integral t^-2(R_Phi^2+R_V^2)/ell^4 dt
                                                   <(8.5e-13)^2,
    D(f_ell)=lambda^-2 gamma0^-1/2 Fdisk(f_ell)
                   <(8.5e-13)(6.1/24)ell^2/g
                   <2.2e-13 ell^2/g.                    (9)

Here C/(lambda B^2)<6.1 and pi>3. These are the actual fixed-source
normalizations; no unit change is used to shrink a Newton product.
The exact projection Pi0 has norm one in Fdisk and uses the actual
adjoint Riesz vector W^-2 chi0. The source conversion is
Fsrc<=124lambda^-2 Fdisk. Equations (5),(6),(9) give respectively

    42000*124*2.2e-13<1.2e-6,
    8e5*124*2.2e-13<2.2e-5,
    4.224e6*124*2.2e-13<1.2e-4,

proving (1). The signed linear residual is still precisely the actual
adjoint direction times its unevaluated moment. Pi0 does not solve it.

## Limits and frozen verification contract

Even the improved weighted gradient bound only yields
||W gradbar tau||/sqrt(gamma0)<3.5e4 ell^2 after composition.
This upper bound does not prove that the gradient is large. The full
fourth/log and local cold divided-jet bounds have not been sharpened
to preserve a nonlinear C5 family; the unweighted metric value does
not imply pointwise control near the cold bolt. A mass upper bound
also does not establish the full exterior/charge return conditions.
No candidate Hessian inverse, nonlinear contraction, projected
nonlinear solution, signed compatibility, stability or physical
wormhole result follows. S167.1-S167.4 remain unchecked.

Freeze this note and its checker before execution. Verify predecessor
receipts and shared hashes, exact noncommuting elliptic regrouping on
a finite-dimensional fixture, rational coercivity/vertex/decay budgets,
the seed polynomial and completed square, low-jet estimates, exact
Ward cancellation, regional geometry/source sums, logarithmic moments,
physical-density conversion and final compositions. The fixture checks
algebra; it is not a physical state or a substitute for the continuum
form-domain argument. Preserve every failed source snapshot and run.
Continuum statements remain ordinary proofs; no Lean campaign is
selected. Do not rerun unchanged predecessor scientific suites.

Remaining to-do list: control the projected correction in the needed
local derivative norms; prove nonlinear family preservation and
remainder bounds; certify the projected full solve and actual signed
compatibility; exterior/state, charge and return; opening, stability,
transport, preparation, separate creation and full-result audit.
