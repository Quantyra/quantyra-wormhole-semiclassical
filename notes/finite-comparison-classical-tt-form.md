# Classical TT form throughout the actual finite comparison

2026-10-20, E003/S167, supporting S139/S074. Ordinary argument completed
before the frozen checks. **Component progress; the full goal is incomplete.**

For the retained fixed action I0=S_classical+Gamma_ref, the tensor
exponential TT Hessian throughout the actual finite path obeys

    Q_exp,ell[h_ell] >= .99 E_ell[c].                         (1)

Here E_ell is the candidate's own TT energy, defined below. This theorem
does not assume the finite comparison is stationary. It does not yet
bound the disk-chart Hessian: the exact off-shell contact in section 7
remains necessary. The canonical LLL disk form also has a lower bound
against E_ell, but forms in different charts cannot simply be added.

## 1. Fixed data and a sharper first derivative

Keep the theory, potential, magnetic flux, temperature, outer clock,
state, full paired tower, renormalization scale, regular bolt and outer
domain of the [centered seed](stationary-centered-covariant-seed.md).
Use the actual path of the [finite scalar comparison](finite-comparison-scalar-resolvent.md),
with lambda=lambda_e, x=lambda B^2, v=r_h/B=1+h and t=r0/B:

    sigma_ell=sigma0+ell u, r_ell=B+(r0-B)exp(ell V), 0<=ell<=1,
    tau=(t-1)/t, nu=tau V, T=1+tau(exp(ell V)-1), S=log T,
    w=ell u+S, g_ell=e^(2w)(d rho^2+a0^2 dtheta^2)+r0^2 e^(2S)dOmega^2,
    y=sqrt(lambda)rho, d rho_ell=e^w d rho,
    dmu_ell=e^(2w+2S)dmu0, |w|,|S|<eps=3e-8.

A0=(log a0)_y, R0=(log r0)_y are dimensionless. The retained seed
bounds are R0^2<=5, 0<=A0 R0<=11 and
.5<=a0,yy/a0=F0,rr/(2lambda)<=10.5. Also
F0>=.95 x(t-v)^2, 5.9<C/x<6.1. The actual linear response has
|u|,|V|,|nu|<1.2e-8 and |b_gap|=|(t-1)V_t|<3e-7.

The [completed metric equation](stationary-completed-metric-response-norm.md)
gives |A_phys^2 delta Rhat_common|<5.2e-8 and
|A_phys^2 Rhat0|<=2, where A_phys=B/sqrt(C). In the seed's
dimensionless canonical proper radius zeta it says exactly

    u_zetazeta+(N_zeta/N)u_zeta=F_u,
    F_u=-A_phys^2 Rhat0 u-A_phys^2 delta Rhat_common/2, |F_u|<5e-8,
    d zeta=sqrt(C) dt/[t sqrt(F0)],
    N=sqrt(C)sqrt(F0)/(kappa B t),
    N d zeta=C dt/(kappa B t^2).

Regularity selects the zero integration constant at the bolt. Thus

    (integral_0^zeta N)/N=sqrt(C)(t-v)/[v sqrt(F0)],
    |u_zeta|<5e-8 sqrt(6.1/.95)<1.27e-7,
    |u_y|<5e-8 (C/x)/[sqrt(.95)vt]<3.2e-7/t,
    u_yy+A0 u_y=C F_u/(x t^2),
    |u_yy+A0 u_y|<3.05e-7/t^2.                       (2)

The last identity uses d zeta/dy=sqrt(C/x)/t and
(log N)_y=A0-R0. It bounds the combination actually present in the
form without requesting a derivative of the source.

For the radius, nu_y=R0(b_gap+V/t), giving
|nu_y|<7e-7 and |A0 nu_y|<3.432e-6. With k=exp(ell V),
B_ell=k-1-ell V k, the exact identity is
T_y=ell k nu_y+B_ell tau_y. The retained exponential remainder
has |B_ell|<2e-16, k<1.001, T>.999, and A0 tau_y=A0 R0/t<=11.
Consequently

    |S_y|<7.1e-7, |A0 S_y|<3.45e-6, |w_y|<1.03e-6,
    |p|=|w_y+S_y|<1.74e-6, |q|=|w_y-S_y|<3.2e-7.   (3)

## 2. Compare genuine TT lifts with the same scalar amplitude

Component order throughout is (proper radial,time,sphere,sphere).
Write h_ell=diag(b_ell,-b_ell-2c,c,c), using the same c(y) on each
geometry. Transversality in the fixed coordinate is

    b_ell,y+2(A0+R0+p)b_ell=-2(A0-R0+q)c.

The [seed lift](stationary-global-tt-reduction.md) has ||b0||<=2||c||.
All norms in this section use dmu0. Set D0=partial_y+2(A0+R0).
Its regular inverse has norm <=sqrt(8)/3<.943. Here are the domain
and range details needed to use this bound. Comparison of a0 with
sinh(y/sqrt(2)) by their Wronskian gives
A0>=coth(y/sqrt(2))/sqrt(2)>1/sqrt(2). On the regular core,

    Re<v,D0v>=integral (3A0/2+R0)|v|^2
                  >=3/(2sqrt(2))||v||^2.

The homogeneous solution (a0^2 r0^2)^-1 is not L2 at the bolt;
the adjoint D0*=-partial_y+A0 has homogeneous solution a0, which
is not L2 at infinity. The closed operator therefore has closed,
dense range, hence is onto, with the stated inverse bound. Its
energy identity is

    ||D0v||^2=||v_y||^2
             +integral [4A0^2+2R0^2-F0,rr/lambda]|v|^2,
    ||v_y||^2+4||A0v||^2<20||D0v||^2.                (4)

The last line uses F0,rr/lambda<=21 and ||v||^2<=8||D0v||^2/9.
Initially take c compact away from the endpoints and the regular
TT lift, with its allowed outer tail. Subtracting the TT equations,

    d_b=b_ell-b0, D0 d_b=f=-2p b_ell-2q c.

A convergent bounded perturbation of D0 gives

    ||b_ell||<2.001||c||,
    ||f||<7.604e-6||c||, ||d_b||<7.2e-6||c||.          (5)

Only d_b, rather than b0 separately, needs to be in Dom D0. For a
regular nonzero bolt value b0(0)=-c(0), the individual term A0 b0
need not be L2. The difference is regular and vanishes to second
order there. Bounds (4)-(5) extend by the established covariant H1
TT closure. No physical bolt datum is set to zero.

## 3. Exact nonstationary exponential-chart form

Keep the fixed reference potential, not a potential recentered at ell:

    V0(r^2)=2+6lambda r^2-2q0/r^2
             +2d(-B^2/r^2+4B^3/r^3-3B^4/r^4),
    q0=B^2(1+3x), x>=1000, |d|<=1.

For ds^2=N^2 d rho^2+a^2 dtheta^2+r^2 dOmega^2, the radial
action after integration by parts is

    I0=-(pi/G) integral [a r'^2/N+2a' r r'/N+(aN/2)V0(r^2)]d rho
          + the fixed-boundary and bolt terms.                 (6)

This is the retained covariant spherical action; equivalently its
areal form is beta/(2G) integral e^alpha[(rF)'-V0/2]dr,
with beta kappa=2pi. Start at N=1. The tensor exponential chart is

    N_e=e^(e b/2), a_e=a e^(-e(b+2c)/2), r_e=r e^(e c/2).

It preserves the four-volume exactly. Expanding (6), its second
variation density, in units dmu/(8pi G), is

    1.5 c'^2+b'c'+2b(A-R)c'-2Rbb'
                -b^2(R^2+2AR)-P c^2,
    P=(V0-X V0'+X^2 V0'')/(2X), X=r^2,              (7)

where A=a'/a, R=r'/r and primes on V0 mean d/dX. Substitute
b'=-2(A+R)b-2(A-R)c and integrate the c c' term. This yields

    Q_exp=1/(8pi G) integral [1.5 c'^2-4Rbc'
       +(E1+2AR-2R^2-P)c^2+(2AR+3R^2)b^2
       +4(AR-R^2)bc]dmu,
    E1=A'-R'+A(A-R).                                (8)

No field equation or stationarity is used in (7)-(8). All coefficients
in (8) are bounded at the bolt. On the endpoint-vanishing core the
integration is direct. For a regular nonzero bolt value the surface
term from (7) is +(pi/G)r_h^2 c(0)^2. The bolt action is
-pi r_h^2/G; its second variation along this chart is
-(pi/G)r_h^2 c(0)^2, and cancels that surface term. The fixed period
is preserved because regular TT data have b(0)=-c(0), so the radial
and time exponential factors coincide at the bolt. At infinity the
permitted TT tail gives a zero surface term. Thus (8) is also the
closed covariant form on regular H1 data, including nonzero bolt value.

Direct differentiation gives

    P=1/r^2-4q0/r^4+d(-4B^2/r^4+25B^3/r^5-27B^4/r^6),
    |P|/lambda<=12+61/x<13,
    |d_logr P|/lambda<=48+321/x<49.                 (9)

The interpolation in r stays above B, since the actual path has
r_ell=B+(r0-B)e^(ell V)>B. At the seed, one can check the
normalization independently in areal coordinates:

    P=F/r^2-F_r/(2r)+r F_rrr/4,
    W_F=3F-rF_r/2+r^2 F_rr/2-r^3 F_rrr/4,
    q_old=1.5r^2 F(c_r+2c/r)^2+W_F c^2-4Fbc,
    r^2 q_exp-q_old=d_r[-rF(b+c)(b+3c)].            (10)

Here q_exp denotes the bracket in (8) evaluated at the seed, and
the seed TT constraint is substituted in the derivative on the
right. Its boundary term vanishes for regular data and the allowed
tail. Therefore the retained seed theorem gives Q_exp,0>=E0[c],
with its original physical normalization.

## 4. Uniform comparison of the forms

Pull (8) back to y, divide by lambda/(8pi G), and denote the
resulting integrand e^(2S) q_ell. Its algebraic form is that of (8),
with c' replaced by c_y and with coefficients

    Rtilde=R0+S_y, ARtilde=(A0+w_y)(R0+S_y),
    jtilde=e^(2w) E1ell/lambda=j0+ell(u_yy+A0 u_y),
    Ptilde=e^(2w) P(r_ell)/lambda,
    j0=A0,y-R0,y+A0(A0-R0).

The exact E1 identity was established in the
[return proof](finite-comparison-disk-return-form.md). Equations (2)-(3),(9)
give the following absolute changes, in the displayed order:

    |Delta Rtilde|<7.1e-7, |Delta Rtilde^2|<3.2e-6,
    |Delta ARtilde|<6e-6, |Delta jtilde|<3.05e-7,
    |Delta Ptilde|<2.3e-6.                          (11)

For example the AR change is at most
3.45e-6+2.24(1.03e-6)+(1.03e-6)(7.1e-7)<6e-6.
Use e^(2eps)-1<6.1e-8 for the potential and measure changes.
The c^2, b^2, bc coefficient changes in (8) are respectively less
than 2.2e-5, 2.2e-5, 3.7e-5. Their seed absolute bounds are
84,37,44, using |j0|<39, A0 R0<=11, R0^2<=5 and (9).

Put a=||c_y||, Cc=||c||, D=||c_y+2R0c|| and
e0=1.5D^2+Cc^2. Then a<=D+4.48Cc. Using (5), the absolute
integrand difference satisfies

    integral |q_ell-q0| <= L a Cc+M Cc^2,
    L=4[(7.1e-7)(2.001)+2.24(7.2e-6)]<7.03e-5,
    M=2.2e-5(1+2.001^2)+3.7e-5(2.001)
        +37(4.001)(7.2e-6)+44(7.2e-6)<.00157,
    integral |q_ell-q0|<.002 e0.                   (12)

Use D Cc<=(D^2+Cc^2)/2 for the last line. Independently,

    integral |q0|<=1.5a^2+17.92a Cc+320Cc^2
                  <=11.96D^2+469.4528Cc^2<500e0.

It follows that integral |q_ell|<501e0. Restoring e^(2S),
the total form difference is less than
(.002+501(6.1e-8))E0<.003E0. Hence

    Q_exp,ell >= .997 E0[c].                       (13)

This comparison keeps both the actual TT lift and the measure.
It does not replace h_ell with a seed TT tensor.

## 5. Candidate energy and a common LLL comparison scale

Define

    E_ell[c]=1/(8pi G) integral
       [1.5(c'_ell+2Rell c)^2+lambda c^2]dmu_ell
      =lambda/(8pi G) integral e^(2S)
       [1.5(c_y+2(R0+S_y)c)^2+e^(2w)c^2]dmu0.

The inequality (a+b)^2<=1.001a^2+1001b^2 with
|2S_y|<1.42e-6, in both directions, gives

    E_ell<=1.002 E0, E0<=1.002 E_ell.               (14)

For the reverse direction use c_y+2R0c=(c_y+2Rtilde c)-2S_y c;
the derivative coefficient is at most 1.001e^(2eps), and the
value coefficient at most e^(4eps)[1+1501.5(1.42e-6)^2].
Both are below 1.002. Since .997/1.002>.99, (13)-(14) prove (1).

For comparison with the previous LLL result, set
E_h=||h'_ell||_ell^2+lambda||h_ell||_ell^2. The
[seed gradient bound](stationary-seed-quantum-form.md) is
||h0,y||^2<=5200(D^2+Cc^2). Equation (4) gives
||h_ell,y-h0,y||^2<40(7.604e-6)^2 Cc^2. Therefore
||h_ell,y||<73 sqrt(D^2+Cc^2); also
||h_ell||_0^2<23Cc^2 by (5). Accounting for the derivative
and value measure factors yields

    E_h<5400lambda(D^2+Cc^2), E_h/(8pi G)<5500E_ell,
    Gamma_hat,LLL,disk''>=-1.3e-27 E_h/(8pi G)>-8e-24 E_ell. (15)

This merely puts the LLL disk bound on the same energy scale as (1).
The chart difference below must still be controlled before combining them.

## 6. Verification scope

The frozen checker pins the preceding return, seed TT and seed quantum
receipts and their input union. It verifies the canonical Volterra
moment, coordinate cancellation, TT lift and inverse identities,
the nonstationary exponential action expansion, integration-by-parts
and seed consistency identities, the explicit fixed potential,
the bolt cancellation and all rational margins in (1)-(15).
These are exact algebra and conservative budgets supporting this
ordinary continuum argument; they are not formal verification or
a numerical solution of the full equations. Every scientific run,
including any failure and its source snapshot, is retained.

## 7. Exact remaining contact and next decision

Let h be the genuine TT tensor, Y=eta partial_rho the actual
outer-fixed return, and k=h-L_Yg the disk tangent. Let H_aff be
the affine tensor Hessian. Diffeomorphism invariance gives the
[retained identity](stationary-seed-quantum-form.md)

    H_aff(k,k)-H_aff(h,h)=I0'[L_Y(2h-L_Yg)].

For a covariant tensor v, write v^2_ab=v_ac g^cd v_db. The tensor
and disk exponential accelerations are respectively h^2 and k^2.
Thus, exactly at every finite comparison,

    I0,disk''[k,k]-Q_exp[h]
       =I0'[L_Y(2h-L_Yg)+k^2-h^2].                 (16)

Only at the seed is I0'=0. This turn supplies no uniform numerical
bound for (16). The classical exponential form is positive; the
candidate classical disk form, full quantum TT block, mixed/density
contacts, full coupled inverse and projected nonlinear solve remain
unproved. Pi_cut remains the fixed seed projection. All S167 tasks
stay unchecked and no physical gate is promoted.

Remaining to-do list: bound (16) and the remaining massive/Weyl/reference
TT and mixed/density terms; complete logarithmic domains and candidate
range; bound cold/higher correction and source jets with the fixed
projection; preserve the nonlinear family and solve its projected
equations; evaluate signed compatibility; establish exterior/state,
charge/return, opening, stability, transport, preparation and creation.
