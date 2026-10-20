# A numerical bound for the actual TT mass coefficient

2026-10-20 HST, E003/S166. Ordinary proof and contract before execution.
**The mass functional now has a numerical upper bound. The global
fast fourth/log remainder and complete S166 inverse remain open.**

Keep the full retained action, analytic seed, physical pairing, clock,
positive temperature, magnetic flux and paired domain. Put
u=B/r, u_h=B/r_h<1, W=u^-2, g=G lambda,
dmu=gamma0 u^-4 du, gamma0=4pi beta B^3, and Dbar=D/lambda^2.
The scalar sphere-component mass of the actual TT inverse satisfies

    |m_c(f)| <=6.6e4 g gamma0^(-1/2) ||Wf||.        (1)

For the actual compatible coupled solution, the
[effective source bound](strong-effective-tt-source.md) therefore gives

    |m_c| <=4.356e12 g gamma0^(-1/2) Fsrc
           <=5.40144e14 g lambda^-2 gamma0^(-1/2) Fdisk. (2)

The physical disk mass cubic has coefficient (3/4)m_c. This is an
upper estimate, not a value of that coefficient, a solution or a
stability result. It does not require a numerical adjoint solve.

## 1. An equivalent auxiliary probe with constant canonical variation

The [mass adjoint theorem](tt-mass-adjoint-certificate.md) gives
psi=1-Dbar_crit^-1 Dbar 1 and m_c(f)=<psi,f>/K, where

    gamma0/(6g)<K<gamma0/(5g),
    ||Dbar_crit^-1||_L2->L2<26g.                  (3)

For any p with p-1 in Dom(Dbar_crit),

    psi=p-Dbar_crit^-1 Dbar p.                    (4)

The correction cancels exactly by the critical inverse. Thus it
suffices to bound one convenient p and Dbar p; the earlier explicit
number Xi=||Dbar 1||/sqrt(gamma0) need not be evaluated.

Choose the fixed disk direction U0=(S,T)=(3/2,-3/2). Its canonical
orbit conformal variation T is constant, while S+T=0 leaves the
physical orbit metric unchanged. Decompose it by the genuine global
York map as

    U0=Rcal(tau,p).

This is an auxiliary linear test direction. It does not add physical
boundary data. The normal TT vector at root zero is exactly U0.

In physical proper radius write A=F_r/(2sqrt(F)), R=sqrt(F)/r,
w=A+2R, and let P be the radial conformal-vector operator. The
[seed projection theorem](stationary-global-tt-reduction.md) gives

    P Y=(3/2)(Y''+wY')-VY,
    V=F_rr/4+F_r/(2r)+3R^2+(3/2)A^2
      >=1.49lambda+3R^2+(3/2)A^2.                 (5)

Primes on Y denote proper derivatives. The projection of U0 solves
P Y=-6R. Introduce the explicitly regular comparison vector

    Y0=(r^3-r_h^3)/(r^2 sqrt(F)), div Y0=3,
    S0=F_rr+2F_r/r, J0=6R-S0 Y0.

Here S0 and J0 are local coefficients, not an action or norm. Direct
differentiation gives P Y0=-S0 Y0. Consequently

    Y=Y0+Z, PZ=-J0,
    tau=-(div Z)/4,
    p=1+2(r_h/r)^3+Z'/2-RZ+AZ/2.                 (6)

Y0=(3/2)rho+O(rho^3) at the regular bolt and approaches
1/sqrt(lambda) at infinity. J0 and Z are critical L2 fields. Solving
only for Z uses the positive critical inverse in (5); the constant
outer vector is never placed in that L2 domain.

The return vector is -Y: its ratio to sqrt(F) tends to zero at
infinity. Thus (6) uses precisely the fixed return in Rcal, with no
extra homogeneous solution. Smooth elliptic regularity and the
two-power outer decay put p-1 and tau in the critical graph domains.
Only domain membership, not an unspecified numerical regularity
constant, is used for this statement.

## 2. Quantitative projection bounds, including the bolt

Here are global point bounds that make (6) quantitative:

    |J0|<=300sqrt(lambda)u^2,
    |nabla(J0 partial_rho)|<=1e4lambda u^2.        (7)

The vector derivative is covariant; it includes J0', AJ0 and the
two sphere entries RJ0. To prove the bounds use the retained guards
F_rr>=lambda, F_rr<=21lambda, R^2<=5lambda and AR<=11lambda.
Convexity at r_h gives
F>=F_r(r_h)(r-r_h)+(lambda/2)(r-r_h)^2. Hence

    sqrt(lambda)|Y0|<=3sqrt(2)<4.3,
    RY0=1-(r_h/r)^3 in [0,1],
    AY0<=63, |Y0'|<=66.

The AY0 bound uses
F_r<=F_r(r_h)+21lambda(r-r_h) and
r[1-(r_h/r)^3]<=3(r-r_h). It stays finite at the bolt.

For an explicit outer calculation put

    nu=1/x, a=(1+2x+delta/2)/x,
    b=(1+3x+d)/x, e=d/x, v=r_h/B,
    f=F/(lambda r^2)
     =1+nu u^2-2a u^3+b u^4-2e u^5+e u^6.

The seed has nu<=.001, a<2.002, 0<b<3.002, |e|<1e-15
and 1<v<1.001. Exactly

    S0/lambda=6+2b u^4-12e u^5+12e u^6=:s(u),
    |S0|<13lambda, |grad S0|<55lambda^(3/2)u^4,
    J0/sqrt(lambda)=N(u)/sqrt(f),
    N=6f-s(1-v^3u^3).

N has powers 2 through 9 and absolute coefficient sum below 62.
On u<=1/2, f>.49, the coefficient sum of f-1 is below 7.01,
and |u f_u|<=43. Thus |J0|<89sqrt(lambda)u^2,
|J0'|<4000lambda u^2, |AJ0|<1500lambda u^2 and
|RJ0|<200lambda u^2.

Globally the elementary bounds above give
|J0|<70sqrt(lambda), |J0'|<1200lambda, |AJ0|<885lambda
and |RJ0|<157lambda. On u>=1/2 these imply (7).
In particular, with j=||J0/sqrt(lambda)||,

    j<=300sqrt(gamma0),
    ||nabla J0||/lambda<=1e4sqrt(gamma0).          (8)

The energy identity (5) yields

    ||sqrt(lambda)Z|| <=(100/149)j <202sqrt(gamma0),
    ||Z'||, ||AZ|| <=sqrt(200/447)j <201sqrt(gamma0),
    ||RZ|| <=sqrt(100/447)j <143sqrt(gamma0),
    ||nabla Z|| <=sqrt(600/447)j <350sqrt(gamma0).

All norms use the fixed physical measure. Equation (6) then gives

    ||W^-1 p||<400sqrt(gamma0),
    ||tau||<180sqrt(gamma0).                      (9)

For the explicit cubic use
||2(r_h/r)^3||=2sqrt(gamma0 v^3/3)<1.2sqrt(gamma0).
The weighted constant has norm sqrt(gamma0 u_h).

A radial vector is locally a gradient, so the exact identity
PZ=(3/2)grad(div Z)+2Ric(Z) gives

    grad tau=(J0 partial_rho+2Ric(Z))/6.

The retained bounds |Ric|<=14lambda and
|grad Ric|<200lambda^(3/2) now imply

    ||gradbar tau||<1000sqrt(gamma0),
    ||gradbar^2 tau||<2e4sqrt(gamma0).             (10)

For the latter the explicit budget is
1e4/6+(200*202+14*350)/3<2e4.
The derivative Ricci guard follows directly from Ric=J+(R4/4)g,
|grad J|<=100lambda^(3/2) and |grad R4|<lambda^(3/2);
it is not a derivative of a curvature value bound.

## 3. Full action on the constant disk direction

Let Abar be the actual fixed-density disk Hessian divided by lambda^2.
We prove

    ||Abar U0||<3g^-1 sqrt(gamma0).                (11)

Use the exact split I_full=I0+W4+Gamma_LLL+Gamma_mass-Gamma_ref.
Only I0 is stationary at its seed. No full-action defect is set to zero.

For I0, global U0 scales X=r^2 by exp(3t) and fixes the orbit metric.
Writing the seed potential as V0(X), stationarity reduces its mixed
variation to the potential

    P0(X)=3[X V0'(X)-V0(X)].

The fixed-density Hessian components are
-(P0+X P0',P0)/(8piG X). With the actual seed potential, they are

    (A_I0 U0)_S=[6-30d u^3+54d u^4]/(8piG r^2),
    (A_I0 U0)_T=[6-12(1+3x+d)u^2+60d u^3-54d u^4]
                                                    /(8piG r^2).

Their normalized coefficient sum is below 37/(8pi g), every power
is at least two, and their physical L2 norm is below
2g^-1 sqrt(gamma0). This uses the complete reference potential.

For the massive determinant, the exact action on a constant
canonical T is the constant-contact row, not a variable-metric
response estimate. Section 5 of
[the curved TT proof](stationary-curved-quantum-tt.md) gives

    |Theta_p0|<=1e20 q^2(1+1/Y), Y>=sqrt(5)u.

This includes the full tower, actual global resolvent, tadpole, finite
contacts and thermal terms. Theta_p0 is B^2 times the source per
canonical area. Since dvol_hat/dmu=B^2/(4pi r^4),

    (Abar_mass U0)_S=0,
    (Abar_mass U0)_T=(-3/2)Theta_p0 u^4/(4pi x^2).

Using q/x<2440 gives norm below 1e27sqrt(gamma0).
The constant canonical LLL row vanishes by its exact Dirichlet
Hessian: one of its two gradients is zero.

The Weyl and reference vacuum terms must be combined before taking
an L2 norm. Their common constant vector is

    lambda^-2 A_W U0 |_(u=0)
     =lambda^-2 A_ref U0 |_(u=0)
     =(11/(80pi^2),11/(160pi^2)).                  (12)

Here is an explicit way to verify the Weyl cancellation and bound its
remainder. Set omega=-log u, Q=q/x, and use dimensionless canonical
coefficients

    Rb=-2nu+12a u-12b u^2+40e u^3-30e u^4,
    Ub=u^-2-nu+4a u-3b u^2+8e u^3-5e u^4,
    Vb=u^-2+nu-2a u+b u^2-2e u^3+e u^4,
    kb=nu^2/15-Q^2/3,
    Lap p=partial_u(f partial_u p).

These are Rhat/lambda, Boxhat(omega)/lambda,
|gradhat omega|^2/lambda and k/lambda^2. In the exact disk Weyl
Hessian from section 8 of that proof define

    Dp=(3/2)(-Rb/15+nu/9)
          -(3/2)[(6omega+1)Rb/45+11Vb/135],
    Lp=-(3/2)[-Rb/15+2Ub/5+7Vb/45],
    Wb=-2Rb/15+14(Ub+Vb)/45,
    Zt=omega Rb^2/30+Rb(Ub+Vb)/15-Ub^2/5
                   -7Ub Vb/45-7Vb^2/90+2kb omega,
    Zst=-Rb^2/30+2kb,
    divW=-partial_u(f Wb/u).

Polarization and integration by parts give

    qS=Lap Lp+(3/2)divW-(3/2)Zst,
    qT=Lap Dp+(1/3)Lap Lp-3Zt+(3/2)Zst,
    Abar_W U0=-u^4(qS,qT)/(16pi^2).               (13)

The only constant terms of -u^4(qS,qT) are (11/5,11/10).
After their removal every term is u^k or u^k log u with k>=2.
Their combined absolute coefficient sum is below 1e8:
the basic sums are |Rb|_coef<61, |Ub|_coef<19,
|Vb|_coef<9, |f|_coef<9 and |kb|<2e6.
They give Lap Lp sum below 6480, Lap Dp below 11000, divW below
1400, and Zst and Zt sums below 4.1e6 and 4.11e6.
The L2 norms of u^k and u^k log u, divided by sqrt(gamma0),
are at most 1 and sqrt(2) for k>=2. Thus the Weyl remainder in
(12) has norm below 1e8sqrt(gamma0).

For the reference row, polarization of its exact disk Hessian gives

    A_ref U0=(-12pS_ref-3r partial_r pS_ref,-6pS_ref).

Its vacuum is p_v=-11lambda^2/(960pi^2), giving (12).
The retained explicit reference profile obeys

    |pS_ref-p_v|<3340(q^2/B^4)u^4,
    |r partial_r pS_ref|<1.4e4(q^2/B^4)u^4.

For example pS_ref-p_v=-(p0A-p_v)u^4+
D_P(u^4-6u^5+6u^6), with |p0A-p_v|<3301q^2/B^4
and |D_P|<3q^2/B^4. Its normalized row norm is below
1e12sqrt(gamma0). Combining (12), not bounding its two constant
vectors separately, proves (11), since g<2000^-10.

Boundary-supported seam contacts are retained with the original
realization. Their pairing here with compact test directions is zero;
the critical/fast extension has zero required boundary trace. The
constant probe does not impose a new finite-cut boundary problem.

## 4. The actual TT probe and numerical mass bound

The fixed-density block identity gives

    Dbar p=R_c^* Abar U0-Bbar* tau.                (14)

The unweighted return bounds ||S_c||<=10||c||,
||T_c||<=sqrt(500)||c|| imply ||R_c^*||<25.
For the mixed transpose use the local estimate already proved in
the strong effective-source note, now without weights. Its unweighted
lift and return adjoints have norms below 5, because
||h_c||<=sqrt(22)||c|| and ||sqrt(lambda)eta_c||<=sqrt(18)||c||.
Consequently, with epslambda=g/(60pi), vmag<4e-27/180,

    ||Bbar* tau|| <=(24g)^-1 {
       5[(7+20vmag+5e5 epslambda)||tau||
           +200epslambda||gradbar^2 tau||
           +8000epslambda||gradbar tau||]
                          +5(6e-22)||gradbar tau|| }
       <300g^-1 sqrt(gamma0).                    (15)

The off-shell return contact has not been omitted.
Equations (11),(14),(15) imply

    ||Dbar p||<400g^-1 sqrt(gamma0).

This is a numerical bound for an equivalent actual adjoint probe,
including the full quantum response. With (3),(4),(9),

    ||W^-1 psi|| <(400+26*400)sqrt(gamma0)
                    <11000sqrt(gamma0).

The exact mass functional and K bound prove (1). The strong
effective-source conversion proves (2). A numerical solve could
sharpen these deliberately conservative constants; none is assumed.

## Frozen verification contract

Freeze this proof and checker before execution. Revalidate the
receipts and exact current inputs for the seed projection, curved TT
form, scalar fourth graph, strong effective source and mass certificate.

Check the radial projection identities, regular-bolt leading term,
covariant divergence identity, seed source polynomials and inner/outer
budgets by exact algebra. Check the projection energy constants,
scalar derivative bounds, critical norm integrals, and all final
mass/source conversions with rational inequalities.

Independently differentiate the full local Weyl density under a constant
U0 and a variable test direction, integrate its test derivatives by the
Euler operator, and compare it with (13). Verify the exact constant
vacuum cancellation, the absence of powers below two in its remainder,
and the full polynomial/log coefficient envelope. Check the I0
potential identity and reference row separately. Preserve every run
and failed scientific source. No actual metric or adjoint fixture is
fitted, and no angular truncation is used.

The critical inverse, contour mass pairing, continuum domain extension
and previously proved full constant-contact bound remain ordinary
proofs. The verification does not certify the still-missing global
fast fourth/log inverse, nonlinear solution or physical stability.

Remaining to-do list: quantify the global fast TT remainder, higher
reconstruction and full inverse/residual composition; nonlinear family
and signed compatibility; exterior/state, charge and return; opening,
stability, transport, preparation, creation and full-result audit.
