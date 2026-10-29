# Complete local/scalar mass return: weighted actual-domain test

2026-10-29 HST, E003/S169. Ordinary argument and numerical contract BEFORE
execution. The preceding goal turn was NO PROGRESS (a status/plan review).
Preserve the closed audit, accepted results, P32, physical norms and full
stable semiclassical wormhole objective. This tests the remaining complete
local/scalar contribution in BOTH mass residuals; it is not a new physical
admission target. The massive/LLL contribution is reused from f0eab59.

## Actual problem and the first domain implication

Use Hmu=L2(z^-4 dz), E^2=ell/N, Es^2=a/N, N=|W_action| beta. In the
following pairings a,b,ell are normalized. The paired first-order mass/log
domain, complementary adjoint and squares, fixed clock/temperature/spin/
scale, common regulator, ultraviolet logarithm, complete angular/thermal
sums, coincidence contacts, finite returns and regular disk bolt are fixed.
No estimate here transfers the seed to a nonlinear metric neighborhood.

The compact canonical mass lift Phi=(nu,T)=-3 chi(n_m,b_m)/2 has the SAME
.1/.2 bump and exact primitive J_m=int_0^z t^3/Q(t)^2 dt as f0eab59.
Its original graph conversion cost 2133--2174 is unchanged. Its York
coordinates (s_m,h_m) obey the original return conditions. Put

    f_m=b h_m-a s_m, k_m=a^-1 f_m,
    a0=alpha A_r k_m-h(Yv0,Phi)+b(s_m,v0),
    a1=b(s_m,v0)+b(k_m,v0-v1)-h(Yv1,Phi).                 (1)

Here v0,v1 are the exact averaged rational trials already fixed, r=1e-110,
and Yv=(Uv,Vv) denotes the canonical image of the pure York shape. The
fixed-domain Ward identity makes the canonical determinant's scalar row
zero; it does NOT remove local anomaly or moving-density terms. Define

    Fnu=-z^4[H_c Phi]_nu, FT=-z^4[H_c Phi]_T,
    Abar k_m=Fnu,
    h_local(Yv,Phi)=(<Fnu,Uv>Hmu+<FT,Vv>Hmu)/beta.       (2)

H_c is the sum of the exact Einstein, Maxwell, R2 and Weyl operators in
signed_local_core.py, with original dz and negative W_action. In particular
it is not the locally cancelling reference operator. The endpoint terms
in its weak-to-strong integration vanish: Phi=O(z^3), the trial canonical
fields have admissible decaying outer rates, and all local products have
positive powers times at most logarithms; Phi vanishes near the bolt.
The C-infinity cutoff is retained through four derivatives. No separate
matter boundary contact is discarded; those are in the accepted massive
and LLL block. Equations (1)-(2) are identities before taking norms.

The first new implication is a WEIGHTED inverse for the actual scalar
source, not qualitative elliptic smoothing. Let Lop=-Delta_bar on the
original radial form realization, Delta_bar=z^2[(Q s')'-2Q s'/z]. Its
nonnegative form is int Q z^-2 |s'|^2. The scalar form domain is D(Lop)
with its graph norm: the positive fourth-order coefficient, bounded
lower-order coefficients and original zero-source core establish this
equivalence. Truncations at the outer end and regular-disk cores at the
bolt give the same closure; no additional fourth-order boundary equation
is introduced by factoring the operator.

Choose p=3/4 and W=[(L+4gamma)/(L-z+4gamma)]^p. W and W^-1 are bounded
for this FIXED gamma, preserve the original outer rates and are smooth
in regular disk radius. Multiplication preserves the scalar graph domain:
the product rule uses bounded W, z sqrt(Q) W', and Delta_bar W. The
same property follows for the York Dirichlet form by its potential bound.
These are auxiliary weights, not new allowed metric variations.

Write sigma2=1.001 p^2. The actual bolt polynomial has
Q<=2gamma(L-z)+1.001(L-z)^2, whence Q(log W)'^2<=sigma2 and also
z^2 Q(log W)'^2<=sigma2. The checker verifies the polynomial signs on
the entire interval, including Q Py>=1 and Q Pxi>=1.9 below.

## Weighted scalar inverse, with the complete remainder

The accepted exact scalar row is

    Abar=6 D P+Arem, P=Lop+4beta, D=Cv+d Lop,
    d=96pi cR C/B2-epsbar, epsbar=C/(60pi B2), Cv=1+11epsbar beta/3.

Arem has gradient coefficient bounded by
KA=(2d+11epsbar/6)w0+22epsbar tau_max/3 and mass coefficient by 2w0,
where w0=1.1e-27 beta, tau_max=3.01beta. These are the COMPLETE retained
off-shell/anomaly coefficients, not a constant-curvature replacement.

Real conjugated second-order forms for W and W^-1 give, with
cp=4beta-sigma2>0 and cd=Cv-d sigma2>0,

    ||W P^-1 W^-1||<=1/cp, ||W D^-1 W^-1||<=1/cd,
    GP=1/sqrt(cp)+sqrt(sigma2)/cp,
    GD=1/sqrt(d cd)+sqrt(sigma2)/cd.

GP,GD bound the ORIGINAL weighted gradients of the two inverse outputs;
the derivatives of W have already been paid. Dual testing with the
opposite weight gives

    R_A=D^-1 Arem P^-1/6,
    ||W R_A W^-1||<=rA=[KA GD GP+2w0/(cd cp)]/6.

If rA<1, the bounded inverse of I+R_A and the original form uniqueness
give k=P^-1(I+R_A)^-1 D^-1 Fnu/6. This factorization holds first as a
form identity and then by closure; it does not assume D and Arem commute.
For Fn=||W Fnu||Hmu and Yc=1/[6cd(1-rA)], put

    Nk0=Yc/cp, NkG=Yc GP, NkL=Yc(1+4beta/cp).

Then weighted L2, gradient and Laplacian norms of k are respectively
<=Nk0 Fn, NkG Fn, NkL Fn. Also Es(k)<=Fn/(beta sqrt(23.99)).

Let b_j be the actual bolt polynomial coefficients, j=2,...,6. An explicit
bound for |Delta_bar log W| is

    Do=p[L^2*1.001+L^2(1/2+sum j|b_j|L^(j-2))
                              +2L*1.001(L+4gamma)].

For u=Wk, define Gu=NkG+sqrt(sigma2)Nk0,
Lu=NkL+2sqrt(sigma2)NkG+(Do+sigma2)Nk0,
Do2=Lu+2sqrt(.17)Gu. Then ||(Q u')'||L2(dz)<=Do2 Fn.
The accepted Q>=gamma(L-z)/2 and regular zero flux imply
|u(L)-u(L-s)|<=4sqrt(s)Do2 Fn/gamma. Average over length
b=gamma Nk0/(4Do2)<L, using ||u||L2<=Nk0 Fn, to obtain

    |A_r k| <= Kobs Fn,
    Kobs=[(4gamma+r)/(L+4gamma)]^p
                 [4sqrt(Nk0 Do2/gamma)+4sqrt(r)Do2/gamma].      (3)

The physical alpha multiplies (3) in (1); it is never absorbed in W.

## Global York estimates and actual mixed pairings

The inverse York for Phi is

    -y''+Py y=S, Py=6/z^2+2(Q''-4Q'/z)/(3Q),
    S=(2/3)(T'/z^2-4T/z^3), y(0)=y(L)=0,
    s_m=u0+z^2 y'/4-z y/2, u0=nu+T/2.

The original mass branch gives y=z^2/4+O(z^3); its leading contribution
to s_m cancels, so s_m has the scalar admissible outer rate. The complete
return is regular at the bolt. With Py>=5/z^2, QPy>=1, ky=1-sigma2,
TH=||WT/z^2||2, weighted energy testing and one integration of T' give

    Ey=(2/3)TH[1+(2+.2p/(L-.2+4gamma))/sqrt(5ky)],
    Yz=Ey/sqrt(5ky), Yp=(1+sqrt(sigma2/ky))Ey, YP=Ey/sqrt(ky).

These bound ||Wy/z||, ||Wy'||, ||W sqrt(Py)y|| respectively. Calculate
the three local quantities U0H=||Wu0/z^2||2,
Gloc=||W sqrt(Q)(u0'-z^2S/4)/z||2,
Lloc=||W[Q u0''+(Q'-2Q/z)u0'-z^2(QS)'/4]||2. Then

    Ns0=U0H+Yp/4+Yz/2,
    NsG=Gloc+sqrt(3.22)YP/4+sqrt(.17)Yz/2,
    NsL=Lloc+3.59Yz+.72Yp.

The last equation follows by differentiating the displayed s_m and using
|z Cy|<=3.59, |By|<=.72, where
Cy=-2Q/z+5Q'/3-2zQ''/3+z^2Q'''/6,
By=Q+z^2Q''/6-2zQ'/3; z^2QPy<=3.22. All are whole-interval guards.
Thus

    Es(s_m)^2<=[6d NsL^2+(6Cv+24beta d+KA)NsG^2
                                      +(24beta Cv+2w0)Ns0^2]/beta. (4)

For forward York on a trial xi, use Pxi=6/z^2+(Q''-4Q'/z)/Q,
-y''+Pxi y=(4xi/z-xi')/z^2, Uxi=xi/2+zy,
Vxi=-xi-zy-z^2y'/2. With kx=1-sigma2/1.9, Pxi>=4.5/z^2,
Q Pxi>=1.9, the SAME weighted test with W^-1 gives

    cY=1+[sqrt(sigma2/1.9)+2/sqrt(4.5)]/sqrt(kx),
    yz=cY/sqrt(4.5kx), yp=cY[1+sqrt((sigma2/1.9)/kx)],
    CUw=1/2+yz, CVw=1+yz+yp/2.

Hence original Hmu norms of W^-1 Uv,W^-1 Vv are bounded by CUw,CVw
times ||W^-1v||Hmu. No finite-domain return is dropped.

The exact mixed adjoint is bbar* s=4 tau_bar D s+brem* s. First test the
COMPLETE refined unweighted constant Kbrem in the accepted scalar-source
receipt: ||brem* s||Hmu<=Kbrem Es(s). Its five terms include density,
curvature, magnetic and global return contacts. A new weighted estimate
for that remainder is NOT assumed. If Vwk=||W^-1 vk||Hmu and V0k=||vk||Hmu,

    Bs0=4tau_max(Cv Ns0+d NsL)Vw0/beta
                                      +Kbrem Es(s_m)V00/beta,
    Bk01=4tau_max(Cv Nk0+d NkL)Fn(Vw0+Vw1)/beta
                                      +Kbrem Es(k)(V00+V01)/beta,
    Hloc=(Fn Knu+Ft KT)(q Nfw0+Nfw1)/beta,
    Elocal=Kg[q alpha Kobs Fn+(q+1)Bs0+Bk01+Hloc].          (5)

Here Kg is the accepted norm of the ACTUAL filtered mass functional; its
adjoint error is already included. q=.338. The two equations are both
charged. No unknown source/boundary/graph comparison is set to zero.

Test one specified sharper estimate in the SAME account. The previous
generic Kbrem uses Es to bound the scalar Hessian, losing d^-1/2. Here
the exact scalar equation and York equations have already bounded that
Laplacian directly. For actual bounds N0,NG,NL set

    NH=sqrt(NL^2+13beta NG^2),
    Rjets=3e6[4400+2*4*2sqrt(2)*20+8], Cjets=100Rjets,
    Brem(N0,NG,NL)=4tau_max(11epsbar w0/36)N0+sqrt(22)d w2 N0
       +epsbar sqrt(22)[100Cjets(NH+2NG+N0)+1e6 N0]
       +[20/(60pi)]4e-27 beta N0+(9e-28 beta)(5/sqrt(beta))NG,
    w2=4.3e-25 beta^2.

This is exactly the accepted tensor/Bochner/mixed-return estimate before
substitution of generic energy-to-derivative constants. The Ricci bound,
core approximation and vanishing endpoint terms justify the unweighted
Bochner estimate NH; W>=1 supplies its unweighted inputs. Evaluate Brem
on (Ns0,NsG,NsL) and Fn(Nk0,NkG,NkL), replacing ONLY Kbrem Es in (5).
The global return remains charged. Report both complete estimators and
their ratio, so this refinement is justified by its full-budget effect.

## Complete trial norms and whole-support arithmetic

For the auxiliary H=-( (2x+a x^2) f')', let Rmu=(H+mu)^-1. This operator
only DEFINES the exact trial; it replaces no physical matter spectrum.
On weight (x+4)^p, real forms give inverse norm <=1/(mu-sigma2).
Truncated weights followed by monotone convergence justify the unbounded
auxiliary weight; the constants are uniform in the truncation.

For nu(nu+1)=muT/a, the exact Green kernel is
Pnu(1+a min(x,y))Qnu(1+a max(x,y)); its flux Wronskian is -1.
The prescribed positive unit density rho is supported in [0,r/gamma].
The positive P solution obeys P<=exp(muT x/2), from its Volterra equation.
Verify exp(muT r/(2gamma))<2, so |RmuT rho|<=2Qnu(1+a x) globally.
The accepted positive-integral Q inequalities give

    BQ^2=(4+1/a)^(2p)[log(3)^2+2log(3)+2]/(4a)
           +pi^2(1+4a)^(2p)/[a^(2p+1)(2nu+1-2p)].

Thus 2BQ bounds the complete weighted first resolvent. For the exact
nominal rational parameters, A0=3alpha/(4gamma_nom A_nom), t=2a_nom,

    Kfw=2BQ A0[1+|muY_nom-t|/(t-sigma2)]
                      [1+2c/(muS-sigma2)],
    h=d1/2, uW=|h^2-d0|/(h-sigma2)^2<1,
    LFw=kappa[(Cv_nom-delta h)/(h-sigma2)^2+delta/(h-sigma2)]/(1-uW),
    Nfw0=gamma^(p+1/2)Kfw/(L+4gamma)^p, Nfw1=LFw Nfw0.

These retain the entire rational functions and auxiliary spectrum. The
Nfw values bound the physical-dz weighted PROFILE, not the physical TT
field. The distinction is essential. Attempt1 omitted the reconstruction
in these two trial norms and is REJECTED despite its numerical pass.
Its immutable receipt preserves that failed estimator implementation.

Retain the exact reconstruction (5) of current-source-adapted-local-return:
r=(H+muY_actual)^-1 f, R=Q/Q2, M=zQ''-4Q', and

    nu=z^4/(3L^5)[-zR f+(zR muY-M)r-(7R-zR_s)Q2 r_s+Wret],
    T=z^4/L^5[zR f+(2M/3-zR muY)r+(5R-zR_s)Q2 r_s-2Wret/3],
    Wret(z)=int_0^z -[3Q''-tQ'''] r(t)dt.

Use 0<R<2, |R_s|<=R1=sum_(j=3)^6 (j-1)|b_j|L^(j-3)/a,
M1=sqrt(2gamma L+aL^2), Dmax=6.3+8.1L, Mmax=2.1L+1.2.
For s0=sigma2 (weighted) or 0 (unweighted), let
gap=muY_actual-s0 and Gr=1/sqrt(gap)+sqrt(s0)/gap. Because W^-1 decreases
with z, the Volterra weight ratio W(t)/W(z)<=1 for t<=z. Both Schur
rows of its weighted kernel are bounded by L Dmax. Consequently

    Knu=[2L+(2L muY+Mmax)/gap+(14+L R1)M1 Gr+L Dmax/gap]/(3L^3),
    KT=[2L+(2Mmax/3+2L muY)/gap+(10+L R1)M1 Gr
                                                    +2L Dmax/(3gap)]/L^3.

These bound the ACTUAL canonical Hmu outputs by the profile's L2(dz)
norm. For the actual inverse V^-1, weighted testing of the Py equation
gives cy=(2/3)[1+(sqrt(s0)+2/sqrt(5))/sqrt(1-s0)] and

    Kinv=1+cy/sqrt(5(1-s0))+cy[1+sqrt(s0/(1-s0))]/2,
    Kxi=Kinv KT.

Use Vwk=Kxi(sigma2)Nfwk and V0k=Kxi(0)Lk, with Lk the accepted COMPLETE
unweighted profile norms. Use Knu,KT directly in Hloc to retain both
explicit canonical rows. These comparisons include the finite Volterra
return, not a flat multiplication (z/L)^4. Global geometry guards are
verified by whole-interval polynomial Bernstein enclosures. Neither the
auxiliary weight nor this reconstruction changes any physical norm.

For exact local rows, factor Phi=z^3(psi_nu,psi_T), with
I=J_m/z^4=int_0^1 t^3/Q(zt)^2 dt,
psi_nu=3beta chi Q I, psi_T=-3beta chi(1/Q+z Q'I)/2.
Taylor coefficients of I are bounded by those of Q^-2 on [0,z_hi],
divided by j+4. This includes z=0 without a singular division.
Assemble all four physical allocations and both directions BEFORE
absolute values. Extract exact polynomial plus log(z) coefficients and
integrate the squared whole-cell bounds with exact moments of
z^2, -z^2 log(z), z^2 log(z)^2. The .1/.2 transition uses the frozen
1024-cell bump enclosure; 512 inner cells cover the remaining support.
The local quantities in (4) are integrated on those same complete cells.
Outward interval arithmetic encloses rounding and coefficient uncertainty.
There is no discretized inverse or omitted scalar application error.

## Complete diagnostic contract and limits

Use the accepted Eshape<2.063278e33 and bridge, with their full precision:

    Emass=Eshape+Elocal,
    Ebridge+sqrt(1+epsilon)[Emass+Efast+Erest+Eapplication]<=1.49e34.

Execute (2)-(5), check the domain/positivity guards, report every debit
and the conditional remaining account. If the resulting estimator is
excessive, quantify its overrun and dominant term before further work.
Efast,Erest,Eapplication remain UNKNOWN. Even a finite affordable Emass
only completes the MASS contribution to this diagnostic. A complete Dj
pass would reject the present signed-trial bolt estimator, as the accepted
ledger explains; it is not physical admission or construction obstruction.
The signed coefficients themselves may remain unenclosed more tightly
than their symmetric absolute intervals.

Remaining to-do list: complete mass-free residual and application errors;
decide the full diagnostic and estimation method; signed Gt, original
coupled/admission budgets, nonlinear existence/contraction, compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
