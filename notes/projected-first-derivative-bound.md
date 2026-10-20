# First derivative bounds for the actual projected correction

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **Small first-order integral
bounds; local nonlinear-family preservation remains unproved.**

Keep the original retained action, analytic seed, positive temperature,
clock, flux, state, paired domain and physical mass column. Use the
actual comparison residual f_ell and exact compatible projection Pi0
of [the value and mass theorem](projected-low-correction-and-mass.md).
Put h_ell=-A_m^-1 Pi0 f_ell=Rcal(tau,c)=(S,T), with 0<=ell<=1.
All norms below use dmu=4pi beta r^2 dr, gamma0=4pi beta B^3,
g=G lambda, W=(r/B)^2 and gradbar=lambda^-1/2 grad. Then

    ||W gradbar tau||/sqrt(gamma0)<2e-4 ell^2,
    ||gradbar(S,T)||/sqrt(gamma0)<1.3e-4 ell^2,
    [||(S,T)||^2+||gradbar(S,T)||^2]^(1/2)/sqrt(gamma0)
                                                  <1.31e-4 ell^2. (1)

The gradient in the last two lines differentiates the two disk scalar
coordinates, not a discarded-mass remainder. It is unweighted and
includes the actual mass and fixed return tail. The scalar estimate
is fast weighted. The scale lambda and physical volume are unchanged.
These are bounds on the seed-preconditioned correction, not values
or errors to a nonlinear solution. No Hessian inverse at the candidate,
signed compatibility, stability or physical solution is asserted.

## 1. A strong total scalar source with the actual compatibility

For arbitrary physical compatible forcing f use the same transformed
source and full coupled scalar equation as the preceding theorem:

    (fbar_tau,fbar_c)=lambda^-2 Rcal* f,
    Fsrc=||W fbar_tau||+||W fbar_c||,
    Tbar=-Box_4/lambda, H0=Tbar+4,
    E=60pi/g+11/3, e0=E+4,
    H0(Tbar-E)tau=q=p+V_T y_T-Cbar tau,
    y_T=V_T*tau, p=80pi^2(fbar_tau-Bbar Dbar^-1 fbar_c).

The reference outgoing inverse is R+=(Tbar-E-i0)^-1 and the
positive elliptic inverse is R_e=H0^-1. The exact scalar solution is

    tau=(R+q-R_e q)/e0.                                (2)

The total q has zero reference moment because f is in the actual
coupled range. No reference compatibility is imposed on p alone.
The previous value theorem gives

    ||y_T||<8644sqrt(g)Fsrc,
    ||W tau||<8e5 g Fsrc,
    ||W Tbar tau||<4.8e7 sqrt(g)Fsrc.

[The TT form/strong-source proof](tt-form-and-scalar-fourth-graph.md)
gives ||W p||<=36000 Fsrc, ||W V_T y||<=250g^-1/2||y|| and
||W Cbar tau||<1e-6 Fsrc. The last estimate is for the actual
local differential expression; it differentiates no square-root or
sign decomposition. These retained bounds give

    Q:=||W q||<2.2e6 Fsrc.                              (3)

Set z=R+q=H0 tau, the exact partial-fraction identity. Consequently

    Z:=||W z||<4.801e7 sqrt(g)Fsrc,                     (4)

using sqrt(g)<4e-17. In particular q,z are strong physical L2
fields, and (Tbar-E)z=q distributionally with the regular bolt and
the known fast tail. The regular scalar realization identifies z
in Dom(Tbar). No derivative of q has been assumed.

## 2. Estimate the derivative of each resolvent term

For z solving (Tbar-E)z=q, a weighted form calculation gives

    G_z^2<=E Z^2+QZ+2sqrt(20)G_z Z,
    G_z=||W gradbar z||.

To justify it, start with truncated weights W_R=min(W,R). Since
|gradbar log W_R|<=sqrt(20), the same estimate has constants
independent of R. The unweighted graph realization just established
permits the test W_R^2 z; the inequality bounds its gradient uniformly.
Pass to the increasing weight limit. The bolt is a smooth interior
point, and the regular physical domain supplies the integration by
parts without a new zero value at the bolt. For complex outgoing
intermediate fields take real parts and absolute values. Solving the
quadratic inequality, and using sqrt(a+b)<=sqrt(a)+sqrt(b), yields

    G_z<=sqrt(E)Z+sqrt(QZ)+2sqrt(20)Z<7.3e8 Fsrc.       (5)

Indeed sqrt(E)<15/sqrt(g), sqrt(g)<4e-17 and g^(1/4)<6e-9.
The leading term is 15(4.801e7)Fsrc=7.2015e8 Fsrc; the other
two terms fit the stated reserve. This avoids interpolating the
large second graph of tau against its elliptic low-frequency value.

For v=R_e q, use the existing weighted positive form estimate from
[S164](global-compatible-scalar-response.md), section 3. With
W_c=max(1,(r/(2B))^2), W_c<=W<=4W_c and
|gradbar log W_c|<3, it states

    [||gradbar(W_c v)||^2+||W_c v||^2]^(1/2)<=20Q.

The product rule and Cauchy--Schwarz therefore give

    ||W gradbar v||<=4[||gradbar(W_c v)||+3||W_c v||]
                    <=80sqrt(10)Q<254Q.              (6)

The same form truncation/limit as S164 justifies the weak product
rule at the corner of W_c. Combine (2),(3),(5),(6), with e0>180/g:

    ||W gradbar tau||
       <[7.3e8+254(2.2e6)]g Fsrc/180
       <7.2e6 g Fsrc.                                 (7)

The old interpolation allowance was 6.3e6 g^(3/4)Fsrc. Equation
(7) gains the missing quarter power of g without any additional
source derivative or change in the physical model.

## 3. A sharper unweighted derivative from the positive TT form

For the full metric reconstruction it helps to retain the positive
TT contribution in the elliptic part. Use precisely the form K
already established in the preceding value theorem:

    B_T=V_T/sqrt(e0), B_s=V_s/sqrt(e0),
    K=H0+B_T B_T*-B_s J_s B_s*,
    K>=(1-b0)H0+B_T B_T*, b0=1e-21,
    tau=z/e0-K^-1 p/e0-K^-1 V J V*z/e0^2,
    V=(V_T,V_s), J=diag(I,-J_s).

Here the a+divbar b representation of p has
A=||W a||<=30000 Fsrc, Bv=||W b||<.2sqrt(g)Fsrc. The following
operator estimates are between the actual form and auxiliary spaces:

    ||H0^1/2 K^-1 p||<=(A/2+Bv)/(1-b0),
    ||H0^1/2 K^-1 B_T||<=1/sqrt(1-b0),
    ||H0^1/2 K^-1 B_s||<=sqrt(b0)/(1-b0).             (8)

For the first, use K>=(1-b0)H0 and
||H0^-1/2 p||<=A/2+Bv. For the second factor through K^-1/2;
K>=B_T B_T* implies ||K^-1/2 B_T||<=1, while
||H0^1/2 K^-1/2||<=1/sqrt(1-b0). For the third factor through
H0^-1/2 B_s, whose norm is at most sqrt(b0). These factorizations
require no commutation with H0 and keep the nonlocal TT block.
Since ||gradbar u||<=||H0^1/2 u||, they also bound derivatives.

Equations (4),(5), and the separate vertex coefficients give

    ||V_T*z||<=201||z||/sqrt(g)+.001||gradbar z||
                 <9.7e9 Fsrc,
    ||V_s*z||<=(1.3e-10)||z||/sqrt(g)+2.24||gradbar z||
                 <1.7e9 Fsrc.

Substitute these and (8) in the exact regrouping. Using
(1-b0)^-1<1.001, sqrt(b0)<1e-10 and 180^(3/2)>2400, the
coefficient of g Fsrc is at most

    7.3e8/180+1.001(15000+.2sqrt(g))/180
       +1.001[9.7e9+1e-10(1.7e9)]sqrt(g)/2400
       <4.1e6.

Thus

    ||gradbar tau||<4.1e6 g Fsrc.                      (9)

All form maps used here were already closed on the actual regular
seed domain. Equation (5) supplies the required z derivative before
the adjoints in this calculation are applied.

## 4. First derivatives of both disk return components

Use the actual strong TT source from the preceding theorem:

    f_TT=fbar_c-Bbar*tau, ||W f_TT||<64Fsrc,
    c=Dbar^-1 f_TT=Dbar^-1/2 y,
    y=Dbar^-1/2 f_TT, ||y||<326.4sqrt(g)Fsrc.           (10)

The last inequality uses ||Dbar^-1/2||<5.1sqrt(g) and W>=1.
The full critical inverse contains the actual mass column. From S165,

    ||c||<5.1sqrt(g)||y||, ||gradbar c||<34sqrt(g)||y||,
    ||b||<=2||c||, ||sqrt(lambda)eta||<22sqrt(g)||y||.

Primes below are physical proper-radius derivatives. The exact
[fixed-disk return](stationary-weighted-tt-inverse.md), with radial b,
is S_c=c/2-R eta, T_c=-(b+3c)/2-(A-R)eta, where

    b'=-2A(b+c)-2R(b-c), eta'=A eta+b+c,
    R'=AR-R^2, A'=F_rr/2-A^2,
    Kgeom=F_rr/2-2AR+R^2.

Differentiate before taking norms. The polar singularity cancels:

    S_c'=c'/2-(2AR-R^2)eta-R(b+c),
    T_c'=-3c'/2+2Rb-Kgeom eta.                         (11)

The retained smooth coefficient guards are |R|<=sqrt(5lambda),
0<=AR<=11lambda, R^2<=5lambda and |Kgeom|<40lambda. Thus

    ||gradbar S_c||<650sqrt(g)||y||,
    ||gradbar T_c||<980sqrt(g)||y||,
    ||gradbar(S_c,T_c)||<1200sqrt(g)||y||
                         <391680 g Fsrc.               (12)

For S_c the budget is 34/2+27*22+3sqrt(5)*5.1<650;
for T_c it is 3*34/2+4sqrt(5)*5.1+40*22<980.
The pair uses 650^2+980^2<1200^2. Start on the regular compact
TT core and extend by the spectral-form approximation of S165:
the bounds on differences make the displayed weak derivatives closed.
The bolt-to-r b integral, mass coefficient and outer-fixed eta are
retained throughout. No norm of the singular A alone is introduced.

Since (S,T)=(tau+S_c,T_c), (9),(12) give the full estimate

    ||gradbar(S,T)||<4.5e6 g Fsrc.                     (13)

## 5. Compose with the actual comparison and retain the limit

The preceding residual theorem gives
D(f_ell)=lambda^-2 gamma0^-1/2 Fdisk(f_ell)<2.2e-13 ell^2/g.
The actual Pi0 has norm one in Fdisk, and the original source
conversion is Fsrc<=124lambda^-2 Fdisk. Consequently

    7.2e6*124*2.2e-13=1.96416e-4<2e-4,
    4.5e6*124*2.2e-13=1.2276e-4<1.3e-4.

The previous metric value bound is 1.2e-6 ell^2. Combining its square
with (1.3e-4 ell^2)^2 proves the last line of (1). No derivative
bound on Pi0 as a map of source Sobolev spaces was required; all
input estimates use precisely its established L2 forcing norm.
The signed actual adjoint obstruction remains unevaluated.

Small H1 norms alone cannot preserve a pointwise cold neighborhood,
even on smooth regular radial fields. Here is an explicit reason.
On a fixed small bolt disk choose rho0>0 and let y=log(rho0/rho).
Let p(v)=2772 integral_0^v t^5(1-t)^5 dt on [0,1], extended
as zero for v<=0 and one for v>=1. This is C5, 0<=p<=1,
p(0)=0, p(1)=1 and |p'|<3. For L>1 set

    psi_L(rho)=p((y-L)/L), psi_L(0)=1.

The field is constant one on rho<=rho0 exp(-2L), zero on
rho>=rho0 exp(-L), and C5 including the regular bolt. Direct radial
integration gives

    integral |psi_L|^2 rho d rho <=rho0^2 exp(-2L)/2,
    integral |psi_L'|^2 rho d rho
          =L^-1 integral_0^1 |p'(v)|^2 dv <9/L.

The actual physical measure is bounded above and below by positive
constant multiples of rho d rho on this fixed disk, and W is bounded
there. The fixed derivative/volume normalizations therefore preserve
the limit: its weighted or unweighted first-order norm tends to zero
while psi_L(0)=1. The same example can be smoothed further by choosing
a smooth step with a bounded derivative. This is a counterexample
to a proposed norm embedding, not a physical correction or solution.
It does not prove failure of our particular projected correction.

Thus (1) advances the derivative estimate but leaves S167.1's actual
local higher-jet requirement open. The full nonlinear source mapping
and its preserved family still need higher spatial regularity and
cold weights; neither a small H1 norm nor the separate mass estimate
supplies them. Further coarse interpolation of the same norms cannot
provide the missing pointwise bolt control. Seek sharper estimates of
the actual higher response, using differentiated residual/source
information with the full compatibility and domain preserved.

## Frozen verification contract

Freeze this note and its checker before execution. Reconcile the
predecessor receipts and all shared input hashes; verify the strong
source budgets, weighted quadratic form algebra, elliptic product-rule
constant, coercive form factorizations on a noncommuting exact matrix
fixture, full return derivative identities, physical compositions,
and the explicit radial non-embedding family. The matrix fixture is
algebra only, not a physical state. Preserve all source snapshots and
terminal runs if a check fails. Continuum domain, weighted limits and
the physical inverse remain ordinary arguments, not formal proofs.
No Lean/cloud campaign is selected. All S167 checkboxes stay open.

Remaining to-do list: local higher derivatives and cold weights for
the actual correction; nonlinear family preservation and remainders;
projected full solve and actual signed compatibility; exterior/state,
charge and return; opening, stability, transport, preparation, separate
creation and full-result audit.
