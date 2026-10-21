# Finite trace gradient and the original fixed projection

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary argument
before frozen verification. **The original fixed projection is
transverse throughout the finite comparison. The projected inverse
exists; its useful full response constant and nonlinear solve remain open.**

Use the full retained model, original comparison, paired state/domain,
temperature, clock, magnetic flux and seed cutoff. Put

    z=sqrt(lambda)rho_ell,
    d_tr=||partial_z(Etr/lambda)||_L2(dz).

We prove d_tr<6e-4. This does not meet the earlier sufficient 3e-4
target. It nevertheless closes the ORIGINAL signed-overlap criterion:
the [full seed error](on-shell-adjoint-error-and-trace-gradient.md)
is below .007, leaving a larger finite-error allowance. In the same
common seed-density norm,

    e_ell=||chi_ell-Phi_ell||_D<.13,
    (theta chi_ell,chi0)_D>.78,
    sigma_ell/sigma0>.048,
    ||L_ell^0||/sigma_ell<21.3.                    (1)

Consequently the [exact projection theorem](finite-comparison-adjoint-and-fixed-projection.md)
makes Pi_cut B_ell an isomorphism on its stated mass-extended physical
domain, uniformly along the comparison. Its inverse norm is below
22.5 times the still-unquantified physical range inverse norm.
The cutoff and density conversion are unchanged. No physical signed
compatibility equation is set to zero by this result.

## 1. Low derivatives in fixed seed coordinates

Write y=sqrt(lambda)rho0, t=r0/B, tau=(t-1)/t, nu=tau V,
k=exp(ell V), T=1+tau(k-1), S=log T, w=ell u+S. All derivatives
in this section are y derivatives. The established low-jet equations give

    |u|,|V|,|nu|<1.2e-8, |u_y|<3.2e-7/t,
    |nu_y|<7e-7, |A0 nu_y|<3.432e-6,
    |L4nu|<6.3e-7<2.2e-6,
    L4=partial_y^2+(A0+2R0)partial_y,
    |w_y|<1.03e-6, |S_y|<7.1e-7.                 (2)

Here R0<2.24, A0R0<=11, F0/x<=5t^2, F0,t/x<=22t,
|F0,tt|/x<=21. In particular |nu_yy|<5.7e-6. No bound on
the singular A0 alone is used at the bolt.

There is a correction to the intermediate coordinate identity in
section 1 of the earlier classical disk-contact note. Let zeta be
canonical orbit proper radius, alpha=d zeta/dy=sqrt(C/x)/t and
Ahat=partial_zeta log N. The exact chain rule is

    partial_zeta^2 nu=alpha^-2(nu_yy+R0 nu_y),
    Ahat=alpha^-1(A0-R0),
    Delta_hat nu=alpha^-2(nu_yy+A0 nu_y).

The earlier pullback omitted the R0 nu_y from the second derivative.
Also alpha^2(Delta_hat phi+|grad_hat phi|^2)=E20+R0^2=2A0R0,
not E20 alone. The CORRECT invariant radius equation is therefore

    L4nu=(C/x)t^-2(A_phys^2 delta Rhat_common)/6
           -delta R4_common/(6lambda)+4u A0R0-(nu/3)R40/lambda.

Its absolute right-hand side is still below 6.3e-7: the old budget
used 22 for the coefficient now exactly equal to 2A0R0<=22.
To repair the subsequent D2nu bound, use the stronger elementary
seed inequality F0<=xt^2+1. Its negative mass/charge terms give
F0/(xt^2)<1.001, hence |R0 nu_y|<1.001(3e-7+1.2e-8).
Thus |D2nu|<1.26e-6<2.3e-6 and
|nu_yy|<1.26e-6+3.432e-6<5.7e-6. The numerical D2S bound
2.31e-6, all subsequent sectional/contact budgets, and the earlier
L4nu bound 2.2e-6 remain valid. The incorrect intermediate identity
and its separate (A0+R0) expression are superseded, not reused.

The exact gap factorization also gives, globally,

    |V_y|<8e-7, |tau_y|<=2.24/t,
    |tau_y|/tau<2.65, |tau_yy|/tau<21,
    |tau V_yy|<1.1e-5,
    |L4tau|<=22/t, |partial_y L4tau|<146/t.       (3)

For the first bound use F0/x<=7(t-1)^2 and |(t-1)V_t|<3e-7.
For the relative tau bounds differentiate tau=1-1/t. The exact
polynomial gives F0,t/x<12.1(t-1): its extremal derivative divided
by t-1 is at most 12+2/x+4|d|/x; the thermal addition
delta/[xt^2(t-1)] is below 7.1h. Thus
|tau_yy|/tau<=6.05/t+14(t-1)/t^2<21. Differentiate
nu=tau V twice, using (2), to obtain the fifth bound in (3).
Finally L4tau=F0,t/(xt^2), whose derivative is bounded by
sqrt(5)(21+44)/t<146/t.

We will use two further consequences of the retained outer low jets:

    |u+nu|<1e-6/t,
    |partial_y(tau V_y^2)|<8e-10/t.              (4)

On t<=40, the second bound follows from
40[2.65(8e-7)^2+2(8e-7)(1.1e-5)]<8e-10.
For t>=4, the independent areal estimates give
|V_t|<1.9e-8/t^3 and |V_tt|<3.6e-7/t^4. Hence
|V_y|<4.3e-8/t^2, |V_yy|<2.1e-6/t^2 and the second
bound in (4) has ample margin. The first uses the core value
bound and |u|<8e-9/t^2, |nu|<2e-9/t^2 outside four.

## 2. One third-derivative estimate, retaining the bolt

The derivative needed in the nonlinear remainder is

    |partial_y L4nu|<500.                         (5)

Here is its coordinate conversion, rather than an unquantified C3
assertion. Let ghat be the seed canonical two-dimensional orbit
metric in A units, phi=log t and alpha=sqrt(C/x)/t. Then

    L4nu=alpha^2[Delta_hat nu+2 grad_hat phi.grad_hat nu]. (6)

The coefficient is TWO, consistently with the corrected chain rule
above and the original invariant radius equation.

On the cold and middle cover, in seed point-normalized orbit
coordinates, the following deliberately enlarged component bounds hold:

    |grad_hat phi|<=t, |Hess_hat phi|<=5t^2,
    |grad_hat nu|<3e-7 t, |Hess_hat nu|<3e-6 t^2,
    |partial_i partial_j partial_k nu|<12t^3.    (7)

The first line is the exact seed log-radius ledger. The gradient
and Hessian of nu follow directly from (2), alpha^2>5.9/t^2,
and the two-dimensional conformal Hessian formula; its correction
has component bound 3R0|nu_y|. The numerical numerator
5.7e-6+3(2.24)(7e-7) is below 5.9(3e-6).

For the last bound use the old low jets in the
[completed metric proof](stationary-completed-metric-response-norm.md).
Its middle canonical radial third derivative is below one, and
its normalized radius third jet is below 1e-6 (coefficient 6).
The triangular radial/normal conversion and its inverse have factor
below two: lower normalized products sum to less than
sum_(j=1)^5 (10/100)^j<.112. On the cold patch the larger radius
contract gives only (1e-5)(.00331)(3!100^3)<.199 before that
conversion. Thus 12t^3 is a common reserve, also covering the
lower connection products when returning to the fixed seed chart.
These are spatial jets of the fixed scalar tangent. No derivative
of an ell-dependent coordinate map is taken.

In these normal coordinates |K_hat|<=2. Differentiating the scalar
Laplacian at the center therefore gives

    |grad_hat Delta_hat nu| <=24t^3+40|grad_hat nu|.

For a radial scalar only its radial derivative enters; the two
diagonal third derivatives account for 24. The reserve 40 includes
all first connection derivatives in dimension two. Differentiate (6),
using |Delta_hat nu|<=2(3e-6)t^2 and (7). After division by
t^3 the bracket is at most

    24+40(3e-7)+2(5)(3e-7)+2(3e-6)
                     +2[2(3e-6)+2(3e-7)]<24.001.

Since alpha^3 t^3<2.5^3, this is below 376, inside (5).
This calculation uses regular normal charts through the bolt.

On t>=40 use the retained COMMON tangent expansion
nu=a3 Y^3+R, |a3|<1e-11 and
|R^(j)|<2e-8 j!100^j Y^(4-j)L(Y), j<=3.
Here Y<3/t<=.075, L=1+log(1/Y),
|Y_y|<1.1Y, |Y_yy|<3Y, |Y_yyy|<24Y.
The last bound follows by differentiating Y_y=-sqrt(C/x)/t.
The third-chain-rule coefficient is at most

    139|a3|Y^3+2e-8[6e6(1.1)^3+3(2e4)(1.1)(3)+2400]Y^4L.

It is below 2e-5. The first two derivative contributions to
partial_y L4nu are smaller, using A0<2, R0<2.24,
|A0,y|<15, |R0,y|<16. Thus (5) holds on the outer end too.

## 3. Differentiate the exact small scalar-curvature remainder

Use the [exact remainder identity](finite-comparison-scalar-curvature-remainder.md),
with R=R40/lambda, n=u+nu, B=k-1-ell V k and
N=R_ell/lambda-R-ell dotR:

    N=R[e^-2w(1+2ell n)-1]+ell(e^-2w-1)dotR
      -6e^-2w[ell(k/T-1)L4nu+(B/T)L4tau
                                      +ell^2(k/T)tau V_y^2]
      +2e^-2w[e^(2ell u)-1-2ell u]/(xt^2).       (8)

The derivative is estimated AFTER this cancellation. In particular
the gradient-square term is retained. Each small exponential or
inverse-T factor is within 1.001 of one. Useful bounds are

    |k/T-1|<1.21e-8/t,
    |partial_y(k/T-1)|<8.5e-7/t,
    |B|<2e-16, |B_y|<1e-14.

For the first term of (8), use the exact cancellation

    partial_y[e^-2w(1+2ell n)-1]
       =2e^-2w[-(S_y-ell nu_y)-2ell n w_y],
    S_y-ell nu_y=ell(k/T-1)nu_y+(B/T)tau_y.

This retains its double zero, rather than separately bounding the
two linear terms. The seed trace identity gives
8piG Arest=Etr0, since the reference theory is exactly stationary
at the seed. Thus dotR=-Etr0/lambda+nu t R_t. The preceding
explicit seed trace derivative bound and |d|/x<1.3e-22 give
|dotR_y|<3e-20/t^4. This is not a derivative inferred from a
value estimate or from an unknown varying quantum source.

The following global bounds are for t times the absolute y
derivative of each line/term in (8), in displayed order:

| Term | Bound |
|---|---:|
| R times the exponential/value remainder | 1e-9 |
| dotR times the exponential change | 1e-23 |
| L4nu term, including derivative of its small coefficient | 3.64e-5 |
| L4tau term | 1e-10 |
| tau V_y^2 term | 5e-9 |
| final u exponential | 1e-20 |

For the dominant row use (5): its leading cost is
6(1.001)(1.21e-8)(500). All other product derivatives fit the
displayed reserve by (2)-(4). Their sum is below 4e-5. The seed
and linear curvature derivatives fit the remaining margin. Since
partial_z=e^-w partial_y, slightly enlarging within that margin gives

    |partial_z(Z/lambda)|<4e-5/t.                (9)

## 4. The remaining fifth-derivative anomaly term

Only a very coarse additional estimate is necessary:

    |grad_ell Box_ell R4|<1e28 lambda^(5/2) t^-1/2. (10)

This is proved using the actual C5 metric, not a sixth derivative.
Here is an explicit extension of the component ledger in
[the mixed-form proof](finite-comparison-weighted-mixed-form.md).
In its point-normalized coordinates use

    (M1,M2,M3,M4,M5)=(30,3000,1e6,4e8,1e16).

The first four bounds are already proved there. For the fifth, the
seed orbit curvature polynomial and its first three covariant jets
are bounded by 2 and 1e6 t^j. The normal-coordinate Jacobi expansion
through order five involves only curvature through three and K grad K.
For example its radial metric factor has fifth coefficient
2K K'/45-K'''/90. Polarization in two dimensions costs at most 5^5;
after point rescaling its normal metric fifth partials are below 1e10.
The sphere fifth normal jet is zero. The log-radius identity
phi'=p, p'=a, a'=pDa yields the radial fifth expression
p(Da)^2+3paD^2a+p^3D^3a. Including angular rotations and
normal-coordinate corrections gives a reserve 1e8 t^5: through
this order there are at most 2^5 5! products, each below 1e4 t^5.

The actual log-direction C5 contract, enlarged to .01 j!100^j as
in the preceding component proof, then permits L5=2e10 along
with L1,...,L4=(10,1000,3e5,1e8). Its fifth exponential Bell
polynomial is below 3e10. Adding the normal-metric product terms,
constant frame conversion and the fifth seed jet stays below 1e16.

For clarity a completely overcounted component recurrence suffices.
Set N0=1 and, for j>=1,

    Nj=4 sum_(p=1)^j binom(j,p) Mp N_(j-p),
    Gj=6 sum_(p=0)^j binom(j,p) Np M_(j+1-p),
    Qj=2G_(j+1)+8 sum_(p=0)^j binom(j,p) Gp G_(j-p),
    Sj=64 sum_(p=0)^j binom(j,p) Np Q_(j-p).

N bounds inverse derivatives, G connection derivatives, Q the
(1,3) curvature derivatives and S scalar-curvature partials. Use
N through four, G through four and Q,S through three. Differentiating
g^ij(R4,ij-Gamma^k_ij R4,k) and converting the four gradient
components gives the bound

    B5=32[N1 S2+S3+4G1 S1+4G0 S2+4N1 G0 S1]
       <2.5e20.

The physical conversion costs 2.5^5 lambda^(5/2). On t<=40,
sqrt(40)(2.5^5)B5<1.6e23, well inside (10).

At the outer end retain the already proved differences of jets
zero through four, |Delta partial^j g|<.001 Mj/t^2. The actual
fifth conormal jet has coefficient at most
1.2e-6(5!100^5)/t^3. The normalized Leibniz/exponential/frame
products through five have coefficient sum below (5!)^2 2^5<1e6,
so it also
fits .001 M5/t^2. Each product contributing to B5 has at most
twelve metric/inverse factors, and at most five further frame
factors. Along the short comparison inverse norms are below two.
The telescoping factor is therefore below 20*2^19<2e7.
Consequently the outer difference is bounded by

    2e4(2.5^5)B5 lambda^(5/2)/t^2.

Multiplication by sqrt(t) and maximization at t=40 gives less
than 2e24. The explicit seed derivative, below 9e4(|d|/x)t^-5
in these units, fits the reserve in (10). Connection and frame
changes are included; no polar coordinate is differentiated at
the bolt and no bounded fifth Y derivative is assumed at Y=0.

## 5. Complete trace gradient and actual projection

Keep every term of the full trace:

    Etr=(1+11epsilon lambda/3)Z+epsilon Box Z
        -11epsilon Z^2/72+(7epsilon/12)C4^2
                         +(11epsilon/6)J^2-20epsilon M_em.

The retained curvature values and first derivatives are
|C4|<50lambda/t, |J|<8lambda/t and
|grad C4|,|grad J|<1e15 lambda^(3/2)t^-1/2.
Also Gq^2/(lambda B^4)<4e-27 and R_ell<2.25sqrt(lambda).
With gamma=G lambda<(4e-17)^2, epsilon lambda<gamma/180,
and integral t^-1 dz<200, equations (9)-(10) give

    d_tr < (4e-5)sqrt(200)
              +(gamma/180)(1e28+1e17)sqrt(200)+1e-24
         <5.7e-4<6e-4.                          (11)

The 1e17 reserve includes both differentiated curvature squares
and Z^2; the magnetic derivative fits the last term. All estimates
use the physical proper radius and actual complete trace.

The preceding adjoint theorem now gives

    e0<.007,
    e_ell<6.7[30.1(6e-4)+1e-8]+.0064<.13.

Use the flexible [cold-overlap inequality](cold-reference-overlap-and-adjoint-budget.md),
not its unnecessarily symmetric .1 corollary:

    (theta chi_ell,chi0)_D
      >1.5-(.4+1.1e0)^2-4(e0+e_ell)-e0 e_ell-1e-5
      >.7848>.78.

Since ||chi0||_D<4.007 and ||chi_ell||_D<4.13, the exact
seed normalization yields sigma_ell/sigma0>.78/(4.007)^2>.048
and ||L_ell^0||/sigma_ell<(4.007)(4.13)/.78<21.3.
The finite physical range inverse composed with the unchanged
projection therefore has norm below 1+1.006(21.3)<22.5 times
that range inverse. This proves (1). It does not supply a useful
numerical value of the range inverse itself.

## Verification and next decision

Freeze this ordinary proof and checker. Reconcile the on-shell,
curvature-remainder, mixed-form and cold-overlap receipts and hashes.
Check the repaired gap/curvature budgets, fixed-coordinate Laplacian conversion,
derivative of the exact remainder, local budgets and fifth-jet
component count. Include independent warped scalar derivative
fixtures, a finite-dimensional fixed projection, and omission controls.
The checks do not replace the continuum jet/domain arguments.
Retain all runs and failures. No Lean run is selected.

The projection question is resolved for the original comparison.
Next obtain useful full range-response bounds in the solver norms,
including cold and higher control, and use them to preserve the
nonlinear family and certify the projected solve. The remaining
signed physical compatibility must then be evaluated. Every S167
acceptance task stays open because none is completed by transversality
alone. No stable wormhole or physical realization is asserted.

Remaining to-do list: useful full range/projected response and
cold/higher control; nonlinear preservation and solve; signed physical
compatibility; exterior/state, charge/return, opening, stability,
transport, preparation, separate creation and full-result audit.
