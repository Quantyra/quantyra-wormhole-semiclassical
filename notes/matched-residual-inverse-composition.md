# Actual comparison residual in the compatible inverse norm

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof
and verification contract before execution. **Finite linear
composition bound; no useful Newton certificate or physical solution.**

The prescribed finite comparison has a forcing bound in exactly
the weighted physical L2 norm used by the completed disk and York
inverse estimates. The actual seed-compatible projection has norm
one in this space. Their composition is therefore quantified, with
the independent signed obstruction retained. This closes the norm
matching obligation without promoting integrated constraint bounds
to pointwise forcing bounds.

## 1. Fixed comparison, full forcing and result

Use [S161's actual comparison](comparison-path-complete-equation-residual.md),
[S160's full source allocation](comparison-path-complete-source-response.md)
and [S162's fixed-volume gradient](general-background-newton-operator.md).
Keep the physical action, clock, temperature, flux, state, paired
domain and subtraction. Write ell for the bookkeeping parameter,
lambda=lambda_e for the fixed outer curvature, and g=G lambda.
The action along the comparison is I_ell=I0+ell Gamma_rest;
I_1 is the complete retained action. Its comparison is

    t=r0/B, s=t-1, tau=s/t, nu=tau V,
    k=exp(ell V), T=1+tau(k-1), t_ell=tT,
    sigma_ell=sigma0+ell u, a=exp(-2ell u).

All differential operators below use the fixed dimensionless
product metric in A units, with C=B^2/A^2. The letters u,V here
are the comparison tangents, not the earlier compactified radius.
For the forcing weight write W=t^2 explicitly. Put

    f_ell=mathfrakF_ell(g_ell)=Jvol F_ell(g_ell),
    Jvol=T^4 exp(2ell u), gamma0=4pi beta B^3,
    Fdisk(f)=||t^2(f_Phi,f_V)||_L2(dmu0),
    dmu0=gamma0 t^2 dt.

Then the actual residual satisfies

    lambda^-2 gamma0^-1/2 Fdisk(f_ell)
                     <(2e-8/g)ell^2 <2.2e25 ell^2.       (1)

Let A_m be the FULL retained Hessian at the analytic seed, with
its one physical mass column, and Pi0 the orthogonal projection
onto Ran A_m in this forcing Hilbert space. With the exact G
and Ynorm of [the higher York theorem](quantitative-higher-york-isomorphism.md),

    h_ell=-A_m^-1 Pi0 f_ell,
    G(h_ell)<(2*10^386012/g)ell^2 <3*10^386045 ell^2,
    Ynorm(Rcal^-1 h_ell)<(2*10^411992/g)ell^2
                                           <3*10^412025 ell^2. (2)

At ell=1 this is the seed-preconditioned correction to the actual
retained candidate residual. A_m is not asserted to be the Hessian
at g_1. No transfer of the inverse to that candidate is implicit.
Equations (1),(2) are upper bounds, not evaluated correction values
or lower bounds on an error. They do not certify a Newton ball.

## 2. Derive the missing orbit equation directly

Let Q(t)=tF0(t)/C and H(t)=t^2 Q'(t). The reduced potential obeys
V0(B^2t^2)=2C Q'(t). The orbit trace of the classical-plus-reference
tensor defect is

    E_orb,0=[Box_(g2) r^2-V0(r^2)]/r^2.

This follows by varying the orbit conformal factor in the reduced
action X R2+(grad X)^2/(2X)+V0(X); the kinetic term is conformally
invariant in two dimensions. It is also the sum of the direct
areal-coordinate Einstein equations. Therefore its normalized
profile is exactly

    O0(ell)=a Delta(t_ell^2)-2H(t_ell),
    E_orb,0=C O0/(B^2 t_ell^4).                         (3)

In particular Delta(t^2)=2H(t), so O0(0)=0. Define
L2=Delta+4 grad(log t).grad and B_ell=k-1-ell V k. Then

    L2 T=ell k L2 nu+B_ell L2 tau
                                 +tau ell^2 k |grad V|^2,
    O_*=(4nu-4u)H+2t^2 L2 nu-2tH'nu.                   (4)

The actual tangent solves the complete first-order equations,
so O_*=O_rest(g0). Substitute (4) into (3) before estimating.
The geometric orbit remainder N_O=O0-ell O_* is

    N_O=ell(aTk-1)O_*
       +2H{a[T^2+2ell Tk(u-nu)]-1}
       +2ell aTk tH'nu-2[H(tT)-H(t)]
       +2a t^2 T B_ell L2 tau
       +2a t^2 T tau ell^2 k |grad V|^2
       +2a t^2 |grad T|^2.                            (5)

The last gradient square has positive coefficient two. No
derivative of a varying quantum stress has been introduced.
The radial identities are

    L2 tau=F0_t/C+2F0/(Ct),
    |grad V|^2=t^2(F0/C)(b/s)^2,
    grad T=sqrt(F0/C)[(k-1)/t+ell k b], b=sV_t.        (6)

All formulas extend smoothly to the regular bolt. They do not
divide a constraint primitive by its vanishing horizon factor.

## 3. Regional bounds with the actual source cancellation

S160's local trace allocation, here denoted scriptT, gives

    O_rest=N_rest-2rho_rest=2(scriptT-A_rest).          (7)

Its seed value and first response imply the following bounds,
where a_x=s[1+log(.01/s)] and L_e=1+log t:

| Region | abs(O_*) | abs(Delta O_rest)/ell | abs(Delta scriptT)/ell |
| --- | ---: | ---: | ---: |
| cold, s<=.001 | 3.401e-6 a_x | 5.201e-10 a_x | 3.2e-16 a_x |
| middle, through t=40 | 5.201e-9 | 4.401e-12 | 3.2e-17 |
| outer, t>=40 | 9.401e-7 L_e^2 | 9.807e-11 L_e^4 | 3.2e-14 L_e^2 |

For estimating (5) at t>=4 one may use the larger outer O_*
bound also on 4<=t<=40. This uses the full retained source and
its common cancellations, not separate absolute canonical/Weyl
pieces. Delta means the endpoint value minus its seed value.

Here is an explicit majorant for (5). Set U=abs(u), Wv=abs(V),
N=abs(nu), and D1,D2 bounds for abs(tH'),abs(t^2H'') on the
displaced segment. With the same scalar exponential reserves as
S161, abs(N_O)/ell^2 is bounded by

    1.01(2U+Wv+1.001N)abs(O_*)
    +2abs(H)[10(U+N)^2+5(U+N)Wv+2tau Wv^2]
    +2.02D1 N(2U+Wv+1.001N)+1.002D1 tau Wv^2+1.004D2 N^2
    +1.01Wv^2 t^2 abs(L2 tau)
    +2.02t^3(F0/C)b^2/s
    +2.02t^2(F0/C)(Wv/t+abs(b))^2.                   (8)

Indeed the H Taylor terms are twice the same scalar Taylor
identity used for S161's P; the Delta u term of that angular
formula is absent here. The coefficient polynomial is

    H=(3x/C)t^4+t^2/C-(1+3x+d)/C+4d/(Ct)-3d/(Ct^2).

The guards x/C<1/5.9, 1/C<.001, (1+3x+d)/C<.501 and
abs(d)/C<1e-24 prove abs(H)<=1.1t^4, D1<=2.1t^4,
D2<=6.3t^4, including the tiny displaced segment. From
L2 tau=Ltau+3F0/(Ct), the S161 drift bounds give
t^2 abs(L2 tau)<=6s in the cold, <=6t^3 through t=4,
and <=11t^3 for t>=4.

Use precisely S161's U,Wv,N,b,F0 rows. In the cold use
s/a_x<=1 and a_x<=.00331 before summing. In the middle
all terms have their maximum at t=4 after writing
t^2(Wv/t+abs(b))^2=(Wv+t abs(b))^2. In the outer region
use t>=4 and t/s<=4/3. Exact rational sums then give

    abs(N_O)<ell^2 times
       3e-12 a_x                 (cold),
       7e-11                     (middle through t=4),
       5e-15 L_e^2               (t>=4).             (9)

Retain S161's geometric angular remainder N_S, bounded by
2e-12 a_x, 6e-11, and 2e-15 L_e^2 in these same regions.
The full normalized orbit residual R_V and full trace R_Phi
are exactly

    R_V=N_O-ell Delta O_rest,
    R_Phi=R_V+2R_S=N_O+2N_S-2ell Delta scriptT.      (10)

Thus the canonical contribution cancels in the full trace
before a norm is taken. Sufficient complete bounds are

| Region | abs(R_Phi)/ell^2 | abs(R_V)/ell^2 |
| --- | ---: | ---: |
| cold | 7.1e-12 a_x | 5.3e-10 a_x |
| middle through t=4 | 2e-10 | 8e-11 |
| t>=4 | 7.4e-14 L_e^2 | 1e-10 L_e^4 |

The last row uses the middle source bounds through t=40 and
the outer bounds thereafter. Consequently, on the entire exterior,

    abs(R_Phi)<=2e-10 ell^2 L_e^2,
    abs(R_V)<=1e-10 ell^2 L_e^4.                    (11)

## 4. Exact fixed-density conversion and integration

S162 identifies F_Phi=(E_orb+2E_theta)/(8piG) and
F_V=E_orb/(8piG). The physical radius factor in (3) is
t_ell^-4. It cancels T^4 in Jvol exactly, giving

    f_ell=C/(8piG B^2) exp(2ell u)t^-4 (R_Phi,R_V). (12)

No mass/lapse primitive Jacobian is inserted in (12). Since
abs(u)<1.2e-8, exp(2ell u)<1.001. Equations (11),(12) and
dmu0=gamma0 t^2 dt give

    Fdisk(f_ell)
      <=C sqrt(gamma0)/(8piG B^2) *1.001 sqrt(5) 1e-10 ell^2
                  *[integral_1^infinity t^-2(1+log t)^8 dt]^1/2.

Extending the lower endpoint 1+h to 1 only increases this
positive bound. With x=log t the last squared integral is

    integral_0^infinity exp(-x)(1+x)^8 dx
                 =sum_(j=0)^8 binom(8,j)j!=109601.

The resulting scalar factor is below 7.5e-8. Now
C/(G B^2 lambda^2)=(C/(lambda B^2))/g<6.1/g;
using pi>3 proves the first bound (1). The retained seed guard
.999*2000^-10<g implies 1/g<1.1e33 and proves its second bound.

## 5. Actual compatible projection and composition

Let H_forc be the real Hilbert space with norm Fdisk. The accepted
adjoint theorem gives the continuous nonzero functional

    L0(f)=integral chi0 dot f dmu0, Ran A_m=ker L0.

Its Riesz representative in H_forc is r0=W^-2 chi0. Set
rhat=r0/||r0||_H_forc and

    eta(f)=<rhat,f>_H_forc,
    Pi0 f=f-rhat eta(f).

Then Pi0 is the actual orthogonal compatible projection,
||Pi0||=1, and abs(eta(f))<=Fdisk(f). This normalizes the
obstruction coordinate without requiring the numerical phase or
a chosen compact test source to have a known nonzero moment.
The Riesz representative is in H_forc by the already established
boundedness of L0; no extra boundary data or physical source is
introduced. The original compatible inverse extends to this closed
range in precisely the forcing norm used here.

Apply the completed disk and higher York estimates to Pi0 f_ell.
Their constants 10^386020 and 10^412000, together with (1), give
(2). The corresponding linear residual is exactly

    f_ell+A_m h_ell=rhat eta(f_ell).

Neither eta(f_1) nor the nonlinear signed equation after a
projected correction has been evaluated. Setting either to zero
would change the theorem. The projection estimate does not assert
that the unprojected forcing lies in Ran A_m.

## 6. Decision and verification contract

The explicit global norm matching and projected residual/inverse
composition complete S166's linear estimate. The constants in (2)
do not certify a small correction, a preserved nonlinear family,
an inverse at g_1 or the product 2M^2 L e<1. Their size also
does not prove failure of the actual Newton step or nonexistence.
Further inflation of these general upper bounds is not a useful
route to the physical goal. The next work must obtain a useful
bound on the projected correction/remainder and control the
admissible nonlinear family, then evaluate signed compatibility.

Freeze this ordinary proof and checker before execution. Reconcile
the exact S160, S161, S162, higher-York and outgoing-inverse receipts
and shared hashes. Check the reduced orbit variation, seed identity,
exact finite decomposition, polynomial/region budgets, source trace
cancellation, fixed-volume conversion, improper integral and final
constants. Verify projection algebra in a nontrivially weighted
finite-dimensional fixture; it is not the physical adjoint solve.

Independently compare the product formula (3), the direct areal
Einstein orbit trace, and (5), at 80 digits, for manufactured
u=.02(t-1)/(1+t)^3, V=.02/(1+t)^3. Use the S161 exact normalized
seed with 1/C=.001, d/C=.00002, x/C=(1-.001-.00002)/6,
h=1e-5 and 1e-30, t=1+2h,1.1,4,50, ell=0,.4,1.
Define O_* from this manufactured tangent; these are not retained
quantum solutions. Require normalized error below
1e-30(1+abs(target)) after division by .02^2 min(t-1,1).
Retain every run, failure and exact failed source snapshot.
Finite checks support the ordinary continuum and Hilbert-space
argument; they do not formally verify it or certify a wormhole.

Remaining to-do list: useful projected correction and nonlinear
family/remainder bounds; actual signed compatibility; exterior/state,
charge and return; opening, stability, transport, preparation,
separate creation and full-result audit.
