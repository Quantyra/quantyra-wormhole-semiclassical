# Nonlinear constraint residual of the centered completion

2026-10-16, S139/S074. Ordinary argument before focused verification.

The explicit common-disk completion of the unique linear correction has
a bounded nonlinear residual for the reduced covariant potential theory
S0. Its signed mass primitive, relative to r F0, is below 5e-13 lambda^2;
its integrated lapse residual is below 1.3e-10 lambda^2. Both residual sources are
absolutely integrable at the outer endpoint. At the bolt the orbit and
angular residual pressures, divided by the prepared thermal force, are
below 1e-17 lambda^2. These statements concern the classical-plus-reference
operator E0 and its particular finite completion. They do not bound the
change of the complete Dirac source on that completion.

## Operator, tangent and source allocation

Use [the centered potential](stationary-centered-covariant-seed.md),
[the complete fixed-seed source](stationary-centered-exterior-source.md)
and [the linear reconstruction and completion](stationary-centered-first-order-geometry.md).
All of their physical parameters, common subtraction scale, paired
domain and prepared HH ensemble are retained. Write

```
t=r0/B, s=t-1, v=1+h, C=B^2/A^2=1+6x+d,
Q(t)=t F0(t)/C, tau=(t-1)/t,
b=(r0-B) V_r=s V_t.
```

Here V=xi/(r0-B), u is the conformal tangent and a is the Eulerian
lapse tangent. The quantities u,V,b are scalar functions on the
fixed physical conformal disk; they are not new independent sources.
For 0<=lambda<=1 put

```
T=1+tau[exp(lambda V)-1],
r_lambda=r0 T=B+(r0-B)exp(lambda V),
sigma_lambda=sigma0+lambda u,
F_lambda/F0=exp(2lambda(V-u))(1+lambda b)^2/T^2.
```

The common disk fixes the original period and boundary clock. The
comparison reference potential V0(X) is held fixed under variation,
including B,A,lambda_e. Define the S0 mass and lapse operators by

```
M_lambda(chi)=[r_lambda F0(r_lambda)-r_lambda F_lambda]/2,
P_M(lambda)=partial_chi M_lambda,
P_A(lambda)=partial_chi alpha_lambda.
```

The exact Einstein constraints with a sphere-integrated source are
P_M=G rho r_chi and P_A=G(rho+p)r_chi/(r F). On the seed both
operators vanish. Their tangents are M1_chi and a_chi, where
M1=-r0 f/2 and M1_r=G rho0. The nonlinear residual in this note is

```
R_M=P_M(lambda)-lambda M1_chi,
R_A=P_A(lambda)-lambda a_chi.
```

These are fixed-chi operator densities. In particular no difference
between r_lambda,chi and r0,chi is inserted into R_M, and no varying
source prefactor is inserted into R_A. In the complete retained
equations all such prefactor changes belong to the complete source
map Gamma_rest(g_lambda)-Gamma_rest(g0). Adding them a second time
here would double count that response.

The mass norm used below is a signed primitive norm, appropriate to
the exact radial inverse of this particular constraint. It is not an
absolute physical energy bound or an operator contraction theorem:

```
N_M(chi)=M_lambda(chi)-lambda M1(chi),
||R_M||prim=sup_{chi>0}|N_M(chi)-N_M(0)|/[r0 F0(r0)],
||R_A||1=integral_0^D0 |R_A|dchi.
```

The limiting ratio at the bolt exists. Absolute integrability of
the mass residual at infinity is proved separately below.

## Low-jet inputs and global variation bounds

The frozen reconstruction establishes

```
|V|<1.2e-8, |b|<3e-7, |u|<2e-7,
|a|<3.3e-10, TV(a)<3.3e-10,
|V_h|<3.4e-10, |u_h|<8.3e-5 h.
```

The total variation of a follows from its absolutely convergent
source integral, not only from its supremum bound. We record the
additional coarse real bounds needed here; no high areal jets enter.

After the outer vacuum has first been cancelled exactly against the
reference, the full seed source obeys

```
|rho0|+|p0|<1e17 q^2/B^2, t>=4,
epsilon_M<4e-11, epsilon_P<4e-10,
epsilon_M=[|kappa R s_R|+G integral_R^infinity |rho0|dr]/(B x),
epsilon_P=1e17 G q^2/(B^2 x).
```

For completeness the canonical massive bound uses
|Theta|<=4K q^2 t and |J|<=K q^2 F(t-v)/(v t), K=1e15.
Its two pressures sum to less than 5K q^2/B^2. The remaining
vacuum-subtracted Weyl, LLL and nonvacuum reference terms are below
1e6 q^2/B^2. The displayed 1e17 bound has ample slack; it does
not involve differentiating the outer logarithmic remainder.

For t>=4 the explicit seed gives

```
.5x t^2<=F0<=5x t^2, |F0_t/F0|<8/t,
|F0_tt/F0|<24.2/t^2,
|f|<=2x epsilon_M/t, |a|<=epsilon_P/t^2.
```

The latter bound follows directly by integrating
|a_t|<=2epsilon_P/t^3. The exact fixed-disk displacement equation
is xi/B=-F0 integral_t^infinity (a+f/F0)/F0 dz. It yields

```
|xi|/B <=(10/3)epsilon_P/t+10epsilon_M/t^2,
|V|<2e-9/t^2, |b|<1.4e-8/t^2, |u|<8e-9/t^2.       (1)
```

For example b=xi_t/B-V and
xi_t/B=a+f/F0+(F0_t/F0)xi/B. Thus
|b|<=(289/9)epsilon_P/t^2+(292/3)epsilon_M/t^3.
Differentiating these exact first-order identities once, using only
the source VALUE bound above and f_t=-(f+2G rho0)/t, gives

```
|V_t|<1.9e-8/t^3, |b_t|<2.5e-7/t^3,
|u_t|<1.3e-7/t^3.                                    (2)
```

One intermediate ledger is
|(a+f/F0)_t|<=6epsilon_P/t^3+36epsilon_M/t^4,
|xi_tt|/B<=522epsilon_P/t^3+1590epsilon_M/t^4.
Consequently |b_t|<=565epsilon_P/t^3+1720epsilon_M/t^4
and |u_t|<=300epsilon_P/t^3+900epsilon_M/t^4.
This is not a claim about three or more conormal derivatives.

Globally TV(V)<2.6e-5. On R<=r0<=4B integrate
|V_r|<=3e-7/(r0-B), obtaining less than 2.5e-5 because
log(3/h)<82. The integral of (2) on t>=4 is below 1e-9.

Here is an independent finite TV(u) bound. In the full common disk
the seed has exp(2sigma0)<60. To check it, put gamma=kappa B/C,
y=t-v, w=h+C gamma integral_t^infinity dz/F0. Convexity gives

```
F0>=2C gamma y+.95x y^2,
F0<=2C gamma y+6.05x y^2,
C gamma integral_t^infinity dz/F0
 <=.5 log[1+2C gamma/(.95x y)].
```

In exp(sigma0)=sqrt(F0/C)sinh(w)/(gamma t), split
sinh(h+z)=exp(h)sinh(z)+exp(-z)sinh(h). The first term
is below sqrt(6.05C/x)/(.95t)<6.86/t. For the second use
F0<=1+x t^2 and h/gamma=1+O(h), obtaining less than .45.
Thus exp(sigma0)<7.31. This also extends continuously to the bolt.

The global fixed-disk linear curvature bound is
|A^2 delta Rhat|<1e-7. Its Eulerian part is

```
delta Rbar/C=-8piG r0^2 Arest/C-6f/C+6G(p0-rho0)/C.
```

The exact vacuum-cancelled local trace polynomial bounds the first
term by 1e-20; the finite mass and the sharper canonical value bound
5.001e15 q^2/B^2 bound the remaining terms together by 1e-10.
The displacement term is bounded by

```
|V| (t-1)|partial_t Rbar0|/C <=4.25|V|<5.1e-8.
```

Indeed Rbar0=-calB''(1/t) and |calB'''|<=17C, giving
(t-1)|Rbar0_t|/C<=17u(1-u)<=17/4. The same explicit seed
has |A^2 Rhat0|<=2. These scalar estimates are from complete
trace and invariant curvature, not separate f''' or a'' estimates.

The radial linear curvature equation now gives

```
(sinh chi u_chi)'=sinh chi[-exp(2sigma0) A^2 Rhat0 u
                           -.5exp(2sigma0) A^2 delta Rhat].
```

Its forcing has magnitude below 2.7e-5; hence |u_chi|<2.7e-5.
To avoid multiplying by the long cold neck, split at s=.01.
The common-disk length outside that seam is below 6: its reciprocal
seed integral is below 175B/C and w_seam/h<211, so
D0-chi_seam<=log(w_seam/h)<log211<6.
On the cold patch the fixed baseline transition to its auxiliary
coordinate leaves total variation unchanged. There
|V_pot|<3, |u|<=1e-6 W, and the curvature forcing is at most
.56e-7 h cosh chi, with
W<=h(Dc+1)cosh chi and h cosh Dc=.01.
Integration from the smooth bolt gives cold TV(u)<1.3e-6.
The outer length contributes less than 6(2.7e-5). Therefore

```
TV(u)<1.8e-4,
TV(b)<4e-4.                                          (3)
```

For the second bound use the exact tangent identity
a=(2tau-1)V+2u-b and TV(tau)<1. It gives
TV(b)<=TV(a)+TV(V)+2|V|_infinity+2TV(u).

## The signed mass primitive

Put delta=s[exp(lambda V)-1] and

```
H=exp(2lambda(V-u))(1+lambda b)^2/T,
H1=2(V-u+b)-tau V,
N_H=H-1-lambda H1,
N_q=Q(t+delta)-Q(t)-lambda Q'(t)sV.
```

The exact mass expression and its tangent imply

```
N_M=(BC/2)[N_q-Q N_H].                               (4)
```

At t=v, Q=0. The thermal term in Q=tF0/C is constant;
therefore Q' is the extremal polynomial derivative and vanishes
at t=1. Direct differentiation of its degree-three rational formula
gives, also on the slightly displaced positive-radius interval,

```
0<=Q'<=3t(t-1), |Q''|<=3t, |Q'''|<=10,
Q(t)>=.14 t(s-h)s, t>=v.                             (5)
```

The last bound follows from F0>=.99x(t-v)(t-1) and x/C>.16.
The first three can be checked term by term from
Q=(x/C)t^3+t/C-2a_seed/C+(1+3x+d)/(Ct)
-2d/(Ct^2)+d/(Ct^3). The constant thermal mass drops out.

Here is the Taylor estimate including its derivative, needed at the
bolt. With E1=exp(lambda V)-1 and E2=E1-lambda V,

```
N_q=Q' s E2+(s E1)^2 integral_0^1(1-z)Q''(t+zsE1)dz.
```

Use |E1|<=1.001lambda|V|, |E2|<=.501lambda^2 V^2,
|E1_t|<=1.001lambda|b|/s and
|E2_t|<=1.001lambda^2|Vb|/s.
Differentiating the displayed integral and applying (5) gives

```
|N_q_t|<=lambda^2 t s[12 Vmax^2+7 Vmax bmax].        (6)
```

The two terms contribute respectively at most lambda^2 times
ts(3.01Vmax^2+3.01Vmax bmax) and
ts(8.04Vmax^2+3.02Vmax bmax); rounding upward gives (6).
Thus, by integration from v and (5),

```
|N_q(t)-N_q(v)|/Q(t)
 <=lambda^2[12Vmax^2+7Vmax bmax]/.14 <2e-13lambda^2. (7)
```

For the other factor |H1|<1.04e-6. The second derivative of
log[1-tau+tau exp(lambda V)] with respect to V is at most 1/4.
Consequently

```
|logH-lambda H1|
 <=lambda^2[bmax^2/(1-bmax)+Vmax^2/8] <9.02e-14lambda^2,
|N_H|<6.5e-13lambda^2.                               (8)
```

Equations (4), (7), (8) prove

```
||R_M||prim <5e-13lambda^2.                          (9)
```

This signed primitive estimate is uniform all the way to the bolt.
No division by h was used in (7). Its subtraction at the actual
fixed-disk bolt is essential.

All these nonlinear bounds carry lambda^2: the Taylor terms contain
either E2 or two factors E1, and the exponential remainder in H
contains two factors lambda H1. At infinity (1) implies
N_q=O(lambda^2 t^-1), Q N_H=O(lambda^2 t^-1), so
N_M(infinity)=0. At the bolt, |N_q(v)|<4lambda^2 h^2 Vh^2 gives
|N_M(0)|/(BC)<2.4e-19lambda^2 h^2. Thus the total signed residual
mass is finite and extremely small. Absolute integrability is
also genuine: using (1)-(2) in (4) gives

```
|N_H|<1.4e-15lambda^2/t^4, |N_H_t|<5e-14lambda^2/t^5,
|partial_t[N_M/(BC)]|<3e-14lambda^2/t^2, t>=4.        (10)
```

For instance |H1|<4.8e-8/t^2 and |H1_t|<8e-7/t^3;
the differentiated logarithmic remainder is below 1e-14/t^5.
The remaining compact interval is smooth. This proves finite
absolute residual mass without replacing the sharper signed norm
(9) by an unproved small absolute-energy norm.

## The integrated lapse residual

The exact reconstructed lapse has a stronger cancellation than an
estimate of separate conformal and radius derivatives:

```
alpha_lambda=lambda a+N_alpha,
N_alpha=2{log[1+tau(exp(lambda V)-1)]-lambda tau V}
                  +lambda b-log(1+lambda b).         (11)
```

It follows by substituting the tangent identity
a=(2tau-1)V+2u-b into the exact lapse reconstruction.
The conformal tangent cancels out of N_alpha itself.
For the first term in (11), its V derivative is bounded by
3|V| and its tau derivative by 3V^2. For the second, its
b derivative is bounded by |b|/(1-|b|). Hence (3) gives

```
||R_A||1=TV(N_alpha)
 <=lambda^2[3Vmax^2+3Vmax TV(V)+bmax TV(b)/(1-bmax)]
 <1.3e-10lambda^2.                                  (12)
```

The endpoint behavior is also explicit: (1)-(2) make N_alpha
O(lambda^2 t^-4), its t derivative O(lambda^2 t^-5). In the actual sharper
mass-correlated cubic expansion it begins still later; that
improvement is not needed for (12).

## Horizon pressures in the thermal ideal

This section uses mixed physical S0 Einstein pressures, with the
base-prefactor linear tangent subtracted at fixed chi. If converting
between a physical pressure and the constraint densities above,
the resulting linear-source prefactor difference is explicitly
bounded below; it is not counted twice in Gamma_rest response.

At the bolt tau_h=h/(1+h), |Vh|<3.4e-10 and
|u_h|<8.3e-5h. Fixed temperature gives

```
b_h=(2tau_h-1)Vh+2u_h-a_h, |b_h|<6.71e-10,
H1_h=3tau_h Vh+2u_h-2a_h, |H1_h|<6.61e-10.
```

The local versions of (6)-(8) now give
|N_q_t(v)|<3.0e-18lambda^2 h, |N_H(v)|<6.8e-19lambda^2,
and Q'(v)=2kappa R/C<3h. Differentiating (4), with Q(v)=0,
proves the orbit constraint remainder

```
|G rho_rem,h|/(kappa R)<5e-18lambda^2.                (13)
```

The harmless conversion between current and base radial prefactors
is below 1e-20 in this ratio: the source ratio is below 5.1e-12
and exp(lambda Vh)(1+lambda b_h)-1 has magnitude below 1e-9.
Smoothness makes p_rem,h=-rho_rem,h for the mixed orbit tensor.

The angular calculation is independent. Put l0=log r0 and
d0=exp(-2sigma0)l0_chichi at the bolt. Smoothness implies
d0=kappa R/C<1.5h. Exact radius differentiation in the completion
gives

```
exp(-2sigma_lambda)l_lambda,chichi=d0 K_h,
K_h=exp(lambda(Vh-2u_h))(1+lambda b_h)/T_h,
K1_h=tau_h Vh-a_h,
|K_h-1-lambda K1_h|<3e-19lambda^2.                   (14)
```

Indeed the logarithmic remainder is at most
b_h^2/[2(1-|b_h|)]+Vh^2/8, with the latter term more sharply
bounded by tau_h Vh^2/2. Together with |K1_h|<3.31e-10
this proves (14).

The exact normalized angular operator at the bolt is

```
E_S,h=-(A^2/2)Rhat_h+4exp(-2sigma_h)l_chichi,h
                         -W_S(r_h),
W_S(r)=r^2 V0'(r^2)/(2C).
```

Its mixed angular pressure is C E_S/(8piG r_h^2).
The product-frame curvature remainder is below 1e-7lambda^2 h^2/A^2
by the frozen curvature identity. The middle term contributes
less than 1.8e-18lambda^2 h by (14). Finally the explicit rational V0
has |partial_logr W_S|<20 and
|partial_logr^2 W_S|<100 on this tiny horizon interval.
The nonlinear change of log r_h is at most
.501lambda^2 tau_h Vh^2, and its linear change at most lambda tau_h|Vh|.
Thus its potential remainder is below 1.2e-18lambda^2 h.
The resulting angular pressure ratio, including conversion to the
base normalization, obeys

```
|8piG pS_rem,h|/(2kappa/R)<1e-17lambda^2.             (15)
```

This uses C/(2kappa R)<.7/h, from kappa B>5xh and C/x<7.
The curvature h^2 term is negligible relative to h at the fixed
preparation. The estimates (13)-(15) retain the quantum center;
there is no constant order-lambda^2 translation or product-curvature
defect whose ratio would diverge as h approaches zero.

## What this establishes

There are also sharper primitive inverse consequences. Both terms
in (11) are nonnegative, by convexity of the log moment-generating
function and log(1+z)<=z. Their values satisfy

```
0<=N_alpha<=lambda^2[Vmax^2/4+bmax^2/(2(1-bmax))]
            <4.6e-14lambda^2,
N_alpha,h<2.4e-19lambda^2.
```

The latter uses the sharper radius term tau_h Vh^2 rather than a
constant Vmax^2 bound, and |b_h|<6.71e-10. Thus the inverse lapse
amplitude of this residual, with its outer value fixed, is below
4.6e-14lambda^2. The radial mass inverse with zero primitive horizon
datum has relative F amplitude below 1e-12lambda^2 by twice (9).
These are inverses of the displayed constraint densities, not a
bounded inverse of the complete retained equations.

For the actual prepared-temperature linear boundary problem, the
residual's horizon datum must instead be fixed by that boundary
condition. The same explicit linear formula in the reconstruction
note gives

```
|s_res|/(Bh) <=[kappa/(D B h)]
                [5e-18+2.4e-19]lambda^2 <1e-17lambda^2,
D=F0''(R)/2+kappa/R.
```

Here kappa/(D B h)<1.5 follows from |D A^2-1|<7h and
the seed gap bounds. This treats the extra prepared-temperature
mass datum separately from the zero-datum primitive inverse.

The exact finite completion has the global constraint residual
(9), (12), finite absolute outer mass/lapse residuals, and both
horizon pressure bounds (13), (15). This provides the nonlinear
E0 part of the order-reduction error for this explicit family in
the stated signed mass/integrated lapse/horizon norms. The mass
norm is intentionally a signed radial primitive; no positivity
claim for this approximation residual is made.

No complete retained-action error estimate follows until the
actual renormalized Dirac source and its HH/domain response on
g_lambda are bounded in compatible norms. The estimate here is
also not a bound on the inverse of the full higher-derivative
retained equations or a proof of a nearby exact stationary
solution. No raw horizon jet or unproved analytic/C5/C6 ball
membership is substituted for those missing statements.

Remaining to-do list: complete focused verification of the accepted
ordinary argument; bound the complete candidate-evaluated source response and combine
both residuals; justify the approximation or a higher-order inverse;
match the corrected exterior/state and complete boundary Hamiltonian.
