# Complete retained-equation defect of the centered candidate

2026-10-16, S139/S074. Ordinary argument before focused verification.

The specified finite centered comparison has a quantitative residual
for the complete retained free-Dirac action. The estimates below combine
the actual metric norms, canonical determinant response, explicit local
sources, stationary Ward reconstruction and geometric source factors.
They concern this candidate's equation defect. They are not a distance
estimate to a solution of the full higher-derivative equations.

## Action split, candidate and exact allocation

Keep the physical q, prepared temperature, outer clock, paired domain,
common scale and finite coefficients of
[the centered seed](stationary-centered-covariant-seed.md). Write

```
S_lambda=S0+lambda Gamma_rest,
DE0(g0)g1+J_rest(g0)=0,
sigma_lambda=sigma0+lambda u,
r_lambda-B=(r0-B)exp(lambda V), 0<=lambda<=1.
```

The tangent and finite comparison are the explicit ones in
[the reconstruction](stationary-centered-first-order-geometry.md).
The original finite comparison, not an auxiliary interpolation used to
estimate the determinant, is evaluated in the equations. Its period,
outer clock and product center are fixed exactly.

The full residual is

```
R_lambda=E0(g_lambda)+lambda J_rest(g_lambda)
 =[E0(g_lambda)-lambda DE0(g0)g1]
       +lambda[J_rest(g_lambda)-J_rest(g0)].        (1)
```

Each bracket is estimated in the same fixed conformal disk. No varying
source Jacobian is added to the first bracket: it belongs to the second.
The first bracket is bounded in
[the E0 calculation](stationary-centered-constraint-residual.md).
The present note supplies the complete second bracket.

The ingredients are ordinary finite estimates in
[the actual metric norm proof](stationary-completed-metric-response-norm.md),
[the canonical response](stationary-centered-canonical-response.md),
[the local response](stationary-centered-local-response.md), and
[the conormal source proof](stationary-centered-conormal-source.md).
They retain the complete angular tower, the HH thermal winding, actual
complementary-adjoint domain, common Mellin finite parts and local
radius source. The differentiated source proof is an input to the
metric norm, not an inference of metric response from a spatial gradient.

For each fixed endpoint lambda, use its fixed proper-coordinate
identification and a fixed-endpoint logarithmic chord to estimate the
source difference. This avoids differentiating the original comparison
in a moving proper gauge, which would request an unavailable sixth
metric derivative. The auxiliary chord has the same endpoint source;
it does not replace the original candidate in (1).

## Finite input constants and units

Let `t=r0/B`, `s=t-1`, `C=B^2/A^2`, `x=lambda_e B^2`. The fixed
benchmark gives

```
G C/B^2 <8.5e-33,
G C^2/(B^2 x)<5e-32,
gamma=Gq^2/(CB^2)<6.8e-28,
L_h=1+log(.01/h)<81.
```

All hats below are dimensionless profiles: for the FULL combined
sphere-integrated source put

```
H_rho=r^2 rho2, H_N=r^2(rho2+p2).
```

Canonical orbit trace is `Theta=B^2 tau_hat`. It is not divided by
q^2 in the following budgets. The actual endpoint field norms are

| Region | Orbit C5 | Radius C5 | Low operator coefficients |
|---|---:|---:|---:|
| physical cold, s<=.001, weight s L_c(s) | `lambda*1e-6` | `lambda*1e-5` | separately small |
| compact non-bolt | `lambda*1e-4` | `lambda*1e-4` | `lambda*1e-6` in the low norm |
| outer, after the mass cubic | `lambda*1e-7` | `lambda*1e-7` | separately small |

The derivative normalization is j!100^j. The compact response uses
two norms: local UV coefficient 1e22 times the high norm, and gluing
coefficient 1e24 times the low norm. Thus its canonical budget is
`2e18 lambda C^2`; the high norm is not used as a global domain
perturbation. The outer pure analytic mass mode has amplitude below
`4e-11 lambda`; its response beyond the exact cubic is retained.

The following enlarged coefficients include the canonical response,
local Weyl/LLL/reference sources, coordinate pullback and radius-source
allocation, before the final mass/lapse Jacobians:

```
T_c=1.1e17, T_osc=2e18,
T_m=3e18, T_o=2e17.                                (2)
```

Every profile bound below has an additional factor lambda. For the
moment set lambda=1. The cold canonical value coefficient is at most
1e23 times 1e-6, whereas the cold oscillation coefficient is at most
1e24 times 1e-6. They must not be replaced by the same larger constant.
The local value coefficient 1e20 times the radius allowance 1e-5
adds at most 1e15 C^2. Its local first spatial response is at most
2e22 times 1e-5 C^2, fitting T_osc. Algebraic q^2 terms and reference
profiles have the stated smaller bounds in the local proof.

For the middle, `1e22*1e-4+1e24*1e-6=2e18`; the local additions
are below 1e16 and fit T_m. For the outer collar,
`1e24*(1e-7+4e-11)<1.001e17`; local additions and finite mass-mode
remainders fit T_o. The position-independent common a0 term has
zero metric response at fixed B and mu0. No product curvature
response is thereby discarded.

## Ward reconstruction and its finite metric term

In the fixed common disk let `W=exp(2sigma)sinh^2 chi`. The canonical
HH regularity condition fixes its pressure, including the integration
constant, by

```
P[Theta,W]=(2W)^(-1) integral_0^chi W' Theta,
rho_hat=P-Theta, N_hat=2P-Theta.
```

The constant part of Theta gives P=Theta/2 for every W and produces
no metric-kernel response in N. Split the exact endpoint difference
into trace variation at W1 and metric variation at the frozen Theta0.
For the latter introduce `W_alpha=exp(2alpha u)W0`. Then exactly

```
Delta N_metric
 =2 integral_0^1 dalpha integral_0^chi dt
    [W_alpha(t)/W_alpha(chi)] [u(chi)-u(t)] Theta0'(t). (3)
```

Only the frozen source gradient occurs in (3). There is no unproved
gradient estimate on the quantum-source response. This is also why
the position-independent a0 value cannot generate a large metric
response here.

For this metric-kernel estimate only, the frozen barrier remains
valid through the auxiliary seam s=.01; it does not extend the
physical cold quantum-response patch beyond s=.001.
Here `|Theta0_t|<=1e15q^2`. The actual conformal field
has its s L_c ideal. For chi>=1, the integrand gives the bound
`|Delta N_metric|<=10*1e-6*1e15q^2 (s^2-h^2)L_c(s)`;
the smooth bolt gives the same or smaller result for chi<=1.
Its lapse image is below 1e-18.

For the remaining exterior, W0(past)/W0(observation)<=1.01.
This does not require a false global monotonicity assumption.
Indeed `calB_z=2z T(z)`, where

```
T(z)=1-3a z+2Qc z^2-5d z^3+3d z^4,
T''=4Qc-30d z+36d z^2>0.
```

T has one zero before the horizon. It lies below .4/a. Beyond the
corresponding maximum toward infinity, `calB>=x` and
`calB_max<=x+2/a^2`, so their ratio is below 1.01. Before that
maximum it is monotone. The small exponential u factor enlarges
the past/current W ratio to less than two.

The frozen endpoint bounds `|u|<2e-7` on 1<=t<=4 and
`|u|<8e-9/t^2` thereafter give `integral |u|dt<6.02e-7`.
Equation (3) therefore gives `|Delta N_metric|<1e10q^2` for
t>=1.01. Its mass-primitive and absolute-infinity images are below
1e-16, and its lapse image below 1e-13. These deliberately loose
allowances are sufficient below. The canonical Ward argument is
not applied to the nonconserved reference source separately;
local/reference stress is added by its explicit variational formula.

## Regional source profiles and the bolt

Write `L_c(s)=1+log(.01/s)`. In the physical cold patch s<=.001 the preceding
trace and local estimates imply

```
|Delta H_rho|, |Delta H_N| <=2 T_c C^2 s L_c(s).
```

At the bolt the full regular mixed source has N=0. For chi<=1,
the canonical oscillation and local first-spatial estimate imply

```
|Delta H_N|<=2 T_osc C^2 h L_h sqrt(chi).            (4)
```

The direct local stress uses C5, not only its value estimate: its
rho+p combination vanishes at the bolt, and the first spatial
bound gives chi, which is smaller than sqrt(chi) here. A mere
constant bound O(h) would not establish lapse integrability.

On the compact non-bolt interval use
`|Delta H_rho|,|Delta H_N|<=2 T_m C^2`. The cold and non-bolt
regions overlap; integrate the cold estimate through s=.001.
The auxiliary H2 comparison seam at s=.01 does not enlarge this
physical response region.

The canonical response near infinity is first expressed in proper
Y and then pulled to fixed chi. Its frozen regular remainder obeys
`|R0_Y|<=2e13q^2/Y`; the coordinate change is
`O(pY^4+1e-7 Y^5 L)`. Its contribution fits (2). The principal
c/Y profile is kept at its actual moving Y, not frozen before this
pullback.

The completion has no cubic areal-lapse term. Its exact geometric
mass coefficient therefore gives opposite canonical and Weyl terms
`+lambda_e m_lambda/(5pi r)` and
`-lambda_e m_lambda/(5pi r)` in the physical density. Add these
before taking an absolute norm. The reference outer vacuum is the
same on both endpoints and cancels identically. This uses the
retained derivative estimates, not just value asymptotics.

Let t_o denote Y=.05. The geometric guards give 40<t_o<60.
For t>=t_o, the Ward memory from the interior is bounded and
`|(log W)_Y|<=2Y`. The outer integral of Y L(Y)^2 is finite.
After the correlated mass cancellation one may use

```
|Delta H_rho|,|Delta H_N|
 <=C^2[2T_m+4T_o L(Y)^2], t>=t_o.                  (5)
```

This estimate includes the pure analytic mass-mode remainder.
The leading cubic is not assigned a separate divergent integral.

## Mass, lapse and horizon norms of the source change

The final norms are norms of fixed-disk densities. With the exact
multipliers (10), define

```
Delta M_density=[J_M H_rho,lambda-H_rho,0] r0,chi/r0^2,
Delta A_density=[J_A H_N,lambda-H_N,0] r0,chi/(r0^3 F0).
```

The zero-horizon-datum primitive norm is

```
Q_M=sup_(chi>0) G|integral_0^chi Delta M_density dt|/[r0 F0(r0)].
```

In the radial calculations immediately below, `Delta rho2 dr` is
shorthand for this density after the hat-profile change and the
separate Jacobian cross term have both been included.

On the cold patch `F0>=.9C(s^2-h^2)/t^2`. Since L_c decreases,
`integral_h^s yL_c(y)dy<=L_h(s^2-h^2)/2`.
The cold primitive is below `2(GC/B^2)T_c L_h<1.6e-13`.

For the middle let s0=.001, t0=1+s0. Anchor the cold primitive;
do not bound it by total exterior mass divided by s0^2. The
source's physical t^-2 factor gives

```
integral_t0^t dv/v^2=(t-t0)/(tt0),
max_(s>=s0)(s-s0)/s^2=1/(4s0).
```

With `rF0>=.9BC s^2/t`, the additional primitive is at most
`600(GC/B^2)T_m<1.53e-11`. The large-radius denominator is
`rF0>=.5Bx t^3`, so (5) cannot create a larger outer maximum.
Including the metric term and the prefactors below,

```
Q_M<2e-11.                                         (6)
```

Absolute infinity mass is a different norm. In units BC/G it is

```
Q_inf=G/(BC) integral_0^D0 |Delta M_density|dchi.
```

The cold integral is negligible. The middle contributes less than
`2T_m`; (5) adds less than `2T_m/40+4T_o`, since
`integral_40^infinity (1+log t)^2/t^2 dt<1`.
The total dimensionless coefficient is below 7e18. Thus

```
Q_inf<1e-13.                                       (7)
```

The lapse norm is `Q_A=G integral_0^D0 |Delta A_density|dchi`;
its radial expression is `G integral |Delta(rho2+p2)|/(rF0)dr`
with the same fixed-density convention.
For chi<=1, `ds/F0<=3dchi/(C h chi)` and (4) give at most
`12(GC/B^2)T_osc L_h<1.66e-11`. For the rest of the cold
patch, `ds/F0<=4ds/(Cs^2)` and
`integral_h^.001 L_c ds/s<(81^2)/2` give less than
`2.63e4(GC/B^2)T_c<2.46e-11`.
For the middle use `integral_t0^infinity dt/F0<1300/C`.
Indeed the retained s=.01-to-infinity estimate is 175/C, while
F0>=.8C s^2 on .001<=s<=.01 gives an additional
`(1000-100)/(.8C)=1125/C`.
Its contribution is below `2600(GC/B^2)T_m<6.63e-11`.
The outer t^-5 integral from (5) is below 1e-17. Therefore,
including (3) and the final prefactors,

```
Q_A<1.2e-10.                                       (8)
```

At the bolt HH regularity gives the orbit pressure equal to minus
the density. The canonical half-trace and the smaller local addition
give `|Delta H_rho,h|<=T_c C^2 h L_h/2`. Since
`kappa B/(xh)>5`,

```
8piG |Delta pA,h|/(2kappa/R)
 <=[GC^2/(B^2x)] T_c L_h/10 <4.5e-14.
```

The angular pressure follows from the FULL retained anomaly,
`2pA+2pS=A4`, with the explicit local radius/reference response.
Its additional normalized force is below 1e-14: the local full-trace
response is at most `1e15 C^2 h L_h` in the hat allocation, so its
half-trace contribution is bounded by
`(2pi/5)[GC^2/(B^2x)]1e15 L_h<5.1e-15`.
The smaller reference/radius terms fit the remaining allowance.
Hence both are
bounded by

```
Q_hor<1e-13.                                       (9)
```

The extra r^-4 frame change at the moving bolt multiplies the
FULL seed horizon source. Its thermal-force ratio is already below
5.1e-12 and `|Delta log R|<3.4e-10 h`; hence this cross term is
below 1e-20. The raw separate product stresses are not used for it.

These are forces at the actual moving bolt pulled to the common
disk, not a same-areal-radius comparison that loses the thermal gap.

## Fixed-disk Jacobians, counted once

For the hat-profile allocation above write

```
T=r_lambda/r0,
r_lambda,chi/r0,chi=exp(lambda V)(1+lambda b),
F_lambda/F0=exp[2lambda(V-u)](1+lambda b)^2/T^2.
```

The exact mass and lapse multipliers relative to their seed kernels are

```
J_M=exp(lambda V)(1+lambda b)/T^2,
J_A=exp[lambda(2u-V)]/[T(1+lambda b)].               (10)
```

They obey `|J_M-1|<4e-7lambda`, `|J_A-1|<8e-7lambda`.
Multiply the changed profile by J once, and apply J-1 once to
the FULL combined seed profile. The frozen seed bounds then give
additional mass image below 3e-18 and lapse below 2.7e-16.
The raw canonical and Weyl pieces are not separately put into an
absolute seed mass norm. Such a separation would destroy their
outer cancellation.

Equivalently one may work with physical rho2,p2 profiles and use
different radius powers in J. Those two allocations cannot be mixed.
Equation (10) specifies the allocation used for (6)-(9).

## Complete equation defect and its interpretation

The original parameter lambda multiplies every source-change bound
(6)-(9). Equation (1) supplies a second lambda. Add the E0 bounds,
which already have lambda^2, to obtain

| Complete retained-equation norm | Bound |
|---|---:|
| signed zero-datum mass primitive | `2.1e-11 lambda^2` |
| integrated lapse constraint density | `2.6e-10 lambda^2` |
| each horizon pressure / prepared thermal force | `2e-13 lambda^2` |

The source's absolute infinity-mass image is below `1e-13lambda^2`.
The E0 mass residual has a finite absolute integral and its signed
endpoint is already controlled separately. A geometric mass statement
is not a renormalized Hamiltonian charge statement.

For the corresponding component Einstein inverse, the zero-datum
relative F image is below `4.2e-11 lambda^2`. The lapse amplitude
image is at most `1.3e-10 lambda^2`, using the source integral
and the sharper E0 amplitude below 4.6e-14, rather than its total
variation. A separately prepared-temperature gap image is below
`2e-10 lambda^2`, using `1.5*(1.2e-10+1e-13)` for the source
piece and the smaller E0 prepared-gap contribution.
Those are defect images under the specified S0 constraint inverse;
they are not an inverse bound for the full retained operator or an
error bar on an unknown exact higher-order solution.

The order-reduction prescription is precise: expand the covariant
action S_lambda regularly in lambda, solve the new coefficient with
the S0 constraints and Ward identity, and evaluate all quantum higher
metric derivatives on already obtained lower-order fields. It does
not admit independently chosen Planck-scale or nonanalytic-in-lambda
homogeneous modes. The explicit finite completion is a representative
of this first-order branch with the stated O(lambda^2) equation defect.
No exact full stationary branch is inferred from that defect alone.

The isolated LLL 1e-8 component target is not used as a complete clock
allowance. The displayed numbers are computed residual norms and
their scope, not a newly fitted physical pass criterion. Corrected
exterior/state matching, the full boundary Hamiltonian, transfer of
the finite-pulse theorem and the remaining physical gates stay open.

Remaining to-do list: independently review and verify this complete
ordinary combination; determine what further inverse or controlled
branch statement is needed for the physical solution; match the
corrected exterior/state and charge, then transfer the driven-source
and capture/tail/area requirements.
