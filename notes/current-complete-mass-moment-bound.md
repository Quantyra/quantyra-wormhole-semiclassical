# Complete source mass moment from value heat bounds

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: PROGRESS. The stable semiclassical wormhole objective
is incomplete; the closed audit and original admission targets remain.

This attempts the COMPLETE source pairing of the regular reference
adjoint, not an isolated contact, a selected mode or an unknown response
replaced by its allowance. The preceding 8.001e15 unsigned mass coefficient
fails that pairing target. Its largest loss came from first estimating a
spatial derivative and then recovering a value. Here use the retained
heat VALUE proof and integrate the whole conserved source directly.
The complete quantum Hessian residual remains a separate requirement.

## 1. Actual realization and the applicable continuum argument

Use the exact current Q, positive gamma, z_h, paired complementary-adjoint
square, original HH period, mu0 and local coefficients. Set t=r/B, z=1/t,
v=1+h, beta=lambda B^2/C, q=2000^4 and r_c=C/q. The smooth source on the
physical half has its one-sided mass cubic and finite-log outer expansion;
no smooth physical double is substituted. The balanced spin trace identity
uses the closed-form even realization, whose coefficient is C2 with a
bounded piecewise third derivative. The regular bolt is unchanged.

The current rational-family hypothesis check is in
fixed-clock-seed-bolt-feasibility.md. Its strict geometry guards transfer
the value estimates in stationary-reference-nonbolt-stress.md, including
the fixed mu0 restoration. They require beta>1/7, C/q<.0035, M1>18,
the actual geodesic boundary, the specified paired projectors and the
positive-square localization; all hold at this seed. The old reference
mass allowance is NOT transferred to the present c_R theory.

For each real thermal momentum the auxiliary cap potentials are
k^2 f^2 +/- k f'. They are continuous, including at the C2 seam, and
bounded below by -1 since |f'/f|<=2. Thus the rescaled Brownian-bridge
formula uses precisely the continuous semibounded case of
[Simon, Theorem 1.1](https://arxiv.org/html/math-ph/9907022), inspected
2026-10-25. It does not require a smooth reflected metric. The caps are
auxiliary only: actual-to-cap transfer uses the original positive
resolvent and its finite-domain return, not a cap boundary condition.

Keep the common angular regulator until the local terms are restored.
The formulas (32)-(33) of the non-bolt proof include a0, a1 cancellation,
a2 and its pole's finite +2. Its current finite-scale a0 value is <167q^2;
|Fang|<100. No mode is renormalized separately. Every thermal winding
and angular level is retained by the bounds below. Their integrable
majorants justify the later spatial integration after this common
restoration. No metric derivative or neighborhood estimate is inferred.

## 2. A boundary heat VALUE remainder, with an explicit constant

Use proper distance Y from the outer boundary in ell=B/sqrt(C) units,
heat time tau, and the same actual even and smooth auxiliary caps. For
0<tau<=.01 define

    R_k=k_even-k_smooth-K3,
    K3=c/(2pi tau^2) integral_0^tau
             E[|Y+X|^3-(Y+X)^3] ds, |c|<1.             (1)

All kernels here include the rank-two spin trace. The retained bridge
has S=sup|sqrt(2) bridge|, P(S>x)<=2exp(-x^2), and probability of crossing
zero exp(-Y^2/tau). Its exact FK factor is J=f(Y) H(Fbar,Gbar).
The already proved pathwise bounds (20)-(23) give, on collar paths,

    |H_even-H_smooth+DeltaF3/2| <=1.3e11(V^4+tau V^2),
    f(Y)<=3, |f(Y)-1|<=12Y^2, |DeltaF3|<=4V^3,
    V=sup|Y+sqrt(tau) bridge|<=.4.

Consequently the complete density-converted remainder is at most
4e11(V^4+tau V^2). On a crossing path V<=2sqrt(tau)S. Using
E S8<=48, E S6<=12, E S4<=4 gives

    E(16S^4+4S^2)^2 <=13888.

Cauchy--Schwarz on the crossing event, followed by the FK divisor
2pi tau, bounds this local contribution by

    [4e11 sqrt(13888)/(2pi)] tau exp(-Y^2/(2tau))
          <8e12 tau exp(-Y^2/(2tau)).                    (2)

For Y<=.1 an exit requires S>=.3/sqrt(tau). The global cap value
factors satisfy |J_even|+|J_smooth|<=5e10 exp(tau). The exit probability
and the cubic Gaussian tail therefore bound the omitted contribution by
2e10 tau^-1 exp(-.09/tau). For the cubic tail one may use, for a>=3,
E[S^3 1_(S>a)]<2.4a^3 exp(-a^2), and on exit
V<=(4/3)sqrt(tau)S. Its extra heat coefficient is <tau^-1 exp(-.09/tau).
Since Y^2/4<=.0025 and

    sup_(tau<=.01) tau^-2 exp(-.0875/tau)<1.585,

this contribution fits below 3.17e10 tau exp(-Y^2/(4tau)).

For Y>=.1, without collar restriction, equality of the cap coefficients
on the positive half gives the direct crossing bound

    |k_even-k_smooth|<=1e10 tau^-1 exp(-Y^2/tau).

Also |K3_Y|<=exp(-Y^2/(2tau)) and K3(infinity)=0 imply
|K3|<=2sqrt(tau)exp(-Y^2/(2tau)). The two ratios to
tau exp(-Y^2/(4tau)) are at most

    1e10 * 64/[9 e^2 (.1)^4] <9.624e13,
    2sqrt(2) exp(-1/2)/.1 <18.

The retained Poisson bound for ALL nonzero windings has minimum cap
circumference >=1e18 and Y<21 on the non-bolt patch. Its value contribution
is <tau^2 exp(-Y^2/(4tau)), absorbed with strict margin. Thus on the
entire non-bolt observation patch

    |R_k|<=C_v tau exp(-Y^2/(4tau)), C_v=2e14.          (3)

This is a value bound derived from the actual FK integrand. It is not a
derivative of an unsigned estimate or a smooth-seam heat expansion.

## 3. Complete common-restored trace bound

Let Atilde=sum_(n>=1) d_n M_n^2 exp(-M_n^2 tau), with BOTH angular signs.
The retained exact sum estimates are

    Atilde<=24C/tau^2,
    |Atilde-2C/tau^2|<=4q/tau+8sqrt(C)/tau^(3/2),
    integral Atilde tau^2 dtau<8C^2/q.                 (4)

Only use the discrepancy on tau<=.01. Subtract the entire universal
cubic profile, whose exact integral is Theta_sing=c C^2/(5pi Y).
The removed long model tail must be added to the error. The short cubic
discrepancy, using the value bound following (2), is at most

    1.6 Cq+16 C^(3/2) E1(25Y^2),

where E1(x)=integral_x^infinity exp(-s) ds/s. The long model tail is
at most 80C^2. The complete short reflected remainder (3) is bounded by
24 C_v C^2 E1(25Y^2). The short smooth remainder and long even-cap value
bounds are respectively 8e15 C^3/q and 1400(1e10)Cq. They hold at EVERY
non-bolt point, not just at the old t=2 anchor. The actual-cap return,
long polynomial subtractions, current a0/a2 restorations, model tail and
the 1.6Cq term fit within an additional 200q^2.

Put R(t)=Theta(t)-Lt, L=-a lambda B^2/(5pi). The actual reciprocal-distance
calculation gives |Theta_sing-Lt|<q^2 for Y<=1. For Y>=1 use t<4,
|L|<q^2/15 and |Theta_sing|<C^2/(5pi), giving the same bound.
It too fits inside 200. Hence for t>=t_c=1.001,

    |R(t)|/q^2 <= B0+A0 E1(25Y(t)^2),
    B0=8e15 r_c^3+1.4e13 r_c+200 <5e10,
    A0=24 C_v r_c^2+16 C^(3/2)/q^2.                   (5)

This is a bound on the COMPLETE massive trace in the original scheme,
including its actual state/domain return. It contains no selected finite
list of angular or thermal modes.

## 4. Integrate the full conserved mass with the correct weight

In the exact original constraint let Bfun=F/t^2=CQ and

    I(t)=integral_v^t Bfun'(s)Theta(s)ds,
    B^2 rho_mass=I/(2F)-Theta/t^2.

Write Theta=Lt+R. Add the retained Weyl L/t BEFORE integration. The
combined source is I_L/(2F)+I_R/(2F)-R/t^2+R_W+Delta R_W+rho_LLL.
For the R part Tonelli with its absolute majorant gives

    integral_v^infinity |I_R/(2F)-R/t^2| dt
      <=integral_0^z_h w(z)|R(1/z)| dz,
    w=1+|Q'(z) I_Q(z)|/2, I_Q=integral_0^z du/Q(u).    (6)

Both original source equations and HH regularity fix the integration
constant. No independent mass or bolt datum is selected to cancel it.

The exact factorization Q=(1-z)^2 R_Q(z)-tau_seed z^3, with R_Q increasing,
R_Q(1)=1 and tau_seed<1e-65, gives on z<=z_c=1/1.001

    I_Q<=1.001 z/[beta(1-z)], |Q'I_Q|<=15z, w<9.      (7)

Here beta>1/7 and 1-z>=1/1001 make the retained thermal correction
<1e-58 relatively. The tiny positive part of Q' is supported in z<.01:
there Q>=.16, Q'_+<5e-11 z, I_Q<=z/.16. Elsewhere the exact Q' formula
is negative. Since Q I_Q tends to zero at the actual simple bolt,

    integral_0^z_h w dz=1.5 z_h+integral_(Q'>0) Q'I_Q dz <1.501. (8)

The current geometry gives Y>=sqrt(5)z. Therefore the entire reflected
term in (5), with its boundary logarithm, has the finite integrated bound

    integral_0^z_c w E1(25Y^2)dz
      <=9 integral_0^infinity E1(125z^2)dz
      =9sqrt(pi)/(5sqrt(5)).                          (9)

The equality follows by Tonelli; it does not exchange divergent signed
sectors. This is the quantitative contribution of (3) to the full moment.

Retain the existing cold derivative bound |Theta'|<3e14q^2 without
optimizing it. At t_c, Y>=sqrt(5)/1.001 makes the E1 term in (5) <q^2.
Over v<=t<=t_c this gives |R|<3.51e11q^2. Equations (7)-(8), now used
at z_c only, imply

    integral_z_c^z_h w dz
      =1.5(z_h-z_c)+Q(z_c)I_Q(z_c)/2 <.006.          (10)

The cold part therefore costs at most 2.106e9q^2. No inverse temperature
or fictitious bolt wall is introduced.

For the linear trace model, its combined Weyl cancellation leaves I_L/(2F).
Its absolute integral is at most
|L| integral_0^z_h |Q'I_Q| dz/(2z). Equations (7),(10) give an integral
<16, so this costs <q^2. The complete Weyl remainder, degree-five/six
perturbation and LLL add <10005q^2 by the retained exact polynomials.
Round all these terms upwards to 10010q^2.

Combining the WHOLE source, not just its boundary term, proves

    B integral |rho_D-rho_vac| dr <= C_mass q^2,
    C_mass=1.501 B0+[9sqrt(pi)/(5sqrt(5))]A0
                     +2.106e9+10010.                 (11)

There is no unknown constant in (11). All full-source pieces and common
finite restorations have specified allowances. The old 8.001e15 bound
remains valid but is no longer the chosen estimator for this moment.

## 5. Original units, execution contract and remaining coupled test

Use the same local and bolt allowances evaluated in
current-reference-clock-adjoint-trial.md, with L_q=Gq^2/(lambda B^4):

    M_local=8pi(5.5e-23)+4pi c_R g(200w0+w0^2+10w2),
    M_bolt=(gamma/beta)(1+h)h(3.33e-10),
    |F(Psi_ref)| <= I_new=2(C_mass L_q+M_local+M_bolt). (12)

This is the same convergent compact-test source moment and original
physical mass normalization. No full-action open-clock boundary one-form
is discarded or added twice. K_m, the original graph conversion and the
cut mass lift are unchanged. Freeze the execution targets

    C_mass<1.7e11, I_new<1.3e-15<3.3e-15.             (13)

If those tests pass, they close the COMPLETE reference-adjoint source
PAIRING allowance by an enclosure about zero. They do not determine its
sign or produce a narrow signed quadrature value. No extra quadrature
error is assumed: (12) is itself the outward analytic error enclosure.
Its contribution to the actual coupled mass certificate is

    |m_actual| <= [I_new+1.1e-10 ||j Psi_ref||E*]/K_m,
    K_m>1.9999999764.                                 (14)

The full residual norm is UNKNOWN. Compute its sufficient ceiling
[K_m(2.29e-15)-I_new]/(1.1e-10) from the new bound, but do not assign
that ceiling to the actual residual. In particular (12) is not a metric
variation estimate; differentiating its envelope would be invalid.

The first unsupported mass step is now the COMPLETE j Psi_ref in the
actual dual-energy norm. The first unsupported full-graph step is still
the differentiated common-regulated core composed with the coupled
inverse. kappa_m<=5e-5, the other coupled columns, the original-norm
comparisons and uniform nonlinear-neighborhood control remain unknown.
The .135973<.1599164 screen still covers only its two-fast-input family.
No actual mass, original graph, nonlinear admission/contraction, exact
compatibility, physical matching, validity/noise or stability follows.

RESEARCH DECISION: retain the construction and the coupled adjoint
method. If (13) passes, stop optimizing this complete source moment and
focus on its actual quantum/contact residual. The improvement comes
from a direct value estimate and integrated boundary weight, not a
normalization change or assumed cancellation of an unknown source.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): apply the complete restored
current Hessian to Psi_ref and certify ||j Psi_ref||E* against the ceiling
computed from (14), including the outer phase, density contacts, scalar
and finite-domain returns. This tests the actual mass certificate; the
independent original graph requirement remains afterward.

Remaining to-do list: full coupled adjoint residual and actual mass;
original coupled graph return; physical-norm/neighborhood comparisons;
nonlinear admission/contraction; exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole objective.
