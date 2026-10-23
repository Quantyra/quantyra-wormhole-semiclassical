# Fixed-clock bolt test on the current finite seed

2026-10-22, E003/S169. Ordinary argument frozen before supporting checks.
This closes the selected seed-source feasibility test. It does not close
the full coupled finite-domain solve or reopen the unsigned estimate chain.
The stable semiclassical wormhole objective remains incomplete.

## Scope and unchanged theory

Use the exact finite seed (7) of
[the current product calculation](recentered-product-and-cold-residual.md),
with its actual product radii B,A, x=lambda B^2, C=B^2/A^2=1+6x+d,
v=1+h, R=Bv, prescribed kappa and alpha0=0. In t=r/B,

    F=1+x t^2-2a/t+(1+3x+d)/t^2-2d/t^3+d/t^4,
    a=1+2x+delta/2, delta=v F_ext(v),
    B^2 R4=-12x+d(4t^-5-6t^-6).

Retain c_R=10^22, G, q=2000^4, g_R=2000^-4, lambda0, mu0=1/b0,
the positive prepared temperature and paired HH domain. R is the bolt
radius; R4 denotes scalar curvature. There is no scale or norm change.

Use the covariant spherical potential V0 of
[the reference action](stationary-centered-covariant-seed.md), now at
these fixed B,A,lambda. This is an exact addition/subtraction:

    S0=S_EH+EM+Gamma_ref,
    Gamma_rest=Gamma_D+c_R integral R4^2-Gamma_ref.

Boundary terms must be repartitioned consistently. The bulk reference
source has rho_ref+p_ref=0. This split defines a FIRST response to the
evaluated complete seed source; it is not an order reduction of the
unknown full equations. An arbitrary off-shell quantum variation cannot
later be substituted here without its metric, density and Ward contacts.

## Transfer of the evaluated quantum lapse, with hypotheses

The [complete exterior source proof](stationary-centered-exterior-source.md)
is an argument on this same rational family, not only a declaration
about its old numerical center. Its required guards and the current
ones are:

| Guard | Required | Current |
| --- | --- | --- |
| x | >=1000 | >8e9 |
| abs(d) | <1e-12 | <8.6e-19 |
| abs(d)/C | <1e-24 | <1.71e-29 |
| h | [2000^-10/8,2000^-10] | same positive-temperature bracket |
| log(1/h) | <80 | unchanged |
| C/q, first angular mass | <.0035, >18 | unchanged |
| abs(log(mu0 B)) | old <1e-21 | <1e-8; recheck finite restoration below |

The degree-six polynomial, derivative guards, complex cold-chart margin,
non-bolt cap derivatives, geodesic outer boundary and its balanced image
coefficient are precisely those in the exterior proof. No hypothesis
requires the old zero-added-c product equations once the geometry guards
hold. The scale exception is substantive and must not be ignored:

    Fang=4 gamma_E-2 H_q+4 log(mu0 B),
    a0 bracket=zeta_h'(-1)+(1+2log(mu0 B))zeta_h(-1).

Since log(q)<32, H_q<33, 0<gamma_E<1, abs(Fang)<71<100.
The existing abs(zeta_h'(-1))<1000q^2 and
abs(zeta_h(-1))<=q^2/6 give the same a0 restoration <167q^2
(rank two, divisor 4pi). The a0 term has zero spatial derivative.
Thus the original local-contact budgets survive without resetting mu0.
The cold/non-bolt derivative subtotal remains <3e14q^2 and its common
K=1e15 bound remains valid for Theta'=partial_t(B^2 tau_hat,mass).

This uses actual-domain localization, not a vacuum replacement. The
[paired realization](stationary-reference-cold-bolt.md), equations
(13)-(14), has Dsharp skew-adjoint, its prescribed complementary square,
||(Dsharp+M)^-1||<=1/M and the separated two-cutoff resolvent identity.
It requires the exact balanced projectors, regular HH bolt, cap agreement
and common angular subtraction. The current rational seed preserves
these hypotheses. The [non-bolt proof](stationary-reference-nonbolt-stress.md)
retains the reflected boundary layer and its Weyl cancellation; a raw
bounded trace at infinity is NOT substituted. These are retained ordinary
continuum arguments. This check is not independent machine certification
of those kernel theorems, nor an extension to a neighborhood of metrics.

For the evaluated massive source, HH regularity and the Ward equation give

    J(t)=integral_v^t F(s)Theta'(s)/s^2 ds,
    partial_t alpha_mass=-(G/B^2)J/(t F^2), alpha_mass(infinity)=0.

F is increasing and F>=.99x(t-v)(t-1). The existing integral estimate is

    ||alpha_mass|| <= [80*1.001/.99] K/2000^8 <3.2e-10.

The scale factor needed here still holds: the current RN reference obeys
lambda B_RN^4=piG q^2/[g_R^2 F_R(3+1/x_RN)], F_R<1; the product
sphere-square shift is <5.5e-23. Therefore

    Gq^2/(lambda B^4)
      <(3+1/(8e9))/(pi*2000^8)*(1-5.5e-23)^-2
      <.956/2000^8 <1.001/2000^8.

The separately retained LLL lapse is <3e-28. The complete Weyl lapse
(15) in the exterior proof, including its d terms, is <2*2000^-10.
The sum is <3.3e-10. Gamma_ref and the vacuum have rho+p=0.
We transfer this evaluated lapse, not the old complete energy-derivative
or global C5 bounds. Their local R2/reference cancellations would need
a separate calculation before those old constants could be reused.

## Exact additional R2 lapse

Let H_ab=2R4 Ric_ab-g_ab R4^2/2+2(g_ab Box-nabla_a nabla_b)R4.
For the alpha0=0 seed, Ric^r_r=Ric^t_t and

    H^r_r-H^t_t=-2F partial_r^2 R4.

The effective stress of +c_R integral R4^2 on the source side is -2c_R H.
After sphere integration,

    rho_R2+p_R2=16pi c_R r^2 F partial_r^2 R4,
    partial_r alpha_R2=16piG c_R r partial_r^2 R4,
    alpha_R2(r)=16piG c_R[r R4'-(R4-R4(infinity))]
               =w_c d[-24t^-5+42t^-6],
    w_c=16piG c_R/B^2.

The boundary condition alpha_R2(infinity)=0 fixes the integration
constant. The large constant-curvature contribution cancels exactly in
this difference; dropping the derivative terms would miss the whole
answer. For f(t)=-24t^-5+42t^-6, t>=1,
f'=(120t-252)/t^7, f(1)=18 and min f=-4/(21/10)^5.
Consequently

    sup abs(f)=18,
    integral_1^infinity abs(f')dt=18+8/(21/10)^5<18.2,
    abs(w_c d)=(64piG c_R/A^2)*abs(d/C)/4<5.045e-38,
    sup abs(alpha_R2)<9.081e-37,
    integral_R^infinity abs(alpha_R2')dr<9.182e-37.

Together with the quantum source, sup abs(a1)<3.3e-10+9.081e-37.
This uses a lapse value bound; no unproved total-variation bound for the
whole quantum source is inferred from it. The outer clock is included, not chosen to cancel
the bolt source. Pure product or bolt values alone cannot control it:
a smooth conserved effective source supported away from the bolt may
have rho(R)=0 but arbitrarily large integral (rho+p)/(rF).

## Executed prepared-temperature boundary functional

The reference constraints for the COMPLETE evaluated joint seed source
are the same exact first-order boundary problem as
[the linear reconstruction](stationary-centered-first-order-geometry.md):

    (r f1)'=-2G rho, a1'=G(rho+p)/(rF), a1(infinity)=0,
    f1(R)+2kappa s_R=0,
    f1'(R)+F''(R)s_R+2kappa a1(R)=0,
    D=F''(R)/2+kappa/R,
    s_R=[G rho(R)/R-kappa a1(R)]/D,
    delta mu_infinity=kappa R s_R+G integral_R^infinity rho dr.

The angular equation follows from the joint orbit/radius Ward identity.
Mass and clock are correlated by these equations, not free seam data.
The quantum massive/Weyl tails must be combined before mass integration;
the extra local R2 and reference pieces have integrable rational tails.
This statement does not compute the full-action Hamiltonian.

The complete current horizon source already satisfies

    sigma_h=abs(G rho(R)/(kappa R))<2e-12.

It includes the exact local R2/product cancellation and actual HH source
from the current finite-bolt proof. The same rational seed gives
abs(D A^2-1)<7h and abs(h/(kappa A^2/B)-1)<12kappa A^2/B;
thus kappa/(D B h)<1.001. Substitution, including the outer clock, gives

    abs(s_R)/(B h)
       <=1.001[2e-12+3.3e-10+9.182e-37]
       <3.33e-10 <2e-9.

This is an executed bound on the actual finite seed's classical-reference
first correction. It is NOT the full current-action Newton correction,
a finite Hessian inverse, a weighted C5 bound or a nonlinear solution.
It closes the bolt VALUE test for that first correction, with an unused
cold value allowance >1.667e-9. An approximate solution/residual and an
appropriate inverse must bound the full correction's DIFFERENCE by that
amount before the 2e-9 physical cold target is certified. No such error
bound is currently proved. Even the first reference correction has not
been shown to satisfy the full physical higher-order boundary phase.

## Feasibility decision and the next calculation

Together with the exact signed massive slope cancellation, this finite
fixed-clock test makes continuation of the current construction plausible
at the source-value level. It removes the need for a 300-fold unsigned
source improvement to decide this particular bolt value. It gives no
factor-of-six control of unknown Hessian feedback merely because such
headroom exists. High spatial frequencies invalidate perturbing the full
quantum Hessian by a second-order inverse on an unrestricted domain.

Retain the construction; park generic unsigned source optimization. The
single next calculation is the full finite-domain matched linear BVP,
using this first correction as a candidate and evaluating its defect
with the ORIGINAL full Hessian. Its boundary form must fix the clock,
temperature, correlated mass, scalar fast phase and paired matter map;
its Ward/range compatibility must be checked before inversion. Obtain
both the bolt difference <1.667e-9 and the physical weighted C5 return
needed for Y0<=5e-9. A product L2 inverse <.51 is not that return.
This BVP/defect calculation remains PROPOSED AND UNCOMPLETED.

Boundary realization and regularity are the first unresolved obstacles.
A generic N2-to-raw-C5 upgrade does not exist. The actual projection is
identity until a range/compatibility theorem defines another one. The
quantum source and its first two metric derivatives on a preserved ball
are still required; spatial derivative estimates do not provide them.
A small stationary Euclidean correction also proves neither all-sector
retarded stability nor physical semiclassical validity. The near-Planck
symbol crossover retains its earlier limitation. No physical mouth
matching, preparation, opening or transport claim is made.

Remaining to-do list: full matched finite-domain BVP and weighted defect;
preserved nonlinear solution, physical matching, semiclassical validity,
all-sector stability and the remaining stable-wormhole requirements.
