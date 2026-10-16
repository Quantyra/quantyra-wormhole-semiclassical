# Prescribed laboratory sources and the clock variation

2026-10-15, E003/S134, supporting S074. An intrinsic source-only contact
has exactly zero direct clock force when the conformal boundary metric and
the source as a function of laboratory time are fixed. Holding its profile
as a function of global time fixed is a different physical protocol. This
corrects the inference in [the allocation note](source-clock-allocation.md)
and [the rule-D window](source-clock-window.md) that varying the clock alone
selects rule D. Their ODE computations retain their conditional meanings;
they do not establish that insertion rule for a prescribed laboratory source.

## Ensemble, comparison and acceptance before execution

Let u be the fixed, normalized laboratory boundary time, with prescribed
metric gamma_b=-du^2+ell^2 dOmega^2. A constant rescaling to the physical
time U only changes the fixed coefficients below. Let t=t(u) be a smooth
strictly increasing bulk-clock map, y=t'>0. Consider just the intrinsic
source-only term

```text
I_s[theta_b,gamma_b]=A integral_J (theta_b'(u))^2 du,
A=2 pi zeta_s ell^2,                                      (1)
```

with finite real A. Neither theta_b nor gamma_b is varied with t in the
prescribed-laboratory-source ensemble. No matter field, extrinsic curvature
or fluctuating controller occurs in (1). The result also holds for an
arbitrary prescribed coefficient A(u). It does not cover an actual metric
variation or the finite induced metric of a moving throat cut.

The hardest distinction is between varying a dynamical embedding at fixed
laboratory data and varying those data with that embedding. For comparison,
retain the same intrinsic action but hold Theta(t) fixed while defining
theta_b(u)=Theta(t(u)). Derive both variations before proposing a clock
equation. A source-only charge accounting convention cannot replace that
action comparison.

Before running the checker, fix exact symbolic cancellation, the Euler
variation, and the rational counterexample below. Numerical checks use
60-decimal arithmetic; action/finite-pullback identities allow 1e-45 absolute
error. Independent action differences must agree with first derivatives
within 1e-14 and pure/mixed second derivatives within 1e-12. The initial
steps were 2e-4 and 1e-4. The first invocation passed the symbolic proof but
failed the first-difference tolerance at 2e-4, with error 1.07235e-13;
pure/mixed errors were already below tolerance. Final steps are refined
to 2e-5 and 1e-5, with no tolerance relaxation. The checker retains the
original-step errors and checks their factor-four reduction on halving
the step, supporting the central difference's O(step^2) truncation diagnosis.
The original failure, runtime and source hashes remain in the raw log.
The nonzero witness must exceed 1e-6 at the declared synthetic
parameters. Retiming derivative differences have tolerance 1e-16 and all
analytic inequalities allow 1e-40 rounding allowance. Examples are smooth
on the closed interval [0,1], with endpoint-preserving variations; they do
not construct a physical controller or a Lorentzian quantum state.

## Exact cancellation at fixed laboratory data

Define the global-coordinate representative for each clock by

```text
Theta_tmap(t)=theta_b(u(t)).
Theta_tmap'(t(u))=theta_b'(u)/y(u),
I_s=A integral_J y^2 [Theta_tmap'(t(u))]^2 du
   =A integral_(t(J)) y(u(t)) [Theta_tmap'(t)]^2 dt.        (2)
```

These are finite changes of variables, with the same action (1). They do
not assert that a physical clock embedding is pure gauge. For a family of
embeddings with the same fixed theta_b and gamma_b, the right hand side of
(2) equals (1) identically. Hence every first and second clock variation is
zero, including along a nonlinear smooth family. This is independent of
the value or sign of zeta_s.

Infinitesimally let delta t=h, delta y=h', and write r=Theta_tmap'. The
required representative variation at fixed global-coordinate argument is

```text
delta Theta_tmap(t)=-r(t) h(u(t)),
delta Theta_tmap'(t)=-r'(t)h(u(t))-r(t)h'(u(t))/y(u(t)).
delta[y r(t(u))]=h' r+y[delta Theta_tmap'+r' h]=0.         (3)
```

Holding Theta fixed instead of making (3) changes theta_b. It is not the
same fixed-source variational problem written in another coordinate.

An intrinsic local source functional can still shift the renormalized
boundary stress, source one-point contact terms and work accounting under
their defining metric/source variations. At prescribed source and metric
it multiplies the matter evolution by a common phase, so normalized matter
states and correlators are unchanged. Its zero direct embedding force is
not a claim that its stress tensor vanishes. A consistent gravitational
charge convention must treat such contact contributions and the boundary
action together; adding a contact energy to an otherwise unchanged bulk
constraint is not derived by (2).

## The different fixed-global-profile problem

If a source implementation enforces theta_b(u)=Theta(t(u)) while the clock
is varied, its intrinsic action is instead

```text
I_D[t;Theta]=A integral_J y^2 r(t)^2 du, r=Theta',
DI_D[h]=2A integral_J [y r^2 h'+y^2 r r' h] du,
       =-2A integral_J [r^2 t''+y^2 r r'] h du              (4)
```

for h vanishing at the endpoints. The local Euler derivative in (4) is
generally nonzero. Its mixed second variation is

```text
D2I_D[h,k]=2A integral_J [r^2 h'k'
          +2y r r'(h k'+k h')
          +y^2((r')^2+r r'')h k] du.                     (5)
```

Equations (4)-(5) include the source composition; merely holding r as an
independent number while differentiating a scalar energy is insufficient.
Even in this clock-dependent implementation they do not establish the
rule-D oscillator constraint. The controller, energy/work exchange and
gravity charges require a complete action and constraint derivation.

An explicit witness uses the same quintic smooth-step polynomial as the
retained switch, on a unit window:

```text
Theta(t)=omega*(10t^3-15t^4+6t^5),
t_s(u)=u+s*h(u), h=u^4*(1-u)^3.
```

Both source endpoint rates vanish. The clock variation and its first two
u derivatives vanish at both endpoints. For small s, y_s>0. At s=0,

```text
I_s[theta_b=Theta(u)]=(10/7)A omega^2,
d_s I_s|_(fixed theta_b)=d_s^2 I_s|_(fixed theta_b)=0,
d_s I_D|_(fixed Theta)=10 A omega^2/1001.                  (6)
```

To verify the last identity directly, use r=30 omega u^2(1-u)^2 and
r'=60 omega u(1-u)(1-2u) in -2A integral r r' h. The beta integral
is 1800 omega^2 [B(8,7)-2B(9,7)]=-5 omega^2/1001, giving (6).
The nonzero result is not supplied by a different value of zeta_s: it
comes from a different source prescription. For the second direction
k=u^3(1-u)^4, the endpoint jets likewise vanish and |h'|,|k'|<=1/4;
therefore t=u+s h+v k remains monotone for |s|,|v|<=1.

## Constructive open-loop continuation and exact retiming

One available laboratory ensemble is concrete: freeze a reference monotone
clock t0(u), a reference profile Theta0, and prescribe once and for all

```text
theta_b(u)=Theta0(t0(u)).                                (7)
```

For an actual monotone clock t(u), put F(t)=t0(u(t)). If both clocks start
at zero, their exact relation and actual global-time source are

```text
F(0)=0, F'(t)=y0(u(t))/y(u(t)),
theta(t)=Theta0(F(t)), theta'(t)=Theta0'(F(t))*F'(t).
```

When y0(u)=eta*z0(t0(u)) and y(u)=eta*z(t(u)) with the same eta>0,

```text
F'=z0(F)/z(t).                                          (8)
```

Thus the reference waveform in global time is recovered on the reference
trajectory; on another trajectory it is retimed. This exact prescription
is not the affine feedback law r0(t)+kappa*z(t). It does not require
observing the actual clock to implement the fixed laboratory waveform (7).
Its construction does not prove that the physical laboratory or boundary
coupling realizes the required ensemble or prepared state.

There is a quantitative finite-interval estimate. Suppose 0<=t<=T, z>=m>0,
z0 is positive and C1 on a common compact interval containing both t and
F(t), |z0'|<=L, z0<=M0, and |z(t)-z0(t)|<=delta. Then

```text
E(t)=delta*(exp(L*t/m)-1)/L, L>0; E(t)=delta*t/m, L=0,
|F(t)-t|<=E(t),
|F'(t)-1|<=(L*E(t)+delta)/m.                            (9)
```

Indeed the difference equation gives |(F-t)'|<=L|F-t|/m+delta/m;
Gronwall proves (9). If |Theta0'|<=R0 and |Theta0''|<=R1 on that
common interval, the source and rate deviations satisfy

```text
|Theta0(F)-Theta0(t)|<=R0*E(t),
|Theta0'(F)F'-Theta0'(t)|
 <=R1*(M0/m)*E(t)+R0*(L*E(t)+delta)/m.                  (10)
```

These are explicit retiming estimates while the clocks remain in the
declared positive domain. They are not clock-solution bounds, and supply
no state, stress or area budget. Lower-bound loss or leaving the reference
profile interval invalidates their hypotheses.

## Correction and physical boundary

The source-only heat coefficients and fixed-source Ward identities are
retained. No finite coefficient is assigned. What changes is the selection
claim: varying y at fixed conformal metric and fixed theta_b does not force
the rule-D insertion. The statement that the zero source-only correction
sets zeta_s to zero or removes it by hand is also incorrect for this
ensemble. The positive/negative rule-D windows remain numerical properties
of their stipulated ODEs, not bounds on a laboratory source counterterm.

This cancellation applies to the intrinsic zeta_s contact (1). It does
not remove an extrinsic-curvature contact c_theta*I_theta, whose geometric
data can vary, or a source contact at a changing finite induced metric.
It does not remove the field-dependent paired fermion bilinear or the
quantum determinant. Its boundary half-density normalization is retained
as in [boundary-vertex-matching.md](boundary-vertex-matching.md). Nor does
it remove genuine controller degrees of freedom, an interacting state,
retarded boundary response, unequal collars or gravitational constraints.

The two-dimensional throat radius and sphere radius also need not agree
for the retained RN-AdS candidate: R2=b/sqrt(C), C=1+6b^2/ell^2. S133's
equal-radius product is explicitly a diagnostic, not that matched metric.
Its cutoff argument supplies no coefficient identification across those
boundaries. Static exterior geometry, boundary-clock normalization and a
common subtraction still require the exterior-to-throat map.

The next well-posed continuation for (7) is to use (8) in the actual
boundary source, solve its field/state response with the matched exterior
return operator, and derive the total action/clock constraints in that
same laboratory ensemble. The fixed-Theta alternative requires its own
implemented controller and work terms. No new physical opening, stability,
transport, preparation or creation is inferred, and no Lean is invoked.

The complementary S134 calculations treat
[static RN-AdS finite-cut geometry](rn-ads-finite-cut-matching.md) and
[the projected LLL exterior delay](exterior-lll-delayed-boundary.md).
Their separate geometric, free-field and clock hypotheses remain in force;
neither turns the present source prescription into a complete interacting
boundary/controller model.

The checker and receipt retain exact symbolic identities, independent
displaced-action and finite-change-of-variable checks, source hashes,
runtime and terminal logs. Only elementary chain rules, integration by
parts and Gronwall are used; no new source theorem or novelty is claimed.

Remaining to-do list: match the actual boundary geometry and source ensemble;
derive the retimed field/state response and complete controller/gravity
constraint if applicable; match extrinsic, quantum and remaining S074 terms;
propagate the full response to the physical clock and area observables.
