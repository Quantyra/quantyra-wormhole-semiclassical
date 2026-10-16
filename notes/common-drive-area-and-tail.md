# Common-drive area control and the residual after source completion

2026-10-15, E003/S136.4, supporting S074/S076. A controlled driving
interval has a precise continuation consequence, but a small bounded
force for arbitrarily long time does not have the same consequence.
Area control also needs more information than the integrated clock
equation. The statements below are ordinary conditional results; they
do not assert the proposed S136 preparation-domain targets are proved.

## Fixed hypotheses and verification target

Retain the single laboratory source and correlated free initial state
of [S135](laboratory-drive-clock-feedback.md), eta=1/2000 and its central
reference w. The original prescribed-global-source theorem gives
0<w<13/10 for that reference. Suppose the actual source completes by
global time 7 and z>=3/200 throughout driving. Suppose its endpoint has
e_f<=-3/200. These are the new S136 targets, not facts supplied by six
sampled trajectories. Set an additive residual to zero after completion
for the first continuation statement. The later tail statement gives a
separate sufficient assumption for a nonzero residual.

Before checker execution require exact symbolic verification of the
spatial Green kernel, the zero-charge area witness and the resonant
oscillator. Use exact rational arithmetic for source and classical-area
constants, with 60-digit evaluations for three smooth tail examples.
Independent quadratures of the tail response and the Green reconstruction
must agree within 1e-35. The examples verify formulas, not physical
stress states or an assumed four-dimensional error bound. No Lean runs.

## Free stress sign from a conditional driving tube

The actual rate q=theta_t=q0(F)w(F)/z is nonnegative and vanishes
outside a source interval of length at most 7. Since q0<=4eta,

```text
q<=4eta*(13/10)/(3/200)=13/75,
D(t)=sum_(j>=0)q(t-j*pi)<=3*(13/75)=13/25<1.
```

At most three nonzero terms overlap because 7<3pi. The same bound
holds after completion: periodization preserves the support-overlap
count even though D need not become constant. The exact chiral stress
p=N(D^2-D)/(8pi) is consequently nonpositive. With
Q=int_(t-pi)^t D^2 and int_(t-pi)^t D=theta,

```text
0<=Q<=13theta/25,
0<=g=(theta-Q)/(4pi*eta)<=1,
g_f>=12/25>0.
```

These deliberately coarse quantum bounds establish signs, not the
sharper capture margin. They do not replace the S136 dynamics test.

## Exact post-source response and a sufficient tail allowance

After completion q=0, so Q and g are constant in the retained free
model. For the clock equation z''+z=g_f/2+r, the unforced solution
from that actual endpoint has

```text
z_free,min=g_f/2-sqrt(g_f^2/4+e_f)
          >=1/2-sqrt(47/200)>3/200.
```

The square root exists because it is the actual oscillator amplitude.
The lower bound uses e_f<=-3/200 and 0<g_f<=1. It is an all-future
statement conditional on certified endpoint data, not a certification
of numerical endpoint samples.

For t>=t_f the exact difference from that free continuation is

```text
delta z(t)=int_(t_f)^t sin(t-s) r(s) ds,
delta z'(t)=int_(t_f)^t cos(t-s) r(s) ds.
```

Thus a sufficient all-future assumption is
int_(t_f)^infinity |r|<=1/200. It retains z>1/100 for all future
times and bounds both coordinate and velocity errors by 1/200.
This is an L1 tail condition; no L1 bound has been derived for the
physical model. A nonzero stationary correction can instead be included
in a matched final equilibrium, with its remaining tail controlled.
Conserved physical interactions can have stronger stability structure
than an arbitrary external forcing class; that structure must be proved.

In contrast, |r|<=R alone, for any R>0, cannot guarantee all-future
capture. Taking r(t_f+v)=R*cos(v) gives the exact perturbation
delta z=R*v*sin(v)/2 from zero initial difference. It has unbounded
amplitude despite uniformly bounded forcing. This rejects that proposed
inference, not a physical interacting wormhole solution: no such
resonant residual is asserted to arise from the matched action.

## A spatial area bound in the same JT normalization

For the symmetric unit AdS2 problem, set W=sin(sigma)*Phi and
B=lim_(sigma->0)W=Neta*z. The exact metric equation gives

```text
(W_sigma/sin(sigma))_sigma=-T_tt,
W_sigma(t,pi/2)=0,
Phi(t,sigma)=B/sin(sigma)
 +int_0^(pi/2) K(sigma,r) T_tt(t,r) dr,
K(sigma,r)=[1-cos(min(sigma,r))]/sin(sigma),
0<=K(sigma,r)<=1, 0<sigma<=pi/2.
```

Reflect to the other half. For the free nonpositive stress,
J=-int_0^(pi/2)T_tt=Neta*g/2, so Phi>=B-J. If the proposed clock
floor holds through driving and the unforced captured continuation,

```text
Phi>=-(97/200)Neta=-(0.485)Neta.
```

This implication applies to a consistent free JT solution with r=0.
Here Phi is the anomaly-shifted free JT dilaton of the constraint note.
The raw physical area field differs by the known constant contact shift;
its reconstruction must retain that shift in the renormalized area ledger.
An arbitrary residual in the clock equation does not specify the bulk
stress or a compatible dilaton. For a consistent corrected unit-JT
equation, write T_eff=T_free+Delta T, and separately establish
B_actual>=Neta*m_B>0. The same spatial identity gives the sufficient bound

```text
Phi_actual >= Neta*(m_B-1/2)
                 -int_0^(pi/2)|Delta T_tt(t,r)|dr.
```

The effective stress must include every geometric term moved to the
reference-JT side. This is not the covariant physical stress alone.
Finite cuts and their physical boundary charge map require the separate
[physical remainder calculation](physical-clock-remainder.md).

For the retained classical RN-AdS candidate, b^2/G=X*k^10,
Neta=2*k^3 and X>=k^3 for k>=2000. Therefore

```text
Phi_ext,class/(Neta)=X*k^7/4>=k^10/4,
Phi_ext,class-(97/200)Neta>0.
```

This evaluates the leading classical constant against the conditional
free bound. It does not perform the anomaly-to-raw-area transfer, fix the
renormalized extremal contribution or
bound physical quantum, state, matching and geometry corrections to area.
Those remain separate terms of the S076 ledger.

## Why the integrated clock error cannot replace an area ledger

There is a simple consistent metric-equation witness. In the unit strip,
add the static reflection-symmetric dilaton
Delta Phi=-A*sin(sigma)^4, A>0, and define Delta T_ab by the JT metric
equation. It is smooth, has zero leading boundary coefficient and zero
boundary charge flux, hence

```text
int_0^pi Delta T_tt d sigma=0,
Delta B=0,
Delta Phi(pi/2)=-A.
```

The clock equation and its integrated residual do not see this change,
while the area at the center changes without bound. The constructed
stress is conserved by the constant-curvature JT identity and generally
has a nonzero trace. No free fermion or physical matter realization is
claimed. It proves that the integrated-error hypothesis alone is
mathematically insufficient: one also needs stress-profile, trace and
boundary information or equivalent control from the full field equations.
The original free traceless covariance has stronger restrictions and
cannot be silently enlarged to include this witness.

Remaining to-do list: certify or reject the S136 driving/endpoint targets;
derive physical stationary/tail control and a spatial stress/metric/area
bound with the matched boundary data; complete the S074-S079 gates.
