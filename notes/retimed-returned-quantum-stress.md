# Exact retimed free-state stress and excitation at finite optical width

2026-10-15, E003/S135.1, supporting S074. The prepared correlated free
fermion state has an exact returned covariance for an arbitrary smooth
laboratory waveform on a positive clock and a fixed complete optical
width P. Point splitting gives its global and laboratory stresses with
an explicit subtraction dictionary. Retiming and width errors have a
finite, all-frequency excitation bound. No interacting or evolving-
geometry matching, physical clock-error bound or opening is inferred.

## Scope and acceptance before execution

Retain N real massless fermions (central charge c=N/2), the initial
correlated NS strip vacuum and reflection phases of
[returned-switch capture](returned-switch-capture.md). Use a fixed
complete optical width P>0, including the retained central interval
and both exterior legs as in [S134.2](exterior-lll-delayed-boundary.md).
The fixed laboratory source is theta_b(u)=Theta0(t0(u)); on an actual
positive clock t=t(u), its global representative is theta(t)=Theta0(F(t)),
F=t0 o u. This is the [S134 laboratory ensemble](prescribed-source-clock-variation.md),
not an affine feedback source. Complete the profile smoothly to constants
in the past/future; assume its global rate has compact support and the
same small-angle spectral branch. All covariance statements concern this
prepared free state, rather than independent reservoirs.

Before running the checker, freeze: independent 65-digit covariance
point splitting at separations 1e-4,1e-5,1e-6 must agree with the
global/laboratory formulas within 2e-8 at the finest separation;
same-map cancellation and nonlinear inverse-clock identities within
1e-45 when evaluated as finite analytic identities; continuum time
quadratures and global work derivatives within 2e-9; Fourier Parseval
partial sums at 64 and 128 harmonics must agree with the time-domain
excitation up to the analytic tail plus 2e-10 rounding. Independent
Fourier quadrature refinements use 512 and 768 Gauss nodes per smooth
source piece and must agree within 2e-10. Proven excitation inequalities
allow relative 2e-8 rounding. Examples use the actual shaped source,
fixed laboratory clock t0(u)=u and actual t(u)=s*u+b*u^3 with
(s,b)=(1.02,0),(1,0.002); P/pi=0.97,1,1.04. These clock and width
changes test the exact state calculation; they are not measured errors
of the RN-AdS candidate. No Lean is used.

## Exact covariance before choosing a clock equation

Write q(t)=theta'(t), q=0 before zero, and
A(t)=sum_(j>=0)theta(t-jP), D(t)=sum_(j>=0)q(t-jP).
At every finite time these causal sums have finitely many terms. With
J=[[0,-1],[1,0]], R(a)=exp(Ja), define

```text
G0,P(d)=I/[4P i sin(pi(d-i0)/(2P))]
                         +iJ/[4P cos(pi(d-i0)/(2P))],
Gout(t1,t2)=R(A(t1)) G0,P(t1-t2) R(A(t2))^T.              (1)
```

Incoming covariance replaces A(t) by A(t-P). The off-diagonal term
is the inherited opposite-endpoint correlation. Successive rotations
and free propagation commute in this one-angle free model, proving
(1) by causal iteration. Its positivity and CAR follow from the
original unitary characteristic construction. It is not an assertion
about a generic mixing interacting exterior.

Point splitting relative to the uncoupled vacuum at the same P gives
the outgoing chiral stress, with N already included,

```text
p_P(t)=N/[8pi] [D(t)^2-(pi/P)D(t)].                       (2)
```

The two bulk chiral stresses are p_P(t-sigma) and p_P(t+sigma-P).
Their spatial integral, relative to that initial vacuum, is

```text
E_P(t)/N=-theta(t)/(4P)+(1/(4pi))int_(t-P)^t D(v)^2 dv,
E'_P(t)/N=[q(t)^2+2q(t)D(t-P)-(pi/P)q(t)]/(4pi).         (3)
```

Indeed D(t)-D(t-P)=q(t), and the last-period integral of D equals
theta(t). No small-angle expansion of the quadratic stress was made.
The returned cross term in (3) is retained for overlapping returns.
For P!=pi, E_P is the full optical-channel energy relative to that
P vacuum. It cannot directly replace the unit-AdS2 JT bulk energy
without the exterior/geometric charge map.

## Laboratory point splitting and the reference dictionary

Let y=t'(u)>0. Canonical endpoint fields transform as half-densities,
so G_lab(u1,u2)=sqrt(y(u1)y(u2))Gout(t(u1),t(u2)). The diagonal
singularity is 1/(2pi i d). Expanding its pullback at separation h,

```text
sqrt(t'(u+h/2)t'(u-h/2))/(t(u+h/2)-t(u-h/2))
 =1/h+h Sch(t,u)/12+O(h^3).
```

The real-fermion energy prescription is -i/2 times the derivative
of the regular covariance. Thus the chiral canonical stress defined
by subtracting the laboratory line singularity is

```text
T_vac,P=-N*pi/(96P^2),
T_lab^can(u)=y^2[T_vac,P+p_P(t(u))]
                                      -N*Sch(t,u)/(48pi). (4)
```

This establishes both the sign and coefficient, c/(24pi)=N/(48pi).
At the same map t(u), subtract the pulled-back uncoupled state at
the same P. Both vacuum and Schwarzian cancel exactly:

```text
Delta_same T_lab=y^2 p_P(t(u)).                          (5)
```

Comparing actual t(u),P and reference t0(u),P0 with one fixed
laboratory normal-ordering is a different subtraction:

```text
T_lab^can-T_lab,0^can
 =y^2 p_P(t)-y0^2 p_P0(t0)
  +y^2 T_vac,P-y0^2 T_vac,P0
  -N[Sch(t,u)-Sch(t0,u)]/(48pi).                         (6)
```

The last line and the vacuum-polarization difference cannot be
discarded by the same-map cancellation. Conversely, covariant bulk
stress includes the geometric Weyl completion. In the standard global
AdS2 frame at P=pi that completion cancels the initial NS Casimir;
its state-independent metric piece is absorbed in the inherited
dilaton convention. Transforming a covariant tensor and transforming
the canonically normal-ordered boundary observable (4) are different
operations. The compensating geometric term restores covariance.
The global JT charge calculation must use its existing renormalized
bulk-frame stress, not add the laboratory Schwarzian a second time.
No geometric completion for an evolving four-dimensional exterior
or a comparison of different physical metrics is supplied here.

At one laboratory boundary, the outgoing and incoming fields have
the same clock map and hence the same Schwarzian. Their difference
is y^2[p_P(t)-p_P(t-P)]. Summing the two ends gives the matter flux
work y^2 E'_P(t). The redshifted charge y E_P instead obeys
d_u(y E_P)=y' E_P+y^2 E'_P; the first term is clock/redshift work.
Neither identity determines an unknown finite source contact.

## Final excitation is an exact nonnegative quadratic functional

Once theta=theta_f and q=0, D is P-periodic with mean theta_f/P.
On the retained branch the final static vacuum shift is

```text
E_vac,f/N=-theta_f/(4P)+theta_f^2/(4pi P).
E_exc/N=(1/(4pi))int_period [D(v)-theta_f/P]^2 dv
       =(1/(2pi P))sum_(n>=1)|qhat(2pi n/P)|^2.          (7)
```

The first equality follows from (3), and the second is Parseval for
the periodized compact rate. This is the continuum quantum-state
excitation; no fermion mode cutoff or all-frequency operator norm
approximation is used. Exact preparation of the final vacuum occurs
precisely when every nonzero period Fourier coefficient vanishes.
After an imperfect retiming, the closed free system retains periodic
local stress and a constant positive excitation energy; it does not
relax to the final vacuum merely because the source has stopped.

For theta=Theta0 o F and inverse clock map T=F^(-1),
qhat(omega)=int Theta0'(v) exp[-i omega T(v)]dv. This
formula retains all frequencies and exposes how clock retiming changes
the spectral zeros even with an unchanged laboratory waveform.

## Quantitative retiming and finite-width bound without a mode truncation

Let P0=pi and let the actual retained smooth source be
q0(t)=4eta[S(t/tau)-S((t-P0)/tau)], tau=1/4, with S the quintic
smooth step. Its support is [0,T0], T0=P0+tau, and theta_f=4eta P0.
Its periodization at P0 is theta_f/P0. For lambda=P/P0 define
q_P(t)=q0(t/lambda)/lambda. This comparison rate has the same final
angle and exactly cancels excitation at width P.

Suppose q=q0(F)F' and q_P are supported within [0,Tc], and on that
interval F'>0 with

```text
|F-t/lambda|<=epsilon0, |F'-1/lambda|<=epsilon1, F'<=B.
R0=||q0||inf=4eta, R1=||q0'||inf=15eta/(2tau),
rho=R1 B epsilon0+R0 epsilon1, M=ceil(Tc/P).
E_exc/N <= M*||q-q_P||2^2/(4pi)
        <= M*Tc*rho^2/(4pi).                            (8)
```

Proof. The source difference is bounded by rho using the global
Lipschitz extension of q0 by zero. Periodization of that difference
has at most M overlapping terms almost everywhere. Cauchy-Schwarz
and integration over one period bound its squared L2 norm by M
times the full-line squared norm. Apply (7). Both retiming and width
enter this one controlled comparison; neither error is assumed to be
the actual RN-AdS residual. S134's positive-clock retiming estimates
can provide epsilon0,epsilon1 only after the actual clock error and
its common interval hypotheses are established.

For a finite observation interval [0,H], put n=floor(H/min(P,P0)).
If ||q-q0||inf<=r and q0 is globally R1-Lipschitz, then
|D_P[q]-D_P0[q0]|<= (n+1)r+R1|P-P0|n(n+1)/2.
Inserting this into (2) gives a direct finite-history null-stress
bound with the linear coefficient pi/P retained. This local estimate
and (8) are useful inputs to a physical remainder calculation, not
already proven clock or area error bounds.

Positivity of the clock alone is insufficient. For F(t)=a t,
with a>T0/P, the entire rate fits inside one return period. Then

```text
E_exc/N=[a ||q0||2^2-theta_f^2/P]/(4pi),
||q0||2^2=16eta^2[P0-50tau/231].                         (9)
```

It diverges with a while every clock is smooth and increasing.
This is an exact retiming of the actual shaped laboratory waveform,
not an arbitrary source substitution. A controlled physical clock
map, rather than monotonicity alone, is the missing quantitative input.

## Certified Fourier tails and checks

For a compact C2 rate, twice integrating by parts gives
|qhat(omega)|<=||q''||1/omega^2. The tail in (7) beyond harmonic K
is at most ||q''||1^2 P^3/(96pi^5 K^3). For the checker clocks
t=T(u)=s u+b u^3, b>=0, the actual rate is q0(u)/T'(u).
Its second global-time derivative is

```text
q''=q0''/(T')^3-3q0'T''/(T')^4-q0 T'''/(T')^4
                                       +3q0(T'')^2/(T')^5.
```

With T'>=s, |T''|<=6bT0, T'''=6b and
int|q0|=theta_f, int|q0'|=8eta, int|q0''|=30eta/tau,
these give an explicit L1 majorant after multiplying by dt=T'du.
The checker uses that majorant, independent source-time Fourier
quadrature, returned time-domain periodization, laboratory covariance
point splitting, the work identity and the analytic retiming bound.
Finite harmonics supplement the exact Parseval/tail proof.

[Maldacena--Qi v3](https://arxiv.org/html/1804.00491v3), appendix C.2,
(C.207)-(C.209), fixes the static real-fermion spectrum, vacuum energy
and c=1/2 normalization; C.1 and (C.210)-(C.211) distinguish the
flat-strip Casimir from the covariant AdS stress and dilaton shift.
The source's HTML C.206 has the previously recorded malformed second
condition; the retained unitary junction convention is used instead.
The time-dependent covariance and its half-density expansion above
are derived directly, rather than imported from that static result.
No novelty or physical creation claim is made.

Remaining to-do list: establish the physical common clock and optical
width/state map with controlled errors; match source contacts and
interactions, evolving geometry and nonzero angular channels; propagate
the complete renormalized stress to the JT constraint and then the
S076 energy/area observables.
