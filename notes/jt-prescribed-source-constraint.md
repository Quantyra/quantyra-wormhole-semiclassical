# Minimal JT constraint and work for a prescribed laboratory source

2026-10-15, E003/S135.2, supporting S074. The symmetric free JT clock
equation follows directly from the bulk metric equation and fixed
metric/dilaton boundary data, even for a time-dependent laboratory source.
The source must be retimed consistently. A midpoint fermion calculation
supplies its work, rather than an independently inserted source energy.
This is a closure of the explicitly minimal two-dimensional model, not
the physical RN-AdS4 action or its complete interacting clock equation.

## Hypotheses and acceptance before execution

Use the unit global AdS2 strip, 0<sigma<pi, the symmetric paired free
Majorana domain and initial correlated strip vacuum of
[the returned-switch calculation](returned-switch-capture.md). N counts
the real fields and is the normalized single-boundary Schwarzian
coefficient after u=N U/phi_r. Work in the leading JT boundary limit with
fixed induced metric and fixed dilaton, and a common positive clock
t=t(u), y=dt/du>0. The source is prescribed in laboratory time:

```text
theta_b(u)=Theta0(t0(u)), h_b=tan(theta_b/2),
theta(t)=theta_b(u(t))=Theta0(F(t)), F=t0 o u(t).
```

The matter action is a free Dirac/Majorana action with that paired
self-adjoint boundary domain. It has no dilaton coupling. Its boundary
interaction is realized by the minimal midpoint contact below. No extra
extrinsic-curvature, dynamical boundary field, boundary dilaton or
unmatched metric functional is included. A fixed-source intrinsic
c-number contact can be retained with its consistent stress/work
accounting as in [S134](prescribed-source-clock-variation.md); it has no
direct embedding variation and is not appended to the clock constraint.

Subtract the same state-independent AdS anomaly as the retained benchmark.
The remaining bulk stress is conserved and traceless. Assume sufficiently
smooth finite-interval stress and the standard asymptotic dilaton expansion
used below. No physical infrared, RN exterior or higher-mode input is
inferred from these hypotheses.

Before checker execution fix: exact symbolic tensor, trace, boundary-series,
Noether and work identities must vanish identically. Numerical checks use
60-decimal arithmetic. Differential tensor, energy and exact work identities
allow absolute error 1e-35. Independent laboratory-energy differences at
steps 2e-5 and 1e-5 allow 5e-8. Midpoint covariance limits at separations
2e-5 and 1e-5 allow 5e-10. The parameter choices are synthetic verification
cases, not physical background or controller assignments.

## Constraint from the metric equation

Let Phi be the anomaly-shifted dynamical dilaton. In the convention of
[MQ section 2.2](https://arxiv.org/html/1804.00491v3), the metric equation
and its trace are

```text
T_ab=g_ab Box Phi-nabla_a nabla_b Phi-g_ab Phi,
(-partial_t^2+partial_sigma^2-2*csc(sigma)^2)Phi=0.
```

For g=csc(sigma)^2 diag(-1,1), direct metric differentiation gives

```text
T_tt=-Phi_ss-cot(sigma)Phi_s+csc(sigma)^2 Phi
    =-partial_sigma[Phi_s+cot(sigma)Phi].                 (1)
```

At the left boundary, the trace equation fixes the first two subleading
terms of the expansion:

```text
Phi=B_L(t)/sigma+A_L(t)+sigma*C_L(t)+O(sigma^2),
A_L=0, C_L=-B_L''/2-B_L/3,
lim_(sigma->0)[Phi_s+cot(sigma)Phi]=-(B_L''+B_L).
```

At the right boundary put v=pi-sigma. The corresponding limit is
+(B_R''+B_R). Therefore the full bulk energy obeys the exact identity

```text
E(t)=integral_0^pi T_tt d sigma
    =-(B_L''+B_L+B_R''+B_R).                            (2)
```

The prescribed source is symmetric, so the two boundary coefficients
agree. The fixed boundary conditions supply their normalization: a cut
with induced ds^2=-du^2/epsilon^2 has sigma=epsilon*y+O(epsilon^3);
Phi|cut=N/epsilon thus requires B_L=B_R=N*y. Consequently

```text
y_tt+y=-E/(2N),
phi_uu=-y^2-E*y/(2N), phi=log y.                        (3)
```

No stationary source or fixed profile in global time was used in (1)-(3).
E is the complete stress integral for the specified free boundary domain,
not merely its positive radiation part. Extra boundary terms that change
the metric/dilaton boundary problem require a fresh derivation.

The same result is the SL(2) constraint. From the two Schwarzian actions,
the pure-gravity global time-shift charge is

```text
Q_0,grav=-2N[y+phi_uu/y]=-2N(y+y_tt)=E.
q_0,matter=-E, Q_0,grav+q_0,matter=0.                    (4)
```

Reflection symmetry makes both nonenergy matter Killing charges zero:
their energy integrands contain cos(sigma), and their momentum integrands
combine the odd T_tsigma with sin(sigma). Identical boundary trajectories
also cancel the two gravitational nonenergy charges. This does not cover
asymmetric perturbations. Adding the interaction term from MQ's already
integrated-out effective action to the full-domain E would double count
that matter response. The static comparison in MQ appendix C.3 agrees
with the small-coupling effective-action charge; it is not a prescription
to add both.

## Laboratory energy and total work

The minimal two-boundary laboratory energy is the Schwarzian expression

```text
H_u=-2N{tan(t/2),u}=N(phi_u^2-2phi_uu-y^2).
```

Using the constraint only after this definition gives

```text
H_u=N(phi_u^2+y^2)+y*E,
dH_u/du=y^2*dE/dt.                                    (5)
```

Indeed differentiating the middle expression gives
2N*phi_u*(phi_uu+y^2)+y_u*E+y*E_u; the first and second terms cancel
using (3) and y_u=y*phi_u. This is the total boundary-energy response
of the minimal constrained system, not an extra matter term to add to
H_u again. The coordinate redshift of E alone is not a total Hamiltonian
when y varies.

For the returned free covariance, define

```text
r(t)=theta'(t), D(t)=sum_(j>=0)r(t-j*pi), D_-=D(t-pi),
p(t)=N[D(t)^2-D(t)]/(8*pi),
E=2 integral_(t-pi)^t p(v)dv
 =N[-theta(t)+integral_(t-pi)^t D(v)^2 dv]/(4*pi),
E_t=N*r*(r+2D_--1)/(4*pi).                             (6)
```

The source is extended by zero into the uncoupled past. Only finitely
many returns occur on a finite interval. Equations (3) and (6), with the
actual theta=Theta0(F), therefore give a causal free mean system while
y remains positive. They do not retain the reference Fourier cancellation
on an altered trajectory without an error estimate.

Here is an independent boundary-interaction work check. For each real
pair of endpoint fields, let J=[[0,-1],[1,0]], R(theta)=exp(J theta),
and use the midpoint convention

```text
gamma=2*tan(theta/2), H_int=-gamma*X,
X=i sum_(species) psi_L,mid psi_R,mid,
psi_mid=M a, M=(1+R(theta))/2=cos(theta/2)R(theta/2).
```

This convention gives the Cayley jump
R=(1-gamma J/2)^(-1)(1+gamma J/2). It fixes a free source normalization,
not a complete RN4 vertex or an interacting controller. The incoming
covariance is the initial correlated kernel rotated by the accumulated
past angle. With d=t1-t2 its unrotated short-distance part is

```text
G0(d)=I/[4*pi*i*sin((d-i0)/2)]+i*J/[4*pi*cos((d-i0)/2)].
```

The midpoint angle has rate D_-+r/2. Its off-diagonal point split has
finite limit

```text
<X_t>=N*cos(theta/2)^2*(1-r-2D_-)/(4*pi),
<X_u>=y*<X_t>.
```

As a static normalization check, the free vacuum energy is
E_vac=N[-theta/(4*pi)+theta^2/(4*pi^2)]. Its derivative with respect
to gamma equals -<X_t> at r=0 and D_-=theta/pi. The midpoint
interaction and the self-adjoint domain describe the same generator.
In a distributional contact representation the kinetic jump and delta
potential cancel by -i(R-1)+i*gamma*J*(1+R)/2=0. Adding -gamma*<X_t>
to the complete-domain stress integral E a second time would double count
the interaction. Its derivative is the work insertion, not an extra energy
to append to (6).

The second relation follows from the canonical half-density on both
endpoint fields. Thus the actual prescribed laboratory coupling does work

```text
P_source=-gamma_b'(u)*<X_u>
 =N*y^2*r*(r+2D_--1)/(4*pi)
 =y^2*E_t=dH_u/du.                                    (7)
```

Both propagation/return work and the prescribed source bilinear are
included. The external source supplies or extracts (7); its apparatus is
not separately modeled. No independent conservation of H_u is asserted
during driving. Once the source is constant, (6)-(7) conserve H_u even
when the returned stress is not a stationary vacuum.

## Covariant charge and normal-ordering terms

E in (2) is a covariant global stress charge after the specified common
AdS subtraction. A coordinate-normal-ordered chiral stress has a Schwarzian
term under t=t(u); the geometric stress restores tensor covariance. In
an incoming-minus-outgoing flux the equal geometric/Schwarzian terms
cancel, leaving the factor y^2. They are not a second source energy.
The Schwarzian H_u in (5) is the gravitational boundary energy with the
retained coefficient N. State-independent renormalizations of that
coefficient or of the extremal dilaton must use the same action/subtraction.

In particular, one cannot infer (3) by redshifting an arbitrary local
source contact as H_contact/y and appending it to E. Such a contact's
metric/source Ward terms must be retained consistently, while the fixed
laboratory intrinsic contact has the zero direct embedding force proved
in S134. Genuine boundary metric, extrinsic, field and controller terms
remain distinct possible changes to this minimal problem.

## Consequence and source scope

For fixed laboratory data, the minimal free oscillator is justified by
the JT field equations rather than selected by rule D. Its source is the
retimed waveform, not the old prescribed global profile or an affine
kappa*z replacement. The total free work is (7). The retained preparation,
vacuum and stability claims require their source/state hypotheses and are
not promoted to the altered protocol merely by this closure.

[MQ](https://arxiv.org/html/1804.00491v3) supplies the JT action and fixed
boundary data in (2.6)-(2.8), the source interaction in (2.12), the
Schwarzian charge conventions in appendix B, and the static free-fermion
stress/dilaton comparison in appendix C. Its small-coupling effective
action (4.26) already integrates out a leading matter response and omits
the general particle production retained here. The direct dynamic
constraint/work calculation above is supplied separately. No new theorem
or physical novelty is claimed from those established ingredients.

The checker derives the tensor from the metric, checks both endpoint
series, tests homogeneous/static/retarded source dilatons, differentiates
the laboratory energy independently, and checks the midpoint covariance
against the physical source-work relation. Its receipt/log retain exact
source identities, runtime, tolerances and terminal status. No Lean or
paid compute is used.

Remaining to-do list: solve and control the retimed free response on the
selected laboratory ensemble; match the physical exterior and boundary
action, state, higher modes, interactions and geometry response; propagate
the complete remainder to the S076 energy, clock and area observables.
