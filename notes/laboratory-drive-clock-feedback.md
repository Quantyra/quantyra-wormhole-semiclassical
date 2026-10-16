# One fixed laboratory waveform and its coupled clock feedback

2026-10-15, E003/S135.3. The laboratory waveform is fixed once using the
central free reference preparation. Different actual preparations then
retime that same waveform through their clocks. This is a stipulated
free minimal-JT coupled system, not a completed four-dimensional source
or quantum-state match. It tests a different question from prescribing
the same global-time angle separately for every preparation.

## Frozen problem and acceptance

Take eta=1/2000, tau=1/4, T_star=pi+tau and the inherited quintic
S(x)=10x^3-15x^4+6x^5 on [0,1], extended by 0 and 1. Put

```
rho0(t)=4*eta*[S(t/tau)-S((t-pi)/tau)],
Theta0(t)=integral_0^t rho0,
d0(t)=4*eta*S(t/tau), Q0(t)=integral_(t-pi)^t d0(s)^2 ds,
g0(t)=[Theta0(t)-Q0(t)]/(4*pi*eta),
w''+w=g0/2, w(0)=53/200, w'(0)=0.
```

The reference is extended through global time 16, including its exact
constant-force continuation after T_star. The laboratory source is
theta_b(u)=Theta0(t_ref(u)), with dt_ref/du=eta*w(t_ref). Actual clocks
have dt/du=eta*z(t), hence F(t)=t_ref(u(t)) obeys F'=w(F)/z.
Use this one source for all six actual initial conditions
z(0) in {1/4,7/25}, z'(0) in {-1/2000,0,1/2000}, F(0)=Q(0)=0.
The full prior source history is zero. It is not an independent-reservoir
or erased correlated-vacuum prescription.

The diagnostic equations, with an optional externally bounded residual r,
are

```
theta(t)=Theta0(F(t)), p=F'=w(F)/z, rho=theta'=rho0(F)*p,
D(t)=rho(t)+D(t-pi), D(t)=0 for t<0,
Q'=rho^2+2*rho*D(t-pi), Q(0)=0,
g=[theta-Q]/(4*pi*eta), z''+z=g/2+r.
```

The constraint insertion has a separate ordinary derivation in the
[minimal-JT review](jt-prescribed-source-constraint.md). The fixed-source
intrinsic-contact cancellation is derived in the
[source-ensemble review](prescribed-source-clock-variation.md).
This note does not extend these results to the missing EFT terms,
nonzero angular modes or full four-dimensional backreaction. The
retiming relation is kinematic once the shared laboratory waveform is
specified; the delayed stress remains in actual global time.

A read-only exploratory diagnostic preceded this acceptance protocol. It
suggested qualitative negative-energy capture at all six corners but
large completion shifts and failure of inherited quantitative margins.
No waveform, preparation corner or stopping rule is changed in response.
Before the new checker, freeze horizon 0<=t<=8, domain stop z=1/200,
and phase stop F=16. Completion is the event F=T_star, not t=T_star.
Run direct variables with DOP853 at rtol=2e-11, atol=2e-13, max step .01;
refine to 2e-12,2e-14,.005. Independently evolve aligned variables with
RK45 at 2e-12,2e-14,.005 and an independently integrated reference.
Uniform common-grid state disagreement must be below 2e-6, completion
time and terminal-energy disagreement below 2e-6. Independent memory
quadrature, energy/work identities and exact-reference control permit
2e-8 absolute error. Analytic defect inequalities allow 1e-8 numerical
rounding allowance. These are numerical diagnostics, not interval-certified
coverage of the preparation box or a certified residual-stability radius.

## Method of steps and continuation

Fix a positive C2 reference w on an open phase domain, a C2 source profile
rho0, and a locally bounded measurable residual r. On the first interval
of length pi the delayed D is prescribed history. The four equations for
(z,z',F,Q) have locally Lipschitz state dependence on z>0 and F inside
that domain. Standard integral Picard iteration, with the measurable
bounded inhomogeneity, gives a unique local absolutely continuous
solution. Its continuation through the next interval uses the already
constructed continuous history D; repeating gives the method of steps.
The alternative identity Q'=D(t)^2-D(t-pi)^2 proves by differentiation
and the initial value that Q is exactly its moving-window integral.

On a finite horizon, a positive lower bound for z, a compact phase range,
and bounded state/source jets exclude the rational singularity and permit
continuation. Blow-up or a domain exit must be reported rather than
continued through z=0. This is a conditional well-posedness statement;
finite numerical trajectories do not prove all preparations avoid exit.

## Exact aligned coordinates and a finite sensitivity criterion

Let a=z-w(F), b=z'-w'(F), q=Q-Q0(F). Then, without linearizing,

```
p=1-a/z,
a'=b+w'(F)*a/z,
b'=[-1+w''(F)/z]*a-q/(8*pi*eta)+r.
```

The leading Theta0(F) forcing cancels. Differentiating a scalar envelope
for Theta0 is unnecessary. Because p>0, the inverse phase s=F(t) is
available while z,w>0. With X=(a,b), it gives

```
dX/ds=A0(s)*X+[a*b/w, -a^2/w]
                    +[0,(1+a/w)*(r-q/(8*pi*eta))],
A0=[[w'/w,1],[-1+w''/w,0]], dt/ds=1+a/w.
```

Thus only a known two-dimensional reference variational equation is
needed for the principal gain. Let Phi be its fundamental propagator.
Variation of constants is an exact Volterra formula for X with the two
displayed nonlinear/input terms. It avoids replacing the entire system
by exp(L*T/min z). No numerical propagator norm is claimed certified here.

For an explicit conditional radius, suppose w>=m0 on phases [0,S],
choose 0<h<m0 and set A=h/(m0-h), K=floor(H/pi)+1,
L=4*K+30*H. The radiation estimate below bounds its normalized force by

```
V_H(A)=(eta/8)*[(8+L*A)*L+32]*A.
G_S=sup_(s<=S) ||Phi(s,0)||,
J_S=sup_(s<=S) integral_0^s ||Phi(s,v)|| dv.
```

If S>=H/(1-h/m0), ||X(0)||<h, |r|<=R, and

```
G_S*||X(0)||+J_S*[h^2/m0+(1+h/m0)*(R+V_H(A))]<h,
```

the first-exit argument in the exact Volterra formula proves ||X||<h
through actual time H. Here the Euclidean norm of [a*b,-a^2] is at
most ||X||^2. A certified bound on G_S,J_S would turn this into a
quantitative residual/preparation theorem; a sampled propagator alone
does not do so. The criterion is not asserted to hold at the six corners.
It also yields |F(t)-t|<=H*A and, if completion is in range,
|t_end-T_star|<=T_star*h/m0 from integrating dt/ds. These distinguish
phase timing, aligned clock error and common-time error.

## Quantitative radiation defect, rather than a value-only bound

Extend F(t)=t before zero for this identity; all source terms there vanish.
Set alpha=|1-p|=|a|/z and B(t)=D(t)-d0(F(t)). Directly from the delay,

```
B(t)=rho0(F(t))*(p(t)-1)+B(t-pi)
               +d0(F(t-pi))-d0(F(t)-pi).
```

The inherited profile has ||rho0||inf<=4*eta and
||d0'||inf<=30*eta. Its last difference is bounded by
30*eta*integral_(max(0,t-pi))^t alpha. Unrolling the disjoint windows gives

```
|B(t)|<=4*eta*sum_(j>=0,t-j*pi>=0) alpha(t-j*pi)
                      +30*eta*integral_0^t alpha,
```

with only finitely many terms. Changing variables v=F(s) in the reference
moving-window integral, including its shifted lower endpoint, then gives

```
|q(t)|<=integral_(t-pi)^t (8*eta+|B(s)|)*|B(s)| ds
                         +32*eta^2*integral_(t-pi)^t alpha(s) ds.
```

In particular alpha<=A on [0,H] implies |B|<=eta*L*A and
|q|/(8*pi*eta)<=V_H(A). The normalized radiation feedback is O(eta)
for a fixed bounded relative retiming defect and horizon. This does not
bound the leading geometry error without the coupled clock evolution.

## Energy, completion and the meaning of capture

For e=z'^2+z^2-g*z, direct differentiation gives

```
e'=2*r*z'-g'*z,
g'=rho*[1-rho-2*D(t-pi)]/(4*pi*eta).
```

These are finite energy/work identities, including the returned cross
term. Consequently |e(t)-e(0)+integral_0^t g'*z|
<=2*R*integral_0^t |z'|. Once F>=T_star, rho=0 even if returned D is
nonzero; Q,theta,g are constant. For r=0, e is then conserved. If
g_final>0 and e_final<0, the autonomous oscillator has positive lower
clock value g_final/2-sqrt(g_final^2/4+e_final). A residual afterward
adds at most 2*R*integral |z'| to the energy ledger. All-future control
with a nonzero residual needs a corresponding all-future assumption.

The diagnostic reports actual completion, its timing shift, finite-horizon
positivity and autonomous capture separately. Common-time errors are
reported both against the central laboratory reference and against each
corner's same-preparation prescribed-global-time baseline
w(t)+(z_init-53/200)*cos(t)+v_init*sin(t). The latter removes the trivial
initial preparation difference. Neither comparison silently reuses a
different completion time as the old fixed-T_star observable.

## Decision scope

The accepted receipt determines whether the exploratory outcome survives
refinement. Qualitative capture, failure of the old -1/25 final-energy
margin, and failure of old common-time correction budgets are distinct
verdicts. No numerical sample proves a whole preparation interval or an
admissible nonzero residual range. The six-corner question is not retuned
if a margin fails. The original prescribed-global-angle theorem retains
its own assumptions. No full S074 pass, physical state preparation,
four-dimensional opening or topology-creation result follows.

The retained acceptance run reproduces all six numerical captures. Actual
source completion is at t=2.42445–6.67360, a shift of -0.96714 to +3.28200
from T_star. The endpoint-dependent autonomous future minimum z is
0.0217529–0.226631. This follows from the exact autonomous formula applied
to numerical endpoints, not interval-certified endpoint enclosures.

All six same-preparation common-time clock correction budgets fail: the
sampled maximum differences on the old [0,T_star] time window are
0.02613–0.05068, above 1/200; on the old [T_star,8] window they are
0.07503–0.16213, above 1/50. The latter window is not the actual
post-source interval for the slow corners. The upper-z three corners
also fail the old completion-energy bound e<-1/25, with
e=-0.02119 to -0.02173. Lower-z corners have e=-0.17298 to -0.17479.
Matched-source-completion energy differences compare each actual endpoint
F=T_star against its same-preparation global-profile baseline at T_star;
all six fail the old |delta e|<1/50 correction budget. The lower-z energy
differences are -0.10844 to -0.10726; upper-z differences are +0.02626
to +0.02635. These compare matched source completion, while the clock
errors above compare explicitly stated common global-time windows.

The largest trajectory-refinement error is 8.54e-10 and independently
aligned RK45 disagreement is 1.53e-10; memory/work residuals are below
3.50e-11 and the same-initial identity control error is 6.22e-15. These
checks support the reported failures, not a uniform rectangle theorem.
The checker's --refresh-report mode reproduces endpoint/baseline reporting
from the accepted receipt without rerunning the six evolutions. Its log
distinguishes this reporting-only review and retains the original
evolution-source hash as well as final source hashes.

Remaining to-do list: certify a useful coupled response/preparation domain
or retain its failed margins; match the source/controller, state, contacts
and nonzero angular sectors; carry the remaining errors into S076.
