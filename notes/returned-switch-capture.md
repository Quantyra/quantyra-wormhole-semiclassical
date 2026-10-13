# Returned fermion stress and capture during a longer switch

S072, E003, 2026-10-13. **The period-averaged switch is compatible with
opening in the symmetric free-Majorana JT benchmark.** Its exact
returned-field covariance reaches the final matter vacuum, and an
explicit redshift and initial-velocity range keeps the boundary clock
positive throughout switching and forever afterward. A conditional bulk
area bound and a geometric probe channel follow from the same stress.

This resolves the duration objection in the previous
[electric return calculation](electric-return-cancellation.md) at the
free-junction level. It does not equate its gauged neutral coupling `h`
to the free-junction angle `theta`, or supply the interacting
four-dimensional state, stress or preparation. The prepared entangled
pair and an existing bridge remain inputs.

## Exact covariance with overlapping returns

Keep the normalization, initial strip vacuum, anomaly subtraction and
real-field count `N` of [switching-radiation-capture.md](switching-radiation-capture.md).
The strip has width `pi`, the reflection is `J=[[0,-1],[1,0]]`, and the
junction is `R(theta)=exp(J theta)`. Extend `theta(t)=0` before zero.
Propagation and junction rotations commute. Iterating the boundary
relation gives the accumulated outgoing rotation and its derivative

```
Theta(t)=sum_(j>=0) theta(t-j pi),
D(t)=Theta'(t)=sum_(j>=0) theta'(t-j pi).
```

Only finitely many terms occur at each finite time. The outgoing
two-point function is

```
G_out(t1,t2)=R(Theta(t1)) G_0(t1-t2) R(Theta(t2))^T,
G_0(d)=I/[4pi i sin((d-i0)/2)]+iJ/[4pi cos((d-i0)/2)]. (1)
```

The incoming covariance has `Theta(t)` replaced by `Theta(t-pi)`.
Constant propagation rotations commute with `G_0` and cancel. The
construction retains the prepared opposite-endpoint correlation; it
does not treat the two boundaries as independent reservoirs.

Applying the same point splitting to (1), the outgoing stress per
chiral component is

```
p(t)=N[D(t)^2-D(t)]/(8pi),
T_--(t,sigma)=p(t-sigma),
T_++(t,sigma)=p(t+sigma-pi),  T_+-=0.                  (2)
```

In particular the *net* boundary work has a returned cross term:

```
E'(t)/N = [theta'(t)^2+2theta'(t)D(t-pi)-theta'(t)]/(4pi),
E(t)/N = [-theta(t)+Q(t)]/(4pi),
Q(t)=integral_(t-pi)^t D(v)^2 dv.                     (3)
```

Indeed `D(t)-D(t-pi)=theta'(t)` and the initial energy is zero.
The spatial integral of (2) gives the same energy. Omitting the cross
term for a switch longer than `pi` can give a spurious negative
excitation above the final vacuum.

For any completed smooth switch, the last period of `D` has mean
`theta_f/pi`. The static vacuum energy from
[Maldacena and Qi, appendix C.2](https://arxiv.org/html/1804.00491v3) is
`E_vac/N=-theta_f/(4pi)+theta_f^2/(4pi^2)` on the small-angle branch.
Consequently the exact final excess energy is

```
E_exc/N = (1/(4pi)) integral_period [D(v)-theta_f/pi]^2 dv
        = (1/(2pi^2)) sum_(n>=1) |F_theta'(2n)|^2 >=0.  (4)
```

This continuum formula uses the common point-split vacuum subtraction.
It is not the unsubtracted energy of an oscillator cutoff. Fourier
uniqueness gives precisely the cancellation-duration condition derived
in the preceding note.

## The shaped switch reaches the matter vacuum

Fix `tau=1/4` and `0<eta=theta_f/(4pi)<=1/2000`. Let

```
S(x)=10x^3-15x^4+6x^5 on [0,1], extended by 0 and 1,
I(t)=integral_0^t S(v/tau)dv, with I(t)=0 for t<=0,
theta(t)=4 eta [I(t)-I(t-pi)],  T_star=pi+tau.          (5)
```

This is the monotone C3 angle obtained by convolving the normalized
quintic derivative with one uniform period. The return sum telescopes:

```
Theta(t)=4 eta I(t),    D(t)=4 eta S(t/tau).            (6)
```

Thus `D=4eta` after the short smoothing interval `tau`. Once this
outgoing state has filled the strip, at `T_star`, (4) gives exactly
zero excess above the final free-fermion vacuum. This is coherent
preparation, not relaxation or discarding quantum noise.

The covariance supplies an independent state identification. For
outgoing times after `tau`, `G_out=R(4eta d)G_0(d)`. In the `J`
eigenchannels its positive frequencies are `1/2+4eta+2n` and
`3/2-4eta+2n`, `n>=0`, with the unchanged vacuum occupations. These
are the small-angle static towers in MQ C.207. They remain positive
in our interval. Equation (1) propagates this stationary vacuum
through the whole strip by `T_star`; no finite frequency truncation
is used in that assertion.

Put `a=theta/theta_f`, `r=Q/(4pi eta)`, and `g=a-r=-E/(N eta)`.
Equations (3) and (6) give

```
r(t)=(4eta/pi) integral_(t-pi)^t S(v/tau)^2 dv,
0<=r(t)<=4eta,  r(T_star)=4eta,
0<=(1-4eta)a<=g<=a<=1,
g_f=1-4eta.                                           (7)
```

Here `r` is monotone: its derivative is the nonnegative difference
of two squared increasing steps. Also `D<=4eta<1`, so the complete
bulk null stress (2) is nonpositive at all times.

## A capture theorem including preparation errors

Use the inherited common boundary time `u`, `y=dt/du=exp(phi)`, and
`y_0=eta z_0`. The full matter charge constraint becomes particularly
simple in global time, because `y_t=phi_u`:

```
y_tt+y=-E(t)/(2N),
z_tt+z=g(t)/2,  z=y/eta,
z(0)=z_0,  z_t(0)=v_0=phi_u(0)/eta.                   (8)
```

The relation to the nonlinear boundary equation follows from
`phi_uu=-y^2-E y/(2N)`. The two nonenergy matter charges cancel by
reflection symmetry as in the earlier proof. The construction uses
global time until positivity of `y` is established, then defines
`u=integral dt/y`; it does not assume that the clock reaches the end
of a long switch.

For every

```
0<eta<=1/2000,  1/4<=z_0<=7/25,  |v_0|<=1/2000,       (9)
```

the switch completes, and

```
z(t)>1/100 during the ramp,
e_f=z_t(T_star)^2+z(T_star)^2-g_f z(T_star)<-1/25,
z(t)>1/25 after the ramp.                             (10)
```

Proof. First omit `r` and the initial velocity and replace `a` by
the unsmoothed ramp `a_0(t)=min(t/pi,1)` for `t>=0`. Its solution on
`[0,pi]` is

```
z_rect(t)=z_0 cos t+(t-sin t)/(2pi).
```

It has one interior minimum at `t=2 atan(2pi z_0)`, of value
`atan(2pi z_0)/pi-z_0`. This minimum decreases on (9), since its
derivative is `2/[1+(2pi z_0)^2]-1<0`. At the upper endpoint,
`atan(14pi/25)>pi/3`, using `pi>25/8` and `sqrt(3)<7/4`.
Hence `z_rect>1/3-7/25=4/75`. On `[pi,pi+tau]`, its continuation
is `1/2-z_0 cos(t-pi)+sin(t-pi)/pi>=1/2-z_0`, also above this bound.

The shaped ramp is the convolution of `a_0` with the normalized
pulse `p_tau=S'(t/tau)/tau`, whose mean delay is `tau/2`.
The Lipschitz bound on `a_0` gives `|a-a_0|<=tau/(2pi)`. For all
`t<=pi+tau`, `integral_0^t |sin v|dv<=2+tau^2/2=65/32`.
Variation of constants therefore bounds the shaping correction by
`65/(512pi)<65/1600`, the `r` correction by `65/32000`, and the
initial-velocity correction by `1/2000`. Thus

```
z > 4/75-65/1600-65/32000-1/2000 > 1/100.              (11)
```

For the final oscillation, integration by parts in (8) gives
`z=g_f/2+Re(C exp(it))`, with

```
C=z_0-i v_0-(1/2) integral_0^T_star g'(t)exp(-it)dt.
```

Symmetry of `p_tau` gives `p_tilde(1)=b exp(-i tau/2)`, `0<b<=1`.
The `a'` contribution has squared amplitude
`z_0^2+2z_0 b sin(tau/2)/pi+b^2/pi^2<=127/625`.
Since `sqrt(127)<113/10`, while `integral r'=4eta`,

```
|C| < 113/250+2eta+|v_0| <=907/2000,
g_f/2 >=499/1000,
z >=(998-907)/2000 >1/25,
e_f=|C|^2-g_f^2/4 < -173355/4000000 < -1/25.          (12)
```

The solution has no finite-time clock endpoint and is bounded above
by `z<13/10` throughout: during switching use
`z<=z_0+|v_0|+(1/2)integral|sin|`, and afterward use (12).
At fixed final matter charge/coupling, a subsequent reduced mechanical
energy increase smaller than `N eta^2/25` retains negative energy.
This tolerance does not cover arbitrary incoming matter or nonspherical
disturbances. The matter is in its vacuum, but the boundary geometry
generally oscillates rather than approaching a static ground state.

## Bulk area and the probe channel

The complete retarded solution from [driven-bulk-geometry.md](driven-bulk-geometry.md)
applies to the new outgoing stress `p` in (2):

```
H(x)=-integral_0^x [1-cos(x-v)]p(v)dv,  H(x)=0 for x<=0,
Phi=N eta[z_0 cos t+v_0 sin t]/sin sigma
    +H'(t-sigma)+H'(t+sigma-pi)
    +cot sigma[H(t-sigma)-H(t+sigma-pi)].               (13)
```

Its boundary coefficient is `B=N y`, since
`B=N eta[z_0 cos t+v_0 sin t]+H(t)-H(t-pi)` and
`B_tt+B=-E/2`. It solves both null constraints and the homogeneous
trace equation. The stress is reflection symmetric and nonpositive.
The accepted spatial-constraint argument therefore gives
`Phi>=B-J`, where `J=-E/2=N eta g/2`. During switching (11) yields
`Phi>=-(49/100)N eta`; afterward (12) gives the stronger lower bound
`Phi>-(907/2000)N eta`. Consequently

```
Phi>=-(49/100)N eta throughout the future strip,
Phi_ext>(49/100)N eta
  => R_4^2 >=2G4[Phi_ext-(49/100)N eta]>0.              (14)
```

This is conditional on the actual renormalized extremal contribution
and the spherical reduction. It does not bound omitted four-dimensional
equations. The dynamical dilaton alone still crosses zero near the center;
its sign is not the total-area condition.

After `T_star` the stress is static, and the geometry is its static
part plus the vacuum homogeneous oscillation with period `2pi`.
There is no secular accumulation from later returns. Since
`eta/100<y<13eta/10`, a negligible-backreaction null probe traverses
the strip in finite boundary time obeying

```
10pi/(13eta) < Delta_u <100pi/eta,
Delta_U=(phi_r/N)Delta_u.                              (15)
```

The interaction joins two externally coupled AdS boundaries. This is
a geometric probe result, not a finite-energy capacity bound or a
locally manufactured asymptotically flat shortcut.

## Verification and the next physical step

[Checker](../checks/returned_switch_capture.py),
[receipt](../evidence/returned-switch-capture.json), and
[source/dependency record](../evidence/returned-switch-source-review.json).
Direct high-precision matrix point splitting checks overlapping returns;
the stationary covariance is independently matched to its two frequency
towers. Time-domain work, spatial stress and a Fourier sum with an
analytic tail agree. Eighteen trajectories include redshift and velocity
endpoints; independent nonlinear boundary evolution and retarded bulk
reconstruction check (8)--(15), with a tighter tolerance repeat.

The initial unweighted oscillatory quadrature exceeded its subdivision
limit and was retained as a development failure. Weighted Fourier
quadrature resolves it without changing a physical formula or suppressing
warnings. All-frequency and all-future-time claims use the ordinary
proof; finite numerical samples test implementation and normalization.

The free-junction duration/capture objection is resolved by this
specialization. The complete interacting electric/angular junction
does not inherit (4) or the final vacuum merely because its second-order
restoring modes cancel. Its normal-ordered vertex, state, constrained
massive interactions, geometry and common stress must be matched with
errors below the capture/area margins. Pair preparation, large-N
fluctuations, broader stability and finite-energy transport remain open.
Two-dimensional opening and Fourier control are established methods;
no novelty acceptance or physical creation claim is made.

Remaining to-do list: bound the interacting junction and its continuum
stress against these margins with matched clocks, vertices, initial state
and geometry; include angular/massive interactions and quantum fluctuations;
establish preparation, sufficient stability and finite-energy transport.
