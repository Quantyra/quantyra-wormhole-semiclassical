# Driven bulk geometry and the extremal area contribution

S072, 2026-10-12. **The smooth-switch solution has an explicit retarded
bulk dilaton, an all-future-time lower bound and a conditional positive
total-area interpretation.** The dynamical dilaton itself crosses zero;
its sign must not be confused with that of the physical area. This resolves
the bulk area gate within the symmetric free-fermion JT model, conditional
on its extremal area contribution. It does not bound four-dimensional
dynamical corrections or establish formation from ordinary matter.

The state, junction, anomaly subtraction, real-field count `N`, clock
`y=dt/du`, and cubic ramp are those of the
[switching-stress calculation](switching-radiation-capture.md). Keep its
entire sufficient range:

```
0<eta_f<=1/2000, 1/4<=delta=eta_f Delta_u<=1/2,
1/4<=z_0=y_0/eta_f<=3/4, initial reduced velocity zero.
```

We use `Phi` for the dynamical bulk dilaton and `Phi_ext` for the constant
extremal contribution, avoiding confusion with the logarithmic boundary
variable `phi=log(y)`. All statements concern expectation values in the
specified large-N semiclassical model and the boundary-cutoff limit.

## The source and its returned pulses

The prior bootstrap gives `z>=11/48` during switching. Since the cubic
has `a_s<=6`, its global-time angular speed obeys

```
0<=theta_dot=4 pi eta_f a_s/z<=1152 pi eta_f/11<1.
```

Consequently the **net** stress emitted at either endpoint is nonpositive,
even after including its positive radiation term. Define the single-pass
pulse, extended by zero outside the switch interval `[0,T]`, by

```
p(v)=N[theta_dot(v)^2-theta_dot(v)]/(8 pi)<=0,
T<pi,
P=integral_0^T p(v)dv = -N eta_f A/2,
A=1-r_f, 3/4<A<=1.
```

This stress sign is a consequence of the proved parameter interval; it
would not hold for arbitrary faster or stronger switches.

In the unfolded two-channel basis, propagation across the strip is
`a(t)=J b(t-pi)` and the junction is `b(t)=R(theta(t))a(t)`, with
`J=[[0,-1],[1,0]]` and `R=exp(J theta)`. These rotations commute.
Iterating the covariance therefore replaces its local angle by
`Theta(t)=sum_(j>=0) theta(t-j pi)`, with only finitely many nonzero
terms at each finite time. Constant propagation phases commute with the
initial covariance and add no stress. Since `T<pi`, at most one term in
`Theta_dot` is nonzero. The same point-split stress formula then gives
the exact outgoing pulse train

```
h(v)=sum_(j>=0) p(v-j pi),
T_--(t,sigma)=h(t-sigma),
T_++(t,sigma)=h(t+sigma-pi),
T_+-=0,  T_tt=T_+++T_--<=0.                                (1)
```

Coordinates are `x^+=t+sigma`, `x^-=t-sigma`, `0<sigma<pi`.
The stress is reflection symmetric. After switching its global energy
is conserved; returning pulses redistribute that energy, without adding
new driving work. The state is not assumed to relax to a static vacuum.

## A retarded reconstruction with no center source

After the state-independent anomaly shift, the JT equations are

```
-(partial_+^2+cot(sigma) partial_+)Phi=T_++,
-(partial_-^2-cot(sigma) partial_-)Phi=T_--,
(-partial_t^2+partial_sigma^2-2 csc(sigma)^2)Phi=0.           (2)
```

These follow from the action and null constraint in
[Maldacena and Qi, section 2.2](https://arxiv.org/html/1804.00491v3).
Let `H(x)=0` for `x<=0` and define

```
H(x)=-integral_0^x [1-cos(x-v)]h(v)dv,  x>0.
```

Then `H'''+H'=-h`, and the required initial-value solution is

```
Phi(t,sigma) = N y_0 cos(t)/sin(sigma)
             + H'(t-sigma)+H'(t+sigma-pi)
             + cot(sigma)[H(t-sigma)-H(t+sigma-pi)].         (3)
```

It equals the initial TFD dilaton and has its initial zero time derivative
on `t=0`. Each source term begins only when the corresponding boundary
signal arrives. To verify (2), write its particular solution as
`f'(x^+)-cot(sigma)f(x^+)+g'(x^-)+cot(sigma)g(x^-)`, with
`f(x)=H(x-pi)` and `g(x)=H(x)`. The two null constraints reduce separately
to `-(f'''+f')` and `-(g'''+g')`; the trace equation vanishes. The formula
is reflection symmetric and smooth across the center. No artificial
localized center stress is introduced. The pulse endpoints are continuous
for the C1 ramp; their higher derivative jumps do not produce a delta
function in (2).

The leading coefficient at either boundary is

```
B(t)=lim_(sigma->0) sin(sigma)Phi(t,sigma)
    =N y_0 cos(t)+H(t)-H(t-pi)
    =N y_0 cos(t)-integral_0^t [1-cos(t-v)]p(v)dv.           (4)
```

The last integral uses the zero-extended single-pass pulse. In particular
`B'''+B'=-p` and `B''+B=-E_t/2`. Since the boundary condition requires
`B=N y`, these are exactly the charge equation and work response of the
earlier coupled switching calculation, now checked from the full bulk.
For `t>=T`, write `C=integral p(v)cos(v)dv` and
`S=integral p(v)sin(v)dv`. Equation (4) becomes

```
B(t)=-P+(N y_0+C)cos(t)+S sin(t)>0.                         (5)
```

Its positivity is the preceding capture theorem. Also
`H(x+2pi)=H(x)-2P` and `H'(x+2pi)=H'(x)` for `x>=0`:
the extra two complete pulses have opposite cosine contributions.
The constant shifts cancel between the two terms of (3). Hence the
complete bulk dilaton is exactly `2pi` periodic for `t>=pi`. This proves
absence of secular growth from arbitrarily many returns in this model,
not merely in a finite numerical time window.

## The dynamical dilaton crosses zero

At the center, no switched stress has arrived before `t=pi/2`. Equation
(3) implies

```
Phi(pi/2,pi/2)=0,
partial_t Phi(pi/2,pi/2)=-N y_0<0.                         (6)
```

Here `H'(0)=H''(0)=0`. Thus the dynamical dilaton becomes negative just
after this time. Demanding `Phi>0` by itself would reject this prepared
TFD protocol. The condition relevant to a near-extremal area is instead
`Phi_ext+Phi>0`. This distinction is explicit in the preceding bulk
transition literature: [Bak, Kim and Yi, sections 2 and
4](https://arxiv.org/html/1901.07679v2) retain the constant extremal term
when identifying an area-zero singularity. Their scalar transition is a
known comparator; it does not prove our free-fermion switching bound or
four-dimensional dynamical matching.

## All-time area bound

Put `W(t,sigma)=sin(sigma)Phi(t,sigma)` and
`k(t,sigma)=-T_tt(t,sigma)>=0`. The spatial JT constraint is

```
partial_sigma[(partial_sigma W)/sin(sigma)]=k.
```

Reflection symmetry sets `partial_sigma W=0` at `sigma=pi/2`.
For `0<sigma<=pi/2`, integration from the center and the boundary gives

```
W=B-integral_0^sigma sin(v) integral_v^(pi/2) k(t,r)dr dv.
```

Let `J(t)=integral_0^(pi/2)k(t,r)dr=-E_t/2`. Since `B>0` and `J>=0`,

```
B-J(1-cos(sigma)) <= W <= B,
B/sin(sigma)-J tan(sigma/2) <= Phi <= B/sin(sigma),
Phi >= B-J.                                                (7)
```

The same inequalities hold by reflection on the other half. During the
ramp, `B>=N eta_f(11/48)` and `J<=N eta_f/2`, so
`Phi>=-(13/48)N eta_f`. After the ramp, `J=N eta_f A/2` and the
conserved reduced energy `e_f<=-3/32` implies

```
B/(N eta_f)=A/2+beta cos(t-t_*),
beta^2=A^2/4+e_f<=5/32.
```

Thus, uniformly over the entire future strip,

```
K=sqrt(5/32),
Phi(t,sigma)>=-K N eta_f,    t>=0, 0<sigma<pi.              (8)
```

In the spherical reduction with this action normalization, the leading
area variable is `R_4^2/(2G4)=Phi_ext+Phi`. The renormalized extremal
constant is retained after the anomaly subtraction; its leading classical
value is `b^2/(2G4)` for extremal radius `b`. A sufficient condition is

```
Phi_ext > K N eta_f
  => R_4^2 >= 2G4[Phi_ext-K N eta_f] > 0.                   (9)
```

This is conditional on the actual extremal coefficient and the validity
of the reduction. It is not an exact solution of the four-dimensional
Einstein equations. Near-extremal control needs substantially more than
bare area positivity. For an inner strip with `sin(sigma)>=rho_c>0`,
the prior capture bounds and (7) also give

```
-K N eta_f <= Phi <= (1/2+K)N eta_f/rho_c.
```

Therefore the necessary small-area-deviation regime has an explicit
candidate hierarchy `N eta_f/(Phi_ext rho_c)<<1`. It has not yet been
converted into a bound on the omitted field equations, gauge response,
exterior transfer or finite-cutoff errors. A positive radius alone is
not that missing estimate.

## Geometric probe channel

Under (9), no area-zero surface truncates the future strip in this
leading description. Both boundary trajectories have unbounded global
time. A zero-backreaction radial null probe sent at `t_i>=0` reaches
the other boundary at `t_i+pi` in the cutoff limit. With
`m=1/2-K`, `M=1/2+K`, the common clock has
`eta_f m<=y<=eta_f M` throughout the future, including the ramp.
Consequently

```
pi/(eta_f M) <= Delta_u=integral_(t_i)^(t_i+pi) dt/y(t)
             <= pi/(eta_f m),
Delta_U=(phi_r/N)Delta_u.                                  (10)
```

This identifies the delayed bulk causal path in the imposed AdS setup.
It does not bound the backreaction or capacity of a finite-energy signal,
provide an asymptotically flat shortcut, or turn the external boundary
interaction into a locally constructed physical device.

## Checks and next physical decision

[Checker](../checks/driven_bulk_geometry.py),
[receipt](../evidence/driven-bulk-geometry.json), and
[source review](../evidence/driven-bulk-source-review.json).
Symbolic checks verify both null equations, the trace equation, spatial
constraint and sufficient sign/area inequalities. Ten parameter cases
compare the retarded reconstruction to independently integrated spatial
constraints and to the original nonlinear boundary equations during and
after switching. They test reflected pulses, late periodicity, the center
zero crossing, the all-time bound on a spacetime grid and probe travel
times. An interpolation-density refinement is retained. The exact proof,
not the grid, supplies the all-time conclusion. Quadrature warnings are
treated as failures; nearly coincident floating-point wavefront cuts are
integrated directly over their tiny interval.

The literature already contains black-hole-to-traversable-wormhole
transitions in two-dimensional gravity. This result is a checked
free-fermion smooth-switch specialization with an explicit sufficient
area margin; no broad novelty claim is accepted. The hardest next task
is the quantitative four-dimensional transfer, particularly the quantum
state and gauge/exterior response under the prescribed coupling. Large-N
fluctuations, initial-pair preparation, additional perturbation sectors
and finite-energy transport remain unproved. No Lean certificate,
publication, outreach or paid computation is selected.

Remaining to-do list: bound four-dimensional state/coupling, gauge/exterior
and finite-cutoff errors against the margins; establish initial-pair
preparation; test quantum fluctuations, other perturbations and signal
backreaction.
