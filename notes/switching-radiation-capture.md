# Smooth switching with correlated fermion radiation

S072, 2026-10-12. **The symmetric free-Majorana JT boundary model has a
finite smooth-switch capture range even after including its switching
stress.** This extends [the reduced capture proof](entangled-pair-capture.md).
It assumes the prepared pair and an existing bridge. Four-dimensional
formation, bulk area positivity and stability in other sectors are open.
This is an ordinary model calculation, not a Lean certificate or a novelty
acceptance.

## Model and conventions

Take global AdS2 with unit curvature radius, strip coordinate
`0<sigma<pi`, and `N` independent **real** massless Majorana fields. Use
the cutoff-boundary limit and an instantaneous, energy-independent
orthogonal junction between the two boundaries. Its prescribed angle is
`theta`, with positive sign chosen to lower the ground-state energy.
Reflection phases are absorbed in the two-channel basis. The initial
state is the uncoupled strip vacuum, with the symmetric boundary trajectory
of a prepared TFD-like pair at its time-reflection slice. This assumes the
appropriate initial entanglement; it is not a preparation protocol.

The two endpoints share global time `t` and rescaled boundary time `u`.
Write `y=dt/du=exp(phi)>0`. Physical boundary time `U` satisfies
`u=N U/phi_r`, where `phi_r` is the single-boundary Schwarzian coefficient.
Our angle variable is **defined** by `eta=theta/(4 pi)`. It agrees with
the leading small-angle MQ coupling; it is not an unproved exact equality
between a four-dimensional interaction coefficient and that coupling.
Subtract the state-independent AdS Weyl-anomaly stress, absorbing its
constant dilaton shift. All energies below use that same subtraction.

The microscopic benchmark is [Maldacena and Qi, sections 2.2, 4 and
appendices B/C](https://arxiv.org/html/1804.00491v3). Their sudden-switch
overlap divergence at dimension 1/2 motivates computing the stress for a
smooth switch. The calculation below uses an explicitly defined rotation,
not the malformed second boundary equation in the HTML rendering of C.206.

## Correlated point splitting before the first return

Unfold the uncoupled strip to an antiperiodic chiral circle of length
`2 pi`. In the chosen basis the two incident channels are separated by
`pi`. Their covariance, with `d=t_1-t_2`, is

```
J = [[0,-1],[1,0]],       R(theta) = exp(J theta),
G_in(d) = I/[4 pi i sin((d-i0)/2)] + i J/[4 pi cos((d-i0)/2)].
```

The off-diagonal term is the prepared state's correlation. It cannot be
discarded by treating the endpoints as independent reservoirs. Until an
outgoing signal first returns at elapsed global time `pi`, incident fields
retain this initial covariance. The outgoing covariance is
`R(theta(t_1)) G_in(d) R(theta(t_2))^T`. Its trace difference is exactly

```
Tr(G_out-G_in) = [cos(Delta theta)-1]/[2 pi i sin(d/2)]
                - i sin(Delta theta)/[2 pi cos(d/2)],
Delta theta = theta(t+d/2)-theta(t-d/2).
```

For a real field, the sum of the two outgoing-minus-incoming energy fluxes
is `-(i/2) lim_(d->0) partial_d Tr(G_out-G_in)`. Expansion gives

```
F_t/N = [theta_dot^2-theta_dot]/(4 pi),
E_t/N = -theta/(4 pi) + R_rad/N,
R_rad/N = (1/(4 pi)) integral theta_dot^2 dt >= 0.             (1)
```

Here `F_t` and `E_t` include all `N` fields; each endpoint supplies half
the total flux. The linear negative term comes from the cross correlation.
The positive quadratic term agrees, with real versus complex field
counting accounted for, with the familiar pump-dissipation coefficient
in [Avron et al., section 5.3](https://arxiv.org/html/math-ph/0305049v1).
That paper's reservoir formula alone would not supply the linear term.
Our point-split derivation is exact for this nondispersive junction before
the first return; it does not use an adiabatic reservoir approximation.

`R_rad` is the integrated positive driving contribution. It is **not**
the exact excitation energy above the final ground state. MQ's exact
static spectrum, with `epsilon=theta/pi=4 eta`, gives

```
E_vac/N = -epsilon(1-epsilon)/4 = -eta+4 eta^2,
E_exc/N = R_rad/N - 4 eta_f^2.                               (2)
```

This avoids counting the quadratic ground-state shift twice. If the
global switch duration is `T<pi`, Cauchy--Schwarz gives
`R_rad/N >= theta_f^2/(4 pi T) = 4 pi eta_f^2/T`, so the excitation
energy in (2) is nonnegative, as required.

## Charges and the common boundary clock

The three global AdS2 Killing fields are

```
K_0 = partial_t,
K_1 = cos(t)cos(sigma) partial_t - sin(t)sin(sigma) partial_sigma,
K_2 = sin(t)cos(sigma) partial_t + cos(t)sin(sigma) partial_sigma.
```

A pair of equal pulses emitted at time `t_e` sits at
`sigma_L=t-t_e`, `sigma_R=pi-(t-t_e)`. The left pulse contributes
coefficients `cos(t_e), sin(t_e)` to its two nonenergy charges, and the
right pulse contributes their negatives. This holds for signed pulse
energy too. The initial subtracted charges vanish. The chosen symmetric
junction preserves equal endpoint stress, including later returns, so
both nonenergy charges remain zero in expectation. This is the symmetry
restriction that permits a single common boundary trajectory; it does
not dispose of arbitrary asymmetric matter perturbations or fluctuations.

The total pure-gravity time-translation charge is
`Q_0,grav=-2N(y+phi_uu/y)`. MQ's convention for matter is
`q_0,matter=-E_t`. Inserting the *full* matter energy (1) into the
constraint `Q_0,grav+q_0,matter=0` yields

```
phi_uu = -y^2 + (eta-R_rad/N)y/2.                            (3)
```

The interaction-induced negative energy is already in (1); adding a
second independent `eta y/2` term would double count it. No independently
variable proper-mouth clock is introduced in the closed flat-space
problem examined by the earlier research cycle.

For relative stress, the common reparametrization has
`F_u=y^2 F_t`; the state-independent Schwarzian anomaly cancels between
the two states being compared. Since `theta_u=4 pi eta_u`, (1) becomes

```
R_rad/N = 4 pi integral (eta_u^2/y) du.                       (4)
```

This inverse redshift is essential. A slow ramp in boundary time need
not be slow in global time. Define

```
s=eta_f u, z=y/eta_f, w=phi_u/eta_f,
a=eta/eta_f, r=R_rad/(N eta_f).
```

Equations (3)--(4), including the radiation response during the ramp,
give the closed system and work identity

```
z_s = z w,
w_s = -z^2 + (a-r)z/2,
r_s = 4 pi eta_f a_s^2/z,
t_s = z,
e = w^2+z^2-(a-r)z,
e_s = -a_s z + 4 pi eta_f a_s^2.                             (5)
```

Thus radiation is coupled to the trajectory, not evaluated on an
unchanged vacuum path and added afterward.

## Capture theorem with radiation

Choose the cubic ramp `a(s)=3(s/delta)^2-2(s/delta)^3` on
`0<=s<=delta`, zero before and one afterward. This extended ramp is C1
and has finite driving cost. Suppose

```
0<eta_f<=1/2000,  1/4<=delta<=1/2,
1/4<=z(0)<=3/4,  w(0)=r(0)=0.
```

Then the switching solution exists throughout the interval, remains
before the first signal return, and obeys

```
r_f < 1/4,  T <= 48/127 < pi,  e_f <= -3/32.                 (6)
```

After switching, the common trajectory has the all-time lower bound

```
y >= eta_f [1/2-sqrt(5/32)] > 0.1047 eta_f.                  (7)
```

Proof. Bootstrap `r<=1/4` up to its first possible failure or the first
return. Since `0<=a<=1` and `r>=0`, (5) implies
`w_s<=-z^2+az/2<=1/16`. Consequently

```
z<=z_0 exp(s^2/32)<= (3/4)exp(1/128)<=96/127=:Z.
```

The last inequality uses `exp(x)<=1/(1-x)` for `0<=x<1`.
On the bootstrap interval,
`w_s>=-Z^2-Z/8= -10740/16129 > -2/3`. Integration gives

```
w>=-2s/3,
z>=z_0 exp(-s^2/3)>=(11/12)z_0>=11/48.                     (8)
```

The ramp cost is `I=integral a_s^2 ds=6/(5 delta)<=24/5`.
Integrating the radiation equation using (8),

```
r_f <= (192 pi/11) eta_f I
    <= (4608 pi/55) eta_f < 1/4.
```

Also `T=integral z ds<=delta Z<=48/127<pi`. Both stopping
conditions therefore remain strictly separated from their thresholds;
the positive lower bound on `z` and bounded right-hand side continue
the solution across the whole ramp. This closes the bootstrap and
justifies using the initial incident covariance throughout switching.

Finally integrate the exact work identity:

```
e_f <= z_0^2-(11/12)z_0+4 pi eta_f I
    <= -1/8+1/32 = -3/32.
```

The quadratic's maximum on `[1/4,3/4]` is `-1/8`; the radiation cost
is below `1/32` at the stated bound on `eta_f`. All strict numeric
comparisons follow rationally from `pi<22/7`.

After the ramp the junction is constant and unitary, hence total global
matter energy is conserved, even when pulses return. Its symmetric
nonenergy charges stay zero. Put `A=1-r_f`, with `3/4<A<=1`.
The autonomous equation has conserved
`e=w^2+z^2-Az=e_f<0`; its turning points are
`z_+-=(A+-sqrt(A^2+4e_f))/2`. The lower turning point decreases with
both `A` and `e_f`, giving (7). The solution stays in a compact
positive `z` interval for all subsequent boundary time, and global
boundary time covers an unbounded interval. This proves the symmetric
boundary capture statement. A subsequent reduced phase-space energy
increase smaller than `(3/32)N eta_f^2`, at fixed matter charge and
coupling, retains negative energy. General injected matter and unequal
endpoint disturbances are not covered. A captured trajectory generally
oscillates; cooling to a ground state has not been proved.

## Physical normalization and the transfer gate

For fixed-charge RN-AdS4 of extremal radius `b`, curvature radius `L4`,
and Newton constant `G4`, write

```
C=1+6 b^2/L4^2,
Q_geom^2=b^2(1+3 b^2/L4^2),
M(R)=[R+R^3/L4^2+Q_geom^2/R]/(2G4),
T_H(R)=[1+3R^2/L4^2-Q_geom^2/R^2]/(4 pi R).
```

Expansion at fixed charge gives
`M-M_ext=2 pi^2 phi_r T_H^2+O(T_H^3)`, where
`phi_r=b^3/(G4 C)`. With `c=N_f q` complex magnetic LLL channels,
`N=2c` real fields. Thus the candidate physical duration is
`Delta_U=(phi_r/N)(delta/eta_f)`. This fixes a thermodynamic clock
normalization; it does not establish the dynamical reduction error.

As a factor check only, the flat limit of the transparent static
Majorana vacuum has `theta=pi/2`, `-E_vac/N=1/16`, and stationary
`y=1/32`. It predicts
`ell=32 phi_r/N=16 b^3/(G4 c)`, the leading
[MMP throat coefficient](https://arxiv.org/html/1807.04726v3).
The transparent angle is outside the theorem's small-angle interval.
This check does not identify the coupled AdS boundary with a physical
exterior or reopen the closed independent-mouth-clock proposal.

The next hard test is reconstructing the bulk dilaton for the driven
stress and determining whether its area interpretation stays admissible.
After that, match the four-dimensional quantum state, boundary interaction,
gauge sector and exterior delay with errors smaller than the capture
margin. A semiclassical use also requires adequate large `N`, negligible
relative quantum fluctuations and a valid near-extremal expansion. No
finite-parameter bound for these omitted effects has been established.
Initial-pair preparation, nonspherical stability and physical signal
transport remain separate uncompleted requirements.

## Reproducible checks

[Checker](../checks/switching_radiation_capture.py),
[receipt](../evidence/switching-radiation-capture.json), and
[source scopes/hashes](../evidence/switching-radiation-source-review.json).
The checker independently multiplies the full two-channel covariance at
65-digit precision in 12 point-split cases, verifies convergence to the
stress, checks the charge cancellation and RN-AdS thermodynamic coefficient,
and symbolically verifies the work identity and rational bootstrap bounds.
It also integrates 27 coupled trajectories and one tighter refinement.
The largest final normalized energy is `-0.1538313554`; largest radiation
charge is `0.120686315`; largest elapsed global time is `0.367603395`.
The work residual is below `5.42e-16`, checked post-switch drift below
`2.89e-15`, and refinement changes the selected final energy by `3.34e-16`.
Finite integrations support the ordinary all-time proof in this model;
they do not establish full spacetime stability.

Remaining to-do list: reconstruct the driven bulk dilaton and area; bound
four-dimensional state, coupling and gauge/exterior matching errors; account
for initial-pair preparation; test other perturbations and physical transport.
