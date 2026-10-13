# A finite-ramp capture range for a prepared entangled pair

S072, 2026-10-12. **The constrained reduced AdS2 dynamics has a finite
range of preparation temperatures and switching times that stays
traversable indefinitely within that model.** The proof has a quantitative
energy margin. Microscopic switching radiation and four-dimensional
transfer are not included, so this is not yet a physical creation method.

## Model and known comparator

Start from the coupled-boundary Schwarzian action and vanishing total
SL(2) charges in
[Maldacena and Qi, section 4](https://arxiv.org/html/1804.00491v3).
Use their dimensionless time u_tilde, call it u below, and their positive
coupling eta. After imposing the constraints and choosing equal left and
right trajectories, phi=log(dt/du) obeys, at Delta=1/2,

    phi'' = -exp(2phi) + eta(u) exp(phi)/2.

For constant eta the reduced physical energy is

    E/N = phi'^2 + exp(2phi) - eta exp(phi).

This equation is obtained from the zero-charge constraint, not merely by
choosing an arbitrary stationary trajectory. N is the action coefficient
in the source's rescaled convention. A four-dimensional identification
must derive it and the time rescaling; it is not automatically q.

The source derives the static minimum and a switch-off comparison with a
thermofield double (TFD), and then shows corrections to the state overlap.
The capture mechanism below specializes that established reduced system;
no novelty verdict is claimed. This AdS boundary problem is distinct from
the closed asymptotically flat proper-mouth clock construction.

We assume a prepared TFD-like pair at its time-symmetric slice, before
turning the coupling on. Its reduced initial data are

    y0 = exp(phi(0)) = 2pi/beta_tilde > 0,     phi'(0)=0.

This assumes the initial entanglement and non-traversable bridge. It does
not create either from ordinary disconnected matter.

## Constant coupling and finite-amplitude criterion

Set y=exp(phi). After the switch is complete, eta=eta_f>0 and

    E/N = phi'^2 + (y-eta_f/2)^2 - eta_f²/4.

If E<0, y is confined to the positive interval

    y_- <= y <= y_+,
    y_+/- = [eta_f +/- sqrt(eta_f²+4 E/N)]/2.

The constant solution is y=eta_f/2. The condition E<0 is stronger than
linear stability: it keeps finite reduced disturbances away from y=0.
Since dt/du=y>=y_->0, global AdS2 time continues without an upper endpoint.
The two boundary trajectories retain communication in the ideal strip;
a global-time null crossing pi takes boundary time at most pi/y_-.
This statement does not include exterior propagation or a matter signal's
backreaction.

For an ideal sudden switch, initial (y0,v0=phi'(0)) is captured exactly
when v0²+(y0-eta_f/2)²<(eta_f/2)². We do not use a sudden microscopic
switch: at Delta=1/2 its particle production is UV divergent, as the
source's section 4.4 explains.

## Finite-ramp theorem

Let eta rise monotonically from zero to eta_f over duration Delta_u.
It may be any absolutely continuous ramp. Normalize

    s=eta_f u, z=y/eta_f, w=phi'/eta_f,
    a=eta/eta_f, delta=eta_f Delta_u.

The constrained equation and reduced energy become

    dz/ds=z w,    dw/ds=-z²+a z/2,
    e=w²+z²-a z,  de/ds=-a' z.

Assume 1/4<=z0<=3/4, w0=0, 0<=a<=1, a'>=0 almost everywhere,
and 0<delta<=1/2. Then at the end of the ramp,

    e_final <= -15/112 < 0.                         (capture margin)

Proof. Completing the square gives w'<=a²/16<=1/16, so
w(s)<=s/16 and z(s)<=z0 exp(s²/32). Conversely w'>=-z² implies

    w(s)>=-z0² s exp(delta²/16),
    z(s)>=z0 exp[-z0² delta² exp(delta²/16)/2].

The exponent's magnitude is at most (9/128)exp(1/64), bounded above by
(9/128)(64/63)=1/14 using exp(x)<=1/(1-x) for 0<=x<1.
Thus z(s)>=13 z0/14, using exp(-x)>=1-x. Integrating the work identity,

    e_final = z0² - integral_0^delta a'(s) z(s) ds
            <= z0² - 13 z0/14.

This convex quadratic reaches its maximum over [1/4,3/4] at an endpoint.
The endpoint values are -19/112 and -15/112, proving the claim. The same
bounds control the solution through the finite ramp, so no assumed
existence beyond a blowup is hidden in the argument.

After switching, energy is conserved and e>=-1/4 automatically. The
uniform turning interval therefore gives

    y >= eta_f [1/2-sqrt(13/112)] > 0.

Any subsequent *reduced phase-space perturbation* whose increase of E is
less than (15/112)N eta_f² still has negative energy. This is not a bound
for arbitrary injected matter: matter carries SL(2) charges and changes
the constraint and potential.

## Verification and the hard physical step

[Checker](../checks/entangled_pair_capture.py),
[receipt](../evidence/entangled-pair-capture.json), and
[source scopes](../evidence/formation-boundary-source-review.json).
Symbolic verification checks the work identity and rational margin.
Independent ODE integration uses 27 combinations of initial redshift,
duration and linear/cubic/cosine ramp. It checks work balance during the
switch, then energy conservation and the analytic turning interval over
150 units of normalized time. These finite integrations supplement the
ordinary argument; they do not establish arbitrary-time physical stability.

The theorem tests an actual finite-time protocol in the specified reduced
model. Its main missing ingredient is the quantum state created by changing
the boundary condition. The source's overlap correction is nonlocal in
time and diverges for an abrupt Delta=1/2 switch. The reduced action omits
that particle production. Replacing a step by a smooth ramp regularizes
UV frequencies but does not alone bound the produced energy below the
capture margin.

The next calculation is therefore the renormalized stress/energy emitted
by a smooth time-dependent massless-fermion boundary junction, with the
physical boundary clock and the coupling normalization fixed. Compare its
work and matter SL(2) charges to the margin above; do not merely compare
an overlap probability. A prospective scale interval
`UV_time << Delta_u <= 1/(2 eta_f)` requires eta_f small, but its existence
and useful energy bound must be derived in the mapped four-dimensional
theory. The classicality condition and throat matching must hold too.

There is no dissipation in the reduced theorem: a generic captured state
oscillates rather than relaxing to the ground state. Full stability would
also require nonsymmetric matter, gauge and nonspherical sectors. The
global research goal remains active and incomplete.

Remaining to-do list: calculate switching radiation and its charges on
the correct boundary clocks, establish controlled four-dimensional matching
and initial-pair preparation, and test the other physical disturbances.
