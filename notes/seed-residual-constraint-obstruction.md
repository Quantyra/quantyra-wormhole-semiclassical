# The seed residual obstructs full homogeneous constraint recovery

2026-10-21, E003/S168. Ordinary argument before focused verification.

The high real-frequency resonances of the reduced charged seed equations
do not lift to homogeneous solutions of the full logarithmic metric-chart
Hessian in the retained auxiliary domain. In fact that full homogeneous
system has zero kernel there, at every sufficiently high real frequency.
This is an off-shell obstruction. It neither removes a physical branch
on an exact solution nor establishes stability or a wormhole no-go theorem.

Keep the original action, contacts, magnetic flux, HH response and
temperature. The [continuous feedback theorem](continuous-causal-feedback-and-resonance.md)
and its reduced resonances remain valid. Here we impose the additional
tensor rows that were explicitly outside that theorem's scope. The seed
is the rational alpha=0 comparison used by
[physical reconstruction](physical-curvature-reconstruction.md).

## 1. Which homogeneous response is being tested

Write g=F(y)(-dt^2+dy^2)+r(y)^2 dOmega^2, r_y=-F, and let

    calE^a_b=diag(E_t,E_y,E_o,E_o),
    Z=E_t+E_y+2E_o=-T, O=E_t+E_y, W=F r^2.

Z here denotes the FULL equation trace, not R+12lambda from some earlier
notes. The background residual is smooth in the open exterior and
conserved in its orbit indices. The magnetic field has only angular
components, so the orbit component of the gauge Ward force F_ba J^a
vanishes for b=t,y, including its variation at fixed magnetic field.
No unproved shape-current cancellation is needed for this necessary test.

Extend the optical chart by independent metric coordinates A,B,p,C:

    g_tt=-F exp(2A), g_yy=F exp(2B),
    g_theta_theta=r^2 exp(2p), g_phi_phi=sin(theta)^2 g_theta_theta,
    g_ty=F C.

At the seed all four coordinates vanish. The directions considered here
have A=B=p+v, C=0. Their volume variation is w0=4p+2v. Define

    S^a_b=delta calE^a_b+w0 calE^a_b.                      (1)

The variations of the diagonal action gradients, divided by their common
nonzero volume factor, are S^t_t, S^y_y and 2S^theta_theta. This includes
the exponential-chart second-derivative contacts. The shift gradient
variation is a nonzero multiple of S^y_t; its background value is zero.
Thus the full homogeneous Hessian rows in this specified chart require
S=0. Off shell, this condition depends on the nonlinear chart. It is not
silently identified with delta calE=0 or with a physical perturbation
equation around an exact solution.

The optical p,v rows are respectively

    S^t_t+S^y_y+2S^theta_theta=delta Z+w0 Z,
    S^t_t+S^y_y=delta O+w0 O.                             (2)

These are precisely the two fixed-density rows retained in
[the shape congruence](charged-curvature-free-shape-form.md). The invertible
curvature reconstruction and shape congruence preserve their homogeneous
zeros. Imposing (2) alone leaves S_t=-e, S_y=e, S_o=0 and S^y_t=j,
S^t_y=-j. Neither e nor j was set to zero by the reduced construction.

## 2. The residual terms in the Ward identity

Varying the mixed-index conservation equation gives

    nabla_a delta calE^a_b
      +delta Gamma^a_ac calE^c_b-delta Gamma^c_ab calE^a_c=0.

Using nabla calE=0 and delta Gamma^a_ab=partial_b w0, (1) yields

    nabla_a S^a_b=sum_a E_a delta Gamma^a_ab
                  =Z partial_b p+O partial_b v, b=t,y.   (3)

For the last equality, the t and y diagonal connection variations are
both partial_b(p+v), and the two angular ones are partial_b p. These
are connection components with the first index fixed, not their sum.
Equation (3) keeps the background residual even when it is small.

For the two-row zero just described, put Ecal=W e and Jcal=W j. Direct
four-dimensional divergence gives

    partial_t Ecal-partial_y Jcal=-W(Z p_t+O v_t),
    partial_t Jcal-partial_y Ecal=-W(Z p_y+O v_y).         (4)

A full homogeneous Hessian zero therefore obeys, for exp(gamma t) with
real gamma>0,

    Zp+Ov=0, Z p_y+O v_y=0.

Eliminating O WITHOUT differentiating the unknown quantum residual gives

    Z(p v_y-v p_y)=0.                                   (5)

No value, derivative bound, analyticity or asymptotic coefficient of O
is needed below. This avoids replacing the actual orbit stress by a
local approximation. Independent conserved state sources cannot change
the geometric residual term in (3); the kernel statement here retains
the original fixed-state homogeneous auxiliary problem.

## 3. The actual seed has a nonzero trace residual near infinity

The complete trace is local by the unchanged anomaly identity:

    T=R+12lambda0+epsilon[7Wc^2/36+22j_R^2/3
                         -11R^2/72-10q_m^2/r^4+Box R],
    Wc=-F_rr+2F_r/r+2(1-F)/r^2,
    j_R=-F_rr/4+(F-1)/(2r^2),
    R=-F_rr-4F_r/r-2(F-1)/r^2,
    Box R=F R_rr+(F_r+2F/r)R_r.

Here epsilon=G/(60pi)>0 and q_m is the fixed nonzero magnetic flux
parameter, not a metric amplitude. With u=B/r, x=lambda B^2,
mu=m/B, e0=epsilon/B^2 and Qhat=1+3x+d, the seed has

    F=x/u^2+1-2mu u+Qhat u^2-2d u^3+d u^4.

The outer vacuum relation lambda0=lambda+11epsilon lambda^2/6 cancels
the constant trace. R+12lambda starts at r^-5, Wc^2 at r^-6 and
j_R^2 at r^-8. Consequently, exactly for a polynomial P,

    B^2 Z=u^4[c+u P(u)], c=10e0 q_m^2>0.                (6)

This is not merely an assumed asymptotic sign. Let |P|_1 be its finite
coefficient norm and u_h=B/r_h. On

    0<u<u_* = min(1/2,u_h/2,c/[2(1+|P|_1)]),

the bracket in (6) is at least c/2. Thus Z is nonzero on an actual
outer collar. The collar may be very small; no useful numerical size
or physical frequency window is claimed. The nonzero coefficient is
the retained magnetic anomaly. Dropping it would change the theory.

## 4. The differentiated metric trace survives the full shape response

Fix gamma above all the previously retained real-frequency thresholds.
For w in the common auxiliary shape energy space, let v=w/r and
p_Q=Q_gamma w/r. The global coercivity bounds control

    J2[v]=||H_s v||^2, J2[p_Q+v/3]=||H_s(p_Q+v/3)||^2,
    H_s=-a^-1 partial_Y(a partial_Y)+s a^-2,
    a=sqrt(CF)/r, dY=a dy, s=gamma^2.

Near the outer seam, a is positive with bounded smooth coefficients,
and Y and y are smoothly equivalent coordinates. The local scalar
graph norm controls H2: integrate by parts for the first derivative,
then solve the displayed equation for the second derivative after a
fixed collar cutoff. The compact shape core has v=0 near the seam and
p_Q=O(y^4). Passing to its form closure therefore gives

    v,p_Q in H2(0,b), v(0)=v_y(0)=p_Q(0)=p_Q,y(0)=0.

This follows from the two retained J2 bounds, not a new boundary axiom.
For any such f, the fundamental theorem and Cauchy--Schwarz give

    |f_y(y)|<=sqrt(y)||f_yy||_(0,y)=o(sqrt(y)),
    |f(y)|<=y^(3/2)||f_yy||_(0,y)/sqrt(3)=o(y^(3/2)).     (7)

For U in the maximal curvature graph domain, the exact identity in
[the metric-trace theorem](coupled-fast-metric-boundary-traces.md) is

    R_gamma U=U/(6d0)+z_U, z_U in D(K_gamma), d0=E+4lambda.

It gives z_U=o(y^(3/2)), z_U,y=o(sqrt(y)). Combining it with (7),
the differentiated maximal trace of U, and r^-1=lambda y+O(y^3), gives
for the FULL p=(R_gamma U+Q_gamma w)/r

    p=c0 y^(3/2) A(y)+o(y^(3/2)),
    p_y=c0 sqrt(y)[3A(y)/2+k B(y)]+o(sqrt(y)),
    v=o(y^(3/2)), v_y=o(sqrt(y)), c0=lambda/(6d0)>0,
    A=A_c cos(k log(y/y0))+A_s sin(k log(y/y0)),
    B=-A_c sin(k log(y/y0))+A_s cos(k log(y/y0)).          (8)

The p_Q term must not be omitted when proving the derivative statement.
Both p and v have C1 representatives on the collar.

## 5. Simple zeros turn the Ward condition into a contradiction

Consider real profiles; reality of the response reduces a complex kernel
to its real and imaginary parts. If (A_c,A_s) is nonzero, (8), in
logarithmic coordinate tau=-log(y/y0), says that p/y^(3/2) and its
tau derivative approach a nonzero sinusoid and its derivative. Choose
small disjoint fixed phase neighborhoods of the sinusoid's zeros.
Outside them the function stays away from zero. Inside each one its
derivative has one sign and a uniform nonzero lower magnitude for all
sufficiently large tau; the endpoint values have opposite signs.
There is exactly one simple zero in each such neighborhood and no
other sufficiently small-y zeros.

On the collar where Z!=0, (5) gives (v/p)_y=0 between consecutive
zeros. Thus v=c_n p on each interval. At a shared simple zero y_n,
p_y(y_n)!=0. C1 continuity implies

    c_n p_y(y_n)=v_y(y_n)=c_(n+1) p_y(y_n),

so ALL constants agree. Evaluating along successive sinusoidal peaks,
the v=o(y^(3/2)) estimate forces that common constant to be zero.
Then Zp+Ov=0 requires p=0 throughout the collar, contradicting its
nonzero fast trace. This proves that no nonzero fast-trace auxiliary
pair can satisfy all homogeneous logarithmic-chart Hessian rows.

Finally, if the fast pair is zero, the earlier
[zero-fast injectivity theorem](high-frequency-zero-fast-injectivity.md)
forces U=0 in the homogeneous reduced problem at sufficiently high
gamma. The positive auxiliary shape inverse then forces w=0 and p=0.
Combining the two cases proves the claimed zero full kernel in this
specified auxiliary sector. It covers every fixed reflecting phase and
the maximal trace domain; no new phase has been selected.

## 6. Consequence, source scope and verification

The reduced resonances are now known to fail a necessary full constraint
test on the off-shell seed. They must not be described as full seed
metric modes awaiting only a boundary-phase choice. No full tensor
inverse or range theorem follows from injectivity.

On an exact self-consistent background calE=0, both Z and O vanish,
and the right side of (3) vanishes. The proportionality argument then
has no premise. A small residual is enough for this exact off-shell
obstruction, so its zero limit cannot be used to infer on-shell mode
exclusion or stability. Chart Hessians and the actual equation
linearization also coincide only after the background equations hold.
This result does not justify order reduction, a source projection,
counterterm changes or deleting the high-frequency branch physically.

[Anderson, Molina-Paris and Mottola, section I around (2)-(3)](https://arxiv.org/html/0907.0823v2)
formulate causal response about a self-consistent background, with local
contacts and independent state variations distinguished. That established
framework motivates this distinction; it does not prove our charged
seed obstruction. Equations (1)-(8) and the simple-zero argument are
derived here using the retained charged estimates.

The focused checker freezes this note, its source review and checker,
then reconciles all predecessor hashes. It computes the connection and
its variation directly from a four-dimensional metric, verifies the
off-shell divergence and both density rows, and checks the exact seed
trace polynomial. Independent finite differences of full metric
Christoffels test the residual contacts at two precisions. Oscillatory
fixtures test the differentiated trace and show that satisfying the
time constraint alone with a logarithmic amplitude ratio fails the
radial constraint; a piecewise-ratio example detects the C1 gluing
requirement. These are mechanism checks. The continuum domain and
zero-kernel theorem remain ordinary proofs, not formal certification
or an evaluation of the charged stress kernel.

The next physical task requires an actual self-consistent background or
a controlled transfer of the FULL constrained problem to one. S167's
signed compatibility and nonlinear solve remain necessary. S168.3-.4
and S077 stay open; no physical fast-branch exclusion is claimed.

Remaining to-do list: self-consistent nonlinear geometry and signed
compatibility; full constrained physical response with boundary/action,
horizon and initial/state data; calibrated validity range and all-sector
stability; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
