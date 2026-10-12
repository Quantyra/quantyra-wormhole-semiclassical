# Complete free matter response of the prescribed proper-mouth family

2026-10-11, E003/S051-S053. The finite circle-state contacts and the
quadratic proper-time trajectory can now be included in the preceding
[covariant anomaly calculation](covariant-mouth-inertia.md). This gives
the entire statically-subtracted free CFT matter kernel for the specified
regularized product-throat family. Its difference from the leading sharp
canonical kernel is uniformly small at damping resolution, for all
positive exterior return lengths.

This closes a free-matter premise, not the interacting semiclassical
wormhole problem. The metric family is prescribed; its gravity constraint,
corrected background, interacting pressure and angular/compact state
with the physical exterior remain to be controlled. The state/trajectory
terms derived here are not adjustable clock coefficients. No physical
damping, finite-disturbance survival, shortcut or novelty is accepted.

## 1. The state part of the metric response has no extra inertia

Keep the metric and notation of the preceding note,

```text
g=e^(2rho_0(z))[-N^2dt^2+B^2(dz+u dt)^2],
log N=alpha(z) phi+O(phi^2), log B=beta(z) phi+O(phi^2),
u=gamma(z) phidot+O(phi^2),
f=alpha-beta, F=f+gamma', k_n=2pi n/L, nu=2pi/L.
```

For now the unspecified quadratic lapse/spatial-scale terms depend only
on phi; the velocity-dependent trajectory term is added in section 2.
Quadratic terms in u have zero first-order expectation in the static
zero-flux state. Static quadratic terms change K(0), not K(w)-K(0).
Time histories are smoothly switched from the circle vacuum in the past.
Its reference stress is T_0^(tt)=T_0^(zz)=e_0, T_0^(tz)=0, where
e_0=-pi c/(6L^2). The Weyl-invariant state part can be computed in the
flat reference metric; the anomaly and rho_0 terms were already included.

Here is a direct covariant derivation of every state contact. Write
U=delta T^(tt), P=delta T^(tz) for the state contribution. The linear
trace equation and both conservation equations give

```text
delta T^(zz)=U+2e_0 f phi,
U_t+P_z=-e_0(2alpha+2beta-gamma') phidot,
P_t+U_z=-e_0[4alpha' phi+gamma phiddot].                   (1)
```

These follow by varying both the connection and the contravariant
stress, rather than treating its components as metric independent.
At temporal frequency w with Im w>0, the unique retarded solution is

```text
U_n=-4e_0 alpha_n+2e_0 w^2 F_n/(w^2-k_n^2),
P_n=i e_0 w gamma_n+2e_0 k_n w F_n/(w^2-k_n^2).            (2)
```

The response of the linear metric vertex, including the determinant
variation, contributes to W'' the frequency-dependent expression

```text
-2e_0 L sum_n w^2 |F_n|^2/(w^2-k_n^2)
                         -e_0 w^2 integral gamma^2 dz.    (3)
```

Terms independent of w are suppressed in (3). The second metric
variation supplies the term omitted by a stress-commutator-only rule:
g_tt contains gamma^2 phidot^2. Its contraction with T_0 yields

```text
W_state,path^(2)=(e_0/2) integral gamma^2 phidot^2 dt dz.   (4)
```

Thus its W'' is +e_0 w^2 integral gamma^2 and cancels the last term of
(3). The restoring kernel is -W''. The spatial zero mode in the first
term is independent of w and disappears in the static subtraction.
The complete remaining state contribution is therefore

```text
Delta K_state(w)=4e_0 L sum_(n>=1) w^2 |F_n|^2/(w^2-k_n^2)
 =-c L nu^2/(6pi) sum_(n>=1) w^2 |F_n|^2/(w^2-k_n^2).     (5)
```

This is precisely the lower-derivative, -nu^2 part of the known circle
stress spectrum. There is no further state inertia for this metric
family. The preceding anomaly kernel gives its k_n^2 part and the
derived local square. Traceless conservation with retarded initial data
has no remaining homogeneous solution in the upper half-plane, so the
state part is fixed. No instantaneous-vacuum assumption is introduced.

## 2. The induced mouth trajectory is exact to quadratic amplitude

Use throat units x=exp(phi), reference x=1, and

```text
A=acos(epsilon), b=sqrt(1-epsilon^2), T=2A,
r=epsilon/b, kappa=-r/A, L=T+D, D>0, 0<epsilon<=1/4.
```

The exact proper-time condition is cos^2(s)=epsilon^2(x^2-sdot^2).
Expanding in the amplitude of any fixed smooth history gives

```text
s=A-r phi-epsilon phi^2/(2b^3)
                  +d_s phidot^2+O(amplitude^3),
d_s=epsilon^3/(2b^3).                                    (6)
```

This is stronger than using a slow-derivative ansatz: at quadratic
amplitude, sdot^2 contains only the square of the first-order trajectory.
There are no missing higher-time-derivative quadratic terms in (6).
The remainder statement is for a fixed smooth history; it is not a
nonlinear bound uniform over arbitrarily increasing frequencies.

For the throat pullback sigma=(s/A)z, an extra displacement delta s at
fixed x changes the logarithms of the relative lapse and spatial scale by

```text
delta log N=(z rho_0'/A) delta s,
delta log B=(1+z rho_0') delta s/A.                        (7)
```

Use the same positive-h metric regularization as the preceding note:
theta_h=G_h*1_th, d_h=G_h*(delta_-A+delta_A), and
Z_h=G_h*(z 1_th), so Z_h'=theta_h-A d_h. The linear profiles are
gamma_h=kappa Z_h, alpha_h=theta_h+rho_0' gamma_h,
beta_h=kappa theta_h+rho_0' gamma_h. Add to their logarithmic lapse and
spatial scale the quadratic trajectory terms

```text
delta log N=d_s a_h phidot^2, a_h=rho_0' Z_h/A,
delta log B=d_s b_h phidot^2, b_h=a_h+theta_h/A.           (8)
```

This is an explicit common metric prescription. Its static quadratic
terms are those of the proper-cut pullback, smoothed at the same seam;
their choice affects only the unsubtracted stiffness. Quadratic shift
terms, including those obtained by differentiating (6), contract with
the zero background momentum and do not affect this quadratic action.
The family is smooth at every h>0, and agrees with (7) off the seam as
h->0. It does not assert a matched Einstein solution.

## 3. Evaluate the trajectory contact with the full circle stress

At the static background, the covariant coordinate components are

```text
T_tt=e_0+c(2rho_0''-(rho_0')^2)/(24pi),
T_zz=e_0-c(rho_0')^2/(24pi), T_tz=0.                      (9)
```

Their trace and conservation give the standard matter anomaly; on the
whole circle the rho_0'' integral vanishes. For a logarithmic variation
(a_h,b_h), the first variation per unit time is

```text
I_h=integral[-T_tt a_h+T_zz b_h] dz
   =2e_0-c/(24pi A) integral
                  [2rho_0'' rho_0' Z_h+(rho_0')^2 theta_h] dz
   =2e_0-c/(24pi) integral (rho_0')^2 d_h dz.              (10)
```

The last equality uses integration by parts on the complete circle and
Z_h'=theta_h-A d_h. It keeps the seam term. Treating the throat as an
isolated patch would miss it and can change the sign of this contact.

Independently, (10) follows by differentiating the exact static functional

```text
W_static/time=pi c/(6J)
       +c/(24pi) integral v[partial_z log(N_physical)]^2 dz,
J=integral dz/v, v=N/B, N_physical=e^(rho_0)N.             (11)
```

This check fixes the trajectory sign without borrowing the earlier
isolated-patch boundary action. Equation (8) adds d_s I_h phidot^2 to W,
so the restoring correction is -2d_s I_h w^2. For every fixed positive
D and smooth static exterior, d_h converges to the two endpoint deltas.
At the product-throat mouths rho_0'(+-A)=+-b/epsilon. Hence

```text
Delta K_trajectory(w)=w^2 B_epsilon(L),
B_epsilon(L)=c epsilon/(12pi b)
                         +pi c epsilon^3/(3L^2 b^3).      (12)
```

It is positive in the restoring kernel, corresponding to a negative
matter contribution to the kinetic mass. Its magnitude is O(c epsilon),
uniformly in D>0, because L>=2acos(1/4). The sign is not an instability
claim; the leading gravitational mass is much larger in the retained
scaling. The finite-h convergence need not be uniform over arbitrary
exterior derivative profiles. The sharp expression (12) is uniform in
the return length once the prescribed mouth jets are fixed.

## 4. The complete subtracted free matter kernel

Let R_0(w;T,D,r) denote the explicit finite part in the
[moving-mouth note](moving-mouth-clock-vertex.md), equation (12). The
preceding covariant calculation proves that it is the sharp limit of
the known vacuum spectrum plus the anomaly contact. Section 1 now
shows that its state contacts are complete. Section 3 adds the remaining
velocity-dependent proper-trajectory contact. Therefore

```text
Delta K_free,proper(w;epsilon,D)
  =R_0(w;2acos(epsilon),D,epsilon/b)+w^2 B_epsilon(L).      (13)
```

This is the entire free matter K(w)-K(0) for the prescribed metric
family, including both real and imaginary parts. The background is the
conformally mapped circle vacuum, with the same fixed NS/holonomy data
as before, and its linear response is causal. The zero-frequency limit
of (13) is zero. To form K(w) one must use the static Hessian of the same
proper-cut energy; it is not a free parameter. The earlier
[static calculation](proper-clock-anomaly.md) controls that quantity in
its stated static join class. No gravitational kinetic coefficient or
four-dimensional correction is obtained merely by completing the matter
response.

For T and Re w in fixed positive compact bands and 0<eta=Im w<=1, the
previous sharp bound and (12) give

```text
|Delta K_free,proper-Delta chi_(0,0)(w;T,D)|
 <= C c epsilon(1+eta^-1), T=2acos(epsilon).               (14)
```

The comparison can also be made to the asymptotic throat length pi.
The canonical formula is

```text
Delta chi_(0,0)=-c pi TD/(3L^3)
  -c[w-4pi^2/(wL^2)] Q(T,D,w)/(6pi),
Q=sin(wT/2)sin(wD/2)/sin(wL/2).
```

Put p=exp(iwT), q_D=exp(iwD). Then

```text
Q=(i/2)(1-p)(1-q_D)/(1-p q_D),
partial_T Q=(w/2)p(1-q_D)^2/(1-p q_D)^2.                 (15)
```

Since |1-p q_D|>=1-exp(-eta L) and L is bounded below, Q and its T
derivative are bounded by C(1+eta^-1) and C(1+eta^-2), respectively.
The rational prefactors and their derivatives are return-uniform.
Using |T-pi|<=C epsilon proves

```text
|Delta K_free,proper(w;epsilon,D)-Delta chi_(0,0)(w;pi,D)|
 <= C c epsilon(1+eta^-2), uniformly over D>0.             (16)
```

For q independent c=1 channels at eta comparable to delta, the error
divided by q delta is O(epsilon/delta^3). It tends to zero for the retained
epsilon=q^-2 and delta=log^-2(q). This supplies the complete free matter
response premise at that resolution, rather than only an absorption or
isolated spectral estimate. It does not extend the canonical singlet
comparison to the interacting mouth-pressure observable.

## 5. Verification, source comparison and decision

Run `python checks/free_proper_mouth_response.py`; the
[receipt](../evidence/free-proper-mouth-response.json) binds this note,
the checker and all reused calculations. It derives the state Ward
system from metric variations, checks the metric-contact cancellation,
and compares independent finite-mode conservation solves with the
circle-state spectrum. It verifies the amplitude expansion of the
proper-cut condition, differentiates the global static energy against
the trajectory stress integral, refines the mouth regularization, and
checks the complete finite-cut response and return-uniform bounds.
The ordinary argument establishes the all-return claims.

The [source supplement](../evidence/free-proper-mouth-source-review.json)
records the scoped primary review. Conservation, trace and state data
are established stress-tensor methods:
[Christensen--Fulling](https://doi.org/10.1103/PhysRevD.15.2088), section
III, treats the two-dimensional case and its state constants. Their
black-hole boundary conditions and Hawking result are not imported into
this circle problem. The anomaly action, Weyl contact and vacuum stress
spectrum use the already reviewed primary dependencies. The Ward and
trajectory calculations here are explicit; this is not a claim that the
full prior-work proof chain or novelty search has been completed.

**Decision:** continue physical long-return damping as provisional primary
and matched finite-duration transfer as fallback. The free product-family
state and trajectory contacts are now closed at linear response. The
next decisive problem is the renormalized interacting proper-pressure
response, together with constraint-compatible geometry and the coupled
angular/compact state/exterior. The remaining physical gaps should no
longer be described as missing free circle-state inertia. No consequential
novelty or physical endpoint is accepted, and no Lean/Lake or GCP run
occurred. Certification remains attached to a suitable completed endpoint.

Remaining to-do list: control interacting pressure and corrected geometry,
construct the coupled state and physical exterior, establish the physical
endpoint and novelty, independently verify, certify suitable mathematics
on GCP, and prepare the reproducible manuscript.
