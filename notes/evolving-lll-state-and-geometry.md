# Evolving spherical LLL state, ray transport and bulk geometry

2026-10-15, E003/S137.2, supporting S074/S076. The finite free
lowest-Landau-level state source is explicit on an evolving spherical
conformal geometry. Its bulk radius source vanishes. Its dependence on
the physical boundary history is an exact causal characteristic map,
including the source's laboratory scattering times, return holonomy and
initial vacuum calibration. The known two-dimensional Weyl completion
has a useful signed charge bound. None of these results supplies the
physical exterior/clock matching or the complete four-dimensional source.

The fixed laboratory waveform, eta=1/2000, initial preparation rectangle
and proposed driving allowance of S135-S137 are unchanged. Internal
matching cuts are transmitting seams, not physical reflecting walls.
The initial full-channel state is the prepared correlated vacuum, not
two independent incoming reservoirs. No new capture certificate follows.

## Exact finite state source in the physical spherical frame

Write

```text
ds4^2=e^(2rho)(-dt^2+d sigma^2)+X(t,sigma) dOmega^2,
R=sqrt(X)>0, x^+=t+sigma, x^-=t-sigma,
g_+-=-e^(2rho)/2, central charge c2=N/2.
```

Fix magnetic flux and retain only the free angular zero modes, with no
electric field and a charge-symmetric prepared state. The exact canonical
rescaling is Psi=R^-1 exp(-rho/2) psi times the normalized sphere zero
mode. The radial spin connection cancels the R derivative and the
two-dimensional spin connection cancels the rho derivative. The zero
angular eigenvalue leaves the flat massless characteristic equations.
This is the reduction of [MMP appendix C](https://arxiv.org/html/1807.04726v3#A3),
also recorded in [the source reduction](dirac-spherical-source-reduction.md).
No expansion in radius, conformal factor or frequency is used.

Consequently the finite state-dependent bulk part of the canonical LLL
functional has no explicit X argument. In a comparison with the same
local subtraction and the same fixed flux-normalized domain data,

```text
Delta J_X=0, Delta j=0,
Delta T_+-=0, Delta T_++=Delta t_+(x^+),
Delta T_--=Delta t_-(x^-),
nabla_g^a Delta T_ab=0, trace_g Delta T=0.              (1)
```

This extends the [fixed-throat state-difference result](lll-state-difference-source.md)
to an arbitrary evolving conformal spherical geometry. A smooth source
change has the same local bulk singularity, so its finite bulk stress
can also be compared by interior point splitting. Source-dependent
boundary counterterms are not thereby canceled. Varying a physical
radius condition may move the domain or change its trace normalization;
those variations give boundary/state response. Equation (1) does not
discard them or set the absolute four-dimensional radius source to zero.
Moving an internal seam produces no extra physical LLL source: its two
transmitting descriptions must be matched.

On a genuinely fixed complete conformal strip of width pi, the actual
source theta(t)=theta_b(u(t)) therefore gives exactly S135's returned
bulk state stress, even when rho and X evolve. The conformal factor and
radius do not create an additional unknown finite LLL bulk state term.
The actual clock t(u), complete domain and local anomaly must still be
matched physically.

## Physical moving endpoints and the calibrated return

For an explicit scalar construction take the mirror-symmetric paired
domain sigma_L=s(t), sigma_R=pi-s(t), with |s_t|<1 and positive width.
These denote the true endpoints of the complete optical problem,
possibly regulated before the asymptotic limit. They are not the S136
internal cuts. Fix the laboratory clock u and a common increasing t(u).
Define the outgoing and incoming null coordinates at a scattering event,

```text
O(u)=t(u)-s(t(u)), I(u)=t(u)+s(t(u)),
Rcal(u)=O^(-1)(I(u)-pi),
v=O(u), b(v)=I(O^(-1)(v))-pi.                          (2)
```

In b's expression O and I are functions of the same laboratory u.
Equivalently O_t(t)=t-s(t), I_t(t)=t+s(t) give
b=I_t o O_t^(-1)-pi. Timelikeness makes these maps increasing and
positive width makes b(v)<v. The exact canonical return in laboratory
half densities contains sqrt(Rcal'). In the inherited real basis it is

```text
psi_out(u)=sqrt(Rcal'(u)) R(theta_b(u)) J psi_out(Rcal(u)),
J=[[0,-1],[1,0]], R(alpha)=exp(J alpha), J^2=-1.         (3)
```

The source is evaluated at the actual physical scattering event u.
Transport through exterior segments must be included before assigning
that event; arrival at an internal cut is a different time. The fixed J
contains the free propagation/reflection holonomy. It cannot be omitted.
The two-end laboratory normalization is the canonical flux normalization
of [the exterior elimination](exterior-lll-delayed-boundary.md).

Assume a stationary past, sufficiently regular timelike data and a positive lower
flight time on the interval. The past vacuum fixes an increasing optical
phase H, and the actual source fixes A, by the causal equations

```text
H(Rcal(u))=H(u)-pi,
A(u)=theta_b(u)+A(Rcal(u)), A=0 before the drive.
h(v)=H(O^(-1)(v)), a(v)=A(O^(-1)(v)),
h(b(v))=h(v)-pi,
a(v)=theta_b(O^(-1)(v))+a(b(v)).                       (4)
```

The method of steps uniquely extends the compatible past data over any
finite interval with finitely many returns. In a past static channel of
physical optical width P_init, h(v)=(pi/P_init)v+constant. Its slope is
set by that actual initial state and time normalization; it is not fitted
to one. For unequal/asymmetric ports, the corresponding channel-dependent
ordered half-density transport is required. A single scalar h is not
asserted for arbitrary unequal physical ports.

The original correlated vacuum kernel is

```text
G0(d)=I/[4pi i sin((d-i0)/2)]+iJ/[4pi cos((d-i0)/2)],
G0(d+pi)=J G0(d).
G_out(v1,v2)=sqrt(h'(v1)h'(v2))
             R(a(v1)) G0(h(v1)-h(v2)) R(a(v2))^T.     (5)
```

Equation (5) follows by iterating (3): the half densities compose, the
rotations commute with J, and the stated quasi-periodicity handles
different numbers of returns on the two arguments. It agrees with the
prepared covariance in the past. Characteristic propagation and unitary
paired scattering preserve positivity and CAR. This derives the state,
rather than postulating independent stress functions.

## Point splitting and covariance of the complete stress

The diagonal singularity of (5) is 1/(2pi i d). Subtract that line
singularity and apply the real-fermion prescription -iN/2 times the
separation derivative. Expanding the half-density pullback and rotations
gives the canonical outgoing stress

```text
t(v)/N=-h'(v)^2/(96pi)-Sch(h,v)/(48pi)
                     +[a'(v)^2-h'(v)a'(v)]/(8pi).      (6)
```

One way to see the last term is to use h as coordinate: its phase rate
is da/dh=a'/h', so h'^2 p(da/dh) gives the displayed polynomial.
The singular kernel supplies the Schwarzian and the NS vacuum term.
An affine h'=pi/P reproduces the complete finite-P coefficient in
[S135](retimed-returned-quantum-stress.md). Fixed endpoints at 0,pi give
h(v)=v, a'=sum_j theta'(v-j*pi), even for a nontrivial actual t(u).
There is then no extra bulk Schwarzian of the laboratory clock.

The covariant stress in the standalone covariant two-dimensional LLL
representative is

```text
T_++=t_+(x^+)-N[(rho_+)^2-rho_++]/(24pi),
T_--=t_-(x^-)-N[(rho_-)^2-rho_--]/(24pi),
T_+-=-N rho_+-/(24pi), J_X,LLL=0.                     (7)
```

Here t_-=t(x^-) and t_+=t(x^+-pi) in the symmetric domain. It obeys
conservation and trace T=N R2/(48pi). This Weyl relation agrees with
[MMP appendix F](https://arxiv.org/html/1807.04726v3#A6).
The complete four-dimensional regulator allocates its local terms
between the full anomaly functional, the remaining determinant and
counterterms. Equations (7) are one explicit LLL representative, not an
additional anomaly to append to the already retained four-dimensional
W+Gamma_hat. Only the finite state difference (1) is independent of
that allocation. The full radius tadpole and other angular modes remain.

The [quantum-reference centering](quantum-reference-clock-centering.md)
retains both stationary stress and radius tadpoles. It is compatible with
this state difference. Equation (16) below explicitly uses bare V and
the un-subtracted representative (7). To use the centered V_eff instead,
the complementary sources must be those of
Gamma-Gamma_ref-Gamma_LLL,2D. Those complementary sources need not
vanish at the product while the absolute LLL contact remains separately
displayed. Subtracting that known LLL reference contact requires its
own consistent transfer. No absolute LLL contact is appended to a total
reference from which it has already been removed. The new reference
note is a scope cross-reference, not a live checker input.

In particular under x^+=f_+(xi^+), x^-=f_-(xi^-),

```text
rho_tilde=rho o f+(1/2)log(f_+' f_-'),
t_tilde_+=f_+'^2 t_+ o f_+-N Sch(f_+,xi^+)/(48pi),
Q_tilde_+=f_+'^2 Q_+ o f_++N Sch(f_+,xi^+)/(48pi),
Q_+=-N[(rho_+)^2-rho_++]/(24pi).                       (8)
```

The analogous minus identities and tensorial mixed component prove
that T transforms covariantly. The two Schwarzians cancel exactly.
Equivalently Sch(h o f)=f'^2 Sch(h) o f+Sch(f) in (6), while the
Weyl term supplies the opposite cocycle. Thus h''' alone is not a
coordinate-invariant physical error or evidence of radiation. Bounds
below apply in the actual matched clock/conformal chart, or to the
complete covariant difference. Merely changing coordinates cannot create
an extra laboratory Schwarzian or a new physical clock force.

Subtracting an undriven state on the same evolving domain cancels the
same h-vacuum, Schwarzian and local Weyl terms, leaving the phase
polynomial in (6). Comparing to a different physical domain/reference
does not cancel them. These are different comparisons.

## A quantitative target for the physical characteristic map

Compare actual (h,a) with a fixed reference (h0,a0). Formula (6) gives

```text
Delta t/N=[a'^2-h'a'-a0'^2+h0'a0']/(8pi)
         -[h'^2-h0'^2]/(96pi)
         -[Sch(h)-Sch(h0)]/(48pi).                     (9)
```

If |a'|<=M, |a0'|<=M0, |h'|<=H1, |h0'|<=H0 and their derivative/
Schwarzian differences are bounded by eps_a,eps_h,eps_S, respectively,
then each chiral error per real field is at most

```text
C_ray=[(M+M0+H0)eps_a+M eps_h]/(8pi)
                    +(H1+H0)eps_h/(96pi)+eps_S/(48pi).
S_state=(1/B) integral_I |Delta T_state,tt|
                       <=2 length(I) C_ray/eta.        (10)
```

Its clock-remainder contribution is at most length(I)*C_ray/eta.
For example, if only eps_S contributes and length(I)<=pi, eps_S<=2.4e-8
would give |r_state|<=1e-6 at eta=1/2000. This is a target in a matched
chart, not an established physical error or a coordinate-invariant
Schwarzian-alone criterion.

The target can be connected constructively to the true boundary rays.
If u_j(v) are the successive physical laboratory scattering times,
then a'=sum_j theta_b'(u_j) u_j'. Comparing to the benchmark u_j0,

```text
eps_a <=sum_j [||theta_b''||inf sup|u_j'| sup|u_j-u_j0|
                  +||theta_b'||inf sup|u_j'-u_j0'|].   (11)
```

The source derivatives are those of the fixed laboratory waveform,
not derivatives of a newly chosen global-time drive. Causal iteration
also gives

```text
h'=h'(b)b',
Sch h=(b')^2 Sch h o b+Sch b,
b'=(1+s_t)/(1-s_t),
Sch b={2s_ttt/(1-s_t^2)+6s_t(s_tt)^2/(1-s_t^2)^2}
                                    /(1-s_t)^2,       (12)
```

where the endpoint derivatives are evaluated at t=O_t^(-1)(v).
For |s_t|<=nu<1, b'<=L=(1+nu)/(1-nu). If at most m returns occur,
an initially affine calibration gives
|Sch h|<=||Sch b||inf sum_(j=0)^(m-1) L^(2j).
A lower flight time d and a finite history of length H allow
m<=1+ceil(H/d). Equations (11)-(12) specify first and third jet
data of the actual full-channel transport. A static delay, small
boundary displacement or monotonic clock alone cannot supply them.

## Known Weyl charge with only weighted first conformal jets

In a fixed matched unit-strip chart write
rho=log(a_ref)-log(sin sigma)+omega. For the NS reference, the local
Weyl completion cancels the flat null Casimir and leaves the known
contact -c_A g0, c_A=N/(48pi). This is the constant anomaly shift of
[the remainder identity](physical-clock-remainder.md), not a change of
physical area or a retuning of the stationary sphere. The difference
of the known covariant two-dimensional completion from this reference is

```text
A_tt/N=-[omega_t^2+omega_sigma^2
              -2cot(sigma)omega_sigma-2omega_sigma_sigma]/(48pi).
integral_I A_tt/N=[omega_sigma+cot(sigma)omega]_L^R/(24pi)
 -integral_I[omega_t^2+omega_sigma^2-2csc(sigma)^2 omega]/(48pi). (13)
```

This instantaneous signed spatial identity is valid for moving internal
cuts too. Differentiating the charge in time would require transport
terms. The endpoint term in (13) must not be dropped or presumed to
cancel the unknown physical boundary defect.

Let cuts stay in [epsilon,pi-epsilon] and suppose

```text
|sin(sigma)omega|<=mu W0, |sin(sigma)omega_t|<=mu Wt,
|sin(sigma)^2 omega_sigma|<=mu W1,
I2=2cot(epsilon), I4=2cot(epsilon)+(2/3)cot(epsilon)^3,
I3=csc(epsilon)cot(epsilon)-log(tan(epsilon/2)).
```

Using the endpoint and bulk terms separately gives

```text
S_A,charge=|(1/B) integral_I A_tt|
 <=2mu(W0+W1)csc(epsilon)^2/(24pi eta)
   +[mu^2(Wt^2 I2+W1^2 I4)+2mu W0 I3]/(48pi eta).      (14)
```

This is not S_T=(1/B) integral|A_tt|. No absolute spatial stress or
area bound follows from it. In the exact S136 clock remainder one may
instead bound the state, rest and geometric terms in absolute norm,
and retain this signed known charge:

```text
|r|<=[S_state+S_A,charge+S_rest+S_G+|beta|]/2.          (15)
```

At eta=1/2000, mu<=3.90625e-33, tan(epsilon)=1e-5 and W0,Wt,W1<=1e5,
(14) is below 5.2e-16, contributing less than 2.6e-16 to |r|.
These are explicit conditional geometry hypotheses, not measured or
proved candidate bounds. They do not impose the S136 raw-area jet
window on a quantum-shifted stationary radius. The corrected stationary
reference, time normalization and internal overlap must be established
consistently; the present state formula does not fit any of them.

## Explicit leading LLL bulk equations and the remaining sources

Take the physical spherical action and source conventions of S136.
In the displayed two-dimensional LLL representative put kappa=GN/(12pi)
and V(X)=2+6X/L4^2-2r_m^2/X. Let T_rest,J_rest denote every physical
source not allocated to (7), including matched local/tadpole terms,
massive angular modes, interactions and finite counterterms. The exact
mixed and angular equations give

```text
X_+-=-exp(2rho)V/4-kappa rho_+-+2G T_rest,+-,
(2-kappa/X)rho_+-=
 exp(2rho)(V/X-V')/4+X_+X_-/(2X^2)
                 -2G T_rest,+-/X-exp(2rho)J_rest/4.    (16)
```

The two null constraints are

```text
-X_(++) +2rho_(+)X_(+) +X_(+)^2/(2X)
 =2G[t_+(x^+)+T_rest,(++)]-kappa[rho_(+)^2-rho_(++)],
```

and the analogous minus equation. Parenthesized plus signs denote null
indices/derivatives; t_+ is the chiral function, not a derivative. The mixed equation
also reads (X+kappa rho)_+-=-exp(2rho)V/4+2G T_rest,+-.
These expressions follow from E_+-=X_+-+exp(2rho)V/4 and
E_(++)=-X_(++)+2rho_(+)X_(+)+X_(+)^2/(2X), with the spherical angular equation
retained. The joint Ward identities of the total matched source propagate
the compatible constraints.

For T_rest=J_rest=0 this closes the leading free-LLL bulk model once
the actual characteristic state functions and compatible initial/
boundary data are given. Its mixed principal system is nondegenerate
for X>kappa/2. If X>=0.99 b_ref^2, the candidate mu bound implies
kappa/X<=mu/(24pi eta*0.99)<1.05e-31. This condition is not a solution
estimate, and nondegeneracy does not provide existence through t=7 or
determine the physical boundary clock. Extra physical higher-derivative
sources cannot be declared harmless by this two-dimensional calculation.

The hardest remaining physical datum is the complete retarded
characteristic/clock map, calibrated to the specified stationary state
and asymptotic laboratory normalization, together with the joint
radius/metric boundary response and remaining sources in (16). Equations
(5)-(12) remove an unspecified free LLL bulk state functional from that
task. They do not determine beta, the full physical initial geometry,
the massive-mode load or controller/surface stress. No complete physical
backreaction, opening, capture, stability, transport or preparation is
claimed.

## Acceptance fixed before checker execution

Use 60 decimal digits. Exact symbolic Ward, trace, coordinate-cocycle,
radius-source, null/mixed/angular-system and signed charge identities
must simplify to zero. Numerical causal return, holonomy, coordinate
covariance, independent signed-charge quadratures and candidate bound
checks use absolute tolerance 1e-35. Independent covariance point
splitting at separations 1e-4,1e-5,1e-6 must agree with (6) within
2e-10 at the finest separation. Wrong omission of the compensating
coordinate Schwarzian must produce discrepancy >1e-7.

The characteristic examples use the actual compact shaped Theta0 with
eta=1/2000 and tau=1/4. For verification their actual clock is chosen
equal to the frozen reference clock, so theta_b(u(t))=Theta0(t) exactly;
this avoids replacing the laboratory waveform. A small compact C3
optical reparametrization constructs a compatible moving endpoint via
(4); its retarded source phase is computed by actual scattering times,
not chosen arbitrarily. These synthetic boundary geometries are not
solutions of the physical four-dimensional initial/boundary problem.
No Lean is invoked.

Artifacts: [checker](../checks/evolving_lll_state_and_geometry.py),
[receipt](../evidence/evolving-lll-state-and-geometry.json),
[log](../evidence/evolving-lll-state-and-geometry.log).

Remaining to-do list: establish the calibrated physical characteristic
map and its source/metric response; solve compatible coupled geometry
and boundary data with the complete matched sources; bound the S136
residual, post-drive tail and physical spatial area.
