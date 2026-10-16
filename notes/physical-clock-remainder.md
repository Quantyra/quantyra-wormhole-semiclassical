# Physical spherical fields and the normalized clock remainder

2026-10-15, E003/S136.3, supporting S074. This derives an exact finite-cut
map from the spherical Einstein equations to the residual in
z_tt+z=g_free/2+r. It identifies a norm of the physical stress, radius,
metric and finite-boundary response that would establish the proposed
|r|<=1e-6 allowance. Candidate scales make one conditional classical
geometric ledger small at a suitable internal seam. They do not provide
the required dynamical stress or boundary-response bounds.

The laboratory source, original preparation rectangle and eta=1/2000
remain those of S135/S136. An internal matching surface is not a new
physical boundary or a change of source. The interval of interest ends
at actual source completion and is at most t=7 only if the separate
robustness calculation establishes that bound. No all-future inference
is made from a bounded driving-interval residual.

## Physical frame, normalization and exact operator identity

Use the physical orbit metric of the spherically symmetric four-metric,
ds4^2=g_ab dx^a dx^b+R^2 dOmega^2. Write X=R^2 and fix magnetic charge.
The classical reduced action and its source conventions are those of
[the nonlinear constraint note](spherical-anomaly-sources.md), (8)-(10):

```text
I0=(1/(4G)) integral sqrt(-g)
       [X R2+(nabla X)^2/(2X)+V(X)],
V(X)=2+6X/L4^2-2 r_m^2/X.
```

For this notation all nonmagnetic sources, including an electric Maxwell
sector if present, are included in Gamma. Define
T_ab=-2 delta Gamma/delta g^ab/sqrt(-g) and
J_X=4G delta Gamma/delta X/sqrt(-g). Thus the metric equation is
E_ab=2G T_ab and the angular equation is
R2-Box X/X+(nabla X)^2/(2X^2)+V'(X)=-J_X.
The magnetic field and cosmological term already in V are not counted
again in T. This is exact for the spherical ansatz; nonspherical
backreaction needs its own projection/error estimate.

Take b to be the extremal reference sphere radius, xi=b^2/L4^2,
C=1+6xi and a=b/sqrt(C), so r_m^2=b^2(1+3xi), V(b^2)=0 and
V'(b^2)=2/a^2. Distinguish the normalized JT area field from X:

```text
Phi=(X-b^2)/(2G),
g=a^2 exp(2omega) g0,
g0=csc(sigma)^2 diag(-1,1), 0<sigma<pi.
```

This is a physical-frame conformal chart on the retained throat patch,
not a removal of the spherical kinetic term by a Weyl redefinition.
Here b is the classical extremal reference of V, and G, L4 and the
finite source/counterterm scheme are fixed throughout the comparison.
A quantum retuning of the extremal radius or a finite boundary
normalization must be supplied as physical matching data; it is not
silently absorbed by fitting b or N to the desired clock.
The chart and its time t must actually be matched to the exterior clocks;
their existence with the bounds below is a hypothesis, not a coordinate
choice that can silently impose a physical result.

The reference normalization is
phi_r=b^3/(G C), u=N U/phi_r, y=dt/du, z=y/eta, B=N eta.
Indeed the static overlap X-b^2 approximately equals 2b^3 rho/L_WH,
while dt/dU=C/L_WH, so the leading area coefficient is
Phi/rho=phi_r*dt/dU=N*y. This fixes the normalization used below.
It does not turn the static dictionary into a dynamical constitutive law;
any renormalized normalization or dynamical matching defect is retained
in the finite-cut response.

Let L0 Phi=g0 Box0 Phi-nabla0 nabla0 Phi-g0 Phi. Direct substitution
in the exact physical metric equation gives

```text
L0_ab Phi=T_ab+G_ab,
G_ab= -G/X [Phi_a Phi_b-g0_ab |dPhi|0^2/2]
      -[omega_a Phi_b+omega_b Phi_a-g0_ab domega.dPhi]
      +g0_ab[(exp(2omega)-1)Phi+a^2 exp(2omega)U(Phi)],
U(Phi)=V(b^2+2G Phi)/(4G)-Phi/a^2
      =-2G r_m^2 Phi^2/[b^4 X].                        (1)
```

All contractions in (1) are with g0. The signs follow from
nabla_g nabla_g Phi=nabla0 nabla0 Phi
-omega_a Phi_b-omega_b Phi_a+g0_ab domega.dPhi.
No small-field expansion is used in (1), except the explicit choice of
the reference linear operator L0. The source T is the sphere-integrated
physical stress: classically T_ab=4pi R^2 T_ab^(4) in the same orbit
coordinates. Transforming a renormalized quantum source to this frame
must also retain the radius source and Ward contacts as in (6)-(7) of
the nonlinear constraint note.

## Exact finite-cut charge, including moving surfaces

Let the internal timelike cuts be sigma=sigma_i(t), i=L,R, with
v_i=sigma_i', |v_i|<1, and s_L=-1,s_R=+1. Put
L_i=a exp(omega_i)/sin(sigma_i), so the cut proper lapse per dt is
L_i sqrt(1-v_i^2). The normal n_i points out of the retained throat;
K_i is its orbit-space scalar extrinsic curvature. Explicitly,

```text
n_i=s_i/[L_i sqrt(1-v_i^2)] (v_i,1),
K_i=s_i/[L_i sqrt(1-v_i^2)]
          [v_i'/(1-v_i^2)+v_i omega_t+omega_sigma-cot(sigma_i)].
```

K_i is the time principal curvature, not the full four-dimensional
timelike-boundary trace. The latter is K_i+2 n_i R/R.
Define the reference-frame charge at the cut by

```text
q_i=s_i[Phi_sigma+cot(sigma_i)Phi]_i
   =L_i sqrt(1-v_i^2)[n_i Phi-(K_i-n_i omega)Phi]
       +s_i[v_i' Phi/(1-v_i^2)-v_i Phi_t]_i.             (2)
```

The physical radius derivative here is n_i Phi=R_i n_i R_i/G.
Equation (2) is an algebraic identity from the normal and curvature,
including acceleration and moving-cut terms. It does not assign the
unknown physical Dirichlet-to-Neumann response. Setting v_i=0 recovers
the static finite-cut expression; discarding the v_i' or n_i omega term
is not valid for a general physical family.

Let T0 be the complete, covariant, anomaly-subtracted free stress of the
unit-width-pi benchmark, evaluated with the actual retimed laboratory
source. It is not the full-P optical-channel energy of S135.1 at P!=pi.
Put E0=integral_0^pi T0_tt, g_free=-E0/B,
I(t)=[sigma_L(t),sigma_R(t)], delta T=T-T0, and

```text
E_miss=integral_([0,pi]\I(t)) T0_tt d sigma,
beta(t)=[q_L+q_R-2N(y+y_tt)-E_miss]/B.                  (3)
```

The elementary identity L0_tt Phi=-partial_sigma(Phi_sigma+cot Phi)
holds without a traceless-source assumption. Integrating (1) over I(t)
at fixed t gives the exact residual formula

```text
z_tt+z=g_free/2+r,
r=-(1/2)[ (1/B) integral_I (delta T_tt+G_tt)d sigma
                                                    +beta(t)]. (4)
```

Beta is a boundary matching defect containing the actual clock jets.
Consequently (4) is not an independently closed predictive equation
until beta is bounded from physical boundary data. Defining r by (4)
does not establish that either r or beta is small.

The moving endpoints do not add terms to this instantaneous spatial
identity. They would add transport terms when differentiating the
finite-cut energy with respect to t. Their instantaneous geometric
effects are already in (2).

E_miss is generally O(epsilon), not O(epsilon^2), even for regular
free stress. It must be retained. For a reference field solving
L0 Phi=T0 with leading clock coefficient N*y,
q_L+q_R=2N(y+y_tt)+E_miss, hence beta=0 at every finite cut.
Subtracting this known free collar energy in (3) specifies the physical
matching target; it does not assume the actual beta vanishes. Surface
stress/actions, exterior stored energy, geometric counterterms and the
asymptotically normalized clock all enter that target together. A fixed
Euclidean source derivative or a static RN mass tangent does not bound it.

### The known AdS anomaly shift must be transferred consistently

The Phi in (1)-(3) is the raw physical area (X-b^2)/(2G). In contrast,
the minimal JT benchmark uses an anomaly-shifted dilaton and a stress
with the same state-independent reference AdS contact subtracted.
They cannot be identified without that known transfer. Let c be the
fixed constant shift specified by that subtraction convention, and set

```text
Phi_J=Phi-c,  L0_ab c=-c g0_ab,
q_c=c[-cot(sigma_L)+cot(sigma_R)],
delta T_J=delta T-L0 c, beta_J=beta-q_c/B.
integral_I (L0 c)_tt=-q_c,
(1/B) integral_I delta T_J,tt+beta_J
                 =(1/B) integral_I delta T_tt+beta.    (4a)
```

Thus (4) is invariant under the consistent reference subtraction. The
physical G in (1), X and its weighted jets still use raw Phi; no area
has been removed from the four-metric. This is comparison bookkeeping;
the physical stress T has not changed. The finite-cut assertion beta=0
for a free JT solution refers to Phi_J solving L0 Phi_J=T0 with the
same subtraction, and consequently to beta_J in this convention.
It need not hold for the raw area charge before the known transfer.

The ledger (5) may use either consistent pair (delta T,beta) or
(delta T_J,beta_J), with the latter appropriate when the inherited
benchmark has absorbed this contact. Their separate absolute-value
bounds differ: a known csc(sigma)^2 contact and its surface charge can
be large and cancel exactly. Such a subtraction must be fixed in
advance by the renormalization convention. It does not permit fitting
c to cancel an unknown state or boundary response. Physical deviations
of the anomaly stress in the evolving geometry, radius-source terms
and other finite counterterms remain in the physical response ledger.

## A sufficient normalized physical ledger

Define on the driving interval

```text
S_T(t)=(1/B) integral_I |delta T_tt| d sigma,
S_G(t)=(1/B) integral_I |G_tt| d sigma.
|r(t)|<=[S_T(t)+S_G(t)+|beta(t)|]/2.                    (5)
```

Thus S_T+S_G+|beta|<=2e-6 is sufficient for the proposed residual
allowance. At eta=1/2000 the combined unnormalized charge/stress error
per real field must be <=2 eta*1e-6=1e-9. This is an expectation-value
stress and canonical surface-response requirement, not a bound on a
determinant value, pointwise gauge variance or optical delay alone.

The geometric term has an explicit weighted-jet estimate. Set

```text
psi=Phi/B, mu=2G B/b^2, vartheta=(1+3xi)/C, 1/2<=vartheta<=1.
```

Suppose both cuts stay inside [epsilon,pi-epsilon], and throughout
the retained strip and driving interval

```text
|sin(sigma) psi|<=M0, |sin(sigma) psi_t|<=Mt,
|sin(sigma)^2 psi_sigma|<=M1,
|sin(sigma) omega|<=mu W0, |sin(sigma) omega_t|<=mu Wt,
|sin(sigma)^2 omega_sigma|<=mu W1.
```

Require mX=1-mu M0/sin(epsilon)>0 and set Omega=mu W0/sin(epsilon).
Then X/b^2=1+mu psi>=mX, |omega|<=Omega. In particular

```text
G_tt/B= -mu[psi_t^2+psi_sigma^2]/[4(1+mu psi)]
         -omega_t psi_t-omega_sigma psi_sigma
         -{[exp(2omega)-1]psi
             -mu vartheta exp(2omega)psi^2/(1+mu psi)}/sin(sigma)^2.
```

Using |exp(2omega)-1|<=2 exp(2Omega)|omega|, obtain

```text
I2=2cot(epsilon), I4=2cot(epsilon)+(2/3)cot(epsilon)^3,
A=Mt^2/(4mX)+Wt Mt,
D=M1^2/(4mX)+W1 M1+2 exp(2Omega)W0 M0
                         +vartheta exp(2Omega)M0^2/mX,
S_G<=mu[A I2+D I4].                                    (6)
```

Proof: the time-derivative products contribute at most mu A csc^2;
the radial, metric-potential and nonlinear-potential terms contribute
at most mu D csc^4. Integrating those two positive weights gives (6).
This is a sufficient separate-term bound. Physical background/charge
cancellations may improve it; it is not a lower bound on the true r.

## Candidate comparison and conditioning of the internal seam

At the retained k=2000 candidate, eta=1/k=1/2000, N=2k^4,
b^2/G=xi k^10 and xi(1+3xi)=pi k^6, with xi>=k^3. Thus

```text
mu=4/(xi k^7)<=4/k^10=3.90625e-33,
lambda0=kappa0*b/C=mu*z_i/2, 1/4<=z_i<=7/25.
```

This comparison uses the frozen eta=1/2000. The relation eta=1/k
is not a proposal to change the laboratory drive in the robustness test.

For the static reference let e=epsilon_cut=lambda0/delta and
epsilon=atan(e), so a cut at |rho|=1/e has precisely that angular
distance from the reference AdS boundary. Then I2=2/e and
I4=2/e+2/(3e^3). The illustrative S134 cut delta=k^-4 gives

```text
e=mu*z_i*k^4/2<=8.75e-21,
mu/e^3=xi^2*k^2/(2*z_i^3)>=5.83e27 at k=2000.
```

Although mu/e is tiny there, the csc^4 ledger in (6) is badly
conditioned unless additional cancellations are used. This does not
make the physical residual large; it shows why relative area smallness
and a tiny static delay cannot certify it by this absolute-value estimate.

An alternative internal comparison seam e=1e-5 gives

```text
delta=lambda0/e<=5.46875e-29,
lambda0/delta=e=1e-5, delta<<1,
mu/e^3<=3.90625e-18.
```

It remains in the reference overlap lambda0<<delta<<1. It changes
neither the source nor the complete propagation channel. The leading
exterior delay replaces the removed throat ends, as in
[RN-AdS finite-cut matching](rn-ads-finite-cut-matching.md); it is not
added to an intact pi strip. The reference thermal mass envelope is
m<=256e^2=2.56e-8 at this seam. That envelope establishes a possible
static comparison, not a bound on the evolving beta or stress.

For illustration only, if all six weighted-jet constants in (6) are
at most one, mX>=0.99 and exp(2Omega)<=1.01, hence A<1.3,D<4.4.
The geometric ledger is then S_G<1.2e-17 and its contribution to |r|
is below 6e-18. These numerical constants are conditional hypotheses,
not established dynamical field bounds. Moving-cut location, velocity,
acceleration and the physical response in (2)-(3) still require control.

More generally, if all six weighted-jet constants are at most K, the
same proof gives A<1.3K^2,D<4.4K^2 whenever mX>=0.99 and
exp(2Omega)<=1.01. At e=1e-5 and k=2000 those last conditions still
hold for K<=1e5. Consequently S_G<1.2e-7 and |r_geom|<6e-8
throughout that conditional jet window. It leaves most of the proposed
allowance for the separate stress and boundary ledgers; it does not
establish that the evolving physical fields stay in the window.

## The angular source is part of the missing physical estimate

The conformal metric jets cannot be assigned independently. In the
same normalized variables, the exact angular field equation is

```text
-2-2 Box0 omega-mu Box0 psi/(1+mu psi)
 +mu^2 |dpsi|0^2/[2(1+mu psi)^2]
 +exp(2omega)[6xi/C+2vartheta/(1+mu psi)^2+a^2 J_X]=0.    (7)
```

Thus a bound on the retarded radius source J_X, compatible stress and
initial/boundary metric data is required to establish the W-jet constants.
The physical-frame Ward identity
nabla^a T_ab=J_X X_b/(4G)+F_ba j^a must hold, including the response
contacts. The existing spherical anomaly functional supplies one explicit
part. It does not determine the conformally invariant state functional,
its massive angular tower, gauge-state response or their boundary data.

Equations (3),(5),(7) specify the useful open datum: a retarded bound on
the physical stress difference and the cut's radius-normal-derivative/
curvature response, in the same laboratory ensemble and clock chart,
through source completion. The finite static delay and small mu do not
give that bound. For example, geometry smallness alone does not exclude
a small symmetric matter excitation supported in the retained patch
with (1/N) integral_I |delta T_tt|>1e-9. That error alone exceeds the
sufficient separate stress allowance in (5), while its area backreaction
can remain tiny. This does not prove the true r is large: boundary or
geometric corrections could compensate it, and energy outside I is not
part of S_T. The specified
physical initial state and its coupled evolution must exclude or bound
it. This is a quantitative state/response target, not another derivative
of an unspecified coupling.

This remainder identity neither solves the state nor proves the proposed
robustness target. Persistent post-completion forcing needs the separate
tail/stationary analysis. No physical opening, stability, transport,
preparation or creation is inferred from this conditional map.

## Predeclared computational checks

The accompanying checker uses 60 decimal digits. Before its first run,
the acceptance tolerances are fixed at 1e-35 for exact numerical tensor,
angular, moving-cut, finite-cut identity and candidate-arithmetic checks.
The symbolic potential and cut-charge identities must simplify to zero.
The weighted estimate is checked with the stated analytic global jet
constants and positive slack; sampled values supplement that elementary
proof. Four-dimensional Ricci/Einstein tensors are built independently
from first and second metric derivatives, with the magnetic stress and
cosmological term subtracted, rather than imported from the reduced
field equation. Moving-cut curvature is also computed from covariant
acceleration. Synthetic fields and a smooth chiral reference stress test
the identities; they do not realize the frozen laboratory experiment or
prove its physical response bound. No physical residual certificate is
inferred from those algebra checks.

The full tensor/cut/quadrature run passed before the explicit contact
transfer subsection was added. The final receipt retains that complete
run and its original hashes. A targeted subsequent run verifies (4a)
symbolically and numerically, records current hashes, and confirms by
AST comparison that the original expensive checks and their constants
are unchanged. The original full-run checker source is preserved in
the append-only log. This provenance does not represent a second full
tensor/quadrature run.

Artifacts: [checker](../checks/physical_clock_remainder.py),
[receipt](../evidence/physical-clock-remainder.json),
[execution log](../evidence/physical-clock-remainder.log).

Remaining to-do list: establish the matched physical chart, stress/radius
source and surface-response bounds in (3),(5)-(7); certify the driving
residual against the common-drive robustness allowance; control the
post-drive tail and complete the remaining physical S074-S079 gates.
