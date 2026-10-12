# The covariant anomaly cancels the sharp pressure inertia

2026-10-11, E003/S051-S053. The divergence found in the
[moving-mouth spectral response](moving-mouth-clock-vertex.md) cancels
against a local term fixed by the covariant anomaly action. It is not
necessary to prescribe a smoothing width to make this particular part
of the free response finite. The cancellation holds in an explicit smooth
metric regularization, includes the curved background's Weyl functional,
and gives a return-uniform sharp limit.

This improves the previous interpretation, not its spectral calculation.
The latter correctly warned that a static subtraction alone was
insufficient. Here the required local term is derived from the same
covariant action, rather than chosen as a counterterm. The result covers
the universal anomaly Hessian and the known free vacuum commutator. It
does not yet include all state-dependent metric contact terms, the
second-order induced-cut trajectory, the interacting pressure response,
or four-dimensional matching. A complete physical damping claim remains
unaccepted.

## 1. Fix the action and distinguish its local and spectral parts

Use the conformal prescription of the earlier notes, signature (-,+),
and matter central charge c, with no gravity/ghost contribution. In the
topologically trivial perturbative sector the Lorentzian anomaly action
and its Weyl difference are

```text
W_an[g]=-c/(96pi) integral sqrt(-g) R Box^-1 R,
W_an[e^(2rho) ghat]-W_an[ghat]
 =c/(24pi) integral sqrt(-ghat) [rho Rhat+(grad rho)^2].      (1)
```

The retarded Hessian is obtained with the retarded inverse, analytic for
Im w>0. This does not mean that a one-branch action with a nonsymmetric
retarded kernel is varied naively: (1) specifies the quadratic form first,
then its causal continuation. Circle-state and global terms are additional
data, not determined by the anomaly alone. They are separated below.

On a flat reference circle of length L take smooth real periodic profiles
alpha,beta,gamma and a small common history phi(t). Prescribe the metric

```text
ghat=-N^2 dt^2+B^2(dz+u dt)^2,
N=exp(alpha phi), B=exp(beta phi), u=gamma phidot.
C=beta-gamma',                 F=alpha-C.                 (2)
```

This is a concrete smooth metric family near phi=0. Additional quadratic
metric terms, such as the induced-cut correction proportional to phidot^2,
have to be evaluated separately against the background one-point stress.
They are not inferred from the linear vertex.

Direct linearization of the curvature gives

```text
Rhat^(1)=2[C phiddot-alpha'' phi].                         (3)
```

Let k_n=2pi n/L and use Fourier coefficients with normalization 1/L.
The quadratic restoring kernel is minus the second variation of W.
For the anomaly action its contribution at each spatial mode is

```text
K_an,n(w)=cL/(12pi)
 (k_n^2 alpha_n-w^2 C_n)(k_n^2 alpha_-n-w^2 C_-n)
 /(w^2-k_n^2).                                            (4)
```

The spatial sum in (4) is over all integers. The n=0 expression is its
continuous algebraic value cL w^2|C_0|^2/(12pi). Subtracting w=0 and
using F=alpha-C yields the exact identity

```text
Delta K_an(w)= cL/(6pi) sum_(n>=1)
                 w^2 k_n^2 |F_n|^2/(w^2-k_n^2)
              +c w^2/(12pi) integral C^2 dz.               (5)
```

For example, at zero frequency (4) gives
-c integral(alpha')^2/(12pi), agreeing with the covariant static vacuum
energy's second variation. For a homogeneous conformal perturbation
alpha=beta, gamma=0, (5) gives cL beta^2 w^2/(12pi), as follows directly
from the negative Lorentzian time-gradient term in (1). These fix the
normalization and sign without matching a desired cancellation.

The first term of (5) is the central, k^3 part of the stress spectrum for
the effective transition profile F. The local second term cannot be
omitted when converting that spectrum into a metric response. For a
circle vacuum the full commutator replaces k_n^2 by k_n^2-nu^2,
nu=2pi/L, in the spectral term. That lower-derivative state contribution
is finite for the mouth profiles below; the anomaly does not fix the
remaining state-dependent metric contacts.

## 2. The fixed curved background is essential

Now let g=e^(2rho_0(z)) ghat, with a fixed smooth periodic rho_0. The
background need not be flat. The Weyl difference in (1) is local, so its
kinetic contribution can be computed without solving another response
problem. Define the connection components

```text
omega_z=[B_t-(Bu)_z]/N,
omega_t=N_z/B+u omega_z,
sqrt(-ghat) Rhat=2(partial_t omega_z-partial_z omega_t).
```

Integration by parts on the circle and in a compactly supported history
gives the exact integrand

```text
W_cocycle=c/(24pi) integral dt dz
 [2rho_0' {N_z/B+u[B_t-(Bu)_z]/N}
                    +(N/B-Bu^2/N)(rho_0')^2].             (6)
```

Its quadratic phidot^2 coefficient is
c[2rho_0' gamma C-(rho_0')^2 gamma^2]/(24pi). Combining its restoring
kernel with the local term in (5) completes a square:

```text
Delta K_an,curved = central spectral term for F
              +c w^2/(12pi) integral [C-rho_0' gamma]^2 dz. (7)
```

Equation (7) is exact for the anomaly Hessian of (2) multiplied by the
fixed background Weyl factor. The other pieces of (6) are frequency
independent and disappear in the static subtraction. No new local
coefficient has been chosen.

## 3. A smooth metric regularization of the proper-mouth variation

Retain the preceding throat notation

```text
A=acos(epsilon), T=2A, r=epsilon/sqrt(1-epsilon^2),
kappa=-r/A, L=T+D, D>0, 0<epsilon<=1/4.
```

On the throat rho_0=log(P sec z), extended smoothly into a fixed exterior.
Only its smoothness is used in the cancellation, not any exterior
Einstein equation. Such extensions are prescribed two-dimensional
metrics, not asserted four-dimensional solutions.
The unsmoothed first-order proper-mouth pullback gives

```text
alpha=(1+kappa z tan z) 1_th,
beta=kappa(1+z tan z) 1_th,
gamma=kappa z 1_th.
```

Inside the throat, beta-gamma'-rho_0' gamma=0. Distributionally
gamma'=kappa 1_th+r(delta_-A+delta_A), so the only remaining term is
-r times the two endpoint delta functions. Squaring distributions is
not permitted. Instead define the whole family at positive h first.

Let G_h be the normalized periodic Gaussian of variance h^2, and set

```text
theta_h=G_h*1_th, d_h=G_h*(delta_-A+delta_A),
gamma_h=G_h*(kappa z 1_th),
alpha_h=theta_h+rho_0' gamma_h,
beta_h=kappa theta_h+rho_0' gamma_h.                       (8)
```

Use these profiles in (2). They are smooth and real, and agree with the
proper-mouth linear variation away from the seam as h->0. The exponential
lapse and spatial scale are positive. Formula (8) regulates the metric,
not a spectral term in isolation. It gives exactly

```text
gamma_h'=kappa theta_h+r d_h,
F_h=theta_h+r d_h,
C_h-rho_0' gamma_h=-r d_h.                                (9)
```

Thus (7) fixes the local contribution c r^2 w^2 integral d_h^2/(12pi).
The static background, contact terms and the profile in the spectrum
use the same metric prescription. Other smooth join families can have
additional finite contributions; no regulator-independence theorem for
arbitrary four-dimensional joins is claimed.

## 4. Exact cancellation and a return-uniform sharp limit

Define the following completed part of the response:

```text
R_h(w)=Delta chi_(r,h)(w)
                    +c r^2 w^2/(12pi) integral d_h^2 dz.    (10)
```

Delta chi is the full known free vacuum commutator of the preceding
note, including its k^3 and state-dependent k terms. Equation (10)
therefore combines that commutator with the derived anomaly contact;
it is not labeled the entire physical influence kernel.

For k_n=2pi n/L the endpoint coefficient is
(d_h)_n=2cos(k_nT/2)exp(-h^2 k_n^2/2)/L. Parseval gives

```text
integral d_h^2 dz=4/L+8S_h/L,
S_h=sum_(n>=1)cos^2(k_nT/2)exp(-h^2 k_n^2).                (11)
```

The spectral endpoint term contained
-2c r^2 w^2 S_h/(3pi L). Equation (11) cancels it exactly at every h,
including the overlapping images when D is small. The finite zero
harmonic 4/L must be retained. The result is

```text
R_h=Delta chi_(0,h)+Cross_h+c r^2 w^2/(3pi L)
 +2c r^2 w^2(w^2-nu^2)/(3pi L)
    sum_(n>=1)cos^2(k_nT/2)exp(-h^2 k_n^2)/(w^2-k_n^2),     (12)
```

where Cross_h is the cross term defined and bounded in the preceding
note. Its conditional sine/k sum is heat convolution of a bounded
sawtooth. The remaining sums are absolutely convergent as h->0.
For every fixed positive T,D,r and Im w>0 with Re w>0, (12) has a finite
limit. It equals the explicit finite part F(w;T,D,r) in equation (12)
of that note. Here that finite part follows from a covariant metric
calculation; it is not selected by discarding a divergence by hand.

With T and Re w in fixed positive compact bands, 0<eta=Im w<=1 and
0<h<=T_min/4, the earlier mode-counting and sawtooth estimates now give

```text
|R_h-Delta chi_(0,0)|
 <= C c [h+h^2/eta+(r+r^2)(1+eta^-1)].                    (13)
|R_0-Delta chi_(0,0)| <= C c(r+r^2)(1+eta^-1).             (14)
```

C is independent of D>0. Unlike the isolated spectral estimate, (13)
has no r^2/h term. Proof: apply the existing canonical Gaussian estimate
to the first term in (12), the sine/k decomposition to Cross_h, and
split the last sum at a fixed k>2|w|. Below that split the denominator
is at least eta Re w and the mode count is O(L); above it the k^-2
tail is O(L). The zero-harmonic term is bounded using L>=T_min.
These arguments also justify passage to the limit and (14).

For q channels, r comparable to epsilon=q^-2 and eta comparable to
delta=log^-2(q), (14) divided by q delta is O(epsilon/delta^2)->0.
The sharp limit can be taken first within this specified response part.
The former h=epsilon^4 spectral example is no longer an obstruction
after (10); its large terms cancel. This says nothing yet about the
size of the omitted interacting and geometric contributions.

## 5. Verification, prior work and next decision

The [checker](../checks/covariant_mouth_inertia.py) and
[receipt](../evidence/covariant-mouth-inertia.json) check curvature,
the spectral/polynomial identity, the Weyl kinetic square, pure-coordinate
and homogeneous conformal limits, smooth-profile identities, direct
curvature versus spectral response for independent Fourier profiles,
and spectral/real-space cancellation with Gaussian refinement. The
ordinary proof, not the finite samples, establishes uniformity in D.

The [source supplement](../evidence/covariant-mouth-source-review.json)
records the targeted searches and scoped primary reads.
[Adam--Rossini](https://arxiv.org/abs/hep-th/9710238v2), equations
(16)-(18), fixes the diffeomorphism-invariant Lorentzian action and
anomaly convention. [Mazur--Mottola](https://arxiv.org/abs/hep-th/0106151v1),
section III, gives the Weyl cocycle used in (1).
[Nguyen](https://arxiv.org/abs/2108.01095v3), section 2, explicitly
distinguishes covariant stress generation from local contacts and
chiral factorization. These are established anomaly mechanisms.
Neither their general formulas nor this application establishes a new
physical wormhole effect. Full source-dependency and novelty review
remain incomplete.

**Decision:** retain physical long-return damping as provisional primary,
with matched finite-duration transfer as fallback. The sharp free pressure
divergence is resolved for the anomaly-completed response part in (10).
Do not carry a physical smoothing-width requirement forward merely from
the earlier spectral divergence. The decisive remaining dynamical work
is the interacting proper-pressure response and the finite state/trajectory
contacts from the matched physical family. Angular/compact sectors,
corrected geometry and the exterior still need joint control.
No Lean/Lake or GCP invocation occurred. This supporting calculation is
not a replacement for the requested consequential result or manuscript.

Remaining to-do list: complete the state/trajectory contacts and
interacting pressure response, match corrected geometry and coupled
state/exterior, establish the physical endpoint and novelty, independently
verify, certify suitable mathematics on GCP, and prepare the manuscript.
