# Transparent returns fix the leading clock response

2026-10-11, E003/S051-S053. The symmetric-clock retarded equation can be
derived from transparent propagation and the JT boundary energy balance.
It is not necessary to guess a finite part of an isolated-mouth action.
In the leading free JT/channel model, this derivation gives exactly the
kernel previously defined by a static subtraction in the
[collective-clock calculation](collective-clock-response.md). It fixes
the real response as well as absorption. It also constructs the causal
matter state for every prescribed smooth finite-time clock history with
positive bounded clock rate.

This closes the conservative-response ambiguity *within this leading
model*. The four-dimensional interacting theory, its quantum-corrected
background and physical channel return map remain to be controlled. The
result is not a full MMP stability or damping theorem. The continuum pole
found earlier is now a consequence of these transparent JT equations;
its interpretation as physical wormhole damping and its novelty remain
unaccepted. No arbitrary kinetic counterterm is selected here.

## Model and the quantities measured at the mouths

Keep the previously derived two-boundary JT coefficient `A=r_e^3/G`,
renormalized in the stated leading theory. The common global throat time
is `tau(t)`, `x=tau'>0`, with `phi=log x`. There are `q` independent
nonchiral CFT channels, each of central charge `c=1`, with fixed NS spin
structure and ground-state holonomy. Channel `j` has a static exterior
optical return time `D_j>=0`. In this section all `D_j` are finite.
The common-clock, reflection-symmetric state makes the two mouths carry
the same outgoing chiral flux `O_j(t)`. No angular or singlet interactions
are included in this free calculation.

Take the sharp asymptotic JT boundary first. A ray exiting the throat at
time `t` entered the opposite mouth at `p(t)`, where

`tau(p(t))=tau(t)-pi`, `p'(t)>0`.                              (1)

It arrived at that mouth from the exterior after time `D_j`. Thus its
previous outgoing-mouth time is `r_j(t)=p(t)-D_j`. The exterior clock is
the physical clock being held fixed; it is not reparametrized together
with the throat. The antiperiodic sign on completing a loop does not
affect the stress tensor.

Let `O_j,I_j` be the state-dependent chiral stresses in exterior optical
coordinates. The exterior is static, so its local geometric anomaly
contributes equally to the two directions and cancels in their energy
flux difference. Current conservation fixes the fermion's amplitude
factor to `sqrt(p')`. Point splitting then gives

`I_j(t)=O_j(t-D_j)`,
`O_j(t)=p'(t)^2 O_j(r_j(t))-c/(24pi) Sch(p,t)`.                (2)

The sign is fixed, for example, by mapping a line vacuum with
`p=exp(k t)`: its outgoing flux is `+c k^2/(48pi)`. This is a statement
about the normal-ordered chiral stress in the exterior frame, not a claim
that the covariant tensor transforms with an anomalous diffeomorphism.
The local geometric part restores covariance in the throat.

For either mouth the renormalized boundary energy and its Ward identity
are

`E_b=-A Sch(tan(tau/2),t)`
`   =-A[phi''-phi'^2/2+x^2/2]`,
`E_b'=sum_j [I_j-O_j]`.                                      (3)

The other boundary uses `-cot(tau/2)`, a Mobius transform of `tan(tau/2)`,
and has the same Schwarzian. Equations (1)-(3), with the JT charge sector
and initial data specified below, are the complete semiclassical mean
flux equations in this symmetric free model. They are established JT
and conformal-stress identities applied to closed delayed returns, not
an action for quantum clock fluctuations or a four-dimensional solution.

The two-mouth factor can be checked without linear response. The
state-dependent exterior energy is

`E_ext(t)=2 sum_j int_(t-D_j)^t O_j(s)ds + constant`.

Consequently `(2E_b+E_ext)'=0` whenever (3) holds. For a prescribed clock
that does not solve (3), its residual is the required external power;
energy conservation must not be asserted for such a drive.

## A finite-history state, not an instantaneous vacuum prescription

Initially take `tau=x_* t` for `t<=0`, with
`T=pi/x_*`, `L_j=T+D_j` and `O_j^0=-pi c/(12L_j^2)`. For a prescribed
smooth history, suppose `0<x_min<=tau'<=x_max` on each finite time
interval, with smooth matching to the static past. Equation (1) gives
`t-p(t)>=pi/x_max`. Recursive backtracking therefore reaches the static
past after finitely many steps on every finite observation interval.

Define the unwrapped initial phase `R_j(t)` by

`R_j(t)=t` for `t<=0`,
`R_j(t)=R_j(r_j(t))+L_j` for `t>0`.                           (4)

The construction is unique, smooth and increasing; the derivatives are
finite products and sums of those of the smooth positive-rate maps.
Smooth matching to the static past makes the definitions agree wherever
the number of backtracking steps changes. Characteristic propagation
extends the same construction from the mouths to the whole channel.

For a complex chiral fermion its two-point function is the pullback of
the initial NS vacuum kernel, with the half-density factors:

`G_j(t,t')=sqrt(R_j'(t) R_j'(t'))`
`          /[2i L_j sin(pi(R_j(t)-R_j(t')-i0)/L_j)]`.          (5)

One may equivalently use its positive-frequency NS mode sum, with the
appropriate field ordering. Positivity follows by pulling back that
positive two-point kernel. On a Cauchy circle the map is an
orientation-preserving diffeomorphism, so the half-density pullback
preserves the fermionic anticommutation relations. The short-distance
singularity is the vacuum one in every smooth spacetime patch. This
does not assert a global smooth-manifold Hadamard condition at an ideal
distributional seam. This argument is finite-time and conditioned on the
positive-rate clock; it does not prove global nonlinear backreaction or
prevent the clock from approaching zero.

Point splitting (5) yields

`O_j=O_j^0 R_j'^2-c/(24pi) Sch(R_j,t)`.                       (6)

The Schwarzian chain rule makes (6) satisfy (2) exactly. For clarity, the
local expansion underlying the sign is

`sqrt(R'(t+h/2)R'(t-h/2))/(R(t+h/2)-R(t-h/2))`
` =1/h + h Sch(R,t)/12 + O(h^3)`.

The cylinder sine adds `h pi^2 R'^2/(6L^2)` to the bracket. The physical
chiral energy is `-1/(2pi)` times that regular coefficient, giving (6).
Thus the radiation is fixed by causal initial data and return propagation;
the free chiral integration constants are not set to zero at each time.

## Linearization fixes the real kernel

Set `tau=x_*[t+psi(t)]`, so `delta phi=psi'`. From (1), to first order,

`delta p=psi(t)-psi(t-T)`,
`delta p'=phi(t)-phi(t-T)`.

Writing `O_j=O_j^0+u_j`, equation (2) becomes

`u_j(t)-u_j(t-L_j)`
` =2O_j^0[phi(t)-phi(t-T)]`
`  -c/(24pi)[phi''(t)-phi''(t-T)]`.                           (7)

For the retarded convention `exp(-i omega t)`, first take `Im omega>0`.
Backtracking is the convergent geometric series, and gives

`u_j/phi = c/(24pi)[omega^2-(2pi/L_j)^2]`
`                   *(1-exp(i omega T))/(1-exp(i omega L_j))`. (8)

Linearizing (3) and using `I_j=O_j(t-D_j)` gives `K(omega) phi=0`, for
nonzero frequency, with

`K=-2A omega^2+2A x_*^2`
`  -sum_j c/(6pi)[omega-4pi^2/(omega L_j^2)]`
`       *sin(omega T/2)sin(omega D_j/2)/sin(omega L_j/2)`.     (9)

There is no undetermined local mass in (9). Dividing the boundary energy
equation by frequency does discard an integration constant; it is fixed
by the initial JT charge sector and equilibrium, not by treating a
zero-frequency solution as an additional oscillating mode. We retain
the same zero total global-time charge as the earlier constrained-clock
reduction. The energy balance alone does not select the throat length.

The static constraint uses the covariant MMP potential established in the
[proper-clock calculation](proper-clock-anomaly.md):

`E_s(x)=A x^2+q c x/24-pi c/6 sum_j (pi/x+D_j)^-1`.

Its positive stationary point is unique. Dividing `D_phi E_s=0` by `x`
gives

`2A x+q c/24-sum_j c pi^2/[6(pi+D_j x)^2]=0`.                (10)

The left side starts at `-q c/8`, is strictly increasing, and tends to
positive infinity. At its root,

`K_stat=D_phi^2 E_s(x_*)`
`       =2A x_*^2+sum_j c pi T D_j/(3L_j^3)>0`.              (11)

The zero-frequency limit of (9) is exactly (11). At equilibrium the
static total-energy check is also exact:

`E_s(x_*)=-A x_*^2-sum_j pi c D_j/(6L_j^2)`
`        =2E_b+E_ext`,                                     (12)

up to the same fixed exterior reference constant. These checks connect
the integration constant and static anomaly to the dynamical equation.

Using the previously derived independent fermion spectrum,

`K(omega)=-2A omega^2+K_stat+sum_j Delta chi_j(omega)`,
`Delta chi_j=sum_(n>=2) 2omega^2 W_jn/[E_jn(omega^2-E_jn^2)]`,
`E_jn=2pi n/L_j`,
`W_jn=2c(n^2-1)sin^2(pi nT/L_j)/(3L_j^2 n) >=0`.             (13)

This converges for `Im omega>0`. The static contact term in (11) cancels
the constant subtraction in the earlier closed formula, leaving (9).
The agreement is a consequence of the flux derivation; the subtraction
is not now an independent postulate of the leading transparent model.

There are no growing linear-response poles in this symmetric free model.
For `z=omega^2` with nonzero imaginary part, (13) implies

`Im[K(omega)/z] = -Im z [K_stat/|z|^2`
`                +sum_jn 2W_jn/(E_jn |z-E_jn^2|^2)]`,       (14)

which is strictly nonzero. On the remaining upper-half-plane axis
`omega=i gamma`, `K(i gamma)>0` for `gamma>0`. These statements exclude
all upper-half-plane zeros, not just purely growing roots. They concern
the symmetric mean clock, not every gravitational or gauge mode.
At finite channel count a passive closed spectrum need not produce
irreversible decay. Linear pole exclusion is not a finite-disturbance
survival theorem.

## Finite proper-time cutoff: an actual ray-vertex error

The sharp ray relation has a controlled finite-cutoff replacement before
any four-dimensional gluing. For the two induced-proper-time curves of
the preceding note, take `sigma=+-s(t)` and

`cos^2 s=epsilon^2[x^2-sdot^2]`,
`tau(t)-tau(p_epsilon(t))=s(t)+s(p_epsilon(t))`.              (15)

At the static rate `x_*`, put
`a=acos(epsilon x_*)`, `b=sqrt(1-epsilon^2 x_*^2)`,
`T_epsilon=2a/x_*`, `L_epsilon,j=T_epsilon+D_j`.
The first-order boundary displacement is
`delta s=-epsilon x_* phi/b`. Therefore

`delta p_epsilon=psi(t)-psi(t-T_epsilon)`
`                +(epsilon/b)[phi(t)+phi(t-T_epsilon)]`.     (16)

The extra endpoint term is required by the induced metric. Repeating the
stress transfer calculation gives (8) with `T,L` replaced by
`T_epsilon,L_epsilon` and its numerator replaced by

`B_epsilon=1-exp(i omega T_epsilon)`
`          -i omega epsilon/b [1+exp(i omega T_epsilon)]`.    (17)

This is a finite-cutoff *matter* response, not the full finite-cutoff
gravitational equation. It includes the static change of the initial
circle vacuum in `L_epsilon`.

A useful uniform bound follows without resolving individual long-return
poles. Fix `m<=|omega|T<=W`, `Im omega>0`, and `epsilon x_*<=1/4`, with
`T=pi/x_*`. Set `d_0=1-exp[-(Im omega)T/2]`. Since

`T_epsilon>=T/2`, `0<=T-T_epsilon<=3epsilon`,
`|B_epsilon-(1-exp(i omega T))|<=6epsilon|omega|`,
`|1-exp(i omega L_epsilon)|>=d_0`,
`|(2pi/L_epsilon)^2-(2pi/L)^2|<=192pi^2 epsilon/T^3`,

subtracting (8) explicitly proves, for each channel's contribution
`C=2i(exp(i omega D)-1)u/(omega phi)` to (9),

`|C_epsilon-C| <= C_(m,W) c epsilon/T^2 (d_0^-1+d_0^-2)`.    (18)

Here `|exp(i omega D)-1|<=2` for every `D>=0`; differentiating the
denominator adds at most `3epsilon|omega|/d_0^2`. The constant is finite
and independent of the return length. The script evaluates the explicit
pre-simplification triangle bound, not a fitted constant.

At damping-scale resolution `(Im omega)T` of order `delta`, (18) is
`O(c epsilon/(T^2 delta^2))`. Its total relative error against
`q delta/T` tends to zero under the sufficient, deliberately conservative
condition `epsilon/T=o(delta^3)`. The previously exhibited choice
`delta=q^-1/16`, `epsilon/ell=q^-1/4` satisfies this stronger condition:
the ratio is `O(q^-1/16)`. It also meets the earlier nominal classical
overlap requirement. This controls the free ray vertex at complex
resolution; it does not bound four-dimensional channel conversion,
interactions, the finite-cutoff gravitational response or the corrected
reference geometry.

## Consequence and contribution comparison

For the inverse-tail continuum `D=kappa d/u`, `0<u<=1`, equations (9)
and (10) are exactly the previously studied closed continuum system.
Hence its locally continued retarded resonance has, with
`x0=qc/(16A)`, `delta=d x0`,

`Re omega_p/x0=1-[5kappa/(3pi)]delta log(1/delta)+O(delta)`,
`-Im omega_p/x0=[4kappa/(3pi^2)]`
`              *[8log 2-13pi^2/24]delta+o(delta)`.            (19)

The earlier continuum calculation supplies the limiting integral and
continuation proof; it is not rerun or silently strengthened here. The
point-dipole absorption comparison and uniform singlet absorption bound
also retain their original scopes. Equation (19) now has a transparent
leading-JT derivation, but neither a finite-channel decay interval nor
the interacting four-dimensional resonance follows from it.

[Engelsoy--Mertens--Verlinde v3](https://arxiv.org/html/1606.03438v3),
sections 3.1, 3.2, 6.3 and appendix A, supplies the boundary dynamics,
energy and global-frame comparison. There is a displayed sign conflict
between its equations 6.15, 6.17, 6.19 and 6.20: substituting the first
three gives the negative of the last right-hand side. The plus sign in
6.19 was confirmed in the PDF on page 24, not inferred from HTML alone.
We use the direct point-split sign above.
[Almheiri--Engelhardt--Marolf--Maxfield v3](https://arxiv.org/html/1905.08762v3),
sections 2.1, 2.2 and 3.1, independently fixes the incoming-minus-outgoing
balance (22), the negative stress Schwarzian (21), and the requirement
to identify the physical times at a transparent interface. These are
established tools, not new results of this project.

[The recent polar-metric/axial-gauge paper](https://arxiv.org/html/2609.33511v1)
is a distinct, close predecessor. Its section 5.4 explicitly discusses
free chiral state data; its conclusion restricts the analysis to the
`l>=2` sector, drops ambient returns, and selects an instantaneous
vacuum. Here the target is the common boundary clock with finite delayed
returns and a causal initial state. That distinction does not by itself
prove novelty. This continuation read the state-choice discussion and
the complete conclusion, and screened the stability criterion and
selected closed-throat mode calculations; the entire paper and its proof
dependencies have not newly been reviewed.

The [independent verification](../checks/transparent_clock_flux.py) checks
point splitting, ray inversion and its derivatives, nonlinear causal
stress against the linear response, the flux/trigonometric/spectral
representations, static energy and stiffness, the pole-exclusion
identity, and the finite-cutoff bound. The
[receipt](../evidence/transparent-clock-flux.json) binds the exact sources;
the [source supplement](../evidence/transparent-flux-source-review.json)
records the scoped citation and novelty searches. No full physical
novelty claim or Lean/GCP run is made.

Remaining to-do list: control the interacting stress and finite-channel
decay window, construct the quantum-corrected four-dimensional background
and channel/state matching, finish primary/fallback novelty and the
consequential physical result, independently verify and GCP-certify
suitable endpoints, and prepare the reproducible manuscript.
