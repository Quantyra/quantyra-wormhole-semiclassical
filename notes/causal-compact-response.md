# Causal charge and rotation response before the opening

S072, 2026-10-13. **The projected Gaussian charge and axial sectors have
a positive constrained wave operator and a quantitative causal response
bound.** For an infinitesimal cubic boundary-transfer drive lasting as
long as the proposed opening switch, the single-species candidate's
charge and angular responses follow their static values to within
`2.33e-6` of the final static response, uniformly in time, including
the slightly nonextremal RN benchmark at the prescribed initial
temperature. This is a response test on a fixed RN exterior. The actual
fermion junction, prepared quantum state and evolving wormhole stress
have not yet been matched to it.

The [static calculation](angular-boundary-response.md) fixes the round
nonrotating AdS boundary metric, magnetic flux and zero electric dipole
source. Those electromagnetic and metric sources remain fixed here.
The new input is an explicitly prescribed phase of the bosonized matter
at the boundary, whose derivative is a transfer chemical potential. It
permits charge or angular momentum to enter through a boundary current.
This boundary-transfer problem must be supplied by the eventual
two-sided junction. An isolated reflecting system instead conserves its
total charges and has a different nonzero-frequency charge response.

## Constrained dynamical sector

Keep `R=r_+`, `X=R^2/L^2`, `Y=Q^2/R^2`,
`tau=1+3X-Y=4pi R T>=0`, `P=q/2` and `Q^2=pi G q^2/g^2`.
The two-dimensional metric is `-f dt^2+dr^2/f`. For one angular axis,
promote the static fields to one-forms `w_a`, `v_a`, with
`w_t=omega`, `v_t=V` in the preceding note. Their radial components
must be retained to impose the constraints. In the electromagnetic
gauge with vanishing dipole `A_theta`, the quadratic Lorentzian action is

```
a=(w,v)^T, n=(1,P)^T, c=(-P,1)^T, d=(1,2)^T/(q+1),
kappa_g=r^4/(6G), kappa_e=4pi r^2/(3g^2), K=diag(kappa_g,kappa_e),
mu=8pi/(3g^2), lambda_J=S_J/pi=N_f q(q^2-1)/(12pi),

S = integral sqrt(-g2) [
    - (kappa_g F_w^2+kappa_e F_v^2)/4
    - mu (c^T a)^2/2 - lambda_J (d^T a-d chi)^2/2 ].        (1)
```

Here `d chi` denotes the derivative of the scalar `chi`; it is distinct
from the column vector `d`. The gauge transformation is
`a -> a+n d alpha`, `chi -> chi+alpha`. In particular `c^T n=0`,
`d^T n=1`. The projected fermion source `d^T a` is exactly
`(w+2v)/(q+1)`, rather than the metric shift alone. The static restriction
`chi=Omega t` recovers the preceding chemical-polarization functional.

The positive Schwinger coefficient comes from the massless Dirac current
sector in [MMP appendix A](https://arxiv.org/html/1807.04726v3). This use is
quadratic in each angular generator. Non-Abelian vertices, higher Landau
levels, the other angular harmonics and background quantum/thermal stress
are outside (1). They remain required for physical transfer.

For `q>1`, `lambda_J>0`, fix unitary gauge locally. Define the constant
matrices

```
M=mu c c^T+lambda_J d d^T=L L^T,
L=(sqrt(mu)c, sqrt(lambda_J)d),  det M=mu lambda_J>0,
p=K(partial_r a_t-partial_t a_r),  p=L Z.
```

The two constrained first-order equations are
`p_x=M a_t`, `p_t=f M a_r`, where `dx/dr=1/f`, `x(infinity)=0`.
Eliminating the potentials gives the canonical two-component equation

```
partial_t^2 Z + H Z = 0,
H=-partial_x^2+V(x),   V=f L^T K^-1 L,                     (2)

V/f = [[ A, B ], [ B, D ]],
A=2/r^2+4Q^2/r^4,
B=sqrt(mu lambda_J)/(q+1) [-P/kappa_g+2/kappa_e],
D=lambda_J/(q+1)^2 [1/kappa_g+4/kappa_e].                  (3)
```

No extra gravitational dipole wave has been introduced. With the matter
polarization removed, the gravitational equation fixes a conserved
angular integration constant. The remaining electromagnetic wave has
potential `f(2/r^2+4Q^2/r^4)`. This agrees with the exceptional-mode
equation (74) of [Kodama and Ishibashi, section 4.2.2](https://arxiv.org/html/hep-th/0308128v3),
at sphere dimension two, after electromagnetic duality. Their equation
allows a cosmological constant. The Schwarzschild limit is the usual
electromagnetic dipole potential `2f/r^2`.

For clarity, the earlier Burko dipole paper was also inspected. Its
equations (26)--(27) and the subsequent physical-mode assignment in the
available v1 text conflict with the constrained equation and the neutral
Maxwell limit used here. That assignment is not used as verification.
The source review records the exact reading scope; no correction or
broader claim about that paper is being published.

Since `K` is positive and `M=L L^T`, `V` is positive semidefinite throughout
the exterior. In fact

```
AD-B^2=mu lambda_J/(kappa_g kappa_e)>0.
```

For homogeneous boundary-transfer data, `Z_x(0)=0`. At the horizon use
the causal ingoing solution. The canonical energy

```
E=(1/2) integral_(horizon)^0
   [|Z_t|^2+|Z_x|^2+Z^T V Z] dx                           (4)
```

is nonnegative. Its only loss for a causal exterior solution is flux
into the horizon; the zero-source boundary contributes no work. The
self-adjoint half-line operator with the homogeneous boundary condition
is nonnegative, and its retarded resolvent has no upper-half-plane poles.
Equivalently, an exponentially growing finite-energy mode is excluded by
the energy identity. This is linear stability of this Gaussian RN
sector, not stability of the subsequently opened wormhole.

## The boundary-transfer response and its static limit

Prescribe `chi_boundary(t)=alpha(t)` and write `Omega=alpha_dot`, while
keeping the physical boundary `w_t=v_t=0`. In unitary gauge the boundary
potentials are `a_t(0)=-n Omega`. Equation (2)'s boundary forcing is thus

```
Z_1,x(0)=0,  Z_2,x(0)=-sqrt(lambda_J) Omega,
J=-n^T p(0)=-sqrt(lambda_J) Z_2(0).                        (5)
```

The radial constraint accounts for `J_dot` through the matter boundary
current. This is an open transfer channel with a specified source, not
an assumption that a conserved isolated charge changes under a pure
time-dependent chemical gauge transformation.

For Laplace variable `s>0`, define the boundary impedance `K_J(s)` by

```
K_J(s)= min_(Z_2(0)=1)
        integral_-infinity^0 [|Z_x|^2+Z^T V Z+s^2|Z|^2] dx. (6)
```

The heavy boundary value is varied, giving `Z_1,x(0)=0`; the solution
decays at the horizon for `s>0`. At `s=0` use the bounded regular
finite-energy solution. The response is

```
J(s)=I_J(s) Omega(s), I_J(s)=lambda_J/K_J(s),
I_J(0)=I_fix.                                              (7)
```

The last equality follows either from (5), or from the static
first-order flux equations and `dF_min/dOmega=J` in the preceding note.
It therefore uses the actual fixed-field coefficient, including the
electromagnetic angular charge, rather than `I_KN`.

The spherical electric current has an analogous single scalar:

```
lambda_e=N_f q/pi, kappa_0=4pi r^2/g^2,
V_e=f lambda_e/kappa_0,
Q(s)=I_e(s) mu_e(s), I_e(0)=K_e,eff.                       (8)
```

Electric-angular cross terms vanish at this quadratic order because
`sum_m m=0`; angular axes are equivalent by rotational symmetry. These
local statements do not specify the compact charge-spin lattice or the
initial TFD's full transfer covariance.

## A bound on the entire causal response

The angular mixing obeys, at every radius,

```
|B|/A <= epsilon,
epsilon^2=4lambda_J/[mu(q+1)^2]
         =N_f g^2 q(q-1)/[8pi^2(q+1)].                    (9)
```

To see this, put `h=kappa_e/kappa_g>=0` in
`B/A=sqrt(lambda_J/mu)(2-P h)/[(q+1)(1+P^2 h)]`.
For `q>=1` the absolute value of its rational factor is at most two.
Assume `epsilon<1`, as in the candidate. The electric case has
`epsilon=0`. In either case put

```
I_0=I(0), lambda=lambda_J or lambda_e,
a_time=I_0/[lambda(1-epsilon^2)].
```

Then

```
I_0/(1+a_time s) <= I(s) <= I_0,     s>0.                  (10)
```

Here is the full ordinary proof. Let `Z0` minimize (6) at zero frequency,
normalized by `Z0_2(0)=1`. Its norm squared is convex because
`(|Z0|^2)''=2|Z0'|^2+2Z0^T V Z0>=0`. Its derivative tends to zero at
the regular horizon, so its maximum is at the boundary. The scalar
maximum principle for `-partial_x^2+f A`, with the natural homogeneous
boundary condition, gives
`sup|Z0_1|<=epsilon sup|Z0_2|`. This also follows by multiplying the
equations for `Z0_1 +/- epsilon sup|Z0_2|` by their positive parts;
the horizon and boundary terms vanish. Therefore

`sup|Z0|^2<=1+|Z0_1(0)|^2<=1/(1-epsilon^2)`.

Adding the positive `s^2` term gives `K(s)>=K(0)`. For the upper bound,
use the admissible trial field `exp(sx)Z0`. Integration by parts with
the zero-frequency equation gives exactly

```
E_s[exp(sx)Z0]=K(0)+2s^2 integral_-infinity^0 exp(2sx)|Z0|^2 dx
             <=K(0)+s/(1-epsilon^2).
```

Combining with `I=lambda/K` proves (10). The single-scalar electric proof
is the same, with `sup|Z0|<=1`. No spectral gap or replacement of the
horizon by a reflecting wall is used in this argument.

To turn this into a real-time statement, use the positive spectral measure
of the Neumann operator in (2). Its boundary Green function is
`I(s)/lambda`. Since `I_0` is finite, the normalized measure `dnu(w)`
obeys

```
integral dnu=1,
I(s)/I_0=integral_0^infinity [w^2/(w^2+s^2)] dnu(w),
integral [s^2/(w^2+s^2)] dnu <= a_time s.                  (11)
```

Boundary evaluation can be defined first at a finite regulator and then
as the half-line Green function; it is not assumed to be a square
integrable vector. Finite static susceptibility excludes a zero-mode
atom in the normalized measure. Equations (10)--(11) bound the
retarded response without replacing a real frequency by its absolute
value in a complex resolvent.

Start with zero perturbation and drive
`Omega(t)=Omega_f a(t/Delta)` with `a(v)=3v^2-2v^3` on `[0,1]`, zero
before and one afterwards. The same argument applies to `mu_e` and `Q`.
The spectral convolution is

```
[J(t)-I_0 Omega(t)]/(I_0 Omega_f)
 = - integral dnu(w) integral_0^min(t/Delta,1)
                    a'(v) cos[w(t-Delta v)] dv.            (12)
```

Monotonicity bounds the inner integral by one. Two integrations by parts,
using `a'(0)=a'(1)=0` and the total variation of the piecewise second
derivative, bound it by `24/(w Delta)^2`, uniformly for all `t>=0`.
Use `s=sqrt(24)/Delta` and
`min(1,s^2/w^2)<=2s^2/(w^2+s^2)` in (11). This proves

```
sup_(t>=0) |J(t)-I_0 Omega(t)|/(I_0 |Omega_f|)
 <= 4 sqrt(6) a_time/Delta.                               (13)
```

The bound is relative to the final static response, so it makes no
claim of a small relative error at the instant the drive starts from
zero. It covers the response after the ramp too. It is a linear mean
response theorem; noise and the charged fermion two-point function are
additional quantities.

## Quantitative opening-timescale test

Use the selected single-species reference family

```
eta_f=k^-1, g=k^-4, L/l_p=k^5, q=k^4, N_f=1, k>=2000,
X(1+3X)=pi k^6, C=1+6X, N=2q, phi_r=b^3/(G C).
```

The prior opening has duration `Delta=delta phi_r/(N eta_f)`,
`1/4<=delta<=1/2`. Let `d_upper` be the proved upper bound for
`G I_fix/b^3`, and `epsilon_pol` the earlier spherical screening bound.
Then the two response-time ratios have the explicit upper bounds

```
a_J/Delta <= 24pi eta_f C d_upper/
             [delta(q^2-1)(1-epsilon^2)],
a_e/Delta <= 8pi eta_f C(1+3X)(1+epsilon_pol)/(delta q^2).
```

Both are asymptotic to `48pi^2/(delta k^3)`. At `k=2000` and the fastest
allowed ramp `delta=1/4`, both are below `2.369e-7`; (13) is below
`2.33e-6`. Slower ramps and larger `k` improve the bound. This uses the
entire static interval, not the smaller trial matter term alone.
For nonzero initial temperature use the actual `R,T,Y`, static coefficient
and physical duration in (9)--(13). The theorem has no small-temperature
assumption. A quantitative extension to the prepared-temperature RN
benchmark is short. Set `T=N eta_f z_0/(2pi phi_r)`, `z_0=1/2`, and
`h=R/b`. At fixed magnetic flux,

```
4pi Tb h^3=h^2-1+3X(h^4-1) >= 2(1+6X)(h-1).
```

For this family `X>=1` and the RN temperature increases with `h>=1`.
The selected temperature is below its value at `h=2`, so
`h-1<=16pi Tb/(1+6X)`. At `k=2000` this is below `7.64e-33`, and it
decreases as `k` grows. The static electric trial bound at `R` is at
most `h` times the bound at `b`. The angular bound is at most `h^3`
times its reference value: its derivative terms grow at most as `R^3`,
while its mixing and matter trial integrals are `R` times functions
decreasing with `X_R`. Using the deliberately loose envelope
`h<=1+1e-12` still gives (13) below `2.33e-6` in both sectors. No floating
point subtraction of nearly equal horizon radii is needed.

This conclusion is uniform for `k>=2000`, not an extrapolation from the
three numerical examples. The exact flux equation gives `X>1e9` and
`X^2<pi k^6/3`; (9) gives `epsilon^2<1e-14`. The static trial bound obeys
`d_upper<X+1.501`, using `D(X)<=1/X`. Substitution bounds both
nonextremal errors by
`4sqrt(6) [48pi^2/(delta k^3)](1+1e-8)`, which is below `2.33e-6`
at the smallest `k` and `delta` and decreases with either parameter.

The earlier single-species scale receipt omitted `z_0` from its displayed
`Tb` field and hence overstated that temperature by two. That field and
its checker are corrected here. Its inertia bounds, phase exponents and
candidate decision are unchanged. The present finite-temperature argument
concerns classical RN coefficients and the projected chemical response;
thermal/background quantum stress corrections remain unbounded.

This is a positive test of the **boundary-transfer response timescale**.
The opening angle `theta(t)` rotates the fermion endpoint junction;
`Omega(t)` and `mu_e(t)` here are infinitesimal conjugate transfer drives.
The physical junction must determine those drives and their fluctuations.
The following are still required before applying the free JT capture
proof: the joint compact-state and fermion matrix elements, the noise
and mixed correlators in the prepared pair, the resulting renormalized
switching stress, and control of the evolving geometry and omitted modes.
The response function is not automatically an off-shell compact-rotor
action or an isolated-system Kubo susceptibility. These identifications
require the boundary variational problem and operator insertions.

## Verification and research decision

[Checker](../checks/causal_compact_response.py),
[receipt](../evidence/causal-compact-response.json), and
[source scopes](../evidence/causal-compact-source-review.json).
Symbolic checks verify the gauge/source projection, determinant, positive
Schur complement and constrained electromagnetic master potential.
Thirty-five direct integrals check the exact cubic-memory expression and
its envelope. Twelve finite passive RN regulators compare independent
linear solves with mode sums, test the Laplace bounds and reconstruct
the real-time response. Their horizon cuts reflect; they check the
spectral proof and discretization, and are not substitutes for the
analytic causal half-line argument. Grid refinement and six candidate
scale bounds are retained. No growing-mode claim relies on a finite
search for eigenvalues; it follows from (4).

The previous turn was progress on static matching. This turn is progress
on the constrained retarded problem and a uniform finite-time diagnostic.
Retain the single-species candidate and continue to its actual quantum
junction and stress. There is no completed formation, full stability,
initial-pair preparation, novelty acceptance, publication, outreach or
Lean certificate.

Remaining to-do list: match the compact state and actual two-sided
fermion-transfer interaction, including noise and switching stress;
bound evolving-geometry and other four-dimensional errors; establish
preparation, sufficient stability and finite-energy transport.
