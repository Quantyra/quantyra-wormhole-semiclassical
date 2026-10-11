# A decay interval with finitely many return channels

2026-10-11, S053. In the leading free transparent JT/channel model,
the symmetric clock has a controlled linear decay interval in a joint
limit of finitely many channels. The sufficient condition below is
`delta -> 0`, `q delta^3 -> infinity`, where `delta=d/ell0`.
The interval contains any fixed number of damping times. Every individual
finite system remains conservative and recurrent. This closes a reduced
model question; it does not prove damping of a four-dimensional wormhole.

The starting response is derived in
[the transparent flux calculation](transparent-clock-flux.md).
Its continuum frequency and width were obtained in
[the collective-clock calculation](collective-clock-response.md).
Here the discrete response, its ultraviolet domain, and the passage from
frequency convergence to real time are treated explicitly.

## 1. Model and statement

Use dimensionless exterior time `t=t_phys/ell0`, `c=1`, `A=q/16`,
`M=2A=q/8`. Set `u_j=(j-1/2)/q`, `D_j=kappa delta/u_j`,
`kappa=1+asinh(sqrt(3))/(2sqrt(3))`, for `1<=j<=q`.
Let `x_q` be the unique positive stationary clock rate, `T_q=pi/x_q`,
`L_j=T_q+D_j`. The equation determining it is

`F_q(x)=x/8+1/24-(pi^2/(6q)) sum_j (pi+D_j x)^(-2)=0`.       (1)

It has `F_q'(x)>=1/8`. Put

`K0_q=M x_q^2+sum_j pi T_q D_j/(3 L_j^3)>0`,
`E_jn=2pi n/L_j`,
`W_jn=2(n^2-1)sin^2(pi n T_q/L_j)/(3 L_j^2 n)`, `n>=2`.

The exact retarded linear kernel is

`K_q(w)=K0_q-M w^2+sum_jn 2w^2 W_jn/[E_jn(w^2-E_jn^2)]`.  (2)

Let `H_q` denote the clock response to unit initial velocity, with all
positions and bath velocities initially zero, in the positive realization
below. Equivalently its retarded frequency transform is `M/K_q(w)`.
Let `K_infinity,delta` use the continuous uniform distribution of `u`
and its own equilibrium `x_delta`. Define `Omega_delta` to be the unique
near-one real zero of `Re K_infinity,delta(w+i0)`. This **exact real
center** is required; an unspecified `O(delta)` frequency error would
accumulate an order-one phase error on the stated time scale.

Write

`S=8log(2)-13pi^2/24>0`, `gamma=4kappa S/(3pi^2)`.

For any fixed `S_time<infinity`, as `delta->0`, `q delta^3->infinity`,

`sup_(0<=s<=S_time) |H_q(s/delta)`
` -exp(-gamma s) sin(Omega_delta s/delta)/Omega_delta| -> 0`. (3)

This is an absolute error theorem for a normalized linear Green function.
It supplies no explicit finite-`q` accuracy threshold. In physical units
the amplitude width is `Gamma=gamma d/ell0^2` and the interval has length
`S_time ell0^2/d`. The working hierarchy `delta=q^(-1/16)` satisfies
the sufficient condition, as does any `delta=q^(-a)`, `0<a<1/3`.
No nonlinear disturbance amplitude is asserted to remain controlled for
this growing time interval.

## 2. Positive realization without an infinite bare counterterm

Index all nonzero weights by `a=(j,n)` and define
`h_a=sqrt(2W_a/(M E_a^3))`. For finite `q`, `h` belongs to `l2`:
each channel's summand is `O(n^-2)`. On `C direct_sum l2` use the form

`Q[u,X]=(K0_q/M)|u|^2+sum_a E_a^2 |X_a-h_a u|^2`,
`Dom Q={(u,X): X-hu belongs to Dom E}`.                     (4)

The bounded invertible map `(u,X)->(u,X-hu)` takes this domain to
`C direct_sum Dom E`. Thus the form is dense, closed and positive.
Since each finite system has `min E_a>0` and `E_a->infinity`, the form
is coercive and its domain embeds compactly into the Hilbert space.
The associated positive self-adjoint operator `L_q` has compact resolvent.
Zero-weight oscillators can be omitted without changing any clock result.

The bare position vector `e0=(1,0)` generally is **not** in the form
domain. Expanding (4) into separate infinite terms is invalid. A velocity
kick at zero position is legitimate: its energy is `M v0^2/2`, and the
wave equation gives `phi(t)=v0 <e0,sin(t sqrt(L_q))/sqrt(L_q)e0>`.
This realizes the linear stress response, not additional physical matter
fields or a nonlinear CFT state construction.

To compute the resolvent, test its weak equation against `(0,Y)` and
against `(1,h)`, which is in the form domain. For `z=w^2` this yields

`X_a=E_a^2 h_a u/(E_a^2-z)`,
`[K0_q/M-z-z sum_a E_a^2 h_a^2/(E_a^2-z)]u=1`.

All sums in these equations converge. Consequently a probability measure
`nu_q` on positive frequencies satisfies

`<e0,(L_q-w^2)^(-1)e0>=M/K_q(w)`
` =int nu_q(domega)/(omega^2-w^2)`,
`int nu_q=1`, `int nu_q(domega)/omega^2=M/K0_q`.             (5)

The last identity also follows by solving the zero-frequency weak
equation: `L_q^(-1)e0=(M/K0_q)(1,h)`.
Compact resolvent makes `nu_q` pure point. Its characteristic function
`C_q(t)=int exp(-i omega t)nu_q(domega)` is a uniformly convergent
sum of exponentials, hence almost periodic. Simultaneous approximation
of any finite set of its frequencies, followed by the small positive
mass of the omitted tail, gives arbitrarily late returns close to one.
Indefinite exponential decay is therefore excluded for fixed `q,delta`.

## 3. A quadrature estimate at the damping resolution

The flux form of (2), including the static contact contribution, is

`K_q(w)=-M w^2+M x_q^2+sum_j C(D_j,w,T_q)`,
`C= -i/(12pi) (w-k_L^2/w)(1-a) F`,
`a=exp(iwT)`, `b=exp(iwD)`, `F=(1-b)/(1-ab)`,
`k_L=2pi/(T+D)`.                                          (6)

For `eta=Im w>0` and `r=1-exp(-eta T)`,

`F_D=-iw b(1-a)/(1-ab)^2`,
`int_0^infinity |F_D|dD <= |w||1-a|/(eta r)`,
`|F|<=2/r`, `int_0^infinity |(k_L^2)'|dD=4pi^2/T^2`.

The integral estimate follows by substituting `v=exp(-eta D)` in
the majorant with denominator `(1-|a|v)^2`. Therefore

`Var_D C <= V(w,T)`
` :=(12pi)^(-1)[8pi^2|1-a|/(|w|T^2 r)`
`       +|1-a|^2(|w|^2+4pi^2/T^2)/(eta r)]`.                 (7)

Monotone reparameterization `D=kappa delta/u` does not increase total
variation. The midpoint discrepancy of a function of bounded variation
is at most its variation divided by `2q`. The summed discrete response
thus differs from the averaged response at the same `T` by at most `V/2`.
For `w=Omega_delta+delta z` in any compact set of `Im z>0`, this is
`O(delta^-2)`, or `O(1/(q delta^3))` relative to `M delta`.

The static integrand in (1) has variation at most `1/6`, so
`|F_q-F_infinity|<=1/(12q)` and

`|x_q-x_delta|<=2/(3q)`.                                  (8)

This equilibrium correction cannot silently be omitted. Both rates tend
to one and their throat lengths differ by `O(1/q)`. For completeness,
at fixed `D` let `H=(1-a)(1-b)/(1-ab)`. Then

`H_T=-iw a(1-b)^2/(1-ab)^2`, `|H|<=4/r`,
`|H_T|<=4|w|/r^2`, `|partial_T k_L^2|<=8pi^2/T^3`.

Differentiating (6) gives the uniform bound

`|C_T| <= (12pi)^(-1)[32pi^2/(|w|T^3 r)`
`                    +4(|w|^2+4pi^2/T^2)/r^2]`.            (9)

The length correction summed over channels is also `O(delta^-2)`;
the direct `M x^2` correction is `O(1)`. Combining these statements,

`[K_q(Omega_delta+delta z)-K_infinity,delta(Omega_delta+delta z)]`
` /(M delta) -> 0`.                                       (10)

Equations (7)-(10) concern the complete complex kernel, not just a
frequency-smoothed positive absorption statistic. The estimate is
sufficient, not asserted sharp.

## 4. The continuum limit in the small frequency window

Here are the analytic details needed from the earlier continuum result.
For `w` near one and `T` near `pi`, the short-return expansion is
`Delta chi(D,w)=-w^2 D/(12pi)+O(D^2)`, uniformly with derivatives in a
fixed complex neighborhood. Integrating `kappa delta/D^2` from
`kappa delta` to a fixed small `D0` gives the real logarithmic term.
The `O(D^2)` remainder gives an analytic `O(delta)` contribution.

For the remaining lengths, the subtracted spectral representation at
`w=i` bounds its weighted positive measure by `O(delta)`: the exact
circle expression is bounded for `D>=D0` there, and `D^-2` is integrable.
Away from a fixed frequency interval about one, the same weighted bound
controls the response and all its derivatives on a smaller neighborhood.
Inside that interval its density is, for small enough `delta`,

`rho_delta(w)/q = kappa delta sin^2(wT/2)/3`
`                 *sum_(n>=2) (1-1/n^2)/(2pi n/w-T)^2`.     (11)

The series and its complex derivatives converge uniformly on a fixed
disk about one; their norms are `O(delta)`. Subtract its value at the
Cauchy pole in the local spectral integral. The remaining integrand is
regular, and the subtracted term is an explicit logarithm with endpoints
outside the smaller disk. Its upper boundary value continues with
imaginary part `-rho_delta`. This gives an analytic `O(delta)` remainder
with derivative `O(delta)` on that disk. It also justifies using the
local continuation in the earlier frequency expansion, without estimating
a pole by a small fixed imaginary regulator.

The static root and the real center consequently obey

`x_delta=1-[8kappa/(3pi)]delta log(1/delta)+O(delta)`,
`Omega_delta=1-[5kappa/(3pi)]delta log(1/delta)+O(delta)`,
`K_infinity,delta'(Omega_delta)/M=-2+o(1)`,
`rho_delta(Omega_delta)/(M delta)->2gamma`.

The derivative statement is for the continued kernel; its imaginary
derivative is `O(delta)`. Strict negative real derivative near one proves
the uniqueness of the real center. Taylor expansion, followed by (10),
now gives, locally uniformly for `Im z>0`,

`K_q(Omega_delta+delta z)/(M delta) -> -2(z+i gamma)`.        (12)

Also `K0_q/M->1`, directly from (1), (8), and the static length integral.
No lower-half-plane continuation of a finite discrete system is used.

## 5. From (12) to the real-time statement

Push the probability measure `nu_q` forward by
`y=(omega-Omega_delta)/delta` and call it `sigma_q`. For
`w=Omega_delta+delta z`, partial fractions in (5) give

`m_q(z):=int sigma_q(dy)/(y-z)`
` =2w delta M/K_q(w)+delta int nu_q(domega)/(omega+w)`
` -> -1/(z+i gamma)`.                                     (13)

The last integral is bounded by `1/Re w`. The limit is the Stieltjes
transform of the probability density `gamma/[pi(y^2+gamma^2)]`.
No mass can escape to infinity: at `z=iR`,
`R Im m_q(iR)=int R^2/(y^2+R^2)sigma_q(dy)->R/(R+gamma)`.
First taking the sequence limit and then `R->infinity` proves tightness
(outside `[-R,R]` the complementary integrand is at least one half).
Every weak subsequential limit has the displayed Stieltjes transform;
Poisson-kernel inversion identifies it uniquely. Thus the measures
converge weakly to this Cauchy law.

Tightness also proves uniform convergence of their characteristic
functions on compact `s` intervals: truncate to `|y|<=R`, use uniform
Lipschitz continuity there and a finite mesh of `s`, then remove the
small tail. Consequently

`sup_(0<=s<=S_time) |exp(i Omega_delta s/delta) C_q(s/delta)`
`                         -exp(-gamma s)| -> 0`.           (14)

Finally `H_q(t)=int sin(omega t)/omega nu_q(domega)`.
Weak convergence above implies that `nu_q` concentrates near one.
The inverse second moment in (5) tends to one, so on a set `B` of small
mass, Cauchy-Schwarz gives
`int_B omega^-1 nu_q <= sqrt(nu_q(B) M/K0_q)`.
On a small fixed neighborhood of one, `omega^-1` is uniformly close
to `Omega_delta^-1`. Hence
`int |omega^-1-Omega_delta^-1|nu_q ->0`.
Replacing that weight is uniformly valid for all `t`; (14) then proves
(3), including `s=0`. Positivity and the inverse moment are what prevent
uncontrolled very soft modes from invalidating this last step.

## 6. Prior work, verification and physical boundary

[Al-Qaiwani--Petrat, 2607.29485v1](https://arxiv.org/abs/2607.29485v1),
Theorem 1.1 and Corollary 1.1, already address exponential decay from
discrete Friedrichs spectra. Sections 2-3 were read completely. Their
regular coupling and quantitative Levy-distance theorem are not directly
assumed for the singular form (4). The Cauchy spectral mechanism is prior
art. Their Lemma 2.2 includes a Lipschitz endpoint estimate whose proof
uses an integral divergent at that endpoint (PDF page 11); our independent
analytic-density argument above does not use it. This localized defect
does not negate their result for strictly smaller Holder exponents.
[Derezinski--De Roeck, math-ph/0604058v2](https://arxiv.org/abs/math-ph/0604058v2),
introduction, is an older weak-coupling and dilation comparator; its full
extended-limit proof is not newly reviewed here. Neither paper is evidence
for the physical validity or novelty of the wormhole application.

An explicit check of the endpoint issue uses the probability density
`rho(y)=(1+|y|)/3` on `[-1,1]`. It is Lipschitz near zero, with
`rho(0)>0`, but its Poisson integral minus `pi rho(0)` is
`[2 atan(1/eta)-pi+eta log((1+eta^2)/eta^2)]/3`.
Dividing by `eta` diverges as `(2/3)log(1/eta)`; the claimed `O(eta)`
endpoint bound cannot hold with those assumptions. This is independently
derived and checked; no conclusion about the paper's other endpoints or
its entire main theorem is inferred from this counterexample.

Run `python checks/finite_channel_decay.py`; the
[receipt](../evidence/finite-channel-decay.json) records algebraic checks,
independent finite oscillator diagonalization, spectral/closed-form tail
checks, complex quadrature at shrinking resolution, equilibrium errors,
and a direct finite-channel spectral time calculation. Numerical plots or
fits are not substituted for the proof. The source review and endpoint
counterexample are recorded separately in the receipt and
[source audit](../evidence/finite-decay-source-review.json).

The finite-spectrum examples deliberately expose the time restriction.
At `delta=0.08`, `q=8192`, the sampled envelope error, including the
positive omitted spectral mass, is below `0.1` through one leading
lifetime, but reaches about `0.50` over three lifetimes. These are sampled
diagnostics, not uniform numerical certificates or evidence of indefinite
decay. The asymptotic theorem requires the joint limit, not merely a large
channel count at fixed `delta`. Both initial numerical failures are kept
in external custody: inadequate quadrature resolution, then an extremely
soft root too close to its bath pole for the chosen double-precision
bracket. The final time check uses a declared fixed frequency band and
bounds all omitted spectral mass and its inverse moment by positivity.

The proved length distribution here is the inverse-tail model. The actual
dipole's leading absorption coefficient alone does not transfer (10), its
real center, or (3). A physical transfer would require a positive clock
spectral representation with bounded inverse second moment and an
`o(M delta)` error in the full response in (12), including its real part.
The prior finite-cutoff matter estimate is compatible with this scale when
`epsilon/ell=o(delta^3)`; it is not a bound on the full matched gravity
and matter kernel. Interactions, the corrected background, full return
propagation, quantum state matching and nonlinear accumulation remain open.
This is a necessary reduced-model closure within the existing primary
damping/fallback pulse program, not a new accepted research target.
No Lean, Lake or GCP execution occurred.

Remaining to-do list: control interacting stress and the full physical
response at the damping resolution; close background, channel and state
matching; finish novelty and the consequential physical result; independently
verify, certify suitable endpoints on GCP and prepare the manuscript.
