# Quadratic electric feedback: positive response and retained quantum noise

> Resumed and reviewed, 2026-10-13. The original checker reproduces its
> saved receipt byte for byte. The bounded stationary-response and
> all-return stiffness arguments are retained. The continuation derives
> [a minimum cancellation duration and a shaped switch](electric-return-cancellation.md);
> physical driven-state and stress matching remain unfinished.

S072, 2026-10-13. **The finite-channel stationary electric response has no
exponentially growing poles within the phase-linear, second-order neutral
influence model.** The proof includes all strip returns and any radial
Gaussian bath with the previously derived positive spectral measure.
For a specified cubic switch, the completed-square stiffness also stays
positive through every return in the selected family. Quantum covariance
and energy checks retain the neutral and radial noises.

These are properties of the resummed quadratic model, not a uniform
higher-order junction theorem. Stationary response, positivity during a
switch and stability of the driven state are different statements. No
full wormhole stability or creation result follows yet.

## The actual projected source, to its accepted order

Return to the electric part of the junction in
[neutral-transfer-current.md](neutral-transfer-current.md), with angular
dynamics held fixed. In common global strip time its interaction is
`-h[X cos(alpha)+Y sin(alpha)]`, with `mu=<X>=q/(2pi)`. The exact neutral
vacuum spectra are `v_E(q)` at positive even energies and `w_E(q)` at
odd `E>=3`. They obey

```
0<=v_E<=q E/pi^2,   0<=w_E<=q(E-1/E)/pi^2,
0<=1-v_E/(q E/pi^2)<=4H_E/q       (even E).              (1)
```

The last inequality is the even-channel proof in
[collective-electric-response.md](collective-electric-response.md),
valid at every energy, though its relative upper bound need not be small
at arbitrarily large energy. The common electric oscillator has already
been removed from these spectra. It is not added a second time.

Expanding the phase through quadratic order, and the neutral influence
through `h^2`, requires both `Y` response and the change of `<X>`:

```
j=-mu h alpha -h alpha(chi_XX*h)+h[chi_YY*(h alpha)]+noise.
chi_AA(t)=2Theta(t)sum_E weight_A(E) sin(Et).             (2)
```

In Laplace space `chi_AA(z)=sum_E 2E weight_A(E)/(z^2+E^2)` before
subtraction. Equation (2) does not retain higher connected neutral
correlators or terms nonlinear in the phase. Resumming (2) does not
prove that those omitted terms are small.

## The finite-channel restoring difference is positive

Use the same Abel regulator in both static susceptibilities. Their
difference has a finite regulator limit for `q>2`:

```
delta_q=lim_(epsilon->0)[chi_XX,epsilon(0)-chi_YY,epsilon(0)]
       =2 integral_0^infinity [C_XX(tau)-C_YY(tau)] d tau.
```

The inherited exact Euclidean kernels give

```
C_XX-C_YY=4q b^2[exp(-(1-1/q)L)+(q-1)exp(L/q)-q],
b=1/(4pi), L=2log coth(tau/2).                         (3)
```

The bracket is nonnegative by convexity, since the weighted linear
exponents sum to zero, and is strictly positive away from `L=0`.
Its short-time singularity is `tau^(-2/q)`; its long-time tail is
integrable. Thus `delta_q>0` for every integer `q>=3`. This is not a
subtraction of two separately assigned infinite static numbers.

Putting `x=tanh(tau/2)` evaluates the same integral:

```
delta_q=q/pi^2 integral_0^1
       [x^(2-2/q)+(q-1)x^(-2/q)-q]/(1-x^2) dx
       =q/pi^2 { (q/2)[psi(1/2)-psi(1/2-1/q)]-q/(q-2) }.
delta_q/mu -> pi/2-2/pi.                               (4)
```

A common finite local contact cancels in this difference. A different
contact in the two kernels would change the model. The `q=2` integral
diverges logarithmically; that case is not included in this theorem.

This calculation also detects a nonuniform cutoff trap. At the candidate
`q=2000^4`, (4) gives `delta_q/mu=0.9341765544...`; keeping 512 even
modes and 512 odd modes gives about `0.29694`. At large `q` such paired
hard cutoffs approach `pi/2-4/pi` before resolving the finite-`q` tail.
They cannot set the continuum restoring coefficient. We keep (4)
explicitly when testing a finite dynamic band.

## Positive oscillator completion, including radial fluctuations

For the neutral current, let `c_E=sqrt(2E w_E)` and introduce canonical
oscillators `(y_E,p_E)`. Their vacuum force
`Y_0=sum c_E y_E^0` has exactly the ordered spectrum `w_E`.
The electric radial sector retains the prior positive representation

```
H_G=Q_0^2/(2I_0)+sum_a (P_a^2+omega_a^2 X_a^2)/2,
alpha=alpha_0+sum_a sqrt(sigma_a)X_a,
G_alpha(z)=1/(I_0 z^2)+integral d sigma(omega)/(z^2+omega^2),
I_0>0, d sigma>=0.                                     (5)
```

Use its physical normalization and a common clock before combining it
with the strip kernels. Positivity is preserved by positive clock/unit
rescalings; actual clock and vertex matching remain required. Equation
(5) can retain the entire positive RN radial measure, not only a static
inertia. At finite regulators a Hamiltonian realizing (2) is

```
H_quad(t)=H_G +k(t) alpha^2/2
          +(1/2)sum_E [p_E^2+E^2(y_E-h(t)c_E alpha/E^2)^2],
k(t)=mu h+h^2 delta_q-h R_X(t),
R_X(t)=sum_even (2v_E/E) integral_0^t cos[E(t-u)] h'(u) du. (6)
```

At a finite common regulator, replace `delta_q` in this exact realization
by `delta_reg=chi_XX,reg(0)-chi_YY,reg(0)` and regulate the kernels
consistently. The displayed expression is the common-regulator limit.
The finite-band demonstrations instead match their local static term
to `delta_q`; they are not a raw truncation of both susceptibilities.

Assume `h(0)=0`. Integrating the `X` response by parts gives
`chi_XX*h=h chi_XX(0)-R_X`. Eliminating each driven `y_E` in (6)
therefore recovers both terms of (2), with the common static contact
canceled. Its neutral noise is `h(t)Y_0(t)`, with symmetrized covariance
`h(t)h(u) Re C_YY(t-u)`. Homogeneous radial fluctuations are also
retained through (5). Neither noise is reset after an interaction.

This oscillator construction is an established method. For an explicit
derivation of the bath elimination, homogeneous noise and propagated
covariance, see [Boyanovsky and Jasnow, section II through equation II.27](https://arxiv.org/html/1707.04135v2).
Their Drude specialization and Born/Markov approximations are not imported.
The finite-`q` difference (4) and its use in (6) are derived here from
the already accepted strip reference.

## Stationary response: no growing poles, with all returns retained

For stationary preparation at constant `h>0`, the restoring coefficient
in the completed square is

```
k_stat=mu h+h^2 delta_q>0.
K_eff(z)=1/G_alpha(z)+k_stat
          +h^2 sum_odd (2w_E/E) z^2/(z^2+E^2).          (7)
```

The inverse phase response is `K_eff`. For `Re z>0`, each
`z/(z^2+omega^2)=(1/(z-iomega)+1/(z+iomega))/2` has positive real part.
Thus `z G_alpha(z)` has strictly positive real part, as does its inverse.
Dividing (7) by `z` gives

```
Re[K_eff(z)/z]
 =Re[1/(z G_alpha(z))+k_stat/z
       +h^2 sum_odd (2w_E/E) z/(z^2+E^2)] > 0.          (8)
```

The sums converge locally uniformly in this half-plane: (1) bounds the
neutral tail by a multiple of `sum E^(-2)`, while the prior radial
measure grows at most linearly. Consequently (8) proves that `K_eff`
has no zero there. This rules out exponentially growing stationary
linear modes, without a Markov approximation or removal of delayed
strip returns. It does not prove irreversible damping; discrete baths
can retain oscillations. At finite cutoff the same result follows from
the positive Hamiltonian and its real normal frequencies.

The continuum statement concerns retarded response. The sum inside a
completed-square Hamiltonian can require a dressed representation, and
bare equal-time phase or momentum covariances need not converge in the
original Fock representation. No continuum quantum state or stress is
constructed merely from (8).

## A single switch has positive stiffness through every return

Take the stated global-time comparison

```
h(t)=h_f[3(t/T)^2-2(t/T)^3], 0<=t<=T; h(t)=h_f after T,
0<T<pi.
```

For `Z_E(t)=int_0^t cos[E(t-u)]h'(u)du`, two integrations by parts give
`|Z_E(t)|<=A/E^2`, with
`A=24h_f/T^2` uniformly for all `t>=0`. This includes the jump of
`h''` at the end of the cubic ramp as part of its total variation.
Equation (1) and `sum_(n>=1) H_n/n^2=2zeta(3)` imply

```
|R_X,q(t)-R_X,infinity(t)| <= (16zeta(3)/pi^2)A.         (9)
```

Here the infinity label denotes the comparison even spectrum
`v_E=q E/pi^2`, not an interchange of the full static cutoff limits.
The weighted difference in (9) is absolutely summable at every `q`.
Its comparison response follows from the cosine Dirac comb:

```
R_X,infinity(t)=mu h'(t)+2mu sum_(m>=1) h'(t-m pi)
                                      -(2mu/pi)h(t).   (10)
```

Terms with negative arguments vanish. Because `T<pi`, returned ramp
packets do not overlap. At most one delayed derivative contributes.
The direct packet has coefficient one; later packets have coefficient
two. Using `delta_q>0` and `max h'=3h_f/(2T)`, (6)--(10) give

```
k(t)/(mu h(t)) >= 1-3h_f/T-epsilon_q    for all t>0,
epsilon_q=(768zeta(3)/pi) h_f/(q T^2).                  (11)
```

Before the first return, replace `3h_f/T` by `3h_f/(2T)`.
For the reference `q=k^4`, `h_f=4pi/k`, `T=0.25`, integer `k>=2000`,
both losses decrease with `k`. At `k=2000`,

```
epsilon_q=1.84636e-12,
pre-return k(t)/(mu h(t)) >= 0.9623008881,
all-time   k(t)/(mu h(t)) >= 0.9246017763.              (12)
```

Thus the instantaneous quadratic Hamiltonian has positive stiffness
through every return for this comparison. This is not a proof of
bounded driven trajectories: its periodic coefficients can in principle
produce parametric amplification. Nor does the switched state become
the stationary preparation of (7) when `h` stops changing. That
distinction must be retained in the full stability calculation.

## Quantum checks and the next physical step

[Checker](../checks/passive_electric_feedback.py),
[receipt](../evidence/passive-electric-feedback.json), and
[source/dependency record](../evidence/passive-electric-feedback-source-review.json).
Six channel counts compare (4) with an independent endpoint-regularized
integral. Three positive oscillator Hamiltonians reproduce (7) by direct
matrix inversion at 20 complex frequencies each, and have strictly
positive squared normal frequencies. Explicit tail bounds cover the
omitted response. Three switches test (9) through ten returns; (11)
supplies the uniform statement beyond those numerical samples.

Two finite-band quantum switches retain 12 neutral bath modes and two
illustrative radial modes, with 128 `X` response modes and the continuum
static difference matched by (4). The full `(alpha,Q_0)` covariance
obeys its uncertainty determinant bound `>=1/4`. Keeping only the
homogeneous initial phase contribution instead gives minima `0.1733`
and `0.05846`, while the full minima are about `0.2655`. This explicitly
repairs the missing-noise test from the earlier deterministic switch.
Direct symplectic evolution also balances the finite-band energy change
against `int <partial_t H_quad>dt`, including mixed phase/bath work.
Time-step tolerance refinement changes the final covariance by less
than `3.03e-9`.

For the `q=1000` finite-band example, the returned Hamiltonian is
periodic with period `pi`. A direct monodromy calculation and tolerance
refinement find its largest multiplier modulus equal to one within
`3e-15`. This is a bounded numerical comparison, not a proof for the
candidate's full spectrum or a uniform Floquet theorem. No finite-band
covariance or unsubtracted oscillator energy is promoted to the
renormalized continuum switching stress.

The Gaussian linear electric sector now has a stationary stability
argument, an all-return stiffness bound and quantum noise accounting.
The next calculation must control the nonstationary coupled state and
higher connected influence, including angular/coset branching and the
massive constrained interactions. Initial dressing, common stress
subtraction and evolving geometry remain necessary. The physical
state, clocks and vertex cannot be supplied by the illustrative bath
or the reference cubic switch. No full stability, formation or novelty
acceptance is claimed.

Remaining to-do list: bound the driven quantum junction and its nonlinear
remainder with matched state, clocks and vertices; include angular and
massive constrained interactions, evolving geometry and complete stress;
establish preparation, sufficient stability and finite-energy transport.
