# Massive current work on global AdS2, including quantum sources

S072, 2026-10-13. **The full massive sphere tower has a curved, positive
current-measure bound in the global AdS2 vacuum.** It controls quadratic
energy deposition by the radial Maxwell energy seminorm, uniformly in
frequency and without a radial supremum penalty at the AdS boundaries.
A second-order extension permits noncommuting quantum sources and
time-dependent radial profiles, under explicit initial-factorization
and finite-moment assumptions. Those assumptions have not been matched
to the actual wormhole source. Evolving geometry and full stability remain
open.

This transfers the integrated-work method of
[the flat longitudinal calculation](massive-bulk-polarization.md) to an
exact curved reference. It does not use a small-curvature expansion.
The proof bounds absorptive excitation energy after a source pulse;
it does not recompute the static Maxwell counterterm or complete stress.

## Geometry, modes and the electric source

Use the fixed product

```
ds^2=a^2 sec^2(sigma)(-dt^2+d sigma^2)+b^2 dOmega_2^2,
-pi/2<sigma<pi/2,        C=b^2/a^2,
M_n^2=n(n+q)/C,          d_n=2(q+2n), n>=1.
```

Assume `q>=1` and `M_1>=4`. The massless LLL is excluded, having been
retained already in the light current sector. The current candidate has
`M_1=18.048...` at `k=2000`. Spinor rescaling gives the one-particle
Hamiltonian and its square

```
H_M=-i sigma_2 partial_sigma+M sec(sigma) sigma_1,
H_M^2=diag[-partial_sigma^2+W^2-W',
           -partial_sigma^2+W^2+W'],   W=M sec(sigma). (1)
```

For these masses the normalizable realization is unique and the scalar
operators in (1) are their Friedrichs realizations. The positive energies
are `E_j=M+1/2+j`, `j>=0`. A normalized positive mode has components

```
u_j^+ = N_j cos^M(sigma) (
 sqrt(1+sin sigma) P_j^(M-1/2,M+1/2)(sin sigma),
 sqrt(1-sin sigma) P_j^(M+1/2,M-1/2)(sin sigma) )^T,
N_j=sqrt[j! Gamma(j+2M+1)]/[2^(M+1/2) Gamma(j+M+1/2)],
u_j^-=sigma_3 u_j^+ .                                (2)
```

See [Serrano Blanco, equations 36, 48, 66, 87–88 and the massive-mode discussion](https://arxiv.org/html/2208.08252v1).
The unique normalizable domain and invariant positive-frequency splitting
specify the vacuum used here. They do not specify the actual RN or
interacting prepared state. Both sphere eigenvalue signs have the same
scalar-source bound and are included in `d_n`.

Start with a smooth bounded scalar electric potential in radial gauge,
`A_t=f(t)V(sigma)`, `A_sigma=0`. Define

```
Q[V]=integral_(−pi/2)^(pi/2) cos^2(sigma)|V'(sigma)|^2 d sigma,
U_E(t)=(2pi C/g^2) f(t)^2 Q[V].                     (3)
```

The second expression is its Maxwell Killing energy in global `t` units,
from the positive coefficient `kappa=4pi b^2/g^2`. A constant `V` is a
pure charge source and does not create pairs. Finite-energy profiles
near the boundaries enter through `Q`, rather than `||V'||_infinity`.
Matching the coefficient `g` and the radial gauge/source to the physical
four-dimensional setup remains necessary.

## A local heat bound with the radial weight retained

Let `delta=pi/2-|sigma|`. For either component of (1),

```
W^2 +/- W' >= M(M-1)sec^2(sigma)
            >= M(M-1)/delta_boundary^2.             (4)
```

Here `delta_boundary` is distance to either chosen endpoint, throughout
the interval. Compare with the half-line inverse-square operator using
the nearer endpoint at the evaluation point. Potential monotonicity and
domain killing in the scalar Feynman–Kac formula bound its diagonal heat
kernel by

```
K_M(tau;delta,delta)
 = delta/(2tau) exp(-z) I_(M-1/2)(z),
z=delta^2/(2tau).                                   (5)
```

This kernel follows directly from the normalized Bessel spectral integral
`delta integral_0^infinity k exp(-tau k^2) J_(M-1/2)(k delta)^2 dk`;
[DLMF 10.22.67](https://dlmf.nist.gov/10.22.E67) gives its closed form.

An elementary bound supplies useful mass decay. In the integral for
`I_nu(z)` from [DLMF 10.32.2](https://dlmf.nist.gov/10.32.E2), set
`v=1-t`. For `nu=M-1/2`, use
`(2-v)^(M-1)<=2^(M-1) exp[-(M-1)v/2]`, then extend the integral
from `[0,2]` to `[0,infinity)`. Its gamma integral gives

```
K_M <= (4pi tau)^(-1/2)
          [1+(M-1)tau/delta^2]^(-M)
    <= (4pi tau)^(-1/2)
          [1+3M^2 tau/(16delta^2)]^(-4).             (6)
```

For the second step, Bernoulli's inequality applied with exponent `M/4`
gives `(1+x)^(M/4)>=1+Mx/4`, and `M-1>=3M/4`. This is why the
theorem assumes `M>=4`; no asymptotic Bessel approximation is used.

The exact sphere counting function has the bound

```
D(Y)=sum_(M_n^2<=Y) d_n <= 2C(1+1/q)Y.              (7)
```

Indeed the first `N` levels have `2N(q+N+1)` states, while
`N(N+q)<=CY` and `N<=CY/q`. Integrating a decreasing function against
this positive measure shows
`sum d_n(1+h M_n^2)^(-4)<=2C(1+1/q)/(3h)`.

The local positive-frequency heat trace is half the full spinor trace,
because `u_j^-=sigma_3 u_j^+` has the same pointwise norm. Applying
(5)–(7) to both scalar components therefore gives

```
sum_n d_n sum_j exp(-tau E_nj^2)|u_nj^+(sigma)|^2
 <= [16/(9sqrt(pi))] C(1+1/q) delta^2 tau^(-3/2).
```

The spectral-projector inequality with `tau=3/(2L^2)` yields

```
N_+(L;sigma) <= A C(1+1/q) delta^2 L^3,
A=16/(9sqrt(pi)) (2e/3)^(3/2)=2.446850... .           (8)
```

All sphere and radial modes are included. The bound is local in `sigma`,
so multiplication by a field profile preserves the useful boundary weight.

## Current conservation turns the projector bound into work

Write `V_jk=<u_j^+,V u_k^->` for a fixed sphere channel and
`E_jk=E_j+E_k`. The exact commutator is

```
E_jk V_jk=<u_j^+,[H_M,V]u_k^->,
[H_M,V]=-i sigma_2 V'.                              (9)
```

This identity keeps the vector Ward cancellation; treating all potential
matrix elements as independent bounded numbers would discard it. Define
the positive energy-weighted pair measure

```
mu_V((0,L])=sum_n d_n sum_(E_jk<=L) E_jk |V_jk|^2.
```

In (9), replace `1/E_jk` by `1/E_j`, extend the hole sum, and use
Parseval. This gives

```
mu_V((0,L]) <= integral |V'|^2
                       sum_(n,j:E_j<=L) d_n |u_nj^+|^2/E_j d sigma.
```

Stieltjes integration of (8) bounds the weighted local sum by
`(3/2)A C(1+1/q) delta^2 L^2`. Finally,
`delta^2<=pi^2 cos^2(sigma)/4`. Hence the central curved estimate is

```
mu_V((0,L]) <= D L^2 Q[V],
D=B C(1+1/q),
B=3pi^2 A/8=9.056041... .                           (10)
```

The response is a discrete positive measure; it is not a pointwise
bounded function on the real frequency axis. In particular a discrete
resonance is not replaced by a flat-space threshold. Start with smooth
bounded potentials so all commutators are legitimate. Equation (10)
extends the quadratic pair form continuously to the completion in `Q`,
modulo spatial constants. This extension does not by itself establish
nonlinear dynamics for an arbitrary unbounded external potential.

For a scalar pulse vanishing with its first derivative at temporal
endpoints, define `F0=||f||_1`, `F2=||f''||_1`. Its quadratic massive
excitation work is `W^(2)=integral |f_tilde(E)|^2 dmu_V(E)`.
The envelope `|f_tilde(E)|^2<=min(F0^2,F2^2/E^4)`, integration by parts
against (10), and positivity give

```
W^(2) <= 2D Q[V] F0 F2.                             (11)
```

For `f(t)=[1-(t/T)^2]^3` on `|t|<=T` and zero outside, with peak one,

```
F0 F2=12288/(875sqrt(5)),
W^(2)/U_peak <= 18.104088... g^2(1+1/q).             (12)
```

This is uniform in duration and radial profile. For the reference
`g=k^-4`, `q=k^4` comparison it is `7.07191e-26` at `k=2000`.
The denominator is the pulse electric energy in (3), not the negative
opening energy. The coefficient is less sharp than in flat space, but
the AdS2 radial operator and its boundaries are retained exactly.

## Nonseparable and quantum sources at second order

The same pair bound also gives a useful statement without assigning a
classical trajectory to the phase. Let the initial massive bath be its
global vacuum, factorized from a source state `rho_s`. In the interaction
picture of the uncoupled bath and the chosen source Hamiltonian, let the
Hermitian radial potential be an operator `A(t,sigma)` on the source.
It may fail to commute at different times. Define ordered norms

```
X_0=integral dt [integral cos^2(sigma)
          Tr(rho_s (partial_sigma A)^dagger partial_sigma A) d sigma]^(1/2),
X_2=integral dt [integral cos^2(sigma)
          Tr(rho_s (partial_t^2 partial_sigma A)^dagger
                       partial_t^2 partial_sigma A) d sigma]^(1/2). (13)
```

Assume these are finite and that temporal integration by parts has no
boundary terms. At finite regulator the pair-creation map `T(V)` has
components `sqrt(E_jk)V_jk`. Equation (10), applied to each source matrix
element, says

```
||P_(<=L) T(A) rho_s^(1/2)||_HS
 <= sqrt(D) L ||partial_sigma A rho_s^(1/2)||_(weighted L2,HS).
```

The coefficient of second-order bath energy is the squared norm of
`integral dt exp(i H_pair t) T(A(t))rho_s^(1/2)`. The frequency projectors
are orthogonal. Its low band `E<=L` is at most `D L^2 X_0^2`.
In each band `2^j L<E<=2^(j+1)L`, integrate twice in time; its norm
is at most `sqrt(D)2^(1-j)X_2/L`. Squaring and summing all high bands,
then optimizing `L`, proves

```
W_bath^(2) <= inf_(L>0) D[L^2 X_0^2+(16/3)X_2^2/L^2]
            = (8/sqrt(3)) D X_0 X_2.                 (14)
```

The minimum uses `L^2=4X_2/(sqrt(3)X_0)`. Zero norms are
handled by a limit. Positive cutoff convergence defines the quadratic
response when the norms in (13) are finite. The proof permits changing
radial profiles and uses the ordered source two-point function, not a
commuting stochastic replacement or a symmetrized covariance.

The assumptions are substantive. A bare quantum field need not have a
finite `X_2`; renormalization or a source cutoff cannot be hidden inside
that norm. The physical source may be initially correlated with dressed
massive modes by gauge constraints. The source's feedback and higher
orders in its coupling to the massive bath also need a remainder bound.
Equation (14) establishes none of those missing conditions. It is an
explicit target for the actual radial quantum source, not its solution.

## Checks and the next physical transfer

[Checker](../checks/curved_massive_response.py),
[receipt](../evidence/curved-massive-response.json), and
[source/dependency record](../evidence/curved-massive-response-source-review.json)
retain the evidence. Exact Jacobi modes check (1), their normalization,
(9), and the vanishing pair amplitude of a spatially constant potential.
Independent Bessel integration fixes (5). Positive mode sums test the
heat comparison, and complete finite-energy pair sums test (10) and (11),
including a profile concentrated toward the radial boundaries.

For the explicit pulse, the exact Fourier transform from the preceding
checkpoint has `|f_tilde(E)|<=D0/(T^3 E^4)` above `E=L`, with
`D0=96[1+6/(LT)+15/(LT)^2+15/(LT)^3]`. Equation (10) bounds every omitted
pair above that cutoff by `4D Q[V]D0^2/(3T^6 L^6)`. This is separate from
spatial quadrature refinement; no finite sphere or radial truncation is
promoted to the whole work.

A two-level quantum source coupled to one pair sector supplies an
independent check of ordering. Its initially excited state samples
`|f_tilde(E-gap)|^2`, while its ground state samples
`|f_tilde(E+gap)|^2`. Exact source-plus-pair evolution approaches these
different second-order predictions with the expected refinement. The
test is not a model of the wormhole source.

The next step must match the actual state and radial potential to (13),
or replace those norms if their ultraviolet limit fails. It must also
control the departure from the fixed normalizable global product: the
RN exterior, boundary transfer, evolving sphere/dilaton and initial
thermal/entangled preparation. The theorem does not cover a persistent
interaction by declaring it a pulse, or turn positive bath excitation
energy into the full renormalized gravitational stress. No physical
creation, complete stability or novelty acceptance is claimed.

Remaining to-do list: establish the actual quantum source moments and
initial-state matching; control radial/exterior transfer, evolving geometry,
feedback and complete stress; establish preparation, sufficient stability
and finite-energy transport.
