# Schwinger singlet and unequal exterior returns

2026-10-10, S052/S053. This note closes a **defined Gaussian graph loading
problem**, not the four-dimensional quantum response. The throat photon
cannot in general be included by replacing `q` independent loops with
`q-1`: exterior propagation with unequal lengths does not preserve the
throat's neutral subspace. We give the exact coupled return matrix, a
finite-resolution bound and its reflecting-singlet limit. None establishes
physical damping or a new wormhole stability theorem.

The starting point is the [constrained clock calculation](collective-clock-response.md)
and [shrinking-sector propagation](long-return-sector.md). All previous
independent-channel response formulas retain their stated model scope.
This is a necessary refinement of that investigation, not a new matter
sector or a selected separate graph-theory paper.

## 1. Model and source applicability

Take `q>=2` canonical real scalar current fields on a common optical throat
interval `[0,T]`, with one exterior edge of length `D_j>0` connecting the
two ends of channel `j`. Let

`v=(1,...,1)/sqrt(q)`, `P=v v*`, `N=I-P`.

The quadratic action is

`S = (1/2) int_throat [|dot phi|^2-|phi'|^2-m(s)^2 |v*phi|^2] dt ds`

`    + (1/2) sum_j int_exterior_j [(dot phi_j)^2-(phi_j')^2] dt ds`.

At each mouth, fields and derivatives match channel by channel. Assume
real nonnegative bounded `m^2`; a singular-endpoint limit would require
its own domain argument. The model retains only the spherical photon
inside the throat. It omits exterior photon propagation, higher angular
gauge exchange, other induced interactions and the dynamical metric
vertex. It specifies linear oscillator propagation, not the fermion
spin structure, compact-boson charge sectors, or a four-dimensional state.

For the convention `-F_ab F^ab/(4g^2)` and unit fermion charge, spherical
reduction of `ds^2=N_lapse(s)^2(-dt^2+ds^2)+R(s)^2 dOmega^2` gives the
electric term `2 pi R^2 E^2/(g^2 N_lapse^2)=E^2/(2 e_2^2)`.
Canonical bosonization couples `E sum_j phi_j/sqrt(pi)`. Completing the
electric square at fixed zero background electric-flux sector gives

`e_2^2=g^2 N_lapse^2/(4 pi R^2)`,
`m(s)^2=q e_2^2/pi=q g^2 N_lapse^2/(4 pi^2 R^2)`.

In the leading AdS throat, `N_lapse/R=sec(s/ell)/ell` in a centered
coordinate. Thus `m>=M=g sqrt(q)/(2 pi ell)` on any truncated interval.
For `k~1/ell`, `k/M=O(1/(g sqrt(q)))`. This establishes a low-frequency
singlet barrier in this reduction; it does not prove the stated matching
conditions for the physical mouths. The ideal `T=pi ell` free-channel
length is used only when comparing scales, not to assert regularity of
the singular AdS endpoints.

The underlying Schwinger mechanism is established in
[Maldacena, Milekhin and Popov v3, Appendix A](https://arxiv.org/html/1807.04726v3#A1):
one common photon gaps the charged combination and leaves neutral current
degrees of freedom. Appendix G treats field-line bosonization for
screening. The more recent [reduction, section 4.3.1](https://arxiv.org/html/2606.12528v2#S4.SS3.SSS1)
also separates the singlet and neutral modes, with boundary conditions and
additional interactions requiring care. Neither source supplies the
unequal-delay stress response used in the present proposed damping test.

The [polar perturbation paper, sections 2.2 and 5.4](https://arxiv.org/html/2609.33511v1)
restricts its master perturbation system to angular momenta at least two;
the symmetric throat-length dynamics are excluded. Its treatment of
radiative state data does not construct our finite-time global state.
These are scoped source checks, not completion of the full novelty audit.

## 2. Exact passive throat and feedback

Let `z=eta-i k`, `eta>0`. Define the scalar singlet Dirichlet-to-Neumann
matrix `Lambda_m(z)` by solving

`f''=(z^2+m^2)f`, `Lambda_m (f(0),f(T))=(-f'(0),f'(T))`.

There is no Dirichlet eigenvalue in this half-plane. Green's identity gives

`Re[b* Lambda_m(z) b/z] = Re(z) int |f|^2`
`                       + Re(1/z) int (|f'|^2+m^2 |f|^2) >=0`.

Write the endpoint value as incoming plus outgoing amplitudes, `b=a+c`,
and the outward derivative as `z(a-c)`. The scattering matrix is

`S_m(z)=(Lambda_m+z I)^(-1)(z I-Lambda_m)`.

Its denominator is invertible and the preceding accretivity identity
makes this Cayley transform contractive. It is analytic for `Re z>0`.
For `m=0`, `S_0(z)=exp(-zT) J`, where `J=[[0,1],[1,0]]`.
For constant mass, an independent closed expression is

`Lambda_m=a [[coth(aT),-csch(aT)],[-csch(aT),coth(aT)]]`,
`a=sqrt(z^2+m^2)` with continuation at removable singularities.

Using port-major ordering on `C^2 tensor C^q`, the full throat and exterior
maps are exactly

`S=S_0 tensor N + S_m tensor P`,
`S_free=S_0 tensor I_q`,
`E=J tensor diag(exp(-z D_j))`.

The sum over any number of returns converges:

`G=(I-E S)^(-1)`, `G_free=(I-E S_free)^(-1)`.

Put `r=exp(-eta D_min)<1`. Both inverse norms are at most `1/(1-r)`.
The resolvent identity, with its order retained, is

`G-G_free=G E (S-S_free) G_free`.

Since `S-S_free=(S_m-S_0) tensor P` has rank at most two and norm at most
two, it follows that

`rank(G-G_free)<=2`,
`||G-G_free||<=2r/(1-r)^2`,
`||G-G_free||_1<=4r/(1-r)^2`.                         (1)

In particular, a bounded boundary readout `B` obeys
`|Tr[B(G-G_free)]|<=||B|| 4r/(1-r)^2`. The estimate is independent of
`q` at fixed `eta,D_min` and includes all repeated returns. It also holds
with any common exterior operator of norm at most `r`; the explicit
secular formula below uses the diagonal delays. The rank statement does
**not** imply uniform smallness on the real axis or a bound for an
unbounded quadratic stress observable.

## 3. Low-energy reflecting limit

Suppose `m(s)>=M>k>0`, and set `a=sqrt(M^2-k^2)`. At `z=-ik`, the real
quadratic form defining `Lambda_m` dominates the constant-`a` form.
Minimizing at fixed endpoints therefore gives

`Lambda_m(-ik)>=a tanh(aT/2) I`.

The matrix is real symmetric. Its scattering matrix is unitary, and
diagonalizing it gives

`||S_m(-ik)+I||<=2k/[a tanh(aT/2)]`.                (2)

Consequently its fixed-frequency large-mass limit is `-I`, not zero.
The singlet reflects separately at the two mouths. Neutral modes still
transmit, so the limiting throat map is

`S_inf=S_0 tensor N-I_2 tensor P`.

In this paragraph the limit means fixed `T,k` with `M` increasing. The
bound (2) is real-frequency local scattering. It is not automatically
a bound on the closed return inverse near its poles; (1) remains the
separate finite-resolution statement.

For equal exterior lengths, one common orthogonal channel rotation
diagonalizes both throat and exterior. The limiting system consists of
`q-1` free circles of length `T+D`, **plus a massless singlet exterior
interval of length `D` with Dirichlet ends**. The massive singlet throat
has decoupled. For unequal lengths the commutator of the delay matrix
with `P` is generically nonzero, so no such global channel decomposition
exists. Reflection creates coherent channel conversion on each return.

## 4. Secular equation and independent eigenvalue bound

In the reflecting limit, exchange of the two mouths has eigenvalue
`sigma=+1` or `-1`. Let `t=exp(-zT)`, `d_j=exp(-zD_j)`. In each parity
sector,

`I-E_sigma S_inf,sigma=I-t D+(t+sigma) D P`.

The determinant lemma gives

`det(I-E_sigma S_inf,sigma)`
` = product_j(1-t d_j) [1+(t+sigma)/q sum_j d_j/(1-t d_j)]`.  (3)

The right side is continued through its apparent poles; individual free
zeros can cancel. Use the full matrix determinant at such points.
For `k>0`, zeros at `z=-ik` are the graph frequencies, including modes
with vanishing endpoint values. At `k=0` the traveling-amplitude
description degenerates; count the kernel by the quadratic form instead.

There is an independent variational check. Let `H_free` be the direct
sum of `q` circle Laplacians on lengths `T+D_j`, initially with no mass.
Restrict its `H^1` form domain by the two scalar endpoint constraints

`sum_j phi_j(0)=sum_j phi_j(T)=0`.

They are independent continuous constraints. The restricted operator
`H_K` splits orthogonally into the reflecting graph `H_inf` and a
Dirichlet singlet interval on the **throat** of length `T`. To see this,
rotate only the common throat fields into `v` and its orthogonal
complement. The singlet has zero endpoints and decouples in both the
energy and `L^2` norm; the remaining variational mouth conditions are
precisely neutral transmission and singlet reflection.

The min-max principle now proves, counting from the first eigenvalue
with multiplicity,

`lambda_j(H_free)<=lambda_j(H_K)<=lambda_(j+2)(H_free)`.

Indeed, restriction raises the infimum. Intersecting the first `j+2`
free eigenspace with the two-constraint subspace leaves dimension at
least `j`, proving the upper bound. This is an application of established
finite-codimension bracketing, explicitly discussed after Theorem 5.3 in
[Berkolaiko and Kuchment, section 5](https://arxiv.org/html/1008.0369v2#S5);
it is not a new interlacing theorem.

For `omega>0`, write `n^+(omega)` for the number of positive frequencies
at most `omega`. The free kernel has dimension `q`, while the limiting
graph kernel has dimension `q-1`. The Dirichlet throat interval has
`floor(omega T/pi)` positive frequencies. Hence

`n_free^+(omega)-floor(omega T/pi)-1 <= n_inf^+(omega)`
` <= n_free^+(omega)-floor(omega T/pi)+1`.           (4)

This retains the extra throat interval in the bracketing comparison.
Dropping it would give an incorrect count. Equation (4) controls mode
counts, not their throat stress matrix elements or excitation weights.

## 5. The oscillator stress observable, including its contact term

For this section define a specific clock-dependent Hamiltonian on fixed
edge coordinates: `H(phi)=exp(phi) V+H_ext`, where

`V=(1/2) int_throat [pi_field^2+|phi_field'|^2+m^2 |v*phi_field|^2]`.

This scales the complete canonical throat energy, including its singlet
mass term. It is a specified reduced-model family, not a derivation of
the total four-dimensional metric vertex. The field called `phi_field`
here is distinct from the scalar clock `phi`. Fix the compact-boson
charge/zero-mode sector and quantize only the strictly positive oscillator
modes. Let their real orthonormal eigenfunctions be `u_a`, with
frequencies `omega_a>0`. Define

`A_ab=int_throat u_a* u_b`,
`B_ab=int_throat [u_a'* u_b'+m^2 (v*u_a)(v*u_b)]`,
`D_ab=B_ab-omega_a omega_b A_ab`.

Expansion of `V` in creation/annihilation operators gives a creation-pair
coefficient `D_ab/(4 sqrt(omega_a omega_b))` in the **ordered** double
sum. The norm of the resulting vacuum excitation, including identical
bosons, therefore gives

`rho_V(Omega)=(pi/8) sum_(a,b) |D_ab|^2/(omega_a omega_b)`
`                         * delta(Omega-omega_a-omega_b)`, `Omega>0`. (5)

The diagonal pair has norm factor two; an off-diagonal unordered pair
has two equal amplitudes. These produce the same ordered-sum formula.
For any compact positive output-frequency band, only finitely many
positive modes of a fixed compact graph contribute. This spectral part
requires no ultraviolet subtraction. It does not include charge-sector
transitions, a thermal population or a time-dependent prepared state.

With a finite common regulator the full retarded susceptibility is

`chi_V(w)=sum_(a,b) |D_ab|^2/(8 omega_a omega_b)`
` * [1/(w-omega_a-omega_b)-1/(w+omega_a+omega_b)]`, `Im w>0`.

The clock force is `exp(phi) V`. Its linear response is
`[<V>+chi_V(w)] phi`, where the first term is the contact term from
`H''(0)=V`. Here

`<V>=(1/4) sum_a [omega_a A_aa+B_aa/omega_a]`.

Removing the regulator in the real part and in this contact term needs
a common vacuum-energy/local-action renormalization prescription. We
have not supplied one for the physical graph/throat system. The contact
term is real and instantaneous, so it does not alter (5) on a strictly
positive frequency band. Omitting it from a stiffness calculation would
nevertheless be wrong.

Two nontrivial checks fix the interpretation. If `V` is the energy of
the entire graph, `A=I`, `B_ab=omega_a^2 delta_ab`, so `D=0`: a change
proportional to the total Hamiltonian cannot create vacuum pairs. For
one free circle of length `L=T+D`, explicit sine/cosine modes instead
give at `Omega_n=2pi n/L`, `n>=2`,

`W_n=2(n^2-1) sin^2(pi nT/L)/(3 L^2 n)`.

One can sum the two real polarizations at fixed positive mode numbers
`r,s` first. Their ordered pair weight is
`4 r s sin^2(pi(r+s)T/L)/(L^2(r+s)^2)`; then use
`sum_(r=1)^(n-1) r(n-r)=n(n^2-1)/6`. Thus (5) independently reproduces
the previously used central-charge-one fermion stress weights, including
the missing first harmonic.

For equal delays in the reflecting limit, the exterior singlet has zero
throat support and therefore `A=B=0`. Its presence changes total spectral
counts but not (5). In this special case the throat spectral response is
exactly `q-1` free-circle responses. With unequal delays, (5) with the
coupled eigenfunctions is required. The check computes this observable
from independent finite-element mass and stiffness forms and verifies
mesh convergence, the equal-delay reduction and failure of the naive
global factor `(q-1)/q` for a specified unequal-delay example. It does
not yet take the large-`q`, shrinking-`delta` response limit.

Concretely, take `q=3`, `T=1.3`, `D=(0.7,1.4,2.8)` and integrate
`rho_V/pi` against `w(Omega)=exp[-x^2/(1-x^2)]` for `|x|<1`, zero
otherwise, where `x=(Omega-4.8)/1.8`. At 96 cells per edge, the
reflecting result is `0.0428331083`; the naive `(q-1)/q` times free
result is `0.0324782303`, lower by 24.17% of the coupled value. The
48-to-96-cell change in the coupled value is 0.65%; successive mesh
changes decrease approximately quadratically. These are convergence
diagnostics, not certified error bars. This example is a finite graph,
not a physical large-flux wormhole parameter point. For the equal-delay
check, the finest mesh differs from the analytic stress value by 0.0111%.

## 6. Consequence for the damping investigation

The [independent-channel model](collective-clock-response.md) has an
absorptive scale `q delta`, `delta=d/ell`, in dimensionless throat units.
If frequency resolution is `eta ell~delta` and `D_min~d`, then
`eta D_min~delta^2`. Bound (1) is `O(delta^(-4))`. Relative to the
reference scale `q delta`, its **bounded boundary-readout** bound is
`O(1/(q delta^5))`. For fixed `g` and `d/r_e=q^nu`, this tends to zero
when `nu>4/5`. At `nu=15/16` the exponent is `-11/16`. This is compatible
with the previously identified scale window; it is not a response error
estimate, since the clock couples to quadratic stress throughout the
throat and also changes the interaction term.

In particular, the low-frequency mass `M ell~g sqrt(q)` does not by
itself delete long-return oscillators or close the absorptive channel.
It also does not establish their coupling to the clock. At a free
resonance the finite-rank difference can diverge as `eta` decreases;
the check retains an explicit unequal-length example. Neither replacing
`q` by `q-1` nor invoking a rank bound at zero resolution justifies the
previous damping coefficient in the interacting graph.

Equation (5) now supplies the exact oscillator spectral formula for a
defined clock family, including the singlet mass vertex, and identifies
its contact term. The decisive open question is its continuum asymptotic
at damping-scale resolution, with control of the real renormalized
response. The physical extension then needs actual throat/exterior
matching, exterior gauge interactions, holonomy/charge-sector and
finite-time state control. A successful graph calculation is necessary
evidence, not sufficient four-dimensional closure. No new novelty,
traversability or shortcut conclusion follows.

## 7. Reproduction and scope

Run `python checks/schwinger_return_feedback.py` with NumPy, SciPy and
SymPy. The receipt records script/note hashes and package versions.
Checks compare a direct scalar boundary-value ODE with the analytic
constant-mass load; verify its energy identity, passivity and reflection
bound for varying mass; compare the full return inverse with its
rank-two reduction; and check (3) against full determinants. Independent
piecewise-linear finite-element forms verify the splitting and min-max
bounds and refine the equal-length analytic spectrum. Sine/cosine
oscillators reproduce the known stress weights independently; convergent
form calculations evaluate the coupled finite-band pair response and
test its equal-delay limit. These are floating
point supplementary checks, not interval proofs or formal certification.
No Lean is invoked; all later Lean development remains GCP-only.

The first complex-valued boundary-value solver run exhausted its mesh;
the same ODE rewritten as a real system passed at the original tolerance,
with the constant-mass analytic comparator and tighter-tolerance checks.
The first stress run missed the declared 2.5% successive-mesh threshold
(2.64% at 48 cells); adding the 96-cell mesh resolved it without relaxing
that threshold. Both failed source versions remain in the external
research cache, and their terminal outcomes are recorded in delivery
verification. Diagnostic runs are not substituted for acceptance.

Remaining to-do list: control the coupled stress response at damping-scale
resolution and its renormalized real part; close physical matching,
gauge-interaction and finite-time state
bounds; finish primary/fallback novelty selection and the consequential
result; independently verify, GCP-certify suitable endpoints and prepare
the reproducible manuscript.
