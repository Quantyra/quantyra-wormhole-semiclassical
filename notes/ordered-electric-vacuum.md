# Ordered electric fluctuations and excitation above a dressed vacuum

S072, 2026-10-13. **The derivative-moment condition of the preceding
quantum-source bound fails for the unregulated Gaussian electric vacuum.
Keeping its ordered spectrum instead gives finite quadratic excitation
bounds, with no electric-mode cutoff.** A smooth-ramp version subtracts
the perturbative endpoint vacuum dressing before measuring excitation.
The calculation is on the fixed global product and does not yet identify
the actual driven junction or its initial state.

## The constrained electric field on the same global product

Retain the geometry, massive fermion vacuum and sphere tower of
[curved-massive-response.md](curved-massive-response.md). Set

```
C=b^2/a^2,  lambda=q/pi,  kappa=4pi b^2/g^2,
s=m_gamma^2=lambda a^2/kappa=q g^2/(4pi^2 C)>0.
```

The electric center-of-mass boson is the one already retained in the
projected LLL action, not an additional oscillator. Its origin is
[MMP appendix A](https://arxiv.org/html/1807.04726v3); the constrained
normalization is derived in [phase-influence-action.md](phase-influence-action.md).
The bosonized global Hamiltonian and vanishing boundary electric field
are also discussed in [Bharadwaj, Isen and Kang, section 3.3](https://arxiv.org/html/2607.19468v1).
Only that action/boundary comparison is imported, not their screening or
lattice results.

Writing electric displacement as `p=sqrt(lambda) phi`, elimination of
the radial constraints gives

```
H_e=(1/2) integral [Pi_phi^2+(partial_sigma phi)^2
                             +s sec^2(sigma) phi^2] d sigma,
partial_sigma A_t = (s/sqrt(lambda)) sec^2(sigma) phi       (1)
```

in radial gauge, up to an irrelevant overall sign convention. A spatial
constant in `A_t` has zero pair matrix element. Choose the Friedrichs,
fast-falloff boundary realization and zero background electric flux.
This is an explicit boundary choice; at small `s` it is not the only
scalar self-adjoint realization. It has not been matched to the initial
RN boundary-transfer problem or to a compact charge rotor.

The normalized scalar modes are

```
Delta=(1+sqrt(1+4s))/2,       omega_j=j+Delta, j>=0,
v_j(sigma)=h_j^(-1/2) cos^Delta(sigma) C_j^Delta(sin sigma),
h_j=pi 2^(1-2Delta) Gamma(j+2Delta)
                  /[j! (j+Delta) Gamma(Delta)^2].          (2)
```

Substitution in `-v''+s sec^2(sigma)v=omega_j^2 v` and Gegenbauer
orthogonality prove (2). Hellmann--Feynman, or a direct weighted integral,
then gives the especially useful identity

```
integral sec^2(sigma) |v_j|^2 d sigma=2omega_j/(2Delta-1).  (3)
```

Let `A_j'=(s/sqrt(lambda))sec^2(sigma)v_j`. In the interaction picture,
the electric potential is

```
A(t,sigma)=sum_j A_j(sigma)/sqrt(2omega_j)
                    [a_j exp(-iomega_j t)+a_j^dagger exp(iomega_j t)].
Q[A_j]/(2omega_j)=K,
K=s^2/[lambda(2Delta-1)],
Q[V]=integral cos^2(sigma)|V'|^2 d sigma.                 (4)
```

Every scalar mode contributes exactly the same positive Maxwell
seminorm. For `N` retained modes the vacuum equal-time squared norm of
`partial_sigma A` is `K N`; that of `partial_t^2 partial_sigma A` is
`K sum_(j<N) omega_j^4 ~ K N^5/5`. Thus both raw norms diverge for any
fixed nonzero coupling. Multiplying by a smooth nonzero pulse does not
remove the divergence of the time-integrated norms in the preceding
theorem. This is a failed hypothesis of that sufficient estimate, not
an infinite physical excitation prediction. Normal ordering a positive
Hilbert-space norm to zero would not repair the argument.

## Ordered transitions give a different, finite estimate

Introduce a real diagnostic envelope `f(t)` multiplying the reduced
cross-sector interaction `V=int A rho_H d sigma`, where `rho_H` is the
massive fermion density. Initially take the source Gaussian vacuum and
the massive vacuum, factorized at zero cross coupling. This pulse is a
specified reduced-Hamiltonian test, not an identification of the
physical Maxwell interaction with a switchable charge.

Let `E` be a massive pair energy and define its positive measure

```
dmu_j(E)=sum_pairs E |<pair|int A_j rho_H|0>|^2/(2omega_j),
mu_j((0,L]) <= D K L^2,
D=B C(1+1/q),  B=9.056041... .                           (5)
```

Equation (5) follows from the preceding curved current bound and (4).
In the vacuum, the source must gain one quantum when a pair is created.
Its ordered amplitude therefore contains `f_tilde(E+omega_j)`, giving

```
W_H^(2)=sum_j integral |f_tilde(E+omega_j)|^2 dmu_j(E),
W_total^(2)=sum_j integral (E+omega_j)/E
                            |f_tilde(E+omega_j)|^2 dmu_j(E). (6)
```

The second expression includes source energy as well as massive energy.
A symmetrized source spectrum would introduce vacuum absorption at
`E-omega_j` and is not (6).

For total transition frequency `U`, define cumulative measures `N_H(U)`
and `N_total(U)` by restricting (6)'s weights to `E+omega_j<=U` before
inserting the Fourier transform. Since `Delta>=1`,

```
N_H(U) <= D K sum_j (U-omega_j)_+^2 <= (D K/3)U^3.       (7)
```

The last sum is bounded by the integral of `(U-x)_+^2` on `[0,infinity)`:
each value at `j+Delta` is at most the integral over `[j,j+1]`.
For the unweighted pair probability, Stieltjes integration of (5) gives
`P_j(R)=int_0^R dmu_j(E)/E<=2D K R`. Hence, with `R=U-omega_j`,

```
N_total(U) <= sum_(omega_j<U) [D K R^2+2D K omega_j R]
           <= D K integral_0^U (U^2-x^2) dx
            = (2D K/3)U^3.                              (8)
```

These are full source, sphere and radial sums. If
`F0=||f||_1`, `F2=||f''||_1` and temporal boundary terms vanish, the
envelope `min(F0^2,F2^2/U^4)` and (7)--(8) prove

```
W_H^(2)     <= (4/3) D K F0^(1/2) F2^(3/2),
W_total^(2) <= (8/3) D K F0^(1/2) F2^(3/2).              (9)
```

For example, for any positive measure with cumulative bound `a U^3`,
integration against this envelope is at most
`4a F2^2/L0`, with `L0=sqrt(F2/F0)`. Positive cutoff convergence
therefore defines the quadratic work although the norms in (4) diverge.
No zero-point energy is inserted into a classical velocity stress.

## Thermal source occupations

For a diagonal Gaussian source with occupations `n_j`, still with the
massive bath in its vacuum, the ordered expression is

```
W_H^(2)=sum_j integral [(n_j+1)|f_tilde(E+omega_j)|^2
                             +n_j|f_tilde(E-omega_j)|^2] dmu_j(E). (10)
```

Define `S_r=sum_j n_j omega_j^r`, `r=0,1,2`. The thermal contribution
above the vacuum answer has the bound

```
Delta W_H^(2) <= D K [4F0 F2 S_0
                +(8/3)F0^(3/2)F2^(1/2) S_1+F0^2 S_2].   (11)
```

To see this, shift the Fourier envelope by `omega`. In integration by
parts, discard the nonpositive contribution from its increasing side.
The decreasing side gives
`4D K F2^2 int_(omega+L0)^infinity E^2/(E-omega)^5 dE`, equal to
`D K[2F0F2+(8/3)omega F0^(3/2)F2^(1/2)+omega^2 F0^2]`.
The additional stimulated-emission part is at most `2D K F0F2`.
Summing with `n_j` proves (11).

For `n_j=(exp(beta omega_j)-1)^(-1)` all three sums converge. They may
be evaluated independently by expanding `n_j=sum_(l>=1)exp(-beta l omega_j)`
and performing the geometric sums over `j`. This covers that stationary
thermal Gaussian source, not a driven or squeezed source, the two-sided
entangled preparation, or a thermal massive bath. Equation (11) bounds
massive energy only; a populated source can lose energy.

## Excitation above the endpoint dressed vacuum

A permanent interaction cannot be treated as a compact pulse from a
factorized vacuum. At finite regulator, take instead
`H(t)=H_0+epsilon h(t)V`, initially in the ground state at `h_i`, and
compare final energy with the ground energy at `h_f`. Let the free
excited channels have energy `U_alpha` and matrix element `v_alpha`.
First-order stationary perturbation theory gives the endpoint ground
correction `-epsilon h v_alpha/U_alpha`. Subtracting that correction
from the evolved state, an integration by parts gives the nonadiabatic
amplitude, up to a phase,

```
c_alpha,ex = epsilon v_alpha/U_alpha
                      integral h'(t) exp(iU_alpha t) dt + O(epsilon^2).
```

Thus the coefficient of excess *total* final energy is

```
W_ex^(2)=integral |h'_tilde(U)|^2/U^2 dN_total(U).         (12)
```

The static second-order vacuum shift is subtracted with the same
regulator at both endpoints. A constant interaction with its own
initial dressed ground state has zero excitation in (12).

Set `G0=||h'||_1`, `G2=||h'''||_1`; assume `h'` and `h''` vanish
outside the ramp. Integrating
`min(G0^2/U^2,G2^2/U^6)` against (8) yields

```
W_ex^(2) <= (8/3)D K G0^(3/2)G2^(1/2).
W_ex^(2) <= (4/3)D K G2^2/U_min^3,
U_min=2M_1+1+Delta.                                     (13)
```

For the first result, split at `L0=sqrt(G2/G0)`; the two derivative
integrals each contribute `2a G0^(3/2)G2^(1/2)`, with
`a=2D K/3`. For the second, use the high-frequency envelope throughout
and the exact lowest joint energy `U_min`. Dropping the nonpositive
lower-end term in the Stieltjes integral is conservative.

A unit quintic ramp `h=10x^3-15x^4+6x^5`, `x=t/T` on `[0,1]`, has
`G0=1`, `G2=40/(sqrt(3)T^2)`. The first bound falls as `1/T`; the
gap-based one falls as `1/T^4`. This is a uniform bound on the quadratic
excitation coefficient after endpoint dressing. It does not construct
a continuum interacting ground state or control a uniform higher-order
remainder. The physical junction changes a boundary interaction, not
necessarily the cross-sector coupling `h`; that transfer is still open.

## Quantitative comparison and verification

At `k=2000`, using the existing reference family and `T=0.25`,

```
s=3.22303287e-26, K=2.03966742e-64, U_min=38.0962606,
W_total,pulse^(2) <= 1.66602802e-50,
W_ex,quintic^(2)  <= 2.98733324e-52.
```

The pulse is `[1-(t/T)^2]^3` on `[-T,T]`; the ramp is the quintic
above on `[0,T]`. These are dimensionless global Killing energies for
the specified envelopes, not bounds on the actual opening error. The
Gaussian electric coefficient and its boundary realization are inputs.
These small values do not match them to the physical preparation.

[Checker](../checks/ordered_electric_vacuum.py),
[receipt](../evidence/ordered-electric-vacuum.json), and
[source/dependency record](../evidence/ordered-electric-vacuum-source-review.json).
The checker verifies the scalar equation, normalization and (3) at four
weights. Complete finite joint-energy sums use the massive Dirac
commutator to evaluate each matrix element directly, with independent
spatial refinement. A separate analytic tail bounds every omitted
source and pair energy, for both the pulse and the dressed ramp.
For a pulse cutoff `L`, its exact transform gives
`|f_tilde(U)|<=D0/(T^3 U^4)`,
`D0=96[1+6/(LT)+15/(LT)^2+15/(LT)^3]`. The omitted total work is at most
`16D K D0^2/(15T^6 L^5)`. The ramp tail is at most
`4D K G2^2/(3L^3)`.

Independent thermal sums and shifted-envelope quadratures check (11).
Three exact two-level evolutions start in the initial dressed ground
state and project onto final dressed excitations. Their excess energies
approach (12) with the expected quadratic-order refinement. They do not
stand in for a wormhole evolution.

The explicit derivative-moment failure selects the ordered spectral
method as the next route. Still required are the driven neutral/compact
source and boundary conditions, initial dressing in the physical theory,
evolving geometry, the other instantaneous massive-charge interactions
generated by Gauss elimination, common counterterms and complete stress.
The current result bounds one quadratic cross-sector channel. It does
not omit those interactions by calling this the whole massive correction.
No physical creation, full stability or novelty acceptance is claimed.

Remaining to-do list: match the driven junction and initial dressed state
to an ordered spectral estimate; include the remaining constrained
interactions, evolving geometry and complete stress; establish preparation,
sufficient stability and finite-energy transport.
