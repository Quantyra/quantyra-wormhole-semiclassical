# Boundary angle and charged-vertex normalization

S072, 2026-10-13. **The weak-branch boundary parameter fixes the free
junction angle, and the initial electric Gaussian midpoint vertex admits
a finite, bounded normalization relative to the charged strip reference.**
For the selected family the matching exponent is below `1/(48 k)`, or
`1.042e-5` at `k=2000`, throughout the new capture preparation interval.
This is a projected-model normalization result. It does not identify the
complete interacting four-dimensional state or bound its switching stress.

## Boundary parameter, contact convention and clock

[Bintanja, Freivogel and Rolph, sections 2.1.1--2.1.3](https://arxiv.org/html/2308.00871v2)
specify a dimensionless boundary parameter `h_b`, normalized spherical
spinors, boundary conditions (2.12), and the frequency shift (2.13)--(2.14).
On `0<=h_b<1`, their `lambda_b(h_b)=atan(h_b)`. Comparison with the
positive towers `1/2+theta/pi+2n` and `3/2-theta/pi+2n` in the
[returned free junction](returned-switch-capture.md) fixes

```
theta = 2 atan(h_b),                 h_b = tan(theta/2).       (1)
```

This branch relation is already present in that source. It is not a new
wormhole construction. We use its boundary conditions and spectrum,
not the apparent extra imaginary factors in HTML equation (2.16), and
not its energy prescription in place of our inherited point-split stress.

An equivalent unfolded free contact makes the coefficient convention
explicit. Put `J=[[0,-1],[1,0]]`, `Psi=(psi_L,psi_R)`, and
`X=i psi_L^dagger psi_R-i psi_R^dagger psi_L=Psi^dagger(-iJ)Psi`.
For `H_contact=-gamma X`, the one-particle delta contact is `M=i gamma J`.
Define it by evaluating the field at the average of its two limits.
Integration of the first-order field equation then gives

```
(1+iM/2) Psi_+ = (1-iM/2) Psi_-,
S = (1-h_b J)^(-1)(1+h_b J) = exp(J theta),
gamma = 2 h_b = 2 tan(theta/2) = theta+theta^3/12+O(theta^5). (2)
```

This defines the midpoint contact convention, rather than deriving an
unregulated interacting Hamiltonian from its free scattering matrix.
The `h` multiplying the normal-ordered neutral vertex in
[neutral-transfer-current.md](neutral-transfer-current.md) includes a
vertex matching factor and has not thereby become `theta` or `gamma`.

For a common boundary clock `t=t(U)` with positive derivative, canonically
normalized massless endpoint fields obey
`psi_U=sqrt(dt/dU) psi_t`. The two factors in the transfer bilinear cancel
`dU` in the classical boundary action. Thus `h_b` has no extra common-clock
power. Unequal endpoint clocks would instead require both square roots.
This classical statement does not omit the quantum stress transformation
or justify using a stationary bath on an evolving geometry.

The free cancellation protocol can consequently be imposed as
`h_b(U)=tan[theta(t(U))/2]`, where `theta` is the shaped angle from the
previous note. Its final values are `theta_f=4pi/k` and
`h_b,f=tan(2pi/k)`. At `k=2000`, `gamma_f/theta_f-1` is about `3.29e-6`.
That small difference still matters to an exact Fourier-zero assertion:
the zeros of `theta'` are not automatically zeros of `gamma'` or of the
matched interacting source. A prescribed laboratory clock must implement
the actual trajectory `t(U)`; assuming the free trajectory in the
interacting problem leaves another error to bound.

## Both boundaries and one regulator

Use the single-boundary projected electric action already derived in
[phase-influence-action.md](phase-influence-action.md):

```
lambda = q/pi,       K(s) = electric Dirichlet-to-Neumann coefficient,
I(s) = lambda/K(s),  I0=I(0),
G_single(s) = 1/[s^2 I(s)] = K(s)/(lambda s^2).              (3)
```

Here there is one massless fermion species, and `K` comes from the
half-line operator `-d_x^2+V_e(x)+s^2` on the initial exterior, with
`V_e>=0` and regular horizon behavior. `I0` is the screened static
susceptibility, not the bare Maxwell coefficient. For two equal copies,
the relative phase is `alpha=chi_R-chi_L`. Its action coefficient is
`I_rel=I/2`, and its covariance is twice (3). With the constant initial
clock `t=kappa U`, `kappa=2pi/beta`, its global-time action coefficient
is `kappa I(kappa |omega|)/2`. This relative factor was already retained
in [collective-electric-response.md](collective-electric-response.md).

For the noncompact Gaussian thermofield comparison define
`W_full=-log <exp(i chi_R)exp(-i chi_L)>`. The midpoint covariance gives

```
W_full(a) = (1/beta) sum_(n!=0) [1-(-1)^n] G_single(s_n) e^(-a s_n)
          = 4/(beta lambda) sum_(n>0 odd) K(s_n)/s_n^2 e^(-a s_n),
s_n=2pi |n|/beta.                                          (4)
```

The zero-frequency constant cancels. For a massless noncompact coordinate,
the thermal trace is defined by first regulating its volume and cancelling
that common volume in the relative midpoint quantity. Formula (4)
includes the resulting free-particle contribution `beta/(8 I0)`.
It is not a sum over the actual compact electric/spin representations.

The reference factor removed from the free strip was
`C_ref=Z_ref^(-1) exp(i phi_c)`, with

```
W_ref(epsilon)=-log Z_ref
    = (2/q) sum_(n>0 odd) e^(-epsilon n)/n
    = (1/q) log coth(epsilon/2).                            (5)
```

Choose the same physical regulator: `epsilon=kappa a`. Setting `K(s)=s`
in (4) gives exactly (5), including the coefficient. Treating only one
boundary as the relative vertex would leave half the required logarithm.
The subtraction below uses the whole reference factor, rather than an
unrenormalized unitary vertex with an unchanged coefficient.

## A finite matching bound

The scalar variational problem gives, for `s>0`,

```
0 <= K(s)-s <= K(0).                                      (6)
```

For completeness, fixing the boundary value to one, dropping `V_e>=0`
from the energy infimum gives `K(s)>=s`. Let `z0` minimize the
zero-frequency form. The scalar maximum principle gives `0<=z0<=1`.
The trial field `exp(sx)z0(x)` has energy
`K(0)+2s^2 integral exp(2sx)z0^2 dx <= K(0)+s`, by integration by
parts using the equation for `z0`. This proves the upper bound. It is
the electric specialization of the existing
[causal-response proof](causal-compact-response.md). At zero potential
the bound follows by its continuous limit. No reflecting horizon or
spectral gap is inserted.

The difference of (4) and (5) therefore has the absolutely convergent limit

```
Delta W = 4/(beta lambda) sum_(n>0 odd) [K(s_n)-s_n]/s_n^2,
0 <= Delta W <= beta K(0)/(8 lambda) = beta/(8 I0),
Z_e,matched = lim_(a->0) exp[-W_full(a)+W_ref(kappa a)]
            = exp(-Delta W),
exp[-beta/(8 I0)] <= Z_e,matched <= 1.                     (7)
```

The sum of positive odd inverse squares is `pi^2/8`, proving the bound.
Dominated convergence also permits any shared bounded cutoff that tends
pointwise to one. If `V_e<=Vmax`, the trial field `exp(sx)` further gives
`K(s)-s<=Vmax/(2s)`. For an odd cutoff index `M>=1`, a simple conservative
tail bound, summing over all omitted integers, is

```
0 <= Delta W-Delta W_M
   <= min[beta K(0)/(pi^2 lambda M),
          beta^2 Vmax/(8 pi^3 lambda M^2)].               (8)
```

Both separate variances diverge when the cutoff is removed; (7) concerns
their ratio in the specified common normalization. No finite microscopic
four-dimensional operator renormalization is fixed by this Gaussian
effective theory alone.

The oscillator-only shift is
`Delta W_osc=Delta W-beta/(8 I0)`, in `[-beta/(8 I0),0]`.
Its normal-ordered factor can exceed one; it is not a probability or the
expectation of an unrenormalized unitary. Adding an exact compact rotor
later must replace the Gaussian free-particle factor, not multiply a
second copy into (7). The charge-spin quotient and neutral/angular
branching still determine the actual compact weights.

## Candidate size and what it controls

Retain `eta=1/k`, `g=k^(-4)`, `L/l_p=k^5`, `q=k^4`, `N=2q`,
`b^2=L^2 X`, `C=1+6X`, and `X(1+3X)=pi k^6`. The actual new capture
interval is `1/4<=z0<=7/25`, rather than the older `z0=1/2` comparison.
With `phi_r=b^3/(G C)` and `beta=2pi phi_r/(N eta z0)`, the static
screening proof gives `I0>=I_bare=4pi b/g^2`. Direct substitution yields

```
x_bare = beta/(2 I_bare) = X/[8 k z0(1+6X)] < 1/(48 k z0),
Delta W <= x_bare/4 < 1/(192 k z0) <= 1/(48 k).            (9)
```

Thus the initial Gaussian matched amplitude loses less than `1.042e-5`
at `k=2000`. This is neither a later-time dephasing bound nor an energy
error of that size. The older compact free-evolution numerical bound was
for a shorter physical duration and must not be transplanted to the
longer returned switch. Likewise (9) does not justify replacing an
operator-valued interacting source by its expectation.

## Independent checks and next step

[The checker](../checks/boundary_vertex_matching.py) verifies the Cayley
map and static frequency branches, and checks its distinction from a
linear contact profile on the actual longer waveform. It compares the
matching sum with an independent continuous spectral integral for
`V=m^2`. For that comparator `K(s)=sqrt(s^2+m^2)` and
`d sigma/dw=2 sqrt(w^2-m^2)/(pi lambda w)` for `w>m`; the zero-mode
weight is `m/lambda`. The integral of the physical-minus-reference
oscillator measures converges without a cutoff and reproduces (7).

A finite positive barrier checks the horizon-asymptotic class:
`V=m^2` on `[-ell,0]` and zero to its left. Its exact coefficient is
`K(s)=k_s(s+k_s tanh(k_s ell))/(k_s+s tanh(k_s ell))`,
`k_s=sqrt(s^2+m^2)`. Direct integration of its Riccati equation checks
this expression and (6). These comparators verify the normalization and
spectral identities; they are not numerical substitutes for the analytic
RN bound. The receipt records regulator convergence, positive tail
bounds and candidate substitutions. No local Lean/Lake run is used.

The next problem is to insert the matched charge vertex into a joint
initial state and a driven current calculation with one regulator and
one subtraction. The exact bath-current identity retains mixed quantum
correlations. The small midpoint exponent cannot discard them.

Remaining to-do list: determine the compact/angular and microscopic
vertex matching; control the interacting driven state and complete
renormalized stress on the evolving geometry below the capture and area
margins; establish preparation, quantum fluctuations, sufficient general
stability and finite-energy transport.
