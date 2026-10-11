# Return phases and lowest Landau level support

2026-10-10. S051/S052 feasibility decision. A unitary return of the lowest
Landau level (LLL) channels need not preserve even the angularly averaged
negative stress in the free-channel approximation. A common electromagnetic
Wilson phase cannot remove relative eigenphases. The exact phase criterion
below identifies the quantity an exterior calculation must control.
The long-cycle suppression mechanism is already discussed by MMP; neither
that mechanism nor the twisted-circle determinant is claimed new.

## Model and observable

Take `q >= 1` free complex massless two-dimensional Dirac channels with a
common optical circumference `L`, zero temperature, and the self-adjoint
vector boundary condition

```text
psi(x+L) = -exp(i phi) U psi(x),        U in U(q).
```

`U` is fixed, energy independent and the same for both components of the
Dirac field. The minus sign fixes the untwisted antiperiodic convention.
There is no reflection, escape, higher-level coupling or interaction.
The electromagnetic variable `phi` shifts all eigenphases together.
Allowing it to minimize the fermionic energy is the most favorable common
phase choice; an actual gauge/geometry action may also contribute.

For the support comparison, additionally take `L = pi ell`: the exterior
optical length is negligible relative to the throat. Compare to `q`
identical reflecting strips on the same conformally related throat, as in
[MMP section 5.2](https://arxiv.org/html/1807.04726v3#S5.SS2).
The local conformal anomaly cancels between states. The reference itself
has zero null stress by its AdS symmetry. This reference differs from the
holonomy average used in the [higher-mode winding bound](warped-winding-bound.md).

An arbitrary `U` can break spherical symmetry. The channel trace below is
the angular integral of the LLL stress, or its spherical average after
division by the area. It is not the pointwise four-dimensional stress.
Nonspherical components cannot be discarded in a physical solution.

## Circle energy and common phase

Write the eigenvalues of `U` as `exp(i theta_j)` and set

```text
a_j(phi) = wrap_to_[-pi,pi](theta_j + phi),
D2(phi) = sum_j a_j(phi)^2.
```

The spatial momenta in eigenchannel `j` are
`2 pi (n + 1/2 + a_j/(2 pi))/L`. For one complex Dirac field the vacuum
energy is minus the sum of absolute momenta. Subtracting the corresponding
infinite-line integral with the same exponential energy regulator gives
`2 pi B2(1/2 + a_j/(2 pi))/L`, where `B2(z)=z^2-z+1/6`.
Equivalently, the two geometric series over nonnegative momenta give this
finite term as the regulator tends to zero. Hence

```text
E_U(phi) = -q pi/(6 L) + D2(phi)/(2 pi L).                 (1)
```

The winding representation gives the same result:

```text
E_U(phi) = (2/(pi L)) sum_{k>=1} (-1)^k
                            Re[exp(i k phi) Tr U^k]/k^2.
```

It also follows from `-(2/pi) integral_0^infty
Re log det(I + exp(-xi L + i phi) U) dxi`. The Fourier series is absolutely
convergent. A change of channel basis conjugates `U` and leaves the energy
unchanged. Off-diagonal matrix entries alone are not a measure of harmful
mixing: the conjugacy class, specifically its phase spread, matters here.

Common-phase minimization is exactly circular least squares. Break the
phase circle at the at most `q` wrap discontinuities. On each remaining
interval `D2` is a quadratic with its vertex at minus the mean of the
unwrapped phases. Comparing these vertices and interval endpoints obtains
the global minimum without a local-optimization assumption.

## Reference subtraction and opening coefficient

The flat reflecting-strip density is `-q pi/(24 L^2)`. The state difference
in the conformal coordinates has equal time and space diagonal entries:

```text
Delta T2_tt = Delta T2_xx = E_U/L + q pi/(24 L^2)
                         = -q pi chi/(8 L^2),
chi(U,phi) = 1 - 4 D2(phi)/(q pi^2).                     (2)
```

Thus negative averaged support occurs precisely when
`D2(phi) < q pi^2/4`. The optimized support uses `min_phi D2(phi)`.
The [reproduced spherical first-order equation](casimir-reference-matching.md)
has `c''(rho)=2 chi/(1+rho^2)^2`: invariant reference pressures do not supply
this opening. For nonscalar `U`, equation (2) specifies only the monopole
projection of the source about that background. It does not solve the
nonspherical equations or prove absence of every possible wormhole.

## Cyclic return and rigid rotation

Let `U` be the permutation sending each channel to the next in a single
`q`-cycle. Its eigenvalues are the `q`th roots of unity, so `Tr U^k=0`
unless `q` divides `k`. Following a field through all channel segments
unfolds the spectrum into one complex Dirac field on a circle of length
`q L`. The common phase can make that long circle antiperiodic. Therefore

```text
min_phi E_U = -pi/(6 q L),
min_phi D2  = pi^2 (q^2-1)/(3q),
chi_best    = -1/3 + 4/(3q^2).                          (3)
```

The first equality proves global optimality; the second also follows by
summing the squares of the centered regular-polygon phases
`2 pi (j-(q-1)/2)/q`. Compared with `q` independent antiperiodic circles,
the magnitude of the circle energy is smaller by `1/q^2`. The best support
coefficient is `1` at `q=1`, zero at `q=2`, and negative at every `q>=3`.
This refutes arbitrary-unitary preservation of free LLL support, even when
the common phase is freely optimized. It is not a physical MMP no-go.

The same eigenphases occur within the actual LLL rotation representation:
MMP appendix C gives spin `j=(q-1)/2`. A prescribed rigid return rotation
`exp(i delta J_z)` with `delta=2 pi/q` has the centered roots in (3), up to
a common phase. This is a test of an allowed rotation representation, not
a derivation of the MMP exterior or its dynamical gauge fields.

It also shows a limitation of a geometric displacement criterion alone.
The displacement `2 pi/q` is asymptotically smaller than the orbital width
scale `q^(-1/2)`, while (3) tends to `-1/3`. For an equatorial spin coherent
state, the squared return overlap is `cos(delta/2)^(2(q-1))`, which tends
to one: the state is the symmetric product of `q-1` equatorial spin halves,
each with overlap `cos(delta/2)`. Small displacement or large wavepacket overlap does not uniformly
control the many-channel vacuum energy. A state-dependent phase can be
nearly invisible to this overlap magnitude.

## Sufficient phase control

Useful quantitative conditions follow directly from (2). If some common
phase makes every `|a_j| <= eta`, then

```text
chi >= 1 - 4 eta^2/pi^2.                               (4)
```

In particular `eta < pi/2` gives strictly negative averaged support. Since
`|a| <= (pi/2)|exp(i a)-1|` for `|a|<=pi`, there is also a sufficient
matrix bound, with the ordinary Frobenius norm:

```text
chi >= 1 - ||exp(i phi) U - I||_F^2/q.                 (5)
```

For a rigid spin rotation, when `|delta|(q-1)/2 <= pi`, taking `phi=0`
gives the exact, unoptimized value

```text
chi(0) = 1 - delta^2 (q^2-1)/(3 pi^2).                 (6)
```

Optimizing the phase can only increase it. In particular `delta=c/q`
gives a uniform positive lower bound for `|c|<sqrt(3) pi`; this range
also satisfies the no-wrap condition. A small coefficient `c` preserves
almost all the free-channel support. The relevant scale is `q delta`,
rather than displacement divided by orbital width.

## Prior work and the physical decision

[MMP section 6.1.1 and appendix G](https://arxiv.org/html/1807.04726v3#S6.SS1.SSS1)
already identify rotation-induced mixing, enlarged return cycles and weaker
Casimir support. They estimate a screened displacement of order `Omega d/q`.
If the return were a rigid rotation with this scaling and a bounded small
`Omega d`, (6) would be consistent with preserved support. That inference
does not determine the actual operator, coefficients, delays or phase
variation. The argument here therefore neither contradicts their stated
slow-rotation regime nor fills its quantitative operator estimate.

[Cherman, Sen, Wagman and Yaffe](https://arxiv.org/html/1612.00403v1)
derive holonomy-trace winding cancellations and enlarged-circle finite-size
behavior in sections I-II. Their adjoint pion theory is a different field
content and dimension. [Cherman, Sen, Unsal, Wagman and Yaffe](https://arxiv.org/html/1706.05385v1)
use the fundamental roots-of-unity twists in equations (3)-(4) and derive
the fermion winding cancellation in (9). These establish the methodological
precedent; their QCD interactions or four-dimensional coefficients are not
imported into (1). Edition and dependency review remains incomplete.

Decision: classify (1)-(6) as a reproducible diagnostic application of known
ingredients. Do not select the long-cycle effect or phase formula as a
standalone novel-paper or Lean target. Drop arbitrary-unitary robustness as
a possible physical bridge for the previous winding bound. Retain a
quantitative exterior-return calculation as the next feasibility test:

1. Derive the projected exterior return operator, including geometric and
   electromagnetic phases, or prove a controlled operator estimate. Begin
   with the static axisymmetric case before its screened slow rotation.
2. Separate spectral phase spread from frequency-dependent delay, unequal
   optical lengths, reflection, leakage and higher-mode transitions. A small
   classical drift estimate alone cannot supply all of these bounds.
3. Account for the interacting singlet and geometric gauge sectors. MMP
   appendix A already discusses their mass generation; freely minimizing a
   prescribed common phase does not integrate those dynamics out.
4. Establish whether the resulting estimate has a consequential improvement
   over MMP's existing regime. If it only reproduces that regime, retain it
   as baseline and reassess the candidate before substantial formalization.

## Verification

[The executable check](../checks/lll_return_holonomy.py) and
[receipt](../evidence/lll-return-holonomy.json) verify seven exact identities,
the global cyclic phase minima, and the rigid-rotation and sufficient-bound
examples. Direct matrix determinants test basis invariance. Independently,
the regulated Hamiltonian spectrum is summed before subtraction of the
same-regulator infinite-line integral. Two successive extrapolations agree
with (1); their residuals and cutoffs are retained. These floating-point
checks do not certify all numerical errors or replace the ordinary proof.
No Lean invocation, physical stability result or novelty acceptance is
claimed. Prior checks and manuscripts are preserved.

Remaining to-do list: derive and bound the physical exterior return, close
the reference/interacting-state comparison, complete the contribution
review, and select a consequential stability, transport or obstruction target.
