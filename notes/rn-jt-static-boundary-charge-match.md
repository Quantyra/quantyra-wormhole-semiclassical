# RN-to-JT boundary charge matching at leading order

2026-10-12. This derives the leading relation between each asymptotically
flat RN mouth's mass parameter and the global throat clock rate in the MMP
overlap. It reproduces the static mass/length matching and identifies the
boundary charge that a dynamical calculation must source. The differentiated
relation is only an adiabatic tangent; it is not a finite-frequency
constitutive law.

## Static overlap matching

At fixed magnetic charge, write the RN mouth metric function as

```text
F(r)=1-2GM/r+r_e^2/r^2
    =(1-r_e/r)^2-2G deltaM/r,
deltaM=M-r_e/G.
```

In the MMP overlap define
`rho=ell(r-r_e)/r_e^2`, `tau=t/ell`, with
`r_e/ell << 1` and `1<<|rho|<<ell/r_e`. Then

```text
F(r)=r_e^2 rho^2/ell^2 - 2G deltaM/r_e + higher overlap orders.
```

The global throat metric has
`g_tt=-r_e^2(1+rho^2)/ell^2`. Matching the constant and `rho^2` terms
therefore fixes

```text
deltaM=-r_e^3/(2G ell^2)=-A x_*^2/2,
A=r_e^3/G,  x_*=d tau/dt=1/ell.                            (1)
```

The same overlap gives
`Phi=R^2/r_e^2-1 ~= 2r_e rho/ell`. MMP's static throat solution
`Phi=alpha(1+rho atan(rho))` has large-`rho` coefficient
`pi alpha/2`, so `alpha=4r_e/(pi ell)`; their equation (5.40) then gives
`ell=16r_e^3/(qG)`. Substitution in (1) gives
`2 deltaM=-q^2G/(256 r_e^3)`, the total binding-energy shift of the two
mouths in their equations (5.31) and (5.42). Thus the static throat profile,
redshift and RN mass shift agree in one normalization.

This is the leading overlap matching already carried out by MMP
[equations (5.39)-(5.42)](https://arxiv.org/html/1807.04726v3#S5.SS3),
rewritten as a boundary charge relation for the present constraint audit.
It is a reproduction, not a new result.

## What this fixes for dynamics

Along a differentiable family of static matched solutions with the same
`r_e` and renormalized coupling, (1) implies

```text
delta(deltaM)=-A x_* delta x.                               (2)
```

Equation (2) is the zero-frequency tangent of the static solution family.
It supplies the mass datum associated with a quasi-static change of the
throat clock. It must not be used as a local-in-time law for a perturbation
with nonzero frequency. In a source-free spherical RN collar, Birkhoff's
theorem instead gives `partial_t deltaM=0`; hence that collar cannot carry
the oscillating monopole load assumed by a freely prescribed `x(t)`.

With charged matter in a collar, the correct condition is the sourced mass
and Gauss balance in `sourced-spherical-mass-charge-balance.md`. For an
oriented timelike cut with unit tangent `w^a`, outward unit normal `n^a`,
the matter energy flux is `T_ab w^a n^b`; integrating it changes the RN
mass between cuts. The two mouths need separate balances, joined by the
retarded fermion energy/current transfer through the throat and exterior.
Their asymptotic time normalizations must also be matched separately before
imposing the symmetric clock condition.

At leading order, equations (1)-(2) provide the boundary charge condition
for the zero-frequency sector of the sourced JT initial-data solution.
They do not determine the three embedding constants at nonzero frequency:
those require the full collar stress/current response and the finite-cut
boundary conditions. Nor do they include the `O(r_e/ell)` overlap terms,
quantum corrections to RN, radiation outside spherical symmetry, or massive
angular modes.

Remaining to-do list: derive the retarded spherical stress/current flux in
the RN collars; match both boundary time normalizations and mass functions
to the throat solution; then solve the full frequency-dependent JT clock
problem and reassess stability.
