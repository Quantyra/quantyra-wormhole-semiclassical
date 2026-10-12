# Long-return enhancement is spatially rank one up to a uniform kernel bound

2026-10-12, S053 primary-response support. For the fixed product-throat
graph used by the interacting Gaussian model, the free resolvent restricted
to the throat admits an exact decomposition into a spatially constant part
and a remainder whose kernel is uniformly bounded independently of the
Euclidean frequency and all positive return lengths. This isolates the
long-return infrared enhancement before the repulsive charged potential is
resummed. It does **not** bound the interacting infrared response, the
shrinking-mouth limit, the moving-metric stress vertices, or the physical
curved exterior load.

## Statement and proof

Let the throat be the interval `I=[-A,A]` of length `T=2A`. Channel `j`
is a circle of length `L_j=T+D_j` with arbitrary `D_j>0`; all channels
are averaged with weights `1/q`. For `Omega>0`, the circle resolvent
kernel at two throat points is

```text
G_(L,Omega)(z,z')
 = cosh(Omega*(L/2-d_L))/(2 Omega sinh(Omega*L/2)),
d_L=min(|z-z'|,L-|z-z'|).
```

Define the restricted averaged kernel and its diagonal constant by

```text
K_Omega(z,z')=(1/q) sum_j G_(L_j,Omega)(z,z'),
c_Omega=(1/q) sum_j G_(L_j,Omega)(z,z).
```

Then, for every `z,z' in I`,

```text
K_Omega(z,z')=c_Omega+Kbar_Omega(z,z'),
|Kbar_Omega(z,z')|<=|z-z'|/2<=T/2.                 (1)
```

To prove this for one circle, set `x=d_L`. On `0<=x<=L/2`,

```text
|dG_(L,Omega)(x)/dx|
 = sinh(Omega*(L/2-x))/(2 sinh(Omega*L/2))<=1/2.
```

The kernel is symmetric under `x -> L-x`, so the same Lipschitz bound
holds as a function of circular distance. The shortest circular separation
between throat points is at most their interval separation, hence
`d_L<=|z-z'|`. The mean-value theorem gives
`|G(z,z')-G(z,z)|<=|z-z'|/2`. Averaging over channels proves (1).
The Schur test also gives

```text
||Kbar_Omega||_(L2(I)->L2(I)) <= T^2/2,             (2)
```

uniformly in `q`, `Omega`, and the complete set of positive `D_j`. The
constant part `c_Omega` may grow without bound in joint low-frequency and
long-return limits; all such unbounded dependence is spatially constant on
the throat.

For the charged rank-one potential `V_0(z)P_c`, with
`P_c=vv*`, `||v||=1`, the throat Birman–Schwinger operator from the
[interacting potential UV split](interacting-potential-uv-split.md) is

```text
B_Omega=V_0^(1/2) K_Omega V_0^(1/2)
       =c_Omega |sqrt(V_0)><sqrt(V_0)|+Bbar_Omega,
||Bbar_Omega||<=||V_0||_infinity T^2/2.             (3)
```

Thus its unbounded long-return/low-frequency component is rank one in the
throat coordinate. The remaining operator is bounded uniformly over all
positive returns at fixed mouth cutoff. Equation (3) is an exact
decomposition, not a perturbative expansion: the bound can grow like
`epsilon^-2` in the shrinking-mouth limit, and `I+Bbar_Omega` need not
be positive. Consequently a scalar Sherman–Morrison resummation cannot
be applied to the remainder without an additional positivity or spectral
argument. The physical contact still requires the full inverse
`(I+B_Omega)^-1` and its retarded continuation.

## Verification and scope

`checks/return_uniform_kernel.ps1` independently evaluates the exact
positive difference in a cancellation-resistant form over deterministic
logarithmic grids of `Omega`, `D`, and throat separations. It checks the
stronger `|G(0)-G(d)|<=d/2` bound and hashes this note and the script into
`evidence/return-uniform-kernel.json`. This numerical check supplements,
but does not replace, the derivative proof above.

The decomposition gives a concrete next step for the primary (l=0)
response: resummation should treat the spatially constant throat mode
nonperturbatively and control the bounded remainder without assuming it is
positive. This step addresses only the Gaussian potential sector. It does
not resolve the interacting joined contact identified in the
[proper-mouth seagull audit](proper-mouth-seagull-audit.md), prove a
physical damping pole, or claim novelty.

Remaining to-do list: obtain a positivity/spectral bound for the bounded
Birman–Schwinger remainder; control the shrinking-cutoff limit; include
kinetic, stress and joined-mouth vertices with one covariant subtraction;
match the curved exterior and state; then decide whether a physical
stability theorem survives.
