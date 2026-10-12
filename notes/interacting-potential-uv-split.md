# UV split for the moving-mouth potential contact

2026-10-11, S053. This isolates the ultraviolet structure of the
potential-profile velocity seagull in the prescribed Gaussian channel
model. It does not bound the full metric-stress response or determine the
finite joined contact.

## Fixed-geometry operator and response identity

Fix `q`, `epsilon>0`, and finite return lengths. Put the throat and returns
on their compact joined one-dimensional graph `G`, with the common
self-adjoint continuity/derivative matching used in the
[interacting pressure-band model](interacting-mouth-pressure-band.md).
On the charged subspace let

```text
H_0=-partial_z^2+V_0(z) P_c,
V_0(z)=a sec^2(z), P_c=v v*, ||v||=1,
R_Omega=(H_0+Omega^2)^(-1).
```

The neutral constant modes are annihilated by `P_c` and are omitted from
these traces. At fixed positive `epsilon`, the mouth potential and each
clock-variation insertion are bounded. Write the one-particle Euclidean
operator expansion as

```text
P[phi]=P_0+phi W_1+...+dot(phi)^2 W_2+...,
W_2(z)=2 mu z V_0(z) tan(z) P_c,
mu=epsilon^3/(2 A b^3), A=acos(epsilon), b=sqrt(1-epsilon^2),
```

where `W_2` is the multiplication-operator convention for the potential
term in equation (2) of the
[proper-mouth seagull audit](proper-mouth-seagull-audit.md). For a real
Gaussian determinant `Gamma=(1/2) Tr log(P)`, its velocity-contact
contribution at external Euclidean frequency `w` is

```text
Pi_contact(w)= (w^2/2) integral_R dOmega/(2 pi)
                         Tr_G[R_Omega W_2].                    (1)
```

The two-potential-vertex response has the standard resolvent form

```text
Pi_bubble(w)=-(1/2) integral_R dOmega/(2 pi)
             Tr_G[R_Omega W_1 R_(Omega+w) W_1],                 (2)
```

up to the fixed Fourier and field-multiplicity conventions. The sign in
(2) is the second variation of `+(1/2)Tr log(P)`; the ultraviolet estimate
below concerns its absolute value. Terms from kinetic/stress vertices,
zero-derivative contacts, anomaly and interface motion are outside (1)-(2).

## Local logarithm and exact profile integral

The one-dimensional large-`|Omega|` diagonal resolvent has leading local
term `R_Omega(z,z)=1/(2|Omega|)+O(|Omega|^-2)` after integration against a
bounded multiplier. The `O(|Omega|^-2)` remainder allows for the potential
step at a mouth vertex. Therefore

```text
Tr_G[R_Omega W_2]
  = (1/(2|Omega|)) integral_G tr(W_2(z)) dz + O(|Omega|^-2).
```

The frequency integral in (1), with a symmetric cutoff `Lambda`, has
logarithmic divergence

```text
Pi_contact(w)/w^2
  = (1/(4 pi)) log(Lambda/mu_R) integral_G tr(W_2) dz
    + finite_renormalized_part + o(1),                            (3)
```

where `mu_R` is fixed by the covariant renormalization prescription. This
is a local contact divergence; the coefficient is independent of the
return lengths. It must be matched with the local counterterm of the same
microscopic action, rather than selected by subtracting the desired
damping response.

Since `tr(P_c)=1` on the charged channel, the throat integral is elementary:

```text
integral_-A^A tr(W_2) dz
 = 2 mu a integral_-A^A z tan(z) sec^2(z) dz
 = 2 mu a (A/epsilon^2-b/epsilon)
 = a epsilon(A-b epsilon)/(A b^3).                              (4)
```

The primitive is `(z sec^2(z)-tan(z))/2`. Formula (4) is positive for
`0<epsilon<=1/4` and is `O(a epsilon)` as `epsilon -> 0`. It gives the
coefficient of the local logarithm in (3), not the finite renormalized
expectation. The seagull-audit quadratic-kernel convention is twice the
effective-action coefficient in (3).

## The nonlocal two-vertex tail is integrable

Let `L_G` be the total graph length and `K=||W_1||_infinity`. For
`|Omega|>=2|w|+1`, positivity and one-dimensional Weyl growth give

```text
||R_Omega|| <= |Omega|^-2,
||R_Omega||_1 <= C_G (L_G/|Omega|+|Omega|^-2),
|Tr(R_Omega W_1 R_(Omega+w) W_1)|
 <= C_G K^2 (L_G/|Omega|^3+|Omega|^-4).                         (5)
```

Indeed, apply `|Tr(AB)|<=||A||_1 ||B||` with
`A=R_Omega` and
`B=W_1 R_(Omega+w) W_1`; the second inequality follows from the compact
one-dimensional eigenvalue counting law. Thus the bubble has a finite
ultraviolet tail, bounded beyond frequency `R>2|w|+1` by
`C_G K^2(L_G/R^2+R^-3)`. This bound is for each fixed graph and fixed
mouth cutoff; its constant is not uniform as return lengths grow or
`epsilon -> 0`.

## Consequence and remaining gap

The potential-sector UV problem separates into a single local logarithm,
whose spatial coefficient is exactly (4), and a convergent two-vertex
response. This is the appropriate starting point for matching the finite
contact to the covariant 1+1 subtraction and then evaluating its retarded
finite part. It does not determine that finite part, bound the kinetic
stress-stress tail, or control the physical exterior boundary as the
geometry moves.

Bilal and Leduc's exact determinant variation for a constant-mass scalar
on a compact boundaryless surface gives the same general local/global
separation through a coincident Green function and a convergent heat-kernel
remainder ([arXiv:1606.01901v2, §§2-3](https://arxiv.org/html/1606.01901v2)).
Their operator and boundary conditions differ from this variable
rank-one potential on a transmitting graph; no determinant formula from
that paper is imported here.

Remaining to-do list: fix the finite local coefficient by the chosen
matched covariant action; evaluate the finite retarded potential response;
bound kinetic/stress vertices and the full UV remainder; then close physical
matching, novelty selection, certification and manuscript preparation.
