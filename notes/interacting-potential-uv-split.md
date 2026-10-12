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
effective-action coefficient in (3). Since `z tan(z) sec^2(z)` increases
on `[0,A]`, the insertion also obeys the exact bounds

```text
0<=W_2(z)<=a/b^2,
I_2<=a epsilon/b^3,
||W_2||_2^2<=a^2 epsilon/b^5.                               (5)
```

Thus its spatial `L2` norm is `O(a sqrt(epsilon))` uniformly over
`epsilon<=1/4`; the local covariant integral `J_2` is
`O(a epsilon[1+|log(lambda_0/epsilon)|])` for fixed positive `lambda_0`.

Define the finite optical resolvent part, excluding the finite local
coefficient `c`, by

```text
R_opt(W_2;mu_R)=1/2 lim_(Lambda->infinity) [
   integral_-Lambda^Lambda dOmega/(2 pi) Tr_G(R_Omega W_2)
   -(1/(2 pi)) log(Lambda/mu_R) I_2 ].                          (6)
```

The subtraction is exactly the divergence in (3). It can equivalently be
defined with the subtracted heat kernel, so it is finite on each fixed
graph and captures the global channel/return dependence.

The already fixed covariant prescription makes the finite **local** part
explicit as well. At the reference clock `phi=0`, the optical lapse is
`N(z)=lambda_0 sec(z)`. In the notation of the proper-mouth local audit,
the local conversion is `c I+J/(4 pi)`, with `I=integral V` and
`J=integral V log(N)`. Applied to this insertion it is

```text
I_2=integral W_2 dz=a epsilon(A-b epsilon)/(A b^3),
J_2=integral W_2 log(N) dz
   =2 mu a [ (A/epsilon^2-b/epsilon) log(lambda_0/epsilon)
             -A/(2 epsilon^2)+3b/(2epsilon)-A ],
C_local=c I_2+J_2/(4 pi).                                    (7)
```

For the second identity, integrate `z tan(z) sec^2(z)` using the primitive
above and integrate the additional `log(sec z)` by parts. Thus the
renormalized potential contact in the specified prescription is exactly
`R_opt+C_local`; only the global finite resolvent part still depends on the
return lengths. This statement varies
the potential insertion only. Variations of the lapse, measure and kinetic
operator contribute separate local terms to the full clock response.

## Return-uniform bound for the two-vertex tail

For the potential vertex, `W_1(z)=w_1(z)P_c` is supported on the throat
and `P_c=v v*` has rank one with `v_j=1/sqrt(q)`. Let
`H_free` be the direct sum of the q free circles of lengths
`L_j=T+D_j`, and `H_0=H_free+V_0P_c`. For `Omega>0`, form ordering gives
`R_Omega<=R_Omega^free`. The diagonal free resolvent on circle j is

```text
G_(L_j,Omega)(z,z)=coth(Omega L_j/2)/(2 Omega).
```

For `epsilon<=1/4`, `T=2 acos(epsilon)>=2 acos(1/4)>2`; hence for
`Omega>=1` the diagonal is at most `C/Omega`, with a universal C. The
normalized charged projector then gives

```text
Tr(R_Omega W_1^2)
 <= Tr(R_Omega^free W_1^2)
 <= (C/Omega) integral_-A^A |w_1(z)|^2 dz.                    (8)
```

This bound is independent of q and of every return length for a fixed
profile `w_1`. For real Euclidean external frequency `w`, positivity also
gives the exact Hilbert-Schmidt identity and bound

```text
Tr(R_Omega W_1 R_(Omega+w) W_1)
 = ||R_Omega^(1/2) W_1 R_(Omega+w)^(1/2)||_2^2
 <= ||R_(Omega+w)|| Tr(R_Omega W_1^2)
 <= C/Omega^3 integral_-A^A |w_1(z)|^2 dz,                  (9)
```

provided `|Omega|>=2|w|+1`. Integrating (8) shows that the bubble tail
beyond `R>2|w|+1` is at most
`C R^-2 integral |w_1|^2`, uniformly in all positive returns and in q.
The bound is not uniform as `epsilon->0` or if the vertex profile norm
grows; those scalings are retained explicitly in its right-hand side.

## Consequence and remaining gap

The potential-sector UV problem separates into a single local logarithm,
whose spatial coefficient is exactly (4), and a convergent two-vertex
response with a return-uniform high-frequency tail. This is the appropriate
starting point for evaluating the global finite part and its retarded
continuation. It does not bound the kinetic stress-stress tail or control
the physical exterior boundary as the geometry moves.

Bilal and Leduc's exact determinant variation for a constant-mass scalar
on a compact boundaryless surface gives the same general local/global
separation through a coincident Green function and a convergent heat-kernel
remainder ([arXiv:1606.01901v2, §§2-3](https://arxiv.org/html/1606.01901v2)).
Their operator and boundary conditions differ from this variable
rank-one potential on a transmitting graph; no determinant formula from
that paper is imported here.

Remaining to-do list: evaluate the return-dependent finite potential
response and its retarded continuation; bound kinetic/stress vertices;
close physical matching, novelty selection, certification and manuscript
preparation.
