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

The two profiles also obey the sharper pointwise comparison

```text
W_2/V_0=2 mu z tan(z)
      =epsilon^3 z tan(z)/(A b^3)
      <=epsilon^2/b^2,                 -A<=z<=A.               (5a)
```

Here `z tan(z)` is even and increasing in `|z|`, and
`tan(A)=b/epsilon`. Hence `0<=W_2<=epsilon^2 b^-2 V_0 P_c` as
operators on the charged throat sector. For every `Omega>0`, trace
monotonicity gives
`0<=Tr(R_Omega W_2)<=epsilon^2 b^-2 Tr(R_Omega V_0 P_c)`.
This exposes an `O(epsilon^2)` relative insertion factor, potentially
useful for the shrinking-mouth family. It does not itself control the
weighted potential resolvent in the infrared uniformly over arbitrarily
long returns; no integrated small-mouth response bound is inferred.

This weighted trace is exactly a coupling derivative of the same finite
graph determinant. Write `V_0=a U`, with `U=sec^2(z)P_c`, and
`B_Omega(a)=a U^(1/2)K_Omega U^(1/2)`. The finite-cutoff identity is

```text
Tr_G(R_Omega V_0 P_c)
   =Tr_th(B_Omega(1+B_Omega)^(-1))
   =partial_(log a) log det(1+B_Omega).                         (5b)
```

The first equality follows from the resolvent identity and cyclicity;
the second is the Fredholm determinant derivative. Combining with (5a)
reduces the positive-frequency contact estimate to a renormalized
`log(a)` derivative of the Gaussian determinant, with a prefactor at
most `epsilon^2/b^2`. The integral over frequency still requires
the derivative of the same local Born subtraction and covariant term.
Existing static determinant bounds concern clock derivatives at fixed
`a`; they do not supply this coupling derivative or its uniformity in
long returns. This is an exact reduction of the missing estimate, not
its resolution.

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

## Exact throat-only resolvent reduction

There is also an exact reduction of the finite-cutoff resolvent that
isolates where the finite return dependence enters. Let `J` embed a throat
function `f` as the channel vector `v f` on the full graph, and let
`R^0_Omega=(H_free+Omega^2)^(-1)`. Define the positive throat operator

```text
K_Omega=J* R^0_Omega J,
B_Omega=V_0^(1/2) K_Omega V_0^(1/2).
```

The resolvent identity for `H_0=H_free+J V_0 J*` gives

```text
R_Omega=R^0_Omega
 -R^0_Omega J V_0^(1/2)(1+B_Omega)^(-1)
     V_0^(1/2)J*R^0_Omega.                                    (12)
```

This follows by multiplication by `H_0+Omega^2` and
`(1+B)^(-1)+B(1+B)^(-1)=1`. For a bounded throat multiplication
insertion `W=w(z)P_c`, cyclicity of the trace reduces the interacting
trace to throat operators only:

```text
Tr_G(R_Omega W)=Tr_th(K_Omega w)
 -Tr_th((1+B_Omega)^(-1) V_0^(1/2) K_Omega w K_Omega
                                      V_0^(1/2)).              (13)
```

For fixed positive `Omega`, finite graph and positive `epsilon`, these
are ordinary resolvent identities. `K_Omega` is the average of the
individual free-circle kernels on the throat, with each return length
entering explicitly as

```text
K_Omega(z,z')=(1/q) sum_j G_(L_j,Omega)(z,z'),
G_(L,Omega)(z,z')=
 cosh(Omega*(L/2-d_L(z,z')))/(2 Omega sinh(Omega*L/2)),
d_L(z,z')=min(|z-z'|, L-|z-z'|),   L_j=T+D_j.
```

The distance formula applies because throat points have
`|z-z'|<=T<L_j`. Thus the finite return-dependent part of the contact
reduces to a one-dimensional integral operator, without assuming equal
returns or replacing the potential by a boundary condition.

For the positive seagull profile `W_2`, this reduction gives a uniform
ultraviolet estimate on the *return-dependent difference* from the free
trace. The resolvent identity and `0<=(1+B_Omega)^(-1)<=1` imply

```text
0<=Tr((R^0_Omega-R_Omega)W_2)
  =Tr((1+B_Omega)^(-1) V_0^(1/2)K_Omega W_2 K_Omega V_0^(1/2))
 <=||V_0||_infinity Tr(W_2 K_Omega^2)
 <=||V_0||_infinity ||K_Omega|| Tr(W_2 K_Omega)
 <=C ||V_0||_infinity I_2 Omega^(-3),   Omega>=1.               (14)
```

Here `||K_Omega||<=Omega^(-2)` because `J` is an isometry and
`||R^0_Omega||<=Omega^(-2)`; the last step is the free diagonal bound
(8). Consequently the integral of this difference over
`|Omega|>=R>=1` is at most
`C ||V_0||_infinity I_2 R^(-2)`, uniformly in `q` and all positive
returns. With `||V_0||_infinity=a/epsilon^2` and
`I_2<=a epsilon/b^3`, the explicit coefficient is
`C a^2/(epsilon b^3)`. The estimate proves convergence of the
return-dependent ultraviolet correction at every fixed mouth cutoff;
its coefficient is not uniform as `epsilon->0`.

The return-dependent ultraviolet correction is now controlled at fixed
`epsilon`, but this identity does not bound the low-frequency contribution
or give a coefficient uniform as `epsilon->0`. The return-uniform tail in
(9) controls the two-vertex bubble, whereas this formula is for the
one-resolvent seagull trace. Kinetic and stress insertions also remain
outside its scope. The local subtraction remains the fixed one in (6)-(7).

## Consequence and remaining gap

The potential-sector UV problem separates into a single local logarithm,
whose spatial coefficient is exactly (4), a return-uniform ultraviolet
correction to the one-vertex seagull trace at fixed mouth cutoff, and a
convergent two-vertex response with a return-uniform high-frequency tail.
This is the appropriate starting point for evaluating the global finite
part and its retarded continuation. It does not bound the kinetic
stress-stress tail or control the physical exterior boundary as the
geometry moves.

Bilal and Leduc's exact determinant variation for a constant-mass scalar
on a compact boundaryless surface gives the same general local/global
separation through a coincident Green function and a convergent heat-kernel
remainder ([arXiv:1606.01901v2, §§2-3](https://arxiv.org/html/1606.01901v2)).
Their operator and boundary conditions differ from this variable
rank-one potential on a transmitting graph; no determinant formula from
that paper is imported here.

Remaining to-do list: control the low-frequency and shrinking-cutoff
finite potential response and its retarded continuation; bound
kinetic/stress vertices;
close physical matching, novelty selection, certification and manuscript
preparation.
