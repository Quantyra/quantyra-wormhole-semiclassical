# Long-return enhancement is spatially rank one up to a uniform kernel bound

2026-10-12, S053 primary-response support. For the fixed product-throat
graph used by the interacting Gaussian model, the free resolvent restricted
to the throat admits an exact decomposition into a spatially constant part
and a remainder whose kernel is uniformly bounded independently of the
Euclidean frequency and all positive return lengths. This isolates the
long-return infrared enhancement before the repulsive charged potential is
resummed. The zero-trace remainder also yields a uniform low-frequency
bound for the positive potential seagull, including the shrinking-mouth
limit. It does **not** bound the full interacting response, moving-metric
stress vertices, or the physical curved exterior load.

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
be positive. However, the compression of `Bbar_Omega` orthogonal to its
rank-one vector is automatically positive. Set

```text
u=sqrt(V_0) 1, e=u/||u||, Q=I-|e><e|,
C_Q=I_Q+Q Bbar_Omega Q,
b_Q=Q Bbar_Omega e,
a_Q=1+c_Omega ||u||^2+<e,Bbar_Omega e>,
s_Q=a_Q-<b_Q,C_Q^(-1)b_Q>.
```

For `v` orthogonal to `u`,
`<v,Bbar_Omega v>=<v,B_Omega v> >= 0`, so `C_Q >= I_Q` and
`||C_Q^(-1)||<=1`. The scalar Schur complement obeys `s_Q>=1` because
`I+B_Omega>=I`. The full inverse therefore has the exact block formula

```text
(I+B_Omega)^(-1) =
  [ s_Q^(-1)                 -s_Q^(-1) b_Q* C_Q^(-1)             ]
  [ -C_Q^(-1)b_Q s_Q^(-1)     C_Q^(-1)+C_Q^(-1)b_Q s_Q^(-1)b_Q*C_Q^(-1) ]. (4)
```

This reduces the return-enhanced inversion to one scalar Schur complement
and an orthogonal block that is positive and bounded at fixed mouth
cutoff. It is a structural reduction of the exact Euclidean resolvent, not
a bound on the frequency-integrated response: `||Bbar_Omega||` may grow
as `epsilon^-2`, the `Q` block still depends on `Omega` and the returns,
and the moving stress vertices and retarded continuation remain to be
controlled.

## Uniform low-frequency bound for the potential seagull

The diagonal of each free circle Green kernel is independent of the
position `z`, so `K_Omega(z,z)=c_Omega` and the weighted remainder satisfies

```text
Tr(Bbar_Omega)=integral_I V_0(z)[K_Omega(z,z)-c_Omega] dz=0.      (5)
```

It is trace class for every fixed `Omega>0` and positive mouth cutoff.
Let `e=u/||u||` and `Q=I-|e><e|` as above. Since
`B_Omega>=0`, its rank-one part vanishes on `Q`, and therefore
`D_Q=Q Bbar_Omega Q=Q B_Omega Q>=0`. Equation (5) gives

```text
Tr(D_Q)=-<e,Bbar_Omega e><=||Bbar_Omega||
       <=||V_0||_infinity T^2/2.                                  (6)
```

Use the Schur block formula (4), with `C_Q=I_Q+D_Q`. Taking the trace of
`B_Omega(I+B_Omega)^(-1)=I-(I+B_Omega)^(-1)` yields

```text
Tr(B_Omega(I+B_Omega)^(-1))
 =1-s_Q^(-1)+Tr(I_Q-C_Q^(-1))-s_Q^(-1)||C_Q^(-1)b_Q||^2
 <=1+Tr(D_Q)
 <=1+||V_0||_infinity T^2/2.                                    (7)
```

The middle inequality uses `s_Q>=1`, `C_Q>=I_Q`, and
`0<=I_Q-C_Q^(-1)<=D_Q`. This is uniform over all positive return
lengths, channel counts, and `Omega>0`; no low-frequency gap is assumed.
The exact coupling-derivative identity identifies the left side with
`Tr_G(R_Omega V_0 P_c)`.

From the pointwise moving-cut inequality
`0<=W_2<=epsilon^2 b^-2 V_0 P_c`, trace monotonicity now gives

```text
0<=Tr_G(R_Omega W_2)
 <=epsilon^2/b^2 * (1+a T^2/(2 epsilon^2))
 =(epsilon^2+a T^2/2)/b^2
 <=(epsilon^2+a pi^2/2)/b^2,                                  (8)
```

because `||V_0||_infinity=a/epsilon^2` and `T=2 acos(epsilon)<pi`.
Consequently, for any fixed `Omega_0`, the positive-frequency contact
integral over `0<Omega<Omega_0` is at most
`Omega_0 (epsilon^2+a pi^2/2)/b^2`, uniformly in the returns and as
`epsilon->0`. This controls the potential seagull's low-frequency
contribution in the prescribed Gaussian graph model. It does not bound
the higher-frequency renormalized finite part, the two-vertex bubble, the
kinetic/stress response, or the retarded physical kernel.

The higher-frequency **interaction correction** can also be bounded
uniformly as the cutoff shrinks. Let
`T_min=2 acos(1/4)` and
`C_0=coth(T_min/2)/2`. For `Omega>=1` and every return circle,

```text
c_Omega(L)=coth(Omega L/2)/(2 Omega)<=C_0/Omega,
Tr(R_Omega^free W_2)=I_2 (1/q) sum_j c_Omega(L_j).
```

Positivity and form ordering give
`0<=Delta_Omega=Tr((R_Omega^free-R_Omega)W_2)
<=Tr(R_Omega^free W_2)`. Split at `Omega=epsilon^-1` (which exceeds
one for the stated cutoff range). On `1<=Omega<=epsilon^-1`, integration
of this free-trace bound gives
`C_0 I_2 log(1/epsilon)`. On `Omega>=epsilon^-1`, the resolvent-identity
tail (14) of the UV-split note gives
`Delta_Omega<=C_0 ||V_0||_infinity I_2 Omega^-3`, so its integral is at
most `(C_0/2)||V_0||_infinity I_2 epsilon^2`. Using
`I_2<=a epsilon/b^3` and `||V_0||_infinity=a/epsilon^2`,

```text
integral_1^infinity Delta_Omega dOmega
 <= C_0 [a epsilon log(1/epsilon)+a^2 epsilon/2]/b^3.           (9)
```

This is uniform in channel count and all positive returns, and vanishes
as `epsilon->0` at fixed `a`. The free high-frequency trace separates
exactly as

```text
c_Omega(L)=1/(2 Omega)+1/[Omega (exp(Omega L)-1)].               (10)
```

The first term is the sole logarithmic divergence with coefficient `I_2`.
The finite return-dependent remainder above `Omega=1` is bounded by
`I_2 E_1(T_min)/(1-exp(-T_min))`, where
`E_1(t)=integral_t^infinity exp(-x)/x dx`; it is therefore `O(a epsilon)`.
Together with the fixed covariant local conversion (7), whose size is
`O(a epsilon[1+|log(lambda_0/epsilon)|])` at fixed renormalization scale,
(9)-(10) control the renormalized high-frequency potential seagull in the
shrinking-cutoff limit. This does not control the low-frequency
interaction correction separately (only the full positive low-frequency
trace is bounded by (8)), and does not extend to the bubble or kinetic and
stress vertices.

## Verification and scope

`checks/return_uniform_kernel.ps1` independently evaluates the exact
positive difference in a cancellation-resistant form over deterministic
logarithmic grids of `Omega`, `D`, and throat separations. It checks the
stronger `|G(0)-G(d)|<=d/2` bound and the high-frequency free-diagonal and
return-remainder inequalities used in (9)-(10), then hashes this note and
the script into `evidence/return-uniform-kernel.json`. These numerical
checks supplement, but do not replace, the analytic operator and trace
proofs above.

The decomposition and (8)-(10) control the low-frequency positive
potential seagull and its renormalized high-frequency interaction part.
The next step is to bound the two-vertex bubble, then the kinetic/stress
response and moving join. This step addresses only the Gaussian potential
sector. It does
not resolve the interacting joined contact identified in the
[proper-mouth seagull audit](proper-mouth-seagull-audit.md), prove a
physical damping pole, or claim novelty.

Remaining to-do list: bound the two-vertex bubble and kinetic/stress
response with the same covariant subtraction; derive the interacting
joined-mouth contact; match the curved exterior and state; then decide
whether a physical stability theorem survives.
