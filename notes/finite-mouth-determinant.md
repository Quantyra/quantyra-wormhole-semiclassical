# Finite-mouth Gaussian determinant at fixed mouth lapse

2026-10-11, E003/S053. The entire renormalized Gaussian electric
determinant has small first and second clock derivatives relative to
q delta in the retained product-throat family, even at the deep cutoff.
This follows without replacing the finite potential by ideal boundary
conditions. Isolating fixed mouth cells controls their interaction
force; the remaining potential has small integrated strength.

The result closes the finite-product-model static premise left by
[the local calculation](proper-mouth-local-stress.md). It does not
identify the proper-mouth clock with the canonical dynamical vertex,
control the corrected physical background's nonlocal determinant,
select compact/angular sectors, or construct the four-dimensional
quantum state. Physical damping and consequential novelty remain open.

## Model and statement

Use the same q separate real current circles, with a common throat and
arbitrary fixed D_j>0. The unit electric singlet is v_j=1/sqrt(q).
On the throat prescribe

```text
x=exp(phi), 1/2<=x<=2, 0<epsilon<=1/32,
b=sqrt(1-epsilon^2 x^2), A=acos(epsilon x), T=2A/x,
V_x(s)=a x^2 sec^2(xs),       -T/2<=s<=T/2.
```

Here a>0, epsilon, D_j and the exterior time normalization lambda_0
are fixed under clock differentiation. The potential matrix is V_x vv*;
it vanishes on the returns. Channelwise continuity is retained at both
mouths. Each finite cutoff is regular. Zero oscillator energies use the
same convention as the previous Gaussian notes, without choosing the
compact charge or holonomy state.

For fixed mu>0 and finite coefficient c, let

```text
E_opt=(1/(2pi)) integral_0^infinity
 [log det(1+V_x^(1/2) G_0,k V_x^(1/2))
                    -I_x/(2 sqrt(k^2+mu^2))] dk + c I_x,
I_x=2a b/epsilon,
E_cov=E_opt+J_x/(4pi),
J_x=2a[(b/epsilon)(log(lambda_0/epsilon)-1)+x A].          (1)
```

G_0,k is the singlet-restricted free circle resolvent. These are
massive-minus-massless energies on the same optical geometry; the
mass-independent Weyl anomaly cancels. The subtraction is the established
[static singlet prescription](singlet-static-curvature.md), not a new
choice of finite counterterm.

Let Lbar=T+average_j D_j. For sufficiently small a, uniformly on this
clock interval, the first two derivatives obey

```text
|D^j E_opt| <= C [1+sqrt(a)+a Lbar
           +a(1+log(1+a/epsilon)+|log mu|+|c|)], j=1,2.  (2)
```

The constant is independent of q, individual D_j, a and epsilon.
The local J derivatives have the already exact bound
O(a+a epsilon[1+|log(lambda_0/epsilon)|]). In particular, for

```text
a=1/(4pi^2 q), epsilon=q^-2, delta=log^-2 q,
average D_j=O(delta log q), lambda_0=1/(16pi q^3),
```

and fixed mu,c (or coefficients/scales with the stated logarithmic
cost controlled), both derivatives of E_cov are O(1)=o(q delta).
No uniform bound on the undifferentiated energy is asserted: it contains
large clock-independent mouth self energies. No claim about a physical
four-dimensional matching coefficient follows from holding c fixed.

## Fixed mouth cells and a genuinely small remainder

Take h=1/8 and L=T-2h>1. Separate a free central interval of length L
from two mouth cells of fixed optical width h. In a cell let d be the
distance inward from its mouth and put

```text
P(d)=a/(epsilon+d)^2,       0<=d<=h,
P=0 on the central interval and the returns,
W=V_x-P.
```

The reference cells move with their respective mouths while their
intrinsic profiles remain fixed. Their outer interfaces are attached
to the fixed returns. The artificial step at d=h is only a decomposition
of the original finite potential, not a physical junction or limit.
Both the P and W forms are nonnegative. Indeed

```text
p_x(d)=sin(asin(epsilon x)+xd)/x
      =epsilon cos(xd)+b sin(xd)/x <= epsilon+d,
V_x(d)=a/p_x(d)^2.
```

The inequality follows from sin(u+t)<=sin u+t for u,t>=0 in this
range. Moreover p_x(d) is bounded below by a fixed positive multiple
of epsilon+d. Write r=epsilon+d. The nonnegative deficit is

```text
r-p_x(d)=x^2[epsilon d^2 A_1(xd)+d^3 A_2(xd)
                   +epsilon^2 d sinc(xd)/(1+b)],
A_1(u)=(1-cos u)/u^2, A_2(u)=(1-sinc u)/u^2.             (3)
```

Extend these functions at zero by 1/2 and 1/6. On |xd|<=1/4,
|epsilon x|<=1/16 they are positive smooth functions with bounded
first two logarithmic clock derivatives. The same holds for p_x/r.
Thus, using W=a(r-p_x)(r+p_x)/(p_x^2 r^2),

```text
0<=W(d)<=C a x^2 d/(epsilon+d),
|D^j W(d)|<=C_j W(d), j=1,2.                            (4)
```

At d=0 these statements extend by continuity. In the central interval
use a fixed normalized coordinate. It stays at least h from either
mouth, so V_x/a and its first two logarithmic derivatives are uniformly
bounded there, including its measure factor L. Consequently the fixed
piecewise form domain has uniform relative derivative bounds for W,
the free form, and the background free-plus-P form. No large derivative
of a translating barrier has been discarded; its translation is absorbed
in the moving cell coordinates and the changing central length.

The strengths are exactly

```text
M=integral P=2a[1/epsilon-1/(epsilon+h)],
eta=integral W=2a[1/(epsilon+h)-epsilon x^2/(1+b)]<=2a/h. (5)
```

M can grow as q; eta tends to zero. This is why the earlier small-I
theorem could not be applied directly to V_x, but can be extended to
the relative W determinant over the strong background.

## The interaction force of the fixed cells

Cut at the two inner cell ends. The outside consists of both fixed
cells and all returns. Its Dirichlet-to-Neumann matrices at k>0 are
N_P(k) and N_0(k), before and after removing P, on C^(2q). They are
positive. The variational characterization of the Dirichlet energy
gives N_P>=N_0. Their difference has rank at most four.

To prove that rank assertion, let u_0(b) and u_P(c) be the outside
solutions for boundary vectors b,c. Integration by parts gives

```text
c* (N_P-N_0) b = integral u_P(c)* P vv* u_0(b).           (6)
```

On each cell the free projection v* u_0 is a linear combination of
exp(kd) and exp(-kd). Therefore P^(1/2)v* u_0(b), as b varies, lies
in a space of dimension at most two per cell. Equation (6) factors
through their direct sum. This uses the fixed singlet on both cells,
not deletion of charged circles or an assumption of equal returns.

The central interval has the matrix

```text
K_L=k [[coth(kL) I,-csch(kL) I],[-csch(kL) I,coth(kL) I]].
F(k,L)=log det(N_P+K_L)-log det(N_0+K_L).                 (7)
```

Gaussian elimination at the cut factors the full determinant into the
outside Dirichlet determinant, the central Dirichlet determinant, and
the determinant of the sum of the two boundary matrices. The central
factors cancel between P and zero; the outside factors and its Born
subtraction depend on k but not L. Equivalently this follows by solving
the Dirichlet problem on each segment, then completing the boundary
quadratic form. Thus all L derivatives of the relative log determinant
are the derivatives of F. No generally divergent constant is integrated
and then silently subtracted.

Put R_i=(N_i+K_L)^(-1), Delta=R_0-R_P. Then
0<=Delta<=K_L^(-1) and rank Delta<=r, r=min(4,2q). The symmetric and
antisymmetric eigenvalues of K_L are k tanh(kL/2) and k coth(kL/2).
Direct differentiation gives the form bounds

```text
-c_1 K_L <= K_L' <= c_1 K_L, c_1=k/sinh(kL),
-c_2 K_L <= K_L''<= c_2 K_L, c_2=k^2/[2sinh^2(kL/2)].
```

Since Tr(Delta K_L)<=r,

```text
F_L=-Tr(Delta K_L'),              |F_L|<=r c_1,
F_LL=-Tr(Delta K_L'')
      +Tr(R_0 K_L' R_0 K_L'-R_P K_L' R_P K_L'),
|F_LL|<=r(c_2+2c_1^2).                                  (8)
```

For the last line, expand the difference into
Tr(Delta K_L' R_0 K_L')+Tr(Delta K_L' R_P K_L'). Each is
nonnegative and at most c_1^2 Tr(Delta K_L), because
K_L' R_i K_L'<=K_L' K_L^(-1) K_L'<=c_1^2 K_L.
These bounds are integrable at zero and infinity. The zero-frequency
oscillators are retained through the k down to zero limit.

Integrating the geometric expansions of csch gives

```text
integral_0^infinity k/sinh(kL) dk=pi^2/(4L^2),
integral_0^infinity (c_2+2c_1^2) dk=pi^2/L^3,
|partial_L E_P|<=r pi/(8L^2),
|partial_L^2 E_P|<=r pi/(2L^3).                          (9)
```

This is uniform in barrier height and all unequal returns. The exterior
can have arbitrarily small gaps; the proof uses positivity and finite
boundary rank instead. Because L=T-2h and

```text
D T=-T-2epsilon/b,
D^2 T=T+2epsilon/b-2epsilon^3 x^2/b^3,
```

its first two clock derivatives are bounded on the fixed interval.
Equation (9) therefore bounds D E_P and D^2 E_P by numerical constants.
The c M counterterm is clock independent in these coordinates.

## Small relative determinant over a strong positive background

Here is the needed extension of the existing small-potential argument.
Let P,W>=0 have integrated strengths M,eta, with eta<=1, on the
same rank-one throat space. Suppose the free, P and W forms obey
uniform first two relative derivative bounds; these imply the bounds
for the free-plus-P form. In this application P'=P''=0.
Let Delta E be the optical energy of P+W minus that of P, with the
same mu,c. Then, for j=0,1,2,

```text
|D^j Delta E| <= C [sqrt(eta)+eta Lbar
       +eta(1+log(1+M)+|log mu|+|c|)].                   (10)
```

For completeness all three frequency ranges can be bounded directly.
Write R_P=(A_k+P)^(-1), ell_k=log det(1+W^(1/2)R_P W^(1/2)).
Positivity gives R_P<=R_0. The free circle zero-mode decomposition
therefore bounds the low-frequency integral by
pi sqrt(eta/T_min)+eta Lbar/12. Applying the exact differentiated
relative-determinant identities to background A_k+P and increment W
bounds |D^j ell_k| by C_j ell_k, j=1,2, as in section 3 of the
earlier static note. It does not require a lower spectral gap.

For 1<=k<=K, with K=C_0(1+M+eta) sufficiently large,
ell_k<=Tr(R_0 W)<=C eta/k. This gives the logarithm in (10).
The same derivative identities apply, and the local subtraction costs
at most C eta[log K+1+|log mu|] over this range and k<=1.

For k>=K expand both log determinants relative to the free operator.
Their first Born difference is Tr(R_0 W). Every term in the remaining
difference contains at least one W and at least two potential factors.
The free trace bounds are Tr(R_0 P)<=C M/k and
Tr(R_0 W)<=C eta/k, while their sum is at most 1/2 by the choice of K.
Schatten Holder with one W trace factor bounds the remainder and its
first two derivatives by

```text
|D^j[ell_k-Tr(R_0 W)]|<=C_j eta(M+eta)/k^2.              (11)
```

To check differentiated convergence, differentiate the resolvent as
R_0'=-R_0 A_k' R_0 and
R_0''=2R_0 A_k' R_0 A_k' R_0-R_0 A_k'' R_0. The relative form
bounds insert bounded factors between the positive square roots.
Each n-th order term has at most a fixed polynomial in n placements,
bounded by C_j n^p b_W(b_P+b_W)^(n-1). Its sum is bounded by
C_j b_W(b_P+b_W); this proves (11) and trace-norm convergence.
Integrating (11) costs C eta(M+eta)/K<=C eta.

Finally the free diagonal is known exactly:

```text
Tr(R_0 W)=eta average_j coth(k(T+D_j)/2)/(2k).
```

Its difference from eta/(2k), including two derivatives, decays
exponentially at a rate fixed by T_min. The difference between
1/k and 1/sqrt(k^2+mu^2) supplies only the displayed scale logarithm.
The finite term c eta has the same derivative control. These estimates
prove (10) without a smallness assumption on M. Substituting (5),
then adding (9) and the exact local J derivatives, proves (2).

## Scope, verification and decision

The determinant, fixed-cell force and Born-subtraction methods are
established. [Harrison and Kirsten](https://arxiv.org/abs/1201.3754)
derive graph vacuum energies and distinguish fixed supported potentials
from length-dependent subtraction terms. The present proof uses an
explicit boundary elimination and a strong-background remainder estimate
for the already selected current model. It does not infer vacuum-energy
convergence from the earlier norm-resolvent boundary limit.

The [checker](../checks/finite_mouth_determinant.py) and
[receipt](../evidence/finite-mouth-determinant.json) test the exact
decomposition, corrected cutoff derivatives, finite-cell boundary rank
and force inequalities, independent determinant assembly, and quadrature
refinement. The [source supplement](../evidence/finite-mouth-determinant-source-review.json)
records the bounded prior-work review. These numerical tests supplement
the ordinary uniform estimates; no Lean certificate is claimed.

Continue physical long-return damping as provisional primary and matched
finite-duration transfer as fallback. The product Gaussian static
premise now holds at the deep proper cut. The actual time-dependent
mouth vertex, nonlocal corrected-background response, angular/compact
state and physical exterior remain the next decisive model questions.
Neither finite-disturbance survival nor a common-observer shortcut follows.
Do not create a separate graph research target or formalization campaign
from this supporting bound; consequential novelty remains unaccepted.

Remaining to-do list: derive the proper-mouth dynamical vertex and
control the corrected-background determinant and coupled state/exterior;
establish the physical endpoint and novelty, independently verify,
certify suitable mathematics on GCP, and prepare the manuscript.
