# Proper mouth time and the local mass contribution

2026-10-11, E003/S053. The order-q local clock curvature in the
[canonical boundary diagnostic](collective-boundary-casimir.md) does not
persist in the prescribed family with fixed mouth proper time. The large
local integral is mainly a clock-independent cutoff term. Its first two
clock derivatives are small, even at the deeper cutoff. The covariant
coincident-subtraction term can also be integrated exactly; its leading
clock response is a x/4, of order 1/q in the retained scaling.

This statement is about the local part of the specified spherical
Gaussian current model. It does not evaluate the nonlocal finite-mouth
determinant or construct the physical quantum state. The earlier
canonical calculation remains correct for its different variation.
The distinction is essential before transferring either calculation to
the physical clock's dynamical response.

## Fixed proper time and the established local prescription

Use reference throat-time units, x=exp(phi), lambda_0=P/ell_0, and
a=q g_c^2/(4pi^2). Prescribe the product throat

```text
ds_2^2=P^2[-x^2(1+r^2)dt^2+dr^2/(1+r^2)],
N=lambda_0 x sqrt(1+r^2),     m_cov^2=a/P^2.
```

Here N is the physical lapse relative to the reference exterior time.
The factor lambda_0 restores that time normalization in local logarithms.
Fix the mouth lapse N_b=lambda_0/epsilon, with epsilon>0 independent
of phi. On a compact clock interval assume z=epsilon x<=1/4. Then

```text
b=sqrt(1-z^2), R=b/z, A=atan R=acos z,
T=2A/x, y=atan(r)/x, V=m_opt^2=a x^2 sec^2(xy).            (1)
```

The ratio 1/epsilon is the leading radial cutoff at x=1; it is not held
fixed in r while changing x. This is the same proper-cutoff convention
as [the earlier clock calculation](proper-clock-anomaly.md).

Retain the finite local coefficient c and the same fixed covariant
subtraction scale as in [the static singlet prescription](singlet-static-curvature.md).
The local difference between optical and covariant subtraction is

```text
E_local=c I+J/(4pi),
I=integral V dy,
J=integral V log N dy.                                   (2)
```

The coefficient of J follows from G_cov,ren-G_opt,ren=log N/(2pi)
and delta E/delta V=G_ren/2. This local statement remains valid for the
spatially varying mass below: the short-distance logarithm is mass
independent. It is not the full massive Weyl action. Bilal and Leduc's
[coincident Green-function derivation](https://arxiv.org/html/1606.01901v2#S3.SS3)
is a predecessor, already reviewed in the singlet calculation. Changing
the subtraction scale also changes the fixed coefficient c; no coefficient
is fitted to the stability calculation.

## Exact product-throat integrals and clock derivatives

Writing L_b=log N_b, direct integration gives

```text
I_0=2a b/epsilon,
J_0=2a[(b/epsilon)(L_b-1)+x A].                           (3)
```

For J substitute r=tan(xy). A primitive of log(1+r^2) is
r log(1+r^2)-2r+2atan r. At the moving endpoint
lambda_0 x sqrt(1+R^2)=N_b; this removes the large variable logarithm.
Let D=partial_phi=x partial_x, with epsilon,lambda_0,a,c fixed. Then

```text
D I_0=-2a epsilon x^2/b,
D^2 I_0=-2a epsilon[2x^2/b+epsilon^2 x^4/b^3],
D J_0=2a x A+L_b D I_0,
D^2 J_0=2a[x A-epsilon x^2/b]+L_b D^2 I_0.               (4)
```

Equations (3)-(4) reproduce the previous proper-cutoff integral I and
extend the calculation to the covariant local logarithm without the
coarse bound proportional to I. For j=1,2,

```text
D^j E_local = a x/4
             +O(a epsilon x^2 [1+|c|+|L_b|]).             (5)
```

The constants are uniform for z<=1/4 and a fixed compact positive x
interval. Indeed |A-pi/2|=asin z<=z/b and b>=sqrt(15)/4;
insert these bounds into (4). The coefficient c multiplies derivatives
of size a epsilon, although its undifferentiated cutoff energy is of
size a/epsilon. Removing an x-independent reference energy does not
change (4), and no such subtraction has been used to cancel a force.

The earlier fixed-r canonical family instead has I_phi=exp(phi)I_0,
and consequently D^2(c I_phi)=c I_0. These are different families of
geometries and joins. One cannot transfer that order-q curvature to the
fixed-proper-time family, or transfer (5) to the canonical retarded kernel
without deriving its moving-mouth vertex and boundary work.

## First-order corrected throat and a cancellation from the field equations

The reference geometry need not be exactly a product for the first-order
local comparison. Reuse the independently tensor-checked solution in
[the reference matching note](casimir-reference-matching.md), with chi=1:

```text
f=1+r^2, S=f+alpha s(r), C=1+alpha c_g(r),
c_g=1+r atan r-U,
s=(V_p+2U)r^2-r^2-r(r^2+3)atan r+log(1+r^2),
ds_4^2=P^2[-x^2 S dt^2+dr^2/S+C dOmega^2].               (6)
```

U,V_p are the fixed invariant pressure parameters of that comparison;
V_p is distinct from the optical mass V. These parameters, a and alpha
are held fixed under this specified clock variation. Their physical
state dependence is not inferred from (6). This is a first-order
background calculation, valid only in its controlled overlap.

Spherical Gauss elimination now gives m_cov^2=a/(P^2 C). Since
dy=dr/(x S), V dy=a x dr/C. The mouth condition is x^2 S(R_alpha)
=epsilon^(-2), so

```text
R_alpha=R+alpha R_1+O(alpha^2), R_1=-s(R)/(2R),
I=2a x integral_0^R_alpha dr/C,
J=2a x integral_0^R_alpha [log(lambda_0 x)+log S/2] dr/C. (7)
```

Define H=log(1+R^2), A=atan R and

```text
K(A)=integral_0^A log(sec^2 theta) dtheta,
B(R)=A-(V_p/2)R-H/(2R),
Q(R)=L_b B+(V_p/2)(R-A)-A H/2+K(A).
```

The exact first-order coefficients are

```text
I=I_0+2a alpha x B+O(alpha^2),
J=J_0+2a alpha x Q+O(alpha^2).                            (8)
```

The cancellation is visible before taking a large-R limit. A primitive is

```text
C_int(R)=integral_0^R c_g(r) dr
        =(1/2-U)R+(1+R^2)atan(R)/2,
R_1-C_int(R)=B(R).                                       (9)
```

The two separate terms in (9) grow as R^2. Their R^2 pieces and all U
dependence cancel. For J, integrate its c_g logarithm by parts to find

```text
Q=L_b B+integral_0^R [r C_int(r)+s(r)/2]/(1+r^2) dr,
r C_int(r)+s(r)/2=(V_p/2)r^2-r atan r+log(1+r^2)/2.        (10)
```

Integrating (10) gives Q in (8). These identities use the corrected
lapse and sphere together with the cutoff shift. Keeping only the sphere
correction would produce a false large remainder. This is a consequence
of the already reproduced linear Einstein-Maxwell solution, not a new
wormhole construction or a general all-orders cancellation theorem.

For R>=sqrt(15), the first two clock derivatives of the terms in (8)
obey, for a numerical constant C independent of the cutoff,

```text
|D^j(2a alpha x B)|+|D^j(2a alpha x Q)|
 <= C a |alpha| x (1+|V_p|)(1+|L_b|+H),     j=1,2.        (11)
```

To see that no power of R is hidden in this bound, use
D R=-R-1/R, D A=-1/R, D H=-2, 0<=K<=pi log 2 and D K=-H/R.
The only apparently growing terms are multiples of xR=b/epsilon, but
D(xR)=-x/R and D^2(xR)=-2x/R-x/R^3. All other terms in xB,xQ and
their first two derivatives are bounded by the right side of (11).
For example D(xB)=x[A+V_p/(2R)-H/R-H/(2R^3)].

The O(alpha^2) symbols in (7)-(8) are at fixed finite cutoff. This note
does not silently promote them to a uniform physical remainder. A full
second-order geometry/state bound is still required. The exact affine
model S=f+alpha s,C=1+alpha c_g can test Taylor convergence, but solving
that prescribed model at finite alpha is not solving the second-order
semiclassical Einstein equations.

## Smooth static joins

The product result is stable under the explicitly controlled static
join class of the earlier proper-clock note. On the optical collar
w=epsilon u, -1<=u<=0, write

```text
N_in/N_b=[cos(zu)-b sin(zu)/z]^(-1),
N_ref/N_b=(1-u)^(-1),
r(z,u)=log(1-u)-log[cos(zu)-b sin(zu)/z],
N_blend=N_ref exp(chi(u) r),
chi=1-10v^3+15v^4-6v^5, v=u+1.                          (12)
```

The sphere radius is fixed in this explicit join test. The exterior
reference and the cutoff lapse are independent of x. The function r
and its first two D derivatives are O(z^2) uniformly on this collar:
it is even analytic in z with r(0,u)=0, and the denominator is bounded
away from zero for |z|<=1/4. The two collars change the local functionals,
relative to N_ref, by exactly

```text
delta I=(2a/epsilon) integral_-1^0 (exp(2chi r)-1)/(1-u)^2 du,
delta J=(2a/epsilon) integral_-1^0
 [(exp(2chi r)-1)(L_b-log(1-u))+exp(2chi r)chi r]/(1-u)^2 du.
```

Therefore their values and first two clock derivatives are
O(a epsilon x^2) and O(a epsilon x^2[1+|L_b|]), respectively.
The same bound controls replacing the corresponding product collar by
the blend, by comparing each to N_ref. These joins show stability of a
local static estimate; they are not a demonstrated Einstein junction.
No claim about nonspherical exterior interactions follows from them.

## Consequence and research decision

In the retained reference family set

```text
g_c^2=q^-2, a=1/(4pi^2 q), alpha=1/(4pi^2 q^3),
lambda_0=1/(16pi q^3), epsilon=q^-2,
delta=log^(-2)q, L_b=-log(16pi q).
```

For fixed c,U,V_p and a compact x interval, (5) is O(q^-1), with
product cutoff error O(q^-3 log q); the known first-order correction is
O(q^-4 log q). Logarithmic growth of the fixed local coefficient retains
these conclusions up to logarithms. Both clock derivatives are o(q delta).
No bound on uncomputed four-dimensional matching coefficients is thereby
created. A coefficient varying with the state requires its own derivative
control.

The preceding turn was progress: it proved the ideal-boundary comparison
and isolated a potentially large local term in the canonical family.
The present calculation resolves that particular size concern for the
fixed-proper-time local comparator, including its known first-order
background correction and a static smoothing class. It leaves the full
finite-mouth optical determinant and its moving-clock response unresolved.
The graph theorem does not supply those by norm-resolvent convergence.

Continue physical long-return damping as provisional primary and matched
finite-duration transfer as fallback. The decisive next task is the
nonlocal finite-mouth determinant and the actual moving-clock vertex,
with angular/compact state and physical exterior matching retained.
Proper-time regulation, covariant subtraction and the source throat
solution are established ingredients. No consequential novelty, physical
stability, finite-disturbance survival or outside shortcut is accepted.

[The verification script](../checks/proper_mouth_local_stress.py),
[receipt](../evidence/proper-mouth-local-stress.json) and
[source supplement](../evidence/proper-mouth-local-source-review.json)
retain exact identities and independent high-precision integration,
clock differentiation and cutoff-shift tests. All Lean/Lake execution
remains reserved for suitable completed endpoints on GCP.

Remaining to-do list: control the nonlocal finite-mouth determinant,
moving-clock vertex and higher-order geometry/state response; close
angular/compact sectors and physical exterior matching; establish the
physical endpoint and novelty, independently verify, certify suitable
mathematics on GCP, and prepare the reproducible manuscript.
