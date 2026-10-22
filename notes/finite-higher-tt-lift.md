# Higher physical TT lift at the original finite comparison

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, actual physical measure, regular
bolt, mass and outer-fixed return. This quantifies the lift factor L8
in the [Schur source split](schur-source-derivative-budget.md). It does
not quantify the nonlocal quantum H4 forcing or the TT H4-to-H8 inverse.

## 1. Intrinsic equations valid off the analytic seed

All derivatives below are for the dimensionless metric lambda g_ell;
L2 still uses the physical measure. Set q=gradbar log r_ell and let
P be the smooth projector onto the orbit plane. Write a radial TT
tensor in its orthonormal frame as h=diag(b,-b-2c,c,c). Put
Y=-sqrt(lambda) eta, identifying Y with its radial vector, and

    Psi=divbar Y+2<q,Y>+c.

For the following radial calculation only, primes and A=a'/a,R=r'/r
are dimensionless. The original equations are

    b'=-2(A+R)b-2(A-R)c,
    Y'=A Y-b-c,
    Psi=2(A+2R)Y-b.

Differentiate the last identity and use
Ric_rr=-A'-A^2-2R'-2R^2. This gives the exact intrinsic system

    gradbar Psi=2q Psi-6q c-2Ricbar(Y)-12q<q,Y>,
    divbar Y=Psi-c-2<q,Y>, curlbar Y=0.             (1)

In particular R'=AR-R^2 has NOT been used. Ricbar is the genuine
finite metric's Ricci tensor, normalized by lambda. The corresponding
unweighted formulas are also obtained by removing the weight from the
[seed higher York system](quantitative-higher-york-isomorphism.md);
the direct derivation above establishes their finite-background scope.

The returned disk fields and original physical tensor are exactly

    S_c=c/2+<q,Y>, T_c=-3c/2+Psi/2-3<q,Y>,
    h=(c+2<q,Y>)gbar+(-3c+Psi-6<q,Y>)P
                                             -2 sym gradbar Y. (2)

The sign follows from subtracting L_eta g from the physical tensor
to return to the disk. The trace of (2) vanishes by (1). These are
geometric tensors through the bolt; no derivatives of a singular
radial unit vector or an isolated coefficient A are estimated.

## 2. Low norms start the higher estimate

The [finite classical disk contact](finite-comparison-classical-disk-contact.md)
gives ||b||<2.002||c||. The [actual return inverse](finite-comparison-disk-return-form.md)
gives ||Y||<=sqrt(2)/(1-.005)||b+c|| and
||A Y||<1.008||b+c||. Since |q|<2.25, safe starting bounds are

    ||Y||<4.3||c||, ||Psi||<47||c||, ||h||<4.7||c||. (3)

For Psi use 2||AY||+4||RY||+||b||. For h use
||h||^2=2||b+c||^2+4||c||^2. The inverse is the already proved
outer-fixed one. Its excluded homogeneous vector a is not reintroduced.

## 3. Full coefficient tensors through eight

Use the frozen improved tangent and finite-geometry arrays from
[the low-gap reconstruction](finite-low-gap-reconstruction.md), with
the covariant reference calculus of [the local coefficient proof](finite-higher-local-coefficients.md).
At each fixed observation D=z0 nabla_gamma and Q=1+Zeta. Let T be
the positive-order majorant for the finite log-radius change, L the
seed -log(1-Zeta) majorant, M the metric array and C the connection
array. Reconstruct M,C by the SAME frozen finite-geometry formula;
their retained entries must reconcile exactly before use.

The reference covector jet array for d log r_ell is (L+T)'. Apply
the old rank-one covariant recurrence, with its extra 2.5 conversion
for this first covector index. This gives full raw bounds q_j for
|gradbar^j q| through eight. Retain q_0=2.25 and q_1=40 if smaller.
The latter follows from the previous |R'|<35, AR<12 and R^2<5.1:
the Hessian log r has components (R',AR,R^2,R^2).

The reference curvature array is exactly bounded by

    RR=B''(Q)+2e-10+2C'+2C^2,
    Rdown=M RR.

Apply the rank-four covariant recurrence to 6.1 Rdown through eight,
retaining the known Rm_0<55 and Rm_1<1e16. Ricci contraction costs
two, giving Ric_j<=2 Rm_j; keep Ric_0<15. Order eight needs metric
order ten, already available. No unavailable T_11 is used. These
are unweighted actual tensors. The newer scalar estimate remains
unchanged and is not used as a substitute for tensor curvature.

Finally the warped-product projector identity

    nabla_a P_bc=(I-P)_ab q_c+(I-P)_ac q_b

gives P_0=2 and P_(j+1)=2 sum_i binom(j,i) Q_i q_(j-i), where
Q_0=2 and Q_i=P_i for i>=1. These bounds include all ordered slots.

## 4. A finite-order curl-div estimate with explicit commutators

For a smooth curl-free one-form Y put y_j=||gradbar^j Y|| and
d_j=||gradbar^j divbar Y||. With the curvature bounds above, for
0<=j<=8 a sufficient integrated estimate is

    y_(j+1)<=2d_j+32(j+1)y_j
       +16(j+1) sum_(i=0)^(j-1) binom(j,i+1) Rm_i y_(j-1-i). (4)

The empty sum at j=0 is zero. Here is a tensor count independent
of polar coordinates. Regard grad^j Y as a one-form in Y's original
index, valued in the bundle of covariant j-tensors. The covariant
curl-div integration by parts has a curvature contraction bounded
by 16(j+1)|Rm||grad^j Y|^2. This includes the Ricci contraction
on the one-form index and the curvature action on all j other
indices. Its square-root term is below 32(j+1)y_j since Rm_0<55.

Commuting either exterior differentiation or divergence through j
ordered derivatives leaves curvature products. One commutator acts
on at most j+1 indices. Differentiating it i times and summing its
possible positions gives sum_(l=i)^(j-1) binom(l,i)=binom(j,i+1).
Antisymmetrization and the four-dimensional trace contraction cost
at most four in total for the two operations; even allowing an
additional factor two leaves less than the coefficient 16(j+1)
in (4). The leading divergence term has coefficient at most two.
Taking square roots of the integrated identity and using triangle
inequalities proves (4), with every commutator retained.

Initially use smooth regular fields and cutoffs at large physical
radius. The low inverse (3), the curl-div equation, local regularity
and this estimate applied to cutoff/regularized fields justify the
successive weak derivatives. Passing to the limit proves the same
bound at each order. Smooth radial H8 approximations then extend
the result to the original L2-defined lift. The bolt is an interior
disk point; no artificial boundary term or new return condition occurs.

## 5. Keep each input derivative separate

Let c_j=||gradbar^j c|| and p_j=||gradbar^j Psi||. Initialize
y_0=4.3 c_0, p_0=47 c_0. Define

    x_j=sum_i binom(j,i) q_i y_(j-i),
    d_j=p_j+c_j+2x_j,
    p_(j+1)=sum_i binom(j,i)[2q_i p_(j-i)+6q_i c_(j-i)
                                              +2Ric_i y_(j-i)]
                 +12 sum_i binom(j,i) q_i x_(j-i).

Use (4) for y_(j+1). Carry this recurrence through j=8, retaining
the coefficients of each independent c_0,...,c_8. Formula (2) gives

    h_j<=2c_j+4x_j+2y_(j+1)
        +sum_i binom(j,i) P_i[3c_(j-i)+p_(j-i)+6x_(j-i)]. (5)

Keep the smaller h_0<=4.7c_0. The H8 tensor norm is the sum of
h_0,...,h_8; the H7 return norm is the sum of y_0,...,y_7.
Each is bounded by its largest matrix column sum times
||c||H8=sum_(j=0)^8 c_j. Their maximum is an explicit L8 in

    ||h_c||H8, ||sqrt(lambda)eta_c||H7<=L8 ||c||H8. (6)

The separate matrices are retained for later estimates; replacing
all input orders by one generic coefficient is unnecessary.

Use (6) with the existing C_loc<1e52 and C_E<1e56 to give an
explicit mixed-row constant 381(C_loc+C_E)L8/(8pi g). Round it
upward and compare with the full Schur target. Neither F4 nor C_D,4
is supplied here. For example an equal split of the target gives
sufficient scalar and tensor budgets, not necessary physical limits.
An oversized upper estimate cannot prove nonexistence. Report the
resulting budget before selecting more coefficient refinements.

## Verification scope

Freeze this note, checker and review before execution. Reconcile
the geometry with the preceding frozen arrays, verify all coefficient
and derivative counts, and compute the exact positive matrices and
their column-sum operator bounds. Independently derive the Ricci,
divergence and return identities in four-dimensional warped coordinates.
Use non-seed regular manufactured metrics with decaying fields to
compare the bolt-fixed TT integral and outer-fixed return integral
through eighth radial derivatives at 60/90 digits, including the bolt.
Omitting the Ricci or quadratic-q terms must fail. These are tests
of the identities and bounds, not actual charged quantum evaluations
or formal verification of the continuum argument.

Remaining to-do list: actual finite quantum/TT H4 source and residual
conversion; H4-to-H8 TT inverse and positive transfer; nonlinear mapping,
neighborhood preservation, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
