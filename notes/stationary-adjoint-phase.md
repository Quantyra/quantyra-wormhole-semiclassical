# The actual adjoint phase left after the mass condition

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument for
review before focused verification. **PROGRESS; physical matching open.**

The [global weighted map](stationary-global-weighted-map.md) leaves a
cokernel of dimension k+1 when its mass-extended kernel has dimension k.
This note identifies those tests more precisely. The mass-column moment
is an explicit nonzero multiple of the adjoint TT source coefficient.
After that coefficient vanishes, the actual adjoint space modulo the
physical kernel has exactly one real oscillatory scalar phase.

Thus solvability requires the k physical-kernel moments and one additional
actual scalar-phase moment. The phase and its moment on the candidate
are not evaluated here. No Ward identity, residual gauge freedom or seed
temperature parameter makes that moment vanish automatically.

## 1. Fixed comparison and finite asymptotic data

Use precisely the real fixed-volume Hessian A, physical measure
dmu=gamma(z)z^-4 dz, gamma=gamma0[1+O(z^2)], and mixed spaces X_delta,
Y_delta of the global theorem. Keep the positive prepared temperature,
fixed flux, clock, actual paired Dirac realization and retained scheme.
All statements are at a fixed admissible nonextremal comparison. They
apply in particular to the accepted finite candidate once its stronger
source tail is checked in section 5.

Write delta=7/2, X_m=X_delta+span{e_m},
e_m=chi(z)z^3 u_m, u_m=(1,-2), with chi=1 at infinity and zero near the
bolt. Set K_-=ker(A:X_(-1/2)->Y_(-1/2)) and K_m=ker(A:X_m->Y_delta).
The adjoint weight -1/2 is an auxiliary linear completion, not an
additional physical phase. Let

```
u_0=(1,-1), e_S=(1,0), s_+=3/2+i omega, s_-=3/2-i omega.
```

The finite Mellin construction already proved in the global theorem
gives every real chi in K_- the unique leading coefficients

```
chi = a0 Psi0 + am Psi_m + A_c Psi_c + A_s Psi_s + v,
v in X_(7/2),                                      (1)
```

where Psi0 starts with u_0, Psi_m starts with z^3 u_m, and Psi_c,Psi_s
start with z^(3/2)cos(omega log z)e_S and
z^(3/2)sin(omega log z)e_S. These are asymptotic representatives, not
separately asserted global null fields. Psi0 includes its forced degree
two and possible degree-three logarithmic terms. The scalar representatives
include their forced degree-7/2 terms. Their coefficients retain every
local and globally dependent contact of the previous construction.
Psi_m-e_m lies in X_(7/2). Leading coefficients do not depend on changes
of cutoff or of the fast representatives.

The [normal symbol](stationary-disk-and-normal-inverse.md) is

```
H(s)=C(3-s)^T diag(J(s),N_T(s)) C(s),
C(s)=[[1,(s+3)/(3(s+1))],[0,2]],
q=s(s-3), d0=1/(480pi^2), epsilon=G/(60pi),
a=1+11epsilon lambda/3,
J(s)=3[a+epsilon lambda q](q-4)/(4pi G lambda),
N_T(s)=d0 q(q+2)[F_M(s)+kappa_M]-A_T q,
F_M(3)=F_M(0)=1-gamma_E.                            (2)
```

Here lambda=lambda_e. The physical-volume normal family is lambda^2 H.
The fixed factor gamma0 includes whichever angular/period integrations
are used in dmu. Keeping gamma0 explicitly avoids silently normalizing
away that factor when comparing boundary coefficients with bulk moments.

## 2. A cutoff Green identity for the actual logarithmic operator

Use the real bilinear pairing <u,f>=integral u dot f dmu. For a smooth
step c_T in x=log z, zero toward infinity and one toward the core,
send its transition to x=T with T tending to minus infinity. Compact
testing of A chi=0, and the actual fixed-volume symmetry, give

```
<chi,A e_m> = -lim_(T->-infinity) <chi,[A,c_T]e_m>. (3)
```

The left integral converges: A e_m=O(z^4 log^N z), while the least
decaying term of chi is constant. Smoothness at the bolt introduces no
other boundary contribution. The same cutoff procedure defines the
antisymmetric boundary pairing of two null fields.

Here is the normal computation, including its applicability to the
nonlocal term. Conjugate the leading physical measure by z^(-3/2).
The normal operator becomes a translation-invariant matrix operator
in x. Write its distribution kernel as K(r) and its exponential symbol
as B(beta)=lambda^2 H(3/2+beta). For complementary powers with root
vectors b,d, the commutator integral uses

```
integral_R [c_T(x-r)-c_T(x)] dx = -r,
integral_R K(r)exp(-beta r)(-r)dr = B'(beta).       (4)
```

Thus the pairing with the second exponent s is the corresponding
matrix entry of lambda^2 H'(s). This expression includes all contact
distributions. It does not replace a finite part by an absolute kernel.
One may justify its diagonal part without a distributional moment by
writing

```
[B(D+beta)-B(beta)]c_T
 = M_beta(D)c_T',
M_beta(i t)=[B(beta+i t)-B(beta)]/(i t),
M_beta(0)=B'(beta).                                (5)
```

The derivative c_T' is smooth and compactly supported with integral one.
The logarithmic symbol has polynomial/logarithmic growth, so Fourier
pairing with this smooth derivative is legitimate. Exponential decay
away from the diagonal, at every root strictly inside -1<Re s<4,
gives the convergent integral and its zero-frequency value in (5).
This also fixes the sign in (3): it is the negative normal derivative.

To pass to A, use its two-power dyadic improvement and finite asymptotic
expansion. A pairing of terms of degrees p,q scales as
exp[(p+q-3)T]; an actual-minus-normal term gains another exp(2T), up
to finite powers of |T|. For source/mass the sum is three. Every forced
source correction paired with mass has larger real sum. For two
source-free fields, the smallest sum is three, from the scalar pair;
their forced corrections and all mass pairings have larger sum.
Therefore all those errors tend to zero. The X_(7/2) remainders give
vanishing terms by weighted Cauchy--Schwarz and the bounded cutoff
commutators. For the source/mass pairing the worst sum with a remainder
is 7/2; for source-free pairings it is five.

The commutator is not confined to the transition annulus. Its separated
tails are controlled by the same positive log-coordinate decay rates
4-Re s and 1+Re s. Split off a fixed collar and then sum the geometric
dyadic tails; compact-to-end contributions vanish with the cutoff.
For example, in the mass pairing the actual-minus-normal crossing
z<epsilon<y is bounded by
integral_0^epsilon dz integral_epsilon^epsilon0 dy (z+y)^-3 y^3,
up to finite logarithms, hence tends to zero as O(epsilon) times such
logarithms. The reversed crossing is O(epsilon^2 log^N epsilon).
For two scalar leading powers both crossings are
O(epsilon^2 log^N epsilon). Thus tail convergence does not require
the annular two-power estimate to hold uniformly on a fixed far patch.
This justifies the limit for the actual nonlocal operator, including
its retained state-dependent terms.

[Witt, section 6, Theorem 6.1](https://arxiv.org/pdf/math/0311028)
expresses the Green form of cone differential operators through Mellin
residues and asymptotic coefficients. That theorem concerns differential
operators. Equations (3)-(5) and the actual kernel estimates above supply
the specific extension needed here for the logarithmic response; a
classical differential theorem is not assumed to cover it automatically.

## 3. Exact mass moment and its nonzero coefficient

The two TT vectors obey C(0)u_0=(0,-2) and C(3)u_m=(0,-4).
Differentiating the congruence in (2), the derivatives of C do not
contribute between these root vectors. Consequently

```
u_0^T H'(3)u_m=8 N_T'(3),
N_T'(3)=3[2d0(1-gamma_E+kappa_M)-A_T],
<chi,A e_m>=B_m a0(chi),
B_m=24 gamma0 lambda^2
          [A_T-2d0(1-gamma_E+kappa_M)] > 0.         (6)
```

This is an actual bulk moment evaluated from boundary data. It is
unchanged by replacing e_m with another cutoff extension. With
L=A_T/d0>30 k^10, k=2000, and 39/5<kappa_M<54/5, the relative loop
correction obeys

```
0 < 2(1-gamma_E+kappa_M)/L < 24/L < 4/(5k^10).    (7)
```

The established bounds 0<gamma_E<1 suffice. Thus B_m is quantitatively
nonzero; at k=2000 the relative correction is below 8e-34. This is not
a bound for the global inverse.

The mass-extended cokernel is now exactly

```
C_m={chi in K_- : a0(chi)=0}.                      (8)
```

Nonzero B_m does not prove that the actual K_- contains a field with
a0 nonzero. In particular this does not assert mass transversality or
that adding mass always removes one actual cokernel direction. The
index increases by one regardless of that rank question.

## 4. One actual oscillatory phase modulo the physical kernel

Define O:C_m -> R^2 by O(chi)=(A_c,A_s). Its kernel is precisely K_m.
Indeed a field in K_m has only mass and fast asymptotics, solves A chi=0
and lies in (8). Conversely a field in C_m with O=0 has, by (1), only
the correlated mass cubic plus an X_(7/2) remainder, hence lies in K_m.
This proves an equality of actual global null spaces, not just a count
of formal normal solutions.

If k=dim K_m, the global index theorem gives dim C_m=k+1. Therefore

```
ker O=K_m, dim O(C_m)=1,
C_m/K_m is one nonzero real scalar phase line.     (9)
```

There is also a direct Green-form check on its geometry. Since
C(s)e_S=e_S, the scalar root derivative is

```
J'(s_+)=-i K, J'(s_-)=+i K,
K=3 omega(a+4epsilon lambda)/(2pi G lambda)
 =omega(omega^2+25/4)/(40pi^2)>0.                  (10)
```

For two real source-free asymptotic fields the oriented outer Green
form is gamma0 lambda^2 K/2 times
A_c,1 A_s,2-A_s,1 A_c,2. Actual null fields regular at the bolt have
zero total Green form, so their scalar-amplitude image is isotropic.
Equation (9) shows that it is a line, not the zero subspace.

Choose a basis k_1,...,k_k for K_m and one actual chi_* in C_m whose
scalar amplitude lies on that line with A_c^2+A_s^2=1. Its sign and
addition by K_m are conventional. The exact source condition becomes

```
A_m u=f is solvable iff
<k_j,f>=0 for all j and <chi_*,f>=0.                (11)
```

Once the kernel moments vanish, the last test is independent of the
choice of chi_* modulo K_m. Its phase angle is determined by the actual
interior problem, including the full Dirac response. It is not selected
by substituting a cosine of the vacuum frequency into the integral.
This theorem does not determine k, that angle or the moment's value.

## 5. The actual candidate is an admissible source, not a zero moment

Let f_app=-F(g_app) be the two-component action gradient at the accepted
finite candidate, expressed per its fixed comparison volume. A mere
bound F=O(z^2 log^N z), sufficient earlier for a density-contact estimate,
would not put this source in Y_(7/2). Use the stronger actual tail.

In [the complete residual](stationary-centered-complete-residual.md),
equation (5) controls the combined physical profiles
H_rho=r^2 rho2 and H_N=r^2(rho2+p2) by finite logarithms after the
correlated mass cancellation. Since rho2 and p2 are sphere-integrated,
their four-dimensional pressures are O(r^-4 log^N r). The seed's full
combined source has the same outer cancellation in
[the source proof](stationary-centered-exterior-source.md). The reference
potential's difference from its fixed outer vacuum also has that
falloff. These are pointwise source statements, not deductions from a
small integrated mass or lapse norm.

The candidate metric has the five differentiated fast jets of
[the fast phase](stationary-fast-phase-variation.md). The Einstein
tensor minus its centered vacuum, apart from the fixed magnetic term,
is O(r^-4 log^N r): the constant and quadratic asymptotic data match
AdS, and the mass cubic solves the vacuum leading equations. The exact
local full trace has this same decay, using the retained anomaly and
its differentiated curvature jets. The angular residual follows by
subtracting the orbit trace from the full trace. Hence both disk
gradients obey

```
f_app=O(z^4 log^N z), f_app in Y_(7/2).             (12)
```

At a fixed positive-temperature bolt they are regular. Thus every
moment in (11) is defined for the actual candidate. This establishes
admissibility; it does not prove any of the moments vanish. The earlier
equation-defect bounds remain equation-defect bounds.

## 6. Why covariance and the thermal family do not solve (11)

There is an explicit Ward witness. In proper radius rho put
U=8piG f_T, H=4piG(f_S-f_T), D=E_rho-E_varphi and
W=a_orb^2 r^2. The spherical Ward identity is equivalent to

```
(W D)'=W[-U'-2(r'/r)U+4(r'/r)H],
E_rho=(U+D)/2, E_varphi=(U-D)/2, E_theta=H.         (13)
```

For any smooth reduced pair compactly supported away from the two
ends, integrate (13) from the bolt with zero constant. Its tensor is
smooth at the bolt and conserved. Beyond the support D=C/W=O(r^-4),
consistent with the fast tensor falloffs. Thus Ward conservation allows
every such reduced pair, including pairs with nonzero moment on chi_*.
Such a pair exists by pairing a bump with the nonzero smooth field chi_*.
This is a symmetry-level witness; it does not replace the Dirac source
or claim that arbitrary conserved tensors arise from that determinant.

An invariant radial diffeomorphism eta(xi)partial_xi preserves the disk
isothermal gauge only when eta'=eta/xi, so eta=c xi. Fixing the outer
endpoint xi=1 forces c=0. If retained, that dilation moves the boundary
and gives a z^-1 metric tangent, outside even X_(-1/2). The actual
scalar phase in (9) is not a surviving gauge direction.

Nor is the seed's thermal parameter a free physical kernel amplitude.
With B and lambda fixed, its disk constants vary as
(delta Phi_0,delta V_0)=(delta kappa/kappa,-delta kappa/kappa), the TT
source root u_0. The seed temperature map is strictly monotone, so a
fixed prepared kappa chooses one member. Changing mass while holding
its old period instead makes a conical bolt; a time rescaling changes
the fixed outer clock. The one allowed corrected mass column is already
included in X_m. No new adjustable physical datum is supplied here.

## Verification contract and next decision

Focused verification may check the congruence/root derivatives,
normalizations and signs, the rational nonzero bound, scalar Green
determinant, Ward reconstruction, gauge equation, thermal source
coefficient and weighted source integrability. It cannot compute the
actual phase or certify the nonlocal cutoff argument and global
quotient theorem. Those remain ordinary proofs. Review the argument
before any scientific execution and retain exact input bytes and runs.

The next decisive task is to determine the actual K_m and phase line,
then evaluate (11) on f_app and obtain a quantitative compatible inverse.
Neither a generic Fredholm theorem nor covariance answers that question.
If controlled spatial reduction is selected instead, its physical
observables and branch error still require proof in the retained model.

Remaining to-do list: determine the actual physical kernel and adjoint
phase, evaluate the candidate moments and quantify the compatible inverse;
close the nonlinear physical logarithmic phase or justify controlled
spatial reduction; match exterior/state and charge, transfer S138, and
complete driven stress, radius, controller, capture, tail, spatial area
and S074-S080.
