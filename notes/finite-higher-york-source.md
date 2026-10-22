# Higher physical York source conversion at the finite comparison

2026-10-22, E003/S167. Ordinary argument before supporting execution.
Keep the original finite comparison, full action, physical volume,
paired realization, mass and outer-fixed return. This quantifies the
H4 source conversion in the [Schur source split](schur-source-derivative-budget.md).
The finite-background quantum source norm and higher TT inverse
remain separate unknowns.

## 1. Convert the density before taking the physical transpose

Let f=(f_S,f_T) be the original fixed-density residual and set

    (s,t)=lambda^-2 Jvol^-1(f_S,f_T),
    f_tau=s, f_c=R_c^*(s,t).

The transpose uses the actual physical measure. All following
derivatives use gbar=lambda g_ell, while L2 uses dmu_ell. For radial
identities put A=a'/a, R=r'/r and mu=a r^2; primes are dimensionless
proper-radius derivatives. Let H be the existing outer-fixed inverse
of partial_rho-A, acting on scalar radial components:

    H d=-a(rho) integral_rho^infinity d(x)/a(x) dx.

The [actual return proof](finite-comparison-disk-return-form.md) gives
the bounded maps

    ||H||<1.422, ||A H||<1.008, ||R H||<3.2.

These bound products, not the singular coefficient A by itself.
For a TT amplitude c, its radial component b and disk return are

    b'+2(A+R)b=-2(A-R)c,
    eta=H(b+c),
    S_c=c/2-R eta, T_c=-(b+3c)/2-(A-R)eta.

## 2. The actual adjoint, with both integral endpoints retained

Put U=a^2 r^2. Integration in the physical pairing gives

    v=U^-1 integral_0^rho U[R s+(A-R)t] dx,
    w=v-t/2, Z=H w,
    f_c=(s-3t)/2+v+2(A-R)Z.                 (1)

For example the formal adjoint of H is
H^*k=-U^-1 integral_0^rho U k dx. The TT integral's adjoint is
2(A-R)H, giving (1). The lower endpoint in v is the regular bolt;
the upper endpoint in Z is the original outer return. Neither is
an arbitrary integration constant.

The potentially singular expression defining v is rigorously the
bounded operator

    v=-(RH)^*s-[(A-R)H]^*t.

Thus it is defined on L2 without asserting A t belongs to L2.
For smooth even inputs, v(0)=t(0)/2 and w=O(rho^2). The return
Z is a smooth odd radial vector. The integrals and boundary pairings
first hold for regular compact inputs and retain their induced tails;
the bounded maps extend them by density.

Define G=s/2-t and Psi=divbar Z. Since Z'=AZ+w,

    Psi=w+2(A+R)Z, f_c=G+Psi-4<q,Z>,
    q=gradbar log r.

Differentiating w gives w'+2(A+R)w=2RG-t'/2. Hence the complete
intrinsic adjoint system is

    gradbar Psi=kappa Z+2qG-gradbar t/2,
    divbar Z=Psi, curlbar Z=0,
    kappa=-tr_P Ricbar-2|q|^2.                  (2)

Here P is the smooth orbit-plane projector. Indeed
kappa=2[A'+A^2+R'+AR] in radial components. The genuine finite
Ricci tensor is retained; no seed-only identity for R' is imposed.
Every object in (2) is a regular tensor at the bolt.

To expose the principal source term, put Upsilon=Psi+t/2. Then

    gradbar Upsilon=kappa Z+2qG,
    divbar Z=Upsilon-t/2, curlbar Z=0,
    f_c=s/2-3t/2+Upsilon-4<q,Z>.               (3)

## 3. Quantitative low data and ordered higher derivatives

The preceding bounded products give

    ||v||<=3.2||s||+4.21||t||,
    ||w||<=3.2||s||+4.71||t||,
    ||Z||<=4.6||s||+6.7||t||,
    ||Psi||<=31||s||+45||t||,
    ||Upsilon||<=31||s||+46||t||.              (4)

For the last estimate use
||Psi||<=[1+2(1.008+3.2)]||w||. Independently the physical return
has norm below 4.6 times the tensor norm, and the TT lift has norm
below 4.7 times ||c||. Thus ||R_c||<22, and duality gives the
smaller output start ||f_c||<=22(||s||+||t||).

Use the full raw coefficient arrays q_j,Rm_j,Ric_j,P_j from the
[finite higher lift](finite-higher-tt-lift.md). They describe the same
finite metric and physical derivatives. A sufficient scalar array is

    kappa_j=sum_i binom(j,i) P_i Ric_(j-i)
                              +2 sum_i binom(j,i) q_i q_(j-i), j<=3.

Keep kappa_0<32: |tr_P Ric|<=sqrt(2)|Ric|, |Ric|<15 and
|q|^2<5.1 imply sqrt(2)*15+10.2<32.

For s_j=||gradbar^j s||, t_j=||gradbar^j t||, define
G_j=s_j/2+t_j, z_j=||gradbar^j Z||, u_j=||gradbar^j Upsilon||.
Starting from (4), use for j=0,1,2,3

    u_(j+1)<=sum_i binom(j,i) kappa_i z_(j-i)
                                  +2 sum_i binom(j,i) q_i G_(j-i),
    z_(j+1)<=2u_j+t_j+32(j+1)z_j
       +16(j+1) sum_(i=0)^(j-1) binom(j,i+1) Rm_i z_(j-1-i).

The second inequality is exactly the ordered curl-div estimate of
the finite lift proof. Its curvature commutators and every tensor
slot remain. Equation (3) gives the output rows

    c_j<=s_j/2+3t_j/2+u_j+4 sum_i binom(j,i) q_i z_(j-i), j<=4,

with c_0 replaced by the smaller 22(s_0+t_0). Carry s_0,...,s_4
and t_0,...,t_4 separately. The maximum column sum of this positive
five-by-ten matrix defines C_Y,4 in

    ||R_c^*(s,t)||H4<=C_Y,4 (||s||H4+||t||H4).   (5)

This is an estimate for the actual transpose, not an inference
from a forward H4 map. On regular compact inputs, local elliptic
regularity for (2) and the cutoff curl-div estimate justify each
derivative successively. Applying the estimates to differences and
using smooth radial H4 approximation extends (5) to the L2-defined
adjoint. The original endpoints and physical measure are unchanged.

## 4. Explicit density matrix and the remaining shape forcing

For users starting with the fixed-density components, keep their
multiplication by Jvol^-1 before (5). The exact volume ratio is

    Jvol=exp(2U+4S), U=ell u, S=log(t_ell/t).

The divided reference-jet majorant for its inverse is
1.001 exp(6Tplus). Convert it by the existing rank-zero covariant
recurrence through four, keeping the smaller zeroth bound 1.001;
call the raw derivative bounds d_j. The product matrix is

    D_(n,l)=binom(n,l)d_(n-l), 0<=l<=n<=4.

Composing both input blocks of (5) with D gives a numerical H4
conversion directly from lambda^-2 f in actual physical norms.
It does not commute the density multiplier through the transpose.
No new higher density norm of the actual residual is presumed.

Alternatively isolate the complete physical shape covector

    b_src=(t-s/3)/2, so t=s/3+2 b_src,
    F_b4=||b_src||H4/sqrt(gamma0).

This follows from theta=S+T/3 and b=2T. It is the complete shape
source, including classical, local and quantum contributions; it
is not identified with the massive determinant alone. Substitute it
in the exact system BEFORE taking absolute values:

    gradbar Upsilon=kappa Z+q(s/3-4 b_src),
    divbar Z=Upsilon-s/6-b_src,
    f_c=-3 b_src+Upsilon-4<q,Z>.                (6)

The direct scalar principal term cancels exactly. Apply the same
recurrence with G_j=s_j/6+2(b_src)_j, divergence allowance
u_j+s_j/6+(b_src)_j, and output allowance
3(b_src)_j+u_j+4(q*z)_j. Its starts are
z_0=(41/6)s_0+13.4(b_src)_0 and
u_0=(139/3)s_0+92(b_src)_0; cap its output value by
(88/3)s_0+44(b_src)_0. In this new five-by-ten matrix the scalar
order-four column is exactly zero: the scalar part gains one
derivative. If a_l,b_l are its two block column sums, set

    C_b4=max_l b_l,
    C_sigma=sum_(l=0)^3 a_l sigma_l,
    sigma_l=(1001/24000) H_l.

H_l are the actual scalar residual's physical derivative bounds
from the [radius/source result](finite-radius-operator-source.md).
Since W>=1 they also bound the unweighted scalar derivatives. Keep
sigma_0<=6e-12 from the earlier complete low residual conversion.
The actual comparison source therefore obeys the explicit split

    F4=||f_c||H4/sqrt(gamma0)
                     <=C_sigma ell^2/g+C_b4 F_b4.          (7)

Report these constants and insert (7) into the unchanged complete
source budget from the [refreshed mixed ledger](finite-refreshed-local-mixed.md):

    N6(p)<=800/g [Cstart ell^2
                         +mhi C_D,4(C_sigma ell^2/g+C_b4 F_b4)].

F_b4 and C_D,4 are still unevaluated. The known low inverse factors,
Cstart<1e29 and mhi<1e13 remain unchanged. At ell=1 the second
equal-split sufficient budget becomes

    (R0 R2 R4) C_D,4(C_sigma/g+C_b4 F_b4)<=6.25e40 g epsilon_c.

It is not a target pass, a necessary physical restriction or an
impossibility statement. At ell=0 the original residual is zero.

## 5. Independent verification

Freeze note, checker and source review. Reconcile all pinned inputs,
the low operator products, curvature identities, complete ordered
matrices, density products and scalar/shape allocation. Independently
derive the four-dimensional warped Ricci contraction and test the
actual adjoint integrals through derivative order four on non-seed
metrics a=rho exp(alpha rho^2), r=exp(beta rho^2), including the bolt.
Use decaying even sources and a decaying odd return at 60/90 digits.
Verify the physical integral pairing against the independently
manufactured forward TT lift; a flat-density substitution must fail.
Retain the outer boundary and test omission of the curvature term.
These checks support an ordinary proof; they do not evaluate the
charged quantum kernel or establish nonlinear existence or stability.

Remaining to-do list: actual finite complete quantum/shape H4 source;
H4-to-H8 TT inverse and positive transfer; nonlinear mapping,
neighborhood preservation, signed compatibility and existence;
physical validity, state/boundary, matching, stability, charge/return,
opening, transport, preparation, separate creation and audit.
