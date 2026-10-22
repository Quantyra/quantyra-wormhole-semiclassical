# Actual bookkeeping scalar residual in weighted H6

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original action split, finite path, physical state/domain,
magnetic flux, mass, subtraction and York return. This estimates the
actual scalar forcing, not the full-action trace coefficient alone.

With g=G lambda, gamma0=4pi beta B^3 and W=t^2, the finite rational
ledger below defines C_tau such that

    ||f_tau||_(H_W6)/sqrt(gamma0) <= C_tau ell^2/g,
    f_tau=lambda^-2 F_ell,Phi=-E_ell/(8pi g).        (1)

Here E_ell is the bookkeeping trace divided by lambda. The receipt
gives its exact rational constant and upward decimal enclosure. The
bound includes actual physical volume and all six ordered covariant
derivatives. It does not give the quantum/TT source H4, TT lift or
inverse, Schur H6 or the positive-transfer certificate. No nonlinear
solution, compatibility zero or stability result is inferred.

## 1. Density conversion and exact quadratic cancellation

The [actual source convention](finite-graph-admissibility.md) is
f_ell=Jvol F_ell(g_ell), I_ell=I0+ell Gamma_rest. Density conversion
precedes the actual York transpose. Since Rcal(tau,0)=(tau,0), the
scalar source is exactly that in (1), with no extra Jvol multiplier.
At ell=1 this is the full retained action; the action along the path
is the original bookkeeping action.

Put z=1/t, d0=d/x, Z0(z)=d0(4z^5-6z^6), and

    Z_ell=(R4(g_ell)+12lambda)/lambda,
    S_ell=log(t_ell/t), z_ell=z exp(-S_ell),
    A_ell=Z0(z_ell), e=epsilon lambda,
    H_ell=e[Boxbar_ell Z_ell+(11/3)Z_ell-(11/72)Z_ell^2
                  +(7/12)C_ell^2+(11/6)J_ell^2]-mu_ell,
    mu_ell=20epsilon M_em(g_ell)/lambda.             (2)

Curvature norms in (2) are divided by lambda^2. The exact reference
potential has U(r)/lambda=12-Z0(B/r). The complete trace therefore is

    E_ell=Z_ell-(1-ell)A_ell+ell H_ell.              (3)

The seed reference equation and the complete tangent trace equation give

    E0=0, dotE0=0,
    dotZ0=-(Z0+H0)+dotA0,
    dotA0=-nu z Z0'(z), nu=(1-z)V.                 (4)

The coordinate Lie term in (4) is retained. Let
N_R=Z_ell-Z0-ell dotZ0 and N_A=A_ell-Z0-ell dotA0. Then exactly

    E_ell=N_R-(1-ell)N_A+ell^2 dotA0
                                      +ell(H_ell-H0). (5)

This identity is used before taking absolute values. It retains the
quadratic ell zero, including at the regular bolt. The constant -12
and the reference potential are not bounded separately on the exterior.

## 2. Fixed-coordinate parameter derivatives with two decay factors

Use the canonical product reference gamma, D=z0 nabla_gamma at each
fixed observation and the positive finite-jet calculus of
[the higher local coefficient proof](finite-higher-local-coefficients.md).
All series coefficients below are divided spatial jets; W=t0^2 stays
outside differentiation. Write T,K for its unweighted and W-weighted
scalar bounds on U=ell u, S_ell and the tangent nu, through ten:

| j | T_j/(100^j) | K_j/(100^j) |
| --- | ---: | ---: |
| 0 | 3e-8 | 1e-4 |
| 1,...,5 | .01 | 10 |
| 6,...,10 | 4e21 | 4e21 |

The same tangent nu allowance follows directly from the retained
radial arrays and their full tensor conversion. For V=nu/(1-z) use
unweighted Vbar and W-weighted KV:

| j | Vbar_j/(100^j) | KV_j/(100^j) |
| --- | ---: | ---: |
| 0 | 1.2e-8 | 4e-5 |
| 1,...,5 | 100 | 4000 |
| 6,...,10 | 2e26 | 2e26 |

These do not assume a lower thermal gap independent of h. On the cold
region the old inverse-(1-z) radial coefficient sum is <2/(1-z),
and a/(1-z)<82. Its low V norm is <.0984 before the tensor factor
914 through five. Its high radial sum is <1.641e18 before the tensor
factor 1.02e7 through ten. On the middle region use inverse sum 1500,
radial nu sum .0006+5e15 and tensor factor 3.12, then W<50.1^2.
The outer inverse sum is <1.07 and the retained remainder/cubic fits
these larger allowances. Thus no inverse power of h enters the table.

Let q=1+Zeta and L=-log(1-Zeta) be the preceding seed majorants. Define
Tplus=T-T0, Vplus=Vbar-Vbar0. Exact parameter derivatives are

    partial_ell S=nu exp(ell V-S),
    partial_ell^2 S=z nu V exp(ell V-2S).          (6)

Consequently their W- and W^2-weighted spatial arrays are bounded by

    S1=1.001 exp(Vplus+Tplus)K,
    S2=1.001 q exp(Vplus+2Tplus) K KV.             (7)

Here every exponential is a finite spatial-jet operation with its
small zeroth value factored off. Large high derivatives are never
exponentiated as metric values. Products of two degree-six high parts
lie beyond the required spatial degree ten.

For the actual metric and its inverse in the fixed canonical frame,
use the following simultaneous parameter arrays:

    M0=3.003 exp(2L+4Tplus),
    M1=2M0(K+S1),
    M2=M0[4(K+S1)^2+2S2].                         (8)

Their meanings are unweighted value jets, W-weighted first parameter
jets and W^2-weighted second parameter jets. The same bounds work for
the inverse, including its second derivative. No moving endpoint
coordinate is differentiated.

## 3. Curvature differentiated in the path parameter

Apply the exact connection, curvature, Ricci, scalar, J and Weyl
formulas in the preceding note to the three arrays (8). For products,
parameter order a uses sum_(b=0)^a binom(a,b)X_b Y_(a-b).
Reference curvature has zero positive parameter derivatives. In
particular let Sc2 be the second parameter array of the normalized
scalar curvature before the C/x conversion. Then

    NR2=6.1 Sc2,
    W^2 |D^j N_R|/j! <= (ell^2/2)(NR2)_j.         (9)

This is Taylor's integral formula applied in fixed coordinates.
It holds through spatial order eight, using metric derivatives through
ten. Let C0,C1 be the unweighted and W-weighted first parameter arrays
of the connection, and J1,CW1 the W-weighted first parameter arrays of
the J/Weyl covariant components, including factor 6.1. The actual
W-weighted component arrays JW0,CW0 and ZW0 are the preceding frozen
coefficient arrays. All indices and metric contraction contacts remain.

## 4. Linear trace and differentiated anomaly differences

For positive polynomial bookkeeping use B from the preceding note and

    Zabs(z)=1.3e-22(4z^5+6z^6),
    jabs(z)=.501z^4+5e-24z^5+4.5e-24z^6,
    cabs(z)=3*.668z^3+6*.501z^4+20e-24z^5+15e-24z^6,
    Box0abs=6.1[z^2 B Zabs''+(z^2B'+2zB)Zabs'],
    J0sqabs=4*(6.1)^2 jabs^2,
    C0sqabs=(4/3)*(6.1)^2 cabs^2,
    em=4e-27/18, emax=2000^-10/180.

The complete seed trace polynomial is bounded by

    Fseed=Zabs+emax[(11/3)Zabs+Box0abs+(11/72)Zabs^2
                               +(7/12)C0sqabs+(11/6)J0sqabs]+em z^4.

Every monomial starts at power four. Substitution z=q bounds its
W^2-weighted divided D jets: each leftover z0^(p-4) is at most one,
while the q exponent remains p. In particular one must not replace
q^p by q^(p-4) when differentiating. With

    Lie=K times (z Zabs') evaluated at q,
    Zd=Fseed(q)+Lie, DZ=Zd+NR2/2,                 (10)

DZ bounds W^2 Delta Z/ell through spatial order eight. This uses the
linear trace equation (4), rather than estimating a generic curvature
difference by its slower t^-2 allowance.

The exact scalar Laplacian in the fixed frame is

    Boxbar_ell f=(C/x) h^-1:(D^2 f-C_ell Df).

Thus a divided-jet bound for W^2 Delta(Boxbar Z)/ell is

    BX=6.1{M0[DZ''+C0 DZ']
                  +M1 Zabs(q)''+(M1 C0+M0 C1)Zabs(q)'}. (11)

Primes on composed formal series in (11) mean spatial series
derivatives, not z derivatives. The second line keeps the variation
of the Laplacian acting on the nonconstant seed trace.

For the curvature-square differences, contractions involve four inverse
metrics for Weyl and two for J. Sufficient W^2-weighted divided arrays
for their parameter derivatives are

    DCsq=2 M0^4 CW0 CW1+4 M0^3 M1 CW0^2,
    DJsq=2 M0^2 JW0 J1+2 M0 M1 JW0^2.             (12)

In each term carry W on two factors and use the remaining weighted
array also as an unweighted upper bound, since W>=1. This is a
two-weight estimate, not multiplication of separately nonintegrable
L2 norms. The magnetic array is

    MU=em*1.001 exp(4Tplus)q^4,
    DMU=4 MU S1.                                 (13)

MU bounds W mu_ell. It uses the physical flux, and DMU bounds
W^2 Delta mu/ell. From (2), (10)-(13), set

    DH=emax[BX+(11/3)DZ+(11/72)(ZW0+Zabs(q))DZ
                                  +(7/12)DCsq+(11/6)DJsq]+DMU.

For the reference remainder use

    RA2=1.001*1.3e-22 sum_(p=5,6) a_p q^p exp(pTplus)
                                       [p^2 S1^2+p S2],
    a5=4, a6=6.

It bounds W^2 partial_ell^2 A_ell. Combining (5) gives the explicit
divided reference-jet bound

    P=NR2/2+RA2/2+Lie+DH,
    W^2 |D^j E_ell|/j! <= ell^2 P_j, j<=6.       (14)

All derivatives in (14) are ordinary ordered tensor derivatives in
the fixed reference. The higher coefficients used by BX are available
because (9)-(10) are established through eight. No eleventh or higher
source derivative is needed here.

## 5. Actual physical derivatives and the noncompact integral

Convert (14) at its fixed endpoint using the scalar instance of the
preceding connection recurrence with C0:

    B(0,k)=P_k,
    B(n+1,k)=(k+1)B(n,k+1)
                          +n sum_(i=0)^k C0_i B(n,k-i),
    H_n=2*(2.5)^n B(n,0), n<=6.                  (15)

It follows that |gradbar^n E_ell|<=ell^2 t^-4 H_n. This is a
pointwise tensor estimate at each observation, so W is not
differentiated. The physical volume ratio Jvol<1.001 and

    dmu0/gamma0=t^2 dt, t>=1+h,
    integral_(1+h)^infinity t^4 t^-8 t^2 dt<1

give the actual weighted Sobolev estimate (1) with

    C_tau=(1001/24000)sum_(n=0)^6 H_n.             (16)

The rational prefactor uses sqrt(1.001)<1.001 and 8pi>24. No finite
outer cutoff, division by the bolt speed or substituted density is
introduced. The original smaller L2 and first-order estimates remain.

## Verification and next step

Freeze the ordinary argument, checker and review before the run.
Verify the reference trace and cancellation identities, parameter
arrays, covariant conversion, physical volume integral and exact
rational enclosure. Independent nonlinear warped-product contours at
60/90 digits test the mixed second-parameter/sixth-spatial curvature
coefficient and detect omission of the radius acceleration. They are
manufactured geometries, not evaluations of the charged Dirac kernel.

This is an actual scalar residual bound, still only one term in the
Schur-source split. The actual quantum/TT H4, higher TT lift/inverse
and transfer constants remain necessary before assessing the full
cancellation target. A large upper bound is not a lower bound on the
residual or a nonexistence result.

Remaining to-do list: actual quantum/TT H4; higher TT lift/inverse and
positive transfer; cancellation, signed compatibility and nonlinear
existence; physical validity, state/boundary, stability, matching,
charge/return, opening, transport, preparation, separate creation and audit.
