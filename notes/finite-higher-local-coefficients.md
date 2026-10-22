# Higher local coefficients of the actual finite mixed row

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, action, state, paired domain,
mass, subtraction and outer-fixed return. This quantifies C_loc and C_E
in [the direct Schur contract](schur-source-derivative-budget.md):

    C_loc=max_(j<=2,k<=6) ||W gradbar^k A_j||_infinity,
    C_E=max_(k<=7) ||W gradbar^k Ebar||_infinity,
    Tbar h=sum_(j=0)^2 A_j dot gradbar^j h,
    Bbar c=-(8pi g)^-1[Tbar h_c-divbar(Ebar eta_bar_c)].

The positive rational algorithm below defines explicit finite constants;
the receipt reports their decimal enclosures and upward powers of ten.
It uses original metric order eleven and requires no new proper-endpoint
pullback. The consequent mixed-row bound is conditional only on the
still-unquantified lift constant L8:

    ||Bbar c||_(H_W6)<=381(C_loc+C_E)L8 ||c||_H8/(8pi g).

This is a coefficient result, not the H6 norm of the actual Schur source.
The quantitative TT inverse, lift, quantum H4 and scalar residual H6
remain necessary. Large upper constants alone do not certify the target
R0 R2 R4 N6(p)<=1e57 epsilon_c or prove that it is impossible.

## 1. A regular covariant reference through the bolt

Write z=1/t and let gamma be the canonical orbit metric in A units
times the round sphere of radius sqrt(C). Its parallel orbit/sphere
projectors are P,Q. With b=calB/C, outward canonical proper radius gives

    z'=-sqrt(b), (sqrt(b))'=-b_z/2,
    Hess_gamma z=(b_z/2)P.                         (1)

The angular component follows from N proportional to sqrt(b).
Equation (1) extends through the bolt by the regular disk geometry.
There is no division by a vanishing radial speed. Sphere derivatives
of z vanish because gamma is a product. Its curvature consists of the
orbit curvature -b_zz/2 and sphere curvature 1/C.

Use the coefficientwise positive polynomial

    B(q)=1/5.9+1e-10 q^2+.668 q^3+.501 q^4
                               +2e-24 q^5+1e-24 q^6.

It majorizes the actual b coefficients; its cubic bounds a negative
actual coefficient. Dropping that negative term gives b<1 on 0<=z<=1,
so |dz|<1. Define a formal series Zeta with zero constant term by

    Zeta'(0)=1, Zeta''=(3/4)B'(1+Zeta).           (2)

All coefficients needed through eleven are determined recursively.
Since |P|=sqrt(2)<3/2, (1), the ordered tensor product rule and the
finite chain rule show |nabla_gamma^j z|/j! <= [q^j]Zeta. The chain
rule estimate is valid for the full ordered covariant tensor, not just
the symmetric radial derivative. Equation (2) is a positive majorant,
not an equation imposed on the physical z.

At an observation z0 fix D=z0 nabla_gamma. Multipliers z0 are constants
while taking these local jets. Since z0<=1, the jets of z/z0-1 under D
are bounded by Zeta. Thus the jets of phi-phi0=log(t/t0) are bounded by

    L=-log(1-Zeta).

The reference curvature coefficients under D, including its two length
dimensions, are bounded by B''(1+Zeta)+2e-10. No analytic continuation
of the actual comparison or nonlocal source is being assumed.

## 2. Small low jets and large high jets remain separate

Let U=ell u and S=log(t_ell/t). The physical metric is

    g_ell=A^2 t0^2 h,
    h=e^(2(phi-phi0+U+S))P+e^(2(phi-phi0+S))Q.    (3)

The [regional metric theorem](seed-nine-metric-eleven.md) and its smaller
low-order inputs give the following global bounds for either scalar,
with coefficients normalized by j!100^j and D as above:

| j | unweighted allowance T_j | W-weighted allowance K_j |
| --- | ---: | ---: |
| 0 | 3e-8 | 1e-4 |
| 1,...,5 | .01 | 10 |
| 6,...,10 | 4e21 | 4e21 |
| 11 | 4e24 | 4e24 |

Here W=t0^2 is placed outside each differentiated coefficient, never
differentiated. These enlarged intermediate allowances do not replace
the smaller original contracts. For the high cold terms use a<=.00331,
t<=1.001; for middle terms t<50.1; for outer remainders use z<Y<3z,
Y<=.1, Y^2L(Y)<=.0331. The cubic is retained. Its normalized full
tensor coefficients have reserve 2Y^3 through eleven, using the existing
outer connection recurrence; both cubic coefficients are below 1e-11.

For clarity the low covariant entries do not follow just from chart
existence. The original radial low jets through five are 1e-4 for u,nu.
The middle finite-log coefficient norm is below .0016 and its tensor
amplification through five below 1.26. The cold finite-log coefficient
norm is below .00073a and its near tensor amplification through five
below 914 (away is smaller). These imply unweighted .01 and weighted
10, including the bolt. The outer low remainder is 1e-7 Y^(4-j)L and
has the retained cubic. Values use the independent 3e-8 and outer t^-2
bounds. In particular large sixth derivatives are not used to bound
the first derivative or the metric inverse at a point.

Below T and K denote the table's formal polynomials with coefficient
T_j 100^j and K_j 100^j. Put Tplus=T-T0 and

    F=exp(2L+4Tplus), M=3.003 F, DeltaM=12.012 F K. (4)

M bounds both h and h^-1 component tensor jets in the gamma frame.
The sum of the two projector norms is below three, and exp(4T0)<1.001.
The weighted differences from the seed are bounded by DeltaM. This
follows by integrating the exponential derivative along the fixed
scalar segment; the weight is carried by one field difference. The
same estimate works for the inverse. Products are finite coefficient
convolutions. Terms of degree at least six still square to zero through
eleven; no large high allowance is exponentiated as a zeroth-order value.

## 3. Curvature and its weighted differences

All quantities in this section are formal positive coefficient bounds
in the D jets. Prime denotes the derivative of a formal polynomial.
The connection difference from gamma and its weighted seed difference
obey

    C=(3/2) M M', DeltaC=(3/2)(DeltaM M'+M DeltaM').

This is the exact connection difference formula with its three terms.
The (1,3) curvature, Ricci, scalar, tracefree Ricci and Weyl bounds are

    R=B''(1+Zeta)+2e-10+2C'+2C^2,
    Ric=2R, Sc=M Ric, J=Ric+Sc M/4,
    Rdown=M R, Weyl=Rdown+2M Ric+(2/3)Sc M^2.     (5)

Ricci contraction costs two in dimension four. The last line follows
from the four terms of the Kulkarni products; using triangle bounds
here is sufficient. Corresponding weighted differences follow from
DeltaR=2DeltaC'+4C DeltaC and the telescoping product rule
Delta(XY)<=DeltaX Y+X DeltaY. Reference curvature itself is unchanged.
Keep coefficients only through the order available: curvature through
nine, connection through ten and metric through eleven.

The seed weighted curvature components can be estimated before the
physical tensor conversion. Define

    pJ=.501(1+Zeta)^2+5e-24(1+Zeta)^3
                                      +4.5e-24(1+Zeta)^4,
    pC=3*.668/(1-Zeta)+6*.501
                              +20e-24(1+Zeta)+15e-24(1+Zeta)^2,
    pZ=1.3e-22[4(1+Zeta)^5+6(1+Zeta)^6].

Then initial weighted jet arrays for J/lambda, C4/lambda and
Zbar=(R4+12lambda)/lambda are, respectively,

    JW=(61/10)DeltaJ+13pJ,
    CW=(61/10)DeltaWeyl+25pC,
    ZW=(61/10)DeltaSc+pZ.                        (6)

To see the seed factors directly, the normalized sectional curvatures
are -(C/x)(b-zb'+z^2b''/2), -(C/x)(b-zb'/2), and
(C/x)(z^2/C-b). Tracefree Ricci has opposite orbit/sphere eigenvalues
proportional to -b4 z^4-(5/2)b5 z^5-(9/2)b6 z^6; Weyl is proportional
to -3b3 z^3-6b4 z^4-10b5 z^5-15b6 z^6. Their tensor factors have norms
below two and four. After W, their conformal lowering factors leave
powers q^(p-2) for J and q^(p-4) for Weyl, q=z/z0. All remaining
powers of z0 are at most one. This explains the sole inverse factor
in pC. The scalar seed tail is exactly (d/x)(4z^5-6z^6), yielding pZ.
The constant -12 is canceled before estimating a weighted scalar jet.

For an initial covariant rank-r tensor array H, convert all its D jets
to actual covariant derivatives using

    H_(0,k)=H_k,
    H_(n+1,k)=(k+1)H_(n,k+1)
                       +(r+n)sum_(i=0)^k C_i H_(n,k-i). (7)

The actual normalized physical norm is below 2*(2.5)^n H_(n,0).
Indeed sqrt(C/x)<2.5 and the pointwise conformal frame changes cost
less than two through the largest rank used. Apply (7) with ranks
2,4,0 to (6), and rank four to (61/10)Rdown for an unweighted Rm array.
Every connection slot is retained, including its zeroth coefficient.

The magnetic term mu=20epsilon M_em/lambda has weighted initial array

    (4e-27/18)*1.001 exp(4Tplus)(1+Zeta)^4.

Its conversion uses rank zero in (7). This keeps the fixed physical
flux; the reference metric charge is not substituted for it.

## 4. The trace and TT local coefficients

Now all arrays denote RAW actual gradbar derivative bounds. Use the
smaller retained values/derivatives whenever available: weighted J,C4
values 500, first derivatives 3e17, J second derivative 1e25; weighted
Zbar value 1.3e-10, first .0017 and second 1e20. Unweighted values are
J<8, C4<50, Rm<55, Zbar<3.1e-12, with |gradbar Rm|<1e16. Parallel
orthogonal projection also gives unweighted J derivatives <=2 Rm
derivatives and C4 derivatives <=Rm derivatives. All larger arrays
remain available where these smaller bounds do not apply.

For raw arrays write (X*Y)_k=sum_i binom(k,i)X_i Y_(k-i), and
(shift_m X)_k=X_(k+m). Let e=(2000^-10)/180, an upper bound for
epsilon lambda, and denote weighted arrays by a W suffix. The full
trace identity gives the explicit bound through k=7

    E_k=(1+11e/3)ZW_k+2e ZW_(k+2)
       +(11e/72)(ZU*ZW)_k+(7e/12)(CU*CW)_k
       +(11e/6)(JU*JW)_k+muW_k.                  (8)

The factor two in the Laplacian is the four-dimensional contraction
norm. Its seventh derivative uses metric eleven exactly.

For a genuine TT tensor, delta R=-J:h and
delta Box R=-Box(J:h)-h:Hess R. Differentiating the curvature
contractions gives sufficient operator-coefficient arrays

    A2=e[(17/3)JW+(7/2)CW],
    A1=4e shift_1 JW,
    A0=(1+11e/3)JW+(11e/36)(ZU*JW)
       +e[2 shift_2 JW+shift_2 ZW
                   +(35/3)(CW*RmU)+(110/3)(JW*RmU)]+2muW. (9)

In deriving (9), the TT Ricci variation has principal bound
|nabla^2 h| and lower bound (|Ric|+2|Rm|)|h|. The Weyl-square
variation has principal bound 6|C4||nabla^2h| and lower bound
20|C4||Rm||h|. The J-square variation is bounded by
2|J||nabla^2h|+20|J||Rm||h|. These larger constants include inverse
metric variations and lowering contacts. The electromagnetic variation
is at most 2mu|h|. All contraction tensors are physically parallel,
so their higher derivatives are bounded by the same product arrays.
No off-shell return or density term is deleted.

Define C_E=max_(k<=7)E_k and C_loc=max_(j<=2,k<=6)(Aj)_k. Equations
(2)-(9), with rational coefficients and finite truncation, specify
these constants completely. The receipt's upward decimal powers are
computed by exact integer comparison, not fitted numerical tolerances.

## 5. Verification and remaining source work

The checker pins all inherited inputs and checks the seed Hessian and
sectional identities, regional allowances, every finite coefficient
array, tensor-slot recursion and final raw Leibniz products. Independent
nonpolynomial scalar jets on a manufactured curved orbit check the
ordered tensor recurrence at the bolt and away from it. A separate
finite conformal metric calculation retains the inverse and curvature
variation contacts. These are supporting controls for the ordinary
argument, not samples of the actual charged quantum kernel or a Lean
certificate.

The comparison residual is f_ell=Jvol F_ell(g_ell), where the subscript
on F is the retained bookkeeping action. Its scalar York component is
lambda^-2 F_ell,Phi because Rcal(tau,0)=(tau,0). Thus its local H6
estimate must retain the homotopy trace and its quadratic ell zero;
the full-action coefficient Ebar in (8) is not silently identified
with that residual. Nor do weighted supremum coefficients alone give
its weighted L2 norm on the noncompact exterior.

Remaining to-do list: actual finite-metric quantum H4 and scalar residual
H6; higher TT lift/inverse and positive transfer; cancellation, signed
compatibility and nonlinear existence; physical validity, state/boundary,
stability, matching, charge/return, opening, transport, preparation,
separate creation and audit.
