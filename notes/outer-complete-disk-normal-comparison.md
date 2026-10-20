# The complete localized disk operator compared with its outer normal family

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Full localized operator
comparison; global matching and the research goal remain open.**

The massive comparison can now be combined with the actual local
Einstein-Maxwell, Weyl and LLL terms. The logarithm cancels before
estimation. The limiting TWO-ROW matrix is the retained normal
family, including its finite contacts and off-diagonal entries.

## 1. Statement and fixed norms

Use the seed, state, clock, full angular tower, paired domain and
finite subtraction of the preceding
[massive comparison](outer-complete-massive-model-comparison.md).
Let A be the complete physical-volume affine disk Hessian, including
the off-shell density contact. Use theta=S+T/3, b=2T, with output
(f_theta,f_b)=(f_S,(f_T-f_S/3)/2). Put delta=7/2, x=log Y and

    ||v||_X^2=||Y^-delta theta||_H4(dx)^2
                 +||Y^-delta b||_H4log(dx)^2,
    ||f||_Z^2=||Y^-delta f_theta||_2^2+||Y^-delta f_b||_2^2.

The weighted inputs are zero extended on the x-line. H4log has
multiplier (1+t^2)^2 log(2+abs(t)). Define M_epsilon as in the
preceding comparison: it equals m=sqrt(a/a(0)) on Y<=epsilon
and one past 2epsilon. The transformed equation has variables
v_hat=M_epsilon v, f_hat=M_epsilon f. Both transformation norms
and their inverses are below 1.01. Write Ahat=M_epsilon A M_epsilon^-1.

For smooth even stationary v_hat supported in 0<Y<epsilon,
0<epsilon<=1e-15, the result is

    ||1_(Y<epsilon) (lambda^-2 Ahat-N0)v_hat||_Z
                          <=1e67 epsilon ||v_hat||_X,          (1)

where lambda^2 N0 is the complete normal operator of
[the normal inverse](weighted-normal-inverse-scale.md).
In particular its inverse on the full log line gives

    ||N0^-1 1_(Y<epsilon)(lambda^-2 Ahat-N0)||_(X->X)
       <8.7e86 epsilon <1e-30 at epsilon=1e-120.               (2)

This is a localized comparison, not an inverse with a new collar
boundary condition. The discontinuous output restriction is only
an L2 restriction, followed by zero extension. The normal inverse
need not preserve support. Mass localization, cutoff commutators,
the intervening annulus and the regular bolt must still be assembled.

## 2. Exact radial density and the half-density change

Set alpha=xi/C, xi=lambda B^2, k=sqrt(alpha), u=B/r and

    w=u/(kY), j=a/a(0)=u'/k, h=log j, v=log w,
    rho=w^4, D=Y partial_Y, H=Dh,
    O=-log(kY), omega=O-v,
    e=Y^2/C, z=Y^2, nu=q/C, chi=alpha Q/(B^2 C),
    ell=lambda0/lambda=1+11g/(360pi), g=G lambda.

The physical charge is Q, not the comparison charge. The exact
off-shell reference-action cancellation from the
[complete local action](complete-compact-disk-graph.md) is retained.
The benchmark bounds give 0<alpha<1, nu<400, abs(chi)<2,
1/C<1 and 1/g<2e33. Also .4<k<.5, so abs(log k)<1.

Define operators on the transformed tangents and rescaled background:

    A1=D-H/2,
    Vhat=(DH-H)/2+H^2/4, B2=D^2-D-Vhat,
    r=-2(DH-H+H^2), o1=-1-Dv, o2=-D^2v,
    U=o2+(H-1)o1, V=o1^2, K=e^2/15-nu^2 z^2/3.

Indeed Rbar=Y^-2 r and Delta omega=Y^-2 U. For a tangent f,

    partial_Y(m^-1 f)=Y^-1 m^-1 A1 f,
    Delta(m^-1 f)=Y^-2 m^-1 B2 f.                         (3)

Each local quadratic density contains two tangent factors.
The canonical measure is j dY (apart from its fixed time factor),
and j m^-2=1. Thus the conjugated canonical local Hessian is the
Euler operator of the following densities against Y^-4 dY.
This is also a direct check of the equation similarity: the
Euler derivative of I''[m^-1 f,m^-1 f] in dY equals m L m^-1,
where L was the original Hessian in j dY.

For shorthand let S=theta-b/6, T=b/2 and

    d=B2 T, f=B2 theta, t=A1 T, s=A1 S, w1=o1 s,
    DD=S(-r/15+e/9)+T[(6omega+1)r/45+11V/135]-11w1/135,
    LL=T[-r/15+2U/5+7V/45]-7w1/45,
    JJ=r/30-7(U+V)/90-e/18,
    WW=-2r/15+14(U+V)/45,
    ZZ=T^2[omega r^2/30+r(U+V)/15-U^2/5
              -7UV/45-7V^2/90+2K omega]+ST[-r^2/30+2K].

The three exact quadratic polynomials are

    Qw=f^2/20-(omega+1/6)d^2/30
                    -(d DD+f LL+JJ s^2+WW T w1-7w1^2/45+ZZ)/2,
    Qc=-1/4 {w^-2[4S^2(r+6V)+48S o1 s+16S o1 t
                       +12s^2+8s t+8e(S+T)^2]
                +6ell w^-4(4S+2T)^2-8chi z^2 T^2},
    Ql=-nu z t^2/12.                                      (4)

The action second variations are respectively C/pi, C/g and
C/pi times the integrals of (4) against Y^-4 dY. For a quadratic
polynomial Q, let

    E_i(Q)=(1/2) sum_(a=0)^2 (3-D)^a partial Q/partial(D^a i).

The factor 1/2 is required because Q represents I'' rather than
I''/2. The adjoint D*=3-D follows from the measure Y^-4 dY.
The complete scaled local physical operator is therefore

    rho {E(Qw)/(4pi^2)+E(Qc)/(4pi g)+E(Ql)/(4pi^2)}.       (5)

The physical factor is C u^4/(4pi xi^2)=rho Y^4/(4pi C).
Both input and output b identifications have already been made
in (4). In particular the massive b entry has the preceding
normalization rho Y^4/(16pi C), with no missing factor four.

## 3. Full limiting matrix and cancellation of the logarithm

Take h=0, v=0, w=1, e=z=0 in (4), keeping DO=-1.
Then r=0, U=V=1, A1=D and B2=D(D-1). Denote the limiting
local operator inside braces in (5) by L0. Its coefficients
are polynomial in D and affine in O. The flat massive term
must be combined with L0 before estimating either one.

On an even Mellin power Y^s the flat massive b response is

    d0 {P4(s)[F_M(s)+kappa_M+O+1/6]-s(s-1)(2s-5)},
    d0=1/(480pi^2), P4=s(s-1)(s-2)(s-3),
    F_M=digamma(s-1)+(pi/2)cot(pi s/2),
    kappa_M=-6/5-log(mu0 ell_A), ell_A=lambda^-1/2.        (6)

To see the derivative term, apply log|D_Y| to partial_Y^2 Y^s
using the accepted even Mellin formula, then apply the outer
two derivatives to the remaining log Y. Also ell_o/k=ell_A.
The O-dependent Weyl term and the derivative term in (6)
cancel exactly. In particular -31/30-1/6=-6/5.

Set q_s=s(s-3), K0=60pi/g+11/3 and F=F_M+kappa_M.
The full remaining matrix is exactly

    N0_11=(q_s-4)(q_s+K0)/(80pi^2),
    N0_12=(s-4)(q_s+K0)/(240pi^2),
    N0_21=-(s+1)(q_s+K0)/(240pi^2),
    N0_22=d0 P4 F
            -[1/(48pi g)+23/(8640pi^2)]q_s
            -(q_s+K0)/(720pi^2).                         (7)

These are all four entries of
Dmat(3-s)^T diag(J,N_T) Dmat(s), with
Dmat(s)=[[1,1/(3(s+1))],[0,1]], as required. The reflection
is 3-s. This checks finite contacts and both cross entries,
not only the highest derivative symbol. For example the scalar
Weyl terms give [P4-7q_s/3]/(80pi^2), while the corrected
classical vacuum supplies 3(q_s-4)/(4pi g)-11/(60pi^2).
Their sum is precisely N0_11.

The Mellin identification uses the accepted even-kernel theorem
on the strip containing Re s=7/2, then weighted completion.
It includes the reflected kernel with no extra half. It does
not assert that individually growing log Y terms are bounded.

## 4. Quantified local error and retained density

The analytic inverse disk established in the preceding actual
remainder gives abs(P/xi-1)<.001 and abs(j-1)<.001 for abs(Y)<.1.
Since u=k integral_0^Y j(t)dt, w is the segment average of j;
abs(w-1)<.001 on that disk as well. Thus h, v, w^-2, w^-4 and
rho are analytic there, with modulus below two. Their deviations
from their values at zero have zero constant and linear terms.
Cauchy on radius 1/20 and the exact positive series sum
sum_(n>=2) n^5(1/2)^(n-2)=4326 imply, through five D derivatives,

    abs(D^a h), abs(D^a v),
    abs(D^a(w^-2-1)), abs(D^a(w^-4-1)),
    abs(D^a(rho-1)) <1e7 Y^2, 0<=a<=5.                  (8)

The first-order vanishing follows from P'(0)=0, j'(0)=0.
These are bounds on the one-sided analytic seed, not a smooth
even continuation assertion across its cubic seam.

Here is a finite coefficient ledger for the COMPLETE local
operator. Apply E to (4), then replace w^-2=1+A and w^-4=1+B.
Use D(w^-2)=-2w^-2 Dv and D(w^-4)=-4w^-4 Dv while
differentiating. Subtract the limiting local operator. Every
remaining monomial has at least one factor among

    A, B, h_1,...,h_4, v_0,...,v_4, e, z.

Assign absolute weight one to these factors and to O and all
tangent jets, weight 400 to nu, two to chi and two to ell.
The sums of absolute weighted monomial coefficients over both
rows, for EACH of E(Qw)-E(Qw0), E(Qc)-E(Qc0), E(Ql),
are below 1e12. This ledger is checked by exact expansion.
It also checks tangent order <=4, background order <=4 and
degree in O <=1. No coefficient is omitted as lower order.

On the collar each small factor is bounded by eta=1e7Y^2<1.
Hence the sum of the coefficients in the local operator
difference, including its prefactors in (5), is below
6e45 eta(1+abs(O)). For a weighted tangent Y^delta f(x),
D^a becomes (partial_x+delta)^a; its H4 or H4log-to-L2 norm
is below 1e3 for a<=4. The latter uses log2>1/2.
Since abs(O)<=1+log(1/Y) and
Y^2[2+log(1/Y)]<4Y for Y<1,

    ||(L-L0)v_hat||_Z <1e58 epsilon ||v_hat||_X.          (9)

The physical rho is common to the local and flat massive
operators. Combine their limiting terms first, giving rho N0.
For the remaining density difference, an explicit safe bound is

    ||N0||_(X->Z)<1e38.                                (10)

Indeed on s=7/2+it, with y=t^2, abs(q_s), abs(q_s-4),
abs(q_s+2) are below 6(1+y); K0<4e35. The accepted digamma
estimates give abs(F)<30log(2+abs(t)). The direct four entries
in (7), or the reflected triangular factors, then prove (10)
using log2>1/2, abs(c_s)<1 and abs(c_(3-s))<1.
No estimate of N0^-1 is used to bound N0.
Consequently the density error is below 1e45 epsilon^2.

Finally the exact decomposition on supported input is

    lambda^-2 Ahat-N0
      = massive_difference +rho(L-L0)+(rho-1)N0.          (11)

The first term has norm <=1e65 epsilon by the preceding proof;
it already includes all actual massive global/constant terms.
Use rho<3/2, (9)-(10), and
1e65+1.5e58+1e45<1e67 to obtain (1). The accepted inverse
norm <8.7e19 then proves (2). No massive term is counted twice.

## 5. Verification contract and remaining work

Freeze this ordinary argument and checker before execution.
Pin the complete massive, complete local disk and quantitative
normal receipts. Verify the radial half-density identities,
weighted Euler adjoint, physical source normalization and all
four normal matrix entries, including the logarithm cancellation.
Expand every local error coefficient, check its vanishing factor,
orders and exact coefficient ledger; verify the analytic series,
weighted derivative, forward normal and final composition budgets.
Preserve every run, failed source snapshot and raw output.

The algebra and rational checks support the ordinary analytic
proof; they do not formally certify the continuum operator,
global inverse, nonlinear branch or stable physical wormhole.
The small localized norm does not control a cutoff commutator,
a noncompact mass input, or propagation through an intervening
annulus by itself.

Remaining to-do list: mass/fast localization and nonlocal cutoffs;
annulus/bolt matching; higher York reconstruction and full
inverse/residual composition; nonlinear family and signed
compatibility; exterior/state, charge, return, opening, stability,
transport, preparation and separate creation.
