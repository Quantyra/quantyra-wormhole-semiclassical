# Quantitative normal inverse and a lower bound for the physical inverse

2026-10-20 HST, E003/S163 supporting S139/S074. Ordinary proof and
protocol before execution. **Quantitative operator bound; no solution
or stability claim.**

In the precise norms below, the complete outer normal inverse satisfies

    8.6e19 < lambda^2 ||N^-1|| < 8.7e19.             (1)

The actual analytic-seed compatible inverse G0 satisfies the lower
bound lambda^2 ||G0||>8.6e19 in the same normalization. There is no
global upper bound in this note. This quantifies a necessary scale
in the inverse problem left by
[S162](general-background-newton-operator.md), using the actual
[normal family](stationary-disk-and-normal-inverse.md). It is not an
inverse of a truncated field equation.

## 1. Fixed norms and retained normal operator

Keep the action, state, domain, subtraction, flux, clock and positive
prepared temperature unchanged. Use delta=7/2 and x=-log z. Normalize
the outer mixed input by

    theta=S+T/3, b=2T,
    ||u||_X^2 =
      ||(1+D_x^2)^2(z^-delta theta)||_L2(dx)^2
      +||(1+D_x^2)^2 log(2+|D_x|)(z^-delta b)||_L2(dx)^2,
    D_x=-i partial_x.

These are the fourth/fourth-log spaces already used by the
[weighted map](stationary-global-weighted-map.md), now with exact
Fourier norm normalization. For the disk action gradient (f_S,f_T),
use the dual transpose output

    f_theta=f_S, f_b=(f_T-f_S/3)/2,
    ||f||_Y^2=||z^-delta f_theta||_2^2+||z^-delta f_b||_2^2. (2)

For the global spaces choose the usual fixed disk/interior patches,
with the outer cutoff equal to one sufficiently far out. Add the single
correlated mass coefficient in quadrature to the remainder norm.
The localized sequences below have zero mass coefficient and vanish
on every interior patch. Thus their numerical lower bound is independent
of how those compact patches are normalized. Different outer norm
normalizations require their explicit conversion constants.

Set s=7/2+it, q=s(s-3), epsilon=G/(60pi), and

    K=E_sp/lambda=60pi/(G lambda)+11/3,
    omega^2=K-9/4=60pi/(G lambda)+17/12,
    J(s)=(q-4)(q+K)/(80pi^2),
    d0=1/(480pi^2),
    H(s)=digamma(s-1)+(pi/2)cot(pi s/2)+kappa_M,
    L=A_T/d0=10pi/(G lambda)+23/18,
    N_T(s)=d0 q[(q+2)H-L].

The existing normal theorem includes the image, Weyl and finite contacts.
In (theta,b) variables its physical-volume symbol is

    N(s)=lambda^2 D(3-s)^T diag(J(s),N_T(s)) D(s),
    D(s)=[[1,c_s],[0,1]], c_s=1/[3(s+1)].           (3)

The reflection is 3-s, not s or conjugate(s). The retained benchmark
gives .999*2000^-10 < G lambda < 2000^-10 and
39/5<kappa_M<54/5. Hence L>800 and omega>10. No physical
temperature is varied in estimating this encompassing parameter box.

The normal problem is the full-line Fourier multiplier in x after
weight conjugation. Its norm is the essential supremum of the
weighted matrix multiplier. It does not impose a new condition at
the actual bolt.

## 2. Scalar inverse: uniform upper and lower bounds

Write C_phi for the scalar normal inverse L2 -> H4 norm. Exactly,

    lambda^2 C_phi=80pi^2 sup_t
      (1+t^2)^2 / [|q-4| |q+K|].

With y=t^2, the elementary identities are

    |q-4|^2=(y+1/4)(y+81/4) >= (1+y)^2,
    |q+K|^2=(y-a)^2+16y, a=omega^2+4.

For a>104, differentiation or completing a square gives

    sup_(y>=0) (1+y)^2/[(y-a)^2+16y]
      =(a^2+2a-15)/[16(a-4)]
      =(omega^2+10+9/omega^2)/16.                 (4)

The maximizer is y=a+8+48/(a-7); it is positive. This yields

    lambda^2 C_phi <=20pi^2 sqrt(omega^2+10+9/omega^2)
                         <20pi^2(omega+1).

For a lower bound evaluate the actual scalar ratio at t=omega:

    lambda^2 C_phi >=
      20pi^2 (1+omega^2)^(3/2)
          /sqrt[(omega^2+1/4)(omega^2+81/4)]
      >20pi^2(omega-11/omega).                    (5)

For the last inequality, the square root is at most omega^2+41/4,
the numerator is larger than omega^3, and
1/(1+v)>=1-v. This is a continuum bound, not a search over sampled
frequencies. It retains the quantum branch near |t|=omega.

## 3. Shape inverse: a uniform fourth-log bound

Put C_T=|| (lambda^2 N_T)^-1 ||_(L2->H4log).
In the normalization of section 1,

    C_T=480pi^2/lambda^2.                           (6)

Here is a direct proof, including all finite frequencies. From the
retained digamma recurrence and reflection identities,

    Im H=-t S_t,
    S_t=1/(t^2+1/4)+1/(t^2+9/4),
    t^2 S_t<2.

For |t|<=2 the existing strip estimate gives |(q+2)H|<400.
Since |q|>=1+t^2 and L>800, this implies

    |(q+2)H-L|>(1+t^2)log(2+|t|).

For t>=2, the positive decreasing series for Re digamma(1/2+it)
is at least its integral, log(1+4t^2)/2, plus digamma(1/2)>-2.
The two recurrence contributions have positive real part. The
real cotangent term is -(pi/2)sech(pi t)>-1/2. Therefore

    Re H >53/10+log(1+4t^2)/2
         >=53/10+log(2+t).

The last inequality holds for t>=2. The needed series, recurrence,
reflection and large-argument expansion are the standard
[DLMF series](https://dlmf.nist.gov/5.7),
[functional relations](https://dlmf.nist.gov/5.5) and
[digamma expansion](https://dlmf.nist.gov/5.11.E2).
Only these special-function identities are imported.

For 2<=t<=3, L-Re[(q+2)H]>L-8>792, exceeding
(1+t^2)log(2+t). For t>=3,

    L-Re[(q+2)H]
      >L-8+(t^2-15/4)[5+log(2+t)]
      >(1+t^2)log(2+t).

Indeed log(2+t)<t there, and
5t^2-75/4-(19/4)t is positive at 3 with positive derivative.
Complex conjugation covers negative t. Also

    |q|^2-(1+t^2)^2=(21/2)t^2+33/16>0.

Together these inequalities prove

    |N_T(7/2+it)| >=
      d0(1+t^2)^2 log(2+|t|).                     (7)

Finally N_T~d0 t^4 log|t|, so the normalized ratio tends to 1/d0.
This proves equality of the supremum in (6), not just an upper bound.
In the scalar sphere-component TT normalization D=9lambda^2 N_T,
the corresponding constant is 160pi^2/(3lambda^2). It must not be
substituted for (6) while keeping the b input/output normalization.

## 4. Complete normal matrix and the benchmark interval

Let w=(1+t^2)^2 and l=log(2+|t|). After applying the input weights,
the inverse of (3), multiplied by lambda^2, is the sum of

    diag(w/J,0)

and the rank-one matrix

    (w l/N_T) [-c_s/l,1]^T [-c_(3-s),1].

One has |c_(3-s)|<=2/3 and |c_s|/l<1/9, using log2>2/3.
The latter rank-one factor has norm at most sqrt(1066)/27<4/3.
Consequently

    lambda^2 ||N^-1|| <=lambda^2 C_phi+640pi^2.     (8)

For a lower bound, apply N to a pure conformal input (theta,b)=(1,0).
Its output is lambda^2 J(1,c_(3-s)). At t=omega this gives

    lambda^2 ||N^-1|| >=
      20pi^2 (1+omega^2)^(3/2)
       /sqrt[(omega^2+1/4)(omega^2+81/4)]
       /sqrt[1+1/(9(omega^2+1/4))]
      >20pi^2(omega-12/omega).                    (9)

For the last step bound the final square root by
1+1/(18omega^2) and combine with (5).
These are the original mixed disk norms, not a favorable change
to unrelated auxiliary source norms.

Using 3.14159<pi<3.14160 in the retained G lambda interval gives
4.39e17<omega<4.40e17. Equations (8)-(9) then imply (1).
All bounds hold throughout that parameter box, not just at one
chosen member of the prepared family.

## 5. Transfer of the lower bound to the actual physical range

Let A0 be the full analytic-seed Hessian with its physical mass column
and let G0 be its inverse on Ran(A0)=ker L0. Existence as a bounded
compatible inverse and injectivity are the earlier
[TT](stationary-weighted-tt-inverse.md) and
[outgoing](stationary-scalar-outgoing-inverse.md) results; no new
surjectivity onto arbitrary fast sources is assumed.

Choose nonzero chi in C_c^infinity(0,1). For fixed frequency t0 and
width H>1 take the pure conformal packet

    theta_(R,H)(x)=exp(-delta x) H^-1/2
                     chi((x-R)/H) exp(-it0 x), b=0. (10)

Its support lies farther down the actual end as R increases, and
its mass coefficient is zero. Its weighted input norm is independent
of R. The accepted actual-minus-normal estimates give

    ||A0 u_(R,H)-N u_(R,H)||_Y ->0 as R->infinity   (11)

for each fixed H. Identify the normal output with the end and extend
by zero through a fixed cutoff. The separated end-to-interior remainder
tends to zero by the same weighted kernel estimates. In seed York
coordinates an alternative check uses the local convergent conformal
coefficients and compactness of Bmix*:Z_(7/2)->Y_(7/2).
Thus (11) does not assume that a full-order error on a fixed compact
annulus is compact. First take R->infinity, then H->infinity.

The Fourier profile of the weighted packet concentrates at t0 as
H increases. Its H4 norm and the polynomial pure-conformal output
converge in ratio to

    (1+t0^2)^2/[lambda^2 |J(s0)|
                         sqrt(1+|c_(3-s0)|^2)].

The second output is polynomial as well:
J(s)/(3(4-s))=-(s+1)(q+K)/(240pi^2). This supplies a direct
dominated-convergence argument for smooth compact chi.

Crucially, define the actual source f_(R,H)=A0 u_(R,H).
It is in the physical range and obeys L0(f_(R,H))=0 exactly,
including all its nonlocal tails. It is not replaced by the normal
source, which need not be compatible. Injectivity gives
G0 f_(R,H)=u_(R,H). Taking the two limits and t0=omega yields

    lambda^2 ||G0|| >= the first expression in (9)
                     >8.6e19.                    (12)

Complex packets are only a convenience: the real operator on the
real Hilbert spaces has the same norm as its complexification.
Equivalently at least one real/imaginary packet has the required
ratio in the limiting argument. No extra physical mass, clock or
source datum is introduced.

## Consequence and limits

The shape normal inverse has a controlled fourth-log constant. The
large normal response is in the retained conformal quantum branch.
Every quantitative upper bound for the global inverse in (2) must
accommodate (12); the actual compact-region matching may increase it.
This lower bound does not show that the specific residual excites
that frequency or that any physical solution is unstable.

The next global estimate must control the interior/bolt and matching
remainders with this frequency dependence, and the nonlinear proof
must bound the source and response in the same norms. The earlier
weak-source outgoing bound and product cosh-weighted response are
different norms, so their numerical constants cannot be substituted.
A normal upper bound is not a global one.

Rescaling norm units alone does not improve a Newton test. If input
and output norms are multiplied by alpha and beta, respectively,
then M'=alpha M/beta, L'=beta L/alpha^2 and e'=beta e, so
2M'^2 L'e'=2M^2 L e. No value of L e is proved here. Frequency-aware
estimates may be useful, but must retain the full equations and
the remaining signed condition; (12) is not a nonexistence theorem.

## Frozen verification protocol

Revalidate the five predecessor receipts and their current input hashes.
Check the scalar identities, exact rational maximizer and positivity,
TT inequalities after their stated elementary majorants, reflected
matrix transformation, pure-conformal polynomial output and benchmark
interval using rational arithmetic.

At 100 digits use both enclosing endpoints of G lambda and kappa_M.
Check the exact digamma imaginary-part identity and the proven scalar,
shape and complete-matrix bounds at t=0,.25,1,2,3,10,1e4, omega/2,
omega-4,omega-1,omega,omega+1,omega+4,2omega,10omega,1000omega.
Use relative tolerance 1e-70 only for identities. Inequalities must
have the stated sign without an added tolerance.

For a full-line normal diagnostic, use a Gaussian Fourier density
centered at t=omega with variance 1/(2H^2). Both squared input and
output norms are degree-eight polynomials, so their Gaussian integrals
are exact finite moment sums. At omega=10,100,1e4 and the midpoint
benchmark, use H=4,8,16,32. The relative error from the point-frequency
ratio must be below .02 and refine by a factor below .35, or lie
below 1e-70. This is not numerical shooting at the quantum frequency
and is not a computed global exterior.

The continuum estimates and localization argument remain ordinary
proofs; sampled multipliers do not establish their suprema. Preserve
all runs and failed sources, and raw UTF8/LF evidence. No claim of a
global upper inverse bound, nonlinear closure, signed compatibility,
an exact solution or physical stability is tested or certified here.

Remaining to-do list: quantitative global compatible inverse and
nonlinear family; signed compatibility and coupled exterior/state,
charge and return; opening, stability, transport, preparation,
creation and full-result audit.
