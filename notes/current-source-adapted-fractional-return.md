# Actual reconstruction-to-flow fractional bound

2026-10-28 HST, E003/S169. Ordinary argument and quantitative contract
BEFORE execution. Previous turn: progress, 1b04eae bounded the entire
local return. Preserve that result, 8dc41ff, P32, original physical norms,
closed audit and the full stable semiclassical wormhole objective.

The remaining first-reference test is

    Kgeo N_9/2(V x_hat)+Eapplication <4.308605e32,
    Kgeo<1.554891e-32.

The actual trial is not H5. Instead bound its LINEAR reconstruction map
on two dense smooth endpoint spaces and interpolate the map. This does
not assume an H5 norm for the actual input, discard its bolt logarithm,
or replace the actual paired flow by the auxiliary geometry.

## Spaces, field and endpoint domains

Use the auxiliary complete hyperbolic disk

    g2=dx^2/q2+q2 dtheta^2, q2=2x+a x^2, theta period 2pi,
    H=-partial_x(q2 partial_x), x>=0.

Its curvature is -a. H is the radial Friedrichs Laplacian with the
regular bolt and decay at infinity already used to DEFINE the trial.
Input spaces have norm ||(1+H)^(n/2) f||L2(dx), n=4,5. They are auxiliary
profile spaces, not selected quantum modes or a new physical domain.

On the ORIGINAL physical interval let s=gamma x=L-z, Q2=2gamma s+a s^2,
R=Q/Q2, q=Q/gamma^2 and r=(H+muY)^-1 f. Put

    D(s)=3Q''(z)-z Q'''(z), W(z)=int_0^z -D(L-t) r(L-t) dt,
    M=z q_xx+4gamma q_x,
    A=z^5 R/L^5,
    B=z^4(2M/3-z R muY)/L^5,
    C=z^4(5R-z R_s)/L^5,
    E=-2z^4/(3L^5).

The complete actual reconstruction from BOTH York rows is

    T f=V x_hat=A f+B r+C <grad_2 s,grad_2 r>+E W.       (1)

This is exactly the previously accepted field; W_s=D r. The physical
shape norm uses original dz=gamma dx and the averaged unit-direction
measure. The actual flow generator P=iX is the complete C2,1 paired
Liouville generator, not H or its square root.

For smooth radial auxiliary inputs, the resolvent r is regular at the
bolt by its regular-flux ODE. Formula (1) is smooth on the physical half,
starts at z^4 at the outer seam and has the original regular bolt.
The paired extension has matching fourth jets. Its fifth weak derivative
is locally square integrable; possible distributional metric coefficients
multiply vanishing lower source jets. Equivalently X^4(Tf) has matching
traces across the paired seam, so its weak flow derivative has no seam
delta. Thus Tf belongs to Dom(P^5) on this smooth core. The actual
mass/log realization, fixed temperature, both York equations and their
boundary data are unchanged. No independent phase or neighborhood claim.

## Quantitative reference covariant derivatives

Write S_j(f)=||nabla_2^j f||L2. On the complete curvature -a disk,

    S_j(f)<=d_j ||(1+H)^(j/2) f||,
    d_0=d_1=1,
    d_j=product_(k=1)^(j-1) sqrt(max(1,4k^2 a)).          (2)

Here is an ordinary justification of the explicit constants. The
commutator of the rough Laplacian with one covariant derivative on a
rank-k covariant tensor contains a Ricci term of norm at most a and
two curvature contractions per tensor index. In dimension two each
curvature action is a permutation minus a metric insertion/trace, of
norm at most 3a. Parallel curvature removes derivative-curvature terms.
Commuting through k derivatives therefore costs at most
a sum_(j=0)^(k-1)(1+6j)<=4k^2 a in the rank-k derivative norm.

The scalar operator A_k=(nabla^k)*nabla^k is a polynomial in H: its
coefficients are invariant, its isotropic principal symbol is a constant
multiple of |xi|^(2k), and successive subtraction of powers of H reduces
the order; odd invariant scalar symbols vanish. Integration by parts on
the complete smooth core and the commutator estimate consequently give

    A_(k+1)<= (H+4k^2 a) A_k,
    A_j<=H product_(k=1)^(j-1)(H+4k^2 a).

Their common scalar spectral calculus gives (2). Closure supplies the
stated domains; no boundary term at infinity or the regular bolt is
introduced. For r=(H+muY)^-1 f and muY>1, (2) bounds S_j(r) by d_j
times the order-n input norm whenever j<=n+1. Indeed the remaining
multiplier (1+p)^((j-n)/2)/(p+muY) is <=1 for p>=0.

The elementary Hessian identity is nabla_2^2 s=(gamma+a s)g2.
Writing b=gamma+a s, the pointwise norms of its next derivatives are

    m1=sqrt(Q2), m2=sqrt(2)b, m3=sqrt(2)a sqrt(Q2),
    m4=2a b, m5=2a^2 sqrt(Q2), m6=2sqrt(2)a^2 b.          (3)

For a scalar coefficient F(s), the ordinary covariant chain/product rule
gives ||nabla^j F||<=sum_k |F^(k)(s)| Bell_(j,k)(m1,...).
The computation uses x jets and (3)/gamma to avoid scaling loss; this is
the identical bound, not a changed norm. Whole interval cells enclose all
derivatives of A,B,C,E,D through the needed orders before maximization.

Products in (1) are estimated by the full multinomial rule. For n>=1,
nabla^n W=nabla^(n-1)(D r ds), while ||W||L2(ds)<=L sup|D| ||r||L2(ds).
This retains the entire finite-domain return. Norms for every derivative
S_0(Tf),...,S_n(Tf) follow using (2)-(3); the derivative term in (1)
requires S_(n+1)(r), already supplied by its two-order resolvent gain.

## Actual metric and connection, including the bolt

In the g2 orthonormal frame the actual metric is diag(1/R,R), with
.1<R<2 proved for the exact finite seed. Both volume forms are dx dtheta.
Let Dconn=nabla_actual-nabla_2. Its only nonzero frame components are

    Dconn^1_11=-b1, Dconn^2_12=Dconn^2_21=b1,
    Dconn^1_22=c1,
    b1=(partial_radial R)/(2R),
    c1=(1-R^2)(q2'/(2sqrt(q2)))-R partial_radial R/2.      (4)

To enclose all three derivatives of this tensor at the bolt, put h=sqrt(x),
S=sqrt(2+a h^2), R=1+h^4 U and

    U=[b3 gamma+b4 gamma^2 h^2+b5 gamma^3 h^4+
                                            b6 gamma^4 h^6]/(2+a h^2).

Then b1=h^3 S(4U+h U')/(4R) and
c1=h^3[-U(2+h^4 U)(1+a h^2)/S-R S(4U+h U')/4].
Every j-th covariant derivative, j<=3, has the explicit factor h^(3-j).
A radial derivative acts by (S/2)partial_h. An angular derivative applies
the frame rotation to every upper/lower tensor index and multiplies by
(1+a h^2)/(h S). Factoring h BEFORE evaluating leaves analytic coefficient
jets even at h=0. This keeps the cancellation needed for finite constants;
no singular polar coefficient is bounded in isolation. Denote the whole
domain tensor bounds by D_j, j=0,...,3.

For T_(k,j)=||nabla_2^j nabla_actual^k(Tf)||L2, the exact connection
difference and product rule give

    T_(k+1,j)<=T_(k,j+1)
             +k sum_(i=0)^j binom(j,i) D_i T_(k,j-i),
    T_(0,j)=S_j(Tf).                                  (5)

Finally |nabla_actual^n(Tf)|actual<=10^(n/2)|...|g2,
and contraction with a unit tangent bounds X^n(Tf). Thus (1)-(5)
give explicit C4,C5 with

    ||P^n Tf||L2(dz dangle/2pi)
       <= Cn sqrt(gamma)||(1+H)^(n/2)f||L2(dx), n=4,5.   (6)

The volume, radius scale and finite-domain costs are retained. These
are quantum-source flow norms, not the original C5 geometric admission.

## Interpolation, application errors and decision test

Apply the spectral interpolation theorem to the LINEAR map T in (6).
One can regularize |P| by a positive spectral floor and truncate both
spectral measures, apply the three-lines argument to the resulting
bounded analytic family, and then take monotone limits. The imaginary
powers of the two self-adjoint generators are unitary. The result is

    N_9/2(Tf)<=sqrt(C4 C5) sqrt(gamma)
                              ||(1+H)^(9/4)f||L2(dx).   (7)

The actual accepted f=v0 belongs to this input space, with norm bounded
by At*Nmodel. It need not belong to the order-five endpoint space.
This addresses the first missing quantitative comparison without asserting
an H5 source estimate or removing the nonzero bolt logarithm.

Execute Cflow=sqrt(C4 C5), Nactual=Cflow*At*sqrt(gamma)*Nmodel and
Egeo=Kgeo*Nactual. Compare Egeo with 4.308605e32 and combine it with all
accepted source/model/local-return entries. If it is excessive, report
the estimator's required improvement; do not infer a construction
obstruction. The common-regulator, full angular/thermal and flow-return
argument for Kgeo remains unchanged, as do all local contacts already
charged in 1b04eae.

An exact analytic trial definition has no additional numerical application
beyond its enclosed coefficients and analytic norm proof. Any implemented
use/evaluation of it in the ordered correction must, however, carry its
actual caller/application errors. The latter are NOT assumed zero by
(7). Distinguish a proved residual for the exact analytic reference trial
from an implemented first actual inverse, later quantum applications,
signed G_t, the full coupled budget or a physical admission certificate.

For the first-stage norm comparison use the already accepted
||A0-I||<=e0, with e0 near .0022. A reference residual Rref implies
||A0^-1 g_hat-x_hat||Et<=Rref/(1-e0). Compare that energy error with the
recorded common stage OUTPUT allowance; do not equate a residual in Et*
with that output error without this factor. The recorded K_first times
a residual is a coefficient product. Applying it to this reference
application requires the above coercivity conversion and, for a source
error repeated in the ordered iteration, its recorded amplification sum.
All such products remain components, not complete bolt/admission bounds.

Also quantify the coarse actual-minus-reference estimate
e0*(Ghat+Rref)/(1-e0). This is only an UPPER REPORT for the uncomputed
actual response difference. If excessive, it shows why merely declaring
the reference trial to be the actual first solution is not certified.
The ordered correction is intended to evaluate and cancel that term;
its amplitude need not fit its numerical ERROR allowance. Preserve this
distinction when selecting the next complete operator application.

Remaining to-do list: execute (7) and its complete residual consequence;
numerical applications and later actual quantum returns; signed G_t,
all coupled/admission requirements, nonlinear existence, compatibility,
matching, semiclassical validity/noise and all-sector stability.
