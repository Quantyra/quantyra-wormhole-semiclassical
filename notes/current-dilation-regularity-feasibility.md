# Complex dilation: a seed resolvent bound and a neighborhood regularity obstruction

2026-10-24, E003/S169. Ordinary argument frozen before acceptance.
Full objective ACTIVE and INCOMPLETE; feasibility audit CLOSED.

The question is whether complex dilation can supply the missing contact
jets without extending the invalid covariant-probe division. Two different
requirements must be separated: a holomorphic first-order inverse at the
seed, and uniform removal of the common regulator in the strong response
throughout the admitted neighborhood. The first can be bounded. The usual
attempt to obtain the second by five Cauchy derivatives in a uniformly
bounded W3,infinity coefficient family fails on the current fast graph.
This is a method obstruction, not a singularity of the quantum operator.

## 1. Exact comparison domain and a finite parameter family

Use ONLY the auxiliary cylinder in current-coincident-domain-transfer.md,
with its already retained actual-minus-cylinder difference <2.012953e-13.
It has Q_c=Q on |z|<=.04, Q_c=beta outside |z|>=.05, and the same even
cubic seam. On the intervening annulus use chi=1-I_t(9,9), t=100|z|-4.
The taper is C8; the seam is C2,1. The original finite geometry and its
regular bolt remain in the comparison difference. This is neither a
smooth physical surrogate nor a new boundary realization.

Let D=z partial_z and define, on REAL signed z,

 Q_lambda=beta+sum_(j=0)^5 lambda^j D^j(Q_c-beta)/j!.       (1)

No taper is evaluated at a complex spatial argument. This is an entire
parameter family into W3,infinity: the cap is C8 and the cubic seam is
retained exactly. Its first five jets at zero agree with real dilation
Q_c(exp(lambda)z), including in W3,infinity. The statement concerns this
SEED. It does not extend a C5 neighborhood to C8.

The actual correlated polynomial is

 Q=beta+u z^2+(-4beta-2u-tau)z^3
                  +(3beta+u+ep)z^4-2ep z^5+ep z^6.

For each j<=5 bound |D^j(Q_c/beta-1)| on the whole line: on the cap by
exact polynomial Bernstein coefficients on 64 subintervals; on the inner
interval by the absolute polynomial coefficient sum at .04; outside by
zero. Take the larger bound, not a sampled cap maximum. Equation (1)
then gives d(r)=sum g_j r^j/j! as a uniform relative deviation bound.
The acceptance records all six guards and both spatial pieces. In
particular d(.20)<.104; d(.25)<.285. The earlier exploratory cap-only
calculation is retained and is not used as a whole-line certificate.

## 2. Complete-mode resolvent estimate on a fixed first-order domain

In the density-conjugated radial transfer equation put

 T_lambda=partial_z-exp(lambda)[M Q_lambda^-1/2 sigma3
                                      +p Q_lambda^-1 sigma1],
 k=sqrt(M^2/beta+p^2/beta^2), M=M_n>0, p=p_j.

The two Pauli matrices here are Hermitian, have norm one and anticommute.
This is the existing first-order equation in a fixed frame, not a square
with a guessed operator domain. For each mode Dom T_lambda=H1(R,C2).
The continuous transmission trace is compulsory: a jump creates a delta
outside L2. The original paired fold is retained. Scalar coefficients
preserve its fixed reflection projector. All modes are included; the
LLL remains in its existing separate allocation.

For d<1 define

 h1=(1-d)^(-1/2)-1, h2=d/(1-d), h=sqrt(h1^2+h2^2).

The convergent binomial series bounds the two coefficient changes by
h1,h2. With alpha=M/(sqrt(beta) k), b=p/(beta k), alpha^2+b^2=1,
the triangle inequality followed by Cauchy--Schwarz bounds the matrix
perturbation by |exp(lambda)| k h. No self-adjointness of a complex
coefficient matrix is assumed.

The free matrix exp(lambda)[M/sqrt(beta) sigma3+p/beta sigma1] has
orthogonal constant projectors and eigenvalues +/-exp(lambda)k. Its
inverse kernel has norm at most exp[-Re(exp(lambda))k|z-z'|], with the
appropriate decaying projector on each side of its identity jump.
Its L2 inverse norm is <=1/[Re(exp(lambda))k], by Fourier transform.
For |lambda|<=r<pi/2 the Neumann ratio in L2 is <=h sec(r). Thus

 ||T_lambda^-1||_2 <= e^r/[cos(r)(1-h sec(r)) k].          (2)

For the pointwise kernel choose decay a=e^-r cos(r)/2. In the weighted
kernel convolution, integration of exp[-(Re(e^lambda)-a)k|s|] costs
2/[(Re(e^lambda)-a)k]. Since a<=Re(e^lambda)/2, its Neumann ratio is
at most 4h sec(r). Consequently, when this is <1,

 |T_lambda^-1(z,z')| <= exp[-a k|z-z'|]/[1-4h sec(r)].    (3)

Bounded perturbations preserve closedness. Equations (2)-(3) establish
the inverse and its holomorphy directly by a normally convergent Neumann
series on the fixed domain. They also hold uniformly for the Hilbert
direct sum with its natural first-order graph domain. At r=.20 use the
OUTWARD guard d=.104. Acceptance verifies a>.40, L2 constant <1.45,
and kernel constant <2.13. Positive and negative thermal frequencies
obey the same bound. There is no angular or thermal cutoff.

This closes a seed first-order domain/resolvent prerequisite. It does
not identify a trace-class difference or a strong renormalized Hessian.
In particular, its bare M^-2 angular envelope has the EXACT partial sum

 sum_(n=1)^N d_n/M_n^2=2C[H_N+H_(N+q)-H_q] -> infinity. (4)

This divergent estimator is not a divergence of the subtracted action.
The common local contacts and principal subtraction must be combined
before this sum and before using Cauchy on the response. Separate
modewise holomorphy cannot replace that missing argument.

## 3. What a certified complex response would have to accomplish

The unchanged real targets are M_A<=3.18e27, M_V<=2.66e30 and
eta_nonlocal<.1599164, in current-symmetric-contact-feasibility.md.
For illustration, a holomorphic extension of the ACTUAL functions
A(exp(lambda)z), V(exp(lambda)z), with uniform bounds S_A L(z), S_V L(z)
on |lambda|<=.20, would give

 M_A<=375000 S_A, M_V<=15000 S_V.

Thus sufficient complex-amplitude targets would be S_A<=8.48e21 and
S_V<=1.773333333333e26. These are thresholds, NOT computed amplitudes.
For the more natural dilated density V_lambda=exp(lambda)V(exp(lambda)z),
converting (partial_lambda-1)^j to D^j gives the fourth-order factor
sum_(k=0)^4 binom(4,k) k!/.2^k=18321, hence S_V<=2.66e30/18321.
This alternative is conditional on the full scaling identity, including
the induced regulator scale and every local anomaly. The resolvent
estimate alone supplies neither identity nor any of these amplitudes.
No finite contact is adjustable to make the target pass.

## 4. Executed continuum test of a uniform neighborhood extension

There is an earlier regularity issue if one tries to repeat section 1
uniformly on the admitted fast family. A C5 bound is not a bound on
three spatial derivatives of five dilation derivatives. This remains
false even in the stronger FULL normal graph used for fast admission.

Here is a continuum counterexample on the actual paired domain, not a
finite physical mode model. Set x=log z, z0=.002, x0=log z0. Fix a real
smooth compact bump psi(t), equal to exp(-t^2) for |t|<=1/2 and zero for
|t|>=1. It has psi^(8)(0)=1680. For N>=100 put

 w_N(x)=N^(-15/2)/log(N) psi(N(x-x0)),
 b_N(z)=z^sigma w_N(log z), sigma=4.005, nu_N=0.         (5)

Each field is smooth, identical on both physical copies, supported
strictly inside the collar, zero near seam and bolt, and has zero mass
and slow coefficients. All boundary and regularity requirements of
these test directions hold exactly. Their admitted weighted C5 norms
tend to zero. No physical source is chosen to cancel a residual.

For completeness the ACTUAL normal multiplier has a global bound

 |B(sigma+it)| <= C_B (1+|t|)^5 log(2+|t|), C_B<2.     (6)

Use the retained B=-K(s-4)s(s-3)[(s-1)(s-2)H(s)-L]. The digamma
series [DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6), split at n=ceil|t|,
gives |psi(3.005+it)|<=log(2+|t|)+4. The cotangent term has modulus
<201 and |kappa|<11, so |H|<320 log(2+|t|). Bounds
|s-4|<=1+|t|, |s|<=5(1+|t|), |s-3|<=2(1+|t|), and
|(s-1)(s-2)|<=12(1+|t|)^2 yield

 C_B=10K[3840+L/log 2]<2.

This uses every frequency, including the large classical coefficient L.
The split-series tail is bounded by |t| sum_(n>ceil|t|)(n+3)^-2<1;
the initial part is bounded by a harmonic integral. The small-|t| case
obeys the same guard. The checker verifies only the finite constants;
the series argument supplies (6).

With the unitary Fourier convention define the finite positive constant

 Cpsi^2=integral (1+eta^2)^3(1+|eta|)^10
       [1+log(2+|eta|)/log 2]^2 |psi_hat(eta)|^2 d eta.

Its finiteness follows from repeated integration by parts on the compact
smooth bump. Rescaling the ENTIRE Fourier integral, (6) gives
||B_sigma w_N||_H3<=2 Cpsi for all N>=2: the squared factors are
N^-15/log^2 N from amplitude, N^-1 from the spatial scaling, N^6
from H3, N^10 from the symbol, and log^2 N. Their product is one.
The coupled scalar source is also bounded: theta=a_op(D)b when nu=0,
and J_c(s)a_op(s) is a polynomial of degree FOUR (the s+1 denominator
cancels). Its H3 norm is O(N^-1/log N), hence tends to zero. Thus a
FIXED rescaling places the whole sequence in any positive-radius
coupled fast graph ball. This is not a frequency-dependent normalization
of that ball or a change of physical source norms.

On the other hand the exact local chain rule is

 partial_z^3 D^5 b_N=z^(sigma-3)
       (sigma+partial_x)^6(sigma+partial_x-1)
                              (sigma+partial_x-2) w_N.

Its degree-eight leading coefficient is one. At z=z0 this equals

 z0^(sigma-3)[1680 sqrt(N)/log(N)+O(N^-1/2/log N)].     (7)

It is unbounded. Lower terms actually improve for this even bump, but
the displayed remainder suffices. Therefore no finite uniform constant
can bound ||D^5 b||_W3,infinity by this coupled fast graph norm, even
on smooth fields with exactly the original boundary traces. Taking a
smaller positive ball or multiplying norms by fixed constants does not
remove the divergence.

This also rules out ANY uniformly W3,infinity-bounded holomorphic
parameter family of fixed radius with these five real dilation jets:
Banach-valued Cauchy would bound its fifth jet, contradicting (7). It
is not just an unfortunate choice of the finite polynomial (1).

The test rejects a proposed UNIFORM PROOF METHOD. It does not reject
the actual first-order domain, weak C2 action on admitted families,
seed contact estimates, or a tame strong response in Sobolev spaces.
The eighth derivative in (7) is harmless in L2 for this sequence:
||partial_x^8 w_N||_2=O(1/log N). Requiring it in L-infinity introduces
the extra loss. Full nonlinear X-ball admission is still unresolved;
if X were strengthened, its inverse and data would require a new audit.

## 5. Decision boundary

EXECUTED: whole-line seed coefficient bounds; a complete-mode fixed-domain
complex resolvent and decay bound; exact divergent bare angular estimator;
conditional contact-amplitude thresholds; and the continuum counterexample
to extending five dilation jets uniformly through W3,infinity on the
current fast graph. Supporting checks verify polynomial, symbol-constant
and scaling identities. They do not replace these continuum arguments.

UNKNOWN: the common-regulated strong contact/nonlocal response and its
neighborhood variation, all actual contact amplitudes, global correction
data, and nonlinear closure. No complete excessive finite upper bound
has been obtained for that operator. The divergent bare estimator and
the W3,infinity embedding obstruction are failures of methods, not
nonexistence theorems. The seed tube is useful but cannot certify them.

RESEARCH DECISION: retain the construction and symmetric subtraction,
but reject a uniform five-jet Cauchy argument through a W3,infinity tube
as the nonlinear bridge. Use a direct Sobolev estimate of the subtracted
response, keeping the highest derivatives in L2. The seed tube can
support the linear calculation; it must not be promoted to Z2 control.
Do not spend another milestone tightening this tube or cold/scalar bounds.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): estimate the first metric
variation of the complete symmetric, common-regulated subtracted quantum
Hessian as a bilinear map from the coupled fast graph into H3. Keep the
eighth derivative in L2, and test whether the common-regulator limit has
a finite uniform norm on the admitted mass/log family. Retain all
contacts, seagulls, reflected and thermal terms and actual-domain return.
This is the decisive test of a tame nonlinear bridge; a seed resolvent
bound or differentiation of a real envelope is insufficient. If finite,
its numerical contribution must fit the remaining total Z2 allocation
<=5e7, alongside the still-open seed targets M_A,M_V and eta_nonlocal.
No part of that allocation is certified or silently assigned twice.

No full signed projection is executed here. Preserve the earlier
sufficient enclosure E_rest<2.123e-42<1.287e-40 with zero as an enclosure
center, not a measured sign. It does not supply an actual solution.
At rho=1e-8 the requirements remain
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1; the allocations
Y0<=5e-9, Z0<=.25, Z2<=5e7 are unproved. Exact compatibility, physical
matching, semiclassical validity/noise and all-sector stability remain
unresolved. No story checkbox is promoted; the closed audit is preserved.

Remaining to-do list: tame common-regulated quantum response, seed contact
and nonlocal bounds, global data and actual controlled correction; full
nonlinear admission and contraction, exact compatibility, physical
matching/causal realization, validity/noise, all-sector stability,
charge/return, opening, transport, preparation and separate creation.
