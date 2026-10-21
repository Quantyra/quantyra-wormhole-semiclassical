# Moving-cut quantum comparison with density and derivative contacts

2026-10-21, E003/S168. Ordinary argument before focused verification.
The complete stable-wormhole goal remains incomplete.

For the retained smooth charged seed and each fixed Euclidean frequency,
the actual canonical shape Hessian differs from the density-conjugated
flat normal Hessian by a form with a finite translated-logarithmic weak
cutoff limit. This includes the density row and derivative contacts.
The statement is about the canonical quantum functional. The original
local Weyl, classical and LLL shape terms still have to be combined with
it before a full pure-shape limit or a conditional branch exclusion.

Keep the paired domain, period, angular tower, mass chart, subtraction
scale and finite coefficients. No finite counterterm is selected here.

## 1. Statement and allocation

Use the actual proper coordinate Y, inward from the limiting boundary,
with metric dY^2+a(Y)^2 dtau^2. On a fixed outer collar put

    a0=a(0)>0, m(Y)=sqrt(a(Y)/a0), u(Y)=m(Y)v(Y).

The seed is smooth as a one-sided function through Y=0, a'(0)=0, and
its fixed smooth continuation has uniformly bounded jets. Cap m to one
outside a slightly larger collar. The resulting far terms belong to
the already controlled global allocation. Localize v inside the collar;
the same allocation controls the localization error. The shape domain
is the actual auxiliary energy space, in particular

    v in H2(0,infinity), v(0)=v'(0)=0,
    v=Y^(3/2)V(-log Y), V in H2 near infinity.

As before, begin with shapes smooth at every positive Y. Limits below
extend their quadratic distributions; no pointwise v''(0) is assumed.
Weak convergence means testing Z=-log d against phi(Z-Z0) as Z0 tends
to infinity, for every compact smooth phi. It is a necessary convergence
test, not a prescription to average the physical action.

Let Qcan,d^gamma be the original canonical Hessian with its actual time
normalization divided out. Let Qflat,d be the complete zero-frequency
flat angular-tower form, with even profile u(d+|x|), defined by the
common meromorphic finite part in the [normal theorem](finite-cut-normal-shape-limit.md).
Then

    Qcan,d^gamma[v] - a0 Qflat,d[u] has a finite weak limit.       (1)

Its value is retained, not asserted to be zero. The established
[constant row](moving-cut-constant-quantum-row.md),
[global terms](moving-cut-global-quantum-locality.md), and
[reflection correction](reflected-boundary-kernel-limit.md) are part
of (1). The new steps are the smooth local kernel, its finite contacts,
and the density row. This argument does not replace the actual smooth
continuation by a reflection-symmetric background.

## 2. The smooth local direct and reflected kernels

For the direct smooth-continuation bubble, center the retained conformal
normal cap at the midpoint of the two physical points. Its first metric
jet is removed before the resolvent comparison. The estimates in
[the curved seed proof, sections 3 and 6](stationary-curved-quantum-tt.md)
give a two-dimensional error bounded by K C/r^4 after the common fourth
allocation, versus a principal kernel of order r^-6. The cap radii and
jet bounds are uniform in small d, since this is the fixed smooth seed
translated through a compact set of centers. The time integral therefore
gives a radial error bounded by K C/|Y-Y'|^3.

Here the density must be kept. Freezing time at the proper midpoint Ym
gives the factor a(Y)a(Y')/a(Ym). Its ratio to sqrt(a(Y)a(Y')) is

    sqrt(a(Y)a(Y'))/a(Ym)=1+O((Y-Y')^2).                       (2)

Taylor's formula for log a proves (2); the odd first jet cancels.
Proper distance and normal distance differ at cubic order. Multiplying
the full flat tower kernel, bounded by K C/|Y-Y'|^5, by either error
again gives order |Y-Y'|^-3. Thus the leading direct kernel is exactly
sqrt(a(Y)a(Y')) times the flat kernel in the proper radial coordinate.

For the crossing kernel, the continuation point is Y=d-x', while the
original positive-side mass vertex remains evaluated at d+x'. Their
relative change, and the corresponding reflected-coordinate change,
are O(|a'(d)|(x+x')+(x+x')^2). This follows by Taylor expansion about d;
a'(d)=O(d). Changing from the isothermal reflection to proper reflection
has displacement O(|a'(d)|(x+x')^2+(x+x')^3), with the same relative
effect on the principal kernel. The crossing remainder is consequently

    |Ecross,d(x,x')| <= K C [d/(x+x')^4+1/(x+x')^3].           (3)

Symmetrize the two crossings as in the preceding reflection theorem.
The original positive-side vertices are unchanged. Equation (3) is
within that theorem's corner envelope (even without its angular log).
After allocation of its constant row, its difference-square form and
its contacts have an ordinary limit. This also covers the bounded
smooth factors introduced by changing the normal coordinate.

The resulting weighted flat difference form is

    Dm,d[v]=-(1/2) FP integral_(x,x'>0) Kflat,ev(x,x')
                m(d+x)m(d+x') [v(d+x)-v(d+x')]^2 dx dx'.      (4)

All finite parts in (4) inherit the original common subtraction. They
are not defined by the separated-kernel bounds alone.

## 3. Why the direct derivative contacts converge

The common fourth contact and its highest boundary derivative are fixed
before comparing kernels. At finite angular cutoff the negative fourth
moment cancels the fourth Taylor subtraction; the actual bubble plus
restoration and the model then have local order at most two. Taking the
same angular finite part preserves their common principal coefficient.
In particular the leading boundary term is the same coefficient of
2a(d)v'(d)v''(d). Terms produced by differentiating m instead have at
most two derivatives on the shapes. The explicit curvature-gradient
contact in the curved seed proof is retained.

For the remaining direct second contact, use the fixed compact coordinate
quadratic at each center from that proof. Its subtracted one-mass
remainder is O(M^-3), uniformly in center and d, and sum d_n/M_n^3 is
finite. The global remainder is exponentially summable. Smooth dependence
on the seed and dominated summation give a bounded continuous full-line
anchor c_d(Y). Restricting this anchor to the positive half-line removes
a second moment with bound

    K integral_x^R ds/s <= K log(R/x), 0<x<R.                (5)

Therefore the actual half-line second coefficient is bounded by
K(1+|log x|), and converges for each x>0. It must not be replaced by a
uniformly bounded coefficient right at x=0.

To check that this is sufficient, let f(x)=v(d+x). Extend f to the line
with a bounded H2 extension. For example 3f(-x)-2f(-2x) on x<0 matches
both traces. Taylor's integral remainder gives

    integral |[f(x+s)-f(x)]^2-s^2 f'(x)^2| dx
       <= |s|^3 ||f'||_2 ||f''||_2 + (s^4/4)||f''||_2^2.     (6)

Restricting the integral to x,x+s>=0 does not worsen the bound obtained
from the extension. Equation (6) is integrable against |s|^-3 ds on a
fixed collar. The coefficient (5) pairs with f'^2: f' is bounded by
its H1 norm and log x is integrable. These bounds give a uniform
continuous quadratic form on H2 of the half-line. Translation is
continuous there. For a fixed smooth profile, kernel convergence off
the diagonal, (6), and the fixed-probe anchor give convergence. Uniform
boundedness and H2 approximation prove ordinary convergence for v.

The remaining corner contact orders are at most two. For a difference
square their only possible nonzero second jet is a multiple of v'(d)^2,
as in the reflection theorem. More general allocation contacts are
smooth bounded coefficients times v(d)^2, v(d)v'(d), v'(d)^2 or
v(d)v''(d). The first three tend to zero. For the last use

    2v v''=(v^2)''-2v'^2

and the retained cubic-test weak lemma. No pointwise second trace is
used. Their coefficients are fixed and continuous: smooth half-domain
heat parametrices and fixed smooth probes, with the same subtraction,
give a compact family of finite anchors. The derivative weights follow
from [Vassilevich, (2.21), (2.28)-(2.32), (5.33)](https://arxiv.org/html/hep-th/0306138v3).
The rough even double is not assigned a smooth-metric heat expansion.

## 4. The density row is a genuine finite contribution

Write B_d(f,g) for the polarized half-line flat form on the even profiles
f(d+|x|), g(d+|x|). Algebra at a finite regulator, followed by the same
meromorphic continuation, gives the exact identity

    Dm,d[v]=B_d(mv,mv)-B_d(m,mv^2).                         (7)

Indeed the difference between the two difference squares is
m_x(m_y-m_x)v_x^2+m_y(m_x-m_y)v_y^2. The second term in (7) is a density
row; it cannot be deleted by invoking constant annihilation. For the
zero-frequency flat form B_d(1,.)=0, so the capped m can be replaced by
m-1 in this row. Its finite limiting value generally is nonzero.

Set U=mv^2. The H2 algebra property gives U in H2; the actual energy
domain gives |U(Y)|<=K Y^3. All caps and subtraction lengths are fixed.
Put h_f(x)=f''(d+x), a_f=f'(d), b_f=f''(d). Polarizing the exact normal
representative gives

    J_d(f,g)=integral_0^infinity h_f h_g dx+a_f b_g+a_g b_f,
    L_d(f,g)=(1/2)<h_f,ev,log|D| h_g,ev>
                         -a_f T[h_g]-a_g T[h_f],             (8)
    T[h]=FP integral_0^infinity h(x) dx/x+gamma_E h(0).

The fixed principal multiplier is c(L_d-ell J_d), c=C/(30pi), with
ell the retained logarithmic scale and constant. Equation (8) follows
from the same Riesz continuation, including its scaleless delta-square
term, as the normal theorem. There is no independent contact choice.

For f=m and g=U, log|D| h_m,ev is uniformly in L2 and converges in L2.
At high frequency this follows from the uniform H1 bound on h_m,ev;
at low frequency its L1 bound and integrability of log^2|p| suffice.
Since h_U,ev converges in L2, the first pairing in (8) converges
ordinarily. Also T[h_m] is bounded and continuous, whereas U'(d)->0.

Only m'(d)T[U''(d+.)] needs a weak argument. With the finite part
defined as the constant Laurent coefficient, integration by parts in
the initially convergent Mellin integral gives

    T[F'']=2 FP integral_0^infinity F(x) dx/x^3
                                      +(gamma_E-3/2)F''(0). (9)

To see the contact, continue
integral x^(z-1)F''=(z-1)(z-2) integral x^(z-3)F.
The residue on the right is F''(0)/2; multiplication by -3z produces
the -3F''(0)/2 term. Formula (9) fixes its sign and normalization.
The Gamma integral and Euler constant are in [DLMF 5.9](https://dlmf.nist.gov/5.9).

The [cubic-test lemma](energy-domain-and-chart-contact.md) now applies
to U(d+x). The finite part has order three and the contact order two.
The extra smooth factor m'(d)=O(d) makes both contributions tend weakly
to zero. In detail, with a log-cut test at scale delta, the near-collar
bound is O(delta^2) up to fixed logarithms; away from that collar it is
O(delta) times the integrable function U(Y)/Y^3. Differentiating the
coefficient while transposing the contact has the same scale bound.
The J term in (8) has the same harmless m'(d)U''(d) contact, plus
ordinary bulk and vanishing U'(d)m''(d) terms.

Thus the principal density row has a finite weak limit, with finite
interior value retained. The full flat remainder E=Kflat-Kprincipal
obeys |E(p)|<=K(1+p^2)log(e+p^2). Its bilinear form is continuous on
H^(5/4)(R). The even extensions of m(d+.)-1 and U(d+.) converge in
that space, by the half-line H2 extension estimate in the normal
theorem. Its density row therefore has an ordinary limit. Together
these prove the needed weak limit of the complete second term in (7).

## 5. Fixed frequency and the consequence

For a fixed Euclidean frequency gamma, cos(gamma tau)-1 adds two powers
of time separation to the local bubble. The direct error is again of
radial order |x-x'|^-3 with the same second-anchor argument. Its corner
kernel has order at most (x+x')^-3. The associated constant row has
order at most x^-2 (possible logarithms included) and boundary contacts
of order at most one. Its pairing with v^2 has a finite weak limit by
the cubic-test lemma; the difference-square part has an ordinary limit.
The compact smooth fixed-probe family fixes these contacts as above.
The flat frequency difference Kflat(sqrt(p^2+gamma^2/a0^2))-Kflat(|p|)
also has order p^2 log p and is continuous on the even H^(5/4) profiles.
The global thermal terms retain the actual state and were controlled
separately. This proves (1) for each fixed gamma, not uniformly for
unbounded gamma or under retarded continuation.

The actual energy-domain transfer puts u in the normal logarithmic
domain. Consequently the canonical Hessian plus the NORMAL compensator

    a0 c [integral_d^infinity log(Y) u''(Y)^2 dY
                                  +2log(d)u'(d)u''(d)]       (10)

has a finite weak cutoff limit, up to its fixed constant J term. This
is a consequence of (1) and the normal theorem, not a replacement of
the original Weyl action by (10). Its original bulk and boundary terms,
and all classical/LLL shape partners, must next be compared with (10).
Only their combined limit could complete the pure-shape calculation
needed by the [nonlinear trace forcing](energy-shape-nonlinear-trace-forcing.md).

Verification below supports the density algebra, bilinear Fourier
continuation and its boundary contacts, the finite-part identity (9),
half-line logarithmic anchors and cutoff examples. It does not compute
the full charged two-point kernel or formally certify the continuum
argument. No nonlinear solution or physical stability follows.

Remaining to-do list: combine the original local Weyl, classical and LLL
shape terms with this canonical normal form; assess the stated candidate
expansion, then physical boundary/source compatibility and continuous
causal constraints; nonlinear existence, all-sector stability, exterior
and state consistency, charge/return, opening, transport, preparation,
separate creation and full audit.
