# Moving-comparison Cauchy test: the coupled return remains unresolved

2026-10-25 HST, E003/S169. Ordinary argument before supporting arithmetic.
Previous goal turn: NO PROGRESS; sources were inspected, with no new bound,
file change or commit. The full stable-wormhole objective remains active.
The closed feasibility audit is unchanged.

OUTCOME: UNRESOLVED complete coupled estimate. The attempted term-by-term
holomorphic estimate for the moving short comparison fails on the ACTUAL
domain. A real finite-jet alternative is quantified below for its cutoff
factor only. No physical Hessian column or strong composite was evaluated.
This is evidence from an incomplete attempt, not completion of the requested
global feasibility milestone or an obstruction to the construction.

## 1. Which missing argument was attempted

The missing numerical step is still equation (7) of
current-short-travel-strong-domain.md: differentiated, common-restored
Dirac-minus-short-geodesic kernels AND coincidence contacts, with the
finite-domain matching and strong inverse. Its auxiliary weight endpoint
is repaired in current-strong-return-weight-repair.md. Retain the weights
3/2 ->13/4 ->4.005 and the original two-row quantum composition factor
<1481.797. Neither that factor nor the H3->L2 bound 1e48 C supplies the
missing higher output derivatives or global matching constants.

An attempted shortcut was to bound the moving remainder in a complex
metric tube, and use Cauchy inequalities for its high derivatives. A
term-by-term version would require the moving short comparison itself to
be holomorphic into fixed Sobolev operator spaces. That hypothesis fails.
This shortcut was proposed in this continuation; it is NOT an assertion
that the previous real finite-jet argument used Cauchy inequalities.

There is also no inherited all-order metric theorem: section 5 of
current-curved-paired-heat-bridge.md proves TWO metric derivatives of the
proper-coordinate identification on the admitted mass/log family. It
explicitly does not complexify the inverse proper-distance map or provide
higher source variations. Analyticity of the heat coefficient in an
independent f cannot be silently transferred through that identification.
We do not claim here that the physical Hessian H is nonanalytic.

## 2. Exact failure on admissible interior variations

Keep the actual paired mass/log realization, original complementary
projectors, common regulator, thermal period and regular bolt. In an
ordinary interior radial chart write the seed metric

    g0=dY^2+A(Y)^2 dt^2.

Choose the chart a positive distance from both ends, containing a radial
interval of radius 4a, a=1/128. Such regular charts are part of the retained
seed construction; short geodesics have curvature <1.5 and no conjugate
point on r<=2a. No extension through the C2,1 seam is used.

Let eps=a/100. Choose nonnegative, stationary smooth probes u and v,
positive inside two radial intervals of radius eps whose centers are
-(a+eps) and +(a+eps). Their supports have distance exactly 2a. They vanish
near the physical boundary and bolt, so are in the original domain with
zero mass, quartic and slow boundary coefficients. They need not be
analytic spatial functions: these fixed smooth domain probes are legitimate
continuous tests on every finite Sobolev space in question.

Choose smooth stationary chi=1 on |Y|<=3.5a and chi=0 for |Y|>=4a, and

    g_t=exp(2t chi)g0.

For sufficiently small real |t| this is an admitted interior family. All
endpoint data, the period and regular-bolt matching are unchanged. The
large but finite derivatives of chi only reduce the allowed interval in
t. This is not a global rescaling of temperature or a replacement seed.

Any path of length <=2a starting in the probes stays inside chi=1 for
small enough |t|: the distance to its edge is at least 2.48a at t=0.
Within that strip all lengths scale by exp(t), and the shortest radial
distance between the supports becomes 2a exp(t). Leaving the strip cannot
produce a shorter path within this travel range. Thus for t>=0 the short
off-diagonal pairing is exactly zero. For every sufficiently small t<0
an open set of pairs in the supports is connected by geodesics of length
<2a; an open set of nearby initial angles has the same property.

Precisely, suppressing the fixed positive action factor, the short
comparison's off-diagonal shape pairing is

    B(t)=-<J_t u, integral_0^(2a) zeta(r) w(r)
                                      cos(rP_t) dr J_t v>.       (1)

The base lift J_t preserves positivity. Each cos(rP_t) is the half sum
of forward and backward geodesic transport; both preserve positivity.
The COMPLETE angular weight is

    w(r)=r sum_(n>=1) d_n M_n^4 [K1(M_n r)^2-K0(M_n r)^2]/(2pi^2)>0.

The positivity and w<=72C/r^5 are the retained entire-tower bounds, not
a finite mode test. Here r stays away from zero on the support of the
pairing, so absolute convergence also follows immediately. The same
argument works if an analytic positive mass/action factor varies with t.
Since zeta>0 for r<2a, (1) is strictly negative for small t<0 and zero
for small t>=0. All identity, differential polynomial and coincidence
contact terms pair to zero between these disjoint supports. No contact
was deleted from the physical operator. Actual thermal and finite-domain
returns remain in H; they are not bounded or cancelled by this test.

Consequently B is not real analytic at zero. Indeed an analytic function
zero on a one-sided open interval is zero near zero, by its convergent
Taylor series, contradicting the strict sign on the other side. The
smooth local flow and smooth cutoff do give finite real derivatives.
If the moving short comparison were holomorphic as a bounded operator
between any of the proposed fixed finite Sobolev spaces, evaluation on
these fixed probes would be analytic. One can use a fixed seed density:
the moving volume factor exp(2t chi) is an analytic bounded multiplier,
so identifying the pairings this way does not alter the contradiction.
Hence no such termwise holomorphic
estimate exists, even on this interior subfamily of the actual domain.

The same cutoff prevents a holomorphic spatial-kernel shortcut across
travel distance 2a at the fixed seed. Smoothness to finitely many orders
is sufficient for the intended real commutator method; analyticity is
not required by that method. This argument says nothing about failure
of the actual strong map or of wormhole existence.

Freezing the comparison at the seed is not a repair of the proposed
uniform lower-order estimate: H(g)-H_s(g0) then retains the variation of
the fourth-log principal symbol. For a conformal interior direction the
principal |xi|_g^4 log|xi|_g factor changes. Without a further moving
principal subtraction it is no longer the order-(2+epsilon) remainder
used in the strong-power bootstrap. At the seed the frozen difference
is the same; this observation only concerns the proposed family shortcut.

## 3. Quantified real finite-jet replacement for the actual cutoff

Keep a fixed. On the transition put x=r/a-1 in [0,1],

    theta(x)=exp(-1/x), theta(x)=0 for x<=0,
    S(x)=theta(x)+theta(1-x), zeta=theta(1-x)/S(x).

For 0<x<=1, theta's elementary complex germ is analytic on the closed
circle |z-x|=x/2; it excludes zero and Re(1/z)>=2/(3x).
[Cauchy's derivative formula, DLMF 1.9.31](https://dlmf.nist.gov/1.9.E31)
therefore gives, for j>=1,

    |theta^(j)(x)| <= j!(2/x)^j exp[-2/(3x)]
                    <= j!(3j/e)^j <= D_j=3^j(j!)^2.             (2)

The last inequality follows from log(j!)>=integral_1^j log s ds
>=j log j-j. Set D_0=1. Smooth extension at zero follows from the same
exponential domination of each finite derivative. This use of a scalar
germ is legitimate; it does not complexify the moving operator.

Since one of x,1-x is >=1/2, S>=exp(-2). Differentiating S zeta=theta(1-x)
gives the explicit positive recurrence

    Z_0=1,
    Z_j=e^2 [D_j+2 sum_(k=1)^j binom(j,k) D_k Z_(j-k)], j>=1.    (3)

It proves sup |partial_x^j zeta|<=Z_j for EVERY transition point, not
only sampled values. Outside the transition every positive derivative
vanishes. Original travel derivatives cost a^(-j) Z_j=128^j Z_j.
For Euler travel derivatives use the Stirling numbers of the second kind:

    T_0=1,
    T_j=sum_(k=1)^j S2(j,k) 2^k Z_k,
    sup |(r partial_r)^j zeta|<=T_j.                         (4)

This follows by expanding ((1+x)partial_x)^j. These two norms are NOT
interchanged. The note does not infer an original physical norm bound
from Euler derivatives without the missing geometric comparisons.

There is a direct quantitative role for (4) in the full cutoff ledger.
Let A_j be the actual full-flow multiplier obtained by replacing zeta
in the [1-cos(rP)] part by (r partial_r)^j zeta, j>=1, with all original
action/York factors retained. Spectral calculus, the complete w bound,
||JT||^2=2pi||T||^2 and the original form conversion give

    ||A_j||_(V->V*) <= delta_a T_j < 3.2e-18 T_j.             (5)

Indeed integrate 72C r^-5 over [a,2a], or more coarsely [a,infinity),
and use |1-cos|<=2. This is the cutoff contribution only: metric/space
differentiation also acts on w, flow, Jacobians, probes and contacts.
Those other terms, and their higher output norms, have not been bounded
by (5). The individual contact terms still have to be combined before
common-regulator removal and differentiation.

The supporting calculation evaluates (3)-(5) through order 16 using
512-bit Arb. Order 16 covers the canonical output order in the selected
eight-power argument, not a claim that no extra kernel difference jets
are needed. Equations (2)-(4) apply to any additional finite order.
All results are upper-bound expressions, not measured derivative norms.
No cutoff-radius optimization or physical operator change is made.

Executed enclosures, rounded UP from the retained Arb expressions:

| Order j | Original travel derivative 128^j Z_j | Euler travel T_j | Cutoff-only form bound 3.2e-18 T_j |
| --- | ---: | ---: | ---: |
| 1 | 8513 | 133.004 | 4.257e-16 |
| 4 | 5.551e16 | 3.358e9 | 1.075e-8 |
| 8 | 1.669e35 | 6.135e20 | 1963 |
| 16 | 2.777e74 | 3.772e45 | 1.207e28 |

The high-order costs are substantial even though the undifferentiated
transfer cost is tiny. They are finite real-jet upper bounds, not lower
bounds on an actual response, and are in different norms from M_graph,4
and M_X,8. It would be invalid to compare 1.207e28 directly with .0012069
or to infer a failure of the coupled budget. Conversely these numbers do
not prove that the missing global constants fit their tail allowances.

## 4. Coupled budgets and what the attempted calculation does not close

The cut reference mass lift remains exactly

    J=int_0^z t^3/Q(t)^2 dt,
    n_m=-2 beta QJ/z, b_m=beta z^3/Q+beta Q'J.

Its uncut reference return is f=2xz, alpha=0, cancelling both local rows.
For the SAME chi_m (one through .1, zero by .2), retain

    f_cut=2xz chi_m-2F n_m D chi_m, alpha_cut=n_m D chi_m,
    s_R=[(G/F_l)rho(R)/R-kappa alpha(R)]/[F''(R)/2+kappa/R],
    delta mu_infinity=kappa R s_R+(G/F_l)int_R^infinity rho dr.

The cut lift has s_R=0 and reference mass -xB; the cut scalar quartic
and correlated slow pair have zero net reference mass. The source and
reference matching are linear, so these identities retain their fast
couplings. They do not evaluate the quantum entries of the coupled
system. The uncut temperature derivative is not kappa_m.

Original-fast conversion remains >.00906808164 |m|. Requiring both
old and lifted graph norms <=1e-11 therefore still requires
|m|<2.20553814952e-9. The conditional coefficient allowance
4.2953e-8+49.462 E_graph is NOT the actual mass coefficient.

If kappa_m<=5e-5 is eventually proved, its conditional budget debits are

    affine: 4.2953e-8 * 5e-5 = 2.14765e-12,
    graph feedback: 49.462 * 5e-5 = .0024731,
    independent energy: 390.478 * 5e-5 = .0195239.

This leaves 5.34235e-12 affine and .0012069 graph. The scalar-quartic,
correlated-slow, fast, interior and bolt returns have NOT been shown to
fit. Even spending the entire remaining graph budget on the energy
term alone would require E_phys/E_graph <.06181654280138; the comparison
is unknown and other entries would reduce that allowance.

Retain eta_seed+1e-11 C3=.135973<.1599164 ONLY in the proved two-fast-input
admitted-family scope. It is not control of independent boundary phases,
an X-ball admission theorem or a global physical-norm comparison.

The full strong constants M_graph,4, M_X,8 and U_feedback are UNKNOWN,
not large evaluated bounds. The conditional 64-step graph allowances
M_graph,4<=6.0612e157 and M_graph,4 U_feedback<=1.5062e156 are NOT actual
coefficients. No 64-step solve, signed coupled center or application-error
certificate exists here. Equations (3)-(5) cannot be inserted for M.

In particular the original preserved-neighborhood requirements remain

    E_center+1e-8 alpha <=1e-11,
    R_X alpha*.028913+Z0_other <=.25,
    R_X alpha^2*1.0706e10+Z2_other <=5e7,
    Y0+1e-8 Z0+5e-17 Z2 <=1e-8,
    Z0+1e-8 Z2 <1.

No new term of these inequalities is supplied by a cutoff form bound.
No quantitative evidence here establishes whether the remaining strong
constants fit the allowances. Their finiteness or a large tail allowance
is not evidence that the complete signed response fits the budgets.

## 5. Research decision

Retain the physical construction. Sharpen the strong-return METHOD by
using real finite-jet kernel/commutator estimates, with (3)-(5) paying the
actual cutoff costs. Abandon the proposed termwise holomorphic moving-
comparison shortcut. Its failure is verified; feasibility of the complete
coupled estimate remains undecided. No estimate of the actual coupled
coefficient has been found to exceed the target.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): evaluate the complete real
differentiated composition H^-1 R_s in the repaired weighted spaces,
including both original rows, common-restored contacts, all angular and
thermal sums, phase/moment rows and finite-domain/bolt reconstruction,
and propagate it to M_graph,4 and M_X,8 with actual column energy norms.
Keep the original quantum prefactor attached. This must produce a numerical
global bound; another cutoff/contact-only estimate would not complete it.
The signed coupled solve and its error certification remain required.

Remaining to-do list: numerical complete strong return and coupled signed
response; original physical-norm comparisons; nonlinear admission and
contraction; exact compatibility; physical matching/causal realization;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole objective.
