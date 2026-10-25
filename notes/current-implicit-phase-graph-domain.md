# Coupled graph attempt: exact phase prerequisite, global estimate unresolved

2026-10-24, E003/S169. Ordinary argument before supporting verification.
Previous turn: NO PROGRESS on the requested strong feasibility requirement.
The full stable-wormhole objective and closed feasibility audit remain intact.

Outcome: the FULL coupled estimate remains UNRESOLVED. An attempted numerical
boundary implementation candidate fails a necessary continuum requirement: substituting
a rounded slow exponent produces an infinite fast-graph defect. A conforming
implicit-root boundary primitive is now implemented. This removes that
implementation defect, not the unknown interior/phase inverse or any source
or contraction requirement. The earlier exact-root mathematics is unchanged.

This is NOT completion of the requested coupled-return milestone. The two
candidate implementations tested below were considered while constructing
a global approximate inverse; no prior computation is alleged to have used
a rounded exponent. No global approximate inverse has yet been assembled.

## 0. Complete system attempted, and where its quantitative argument stops

Retain the exact reference mass lift, with J=int_0^z t^3/Q(t)^2 dt,

    n_m=-2 beta QJ/z, b_m=beta z^3/Q+beta Q'J.

Its uncut areal fields f=2xz and alpha=0 cancel BOTH reference constraints
locally. Its actual cutoff is one through .1 and zero by .2, and gives

    f_cut=2xz chi_m-2F n_m Dchi_m, alpha_cut=n_m Dchi_m.

The transition cannot be discarded. It has zero reference bolt displacement
and fixed temperature because it vanishes near the regular bolt; its net
reference mass is -xB. The original reference matching is

    rho=-(F_l/2G)(rf)', rho+p=(F_l/G)rF alpha',
    s_R=[(G/F_l)rho(R)/R-kappa alpha(R)]/[F''(R)/2+kappa/R],
    delta_mu_infinity=kappa R s_R+(G/F_l)int_R^infinity rho dr.

These identities, including cutoff terms and the regular bolt, are retained
from current-coupled-global-return-attempt.md. They give no quantum zero
mode. Changing to this mass basis costs ||d_chi||H3>.00906808164 per unit
mass in the ORIGINAL fast graph, so retaining both 1e-11 graph ceilings
requires |m|<2.20553814952e-9. The conditional affine mass allowance
4.2953e-8 is neither the actual mass nor proof of that comparison.

The full attempted reduction, with Hbar=8piG A_r^2 H, is still

    Hbar=[[-Abar,-Cbar],[-Cbar*,Dbar]],
    K=Cbar^-1 Abar(Cbar*)^-1, S=I+K Dbar,
    S T=-Cbar^-1 f+K g, nu=(Cbar*)^-1(Dbar T-g).

Dbar includes the complete common-restored shape response, every angular
and thermal mode, contacts and finite returns. It is not replaced by B(D).
The known mixed inverse <12.500001 is a physical L2 estimate, not the
required graph/phase estimate. The FIRST UNSUPPORTED STEP in this ordinary
global argument is quantitative control of S^-1 and reconstruction in the
original strong graph, including the independent mass, scalar quartic,
slow phase, interior/bolt and fast inputs. Qualitative fixed-seed existence
does not supply it. The reviewed local symbol and smooth-boundary logarithm
theorems in current-coupled-schur-graph-attempt.md do not apply as numerical
global inverse theorems for these data. No new application is made here.

Attempting a conforming P first requires an exact slow border. Sections
1-3 execute that necessary representation calculation. They do not address
the missing strong norm of Dbar S^-1 or its original finite-part rows.
The actual nonlinear-family domain and operator comparisons remain separate
unproved requirements; this note only uses the retained fixed-seed scope.

## 1. The actual domain test preceding a numerical inverse

Keep both original rows, the actual paired mass/log domain, the common
regulator and fixed finite seed. The complete angular/thermal response,
coincidence contacts and finite returns remain in H. Its retained outer
normal operator has, with s=4+e and q=s(s-3),

    B(s)=-K(s-4)q[(q+2)H(s)-L],
    H(s)=psi(s-1)+(pi/2)cot(pi s/2)+kappa_M,
    K=G lambda beta/(15pi), L=10pi Fstar/(G lambda).

These are the full common-restored normal coefficients, not a selection
of modes. Set h(e)=psi(3+e)+kappa_M+(pi/2)cot(pi e/2)-1/e. Then exactly

    R(e)=6-L e+5e+e^2+(6+5e+e^2)e h(e),
    B(4+e)=-K(4+e)(1+e)R(e), B(4)=-24K.              (1)

The implicit e_* is the unique small positive zero of R. In the present
parameter intervals, L>3e34, 7.8<kappa_M<10.8. The convergent
[digamma series](https://dlmf.nist.gov/5.7.E6) and
[cotangent partial fractions](https://dlmf.nist.gov/4.22.E3) give, for
0<=e<1e-30,

    h(e)=psi(3+e)+kappa_M-pi^2 e/12+r(e), |r(e)|<=e^3,
    |h(e)|<12, |h'(e)|<2, |R'(e)+L|<120.             (2)

For the derivative bound use psi'(3+e)<1/2 and the absolutely convergent
cotangent remainder derivative <3e^2. No boundary regularity theorem is
being inferred from these scalar series. With eta=6/L, the correlated
endpoints eta(1 +/-300/L) have opposite signs in (1); (2) gives uniqueness.
The implementation retains THIS exact root, with an enclosing interval.

In the original source graph sigma=4.005, with the KNOWN affine g0 z^4
already subtracted. This is the fixed split used for the proposed return,
not an augmented source space with a freely variable leading coefficient.
A nonzero residual coefficient
r_hat z^s_hat, with real s_hat<sigma, has infinite weighted L2 norm and
therefore infinite H3 norm:

    int_0^a |z^-sigma r_hat z^s_hat|^2 dz/z=infinity.  (3)

Making r_hat very small does not change (3). This is a domain condition,
not a numerical error allowance.

For a cut slow column (nu,b)=chi(-a(s_hat),1)z^s_hat,
a(s)=(s+3)/(3(s+1)), chi=1 near zero and zero before the bolt, the trace
normal row cancels exactly, while the differentiated shape row has leading
coefficient B(s_hat). The retained actual-minus-normal kernel gains two
powers locally; the original remote z^4 finite part is killed by the
differentiated row and must still be kept in the original matching row.
The remaining terms cannot cancel the nonzero leading power in (3).
Specifically, current-global-boundary-moment-return.md equations (1)-(2)
retain A(m e_m), Rw and [A0,chi]w and give a differentiated remainder
O_b(z^5 log^N z). Since s_hat<4.005<5, this is o(z^s_hat).
For the explicit cut-power input no weak-solution bootstrap is needed:
its jets are explicit; the same retained separated-kernel expansion and
local two-power gain apply directly. The separated commutator is kept,
including its original z^4 moment. This uses the fixed one-sided analytic
seed and completed common contacts, not a smooth double or a theorem on
an arbitrary mass/log nonlinear neighborhood. The remainder constant is
not needed to reject the leading exponent; it IS needed for a finite
global defect bound and remains unknown here.

Consequently a preconditioner that inserts such an approximate exponent in
an independent slow input does not have a bounded defect in the ORIGINAL
fast source graph. This rejects the whole proposed operator-norm certificate,
even if every other column were exact. It is not a finite-mode substitute
for a continuum estimate: one allowed input suffices to disprove boundedness.
The field still belongs to the broader physical mass/log variation domain;
failure of this numerical graph construction is not failure of the action.

## 2. Executed tests of two approximations

The actual offset is e_*=1.8650970113308049749...e-34. Binary64 addition
4+e_* is exactly 4. The independent rounded slow column then has

    B(4)=-8.289319952534124466...e-35,                 (4)

which is nonzero. For the forced correlated pair
B0(z^4-z^(4+e_*)), B0=530514.97074458..., this rounding collapses the
pair and loses its known differentiated source coefficient

    g0=B0 B(4)=-4.3976083321111487851...e-29.

The leading missing term alone, on a finite interval z_min<=z<=.05,
has the exact derivative norm

    |g0|(1+delta^2)^(3/2)
       sqrt[(z_min^(-2delta)-.05^(-2delta))/(2delta)], delta=.005.

It is <1.3900e-26 at z_min=10^-300, yet >1.3906e-11 at z_min=10^-3300.
It reaches 1e-11 near log10(1/z_min)=3271.35343347. These are norms of the
mandatory LEADING TERM, not enclosures of the full finite-cut residual.
They show why a finite endpoint test can conceal the infinite-domain defect.
Equation (3), not the diagnostic cutoff, proves the continuum failure.

Even retaining the analytic leading approximation eta=6/L is insufficient.
Preserve 6-L eta=0 before interval evaluation. Then

    R(eta)=5eta+eta^2+(6+5eta+eta^2)eta h(eta)>0,
    B(4+eta)=-1.850453702...e-67,
    e_*-eta in (4.16352087e-67,4.16352090e-67).        (5)

The last bound follows from R(eta)/(L+120)<=e_*-eta<=R(eta)/(L-120).
This calculation preserves correlations even though separately stored
intervals for e_* and eta overlap. (5) again produces (3). Increasing
precision or taking a finite root-iteration count is not by itself a domain
proof. For any approximate root with certified nonzero B, the permitted
leading-coefficient error in this fast graph is EXACTLY ZERO.

## 3. Executed replacement, with both canonical components

The new operators/exact_slow_phase.py defines e_* implicitly by (1)-(2)
and evaluates the correlated pair without subtracting large rounded fields.
For x=log z, p=B0 e_*, E=exp(e_* x), and

    Q_j(e)=sum_(k=1)^j binom(j,k)4^(j-k)e^(k-1),
    z^-4 D^j b=p[-4^j expm1(e_* x)/e_*-Q_j(e_*)E],
    z^-4 D^j nu=-(7/15)z^-4 D^j b
                    -2p(4+e_*)^j E/[15(5+e_*)].     (6)

These are identities, not a truncation in e_*x. In particular they remain
valid at arbitrarily small z, where a fixed log-Taylor approximation would
also fail the fast-graph requirement. The interval is an enclosure of an
exact root, never a point exponent substituted for it. The exact identities
theta=nu+a(D)b=0 and B(D)b=g0 z^4 are used symbolically for leading rows.
This use of B(4+e_*)=0 follows from the proved unique-root definition; it is
not a small interval residual reset to zero.

The implementation evaluates (6) and the independent column
(-a(4+e_*),1)z^(4+e_*) through eight conormal derivatives,
including x=-1e35; it avoids forming 4+e_* in binary64 and avoids cancellation
of B0-sized amplitudes. It deliberately leaves the ORIGINAL shape finite
part and cutoff/global returns unevaluated, rather than returning zero.

All other original rows remain necessary. The exact reference mass lift,
its cutoff transition and regular-bolt matching are unchanged; it is not a
global fixed-temperature zero mode. A scalar quartic still contributes
20 beta Fstar c4 z^4 to the original second row. Neither it, the mass row,
nor a weight-sigma fast remainder cancels the independent erroneous slow
power in (3). The corrected phase primitive also does not evaluate their
couplings. Multiplying (6) by its retained cutoff requires ALL product-rule
and nonlocal commutator terms in the actual global defect. None is charged
as zero by the leading identity.

The quantitative contribution to the full preconditioner is specific:
the inadmissible leading slow-column defect changes from an infinite graph
norm to exactly zero. Its regular remainder, every global column and the
complete defect norm are still UNKNOWN. The primitive is not a full P.

## 4. Decision and remaining feasibility conditions

Retain the construction and exact mixed Schur algebra. Reject rounded-power
and finite log-Taylor boundary implementations; use the implicit-root pair
(6) with the original finite-part matching equations. This is a change in
numerical representation, not in the physical domain or admission target.
No global or nonlinear requirement is marked complete.

The next SINGLE calculation is the complete graph defect of a global P
using this conforming border, including cutoff/interior/bolt return and
both original matching rows. It must keep all angular/thermal tails and
common contacts, and establish delta=||I-HP||<=1/2 together with the
physical output maps and joint source/feedback errors. That full calculation
is PROPOSED, NOT EXECUTED. The boundary primitive and failures (3)-(6) ARE
executed/proved within the stated seed scope.

No kappa_m or other post-elimination column is measured. At conditional
kappa_m=5e-5, residual budgets remain 5.34235e-12 affine and .0012069 graph
feedback before independent energy feedback and every other entry. The
independent mass-energy debit is .0195239 times the energy/graph increment
ratio: allocating ALL .0012069 to it would require a ratio below
.06181654280138. Other entries make that restriction stricter. The scalar
quartic, independent slow phase and fast/interior cross entries have no
certified values; their budgets cannot be presumed available twice. The
original mass conversion still requires |m|<2.20553815e-9 to retain both
graph ceilings. Actual m is unknown. The fast screen
.028913+1e-11(1.0706e10)=.135973<.1599164 retains only its two-fast-input
admitted-family scope. It supplies no independent phase derivative.

For the proposed global P, E=I-H P, delta<=1/2 and M_T=||T P|| must be
certified on the ORIGINAL spaces. The complete requirements remain

    k_m_trial+2 M_graph ||E b_m||<=5e-5,
    S_trial+2 M_graph epsilon_aff<=7.49e-12,
    eta_trial+2 M_graph epsilon_fb<=.00368,
    |m_trial|+2 M_m epsilon_aff<2.2055381495e-9.

Here epsilon_fb includes all independent phase and interior inputs, and
T is an original physical output map. No delta, M_graph, M_m, trial column
or aggregate defect has been computed. Thus the conditional residual
targets in these inequalities are not measured coefficients. The new
boundary representation removes one infinite term from E, but gives no
finite upper bound on E or any of those rows. The already known loose L2
upper estimate is not an upper estimate of these different graph maps.

The original comparisons and nonlinear inequalities remain unproved:

    E_center+1e-8 alpha<=1e-11,
    R_X alpha*.028913+Z0_other<=.25,
    R_X alpha^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

No actual controlled solution, exact compatibility, physical matching,
semiclassical validity/noise or all-sector stability follows. The complete
graph estimate remains unresolved; the failed numerical realization is
not an obstruction to the construction. Supporting arithmetic checks do
not substitute for (1)-(3) or the missing global argument.

Remaining to-do list: complete coupled graph defect and global phase/bolt
returns; actual correction and original physical-norm comparisons; nonlinear
admission/contraction and exact compatibility; physical matching/causal
realization; semiclassical validity/noise; all-sector stability; charge/return,
opening, transport, preparation, separate creation and the full objective.
