# Paired boundary cancellation and the unresolved regulator identification

2026-10-23, E003/S169. Ordinary argument and execution protocol.
Previous turn: progress, finite-mode geometric conversion/contact cancellation.
The full stable-wormhole objective is active and incomplete. The feasibility
audit remains closed. This milestone addresses the domain prerequisite of
the remaining quantum/contact projection, before interpreting a number for it.

## 1. The first unsupported step, stated precisely

The retained actual-kernel construction gives a candidate weak quantum form
on the physical core and its form closure. Its separated kernel, constant
and second contacts, geodesic fourth restoration, global returns and paired
normalization are essential; the existing form inverse is not being erased.
See current-geodesic-two-sided-remainder.md and the ordinary arguments it
cites. Those estimates are much stronger than a finite-mode check.

There is nevertheless a distinction between (i) extending that form by
continuity and (ii) identifying its extension with two derivatives of the
original commonly regulated action on a family of mass/log metrics. For
the latter one must justify the simultaneous approximation/regulator limit
and metric differentiation, including the contact terms. A bound already
stated for the limiting form does not, on its own, bound the regulated
forms uniformly. Nor does convergence on a smooth core imply uniform
convergence on its closure. The raw thermal ODE representation additionally
needs its full matching to the auxiliary mass modulus/common subtraction.

The first unsupported step in the proposed numerical implementation is
this UNIFORM identification. The current-canonical-heat-restoration-assembly
note explicitly leaves its differentiated paired heat remainder unknown.
The finite-M conversion of the preceding milestone does not fill that gap.
No complete projection value is interpreted here. In particular, this note
does not retract the scoped weak-form estimates or promote them to a C2
action theorem on a nonlinear neighborhood.

## 2. The actual domain and a boundary obstruction test

Keep the original folded pair, complementary first-order adjoint domains,
regular thermal bolt and scalar-identical deformations. The original D is
not self-adjoint; its positive modulus uses both partner squares. Scalar
multiplication preserves the first-order projector. The induced square has
the retained mixed condition, with S=-K Pi_+/2 in inward heat conventions.
No independent Robin parameter or wall is introduced.

On the limiting physical side, Q'(0)=0 and b(0)=b_z(0)=0. Thus K=0,
S=0 and the tangential covariant derivative of chi=Pi_+-Pi_- is zero.
These statements hold for the mass/log phase, not just b identically zero.
For the rank-two geometric orbit square,

    tr chi=0, chi^2=1, E=-R I/4,
    Omega_an=(R/4) gamma_a gamma_n,
    chi Omega_an chi=-Omega_an.                     (1)

There is one tangential direction, so Omega_ab=0. At K=S=0 the a1 and
a3 boundary traces vanish. At a4 the bulk divergence has fiber trace
-6 Delta R. Its integral is +6 int_boundary R_;n for the INWARD normal.
The boundary derivative trace is

    tr[(240Pi_+-120Pi_-) E_;n
                    +(42Pi_+-18Pi_-) R_;n] = -6 R_;n. (2)

The two terms cancel before a limit or a variation. The remaining bulk
coefficient is -R^2/120, giving exactly the retained integrated A2.
Dropping either divergence term separately would be wrong.

Execute the next odd-coefficient test using the complete mixed a5 formula.
All surviving terms in its A1 part contain chi multiplying a scalar and
have zero trace. Its A2 part vanishes with the second fundamental form.
After S, chi_:a and Omega_ab are set to their actual zero values, the
remaining A3 trace is

    tr[-180 E^2+180 chi E chi E
       -45 Omega_an^2+180 chi Omega_an^2
                         -45 chi Omega_an chi Omega_an] = 0.  (3)

The first pair cancels because E is scalar; the last three cancel by
(1). In particular the potentially ill-defined R_;nn term is multiplied
by an identically zero fiber trace. Cancel this algebraically before
evaluating any fourth metric jet. This is not assigning a value to infinity.

The coefficient formulas are those of Branson--Gilkey--Kirsten--Vassilevich,
[Theorems 2.1 and 2.4](https://arxiv.org/pdf/hep-th/9906144). Their theorem
assumes smooth coefficients and boundary. We use its local algebra to test
candidate counterterms, NOT as a heat-remainder theorem for our log phase.
For every angular mass, multiplication by exp(-M^2 t) gives

    a5,M=a5-M^2 a3+(M^4/2)a1=0.                    (4)

This is an ALL-angular identity. No selection of low modes was substituted.
It removes this candidate M^-3 boundary term from the smooth local
coefficient ledger. It does not exclude nonclassical powers/logarithms in
the actual remainder, prove M^-4 decay, or justify regulator removal.

## 3. Execute the mass/log regularity and cutoff test

Write the leading physical-side shape as

    h(z)=c z^3+d z^4 log z,       v(z)=-z^4.        (5)

The actual phase also admits finite log polynomials and a remainder with
the retained open-collar conormal derivative bounds. Formula (5) isolates
its two leading mechanisms; it is not a replacement for the global h1.
Since Q(0)>0, change to proper distance preserves these local regularity
orders. On the physical side

    h'''=6c+d z(24 log z+26),
    h''''=d(24 log z+50).

For d nonzero the fourth derivative is unbounded, although square
integrable. On the actual even double, the cubic has fourth distributional
derivative 12c delta_0. Thus a nonzero mass cubic is in H3 locally but
not H4 across that seam. The log term alone is not a bounded C4 field.
Neither fact is a divergence of the renormalized weak quantum pairing.

For any 0<delta<1 the physical-side leading moments are EXACT:

    int_0^delta (h''')^2 = 36c^2 delta
       +cd delta^2(144 log delta+84)
       +d^2 delta^3[192 log^2 delta+288 log delta+388/3],
    int_0^delta (v''')^2 = 192 delta^3,
    int_0^delta h''' v''' = -72c delta^2
                              -d delta^3(192 log delta+144). (6)

Finite powers of log have the same integrability after three derivatives.
Consequently this endpoint supplies no failure of the weak H3 requirement.
At the bolt use the existing regular disk fields; there is no new cutoff.

Test the often implicit approximation argument quantitatively. Define
theta(t)=B(8,8)^-1 int_0^t x^7(1-x)^7 dx on [0,1], extended by 0 and 1.
It is C7, monotone and has seven matching derivatives at the joins.
On the physical half let h_e=z^3 theta(z/e-1). It is zero through e,
equals the cubic from 2e onward, and obeys the same first-order boundary
projector. It is only an APPROXIMATING TEST, not a substituted geometry.
With x=z/e,

    ||h_e'''-6||^2_(0,2e)=K3 e,
    ||h_e''''||^2_(0,2e)=K4/e,
    K3=36+int_1^2 [d_x^3(x^3 theta(x-1))-6]^2 dx,
    K4=int_1^2 [d_x^4(x^3 theta(x-1))]^2 dx >0.     (7)

For the quartic test the analogous third-derivative error is K3v e^3.
Compute K3,K4,K3v exactly and verify their scaling independently. Equation
(7) proves H3 convergence but an UNBOUNDED higher-norm approximation cost.
A heat theorem whose remainder constant requires that uniform H4/C4
control cannot be imported by this approximation. This is a failure of
that proof strategy, not an obstruction to the physical construction.

The candidate geodesic form has its sixth-moment bound proportional to
sum d_n/M_n^4 <=4C^2/q; its lower-order absolute kernel/contact estimates
also give continuity in the retained stronger form norms. This explains
why a weak extension can exist despite (7). To identify it with the
original determinant derivative, a regulator-uniform bound in a suitable
weak norm is still required. Hart Smith's
[finite-regularity heat theorem](https://arxiv.org/pdf/1809.05614) concerns
potential perturbations on a fixed complete manifold. It does not cover
our changing metric, paired boundary and uniform angular/metric limits.

## 4. Quantitative consequences and classifications

Let r_b(t)=gamma R_b(t)/(2pi) be the time-reduced geometric heat remainder
AFTER the same common a0/a1/a2 completion. Let r_hv^prod(t) include its
product-gradient contact. A sufficient uniform endpoint estimate is

    |r_hv^prod(t)| <= t omega_hv(t),
    D_hv(t0):=int_0^t0 omega_hv(t) dt/t < infinity. (8)

Uniformity means independence of the allowed field approximation and
common ultraviolet regulator, for the required variations. Together with
actual positive-time trace convergence and large-time massive decay,
this would justify dominated removal. For C2 differentiability of an
action family the corresponding 0th/1st/2nd estimates must be uniform
on that family; an estimate only at the selected pair is not enough.

Using the complete Theta(t)<=4C/t and the SAME physical normalization,

    E_cap <=4GC D_hv(t0)/(B^2 F_lambda),
    D_hv(t0)<=1.0666666e-11                       (9)

is sufficient for the reserved 2.5e-43 cap error. Evaluate with the actual
interval parameters. For t0=.001, substitute both

    omega=A t^alpha, D=A t0^alpha/alpha,
    omega=A_log t(1+|log t|),
           D=A_log t0(2-log t0).                 (10)

The latter explicitly allows a logarithmic remainder. These constants
are UNKNOWN, not measured failures. Finiteness in (8), rather than the
small precision target (9), is what first matters for the domain argument.
The a5 algebra does not supply either constant or imply alpha=1.

Also quantify a KNOWN excessive estimator. The existing row-subtracted
global gluing allowance is 6e24 C ||h||_2 ||v||_2 in disk action units.
The retained global |h1|<1.2e-8, z_h<1 and v=-z^4 give a sufficient
product bound 4e-9. The normalized allowance therefore is

    U_glue=12e24 G C (4e-9)/(B^2 F_lambda).       (11)

Compare it to 1.287e-40. It is an UPPER ALLOWANCE for one estimator,
not a value or lower bound on the remaining projection. The maximal
gluing coefficient K in C K ||h||_2||v||_2 that this crude norm pairing
could certify is 1.287e-40 B^2 F_lambda/(8e-9 G C). This coefficient
test changes neither normalization nor construction. Even improving
the entire cold/scalar allowance cannot repair that estimator loss.

## 5. Decision boundary

Execute (1)-(7) and (9)-(11), with frozen sources, all attempts retained,
and independent quadrature/arithmetic checks. Do not interpret successful
checks as a proof of the unsupported uniform continuum step.

Result classification fixed before execution:
- No boundary counterterm obstruction is detected through a5: its
  all-angular candidate coefficient cancels before any sum.
- A higher smoothness approximation is INAPPLICABLE to the actual mass/log
  class. The weak endpoint moments are finite; they do not close (8).
- The existing gluing estimator is an EXCESSIVE KNOWN upper allowance.
- Regulator-uniform mixed differentiation, the thermal representation
  match, full signed remainder and certified error remain UNKNOWN.

Research decision: retain the construction and sharpen the method at its
continuum prerequisite. Use the actual paired geometric Duhamel formula,
including the second vertex, and keep the common local contacts inside
the calculation. Park further adiabatic/finite-mode numerical projection
work until the uniform identification has been proved. No new cold or
scalar optimization is selected. No construction switch is justified by
the present evidence, and no controlled solution is established.

The SINGLE next calculation is a regulator-uniform estimate of the
product-subtracted mixed Duhamel remainder on the actual h1,-z^4 pair,
with the mass/log endpoint terms retained. Test a bound of the logarithmic
form (10), computing A_log rather than assuming it. The first question
is finiteness/uniformity; then test (9). If it holds, evaluate the complete
remaining signed projection with all thermal/angular/spatial errors and
test |P_hat_rest|+E_rest<1.287e-40. That projection was NOT executed here.
The target is sufficient screening, never a necessary existence condition.

Preserve the already controlled magnetic <1.284506e-40, complete cold
<1.693010e-43 and scalar-trial <2.331e-44 components. They do not establish
the remaining source or an actual corrected solution. The exact nonlinear
requirements at rho=1e-8 remain

    Y0+1e-8 Z0+5e-17 Z2<=1e-8,
    Z0+1e-8 Z2<1,

with sufficient Y0<=5e-9, Z0<=.25, Z2<=5e7 still unproved in the actual
admission norm. Exact compatibility, physical boundary/matching and causal
energy balance, semiclassical validity/noise and all-sector stability stay
open. The stationary boundary-work cancellation is not a retarded boundary
condition. The full stable-wormhole objective is unchanged.

Remaining to-do list: establish the actual-domain regulator-uniform mixed
remainder; then certify the full quantum/contact projection and global
forward/adjoint correction with admission; nonlinear closure and exact
compatibility; physical matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
