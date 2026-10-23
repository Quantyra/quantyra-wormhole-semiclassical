# Vacuum-contact budget and feasibility closeout

2026-10-23, E003/S169. Ordinary argument before execution.
The full stable semiclassical wormhole objective remains incomplete.

This is the last budget substitution for the closed audit, not another
source-estimate chain. Decide whether computing the angular vacuum finite
part more accurately could change the present local integration decision.
Keep the current action, original G, q=2000^4, c_R=1e22, mu0=1/b0,
temperature, paired Dirac domain and physical observable unchanged.

## 1. Use the available bound before commissioning a spectral calculation

The [common angular subtraction](stationary-reference-nonbolt-stress.md),
equations (32)-(36), proves

    zeta_h(-1)=1/30-q^2/6,
    |zeta_h'(-1)|<604 q^2 <1000 q^2.

This is a known upper bound, not an unknown constant. The argument uses
Euler--Maclaurin through B6, with remainder analytic for Re(s)>-2.
At s=-1 the differentiated remainder is bounded by
8/720 times integral_1^infinity 31968 q^2 x^-3 dx. Endpoint terms
are bounded separately. It is an ordinary summation argument; the
supporting arithmetic checker is not its proof. No split logarithmic
regularization or new finite counterterm is used here.

The [retained action allocation](stationary-curved-quantum-tt.md), section 5,
has, per physical orbit area,

    J=zeta_h'(-1)+(1+2log(mu0 B))zeta_h(-1),
    A0=J/(2pi B^2),   Gamma_a0=-A0/2.

The canonical seed orbit metric is physical orbit metric divided by
A^2=B^2/C. For the affine shape variable b the complete volume action is

    Gamma_a0[b]=-(A^2 A0/2) integral_gamma exp(2b).

Its Hessian at b=0 is multiplication by -2 A^2 A0. Consequently its
normalized fixed-area row, with the SAME -2G/B^2 factor as the complete
cold action defect, is

    H_a0,norm u=c0 u,
    c0=4G A0/C=2G J/(pi C B^2).                       (1)

The factor A^2 must not be omitted and no extra angular sphere volume
may be inserted. A0 already sums the angular multiplicities. The
fourth-order restoration converted with this same A^2 gives
-C F_eff/(960pi), consistent with the previous diagonal calculation.
The canonical massive action depends on b, not on the separate total
four-dimensional conformal variable phi. The latter's anomaly remains
in the Weyl action. Equation (1) is just this allocated contact, not the
complete constant row or the full Hessian.

Let L=log(mu0 B)=log(B^2/b0^2)/2 and use the current interval for B^2,C.
In G=1 units an explicit positive upper allowance is

    U0=2[604 q^2+(1+2|L|)(q^2/6-1/30)]/(pi C B^2).   (2)

Evaluate (2), also retaining the older 1000 version as a comparator.
The exact sign and value of c0 remain unevaluated. The bound is retained
as an error allowance; no physical term is dropped or set to zero.

On the current reference direction in the interior band J=[.49,.88],
the established |u|<1.2e-8 and |D_rho u|<7.7e-6 give

    |H_a0,norm u|<=1.2e-8 U0,
    |D_rho H_a0,norm u|<=7.7e-6 U0.                  (3)

Here c0 is spatially constant in the fixed canonical density. Test
both against the already selected local 1e-28 integration budget.
This does not claim that all remaining components fit that budget.
On the cold region, the already established |u|<=1e-6 a_s gives

    |H_a0,norm u|<=1e-6 U0 a_s,
    contact fraction of massive allowance <=1e-6 U0/(1.72e-15). (4)

This contribution is ALREADY INCLUDED in the complete cold bound.
Do not add it a second time or subtract its upper bound from that bound.
If (3) is below budget and (4) is tiny, resolving zeta_h'(-1) more
accurately is not the next feasibility calculation. This tests the
materiality of that proposed calculation, not the full solution.
It makes no global outer estimate: vacuum cancellation there must be
performed with the physical reference before taking absolute values.

## 2. Retain the physical bolt tolerance, including its h dependence

The preceding [bolt-domain calculation](current-bolt-adjoint-domain.md)
has h in its certified current interval and

    L_b e=ell_b [e_s(0)+e_xi(0)/2], ell_b=(1+h)/h.

Thus the cold difference target eta=1.667e-9 is exactly equivalent to

    |e_s(0)+e_xi(0)/2|<eta h/(1+h).                  (5)

Substitute the retained root enclosure and report the resulting interval
for this threshold, not a fabricated interval for the correction. An
implementation must keep exact h-dependent cancellations before numerical
subtraction. Equation (5) does not require blindly computing two unrelated
order-one fields to 43 decimal places; the actual fields are not known to
be order one. Rescaling the unknown alone cannot improve the certificate.
The combined physical response and its error remain uncomputed.

For e_app in the actual forward domain, d=H e_app+r1 in V*, the retained
unit form inverse gives e-e_app=-H^-1 d. With the exact bolt adjoint psi,

    L_b(e-e_app)=-<psi,d>.

For any psi_app in V, d_b=L_b-H psi_app, symmetry gives the sharper
two-residual identity and bound

    L_b e=L_b e_app-<psi_app,d>-<H^-1 d_b,d>,
    |L_b e-[L_b e_app-<psi_app,d>]|<=||d_b||_V*||d||_V*. (6)

This can exploit a small FORWARD residual instead of multiplying the
adjoint error by the full original defect. It is a specified method, not
an evaluated error. A conforming weak approximation need not reproduce
the two delta sources pointwise, but any missing sources remain in d_b
and MUST be bounded in V*. In an L2 strong-residual implementation, the
complete singular parametrix is needed before an L2 residual is asserted.
Smooth off-bolt samples cannot measure a missed delta. This qualifies the
previous sufficient implementation contract; it does not remove either
physical point source or change the forward/adjoint boundary realization.

If half the eta budget is allocated to the corrected observable in (6),
the remaining residual product must be <8.335e-10. For certified forward
residual D>0 the adjoint residual must be <8.335e-10/D; conversely for
certified adjoint residual E>0, D<8.335e-10/E. Neither D nor E is known.
If either is zero its residual contribution vanishes. The exact
complementary equations and the full X admission still apply.

## 3. Close the audit; perform the actual solve next

The completed cold value/oscillation, trace-only response and current
interior modulus have the scopes recorded in the closed planning audit.
None proves the shape/bolt correction, the full differentiated quantum
response, physical return constants, full-ball embeddings or nonlinear
Y0/Z0/Z2. The latest actual-reference H5 and H6 close a local principal
integration gap, not regularity of every nonlinear direction. The old
huge upper bounds fail their certificates and are not nonexistence bounds.

Boundary/domain requirements come first: use smooth forward bolt data,
the actual paired matter domain, retained mass plus weighted/log outer
phase and all action contacts. Do not impose zero mass, use the old
zero-c projection, or treat the C2,1 reflected seam as globally C5.
Stationary boundary work cancels in the selected phase. Retarded
gravitational constraints and renormalized energy balance require their
own physical realization; they are not supplied by the matter projector.
Product and compact-interior parametrices do not give the quantitative
finite-geometry inverse. Semiclassical validity and all-sector stability
remain open even if the stationary solve succeeds.

Research decision: pursue the direct signed physical finite-domain
H e=-r1 calculation with full quantum moments and a certified weighted
remainder. Use (6) if it reduces the cost of its bolt certificate. The
single next calculation must return the actual coupled correction,
its <1.667e-9 cold error, exterior return and X/admission bounds together.
Do not commission another separated-kernel sample or more accurate a0
value in place of it. No construction switch follows from the present
upper bounds. An actual response incompatibility would justify a dated
revision; another excessive upper bound would only reject that estimator.

The accompanying execution tests (1)-(5), the residual identity (6),
and the available-contact budget. It does NOT execute H e=-r1, compute
either residual norm, or certify the full construction. The next complete
solve remains proposed; the audit must state this without ambiguity.

Remaining to-do list: execute the complete physical coupled correction
with quantified weighted error and admission; nonlinear closure and exact
compatibility; physical matching, causal realization, semiclassical
validity, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
