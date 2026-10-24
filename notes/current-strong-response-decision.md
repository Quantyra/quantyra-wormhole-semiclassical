# Direct strong-response attempt: a cutoff derivative gap and the remaining test

2026-10-24, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress. Full stable semiclassical wormhole objective
ACTIVE/INCOMPLETE. The feasibility audit remains CLOSED.

Outcome: the COMPLETE strong estimate remains UNRESOLVED. This attempt
does not certify eta_seed or C3. It identifies a failed step in a precise
candidate proof of that estimate, quantifies a replacement comparison,
and states the missing combined estimate and original-norm requirements.
It is not a new contact milestone or evidence against existence.

## 1. The first unsupported step, and the direct attempt

The existing weak action, closed first-order paired realization and
relative/thermal estimates do not establish the local zero-winding
H3 response uniformly on the nonlinear graph/mass-log family. The
principal-coefficient estimate in current-uniform-principal-coefficient.md
is a multiplier estimate only. The first unsupported step is promoting
the COMMON, combined, differentiated heat finite part to a bounded
bilinear graph-to-H3 map on that family. In particular differentiating
a matched *linear* principal subtraction does not supply a matched
third-response subtraction in every frequency orientation.

Here is the direct expression attempted. On the original fixed paired
first-order domain use A_b=Dsharp+M exp(b) and its composed positive
square P_b. Write K_t=exp(-t P_b) and

 P_I=2^|I| M^2 exp(2b) product_(i in I) u_i
                          + M[exp(b) product_(i in I)u_i,Dsharp].

At positive COMMON heat cutoff epsilon the third heat variation is

 -t Tr(P_uvw K_t)
 +t^2 sum_(uv,w;uw,v;vw,u) int_0^1
                       Tr(P_uv K_((1-s)t) P_w K_(s t)) ds
 -t^3 int_(alpha,beta>=0,alpha+beta<=1)
   {Tr(P_u K_(alpha t) P_v K_(beta t) P_w K_((1-alpha-beta)t))
    +Tr(P_u K_(alpha t) P_w K_(beta t) P_v K_((1-alpha-beta)t))}
                                                        d alpha d beta. (1)

Cyclicity combines the six ordered insertions into the two triangle
orientations in (1). The simplex area is 1/2, so the commuting limit
is -t P_uvw+3t^2 P_uv P_w-t^3 P_u P_v P_w, with the three distinct
pair terms understood. Integrate (1) against dt/(2t), sum EVERY angular
level with its prescribed multiplicity, and retain the FULL thermal
trace. Differentiate the same local subtractions, finite restoration
and geometric/mass conversion. No separately renormalized word occurs.
Its finite-domain kernels still have the original paired seam and regular
bolt. The actual-to-cylinder split retains the full relative return.

Apply both existing normal inverses, the sharp input split, the original
shape row, and the H3 output. Trying to obtain an operator-norm Cauchy
limit of this combined expression after the latest LEFT-coefficient
principal subtraction fails as shown below. This tests the complete
actual expression, rather than an isolated triangle or contact. The
previous small coefficient and flat-comparison bounds remain valid.

## 2. Actual-domain failure of this particular differentiated cutoff limit

Set x=log z, sigma=4.005, P=K D^5 log(2+|D|),
K=GC beta^2/(15pi B^2)>0, a0=(Q/beta)^2, rho=a0-1.
The selected left principal operator is -a0 exp(-2b)P. After the
already allocated seed curved and flat linear pieces its coefficient is

 r_b=a0(exp(-2b)-1)+2b,       D_b r_0[v]=-2 rho v.    (2)

Include chi in every output. Differentiating (2) gives a term proportional
to rho v P u. It places the logarithm on the *u* input. It does not
include the equally required high-frequency *v* input contribution.

Use the actual compact interval I=(.0019,.0021), strictly inside the
physical half and inside chi=1. Let u be a fixed smooth admitted field,
equal to 1 near I and supported in (0,.003), and let

 v_N=z^sigma N^-6 psi(log z) exp(iN log z),

with psi smooth, nonzero, supported in I. Set the scalar metric component
to zero, retaining its induced triangular comparison datum. Real/imaginary
parts give real fields. Both directions vanish near the ACTUAL seam and
bolt, satisfy the paired traces, have zero mass/slow coefficients and
finite log degree zero. They are tests of the original operator. No
physical mode truncation or smooth replacement of the domain is made.

The complete local symbol and its soft variation were derived in
current-zero-winding-regularity-test.md from the actual integrated heat
coefficient, with all coincidence vertices. In the present normalization,
on the patch where u=1,

 S D H_0[v_N]u = 2K a0 D^5 log|D| (z^-sigma v_N)
                                      + terms of order at most five. (3)

The complete seed flat comparison has the same expression with a0=1.
Allow its entire third comparison to be cut off with the SAME epsilon;
this is more cancellation than merely subtracting an unregulated flat
response. The remaining leading coefficient is still 2K rho.
All finite polynomial contacts have no frequency logarithm. Separated
seam/bolt returns and nonzero thermal windings cannot change this local
diagonal symbol. They remain in the operator, with their established
locality/gap bounds, rather than being assigned zero.

At fixed epsilon>0, the full heat expression (1) and its common local
restorations have order at most four on a compact smooth interior, hence
at most five after the shape row. To check the needed packet bound, split
each Duhamel simplex according to its longest heat segment. Its length
is at least t/3; its smooth kernel and massive Gaussian dominate the
entire angular/thermal sums for t>=epsilon. Place the finitely many
metric derivatives on the packets. The resulting bound is
O_epsilon(N^8 N^-6) in H3. The constant may diverge as epsilon decreases.
The same holds for the matched flat heat expression. It is NOT a bound
uniform in epsilon. There is no N^8 log N term at fixed positive cutoff.

Any addback obtained ONLY by differentiating the left coefficient in
(2) acts as v_N times the fixed smooth P u (or its regulated version).
Its localized H3 norm is O_epsilon(N^-3), and cannot restore (3)'s
missing logarithm. Lower local row terms cannot do so either.

Let T_epsilon^left denote this attempted fully regulated, flat-matched,
left-principal matched third remainder, and T^left its distributional
finite part. The existing weak theorem defines the latter on these
smooth tests. If it admits the desired strong extension, the packet test
gives, with a fixed localization zeta inside u=chi=1,

 lim_(N->infinity)
 ||zeta (T^left-T_epsilon^left)(u,v_N)||H3 / ||data(v_N)||H3
       =2 ||rho psi||2/||psi||2 >0,   every epsilon>0. (4)

Here ||data(v_N)||H3=K N^2 log N ||psi||2(1+o(1)); the scalar datum
is lower order. Divide also by the finite fixed ||data(u)||H3 to get
a bilinear operator-norm lower bound. This last datum is not presumed
to be one or numerically measured. Multiplication by zeta is bounded,
so (4) rules out norm-Cauchy convergence of the *global* bilinear map.
If the strong extension does not exist, the proposed inference fails
earlier. No finite packet frequency is claimed to violate a physical
budget: its asymptotic onset has not been numerically bounded.

Using the actual correlated Q, the complete interval I satisfies

      1.07e-7 < 2|rho(z)| < 1.49e-7.               (5)

Thus this is a nonzero obstruction to the SPECIFIED norm-limit argument,
not an oversized upper bound. Its small size does not repair norm
convergence. It also does not show that the desired C3 bound is too large.
The completed weak projection and all earlier component bounds stand.

## 3. A quantified replacement in the COMPLETE bilinear comparison

Use the already proved COMPLETE T_ge4 comparison, including the three
pair contacts, both triangle orientations, its fixed fourth restoration,
all angular levels, original sharp input split and every mixed tail.
Let R_flat be its weighted graph-to-H3 map. Its proof gives
C_flat<1.295e10, with an explicit final cutoff multiplier allowance 100.
It retains the highest logarithmic derivative in L2.

For every admitted real b with ||b||A<=.001 define

 delta_b=a0 exp(-2b)-1,
 R_cmp,b=(1+delta_b) R_flat.                         (6)

This is a comparison at the BILINEAR RESPONSE level; it is not falsely
identified with the derivative of a newly chosen physical action. The
exact remainder is DEFINED by subtracting (6) from the actual combined
third response. It retains all lower quantization terms, low contacts,
conversion, boundary completions and non-fast background contributions.
In particular multiplying the output coefficient need not preserve the
three-slot action symmetry; its compensating terms stay in that exact
remainder. The physical action/source is unchanged.

Unlike (2) alone, (6) includes the high-frequency orientation in BOTH
inputs and the complete mixed-frequency region. At b=0 it restores the
2K rho term in (4). It is not a proof that its uncomputed actual remainder
has a uniform strong bound. Its use with the same cutoff must use the
ENTIRE regulated comparison, not just its high-input coefficient.

The existing annular proof applies to the product chi delta_b by replacing
its final factor 100. Write a=.01, B_j=.001*3^j*a^3/6, j=0,...,3.
For a0=sum a_k z^k, a_0=1, put

 R_j=sum_(k>=1)|a_k| k^j a^k, e=exp(2B_0),
 E=[e-1,2e B_1,e(2B_2+4B_1^2),
                         e(2B_3+12B_1 B_2+8B_1^3)],
 F=[e,E_1,E_2,E_3],
 A_j=sum_(i=0)^j binom(j,i)R_i F_(j-i)+E_j,
 M_delta=100A_0+51A_1+12A_2+A_3.                    (7)

The actual vanishing traces and bounded third physical derivative give
the B_j caps; no bounded fourth logarithmic derivative is assumed.
The Leibniz rule gives |D^j delta_b|<=A_j. The same argument for the
ENTIRE near and remote output then proves

 ||delta_b R_flat||_(H3 x H3 -> H3)
              <= (M_delta/100) C_flat =:C_delta.    (8)

This is not multiplication of two differently normalized norms. The
factor 100 is removed and replaced within the SAME explicit proof.
No new estimate of a small cold, scalar, or isolated contact is used.
It bounds the complete comparison correction required by the failed
orientation test, including its mixed-frequency and separated inputs.

Development evaluation gives M_delta=0.003137746541395028... and
C_delta=406338.177110656... . Freeze the safe target

                         C_delta < 407000.          (9)

Its feedback cost at radius 1e-11 is <4.07e-6. That quantifies the
replacement's affordability, not the unknown actual remainder.

## 4. What regulator conclusion is actually needed

The failed operator-norm cutoff convergence is STRONGER than necessary
for a contraction argument. Do not add another impossible gate. At a
fixed background a uniform bound of the COMBINED cutoff responses,
followed by their established distributional common limit, gives a
bounded H3-valued limit by weak compactness and lower semicontinuity.
Uniformity on a line segment of backgrounds then permits testing the
weak fundamental-theorem-of-calculus identity against every H3 dual
vector. This yields the required operator Lipschitz bound for DF.
It does not require operator-norm convergence of D H_epsilon.

Pointwise strong convergence for each graph input would additionally
follow from convergence on a dense core IN H3 and that uniform bound.
Neither distributional convergence alone nor the current weak heat
remainder proves this dense-core H3 convergence at the mass/log seam.
The finite-dimensional mass/slow components must be retained separately;
they cannot be removed by taking a compact fast core.

The direct attempt therefore stops at the following PRECISE missing
argument: uniformly on the stated graph/mass-log neighborhood, bound
the combined actual response (1), with all its common subtractions and
returns, minus (6), after both normal inverses, in H3. One needs tame
simultaneous output/input derivative estimates of this ONE combined
kernel/finite-part expression. Existing gain-two VALUE/form envelopes
do not bound those derivatives. The scalar left coefficient does not
bound the high-metric-frequency terms, as (4) verifies. No numerical
constant for this remainder has been obtained in the present attempt.

The symbol shortcut also fails its hypotheses: the cited rough-symbol
mapping theorem requires its actual normalized symbol seminorms and
|s|<tau. A reflected C2,1 coefficient is not a C^tau symbol with tau>3
on which one can simply assert an H3 estimate. Weighting/localizing may
repair the problem only with the very derivative estimates just named.

## 5. Full feasibility and ORIGINAL nonlinear norm ledger

Keep the user's target eta_seed+1e-11 C3<.1599164. Retain the prior
4.001e-9 reserve, giving .159916395999 before this comparison. With
the exact repartition (6), write C3<=1.295e10+407000+C_rem. The previous
2e-11 direct principal-family debit is NOT charged again: its nonlinear
principal component has been reassigned to (6) and the exact remainder.
The old multiplier bound remains true and historical; it is not a
second allocation. A sufficient revised test is

             eta_seed+1e-11 C_rem < .030412325999.   (10)

Both eta_seed and C_rem remain UNKNOWN. Best case eta_seed=0 requires
C_rem<3.0412325999e9. Equal splitting requires
eta_seed<.0152061629995 and C_rem<1.52061629995e9.
The replacement uses about 0.0134 percent of the previously remaining
budget. It is quantitatively affordable; no evidence yet bounds the
actual combined remainder by the available amount.

The graph test is not the original X-space contraction. Let alpha be
the norm of the ORIGINAL X-to-fast-data map, beta the ORIGINAL
preconditioned H3-to-X output return, and E_center its graph center
offset. Required comparisons include

 E_center+1e-8 alpha<=1e-11,
 beta alpha eta_seed+Z0_other<=.25,
 beta alpha^2 C3+Z2_other<=5e7.                     (11)

Z0_other and Z2_other must include all other rows, finite-dimensional
mass/slow directions and global matching; no values are assigned.
Even with E_center=0 the first requires alpha<=.001. Spending the
ENTIRE Z2 budget just on the guarded comparison (6) would require
beta alpha^2<=5e7/(1.295e10+407000), approximately .003861.
If alpha=.001 and every other allowance were zero, beta would have
to be below approximately 3861. These are conditional thresholds,
not established comparison constants. Lower alpha also changes its
actual return comparison; rescaling alone does not improve the product.

The known fast-to-A embedding 5.681e6 gives A_fast<=5.681e-5 at this
radius and requires the rest <=.00094319 for A<=.001. The H6 inverse
1801 and the fast weighted-C5 embeddings are one-way statements. They
neither prove alpha nor identify beta with 1801 or 5.681e6. Trying the
reverse C5-to-full-graph comparison fails on the retained oscillatory
packets: the graph carries an eighth derivative and a logarithm.
The frozen original audit specifies X and its embeddings as obligations;
it does not supply these quantitative maps. This attempt therefore does
not establish the original norm comparisons requested in (11).

The exact original requirements remain
 Y0+1e-8 Z0+5e-17 Z2<=1e-8,        Z0+1e-8 Z2<1,
with unproved allocations Y0<=5e-9, Z0<=.25, Z2<=5e7. The ceilings
consume all preservation budget. Actual residual/correction, global
admission, contact targets M_A<=3.18e27 and M_V<=2.66e30, original
range compatibility and complementary equations remain unresolved.

## 6. Source hypotheses, verification boundary, and decision

The source checks used in this attempt are deliberately limited:

* [McIntosh--Morris, Theorem 3.1](https://arxiv.org/pdf/1201.1818): a C0
  group, domain-preserving Lipschitz multipliers, bounded first and zero
  second commutators give finite propagation. The existing doubled
  self-adjoint first-order realization verifies these on the paired
  metric space. This controls separated returns; it gives no coincident
  third-response constant or derivative bound.
* [Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf): mapping
  for actual C^tau S^d_(1,0) symbols with |s|<tau and N>m/2+1. The
  required symbol seminorms and seam hypotheses for the present H3 map
  have not been proved. This theorem is NOT used to assert the bound.
* [Moller, introduction and Theorem 1](https://arxiv.org/html/0705.3857v3):
  smooth closed spin geometry and smooth metric variations. It does not
  apply globally to the paired mass/log domain. The local packet argument
  instead uses the retained actual-domain compact-interior heat-symbol
  proof, whose smooth nested patches avoid the seam without replacing it.

The arithmetic verification checks (5), (7)-(11) and the complete prior
comparison inputs. It does not prove an uncomputed kernel estimate.
All new proof scope and guards are frozen before acceptance; development
and failed verification attempts, if any, are retained. No Lean is used.

RESEARCH DECISION: CHANGE THE ESTIMATION METHOD, retain the construction.
Do not seek a bilinear operator-norm cutoff limit by differentiating
only the left principal coefficient. Use the complete mixed-frequency
comparison and a uniform bound for the combined renormalized response,
then the common weak limit and the pathwise identity. The explicit
comparison cost fits the numerical ledger, but full attainability is
still undecided. No actual construction obstruction has been found.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the combined
actual-minus-(6) response uniformly, with its seed partner, in the
original weighted graph, testing (10) and the actual X returns (11).
In concrete terms this is the simultaneous-derivative/finite-part
estimate for (1), including low contacts, paired seam and finite-domain
returns; it is not another principal coefficient calculation. The
present attempt establishes the failed cutoff inference and the cost
of its complete comparison replacement. It does NOT execute that
remaining estimate or a controlled solution.

The certified selected projection E_rest<2.123e-42<1.287e-40 remains
unchanged; its sign is unmeasured and it is a sufficient screen only.
Exact compatibility, physical matching/causal realization, semiclassical
validity/noise and all-sector stability remain unproved. No claim of
existence, nonexistence, stability or completed full objective is made.

Remaining to-do list: combined strong quantum remainder and original
norm comparisons; actual correction and nonlinear admission/contraction;
exact compatibility; physical matching/causal realization; validity/noise;
all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
