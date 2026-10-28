# Complete scalar-log solver attempt: norm transfer and an unresolved complement

2026-10-28 HST, E003/S169. The preceding physical-log reference rejection
was progress. The current full-range solve remains UNRESOLVED: no implemented
inverse or complete defect upper bound was obtained. The calculations below
test a specified possible finite correction and fix an internal norm transfer;
neither substitutes for the missing complementary operator estimate.

## Actual operator, domain and the first missing implication

Keep A_sc=(ell_sc+B)/N on E_t->E_t*, B=b* a^-1 b, N=|W_action| beta.
The actual paired Dirac/adjoint/square domains, common regulator, UV log,
full angular/thermal sums, contacts, physical mass/log completion and regular
bolt are those in current-global-scalar-log-comparison.md. Both original
equations remain in B. The physical scalar Laplacian reference f(L) and its
domain are fixed in current-physical-log-preconditioner-test.md; the original
admission norms and P32 do not change.

Attempt a full-range reference corrected in the admissible quartic direction:

    F_gamma=f(L)+gamma |psi><psi|_H, psi=z^4, gamma>=0,
    P_gamma=F_gamma^-1 : E_t* -> E_t.                   (1)

The finite term is bounded on E_t because E_t embeds in H. It preserves
coercivity, the mass/log form domain, and regular-bolt realization. The exact
rank-one inverse formula makes (1) explicit in terms of f(L)^-1 and its
quartic application, but these applications are NOT numerically implemented.
Unlike a finite output alone, (1) retains an infinite-dimensional range.

The first missing quantitative implication is an upper bound on the complete
defect on the complementary space, together with its coupling to the finite
part, in the ORIGINAL dual norm. No bound for A_sc-A_Q supplies a bound for
A_sc-F_gamma. The qualitative graph equivalence for f(L) gives no numerical
comparison for this difference. The actual-to-normal parametrix remainder
C h^2(1+|log h|)^N still has an uncalibrated C; its cutoff commutators and
interior/boundary returns cannot be assigned zero or a conditional ceiling.
This is an unknown estimate, not a newly undefined seed form.

## Internal energy comparison needed before a numerical defect is interpreted

The accepted complete estimates give, without a physical normalization change,

    .997216 R_t <= A_sc <= 1.341 R_t, R_t=ell_t/N.      (2)

Indeed |ell_sc-ell_t|/N<.002784 R_t and B/N<=q_t R_t,
q_t=.338/(1-eps), eps<1.227e-74. Set ||f||sc,*^2=<f,A_sc^-1 f>.
Inversion of (2), followed by the two norm comparisons, proves for ANY
bounded defect D on the common dual

    ||D||Et* <= sqrt(1.341/.997216) ||D||sc,*
              <1.159631 ||D||sc,*.                    (3)

Thus a defect <=.49 measured in the compensated operator's dual energy
norm is insufficient by itself to claim <=.49 in the original dual.
The sufficient internal target is

    ||D||sc,* <= .49/sqrt(1.341/.997216)
                 =.42254825263115... .                (4)

This is an established comparison between two specified INTERNAL norms.
It is not a comparison with the original graph, radius or physical output.
The .49 original-dual target and its <.490978 transfer to A_Q are unchanged.

## Executed necessary test of the finite correction

Use the two actual fields already admitted in the full quartic test:

    xi_1=4[z^3 Q+z^4 sum_(k>=2) k Q_k z^(k-1)/(k-1)], J_1=z^4,
    xi_2=psi=z^4, J_2=0.

Their exact York equations, outer conditions, cutoff approximation and regular
bolt are retained. xi_1 carries a nonzero admissible mass cubic; neither is
a replacement matter mode. Let K be the accepted complete local/LLL/vacuum
energy matrix, and B_low the accepted scalar variational LOWER matrix.
The positive complete geodesic remainder obeys a matrix upper bound
G<=1e-35 I: Cauchy-Schwarz gives G<=2 diag(e_geo,1,e_geo,2), and both
accepted diagonal upper bounds are below 5e-36. This includes all returns.

With e=.002784, d=.0022 and eps=1.227e-74, the complete form compression has

    A_lower=(1-e)(1-eps)/(1+d) K+B_low,
    A_upper=(1+e+.338000000000001)(1+eps)/(1-d) (K+1e-35 I).

The upper B allowance uses the complete mixed-form bound, not B_low as
an approximation to the full inverse. The lower bound uses both scalar
tests in the variational maximum. Hence A_lower<=A_sc|U<=A_upper.

Set M1_ij=<xi_i,L xi_j>H and M3_ij=<xi_i,L^3 xi_j>H. Both are exact
polynomial Dirichlet integrals on the correlated finite seed. The whole-
spectrum bounds from the preceding physical-log test give

    F_lower=[1.5(1-1.179e-8)-60e-33/pi] M1,
    F_upper=1.5(1+1e-30) M1+(3e-33/(20pi))(400M1+14M3).

These are matrix inequalities, not selected spectral-mode replacements.
Write v_i=<psi,xi_i>H, so the finite update is gamma vv^T. If the full
original-dual defect in (1) were <=.49, spectral necessity would give

    1.49(F_upper+gamma vv^T)-A_lower >=0,
    A_upper-.51(F_lower+gamma vv^T) >=0.               (5)

The original-dual spectral implication follows from the same positive
similarity as the physical-log rejection. Conditions (5) are necessary,
not sufficient; even all finite compressions require an additional uniform
argument before they can certify the whole space. In two dimensions the
determinants in (5) are affine in gamma. The outward calculation gives
the necessary window, rounded OUTWARD,

    66.7925 < gamma < 90.2951.                         (6)

This NONEMPTY window does not reject the candidate family. Nor does it
certify any gamma: the actual scalar form has not been substituted by
either envelope, and the complementary space remains uncontrolled.

The current separated upper-envelope estimator cannot certify the desired
upper spectral inequality for ANY gamma. Take w=xi_1-c psi,
c=<psi,xi_1>H/<psi,psi>H=-1.01190476185265... . Then <psi,w>H=0, so the
finite update vanishes. The recorded calculation evaluates

    A_upper[w,w]/F_upper[w,w] >1.69,

already above the necessary upper endpoint 1.49. This is an EXCESSIVE
KNOWN UPPER ENVELOPE, not a lower bound on the actual Rayleigh quotient.
The corresponding A_lower/F_upper quotient is below 1.31, so this trial
does not establish an actual violation. The receipt also records the
scalar B allowance needed to make this particular row fit, its current
upper allowance, and their ratio. Refining that one row would still not
control the complementary map, so no isolated scalar optimization is selected.

Verification provenance: the first arithmetic run incorrectly asserted
the stronger threshold 1.72 and failed. Its immutable receipt is retained.
Inspection gives 1.69183290450261213..., so the final outward assertion
uses 1.69. The separated envelope still exceeds 1.49; no failed run is
counted as acceptance. Use the separately named accepted receipt.

## Decision and consequence for the complete goal

Outcome: UNRESOLVED complete solver. Executed evidence establishes the
internal conversion (3), a nonempty necessary finite-correction window,
and failure of the specified separated upper-envelope certification.
It supplies no complete defect upper bound, application error, signed G_t,
source allowance or physical margin. The 4.9887427328e-32 source screen
remains conditional, and all actual S/K errors stay unknown.

Retain the compensated scalar-log method and construction; do not infer a
working inverse from a small finite matrix or from an exact inverse formula.
The decisive next calculation remains the complete complementary remainder
of an IMPLEMENTED full-range action, with numerical bounds for the joint
logarithm, cutoff/interface terms and paired boundary/bolt returns. A
corrected reference may be used only if its entire source action and these
errors are enclosed. This calculation is proposed, not executed here.

Remaining to-do list: complete implemented inverse/source/application
certificate; signed G_t and all original graph, physical/radius, mass/phase
and coupled budgets; nonlinear admission/contraction, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
