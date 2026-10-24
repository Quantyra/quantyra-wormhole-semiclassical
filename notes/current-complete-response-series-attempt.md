# Complete response attempt: convergent resolvents do not give the strong trace estimate

2026-10-24, E003/S169. Ordinary argument frozen before arithmetic acceptance.
Full stable semiclassical wormhole objective ACTIVE/INCOMPLETE; audit CLOSED.

Outcome: UNRESOLVED ESTIMATE, with an executed failure of the specified
global Sobolev Born majorant below. Neither eta_seed nor the actual C3
is certified. The complete comparison allowance can, however, be reduced
by a factor of almost 36 in the SAME physical coordinates and norms.
This changes the remaining numerical target substantially. It does not
establish that the unknown actual remainder meets it.

## 1. First unsupported step and the complete expression being attempted

The unsupported step remains a UNIFORM weighted graph-to-H3 bound on
the combined differentiated common finite part. Weak mass-action
existence, seed conversion and a bounded massive resolvent do not give
this bound throughout the nonlinear mass/log family.

Use exactly the full third heat variation (1) of
current-strong-response-decision.md: the tadpole, all three pair contacts,
both triangle orientations, every angular level d_n=2(q+2n), the entire
antiperiodic thermal trace, and the original local restorations. Keep
the composed square of A_b=Dsharp+M exp(b), the paired first-order domain
and its complementary adjoint, and the regular finite-geometry bolt.
No separate finite part is assigned to any individual word.

At a common positive cutoff this expression is well defined. The weak
integrated mass variations and their seed geometric conversion have the
scope proved in current-third-metric-conversion.md. The finite-domain
minus cylinder return and the nonzero thermal windings retain their
existing bounds. They are not recomputed, suppressed or substituted for
the remaining zero-winding expression. The relative split uses Q_c=Q
through .04 with the retained taper to beta outside .05; b_c is the
same physical b through .04, with its retained remote taper. It leaves
the mass cubic at the paired seam intact.

The direct attempt in sections 2-4 treats the ENTIRE variable-mass
resolvent series, not a selected contact or a finite set of modes.
The attempted estimator fails at the strong summation step. Therefore
there is no justification here for exchanging that series, all angular
levels, the strong output derivatives and common regulator removal.

## 2. Exact cylinder reduction on the actual paired domain

Set j=(Q_c/beta)^(1/2), ds/dz=beta/Q_c=j^-2 and s(0)=0. Keep the
physical thermal coordinate tau and its period. The cylinder metric is

 g_c=j^2(ds^2/beta+beta d tau^2).

On paired spinors the unitary density change is U psi=j psi after
coordinate pullback. Directly conjugating

 D_Q=gamma1(sqrt(Q) partial_z+Q'/(4sqrt(Q)))
                                      +gamma2 Q^-1/2 partial_tau

gives, with D_beta=gamma1 sqrt(beta) partial_s
+gamma2 beta^-1/2 partial_tau,

 U(D_Q+M exp(b_c))U^-1
       =j^-1/2 [D_beta+M F] j^-1/2, F=j exp(b_c).       (1)

The coefficients j and j^-1 are bounded Lipschitz, equal on paired
copies, j(0)=1 and j'(0)=0. The coordinate map is bi-Lipschitz. Their
product rules hold in the first-order graph norm; their traces commute
with the original constant paired projector. The same calculation on
the complementary projector preserves the adjoint domain. Close the
identity from the original piecewise smooth domain core. This proves
(1) for the retained first-order realizations, with squares still defined
by composition. It does not invoke a smooth doubled metric. The actual
bolt is carried by the retained relative return, not flattened through
its zero lapse. The thermal period and spin structure are unchanged.

The retained collar/taper guards give |Q_c/beta-1|<.001 and
||b_c||infinity<=.001. Thus f=F-1 satisfies

 ||f||infinity <=sqrt(1.001) exp(.001)-1
          <.001500876<.0016=:r0.                     (2)

The smaller negative endpoint is covered by the same bound. For every
M and the FULL thermal operator, skew-adjointness of D_beta gives
||(D_beta+M)^-1||<=1/M. Consequently the Neumann series for
(D_beta+M F)^-1 converges in operator norm with ratio at most r0.
This is an operator statement, NOT a trace-ideal or renormalization
statement about its full angular direct sum.

Nor does (1) identify the geometric heat determinant with a flat mass
determinant without a contact. With omega=log j and
A_r=j^-r/2(D_beta+M F)j^-r/2, the positive-cutoff conversion is the
endpoint heat term sum_n d_n int_0^1 Tr_bar(omega exp(-epsilon A_r* A_r))dr,
with the original sign convention. It must be combined with the previous
mass/geometric conversion and its paired boundary completion. Omega has
a nonzero second seam jet. Boundary terms cannot be dropped by importing
the special zero-second-jet smearing cancellation for physical b.
We do not claim a numerical strong bound for this endpoint term.

## 3. All-order, all-angular loop estimate attempted before strong summation

Here is the ultraviolet estimate actually derived in this attempt for
the direct words in the unfolded paired resolvent. Reflected words and
their boundary completions remain in the exact inventory; (5) is not
an estimate of those additional terms. Failure of the direct-word
positive majorant already stops the proposed absolute-value proof of
the full expression. In proper flat momenta let
R_M(l)=(i gamma.l+M)^-1. For a cyclic k-vertex
word, k>=3, with cumulative external shifts q_i and P=sum |p_i|, define

 L_k,M=M^k int tr product_i R_M(l+q_i) d^2l/(2pi)^2.

The following bounds hold also for its symmetrization. Each |q_i|<=P,
||R_M(l)||=(M^2+|l|^2)^-1/2, and the spin trace costs at most 2.
The product rule and d^m R/dt^m=(-1)^m m! R(i gamma.q R)^m give

 |d^m L_k,M(t p)/dt^m|
 <=m! binom(k+m-1,m) P^m M^(2-m)/[pi(k+m-2)].        (3)

To verify integrability rather than assume it, apply Holder to the
k factors with exponents (k+m)/(1+alpha_i), sum alpha_i=m. Translation
does not change each integral. The common integral is
int (M^2+|l|^2)^(-(k+m)/2)d^2l/(2pi)^2
=M^(2-k-m)/[2pi(k+m-2)]. There are binom(k+m-1,m)
compositions. All arguments are real; no spatial analyticity is used.

Conjugation by the grading and l->-l show L_k,M(-p)=L_k,M(p).
Subtract its momentum Taylor polynomial through degree four. The
sixth-order integral remainder, using (3), is bounded by

 |L_k,M-Taylor_4 L_k,M|
       <=binom(k+5,6) P^6/[pi(k+4) M^4].              (4)

This subtraction is an estimator decomposition. Its zero/two/four
polynomials still require their common-regulator restoration; they
are NOT discarded physical contacts. For the full massive tower use
N(M)<=4 C M^2, M1>m=18. Split at M=P. On M>=P, (4) and
sum_(M>=P) d/M^4<=8C/P^2 apply. On M<P use (3) for the raw word and
Taylor coefficients. Stieltjes integration gives
sum_(M<P) d/M^2<=4C+8C log(P/m). Thus, for EVERY k>=3,

 sum_n d_n |L_k,Mn-Taylor_4 L_k,Mn|
 <=C P^4/pi {8/(k-2)+2(k+1)
       +binom(k+3,4)/(k+2) [4+8 log_+(P/m)]
       +8 binom(k+5,6)/(k+4)}.                       (5)

For P<m the high-mass estimate alone gives (5), using P^2/m^2<=1;
P=0 has zero remainder. No angular level is sampled or omitted.
For example the k=4 brace is 364/3+(140/3)log_+(P/m).
The coefficients grow polynomially in k, not exponentially.

In the formal log series the k-word coefficient is (-1)^k/k. Third
mass-field differentiation includes all vertex choices, giving the
three chain-rule weights

 (k-1)(k-2) F^3 f^(k-3),
 3(k-1) F^2 f^(k-2), and F f^(k-1),                 (6)

where products denote the distinct slots of the multilinear kernel.
They represent the triangle, three pair contacts and single contact;
k=1,2 must additionally retain their original common finite parts.
Equation (6) is not permission to multiply divergent mode traces.
Equations (3)-(5) justify the subtracted word on bounded momentum sets.
They do not yet identify the all-order strong physical response.

To try that identification directly, the global Fourier-Young estimator
allocates P^(4+r) log(2+P), r=1,...,4, among ALL vertex slots. The highest
allocation requires H8,log of one slot and Fourier L1 of the others.
At k>=4 there are background slots as well as the two corrected fields.
Convergence in (2), which is an L-infinity operator perturbation bound,
does not control these Fourier norms. The failure is explicit below.

## 4. Actual seam calculation: the proposed global majorant is infinite

Write the actual even seam expansion

 Q_c/beta=1+a2 z^2+a3 |z|^3+O(z^4),
 a2=u/beta, a3=(-4beta-2u-tau)/beta.

The physical admitted b has cubic coefficient m_b, |m_b|<=.001/6;
its quartic logs and fast graph remainder are retained. Since
z=s+O(s^3), the exact flattened f has

 f(s)=(a2/2)s^2+c3 |s|^3+O(s^4 log^d|s|),
 c3=a3/2+m_b,
             -2.000167<c3<-1.999833.                (7)

Here d is the retained finite log degree. Choose an analysis cutoff equal
to one near the seam and supported inside the common collar; it does
not change the physical domain. Its fourth distributional derivative is

 D_s^4(cutoff*f)=12 c3 delta_0+g, g in L1.

The remainder assertion follows from the actual finite mass/log phase
and the fast graph embedding, not from a smooth extension of the cubic.
The unitary Fourier convention and the Riemann-Lebesgue lemma yield

 p^4 Fourier(cutoff*f)(p) ->12 c3/sqrt(2pi) !=0.      (8)

Hence this coefficient is not in H8,log; it is not even in H^(7/2).
Its Fourier moment int |p|^a |f_hat|dp diverges for a>=3. The specified
global positive majorant of (5)-(6) is therefore infinite, already at
the first background-inserted word, even at b=0. Summing further words,
thermal tails or positive return allowances cannot make that estimator
finite. Smoothing the seam would change this finding by changing the
domain and is excluded. This is an executed CONTINUUM failure of the
candidate estimate, not a finite-mode diagnostic.

This does NOT prove that the combined actual response is unbounded.
Both fast physical inputs vanish as z^sigma, sigma=4.005. Taking the
global H8 norm of a background before combining its kernels with those
inputs loses precisely that information. Even the local product
|s|^3 u v starts at order 3+2sigma>11; treating its differentiated
background and its two inputs separately creates spurious seam terms.
Nonlocal terms need their own simultaneous-derivative argument; this
observation alone is not its bound. The missing argument is an estimate
on the COMBINED weighted kernel, its conversions and local restorations,
with high background frequencies retained until the vanishing inputs
and actual paired traces have been used. No such finite C_rem is obtained.

## 5. Quantitatively material correction to the complete known allowance

The previous mixed-frequency note used proper flat p. In the ORIGINAL
z coordinate, for g_beta=dz^2/beta+beta d tau^2, the exact stationary
comparison is T_beta(k_i)=T(sqrt(beta) k_i). Area density is one.
Proper time contributes sqrt(beta) and proper spatial measure
1/sqrt(beta); they cancel in the action per fixed physical tau. No
change of time period, physical units, graph data, inverse or norm occurs.

Consequently the full fourth-and-higher kernel bound is

 |T_beta|<=7 beta^2 C P^4 log(2+P), beta<1.          (9)

This includes the fixed fourth contact, which scales by beta^2 exactly.
For the separated kernel,
k_A,beta(s)=beta^-1/2 k_A(s/sqrt(beta)), so its r-th derivative guard
16 r!4^r C |s|^(-5-r) also acquires beta^2. In the mixed kernel take
B_beta(k)=beta^-1 B(sqrt(beta)k); then (k B_beta)'=-k^3 A_beta'.
The primitive identity, its remote derivatives and ALL mixed tails
therefore acquire the SAME factor. These identities verify both parts
of the existing infinite-annulus proof; no local-to-global norm is guessed.

The proof of current-weighted-third-comparison.md explicitly replaced
beta^2 by 1. Retain it in both its near and far constants, keeping K,
W, chi, weights, cutoffs and the actual physical fields unchanged. The
complete coefficient repair of current-strong-response-decision.md
inherits the same factor. At the correlated physical parameters,

 C_flat<1.295e10 beta^2<359723000,
 C_delta<407000 beta^2<11306.                       (10)

This closes these TWO comparison allowances with stronger constants,
not the actual variable-geometry remainder. They formerly cost .12950407
at radius 1e-11 and now cost <.00359734306. Retain the previous 4.001e-9
reserve, without recharging the already reassigned 2e-11 principal debit.
With C3<=359734306+C_rem, a sufficient version of the requested test is

           eta_seed+1e-11 C_rem<.156319052939.       (11)

Both quantities on the left remain UNKNOWN. In the unattained best
case eta_seed=0, C_rem must be <1.56319052939e10. Equal allocation asks
eta_seed<.0781595264695 and C_rem<7.81595264695e9. These thresholds are
about 5.14 times the previous remaining thresholds. The discarded beta
factor was a real estimator loss; restoring it is not renormalization.
There is still no evidence placing the actual remainder below (11).

## 6. Original admission, theorem scope and research decision

Use alpha for ORIGINAL X-to-fast-data norm and R_X for the original
preconditioned H3-to-X return (avoiding confusion with physical beta).
Exactly as before one needs

 E_center+1e-8 alpha<=1e-11,
 R_X alpha eta_seed+Z0_other<=.25,
 R_X alpha^2 C3+Z2_other<=5e7.                     (12)

Even E_center=0 requires alpha<=.001. Using the entire Z2 allowance for
the comparison (10) alone requires R_X alpha^2<=5e7/359734306,
approximately .13899147. If alpha=.001, R_X must be below about 138991
in that optimistic comparison-only calculation. These are thresholds,
not established norm comparisons. If the actual C3 used the full
best-case user budget about 1.599e10, the corresponding product ceiling
would be only about .00313. The improved comparison-only number must
not be advertised as the allowable return for the complete response.

The original preservation and contraction requirements remain
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1. The allocations
(Y0,Z0,Z2)<=(5e-9,.25,5e7) remain unproved and leave no spare preservation
budget at their ceilings. The fast A contribution <5.681e-5 leaves
.00094319 for the admitted background/other components, also unproved.
The one-way inverse/embedding constants do not establish (12).

No external smooth-domain theorem is used to fill this gap. Holder,
the explicit two-dimensional integral and Taylor's integral remainder
have their hypotheses verified in (3)-(5). The Neumann theorem applies
to the bounded operator perturbation (2), not to its angular trace.
The Riemann-Lebesgue hypothesis in (8) is g in L1, checked from the
finite mass/log phase. The first-order domain transport is proved by
the actual product rule and paired traces in section 2.
[Grubb Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf) requires
actual symbol seminorms and |s|<tau. Its hypotheses do not give H3
from this C2,1 seam, and no numerical constant is extracted from it.
The existing weak conversion theorem is used only at its stated scope;
neither its seed result nor finite propagation gives the missing uniform
strong bound. No determinant Schatten assumption is inferred from (2).

The prior COMPLETE selected projection screen E_rest<2.123e-42<1.287e-40
stands; it is not newly evaluated and its sign is unmeasured. No actual
controlled correction follows. Exact compatibility and complementary
equations, physical matching/causal realization, semiclassical validity
and noise, and all-sector stability remain unresolved. The old audit is
not reopened. There is no demonstrated obstruction to the construction.

DECISION: change the estimation method, retain the construction. The
global Sobolev Born majorant fails on the actual seam even though its
resolvent series converges. Use the paired weighted kernel, keeping each
background insertion together with the two vanishing graph inputs until
the common local subtraction and metric conversion have been combined.
The recovered comparison budget makes that method numerically worth
testing; it does not establish its quantitative attainability.

SINGLE NEXT CALCULATION, PROPOSED AND NOT EXECUTED: bound the complete
actual-minus-comparison zero-winding response and its seed partner on
the original infinite paired annuli, retaining all background insertions,
local contacts, conversion and returns, to test (11). It must supply
simultaneous input/output derivative bounds with an all-order summable
remainder, not another isolated contact. Report its bound in the original
norms (12); a finite graph bound alone will not close original admission.

Remaining to-do list: the combined strong estimate and original norm
comparisons; actual correction and nonlinear admission/contraction;
exact compatibility; physical matching/causal realization; validity/noise;
all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
