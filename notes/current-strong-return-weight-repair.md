# Strong return attempt: logarithmic weight loss and the original two-row normal inverse

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: progress in the comparison method, not in quantitative
coupled admission. Full stable-wormhole objective ACTIVE/INCOMPLETE; audit CLOSED.

OUTCOME: the actual quantitative coupled return is UNRESOLVED. Attempting
to assign its constants reveals an unsupported endpoint in the preceding
short-travel proof. We repair that endpoint and numerically bound the
complete TWO-ROW NORMAL inverse needed by the repair. This is neither a
finite-geometry inverse nor an evaluated quantum response. Its role in the
complete calculation is specified below; no physical operator is replaced.

## 1. Correction before using the preceding strong-domain argument

Section 4 of current-short-travel-strong-domain.md passes from weight 3/2
to EXACTLY 7/2 after one application of J_s=H^-1 R_s. Its preceding
parametrix estimates allow h^2(1+|log h|)^N. With N>0 that allowance does
not prove a bounded gain of exactly two powers in the unlogged weighted
Sobolev spaces. Finite differentiability of coefficients is not enough
to remove this loss. The old evidence is retained, not overwritten.

Here is a complete test of this inference, not a substitute quantum
operator. Set x=-log z, delta=3/2 and let psi be a nonzero smooth function
supported in (0,1). For L>log 10 let u_L=z^delta psi(x-L). Its
z^delta H_b^k norm is independent of L. Multiplication by

    A_N(z)=z^2(1-log z)^N

satisfies the displayed coefficient envelope and every finite conormal
derivative bound of the same type. Nevertheless

    ||z^(-delta-2) A_N u_L||_L2 >=(1+L)^N ||psi||_L2.

The unlogged endpoint estimator is therefore unbounded for every N>0.
This does NOT prove the actual R_s has such a noncancelling coefficient,
or that its endpoint map is unbounded. It proves that the stated allowance
cannot justify that step without an additional cancellation or logarithmic
weight. No physical source has been selected to equal A_N u_L.

Use instead the auxiliary first weight delta_1=13/4. The first gain is
7/4, leaving 1/4 for outer logarithms. For x>=log 10 define

    m_N=sup exp(-x/4)(1+x)^N,
    m_0=10^(-1/4),
    m_N=(4N)^N exp(1/4-N), N>=1.                     (1)

Differentiating the positive scalar function proves (1): its only
interior maximum is x=4N-1, which exceeds log 10 for every integer N>=1.
For integer k>=0 the full H_b^k multiplication bound is

    B_N,k=sum_(l=0)^min(k,N) binom(k,l)(5/4)^(k-l)
                                      (N)_l m_(N-l).           (2)

Indeed bound each derivative of exp(-x/4)(1+x)^N by the product rule,
then sum the H^k product estimate. The binomial identity
sum_(j=l)^k binom(k,j)binom(j,l)(1/4)^(j-l)
=binom(k,l)(5/4)^(k-l) proves (2). Thus every ACTUALLY established
finite log degree N has an explicit absorption price; N is not assigned
an invented value. The unknown kernel seminorm multiplies this price.

The second gain 7/4 reaches weight five, beyond sigma=4.005. The exact
mass root at three is retained on the first strip; the original quartic
and exact slow root 4+e_* on the second. Subsequent steps stay at sigma
with those phases retained. The earlier derivative-order loss 1/4 is
separate: spatial log absorption and frequency log absorption must both
be paid. With the asserted finite differentiated kernel bounds, the
four/eight POWER counts are unchanged. The exact 3/2-to-7/2 assertion
is superseded. No original admission norm or physical boundary condition
has been changed by these auxiliary proof weights.

## 2. All-frequency inverse on the repaired line, keeping both original rows

Use s=13/4+it, q=s(s-3), r=D_R lambda, beta=lambda B^2/C and

    Jc(s)=-6 beta [Cv-r q](s-4)(s+1),
    B(s)=-K(s-4)q[(q+2)H(s)-L],
    H(s)=psi(s-1)+(pi/2)cot(pi s/2)+kappa_M,
    K=G lambda beta/(15pi), L=10pi Fstar/(G lambda),
    a(s)=(s+3)/(3(s+1)), theta=nu+a(s)b.

These are the SAME complete normal symbols and finite contacts as in
current-coupled-outer-return-test.md. Do not replace H by its logarithmic
principal part. For the original normalized row data (f,g), inversion is

    theta=f/Jc(s),
    b=[(s-4)g-(s-6)f/3]/B(s), nu=theta-a(s)b.        (3)

Both original equations are recovered on this line by substitution.
No differentiated shape row is substituted for the scalar finite-part
matching equation in the physical problem. There is no pole on this
line. The previously established normal strip to weight 7/2 crosses
only the mass root, so it also covers the new line.

The series used in the following elementary estimates are
[DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6) and
[DLMF 4.22.3](https://dlmf.nist.gov/4.22.E3), inspected this turn.
They provide special-function identities, not an operator-domain theorem.

Write a0=9/4=s.real-1. The actual contact enclosure has
7.8<kappa_M<10.8. For |t|<=2, the digamma series gives

    Re psi(a0+it)<=psi(a0)+4[a0^-3+1/(2a0^2)]<2,
    Re psi(a0+it)>=psi(a0)>0, |Im psi|<pi/2.

The last inequality follows from decreasing-series comparison:
integral_0^infinity t/((a0+y)^2+t^2)dy <=Im psi
<=that integral+t/(a0^2+t^2)<pi/2, since a0>1.
Hence |psi|<4. On the entire vertical line |cot(pi s/2)|<=1,
because cos^2(13pi/8)<sin^2(13pi/8). Thus |H|<17 and
|q+2|<14 for |t|<=2; their product is <238<300.

For t>=2 the same decreasing-series comparison gives
0<pi/2-Im psi(a0+it)<=arctan(a0/t). Also
|cot(pi s/2)+i|<=2/(exp(pi t)-1). Consequently

    |t Im H|<=9/4+pi t/(exp(pi t)-1)<3, Re H>5.

For the real part use Re psi>=psi(a0)>0 and
Re cot bounded below by -1; 7.8-pi/2>5. Since
Re(q+2)=45/16-t^2<0 on this range,

    Re[(q+2)H]<=393/16-5t^2,
    |L-(q+2)H|>=L/2+4t^2.                          (4)

Here L>3e34 supplies the last inequality. Conjugation covers t<0.
For all t, |s-4||s||s-3|>=(3/16)(1+t^2)^(3/2).
Together with the low-frequency bound this proves

    |B(s)|>=c3 (1+t^2)^(3/2), c3=3KL/32.

For |t|>=2 use instead |(s-4)s(s-3)|>=|t|^3 and
(1+t^2)^2<=(25/16)t^4. Equation (4) and AM-GM then give

    |B(s)|>=c4 (1+t^2)^2,
    c4=min(c3/sqrt(5), (32/25)K sqrt(2L)).           (5)

That bound alone loses one derivative when applied to (3). The complete
four-derivative return instead needs a fifth-order bound for B. For
|t|<=2 divide the c3 bound by 1+t^2<=5. For |t|>=2 retain the 4t^2
part in (4), and use (1+t^2)^(5/2)<=(5/4)^(5/2)|t|^5. Thus

    |B(s)|>=c5 (1+t^2)^(5/2),
    c5=min(c3/5, 128K/(25sqrt(5))).                 (5a)

Put A=Cv-13r/16. Since A>r>0, the scalar bound is

    |Jc(s)|>=(9/2) beta r (1+t^2)^2,
    C_theta=2/(9 beta r), C_b=1/c5.                 (6)

These are entire-frequency inequalities, not sampled mode tests. For
the usual Fourier H^k(dx) norms, (3)-(6) yield the complete entrywise
bound from the ORIGINAL (f,g) H^k norms to (nu,b) H^(k+4):

       [ C_theta+(275/612)C_b     (25/51)C_b ]
       [          (11/12)C_b                  C_b ].          (7)

This follows from |a(s)|<=25/51,
|(s-4)/<t>|<=1 and |(s-6)/(3<t>)|<=11/12.
Replacing C_b by 1/c4 gives a smaller certificate ONLY from H^(k+1)
inputs to H^(k+4) outputs. Record both, and the extra derivative's price;
the smaller number cannot support the claimed four-derivative bootstrap
with the original input norm.
No datum is relabelled or rescaled to improve the bound. A comparison
of (7) to a graph target is legal ONLY with the required localized source
norms, phase extraction, finite-domain correction and original output map.

The small physical quantum coefficient must remain attached. The physical
canonical area is 2pi N dz/gamma, whereas the original normalized action
pairing has W=-N pi B^2/(G gamma). Thus a fixed-area remainder C E_s has
original canonical source (0,-2GC E_s/B^2). This is exactly the conversion
in current-physical-york-dual-assembly.md, not an adjustable normalization.
Using the shape column of (7) gives the composed NORMAL coefficient

    C_joint=(2GC/B^2) C_b
           =750pi sqrt(5)/(128 beta^2),             (8)

because K=GC beta^2/(15pi B^2) and c5's high-frequency branch is the
smaller one at the actual parameters. The nu coefficient is (25/51) times
(8); the shape coefficient is (8). No source-independent scalar inverse
is multiplied in again. The quantum remainder has zero canonical trace
row, while the actual H inverse and original second row remain coupled.
This exact product has not improved by renaming E_s=R_s/C: its norm and
every contact are still required. Equation (8) exposes the relevant
factor and avoids interpreting the bare 1/K inverse as a physical loss.
The preceding interior 1e48 C bound is only H3->L2 and cannot be inserted
as the missing differentiated E_s bound.

## 3. Quantitative consequence and exact boundary of execution

Execute (1)-(8) on the current parameter enclosure. Record C_theta, both
candidates in (5) and (5a), C_b, every entry of (7), its max-row norm and
the original-scale product (8). These
numbers quantify one necessary factor in the repaired FIRST weight
transition. They are not M_graph,4, M_X,8, kappa_m, or an estimate of the
actual coefficient. The full normal symbol already retains its entire
angular contact; actual thermal and finite-domain returns remain in R_s
and have NOT been evaluated by this normal-line calculation.

To turn this into the requested full bound requires numeric C_jk and
their actual log degrees N_jk for the differentiated, common-restored,
moment-subtracted R_s, plus the actual-H interior/normal gluing constants.
For a coefficient with a known log degree, its factor in that ledger is
C_jk B_(N_jk),k. Neither the earlier unsigned H3->L2 bound nor finite
existence of C_jk gives this number. Normal inversion (7) does not absorb
the compact/interior coupling. Applying it alone would repeat the
already disproved outer-data-only return estimator.

The exact short-travel form cost <3.2e-18 and energy q_s are unchanged.
If the repaired composite constants are bounded, the conditional tail
still has error_ij<=M_i,p q_s^(N-p) U_j. At N=64, the existing half-budget
allowances M_graph,4<=6.0612e157 (affine) and
M_graph,4 U_feedback<=1.5062e156 (feedback) retain their stated scope.
No new value of either left side or of the signed center is computed.

The cut reference mass lift, its local cancellation, transition and
regular-bolt matching are unchanged; it is not a fixed-temperature zero
mode. The original conversion >.00906808164|m| is retained. The mass
target kappa_m<=5e-5 is still UNKNOWN, as are the actual m and other
post-elimination columns. Its conditional residual budgets remain
5.34235e-12 affine and .0012069 graph, with independent energy cost
.0195239 requiring its own comparison. The fast screen
.028913+1e-11(1.0706e10)=.135973<.1599164 retains only its two-fast-input
admitted-family scope. No phase or nonlinear-neighborhood estimate follows.

RESEARCH DECISION: keep the physical construction and short-travel method,
but use the repaired weight path 3/2 ->13/4 ->4.005 and pay (2).
The endpoint estimator fails; this is not a construction obstruction.
The full quantitative estimate remains unresolved, rather than having
an excessive certified kappa_m bound. The coupled numerical solve has
NOT been executed and this does not complete the requested milestone.

ONE NEXT CALCULATION, still proposed: bound the complete original-scale
composition H^-1 R_s, keeping the factor (8) attached to the differentiated
kernel and its explicit log-absorption costs, then certify its strong powers
for all columns. The original finite-part moments, compact/interior gluing
and finite-domain returns must enter that same bound. This is the missing
numerical composite-return calculation, not another normal-symbol or
isolated-contact calculation.

Remaining to-do list: numerical complete coupled return and signed response;
original physical-norm comparisons and nonlinear admission/contraction;
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability and the full stable semiclassical wormhole objective.
