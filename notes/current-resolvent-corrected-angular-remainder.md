# Exact resolvent corrections: a complete angular/thermal remainder bound

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
The full stable semiclassical wormhole objective remains ACTIVE/INCOMPLETE.
The closed feasibility audit is preserved.

OUTCOME: the complete physical coupled return remains UNRESOLVED. Two
EXACT resolvent corrections repair the previously divergent auxiliary
split: its remaining operator and first conormal commutator now have
absolutely convergent complete angular AND thermal sums in L2. The
retained common-regulated core and original strong coupled observation
still need evaluation. No kappa_m, actual correction or solution gate is
certified. This is a bound for an entire operator remainder, not an
isolated contact or a sampled spectral calculation.

## 1. First unsupported step and the scope of this repair

The predecessor proved that separately summing H_M-H_* fails even after
the actual weak coupled inverse. Its complete thermal large-mass limit
C_aux=-sqrt(Q)Q''/120-23(Q')^2/(1680sqrt(Q)) is nonzero. Subtracting a local
large-mass coefficient alone does not supply an operator-norm remainder:
that asymptotic was proved on fixed smooth test directions, without
uniformity in arbitrarily short radial wavelengths.

We therefore retain the full frequency dependence in the exact scalar
resolvents. This changes only the auxiliary decomposition. Keep the actual
paired mass/log domain, complementary adjoint and composed Dirac square,
fixed period, both original gravitational equations, regular bolt,
original contacts and common physical heat prescription. The first
unsupported physical step remains quantitative evaluation of the
common-regulated retained core followed by the strong coupled inverse
and original matching rows. This note does not infer that result from
the new L2 bound.

## 2. Exact decomposition on the actual endpoint realization

Use the same notation and proved hypotheses as
[the paired response comparison](current-paired-response-parametrix-bound.md):

    H_M,p=-Mt-4M^2 R, R=L^-1,
    L=-partial_x(P partial_x)+V,
    P=1/(Wa), V=4W/a+4Mt/a^2,
    dz/dx=f=sqrt(Q), W=p/f, K^2=M^2+p^2/.17.

The original graph cutoff chi is one through .05 and zero by .1, and

    P_*=P+chi(P0-P), V_*=V+chi(V0-V),
    P0=k/W^2, V0=4k^3/W^2, k^2=M^2+W^2,
    Delta=L_*-L, S=L_*^-1.

Both scalar forms have the same natural Neumann seam condition with the
value free, and the actual finite-form regular-bolt branch. Coefficients
coincide outside the collar. The cutoff transition is retained in Delta.
It is NOT the independently retained mass cutoff chi_m. No wall or
smooth extension across the cubic seam is introduced.

For every actual M_n>18 and p_j=gamma(j+1/2)>0, retain

    ||R||<=c0/K, c0=sqrt(2)/[4(1-sqrt(1.05)/36)],
    ||L^-1/2 Delta L^-1/2||<=delta=1.85/K^2<1.       (1)

Write E=L^-1/2 Delta L^-1/2 and F=(I+E)^-1. These are bounded form
operators; E is self-adjoint and F is its functional calculus. Products
such as S Delta S denote their bounded form extensions, not composition
of unbounded differential expressions on unspecified domains.

Define the retained core and the exact error per mode by

    H_core=-Mt-4M^2(S+S Delta S+S Delta S Delta S),
    T_M,p=H_M,p-H_core.

The resolvent identity, iterated twice, gives

    T_M,p=-4M^2 S Delta S Delta S Delta R
          =-4M^2 L^-1/2 E^3 F^3 L^-1/2.             (2)

Indeed F+F^2 E+F^3 E^2=I-E^3 F^3. E and F commute in this bounded
functional calculus; no commutation with spatial differentiation is
assumed. The original tadpole -Mt is kept EXACTLY in H_core. In terms
of the old H_* the core includes its missing tadpole difference and
both displayed resolvent corrections. No coincidence term is erased.
All propagation through the finite paired geometry occurs in R and S.

Equation (2) has the uniform all-frequency bound

    ||T_M,p|| <= A M^2/K^7,
    A=4c0(1.85)^3/(1-1.85/324)^3 <9.377.           (3)

This is not a fixed-radial-frequency expansion. The full resolvents are
kept, and (1) bounds the remainder on all proper-radial L2 inputs.

## 3. First conormal commutator, including the full cutoff cost

Use the SAME X=chi z partial_z. It preserves the seam condition and
vanishes before the bolt. The predecessor established the complete form
commutator estimates

    B=L^-1/2[X,L]L^-1/2, ||B||<996,
    C_X=L^-1/2[X,Delta]L^-1/2, ||C_X||<1070/K^2.     (4)

They include the first-order term -P h'' partial_x for X=h partial_x
and all derivatives of the original cutoff. No derivative is obtained
by differentiating an unsigned numerical envelope.

Commute the first expression in (2) on the original regular core.
The three S factors each give -S([X,L]+[X,Delta])S; the three Delta
factors give [X,Delta]; the last R gives -R[X,L]R. Conjugating each
product by L^1/2 and using (1),(4) yields

    ||[X,T_M,p]|| <= A_X M^2/K^7,
    A_X=4c0{3a^3(996+1070/324)/(1-a/324)^4
             +3a^2(1070)/(1-a/324)^3
             +a^3(996)/(1-a/324)^3}<53879,
    a=1.85.                                       (5)

All seven differentiated terms are included. This argument requires only
the previously proved first commutators, not unknown higher coefficient
jets or a smooth-boundary theorem. The bounded form identities extend
the core result to the indicated operators. It supplies ONE commutator;
it is not the H3 graph estimate or a metric-neighborhood variation bound.

## 4. Complete thermal and angular sums

For the actual fixed positive gamma, the decreasing-function bound gives

    gamma sum_(j>=0) M^2/[M^2+p_j^2/.17]^(7/2)
       <= M^-4[sqrt(.17)*8/15+gamma/M].             (6)

The integral is evaluated over the entire half-line; its dimensionless
value is int_0^infinity(1+s^2)^(-7/2)ds=8/15. The actual mesh endpoint
correction is retained, with gamma<1e-33 and M>18. Thus

    ||gamma sum_j T_M,pj|| <2.062/M^4,
    ||gamma sum_j [X,T_M,pj]|| <11848/M^4.           (7)

Both thermal signs are included by the original action factor in section
5. No physical zero-temperature replacement or finite mode list is used.

For M_n^2=n(n+q)/C, d_n=2(q+2n), the exact entire-tower estimate is

    sum_(n>=1)d_n/M_n^4 <4C^2/q.                    (8)

One direct proof uses the decreasing function
F(x)=(q+2x)/[x^2(x+q)^2], with integral_1^infinity F=1/(q+1), and
F(1)+1/(q+1)<2/q. Multiplication by 2C^2 gives (8).

Consequently the double series for T=sum_n d_n gamma sum_j T_Mn,pj
and its FIRST commutator converges absolutely in operator norm. Tail
control covers every n and j. Closedness of X identifies the summed
commutator on Dom X. The actual bounded mass/log directions are included
in L2; this extension is not a nonlinear family of L2 metrics.

## 5. Original source units, common regulator and the remaining gap

For near input u=1_[0,.2]T_phys the established comparison is
||u||L2(dx)<35.534 E(u_phys). On the original output z<=.1, f>.4 and
the conversion of f^-1 output from L2(dx) to L2(dz) costs sqrt(2.5).
The physical pair-of-thermal-signs source factor is -2G/(pi B^2 f).
Using the actual C,B,q,

    A4=(2G/(pi B^2))(4C^2/q)<4.583e-35,
    ||source_tail||L2(dz) <5.31e-33 E(u_phys),
    ||[X,source_tail]|| <3.052e-29 E(u_phys).        (9)

For the last operator bound include the density commutator:
|X log f|<.0625, so its coefficient is bounded by
A4 sqrt(2.5)*35.534*(11848+.0625*2.062). The input indicator commutes
with X, since X vanishes before .1; it does not create a metric wall.
The last line denotes the induced commutator on the near-input operator,
not a bound for X applied to arbitrary L2 inputs without their X norm.

These are original canonical L2 source bounds. In particular neither
2.062 nor 5.31e-33 is the missing composed C4 or a kappa_m value. Passing
from these bounds to T_old H^-1 and its H3/phase/bolt observations still
requires the actual structured smoothing and higher commutators. A bare
bounded L2-to-original-graph inverse is not assumed. The previously
controlled far-input graph bound <4.583e-13 E remains in its different
norm; it has not been added to (9) as if the norms agreed.

Common-regulator scope is precise. T is an absolutely convergent bounded
operator at the fixed seed. It can be subtracted from the original
regulated physical form H_e, defining the companion H_core,e=H_e-T
with the original action factors and all other sectors retained. This is
an exact bookkeeping identity, so the established admitted-family WEAK
limit gives H_core,ren=H_ren-T without changing any finite contact.
For any finite set of separated labels, removing the original regulator
recovers the resolvent identity (2). We do not replace that regulator
by exp(-e L_*) or claim the bare infinite sum of H_core is already
renormalized. The common finite-part/restoration terms remain in its
companion. No new determinant or quantum state is assigned to L_*.

This establishes a legitimate finite bounded subtraction from the weak
physical operator. It does NOT execute the common-regulated core or
prove uniform metric variations of this split on a nonlinear ball.
The earlier C2 admitted-family result retains only its proved scope.

## 6. Coupled decision and unchanged original requirements

The full physical solve still includes both original equations, the
scalar quartic's undifferentiated row, exact mass compatibility, the
correlated slow pair, fast remainder, interior and regular-bolt matching.
The original mass lift remains J_m=int_0^z t^3/Q^2,
n_m=-2beta QJ_m/z, b_m=beta z^3/Q+beta Q'J_m. Its local reference rows
cancel; its original transition is f_cut=2xz chi_m-2F n_m Dchi_m and
alpha_cut=n_m Dchi_m. Its reference bolt shift is zero and net mass -xB.
It is not a global fixed-temperature zero mode. The original-fast cost
>.00906808164|m| and required actual |m|<2.20553814952e-9 remain unchanged.

Actual m and kappa_m remain UNKNOWN. IF kappa_m<=5e-5, its affine, graph
and energy debits are 2.14765e-12, .0024731 and .0195239. After the
controlled far entry the remaining budgets are still

    Y_other<=5.34234999994958e-12,
    eta_other+(.0195239+4.583e-13)gamma_E<.0012069.

The source-family gamma_E and all remaining coupled entries are unknown.
No bound on them is inferred from (9). The previous composed C4<=1e33
allowance, conditional on the ORIGINAL composed norm, has not been met
by renaming the new L2 constant. Its stricter conditional gamma_E<.01846
alternative remains an allowance only. The two-fast-input screen
.028913+1e-11(1.0706e10)=.135973<.1599164 keeps its admitted-family scope.
Original physical-norm/neighborhood comparisons and nonlinear conditions
Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1 remain unproved.

EXECUTED: exact full-resolvent decomposition, actual-domain all-frequency
error and first commutator, COMPLETE norm-convergent angular and thermal
sums, original source conversion and a weak-regulator-compatible bounded
subtraction. No physical signed coupled center, strong composed bound,
independent phase coefficient or nonlinear solution is computed.

RESEARCH DECISION: retain the construction and replace the failed
separate auxiliary angular estimate by this two-resolvent-correction
decomposition. Keep its exact tadpole and corrected core in the original
common subtraction. The L2 angular remainder is now controlled; further
optimization of that value bound is not the next feasibility task.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): evaluate the common-regulated
retained core together with its structured remainder in the ORIGINAL
coupled strong return T_old(H_s^-1 R_s)^4, including higher commutators,
mass/quartic/slow/bolt rows, signed columns and certified regulator error.
Test those actual entries against the unchanged affine/feedback budgets.
The unknown core and strong composition, not the now-summable L2 tail,
are decisive for whether an actual controlled solution is attainable.

Remaining to-do list: complete common-subtracted strong coupled return and
signed response; original physical-norm/neighborhood comparisons;
nonlinear admission/contraction and exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable semiclassical
wormhole objective.
