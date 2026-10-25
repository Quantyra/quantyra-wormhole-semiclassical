# Global comparison attempt: an artificial interior degeneracy

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress (failure of the original-cutoff positive
estimator). Full stable-wormhole objective ACTIVE/INCOMPLETE; audit CLOSED.

OUTCOME: the complete coupled return is still UNRESOLVED. A specified
candidate global preconditioner fails an ACTUAL continuum defect test,
not merely an oversized upper bound. This is an obstruction to that
comparison, not to the physical Hessian or wormhole construction.

## 1. Construct the candidate before assuming an inverse

The preceding .002555496 remainder bound applies only at the original
.05/.1 graph cutoff, relative to the normal inverse. It did not construct
a finite-geometry inverse. Test the literal global extension of its
suggestion to absorb the cubic and quartic local matrices:

    A34=A0+(Ac)_3+(Ac)_4, A=z^4 H, Ac=z^4 Hc.        (1)

Here (Ac)_p is the complete coefficient of z^p in BOTH original local
rows, before the differentiated shape combination. A0 is the complete
original normal operator, including its logarithmic/contact shape part.
The proposed inverse would retain the original outer and regular-bolt
conditions. Formula (1) is an explicitly tested comparison, not the
physical equations. A collar parametrix patched to a different interior
operator is a different candidate, not ruled out by this test.

The first unsupported step is asserting that (1) has a conforming global
strong inverse. Its highest-derivative coefficient must be checked before
bounding quantum remainders or using the collar allowance globally.

Keep the exact reference mass lift and its .1/.2 cutoff, original boundary
data and bolt matching. Its two reference rows cancel locally, not globally;
its original graph conversion still costs 2133..2174 per unit coefficient.
None of these quantities or the original admission radius is changed here.

## 2. Both original equations and the full conformal column

Work in the retained affine canonical action chart (nu,b). Its exact
common-restored action is the local action plus Gamma_hat[b]. The latter
contains the complete angular/thermal tower, LLL and finite-domain returns.
It has no nu dependence. This is the action identity proved in
current-finite-local-action-assembly.md, not an assertion that physical
quantum stress has zero trace. The FULL Weyl anomaly remains local and
must be included in the nu column.

For the tests below b=0 and nu is compactly supported in .6<z<.63. It
vanishes near both the paired boundary and bolt. Thus the original paired
projectors, complementary adjoint, mass/log phase and regular-disk traces
are unchanged. The canonical massive operator is unchanged; the common
regulator's conformal conversion contributes exactly the retained local
anomaly. Boundary representatives vanish on this compact support. This
uses the established action identity on a bounded interior C5 family;
it does not extend a seed-domain theorem to the unknown nonlinear ball.
No angular or thermal sum is truncated, and no quantum entry is discarded:
its nonlocal nu column is identically zero, while its anomaly column is kept.

Put c=8piG c_R C/B^2 and p=GC/(2pi B^2), in the unchanged normalized
action units. The principal second-variation densities contain

    c Q^2[72 nu_i''nu_j''+24(nu_i''b_j''+b_i''nu_j'')+8b_i''b_j''],
    p Q^2[-nu_i''nu_j''/5-(nu_i''b_j''+b_i''nu_j'')/15+...].

The second line is the Weyl contribution. Consequently the leading D^4
coefficient of the COMPLETE nu column, D=z partial_z, is

    a4(z)=(72c-p/5)Q(z)^2,  second original row = a4(z)/3. (2)

The positive scalar coefficient agrees independently with
6D_R C/B^2 times Q^2 in the covariant conformal density. This verifies
the anomaly sign and factor, not just the symmetry of a matrix.

For completeness the full Weyl nu-nu strong row used in the calculation is

    W_nu nu=p[-Delta^2 nu/5-partial_z(d_W partial_z nu)],
    Delta=partial_z Q partial_z,
    d_W=-Q''Q/15+7Q'Q/(45z)-7Q^2/(15z^2)-uQ/9, u=1/C.

This follows directly from the local action: the mixed Delta nu/nu'
bilinear combines to -partial_z[7Q^2 nu'/(45z^2)]. All other local
nu-column terms are the exact Einstein/cosmological and R2 rows from
that action. No independently adjustable quantum source is introduced.

## 3. The global truncation loses an interior principal coefficient

Write Q=beta q, q=1+q2 z^2+q3 z^3+q4 z^4+q5 z^5+q6 z^6,
with the ACTUAL finite-temperature polynomial. Set

    k=72c beta^2=576pi G c_R lambda beta,
    d_a=1/(5760pi^2 c_R)  (G=1 units).

Since (1) restores only the classical p=3,4 rows, its principal coefficient is

    a34_4/k = 1-d_a+2q3 z^3+(2q4+q2^2)z^4 =:v(z),
    a4/k=(1-d_a)q(z)^2.                              (3)

Both original comparison rows have the SAME zero, with ratio 1:1/3.
The complete normal quantum term acts only on b and cannot repair it.

Whole-interval arithmetic using the retained parameter enclosures gives

    .61427243184514 < z0 < .61427243184516,
    v(z0)=0, -3.494<v'(z0)<-3.493,
    .2499 < (1-d_a)q(z0)^2 < .2501,
    7.36e-10 < a4(z0) < 7.37e-10.                    (4)

The zero is unique between .61 and .62: v has opposite endpoint signs
and its derivative is negative throughout that interval. It is far from
the regular bolt. On the outer cutoff v(.1)>.9925, but v(.7)<-.3033.
Thus a collar bound alone conceals this global degeneracy. The omitted
p=6,7,8 square terms are essential: even at q2=q5=q6=0 they complete
(1-4z^3+3z^4)^2. They cannot be discarded globally because their
localized norm was small. Reinstating the ENTIRE q^2, with the anomaly,
removes the artificial zero; it does not yet prove a global inverse.

## 4. Execute a continuum defect test, including both rows

Let x0=log z0, psi(t)=(1-t^2)^10 for |t|<1 and zero otherwise, normalized
in L2(dt). This is C9 and H10 and lies in the completed strong test
domain; smooth compact approximations give the same strict estimates.
Define the ONE specified continuum field, not a finite physical mode sum,

    nu_N(x)=N^(-27/4) psi(sqrt(N)(x-x0)) exp(iN(x-x0)), b_N=0,
    N=10^8.                                          (5)

Its whole support lies in .6<z<.63. The real operator can be complexified
without changing its Hilbert-norm bound. All outer and bolt jets vanish.
Finite C5 norms permit a small path in the same admitted action domain;
the field is a test direction, not a proposed physical excitation or
chosen residual. No smooth surrogate boundary is used.

Here is an explicit error bound for (5), without spatial sampling.
Let n_j=||psi^(j)||2/||psi||2, j=0..7, calculated by exact polynomial
integration, and

    U_l=N^(l-7) sum_(j=0)^l binom(l,j)N^(-j/2)n_j.

Then ||D^l nu_N||2<=U_l and ||D^7 nu_N||2>=2-U_7. For a coefficient
f(z)=sum f_p z^p on z<=.63, use the WHOLE-interval jet guard
C_j(f)=sum |f_p| p^j .63^p. All coefficients here are actual finite
polynomials. For a comparison leading coefficient use C_1(f)/sqrt(N)
in place of its value guard because f(z0)=0 EXACTLY. The center is the
implicit root, not a rounded decimal substituted for it.

For each original comparison row sum the product-rule bounds

    L_h=sum_(l=0)^4 sum_(j=0)^h binom(h,j) C_j(a_l) U_(l+h-j),
    ||row||H3 <= sqrt(sum_(h=0)^3 binom(3,h)L_h^2).    (6)

Use Y3=H3(dx) direct_sum H3(dx) with its Hilbert norm; the sum of the
two individual row bounds is an upper bound for that norm. No row is
rescaled. For the actual scalar row, take the D^3 output and retain its leading
a4(z0)D^7 nu_N. Subtract the variation C_1(a4)U_7/sqrt(N) and EVERY
remaining Leibniz term in (6). Exact polynomial integration and 512-bit
outward arithmetic prove

    ||A34(nu_N,0)||_Y3 <3.853e-12,
    ||[A(nu_N,0)]_nu||H3 >7.27e-10,
    ||(A-A34)(nu_N,0)|| / ||A34(nu_N,0)|| >187.        (7)

The numerator lower bound uses its scalar row; the denominator upper bound
includes BOTH original rows. These are bounds on this actual operator
application, unlike the preceding positive-matrix upper-bound test.
If P34 were a two-sided inverse of A34 on this domain, (7) would imply
||I-A P34||>187 in the joint original-row H3 norm, failing delta<=1/2.
If that inverse does not exist, it is already unavailable as a preconditioner.
This tests the proposed unrestricted source-space norm, not the actual
restricted affine/feedback source family. It does not assert that the
prescribed quantum stress can be chosen to equal A34(nu_N,0).

More generally the SAME proof with N tending to infinity gives
||A34 u_N||=O(N^-1/2), while ||A u_N|| has a positive lower limit. Thus
the relative strong defect is unbounded for this comparison. Any fixed
smooth positive physical weights on this compact interval are equivalent,
so a normalization change cannot restore its boundedness. This statement
concerns a global strong-data defect. It is NOT a bound for the localized
outer mass column kappa_m, and does not prove an actual physical instability.

## 5. Repair and remaining complete calculation

The repair is to retain the FULL geometric local operator in a global
comparison, not just the powers that dominated one collar estimate.
The already defined short-travel geodesic comparison H_s does this,
and keeps the complete logarithmic principal part, original fourth
contact, paired completion and regular bolt. Its definition in
current-short-travel-strong-domain.md is unchanged. The established bounds

    ||H-H_s||_(V->V*)<=eta_s=.0022000000000000032,
    ||H_s^-1||_(V*->V)<=1/(1-eta_s)<1.002205

are valid FORM bounds, not the missing strong graph/phase bounds.
They justify retaining H_s as the conforming comparison candidate;
they do not certify the requested complete defect or signed correction.

The first remaining strong step is a numerical estimate for the complete
restored near kernel/contact operator composed with this geometric
inverse and original phase/matching maps. Its entire angular/thermal
sums and finite returns remain uncomputed. No theorem is imported here
to assign those unknown seminorms or extend them to a nonlinear ball.
The exact mass lift, scalar quartic's undifferentiated row, correlated
slow pair and independent boundary-phase columns remain in that solve.

The full linear mass <1.627532e-15 and distant shape <1.24e-20 remain
established. The >6.4617454e-12 near allowance is SHAPE ONLY. kappa_m<=5e-5
is still conditional; its remaining affine and feedback budgets are
5.34234999994958e-12 and .0012069, including the independent energy debit.
The .135973<.1599164 quantum screen remains restricted to its two-fast
admitted-family scope. No physical-norm/neighborhood comparison follows.

DECISION: retain the construction; reject the globally truncated A34
inverse. Use the full geometric comparison already available, preserving
all local polynomial powers and contacts. The failure (7) is specific
and verified; the complete physical return remains unresolved.

ONE NEXT CALCULATION, PROPOSED NOT EXECUTED: apply the full H_s inverse
in the original bordered graph and certify the COMPLETE common-restored
strong defect and output/phase maps, including the actual near quantum
response. Test its signed affine and feedback errors against the unchanged
budgets. No further cold/scalar-source or isolated contact improvement is
selected. This turn executed (7), not that full inverse calculation.

Remaining to-do list: complete coupled inverse/defect and independent phase
columns; original scalar graph and physical-norm/neighborhood comparisons;
nonlinear admission/contraction; exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole objective.
