# Coupled local graph feedback: exact cancellation and a quantitative bound

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The full stable semiclassical wormhole objective remains active and incomplete.
The feasibility audit stays closed. This is an actual finite-seed operator
component estimate, not a solution, a nonlinear-family theorem or a source
chosen to cancel the geometric residual.

## 1. Scope and the domain that is retained

Keep the paired boundary realization, common regulator, prescribed field
theory/state, finite temperature, mass/log fields and regular bolt. The action
and its first two variations are available on bounded ADMITTED stationary
weighted C5 families, as specified in current-curved-paired-heat-bridge.md.
Admission of the proposed nonlinear X ball is still unproved. The fixed-seed
strong graph estimate is also incomplete. No theorem from the Navier--Stokes
or Jacobian--Weyl projects supplies either missing implication.

Here the local Einstein/cosmological, magnetic and additional R2 Hessians
are taken from current-finite-local-action-assembly.md. Call their sum Hc.
It includes the second curvature variation and the off-shell density terms.
All Weyl anomaly and paired determinant terms remain in the complementary
operator; they are NOT dropped. In particular it would be wrong to bound the
explicit outer log z Weyl contact alone: its cancellation with the mass
determinant is required first. The earlier .0022 weak-form bound is not a
strong graph estimate.

Let A_c=z4 Hc, and A_c0 be its z-independent Euler normal operator. This
definition fixes the normal component of THIS action sum, without changing
any coefficient of the full normal operator A0. On x=log z use

    sigma=4.005, D=z d_z, a(s)=(s+3)/(3(s+1)),
    theta=nu+a(D)b,
    T(V)=(D-4)V_b-(D-6)V_nu/3.

The input normal responses are theta=Jc_sigma^-1 f and b=B_sigma^-1 g,
where f and g are the original normalized scalar and differentiated shape
rows, respectively; there is NO additional row rescaling. The complete
normal multipliers are

    Jc(s)=-6 beta [Cv-r q(s)](s-4)(s+1),
    r=D_R lambda, q=s(s-3), Cv=1+11G lambda/(180pi),
    B(s)=-K(s-4)q[(q+2)H(s)-L].

All parameters, finite contacts and H(s) are those of the preceding full
normal proof. In particular K=G C beta^2/(15pi B^2) and L=10pi Fstar/(G lambda).
The relation lambda B^2=C beta fixes the normalization of Jc.

We bound the four entries of

    chi ( [1,0; -(D-6)/3,D-4](A_c-A_c0)
                    [1,-a(D);0,1] ) diag(Jc^-1,B^-1),             (1)

after conjugation by z^sigma. The row differential operator acts BEFORE
the output cutoff. This is the same a=.01 output cutoff used previously:
chi=1 for x<=log(a)-1, chi=1-S(v) for v=x-log(a)+1 in [0,1], zero afterwards,
S=35v4-84v5+70v6-20v7. It is not an input cutoff, a physical wall or a
replacement for the original domain. Normal inverses on the full comparison
line are estimators; original finite-dimensional, cutoff, remote and bolt
data must still be bounded in the actual global problem.

## 2. Exact differential cancellations, before taking absolute values

Write c=8piG c_R C/B^2, u=1/C, ell=lambda0 B^2/C and
mag=piG q_flux^2/(g_R^2 C B^2). For a direction (n,b), let

    Delta v=(Q v')', D0=-Q''+6Q'/z-12Q/z2,
    d=-2Delta b-6Delta n+12Q n'/z, v=d-2D0 b.

The additional R2 rows, divided by c, follow by two ordinary integrations
by parts in the bilinear density (5) of the local-action note:

    F_n=-12(Qv)''-[(24Q/z-12Q')v-12Q(4u+2D0)n']',
    F_b=(16u2+4D0^2)b-4D0 d+4(Q'v)'-4(Qv)''.

These formulas apply to arbitrary direction jets, not just powers. In
T(z4 F), the fifth AND fourth derivatives cancel identically for arbitrary
Q. The surviving coefficients of n''' and b''' are, respectively,

    -16 z4 Q [u z+z Q''-3Q'],
    -8 z3 Q [z2 Q''-4z Q'+4Q].                       (2)

Thus the scalar-to-shape R2 row is at most third order. In the input
combination n=-a(D)b, the first row is at most order two (a rational Euler
operator with a regular D+1 denominator). This cancels another spurious
high-order loss. The Einstein and magnetic rows are included with their
exact signs before collecting terms. A direct differentiation of the
original weak densities supplies an independent check of these identities.

Use the full actual polynomial, including the correlated temperature:

 Q=beta+u z2+(-4beta-2u-tau)z3+(3beta+u+epsilon)z4
                                             -2epsilon z5+epsilon z6.

After subtraction of A_c0, every entry in (1) before the normal inverse is
sum_(p=2)^12 z^p c_ij,p(D). The scalar-input column has polynomial degrees
at most 4 and 3; the shape-input column has degrees at most 2 and 3 plus
a remainder divided by D+1. These exact coefficients are recorded in the
executed evidence. No angular or thermal spectrum is truncated by this
finite polynomial representation of the LOCAL action.

The boundary representatives of the local Hessian are the original ones.
The paired mass/fast tangents have vanishing second-variation endpoint
flux, as established in the local-action note; the first-variation boundary
terms need not vanish and are not omitted from the global source. Local
coefficient estimation supplies no new boundary condition or adjoint domain.

## 3. Whole-frequency multiplier bounds

Put delta=.005 and A=Cv-r sigma(sigma-3)>0. On s=sigma+it,

    |Jc(s)| >= 6beta(A+r t2)|s-4||s+1|,
    |B(s)| >= (K L/2)|s-4||s||s-3|.                 (3)

The second inequality follows from the previous low-frequency bound
|(q+2)H|<5000 for |t|<=4 and the high-frequency real-part bound for |t|>=4.
It holds at EVERY frequency, including all contacts. For |t|<=1, (3) gives

 |Jc^-1| <= [6beta A(sigma+1)]^-1 /(delta2+t2)^1/2,
 |B^-1| <= [K L sigma(sigma-3)/2]^-1 /(delta2+t2)^1/2.

For |t|>=1, the bounds on |t|^j times the inverse, j=0,...,4 for Jc, are

 [1/(6beta A), 1/(6beta A), 1/(6beta A),
                         1/(12beta sqrt(A r)), 1/(6beta r)].       (4)

For B and j=0,...,3, every entry can be bounded by 2/(K L).
Use |s+1|>=sigma+1 for the remaining rational term. Expand each polynomial
as c(sigma+it) BEFORE absolute values. Sum its coefficient magnitudes times
(4) for a certified high-frequency multiplier bound h_p. On |t|<=1,
the sum of coefficient magnitudes, including the rational remainder,
times the corresponding low-frequency constant is a bound l_p for
|c(sigma+it)/Jc(s)| sqrt(delta2+t2), or the analogous B expression.

These are analytic whole-frequency estimates. Frequency samples are neither
the proof nor a finite selection of physical modes.

## 4. Localization resolves the excessive crude bound

For alpha_p=chi exp(p x), let

 M_p=sum_(j=0)^3 binom(3,j)||alpha_p^(j)||_infinity,
 L_p=sum_(j=0)^3 binom(3,j)||alpha_p^(j)||_L2(dx).

The usual H3 product rule gives ||alpha_p v||_H3<=M_p||v||_H3.
A useful stronger statement for the low-frequency multiplier R is

 ||alpha_p R 1_(|t|<=1)||_(H3->H3)
        <= L_p ||R 1_(|t|<=1)||_L2(dt)/sqrt(2pi)
        <= L_p l_p/sqrt(2delta).                    (5)

To prove (5), insert the input weight <t>^-3 into the Fourier kernel and
differentiate the output k times. For l derivatives on alpha, the squared
derivative weights are
sum_(k=l)^3 binom(3,k)binom(k,l)^2 t^(2(k-l))/<t>6
<=binom(3,l)^2/<t>^(2l). The Hilbert--Schmidt kernel bound and the triangle
inequality yield (5), using integral_R dt/(delta2+t2)=pi/delta.
This controls a narrow low-frequency response by its integrated size,
rather than the crude maximum 1/delta. It changes neither units nor domain.

For each matrix entry sum the low bounds (5) over p and sum M_p h_p for
the high-frequency part. Orthogonality of the two input Fourier ranges
then gives sqrt(low_total^2+high_total^2) as a bound on that entry.

All cutoff constants are certified on the whole transition interval.
Let R_j(v)=(d_v+p)^j(1-S(v)). On each of 128 rational subintervals,
the Bernstein convex hull bounds |R_j|; multiply by the exponential at
the RIGHT endpoint to bound exp(p(v-1))|R_j(v)|. Let this maximum be m_j.
The left tail contributes p^j exp(-p) to the supremum and
p^(2j)exp(-2p)/(2p) to the squared L2 norm. Thus the respective bounds are

 a^p max(m_j,p^j exp(-p)),
 a^p sqrt(m_j^2+p^(2j)exp(-2p)/(2p)).                (6)

Equations (1)--(6) constitute the ordinary argument. Exact algebra and Arb
interval evaluation verify its hypotheses and constants; passing checks
cannot establish the missing full operator/domain theorem.

## 5. Decision and quantitative requirements

The crude all-frequency supremum/product estimate is also recorded. Its
failure is an EXCESSIVE UPPER BOUND, not an actual local obstruction.
The frozen acceptance target for the sharper local block is maximum row
sum <.06 on the SAME .01 collar with unchanged physical normalization.
Report all four entries, both estimates and their difference in the evidence.

To use this block in a complete closure, the remainder must control BOTH
original rows in E=max(||f_R||_H3,||G_R||_H3), after the same prescribed mass,
quartic and slow-mode decomposition. It is not legitimate to treat the
unknown scalar residual as zero in a shape-only inequality. In particular
E3=||G_R||_H3<=E, but its preceding target alone does not bound f_R.

If eta_c is the certified local row-sum bound, the remaining operator
target under the previous allocation is eta_remaining<=.24-eta_c.
Here eta_remaining includes every quantum/nonlog contact, variable-
coefficient/curved and separated kernel contribution, and the feedback
of original cutoff, mass, scalar, slow and bolt/global sectors. The known
log-principal component eta_pr<.007302177 is separately retained. A complete
coupled estimate must establish

 E <= S_coupled+(eta_pr+eta_c+eta_remaining)(E+k_coupled).            (7)

Both the joint affine bound S_coupled and the known-data norm k_coupled
must be evaluated from the original equations. The previous k3<1.613e-27
is ONLY the known shape-source norm; it is not automatically k_coupled.
For explicit conditional targets choose S_coupled<=7.49e-12 and
k_coupled<=1e-14. Substitute the certified eta_c and eta_remaining<=.187
in (7), and check whether it implies E<1e-11. None of these still-missing
joint source/global estimates is supplied by the local calculation.
The independent Efar<=1e-10 requirement remains unknown.

RESEARCH DECISION: retain the construction and use the exact coupled
cancellations and the low-frequency integral bound when completing the
strong estimate. Do not discard a route because its crude bound exceeds
one; the sharper, executed test must decide this local component. Do not
optimize the already controlled cold or scalar source.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): derive and bound the
remaining Weyl-plus-paired-determinant operator in this same coupled graph,
with explicit outer log cancellation, full angular/thermal tails and the
original-domain remote/boundary terms, against the .187 remaining feedback
allocation. This is the outstanding quantum strong-operator calculation,
not another normal-mode or weak-form screen. Its failure must be classified
as an unknown estimate, excessive upper bound or actual obstruction.

The earlier sufficient quantum/contact screen E_rest<2.123e-42<1.287e-40
retains its scope and unmeasured sign; it is not necessary for existence.
No actual E, E3, Efar or complete correction has been certified here.
At rho=1e-8 the unresolved preserved-ball and contraction conditions remain
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with unproved allocations
Y0<=5e-9, Z0<=.25 and Z2<=5e7. Exact compatibility, physical matching and
causal realization, semiclassical validity/noise and all-sector stability
remain open. The full stable-wormhole objective is preserved.

Remaining to-do list: remaining quantum/global graph estimate and joint
source return; actual controlled correction and nonlinear admission; exact
compatibility; physical matching/causal realization; semiclassical validity/
noise and all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
