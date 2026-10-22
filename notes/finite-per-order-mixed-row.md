# Mixed row by derivative order and its low inverse composition

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the full physical Hessian at the original finite comparison, the
actual paired realization and physical measure, mass and fixed return.
This combines existing coefficient and lift arrays without changing
the operator or enlarging every derivative to its highest-order value.

## 1. Complete operator and fixed weights

Use [the exact Schur split](schur-source-derivative-budget.md):

    Bbar c=-(8pi g)^-1[Tbar h_c-divbar(Ebar eta_bar_c)],
    Tbar h=sum_(j=0)^2 A_j dot gradbar^j h,
    eta_bar=sqrt(lambda)eta, g=G lambda, W=t_seed^2.

The [local coefficient ledger](finite-higher-local-coefficients.md)
supplies raw full tensor bounds a_(j,k)>=||W gradbar^k A_j||inf,
0<=j<=2,0<=k<=6, and e_k>=||W gradbar^k Ebar||inf,0<=k<=7.
All metric contractions, local Weyl/Ricci and magnetic contacts are
already included. Ebar is the full off-shell trace coefficient. It
is not identified with the bookkeeping scalar residual along the path.

The [finite higher lift](finite-higher-tt-lift.md) supplies matrices

    ||gradbar^r h_c|| <=sum_l H_(r,l)c_l, r<=8,
    ||gradbar^r eta_bar_c|| <=sum_l Y_(r,l)c_l, r<=7,
    c_l=||gradbar^l c||, 0<=l<=8.

The sign Y=-eta_bar in that construction does not affect these norms.
Every norm uses the actual physical volume. W is outside the
differentiated output; in a product it is assigned to the coefficient.
No derivative of W, new end condition or weighted TT inverse is used.

## 2. Ordered product bounds and the divergence contact

At each output derivative order 0<=n<=6 define positive row vectors

    T_(j,n,l)=sum_(k=0)^n binom(n,k) a_(j,k) H_(n-k+j,l),
    D1_(n,l)=sum_(k=0)^n binom(n,k) e_(k+1) Y_(n-k,l),
    D2_(n,l)=2 sum_(k=0)^n binom(n,k) e_k Y_(n-k+1,l).

These follow by differentiating ordered covariant tensors. A metric
contraction is parallel. Permuting tensor slots to identify a product
norm does not commute covariant derivatives or omit a curvature term.
The product rule gives the stated binomial counts even though higher
covariant derivatives need not be symmetric.

For the last two rows expand first

    divbar(Ebar eta_bar)=<gradbar Ebar,eta_bar>
                                         +Ebar divbar eta_bar.

The first contraction costs one; the four-dimensional divergence
costs two. Its derivatives are contractions of ordered derivatives
of eta_bar, so no commuting of divergence past covariant derivatives
is needed. D1 retains e_7 Y_0 at n=6, including the eleventh metric
derivative in the original trace ledger. D2 separately retains e_0 Y_7.

Since 8pi>24, set

    M_(n,l)=[sum_(j=0)^2 T_(j,n,l)+D1_(n,l)+D2_(n,l)]/24,
    m_l=sum_(n=0)^6 M_(n,l), Cmix=max_(l=0)^8 m_l.

Then the complete row satisfies

    ||Bbar c||H_W6 <=g^-1 sum_l m_l c_l
                                         <=(Cmix/g)||c||H8. (1)

There are no omitted source sectors in (1). The derivative-sum norm
makes the maximum COLUMN sum the relevant matrix norm. Keep every
column and the five channel contributions. Compare Cmix with the
previous valid enclosure 1e102 rather than promoting a smaller
individual channel into the full answer.

For comparison, replacing all a_(j,k) by C_loc, all e_k by C_E and
using the shared L8 recovers the valid uniform enclosure
381(C_loc+C_E)L8/24. The new array retains both the different tensor
and return sizes and the derivative allocations that this bound loses.

## 3. Preserve the established low inverse factors

Let c=Dbar^-1 f_c be the original positive critical TT solution and
F0=||f_c||L2/sqrt(gamma0), F4=||f_c||H4/sqrt(gamma0). The
[actual finite realization](finite-comparison-tt-realization-and-range.md)
and [strong-source form estimates](finite-strong-schur-source.md) give

    c_0/sqrt(gamma0)<28g F0,
    c_1/sqrt(gamma0)<265g F0,
    c_2/sqrt(gamma0)<531sqrt(g) F0.                (2)

For the first inequality 8pi/.90<28. For the others put
y=Dbar^-1/2 f_c. Its norm is below 5.3sqrt(g)||f_c||. The existing
form estimates give ||gradbar Dbar^-1/2 y||<50sqrt(g)||y|| and
||Boxbar Dbar^-1/2 y||<100||y||. The scalar Bochner estimate with
|Ric|<15 gives ||Hessbar c||<=||Boxbar c||+4||gradbar c||.
Thus c_2<(530sqrt(g)+1060g)||f_c||<531sqrt(g)||f_c|| because
sqrt(g)<4e-17. These are full scalar covariant derivatives.

Let mhi=max_(3<=l<=8)m_l and retain the unevaluated actual inverse
constant C_D,4 for H4 to H8. Applying (1),(2) yields

    ||Bbar Dbar^-1 f_c||H_W6/sqrt(gamma0)
       <=Klow(g)F0+(mhi/g)C_D,4 F4,
    Klow(g)=28m_0+265m_1+531m_2/sqrt(g).            (3)

This uses the low estimates only on their proved orders. It does
not assign 28g to the higher inverse. Since sqrt(g)>3e-17, use the
explicit upper coefficient

    Klow_bar=28m_0+265m_1+[531/(3e-17)]m_2.

The actual unprojected forcing already has F0<6e-12 ell^2/g, by
27(2.2e-13)<6e-12 and [the finite source conversion](on-shell-finite-range-response.md).
The scalar coefficient stays below 1e20. Hence the SAME full Schur
source obeys the conditional estimate

    N6(p)<=800/g [Cstart ell^2+mhi C_D,4 F4],
    Cstart=1e20+6e-12 Klow_bar.                    (4)

Retain F4 as an unknown physical source norm. At ell=0 the actual
residual vanishes; (4) does not create a nonzero source there. At
ell=1 compare both terms with the retained target
R0 R2 R4 N6(p)<=1e57 epsilon_c. Upward decimal enclosures for Cstart
and mhi give explicit sufficient equal-split budgets. They are not
necessary physical bounds or a proof of nonexistence. The forcing
H4 norm, its conversion, C_D,4 and higher positive transfers still
need numerical estimates before a useful certificate is established.

## 4. Verification and next decision

Freeze this argument, checker and source review. Pin the complete
local coefficients, physical lift matrices and prior scalar evidence.
Check exact ordered Leibniz sums by a separate divided-coefficient
calculation, retain all five channels, verify column norms and
negative controls for the seventh trace derivative and divergence.
Test independently the four-dimensional physical-volume divergence
and its product expansion on nonconstant conformal metrics, including
mixed spatial derivatives through six at 60/90 digits. A polynomial
operator test separately includes A0,A1,A2 and both contact terms.
These are supporting algebra/numerical checks for an ordinary proof;
they do not sample the actual charged quantum kernel or prove a
nonlinear solve. Report the dominant coefficient and resulting
source budget before selecting further refinement or inverse work.

Remaining to-do list: actual finite quantum/TT H4 source and residual
conversion; useful mixed/inverse and positive-transfer estimates;
nonlinear mapping, neighborhood preservation, signed compatibility
and existence; physical validity, state/boundary, matching, stability,
charge/return, opening, transport, preparation, separate creation and audit.
