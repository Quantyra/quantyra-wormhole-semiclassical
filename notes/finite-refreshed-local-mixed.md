# Refreshed finite local coefficients and complete mixed source budget

2026-10-22, E003/S167. Ordinary argument before supporting execution.
Keep the original finite comparison, full action, physical volume,
paired domain, angular tower, mass and outer-fixed return. This updates
the local coefficient inputs to the [complete per-order mixed row](finite-per-order-mixed-row.md).
It does not assume a nonlinear solution or a quantum H4 source bound.

## 1. Finite fields with an explicit order-eleven reserve

Use the reference derivatives D=z0 nabla_gamma and W=t0^2, fixed
at each observation as in the [original local ledger](finite-higher-local-coefficients.md).
All polynomials below contain divided D-derivative bounds, with the
100^j factor already included. The physical finite fields are

    U_a=a u, S_a=log(1+tau(exp(a V)-1)), 0<=a<=1,
    partial_a S_a=nu exp(a V-S_a), S_0=0.

The [low-gap reconstruction](finite-low-gap-reconstruction.md) supplies
unweighted/weighted tangent arrays U,K and the unweighted Vbar through
order ten. Its six finite integrations supply Tplus, bounding positive
derivatives of both U_a,S_a through ten. The separate finite value
bound remains 3e-8. Put Vplus_0=0, Vplus_j=Vbar_j for j>=1.
Through order ten only, form

    S1W=1.001 exp(Vplus+Tplus) K.

This bounds W times every divided spatial derivative of partial_a S_a,
uniformly in a. The small value estimates justify the 1.001 prefactor.
Integration from zero to a commutes with the fixed reference
derivatives. Hence a valid simultaneous weighted bound for U_a,S_a is

    Kfield_j=min(Kold_j,max(K_j,S1W_j)), 0<=j<=10.

Here Kold is the original finite-field table. This is a bound on
the finite fields, not a new order-eleven tangent theorem. Complete
the arrays by retaining the established finite-field bounds

    Tplus_11=Kfield_11=4e24*100^11.

Neither Vbar_11 nor a parameter derivative at order eleven is needed
or assigned zero. Polynomial padding carries no derivative claim.

## 2. Local scalar segment through the full required order

For the local metric expression, integrate along the scalar segment
(U_a,S_a) -> b(U_a,S_a), 0<=b<=1. This is an algebraic device for
the finite coefficient difference, not a change to the physical path,
state, source or domain. The same finite Tplus and Kfield bounds hold
on the segment. The original formulas therefore apply through eleven:

    M=3.003 exp(2L+4Tplus), DeltaM=4 M Kfield,
    C=(3/2) M M', DeltaC=(3/2)(DeltaM M'+M DeltaM'),
    R=Rref+2C'+2C^2, DeltaR=2DeltaC'+4C DeltaC,
    Ric=2R, DeltaRic=2DeltaR,
    Sc=M Ric, DeltaSc=DeltaM Ric+M DeltaRic,
    J=Ric+Sc M/4,
    DeltaJ=DeltaRic+(DeltaSc M+Sc DeltaM)/4,
    Weyl=M R+2M Ric+(2/3)Sc M^2,
    DeltaWeyl=DeltaM R+M DeltaR
      +2(DeltaM Ric+M DeltaRic)
      +(2/3)(DeltaSc M^2+2Sc M DeltaM).

L, Rref and the seed arrays pJ,pC,pZ are unchanged. The initial
weighted component arrays are

    JW=6.1 DeltaJ+13pJ,
    CW=6.1 DeltaWeyl+25pC,
    ZW=6.1 DeltaSc+pZ.

The metric/connection/curvature orders are eleven/ten/nine. Every
product has nonnegative coefficients; sharper valid inputs preserve
the inequalities. The old inputs must reproduce the old ledger
exactly before the new result is used.

Through order eight, take minima with the already proved integrated
first-path curvature arrays. Also cap ZW there by pZ+DZ from the
[radius/source result](finite-radius-operator-source.md), with

    DZ=complete_seed_polynomial(1+Zeta)+Lie+N.

Retain the degree-twelve seed polynomial and the complete selected
curvature remainder N through eight. That earlier theorem bounds
W^2(Z_a-Z0)/a, so W>=1 and a<=1 justify its use here. This does
not identify the full-action trace with the bookkeeping residual.

Convert the component arrays to raw actual covariant derivatives
using the unchanged ordered connection recurrence. Use ranks two,
four and zero for J,C4,Z; use the full unweighted Rm array and fixed
physical magnetic flux. Preserve all established smaller low bounds.
Then apply the original complete A0,A1,A2 and E formulas, with A_j
through six and E through seven. E7 still contains the Z ninth
derivative, supplied by the retained finite metric order eleven.
No inverse metric, magnetic, curvature or off-shell contact is removed.

## 3. Recompose the complete mixed row and known low inverse

Keep the physical H and Y lift matrices from the
[finite higher TT lift](finite-higher-tt-lift.md) unchanged. Feed the
new coefficients to the same seven-by-nine, five-channel ordered
Leibniz matrix from the preceding mixed-row proof. Its column sums
m_l give

    ||Bbar c||H_W6 <=g^-1 sum_(l=0)^8 m_l ||gradbar^l c||,
    Cmix=max_l m_l.

In particular div(E eta)=<grad E,eta>+E div eta has two separate
contraction costs, one and two. Keep the seventh derivative of E.
The weight remains outside differentiated outputs. Physical volume
and the fixed return are unchanged.

For c=Dbar^-1 f_c, retain only the proved low inverse factors:

    c0/sqrt(gamma0)<28g F0,
    c1/sqrt(gamma0)<265g F0,
    c2/sqrt(gamma0)<531sqrt(g) F0,
    F0=||f_c||L2/sqrt(gamma0)<6e-12 ell^2/g.

Since sqrt(g)>3e-17, define

    Klow_bar=28m0+265m1+[531/(3e-17)]m2,
    mhi=max_(3<=l<=8)m_l,
    Cstart=1e20+6e-12 Klow_bar.

The scalar coefficient remains below 1e20. The full source retains
the conditional bound

    N6(p)<=800/g [Cstart ell^2+mhi C_D,4 F4],
    F4=||f_c||H4/sqrt(gamma0).

Report exact rational constants, decimal enclosures, improvement
against the previous complete row, and both sufficient equal-split
budgets for R0 R2 R4 N6(p)<=1e57 epsilon_c at ell=1. F4, its
residual conversion, C_D,4 and higher positive transfers are still
unknown. An improved upper bound is not a target pass, a necessary
physical limit, a lower bound or a nonexistence theorem.

## 4. Verification and decision

Freeze note, checker and review before execution. Pin all inherited
inputs. Reproduce the old finite geometry exactly, recover the new
metric/connection at available orders, and verify componentwise
nonincrease through every valid order. Record the explicit validity
limits for each input. Retain the full updated local arrays and
all mixed channels, column sums, low inverse factors and budgets.
Identify the largest terms in columns zero and two after refinement.

An independent finite exponential expansion verifies the scalar-
segment identity through order eleven. Manufactured finite-log
fields at two precisions check the integrated weighted bound through
order ten. An order-eleven omission control must change Z9/E7;
the omitted calculation is never accepted as a bound. Earlier
independent tensor, curvature, volume-divergence and operator tests
remain available with their exact source hashes. These are ordinary
proof support, not formal certification or charged-kernel sampling.

Select the next step from the resulting complete source budget.
Nonlinear existence, compatibility, matching and stability remain
unproved regardless of whether a coefficient upper bound improves.

Remaining to-do list: useful complete mixed/inverse and positive-transfer
estimates; actual finite quantum/TT H4 source and residual conversion;
nonlinear mapping, neighborhood preservation, signed compatibility
and existence; physical validity, state/boundary, matching, stability,
charge/return, opening, transport, preparation, separate creation and audit.
