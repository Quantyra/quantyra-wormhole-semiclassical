# Complete coupled transport error in original energy and outer observations

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: PROGRESS. Full stable-wormhole objective incomplete;
closed feasibility audit and original physical targets unchanged.

OUTCOME: the signed coupled estimate remains UNRESOLVED. This calculation
certifies its transport approximation error, including the perturbed
interior inverse and EVERY Schur product. It does not evaluate the retained
signed center or solve the complementary transition/interface equations.

## 1. First missing implication and why the truncation depth matters

The previous original-outer tail has input norms U3,W4; the global weak
tail is Linfinity -> L1. Neither by itself bounds a perturbation of the
full coupled energy operator or its interior inverse. Assigning the
conditional inverse/error constants one would be unjustified. We supply
the missing global shape-supremum comparison on the actual form domain.
It has a large finite-temperature cost, which is paid explicitly below.

Keep H, the original closed positive-diagonal space V, its canonical
return (nu,T), and the original energy E^2=(a+l)/(|W_action| beta),
|W_action|=N_halves pi B^2/(G gamma). The corresponding normalized dual
norm uses this SAME action factor. G=1 in the retained units.
Keep the canonical interior subspace W_I (zero on 0<z<.2) and the FULL
quotient V/W_I of current-canonical-interior-schur-bound.md. In particular
the annular equations on .1/.2 and their matching data are not replaced
by the original outer graph and three phase observations.

## 2. Global shape-supremum bound, including the regular bolt

Use the joint lower bound already retained in
current-geodesic-two-sided-remainder.md, equations (5)-(7):

    l>=.99 F_l E0+C S2/(100pi),
    S2=sum_(j=o,i) J2[chi_j T],
    J2[T]<=2S2+2(500)^4 J0,outer[T].                         (1)

These are the ACTUAL commonly restored forms, not only the local a2
contact. The positive fourth reserve is left after the full angular,
thermal, curved, reflected, off-shell and chart/contact costs. The
original York reconstruction and the same paired domain are used.
The current finite-seed transfer of (1) is explicit in that source; no
RN-to-neighborhood constant is imported. The original partition in (1)
is used only for this norm comparison, not to replace the output cutoff.

For stationary fields J2=N_halves(2pi/gamma) int |Delta T|^2 dz,
Delta T=(QT')', and similarly for J0. The prior energy-to-L2 result is

    ||T||2<=C_L E, C_L=sqrt(500/(.99 F_l)).

Divide (1) by its actual action/period factors and conservatively use
J0,outer<=J0. This gives

    ||Delta T||2<=C_D E,
    C_D^2=100pi B^2 beta/(G C)+2(500)^4 C_L^2.                (2)

No term is dropped from the fourth reserve or counted twice.

An elementary radial estimate avoids an unspecified Sobolev-chart
constant. Write L=z_h, v=1+h, t=1/z, gamma=hK(v). The exact finite-bolt
factorization of current-coincident-dirac-action-assembly.md is

    Q=(t-v)[2gamma v+(t-v)H1(t,v)]/t^3, H1>=.15t.

For z>=.5 it implies Q>=2gamma v^2 z^2(L-z)>=gamma(L-z)/2.
For z<=.5 the actual polynomial gives Q>=beta/4-|eps|-tau>.04,
which is larger than gamma(L-z)/2. The factored Q(.99)>0, Q(1)=-tau<0
and Q'<0 on [.5,1] retain the actual root .99<L<1. For the latter sign,
the leading derivative is -12beta z^2(1-z), the u term is nonpositive,
and the eps term is dominated by 4|eps|z^2(1-z); tau contributes
-3tau z^2. Hence on the WHOLE
physical interval

    Q(z)>=gamma(L-z)/2.                                    (3)

For an original smooth regular-bolt core field, QT'(L)=0. With d=Delta T,

    QT'(z)=-int_z^L d(t)dt,
    |T'(z)|<=2||d||2/[gamma sqrt(L-z)].                      (4)

There exists z0 with |T(z0)|<=||T||2/sqrt(L). Integrating (4) yields

    ||T||infinity <= ||T||2/sqrt(L)+4sqrt(L)||d||2/gamma
                  <=C_infinity E,
    C_infinity=C_L/sqrt(.99)+4C_D/gamma <7.477e50.            (5)

The large 1/gamma cost has not been renormalized away. Inequalities
(1)-(2) apply to differences of the original form-core approximants.
Equations (4)-(5) then give uniform convergence to a continuous T and
convergence of QT'=-int_z^L d. The distributional limit is the original
York return. Thus (5) holds on the actual energy completion with the
regular-bolt flux; no puncture condition or new scalar boundary is used.
It concerns T, not the physical bolt radius s+xi/2 or an entire C5 norm.

## 3. An entire-spectrum depth that controls the coupled energy operator

Retain the exact Dirac transport and the full commonly regulated low-word
part, all geometric-to-mass conversion terms, LLL/Weyl/local contacts and
boundary representatives. Only the trace-class transport tail is removed
to define tilde H. For EVERY angular/thermal label choose

    N(M,p)=256+ceil{[8log(M/18)+6log(1+|p|)]/log(1/.401)}.    (6)

There is no omitted physical mode. N depends only on seed spectral labels
and is held fixed in variations. Compared with the earlier strong rule
it increases every order by 128. It exceeds the earlier weak rule by
at least 128. Consequently set

    r_w=.3961^128, r_s=.401^128,
    e1=9.204e-79 r_w, e2=9.204e-77 r_w.                     (7)

The old trace-ideal/common-regulator proof supplies the ENTIRE normalized
source L1 bound e1 and Hessian Linfinity -> L1 bound e2. Higher retained
word sums are still to be evaluated; no independently renormalized
low-word determinant is asserted. The tail is the second derivative of
its absolutely convergent action series, so its Hessian is symmetric.

For delta H=H-tilde H and delta f=f-tilde f, (5) restores the original
energy/dual units:

    ||delta H||_(V->V*)<=epsilon=e2 C_infinity^2/beta<1.020e-25,
    ||delta f||V*<=e_f=e1 C_infinity/beta<1.365e-78.          (8)

The beta divisor is essential: normalization alone is not an improvement.
At the old base depth 128, this SAME valid upper-bound expression is
between 3.0865e26 and 3.0866e26. It cannot certify a small energy defect.
That is an excessive upper certificate, NOT a lower bound on the true
error. With this estimator, 199 is the first integer base depth giving
epsilon<.01; (6) gives the stated much smaller complete error.
Thus increasing the depth is justified by a quantitative sensitivity
test, not by optimizing an already useful error in its old norm.

The perturbation is bounded on the existing form space; no assertion
about a nonlinear action on arbitrary energy fields is needed. The
original mass/log core identity extends by (5) and (8). Regulator removal
occurs first in the established common construction. We do not interchange
it with an unsupported strong convergence of the raw heat cutoff.

## 4. Both exact and approximate interior inverses are controlled

The full fixed-seed operator obeys ||H^-1||<=1 and
||H||<=h_H=sqrt(1.338) in these energy units. Its actual restriction has
||P_I||<=K=7. Since epsilon<1 and K epsilon<1, the convergent Neumann
identities on these Hilbert spaces give

    ||tilde H^-1||<=1/(1-epsilon),
    ||tilde P_I||<=K_t=7/(1-7epsilon)<7.000001.              (9)

For example tilde H=H(I-H^-1 delta H); every factor is a bounded map
between the stated spaces. The same calculation is valid on the CLOSED
W_I, where the inverse was proved by its constrained signed form. This
does not infer an interior inverse from a compression of an indefinite
full inverse. Original interface matching remains in the quotient.

## 5. Every interior Schur-error term, in the original output rows

Let i be one of the four ORIGINAL observations. For an interior field w
in W_I, its local outer jets vanish, so U3(T_w)=||T_w||infinity and
W4(T_w)=0. Use the earlier outer-tail coefficients a_i and actual far
return coefficients c_i:

| i | a_i before multiplying r_s | c_i per energy |
| --- | ---: | ---: |
| Original H3 graph | 6.861e-45 | 4.583e-13 |
| Quartic finite part | 1.571e-51 | 7.949e-26 |
| Mass moment | 1.571e-52 | 7.949e-27 |
| Exact slow moment | 1.896e-51 | 1.141e-22 |

Write A_i=a_i r_s, d_i=A_i C_infinity. Then

    ||delta H_iI||<=d_i, ||tilde H_iI||<=c_i+d_i.             (10)

This is the actual differentiated kernel bound, with full original
.05/.1 transition, both paired paths and the exact slow-root subtraction.
Local contacts vanish only across these separated supports. They remain
unchanged in both complete centers. The bare mass-tail nu row is zero,
but neither original coupled gravitational equation is removed.

For B=H_iI, A=H_Iy and their tilde versions, the exact identity is

    BPA-tilde B tilde P tilde A
      =delta B P A+tilde B(P-tilde P)A+tilde B tilde P delta A,
    P-tilde P=-P delta H_II tilde P.                         (11)

Using (8)-(10) and ||H_Iy v||<=h_H E(v) gives the COMPLETE interior
product error, for EVERY actual column v,

    |delta(interior product)_i[v]|<=b_i E(v),
    b_i=d_i K h_H+(c_i+d_i)K_t epsilon(1+K h_H).             (12)

The last term includes both the perturbed input block and the perturbed
interior inverse. No direct-error-only approximation is used. Combining
(12) with the direct outer tail gives

    ||delta S_graph[v]||H3
      <1.094e-95 U3(T_v)+2.159e-104 W4(T_v)+2.977e-36 E(v),
    |delta S_4[v]|<2.505e-102 U3(T_v)+5.314e-49 E(v),
    |delta S_m[v]|<2.505e-103 U3(T_v)+5.314e-50 E(v),
    |delta S_slow[v]|<3.023e-102 U3(T_v)+7.410e-46 E(v).     (13)

S here is the actual Schur return BEFORE solving the remaining original
quotient/border system. All mass/quartic/correlated-slow/fast columns use
the same bounds. Their U3,W4,E values are NOT assigned one. For a bare
canonical scalar input T=0, the direct quantum-tail terms vanish, while
the interior mixing term in (13) remains.

For the affine right-hand side, ||f||V*<=h_H(1.1e-10)=F from the known
actual full linear response. The source has the preceding C4 kernel
bound without a differentiated mass vertex, so the conservative A_i
also bounds its direct output error. Substitution in (11) with f gives

    affine_error_i <= A_i+d_i K F
                         +(c_i+d_i)K_t(e_f+K epsilon F).     (14)

Executed full affine-error bounds in the graph, quartic, mass and exact
slow rows are respectively

    2.914e-46, 5.221e-59, 5.221e-60, 7.255e-56.              (15)

These are approximation errors, not values of the physical source or
signed solution. In particular (15) is not a measurement of P_rest in
its different projected units. Subsequent original phase-to-graph and
physical maps still need their actual norms.

## 6. Decision and unchanged feasibility requirements

The first previously missing numerical comparison for a complete
transport-error/inverse calculation is now supplied by (5)-(15).
The remaining unsupported step is the signed retained complementary
operator in the ORIGINAL strong quotient, including coincidence
restorations, the .1/.2 transition equations and exact phase/interface
matching. It has not been evaluated, nor has its strong inverse/border
amplifier. These are unknown estimates, not demonstrated violations.

The exact mass lift keeps its local cancellation in BOTH reference
equations, its original cutoff transition and fixed-temperature regular
bolt matching. It is not a global zero mode. The conversion
2133<||d_chi||H3<2174 and LINEAR |m|<1.627532e-15 still cost at most
3.538254568e-12, leaving near SHAPE >6.4617454e-12 after the known
distant return. Scalar quartic and correlated exact slow terms remain.

kappa_m<=5e-5 remains UNKNOWN. If certified, its original conditional
affine/feedback/energy debits remain 2.14765e-12, .0024731 and .0195239.
The residual ledger of the preceding interior result still requires

    Y_other+A_o 4.084e-22+Y_phase_return<=5.34235e-12,
    eta_other+.0195239 gamma_E+A_o 3.712e-12 J_E
                                      +eta_phase_return<.0012069.

When replacing exact entries by computed retained centers, add (13)-(15)
through the SAME proved remaining maps, without counting them twice.
For example the raw graph error for a mass lift must be included with
its actual U3,W4,E before comparison with 5e-5. No retained center is
known here, so no original mass or complete feedback target closes.

The conditional fast screen .135973<.1599164 retains its two-fast
admitted-family scope. Original nonlinear requirements remain

    E_center+1e-8 alpha_embed<=1e-11,
    R_X alpha_embed*.028913+Z0_other<=.25,
    R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

Equation (5) is a seed energy-to-shape-value comparison, not any of
these missing embeddings or a nonlinear-neighborhood theorem. Exact
energy E(h_lin)<1.1e-10 inserted in this deliberately coarse comparison
only certifies ||T_lin||infinity<8.225e40. It cannot certify the admitted
.001 value margin; using this comparison alone would require the much
stronger sufficient E<.001/(7.477e50). This is failure of that admission
upper certificate, not evidence that the actual correction is large.
The original physical-X comparisons or direct correction bounds remain
essential. Exact
compatibility, physical matching, validity/noise and all-sector stability
remain unresolved. No actual controlled wormhole or obstruction is proved.

RESEARCH DECISION: retain the construction and exact transport method,
using base depth 256 for the complete coupled error. The old depth's
generic energy certificate is excessive after its actual bolt conversion;
the stronger executed certificate closes that error requirement. Stop
refining it. The signed complementary center and original matching maps
now decide whether the sufficient global budgets are attainable.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): evaluate the complete
commonly restored signed complementary quotient/border for the affine
source and all coupled columns, including the independent transition
equations. Combine its actual input/output maps with (13)-(15), then
test the original mass and residual budgets simultaneously. This turn
executes a full coupled approximation-error certificate, not that center.

Remaining to-do list: signed complete coupled quotient/border and original
physical-X/neighborhood comparisons; nonlinear admission/contraction;
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
