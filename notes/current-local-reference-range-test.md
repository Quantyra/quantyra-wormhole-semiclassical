# Full-range residual attempt: reject a positive local reference

2026-10-27 HST, E003/S169. First turn of the newly active goal to calculate
the full-domain bound on T_t-(I+T_t)F_n, including application errors.
The full stable-wormhole objective, closed audit, accepted results, P32,
physical norms and all unresolved requirements remain unchanged.

OUTCOME: the requested numerical residual certificate is UNRESOLVED.
This attempt tests one proposed way to avoid the missing quantum solve:
use the complete local coupled Hessian as a positive reference. That
method fails its coercivity hypothesis on the actual domain. The complete
geodesic reference avoids this particular failure, but its remaining
signed action and numerical application errors have not been evaluated.
This is not a completed feasibility result or a construction obstruction.

## 1. Actual spaces and the numerical implication still missing

Use E_t, norm squared ell_t/N, N=|W_action| beta, with the retained paired
first-order domain, complementary adjoint and both composed squares.
Keep the one-sided mass cubic and fixed finite logarithmic phase, original
cutoff transitions, clock/period/spin/scale and smooth regular-disk bolt.
York has J(0)=0 and (Q^2 J'/z^4)(L)=0. No smooth reflected surrogate,
inner wall, selected angular/thermal modes or nonlinear-ball extension
is introduced. The exact fixed-seed form defines

    T_t=ell_t^-1 B, B=b* a^-1 b, 0<=T_t<=q_t,
    q_t=.338/(1-epsilon), F_t=(I+T_t)^-1 T_t.

For a conforming finite trial space V_n in this domain, exact Ritz assembly
would give (ell_t+B)[F_n u,v]=B[u,v] for v in V_n. Its complete residual is

    r_n(u)[v]=B[u,v]/N-(ell_t+B)[F_n u,v]/N,
    ||R_n||=sup_{||u||E_t=1} ||r_n(u)||E_t*.

The first unexecuted numerical implication is an enclosure of the full
functionals (ell_t+B)[psi_j,.]/N, including the complementary test space,
and of their numerical source/solve applications. The available action
code still covers compact directions and two angular levels. It supplies
neither these columns nor their full-domain errors. No finite trial matrix,
application error or omitted-space constant is declared computed here.

## 2. Test the proposed local positive reference before using it

Define ell_loc to retain the exact local TT Hessian, all fixed local
restorations and the same York return, while omitting the positive
nonlocal geodesic remainder and the actual-minus-geodesic remainder.
The intended positive-reference solve would use ell_loc+B. This is a
test of a computational reference; ell_t in the target remains complete.

The exact local generator in current-finite-local-action-assembly.md is
in the fixed affine chart. Put phi0=-log z and d=Delta T. At highest
derivative order the exact York return has Delta nu=-d/3. The fourth-order
part of its normalized Weyl Hessian is obtained from

    phi L^2/60+L M/30-M^2/10,
    delta L=-2d, delta M=-d/3.

Its second variation is (2/15)(phi0+1/6)d^2. Multiplication by the
original negative action factor, and division by N, give

    -G C/(15 pi B^2 beta) (phi0+1/6)d^2.

The R2 fourth term cancels since -2d-6(-d/3)=0. Einstein, Maxwell,
LLL and the other contacts have lower derivative order on this TT
principal direction. The accepted common fourth restoration is exactly
G C kappa4/(15 pi B^2 beta) times d^2, where

    kappa4=H_q/2-EulerGamma-log(mu0 B), q=2000^4.

Consequently the complete local TT fourth coefficient is

    p4(z)=G C/(15 pi B^2 beta) a4(z),
    a4(z)=kappa4+log z-1/6.                         (1)

The 1/6 in (1) is the exact principal coefficient, not the 19/60 from a
later lower-bound estimate that also spends a fourth-order reserve.
The exact finite angular restoration uses the harmonic number H_q;
it is not a selection of physical angular modes. Keep the actual scale
enclosure |log(mu0 B)|<1e-8, rather than resetting that scale.

The interval calculation gives, on I=[1e-10,2e-10],

    a4(z)<-7.58, p4(z)<-5.65e-33.                  (2)

The second bound uses the original E_t normalization and inherited G=1
units. It is a coefficient of integral (Delta T)^2 dz, not a physical
output constant or an available admission margin.

## 3. Admissible full-domain test, including both equations

Here is an explicit ordinary test of positivity. Choose a nonzero smooth
bump chi supported inside I, for example the usual exponential bump,
and put m=(Q'/z^4)', rho=m chi^2. The actual nonconstant polynomial Q
makes integral m rho=integral m^2 chi^2 positive. Set

    c_k=integral m chi sin(kz) / integral m rho,
    J_k=(chi sin(kz)-c_k rho)/k,
    xi_k=(Q J_k)'-z^4 integral_0^z m(t)J_k(t)dt,
    nu_k=xi_k/2+Q J_k/z, T_k=-xi_k-(Q J_k)'/2.       (3)

The zero moment in (3) makes xi_k compactly supported along with J_k.
Differentiation verifies the exact York equation

    -(Q^2 J_k'/z^4)'=Q(4xi_k/z-xi_k')/z^4.

All fields vanish near the paired seam, original cutoffs and regular bolt;
their mass/log coefficients and endpoint returns are zero. They belong to
the actual core; this test does not replace the domain or truncate matter
modes. Integration by parts gives c_k=O(k^-m) for every fixed m.
On I, xi_k=Q chi cos(kz)+O(k^-1), T_k=-3xi_k/2+O(k^-1), and
nu_k+T_k/3=O(k^-1), with the corresponding differentiated asymptotics.
Thus integral (Delta T_k)^2 is a positive constant times k^4 to leading
order. Equation (2) makes ell_loc[xi_k,xi_k] tend to minus infinity.
The remaining terms have lower derivative order on this fixed compact
interval and cannot cancel that leading coefficient.

Both original equations give the exact Schur form ell_loc+B. The retained
mixed-row identity has order at most two, a has order four, and B=b*a^-1b
has local order zero. Its complete finite-domain returns are smooth on I.
Therefore its contribution to (3) is lower order as well. It cannot make
ell_loc+B positive. This disproves the positive/coercive reference
hypothesis before any numerical inverse or quantum error is assigned.

No finite change of the local fourth contact repairs this on the whole
domain, because log z is unbounded below. The actual nonlocal logarithm
is essential to the accepted positive full form. This argument does NOT
prove that every indefinite local inverse fails, or that the composed
residual on the structured B range is unbounded. Those stronger claims
would require separate arguments and are not made.

## 4. Quantify the complete-reference alternative without hiding errors

Retain instead the complete geodesic comparison ell0, including its full
logarithm, common restoration and actual domain. All actual angular,
thermal, coincidence and finite-return differences remain in the accepted
bound ||Q0||<=e0, Q0=ell_t^-1(ell0-ell_t), with

    e0=(.0022+epsilon)/(1-epsilon),
    F0=(ell0+B)^-1 B, f_t=q_t/(1+q_t).

This is an analytic reference operator, not an implemented numerical F_n.
Writing L=I+T_t gives (L+Q0)F0=T_t, hence EXACTLY

    R0=T_t-L F0=Q0 F0,
    ||R0||<=e0 f_t/(1-e0)<5.569803e-4.              (4)

No commutation is used: F0=(I+L^-1 Q0)^-1 F_t. The original domain
and coercivity used here are the accepted geodesic comparison, not the
rejected local reference. For a numerical approximation with a proved
a_app=||F_n-F0||, the complete bound would be

    ||R_n||<=e0 f_t/(1-e0)+(1+q_t)a_app.             (5)

The actual a_app is UNKNOWN. Equation (5) is not an attained numerical
certificate. In particular (4) cannot be relabelled as the bound for F_n.

The accepted original-output weights give the exact-reference error
estimator Gamma(e0/(1-e0)), with

    Gamma(z)=W0(3z+3z^2+z^3)+W1(2z+z^2).

Its upper expression is about 1.69436e14, or 1.01642e23 times the ENTIRE
1.667e-9 bolt allowance, even granting every other debit zero. That is
a quantified failure of this separated reference-only estimator; it is
not a lower bound on the actual error or a nonexistence theorem. No
conditional margin is treated as available. The already accepted internal
Neumann and P32 transport bounds are not refined or replaced.

## 5. Decision and remaining calculation

Reject the proposed positive local-only reference. Retain the complete
geodesic principal operator and the composed Ritz/direct-adjoint route.
The executed result is an actual-domain failure of that specific reference
and the quantitative comparison (4)-(5). It supplies NO implemented F_n,
signed G_t value, complete numerical residual or admission improvement.

The single next calculation remains a certified full quantum-action
column on an explicit admissible trial, as an E_t* functional over ALL
test fields, followed by the full range residual and application errors.
It must include the signed actual-minus-geodesic kernel and its matched
contacts; the .0022 envelope cannot substitute for those values. This
calculation is PROPOSED, not executed. The goal remains active/incomplete.

Remaining to-do list: full-domain residual with actual quantum/source/
solve/arithmetic errors; signed G_t; every original coupled admission and
mass/phase budget; nonlinear existence/contraction, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
transport, preparation, separate creation and the full objective.
