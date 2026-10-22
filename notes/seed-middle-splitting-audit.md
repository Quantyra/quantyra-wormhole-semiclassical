# Correction: fixed-cap jets and a moving-center UV subtraction

2026-10-21, E003/S167. Ordinary audit before supporting verification.

The complete middle C3 bound 1e22 asserted in
[the preceding note](seed-middle-coercive-jets.md) is not established by
its proof. Its fixed-cap gluing jets were combined with local UV jets
obtained by varying the center of the normal cap. An unbounded difference
between those two differentiations was omitted. Here "unbounded" means
not yet quantitatively bounded in this argument, not an infinite norm.

Withdraw the complete C3 improvement and its proposed derivative-based
metric consequences from current certified inputs. Keep the successful
arithmetic run and original proof immutable as historical evidence.
The fixed-cap row/gluing estimates themselves remain available. The
zeroth complete source bound also survives: |Theta_seed|<1e14 q^2 in
the stated middle region. This audit does not refute the numerical
inequality 1e22 for the actual source; it identifies its missing proof.

## 1. The two operations are different

The retained [local UV construction](stationary-reference-cold-bolt.md),
"Analytic normal charts and the complete local UV remainder", solves
the Goursat equations in (u0,p,pbar), caps v_N in THAT normal chart,
and takes the complete Born remainder. Its complex parameter estimate
and the [middle source adaptation](stationary-centered-canonical-response.md)
bound derivatives of that centered construction. The full remainder
includes all Born orders; it is not just a finite local heat polynomial.
Changing the cap away from the center can change this remainder.

In contrast, equations (8)-(11) of the preceding note explicitly fix
the normal chart, cutoffs and middle inverse while differentiating the
two observation rows. They bound the gluing function for one fixed cap.
The missing comparison is not supplied by a bound on the same zeroth
value, by locality of the heat counterterms, or by regularity alone.

Identify a small family of centered UV caps in common PHYSICAL
coordinates. Let x be the observation position, b the cap center, and
let U(x,b) be the renormalized UV-cap scalar trace, with the original
common finite restoration. Work on a common neighborhood where all
these caps agree with the seed. Let Theta(x) be the unchanged actual
trace, and define

    G(x,b)=Theta(x)-U(x,b),
    U_mov(x)=U(x,x), G_fix(x)=G(x,0),
    C(x)=U(x,x)-U(x,0).

At x=0 the caps are identical. The exact identity is

    Theta(x)=U_mov(x)+G_fix(x)-C(x).               (1)

The preceding combination used the first two terms and did not bound
C. The local finite counterterms are the same near the observation;
their b-variation there is zero. This does not make the nonlocal
renormalized UV-cap remainder independent of b.

For directional derivatives, with all terms evaluated at (0,0),

    C(0)=0,
    C^(j)(0)=sum_(r=1)^j binom(j,r)
                         partial_x^(j-r) partial_b^r U, j>=1.    (2)

For several observation/center coordinates, use the corresponding
symmetric derivative tensors and sum over ordered placements. The
Hilbert tensor norm is bounded by the same binomial sum of mixed norms.
The formula is invariant under choosing a common physical chart. It
does not identify a cap variation with a change of the physical state.

## 2. Resolvent terms that must be retained

Let S_b=(A_b+M)^-1 for the cap's fixed paired realization, expressed
in a common half-density space. The vertices A_b^(r) include all
coefficient, spin-connection and density terms of that identification.
They are not assumed to be scalar multiplication operators. On a
regular family with the required common graph domain, the resolvent
identity gives

    S'=-S A' S,
    S''=2S A' S A' S-S A'' S,
    S'''=-6S A' S A' S A' S
          +3S A'' S A' S+3S A' S A'' S-S A''' S.    (3)

The two order-three mixed products are distinct. These identities
follow by differentiating (A_b+M)S_b=I; see also
[Teschl, Mathematical Methods in Quantum Mechanics, Lemma 6.5](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf),
equation (6.5), printed page 159, PDF index 170. That identity supplies
no quantitative bound for the actual cap vertices.

In common physical coordinates the family agrees with the seed on
the observation neighborhood, so its genuine cap variation is remote.
Separated row estimates are relevant to (3), but each insertion and
all mixed observation derivatives must still be bounded. A proof may
instead bound the UV term differentiated with the SAME cap fixed.
Either method must retain the full angular sum and common subtraction.
The preceding c_j control observation rows; they are not bounds for
A', A'', A''' or for the mixed derivatives in (2).

## 3. A continuum control detecting the omitted derivative

This is a mechanism fixture, not a charged geometry or an evaluation
of the physical source. On L2(R), let L0=-d^2/dx^2+M^2, M>=2, and
v=1_[1,2]. The cap family L_b=L0+b|v><v| has the same local operator
near x=0 for every b. For |b|<1/10 it is coercive with lower bound
M^2-1/10. Its free Green kernel is exp(-M|x-y|)/(2M).

For x<1 define

    B=(exp(-M)-exp(-2M))/(2M^2),
    a=<v,L0^-1 v>=1/M^2-(1-exp(-M))/M^3>0.

Direct rank-one inversion gives

    U(x,b)=1/(2M)-b B^2 exp(2Mx)/(1+a b).

Take the actual operator to remain L0. Then Theta(x)=1/(2M) and
G_fix(x)=0, whereas

    d U(x,x)/dx at zero = -B^2 !=0.                (4)

All fixed-cap observation derivatives of G_fix vanish, and every
cap perturbation is separated from the observation. Nonetheless the
naive sum has a nonzero derivative while the actual trace is constant.
Equation (1) cancels it exactly. Higher derivative discrepancies are
the derivatives of -x B^2 exp(2Mx)/(1+a x), not zero by positivity.
This also shows why point agreement cannot transfer derivative bounds.

The issue does not require a nonlocal cap perturbation. For the local
family L_b=L0+b V with V(y)=(y-1)^4(2-y)^4 on [1,2] and zero elsewhere,
the common domain is H2 and V is C3, bounded and remote. The resolvent
identity gives

    partial_b U(0,0)
       =-integral_1^2 exp(-2My)V(y)dy/(4M^2)<0.     (5)

No finite-sampling inference is needed for its strict sign. Higher
cutoff smoothness can be chosen without changing that sign argument.

## 4. Retained bounds and a concrete repair budget

Let G_j be the rational fixed-cap bounds from the preceding frozen
ledger, normalized by j!100^j q^2. These remain component bounds in a
fixed normal chart through order nine. At order three the covariant
conversion gives G_3+G_1/200 as before. The local value bounds also
apply to the identical cap at x=b=0, hence

    |Theta(0)|/q^2 <=G_0+2e13+1e12 <1e14.          (6)

This pointwise statement does not differentiate an upper bound.

For a proposed repair define kappa_j as a uniform bound on the full
Cartesian derivative tensor |partial^j C(0)|/(j!100^j q^2), j=1,2,3.
At the normal center Gamma=0 and |partial Gamma|<300, so

    |nabla^3 C|/(3!100^3 q^2) <=kappa_3+kappa_1/200.

A sufficient repair condition for the PREVIOUS 1e22 target is

    kappa_1<=4e21, kappa_2<=4e21, kappa_3<=4e21.    (7)

Indeed the largest previously computed combined ledger is below
5.955e21, and 5.955e21+4e21+4e21/200<1e22. Each lower-order target
has more margin. Condition (7) is not yet verified for the actual caps.
The direct alternative is a complete same-fixed-cap UV derivative bound.
Both are meaningful next actions; no user input or external change is
needed merely because the estimates are unresolved.

The earlier 1e26 source-three ledger also combines a center-Cauchy
UV estimate with a spatial gluing estimate. This audit does not certify
that earlier combination as a fallback. Its cap-identification step
needs the same scope check. Preserve the earlier records and numerical
allowances, but do not infer a newly validated complete C3 input from
this correction. Descendant estimates that require that source input
remain conditional until the splitting comparison is supplied.

The source-dependent metric improvements 4.08e-7 and 4.488e-5 quoted
in the preceding note consequently remain CONDITIONAL on restoring
the complete C3 bound. No endpoint C5 contract had been changed, and
none is changed now. The actual cancellation, signed zero, nonlinear
existence and physical stability requirements remain unfulfilled.

## Verification contract

Freeze this audit, checker and review before execution. Verify the
fixed/moving split and mixed derivative identities through nine, the
noncommuting resolvent formulas through three, and the rank-one
continuum formulas. Independently integrate the rank-one coefficients
and the local-potential derivative at 60 and 90 digits. Reconcile the
retained point-value and conditional repair budgets with the exact
previous G_j values. The fixtures detect a scope error; they do not
disprove the actual source's 1e22 bound or evaluate its cap variation.

The preceding 101 checks remain successful arithmetic/mechanism tests;
they did not cover this decomposition compatibility. Do not overwrite
their result or claim that passing tests verified the whole proof.
No Lean invocation is selected and no physical gate is completed.

Remaining to-do list: bound the actual moving/fixed cap correction or
the UV jets with the same fixed cap; restore a complete source-three
input and audit its use in metric bounds; complete higher source and
transfer estimates, cancellation, signed compatibility and nonlinear
existence; physical validity, state/boundary, stability, matching,
charge/return, opening, transport, preparation, separate creation and audit.
