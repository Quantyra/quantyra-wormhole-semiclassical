# Actual affine regulator attempt: the quartic certificate does not transfer

2026-10-28 HST, E003/S169. Second turn of the NEW "execute the plan"
goal. Previous turn: limited progress, rejecting specified residual shortcuts.
This attempt supplies NO NEW complete quantum, residual or output bound.
The goal is active and incomplete; previous blocked runs are not counted.
Preserve the closed audit, P32, all accepted results and original norms.

## Attempt and first missing implication

The selected calculation remains the complete shape part of
r1=Gamma_rest''(g0)h1, on the actual paired mass/log realization. The
fixed-seed weak form defines it. Numerical h1 quadratures, the complete
signed action and the phase-aware strong regulator remainder are not
implemented. Both original rows, the exact UV logarithm, common contacts,
complete angular/thermal sums, cutoff transitions and regular-bolt returns
are required. None is replaced by the polynomial calculation below.

Inspected the ordinary arguments and implementation in
current-quartic-column-energy-dual.md,
current-quartic-resolvent-dual-test.md,
current-geometric-mass-conversion.md and
current-seed-separated-quantum-response.md. The accepted quartic result
uses the explicit direction -z^4 and polynomial coefficient enclosures.
The separated kernel evaluator does not evaluate the coincident action
or its integral against h1. Neither supplies the needed transfer to h1.

For the quartic heat argument, Delta=(Q partial_z)' and Hc=Delta+R.
Its sixth-order coefficient has leading part

    O6[v]=-(1/105) Hc Delta Hc v + terms of order at most four.

Here Q=beta+O(z^2) at the physical outer face. For an allowed cutoff
direction v=b_log z^4 log z near that face, exact differentiation gives

    (z^4 log z)''''=24 log z+50,
    (z^4 log z)'''''=24/z,
    (z^4 log z)''''''=-24/z^2,
    O6[v]=(24/105) beta^3 b_log z^-2+O(log^k z).       (1)

The cutoff equals one near the outer face and vanishes before the bolt;
its transition remains in the full expression. The displayed leading
term is unaffected by that transition. If b_log is nonzero, the separated
L2(dz) coefficient bound A6 used for the quartic proof is infinite.
The derivative identities were executed with exact rational arithmetic;
they are also obtained directly by six differentiations. They do not
show that the complete commonly restored quantum response diverges.
The heat expansion is not uniform in the unresolved boundary layer;
its joint contribution can cancel this separated coefficient.

This is a hypothesis test on an admitted log direction, NOT a counterexample
to the actual structured correction. For actual h1 the shape-log coefficient
has not been enclosed. The nonzero trace-log information alone is insufficient
to declare it nonzero, as the following exact reconstruction demonstrates.

## The actual coefficient that must be determined

Use the original reference constraints and fixed-disk return, with
F=lambda r^2+1+O(r^-1). If the actual asymptotic expansions contain

    f1=-2 delta_m/r+F_L r^-2 log(r/B)+...,
    a1=A_L r^-4 log(r/B)+...,

then xi'=a1+f1/F+xi F'/F and u1=(a1+xi')/2-xi/r give

    [xi]_(r^-3 log(r/B))=-(A_L+F_L/lambda)/5,
    [u1]_(r^-4 log(r/B))=A_L+F_L/(2lambda).

In z=B/r units, b_log=-(A_L+F_L/(2lambda))/B^4.
The accepted full trace identity says F_L+8lambda A_L=a_ref, so

    b_log=-(a_ref/(2lambda)-3A_L)/B^4.               (2)

In particular b_log=0 would require A_L=a_ref/(6lambda) and
F_L=-a_ref/3. The trace identity neither proves nor excludes this
correlation. The complete source/shape calculation must determine it.
No coefficient is selected to manufacture cancellation. Even a verified
zero would leave the mass seam, higher remainders and global returns.
There is no numerical smallness threshold that makes a nonzero b_log
acceptable in the separated A6 norm; the representation must change.

## Specified alternative theorem checked

[van den Berg, Gilkey and Kirsten, arXiv:1006.5644v1](https://arxiv.org/pdf/1006.5644v1),
Sections 1.1-1.9, Theorems 1.4 and 1.6, concerns a smooth compact
Laplace-type realization with mixed boundary conditions and weights
r^-alpha times a smooth factor, alpha<1. This is relevant to the smooth
physical side of our fixed seed. It is not a theorem on its rough double
or a rough metric neighborhood. Theorem 1.4 gives a trace expansion;
Theorem 1.6 computes its first boundary coefficients. Analytic dependence
of those coefficients does not provide a differentiated numerical remainder
for a logarithmic metric vertex. The theorem therefore does not supply
the complete determinant Hessian or our weighted H3 regulator estimate.
Differentiating in alpha to represent log weights remains an additional
uniform-remainder argument, not a certified bound obtained from this paper.

## Exact unresolved estimate and method decision

For the actual h1, form the complete regulated residual r1(tau) in the
original normalized rows, retaining source and response contacts together.
What is still missing is an explicit end polynomial and numerical A,alpha,m,
with alpha>0, such that the remainder after that same phase extraction has

    ||r1(0)-r1(tau)||Y_phase <= A tau^alpha(1+|log tau|)^m,

including the regular-atlas H3 observations, full finite returns and
actual numerical applications. Y_phase denotes the explicitly proposed
strong residual topology in current-affine-certificate-domain-test.md,
not an already identified original physical/source norm. The necessary
output comparisons to the original geometric observations remain unknown.
No source target or physical margin can be assigned from this theorem review.

Decision: do not transfer the polynomial A6 certificate to h1. Retain
the joint paired boundary-layer calculation with the actual log phases;
the alternative compensated inverse still faces its unknown complete
complementary estimate. No quantitative advantage or reason to change
construction is established. This is an unresolved research attempt,
not completion of the selected feasibility calculation. No further
scalar, degree, tail or checker optimization is selected.

One next calculation remains PROPOSED: construct the complete paired
boundary-layer remainder for the actual h1, with the above strong norm
and numerical constants, before evaluating its signed original outputs.
All original coupled budgets are unchanged. Exact compatibility, nonlinear
admission/contraction, physical matching, validity/noise and all-sector
stability remain unproved.

Remaining to-do list: complete quantum/strong-regulator and original-output
certificate; all coupled entries, nonlinear existence and compatibility;
physical matching/opening, semiclassical validity/noise, all-sector stability
and the full stable semiclassical wormhole objective.
