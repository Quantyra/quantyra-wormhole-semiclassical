# Continuum and regulator leads: a hypothesis comparison

2026-10-26 HST, E003/S169. Requested literature and method assessment.
This is not a completed signed-response feasibility milestone. The full
stable-wormhole objective, accepted results, closed audit, original norms
and blocked research checkpoint are unchanged.

## The implication still needed

The target remains a computable enclosure in the original energy dual
E_t* of L_t[phi,.], phi=chi_N,r-alpha*z^4, and the complete correlated
three-residual transfer Zr. A defined weak form is not a certified evaluator
for that form. The selected-pair heat theorem and the accepted qualitative
C2 bridge on admitted stationary scalar families remain accepted within
their scopes. Neither supplies the required numerical E_t* error on the
entire test ball or establishes admission of the unknown correction.

Retain the actual first-order paired realization and its complementary
adjoint; do not declare that first-order operator self-adjoint. Both
positive partner squares, the induced mixed condition S=-K Pi_+/2, the
original projector, period and spin structure, and the regular disk bolt
must enter. The one-sided cubic mass and finite quartic logarithms cannot
be replaced by a smooth reflected metric. Cutoff transitions, common
angular/thermal regulation, coincidence contacts and finite-domain returns
remain part of the action being bounded.

## Ranked leads and their transfer limits

1. **Nonsmooth boundary resolvents: the strongest boundary-specific lead.**
   Abels, Grubb and Wood's
   [Extension Theory and Krein-type Resolvent Formulas for Nonsmooth
   Boundary Value Problems](https://arxiv.org/abs/1008.3281v2) treats
   second-order elliptic realizations through Poisson and boundary
   operators, including nonselfadjoint cases. Assumption 2.18 requires
   H_q^1 coefficients and B_(p,2)^(3/2) boundary regularity with
   1-n/q >= 1/2-(n-1)/p > 0, p>2, q>2. Theorem 6.5 expresses a resolvent
   as a reference resolvent plus a boundary correction. It requires the
   relevant resolvent sets and boundary inverse; those are not free gains.
   This offers a way to retain boundary returns explicitly while estimating
   the complete action. Applying it here still requires the actual matrix
   pair, induced mixed projector, regular bolt, parameter-dependent bounds
   and metric derivatives in the original energy topology. Its scalar
   boundary realization formulas are not already that extension. The paper
   does not directly treat our fourth-order response as a second-order
   operator. Apply such machinery, if justified, to the partner squares.

2. **Finite-regularity relative heat estimates: the complementary lead.**
   Hart Smith's [On the Trace of Schrodinger Heat Kernels and Regularity of
   Potentials](https://arxiv.org/abs/1809.05614) relates finite relative
   heat expansions to Sobolev regularity of a bounded compactly supported
   potential. Theorem 1.1 works on a fixed complete smooth manifold;
   in dimensions at least four it also assumes a localized trace bound.
   Its Duhamel/remainder method is relevant to avoiding unavailable smooth
   jets. It does not cover metric principal-coefficient variations,
   our paired boundary or uniform energy-dual response. This source was
   already identified in current-paired-regulator-domain-audit.md; it is
   not a newly discovered theorem closing that audit. The potentially
   useful adaptation is to combine finite-regularity remainder accounting
   with the actual boundary resolvent, preserving cancellations before
   estimating the differentiated trace.

3. **Certified primal/adjoint weak residuals: useful after source control.**
   Mallik's [equilibrated-moment goal estimator](https://arxiv.org/html/2107.06614v1)
   provides a guaranteed output-error framework for a clamped biharmonic
   problem on a polygonal domain. Its weak space is H_0^2, its load may
   lie in H^-2, and the displayed goal uses an L2 weight. This is closer
   in differential order to our weak response than a generic second-order
   finite-element estimator, but its boundary conditions and local PDE
   are different. The transferable idea is a joint primal/adjoint residual
   representation with certified reconstruction error, not its constants.
   Ishizaka's [2026 equilibrated-flux paper](https://arxiv.org/html/2602.01636v2)
   gives a useful residual/stability/Lipschitz verification pattern, but
   Remark 2.1 explicitly restricts its computable residual estimator to
   L2 data. Its abstract H^-1 formulation does not remove that restriction.
   Neither paper certifies the nonlocal renormalized quantum source.

4. **Filtered relative determinants: a numerical alternative with an
   analytic entry requirement.** Bornemann's
   [Fredholm determinant evaluation](https://arxiv.org/abs/0804.2543v2)
   supplies projection/Nystrom error analysis once the requisite operator
   and kernel regularity are established. The already accepted failure of
   compact approximation for the unfiltered response prevents using this
   directly on the full action. A genuinely trace-class relative remainder,
   or a suitable filtered representation, must first be proved, together
   with its differentiated error in our norm. Mass/log regularity does not
   license the analytic-kernel exponential rate. The author's
   [erratum](https://www-m3.ma.tum.de/bornemann/errata/FredDetErratum.pdf)
   also matters: Lemma 4.1 contains a quadratic remainder, not an exact
   global one-Lipschitz determinant estimate.

The [transfer-boundary heat-coefficient paper](https://arxiv.org/abs/hep-th/0208130)
by Gilkey, Kirsten and Vassilevich is useful for checking interface/contact
bookkeeping. Its smooth pieces, smooth smearing and Robin-type transfer
operator are not our mass/log mixed projector, and finitely many heat
coefficients are not the needed remainder bound. The previously reviewed
[Dunne-Kirsten radial determinant method](https://arxiv.org/abs/hep-th/0607066)
likewise needs a new curved paired-domain, whole-sum and differentiated
error argument. Neither should become another isolated coefficient study.

## One decisive proposed test

Use the actual structured phi above at the retained seed. Seek a complete
subtracted heat/resolvent representation in which the exact ultraviolet
logarithm and all finite contact restorations are explicitly retained.
Subtraction must be an identity with restoration, not a change of theory.
The first new implication to prove is a uniform energy-dual remainder
estimate for its metric vertices and boundary correction. Fixed-direction
continuity is insufficient. No seed-to-neighborhood transfer is assumed.

Here is a sufficient error specification, not an attained estimate. Use
tau for heat time, distinct from the family label t. Define H_phi(tau) to
be the FULLY angular/thermal-summed remainder of this proposed
representation as an E_t*-valued integrand, with all physical normalization
factors included, so that its contribution is int H_phi(tau) d tau/tau.
Keep the separately evaluated principal and restored terms with their
own certified errors. Prove, uniformly in admissible core approximations
and the common regulator,

    ||H_phi(tau)||_(E_t*) <= B tau^epsilon,
    epsilon>0, 0<tau<=tau_*,                     (1)

with explicit B, plus positive-time E_t* convergence and large-time
control. Then the omitted short-time contribution is at most

    E_UV <= B tau_*^epsilon/epsilon.             (2)

Equation (2) follows by elementary integration; (1) has NOT been proved
or evaluated for this complete structured action. It may be replaced by
another explicit integrable majorant. No angular weight may be counted
again in (2), since H_phi already includes the whole sum. An implementation
using a time-reduced radial integrand must restore that weight and all
original physical factors before comparing its error with (2).

The certificate must also account for the retained principal/contact
evaluation, positive-time quadrature, angular/thermal tail representation,
boundary and bolt returns, inverse applications and numerical arithmetic.
Only after all of these are enclosed can the correlated transfer and
original physical/source output maps be tested. A finite interval of
heat times or a finite collection of modes alone cannot settle (1).

The actual mechanism sought is cancellation of the nonintegrable
short-time pieces in the common, fully restored expression, leaving a
quantified integrable remainder on the whole E_t test ball. The old
physical-side cubic approximation has an H4 cost growing as the inverse
square root of the collar width. A bound depending on that uncontrolled
norm is unsuitable. The alternative determinant method faces the same
obstacle in trace-class norm and metric differentiability; faster finite
matrix contraction does not supply either property.

## Decision and feasibility consequence

Prioritize the nonsmooth paired-resolvent/finite-regularity weak-heat
approach for a bounded next attempt. Keep the construction and P32. The
reason is its direct treatment of the two missing ingredients, boundary
returns and regulator-uniform weak regularity. This is a method-selection
judgment, not a claim that the required constants are attainable.

Executed here: primary-source hypothesis comparison and identification of
the precise quantitative test (1). Not executed: its B or epsilon, a signed
center, a whole-action error enclosure or any original-norm improvement.
The contribution of this review to established admission margins is zero.
The accepted complete positive bolt estimator still exceeds its allowance
by approximately 1.533249982e25; this is an estimator failure, not a lower
bound on the correction. This comparison supplies no reason to change
construction and no evidence of wormhole existence.

The single proposed next calculation is the uniform E_t* remainder bound
for this actual structured input, including the original paired boundary
return and complete common-regulator sum. Its quantitative contribution is
the controlled continuum error in the signed certificate. If unavailable,
the blocked checkpoint remains; do not substitute more coefficient,
transport-degree or contraction-planner optimization.

## Provenance

Read the primary sources above and the current impasse, paired-regulator,
selected-smearing and curved-domain notes. The qualitative curved C2 result
is preserved; the missing numerical energy-dual estimate is distinguished
from it. No simulation, determinant evaluation or Lean run was performed.
Downloaded public PDFs are retained outside Git; SHA-256 identifiers:

- Abels-Grubb-Wood v2: `4d6de512a8e6a19676fb802f5734bc55ffa3cdb1165807e506e4cff9fff2c777`
- Smith author PDF: `a3afc11e46cfd382f23f1e3311302b454869628e1963460dec471de80de99af3`
- Transfer-boundary v1: `1b3135d3ba64f86a3d98ecb700bf604688260c3728b8c967e6fad0d5c7d534d9`

Remaining to-do list: complete signed massive response and original output/source maps; every coupled entry and application error; nonlinear admission/contraction; exact compatibility; physical matching; semiclassical validity/noise; all-sector stability; opening, transport, preparation, separate creation and the full stable-wormhole objective.
