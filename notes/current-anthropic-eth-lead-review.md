# Anthropic ETH paper: computational lead and transfer limits

2026-10-26 HST, E003/S169. User-requested review of an emailed public
paper. Preserve the construction, closed audit, accepted results, original
norms and blocked research checkpoint. This is a literature assessment,
not execution of the missing signed-response certificate.

## Source and theorem scope

The supplied public post resolves to Claude (Anthropic), *A 2^{o(n)} Time
Algorithm for 3-SAT*, dated October 26, 2026, at repository commit
`16a5ce41304fcdaea4e0be7c132e1c361dc26da5`:
[paper](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/paper.pdf),
[README](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/README.md),
[formal statements](https://github.com/anthropics/formal-math/blob/16a5ce41304fcdaea4e0be7c132e1c361dc26da5/eth/MainResults.lean).
ETH here is the Exponential Time Hypothesis in complexity theory. It is
not the Eigenstate Thermalization Hypothesis in quantum statistical physics.

Theorem 1, pages 1-2, states a deterministic exact counting algorithm for
k-CNF formulas, with n occurring variables and m clauses, of time

    2^{O(min{k^2 n log(log n)/sqrt(log n), k m/sqrt(log m)})}.

Logs have base two; the O constant is absolute. For fixed k this is
subexponential. The stated theorem refutes ETH and SETH, not P != NP;
it does not establish a polynomial-time SAT algorithm. The paper explicitly
disclaims practical SAT performance. For arbitrary 3-CNF with m=5n, its
proved bound does not improve on exhaustive search unless n>10^10000;
for its bound in n alone the stated necessary threshold is n>10^(10^8).
These are limits of the proved bounds, not lower bounds for all algorithms
or a claim that performance becomes useful beyond the thresholds.

The repository reports Lean proofs on finite-alphabet multi-stack Turing
machines and checks with additional kernels. Appendix E, page 40, states
that correspondence of the formal definitions with the intended theorem
requires human reading; it also says no person was known to have checked
the development at writing. MainResults.lean explicitly leaves the standard
machine simulation and input-rewriting comparison outside its statements.
This review read the theorem, algorithm overview, Appendix E and formal
definitions; it did not independently audit the full proof or run Lean.
Reported kernel acceptance is not represented as our own verification.

## The useful idea, and the missing transfer

Sections 2.1-2.6 rewrite an exact finite tensor contraction as a SIGNED sum
of contractions. An eight-square identity rewrites crossing equality
constraints; terms can then use different row or column partitions.
This escapes a restriction to a single contraction tree. The recursion
tracks integer value bit lengths as well as operation counts. The README
also gives a condensed SWAP/singlet version of the rank idea.

This is a plausible method lead for evaluating a finite, commonly restored
representation while preserving cancellations. It is not a norm bound on
the continuum response. A transfer would first require a conforming finite
representation of the actual paired mass/log action and a certified error
in the original energy dual and physical output norms. The exact ultraviolet
logarithm, full angular/thermal tails, coincidence contacts, cutoff
transitions and regular-bolt returns must remain in that error accounting.
Finite rationalization and contraction errors also need bounds; exact
integer computation does not provide those automatically. The already
controlled P32 transport error does not cover a new discretization.

The current obstruction to progress is that representation/error bridge,
not an assumed ETH lower bound. Faster evaluation after discretization
does not justify the discretization or common-regulator limit. Further,
the accepted unfiltered operator has distance at least 1 from compact
operators in its original energy normalization. Finite tensor identities
do not overturn that result or produce a whole-action compact factorization.
This does not rule out compression of a retained finite representation,
the filtered F^3 map, or the actual affine correction.

## Decision and verification scope

Keep this as an algorithm-design lead behind the missing continuum estimate.
It is more specific than a general appeal to faster computation: its proposed
transfer is signed contraction with term-dependent partitions. No applicable
global error estimate, practical speedup, smaller physical-output constant,
or new available margin has been established. Current budget improvement
is zero; future applicability remains unknown. Do not change the action,
state, boundary realization or blocked status on this evidence.

The next decisive calculation is still the complete commonly restored
massive weak action on chi_N,r-alpha*z^4 and its correlated three-residual
transfer, with a certified global error. That calculation is NOT executed
by this review. The paper supplies neither negative energy nor a statement
about thermalization, physical matching or wormhole stability.

Public source bytes were retained outside Git and identified by SHA-256:

- paper.pdf: `7af518b09342c967ea371c74b8989806684c5d02aabe5763ca08554aff5025a0`
- README.md: `844e0dd78593f728e6a6a1f52b7b3d1e162020c1746ac1f4b3c4e5f2ff7f494e`
- MainResults.lean: `47f78a9550d95b37d8718dd1ee16ae824da893a5e17c319f737e6d0a5ed24944`

The large-exponent typography on page 2 was checked against a rendered
page, not inferred from flattened text. Private mail remains outside Git.
No independent theorem certification or scientific checker campaign was run.

Remaining to-do list: signed massive remainder and complete original output/source maps; every coupled entry and application error; nonlinear admission/contraction; exact compatibility; physical matching; semiclassical validity/noise; all-sector stability; opening, transport, preparation, separate creation and the full stable-wormhole objective.
