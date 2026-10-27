# ArcTN: practical contraction planning, no continuum certificate

2026-10-26 HST, E003/S169. Requested review of the latest emailed public
paper. The closed audit, accepted construction and original physical norms
remain unchanged. The signed-response research checkpoint stays blocked.

## Primary evidence and scope

Liu et al., [ArcTN: Efficient Tensor-Network Contraction Planning for Quantum
Circuit Simulation](https://arxiv.org/html/2610.31554v2), v2, October 23, 2026.
The paper concerns finite quantum-circuit networks. It supplies contraction
planning and slicing methods, with implementations and benchmark results.

Section 6.3 reports total planning/execution geometric means of 67.80 seconds
for Heavy and 146.34 for Cotengra AutoHQ: a 53.67% reduction on 134 runs
completed by all five compared configurations, across 29 instances. This
common-success subset is not a universal speed guarantee. Section 6.4's
slicing comparison measures plans and modeled costs without numerical
execution of those sliced contractions. Section 5's exact integer arithmetic
counts costs; Section 6.3 evaluates tensor values in complex floating point.

The [public interface](https://github.com/Quill-ArcLight/ArcTN/blob/3352ed53d5863ed6372d4176d90f94c872173431/README.en.md)
accepts indices, output labels and dimensions for path search, without tensor
values. It can save and reuse plans. Its Light/Heavy implementation is
proprietary and supplied under a
[noncommercial binary license](https://github.com/Quill-ArcLight/ArcTN/blob/3352ed53d5863ed6372d4176d90f94c872173431/LICENSE-ARCTN-ENGINE.txt);
the repository is not a complete openly licensed engine implementation.
No installation or license-dependent integration was attempted here.

## Consequence for the actual construction

This is a more concrete candidate for finite-evaluator acceleration than
the asymptotic ETH algorithm. That assessment follows from the reported
executed benchmarks and usable path interface, not a local benchmark of
our problem. The methods are different: the ETH paper rewrites a network
as a signed sum with term-dependent partitions; ArcTN searches contraction
orders and slicing for a supplied network. The path-only interface cannot
discover value-dependent quantum cancellations from dimensions alone.

A returned order may help execute an already justified finite representation
of our response. All slices would have to be accumulated, including their
roundoff. Finite slicing is not permission to omit angular/thermal tails.
The actual paired mass/log domain, common regulator, ultraviolet logarithm,
coincidence contacts, transitions and bolt returns still require their
continuum argument. No bound for that argument, the original energy-dual
error, or the original physical output is supplied by this paper.

An opaque planner would not by itself prevent a certificate: one could
independently validate its returned order and execute it with rigorously
bounded arithmetic. That does not supply the missing continuum-to-finite
comparison. Likewise, fewer operations can reduce a particular rounding
bound after tensor magnitudes and arithmetic are specified, but cannot be
translated into an admission improvement from these timing ratios alone.

Decision: retain ArcTN as an evaluator-engineering candidate. First obtain
the complete signed-action representation and its global certified error;
only then compare planners on that SAME expression, with the SAME precision,
hardware, error target and full planning-plus-execution cost. That benchmark
is proposed, not executed. No actual correction, source comparison, physical
margin or nonlinear requirement has been certified by this review.

The decisive research calculation remains the commonly restored massive
weak action on chi_N,r-alpha*z^4 and its correlated three-residual transfer.
This review does not perform it or change the construction. Current
established feasibility-budget improvement from ArcTN is zero.

## Provenance and completed work

Read the primary paper's method, cost model, execution and slicing sections,
and the pinned repository interface/license. No benchmark was independently
reproduced, no new scientific checker was run, and no Lean was invoked.
The reviewed email was marked read under the standing mailbox protocol;
private mail and downloaded public sources remain outside Git.

Repository snapshot: `3352ed53d5863ed6372d4176d90f94c872173431`.
Retained source SHA-256 values:

- Paper v2 HTML: `cda8633799cdfa0e814daffd7a5260d5fc5dbf884643a110c65b0497aa688752`
- README.en.md: `0e0ea44287be4417c441351367b3d06374f939b47852b074c42b218ad9720f95`
- Engine license: `8f96a6fa83a4ed41575a640bb0b39149db64f83e2b57c85f8d405c27b8cb1a76`

Remaining to-do list: complete signed massive response and original output/source maps; every coupled entry and application error; nonlinear admission/contraction; exact compatibility; physical matching; semiclassical validity/noise; all-sector stability; opening, transport, preparation, separate creation and the full stable-wormhole objective.
