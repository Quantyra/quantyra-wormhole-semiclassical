# Noise-assisted entanglement: transfer test for S169

2026-10-26 HST. Requested literature lead, subordinate to the active
holographic comparison. No physical model, state, boundary condition or
closed audit is changed. The private mail remains outside Git.

The linked public post was resolved through its public oEmbed metadata
to Wang, Ding, Jing and Yu, arXiv:2610.31947v1, submitted 2026-10-22:
https://arxiv.org/html/2610.31947v1 . Read the model, noise assumptions,
Gaussian average, dark-state argument and limitations. This is not a
full independent audit of every concurrence formula or numerical figure.

The paper studies two resonant, identically coupled qubits and an ideal
lossless single cavity mode. A classical Gaussian process multiplies
one fixed interaction operator. It reports entanglement preservation
and revivals, with an exactly protected antisymmetric state. These are
quantum-information results; the work supplies no curved-spacetime
stress tensor, renormalized Dirac Hessian or gravitational stability bound.

## Executed implication 1: what the soluble noise mechanism actually gives

Write H_I(t)=(g0+eta(t))K, where
K=a_dagger(sigma_A^-+sigma_B^-)+a(sigma_A^++sigma_B^+).
All H_I(t) commute on their common excitation-sector domains. For
R_eta(s)=Gamma gamma exp(-gamma|s|)/2, direct double integration gives

    V(t)=Var integral_0^t eta(s) ds
        =Gamma[t-(1-exp(-gamma t))/gamma],
    V'(t)=Gamma(1-exp(-gamma t)) >= 0.                 (1)

In a K eigenbasis the exact averaged joint channel multiplies a coherence
with eigenvalue difference k by

    exp[-i k g0 t] exp[-k^2 V(t)/2].                  (2)

Thus its damping envelope NEVER increases for positive Gamma,gamma.
This directly conflicts with the statement in the paper's conclusion
that this OU decoherence factor recovers between oscillation cycles.
Reduced qubit concurrence can nevertheless revive because it is a
different quantity, obtained after coherent qubit-cavity evolution and
tracing the cavity. Equation (1) does not disprove those revivals.

Also R_eta(0)=Gamma gamma/2. Comparing smaller gamma at fixed Gamma
changes both correlation time and instantaneous variance; it is not a
test of memory alone at fixed noise variance. At fixed variance sigma^2,
V(t)=2sigma^2[t/gamma-(1-exp(-gamma t))/gamma^2], tending to sigma^2 t^2
as gamma tends to zero. This is a precisely stated comparator correction,
not an assertion that all colored noise is ineffective.

For psi_minus=(|01>-|10>)/sqrt(2), both collective raising and lowering
operators annihilate psi_minus. Hence K|psi_minus,n>=0 for EVERY n and
every realization, including zero noise. Protection here follows from
the common-coupling symmetry; noise does not create that symmetry.
There is no demonstrated analogous protected sector with the required
nonzero stress and full gravitational response in our construction.

The paper's rotating-wave model requires a controlled weak-coupling
regime. Its examples g0/omega=.5 and 5 are not evidence of such control;
the text itself identifies 5 as outside the standard regime. An
ultrastrong physical implementation requires retaining the omitted
counter-rotating interactions and testing the protection again.

## Executed implication 2: injecting this noise at our actual boundary

The existing coupled-fermion-domain.md derives, for a differentiable
boundary phase theta(t), an explicit boundary work term proportional to
theta_dot and a source-kinetic counterterm proportional to
integral theta_dot^2. Its finite coefficient requires physical matching.
The identical static angle can be conjugated away, but a time-dependent
conjugation adds the connection -i V_theta* partial_t V_theta. This is
not the paper's fixed scalar-multiple Hamiltonian. No domain-preserving
commuting-generator theorem has been proved for that driven problem.

Test the direct proposal theta(t)-theta0=eta(t), with the same OU
correlation as above. Its two-sided spectral density is

    S_theta(omega)=Gamma gamma^2/(gamma^2+omega^2).

With a PHYSICAL noise bandwidth Omega_c, distinct from the common matter
ultraviolet regulator, direct integration gives

    E[theta_dot^2]_(Omega_c)
       = integral_(-Omega_c)^Omega_c omega^2 S_theta(omega) d omega/(2pi)
       = Gamma gamma^2/pi [Omega_c-gamma atan(Omega_c/gamma)].  (3)

It grows linearly as Omega_c tends to infinity. Equivalently the
mean-square difference quotient is Gamma gamma(1-exp(-gamma dt))/dt^2,
which diverges as Gamma gamma^2/dt. Ideal OU phase noise therefore does
not satisfy the differentiable/finite-source-action hypotheses of the
existing boundary argument. White noise is no improvement. This is a
failure of a DIRECT transfer, not an obstruction to a smooth driven model
or to the wormhole construction. A finite bandwidth, full domain evolution,
renormalized boundary work and stress/noise response must all be included.
No finite source coefficient or useful bandwidth is invented here.

## Decision and budget consequence

Keep the paper as a possible later state-preparation/noise-control lead.
Its actionable idea is to test a genuine symmetry-protected sector or a
physically bandwidth-limited controller. Neither supplies the missing
static signed massive action or any original geometric-admission constant.
The contribution to the current coupled feasibility budget is ZERO.

Do not switch the static construction to ideal OU boundary noise. The
active holographic subtraction result and its excessive bolt estimator
remain as calculated. The single next selected calculation remains the
complete common-restored massive weak action on the actual Green-return
input, with its correlated physical-output error. A boundary-noise program
is a proposal for later, not an executed feasibility improvement.

Remaining to-do list: signed massive response and original output/source
maps; nonlinear admission and contraction; exact compatibility; physical
matching; semiclassical validity/noise; all-sector stability and the full
stable-wormhole objective.
