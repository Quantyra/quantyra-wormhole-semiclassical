# Free magnetic channels do not supply the SYK entropy flow

S071, 2026-10-12. **The direct free-channel transfer fails: its bilinear is
marginal and the specified contact determinant has no infrared constant
entropy change.** The interacting gauge theory is not settled by this
calculation. S071 is parked with that limitation; the distinct prepared-pair
route is developed in [entangled-pair-capture.md](entangled-pair-capture.md).

## Dimension and contact response

The massless lowest Landau channels of MMP are free two-dimensional Dirac
fields at leading order with the gauge background held fixed. Near an
AdS2 boundary, ds²=L²(-dt²+dz²)/z², their equation reduces to

    [gamma^z(z partial_z - 1/2) + z gamma^t partial_t] psi = 0.

At zero frequency psi=z^(1/2) psi_0. Equivalently, rescaling psi by
sqrt(L/z) gives a flat massless spinor. The boundary weight is Delta=1/2;
the left-right bilinear has dimension one in boundary time. Its coefficient
has zero engineering dimension. This is the near-throat statement; the
massless AdS4 field has boundary dimension 3/2 and its three-dimensional
bilinear is also marginal at this order. The marginality and the need to
check higher corrections are already noted in
[Bintanja et al., discussion](https://arxiv.org/html/2102.06628v2).

For an explicit free local junction, take incoming and outgoing two-channel
amplitudes related by

    U(h) = (1-i h sigma_x)(1+i h sigma_x)^(-1),
    |U_12|² = 4 h²/(1+h²)².

U is unitary for real h. Its eigenphases are -/+ 2 arctan h and are
independent of frequency. This exactly quadratic contact problem generates
no running energy scale. The h here defines this contact convention; no
unverified equality with the full four-dimensional coupling is assumed.
Finite propagation through an exterior adds a frequency-dependent delay
and possibly other channels. A regular low-frequency delay is not by
itself a relevant power-law coupling; zero modes, resonances and interacting
exterior response must be checked separately.

This conclusion does not forbid negative Casimir energy on a connected
strip. [Maldacena and Qi, appendix C.2](https://arxiv.org/html/1804.00491v3)
explicitly exhibit such energy for massless fermions. Energy supporting
an already-connected configuration is a different quantity from changing
the extremal entropy of two thermal disks.

## A regulated thermal determinant, with its normalization stated

Choose two identical real boundary modes with Green function
G(i omega)=i C sign(omega), no electric chemical potential, no extra
boundary zero mode, and antiperiodic thermal frequencies
omega_n=(2n+1)pi/beta. For a Hermitian inter-side bilinear, the determinant
ratio of the two-by-two inverse kernel is 1+lambda², lambda=C h. One pair
of real modes gives one log(1+lambda²) for each positive Matsubara mode.
A complex pair doubles this expression without changing the limit below.

Use the same smooth physical-frequency cutoff exp(-|omega|/Lambda) for
coupled and uncoupled theories, and set A=log(1+lambda²), a=pi/(beta Lambda).
The relative logarithm is

    log[Z_h/Z_0]_Lambda = A sum_(n>=0) exp[-(2n+1)a]
                       = A/(2 sinh a)
                       = A beta Lambda/(2pi) - A pi/(12 beta Lambda) + ... .

The first term is a local ground-energy counterterm, removed by the
relative subtraction. It is also killed by the entropy operator. Thus

    DeltaS_Lambda = (1-beta partial_beta) log[Z_h/Z_0]_Lambda
                 = A [1-a coth(a)]/(2 sinh a)
                 = -A pi/(6 beta Lambda) + O((beta Lambda)^(-3)).

There is no beta-independent entropy term: DeltaS tends to zero in this
specified free contact problem. The small negative finite-cutoff value is
a relative regulator effect, not negative absolute thermodynamic entropy.
No singular continuation from a relevant Delta<1/2 flow was taken.

An arbitrary multiplicative boundary partition-function normalization or
an added boundary state could supply an entropy constant. Those are extra
physical/UV data; this determinant does not determine them. Consequently
the result is not an absolute computation of the MMP gravitational S0.
It shows that the proposed entropy-removing flow is absent in this explicit
free completion of the light channels. The previous
[entropy/throat hierarchy argument](formation-entropy-transfer.md) still
applies conditionally to any bounded additive per-field entropy shift.

## Gauge field and applicability limit

In a common gauge theory the local-looking charged transfer operator needs
a Wilson line: psi_L^dagger W_(L,R) psi_R + h.c., where
W_(L,R)=exp(i integral_R^L A) transforms as
exp(i alpha_L) W exp(-i alpha_R). It preserves gauge invariance while
transferring charge. A separately imposed coupling of two holographic CFTs
instead specifies which boundary symmetry is retained. Neither permits
silently discarding the charge/flux dynamics.

[MMP appendix A](https://arxiv.org/html/1807.04726v3) gives the relevant
warning explicitly. For n=Nf q channels and the spherical Maxwell mode,
g_2²=g²/(4pi r_e²), so the Schwinger mass satisfies

    m_S² r_e² = n g²/(4pi²).

Even small four-dimensional g does not ensure n g² is small. The gauged
charge sector is massive; neutral modes remain. A gauge-variant single
fermion dimension cannot simply be substituted for the dimension of the
Wilson-line-dressed transfer operator. The latter's infrared response,
including the exterior, has not been derived here. We do not claim that
all interactions remain exactly marginal or that cooling can never work.

## Decision and checks

[Checker](../checks/formation_boundary_flow.py),
[receipt](../evidence/formation-boundary-flow.json), and
[source scopes](../evidence/formation-boundary-source-review.json).
The checker verifies unitarity, transmission and the kernel determinant,
derives the entropy limit symbolically, and compares twelve finite
Matsubara sums and entropy derivatives with the closed expression.

The free approximation supplies no controlled SYK-like route through the
extremal entropy barrier. The full interacting continuation is parked,
not declared impossible. It would require a gauge-invariant infrared
kernel and renormalized boundary entropy with the four-dimensional state
and exterior fixed. Rather than assume that missing enhancement, active
work tests opening a specified entangled pair. Initial-pair preparation
remains part of the broad objective.

Remaining to-do list: S071 reopening needs the dressed interacting kernel
and physical entropy normalization; active S072 must control finite-time
switching radiation, four-dimensional matching, initial preparation and
disturbances before a physical formation or stability claim.
