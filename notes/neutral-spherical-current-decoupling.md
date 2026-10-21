# Neutral spherical metric response decouples from electric current

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

For the selected free vectorlike Dirac field, fixed magnetic flux,
zero-source radial projector and neutral HH restriction, spherical
metric perturbations do not induce electric charge or radial current.
The mixed retarded metric/current kernels vanish, including the fixed
local contacts. Conversely an electric perturbation has no linear
metric-source response at this background. Thus delta A=0 and zero
electric charge form a consistent neutral spherical linear sector.

This removes the unevaluated electric response as a prerequisite for
testing the fast scalar in that sector. It neither proves stability
of the independent electric sector nor constructs a metric mode.
The physical metric boundary realization, charged curvature mixing
and retarded two-row metric solve remain unresolved.

## 1. State, domain and angular signs are essential

Use the complete angular reduction and canonical fields in
[the source reduction](dirac-spherical-source-reduction.md), with the
[actual zero-source domain](coupled-fermion-domain.md). On each open
exterior the optical Hamiltonian is

    h_0[V]=-i alpha partial_x+beta m(x,t) D_q,
    alpha=sigma_2, beta=sigma_1, m=exp(V)/B.              (1)

The radius Weyl field Phi contributes its already retained local
Weyl functional. It does not introduce another one-particle vertex
in (1). For real spherical electric potentials a_t,a_x, choose the
minimal-coupling convention

    h[V,a]=h_0[V]-a_t I-a_x alpha.                      (2)

Reversing the common convention for a reverses both vertices without
changing any zero-response statement below. At the zero-source outer
boundary, each exterior allows one beta eigenspace. In the common
folded notation the allowed projector is

    Pi_0=(1+tau_3 tensor beta tensor I)/2.              (3)

The spin and monopole bundle identification is the fixed one in the
domain note. No new projector is substituted. The original radial
projectors and h=0 specialization of the coupled condition are in
[Bintanja, Freivogel and Rolph, (2.8)-(2.12)](https://arxiv.org/html/2308.00871v2).
That source supplies the boundary input, not the symmetry proof here.

The state is each exterior's HH KMS restriction at the prepared beta_T,
with ZERO electric chemical potential. The fixed magnetic-flux sector
or a chemical potential selecting that flux is not an electric
occupation bias. At zero boundary source the exteriors are decoupled;
their thermofield correlations do not alter their individual retarded
responses, as established in the thermal response note. Initially
there is no independent charge-odd state perturbation or electric drive.

## 2. Constructing an antiunitary on the fixed magnetic Hilbert space

Let Gamma_S be sphere chirality. It is parallel, Hermitian, squares
to one, and anticommutes with the twisted sphere operator D_q.
For each positive eigenvalue lambda choose an orthonormal basis u_A
and choose Gamma_S u_A as the basis of eigenvalue -lambda. In this
pair D_q=lambda sigma_3 and Gamma_S=sigma_1 are real matrices. On
the zero-mode space Gamma_S has its definite sign and D_q=0.

Let K_q mean conjugation of coefficients in this chosen complete
angular basis and of radial coefficients. Then

    K_q D_q K_q=D_q, K_q Gamma_S K_q=Gamma_S.

This conjugation is defined on the fixed q Hilbert space. It is not
naive complex conjugation of a monopole section in one gauge patch,
which would change its bundle. It need not be a local four-dimensional
charge-conjugation operation. The restricted spherical Hamiltonians
and angularly summed orbit currents are all that are used here.

Define the antiunitary

    C=(beta tensor Gamma_S) K_q, C^2=1.                (4)

It does not exchange exteriors; put the flavor identity in (4) when
folding. Direct coefficient algebra gives

    C alpha C^-1=alpha, C(beta D_q)C^-1=-beta D_q,
    C(-i alpha partial_x)C^-1=+i alpha partial_x,
    C h[V,a] C^-1=-h[V,-a].                            (5)

The sign of i is important in the kinetic term. The actual projector
is real in this basis and commutes with beta Gamma_S, hence
C Pi_0 C^-1=Pi_0. Radial L2/horizon domains are preserved by this
bounded coefficient map. For compact metric histories no endpoint
domain is changed. The same identity holds on any admitted continuation
with the same projector and real scalar canonical rescalings.

The nonzero angular signs are both essential to (5); charge
conjugation within a single massive sign would generally change the
allowed boundary subspace. In the LLL the mass is zero and C reduces
to plus or minus beta K_q, so no missing negative-chirality zero modes
are invented. This uses a four-component vectorlike Dirac field.
A single four-dimensional Weyl field has a different radial LLL and
is outside this theorem.

## 3. The actual neutral thermal state preserves the symmetry

The canonical CAR algebra admits the complex-linear particle-hole
automorphism

    P(a(f))=a^dagger(C f), P(a^dagger(f))=a(C f).        (6)

Here a(f) is antilinear and a^dagger(f) linear in f, so the
antilinearity of C is exactly what (6) requires. At a finite regulator
the centered bilinear O_B=dGamma(B)-Tr(B)/2 transforms as

    P(O_B)=-O_(C B C^-1).                              (7)

Thus the many-body Hamiltonian is even, while the smeared electric
density (B=w I) and radial current (B=w alpha) are odd for real radial
w. Metric vertices obtained by differentiating h_0[V] are even by
(5),(7). The symmetry is an automorphism of observables, not an
antiunitary time-reversal operation on their expectation values.

The neutral thermal covariance is F_beta=(1+exp(beta_T h_0))^-1.
The exact functional-calculus relation

    C F_beta C^-1=1-F_beta                              (8)

shows that its quasifree state is invariant under (6). At a zero
eigenvalue F_beta=1/2; no zero-mode occupation is assigned by hand.
One can verify invariance of all observables by the CAR and Wick
rule, not just equality of the total charge. The same construction
defines the automorphism and state directly in the continuum CAR
algebra; no trace-class thermal density matrix on the horizon half-line
is required.

For any real spherical metric history at a=0, (5) holds at each time.
If S(t,t0) solves i partial_t S=h_0(t)S, antiunitarity and the minus
sign in (5) give C S C^-1=S. The evolved covariance therefore still
satisfies (8). This is evolution from the initial HH state, not
instantaneous rethermalization. Every smeared orbit current remains
zero. Rotational invariance then makes the four-dimensional current
zero pointwise in angle; an invariant tangential vector on the sphere
is also zero. This conclusion is independent of the metric being
on shell.

## 4. Retarded mixed response and the common subtraction

At the stationary background, if J is a current bilinear and O is a
metric vertex, then P([J(t),O(t')])=-[J(t),O(t')]. Invariance of the
state makes the expectation of this commutator zero. Equivalently,
the full regulated in-in generating functional for two metric and
electric source histories obeys

    W[g_+,g_-,a_+,a_-]=W[g_+,g_-,-a_+,-a_-].            (9)

Differentiate once in an electric source and arbitrarily in metric
sources at a_+=a_-=0. This proves vanishing of mixed responses and
their coincident-point contacts, rather than discarding such contacts.
For the bulk linearized source map in particular,

    delta_g j^t=delta_g j^x=0,
    delta_a T_ab=0, delta_a J_X=0.                     (10)

The reverse response in (10) follows either from (9) or from the
odd commutator of an even metric observable with an odd electric
vertex. It is not inferred by interchanging a retarded and advanced
kernel. Ordinary metric factors multiplying a current also have zero
variation because its unperturbed expectation is zero.

To remove the regulator, retain both signs in every included massive
angular level, all zero modes, and radial regulators invariant under C.
At each finite regulator (9) is an exact CAR identity. Sum the angular
contributions before the common four-dimensional subtraction. The
retained parity-even prescription and continuously fixed phase are
specified in the uniform-reference and paired-phase notes. No finite
coefficient is altered here. The available local gauge term and the
Weyl gauge contact depend on

    F_4^2=F_2^2+q^2/(2r^4),                            (11)

which is even under a->-a. Their metric/current cross derivative at
F_2=0 vanishes. Curvature-only terms have no current derivative.
Boundary subtractions inherited from the same invariant regulated
domain obey (9) as well. No independent parity-odd boundary phase,
theta F wedge F term or electric boundary source is added. Such a
change would require a new domain/state analysis, not an automatic
use of this theorem.

Consequently the common renormalized distributions retain (10).
No separate two-dimensional finite renormalization is used, and an
unknown current kernel is not set to zero by a large-mass estimate.
The current/current response is even and is generally nonzero;
nothing here computes or bounds it. The finite fixtures explicitly
retain a nonzero electric susceptibility.

## 5. Consequence for the actual coupled equations

The previous [Ward recovery](spherical-retarded-constraint-recovery.md)
derived the spherical electric equations, in its current convention,

    partial_x Q_e=U delta j^t,
    -partial_t Q_e=U delta j^x,
    Q_e=r^2 delta F_tx/(g_c^2 N^2), U=N^2 r^2.           (12)

For a neutral metric disturbance with delta A=0, (10) makes both
right sides zero. Zero initial electric charge and no electric
boundary drive select Q_e=0. This is a solution of the full linear
Maxwell/current subsystem. It does not require invertibility or
stability of the separate electric response operator. Conversely,
the Maxwell stress cross variation about a purely magnetic background
is zero, and (10) eliminates the linear Dirac stress/current mixing.
Thus the linear source/equation map separates into neutral spherical
metric and electric sectors under the stated preparation and domain.

There is no electric Schur-complement correction to the neutral metric
kernel. The two-row retarded metric operator already identified in the
preceding notes can be analyzed with delta A=0 consistently. This is
stronger than merely noticing the vanishing magnetic orbit Lorentz
force. It resolves the current question for this selected neutral
sector on the actual variable-radius exterior, without replacing that
exterior by the AdS normal limit.

The Ward theorem's gravitational incoming boundary datum is still not
fixed by the fermion projector. The seed is still off shell, so both
its residual contacts remain. No metric inverse, constrained fast mode,
nonlinear exterior or physical disturbance tolerance is obtained.
If a neutral instability is eventually proved, electric decoupling
will make it a consistent sector of the coupled problem; stability
of this sector alone would not establish stability of all disturbances.
Nonzero electric chemical potential, charge-odd initial data, changed
boundary projectors, nonspherical mixing or a different field theory
lie outside the result. No fast branch is removed or order-reduced.

## 6. Verification contract

The checker verifies the exact antiunitary coefficient identities,
the actual folded projector, the LLL and both massive signs, current
and metric vertex parities, and failure controls for a changed
projector or an omitted angular sign. It also constructs a finite
Galerkin Hamiltonian using functions that obey the reflecting beta
projector, then compares one-particle thermal/response formulas with
an independent full CAR Fock-space calculation. Thermal currents,
mixed retarded commutators and metric-quench currents vanish, while
electric susceptibility remains nonzero. A nonzero electric chemical
potential breaks the neutral-state identity. These finite fixtures
test the symmetry and its limits, not a physical exterior spectrum.
The continuum result is the ordinary argument above; no Lean run or
formal certification is claimed.

Remaining to-do list: physical metric boundary and two-row retarded
stability transfer; S167 signed coefficient, cold/higher response,
nonlinear solve and compatibility; exterior/state, charge/return,
opening, full-sector stability, transport, preparation, separate
creation and full audit.
