# Gauging the complete mouth interaction

2026-10-11, S053. A return map at zero photon source does not specify the
vacuum interaction. Integrating out the exterior also generates a quadratic
photon vertex and an exterior determinant. This note constructs all these
terms at a gauge-covariant finite regulator, proves their Ward cancellation,
and bounds the first two bounded-potential insertions in a conditional
continuum absorbing-port problem. These are matching requirements using
established methods. They do not establish the interacting physical
wormhole, its renormalized stress, or the proposed damping result.

## 1. Exact regulated action, including the exterior

Let a finite fermion regulator and an orthogonal core/exterior partition
be fixed independently of the source. Write its C2 source-dependent matrix
as

    D[a] = [ A[a]  B[a] ],   R[a]=E[a]^-1,
           [ C[a]  E[a] ]    S[a]=A[a]-B[a] R[a] C[a].       (1)

Assume E and D are invertible in a neighborhood of zero source. Block
Gaussian elimination gives det D=det E det S. On continuous local logarithm
branches, the effective action and its Hessian are

    Gamma[a]=-log det D[a]=-log det E[a]-log det S[a],
    Pi_ab = Tr(R E_b R E_a - R E_ab)
          + Tr(S^-1 S_b S^-1 S_a - S^-1 S_ab).              (2)

Indices a,b mean directional derivatives in two specified classical photon
profiles at zero source, not spacetime or angular indices. No sign of Pi is
asserted: Euclidean source conventions can include factors of i. Equation
(2) is identical to Tr(D^-1 D_b D^-1 D_a-D^-1 D_ab).

To prove (2), differentiate log det X using Tr(X^-1 X_a), and use
(X^-1)_b=-X^-1 X_b X^-1. The vertices required by (1) are

    R_a = -R E_a R,
    R_ab = R E_a R E_b R + R E_b R E_a R - R E_ab R,
    S_a = A_a - B_a R C - B R_a C - B R C_a,
    S_ab = A_ab - B_ab R C - B_a R_b C - B_a R C_b
                    - B_b R_a C - B R_ab C - B R_a C_b
                    - B_b R C_a - B R_b C_a - B R C_ab.   (3)

Thus no derivatives of boundary hopping/identifications can be discarded
when the source acts there. In the useful special case of constant B,C and
bare linear source dependence in A,E,

    S_a = A_a + B R E_a R C,
    S_ab = -B (R E_a R E_b R + R E_b R E_a R) C.             (4)

The second term is an induced quadratic vertex, despite D_ab=0. Keeping
only the bubble with the effective linear vertex does not give (2).
It is also necessary to retain -log det E. Subtracting a reference removes
it only if its full source dependence and regulator match. Equality at
zero source is insufficient.

This is finite-dimensional algebra, not a determinant gluing theorem for
sharp continuum Dirac domains. Invertibility at a finite mass regulator
does not prove massless, infinite-volume or ultraviolet limits. A physical
calculation must control those limits and fix local renormalization terms.

## 2. Gauge cancellation and the information that is missing

Suppose the regulator obeys vector gauge covariance

    D[a+d chi]=G_chi D[a] G_chi^-1,

and the partition is by sites, so G=diag(G_core,G_ext). Then E transforms
by G_ext and S by G_core. Each determinant is invariant. Along a pure gauge
path a=t d chi, both terms of (2) vanish separately. Explicitly, for
X=E or S and Q its site gauge generator,

    X_chi=i[Q,X],  X_chichi=-[Q,[Q,X]],
    Tr(X^-1 X_chi X^-1 X_chi-X^-1 X_chichi)=0.              (5)

This follows by differentiating the exact similarity transformation; it
can also be checked by cyclicity of finite traces. The bubble alone need
not vanish. The statement concerns a vector gauge-covariant finite
regulator, not an unregulated anomalous chiral determinant.

There are two distinct information requirements:

* S[0], or a one-particle return at zero photon source, does not determine
  S_a or S_ab. Fixing D[0] while varying an exterior source profile changes
  (2). In a specified physical theory minimal coupling fixes the profile;
  it still has to be calculated rather than inferred from S[0].
* Even S[a] for all a does not determine det E[a]. Adding a decoupled
  exterior charged sector leaves S[a] unchanged but changes Gamma[a].
  This changes the model, so it is an information counterexample, not a
  proposed extra sector of the fixed wormhole theory.

For a normalized finite Gaussian photon covariance K, with zero one-point
term Gamma_a at the expansion point, the order-g2 vacuum correction is
(g2/2) Tr(K Pi). If a one-point term is present, completing the Gaussian
also gives -(g2/2) Gamma_a K_ab Gamma_b. Physical background stationarity
and counterterms therefore have to be specified. Neither formula licenses
contracting only the throat part of (2), nor combining the spherical
reference with a finite-return term on an incompatible exterior.

## 3. Continuum source insertions at positive resolution

Here the assumptions are explicit and independent of the finite regulator.
Let H be a magnetic Weyl/Dirac operator in an exterior with incoming and
outgoing flux traces. Assume a common graph domain, a well-posed absorbing
incoming boundary problem, and the Green current identity below. At
z=omega+i eta, eta>0, define R_z as the inverse of H-z with homogeneous
incoming data; P_z maps incoming port data to solutions of (H-z)u=0;
T_z=Gamma_out P_z. The trace norm includes the current's normal velocity.
For (H-z)u=f and incoming data h, assume

    2 eta ||u||^2 + ||Gamma_out u||^2
       = ||h||^2 - 2 Im <u,f>.                            (6)

The inner product is conjugate-linear in its first argument. The sign of
the last term will not matter for the bounds. Setting f=0 or h=0 gives

    ||R_z|| <= 1/eta,  ||P_z|| <= 1/sqrt(2 eta),
    ||Gamma_out R_z|| <= 1/sqrt(2 eta),  ||T_z|| <= 1.      (7)

For the third inequality maximize 2||u||||f||-2 eta||u||2
over ||u||. This avoids replacing a boundary estimate by a formal bulk
resolvent norm. No spectral gap or maximum return length is used.

Let H(t,s)=H+t V_a+s V_b, with bounded Hermitian multiplication potentials
on the same domain and source-independent traces. In particular a static
vector potential gives V=-sigma dot a. Differentiate the equation with
fixed incoming data. The first derivative solves a zero-incoming problem
with source -V_a P_z; the second solves one with the two ordered sources.
Consequently,

    T_a = -Gamma_out R_z V_a P_z,
    T_ab = Gamma_out R_z (V_a R_z V_b+V_b R_z V_a) P_z,
    ||T_a|| <= ||V_a||/(2 eta),
    ||T_ab|| <= ||V_a|| ||V_b||/eta^2.                    (8)

These derivatives exist in operator norm near zero by the bounded
perturbation resolvent identity and its Neumann series. More generally
an ordered n-insertion argument gives
||D^n T[V,...,V]|| <= n! ||V||^n/(2 eta^n) for n>=1.
Nonlinear H(t,s) adds -Gamma_out R_z H_ab P_z to T_ab.
Moving mouths or source-dependent trace identifications add further
boundary terms. None are assumed absent in the eventual physical problem.

Gauge covariance imposes an additional independent check. For a smooth
bounded chi with allowed graph-domain action,

    T[A+grad chi]=exp(i chi_out) T[A] exp(-i chi_in),
    T_chi=i chi_out T-i T chi_in,
    T_chichi=-chi_out^2 T+2 chi_out T chi_in-T chi_in^2.     (9)

Both derivatives vanish if chi vanishes at both mouths. This follows from
the continuum similarity transformation, and tests the endpoint phases
as well as the volume insertion. Chi can act as a port multiplication
operator; it need not be constant on a mouth.

## 4. Why these estimates alone do not close vacuum matching

At the working resolution eta=q^-1/16, (8) permits first and second
bounded-potential responses of orders q^1/16 and q^1/8. This is a finite
resolution control statement, not an interaction error estimate. The
photon vacuum diagram integrates frequencies down to zero and sums
spatial/angular photon profiles, whose relevant norms need not be
uniform. The bounds eta^-1 and eta^-2 are not integrable ultraviolet or
infrared estimates for that diagram. Squaring (8) cannot establish a
finite, renormalized polarization trace.

Likewise a prepared-packet estimate for T[0] does not bound the insertion
R_z V P_z: V can excite complementary angular states and the gapless
exterior continuum. A bound on T-T_model is not a bound on their source
derivatives without additional uniform source-dependent control. Existing
fixed-trim and shrinking-tail wave-packet residual estimates must be
extended to these inhomogeneous sources, including the induced quadratic
vertex and exterior determinant. Clock variations are generally
unbounded differential perturbations and are not covered by (8).

The next physical calculation is therefore the matched source-dependent
exterior Dirac response, with its trace/renormalization estimates and a
complete finite-return subtraction. The complete spherical unwrapped
reference remains a useful check, not the exterior of the actual two-mouth
geometry. Primary long-return damping and fallback finite-duration band
transfer remain candidates; this note establishes neither novelty nor
their interacting validity.

## 5. Verification and provenance

[The executable check](../checks/gauged_mouth_interaction.py) independently
compares (2)-(3) with full regulated Hessians and direct log-determinant
differences. A massive one-dimensional Wilson ring supplies gauge-covariant
matrices, source directions on core/exterior/crossing links, different
partitions, pure gauge cancellations and nonzero holonomy response. It is
an algebra/identity regulator, not a lattice simulation of the wormhole.
Linear exterior-only jets test (4) with no bare quadratic vertex. An
exterior spectator tests the determinant information requirement.

Continuum tests use a two-component Dirac interval with noncommuting
bounded piecewise constant sources. Block matrix exponentials compute
ordered transfer derivatives; direct source differences and an adaptive
variational ODE provide separate comparisons. Incoming/outgoing scattering
norms test (7)-(8). Gauge sources test both zero-mouth phases and nonzero
endpoint phases in (9). An exactly solvable chiral interval is a further
check of length and eta dependence. Numeric tests supplement the above
proofs; they do not prove the curved three-dimensional domain assumptions.

Results and source hashes are in [the receipt](../evidence/gauged-mouth-interaction.json).
Tolerances and cases are declared in the executable before calculation.
No Lean/Lake invocation is made; any later certification remains GCP-only.

The determinant factorization is established methodology. The primary
source reviewed here is Martin Luescher, *Schwarz-preconditioned HMC
algorithm for two-flavour lattice QCD*,
[hep-lat/0409106v1](https://arxiv.org/abs/hep-lat/0409106v1),
especially equations (2.1)-(2.3), (2.8)-(2.10), (3.1)-(3.8), and the
boundary projection explanation preceding and in Appendix B, (B.1)-(B.3).
It provides regulated block factorization and gauge-link differentiation,
not this wormhole's continuum determinant or a physical damping theorem.
Our normalized/un-normalized Schur conventions differ by det A; both
retain all determinant factors. This note's second derivatives and port
bounds are derived above rather than attributed to that paper.
[The source review receipt](../evidence/gauged-mouth-source-review.json)
records the version, cached hash, reviewed scope and narrow novelty verdict.

Remaining to-do list: match the full photon-source-dependent exterior and
finite-return subtraction, establish physical renormalization and static
and dynamic stress bounds, finish the interacting geometry/state and
primary/fallback novelty decision, then independently verify, GCP-certify
the consequential endpoint and prepare the reproducible manuscript.
