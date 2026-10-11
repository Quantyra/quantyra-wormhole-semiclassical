# Photon response with a nonnegative exterior extension

2026-10-11, S053. The projected quadratic photon estimate survives a
radiating or reflecting exterior without a global photon gap, under the
explicit energy hypothesis below. The price is one more inverse power of
the imaginary frequency resolution. This removes the dual-Dirichlet
boundary requirement from this estimate. It does not derive the physical
coupled metric/photon load or the renormalized clock stress.

## 1. The boundary issue and the model

The [preceding calculation](propagating-photon-polarization.md) used a
fixed-metric Maxwell throat with reflecting dual-scalar ends. Its lowest
frequency is at least `2x`. An asymptotically flat radiating exterior has
continuous spectrum down to zero, even at fixed angular momentum. Its
full operator cannot inherit that gap. Replacing the ends while retaining
the gapped resolvent estimate would therefore be unjustified.

Keep the same q free Dirac channels, unequal circle lengths, angular
vertices, contact term and exact free susceptibility `Sigma(w)`. The
new statement changes only the photon continuation outside the throat.
Let `I=(-a/x,a/x)`, `0<a<pi/2`, `x>=x_min>0`, and let
`M=diag(l(l+1))`, `1<=l<q`, with angular multiplicities. Use canonically
normalized electric dual scalars `u_A`. Suppose the entire photon
extension has a densely defined closed nonnegative quadratic form h and
self-adjoint operator H on a Hilbert space containing the core fields,
with their ordinary L2 norm. Require

    h[u] >= sum_A int_I (|u_A'|^2
                   + l_A(l_A+1) x^2 sec^2(xs) |u_A|^2) ds.       (1)

The core restriction is defined on the form domain. The exterior may be
infinite, have arbitrarily low frequencies, mix angular channels, or
contain conservative auxiliary degrees of freedom. A nonnegative
exterior form joined by continuity, or by positive spring forms, satisfies
(1). Eliminating such an infinite exterior gives its retarded load at
`Im w>0`. No real-frequency outgoing limit is used in the proof.

This is a specified class of nonnegative extensions. The word passive
alone is insufficient: positivity of H without the local domination
(1), an arbitrary fitted frequency-dependent impedance, an active load,
or a non-self-adjoint generator does not establish these hypotheses.
The frequency-dependent boundary operator must come from the *same*
positive realization for all w. Arbitrary pointwise Robin data do not.

For an explicit physical test-field example, prescribe a static spherical
metric `ds^2=N(s)^2(-dt^2+ds^2)+R(s)^2 dOmega^2`, with N,R positive.
The Maxwell dual operator for either radiative parity is

    H_l=-partial_s^2+l(l+1) N(s)^2/R(s)^2.                       (2)

When `N/R=x sec(xs)` on the core, its standard full-line form satisfies
(1). This includes arbitrary positive spherical exterior profiles, not
only a hard reflecting wall. It is a probe Maxwell field on a prescribed
metric. Holding the metric fixed in a magnetically charged background
does not give the coupled Einstein-Maxwell perturbation problem.

Equation (2) is established separation of Maxwell's equations. Compare
[Cardoso et al., 1901.01265v2, Appendix A](https://arxiv.org/html/1901.01265v2#A1),
equations (43)-(46). In optical coordinates its electric-parity reduced
action has coefficient `alpha=R^2/N^2` for the two-dimensional field
strength and mass coefficient `kappa=l(l+1)` for the gauge-invariant
potential. Eliminating that potential gives the scalar action with
potential `kappa/alpha`; the magnetic parity gives the same scalar
operator directly. This derivation does not assume a Schwarzschild
metric, an Einstein equation or a photon gap.

## 2. Local energy gives the resolvent estimate

Set `m0^2=2 x_min^2`. Equation (1) implies

    ||u'||_core^2+m0^2 ||u||_core^2 <= h[u].                    (3)

Let `L(w)u=(u'|_I,-w u|_I)`, with potential-component signature
`J=diag(1,-1)`. Define `B(w)=L(conjugate(w))^* J` in the form-dual
sense. Thus B is the distribution `-partial_s f_0+w f_1`, including
the correct endpoint distributions for a core-supported source. It does
not drop surface terms or require the source to vanish at the mouths.

The weighted potential propagator of this reduced model is

    M^(1/2) D_ext(w) M^(1/2)
      = U^(-1) [1-L(w)(H-w^2)^(-1) B(w)] U^(-1),             (4)
    U=diag(1,-i).

This follows by eliminating the electric scalar, including its exterior
load. It also defines the reduced extension for general mixed scalar
loads. On a diagonal spherical extension it is the Maxwell potential
response, with the same Wick-rotation and contact conventions as the
preceding note. On a general nonspherical physical exterior, deriving
this representation remains a matching requirement, not an assumption
hidden in the word exterior. U is unitary; (4) is analytic in w and
uses a transpose-type analytic companion, not an adjoint at complex w.

For `w=v+i eta`, `eta>0`, define

    r(w)=sup_(lambda>=0) lambda/|lambda-w^2|
        = 1,                                  v^2<=eta^2,
        = (v^2+eta^2)/(2 |v| eta),             v^2>eta^2,
    c_B(w)=max(1,|w|^2/m0^2),
    C_ext(w)=1+c_B(w) r(w).                                  (5)

Then

    ||M^(1/2) D_ext(w) M^(1/2)|| <= C_ext(w).                 (6)

**Proof.** By (3), `L(w) H^(-1/2)` extends to a bounded operator of
norm at most `sqrt(c_B)`, on the closure of the range of `H^(1/2)`.
One can first use `H+e`, `e>0`, and pass to the form limit; any exact
kernel of H has zero core restriction and derivative by (3). The
same estimate holds for `H^(-1/2) B(w)`. The spectral theorem bounds
`H^(1/2)(H-w^2)^(-1)H^(1/2)` by r(w). Factoring the middle term in
(4) through these three operators proves (6). Differentiating
`lambda^2/[(lambda-v^2+eta^2)^2+4v^2 eta^2]` gives (5): when
`v^2>eta^2` its maximizer is `|w|^4/(v^2-eta^2)`; otherwise its
supremum is its limit one at infinity. In particular
`r(w)<=1+|w|/eta`. No smallest eigenvalue, upper bound on exterior
size, or diagonalization with M enters the argument.

For bounded |w| and fixed x_min, `C_ext=O(1+eta^-1)`. Even a resonance
at the clock frequency is included. This is a resolvent bound at positive
imaginary resolution, not the existence of a pointwise real-axis limit.

The extra resolution factor cannot be removed uniformly over this class.
For a continuum example attach a finite zero-potential interval of
length L to one core end, with ordinary continuity, a Neumann outer end
and a Neumann condition at the other core end. The full form satisfies
(1). A sine test supported on that interval has Rayleigh quotient
`pi^2/L^2`, so a sufficiently long exterior supplies an eigenvalue below
the reflecting throat gap. The lowest eigenvalue is strictly positive:
zero form energy would force a constant field vanishing in the core.
An eigenfunction cannot vanish throughout the core, since its value and
derivative at the join would then both vanish, forcing the exterior
solution to vanish by ODE uniqueness. Consequently L(w) sees this
eigenmode and (4) has a nonzero pole at its positive frequency. Its norm
grows as `1/eta` when approached from the upper half-plane.

A finite-dimensional check of this mechanism takes a one-coordinate
core with `m0^2=2` and one exterior coordinate, with form
`2|u|^2+|u-v|^2`. Its operator is
`[[3,-1],[-1,1]]`, with eigenvalues `2+-sqrt(2)`. The lower mode has
nonzero core weight `(2-sqrt(2))/4`. A nonzero spatial potential probe
therefore has a pole as w approaches `sqrt(2-sqrt(2))` from above:
its magnitude grows as `1/eta`. This positive extension satisfies (3).
It refutes only a uniform gap-based improvement for this class, not
physical wormhole stability.

## 3. Quadratic fermion resummation remains small

The preceding proof of the exact free current bound gives

    ||P^(-1/2) Sigma(w) P^(-1/2)|| <= B_J(w),
    P=diag(tau_l), tau_l=q beta_l/(4pi),
    beta_l=(q-1)! q!/[(q-1-l)! (q+l)!].

Let `s_q=max(tau_l/kappa_l)=q beta_1/(8pi)`. Define the dressed
propagator using its bounded resolvent identity, without requiring
invertibility of the compressed D_ext:

    D_g = D_ext + g^2 D_ext Sigma D_g,
    epsilon_ext = g^2 C_ext B_J s_q.

For `epsilon_ext<1`, a Neumann series on the M-weighted core spaces
gives

    ||M^(1/2)(D_g-D_ext)M^(1/2)||
       <= C_ext epsilon_ext/(1-epsilon_ext),                 (7)
    ||P^(-1/2)(Sigma+g^2 Sigma D_g Sigma-Sigma)P^(-1/2)||
       <= B_J epsilon_ext/(1-epsilon_ext).                   (8)

No angular commutation is assumed, either for the load or for Sigma.
At bounded real frequency and `eta` comparable to `delta`, both
`B_J` and `C_ext` are `O(delta^-1)`. For the working family
`g^2=q^-2`, `delta=q^-1/16`, this gives

    epsilon_ext=O(g^2 q/delta^2)=O(q^-7/8),
    bounds in (7),(8)=O(g^2 q/delta^3)=O(q^-13/16).           (9)

Thus the positive exterior load, including resonances and angular
mixing, does not by itself defeat this projected quadratic response
estimate in the stated joint limit. These are different norms from a
metric stress trace. Neither (7) nor (8) is a clock self-energy estimate.
For `g^2=q^-b`, `delta=q^-a`, sufficient smallness of epsilon is
`b>1+2a`; vanishing of the displayed absolute normalized corrections
uses `b>1+3a`. The earlier working window `b>1+5a` meets both.

## 4. Reproducible checks and remaining physical transfer

The [script](../checks/open_photon_load.py) verifies the scalar spectral
maximum and resonant counterexample independently, compares full matrix
inversion with spectral and Schur reductions, and tests angularly mixed
nonnegative loads. It also implements exact outgoing flat radial loads
on both ends of an AdS2 core. For each l the exterior potential is
`l(l+1)/(r_b+t)^2`, `t>=0`, with `r_b=cos(a)/x`. Its outgoing solution
is `r h_l^(1)(w r)` and its logarithmic derivative is

    zeta_l(w)=1/r_b+w h_l^(1)'(w r_b)/h_l^(1)(w r_b).        (10)

The outward core derivative equals zeta times the boundary value.
In a weak core resolvent the boundary contribution is therefore
`-zeta`, at both ends. At `w=i k`, it is real nonnegative. The l=1
identity `zeta=iw-i/[r_b(w r_b+i)]` checks the sign and normalization.
This removes an exterior truncation from the test. A finite-element
core calculation is independently compared with continuum ODE shooting;
its predeclared final probe tolerance is 0.1%. Finite fermion frequency
and core mesh refinements use the previous 1% probe tolerance. The
ordinary proof, not these finite samples, removes the cutoffs.

These radiative radial test loads are prescribed test-field geometries;
they are not the nonspherical two-mouth exterior of MMP. The abstract
estimate allows mixing only when its full positive realization and
local energy inequality have been established. In particular the
coupled gravito-electromagnetic system, rotational isometry gauge fields,
and moving mouths do not follow from fixed-metric Maxwell positivity.
[MMP v3, Appendix A.1](https://arxiv.org/html/1807.04726v3#A1.SS1)
already identifies a separate classically massless sphere-isometry
SU(2) sector. Assigning it the Maxwell multipole mass in (1) would be
incorrect. Rotation or an ergoregion also requires a separate energy
analysis. No full physical load is claimed here.

**Decision:** continue the physical long-return damping primary and
matched finite-duration energy-transfer fallback. Nonnegative photon
continuations are controlled at quadratic order with the revised
resolution loss. The estimate uses standard spectral/form and Maxwell
methods; it is a matching tool, not an accepted novel physics endpoint.
The decisive outstanding step is still the renormalized stress and
higher-current response with the actual coupled geometry and state.
Do not commission a separate formalization campaign for this bound or
substitute it for the consequential result. No Lean/Lake invocation or
cloud computation occurred.

Remaining to-do list: derive renormalized clock stress and higher-current
control, match the coupled metric/photon/fermion load and state, complete
the novelty decision and consequential result, then independently verify,
GCP-certify suitable endpoints and prepare the reproducible manuscript.
