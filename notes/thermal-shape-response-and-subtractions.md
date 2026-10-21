# Thermal shape response and frequency subtractions

2026-10-20, S168. Ordinary proof before focused verification.

The nonlocal shape row has an exact thermal spectral representation.
For a fixed compact radial probe, its third divided difference in the
squared imaginary frequency is nonnegative, with all retained local
terms included. Three positive-frequency values consequently give a
quadratic lower bound at larger imaginary frequency. This provides a
specific possible route to the missing causal estimate; the three
operator values and their needed uniform bounds are not evaluated here.

Neither thermal passivity nor a positive static Hessian alone supplies
that estimate. Local frequency terms survive renormalization and can
change its sign. This note derives the representation, the subtraction
requirement and an explicit counterexample to that shortcut. It does
not prove a causal inverse, a physical mode or an instability transfer.

## 1. Actual channel and state

Use the selected free vectorlike Dirac field, full magnetic angular
tower, original local coefficients, zero-source stationary paired domain
and HH restriction on each exterior. The [stationary state allocation](stationary-lll-quantum-exterior.md)
specifies reflecting exteriors at zero source and their thermofield
purification. Each exterior's local response uses its KMS restriction
at beta=2pi/kappa in the prepared Killing time U. Even operators on
the decoupled exteriors commute, so cross-purification correlations do
not create an off-diagonal retarded stress response at zero coupling.
No thermal ensemble for the sum of two forward-time Hamiltonians is
substituted for that purification. Driven coupling is outside this note.

On the open exterior choose the fixed optical chart and write

    g=e^(2Phi)[e^(2V)(-dU^2+dx^2)+B^2 dOmega^2],
    e^Phi=r/B, e^V=B N/r.                               (1)

Here N is the physical lapse; x runs to the horizon at minus infinity
and has a finite outer endpoint. Fix a real smooth radial test v with
compact support strictly inside this interval. V changes by a(U)v(x),
Phi is held fixed, and the boundary projector and initial state remain
fixed. All results below use the same radial probe at each frequency.
They are not estimates uniform as its support approaches an endpoint.

The [exact spherical reduction](dirac-spherical-source-reduction.md)
and the two-dimensional Weyl spinor rescaling give

    h_nu=-i sigma_2 partial_x+sigma_1 m_nu(x),
    m_nu=lambda_nu e^V/B,
    B_nu[v]=sigma_1 m_nu v,
    W_nu[v,w]=sigma_1 m_nu v w.                         (2)

lambda_nu runs over the q zero modes and BOTH signs of
sqrt(n(n+q)), with multiplicity q+2n for each sign. Keep the actual
domain in each channel; do not identify opposite-sign domains without
their intertwiner. The one-particle Hamiltonian is self-adjoint with
the original reflecting/paired prescription and horizon realization.
The matrices and the zero-source radial projector can be chosen real
in each angular block: -i sigma_2, sigma_1 and the projector are real.
Complex conjugation therefore supplies reciprocity for these real
shape probes. This is not an assumption of time reversal for arbitrary
electromagnetic perturbations of a magnetic four-dimensional state.

Equation (2) follows by differentiating m exp(a v); W is the second
vertex, not an optional contact. There is no time derivative in this
canonical Hamiltonian vertex. The zero modes have B=W=0 in this
mass representation. Their local two-dimensional anomaly and all
four-dimensional Weyl/measure terms must still be restored once.
The angular determinants are NOT independently renormalized.

At finite regulator, second-quantize (2). Write O_v=dGamma(B[v]),
W_vw=dGamma(W[v,w]), and use H(a)=H+a O_v+a^2 W_vv/2+O(a^3).
With the convention exp(-i omega U), define the causal shape response

    K_vw^R(omega)=<W_vw>
          -i integral_0^infinity e^(i omega t)<[O_v(t),O_w(0)]>dt
          +K_local,vw(omega), Im omega>0.                (3)

K_local contains the selected Weyl and measure restoration, common
counterterms, and the local gravitational and magnetic action if the
complete two-row Hessian is wanted. Its fixed-probe value at omega=i
sqrt(s) is a polynomial of degree at most two in s. All these actions
have at most four derivatives; integrations by parts and the static
background remove odd time derivatives in the reciprocal scalar form.
The definition includes the affine/exponential chart and density
contacts appropriate to the chosen gradient. None changes that degree.

Thus in (Phi,V) variables the nonlocal part of the uneliminated
two-row retarded Hessian occurs only in its V,V entry. The Phi,Phi
and mixed entries are local Weyl/gravitational variations, with the
[previous conformal column](charged-conformal-response-and-constraints.md).
Independent initial-state variations supply a separate source; their
Phi row has zero bulk trace. The spherical Ward constraints must also
be imposed with their initial and incoming data. This is not yet a
constraint-reduced operator. In particular the V,V entry is not the
stationary York TT block used by S167. A Schur complement with
frequency-dependent coefficients need not inherit the signs below.
Any electric/current response needed by the joint Maxwell constraints
is not evaluated by this metric-only row.

## 2. Lehmann measure and the zero-frequency distinction

First work with a finite spectral regulator solely to justify algebra.
For each channel let h|a>=e_a|a>, f_a=(1+exp(beta e_a))^-1 and
Delta_ab=e_b-e_a>0. Wick contraction and the causal commutator give

    D_v(s):=K_vv^R(i sqrt(s))
      =C_v(s)-sum_(nu,a<b) 2Delta_ab(f_a-f_b)
                                  |B_nu[v]_ab|^2/(Delta_ab^2+s),
    dmu_v(t)=sum_(nu,a<b) 2Delta_ab(f_a-f_b)|B_nu[v]_ab|^2
                                      delta(t-Delta_ab^2)dt >=0. (4)

The ordering a<b here means e_a<e_b, not an arbitrary index ordering.
Include multiplicities exactly once. C_v contains <W_vv> and the
local polynomial. In continuous spectrum, use the corresponding
two-energy spectral projectors of h; the sums in (4) become that
positive measure. The formula retains every angular channel and the
actual optical potential m_nu(x). It does not replace that potential
by its proper mass, which remains nonzero at the horizon while the
optical mass tends to zero. No optical gap is asserted.

For reciprocal real probes, polarization of mu_v gives a positive
matrix-valued measure: mu_vw uses Re(B[v]_ab conjugate(B[w]_ab)).
This supplies inequalities of quadratic forms, not just one diagonal
numerical sample. Causal evolution is analytic in Im omega>0 before
coupling to dynamical gravity; this matter spectral statement does
not establish analyticity of the inverse gravitational operator.

The static equilibrium Hessian is not automatically the zero-frequency
causal Hessian. At finite regulator, direct differentiation of
-beta^-1 Tr log(1+exp(-beta h(a))) gives

    D_eq,v(0)=<W_vv>-integral dmu_v(t)/t-Q_v+K_local,v(0),
    Q_v=beta sum_(e_a=e_b) f_a(1-f_a)|B[v]_ab|^2 >=0.     (5)

Degenerate blocks, including diagonal entries, belong in Q. The
retarded commutator has no equal-energy contribution. If its limit
exists, D_v(0+)=D_eq,v(0)+Q_v. This finite-regulator identity does
not assign a nonzero Q to the continuous horizon spectrum, nor does
it prove that either zero-frequency limit exists there. Positive
frequency subtractions below avoid those infrared limit questions.
For nonzero bosonic Matsubara frequencies omega_l=2pi l/beta, the
Euclidean connected covariance has the same denominator as (4);
the static extra term in (5) integrates to zero when l is nonzero.

## 3. Common ultraviolet subtraction and the positive remainder

The unrenormalized integral in (4) is divergent in four dimensions.
It cannot be used to infer D(s)>=D(0) after removing the cutoff:
C_v(s)-C_v(0) also contains frequency-dependent counterterms. Apply
the selected four-dimensional prescription after the angular sum.
There is no new adjustable polynomial in the model. Changing the
notation for its subtractions does not change its prescribed value.

For a fixed compact probe, the third subtracted integral is convergent.
Here is the local ultraviolet argument rather than a mass-gap claim.
The spatially smeared shape operator has four-dimensional engineering
dimension at most four: its vertex contains at most one angular/spatial
derivative of the Dirac field. The free Hadamard two-point singularity
and its differentiated products bound its separated Euclidean covariance
by C tau^-5(1+|log tau|)^p for small positive tau. Indeed, before
smearing the stress product has order distance^-8; integrating the
relative three spatial coordinates gives order tau^-5. The metric and
vertex coefficients are bounded on the fixed support. Image paths
have positive length there and do not change this coincident bound.
This uses the selected regular HH/Hadamard state, not a vacuum
substitution or a bound uniform on escaping supports.

More explicitly, let dW_v(Delta) be the nonnegative upward-transition
weight f_a(1-f_b)|B_ab|^2, including all channels. The connected
Euclidean covariance is positive and is at least
int_0^Omega exp(-Delta/Omega)dW_v(Delta) at tau=1/Omega.
It follows that W_v([0,Omega])<=C Omega^5(1+log Omega)^p.
Since f_a-f_b=(1-exp(-beta Delta))f_a(1-f_b), (4) gives

    mu_v([0,T])<=C' T^3(1+log T)^p,
    integral_0^infinity dmu_v(t)/(t+s_*)^4 < infinity,
                                      for every s_*>0.   (6)

Local finiteness at t=0 follows from the same positive Euclidean
covariance at a fixed interior thermal time. Integration by parts
against the cumulative bound proves convergence at infinity. No
numerical value of C or C', and no endpoint-uniform estimate, is
supplied. At finite spatial/energy regulator the same formulas hold
before removal; (6) permits the positive subtracted integral limit.

At s_*>0, the renormalized response therefore has the exact form

    D_v(s)=D_v(s_*)+(s-s_*)D_v'(s_*)
                +(s-s_*)^2 D_v''(s_*)/2
       +(s-s_*)^3 integral dmu_v(t)/[(t+s)(t+s_*)^3].     (7)

The first three coefficients are the values/derivatives of the actual
renormalized response, including its fixed local terms. They are not
free choices and are not determined by the spectral discontinuity.
Equation (7) follows by subtracting the first three terms of the
resolvent's geometric identity. Local terms of degree at most two
are preserved in the displayed quadratic and cancel from its remainder.
Uniform domination on s in a compact positive interval then gives

    D_v'''(s)=6 integral dmu_v(t)/(t+s)^4 >=0.             (8)

This is a sign for the third frequency derivative of the direct shape
quadratic form. It is NOT positivity of D_v, its first derivative,
its inverse, the gravitational energy, or the reduced physical scalar.

For distinct positive a<b<c define the quadratic interpolation of the
ACTUAL values D(a),D(b),D(c):

    I_2D(s)=sum_(u in {a,b,c}) D(u)
                            product_(w != u)(s-w)/(u-w).

The finite resolvent identity gives the stronger usable formula

    D_v(s)-I_2D_v(s)
     =(s-a)(s-b)(s-c) integral dmu_v(t)
                         /[(t+s)(t+a)(t+b)(t+c)].         (9)

Consequently D(s)>=I_2D(s) for s>c as forms on the common compact
core. The third divided difference at ANY four distinct positive
arguments equals the nonnegative integral of the four denominators.
Repeated arguments recover (8). Polarization preserves these form
statements; no simultaneous diagonalization of spatial operators is
needed. Nonzero Matsubara values, a=kappa^2 l_1^2, etc., are possible
anchors. Three scalar evaluations for one v do not certify the
operator inequality uniformly over all physical perturbations.

If a closed, appropriately realized interpolation form is coercive
on a specified physical shape space uniformly in the required s range,
(9) is a sufficient route to its direct-channel lower bound there.
Recovering an inverse also needs a closed realization, constraints,
couplings and surjectivity. None of these follows from (9) alone.
In particular no lower bound on I_2D for the retained exterior is
evaluated in this note. The previously established static TT gap
is not three positive-frequency values of this uneliminated shape row.

The spectral and subtraction framework is established work; compare
[Romatschke--Son, section II, (4)-(17)](https://arxiv.org/pdf/0903.3946v1).
Their stress-channel discussion explains why ultraviolet subtractions
and response contacts must accompany dispersion relations. No QGP
sum rule or flat-space estimate is imported into this exterior. The
retained optical vertex and frequency-interpolation application are
the present contribution; the underlying Lehmann principle is not new.

## 4. A static gap and passivity do not fix the frequency bound

Here is an exact diagnostic with positive measure delta(t-2)dt:

    D_+(s)=1+s+s^2+1/2-1/(2+s),
    D_-(s)=1-10s+s^2+1/2-1/(2+s).                       (10)

Both have D(0)=1, the same nonnegative spectral measure, the same
positive leading quadratic coefficient and D'''(s)=6/(2+s)^4>0.
Yet D_+(1)=19/6>0 while D_-(1)=-47/6<0. The difference is the local
linear term. This refutes inference of a positive imaginary-frequency
bound solely from a static positive value and spectral passivity.
It does not exhibit a second allowed choice of counterterms in the
retained model: those coefficients are fixed and must be evaluated.
Nor is this finite diagnostic a spectrum of the charged exterior.

The next required calculation is thus concrete. Bound the actual
positive-frequency shape forms (or the three coefficients in (7))
in the same gradient and physical domain, retain the local conformal
and mixed rows, and impose the initial Ward constraints. Establish
uniform control on the escaping annuli before applying the vacuum
scalar-tail criterion. This turn supplies the nonlocal spectral row
and a contact-preserving sufficient bound, but not that control.
S168.3-S168.4 and the signed/nonlinear tasks of S167 remain open.

## Verification contract

Freeze note, checker and source review before execution. Reconcile
the preceding conformal/constraint receipt and its transitive hashes.
Check the optical mass vertices, KMS occupation identity, seagull and
static degenerate-block terms, resolvent subtraction/interpolation
identities, signs and counterexample exactly.

Use a fixed three-level free-fermion fixture, beta=13/10 and energies
(-2,1/2,3), with the real matrices specified in the checker. At 70
digits compare the spectral causal formula against direct matrix
commutator Laplace integration at gamma=.4,1.1,2.7. Integrate through
time 256, with an explicit exponential tail bound; require scaled
error below 1e-30. Independently integrate the connected Euclidean
matrix covariance for Matsubara l=0,1,2,3 and require error below
1e-50; include Q only for l=0. Check the equilibrium Hessian using
eigenvalues of h+aB+a^2 W/2 at steps 1e-4,5e-5, requiring scaled
error below 1e-6 and error refinement factor below .35. Check (9)
against independent quadrature of a continuous positive measure and
four positive anchors; detect wrong bubble sign and omitted contact.

These are exact identities and independent finite fixtures, not
numerical exterior spectral data. The Hadamard bound, spectral
limits and core form inequalities remain ordinary arguments. No
Lean campaign is selected. Retain every run and failed exact source;
do not edit scientific sources after a successful frozen run.

Remaining to-do list: actual positive-frequency shape estimates,
constraint/domain control and instability transfer or justified
exclusion; S167 signed coefficient, nonlinear solve and compatibility;
exterior/state, charge/return, opening, stability, transport,
preparation, separate creation and full-result audit.
