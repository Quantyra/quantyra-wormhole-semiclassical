# Fast scalar stability preflight in the exact AdS vacuum comparison

2026-10-20 HST, E003/S168 supporting S139/S074 and the unresolved
S077 stability requirement. Ordinary argument before frozen verification.
**The full retained coefficients have a below-BF scalar branch in the
exact conformal AdS vacuum comparison. Its radial operator is unbounded
below for every self-adjoint boundary extension. Transfer to the actual
charged exterior and its state remains open.**

The [preceding cold-source result](cold-projection-source-lower-bound.md)
makes a cancellation-sensitive signed coefficient estimate necessary for
the proposed cold source control. Before extending that static calculation,
this note tests whether its fast branch can support the required dynamics.
The earlier Euclidean factor alone is not used as an instability proof.

## 1. Scope and prior result

Use the fixed curvature/anomaly coefficients of
[the retained general-background operator](general-background-newton-operator.md).
The comparison is exact AdS4 with R_ab=-3lambda g_ab, lambda>0,
zero background electromagnetic field and an AdS-invariant conformal
vacuum. Its state is transported by conformal covariance when the metric
is varied conformally. Variations are first taken compactly supported
away from the conformal boundary. The Weyl-invariant part of the quantum
action retains its full field content; it is not approximated by an
angular truncation. Its vacuum stress is zero by maximal symmetry and
tracelessness, so it has no conformal response at this comparison.
The local anomaly and the fixed finite contacts remain.

These comparison assumptions do NOT replace the nonzero magnetic flux,
prepared state, finite-temperature exterior or boundary data of S167.
The charged solution has not been constructed. The comparison identifies
a potential obstruction that must be tested against that actual problem.

The general scalar AdS boundary-extension result is already known:
[Ishibashi--Wald, Theorem 3.2(iv)](https://arxiv.org/pdf/hep-th/0402184v4)
places nu^2<0 in the regime where every self-adjoint extension is
unbounded below. The retained calculation here identifies its parameter
and supplies a direct compact witness for this benchmark. It does not
claim a new general Breitenlohner--Freedman theorem.

[Flanagan--Wald, section IV](https://arxiv.org/pdf/gr-qc/9602052v2)
distinguishes higher-derivative runaway solutions from prescriptions
that select a perturbative effective branch. No order reduction or new
finite counterterm is adopted in this note. The result concerns the
FULL retained equations; it is not a statement that an ultraviolet
completion or a differently specified effective prescription is unstable.
Source versions, reading scope and downloaded-byte hashes are in the
[source review](../evidence/ads-fast-scalar-stability-source-review.json).

## 2. Lorentzian trace and an actual conformal metric lift

Use signature (-,+,+,+), Box=g^ab nabla_a nabla_b and

    epsilon=G/(60pi), lambda0=lambda+11epsilon lambda^2/6,
    Etr=R+12lambda0+epsilon[(3/2)C^2-(11/12)Euler4
                                             -20M_em+Box R].

This is the covariant trace identity in the retained scheme, with
the same sign of the local Box R term. Set

    J_ab=Ric_ab-(R/4)g_ab,
    Euler4=C^2-2J^2+R^2/6.

At the comparison C=J=M_em=0, R=-12lambda and Euler4=24lambda^2.
The displayed lambda0 makes Etr=0 exactly. For ANY metric variation,
psi=delta R obeys

    delta Euler4=-4lambda psi,
    delta(Box R)=Box psi,
    delta Etr=epsilon(Box+E)psi,
    E=60pi/G+11lambda/3.                              (1)

Thus the curvature branch has m_eff^2=-E in the convention
(Box-m_eff^2)psi=0. The sign in (1) is Lorentzian and follows
from the covariant trace, not from treating a Euclidean negative
factor as a time-dependent eigenvalue. Since R is constant at the
comparison, psi is invariant under linearized diffeomorphisms.

Here the trace branch can be lifted to the full tensor equation under
the comparison's stated conformal-vacuum assumption. Let h_ab=2s g_ab
and let calE_ab denote the full gravitational equation tensor, normalized
so its trace is -Etr and its classical part is G_ab-3lambda g_ab.
Its background value vanishes. Define

    Q_ab f=nabla_a nabla_b f-g_ab Box f+3lambda g_ab f.

On AdS, divergence Q=0 and tr Q=-3(Box-4lambda). The pure conformal
response of the Weyl-invariant action vanishes as stated in section 1;
the remaining tensor response is local of differential order at most
four. Maximal symmetry puts every such scalar-to-symmetric-tensor
operator in the form nabla_a nabla_b A(Box)s+g_ab B(Box)s,
with A of degree at most one and B of degree at most two. Commuting
covariant derivatives only adds powers of lambda. Conservation of
the linearized equation, since calE_ab=0 at the basepoint, requires
B=-(Box-3lambda)A. It is therefore Q_ab A(Box)s.

The conformal curvature variation is

    psi=-6(Box-4lambda)s.

Equation (1) fixes the polynomial A uniquely by its trace. A nonzero
polynomial cannot be annihilated identically by Box-4lambda. Hence

    delta calE_ab[2s g]=-2epsilon Q_ab(Box+E)s.         (2)

Its epsilon->0 classical part is -2Q_ab s, providing a sign check.
Every solution of (Box+E)s=0 solves all components of (2), including
the linear constraints, and has

    psi=6(E+4lambda)s !=0 for a nonzero s.              (3)

This is a constrained metric branch in the exact vacuum comparison,
not just an independently postulated test scalar. It is not a lift
theorem on the non-Einstein charged background. No independently
varied quantum state or general retarded shape response is removed
from that still-unresolved problem.

## 3. Global AdS radial operator and its BF parameter

Use dimensionless global coordinates

    ds^2=lambda^-1 cos(x)^-2[-d tau^2+dx^2+sin(x)^2 dOmega^2],
    0<x<pi/2, z=pi/2-x.

For a spherical curvature perturbation psi, put u=tan(x)psi.
The scalar equation in (1) becomes

    partial_tau^2 u=-A u,
    A=-partial_z^2+(2-E/lambda)/sin(z)^2
      =-partial_z^2-a/sin(z)^2,
    a=E/lambda-2=60pi/g+5/3, g=G lambda.               (4)

The radial Hilbert measure is dz. It follows either from the
Klein--Gordon static measure or directly from the displayed wave
equation. The center is at z=pi/2 and the conformal boundary at z=0.
The BF parameter and static exponents are

    nu^2=9/4+m_eff^2/lambda
        =-60pi/g-17/12=-kbar^2<0,
    psi ~ z^(3/2+i kbar), z^(3/2-i kbar).               (5)

This is precisely the spectral frequency of the preceding static
analysis. Equation (5) now has a Lorentzian interpretation because
of (1)-(4). The physical coefficient guard g<(4e-17)^2 and pi>3 give

    a>1.125e35>200.                                    (6)

No claim that these frequencies lie inside a controlled low-energy
effective approximation follows. Excluding them requires a specified
physical prescription; it cannot be inferred from the static range
projection, which projects equation sources rather than initial data.

## 4. A compact witness, independent of the outer boundary extension

Take the fixed real profile

    eta(v)=(v-1)^2(2-v)^2 for 1<v<2, zero otherwise.

It has zero value and first derivative at both endpoints and belongs
to H2 with compact support. Direct integration gives

    integral eta^2 dv=1/630,
    integral (eta')^2 dv=2/105,
    integral eta^2/v^2 dv >=(1/4) integral eta^2 dv.    (7)

For 0<d<pi/4 put u_d(z)=d^-1/2 eta(z/d). It is supported in
[d,2d], so it vanishes on a whole neighborhood of BOTH radial
endpoints. Using sin(z)<=z and (7),

    <u_d,A u_d>/||u_d||_2^2
      <=d^-2[12-a/4]<-a/(8d^2).                     (8)

The last strict inequality uses a>96, amply implied by (6).
For each fixed d the potential is smooth and bounded on a slightly
larger compact annulus. Approximation of eta by compact smooth
profiles in H2 therefore gives C_c^infinity test fields with

    Rayleigh quotient <-a/(16d^2)<-7e33/d^2.          (9)

The factor two is reserved for this approximation; it is not a
claim that the polynomial profile itself is infinitely smooth.
Let d tend to zero. Every self-adjoint extension of the minimal
radial differential operator agrees on these compact smooth tests,
so EVERY such extension is unbounded below. Imposing a boundary
condition cannot change their quadratic form. Each test already
preserves every outer source/metric germ at its initial time.

For any self-adjoint extension, spectral subspaces at arbitrarily
negative energies have arbitrarily large exponential wave growth.
For fixed tau>0, cosh(tau sqrt(-A)) on the negative spectrum is
unbounded. Thus this branch has no bounded wave evolution operator
on the usual scalar L2 data space, or on a fixed finite graph-norm
scale that supplies only polynomial spectral weights. In particular
there is no positive self-adjoint dynamics or positive-energy
stability estimate on that class. The branch lift (2)-(3) makes
this a vacuum-comparison metric obstruction under the stated state
assumptions; the scalar quadratic form is not identified with the
complete gravitational Hamiltonian.

Vanishing fast coefficients imposed on a stationary inverse is not
a proof that arbitrary compact initial disturbances avoid the branch.
If a proposed dynamical prescription excludes it, its initial-data
space and its evolution must be stated and checked separately.

## 5. What a transfer to the actual exterior would have to prove

There is a useful elementary robustness criterion. Suppose a physical
scalar reduction on the actual tail, in its correctly normalized L2
measure, has the self-adjoint local form

    A_tail=-partial_z[p(z)partial_z]-a/z^2+r(z),
    p(z)->1, z^2 r(z)->0 as z->0,                     (10)

and includes compact smooth test data of this sector. On sufficiently
small test annuli p<=2 and |z^2 r|<=1. The same profile gives

    Rayleigh(A_tail)<=d^-2[24-a/4+1]<-a/(8d^2),       (11)

since a>200. A further form term that is o(d^-2)||u_d||^2 on this
family cannot restore a lower bound. The corresponding smooth
approximation and unbounded-below conclusion then follow.

Equation (10) is an EXPLICIT TRANSFER HYPOTHESIS, not a theorem about
the actual charged exterior. The actual problem must establish the
constrained gauge-invariant variable, retarded state response,
boundary realization and permitted compact initial disturbances.
A nonlocal or frequency-dependent coupled reduction cannot simply
be inserted into (10). Background curvature/flux falloff alone does
not prove the required response or constraint estimate.

This preflight therefore completes the vacuum comparison, not the
physical instability verdict. It warrants checking the transfer or
a justified exclusion before investing further in a static solution
that may not meet the stability objective. S167's signed coefficient,
cold/higher inverse, nonlinear solve and compatibility remain open.
S077 still requires an accepted physical geometry and a specified
disturbance class; none of its acceptance tasks is completed here.

## Frozen verification contract

Freeze this note, source review and checker before execution. Reconcile
the current S167 receipt and retained action identities. Check vacuum
trace cancellation, Lorentzian mass/sign, conformal curvature/tensor
factorization, exact global radial conjugation, BF parameter, compact
profile integrals, benchmark reserves and perturbation criterion.
Use an independent finite-difference Rayleigh calculation for the
explicit polynomial profile at successively smaller annuli, compared
with its direct integral; this is a comparison-operator fixture and
does not simulate a wormhole or the retarded quantum response.
Retain all runs and exact failed sources. Continuum state/conservation,
domain, approximation and spectral arguments remain ordinary proofs.
No Lean campaign or independent specialist review is selected.

Remaining to-do list: actual charged exterior's constrained retarded
scalar problem and transfer or justified exclusion; S167 signed
coefficient, cold/higher response, nonlinear preservation and solve;
signed compatibility; exterior/state, charge/return, opening,
stability, transport, preparation, separate creation and full audit.
