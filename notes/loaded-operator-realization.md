# Realizing the higher-mode load and smoothing its gauge seam

2026-10-11, S053. The prescribed spherical higher-mode load can be attached
to the flat dipole as a closed operator with a resolvent and genuine L2
interface traces. Its existence is no longer an assumption for the collar
class specified below. Keeping the throat in the Hilbert space is essential:
its arbitrarily high angular modes do not have a uniform scattering loss.
An explicit smooth gauge completion also preserves the asymptotic matching
estimate after a compatible choice of the auxiliary port radius.

This constructs a prescribed free-field background, not an Einstein-Maxwell
solution. The smooth gauge interpolation generally requires a collar current;
its physical admissibility, metric backreaction, state and renormalized
stress have not been established. There is no stability or novelty acceptance.

## 1. Geometry, Hilbert space and the exact boundary condition

Use d=1 units. Let E be the flat exterior of the two radius-p balls with
the exact point-dipole magnetic connection. Join their boundaries by the
prescribed spherical throat of [the load note](throat-load-feedback.md),

    ds^2=N(x)^2(-dt^2+dx^2)+R(x)^2 dOmega^2,  0<=x<=L.

In this realization require Euclidean end collars:

    N=1, R=p-x                 for 0<=x<=p/4,
    N=1, R=p-(L-x)             for L-p/4<=x<=L.              (1)

Take L>p/2 and smooth positive N,R in between. Positive interpolants exist;
no field equation is imposed. These conditions ensure a smooth principal
Dirac symbol across each seam and give the previous collar constants
ell_c=p/4, beta=2/p. The theorem is not asserted for every unspecified
throat profile with a collar bound alone. A physical near-RN throat must
be compared to this class, or treated by its own transmission analysis.

Identify the magnetic bundles and spin frames at the seams with the
flux-normalized conventions of the mouth calculation. The charges at
the two exterior boundaries have the opposite orientations required by
one flux-q throat. A fixed unit phase at a seam is allowed; it specifies
a holonomy convention and is kept fixed. It does not follow from the
existence proof. Gauge strings are patch choices, not extra boundaries.

On the throat, let G be the kernel of the unit-sphere magnetic Dirac
operator and P_H=1-P_G. G has dimension q per Weyl chirality. It is smooth,
finite-dimensional, invariant under the radial current matrix, and its
projection commutes with the separated throat operator. Work on

    Hcal = L2(E) direct_sum P_H L2([0,L] x S2).               (2)

The throat measure is the optical, flux-normalized measure after the
known spinor rescaling. Equivalent smooth positive measures may be used
in the local regularity argument. All higher angular modes, signs and
degeneracies are present; no angular cutoff is taken.

Let H_E be the flat magnetic Weyl expression and H_T the spherical one.
Write K=i(H_E direct_sum H_T). Start with compactly supported, piecewise
smooth sections satisfying, at each seam,

    P_H trace_E u = U trace_T v,
    P_G P_in,E trace_E u = 0.                                (3)

Here U comes from the geometric bundle identification and smooth flux-density
rescaling, including orientation; an arbitrary nonlocal high-mode unitary
is not included in this claim. The incoming
and outgoing ground spaces mean the actual current eigenspaces at the
two mouths; a single Weyl chirality need not have ground inputs at both
mouths. Vectorlike Dirac fields are direct sums of the chirality results.
Let K_out be the graph closure of this core. The other ground traces are
free. High-mode transmission is imposed on the full spinor, so its two
current components match, with opposite outward normals on the two sides.

The sharp connection is spherical on T and dipolar on E. Their difference
at a seam is a smooth tangential one-form. A smooth reference connection
on the joined bundle exists; the difference between its Dirac expression
and the sharp one is a bounded Hermitian multiplication operator for each
fixed q,p. Its bound need not be uniform in q. The principal symbol is
smooth because of (1). A zero-order jump does not create a distributional
delta in the first-order operator; a spinor jump does.

## 2. Domain and surjectivity without a high-mode contraction constant

We give the domain argument, since passivity of a boundary matrix with
norm one would not prove this result. Distributional graph sections have
weak H^-1/2 traces at each smooth noncharacteristic interface: define them
by Green's formula against H1 extensions of H^1/2 boundary tests. The
normal symbol is invertible, so this defines the whole trace. This argument
works on either side before any L2 trace is assumed.

Suppose a distributional graph pair in (2) satisfies the first condition
of (3) weakly. Temporarily regard v as a full throat spinor with zero ground
component. On the joined manifold the jump of (u,v) is then entirely in G.
Its coefficients are well-defined numbers, obtained by pairing the weak
trace with the finitely many smooth ground modes. Thus the jump is smooth
in the angular variables, even though the original trace was only H^-1/2.
Subtract a smooth exterior collar lift of that finite-dimensional jump.
The resulting pair has zero jump and its distributional Dirac image is L2
on the joined manifold: the boundary delta cancels exactly.

Interior elliptic regularity for the smooth reference principal operator
now puts that corrected section in H1 locally through the seam. The bounded
zero-order connection difference is absorbed in the L2 right-hand side.
Adding the smooth lift back shows that u,v are piecewise H1 up to each
seam. In particular, **all interface traces are H^1/2 and L2**. No uniform
in-q H1 estimate is needed for this qualitative step.

Here and below interior elliptic regularity and the H1 trace/extension
theorem are standard analytic inputs. The transmission example and trace
proof in the source audit independently corroborate their use. Neither a
global spectral gap nor compactness of the exterior is assumed.

This also proves the converse graph-domain characterization. Localize an
H1 pair satisfying (3) near the seams and approximate it in piecewise H1
by smooth pairs. Project its throat approximation with P_H, which preserves
every Sobolev order and commutes with the spherical operator. The remaining
trace mismatch and forbidden ground trace tend to zero in H^1/2. Correct
them on the exterior side using bounded H^1/2-to-H1 collar extensions;
smooth input traces have smooth extensions. The corrected pairs satisfy
(3) exactly and converge in graph norm. Away from the seams use ordinary
interior smoothing. At Euclidean infinity, cutoffs with |grad chi_R|<=C/R
have commutators tending to zero in L2. A diagonal sequence gives the core
approximation. The same argument applies with incoming/outgoing exchanged.

On the core the internal flux cancels, leaving

    2 Re <w,(s+K_out)w>
       = 2 nu ||w||^2 + ||Gamma_out,G w||^2,  Re s=nu>0.      (4)

The identity extends by graph closure; it also controls the ground trace
of graph-Cauchy sequences. Arbitrary smooth matched high traces and arbitrary
allowed ground traces in Green's formula show that the Hilbert adjoint is
the formal expression -K with the same high-mode transmission and the
complementary condition Gamma_out,G=0. The preceding finite-jump argument
and density characterize this adjoint domain as well. Its analogous identity
has positive incoming ground flux. Consequently

    ||(s+K_out)w|| >= nu ||w||,
    ||(conj(s)+K_out^*)w|| >= nu ||w||.                       (5)

The first inequality makes the range closed. The second makes its orthogonal
complement, ker(conj(s)+K_out^*), zero. Hence the range is all of Hcal and

    ||(s+K_out)^-1|| <= 1/nu.                               (6)

This proves maximal accretivity in the needed half-plane directly. It does
not invert I-K_load S by a purported uniform Neumann contraction. Indeed
the exact high-mass defect tends to zero, as the earlier load proof shows.

For incoming ground data g, choose a smooth finite-dimensional exterior
collar lift Eg satisfying the high transmission condition and set

    w=Eg+(s+K_out)^-1[f-(s+K)Eg].                           (7)

Its value is independent of the lift. Every g is smooth angularly because
the input space is finite-dimensional for fixed q. The exact identity is

    ||Gamma_out,G w||^2 + 2 nu ||w||^2
       = ||g||^2 + 2 Re <w,f>.                              (8)

Thus the ground scattering map is a contraction and its Poisson map P
has norm <=(2nu)^-1/2. For zero incoming data, square completion also gives
||Gamma_out,G (s+K_out)^-1|| <=(2nu)^-1/2. These constants do not grow with
q or the throat length; the domain regularity constants may do so.

With zero throat forcing, the unique spherical throat solution with its
actual incoming traces is precisely the all-mode radial solution of the load
note. Eliminating it therefore gives in_H=K_load(s) out_H, including its
cross-mouth transmission. Separate integration on E gives the original
exterior current identity, now justified by the proved L2 traces. The
resolvent, current and trace assumptions in the joint and logarithmic
matching arguments are established for this specified realization.

Restoring the ground throat modes and imposing transmission on every mode
gives a conservative operator. In this case the jump is zero, the same
domain argument identifies the adjoint with itself for H, and the two
imaginary shifts are onto. The full H is self-adjoint. This is an operator
statement on a prescribed background, not a selected physical vacuum or a
calculation of its renormalized stress.

## 3. Exact bounded-source derivatives

The domain does not depend on a bounded Hermitian multiplication source V
on E, extended by zero on T. Real amplitudes preserve the identities; small
complex amplitudes preserve the resolvent by its Neumann test. In z=is
notation, write R_z=(H_out-z)^-1. For fixed incoming g the exact return T_z
therefore has the insertion formulas

    partial_a T_z = -Gamma_out R_z V_a P_z,
    partial_a partial_b T_z
       = Gamma_out R_z (V_a R_z V_b+V_b R_z V_a) P_z.          (9)

Formula (9) uses H-z throughout. In the s+K convention the source is iV
and the resolvent changes by the corresponding scalar factor. The norms are

    ||partial_a T_z|| <= ||V_a||/(2nu),
    ||partial_a partial_b T_z|| <= ||V_a||||V_b||/nu^2.        (10)

Unlike the propagation comparison, this realization and insertion statement
does not require axisymmetric sources. It does not bound a continuum vacuum
trace or an unbounded metric/clock insertion. The earlier regulated full
Hessian still requires the exterior determinant and induced quadratic vertex.

## 4. Smooth gauge completion and its quantitative cost

The sharp connection can be replaced by a smooth prescribed connection
without losing the long-return comparison. At either mouth write
A_dipole=A_mono+a_remote in its regular patch. From the exact port formula,

    |a_remote| = q rho sin(theta)/[2 S(S+1-rho cos(theta))]
               <= q rho/[4(1-rho)^2],
    S=sqrt(1+rho^2-2rho cos(theta)).                         (11)

This is the physical orthonormal one-form norm in d=1 units; a_remote has
zero additional flux and is globally regular on this sphere. Take a fixed
smooth chi with chi(t)=0 near t=1 and chi(t)=1 for t>=2, and set in E

    A_smooth=A_mono+chi(rho/p) a_remote.                     (12)

The old monopole throat connection then joins smoothly. Only the two shells
p<=rho<=2p are changed. The operator difference W is Hermitian multiplication,
supported there, and ||W||<=C q p. It is not small in operator norm when
p=q^-3/4, so that norm alone cannot justify a resolvent comparison.

Instead apply W to the existing matched trial Phi. On these shells the
trial is exactly its normalized spherical angular collar with the integral
lower-spin corrector. The latter has norm <=27 rho^2/sqrt(160), uniformly
over the angular ground sectors. Its radial norm is consequently bounded
by C sqrt(p). The source profiles vanish near the mouths; on the complex
source polydisc the longitudinal phase is bounded as in the long-return
proof. Therefore

    ||W Phi|| <= C q p^(3/2).                               (13)

This bound holds on all retained superpositions by azimuthal orthogonality.
The incoming/outgoing projectors and comparison bases at the port are the
same isolated-monopole bases. Thus the old boundary mismatch estimates
remain valid. Apply the now-realized current estimate to the old trial with
its additional residual W Phi. This compares the exact *smooth* return to
the same Wilson/length trial, rather than bounding R_smooth-R_sharp by a
large global operator perturbation norm.

At the previous p=q^-3/4, (13) is O(q^-1/8) and still vanishes faster than
any power of delta=(log q)^-2 after the polynomial trim/resolution losses.
A stronger compatible auxiliary choice is

    p=q^-1, a=c_a q^-1/4, g^2=q^-2,
    delta=(log q)^-2, nu=delta^2.                            (14)

Then (13) is O(q^-1/2), smaller than the O(q^-3/8) bulk residual. The old
weighted boundary cost only improves: its largest port term is
nu^-1 q^(3/2) p^3=nu^-1 q^-3/2. Inserting p=q^-1 into every displayed
port term of joint mouth matching verifies this claim. The proof of that
note uses p<a and the listed smallness conditions, all satisfied here.
At the leading MMP scale ell/r_e=16pi q^3, the physical port has
rho_p/r_e=16pi delta q^2 -> infinity. This is a scale check, not a metric
perturbation theorem.

For (14) the full logarithmic return and its first two source derivatives
retain their O(q^-3/8 times powers of log q) asymptotic bounds for any fixed
finite resonant harmonic set. The loaded-existence qualifier is removed
for (1),(12); the geometry is prescribed and the finite-q threshold is still
uncomputed. The earlier large sampled bulk residuals are not revised away.

There is a precise obstruction to making this particular interpolation
source-free. In a connected Euclidean shell with static vacuum Maxwell
equations, div B=curl B=0, every Cartesian component of B is harmonic.
If B equals the isolated monopole field on an open inner subannulus,
analytic continuation of their harmonic difference makes them equal on
the whole shell. They therefore cannot also equal the point-dipole field
on an open outer subannulus, where the remote monopole field is nonzero.
This concerns the imposed flat, exactly spherical collar; it is not an
obstruction to a nonspherical Einstein-Maxwell throat.

One can see the same issue without a vector unique-continuation theorem.
The nonzero uniform-field multipole of the remote field has potential
a_phi=(B0/2) rho^2 sin^2(theta). Write chi_p(rho)=chi(rho/p).
The radial cutoff changes its azimuthal curl to

    (curl B)_phi=-(B0/2)[rho chi_p''(rho)+4 chi_p'(rho)] sin(theta). (15)

Zero current would imply rho chi_p''+4chi_p'=0. Its solutions are a+b/rho^3;
none is zero on an inner open interval and one on an outer open interval.
Other angular multipoles cannot cancel this multipole equation. Equation
(15) is checked independently by differentiating the spherical curl formula.
For a transition of width p its current is of order |B0|/p, even though
its field is small relative to the monopole field. A small relative field
does not remove this source-free consistency requirement.

The next physical matching step must thus permit the remote multipoles
inside the throat or a controlled metric/matter response. An arbitrary
flat spherical cutoff cannot supply that step. Neither its extra current
nor its stress/backreaction is accepted as matter in the selected model.

## 5. Verification and source boundary

Run `checks/loaded_operator_realization.py --output
evidence/loaded-operator-realization.json`. Its coupled-channel networks
retain reflected high modes and compare a direct real-component boundary
value solve with independent scattering-matrix feedback. They verify the
total current, Poisson bound, inhomogeneous resolvent bound, and the defect
tending to zero at large angular mass. Two BVP tolerances/meshes are used.
The gauge shell checks evaluate (11) and an all-angular-profile upper
integral for (13), with exact rational port exponents. These are supplementary
component checks, not a three-dimensional dipole or quantum-field simulation.

The analytic tools are established. Bär and Ballmann's
[2011 primary preprint](https://archive.mpim-bonn.mpg.de/id/eprint/327/1/preprint_2011_2.pdf),
Theorem 6.7 and Example 7.28, give the relevant trace and transmission
framework. The proof above specializes finite-dimensional ground jumps
and dissipative ground ports to this loaded geometry. The source receipt
records the version, pages and dependencies actually inspected; the full
paper and its entire regularity chain are not declared newly reviewed.
No claim that general transmission existence is novel is made.

Remaining to-do list: replace prescribed collars by a controlled physical
Einstein-Maxwell interface and state; derive full angular/determinant and
discarded-tail response and renormalized stress; complete novelty review
and the consequential endpoint; independently verify it, certify suitable
mathematics exclusively on GCP, and prepare the reproducible manuscript.
