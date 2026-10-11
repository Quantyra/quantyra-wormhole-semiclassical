# Source-free Maxwell matching through the two mouths

2026-10-11, S053. Solve the exterior and throat Maxwell boundary problem
together. For the prescribed smooth geometry with Euclidean end collars,
there is a source-free magnetostatic completion of the fixed flux. Both
the tangential potential and its conjugate normal flux match. Its squared
field correction is O(q^2 p^3), uniformly in the intervening spherical
throat profile. The classical energy correction is subleading to q/ell
in a compatible auxiliary-port family, p=q^-3/2.

This resolves the artificial collar current of the previous cutoff
construction on a prescribed metric. It does not solve the Einstein
equations, bound the local stress response, select a quantum state, or
transfer the previous Dirac-return estimates to the corrected connection.
An integrated classical energy bound does not supply those conclusions.

## 1. Fixed geometry, flux and normalization

Lengths are in separation units d=1 until section 6. Use the connected
flat exterior E of two radius-p balls centered at z=+/-1/2, 0<p<=1/10.
Join their boundaries by the spherical segment

    ds^2=N(x)^2(-dt^2+dx^2)+R(x)^2 dOmega^2, 0<=x<=L,
    w(x)=N(x)/R(x)>0.

Keep the smooth Euclidean end collars of
[the loaded realization](loaded-operator-realization.md):
N=1, R=p-x on 0<=x<=p/4 and the reflected condition at the other end.
The collars are disjoint. The middle N,R are any smooth positive functions.
These conditions give w<=2/p on the collars. Constants below do not depend
on the length or the middle profile. They may depend on fixed collar ratios.

The magnetic flux through the handle is 2pi q in the established connection
normalization. Work in the static, axisymmetric, purely magnetic class.
Regular gauge patches carry the monopole connection; each zero-flux
perturbation is a globally regular azimuthal one-form a dphi. At an axis,
regularity means a vanishes quadratically with cylindrical radius.
At infinity use the finite-energy, no-applied-field completion of compact
perturbations. A flat Wilson-line phase can be fixed separately; it has no
Maxwell field energy and is not determined by this construction.

On a spherical region the Maxwell energy, with its common factor pi/g^2
removed, is

    Q_T[a]=integral_0^L integral_0^pi
              [|a_x|^2+w(x)^2 |a_theta|^2] dtheta dx/sin(theta). (1)

On E, the same normalized form is the Euclidean magnetic energy. In a
spherical annulus it is (1) with x=rho and w=1/rho. Equivalently, in
cylindrical meridional coordinates it is

    Q_E[a]=integral_E,meridian (|a_r|^2+|a_z|^2) dr dz/r.       (2)

These follow directly from (2g^2)^-1 integral N |B|^2 dvol. Varying (1)
gives a_xx+w^2 sin(theta) partial_theta(a_theta/sin(theta))=0.
All other components of the static Maxwell divergence vanish for this
axisymmetric azimuthal potential, so this scalar equation is the full
source-free Maxwell condition in the chosen class.

## 2. Boundary space and the two Dirichlet-to-Neumann maps

Let B=L2((0,pi),dtheta/sin(theta)) for zero-flux potential traces. The
positive angular operator is

    A=-sin(theta) partial_theta [sin(theta)^-1 partial_theta],
    Y_l(theta)=(1-cos(theta)^2) P_l'(cos(theta)), l>=1,
    A Y_l=l(l+1)Y_l,
    ||Y_l||_B^2=2l(l+1)/(2l+1).                               (3)

The regular endpoint realization is understood. There is no l=0 zero-flux
mode. Put X=D(A^1/4), with squared norm sum_l sqrt(l(l+1)) |h_l|^2 in the
normalized basis. At the two ports use X direct_sum X. Orientation and
gauge identifications act isometrically; at the reflected mouth the
global-axis l basis changes by (-1)^(l+1).

For h in this two-port space, minimize (2) among finite-energy regular
exterior perturbations with trace h. Its minimum defines
Q_E(h)=<h,Lambda_E h>. The exterior here is the actual connected
two-ball exterior; Lambda_E is not a sum of two isolated-sphere maps.
Likewise minimize (1) with endpoint values h to define
Q_T(h)=<h,Lambda_T h>. The throat multipoles obey

    a_l''=l(l+1) w^2 a_l.                                    (4)

Its DtN block maps (a_l(0),a_l(L)) to (-a_l'(0),a_l'(L)).
The maps are defined as quadratic forms X to X*, not presumed bounded
operators on B. All angular modes enter these definitions.

Here is the uniform coercivity estimate needed for the matching. Restrict
an exterior field to either annulus p<=rho<=2p. The lth contribution is
at least the one-dimensional energy with potential l(l+1)/(4p^2), fixed
inner value and free outer endpoint. Minimization gives

    Q_E(h)>= (c_*/p) ||h||_X^2,
    c_*=(1/2)tanh(sqrt(2)/2)>0.                              (5)

The two annuli are disjoint. The estimate also follows for limits in the
energy space. Conversely, exponential radial extensions of each h_l,
cut off inside these annuli, give Q_E(h)<=C ||h||_X^2/p. Extensions into
the two throat collars, cut off before their far endpoints, give

    0<=Q_T(h)<=C ||h||_X^2/p.                               (6)

For example exp[-sqrt(l(l+1)) x/p] times a fixed collar cutoff has the
required order. Its cutoff derivative costs the same order since l>=1.
No property of w in the throat middle is used.

These estimates justify the minimizers and the function spaces. Smooth
regular azimuthal fields are completed in the magnetic-energy norm;
fixed traces form closed affine subspaces by (5), and collar extensions
show they are nonempty. Orthogonal projection gives the unique minimizer,
modulo pure gauges, which are fixed by axis regularity and the condition
at infinity in this azimuthal class. In the zero-trace exterior subspace,
extension by zero through the balls and the Coulomb Sobolev estimate
also make the energy completion concrete. The azimuthal potential is
divergence-free. The throat argument is the corresponding positive
one-dimensional variational problem, summed in the normalized basis.

## 3. The actual coupled transmission equation

The piecewise reference field F_ref is the exact point-dipole field on E
and the isolated flux-q monopole on T. Its boundary potential mismatch is
the remote connection v. At the first sphere, in its regular patch,

    v(rho,theta)=q rho^2 sin(theta)^2/[2 S(S+1-rho cos(theta))],
    S=sqrt(1+rho^2-2rho cos(theta)).                         (7)

At the second sphere apply the orientation reflection. This is the exact
formula from [mouth modes](mouth-mode-matching.md), not a leading dipole
Taylor truncation. Let f be the outward exterior normal flux of this
background for a potential variation. In local spherical coordinates
f=-partial_rho v at rho=p; the isolated monopole contributes no such
radial derivative. The remote field is smooth one separation away, hence

    ||v||_X<=C q p^2,       ||f||_(B direct_sum B)<=C q p.     (8)

For (8), differentiate the rationalized expression (7); S>=1-p and its
denominator stays away from zero. v and its first angular derivative are
O(qp^2) with the regular sin(theta)^2 endpoint behavior. Interpolation
between the B and A^1/2 norms gives the stated X bound. Radial
differentiation gives the second estimate.

Let h be the actual throat perturbation trace relative to its monopole.
The exterior correction then has trace k=h-v. Its energy cross term with
the exterior reference is <f,k>, by Green's formula. The throat monopole
has zero cross term: its angular derivative is q sin(theta)/2, so the
angular integral is proportional to a(pi)-a(0)=0. Thus the change in
normalized Maxwell energy is

    J(h)=2<f,h-v>+Q_E(h-v)+Q_T(h).                          (9)

All fields and forms here are real; complexification only replaces the
cross term by its real part. Stationarity of (9) is exactly

    (Lambda_E+Lambda_T)h=Lambda_E v-f.                     (10)

By (5)-(6), this equation has a unique solution in X direct_sum X by
the coercive variational theorem. Importantly, (10) enforces the normal
flux condition as well as the already matched potential trace. Prescribing
only h=v and solving a throat Dirichlet problem would not do so.

The reconstructed field is closed and has the original flux. Equation
(10) cancels the interface term in the global Maxwell variation, hence
there is no surface current or charge. It solves the full magnetostatic
Maxwell equations weakly. The metric and lapse are smooth through the
Euclidean seams; local elliptic regularity then makes the field smooth
there and in the interior. Regular angular domains remove axis artifacts.
This is a source-free field on the prescribed geometry, rather than a
cutoff maintained by an additional current.

## 4. Uniform field and energy correction

Test (10) against h. Use Cauchy-Schwarz for Q_E and (5) for the f term:

    Q_E(h)+Q_T(h)
       <= sqrt(Q_E(v)) sqrt(Q_E(h)) + ||f||_B ||h||_B
       <= C q p^(3/2) sqrt[Q_E(h)+Q_T(h)].                 (11)

In the last inequality, (8) and Q_E(v)<=C||v||_X^2/p were used. Therefore

    Q_E(h)+Q_T(h)<=C q^2 p^3,
    ||h||_X+||h||_B<=C q p^2,
    Q_E(h-v)+Q_T(h)<=C q^2 p^3.                           (12)

The final line is the squared magnetic-field correction, measured
piecewise against F_ref in the positive static energy norm. The reference
itself need not be a globally matched Maxwell field for this comparison
to be well defined. Equations (8),(9),(12) further give

    |J(h)|<=C q^2 p^3.                                    (13)

No sign of J is asserted. The total energy is positive, but changing the
reference field can increase or decrease it. The proof uses all modes and
does not infer the bound from the numerical multipole truncation.

## 5. What changes for the matter calculation

The artificial flat collar current is removed. The resulting smooth
connection, with a separately fixed flat holonomy, also admits the usual
self-adjoint full free Dirac realization on this complete prescribed
static background. This qualitative fact follows from the transmission
argument of the previous note; it is not a quantitative return comparison.

The field in E is now corrected as well as the field in T. The latter is
not exactly spherical, so the isolated angular ground projector need not
commute with the throat Dirac operator. Its former diagonal higher-mode
load cannot simply be reused. Bounds (12) are energy norms of fields,
not uniform multiplication-operator norms of potentials or bounds on a
vacuum trace. Controlling the actual channel decomposition, induced mode
mixing and return/source derivatives is the next matter step.

Nor does (13) bound the local stress tensor difference, its length
derivatives, or the Einstein response. Angular cancellations enter the
integrated energy. A metric perturbation can couple to stress components
that these cancellations do not control. This prevents promoting the
following scale comparison to a solved semiclassical background.

## 6. Compatible physical scale and rejection of p=q^-1 for this bound

Restore d. The classical Maxwell energy correction is

    |Delta E_Maxwell|<=C q^2 p^3/(g^2 d),                 (14)

where p remains the dimensionless radius rho_p/d. The relevant quantum
energy scale is q/ell=q delta/d. With g^2=q^-2,

    |Delta E_Maxwell|/(q/ell)<=C q^3 p^3/delta.            (15)

For p=q^-1, the right side is C/delta and does not prove smallness. This
does not prove that the actual correction is large; it means that the
previous port choice does not close this energy estimate. Choose instead

    p=q^-3/2, a=c_a q^-1/4, delta=(log q)^-2,
    g^2=q^-2, nu=delta^2.                                 (16)

Then (15) is O(q^-3/2 (log q)^2)->0. The leading scale
ell/r_e=16pi q^3 gives rho_p/r_e=16pi delta q^(3/2)->infinity,
so p<a and the port stays outside the gravitational radius in this
comparison. The earlier *reference-connection* collar and boundary error
orders improve as p decreases; the Maxwell-corrected Dirac errors still
need their own estimates as explained above. This is the same weak-coupling
theory family, with a different auxiliary cut, not a fixed-coupling theorem.

## 7. Verification and prior work

Run `checks/coupled_maxwell_matching.py --output
evidence/coupled-maxwell-matching.json`. The check uses actual exterior
two-sphere multipoles, including translations between the mouths, and
continuous boundary-value throat DtN blocks independently checked against
conforming finite elements.
It varies the retained angular order, angular quadrature and throat mesh,
and tests the unprojected boundary/current mismatch as well as the
projected stationarity equation. The throat tests include a length/p ratio
of 1000. They are convergence checks, not validated numerics or a proof
of the all-mode theorem. Exact angular and rational scale identities are
checked separately. No Einstein or quantum-field simulation is performed.

Source-free wormhole flux and electromagnetic multipoles are established
subjects. The targeted source review includes
[Krasnikov, 0802.1358v1](https://arxiv.org/abs/0802.1358v1), whose Section 2
proves uniqueness up to the flux sector for a spherically symmetric
two-end electrostatic problem. Our two mouths share one exterior, and
their nonspherical fields are fixed by the same flux and transmission
geometry. They are not additional freely specifiable source-free
multipole charges. The existing MMP construction remains the physical
baseline. Neither this matching method nor flux without local sources
is accepted as a new consequential wormhole result.

Remaining to-do list: control the corrected connection in the full Dirac
return and channel mixing; derive metric and local-stress backreaction,
state and renormalized response with exterior determinant and discarded
tail; complete literature/proof/novelty review and the consequential
endpoint; independently verify it, certify suitable mathematics only
on GCP, and prepare the reproducible manuscript.
