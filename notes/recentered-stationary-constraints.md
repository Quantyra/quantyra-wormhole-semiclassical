# Stationary current, off-shell constraints and the weight-four gate

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The current equations and static radial recovery hold throughout the
selected identical magnetic disk family. This supplies a missing
hypothesis of the reduced solver, not its nonlinear solution. We derive
the off-shell reconstruction with both density and residual contacts.
The remaining endpoint problem is real: the shape normal family has a
nonzero pole at weight four. Its finite-part response contains a global
moment; merely continuing the normal symbol does not solve that problem.

## 1. Fixed construction and scope

Keep c_R=10^22, the full vectorlike Dirac field, fixed magnetic bundle,
zero electric source/chemical potential, positive fixed thermal period,
identical clocks and paired theta=0 domain. All other additional finite
coefficients remain zero at the retained mu0. The reference and physical
fast phase are those of [the mass-boundary audit](recentered-mass-boundary.md).
Work on a smooth nonextremal disk times the round sphere, initially at a
finite spherical cut, then in the retained common-subtraction limit.
No assertion is made about a driven join, an electric susceptibility,
a different state, an extremal disk or Lorentzian stability.

The domain is the actual one in [the paired-phase derivation](stationary-paired-dirac-phase.md):
chi=tau3 gamma_n, Dom D has Pi_- psi=0 and Dom Ddagger has Pi_+ psi=0.
The squares have complementary mixed data, including the K/2 Robin
coefficient. D itself is not self-adjoint on its first-order domain.
The full positive singular-value determinant, including the paired phase
convention, is used. These domain qualifications cannot be replaced by
an unspecified bag boundary or by an LLL determinant.

The source boundary condition is consistent with
[Bintanja--Freivogel--Rolph, (2.8)-(2.12)](https://arxiv.org/html/2308.00871v2#S2.SS1.SSS2).
The spin connection and mixed-operator conventions are compared with
[Vassilevich, sections 3.3 and 5.3](https://arxiv.org/html/hep-th/0306138v3).
Neither reference proves the current cancellation or nonlinear theorem
for the present construction; the domain intertwiner below is needed.

## 2. An actual symmetry at fixed magnetic flux

Let c_S be the sphere chirality, globally defined on the spinor bundle
twisted by the monopole line bundle. It is parallel, has square one,
and anticommutes with the magnetic sphere Dirac operator D_q, including
its zero modes. It does not conjugate or change the monopole bundle.
In the orientation-aware frame of the paired-phase note put

    gamma_tau=-tau3 tensor sigma3 tensor 1,
    gamma_n=1 tensor sigma1 tensor 1,
    gamma_A=-tau3 tensor sigma2 tensor gamma_A,S,
    U=1 tensor sigma1 tensor c_S.                         (1)

Then U is unitary and Hermitian, U^2=1,

    U gamma_tau U=-gamma_tau,
    U gamma_n U=gamma_n, U gamma_A U=gamma_A,
    [U,chi]=0.                                           (2)

Thus time reflection tau->-tau has a lift R to the full Dirac bundle:
in this radial frame it is U followed by pullback. This is a reflection
of the Euclidean disk, with no reflection of the sphere or magnetic
connection. The full vectorlike spinor, rather than a Weyl subbundle,
admits this orientation-reversing lift. The spin connection transforms
by its Clifford bivectors under (2). One can also check smoothness at
the bolt directly: the polar spin rotation S(phi) satisfies
U S(-phi)=S(phi) U. The polar-frame expression therefore becomes a
constant reflection matrix in the smooth Cartesian disk frame. It
preserves the regular antiperiodic polar spin structure and extends
over the center. There is no added condition at the bolt.

For a static spherical orbit gauge field A=A_tau(r)d tau+A_r(r)dr,
with A_tau=O(rho^2) at the bolt, covariance gives as closed operators

    R D[A_tau,A_r] R^-1=D[-A_tau,A_r].                    (3)

Both first-order domains are preserved by (2); composition then
intertwines both squares, their Robin data and their heat kernels.
This is equality of the finite-cut operators and their spectra, not
only an equality of a few local heat coefficients. The sphere angular
tower and all contact vertices are retained. The inherited kernel-free
disk result permits differentiation locally about A_tau=0. A numerical
uniform gap is not needed to differentiate this exact symmetry, and
is not supplied here.

The determinant modulus is even in A_tau. Its specified paired phase
introduces no varying parity-odd term. The full common subtraction is
also even: the local F^2 terms have no magnetic/electric cross term,
and no additional parity-odd gauge or gravitational term was selected.
Consequently the renormalized functional is even, and its derivative
against every compact static A_tau test is zero. This proves j^tau=0
pointwise in the interior. Symmetry is exact before cut removal, so
no estimate of a small residual odd contribution is needed. It also
annuls any allowed boundary charge response odd in the electric source.

SO(3) invariance forbids an angular one-form current. Vector gauge
invariance of this nonchiral determinant gives current conservation.
For a static radial current, a_orb r^2 j^rho is constant. A smooth
rotation-invariant vector at the disk center has j^rho=O(rho), while
a_orb=O(rho); hence that constant is zero. Thus the complete stationary
current vanishes. Reflection likewise sets the orbit off-diagonal
stress to zero. The magnetic field q sin(theta)dtheta wedge dphi
solves the source-free Maxwell equation for every metric in this warped
family. All stationary metric derivatives of the zero-current identity
vanish wherever the corresponding derivative exists, even off shell.
This statement does not set the electric/electric response to zero.

An independent algebraic check is the folded Hamiltonian:
H0=-i tau3 sigma2 partial_x+sigma1 a(x)D_q. The same U anticommutes
with H0 and preserves chi. It pairs opposite energies without changing
q or exchanging the copies. The covariant disk proof above is what
justifies the bolt and thermal domains; a formal Hamiltonian symmetry
alone would not suffice there. Charge conjugation alone would change
the magnetic bundle and is not the argument being used.

## 3. Exact recovery in the fixed disk chart

Use a cylindrical disk coordinate x=log xi, with bolt x=-infinity,
azimuth phi of period 2pi, and write

    g=exp(2L)(dx^2+dphi^2)+r^2 dOmega^2,
    r=B exp(p), L=p+v, v=V+x,
    W=exp(2L) r^2=B^2 exp(4p+2v).                       (4)

The fixed period is already absorbed in this chart. Tangents are
s=delta p, h=delta v=delta V and d=delta log W=4s+2h.
Let E_x,E_phi,E_o be the mixed components of the complete tensor
residual, with the sphere component repeated. Suppress only the common
nonzero action normalization in defining

    t=E_x+E_phi+2E_o, o=E_x+E_phi,
    f1=W t, f2=W o.                                      (5)

These are the density-correct reduced rows; they have the same zeros
as the two physical-volume rows used in the reference inverse.
Gauge/diffeomorphism invariance and the current result give exactly

    E_x'+L'(E_x-E_phi)+2p'(E_x-E_o)=0,
    (W E_x)'=v' f2+p' f1.                               (6)

All quantum and R^2 contributions obey this identity as variations of
the complete action with the specified state/domain. The curvature
auxiliary variable, if used, is eliminated by its defining equation;
an independent off-shell auxiliary field would add its own Euler term.
Boundary sources are fixed, and compact diffeomorphisms suffice for
the local identity. Smoothness at the bolt implies W E_x->0. Therefore

    E_x(x)=W(x)^-1 integral_(-infinity)^x
                              [v' f2+p' f1](y)dy.       (7)

This holds for actual residuals of arbitrary admissible metrics in the
family, not arbitrary independently prescribed tensor components.
If f1=f2=0, (7) forces E_x=0, hence E_phi=E_o=0. There is no separate
static radial integration constant left. This reuses the exact-solution
Ward recovery already present in the disk note; the new ingredients
are its verified current hypothesis and the off-shell formula below.

Let b_x,t0,o0,f10,f20 denote reference residuals, a dot a directional
derivative, and e_x=dot E_x. Differentiate (6) BEFORE imposing any
equation on the reference:

    dot f1=W(dot t+d t0), dot f2=W(dot o+d o0),
    [W(e_x+d b_x)]'=v' dot f2+p' dot f1
                                      +h' f20+s' f10.  (8)

Thus, with the same regular bolt constant,

    e_x=-d b_x+W^-1 integral_(-infinity)^x
                [v' dot f2+p' dot f1+h' f20+s' f10]dy.  (9)

This is the off-shell constraint reconstruction for the coordinate
Hessian. It includes both the moving-density contact and the connection
times reference-residual contribution. In particular dot f1=dot f2=0
does NOT imply e_x=0 on the recentered RN reference. For a Newton step
dot f_i=-f_i0, (8) yields

    [W(b_x+e_x+d b_x)]'=h' f20+s' f10.                  (10)

The right side is of order correction times reference defect; it must
stay in the Newton remainder. It is not a new freely adjustable source.
A useful necessary bound follows directly from (9), replacing each
integrand by its absolute value. No complete admission estimate for
that integral is claimed by writing this identity. At an exact solution
all b and f0 vanish, and the usual homogeneous linear recovery follows.

The Lorentzian constraint propagation equations and counterexamples in
[the retarded recovery note](spherical-retarded-constraint-recovery.md)
remain in force. A regular Euclidean bolt cannot replace Lorentzian
initial constraints and incoming gravitational boundary data.

## 4. Execute the endpoint diagnostic before extending estimates

The reference inverse only gives a mass cubic plus W_(7/2). Even the
cut-off power z^(15/4) lies in that remainder space and fails every
O(z^4 log^p z) bound with finite p. The inverse must therefore be
upgraded for its actual source before the mass-boundary theorem applies
to its output. This is an unknown endpoint theorem, not a measured large
constant and not a license to change the phase.

Use precisely the retained TT normal symbol (normal TT amplitude,
before the fixed scalar/TT coordinate conversion):

    q(s)=s(s-3), d0=1/(480pi^2),
    F_M(s)=psi(s-1)+(pi/2)cot(pi s/2),
    kappa_M=-6/5-log(mu0 ell),
    A_T=F_lambda/(48pi G lambda)+23/(8640pi^2),
    N_T(s)=d0 q(s)[q(s)+2][F_M(s)+kappa_M]-A_T q(s).

Writing e=s-4, direct Laurent expansion gives

    N_T(4+e)=24d0/e+C4+O(e),
    C4=d0[86-24 gamma_E+24 kappa_M]-4A_T,
    1/N_T(4+e)=e/(24d0)+O(e^2).                        (11)

The residue is 1/(20pi^2), about .00506605918212; it is nonzero.
The reciprocal slope is 20pi^2, about 197.392088022. For a cut-off z^4
tail, finite-part continuation therefore contributes

    24d0 z^4 log z + (C4+C_remote) z^4.                 (12)

The logarithmic coefficient is fixed. C_remote depends on the cut-off
and the complete interior kernel; (11) does not specify it. The leading
homogeneous remote kernel is -24d0 y^-5, so its integral against the
tail produces exactly this logarithm. The actual kernel adds its own
finite global moment and lower jets. Mass-cubic subtraction and all
local contacts must precede defining that moment. It cannot be discarded
as a harmless faster remainder at the endpoint.

For perspective, the meromorphic extension of N_T has a zero just to
the right of the pole. Put e0=6d0/A_T. At the retained benchmark,

    1.86e-34<e0<1.87e-34,
    N_T(4+e_*)=0, |e_*/e0-1|<1e-29,
    5.34e33<1/e_*<5.38e33.                              (13)

The checker executes the scaled root calculation at 100/150 digits;
the existence of this local meromorphic zero also follows by opposite
signs and strict negative derivative near e0. The derivative of its pole
dominates the bounded analytic remainder there. Crucially the actual
homogeneous kernel was justified only for -1<Re s<4. This zero is NOT
an established physical mode, instability, new boundary datum or
eigenvalue of the complete operator. A power with that exponent would
remain indistinguishable from z^4 on log-radius intervals much shorter
than 1/e_*. This is a warning against inferring endpoint estimates from
finite-radius fits, not an obstruction theorem.

The decisive next endpoint calculation is the finite-part remote moment
of the actual Green response, coupled to the forced scalar z^4 log z
coefficient and the actual source. It must prove five differentiated
fast-tail bounds for the solution, retaining the mass cubic, and bound
their constants in the physical admission norm. Computing more normal
roots or improving an interior constant cannot substitute for it.
This turn executes the residue/logarithm and scale diagnostic, not the
unproved actual-kernel moment theorem.

## 5. Quantitative status and research decision

Continue this specified construction with the sharper method: recover
the actual endpoint expansion and estimate the actual preconditioned
correction, with the complete remote moment. The current and static
constraint requirements need exact identities; they are now supplied
on the selected stationary family. They consume no small-error budget.
Keep the existing exact mass boundary cancellation and reference energy
inverse bound 1. The finite admission inverse is still unquantified.

For a fixed physical Banach ball ||u||_X<=rho and preconditioner P, let
Y0=||P F(0)||, Z0=||I-P DF(0)|| and
||P(DF(u)-DF(0))||<=Z2 ||u||. Preservation and contraction require

    Y0+Z0 rho+(Z2/2)rho^2<=rho,
    Z0+Z2 rho<1.                                        (14)

P must be injective on the full required output, compatibility must be
met, and each physical admission margin a_i must obey
|delta a_i(u)|<=a_i(0) with strict positivity where required. At
rho=1e-8, the sufficient targets remain Y0<=5e-9, Z0<=.25,
Z2<=5e7: the first bound is exactly 1e-8 and the second is at most .75.
If bounds E_i rho control the admission changes, each missing E_i must
satisfy E_i<a_i(0)/1e-8. Unknown a_i and E_i need actual evaluation.
An energy-space defect <=.25 would give inverse <=4/3 in those SAME
spaces; it does not establish these admission-norm targets.

The original C_I=9.1591500000000000015e103 and its failed uniform
certificate remain parked. They are known upper bounds that fail to
certify the old target, not lower bounds on a solution correction. The
new Y0,Z0,Z2 and endpoint constants are unknown. Changing normalization
does not improve any invariant inclusion or contraction inequality.
Fixed-flux/temperature compatibility, matching about the new full
reference, omission/noise control and all-sector retarded stability
remain independent earlier-to-later gates. No finite-state solution or
controlled wormhole has been constructed here, and S169.2 stays open.

Before execution freeze this argument, checker and source review. Check
the Clifford/domain intertwiner, polar lift and full angular-pair
algebra exactly. Independently derive Ward recovery using the Christoffel
symbols of a general spherical metric, then differentiate an actual
off-shell Einstein residual along two analytic disk deformations; keep
both contacts. Check the Laurent coefficients, the scaled normal root
and the nonlinear target arithmetic. Preserve every run and frozen
snapshot. These checks support the argument; they do not verify a
continuum kernel limit or the missing endpoint theorem. No Lean run.

Remaining to-do list: actual-source weight-four moment and physical
fast-tail/admission bounds; finite geometry/state, residual and nonlinear
solution; exact compatibility, same-source matching, semiclassical
validity, all-sector stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
