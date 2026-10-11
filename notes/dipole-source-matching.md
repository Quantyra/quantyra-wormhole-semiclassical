# Photon-source derivatives of the trimmed dipole return

2026-10-11, S053. The actual flat point-dipole operator admits a gauged
extension of the existing global mouth trial. A line-integral phase and
one transverse spin correction control its first two source derivatives,
including a shrinking frequency resolution. This replaces a zero-source
comparison by a source-dependent one for the stated axisymmetric profiles.
It remains conditional on the specified loaded boundary realization.
The resonant shrinking flux sector, general photon angular mixing,
continuum determinant and renormalized physical stress remain open.

## 1. Model, units and statement

Lengths in this note are in units of the mouth separation d. Thus
`z=d omega_phys+i nu`, `nu=d eta_phys>0`, and a physical potential is
rescaled as `V=d V_phys`. Use the flat two-mouth point-dipole Weyl operator
and spherical higher-mode load of [joint mouth matching](joint-mouth-matching.md).
The load is left attached on all higher angular modes; the ground input
and output ports are open. Keep flux labels `u in [u_*,1-u_*]`, where
`u_*>0` is fixed, and physical sectors `q(1-u) in Z+1/2`. Set

    p=q^-3/4, a=c_a q^-1/4, |Re z|<=Omega, 0<nu<=1.          (1)

The port radius is p and the collar/tube overlap radius is a. The load
has the two collar bounds already stated in the linked note, with
constants independent of q. Assume its coupled exterior problem has
the current/trace identity and a zero-incoming resolvent with norm
`<=1/nu`. This is a well-posedness assumption, not a consequence of the
trial. In particular, the uncommitted flat-boundary argument is not used
as an accepted construction of the loaded physical operator.

Let `V_j` be finitely many fixed, bounded Hermitian axisymmetric Pauli
multiplication potentials with bounded C2 coefficients on the retained
tubes. They vanish on fixed neighborhoods of both mouths. The profiles,
their regularity bounds and the neighborhoods are independent of q.
These are static external photon/electric sources, not dynamical photon
fields. Axisymmetry preserves the azimuthal sector; angular mixing beyond
that class is not claimed. Put `H(t)=H+sum_j t_j V_j`.

In the field frame write `V=V_0 I+V_x sigma_x+V_y sigma_y+V_z sigma_z`,
and let `v_j(s)=(V_j)_{++}(s,0)`. Define the comparison return

    R_W(t;u,z)=exp[i z L_p(u)-i sum_j t_j I_j(u)],
    I_j(u)=integral_0^L_p v_j(s) ds.                        (2)

For a physical vector/electric source, `v=V_0-A.tangent` in the chosen
units. Equation (2) has the minimally coupled Wilson phase. The gauge
and spin-frame phases of the ports are fixed as in the joint trial;
no closed-loop state or holonomy is selected here.

**Conditional return theorem.** With the above assumptions and fixed
source profiles, at t=0,

    ||R-R_W|| <= C q^-3/8 nu^-1/2,
    ||partial_j(R-R_W)|| <= C_j q^-3/8 nu^-3/2,
    ||partial_j partial_k(R-R_W)|| <= C_jk q^-3/8 nu^-5/2.  (3)

The constants may depend on the fixed trim, frequency bound, sources and
collar constants, but not on q or nu. The norm is on all retained ground
inputs, with all higher modes present in the error problem. Repeated
derivatives are included; their Cauchy factorial changes only C_jj.
In particular,

    partial_j R_W=-i I_j R_W,
    partial_j partial_k R_W=-I_j I_k R_W.                   (4)

These are approximations to the derivatives of the exact return, not
derivatives inferred solely from a small zero-source error.

## 2. The transverse source corrector

In the exact field-aligned frame, with the prior notation,

    H_F=i sigma_x partial_y-i sigma_z h^-1/2 partial_s h^-1/2
                                       +q sigma_y m(s,y),
    h=1-k(s)y, b=m_y(s,0)>0, c=m_yy(s,0),
    g=(q b/pi)^1/4 exp(-q b y^2/2), omega=b_s/(4b),
    P_0=1-q c y^3/6, Q=omega y.

For the combined source V(t), let

    v(s)=V_++(s,0), w(s)=V_-+(s,0),
    P=P_0+i w(s)y,
    F(s)=exp[i z s-i integral_0^s v(s')ds'],
    Phi_F=F(s) chi(s,y) g(s,y) (P,Q).                      (5)

The tube cutoff and both collar pieces are those of the joint trial.
The source is zero near both overlaps for sufficiently large q, so the
extra transverse term is zero there. The phase is constant relative to
the zero-source phase on each collar. It therefore multiplies the
receiving collar by exactly the factor in (2), without an exposed
longitudinal cutoff or new overlap discrepancy.

The upper projected source is canceled by `F'=i(z-v)F`. The remaining
order-one lower source `w g` is canceled because

    i(partial_y+q b y)(i w y g)=-w g.                      (6)

The normal and azimuthal Pauli components enter w. Discarding them while
retaining only the phase leaves this order-one residual. No positive gap
for the complete noncompact exterior is assumed in (6).

Here is an exact residual identity for independent verification. Set
`G=g_s/g=omega(1-2q b y^2)`, `h_s=-k_s y`. Before the cutoff, the
two components of `(H_F+V-z)Phi_F/(F g)` are

    U=i[Q_y-q(b y+m)Q] - i(G P+P_s)/h + i h_s P/(2h^2)
                       +[(z-v)/h+V_++-z]P + V_+- Q,
    W=i[P_y+q(m-b y)P] + i(G Q+Q_s)/h - i h_s Q/(2h^2)
                       +[-(z-v)/h+V_---z]Q + V_-+ P.       (7)

The term `V_---z` means `V_{--} - z`. The cutoff adds
`i sigma_x chi_y g(P,Q)-i sigma_z chi_s g(P,Q)/h`, with the
same longitudinal phase. Formula (7) uses the exact dipole mass m.

On the fixed compact part of the retained tubes where the source is
nonzero, geometry and C2 source bounds are uniform. Relative to the old
residual, the upper new terms contain `V_++(s,y)-v(s)/h`, `V_+- Q`,
and derivatives or bounded coefficients times `i w y`. Their Gaussian
norms are O(q^-1/2). In the lower component, after (6), the new terms are
`V_-+(s,y)-w(s)`, `V_-+(P_0-1)`, `V_-+ i w y`, a bounded coefficient
times Q, and `-q w y[m-b y]`. Their norms are O(q^-1/2) because
`m-b y=O(y^2)` and `||q^r y^k g||=O(q^(r-k/2))`. Source-dependent
cutoff terms are exponentially small. These bounds are uniform on any
bounded complex polydisc of source amplitudes.

Thus the previous full residual bound becomes

    ||f(t)|| <= C [(qa)^-1/2+a^3/2+sqrt(q)a^7/2]
                +C_V q^-1/2+C a^-1/2 exp(-c q),           (8)

up to a bounded modulus of F. The boundary norms remain those of the
joint trial, multiplied by that same bounded phase modulus.

## 3. Uniform positive-resolution estimate and analytic differentiation

First extend the old imaginary-frequency estimate to `z=omega+i nu`.
In a throat collar the real part omega is an anti-Hermitian part of the
radial evolution generator. Its Hermitian part has norm at most
`nu+lambda_n beta`, so the lower Gronwall bound used for collar loss
is unchanged. The current identity still gives

    I-K_n(z)^*K_n(z) >= d_n I,
    d_n=nu/(nu+lambda_n beta)
                  [1-exp(-2(nu+lambda_n beta) ell_c)].     (9)

This is a bound for the prescribed spherical load, not the actual
nonspherical throat. For `beta=C_beta/p`, `ell_c=c_ell p` and `nu<=1`,
the exponential bracket is bounded below uniformly for large q, giving
`d_n^-1<=C[1+lambda_n/(nu p)]`. Summing the same complete angular
boundary estimates as before now gives a bound `C nu^-1 q^-3/4`
at p in (1). The bulk cost from (8) is also at most
`C nu^-1 q^-3/4`. Ground trace errors are O(q p^4)=O(q^-2).
The bounded real frequency only changes constants in the trial.

For complex source amplitudes let `M_j=||V_j||` (discard identically zero
directions). For one or two active variables choose

    |t_j|<=r_j=nu/(8 M_j).                                (10)

Then `||Im sum t_j V_j||<=nu/4`. The exterior current estimate has
bulk coercivity at least `3nu/4`; the throat loss (9) is unchanged since
the source vanishes there. The phase modulus in (5) is at most
`exp(nu L_p/4)`, uniformly bounded on a fixed trimmed contour family.
All coefficient bounds in (8) are uniform on (10). Consequently

    sup_polydisc ||R(t)-R_W(t)|| <= C q^-3/8 nu^-1/2.       (11)

The assumed zero-incoming resolvent and a bounded-perturbation Neumann
series construct a holomorphic solution throughout a slightly larger
polydisc, including its port trace through the current estimate.
Cauchy's integral formula applied to this operator-valued difference
gives (3). For distinct variables the extra factors are `1/(r_j r_k)`;
for a repeated derivative they are `2/r_j^2`. Axisymmetric sectors are
orthogonal, so taking their supremum loses no mode-count factor.

At the working clock scale `delta=d/ell=q^-1/16` and physical damping
resolution `eta_phys=delta/ell`, the dimensionless quantity in this
note is **nu=delta^2=q^-1/8**, not delta. Thus (3) gives respectively

    O(q^-5/16), O(q^-3/16), O(q^-1/16).                   (12)

These orders concern derivatives in fixed dimensionless profiles `d V_phys`.
Changing the units or scaling of the source requires rescaling those
derivatives. Most importantly, the constants still require fixed u_*.
The resonant labels `u~delta/n` are not included by (12).

## 4. Gauge test and the unresolved physical terms

For a smooth axisymmetric pure gauge chi constant near each mouth,
`V_x=chi_y`, `V_z=-chi_s/h`, `V_y=V_0=0`. Thus `v=-chi_s(s,0)`
and (2) transforms by `exp[i(chi_out-chi_in)]`, exactly as the port
gauge law requires. Locally the corrector `i w y=i chi_y(s,0)y`
is the first normal term in `exp[i(chi(s,y)-chi(s,0))]`.
Multiplying the whole ungauged trial by `exp(i chi)` is an exact
alternative for a pure gauge and directly proves residual covariance.
Its difference from (5) is O(q^-1) in sectional norm. This tests the
signs of the phase and transverse correction independently.

This theorem does not cover sources with azimuthal mixing, derivatives
growing with the magnetic scale, source support at a moving physical
mouth, or the shrinking long-return sector. It does not determine the
exterior determinant required by [the regulated mouth action](gauged-mouth-interaction.md).
There is no frequency/angular trace estimate or fixed physical
renormalization prescription in (3). These restrictions preclude
claiming a complete vacuum-stress or physical damping result.

An exact transverse zero mode removes the nonlinear transverse-mass
residual, but does not by itself solve the long-path problem: its width
derivative and the longitudinal transport solvability condition still
vary near the rescaled endpoints. No uniform endpoint theorem is inferred
from that replacement. The next propagation calculation must track the
shrinking labels with `Q=q u` and match both endpoints, while extending
the source class and retaining the exterior determinant.

## 5. Checks and source review

[The script](../checks/dipole_source_matching.py) evaluates (7) in the
actual dipole geometry, compares it with direct Pauli multiplication and
with finite differences of the full trial, and checks the residual's
large-flux order at two geometry resolutions. Nonzero scalar, radial,
azimuthal and axial source components are retained, including complex
amplitudes. A phase-only trial is compared to (5). Separate gauge checks,
complex-frequency collar loss checks and exact scale exponents test
other parts of the argument. These are supplementary calculations, not
a solution of the full loaded exterior PDE. The [receipt](../evidence/dipole-source-matching.json)
pins the source files and prior joint-matching dependencies.

[Bal et al., 2106.00729v1](https://arxiv.org/abs/2106.00729v1) is an
established predecessor for corrected Dirac packets at curved interfaces.
The kernel and complement analysis in Lemmas 3.2-3.3 was read, including
the original PDF's page 18. Two displayed bookkeeping issues are recorded
in the [source audit](../evidence/dipole-source-review.json): the displayed
inverse after PDF (3.19) has the opposite overall sign, and the sequence
space's additional condition at n=0 conflicts with its shifted Hermite
transform. Negating the displayed matrix and allowing all rapidly
decreasing shifted sequences fixes this complement calculation; explicit
checks are included. The note above derives its own source cancellation
and does not import either displayed error. This observation is not a
new physical result or a claim that the paper's main theorem fails.

Remaining to-do list: prove uniform source-dependent long-return endpoint
matching and the full angular response, match the exterior determinant
and physical renormalization, control static/dynamical stress and the
coupled geometry/state, finish novelty and the consequential result,
then independently verify, GCP-certify and prepare the manuscript.
