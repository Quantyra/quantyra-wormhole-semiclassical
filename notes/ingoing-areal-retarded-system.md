# Ingoing areal coordinates and the constrained retarded system

2026-10-21, E003/S168. **The coordinate clock introduced by the optical
gauge is a source identification, not an additional dynamical unknown.**
An ingoing areal gauge represents the exterior perturbation with the
prepared boundary time fixed. It requires one incoming transport and
does not impose the optical outgoing moment condition. The full spherical
equations then take the three mass/redshift response rows below. This
is a gauge and equation-reduction result; the quantum kernel and physical
metric boundary domain remain to be solved.

This corrects the interpretation of the outstanding S168 "physical clock
equation" where it meant this coordinate clock. It does not erase physical
relative mouth times, a controller degree of freedom, or boundary equations
from a separately specified boundary action. The retained paired source,
HH preparation, magnetic bundle, all angular modes, local action and
renormalization scheme stay fixed. Transform the state, frame, source and
initial data together, as in the [transport identity](transported-clock-dirac-response.md).

## 1. Gauge coverage with the prepared time fixed

On one static exterior, use the background conventions

    g=F(-dt^2+dy^2)+r^2 dOmega^2, r_y=-F,
    v=t-y, g=-F dv^2+2 dv dr+r^2 dOmega^2.

Write h_tt=-2F A, h_yy=2F B, h_ty=F C and delta r=r p=rho.
Take h_new=h-Lie_xi g. Areal gauge fixes

    xi^r=rho, b=xi^y=-rho/F, kappa=xi^v=xi^t-xi^y,
    xi^t=kappa+b.

The remaining condition h_new,rr=0 is exactly

    (partial_t+partial_y)kappa=A-B-C=:f_-.                 (1)

For example the old-coordinate rr component is 2(-A+B+C)/F,
and (Lie_xi g)_rr=2 partial_r kappa, which proves (1) because
partial_r|v=-(partial_t+partial_y)/F. There is no outgoing transport.
For bounded outer p, b=O(y) and therefore b(t,0)=0. If the prepared
source clock in the original chart is T=t+c_b(t), choose
kappa(t,0)=c_b(t). Then delta T_new=c_b-xi^t(t,0)=0. If the original
chart already uses the prepared time, c_b=0.

With initial identification kappa_i(y), the unique retarded solution is

    y<t: kappa=c_b(t-y)+integral_0^y f_-(t-y+s,s) ds,
    y>=t: kappa=kappa_i(y-t)+integral_0^t f_-(s,y-t+s) ds.  (2)

The corner data must agree, with higher corner compatibility when smooth
solutions are required. Every source time in (2) is at most t. The initial
identification is retained; setting it to zero is a choice of data, not an
identity. If an initial hypersurface is to be fixed, its time displacement
also constrains this identification. All induced initial-state contacts
must accompany that choice.

For Re gamma=a>0 the one-sided Laplace solution is

    kappa_hat=e^(-gamma y)c_b_hat
       +integral_0^y e^[-gamma(y-s)](f_-hat+kappa_i)(s) ds,
    ||kappa_hat||_2 <= |c_b_hat|/sqrt(2a)
                           +||f_-hat+kappa_i||_2/a.       (3)

Young's inequality gives the bound and (1) controls its H1 graph norm.
The theorem applies when the displayed input is in L2. It does not
identify optical H1 with the full physical horizon domain. In particular
b=-rho/F need not be optical H1 even when xi is a smooth vector in (v,r).

For a boundary-reaching ingoing characteristic the same solution is

    kappa(v,r)=c_b(v)-1/2 integral_r^infinity h_rr(v,R) dR. (4)

For a characteristic meeting the initial surface, replace the boundary
anchor by that initial value and integrate from its intersection. Thus
(4) is not an instruction to invent a pre-initial boundary history.
If h_rr and rho are smooth through a simple future horizon and the
outer tail in (4) is integrable, xi is smooth there too. No factor 1/F
appears in (4), and no exception at gamma=kappa_h is needed. The assertion
is conditional on these regularity and anchoring hypotheses, not on all
optical L2 probes. It is an exterior linear gauge statement, not a global
nonlinear foliation through a throat where the areal coordinate can fail.

A residual gauge has xi^r=0 and partial_r xi^v=0. Its boundary and initial
identifications fix the remaining function of v. For a full prepared
history, its boundary value alone fixes that function. Apply the construction
on the paired exteriors with the same prescribed time identification;
there is no independently adjustable relative clock or change of pairing.

The compact probe A=q(t)k(y), B=C=p=0 used in the
[optical coverage result](global-optical-gauge-and-boundary-clock.md)
has a nonzero outgoing moment Theta, but (2) still exists with c_b=0.
Hence the optical codimension-one condition is not an obstruction to this
gauge. Neither this example nor the gauge theorem is an on-shell mode.

## 2. The boundary condition and the two metric amplitudes

Use the exact areal family

    ds^2=-e^(2psi) Fcal dv^2+2e^psi dv dr+r^2 dOmega^2.

For the fixed vacuum curvature Fcal/r^2 -> lambda and conformal frame
Omega=B_ref/r, its boundary metric is

    gbar=-lambda B_ref^2 e^(2psi_infinity)dv^2+B_ref^2dOmega^2.

Matching the prepared metric in T gives e^psi_infinity=dT/dv.
Thus psi_infinity=c_b' to first order in a chart T=v+c_b, and
psi_infinity=0 in the prepared chart v=T. This condition fixes the lapse
source. It does not select all higher-curvature boundary traces or a fast
scalar self-adjoint phase. The fixed source is consistent with the earlier
[centered-source normalization](stationary-centered-exterior-source.md).
We do not vary this coordinate label as a new physical oscillator.

Put delta psi=zeta and delta Fcal=f at fixed areal r. Direct subtraction
of the Lie derivative gives

    zeta=A+B-partial_r rho-partial_v kappa,
    f=-2FB+2F partial_r rho-F_r rho+2 partial_v rho
     =-2FB-2 partial_y rho-F_r rho.                      (5)

The derivatives in the first line are at fixed v; the last y derivative
is at fixed t. The second amplitude equals
delta[(nabla r)^2]-F_r delta r, the spherical invariant on the static seed.
The new components are h_vv=-f-2Fzeta, h_vr=zeta, h_rr=0 and h_AB=0.

## 3. Exact constraints and their retarded linearization

Let Q^2=pi G q^2/g_c^2, Lambda_0=-3lambda_0, and T^Q be the complete
renormalized four-dimensional stress after moving the retained local
curvature action to the right-hand side. This move is bookkeeping only:
T^Q still depends on the full metric and its higher derivatives. It is
not a prescribed load and no order reduction is made. The fixed magnetic
stress has T^v_v=-Q^2/(8pi G r^4), T_rr=T^r_v=0. The retained
[neutral spherical sector](neutral-spherical-current-decoupling.md)
allows zero electric variation; independent electric sectors remain open.

Direct four-dimensional curvature calculation yields, for arbitrary psi
and Fcal,

    G_rr=2psi_r/r,
    G^v_v=(r Fcal_r+Fcal-1)/r^2,
    G^r_v=-Fcal_v/r.                                    (6)

Define the charged AdS mass aspect

    M=r[1+Q^2/r^2+lambda_0 r^2-Fcal]/(2G).

The three equations calE=G-3lambda_0 g-8pi G(T_EM+T^Q)=0 are exactly

    psi_r=4pi G r T^Q_rr,
    M_r=-4pi r^2 (T^Q)^v_v,
    M_v= 4pi r^2 (T^Q)^r_v.                             (7)

These specialize the established spherical mass balance; the magnetic
energy is absorbed once. The older [mass-balance note](sourced-spherical-mass-charge-balance.md)
uses sphere-integrated stress, 4pi r^2 times the stress used here.
M is not automatically the finite renormalized boundary energy: if the
physical lambda differs from lambda_0, even the seed M has an r^3 term.

At psi_0=0 set delta M=mu, so f=-2Gmu/r. The actual retarded Jacobian has
the three rows

    zeta_r-4pi G r delta T^Q_rr[mu,zeta]=j_1,
    mu_r+4pi r^2 delta (T^Q)^v_v[mu,zeta]=j_2,
    mu_v-4pi r^2 delta (T^Q)^r_v[mu,zeta]=j_3.            (8)

Each stress variation includes its inverse-metric, state, domain, second
vertex, local and initial contacts. Independent initial-state data add
their source to (8); mixed tensor variations include variation of index
raising. The notation does not assert that the kernels have been evaluated.
For a Newton correction about an off-shell seed with residual calE_0,

    j_1=-r calE_0,rr/2,
    j_2= r^2 calE_0^v_v/(2G),
    j_3=-r^2 calE_0^r_v/(2G).                            (9)

A homogeneous mode test represents physical stability only on an exact
background (or with a proved transfer from an approximate one).

At an exact solution the three rows and the full Ward identity imply the
angular equation. Symmetry gives calE^r_r=calE^v_v+Fcal calE_rr. Once the
orbit block vanishes, the radial divergence is -2calE^theta_theta/r,
without division by Fcal. Off shell, the corresponding linear relation is

    delta calE^r_r=delta calE^v_v+F delta calE_rr
                                      +f calE_0,rr.     (10)

The linearized Ward identity also contains variation of the connection
acting on calE_0. Discarding these residual terms would falsely promote
three homogeneous seed rows to the full physical mode equation.

Only if the outer integral converges may the homogeneous first row be
written zeta(v,r)=-4pi G integral_r^infinity R delta T_rr(v,R)dR,
using zeta_infinity=0. Otherwise retain a finite outer anchor and the
renormalized boundary trace problem. Likewise at a finite R,

    mu(v,r)=mu(v,R)+4pi integral_r^R s^2 delta T^v_v(v,s)ds
                                      -integral_r^R j_2(v,s)ds.

The mass anchor and initial flux data are not discarded. Cross-derivative
compatibility of the two mass rows is part of the actual response problem.

## 4. Connection to the retained Dirac response

On any fixed compact exterior region F>0, use the background transformation
dv=dt-dy, dr=-Fdy. With k=f/F, (5) becomes the ADM tangent

    A=zeta+k/2, B=-k/2, C=zeta+k, p=0,
    d=A-B=C=zeta+k, w=A-p=zeta+k/2.

The [full ADM Hamiltonian](adm-dirac-boundary-response.md) therefore gives

    B_i=i(I-alpha)(d_i partial_y+d_i'/2)+beta m_nu w_i.   (11)

The kinetic and shift terms remain even for the angular zero mode.
To fix all second contacts, use affine psi and k coordinates in the exact
family, E=e^psi and a=2-E(1+k). Direct lapse/shift extraction gives

    L^2=F E a, N^2=F E/a, speed=N/L=1/a,
    S=speed-1, M_nu/m_nu=sqrt(E/a).                      (12)

For two tangents the mixed second derivatives are

    speed_ij=3zeta_i zeta_j+3(zeta_i k_j+zeta_j k_i)+2k_i k_j,
    mass_ij/m_nu=2zeta_i zeta_j
                 +3(zeta_i k_j+zeta_j k_i)/2+3k_i k_j/4.

Hence W_ij=i(I-alpha)(speed_ij partial_y+speed_ij'/2)
+beta m_nu(mass_ij/m_nu). At the fixed boundary source zeta_infinity=0
and k_infinity=0, S_infinity=0 and the paired projector is the retained
Pi_0. This does not discard any bulk/state/observable/local contact.

There is a material horizon limitation. Smooth nonzero f at a shifted
horizon can make k=f/F diverge. The old static slices then need not remain
spacelike: (12) requires a>0. Thus (11)-(12) are compact-region response
vertices, not a global bounded static-slice perturbation theorem. The
stress must be matched covariantly through the horizon with the prepared
state. No null-slice self-adjoint Hamiltonian or unitary map from one
static exterior to a truncated horizon-crossing exterior is asserted.

## 5. Source review, verification and decision

[Maedler and Winicour, sections 2-3, equation (28)](https://arxiv.org/pdf/1609.01731v3)
give the characteristic areal construction and radial lapse equation.
Their outgoing convention has psi=2beta after the explicit ingoing
conversion; in spherical symmetry (28) gives the first row of (7).
[Hayward, section II, equation (2.7)](https://arxiv.org/pdf/gr-qc/9710089)
gives the spherical energy balance. These established identities are not
claimed as new, and neither source evaluates this charged quantum kernel.
The present result identifies their use with the fixed S168 preparation,
the exact optical-gauge correction, and the retained full response.

The ordinary argument precedes the frozen focused checker. Verification
derives the four-dimensional Einstein tensor from Christoffels, checks
every displayed gauge/metric/contact coefficient componentwise, and uses
independent characteristic quadrature, complex-frequency ODE refinement
and finite metric differences. Negative controls detect lost initial data,
off-shell residuals and the horizon failure of the old static slices.
These checks do not replace the analytic domain hypotheses or evaluate
the actual response. No Lean campaign is selected.

The next calculation is the renormalized retarded stress in (8), with its
physical outer traces, horizon regularity and initial/Ward compatibility.
There is no additional force equation for the coordinate label c_b.
Physical boundary degrees of freedom, if present in the retained physical
model, remain separate. S168.3-S168.4 and all S167 tasks stay unchecked.

Remaining to-do list: evaluate the coupled charged retarded response and
its physical domain; decide global stability; establish nonlinear existence
and compatibility; complete exterior/state, charge/return, opening,
all-sector stability, transport, preparation, separate creation and audit.
