# Spherical metric vertices and the moving Dirac boundary domain

2026-10-21, E003/S168. The complete canonical spherical Dirac Hamiltonian,
its first and mixed second metric vertices, and the geometric continuation
of the retained zero-source projector are explicit below. A fixed-domain
representation has a time-derivative source vertex. Its retarded readout
therefore differs from its driving vertex. These formulas provide the
response to be evaluated; they do not evaluate the charged boundary kernel,
select a physical metric boundary condition, or prove stability.

Keep the prepared HH state, fixed magnetic bundle, both angular eigenvalue
signs and zero modes, and common four-dimensional subtraction. The prior
[transported-clock identity](transported-clock-dirac-response.md) and
[curvature reconstruction](boundary-clock-curvature-reconstruction.md)
remain in force. The local Weyl, gravitational and boundary actions are
additional fixed parts of the total response, not removed by this reduction.

## 1. General spherical Hamiltonian

Write the Lorentzian metric and canonical field as

    ds^2=-N^2 dt^2+L^2(dy+S dt)^2+r^2 dOmega^2,
    phi=r sqrt(L) Psi, v=N/L, M_nu=lambda_nu N/r.

Here N,L,r are positive. Angular normalization is the retained one.
The warped-product rescaling chi=r Psi removes the radius connection.
For e^0=N dt, e^1=L(dy+S dt), Cartan's equations give

    omega^0_1=omega^1_0=a e^0+k e^1,
    a=N_y/(NL), k=[L_t-(LS)_y]/(NL).

Choose gamma^0=beta alpha, gamma^1=beta, with Hermitian involutions
alpha,beta anticommuting. Then gamma^0 gamma^1=-alpha. In this radial
Clifford frame the angular term is a pseudoscalar: the reduced equation
is (gamma^a nabla_a-i alpha lambda_nu/r)chi=0. This keeps gamma^1=beta
consistent with the retained radial projector. Its evolution is

    chi_t=(S-v alpha)chi_y
          -N_y alpha chi/(2L)-[L_t-(LS)_y]chi/(2L)-i beta M_nu chi.

Substituting phi=sqrt(L)chi cancels L_t and gives exactly

    H=-i alpha(v partial_y+v_y/2)+beta M_nu
                         +i(S partial_y+S_y/2) I.                 (1)

The canonical norm is integral phi^dagger phi dy. With
K=v alpha-S I, the polarized current is phi^dagger K psi, so (1)
has the Green form determined by -i K at an oriented endpoint.
The sign of the endpoint convention does not affect maximal isotropy.
No N_t or r_t potential remains. This cancellation is an exact field
normalization, not an adiabatic approximation.

At the unchanged static seed N_0=L_0=sqrt(F), set

    N=N_0 exp(a), L=L_0 exp(b), r=r_0 exp(p), S=s,
    d=a-b, w=a-p, m_nu=lambda_nu sqrt(F)/r_0.

For tangent profiles i,j the differential-expression vertices are

    B_i=-i alpha(d_i partial_y+d_i'/2)+beta m_nu w_i
                                  +i(s_i partial_y+s_i'/2) I,
    W_ij=-i alpha(d_i d_j partial_y+(d_i d_j)'/2)
                                  +beta m_nu w_i w_j.             (2)

The shift chart is affine; its mixed second derivative is zero in (2).
A common conformal tangent a=b=p has zero canonical vertex. The retained
local Weyl response still contributes. Optical a=b and s=0 recover the
old mass vertex and its second contact. At lambda_nu=0 the kinetic and
shift terms still contribute. Thus an LLL mass-only treatment misses
metric response outside optical gauge.

## 2. The projector must follow the boundary normal

For the folded paired collars use

    A=tau_3 tensor sigma_2, beta_f=I tensor sigma_1,
    chi_0=tau_3 tensor sigma_1, Pi_0=(I+chi_0)/2.

An angular identity is understood. The following formulas assume matched
collar metric functions on the paired copies, as in the clock calculation;
at zero source they also apply separately to each decoupled exterior.
Keep the same geometric radial projector in the boundary rest frame,
the same spin/bundle identification and zero coupling angle. The shift
does not introduce an independent boundary coupling.

At y=0 put u=S/v and require |u|<1, exactly the condition that the
boundary is timelike. In an ADM slice frame its allowed orthogonal
projector is

    Pi_u=[I+sqrt(1-u^2) chi_0+u A]/2.                            (3)

To see its geometric origin, set eta=atanh(u), D=exp(eta A/2).
The boundary rest-frame boost satisfies

    D^dagger(v A-S I)D=sqrt(v^2-S^2) A.

The range of D Pi_0 is the range of (3). Also Pi_u^2=Pi_u,
Pi_u^dagger=Pi_u, rank Pi_u=2 per angular component, and

    Pi_u A Pi_u=u Pi_u,
    Pi_u(v A-S I)Pi_u=0.                                      (4)

This establishes the maximal isotropy of the geometrically transported
domain. Flux conservation alone would allow other boundary couplings;
it is the retained rest-frame projector that fixes (3). Keeping Pi_0
in the ADM frame instead gives Pi_0(v A-S I)Pi_0=-S Pi_0 and fails
when the boundary shift is nonzero. At conformal infinity these equations
refer to canonical data in the same conformal frame. A finite regulator
does not become a new physical reflecting wall.

## 3. All contacts on a fixed domain

Let J=A chi_0, so J^dagger=-J, J^2=-I. Choose a smooth collar function
rho equal to one at the boundary and zero toward the unchanged horizon.
Extend only the boundary angle, not the entire metric, by

    theta(t,y)=rho(y) arcsin(S(t,0)/v(t,0)),
    R=exp(theta J/2), phi=R zeta.

Then R Pi_0 R^dagger=Pi_u at the boundary. The exact fixed-domain
Hamiltonian and its principal boundary coefficient are

    Hhat=R^dagger H R-i R^dagger R_t,
    Khat=R^dagger(v A-S I)R
         =sqrt(v^2-S^2) A+S(chi_0-I)             at y=0.        (5)

Thus Pi_0 Khat Pi_0=0 for the entire family. Smooth changes of rho give
a unitary identification equal to the identity at the boundary and
near the horizon. With the state, observables and local terms transformed,
they change representation, not the boundary drive. No global physical
metric inverse or new horizon-domain theorem is inferred from this fact.

For static spatial tangent profiles with amplitudes lambda_i(t), define

    X_i=rho s_i(0) J/2,
    X_ij=-rho [s_i(0)d_j(0)+s_j(0)d_i(0)] J/2.

These are the first and mixed second derivatives of log R. All X commute.
The expansion through second order is

    Hhat=H_0+sum_i(lambda_i A_i+dot(lambda_i) D_i)
       +1/2 sum_ij lambda_i lambda_j T_ij
       +sum_ij lambda_j dot(lambda_i) E_ij+O(lambda^3),
    A_i=B_i+[H_0,X_i], D_i=-i X_i, E_ij=-i X_ij,
    T_ij=W_ij+[B_i,X_j]+[B_j,X_i]
       +1/2([[H_0,X_i],X_j]+[[H_0,X_j],X_i])+[H_0,X_ij].       (6)

Here the remainder counts amplitudes and their time derivatives together.
For arbitrary time-dependent tangent directions the mixed second vertex
is T_ij-i partial_t X_ij, with the time dependence also in B_i,X_i.
No extra connection commutators occur because these X commute.
The boundary principal coefficients provide a useful independent check:

    delta_i Khat=d_i A+s_i(chi_0-I),
    delta_i delta_j Khat=(d_i d_j-s_i s_j) A.                  (7)

Omitting X_ij leaves an erroneous mixed d,s boundary term. Omitting the
projector transport altogether already violates the first flux identity.
Although (6) contains commutators, X_i generally does not preserve the
old boundary domain. B_i and [H_0,X_i] must be combined before treating
them as a symmetric fixed-domain vertex. One cannot discard either part
as an independent pure gauge perturbation.

## 4. Retarded driving and readout are different

The derivative coupling in (6) matters when differentiating the action.
Define the canonical generalized force, with this fixed sign convention,

    F_i=Tr(C partial Hhat/partial lambda_i)
                      -partial_t Tr(C partial Hhat/partial dot(lambda_i)).

For a stationary prepared C=f_beta(H_0), the first variation is

    delta F_i=Tr(A_i delta C)-Tr(D_i delta C_t)
                             +sum_j lambda_j Tr(C T_ij).       (8)

The E_ij dot(lambda_j) contacts cancel under the time integration by
parts, since E_ij=E_ji and C_t=0 at the background. They have not been
dropped from the Hamiltonian. Initial/final source endpoints must also
be retained; take sources zero near the initial surface here.

Let L_gamma=gamma+i ad(H_0), Re gamma>0. For zero initial covariance
variation, the canonical retarded response to lambda_j is

    delta C_j=-i L_gamma^-1 [A_j+gamma D_j,C],
    K_ij(gamma)=Tr(C T_ij)
       -i Tr[(A_i-gamma D_i)L_gamma^-1[A_j+gamma D_j,C]].       (9)

There is no complex conjugation on gamma in this bilinear Laplace
formula. Sesquilinear versions require the corresponding adjoint test
frequency. Using A_i+gamma D_i as the readout is incorrect. For an
independent initial covariance c_0 add

    Tr[(A_i-gamma D_i)L_gamma^-1 c_0]+Tr(D_i c_0)              (10)

to the Laplace force. This keeps the initial derivative contact.
There is no instantaneous replacement of C by f_beta(Hhat).

Equations (8)-(10) are exact at a common finite regulator and as retarded
distribution/form identities wherever the corresponding renormalized
response is defined. At continuum level sum all angular channels and
perform the same four-dimensional subtraction, adding the retained local
metric, Weyl and boundary variations. Separate displayed traces need not
be finite. These formulas do not prove a new trace-class estimate, remove
boundary divergences, or establish endpoint bounds for the genuine drive.

For a finite matrix or a domain-preserving R, (9) agrees exactly with
the unrotated response using B_i,W_ij; this follows by cyclicity and
L_gamma C=gamma C. This is a useful check on both frequency signs and
all second contacts. It cannot justify cancelling the actual boundary
domain contribution: X_i mixes the boundary eigenspaces there. Only a
change of collar extension, whose relative map preserves the domain,
is covered by that representation-independence argument.

## 5. Evidence, sources and next decision

The accompanying checker independently derives (1) from the Cartan
connection, checks its metric Taylor coefficients and the exact folded
projector identities component by component, and tests (6)-(10) against
finite matrix evolution and exact unitary conjugation. Negative controls
remove the projector transport, X_ij, the time connection, the readout
sign, and the initial contact. These are mathematical identity fixtures,
not a numerical charged Dirac spectrum or physical stability calculation.

[Maldacena, Milekhin and Popov, Appendix C](https://arxiv.org/pdf/1807.04726v3)
provides the established spherical spinor rescaling and conformal special
case. [Bintanja, Freivogel and Rolph, (2.8)-(2.12)](https://arxiv.org/html/2308.00871v2)
provides the retained radial projector and coupled boundary action.
The ADM reduction, boost, fixed-domain contacts and readout formulas above
are derived here. Neither source supplies the charged physical response
or a global stability theorem. The independent stationary framework of
[Finster and Roken](https://arxiv.org/pdf/1512.00761v3) uses additional
geometric hypotheses; it is not a time-dependent theorem for this model.

Next evaluate the genuine boundary column with (9) and the fixed local
terms, obtain the physical clock compatibility, and control the global
coupled problem. The preceding curvature-free boundary column is not
declared pure gauge. Keep all off-shell residuals and independent initial
data. S168.3-S168.4, S167 and S077 remain open.

Remaining to-do list: evaluate the renormalized boundary response and
clock compatibility; global coupled charged stability; nonlinear existence
and compatibility; exterior/state, charge/return, opening, all-sector
stability, transport, preparation, separate creation and final audit.
