# The transported clock's Dirac response and off-shell contacts

2026-10-21, E003/S168. The coordinate part of the boundary-clock
response is fixed by covariance: the prepared Dirac state is transported,
not rethermalized. Its Hamiltonian vertex, state variation and observable
contacts cancel the apparent retarded spectral denominators. The full
equation tensor transforms by its Lie derivative, including the seed's
nonzero residual. This determines the coordinated clock/frame response;
it does not evaluate a genuine boundary drive or the full shape operator.

Use the unchanged vectorlike Dirac theory, magnetic bundle, full angular
tower, common subtraction, prepared HH state and zero-source paired
projector. The [boundary reconstruction](boundary-clock-curvature-reconstruction.md)
remains in force. Its curvature-free boundary column is NOT presumed to
be a pure diffeomorphism, and its physical response remains to be solved.

## 1. Coordinate transport on the actual one-particle domain

In a canonical radial representation on each exterior write

    h=-i alpha partial_y+beta m(y),
    alpha^2=beta^2=1, {alpha,beta}=0,
    m=lambda_nu sqrt(F)/r.                            (1)

The sign of alpha incorporates the chosen inward coordinate and folded
mouth orientation. Keep both angular signs and all zero modes. The
allowed zero-source boundary projector is the previously fixed Pi_0;
its maximal isotropy gives Pi_0 alpha_normal Pi_0=0. No replacement
boundary condition or angular truncation is used in the theorem.

Consider an infinitesimal pullback with positive sign

    T=t+epsilon a(t,y), Y=y+epsilon b(t,y),
    xi=a partial_t+b partial_y,
    delta g=Lie_xi g.                                 (2)

Take a,b smooth, identical on the paired copies, zero near the initial
surface and at the horizon, with b(t,0)=0. Transport the boundary clock,
frame, regulator surface, projector and state along with the metric.
For explicit clock/frame extensions one can use

    a=c(t)chi(y), b=partial_t c(t) y chi(y),            (3)

where chi is one near zero and zero in the interior past a finite collar.
Then a_y(t,0)=b_t(t,0)=0; the scalar boundary normalization commutes with
Pi_0 and its constant algebraic matrix is unchanged in this frame. The
time parametrization and radial defining function do transform. This is
a common coordinate change, not a change of the boundary coupling angle
or an independently selected state perturbation. General spin frames
require the corresponding transported spin lift of Pi_0.

The canonical conserved current is j^t=psi^dagger psi,
j^y=psi^dagger alpha psi. On the pulled-back t slice its norm density is

    psi(T,Y)^dagger [Y_y I-T_y alpha] psi(T,Y).

For small epsilon the slices remain spacelike and this matrix is positive.
The canonical identification therefore includes the square root of this
matrix, not just evaluation at (T,Y). Differentiating and using
i partial_t psi=h psi gives

    delta psi=-i G_xi psi,
    G_xi=Q_a-P_b,
    Q_a={a,h}/2=-i alpha(a partial_y+a_y/2)+beta a m,
    P_b=-i(b partial_y+b_y/2).                         (4)

In particular the terms a_y/2 and b_y/2 are required by the current norm.
They are not optional ordering conventions.

This identification is an isometry on the actual Cauchy solution space:
integrate the polarized conserved current between the two slices. The
outer flux vanishes by the transported maximal-isotropic projector; the
slices coincide near the horizon. Solving backwards between the slices
gives the inverse. The existing self-adjoint Dirac evolution supplies
these smooth spacelike Cauchy solutions. Thus the geometric map is unitary
on one-particle data. Formula (4) is its derivative on the smooth common
core; it does not require treating G_xi as a bounded operator.

At the boundary, smooth solutions and all their time derivatives obey
the same projector. For (3), the infinitesimal half-density factor there
is scalar and the shifted time evaluation preserves the domain. On the
full angular Hilbert space the current isometry is the direct sum of
the channel maps. The CAR state can be pulled back by this isometry
without asserting existence of a unitary implementer in a chosen
infinite-volume Fock representation.

## 2. Exact canonical clock vertex

Let U_epsilon(t)=I-i epsilon G_xi(t)+O(epsilon^2) be that Cauchy-data
identification. Its transformed Hamiltonian is
U_epsilon h U_epsilon^dagger+i partial_t U_epsilon U_epsilon^dagger.
Consequently

    V_xi=partial_t G_xi-i[G_xi,h].                    (5)

Direct differential-operator algebra gives the equivalent local formula

    V_xi=-i alpha[(a_t-b_y)partial_y+(a_ty-b_yy)/2]
          +beta(a_t m+b m_y)
          +i[(b_t-a_y)partial_y+(b_ty-a_yy)/2]I.       (6)

It agrees with the lapse, radial-metric and shift changes of (2).
For b=0 and a=c(t)chi(y), (5) becomes

    V=partial_t c Q_chi-i c[Q_chi,h],
    [Q_chi,h]=chi_y partial_y+chi_yy/2.               (7)

The mass terms cancel in this commutator because [Q_chi,h]=[chi,h^2]/2.
Dropping the shift term in (7) changes the physical perturbation.

On any region where the diffeomorphism preserves optical gauge,
a_t=b_y and b_t=a_y. There (6) reduces exactly to the usual mass vertex

    V_xi=beta m v_xi,
    v_xi=a_t+(m_y/m)b,
    p_xi=(r_y/r)b.                                   (8)

For a zero angular eigenvalue, use the nonsingular form
V_xi=beta(a_t m+b m_y)=0 in this optical region. The kinetic/shift
vertices outside that region and all anomaly terms remain. A cutoff
extension of a boundary clock generally has those extra vertices.
Equation (8) is not permission to apply a shape-only kernel to the
entire cutoff diffeomorphism.

## 3. The retarded state is transported, with all initial contacts

Let C_beta=(1+exp(beta_T h))^-1 be the initial one-particle covariance.
It commutes with h. In the canonical identification its coordinate
variation is exactly

    delta C(t)=-i[G_xi(t),C_beta].                    (9)

Indeed the linearized evolution is
i partial_t delta C=[h,delta C]+[V_xi,C_beta].
Substitute (5),(9) and apply the Jacobi identity. With arbitrary initial
covariance perturbation delta C_0 the complete solution is

    delta C(t)=-i[G_xi(t),C_beta]
       +exp(-iht)[delta C_0+i[G_xi(0),C_beta]]exp(iht). (10)

For the zero-initial coordinate change used here the second term vanishes.
It is retained when independent initial state data are present. Equations
(9)-(10) hold as forms on the common core and for smeared CAR observables;
no trace-class thermal density matrix on the infinite exterior is assumed.

For a separated right-half-plane frequency gamma, the generator amplitude
G_gamma may itself contain gamma, as in (3). For h eigenvalues e_a,e_b,
Delta=e_b-e_a, its vertex and covariance response obey

    (V_gamma)_ab=(gamma-i Delta)(G_gamma)_ab,
    (gamma-i Delta)(delta C_gamma)_ab
                    =-i(V_gamma)_ab(f_b-f_a),
    (delta C_gamma)_ab=-i(G_gamma)_ab(f_b-f_a).         (11)

The apparent spectral denominator cancels. Degenerate blocks have
f_a=f_b and no extra thermal term. This cancellation concerns the
coordinated diffeomorphism column, not the inverse gravitational response.
It persists channel by channel before common subtraction. It does not
remove poles belonging to independent state data in (10).

A transported smeared observable also changes:
delta O=-i[G_xi,O]. At finite regulator cyclicity gives

    Tr(delta C O)+Tr(C_beta delta O)=0                 (12)

for the coordinate part. Tensor components and coordinate densities,
rather than a fully transported smearing, instead acquire their ordinary
Lie derivative. Missing the second term in (12) produces a spurious
clock response. Covariant subtraction carries this identity to local
renormalized observables, as discussed below.

For a uniform time reparametrization, a=c(t),b=0, one has G=c h,
V=partial_t c h and delta C=0. In contrast, replacing the state by
f_beta((1+partial_t c)h) would give
delta C_eq=partial_t c h f_beta'(h), generally nonzero. That replacement
is an instantaneous thermal prescription, not the evolved prepared state.
No temperature or thermal occupation change is adopted here.

## 4. Full covariant response, boundary contacts and the off-shell seed

The vectorlike four-dimensional theory admits a covariant conserved
current and stress prescription. The selected local gravitational,
Maxwell, Weyl and boundary terms are transported in their fixed scheme;
no finite counterterm is altered. One must transport the cutoff surface
and boundary geometry as well. Keeping a coordinate cutoff fixed while
changing its physical embedding is a different variation with extra
boundary terms. The Weyl anomaly is retained; it is not a failure of
diffeomorphism covariance.

It follows, including the induced state/domain variation, that

    delta_xi calE_ab=Lie_xi calE_ab.                  (13)

For the stationary mixed residual diag(b_t,b_r,b_o,b_o), direct tensor
transformation gives

    delta calE^t_t=b partial_y b_t,
    delta calE^y_y=b partial_y b_r,
    delta calE^theta_theta=delta calE^phi_phi=b partial_y b_o,
    delta calE^y_t=(b_r-b_t)b_tcoord,
    delta calE^t_y=(b_t-b_r)a_y.                      (14)

Here b_tcoord means partial_t b, to distinguish it from the residual b_t.
In particular the two off-diagonal mixed components are not generally
opposites off shell: their sum is (b_r-b_t)(b_tcoord-a_y).

Let w=F r^2, tau_0=b_t+b_r+2b_o and o_0=b_t+b_r. Since the volume
variation is div xi=a_t+b_y+2(n_y+l_y)b, the two reduced density rows obey

    delta_xi(w tau_0)=partial_t(a w tau_0)+partial_y(b w tau_0),
    delta_xi(w o_0)=partial_t(a w o_0)+partial_y(b w o_0). (15)

The fixed action normalization and sphere factor multiply both sides
unchanged. These are the required coordinate/density contacts. The trace
residual is known locally; the orbit residual still contains the actual
stationary quantum stress. Formula (15) does not claim its numerical value.
Only on an exact solution does the coordinated gauge column vanish.
Neither the seed residual nor its boundary contact may be discarded.

For complete perturbations, the combined response therefore satisfies

    R[h+Lie_xi g, delta boundary+delta_xi boundary,
                       delta state+delta_xi state]
       =R[h,delta boundary,delta state]+Lie_xi calE.  (16)

Use compact smearings first, or include the specified boundary densities
when the smearing reaches the boundary. A boundary-only response is not
asserted zero separately. Equation (16) fixes the coordinated sum and
explains how changes of coordinate extension are compensated.

[Sanders, section 4.3 and Theorem 4.18](https://arxiv.org/pdf/0911.1304v2)
provides the established Dirac relative-Cauchy/stress-commutator framework.
[Zahn, sections 4.1-4.2](https://arxiv.org/pdf/1210.4031v3) provides the
locally covariant conservation framework. These bulk results do not
supply our paired timelike boundary condition. The boundary current
isometry and transported-domain argument above use the actual retained
projector. Nor do those sources compute the charged shape response.

## 5. Consequence for the remaining physical clock column

The coordinate clock/frame response now has no independently adjustable
thermal susceptibility: (5),(10),(15) fix it. The previous reconstruction
column with U=0 and V=chi is more general. Curvature-free is not the same
as pure gauge. A nontrivial global optical diffeomorphism with preserved
endpoint and ordinary horizon behavior is not available at generic large
real frequency, as the earlier gauge theorem shows. Thus neither (11)
nor (13) sets that entire boundary column's shape response to zero.

The physical next step is to separate the coordinated contribution using
(16) and evaluate the remaining genuine metric/boundary response with
the prepared state and domain fixed as required. Its scalar clock value
still has to satisfy the original boundary data and coupled curvature
condition. No new state knob, order reduction, global metric inverse,
charged mode or stable nonlinear wormhole is introduced by this result.

Freeze this note, checker and source review before acceptance. Verify
the differential generators and local vertex, actual projector flux,
matrix retarded evolution with and without independent initial data,
spectral cancellation, observable contacts, and all four-dimensional
off-shell density identities. Independent coordinate pullbacks of exact
Dirac waves test the half-density generator. Finite matrix ODE fixtures
test the response identities, not the actual charged spectrum. Preserve
every failure and exact source; no Lean or Fock implementability claim.

Remaining to-do list: genuine boundary shape response and physical clock
compatibility after coordinate transport, global coupling and charged
stability; S167 signed source, higher response, nonlinear existence and
compatibility; exterior/state, charge/return, opening, all-sector stability,
transport, preparation, separate creation and final audit.
