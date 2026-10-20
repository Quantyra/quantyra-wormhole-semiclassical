# The retained Newton operator beyond the analytic seed

2026-10-20 HST, E003/S162 supporting S139/S074. Ordinary argument and
verification contract written before execution. **Operator identification;
quantitative compatible inverse and full research goal remain open.**

The pure conformal Hessian has the same exact curvature factorization
at every admissible background, without using a background field equation.
The coupled operator still contains the actual nonlocal shape response.
Its physical fast range has negative index. We give an exact projected
formulation that retains the missing signed equation.

This extends the seed block calculation in the
[outgoing construction](stationary-scalar-outgoing-inverse.md), and
uses the already proved [weighted map](stationary-global-weighted-map.md)
and [adjoint quotient](stationary-adjoint-phase.md). It does not re-prove
those global estimates or turn their qualitative bounds into numerical
Newton constants. The [complete candidate defect](comparison-path-complete-equation-residual.md)
remains an equation defect.

## 1. Coordinates, equations and the density contact

Keep the retained action, fixed flux, outer curvature lambda=lambda_e,
clock, positive prepared temperature, subtraction scale and paired
kernel-free Dirac realization. Use the
[fixed disk](stationary-disk-and-normal-inverse.md):

    g=exp(2Phi)[exp(2V)(d xi^2+xi^2 dvarphi^2)+B^2 dOmega^2],
    period(varphi)=2pi, 0<=xi<=1, r=B exp(Phi).

Admissible fixed comparisons are smooth on the open exterior and the
regular disk bolt; at infinity they satisfy the five differentiated
conormal jets, fixed constant/quadratic data, correlated mass cubic and
fast remainder required by the weighted-map theorem. This statement
does not establish a nonlinear neighborhood preserved by a solver.

For affine disk directions h=(S,T), delta log dvol=4S+2T. Write the
moving-volume action gradient as F=(F_Phi,F_V). If E_rho,E_varphi,E_theta
are the full mixed tensor defects, then

    F_Phi=(E_rho+E_varphi+2E_theta)/(8piG),
    F_V=(E_rho+E_varphi)/(8piG),
    E_theta=4piG(F_Phi-F_V).                         (1)

Thus the independent radius equation is present: keeping the physical
orbit metric fixed means (S,T)=(s,-s), and delta(r^2)=2r^2 s.
The corresponding radius Euler derivative per physical volume is
(F_Phi-F_V)/(2r^2). It cannot be dropped in favor of the two constraints
from a prescribed-source Einstein solve. At F=0 the Ward identity and
regular bolt give (a_orb^2 r^2 E_rho)'=0 with zero integration constant,
so (1) recovers the full static spherical tensor equation.

Define the genuine affine Hessian A at the current comparison by
I''(h,k)=integral k dot A h dvol_g. Then

    A h=DF[h]+F(4S+2T).                              (2)

This is a symmetric Hessian, including the boundary completion. DF
alone generally is not symmetric off shell. In particular
DF*=A-l^T F^T, l=(4,2), in the physical pairing. Its adjoint moment
cannot be silently identified with that of A.

For one fixed reference volume dvol_0, put

    Jvol=exp[4(Phi-Phi_0)+2(V-V_0)],
    mathfrakF=Jvol F, D mathfrakF[h]=Jvol A h.          (3)

This is the action gradient in a fixed output identification.
Newton's equation D mathfrakF[h]=-mathfrakF is exactly A h=-F.
The positive volume ratio is bounded above and below for each admitted
comparison pair. The leading mass cubic cancels in log Jvol, since
4S_3+2T_3=0. No smallness or Lipschitz bound for (3) is inferred.

## 2. Exact conformal block at a general background

Use physical four-dimensional curvature, positive T_g=-Box_g, and

    epsilon=G/(60pi), lambda0=lambda+11epsilon lambda^2/6,
    J_ab=Ric_ab-(R/4)g_ab, M_em=q^2/(2r^4),
    B_an=(3/2)C^2-(11/12)Euler4-20M_em+Box R,
    E_tr=R+12lambda0+epsilon B_an, F_Phi=-E_tr/(8piG).

The electromagnetic normalization is that of the
[retained reduction](exact-spherical-em-reduction.md); the comparison
metric's effective charge coefficient is not substituted for q.
For delta g=2s g, delta r=r s, the local variations are

    delta R=-6Box s-2Rs,
    delta C^2=-4s C^2, delta M_em=-4s M_em,
    delta Euler4=-4s Euler4+8J:Hess s-2R Box s,
    delta(Box R)=-6Box^2 s-2R Box s-2grad R.grad s-4s Box R.

The last identity includes variation of Box itself. Combining them,

    delta E_tr=-6epsilon Box^2 s-(6+epsilon R/6)Box s
       -(22epsilon/3)J:Hess s-2epsilon grad R.grad s
       +(-2R-4epsilon B_an)s.

The density term 4s E_tr cancels the undifferentiated anomaly:

    A_PhiPhi s=(8piG)^-1 [
      6epsilon Box^2 s+div(K grad s)-(2R+48lambda0)s ],
    K=(6+epsilon R/6)g+(22epsilon/3)J,
    div K=2epsilon grad R.                           (4)

Here div J=grad R/4 by contracted Bianchi. Formula (4) is formally
symmetric on compact smooth directions, and the accepted boundary
realization supplies its extension. No background Einstein or
semiclassical equation was used.

Set Z=R+12lambda and E_sp=1/epsilon+11lambda/3. Since
80pi^2/(8piG)=1/(6epsilon), (4) is exactly

    80pi^2 A_PhiPhi =
      (T_g+4lambda)(T_g-E_sp)
      +div[(Zg/36+11J/9)grad]-Z/(3epsilon).           (5)

This is a local scalar BLOCK, not the coupled scalar Schur complement.
If the actual shape block D has an appropriate inverse, that complement
is A_PhiPhi-Bmix D^-1 Bmix*. Positivity, invertibility or norm bounds
proved for the analytic seed's TT block do not automatically hold at
another comparison. Nor is the Euclidean factor T_g-E_sp a Lorentzian
instability theorem.

As a bookkeeping check, for
E_tr,ell=R+(1-ell)U(r)+12ell lambda0+ell epsilon B_an,
replace epsilon by ell epsilon in K and the highest derivative, and
replace the potential in (4) by

    -2R-(1-ell)(r U'+4U)-48ell lambda0.               (6)

At ell=0, R=-U this gives
-(6/(8piG))(T_g+V0), V0=(2U+rU')/6, the earlier seed
potential-theory Hessian. Bookkeeping ell is not the outer lambda.

## 3. Mixed row and the actual shape response

An orbit-only direction v has h_ab=2v P_ab, where P projects onto
the physical orbit plane, delta r=0 and delta log dvol=2v. Thus

    A_PhiV v=-(delta_v E_tr+2v E_tr)/(8piG),
    A_VPhi=(A_PhiV)*.                                (7)

This specifies a local mixed row at any fixed background. It can be
evaluated without guessing a quantum stress derivative: for H=tr_g h,

    delta R=nabla_a nabla_b h^{ab}-Box H-Ric_ab h^{ab},
    delta(Box R)=Box delta R-h^{ab}nabla_a nabla_b R
                    -(nabla_a h^{ab}-(1/2)nabla^b H)nabla_b R,
    delta M_em=0.

For the remaining invariants use C^2=Rm^2-2Ric^2+R^2/3 and
Euler4=Rm^2-4Ric^2+R^2, differentiating contractions as well as
curvature. Explicitly
delta Gamma^a_bc=g^{ad}(nabla_b h_cd+nabla_c h_bd-nabla_d h_bc)/2
and delta R^a_bcd=nabla_c delta Gamma^a_db-nabla_d delta Gamma^a_cb.
These formulas fix all local contacts. The seed's order-two bound
for the actual York/TT mixed block is not an order claim for (7) in
the fixed disk variables.

The product metric depends only on V. Identify orbit spinors with
one fixed spinor-density Hilbert space. In two dimensions the
conformal Dirac law becomes

    D(V+v)=exp(-v/2) D(V) exp(-v/2),
    D_i=-(v_i D+D v_i)/2,
    D_12=(v1 v2 D+v1 D v2+v2 D v1+D v1 v2)/4.         (8)

Multiplication by v_i commutes with v_j, not with D. For Q=D^dagger D,

    Q_i=D_i^dagger D+D^dagger D_i,
    Q_12=D_12^dagger D+D^dagger D_12
                   +D_1^dagger D_2+D_2^dagger D_1.

The outer complementary first-order projectors remain the retained
ones. Allowed fast and correlated-mass directions have v=partial_n v=0
at the product boundary, so multiplication in (8) preserves the first
domains and the induced squared mixed data, including their normal
derivative. The bolt is a smooth disk origin, not a new boundary.
For compact smooth directions the differentiation is immediate on
the common domains. The established weighted kernel estimates extend
the resulting linear Hessian to the specified remainder space; this
does not prove nonlinear C1 dependence on an arbitrary C5 family.

With the complete magnetic angular heat trace Theta_q, retain exactly
the [heat Hessian](stationary-shape-logarithmic-response.md):

    K12(t)=-t Tr(Q_12 exp(-tQ))
       +t^2 integral_0^1 Tr[Q_1 exp(-t(1-a)Q) Q_2 exp(-taQ)] da,
    Gamma_can''=(1/2) partial_s|0 {
       mu0^(2s)/Gamma(s) integral_0^infinity
       t^(s-1) Theta_q(t/B^2) K12(t) dt }.             (9)

The braces denote common four-dimensional Mellin continuation with
the retained subtraction. The full angular tower, lowest level and
contact vertex remain. The local part of A_VV is the second variation
of the selected gravitational/magnetic action and the complete bulk
and boundary Weyl restoration. Its nonlocal part is (9). Phi enters
the determinant through that restoration, so a second copy must not
be added. This identifies the full coupled Newton operator.

The original local-symbol input is
[Moller](https://arxiv.org/html/0705.3857v3), used only in the local
scope justified in the shape-response note. The retained heat
conventions agree with [Vassilevich](https://arxiv.org/html/hep-th/0306138v3).
The [Minkowski linear-response study](https://arxiv.org/html/gr-qc/0209075)
does not supply a global inverse or a stability result for this
curved thermal boundary problem.

## 4. Weighted consequence and an exact compatible formulation

The existing mixed principal form in theta=S+T/3, b=2T is

    p^4 theta^2/(80pi^2)
      +p^4 T^2[log(p/reference)-6/5]/(120pi^2),

with the retained ordinary shape term and all lower terms understood.
Use the weighted-map theorem's disk/interior patches and outer spaces

    X_delta: z^-delta theta in H_b^4,
             z^-delta b in H_b^(4,log),
    Y_delta: z^-delta outputs in L_b^2 x L_b^2,
    dvol_b=dz/z, dmu_g=gamma(z)z^-3 dz/z,
    X_delta^m=X_delta plus span{chi(z)z^3(1,-2)}, delta=7/2. (10)

Use the dual transpose when changing output variables. These are b-volume
weights, not the same exponents in physical L2. Under that theorem's
admissibility hypotheses, ind(A)=-2 and ind(A_m)=-1. This is reuse of
an existing theorem, not a new index calculation.

The density correction F tensor l is compact X_delta^m -> Y_delta.
On compact regions use Rellich for the fourth-order input. On a short
outer end its operator norm tends to zero because F=O(z^2 log^N z)
at worst; for the actual accepted candidate F=O(z^4 log^N z).
The mass column is finite dimensional and has a fast product with F.
Thus DF has the same index, although its adjoint differs. Multiplication
by Jvol in (3) is an output isomorphism and also preserves the index.
If the physical kernel has dimension k, the cokernel has dimension k+1.
There is no inverse onto all Y_delta in these physical spaces.

At the analytic seed the earlier outgoing/TT proof gives k=0. Let
A0 denote that full Hessian with its mass column and chi0 its actual
adjoint mode, normalized in any fixed nonzero way. The
[phase moment](stationary-scalar-phase-moment.md) gives

    L0(f)=integral chi0 dot f dmu_0,
    Ran(A0)=ker L0, A0 injective.

The pairing is bounded on Y_delta: the slow scalar amplitude is
O(z^(3/2)) with the proved subleading corrections, while delta=7/2
and the physical volume has exponent -3. The remaining branches
are controlled by the accepted adjoint theorem. Choose a smooth
compact source q with L0(q)=1 and set Pi f=f-q L0(f). Such q exists
by density and the nonzero adjoint. The bounded inverse theorem
gives G0:ker L0 -> X_delta^m. Its norm is not quantified here.

There is then an exact augmented isomorphism:

    [A0 q](u,eta)=f,
    eta=L0(f), u=G0 Pi f.                            (11)

For sum norms its inverse norm is at most
||G0||(1+||q|| ||L0||)+||L0||. The auxiliary eta is a computational
obstruction coordinate, not a new physical source, mass or temperature.
An augmented solve is a physical solve only when eta=0.

For the fixed-density nonlinear gradient (3), on any family where it
is defined into Y_delta, the exact equations are equivalently

    Pi mathfrakF(u)=0, L0(mathfrakF(u))=0.             (12)

The seed projection is fixed. At u=0, Pi D mathfrakF=A0 is an
isomorphism onto ker L0. A local projected construction still needs
the nonlinear mapping theorem and quantitative inverse/remainder
bounds, together with a defect inside their radius. Even if it is
constructed, its remaining scalar value in (12) must be evaluated.
No adjustable gauge, second mass datum or variable prepared clock
has been introduced to force that value to zero.

The signed quantity is determined by the retained model; its value
has not been proved. It is not a freely selectable missing physical
input. Neither negative index nor an unevaluated moment proves
nonexistence. Conversely the small complete defect does not prove
compatibility or a nearby exact physical solution.

## Verification contract and next proof

Before execution, freeze the note and checker. Revalidate the exact
current inputs of the seven cited component receipts. Check the
curvature/density identities, arbitrary-background factorization,
bookkeeping specialization, and block/range algebra exactly.

Use an independent warped-metric fixture
a=rho(1+rho^2/7+rho^4/101), r=2+rho^2/5+rho^4/103 and
s(rho)=1/13+rho^2/17+rho^4/109+rho^6/601.
Calculate sectional curvatures directly after the metric deformation,
including the radial lapse, at rho=1/5,2/3,7/4.
At 80 digits compare the predicted Hessian with 64/128-node complex
contours of radius .005 (relative tolerance 1e-45) and central real
steps .0001/.00005 (tolerance 1e-5, refinement factor .35).
Detect the error from omitting the density contact.

For a fixed nonnormal 3-by-3 Dirac matrix and two commuting diagonal
multipliers, check the square vertices and heat Hessian against
independent real mixed differences at steps .0002/.0001 and heat
times .07/.4. Use tolerance 1e-5 and refinement factor .35, allowing
an absolute roundoff floor 1e-50. Check cancellation of the finite
matrix log-determinant Hessian and detect omission of its contact.
These fixtures test formulas; they are not physical fields or a
spectral approximation to the actual exterior.

No test is claimed to certify the analytic continuation, the global
Fredholm theorem, a numerical inverse bound or source compatibility.
Preserve all runs and failed/superseded exact sources. The next proof
is a quantitative bound for G0 in the actual fast fourth/log norms
and nonlinear family control for (3), followed by evaluation of
the scalar in (12). An outgoing inverse alone does not supply them.

Remaining to-do list: quantitative compatible inverse and nonlinear
family; actual signed compatibility and coupled solution; exterior/state,
charge and return; opening, stability, transport, preparation, creation
and full-result audit.
