# Recentered canonical variables, mixed boundary limit and TT budget

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The preceding goal turn was progress (technical e4e94a2, planning 6416b1f).
The full stable semiclassical wormhole objective and S169.2 remain open.

This audit continues the SAME additional Lorentzian c_R=10^22 at the fixed
mu0, with all other added finite coefficients zero. It supplies a local
TT comparison budget and a mixed boundary limit, corrects an overbroad
regularity warning, and specifies the quantum calculation still needed.
It does not transfer the old quantum constants or construct a solution.

## 1. Exact conformal quotient and a dated regularity correction

Fix a reference sphere radius b_ref, the thermal period, flux, paired
Dirac domain and subtraction scale. For a spherical four-metric write

    g4=exp(2sigma)(h2+b_ref^2 dOmega^2), r=b_ref exp(sigma),
    hhat4=(b_ref/r)^2 g4=h2+b_ref^2 dOmega^2.                 (1)

A pure four-dimensional Weyl change g4->exp(2u)g4 changes r by exp(u).
It therefore leaves hhat4 EXACTLY unchanged, at fixed coordinates.
The conformal covariance and paired cocycle identity retained in
[the shape response](stationary-shape-logarithmic-response.md) give

    Gamma_D[g4]=Gamma_hat[h2]+W[sigma;h2],
    D_sigma Gamma_hat=D_h D_sigma Gamma_hat=0.               (2)

W includes the full local boundary cocycle and common finite restoration.
Its mixed derivatives do not vanish. Shape-shape derivatives of Gamma_hat
remain nonlocal, with the full angular tower and contact vertex. In an
affine tensor chart the corresponding pullback/gradient contacts must
still be included. Equation (2) is not permission to drop them. It is
an existing structural identity, not a newly discovered cancellation.
Neither b_ref nor mu0 is reset as the metric varies.

The [preceding scalar note](recentered-r2-scalar-boundary.md) correctly
finds u~r^-4 log r and failure of closed-compactification C4. Its inference
that this alone requires a new rough-collar determinant framework is too
broad. The [retained canonical theorem](stationary-centered-canonical-response.md)
explicitly uses C5 on OPEN charts with the weighted collar condition

    |v^(j)(Y)| <= e j!100^j Y^(4-j)(1+|log Y|), j<=5.        (3)

The [actual-domain outer neighborhood](finite-outer-admission-and-response.md)
uses exactly N_(4,5) from (3), the cubic mass coefficient and global C1
control. It already allows Y^4 log Y; it does not assume closed-boundary
C5 or spatial analyticity of the unknown. A pure conformal correction
does not change its canonical input at all by (1). This is a dated
correction, not a modification of the frozen preceding note or receipt.

For the ORIGINAL zero-added-c finite comparison, the outer theorem gives
background N<1.01e-7, allowed N<=1.5e-7, margin >4.9e-8,
first response <1.8e-7 and mixed response <.02 in its specified norms.
These numbers are NOT transferred to the new c_R geometry. Its actual
canonical coefficient bounds, state/gap, cap attachment/separation and
full two-direction response still need verification. The finite-cut
normal heat rows found in the previous audit are also already included
in [the complete boundary Weyl formula](stationary-exterior-boundary-weyl.md);
they must not be counted as a missing term newly discovered here.

## 2. Local action in the same variables

The auxiliary local Lorentzian density is

    (R-2Lambda)/(16piG)+c_R(2phi R-phi^2)-I_EM/(4g_R^2),
    F_s=1+32piG c_R phi.

With H=R[h2]+2/b_ref^2, integration of the Weyl Laplacian gives

    S_loc=b_ref^2/(4G) int sqrt(-h2) {
      exp(2sigma)F_s(H+6|grad sigma|^2)
      +6exp(2sigma) grad F_s.grad sigma
      -(2Lambda+16piG c_R phi^2)exp(4sigma)}
      -pi b_ref^2/g_R^2 int sqrt(-h2) I_hat.               (4)

The local GHY term after that integration is
b_ref^2/(2G) int sqrt(|gamma_h|) exp(2sigma)F_s K_h.
Indeed the original four-dimensional K transforms as
exp(-sigma)(K_h+3n sigma); its normal-sigma term cancels the integrated
Laplacian flux. This step does not supply the quantum boundary action.

For F_s>0 set alpha=sigma+log(F_s)/2. The braces in (4) become

    exp(2alpha){H+6|grad alpha|^2-(3/2)|grad log F_s|^2}
      -exp(4alpha)(2Lambda+16piG c_R phi^2)/F_s^2.          (5)

The local boundary term is now b_ref^2/(2G) int exp(2alpha)K_h.
The canonical quotient hhat4 is unchanged even in the Einstein frame.
The complete quantum action is still Gamma_hat[h2]+W[alpha-log(F_s)/2;h2].
Thus this change of variables exposes the local scalar/metric separation
without deleting the anomaly or changing any physical normalization.

## 3. Mixed scalar/shape boundary limit on a specified class

Work first at smooth finite cuts, using the complete local Weyl boundary
functional that accompanies the physical-radius bulk representative in
stationary-exterior-boundary-weyl.md. In a fixed outward coordinate x,
write hhat4=M^2 dx^2+N^2 dU^2+b_ref^2 dOmega^2,
D=M^-1 partial_x, a=D log N, s=D sigma and b2=b_ref^-2. Then

    W_boundary=-4pi b_ref^2 P_f int dU N Gcal,
    P_f=4/[360(4pi)^2] per exterior,
    Gcal=sigma(12a Da+44a^3/7+10a b2)
       -6(Da)s-(24/7)a^2s-6a Ds-(96/7)as^2
       -18s Ds-(85/21)s^3.                              (6)

Use an inward smooth canonical Y=-x with background M=1,
N=N0[1+O(Y^2)], sigma=-log Y+O(1). For fixed-source tangents let

    delta_sigma=u=O(Y^4 log^p Y),
    delta log N=j=O(Y^3)+O(Y^4 log^p Y),
    delta log M=h=O(Y^3)+O(Y^4 log^p Y),                 (7)

with differentiated bounds through the orders in (6). The cubic
coefficients are finite but not set to zero. These are selected AdS
tangents, not arbitrary finite-wall Dirichlet data. Proper derivatives
must vary: delta a=Dj-a h, delta s=Du-s h,
delta(Da)=D^2j-2(Da)h-a Dh, and similarly for Ds.

The leading mixed term can be given explicitly. Put N0=1 locally,
j=j3Y^3+..., h=h3Y^3+..., u=Y^4(u4 log Y+v4)+... . For the two-parameter
family N->N exp(tj), M->M exp(th), sigma->sigma+e u,

    partial_e partial_t (N Gcal)|0
      =Y^4/7 { (5028h3-452j3)(u4 log Y+v4)
                        +(2769h3-113j3)u4 }+o(Y^4).     (8)

The remainder may contain higher powers times finite logarithms. The
N0 factor is restored for a general fixed clock. Equation (8) follows
by differentiating D as well as N; setting h=0 is not the full test.
For general (7) every term is O(Y^4 log^q Y), hence tends to zero.
For example s=O(Y^-1), Ds=O(Y^-2), delta_u s=O(Y^3 log^p),
delta_u Ds=O(Y^2 log^p); inserting j,h or their derivatives in the
cubic terms of (6) leaves precisely that weight. The other terms are
higher weight. This is a finite-jet estimate before taking the cut limit.

The local EH+c_R R^2 mixed boundary variation has the same or better
weight: physical boundary measure is O(Y^-3), scalar field/curvature
jets have O(Y^4 log^p), shape jets O(Y^3), and physical normal derivatives
preserve these powers. The intrinsic retained counterterms also vanish
in this mixed limit. The canonical determinant's sigma-h mixed derivative
is identically zero in (2). Thus the mixed scalar/shape boundary pairing
vanishes on (7), including the radial metric variation and fixed clock.

This does not prove the shape-shape boundary realization, all off-shell
boundary momenta, or an open physical metric domain. The finite-cut
Navier obstruction from the preceding note remains valid. The old
mass-channel boundary cancellation is useful but is not a full Dirichlet
completion. No extra wall or compensating action has been selected.

## 4. A quantitative TT comparison at the recentered RN reference

Use g=G lambda, g+11g^2/(360pi)=2000^-10,
F_lambda=1-384pi c_R g, Z=piG q^2/(g_R^2 F_lambda),
f=1-2M/r+Z/r^2+lambda r^2 and R0=-12lambda. Here

    0<g<9.765625e-34, 1-1.179e-8<F_lambda<1,
    Z=B^2+3lambda B^4, B^-2/lambda<beta0=1.223e-10.       (9)

There is an exact local polynomial identity, with Lambda=-3lambda0:

    (R-2Lambda)/(16piG)+c_R R^2
      =F_lambda(R-R0/2)/(16piG)+c_R(R-R0)^2
                              +(R0-4Lambda)/(32piG).   (10)

In the Euclidean action the gravitational signs are reversed. Define
Gamma_rest by subtracting the constant vacuum-volume functional from
Gamma_D so that it cancels the last term in (10). This is an exact add/
subtract split, not an assertion that the nonconstant quantum stress is
known. The comparison is Einstein-Maxwell with G_eff=G/F_lambda and
cosmological constant -3lambda; it is stationary at this RN geometry.
The remaining local square has zero first variation there and Hessian

    I_square''[k,k]=-2c_R int (delta R[k])^2 dvol.         (11)

For genuine spherical TT components k=diag(-b-2xi,b,xi,xi),
the constraint and its unique bolt-regular solution are

    (r^2 f b)'=-(r^2 f'-2rf)xi,
    b=-(r^2 f)^-1 int_(r_+)^r (y^2 f'-2yf)xi dy.        (12)

The end tail in (12) is retained. Since tr k=div k=0,
delta R=-J:k. Here J=diag(-t,-t,t,t), t=Z/r^4, so

    delta R=-4t xi,
    I_square''=-128pi beta_T c_R int r^2 t^2 xi^2 dr.    (13)

The symbol beta_T is the physical Euclidean time period, not beta0.
The extra R^2 term is thus zeroth order on this TT space; it is not
estimated by a fourth-derivative norm or treated as a small scalar term.

The positive classical argument of
[the global TT reduction](stationary-global-tt-reduction.md), sections
3-4, applies with d=0 and G_eff. For clarity its quantitative hypotheses
are checked, not merely transferred. Write x=lambda B^2>1000,
M=B(1+2x+delta/2), delta>0. The retained preparation has
gamma_old=kappa B_old/C_old<1e-30, C_old<8k^3.
The old product has |B_old^2/B_class^2-1|<1.4e-22. Equations (9)
give B/B_old<2 and C_RN=1+6x>49119000000. Thus at the SAME kappa,
gamma_new=kappa B/C_RN<3e-30. The d=0 thermal relation hK(h)=gamma_new,
|K'|<3 on h<=.001, gives h=(r_+-B)/B<6e-30 and
delta=(h^2/(1+h))[1+x((1+h)^2+2(1+h)+3)]<1.
These loose guards suffice and do not change the preparation interval.

In particular the York projection potential obeys

    V=f''/4+f'/(2r)+3f/r^2+3(f')^2/(8f)
      >=3lambda/2+Z/(2r^4)>=3lambda/2,
    ||(-P)^-1||_(L2->L2)<=2/(3lambda).                 (14)

This is a projection bound, not a physical inverse. For the stationary
classical comparison, integration by parts using (12) gives

    I_EM,eff''=beta_T F_lambda/(2G) int {
       (3/2)r^2 f(xi'+2xi/r)^2+W_f xi^2-4fb xi}dr,
    W_f=3-12M/r+13Z/r^2+3lambda r^2.

If H(r)=int_r^infinity xi(y)/y^2 dy and P0=r^4f'-2r^3f,
the cross term is 2int P0'H^2+2r_+^4f'(r_+)H(r_+)^2.
Here P0'>=-6r^2 and int r^2H^2<=(4/9)int xi^2. Consequently
it is >=-(16/3)int xi^2. As in the cited proof,
W_f-16/3>=lambda r^2 for x>=1000, delta<=1, d=0. This uses
max_(t>=1)(24t-t^4)=18*6^(1/3)<36<39 and
3-12/t+13/t^2>=3/13. Therefore

    I_EM,eff''>=E_eff[xi]
      :=beta_T F_lambda/(2G) int {
        (3/2)r^2 f(xi'+2xi/r)^2+lambda r^2xi^2}dr.      (15)

The same bolt cutoff and infinity closure proof applies to this analytic
RN geometry. The Schur estimate for (12) gives ||b||<=2||xi||,
||-b-2xi||<=4||xi|| in L2(r^2 dr). Hence (15) has physical tensor
gap F_lambda lambda/(176piG). No unprojected TT invariance is asserted.

Since t<=(3+beta0)lambda on this exterior, (13)-(15) imply

    |I_square''| <= eta_R E_eff,
    eta_R=256pi c_R g(3+beta0)^2/F_lambda <7.07e-8,
    I_EM,eff''+I_square'' >=(1-7.07e-8)E_eff.            (16)

The bound is relative to the displayed fixed physical energy, not a
normalization gain. It is a new local comparison result on the selected
reference. It does not include the nonlocal quantum remainder. In
particular (16) is not the full inverse or an all-sector stability result.

## 5. Feasibility targets and the research decision

For the full pulled-back TT row, including off-shell gradient/chart/Lie
contacts and the physical boundary action, the decisive remaining form
estimate is

    Gamma_rest''[TT,TT] >=-eta_Q E_eff,
    eta_Q <1-7.07e-8; select eta_Q<=1/4.                (17)

The sufficient chosen target leaves >.7499999293 E_eff. Eta_Q is UNKNOWN
on this new geometry. The old .005 quantum loss was for another seed;
it is not a value for eta_Q. A C5-to-C0 source estimate is not (17).
The positive p^4 log p term is unbounded relative to classical energy,
so a two-sided small-perturbation estimate is inappropriate. Prove a
one-sided bound on the negative finite part and contacts, retaining the
positive high-frequency form. Closedness, the weighted outer/bolt
realization and the quantitative inverse in the solver norm are separate
requirements. Equation (16) makes the added local square a negligible
loss in this channel; it does not settle those continuum steps.

For the properly realized coupled Jacobian [A+E,Bmix;Cmix,L], retain
the scalar energy comparison ||A^-1||<=1 and require

    ||E||<=1/4,  ||L^-1 Cmix|| ||Bmix||<=3/8.            (18)

Then the Schur scalar inverse is <=8/3, since ||(A+E)^-1||<=4/3
and the loop is <=1/2. The shape/current/constraint row L, its actual
inverse norm, the two mixed norms and transfer E are UNKNOWN. The exact
canonical cancellation (2) removes one source of mixed nonlocal terms;
the local anomaly and reconstruction still have to meet (18).

For a complete physical-admission norm X and projected residual F,
an injective preconditioner Bp on the projected range must satisfy

    Y0=||Bp Pi F(0)||_X,
    Z0=||I-Bp Pi F'(0)||,
    ||Bp Pi(F'(u)-F'(v))||<=Z2 ||u-v||_X,
    Y0+Z0 rho+(Z2/2)rho^2<=rho,  Z0+Z2 rho<1,
    N_i(background)+C_i rho<=R_i for EVERY admission norm i. (19)

For rho=1e-8 the sufficient allocation is Y0<=5e-9, Z0<=.25,
Z2<=5e7: the radius sum is exactly 1e-8 and contraction <=.75.
Each unknown embedding must meet C_i<=(R_i-N_i(background))/1e-8;
a nonpositive margin cannot be repaired by a smaller inverse. The old
outer margin 4.9e-8 would allow C_outer<=4.9 if it applied; its new
background margin is unknown. No numerical solver ball is certified.

The available scalar-row particular has
||u_tr||_E,rad/sqrt(lambda r_+^3)<1.82e-29 and forced coefficient
a/r_+^4<6.60e-30. On its own it could fit a 5e-9 admission allocation
if the corresponding embedding constant were <=2.74725e20. Allowing
the Schur amplification 8/3 lowers that scalar-only target to
1.03021e20. Both are explicit thresholds for UNKNOWN embeddings, not
their estimated values. The shape-forced contribution
(8/3)||Bmix L^-1 F_shape||_E* and reconstruction must fit the remaining
budget; neither is set to zero. No scalar energy bound substitutes for
cold jets, canonical admission, positivity of F_s or the complete metric.

The original construction's latest complete interior allowance is still
C_I=9.1591500000000000015e103; its global C_F=1.069626977892751e201,
interior-only C_F,I=4.141944868545865e118, Cstart=3.123143830205975e28,
mhi=2.577600186401461e12. In the original fixed norms, the actual lower
bounds R>.97 and C_D4>2 make the current global upper certificate miss
the conservative target by >4.4868e225/epsilon_c, and the interior-only
one by >1.7374e143/epsilon_c. Even deleting the high term would require
R<=3.90857e-8 epsilon_c. These are KNOWN upper allowances that FAIL to
certify, not unknown constants or lower bounds on the actual correction.
Their dominant 1e160 outer and 1e112 interior Cauchy-derivative losses
remain a reason to park that uniform chain. They are not estimates for
the new c_R solution, nor a nonexistence theorem. The dated
[full feasibility audit](construction-feasibility-audit.md) retains all
exact thresholds and signs. Existing cancellations are not credited twice.

DECISION: continue the specified c_R route using the canonical quotient
and a one-sided, frequency-resolved TT form calculation. The calculation
executed here shows the added local TT loss <7.07e-8, rather than a
large fourth-derivative loss; the mixed boundary limit also passes on
(7). Next bound the actual negative finite quantum/contact form against
(17), then the local mixed Schur product (18). Work on the actual state,
paired boundary realization and reference-to-finite admission in parallel
as mathematical requirements, not by importing old receipt successes.
This is a specified sharper method and a quantitative next decision,
not a conclusion that the full construction is feasible.

The separate radial/Maxwell constraints and fixed charge must be retained:
trace and one shape row alone do not determine all integration constants.
Signed compatibility remains an exact equation after projection, with no
spare temperature, flux or asymptotic-curvature knob. The old signed moment
cannot be reused on this theory. Physical matching and the nonlinear
source/metric domain precede any claim of a wormhole. The selected scalar
mass m^2G in (3.31e-25,3.32e-25) avoids the old trans-Planckian scalar
scale, but does not control omitted operators, noise, the quantum state,
all modes of the coupled system or causal stability. The current RN
comparison itself is not a traversable wormhole. No conclusion of
stability, construction feasibility or nonexistence follows from this audit.

## 6. Frozen supporting verification contract

Freeze this note, checker and review before the first scientific run.
Retain every run, source snapshot, failure and raw row. Check (4)-(5),
the exact local split and TT contraction independently by symbolic
variation. Check (8) by varying the unexpanded finite-cut expression,
including radial M, and by 60/90-digit finite differences of that
expression. Use fixed central steps 1e-5 and 5e-6 with expected quadratic
refinement; the steps change tangent parameters, not geometry or theory.
Use rational upper/lower arithmetic for (9), the thermal guard and (16).
Support the classical TT identity with independent RN radial quadrature
of the original areal variation and the completed-square identity, keeping
the nonlocal b tail and an explicit positive-temperature bolt.
These are proof diagnostics, not quantum spectral data or a solve of the
coupled geometry. They cannot establish the continuum hypotheses in
(17)-(19). No Lean invocation is selected.

Remaining to-do list: actual quantum TT/contact form and mixed Schur
bound; full physical boundary domain and coupled inverse; finite state,
residual, admission and nonlinear solution; exact compatibility, matching,
validity, all-sector stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
