# Actual global mixed inverse and the remaining logarithmic Schur operator

2026-10-24, E003/S169. Ordinary argument before frozen verification.
Full stable semiclassical wormhole objective ACTIVE/INCOMPLETE. The closed
audit, physical action, paired domain, state, scale and admission targets stay.

OUTCOME: complete coupled estimate UNRESOLVED, with a certified global
mixed-row inverse and an exact coupled reduction. In dimensionless physical
L2 units its weighted inverse and weighted adjoint inverse are <12.500001.
The preconditioned mixed remainder is <1.556e-23. This closes an actual
global operator ingredient, not a mass column or the original fast graph.
The remaining Schur operator has a logarithmic principal part; treating it
as identity plus a small bounded error fails in the continuum.

## 1. First unsupported step and the actual mixed row

The [preceding two-row test](current-canonical-weighted-inf-sup-attempt.md)
rejects unchanged York-sign coercivity after canonical weighting. Keep the
same W=exp(omega), exponent s=1.4, and original global two-equation Hessian.
Instead of estimating a diagonal alone, construct its mixed inverse first
and carry both equations into the exact Schur reduction below.

Let epsilon=G/(60pi), D=96piG c_R-epsilon, c_R=1e22. Use A_r for the
physical orbit scale, to distinguish it from the conformal operator. Put
beta=lambda A_r^2. The exact conformal trace density is

    T=R+12lambda0+epsilon Q_anom-D Delta4 R,
    Q_anom=7 C4^2/12+11 J4^2/6-20 I-11R^2/72.

For a canonical orbit variation T_c at fixed physical sphere radius r,

    delta R=2[P-t0]T_c,
    P=-Delta2+2/r^2+12lambda0, t0=R+12lambda0,
    delta(Delta4 R)+2T_c Delta4 R=Delta4(delta R).

Therefore the COMPLETE mixed trace-density operator C is exactly

    C=2 L_D P+E,
    L_D=1-D Delta4,
    E=2D Delta4(t0 .)+epsilon delta_density Q_anom.       (1)

The action mixed row is -C/(8piG); its other row is the physical L2 adjoint.
The density contact in (1) has been cancelled exactly, not omitted by assuming
stationarity. The full conformal Ward identity is used at a COMMON regulator
with its fixed local restoration. The canonical nonlocal determinant has
zero conformal and mixed derivatives in these quotient coordinates. Thus
all angular/thermal sums and finite returns are accounted for in this row
by that identity and the full anomaly. The independent shape-shape operator
still contains their complete nonlocal response and is retained in section 6.

## 2. Domain and both weighted second-order inverses

Use physical radial L2 with measure r^2 dr, including the unchanged common
period/sphere factor when desired. The orbital measure for Delta2 is dr.
The map U:u->r u is UNITARY between these two L2 spaces, not a rescaling of
the physical admission norm. The physical adjoint is P*=r^-2 P r^2.

Keep the exact finite seed

    z=B/r, y=1-z, z_h=1/(1+h),
    Q=y^2 P_seed(z)-tau z^3,
    P_seed=beta(1+2z+3z^2)+u z^2+epsilon_seed z^4,
    6beta+u+epsilon_seed=1, u=1/C_geom, epsilon_seed<0.

P_seed is positive and increasing on [0,1]: P_seed'>=2beta-4|epsilon_seed|>0.
Its endpoint is exactly one. Hence 0<Q<=y^2 on the actual exterior; Q>0
uses the retained regular finite-temperature root and monotonic lapse.
There is no h=0 replacement. The original weight obeys

    omega_z=s theta/(h+y), 0<=theta<=1,
    theta=0 through .1, theta=1 from .2,
    theta=I_(10z-1)(11,11) in between.

It is one at the paired seam and smooth even at the regular bolt. In physical
units A_r^2|grad omega|^2=z^2 Q omega_z^2<=s^2. More strongly, for BOTH signs,

    A_r^2 |grad(omega +/- log r)|^2
       =Q(z omega_z -/+1)^2
       <=[s z+y]^2<=s^2.                              (2)

The same bound holds with -omega. This uses s>=1 and Q<=y^2; it includes
the radial-measure conversion. Bounding the two gradients separately would
lose this cancellation and would not prove the needed physical L2 estimate.

Under U, W P W^-1 and W P* W^-1 become the conjugates of the self-adjoint
orbital P by exp(omega+log r) and exp(omega-log r). Their real quadratic
forms are

    ||grad v||^2+integral [2/r^2+12lambda0-|grad psi|^2]|v|^2 dr.

Consequently, with p0=12beta-s^2,

    p0> .03999999995,
    ||W (A_r^2 P)^-1 W^-1||,
    ||W (A_r^2 P*)^-1 W^-1|| <=1/p0<25.0000001.       (3)

The same estimates hold with W replaced by W^-1. The real part also controls
the orbital gradient. The first-order imaginary/skew form terms are bounded
by (2); the positive potential and bounded seed coefficients make the forms
continuous and coercive on the corresponding H1 space. The variational
inverse therefore exists, and the adjoint estimate holds on the same domain.

Use smooth bolt-regular functions and the physical decaying end, then take
their form closure. A punctured-bolt boundary condition is NOT used. Radial
regularity is the rotational smooth-disk condition. At infinity the allowed
mass cubic, quartic/log and correlated slow data are in this physical L2
class. Boundary products for these local second/fourth-order forms vanish
as O(r^-3 log^N r) or faster. The retained fixed-source boundary one-form
has zero linearization on this phase; see
[the actual boundary check](current-seed-boundary-phase-and-trace-defect.md).
No smooth reflected mass cubic is presumed. No quartic coefficient is set
to zero: the inverse generates its coefficient from the full source and
regular bolt. Its numerical fast/finite-part trace is still a separate map.

The endpoint can be checked directly rather than inferred from an L2 bound.
For a power z^a the exact outer normal mixed factors are

    C_normal(a)=2lambda[C_v-Dlambda a(a-3)](4-a)(a+3),
    C*_normal(a)=2lambda[C_v-Dlambda a(a-3)](6-a)(a+1),
    C_v=1+11epsilon lambda/3.                         (3a)

They follow from (1), including the R^2 anomaly and the exact lambda0 relation.
The first factor's decaying root is the heavy scalar root, well above six.
The slow local mixed roots are four and six respectively. The fourth-order
ODE has analytic one-sided seed coefficients, so its ordinary regular-singular
recurrence retains the mass power and can generate resonant quartic/log data.
Intermediate inverses are taken in this ambient physical L2/ODE domain; they
are NOT asserted to preserve the restricted forced-log coefficients or the
fixed fast-data ball. Generated coefficients, including any extra resonant
logarithm, must be carried through (11) and matched to the original phase.
This auxiliary domain does not replace the paired quantum realization or
the original final boundary conditions. In particular, a large separated
quartic/slow coefficient can cancel in L2; (3) does not control it separately.

For L_D, the physical measure is already correct. Its conjugated real form
gives m_D=1-Ds^2/A_r^2, with

    m_D>.99999996536,
    ||W L_D^-1 W^-1||<=1/m_D<1.000000035.             (4)

Again both signs of omega work. These elementary form arguments use the
actual smooth physical interior and original asymptotics. They invoke no
smooth-double pseudodifferential theorem for rough reflected boundary fields.

## 3. Quantitative physical derivatives needed by the entire anomaly

Here f is the dimensionless datum in A_r^2 P u=f. Write v=Wu and measure
all L2 norms physically. The unitary form estimate in section 2 gives

    ||v||<25.0000001||Wf||,
    A_r||grad v|| <=(p0^-1/2+p0^-1)||Wf||<30.000001||Wf||.

The beta-polynomial cutoff has theta_z<=37.001381<38. On its transition,
z^2<=.04 and y<=.9. The actual polynomial gives |Q'|<=4y, including the
thermal term 3tau/y<1e-30. Direct differentiation yields

    A_r^2|Delta4 omega|<9.6152<10,
    A_r^2|Hess omega|<7.209<8.                       (5)

For the Hessian, its radial, time and sphere entries are bounded by
s(1.368+1+2+.25), s(2+.25), and s/4, respectively, with two sphere entries.
No singular bolt connection is bounded separately.

Use P=-Delta4+2grad(log r).grad+2/r^2+12lambda0, where
A_r|grad log r|=sqrt(Q)<=1 and A_r^2(2/r^2+12lambda0)<2.001. The exact
conjugated equation and (5) give

    A_r^2||Delta4 v||<565||Wf||.

Integrated Bochner, with Ric>=-13lambda g and vanishing endpoint flux, then
gives A_r^2||Hess v||<610||Wf||. This follows directly by integrating the
scalar commutation identity; no global coordinate Schauder constant is
assumed. Returning derivatives through W gives the useful guards

    ||Wu||<26||Wf||,
    A_r||W grad u||<66||Wf||,
    A_r^2||W Hess u||<950||Wf||.                     (6)

All estimates also hold for W^-1. They are physical L2 derivative bounds,
not the original fast H3 data norm, C5 norm, or a nonlinear neighborhood map.

To bound the complete anomaly on u=(A_r^2P)^-1 f, let Pi be the orbit
projector and h=2u Pi. The warped metric identities give
A_r|grad Pi|<20 and A_r^2|Hess Pi|<100; these follow from
F/r^2<=5lambda, F'/(2r)<=11lambda and lambda<=F''<=21lambda.
Using (6),

    ||Wh||<78||Wf||,
    A_r^2||W nabla^2 h||<14000||Wf||.

The connection-difference formula bounds
|delta Rm|<=6|nabla^2 h|+10|Rm||h|. Since A_r^2|Rm|<10, it gives
A_r^2||W delta Rm||<1e5||Wf||. Contractions, including inverse-metric
variations, give respective guards 1e6,1e7,1e7,1e8 for delta Ric, delta R,
delta J and delta C in the same units. Together with |C|<35lambda,
|J|<7lambda, the exact curvature-quadratic coefficients and density variation,
the ENTIRE curvature anomaly obeys the generous bound

    ||W A_r^4 delta_density Q_curv u||<1e12||Wf||.    (7)

For the magnetic term I=q^2/(2r^4), delta I=0 in an orbit variation but
the density contributes -40 I u. Its contribution below is bounded using
Gq^2/(lambda B^4)<4e-27 and the first bound in (6). It is not dropped.

## 4. Execute the full mixed inverse with its remainder

Define Pbar=A_r^2 P and Cbar=A_r^2 C. Factor (1) as

    Cbar=2 L_D (I+R_mix) Pbar,
    R_mix=(1/2)L_D^-1 (A_r^2 E) Pbar^-1.             (8)

The potentially dangerous derivative of t0 cancels at the operator level:

    L_D^-1 D Delta4=L_D^-1-I.

Thus the first part of R_mix is
(L_D^-1-I)(A_r^2 t0)Pbar^-1. This identity is first on the regular core,
then on the closed composed domain; it does not commute t0 through Delta4.
The current seed bounds give

    A_r^2||t0||infinity
      <=beta[1.1e-27+12*11G lambda/(360pi)].

With the exact vacuum relation and retained parameter intervals, its
preconditioned norm is <9.167e-27. Equations (4),(7) and the full magnetic
term give, for either sign of the weight,

    ||W R_mix W^-1||
       <=(1+1/m_D) A_r^2||t0||/p0
         +[epsilon/A_r^2 *1e12
              +40beta(4e-27)*26/(60pi)]/(2m_D)
       <1.556e-23.                                  (9)

These are actual global row estimates, not comparison-model values. The
Neumann inverse in (8) is now justified on physical L2, and

    ||W Cbar^-1 W^-1||
       <=1/[2p0 m_D(1-1.556e-23)]<12.500001.         (10)

Taking the adjoint of the estimate for the opposite weight proves the same
bound for W(Cbar*)^-1 W^-1. In physical action units the mixed inverse
constant is (8piG A_r^2) times (10); that factor has not disappeared.
The composed domain is {u in D(Pbar): Pbar u in D(L_D)}, with the original
regular-bolt and normalizable endpoint realization, and its small perturbation
from (9). Local regularity identifies the differential expression. This
constructs the seed mixed inverse; it does not prescribe independent physical
boundary phases or transfer it to an unknown nonlinear metric family.

## 5. Both original equations and the exact coupled reduction

Normalize the full action Hessian by 8piG A_r^2, explicitly retaining that
physical row conversion. In canonical variables write

    Hbar=[[-Abar,-Cbar],[-Cbar*,Dbar]],
    Hbar[nu,T]^T=[f,g]^T.

Abar is the complete scalar row (including its fourth order), Cbar is the
actual mixed row just controlled, and Dbar is the COMPLETE shape row. Define

    K=Cbar^-1 Abar(Cbar*)^-1,
    S=I+K Dbar.

The exact elimination, initially on the common strong core, is

    S T=-Cbar^-1 f+K g,
    nu=(Cbar*)^-1(Dbar T-g).                         (11)

This retains both original equations and their signs. Every mass, scalar
quartic, correlated slow and fast source column enters the same (11).
The actual cut mass lift, including its transition and regular bolt, is
unchanged. It is not promoted to a global fixed-temperature zero mode.

The inverse (10) supplies global L2 bounds for two factors in (11). It does
NOT supply numerical strong graph bounds for K, S^-1, Dbar S^-1 or their
original finite-part rows. Those compositions must be interpreted in the
actual coupled mass/log graph, carrying generated quartic data explicitly.
No boundary equation is removed by the algebra, and no unknown phase map
is assigned zero. The previous TT phase defect is not reintroduced: W=1
on the entire original outer localization.

## 6. A principal cancellation the preconditioner must preserve

There is a further definite restriction on how to estimate (11). In a
smooth compact physical interior patch let k be proper frequency. The
retained full principal symbols, in unscaled trace units, are

    A ~6D k^4,
    C ~2D k^4,
    D_shape ~-(2/3)D k^4
                 +G k^4[log(k/mu0)-6/5]/(15pi).

The scalar fourth coefficient includes the same anomaly repartition as
[the physical metric block](physical-metric-block-and-preconditioner.md).
The quantum fourth/log coefficient is the complete four-dimensional local
coefficient, not one angular mode. The cited
[Moller Theorem 1 and Lemma 6.1](https://arxiv.org/html/0705.3857v3) concern
closed smooth spin manifolds and stable kernel dimension. Here only the
already derived local interior symbol is retained. That theorem is NOT
used for the paired global realization, its boundary inverse or neighborhood.

Since K~3k^-4/(2D), the leading constant parts in S cancel EXACTLY:

    1+[3/(2D)](-2D/3)=0,
    S ~ [G/(10piD)] [log(k/mu0)-6/5]+lower order,
    G/(10piD)=1.05542899627435...e-26.               (12)

This cancellation is structural, not a change of normalization. It makes
an estimate of K and Dbar separately followed by an identity-plus-small-error
argument inappropriate. Compact oscillatory functions on this actual smooth
interior patch see the nonzero log coefficient in (12). Consequently S-I
is not bounded on ordinary same-order L2/Sobolev space; its high-frequency
multiplier grows without bound. No frequency cutoff is introduced. This
rejects that proposed Neumann estimator, not S on its logarithmic graph.

A valid remaining preconditioner must preserve (12), retain the full
quantum principal order, supply low-frequency completion, and incorporate
the original global phase/bolt rows. The number in (12) alone supplies none
of those inverse constants. A small principal coefficient is not a small
operator over all frequencies, nor a certificate about physical fluctuation
power or semiclassical validity.

## 7. Feasibility decision and the one next calculation

CONTINUE WITH THE EXACT MIXED-INVERSE SCHUR METHOD, sharpening the remaining
logarithmic graph estimate. The seed mixed inverse is quantitatively
attainable without weakening s=1.4 or changing the boundary realization.
Do not optimize its already negligible 1.556e-23 remainder further.

The FIRST REMAINING UNSUPPORTED CONTINUUM STEP is an inverse/defect estimate
for the actual S on the mass/log graph with logarithmic principal part and
original finite-part rows, followed by the reconstruction in (11). This
includes the complete shape-shape quantum operator, all angular/thermal
tails, coincidence contacts and finite-domain returns. Neither the scalar
Ward reduction nor the mixed inverse evaluates that independent operator.

No complete post-elimination kappa_m, other return column, signed correction
or physical X-norm comparison is certified. At the proposed kappa_m=5e-5,
the affine debit remains 2.14765e-12, residual 5.34235e-12; the graph debit
remains .0024731, residual .0012069, before other entries and independent
.019523868 energy feedback. No claim is made that those entries fit.
The original mass conversion >.00906808164 and two-graph ceiling
|m|<2.205539e-9 remain. (10) is not interchangeable with kappa_m.

The conditional two-fast-input screen remains .135973<.1599164, reserve
.0239434, solely in its admitted-family scope. The original comparisons
and nonlinear inequalities remain unproved:

    E_center+1e-8 alpha<=1e-11,
    R_X alpha*.028913+Z0_other<=.25,
    R_X alpha^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the actual logarithmic
Schur operator S and the reconstruction (11) in the original coupled graph,
including its global finite-part/phase rows and all four return columns.
Test the resulting affine, feedback and physical admission constants against
the unchanged targets. Another mixed-row or local-symbol calculation is not
that calculation.

Executed here: the actual global mixed inverse and its adjoint, complete
trace/anomaly remainder, physical derivative conversions, exact coupled
reduction and its principal cancellation. Full coupled solution, nonlinear
admission/contraction, exact compatibility, physical matching, validity/noise
and all-sector stability remain unresolved. Supporting checks verify the
ordinary identities and enclosures; they do not replace the remaining argument.

Remaining to-do list: complete logarithmic Schur inverse and global phase returns;
actual correction and original physical-norm comparisons; nonlinear admission/
contraction and exact compatibility; physical matching/causal realization;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole objective.
