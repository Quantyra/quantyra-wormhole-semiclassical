# A common cap family restores the middle source-three allowance

2026-10-21, E003/S167. Ordinary argument before supporting verification.

With the original seed, state, angular tower and common subtraction,

    |nabla^j Theta_seed|/(j!100^j q^2)<1e26,
    0<=j<=3, s>=.001, Y>=.05.                       (1)

The proof below includes the fixed/moving cap correction identified in
[the splitting audit](seed-middle-splitting-audit.md). It restores the
OLDER 1e26 allowance; it does not restore the withdrawn 1e22 improvement.
The auxiliary cap is specified explicitly, transported to common physical
coordinates and differentiated with its density term intact. No physical
boundary condition, state or action coefficient is changed.

This resolves this particular middle-source input. The original metric
ledger used 1e26, so its numerical input is recovered without improvement.
Dependent estimates still require a scope review; the actual finite-metric
H6 source, transfer constants, cancellation and nonlinear solution are not
supplied here.

## 1. Analytic seed neighborhood and explicit auxiliary cutoffs

Use R0=1/32 and the retained middle normal chart. The Goursat equations
from [the source construction](stationary-reference-cold-bolt.md) are

    sigma_zw=Bcal''(u)exp(2sigma)/8,
    u_zw=Bcal'(u)exp(2sigma)/4,
    sigma(z,0)=sigma(0,w)=0.

The same iteration extends this chart to R*=1/30. Indeed its u change
is bounded by 2(.224)/30+7.08 exp(.004)/(4*30^2)<.02, and its sigma
bound is 8 exp(.004)/(8*30^2)<.002. Thus the original |Bcal'|<7.08,
|Bcal''|<8 guards still apply. The exact normalized polynomial

    Bcal(u)=x/C+u^2/C-(2+4x+delta)u^3/C
             +(1+3x+d)u^4/C-2d u^5/C+d u^6/C

also gives |Bcal'''|<17 and |Bcal''''|<13 in that neighborhood, using
x/C<1/5.9, 1/C<1e-10, |d|/C<1e-24, delta<=1 and |u|<1.02.
The Goursat Picard row sums are below .02 there. In particular
|sigma_zw|<1.005, |u_z|<.284, |sigma_zz|<.7|w|, and likewise with
z,w exchanged. On the slightly enlarged real/complex observation disk,

    |sigma|<.002, |partial sigma|<.067,
    |partial^2 sigma|<3.                            (2)

Use the degree-nine step

    eta(t)=sum_(i=0)^4 binom(9,i)t^i(1-t)^(9-i), 0<=t<=1,

extended by one to the left and zero to the right. It is C4, lies in
[0,1], and has zero derivatives through four at both endpoints. Set

    chi(r)=eta(128r-3),       transition [3/128,1/32],
    zeta(r)=eta(64r-1),       transition [1/64,1/32].

Here chi caps the normal conformal factor, while zeta extends the
normalizing coordinate map. The Born estimate allowed ANY such fixed
bounded cap cutoff; this selects one without changing the actual source.
For Cartesian derivative tensor norms the following bounds hold:

| cutoff | first | second | third | fourth |
| --- | ---: | ---: | ---: | ---: |
| chi | 315 | 7e5 | 3e9 | 4e12 |
| zeta | 158 | 1.8e5 | 4e8 | 4e11 |

For the one-dimensional step, |eta'|<=630/256, |eta''|<40,
|eta'''|<=1008 and |eta''''|<=9072. The last two follow from the
Bernstein finite differences. Radial tensor norms are bounded by
|f''|+|f'|/r at order two, |f'''|+2|f''|/r+2|f'|/r^2 at order three,
and |f''''|+3|f'''|/r+9|f''|/r^2+9|f'|/r^3 at order four.

## 2. Moving normal charts in common physical coordinates

Let b=(b1,b2) denote the physical center, beta=b1+i b2 and
betabar=b1-i b2. These are independent analytic expressions when b
is complexified; no complex conjugation of a varying parameter is used.
The local normalizer solves

    F_b(0)=beta,
    F_b'(z)=exp[sigma(beta,betabar)-2sigma(F_b(z),betabar)],          (3)

with the corresponding equation for the other coordinate. Its pulled
conformal factor is

    sigma_b(z,w)=sigma(F_b(z),Fbar_b(w))
                         +(log F_b'(z)+log Fbar_b'(w))/2.

It vanishes on both axes, as required. On the complex parameter bidisk
|b1|,|b2|<=rho=1e-5, for spatial arguments within 2rho of B_(R0), (2)-(3)
give the following bounds in Cartesian operator/tensor norms:

    |F_b-id|<=1.5rho, |DF_b-I|<=.2rho, |D^2F_b|<=10rho,
    |sigma_b-sigma|<=.3rho,
    |partial(sigma_b-sigma)|<=10rho.                (4)

For example |sigma_z(F,betabar)|<=1.005|betabar| and
|log F_b'|<.127rho. Integrating F_b'-1 costs less than .004rho,
in addition to the center vector of norm at most sqrt(2)rho.
The second derivative follows from F_b''=-2sigma_z(F,betabar)(F_b')^2.
The last bound follows from the Hessian bound 3, the displacement
1.5rho, the derivative change .2rho and half of the second-map bound.
All arguments stay inside R*=1/30. No analyticity of a cutoff is used.

Define the real coordinate extension and cap metric by

    G_b(w)=w+zeta(|w|)(F_b(w)-w),
    u_b(w)=chi(|w|)sigma_b(w),
    cap_b=(G_b)_*[exp(2u_b)|dw|^2].                 (5)

For small real b, G_b is a global diffeomorphism equal to the identity
outside R0. The cap agrees with the seed wherever zeta=chi=1. Align
the pushed spin frame with the seed frame there, using the spin rotation
of F_b', extended by zeta. Write this rotation Q_b. At b=0, G_0=id,
Q_0=I and the cap is exactly the fixed normal cap used for the previous
row estimate. Every cap derivative at b=0 vanishes on |x|<1/64.

This common physical family represents the original point-centered UV
cap at its own center: the embedding is an isometry and sends 0 to b.
Changing its embedding in the transition does not change that centered
trace. It does change U(x,b) away from x=b, and those derivatives are
precisely what is included below.

## 3. Finite-jet coefficient control, including the density term

Represent spinors as half-densities in the fixed physical Lebesgue
space. The massless cap Dirac operator has the two-dimensional form

    A_b=P_b^i partial_i+(partial_i P_b^i)/2,
    P_b^i(x)=Q_b(x) exp(-u_b(w)) gamma^a
                          (DG_b(w))^i_a Q_b(x)^-1,
    w=G_b^-1(x).                                   (6)

The unitary coordinate transport contributes half the derivative of
log det DG. Direct differentiation of (6) retains it. Conjugating by
the Spin(2) rotation introduces no extra zero-order term: its generator
anticommutes with both planar Clifford matrices, so the additional
anticommutator is zero. This also verifies that A_b equals the original
seed operator throughout the common inner neighborhood. Conformal
covariance is the dimension-two case of
[Walpuski, Spin Geometry, Proposition 15.5](https://walpu.ski/Teaching/SpinGeometry.pdf),
printed page 79 / PDF index 78; the transport formula (6) is derived
directly here, not imported as a numerical estimate.

The cutoffs are C4, not analytic. To bound only the needed parameter
jets, fix a real spatial x and replace each cutoff by its Cartesian
Taylor polynomial through degree four at x. Use this polynomial in
(5)-(6) for the local complex calculation. The resulting coefficient
jets through parameter order three, including one spatial derivative,
are EXACTLY the true jets: cutoff remainders in u_b are of order five,
and a cutoff remainder in G_b has the additional factor F_b-id=O(b).
Even after two spatial derivatives it cannot affect the third b jet.
This is a finite-jet device, not a holomorphic extension of the actual cap.

On |w-x|<=2rho, the cutoff tensor table gives

    |chi|<1.007, |Dchi|<330, |D^2chi|<7.7e5,
    |zeta|<1.004, |Dzeta|<162, |D^2zeta|<1.9e5.

For the local polynomial extension, (4)-(5) consequently give

    |DG_b-I|<250rho, |D^2G_b|<286000rho,
    |DG_b^-1|<1.003, |G_b^-1(x)-x|<1.6rho.          (7)

The inverse exists by contraction on the 2rho disk: its displacement
is below 1.004(1.5)rho and its derivative error below 250rho<.003.
The complex cutoff derivatives used above are Taylor sums through
order four, so these bounds remain valid for that inverse construction.

At real b=0, |u_0|<.002 and |Du_0|<.7. On the complex neighborhood
one may use |Du_0|<.75, |D^2u_0|<1700. Equation (4) gives
|u_b-u_0|<.303rho and |D(u_b-u_0)|<110rho before the inverse map.
After composition with G_b^-1, the scalar principal multiplier obeys

    |exp(-u_b(w))-exp(-u_0(x))|<1.6rho,
    |D_x[exp(-u_b(w))-exp(-u_0(x))]|<3200rho.

The extended spin rotation has |Q_b-I|<.21rho and |D_x Q_b|<45rho.
Using (6), including all inverse-map and rotation products, yields

    |P_b-P_0|<600rho,
    |div P_b-div P_0|<1.2e6rho.                    (8)

For the largest term in the second bound use
div_x DG=grad_w log det DG. Contracting D^2G against DG^-1 and the
two Clifford matrices costs at most 2(1.003)(1.003)(286000)rho.
The scalar-gradient difference contributes below 10000rho, and the
remaining scalar/matrix and rotation products below 10000rho. Their
sum is below 1.2e6rho. The principal bound follows by twice the sum
of 1.6(1.003), 1.003(250), and 2(.21)(1.003)^2, all times rho.

Let ||v||_EM=||partial v||_2+M||v||_2. Since M>=18, (8) bounds the
local polynomial operator difference by

    [600rho+1.2e6rho/(2*18)]||v||_EM <.35||v||_EM.

Cauchy in the two b coordinates, and the Hilbert tensor norm over
ordered derivatives, therefore prove for the TRUE cap vertices

    ||partial_b^j A_b|_0 v||_2 <=K_j||v||_EM,
    K_j=j!10^(5j), j=1,2,3.                        (9)

Indeed .35(sqrt(2))^j<1 for j<=3. The cutoff replacement was made
only to bound the exact finite jets pointwise; it does not define a
different quantum state or require a global complex inverse.

## 4. Resolvent graph rows and the complete angular insertion sum

At b=0, A_0=exp(-u_0/2)D exp(-u_0/2) is skew-adjoint on its cap
domain. Thus ||S||<=M^-1 for S=(A_0+M)^-1, and

    ||partial S||<2.1, ||S||_(L2->EM)<4.            (10)

To see the derivative bound, Dv=exp(u_0)A_0v+(Du_0)v/2;
exp(.002)*2+.7/(2*18)<2.1. This uses no curvature bound for the
transition. The same statements hold for the adjoint.

Keep the fixed physical row constants from
[the coercive interior estimate](seed-middle-coercive-jets.md):

    c0=30000,
    c_j=520000(j+3)[600(j+3)]^j, 1<=j<=3.

Each cap vertex in (9) is supported outside |x|=1/64. The radius-1/256
observation ball therefore leaves the same separation >1/88 used
there. Plain rows have bound c_j exp(-M/176). Their input-gradient
rows have bound 40M c_j exp(-M/176). Here is a graph-row justification:
conjugation by a distance weight with Euclidean gradient at most M/2
adds norm at most exp(.002)M/2. The weighted inverse is below 2.01/M,
and its flat gradient below 4.1. Duality then bounds the weighted
inverse followed by divergence below 5.2. The original local Fourier
evaluation calculation has prefactor below 12000M, including its
small-ball weight and endpoint factors. This is below 40M times its
plain-row allowance. Apply the same finite interior recurrence to this
Hilbert-valued incoming-gradient row to obtain all j<=3. No derivative
of a distant cutoff is assigned to a pointwise kernel bound.

The ordered inverse derivatives from the splitting audit, (9) and
(10) give the insertion factors

    V1=K1,
    V2=K2+8K1^2,
    V3=K3+24K1K2+96K1^3.                           (11)

There are two endpoint rows in every term. Intermediate graph inverses
cost at most four. The final vertex acts on the right row, costing at
most 41M times its plain constant, and the trace itself contributes M.
The spin trace costs two. This accounts for the factor 82M^2 below;
no unsupported cancellation of these two mass powers is used.

For the exact paired spectrum, an integral plus twice the maximum of
the unimodal summand gives

    sum d_n M_n^2 exp(-M_n/88) <25 C 88^4.          (12)

The integral is exactly 24C*88^4. An upper bound for the error is
4q(176/e)^2+8sqrt(C)(264/e)^3. Using q/C<361, 1/sqrt(C)<1/200000
and e>27/10 makes the error less than C*88^4. The first two guards
follow from 18<M1<19 and q=2000^4. Unimodality follows because the
derivative of 4CM/(q^2+4CM^2)+2/M-1/88 is negative: after multiplying
by M^2 its first contribution is at most 1/8, whereas the second is -2.
The whole angular sum is retained and absolute convergence justifies
the differentiated cap insertions and common subtraction.

## 5. Actual cap correction and a consistent complete source estimate

Use U(x,b), G_fix and C(x)=U(x,x)-U(x,0) in common physical coordinates
as in the splitting audit. Every local counterterm has zero b derivative
on the observation neighborhood. All first through third cap derivatives
of the finite remainder are thus given by the separated resolvent
insertions just bounded. Set

    P_k=sum_(a=0)^k binom(k,a)c_a c_(k-a),
    kappa_j=82*25*88^4*(31/2500000)/(j!100^j)
                *sum_(r=1)^j binom(j,r)V_r P_(j-r), j=1,2,3.       (13)

Equations (9)-(12) and the mixed derivative identity prove
|partial^j C|/(j!100^j q^2)<=kappa_j. These are now bounds for the
ACTUAL selected centered cap family, rather than hypothetical budgets.
They need not meet the earlier proposed 4e21 budget for the sharper
1e22 target.

For completeness do not simply reuse parameter Cauchy as a physical
derivative estimate for U_mov. The normal-cap family in its own w chart
is genuinely holomorphic in the physical b bidisk of radius rho: chi is
fixed, while (3) is analytic. Its complex normal factor is below .003,
vanishes on both axes, and has quadratic reserve c2=8. The retained
flat Born argument therefore has uniform remainder below 1e20 M^-4.
One may check this directly with degree-six remainder coefficients
94(.003)32^6 and 72(.003)^2*32^6 for Born one/two, and the full higher
sum 2949120000*8^3/(1-.015); their combined coefficient is below 1e20.

After the complete angular sum, the analytic local remainder is below
2e13 q^2. The common finite polynomial restoration has an analytic
variation allowance below 1e12 q^2: R=-Bcal'' has magnitude below eight,
Delta R=-Bcal Bcal''''-Bcal'Bcal''' below 150, and its fixed zeta/scale
factor below 1000C. The LLL curvature term fits this reserve. The common
a0 vacuum term is position independent and has zero derivative; keep
its previously bounded value at order zero. No scale is reset to B.

Consequently for j=1,2,3 a sufficient coordinate bound is

    |partial^j U_mov|/(j!100^j q^2)<=L_j,
    L_j=3e13*1500^j.                               (14)

The factor 1500 is larger than sqrt(2)/(100rho). This explicit center
transport replaces the old unexamined inference between parameter and
physical derivative bounds. Let G_j denote the previous fixed-cap
gluing ledger. At the normal center the third covariant conversion
adds the first-coordinate bound divided by 200. Therefore

    B1=G1+L1+kappa1, B2=G2+L2+kappa2,
    B3=G3+L3+kappa3+B1/200,
    max(1e14,B1,B2,B3)<1e26.                       (15)

The zeroth bound is the retained complete point-value estimate. Exact
rational verification of (13)-(15) completes (1). No sign of C is
assumed; it is subtracted exactly before absolute values are taken.

## 6. Consequence and verification scope

This supplies a complete middle source-three input at the original
1e26 allowance with one consistent cap family. It closes the specific
scope gap in the earlier middle ledger. The withdrawn 1e22 claim and
its proposed metric improvement are not reinstated. The existing
middle metric ledger may reuse its old K_m once its other hypotheses
are checked; this is not a new endpoint C5 or finite-metric H6 estimate.

Freeze the proof, checker and review before execution. Verify the
cutoff endpoint/radial bounds, Goursat and map budgets, finite-jet
replacement order, exact half-density transport and spin-rotation
anticommutator, all vertex/graph/insertion/angular budgets, and the
actual source aggregation. Use an independent constant-curvature
normalizer to check the coordinate formula, and a nonlinear compact
coordinate transport to check its density term. Independently compare
the angular moment integral and a guard-point complete-tail enclosure.
Fixtures do not evaluate the charged source or replace the ordinary
continuum domain, graph-row and coefficient arguments.

No Lean run is selected. Retain every scientific run, including failed
snapshots if any. No cancellation target, signed zero, nonlinear
existence or physical stability is established by (1).

Remaining to-do list: verify source-dependent metric estimates against
the repaired input; complete higher local/cold/outer and finite-metric
source/transfer bounds; cancellation, signed compatibility and nonlinear
existence; physical validity, state/boundary, stability, matching,
charge/return, opening, transport, preparation, separate creation and audit.
