# Finite cold admission: geometric bounds and an unresolved source transfer

2026-10-22, E003/S167. Ordinary argument before frozen supporting checks.

The actual-map and cap calculation advances the geometric inputs of
[the cold transport audit](finite-cold-neighborhood-transfer.md). It does
NOT close cold admission. In particular, the proposed enlarged endpoint
source proof has an unresolved separated-row and center-transport step.
Sections 4-6 isolate that step, calculate its exact numerical targets,
and select a fixed-domain method to resolve it. No conditional response
constant below is adopted as a solver certificate.

Keep the original B,A,C,q,mu0, temperature, paired domain and source units.
The actual geometry is sigma_ell=sigma0+ell u, 0<=ell<=1. The retained
finite endpoint contracts and P5_40(u)<1.3e-4 are unchanged. Write
epsilon0=1e-6 ell, a_x=s[1+log(.01/s)] and L_h=1+log(.01/h). We use
a_x<.00331, a_x<=2hL_h on chi<=1 and exp(-80)<h<1e-30.

## 1. An explicit fixed linear direction norm

All physical fields below live on the original common manifold. They are
real C5 on the cold/interior charts, stationary and spherical, regular
even fields at the bolt, and retain the outer clock and fixed projector.
At the outer boundary retain the existing lower regularity/weighted
expansion conditions. No compact-support restriction on a correction is
introduced. Independent radius and correlated mass conditions remain part
of the full metric problem, not of this scalar theorem.

For a fixed ell and each cold center x, freeze the seed normal chart F0,x
and the actual finite endpoint chart F*,x. The latter is obtained by the
real quadratic normalizer of ell u. Neither map is varied when defining
the norm of a field v. Use the original j!100^j derivative units. Define
D_ell(v) as the maximum of these LINEAR-field seminorms:

* P5_40(v), with the original fixed interior charts, local caps and actual
  global low norms of [the fixed-chart theorem](comparison-path-interior-response.md).
* Cold chart jets of v pulled back by F0,x and by F*,x, through order five,
  divided by a_x j!100^j. Include their one-center derivatives through
  order four divided by a_x j!100^(j+1), where center differentiation
  means the fixed seed H2 transvection used in the retained proof.
* |v|/a_x and |v_chi|/(8a_x) on the cold region, for the Ward kernel.
* Actual global value/gradient norms and the global value/gradient norms
  of the following fixed, explicit H2 extension Ev.
* The weighted suprema of |v|, |grad v| on the actual manifold and of
  |Ev|, |grad Ev|, |X_x Ev| on the H2 cap, divided by
  a_x exp(d_x(y)), uniformly over cold centers x and points y. Here d_x
  is the fixed seed distance used on each side of the retained gluing
  proof, and X_x is its unit radial transvection. These are real fixed
  coefficient functions, not resolvents or desired source estimates.

For Ev use the fixed seed H2 identification through the auxiliary seam
chi1 corresponding to s=.01. Set Ev=v there. For t=chi-chi1 in [0,1],

    Ev=eta5(t)[v1+v1' t+v1'' t^2/2],
    eta5=1-10t^3+15t^4-6t^5,

and set Ev=0 beyond chi1+1. The primes are fixed chi derivatives at the
seam. This linear extension is C2 at both joins and has the SAME clock.
The distant inverse requires C1, while the local C5 fields are unchanged.
It is the width-one C2 cap, not the narrow C3 cap used for seed jets.

This defines D explicitly without an unknown quantum condition. For
radial fields X_x Ev=(Ev)' cos(phi), so its weighted bound is directly
checkable. All charts, distances and extension maps are fixed once for
ell. At the fixed positive h the weighted suprema are finite for smooth
compatible fields with the stated endpoint behavior; no extremal limit
is asserted. The extra weighted/cap entries are genuine correction
requirements. They are not inferred from P5_40 or discarded by a norm
rescaling. An inverse-to-D embedding remains to be proved.

## 2. Actual finite maps and their center derivatives

The retained seed normal map satisfies, on |z|<=7/32,

    |F0'|<1.05, |F0''|<2.5, |F0'''|<300,
    |F0,x|<=260a_x, |F0,x'|<=1400a_x,
    |F0,x''|<=1.2e5a_x.

These are the actual analytic Goursat-map bounds in
[the cold source proof](stationary-reference-cold-bolt.md), retained for
the centered seed by [its source transfer](stationary-centered-exterior-source.md).
Extend a map F by G=id+eta(F-id), with eta one through 5/32 and zero
from 7/32. Its radial derivative operator bounds are 30,1500,3e5.
For the seed, ||DG0-I||<.439, ||DG0^-1||<=2 and ||D2G0||<22.

The actual low metric bounds supply |u|<1.2e-8, |grad u|<7.7e-6.
Their cold refinements follow from the real barrier in
[the finite metric proof](stationary-completed-metric-response-norm.md):
W<=2a_x, |du|<=1.2e-6 s(L_c+2), |Hess u|<=3.6e-6 s(L_c+2).
Since L_c>=1+log 10, fixed H2/frame conversion gives at the center

    |u_x|<4e-6a_x, |u_z,x|<2e-5a_x.               (1)

For the second bound include the derivative of the seed frame: its
coefficient is bounded by 200a_x; multiplying by the preceding gradient
adds less than 2.7e-6a_x. The Hessian contribution is below 6e-6a_x.
Thus (1) has slack. These require two derivatives, not an unproved sixth.

Set H_ell(z)=a_ell z+c_ell z^2 with
a_ell=exp[-ell u(0)] and c_ell=-a_ell^2 ell u_z(0). Then

    |a_ell-1|<2e-8, |c_ell|<4e-6,
    |a_ell,x|<5e-6a_x, |c_ell,x|<3e-5a_x.

The finite map is F*=F0 composed with H_ell. The image retains the
original analytic seed margin. Finite chain/product rules give

    |F*'|<1.052, |F*''|<2.7, |F*'''|<310,
    |F*,x|<261a_x, |F*,x'|<1401a_x,
    |F*,x''|<1.21e5a_x.                            (2)

In obtaining the last bound, the derivative F0'''_x is needed only
multiplied by H_ell-id if a difference is used. Avoid that route: compute
F*,x'' from F0,x''(H_ell), F0'''(H_ell) H_ell,x and the displayed
quadratic derivatives directly. The seed center bounds hold on the
slightly enlarged image by their strict .018 margin. No unknown-field
third center derivative enters.

For a real endpoint perturbation v with d=D_ell(v)<=5e-10, normalize it
in the F* chart by H_v=a_v z+c_v z^2. The new local map is F=F* composed
with H_v, and its global extension is G=id+eta(F-id). This composition
is used consistently on the actual and auxiliary operators. We do not
replace it by a different quadratic map in the seed chart.

Put t=d a_x<=1.655e-12. The normalizer has |a_v-1|<=2t, |c_v|<=100t.
For a direction w of norm U, differentiating its explicit coefficients
at any such real v gives, after extracting U a_x,

    |delta H|<6, |delta H'|<46, |delta H''|<200,
    |delta H_x|<600, |delta H_x'|<5000.             (3)

For example |delta a_x|<101Ua_x and |delta c_x|<11000Ua_x;
at R=7/32 these give 101R+11000R^2<600 and 101+22000R<5000.
The products involving v itself are below 1e-8 of these reserves.

Using (2),(3), the actual map-direction arrays BEFORE the cutoff obey

    |delta F|<6.4Ua_x, |delta F'|<66Ua_x,
    |delta F''|<2400Ua_x,
    |delta F_x|<670Ua_x, |delta F_x'|<10000Ua_x.

For the last line, the largest whole-cold terms are
121000(.00331)(6), 2.7(600), 1401(.00331)(46), and 1.052(5000);
their sum is below 10000. The H_v,x cross terms use t<=1.655e-12.
Applying the fixed cutoff gives

    |delta G|<11Ua_x, |delta DG|<420Ua_x,
    |delta D2G|<3e4Ua_x,
    |delta G_x|<1100Ua_x, |delta DG_x|<5e4Ua_x.   (4)

For the unvaried finite maps, direct substitution of |a_ell-1| and
|c_ell| into F0(H_ell)-F0 gives changes below 1e-5 in DG and .001 in
D2G. The extra v is much smaller. Hence, uniformly on the real tube,

    ||DG||<1.5, ||K||<3, ||D2G||<23,
    K=(DG)^-1, |grad K|<=9(23)<210,
    |G_x|<262a_x, |DG_x|<9400a_x.                (5)

The global map is a diffeomorphism: its extension differs from the
identity by a map of Lipschitz constant <.45, and is the identity outside
the fixed disk. It preserves the physical outer realization. At the bolt
the physical conformal family is even and radial; coordinate transport
does not introduce a new boundary or change the period.

## 3. Geometric cap matching and local transported-coefficient budgets

Use sigma_cap=sigma_seed,cap+ell Eu+Ev. This equals the ACTUAL geometry
through chi1 for every endpoint, including its local density. This local
agreement does not identify the two global spinor domains.
For the original finite u the seam bounds are

    |u1|<3e-9, |u1'|<4e-8, |u1''|<2e-7.

They follow from the same low barrier/curvature equation at s=.01.
The width-one polynomial gives |Eu|<3e-8, |(Eu)'|<6e-7. For v the norm
explicitly supplies |Ev|,|grad Ev|<=d. The seed cap has value <.042 and
gradient <.4. Thus on the ENTIRE cap

    |sigma_cap|<.043<.05, |sigma_cap'|<.401<.5.    (6)

On the normalizer image, the seed has |grad Phi|<=4 and |Hess Phi|<=400,
where Phi includes the fixed H2 conformal factor. The actual low u jets
and D-controlled v jets enlarge these to <5 and <500. Consequently
nu=Phi composed with G satisfies exp(-nu)<1.1, |grad nu|<8.
Freeze the explicit conformal variation and handle it separately by
S_delta=(S f+f S)/2-M S f S. The COORDINATE variation therefore satisfies

    |delta nu|<100Ua_x,
    |grad delta nu|<21000Ua_x.                   (7)

Indeed 5(11)<100 and 500(1.5)(11)+5(420)<21000.
On the local image |grad Phi_x|<21000a_x: the retained seed bound is
2e4a_x and the actual u/v second-jet additions fit its spare 1000a_x.
On chi<=1 the mixed products carrying two a_x factors are bounded with
a_x<1e-27. The largest remaining term in delta nu_x is
5|delta G_x|<5500Ua_x. The explicit cap multiplier's center derivative
is controlled by the fixed weighted norm, and is not counted again in
nu_delta. Equations (2)-(4) give |delta nu_x|<10000Ua_x and
|delta DG_x|<1e5Ua_x.

The background center derivative satisfies |nu_x|<10: a radial cap has
X_x sigma_cap=sigma_cap' cos(phi), bounded by .5, and its local coordinate
term is <5(262)a_x near the bolt. Also |K_x|<=9(9400)a_x<1 there.
These are cap/local coefficient inputs. Extending them to a global
center-dependent ACTUAL-domain realization requires the argument in
section 4; local coefficient bounds alone do not verify that argument.
The fixed weighted entries of D supply the global exp(distance)
envelopes for the explicit multipliers on BOTH actual and capped domains.
All remaining coordinate coefficients have compact support and (4).

The exact symmetric-divergence formulas in the preceding audit now apply
to these local/cap arrays. They give the following algebraic budgets,
with U a_x extracted,

    |delta B|<10230,
    |grad delta B|<21360240,
    |delta B_x|<1132010, chi<=1.                 (8)

Each is less than twice its old reserve. The symmetric form includes the
Jacobian connection. A local spin-frame rotation expresses the cap kernel
with constant Clifford matrices for the Born calculation. It is unitary
and retains the traced diagonal. Undo it before common-frame gluing:
the graph estimates above use the unrotated half-density frame on BOTH
operators. No extra derivative of a spin rotation is silently charged
to a different graph norm. The spin connection is not omitted or added
as a second source.

## 4. Where the proposed closure stops

The local endpoint calculation retains the preceding coefficient 122.
The explicit D norm also makes its cap-field requirements inspectable.
These are genuine added hypotheses on a future correction. They do not
prove that an inverse maps into D, or improve a bound by rescaling units.

A proposed continuation would use the global coefficient interpolation

    B_tau=B_0 exp[tau log(B_0^-1 B_1)], 0<=tau<=1.       (9)

Near the observation the two frame matrices are scalar times rotations;
there (9) interpolates the normal conformal factors and rotation angles.
Outside the compact transport they are conformal in the fixed chart.
In the transition annulus the matrices need not commute. A small relative
matrix norm, orientation and an inverse bound alone are NOT a proof of
the complete differentiated point-row theorem for this interpolation.
The physical paired domain, local spin rotation, divergence connection,
spatial evaluation, and parameter/center derivatives must be kept together.
Undo any local spin rotation before gluing. In two dimensions the identity
{gamma_i,gamma_1 gamma_2}=0 explains its symmetric-divergence cancellation;
it does not supply a bound for derivatives of a separated resolvent row.

The proposed constants 800 exp(-M/32) for the local commutator row and
5e20,1e30 for the far-seam value/gradient prefactors have NOT been derived
for this entire actual/cap construction. The first is not justified by
doubling a previous 400 bound after changing the interpolating operator.
Likewise, a weighted L2 inverse does not imply a point evaluation bound
without local regularity and cutoff estimates. The old one-path theorem
cannot simply be applied to every member of a new D ball.

There is a concrete domain distinction to enforce. In the unit disk let
the physical boundary have radius R=exp(-h)<1 and consider the H2
transvection T_b(z)=(z+b)/(1+bz), for real b. At b=0 its normal velocity
on z=R exp(i theta) is

    Re[exp(-i theta) partial_b T_b(z)]
                      =(1-R^2) cos(theta).             (10)

This is nonzero for finite h. Thus this isometry of the COMPLETE H2
disk does not preserve the finite physical disk or its boundary projector
without a further identification. Its small boundary velocity is not zero.
Equation (10) is a counterexample to that automatic identification, not
an obstruction to a correct fixed-domain proof. A compactly supported
normalizer G does preserve the boundary; a complete center transvection
and G are different maps. Any proof using both must account for both.

The decisive domain calculation (10) and the threshold calculation below
are executed in the frozen checker. They reject promoting the current
local/cap calculation to cold admission. They do not reject the earlier
outer component, or prove that the desired cold estimate is false.

## 5. Exact missing thresholds, including the far seam

The proposed endpoint tube still needs

    epsilon_end <=1.000001e-6 + K_tr d <=2e-6,
    d<=5e-10, hence K_tr<=1999.998.                     (11)

K_tr is an UNKNOWN complete transport/source-contract constant. The
local bound 122 is known; it is not a known upper bound for K_tr. A
candidate 1000 fits (11), but has not been established for the full
contract. Local potential and curvature tests alone do pass:
exp(.025+2e-6*.00331)-1<.026 and
2!100^2(2e-6)(.00331)<.000133.

The proposed auxiliary graph test is also only part of the argument.
If the chord coefficient differences obey the specified original arrays,
b0=1e4(2e-6)(.00331), b1=2e7(2e-6)(.00331), its ratio would be

    13 b0 +(b1/2)(1.5/18)+(3/4)b0(1.5)<.007.           (12)

This includes the divergence connection. It is a passing CONDITIONAL
budget, not proof of the actual weighted point-row or common-domain
hypotheses. An L2 Neumann pass is not a substitute for those hypotheses.

The explicitly matched cap gives useful separation. The cold lower guard
.946 log 10-.25>1.928 leaves D>=1.9 after the finite value perturbation.
At chi<=1, chi1=arcosh(.01/h), h<1e-30, and the quarter-unit seam cutoff
give D>=60. Do not assert a long separation on the whole cold patch.

Define UNKNOWN prefactors T0,T1 by sufficient per-mode far-seam estimates

    |delta E_M,far| <=T0 epsilon a_x exp[-(M/2-1)D],
    |grad delta E_M,far| <=T1 epsilon a_x M
                                      exp[-(M/2-1)D].  (13)

The exponent in (13) is itself a target: paired separated rows must
justify it, including the one exp(D) insertion weight. It is not inferred
from geometric distance alone. To allocate at most 1e20 epsilon a_x M^-3
and 1e20 epsilon a_x M^-2 respectively, the EXACT sufficient thresholds are

    T0 <= 1e20 exp(15.2)/18^3,
    T1 <= 1e20 exp(480)/18^3.                           (14)

Indeed M^3 exp[-(M/2-1)D] decreases for M>=18,D>=1.9,
since its logarithmic M derivative is 3/M-D/2<0, and it also
decreases in D. The maxima are at the endpoints (18,1.9) and (18,60).
The suggested 5e20 and 1e30 satisfy (14) arithmetically. They remain
unproved upper estimates, rather than known estimates that fail.

For the LOCAL flat-to-H2-cap correction, suppose the full mixed ledger
is C_glue epsilon a_x R^2 exp(-M/16), with C_glue=1e10 as a target.
The exact maximum M^3 exp(-M/16) on M>=18 is 48^3 exp(-3).
Allocating 1e20 of the total 1e22 reserve therefore requires

    C_glue R^2 <=1e20 exp(3)/48^3.                     (15)

At C_glue=1e10 this allows R<=sqrt(1e10 exp(3)/48^3).
R=800 passes. Neither C_glue nor R has been proved on the new chord.
This reveals useful slack; obtaining another very small L2 inverse is
not the decisive missing task. Derive the differentiated row argument.

The full endpoint theorem must then provide, including the common finite
restoration, all angular modes, and the actual paired inverse,

    |S_endpoint-S_seed|<=1e23(2 eta_G)epsilon a_x,
    |delta S(chi)-delta S(0)|
           <=1e24(2 eta_G)epsilon hL_h sqrt(chi).       (16)

The counting bound N(E)<=4CE^2 implies sum d_n/M_n^3<=2C/3 when
M1>=18; this summation is legitimate only after a uniform per-mode
estimate. Finite-mode numerical samples cannot establish (16).

IF (11) holds with K_tr<=1000 and (16) holds through epsilon<=2e-6,
rho=2.5e-10 would admit a real background radius rho and parameter
radius rho. With the independent fixed-chart complex bound 5e-11 and
real amplitudes 3.4e-15 a_x, 3.4e-14 hL_h sqrt(chi), the earlier
half-disk calculation would give first/mixed constants .007,2e10,
and bolt oscillation constants .07,2e11. These remain CONDITIONAL.
The Ward/lapse mixed allowance 8e14 would give remainder 4e-6 at
D(h)=1e-10. That is not a preconditioned Newton remainder.

## 6. Full feasibility ledger and research decision

Retain [the quantitative construction audit](construction-feasibility-audit.md)
and its exact unrounded input ledger. The latest complete interior shape
coefficient is C_I=9.1591500000000000015e103 (upper enclosure 1e104);
the global outer-dominated coefficient is C_b=2.365283519093731664e186.
With C_sigma=1421.40132786776365 and C_b4=4.52219350981899516e14,
the tight forcing recombination is about 1.069626977892750664e201;
the interior-only counterpart is about 4.14194486854586495e118.
The retained physical g is approximately 9.765625e-34. Uniform norm
lower bounds R0 R2 R4>.97 and C_D4>2 mean those known upper allowances
miss the old source budget by more than 4.4868e225/epsilon_c globally
and 1.7374e143/epsilon_c even in the interior diagnostic. The low mixed
allowance independently fails. These are certification failures, not
lower bounds on the actual correction and not nonexistence theorems.

For a COMPLETE solver space X the preserved-ball and projected-solution
conditions remain

    Y0+Z0 r+Z2 r^2/2<=r,  Z0+Z2 r<1,
    Y0=||B Pi F(0)||_X, Z0=||I-B Pi F'(0)||,
    ||B Pi(F'(h)-F'(k))||<=Z2 ||h-k||_X.               (17)

B must be injective on the projected source range. F must be C1 on the
whole admitted ball with the correct physical domains. For r=1e-8,
sufficient numerical targets are Y0<=5e-9,Z0<=.25,Z2<=5e7. All three
are still unknown in such a complete X. The established outer canonical
margin is >4.9e-8; the interior allowance leaves 8.7e-4. They require
E_o r<=4.9e-8 and E_45 r<=8.7e-4. A conditional cold radius 2.5e-10
would additionally require E_c<=.025 at that r. E_c is unknown, and the
cold radius itself remains unproved. Radius positivity, the common metric
boundary action, remaining source sectors and all full-X embeddings stay
separate. Changing units leaves response-times-radius products invariant.

Decision: pursue a SPECIFIED sharper method, keeping this construction
under investigation and the failed uniform-product route parked. Work in
the original physical common chart for the actual-minus-H2-cap identity,
with a fixed physical seam cutoff. Keep the actual inverse and its paired
domain fixed under observation differentiation; differentiate the separated
cap rows and their evaluation point. Treat coefficient variations on that
same domain. Do not transport the actual boundary by a complete H2
transvection. Use the explicit matched cap above, retain its weighted
insertions, and prove (13)-(16), with all finite restoration terms.
Then combine with the local source argument or derive weighted parameter
responses directly, without claiming that a changed normalization is a gain.

After admission, enclose the ACTUAL correction and signed moment jointly,
using low physical frequencies and a certified fourth-log tail. For the
illustrative epsilon_c=1e-4 allocation, targets remain |beta|<=1e-53 or
N(q3)<=1e53, including all residual errors. This is distinct from metric
admission and from solving the exact nonlinear compatibility equation.

Compatibility and physical filters can change the ordering: on a projected
branch z_i in [.25,.28], a derivative of fixed sign with lower magnitude m
and midpoint error e requires m>(|eta(mid)|+e)/.015 for a sign bracket.
No signed data are supplied. The retained fast frequency obeys
k sqrt(G)=sqrt(60pi+17g/12) approximately 13.7294. Small background
curvature does not establish semiclassical validity there. Specify the
metric boundary/action and validity errors alongside the direct solve,
with ||B Pi delta F||_X<=r/4 as an explicit omission-error target.
The off-shell constraint obstruction is not physical stability; the
on-shell constrained causal problem remains open. A demonstrated boundary
or validity incompatibility would justify switching construction. The
present oversized estimates do not.

The audit executes the decisive tests for this proposed admission proof:
the maps and cap have plausible numerical room, while the full source
claim is not supported by its current continuum argument. Feasibility
of an actual controlled solution remains undecided. The next calculation
is now the fixed-domain differentiated resolvent identity, not another
unrelated estimate or another checker pass.

## 7. Frozen verification scope

Freeze this note, checker and review before execution. Check finite-map,
cutoff and cap arithmetic; exact admission and source thresholds; the
finite-disk domain diagnostic; matrix-order and spin-frame identities;
and the conditional graph/response budgets. Use independent 60/90-digit
evaluations where applicable. Preserve prior hashes and every run. These
checks are neither actual quantum-state samples nor machine verification
of the continuum proof. No Lean is invoked. No full cold-admission flag
or nonlinear acceptance gate is promoted by a passing run.

Remaining to-do list: fixed-domain cold quantum-source proof and admission;
actual correction in admitted norms and signed moment; complete nonlinear
map, projected solution and exact compatibility; physical boundary/state,
validity, matching, stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
