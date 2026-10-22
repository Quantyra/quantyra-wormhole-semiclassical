# Local cold gluing and finite-centered canonical admission

2026-10-22, E003/S167. Ordinary argument before frozen supporting checks.

The [distant contribution](fixed-domain-cold-far-seam.md) is already
controlled on the fixed physical domain. Here the remaining local
flat-to-H2-cap gluing is estimated along the actual REAL conformal family
on complete H2. The flat local source is compared separately at its
normalized endpoints. There is no global matrix-log interpolation and no
parameter differentiation of five spatial jets of a moving C5 field.

Together these estimates give a massive-canonical cold neighborhood in
the explicit D_ell space of [the audit](finite-cold-admission-and-response.md).
They do not put the unknown correction in D_ell or solve the full
metric/radius equations. The physical metric boundary/action and stability
requirements remain separate.

## 1. Three contributions and two different comparison paths

Fix ell in [0,1], d=D_ell(v)<=5e-10, epsilon0=1e-6 ell and
a_x=s[1+log(.01/s)]. The physical endpoint is sigma0+ell u+v. Its
explicit H2 cap is sigma0,cap+ell Eu+Ev, with the width-one C2 extension
and the same clock. The direction norm, fixed charts, weights and all
physical constants retain their previous definitions.

At each real endpoint decompose its canonical scalar into the flat normal
cap source, the H2-minus-flat gluing term and the actual-minus-H2 term.
Subtract the same three quantities at the seed. The last difference is
already proved. The flat source is compared along a LINEAR chord of its
two fixed normalized conformal factors. Only its endpoint C5 difference
and one-center C4 difference are used. Its complete local Born argument
does not require this chord to be a global physical interpolation.

For the gluing term take the real family sigma0,cap+ell Eu+t Ev,
0<=t<=1, and normalize each member. Let F*=F0 composed with H_ell be
the frozen finite endpoint normal chart. In that chart put

    H_t(z)=a_t z+c_t z^2,
    a_t=exp[-t v(0)], c_t=-a_t^2 t v_z(0),
    F_t=F* composed with H_t.

Extend F_t by G_t=id+eta(F_t-id), with eta one through 5/32 and zero
from 7/32. These are the actual maps already bounded in the audit.
They are used consistently on the complete H2 cap. Their spin-frame
rotation is undone before gluing to the transported H2 inverse.
The observation center is moved by an H2 isometry ONLY on this complete
domain. No such map acts on the finite physical boundary.

The local flat potential is V=kap(exp(nu)-1), where nu is the normalized
factor and kap is fixed radial, one through 3/16 and zero by 7/32.
The cutoff chi for the gluing identity is one through 1/8 and zero by
5/32. Thus both operators agree on its whole commutator support.
Earlier Taylor-tail arguments allowing a cap out to 1/4 also apply to
this smaller support. The unknown field is never evaluated outside its
proved chart.

## 2. Actual low derivatives of the normalized potential

Write a dot for t derivative and x for one spatial-center derivative.
The audit's map arrays give, with d a_x extracted,

    |dot H|<6, |dot H'|<46, |dot H''|<200,
    |dot H_x|<600, |dot H_x'|<5000.

On |z|<=1/8 use the sharper bounds 2,27,200 and 200,3000,22000
for dot H, dot H', dot H'', dot H_x, dot H_x', dot H_x''.
They follow directly from a_t,c_t and their center derivatives; the
complex first derivative v_z has modulus at most 100d a_x/sqrt(2),
leaving slack even after the tiny exponential factors.

The normalized finite background on this inner disk has successive
spatial bounds 1,16,400,2e4 through order four. For the analytic seed
these follow from .027/(1-8T)^2; the finite comparison adds
epsilon0 a_x j!100^j. Its largest fourth-jet addition is below eight.
The tiny H_t displacement and t v additions fit the same reserves.
On the larger support retain the earlier gradient bound five.

Differentiate nu_t=(nu_*+t v)(H_t)+log|H_t'| only to spatial order
two, and once in the center. The inner second-derivative ledger for
dot nu is bounded by

    20001+400(2)+2(16)(27)+200+1 < 24000

times d a_x. The unit at the end covers the log term and products
containing d a_x. At zero both dot nu and its gradient vanish because
each endpoint is normalized. For |z|>=1/8 the direct value ledger is
below 140d a_x; division by |z|^2 costs at most 64.
Exponentiation and the fixed cap consequently give

    |dot V(y)|<=alpha d a_x |y|^2, alpha=2e4.           (1)

For the mixed derivative the principal inner terms are the direction's
mixed second jet 2e6, and
400(200)+2(16)(3000)+22000=198000. After Taylor division by two,
exponentiation and the small products, use

    |dot V_x(y)|<=beta d a_x |y|^2, beta=3e6.          (2)

Outside 1/8 its value ledger is below 8000d a_x; multiplication by
64 fits (2). The fixed-seed analytic center derivative is O(a_x):
Cauchy on its Goursat parameter margin .0089 and spatial margin 1/8,
with physical parameter derivative sum <=3s, bounds its second jet by
1e5 a_x. The finite e_ell mixed second jet adds at most 2a_x.
The perturbation adds at most beta d a_x. Thus the generous background
bound

    |V_x(y)|<=b a_x |y|^2, b=1e7                     (3)

holds near the bolt. There a_x<1e-27, so products of two such factors
are explicitly negligible below; they are not discarded elsewhere.
Equations (1)-(3) use at most four derivatives of v. They do not assert
C5 Lipschitz dependence under moving coordinates.

The whole cap has |V|<delta=.026 and |V(y)|<=10|y|^2. The inner
quadratic coefficient follows from the second-jet bound 16 and
normalization; outside 1/8 it follows from delta/|y|^2<=1.664.

## 3. Flat resolvent rows, including the missing derivation

Let G_M=(D_flat+M)^-1 and T=D_flat+M(1+V). Use the weight
w(r)=(M/4+2)r, M>=18. The radial operator-norm envelope for G_M is
k_M=M(K0+K1)/(2pi). Its weighted L1 norm is <7/M. To verify this
without a sampled mass cutoff, put a=1/4+2/M<=13/36. The positive
[Bessel integral representation](https://dlmf.nist.gov/10.32.E9) gives

    M ||e^w k_M||_1
      = integral_0^infinity (1+cosh t)/(cosh t-a)^2 dt
      = integral_0^1 4/[(1-a)+(1+a)y^2]^2 dy
      <=36/23+(22/7)1296/(23*33)<7.                  (4)

Here atan is bounded by pi/2, pi<22/7 and sqrt(1127)>33.
The same positive representation gives K0<=K1 and
K1(s)<=sqrt(2)exp(-s/sqrt(2))/s. With Gamma(2/3)<2 and
1/sqrt(2)-13/36>1/3 this implies weighted L^(4/3) norm squared
<22/M<56/M. It also gives

    ||r^2 G_M e^w||_2<10/M^2.                        (5)

Indeed its squared norm is bounded by 24/[pi M^4 b_M^4], where
b_M=sqrt(2)-1/2-4/M>2/3. The free far row at r>=1/8 has weighted
norm <1.5: its square is <=(4/pi)E1(M b_M/8), with
M b_M/8>1.5 and E1(q)<=exp(-q)/q. For the gamma bound, split
its defining integral at one to obtain Gamma(2/3)<=3/2+exp(-1)<2.

The full Neumann inverse therefore has weighted norm <9/M since
7delta<1. Its nonfree row has norm
<=56delta/(1-7delta)<1.8, by Young's inequality with the two
L^(4/3) kernels. This includes EVERY potential insertion.

For the canonical commutator row r_M=ev_0 S_flat C_chi, normality gives
nu(0)=0 and

    r_M=ev_0 T^-1 J, J=exp(-nu/2) gamma.grad chi.

The commutator slope is <33. The retained distance conversion on its
annulus has exp[2(d_H2-|y|)]<1.02. All endpoint/frame allowances fit
an extra factor two. Extract exp(-M/32) from r>=1/8. Then

    ||W^-1 r_M||_2 <= R_M=400 exp(-M/32).             (6)

The explicit product 2(1.02)(33)(1.5+1.8) is below 400.
The same holds for the adjoint column. This is the flat resolvent proof
of the row bound; no row constant is assigned by doubling an old result.

For any quadratic insertion |q(y)|<=A|y|^2, the exact identity
T^-1=G_M-G_M M V T^-1 and (5) give

    ||ev_0 T^-1 q e^w||_2
       <=[10+(100)(9)/16] A/M^2 <70 A/M^2.           (7)

This is what controls pointwise parameter differentiation. Use (7),
the inverse 9/M, J<33 and
dot J=-(dot V/[2(1+V)])J on the commutator annulus. The normalized
row derivative costs at most .4 A R_M. For example its resolvent part
is <=2(1.02)(70)(9)(33)A/(400 M^2), and its endpoint part is
<=A(5/32)^2/[2(1-delta)]; their sum is <.4 A at M>=18.

Differentiating again with two quadratic insertions bounds the bilinear
row variation by 4AB R_M. The two ordered inverse products cost at most
twice the preceding resolvent coefficient times 9B/16; the two endpoint
products and the second derivative of (1+V)^(-1/2) fit the remaining
reserve. The direct mixed potential insertion costs .4 beta.

Undo the local spin rotation before gluing. Its first, center and mixed
derivatives are bounded by 100d a_x, 2000a_x and 1e5d a_x respectively,
from F_t', dot F_t', F_t,x' and dot F_t,x'. These are included below.
The Clifford identity {gamma_i,gamma_1 gamma_2}=0 gives the corresponding
symmetric-divergence expression; it is not permission to omit derivatives
of the rotation. The complete common-frame row bounds are

    ||W^-1 dot r|| <=1e4 d a_x R_M,
    ||W^-1 r_x|| <=1e7 a_x R_M,
    ||W^-1 dot r_x|| <=2e6 d a_x R_M, chi<=1.          (8)

For the last line .4 beta+1e5<2e6; the remaining products are bounded
by 4 alpha b a_x and the displayed frame products, all <1 in these
coefficient units since a_x<1e-27. The factor a_x is retained.

## 4. The actual transported H2 inverse and its mixed insertion

Use a FIXED center-coordinate weight W=exp(-2d_H2(0,y)). Pulling it
through G_t gives physical slope at most six. More precisely,
||DG_t-I||<.45 and the local H2 conversion <1.05 give proper inverse
differential <1.05/.55<2, so using three leaves a strict reserve.
Thus the cap inverse obeys

    M ||e^w S_H2 e^-w||
       <=exp(.043)/[exp(-.043)-6/18]<2.               (9)

This accounts for the changed weight gradient; a slope-two bound is not
silently reused after transport. Dirac energy follows directly from
D_H2 T^-1=I-M exp(sigma)T^-1, with the weight commutator included.
H2 Lichnerowicz adds only curvature 1/2. Canonical/local flat conversion
gives energy <4.2 before G_t. The half-density pullback then costs
at most 1.5(4.2)+(3)(23)(2)/18<20. No curvature bound for an unrelated
narrow cap is used.

The real G_t are global diffeomorphisms. The cap has bounded C1
coefficients and a complete, fixed H1 domain. Its first-order expression
in the unrotated half-density frame is {B^a,partial_a}/2. The audit's
actual matrix arrays, including the Clifford factor sqrt(2), give

    |dot B|<=2e4 d a_x, |grad dot B|<=4e7 d a_x,
    |dot B_x|<=2e6 d a_x, chi<=1.                    (10)

The explicit conformal multiplier Ev is handled separately, once, by
S_dot=(S Ev+Ev S)/2-M S Ev S. Its value and center derivative have
the weighted envelopes in D_ell. Compact coordinate contributions in
(10) are not counted again as a conformal source.

For the background center variation, X sigma_cap=sigma_cap' cos(phi)
has value <.5 and gradient <20. The latter follows from the seed bound
|sigma_tt|<4 and the finite C2 cap extension: P5_40 bounds seam
|v_chichi|<3e4d, and differentiating the width-one polynomial gives
|Ev_chichi|<3e5d. At d<=5e-10 this adds <.00015. Regular evenness
controls sigma'/sinh(chi) at the bolt. These are coefficient bounds,
not unproved resolvent assumptions.

The coordinate-only center arrays have |B_x|<2e5 a_x and
|grad B_x|<1e9 a_x. They follow from G_x<262a_x,
DG_x<9400a_x,D2G_x<1e6a_x and the displayed inverse/product rules.
Their graph contribution is negligible for a_x<1e-27. The multiplier
identity and (9) therefore give uniform center inverse/energy bounds
10/M and 40. No derivative of a moving finite boundary occurs.

In the divergence-form bilinear insertion, half of dot B multiplies
each first-derivative energy term. The weighted inverse/energy constants
2 and 20 give 8e5 d a_x/M. The explicit conformal terms fit the spare
reserve. For the mixed insertion, differentiating both inverse factors
costs 5.6e6, and the direct dot B_x term costs 8e7, in d a_x/M units.
The conformal mixed terms fit below 1e3. Hence

    M||W S_dot W||<=1e6 d a_x,
    M||W S_dot,x W||<=1e8 d a_x, chi<=1.              (11)

Parameter differentiability here is H1-to-L2 coefficient differentiability
followed by the resolvent identity; it uses the C1 operator arrays and
common complete domain. It is not inferred from a manufactured matrix.

## 5. Close the local gluing budget

The exact diagonal identity is -r_L S_H2 r_R, with both commutators
already included in the rows. Retain rank two and the external M.
For the t derivative its coefficient is
2[2(1e4)(2)+1e6]<3e6. For the mixed center derivative the nine product
terms have coefficient, before rank/physical-center conversion,

    2(2e6)(2)+2(1e4)(10)+1e8
              +[2(1e7)(1e6)+2(1e4)(1e7)(2)]a_x.

A factor four covers rank and proper-center conversion. The result is
<1e9. Thus the finite-to-perturbed local gluing difference obeys

    |delta E_M,glue|,
    |grad delta E_M,glue| (chi<=1)
          <=1e9 d a_x (400)^2 exp(-M/16).             (12)

The maximum M^3 exp(-M/16) is 48^3 exp(-3)<5600, and the M^2
maximum is smaller. The contribution is therefore <9e17 d a_x in
the required M^-3/M^-2 units. It fits the 1e20 allocation with a
factor greater than 100, even before assigning 1000d to the endpoint
contract. The seed-to-finite local term retains the original real-path
bound 1e10 epsilon0 a_x R_M^2 from the complete-H2 argument; (6)-(7)
supply its previously implicit row mechanism. Its hypotheses are the
retained finite metric contract, and it costs <9e18 epsilon0 in the
same units. This use is restricted to that original path.

## 6. Finish the endpoint contract without a sixth derivative

The pure endpoint divided-jet bound is already
1.000001 epsilon0+122d. For the one-center bound through order four,
differentiate the endpoint expression only once in x:

    e_ell,x(H)+e_ell,z(H)H_x,
    v_x(H)+v_z(H)H_x,
    sigma_seed,x(H)-sigma_seed,x+sigma_seed,z(H)H_x,
    partial_x log|H'|.

Divide by the additional 100 in the mixed norm. The first term retains
the pure composition multiplier <1.000001. The e_ell,z term costs
<.001d, since its finite coefficient multiplier is <1e4 and
epsilon0 a_x<=3.31e-9. The two v terms cost <1.001d.
For the seed, Cauchy on the fixed analytic center disk of radius .01
gives |sigma_seed,x|<5. Its full rational shift majorant therefore uses
.05 in place of .027; together with the seed-z term it costs <60d.
For the log term use

    (50t+220t T)/(1-100t-200t T), t=d a_x<=1.655e-12,

after dividing by the extra 100. Its divided coefficients through order
four are <100t. These finite coefficient calculations give the generous
mixed allowance 1.000001 epsilon0+500d. They never bound a moving
C5 difference by a sixth derivative. The entire composed e_ell(H)
is bounded, rather than asserting C5 Lipschitz translation.

Set

    epsilon_end=1.000001 epsilon0+1000d
                              <=1.500001e-6<2e-6.     (13)

The normal potential remains <.026 and its raw second-jet addition is
<.000133. The retained real-C5 flat Born proof therefore applies to the
fixed normalized endpoint chord with (13). Its one/two-vertex remainders
and full tagged Neumann tail cost <1.2e21 in M^-3/M^-2 units. The
original quadratic reserve eight remains valid on that chord: the analytic
seed coefficient is <1.664 and its endpoint increment is <.000133.
The looser coefficient ten used in the row proof does not enlarge the
Born tail. Its denominator remains 1-5(.026)>.87. The
general-Hessian cancellation, odd polynomial cancellations and the same
finite a2 restoration are retained. Local interpolation is used only for
this flat source, not to create a fictitious global physical metric path.

Adding (12), the original path's local term and the proved distant
contribution fits the original per-mode reserves K=K'=1e22. The
complete N(E)<=4CE^2 angular summation and its low/high split then give

    |S_endpoint-S_seed|<=1e23(2 eta_G)epsilon_end a_x,
    |delta S(chi)-delta S(0)|
           <=1e24(2 eta_G)epsilon_end hL_h sqrt(chi). (14)

The local finite restoration fits its retained 1e13 epsilon_end C^2
allowance and is included. The homogeneous value is not assigned zero.
This direct decomposition establishes the needed enlarged endpoint source
contract. It does not claim every proposed matrix-log chord estimate.

## 7. Preserved cold canonical neighborhood and its exact limit

For each fixed ell, take rho=2.5e-10 in D_ell. Every real background
D_ell(v)<=rho and another real parameter radius rho lie in (13).
The independent fixed-chart COMPLEX source theorem applies because
P5_40(ell u+v)<1.3e-4+5e-10<.001<.011. Its two-point bound is 5e-11.
The real amplitudes from (14) are 3.4e-15 a_x and
3.4e-14 hL_h sqrt(chi). The proved half-disk smallness argument, with
logarithmic bound 91 and polarization, now gives actual massive-canonical
response bounds on this neighborhood:

    |D S[w]|<.007 D_ell(w) a_x,
    |D^2 S[w,z]|<2e10 D_ell(w)D_ell(z) a_x,
    bolt-oscillation constants .07 and 2e11
                       in hL_h chi^(1/4).             (15)

The sharper logarithmic tail gives C2 continuity in the quarter-power
weighted target, not just pointwise derivatives. The canonical deviation
Ward constants .024,6.9e10 and bolt constants .24,6.9e11 follow as in the
earlier conditional ledger. Independent Weyl, LLL, radius and higher
spatial source requirements are not supplied by this scalar theorem.

A correction with D_ell(h)<=1e-10 would leave margin 1.5e-10. No such
actual correction bound is proved. The conditional radius geometry is
still needed even to use the deviation lapse mixed allowance 8e14;
its 4e-6 Taylor allowance is not the full Newton constant. The old
r=1e-8 design needs the unknown solver embedding E_c<=.025.

Decision: retain this cold canonical neighborhood and the earlier outer
one, and move to the ACTUAL correction and signed moment in the complete
solver norms. Keep the failed uniform source-times-inverse route parked.
Full nonlinear mapping, projected existence, exact compatibility, physical
boundary/state, semiclassical validity and stability remain open. No S167
acceptance task or full stable-wormhole objective is marked complete.

## 8. Frozen verification

Freeze the note, checker and review before execution. Check the flat
kernel integrals, quadratic insertion and row ledgers, fixed-weight graph
budgets, nine-term mixed product and mass extrema. Independently compare
Bessel quadrature with the positive-integral transformation, finite mixed
normalizer derivatives on real C5 non-C6 fields, and noncommuting inverse
product identities. Verify the full finite mixed-jet majorants, response
constants and unchanged completion gates. Use normalized 60/90-digit
diagnostics, preserve every run and inherited hash, and invoke no Lean.
These checks support the ordinary continuum proof; they are not samples
of the actual quantum state or machine verification of the theorem.

Remaining to-do list: actual correction in the admitted coefficient norms
and signed moment; complete nonlinear map, projected solution and exact
compatibility; physical boundary/state, validity, matching, stability,
charge/return, opening, transport, preparation, separate creation and
full-result audit.
