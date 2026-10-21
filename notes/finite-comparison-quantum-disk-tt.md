# Full disk-chart TT positivity at the finite comparison

2026-10-20, E003/S167 supporting S139/S074. Ordinary proof before
frozen verification. **The full retained disk-chart action has a
positive TT-return form. Physical tensor/density contacts, mixed
blocks, coupled domains/range and the nonlinear goal remain open.**

For every original finite endpoint, let h be a genuine spherical TT
input in the smooth bolt-regular, compact-outer core, and (S,T) its
original outer-fixed disk return. Then

    I_full,disk''[(S,T),(S,T)] >= .90 E_ell[c].          (1)

The stronger estimate below retains positive fourth-derivative terms.
Equation (1) uses the complete paired angular tower, original common
scale and coefficients, fixed reference subtraction, Weyl completion
and canonical LLL. It is a disk-chart Hessian, evaluated on a returned
TT input. At a nonstationary finite endpoint it is not automatically
the physical exponential-tensor Hessian or fixed-density gradient's
TT block. Their nonzero-gradient contacts are still required.

## 1. Retained geometry, return and allocation

Use lambda=lambda_e, x=lambda B^2, ell_o=B/sqrt(C),
gbar=ghat/ell_o^2, t=r_ell/B and omega=log t. The
[massive comparison](finite-comparison-massive-geodesic-form.md) gives

    Gamma_hat,mass,disk'' >= Q_geo-.001 C J2[T]-.005 E_ell,
    Q_geo=sum d_n P_Mn[T]+C kappa4 J2[T]/(30pi)
                  -C kappa4 integral Rbar |grad T|^2/(90pi),
    13.9<kappa4<17.5, |Rbar|<2.001.                    (2)

Jj are the actual physical-half orbit norms defined there. Its global
correction and p0 contact are already included; do not add either
again. Keep the positive spectral tower in Q_geo for the outer end.

The [classical disk theorem](finite-comparison-classical-disk-contact.md)
gives (I0+Gamma_hat,LLL),disk''>=.92 E_ell, where
I0=S_classical+Gamma_ref. The remaining allocation is therefore
Gamma_hat,mass+W4-Gamma_ref. This is the original fixed potential,
not a recentered reference functional.

Put e=8pi G E_ell, Cc=||c|| in actual physical four-volume.
The same theorem and differentiated disk return give

    ||S||<12 Cc, ||T||<18 Cc, ||S+T||<30 Cc,
    ||S'||<116 sqrt(e), ||T'||<210 sqrt(e),
    |Rphysical|<2.25 sqrt(lambda), |E2|<22.001lambda.

For example S=c/2-Rphysical eta, ||eta||<4.83 Cc/sqrt(lambda),
||b||<2.002 Cc, ||A eta||<1.008||b+c|| and
S'=c'/2-Rphysical(b+c)-E2 eta, ||c'||<5.32sqrt(e).
No separate singular A bound is used at the bolt.

In physical canonical orbit units the finite geometry has

    |Rhat|<13lambda, |Delta_hat omega|<23lambda t^2,
    |grad_hat omega|^2<5.1lambda t^2, B^-2<=lambda,
    0<=omega<=t^2, dmu_hat=dmu_4/(4pi B^2 t^4).       (3)

The first bound follows from (2) and C/x<6.1. The next two are
the exact identities Delta_hat omega=t^2 E2 and
|grad_hat omega|^2=t^2 Rphysical^2. Thus these are actual finite
coefficient bounds, not substituted seed equations.

## 2. Exact Weyl Hessian and finite-background bound

Use the compact bulk polynomial and retained boundary representative
of the [seed curved proof](stationary-curved-quantum-tt.md), with
actual endpoint coefficients throughout. Write R=Rhat,
U=Delta_hat omega, V=|grad_hat omega|^2, y=B^-2,
k=y^2(1/15-q^2/3), d=Delta_hat T, f=Delta_hat(S+T/3),
w=grad_hat omega.grad_hat S and j=|grad_hat S|^2. Direct variation
of the fixed anomaly polynomial gives exactly

    W4''= B^2/(20pi) integral f^2
        -B^2/(30pi) integral(omega+1/6)d^2
        -B^2/(2pi) integral[d D+f L+J j+W T w-7w^2/45+Z],
    D=S(-R/15+y/9)+T[(6omega+1)R/45+11V/135]-11w/135,
    L=T[-R/15+2U/5+7V/45]-7w/45,
    J=R/30-7(U+V)/90-y/18,
    W=-2R/15+14(U+V)/45,
    Z=T^2[omega R^2/30+R(U+V)/15-U^2/5
                         -7UV/45-7V^2/90+2k omega]
                  +ST[-R^2/30+2k].                  (4)

The checker differentiates the original density polynomial, rather
than merely rechecking this displayed completion. U cancels from D;
the ordinary 1/6 term remains.

Here is a ledger in the single measure dmu_4. Divide each integrand
by t^4, and divide its integral by lambda e:

| Integrand | Upper bound |
| --- | ---: |
| D^2 | 1e4 |
| L^2 | 6e4 |
| abs(J) j | 4e4 |
| abs(W T w) | 6e4 |
| abs(Z without k) | 6e4 |
| w^2, for absolute continuity | 7e4 |

Indeed ||D/t^2||<lambda||S||+2.5lambda||T||
+(11/60)sqrt(lambda)||S'||; and
||L/t^2||<11lambda||T||+.35sqrt(lambda)||S'||.
Also |J|/t^2<2.7lambda, |W|/t^2<11lambda, and the non-k
T^2 and ST coefficients of Z/t^4 are below 160lambda^2 and
6lambda^2. These bounds include all radius and measure weights.

Use |dD|<=d^2/100+25D^2 and
f^2/10-|fL|>=f^2/20-5L^2. The w^2 term is favorable.
The combined non-k ledger is below 710001lambda e; multiplying
by 1/(8pi^2) gives a cost below 1e6 G lambda E_ell.
Since |k|<q^2/(3B^4), the k terms cost at most
(360/pi)Gq^2/(lambda B^4)E_ell<1e-24 E_ell. Consequently

    W4'' >= -C/(30pi) integral(omega+19/60)(Delta_bar T)^2
        +C/(40pi) J2[S+T/3]
        -[1e6 G lambda+1e-24]E_ell.                 (5)

The retained finite-boundary representative adds no Hessian on this
compact-outer core. No claim about general outer traces follows;
their stronger approximation and boundary flux remain separate.

## 3. Actual outer collar and the complete scalar model

Let Y be actual normalized proper distance from the outer seam,
gbar=dY^2+a(Y)^2dtau^2, h=log(a/a(0)). The even seam has h'(0)=0.
The retained proper chord, after its cubic, satisfies

    |partial_Y^j delta h| <=1e-7 j!100^j Y^(4-j)L(Y)
                 +|partial_Y^j(c3 Y^3)|, j<=5,
    L(Y)=1+|log Y|, |c3|<1e-11.                      (6)

Here delta h is the actual endpoint minus the fixed seed in each
endpoint's proper coordinate. Thus on Y<=.1 the new value, first,
second and third jets are below 4e-11, 4e-8, 7e-5 and .2.

Use the seed's complete even extension h0,ext. Extend delta h from
.1 to .2 by its quadratic Taylor polynomial at .1 times
1-10s^3+15s^4-6s^5, s=10(Y-.1), and zero beyond .2; reflect evenly.
This matches its first two jets at .1 and has a Lipschitz second
derivative. Its first three derivative increments are below
2e-5, .001 and .2. Adding the retained seed extension bounds gives

    |h_ext'|<1.201, |h_ext''|<42.001, |h_ext'''|<3271,
    H=-partial_Y^2+V_H, V_H=h_ext''/2+h_ext'^2/4,
    |V_H|<22, Lip(V_H)<2000,
    H=(partial_Y-h_ext'/2)^*(partial_Y-h_ext'/2)>=0.    (7)

This is a complete even line, agreeing with the actual metric on
the later outer support. No artificial Dirichlet boundary is imposed.
Set v=sqrt(a_ext)T and w_H=Hv. The unchanged all-frequency flat
inequality, with the same angular masses and finite part, is

    Kflat(p)>=C p^4/(30pi)[log(p/(mu0 ell_o))-38/15].  (8)

The seed proof derives (8) at every frequency by differentiating its
convergent remainder and using the digamma integral; it is not an
asymptotic replacement. From H+24>=-partial_Y^2 and operator
monotonicity of log, positive regularization gives

    <w_H,log sqrt(H)w_H> >= <w_H,log|D|w_H>-12<v,Hv>.

The inequality log(1+24/H)<=24/H is legitimate on w_H=Hv and
needs no gap. [Beckner, Corollary 2, equation (10)](https://arxiv.org/pdf/math/0612486)
at n=1,p=2 gives the logarithmic uncertainty constant
D_H=log2+psi(1/4)=-gamma_E-2log2-pi/2>-4. Therefore

    <w_H,log sqrt(H)w_H> >= integral -log|Y| |w_H|^2
                              -4||w_H||^2-12<v,Hv>.  (9)

The full-line inequality applies to the even extension with the
physical-half normalization. Approximation starts on the compact
core; V_H is bounded Lipschitz and w_H has sufficient regularity.

The scalar comparison retains its exact positive heat weight

    beta(t)=sum d M^2/(2pi t) integral_0^1 exp[-tM^2/(alpha(1-alpha))]d alpha,
    beta(t)<=.8C pi^-1 t^-3 exp(-648t).

[Simon's Feynman-Kac formula, theorem 1.1](https://arxiv.org/html/math-ph/9907022),
with variance two for H=-partial^2+V_H, bounds the heat-kernel
difference by (exp(24t)-1)(4pi t)^-1/2 exp(-s^2/(4t)).
Its integrated bound remains 13C/|s|^3. Actual midpoint distance
and half-density comparisons retain |b'|<3, |b''|<6 for b=a^2
normalized at the midpoint, |d-r|<=17r^3 and density error <.31s^2.
Together with the finite normal-cap bound in (2), the radial
comparison fits 2e9 C/|s|^3 locally.

Fix its second contact with the same compact quadratic. For the scalar
test q(Y)=Y^2 at its center, Hq=-2, H^2q=-4V_H(0), and

    |(exp(-tH)-exp(-tD^2))q(0)+2V_H(0)t^2|
         <=4 Lip(V_H)t^(5/2)+484t^3 exp(22t).

The complete massive integral costs <100/M^3+3/M^4; the tower
is <70C. Compactification costs <1e4C. The real normal and
half-density conversion, including their cubic/quartic jets, fits
the finite second-contact reserve 3e13 C in (2); take 1e14 C
for this model comparison including the known curvature-gradient
contact. Only second and zero contacts remain after the common
finite-mode fourth cancellation. Both complete models annihilate
constant base T: H sqrt(a_ext)=0 and Kflat(0)=0. This fixes the
half-density zero contact without deleting restricted-collar escape.
Retain both completions' bounded tails: the auxiliary tail at
R=1/400 costs <2e11 C J0 and the geodesic cos-1 tail <2e12 C J0.
Their derivative polynomials are already in the contact calculation.

## 4. Localization with the finite fifth jets

Use a_cut=1/100 and the same C5 partition as the seed proof:
theta=(pi/2)P((Y-a_cut)/a_cut),
P(s)=462s^6-1980s^7+3465s^8-3080s^9+1386s^10-252s^11,
extended constantly, chi_o=cos theta, chi_i=sin theta.
The derivative reserves for chi, orders one through six, are
a_cut^-j(5,375,12000,650000,3.8e7,3e9), the sixth a.e.

On the transition's R=1/400 neighborhood, .0075<=Y<=.0225,
(6) has L<6. Combining it with the seed gives actual h derivative
bounds through five

    (.05,3,11,1500,1e8).

The actual geodesic equation Y''=h'(Y)(1-Y'^2) then gives the
successive reserves

    (|Y'|,...,|Y6|)<=(1,1,10,100,2000,2e8).

In particular the larger fifth metric derivative is retained. The
sixth chain-rule coefficients are (2e8,14500,2115,245,15,1).
They still yield |X6 chi|<1e22 and |X4 chi|<1e14. For the same
anchored q_x=(Y-Y_x)^2 xi, radius R, Leibniz and Bell composition
give |X6(chi q_x)|<1e20 and |X4(chi q_x)|<1e12.
These estimates require no sixth metric derivative and are used only
away from the seam. When the probe misses the transition the IMS
operator on that probe is zero.

The exact spectral point-action coefficient is <2e-5 C sup|X6 u|,
and the common fourth action is <=C sup|X4 u|. Hence both finite
IMS contacts stay below 1e22 C. The gain-two separated IMS kernel,
the scalar comparison and its anchor all fit K=1e14 in the
unitary-transport split of (2). The full-flow cosine IMS tail costs
<2e13 C J0, including returns; polynomial contacts are counted in
the point actions. No global pointwise geodesic-kernel assertion is made.

For the Weyl weight w0=omega+19/60 the exact integrated fourth IMS
difference, before multiplying C/(30pi), is

    6w0 theta'^2 T'^2
    +[8w0 theta'theta''+4w0 h'theta'^2+2w0' theta'^2]TT'
    +w0[theta'^4+(theta''+h'theta')^2]T^2.

All coefficients are supported in the transition, where
w0<7, |w0'|<2/a_cut, |h'|<5a_cut. Its absolute cost is bounded by
6e4 a_cut^-2 J1+1e6 a_cut^-4 J0,transition. The same identity with
w0=1 converts the global unweighted J2 into the localized sum
without a factor-two loss; its lower terms are included below.

## 5. Cancel the outer logarithm and close the disk form

The original fixed-scale benchmark gives
-log(mu0 ell_o)+.5log(x/C)>9. Actual proper/radius changes are
below 1.2e-8 and 3e-8, so on the outer support

    omega<=-log Y-.5log(x/C)+.5log2+1e-6.

The seed polynomial bound used here is valid on the larger .1
collar, covering the small proper-coordinate shift. Combining
(5),(8),(9) leaves outer coefficient, in C/(30pi) units,

    9-38/15-19/60-4-(log2)/2-1e-6 >8/5.             (10)

The Hardy gradient cost is 12C J1/(30pi). On the interior support
monotone radius and the same comparison give omega<6.001; hence
kappa4-omega-19/60>7.5. The curvature-gradient part of Q_geo
is bounded separately. Thus the unbounded position logarithm has
been cancelled, not capped by a constant on the entire geometry.

Apply the gain-two split with epsilon=.001, K=1e14, R=1/400.
Its gradient coefficient is <1e16 C. The massive comparison in
(2) costs .001 C J2 globally; model IMS costs at most
epsilon*pi*C J2/4 plus a curvature-gradient term; the scalar
comparison costs at most 2epsilon*C J2[chi_oT]/3 plus
968epsilon*C J0[chi_oT]/3. The latter uses
||v''||^2<=2||Hv||^2+968||v||^2, from |V_H|<22.
After unweighted IMS conversion the total fourth loss is below
.00246 C times the localized J2 sum. The margins in (10) and
the interior therefore retain at least

    C/(100pi) sum_(j=o,i) J2[chi_jT]
                +C/(40pi) J2[S+T/3].               (11)

For the lower terms use the literal aggregate J1* of T, both
localized fields and auxiliary v, and the corresponding J0*.
The physical weights from (2) give C J0/E_ell<5.2e-29 and
C J1/E_ell<8.5e-28. On the cutoff support Y<=2a_cut,
Y t>2.23 gives C J0,outer/E_ell<3.4e-37. Thus

    C sum J1[chi_jT]/E_ell<8.6e-28,
    C J1*/E_ell<1e-26, C J0*/E_ell<2e-28.           (12)

The auxiliary half-density derivative costs at most twice the
localized Dirichlet norm plus .001 J0,outer. All new finite
second contacts, IMS, curvature, Hardy and split gradient terms
fit 1e23 C J1*. All new bounded zero terms and model/IMS escape
terms fit 1e23 C J0*. The actual massive gluing is already in
the .005 E_ell term of (2). The new costs are consequently below
.001 E_ell and .00002 E_ell. The weighted Weyl IMS cost is below
1e-12 E_ell by its explicit coefficients and (12).

Finally, with the same fixed reference profiles evaluated at the
actual r>=B, its exact Euclidean disk Hessian is

    Gamma_ref''=-4 integral_4[(S+T)^2 pA+2S(S+T)pS
                                +S^2pS+(r/2)S^2 partial_r pS].

The second variation of X=e^(2S)r^2 supplies S^2pS. From the
explicit inverse-power profiles, |pA|,|pS|<3341q^2/B^4 and
|r partial_r pS|<1.4e4q^2/B^4 for every r>=B. Their coefficient
sum using (12Cc,18Cc,30Cc) is below 7e6. Therefore
|Gamma_ref''|<1e-16 E_ell. This subtracts the reference added in
I0 and introduces no new action coefficient.

Equations (2),(5),(10)-(12) prove

    (Gamma_hat,mass+W4-Gamma_ref),disk''
       >= -.01 E_ell + the positive terms (11).

Combining with the .92 classical/LLL margin gives .91 E_ell
plus (11), proving the conservative claim (1).

## 6. Completion and unresolved physical contacts

The estimate extends to the explicit completion of the compact-outer
regular TT core in

    N_ell^2=E_ell+sum d_n P_Mn[T]
          +C[J2[T]+J2[S+T/3]+integral omega(Delta_bar T)^2].

Positive spectral/weak derivative components are closed after the
fixed actual return, whose L2 norm is controlled by E_ell. The
displayed estimates make the complete disk Hessian continuous in
this norm; its lower bound passes to the completion. This does not
identify N_ell with the full physical logarithmic operator domain,
prove density for all prescribed mass/fast directions, or close the
candidate range/adjoint theorem. In particular outer boundary traces
and approximation must be proved in the norm retaining the logarithm.

At this nonstationary endpoint the full gradient's tensor-chart,
gauge and fixed-density contacts must still be bounded before
identifying (1) with the block used by a coupled physical inverse.
The earlier classical contact theorem controls I0's contact only.
Mixed scalar/TT blocks and nonlinear regularity are not consequences
of disk positivity. Pi_cut stays fixed at the seed.

The frozen checker pins (2), the classical/return geometry and seed
spectral comparison. It checks (4) directly, actual low Weyl budgets,
the auxiliary extension and Hardy constant, finite fifth-jet IMS
chain rules, model/contact budgets, exact weighted IMS identities,
physical weights, fixed reference variation and the final margin.
Continuum heat, flow, common finite-part and domain arguments remain
ordinary proofs. This is not formal verification or a numerical
semiclassical solution. Every scientific run is retained.

Remaining to-do list: physical tensor/gauge/density contacts and mixed
blocks; full logarithmic domains and candidate coupled range; fixed
projection and cold/higher response; preserved nonlinear family and
projected solve; signed compatibility; exterior/state, charge/return,
opening, stability, transport, preparation, creation and full-result audit.
