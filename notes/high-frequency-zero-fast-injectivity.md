# High-frequency injectivity on the auxiliary zero-fast domain

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

At every sufficiently large admitted positive thermal frequency, the full
coupled auxiliary operator with BOTH fast curvature traces removed is
injective. Its range has codimension exactly one, and the inverse on that
range is O(gamma^-2). Its single adjoint compatibility vector has a
convergent bounded-perturbation formula. This improves the earlier index
statement, which allowed a nonzero kernel and more than one condition.

The threshold is finite but not numerically calibrated for the retained
benchmark. The domain is still an auxiliary choice. No physical boundary
selection, continuous causal inverse, nonlinear background or stability
theorem follows. In particular a codimension-one range is not an inverse
on freely prescribed sources.

## 1. A uniform bound for the actual coupled remainder

Keep the original action, HH state, paired domain, clock, flux, angular
tower and all finite/residual contacts. Let gamma be a positive Matsubara
frequency at least Gamma_*=max(Gamma_T,Gamma_R)>0, from the retained
[shape inverse](global-thermal-shape-inverse.md) and
[physical reconstruction](physical-curvature-reconstruction.md).
The [strong coupled operator](coupled-thermal-boundary-index.md) is

    L_gamma=Atilde+gamma^2+B_gamma,
    B_gamma=C_gamma-epsilon^-1 P R_gamma
                         +epsilon^-1 D_gamma K_cross^dagger,
    C_gamma=A1/6+(B1 partial_y+B0)R_gamma,
    epsilon=G/(60pi), P=F(r calT_r+4calT).                 (1)

All coefficients here are the retained fixed seed coefficients. In
particular P and the density contacts inside D are not set to zero.
The strong feedback is D K_cross^dagger; it need not be self-adjoint.

Set C_R=1/(18lambda), rho=ell y/(ell+y), ell=lambda^-1/2,
M0=sup rho^2 F and

    C_d=[M0/(108lambda)+1/(324lambda^2)]^(1/2),
    C_C=||A1||_inf/6+||B1/rho||_inf C_d+||B0||_inf C_R.

Then ||R_gamma||<=C_R, ||rho partial_y R_gamma||<=C_d and
||C_gamma||<=C_C uniformly along this ray. The finite constants K2,K1,K0
of the [mixed-form theorem](thermal-curvature-shape-schur.md) give

    D_*=epsilon K2 sqrt(200pi)+epsilon K1 sqrt(33G)
                                      +K0 sqrt(40G)/Gamma_*,
    ||D_gamma||_(E_gamma->H)<=D_*,
    ||K_cross^dagger||_(H->E_gamma)<=K_*=C_R D_*/(2G).

Consequently

    ||B_gamma||<=B_*=C_C+epsilon^-1 ||P||_inf C_R
                                      +epsilon^-1 D_* K_*. (2)

This estimate compares norms separately in each E_gamma; it does not
identify different frequency energy spaces or assume uniform compactness.
It follows from the existing energy estimates, not from ignoring the
full nonlocal Dirac response. That response is retained in E_gamma.

Put a=E/lambda-2>1/4 and k=sqrt(a-1/4). The exact seed potential gives

    Atilde=T_min+W, T_min=-partial_y^2-a/y^2,
    W=V_R-F c_R+a/y^2.                                   (3)

W is globally bounded: the scalar expansion cancels the y^-2 term at
zero; F c_R is bounded and tends to zero at both ends; and V_R and a/y^2
tend to zero at the horizon. The middle interval is smooth. Bounded
multiplication leaves the maximal graph domain and both outer traces
unchanged, so the two operators have exactly the same D_min. Define

    V_gamma=W+B_gamma, M_*=||W||_inf+B_*,
    L_gamma=T_min+gamma^2+V_gamma, ||V_gamma||<=M_*.         (4)

M_* is finite and independent of gamma on the admitted thermal ray.
The numerical benchmark value is not evaluated here.

## 2. A scaled inverse-square inverse on its compatible range

Use H=L2(0,infinity) and the minimal outer domain A_c=A_s=0, with the
limit-point horizon. At unit frequency define

    h(y)=sqrt(y) I_(ik)(y), b(y)=sqrt(y) K_(ik)(y),
    G_a(y,s)=h(min(y,s)) b(max(y,s)).                      (5)

Both solve (T+1)u=0, and W(h,b)=-1. The derivative jump of G is
therefore -1, as required for -partial_y^2. These facts follow from
[DLMF 10.28.2](https://dlmf.nist.gov/10.28).
The [imaginary-order small-argument formulas](https://dlmf.nist.gov/10.45)
give h,b=O(sqrt(y)) at zero for this fixed k>0. The
[large-argument formulas](https://dlmf.nist.gov/10.40) give h=O(exp y)
and b=O(exp(-y)). In particular b is real, square integrable and nonzero.

For an explicit finite bound define

    H_a=max{sup_(0<y<=1)|h(y)|/sqrt(y),
                         sup_(y>=1)exp(-y)|h(y)|},
    K_a=max{sup_(0<y<=1)|b(y)|/sqrt(y),
                         sup_(y>=1)exp(y)|b(y)|},
    C_a=2 H_a K_a.                                       (6)

These constants are finite for each fixed a. The kernel is bounded by
H_a K_a times sqrt(ys) if y,s<=1, exp(-|y-s|) if y,s>=1,
and sqrt(min(y,s)) exp(-max(y,s)) otherwise. Its absolute row integral
is at most C_a; symmetry gives the same column bound. Indeed for y<=1
the bound is H_a K_a sqrt(y)(2/3+exp(-1))<2H_a K_a, and for y>=1
it is H_a K_a[2+(2/3-e)exp(-y)]<2H_a K_a. Schur's estimate gives
||G_a||<=C_a. This is a deliberately coarse bound, not a numerical
large-order estimate for the enormous benchmark k.

For f orthogonal to b, the Green output has both outer traces zero.
Near zero its expression is

    b(y) integral_0^y h(s)f(s)ds
                     -h(y) integral_0^y b(s)f(s)ds.

Cauchy--Schwarz makes this o(y^(3/2)) with derivative o(y^(1/2)).
It is in H, solves (T+1)u=f and belongs to D_min. Conversely Green's
identity imposes <b,f>=0, and the only horizon-admissible homogeneous
solution is b, whose outer trace pair is nonzero. Thus T_min+1 is
injective with range b^perp and compatible inverse norm at most C_a.

Unitary dilation S_gamma u(y)=sqrt(gamma)u(gamma y) preserves D_min:
it rotates and rescales the two fast traces and preserves their joint
zero value. It obeys

    (T_min+gamma^2) S_gamma=gamma^2 S_gamma(T_min+1).       (7)

Let b_gamma=S_gamma b/||b|| and P_gamma=I-|b_gamma><b_gamma|.
The scaled Green kernel is gamma^-1 G_a(gamma y,gamma s). Define its
composition with P_gamma as T_gamma. Then

    ||T_gamma||<=C_a/gamma^2,
    T_gamma A_gamma=I on D_min,
    A_gamma T_gamma=P_gamma,
    A_gamma=T_min+gamma^2.                                (8)

In particular ||u||<=C_a ||A_gamma u||/gamma^2 on D_min. The symbol
T_gamma in (8) denotes this compatible scalar inverse, not the prior
shape reconstruction map Q_gamma.

## 3. Injectivity, a single adjoint vector and exact range

Suppose

    eta_gamma=C_a M_*/gamma^2<1.                          (9)

For example gamma^2>=2 C_a M_* gives eta<=1/2 when M_*>0.
If M_*=0 the conclusion holds at every admitted gamma. Combining
(4) and (8) proves

    ||u||<= C_a/(gamma^2-C_a M_*) ||L_gamma u||.           (10)

Hence L_gamma is injective and has closed range. Closedness follows
from the unchanged domain and bounded perturbation; (10) then turns a
Cauchy sequence of outputs into a Cauchy sequence of inputs. The earlier
index -1 theorem already implies exactly one compatibility condition.
The following construction proves that fact directly and identifies it.

The identities (8) imply A_gamma^* T_gamma^*=I and
T_gamma^* H subset D(A_gamma^*) intersect b_gamma^perp. The first follows
by testing against D_min; the orthogonality follows from T_gamma b_gamma=0.
Define the convergent Neumann formula

    z_gamma=(I+T_gamma^* V_gamma^*)^-1 b_gamma.            (11)

Its defining equation puts z in D(A_gamma^*) and gives L_gamma^* z=0.
Also <b_gamma,z_gamma>=1. Any other adjoint null vector satisfies the
same equation multiplied by its b_gamma coefficient, because
T_gamma^* A_gamma^* z=P_gamma z. Thus

    ker L_gamma^*=span{z_gamma},
    ||z_gamma||<=1/(1-eta_gamma),
    ||z_gamma-b_gamma||<=eta_gamma/(1-eta_gamma),
    Ran L_gamma={f:<z_gamma,f>=0}.                        (12)

The full coupled auxiliary pair therefore has a UNIQUE solution for
each compatible source and no solution for an incompatible one. For
its strong source F=f+epsilon^-1 D_gamma w_ell,

    ||U||<=C_a ||F||/(gamma^2-C_a M_*),
    w=w_ell-K_cross^dagger U.                            (13)

This is a norm inverse estimate on the compatible range, not positivity
of the curvature quadratic form or a right inverse on all sources.

For completeness the bounded formula

    J_gamma=(I+T_gamma V_gamma)^-1 T_gamma

satisfies J_gamma L_gamma=I and

    L_gamma J_gamma f=f-b_gamma<z_gamma,f>.               (14)

Indeed its output is in D_min since u=T_gamma(f-V_gamma u), and the
residual lies along b_gamma by (8); pairing with z fixes its coefficient.
Equation (14) exposes the discarded source component. It does not
authorize replacing the physical source by a projection.

## 4. Consequences for the original source test and physical scope

The [original trace/orbit source test](ward-sources-and-fast-compatibility.md)
now has exactly one condition above (9). With alpha=epsilon/(2G), set

    nu_gamma=epsilon^-1 D_gamma^dagger z_gamma,
    mu_gamma=z_gamma+alpha Q_gamma nu_gamma.

The condition is <mu_gamma,f>+<j nu_gamma,g>=0, with the shape inclusion
j:E_gamma->H. Its existing coercivity and ||Q_gamma||<4/5 imply

    ||j nu_gamma|| <= sqrt(40G) D_*/
                     [epsilon gamma(1-eta_gamma)],
    ||mu_gamma-b_gamma|| <= eta_gamma/(1-eta_gamma)
                                  +(4alpha/5)||j nu_gamma||. (15)

The adjoint vector is close to the normalized inverse-square horizon
solution at O(gamma^-2); the original shape-source pairing is bounded
at O(gamma^-1). These estimates do not make either overlap exactly zero.
Initial metric source polynomials can grow with gamma, so their physical
compatibility does not follow by multiplying these bounds informally.
Whether nu_gamma vanishes identically is not decided.

The [fixed-HH exponential preflight](causal-preparation-source-map.md)
has no independent bounded growing covariance source. On the auxiliary
two-zero domain considered here, it now has no nonzero homogeneous pair
at all sufficiently high admitted thermal frequencies. Finite frequencies
below the threshold, other boundary domains and continuous causal
frequencies remain untreated. A one-phase gravitational domain is not
D_min, and cannot inherit (10) merely because its scalar endpoint is
reflecting. No physical fast mode is excluded by this auxiliary theorem.

Neither C_a nor M_* has been calibrated on the benchmark. Their present
bounds can be very coarse, so no frequency within a physically controlled
semiclassical regime is certified. The centered seed remains off shell.
The original conditional fast-family obstruction and the static signed
condition remain in scope. No new boundary prescription, finite action
term, source projection or order reduction is adopted.

Focused checks use Bessel Green reconstruction and moment identities at
two precisions, exact dilation algebra, and finite-rank perturbations of
an infinite unilateral shift to test the adjoint Neumann formula and
unprojected residual. The shift tail remains infinite and unchanged;
the fixture is not an approximation to the charged Dirac kernel. The
ordinary operator proof is not inferred from a finite matrix rank test.
No Lean campaign is selected.

Remaining to-do list: physical gravitational boundary domain and
constrained initial/state source overlaps; finite/continuous causal
spectrum and quantitative thresholds; static signed condition, nonlinear
existence and all-sector stability; exterior/state, charge/return,
opening, transport, preparation, separate creation and full audit.
