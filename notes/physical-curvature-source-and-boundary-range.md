# Physical curvature source norm and the charged boundary range condition

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full wormhole goal remains incomplete.

The gauge-invariant curvature row has an exact wave formulation on the
retained charged seed. Its natural scalar Hilbert space, source weight,
diagonal potential and two outer traces can all be specified. If both
fast traces are removed, the horizon-admissible scalar equation has a
closed range of codimension one at every Laplace frequency Re gamma>0.
The missing equation is an explicit signed moment of the TOTAL forcing,
including its metric coupling and initial-data terms.

This supplies a physical norm and an actual seed boundary test beyond
the preceding volume-norm obstruction. It is not a decoupling of the
curvature row, a physical metric boundary choice or a charged instability
theorem. The other metric row and constraints must determine the coupling.

## 1. Gauge-invariant row, including the off-shell contact

Use the full retained trace T=-tr calE and constants from
[the charged conformal identity](charged-conformal-response-and-constraints.md):

    epsilon=G/(60pi), E=60pi/G+11lambda/3, Z=R+12lambda.

In a monotone static exterior, let delta r be the area-radius variation,
xi^r=delta r, h_A=h-Lie_xi g and

    Psi=delta R-R_r delta r,
    f_A=(DT[h]-T_r delta r)/epsilon,
    B[h_A]=K_R[h_A]+7 delta_A(C^2)/12+11 delta_A(J^2)/6. (1)

Here R_r,T_r are derivatives with respect to background areal radius.
The areal gauge may have a time-radial component; it is not set to zero.
Both Psi and f_A are invariant under spherical infinitesimal coordinate
changes. Scalar covariance gives DT[Lie_xi g]=delta r T_r even when the
seed is off shell. At fixed magnetic flux delta_A M=0. Thus the exact
row on the seed, or on an exact background in this model, is

    (Box+E-11Z/36)Psi+B[h_A]=f_A.                       (2)

On an exact solution T=0 the contact in f_A vanishes. It does not vanish
by gauge choice on the seed. B is the complete local metric coupling
already specified in the charged note. It contains curvature variations
and the variation of Box acting on R. No independent state perturbation
changes the bulk anomaly trace, but state response still enters the other
metric equation. Equation (2) is not a closed scalar evolution law.

## 2. Exact wave coordinate and source norm

Write the static metric using inward tortoise coordinate y as

    ds^2=N(y)^2(-dt^2+dy^2)+r(y)^2 dOmega^2,
    0<y<infinity,

with the outer boundary at zero and the nonextremal horizon at infinity.
For spherical fields,

    Box Psi=N^-2[-Psi_tt+Psi_yy+2(r_y/r)Psi_y].

Set u=r Psi. Multiplication of (2) by r N^2 gives EXACTLY

    u_tt+A_R u = g,
    A_R=-partial_y^2+V_R,
    V_R=r_yy/r-N^2(E-11Z/36),
    g=r N^2(B[h_A]-f_A).                               (3)

The scalar kinetic Hilbert norm and the corresponding row-source norm are

    ||Psi||_kin^2=integral r^2 |Psi|^2 dy=||u||_L2(dy)^2,
    ||f||_src^2=integral r^2 N^4 |f|^2 dy
                                      =||r N^2 f||_L2(dy)^2. (4)

Constant angular normalizations are suppressed equally on both sides.
In an areal chart ds^2=-N(r)^2dt^2+A(r)^2dr^2+r^2dOmega^2,
dy/dr=-A/N, so the two weights are r^2 A/N and r^2 A N^3.
This is the scalar Hilbert measure N^-1 dSigma. It is not spacetime
volume N dSigma, and positivity of this Hilbert norm does not imply
positivity of the wave energy involving A_R.
[Ishibashi--Wald, section 1, equations (2)-(4)](https://arxiv.org/pdf/hep-th/0402184v4)
states the established static scalar Hilbert prescription. The metric
coupling and off-shell source in (1)-(4) are retained here explicitly.

These weights specify a necessary curvature-row component of a future
physical norm. They do not by themselves define a complete metric/state
energy, an admissible gauge at the horizon or all constrained initial data.
In particular B[h_A] belonging to the source space is an obligation,
not a consequence of assuming an areal gauge on the open exterior.

## 3. Actual seed potential at both ends

For the actual centered seed N^2=F(r), dr/dy=-F and
y(r)=integral_r^infinity F(s)^-1 ds. Its exact rational function is

    F=lambda r^2+1-2m/r+Qhat/r^2+c3/r^3+c4/r^4,
    c3=-2d B^3, c4=d B^4,
    Qhat=B^2(1+3lambda B^2+d).                         (5)

The thermal mass m includes the prescribed positive mass increment.
Qhat is the seed metric coefficient, not a changed physical magnetic
flux. The preceding centered-seed theorem proves a unique simple
outer root r_h with F_r(r_h)=2kappa>0 and F>0 outside it.

Direct curvature differentiation yields

    Z=-2c3/r^5-6c4/r^6,
    V_R=F[F_r/r-E+11Z/36].                            (6)

At infinity the exact tortoise integral gives

    y=1/(lambda r)-1/(3lambda^2 r^3)+O(r^-4),
    r=1/(lambda y)-y/3+O(y^2),
    F=1/(lambda y^2)+1/3+O(y),
    r_yy/r=2/y^2+2lambda/3+O(y).

It follows that, with a=E/lambda-2 and k^2=a-1/4,

    V_R=-a/y^2+W(y),
    W(y)->(2lambda-E)/3, y->0,
    a=60pi/(G lambda)+5/3,
    k^2=60pi/(G lambda)+17/12.                         (7)

In particular W is bounded on a sufficiently small fixed outer collar.
At the horizon, r-r_h=O(exp(-2kappa y)), Z and F_r/r remain bounded,
and V_R=O(exp(-2kappa y)). Thus the horizon potential is integrable
and the endpoint is limit point. These are statements about each
retained positive-temperature member, not a zero-temperature limit.
No exact semiclassical solution is substituted for the seed.

## 4. Outer traces with a controlled bounded remainder

For (A_R+gamma^2)phi=0 put q_gamma=gamma^2+W and k=sqrt(a-1/4)>0.
Fix a reference length y0 for the logarithmic phase. The normal modes
are phi_+^0=sqrt(y)exp(ik log(y/y0)) and its minus partner. The Volterra
construction is

    phi_+(y)=phi_+^0(y)+integral_0^y
       [sqrt(y s)/k] sin(k log(y/s)) q_gamma(s) phi_+(s) ds. (8)

The same construction with the minus sign gives phi_-. If |q_gamma|<=M,
division by sqrt(y), iteration and the elementary integral of s give

    |phi_+(y)/sqrt(y)| <= exp(M y^2/(2k)),
    |phi_+(y)/sqrt(y)-exp(ik log(y/y0))|
                                      <=exp(M y^2/(2k))-1. (9)

The conormal derivative of the difference is bounded by
M y^2 exp(M y^2/(2k))/2. The bounds also hold for phi_- and uniformly
for gamma in a fixed compact causal-frequency set. Let
phi_c=(phi_++phi_-)/2 and phi_s=(phi_+-phi_-)/(2i).
Their Wronskian phi_c phi_s'-phi_c' phi_s is exactly k.
For real gamma these basis functions are real.

Every maximal-domain function has well-defined leading coefficients
(A_c,A_s) in this basis. For L2 forcing, variation of constants gives
a remainder o(y^(3/2)); its derivative is o(y^(1/2)). The leading
coefficients are independent of the gamma used to define the basis.
Both homogeneous solutions are L2 near zero, and the boundary Green form is

    k(conj(A_c)B_s-conj(A_s)B_c).                      (10)

The one-phase scalar self-adjoint domains are the maximal isotropic lines
cos(theta)A_c+sin(theta)A_s=0. Setting both traces to zero is the minimal
outer domain, not such a line. The established endpoint framework is
[Teschl, section 9.2, Theorem 9.6](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf).
Neither this scalar phase nor the two-zero prescription is selected
as a physical metric boundary condition here.

The actual diagonal A_R is unbounded below for every self-adjoint
extension. The retained compact test eta=(z-1)^2(2-z)^2 on (1,2)
has ||eta'||^2/||eta||^2=12. Its dilation to (d,2d) gives

    Rayleigh quotient <=(12-a/4)/d^2+||W||_infinity.

Since the retained a>96, this tends to minus infinity, independently
of the boundary phase. Smooth compact approximations preserve the
strict negative margin. This conclusion concerns the diagonal scalar
block in (3); the B[h_A] coupling is not set to zero in the physical
metric problem. It is not an actual charged growing-mode theorem.

## 5. Horizon-admissible solution and exact range at every causal frequency

For Re gamma>0 there is a unique horizon solution w_gamma normalized by
w_gamma(y) exp(gamma y)->1. On a sufficiently far tail, writing
m_gamma=exp(gamma y)w_gamma gives

    m_gamma(y)=1+integral_y^infinity
       [1-exp(-2gamma(s-y))]/(2gamma) V_R(s)m_gamma(s) ds. (11)

The integrable horizon potential and Volterra iteration prove existence,
uniqueness and holomorphy in gamma. Extend by the ordinary differential
equation through every compact interval. By (8)-(9) w_gamma is L2 also
at zero, with nonzero trace pair. It is therefore in L2(0,infinity).
Its L2 norm is locally bounded in the open causal half-plane. The
second horizon solution grows as exp(gamma y).

Let D_min consist of U in the maximal L2 domain of A_R with A_c=A_s=0
at zero, and the limit-point domain at the horizon. Define the bounded
complex-linear functional

    M_gamma(f)=integral_0^infinity w_gamma(y) f(y)dy,
    ||M_gamma||=||w_gamma||_L2.                         (12)

The integral is bilinear, not an inner product with w_gamma. Its Riesz
representative is conjugate(w_gamma)=w_conjugate(gamma). Then

    Ran(A_R|D_min+gamma^2)=ker M_gamma,
    ker(A_R|D_min+gamma^2)={0}.                        (13)

In particular the range is closed of codimension one and the operator,
with its graph domain, has Fredholm index -1. This holds for the actual
seed diagonal, including complex gamma with positive real part.

For completeness, choose a second solution phi with W(phi,w_gamma)=1.
The Green kernel G(y,s)=-phi(min(y,s))w_gamma(max(y,s)) has derivative
jump -1 and gives a bounded L2 map. Near zero its kernel is bounded by
C sqrt(y s). On the horizon tail it is bounded by C exp(-c|y-s|), plus
decaying products, for c=Re gamma>0; the mixed regions have the same
integrable bounds. Schur estimates prove boundedness. Its output has
leading trace pair equal to -M_gamma(f) times the nonzero trace pair
of phi. The other variation-of-constants terms are o(y^(3/2)). Thus
M_gamma(f)=0 is sufficient for both zero traces and horizon admission.
Green's identity proves necessity. The homogeneous solution with both
zero traces is zero, proving uniqueness. This also constructs a bounded
inverse on the compatible subspace, without choosing an inverse on
all sources. For real gamma, a compact f=chi w_gamma with chi>=0,
nonzero where w_gamma is nonzero, has M_gamma(f)>0. The functional is
therefore nontrivial. The same follows for complex gamma using
f=chi conjugate(w_gamma).

For compact f, solving with both outer traces zero produces a horizon
growing coefficient exactly -M_gamma(f)/(2gamma), when the growing
solution is normalized as exp(gamma y). This fixes the sign and
normalization of the range condition. Arbitrarily small incompatible
scalar sources exist, but they are not automatically admissible
initial data of the full charged metric problem.

## 6. What the coupled metric problem must satisfy

Assume the time transforms exist in the spaces (4), and let u0,u1 be
the initial curvature data. The Laplace transform of (3) obeys

    (A_R+gamma^2)U=gamma u0+u1+g_hat(gamma).

Consequently the two-zero-trace proposal requires at EVERY such frequency

    M_gamma[gamma u0+u1
                   +r N^2 (B[h_A]-f_A)_hat(gamma)]=0. (14)

The hat means the transform of the complete time-dependent expression.
If B contains time derivatives, their initial contacts must be retained;
one cannot use a zero-initial-data frequency symbol on nonzero initial
jets. The gauge contact T_r delta r is likewise kept off shell.

Equation (14) is a signed, frequency-dependent condition on the coupled
metric problem. The other retarded row and physical initial constraints
may change its value; no cancellation or noncancellation is established
by (13). The metric coupling is unknown, not zero. Independent state
response can influence that other row despite its zero trace. Ward
recovery and neutral current decoupling are available but do not supply
the required metric solution or boundary phase.

The next concrete task is to estimate the mixed curvature/shape moment
in (14) in these physical source norms, with admitted initial data and
metric boundary traces. The previous volume-norm sequence does not
answer this question. No nonlinear exterior, physical mode, cutoff,
order reduction or S077 survival tolerance is obtained here.

## Verification scope

The focused checker derives the scalar wave transform, physical weights,
seed curvature and potential asymptotics exactly. It tests the bounded
remainder Volterra estimate with imaginary-order Bessel solutions, and
the horizon moment sign against independently integrated forced ODEs
on two grids for a solvable short-range model. Compatible compact sources
are separate controls. These fixtures verify the range mechanism, not
the physical benchmark spectrum or a full metric evolution. The ordinary
continuum argument is not a Lean certificate; no Lean run is selected.

Remaining to-do list: coupled physical source moment, metric initial and
boundary data, and constrained stability transfer or justified exclusion;
S167 signed source, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, full-sector stability, transport,
preparation, separate creation and final audit.
