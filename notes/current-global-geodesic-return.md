# Global geodesic return: actual orbit decomposition and the mass seam

2026-10-24, E003/S169. Ordinary argument before supporting verification.
The full stable-wormhole objective is active and incomplete; the closed
feasibility audit and original physical targets are preserved.

OUTCOME: UNRESOLVED COMPLETE COUPLED ESTIMATE. This turn makes the global
geodesic comparison explicitly evaluable on radial fields, retaining the
actual finite geometry and its trapped trajectories. It also rejects a
smooth-periodic H4 tail estimator on the actual mass lift and derives an
analytic replacement for its entire cubic-jump tail. Neither the global
comparison matrix nor the physical coupled inverse has been evaluated.
No mass-column, admission or nonlinear requirement is marked complete.

## 1. Complete problem and the first unsupported step

The first missing quantitative implication remains the actual logarithmic
Schur inverse and reconstruction in the original strong graph:

    Hbar=[[-Abar,-Cbar],[-Cbar*,Dbar]],
    K=Cbar^-1 Abar(Cbar*)^-1, S=I+K Dbar,
    S T=-Cbar^-1 f+K g, nu=(Cbar*)^-1(Dbar T-g).

Dbar is the complete common-restored shape operator. Both original rows,
mass/slow moments, scalar finite part and regular-bolt matching remain.
The fixed-seed comparison H0 in current-geodesic-two-sided-remainder.md
retains the exact local rows, LLL, local contacts and the complete geodesic
massive form. Its actual H-H0 form error is <.0022; the contact-preserving
multiplier implementation also fits that guard. This is only V->V*, not
the required strong graph estimate. No physical thermal return is set to
zero. The finite-temperature/state-dependent difference is in H-H0 with
its existing complete form bound; its strong return remains uncomputed.

The present attempt develops the ACTUAL global flow in that H0, rather
than replacing it by a scalar Laplacian or a flat finite interval. The
obstacle encountered before numerical matrix assembly is concrete:
mass/log fields are not smooth periodic functions along seam-crossing
geodesics. A tail bound assuming a fourth L2 flow derivative fails even
for the exact cut mass lift. Sections 2-5 give the full orbit decomposition
and address this obstacle without modifying that lift or its physical norm.

## 2. All orbits of the finite seed, including the narrow trapped family

The canonical orbit metric is dz^2/Q+Q d tau^2, tau of ORIGINAL period
calB=2pi/gamma, 0<=z<=z_h. Its even double is C2,1 at z=0 and regular at
both bolts; it is the existing comparison geometry, not a new physical
Dirac realization or a smooth replacement. Write

    Q=beta+u z^2-a z^3+b4 z^4-2eps z^5+eps z^6,
    a=4beta+2u+tau_seed, b4=3beta+u+eps,
    6beta+u+eps=1, Q(z_h)=0, Q'(z_h)=-2gamma.

For a unit-speed geodesic let j be its conserved angular momentum.
Then dot z^2=Q-j^2, ddot z=Q'/2. For 0<z<z_h,

    Q'/[2z(1-z)]=u-(6beta+2u)z+eps z^2(2-3z)
                                      -3tau_seed z/[2(1-z)].

The derivative of this bracket is less than -.99, using eps<0,
|eps|<1e-28, tau_seed>0 and 6beta+2u=1+u-eps. Thus there is exactly one
maximum z_c in each physical half. Q first increases ABOVE beta and then
decreases to zero. The implementation retains this small increase exactly:

    Q(u x)-beta=u^3 F(x),
    F=x^2-a x^3+b4 u x^4-2eps u^2 x^5+eps u^3 x^6.

It never obtains this quantity by subtracting independently rounded values
of Q and beta. The parameter enclosures give

    u=2.03584378760339...e-11,
    z_c=2.03584378756194...e-11,
    F(z_c/u)=.333333333316368...,
    z_0=3.05376568135069...e-11, Q(z_0)=beta.

For 0<=j^2<beta the trajectory crosses the paired seam. Let z_+(j) be its
outer turning point. A radial EVEN base profile has period

    T(j)=2 int_0^z_+ dz/sqrt(Q-j^2).                       (1)

Its full signed radial orbit has twice this period; the base profile is
identical on the paired halves. For beta<j^2<Q(z_c), a trajectory is
trapped in one half between two turning points z_- and z_+, and

    T(j)=2 int_z_-^z_+ dz/sqrt(Q-j^2).                     (2)

Both copies and both angular directions are counted below. The separatrix
j^2=beta approaches the seam in infinite flow time; the limiting circular
orbits and this separatrix have zero Liouville measure. They are treated
by the endpoint limits of the integral, not assigned an artificial wall.

Validated quadrature, with the exact implicit turning-point identities
factored BEFORE square roots, gives these finite-geometry intervals:

    T(0) in (156.704452148,156.704452150),
    T(j^2=beta/2) in (6.43647293852,6.43647293856),
    T(j^2=beta+u^3/6) in (1529608.55545,1529608.55547).

T(0) uses Q(z_h)=0 and Q'(z_h)=-2gamma explicitly and a hyperbolic change
of variable at the regular bolt. The other two values verify the two
branches of the implementation; sampled periods are not a global estimate.

The ENTIRE trapped Liouville measure, divided by calB, is

    M_trap=4 int_0^z_0 atan sqrt[(Q-beta)/beta] dz
          in (1.0993804115e-26,1.0993804117e-26).           (3)

For an analytic enclosure put xcap=1.6 and
Delta=xcap |b4|u+2|eps|u^2 xcap^2+|eps|u^3 xcap^3. The positive part of
F is squeezed between x^2[1-(a+Delta)x]_+ and
x^2[1-(a-Delta)x]_+. Their square-root integrals are
4/[15(a+/-Delta)^2]. Use atan t<=t and atan t>=t(1-tmax^2/3), with
tmax^2<=u^3*4/[27(a-Delta)^2 beta]. This covers the separatrix and all
trapped trajectories, not sampled momenta. Small measure is NOT a small
operator-norm bound on arbitrary localized inputs. Nothing is discarded.

## 3. Exact global compression of the complete multiplier

Let f_j(t)=T_field(z_j(t)) be a physical radial shape profile along the
orbit, f_(j,k)=T(j)^-1 int_0^T(j) f_j(t) exp(-2pi i k t/T(j))dt.
Let A(w) be the COMPLETE flow multiplier of
current-geodesic-flow-multiplier.md, including its restored fourth contact:

    A(w)=C w^4 chi(w)/(90pi^2).

The exact global form, divided by the unchanged period calB, is

    Qgeo[T_field,S_field]/calB
       =int_0^sqrt(Qmax) 2T(j) sum_(k in Z)
          A(2pi k/T(j)) conjugate(f_(j,k)) g_(j,k) dj.       (4)

This is a direct-integral identity, not a finite selection of flow modes.
The factor follows directly from dmu_SM=dz d tau dtheta and
dtheta=dj/sqrt(Q-j^2) on each radial branch. In particular the identity
operator gives 2 int T(j)dj=2pi z_h, the full physical-half angular measure
divided by calB. No physical normalization or admission target changes.

The flow Fourier integer k is NOT a thermal Matsubara mode. A(w) already
contains the entire massive angular tower. The original thermal sums and
state-dependent finite returns are retained in the physical H and its
comparison remainder. Equation (4) does not replace those sums by k.
All other local, LLL, scalar and mixed terms in H0 remain to be assembled.
A(0)=0 removes only this flow multiplier's constant mode, not its separately
retained constant contacts or either original matching equation.

There is an essential assembly identity: A already contains Q4geo. If H0
is assembled from the retained local Ha2 Hessian, add Ha2-Q4geo to (4),
not a second copy of Ha2's fourth term. The exact remainder is

    Ha2-Q4geo=C F_eff/(72pi) int Rbar |grad T|^2
             -C F_eff/(960pi) int (4Rbar^2+8Delta Rbar)T^2.

Its paired seam distribution and original boundary completion are retained,
as are the shared constant-contact uncertainty and every other local row.
This is the prior complete-restoration identity, not a new subtraction
choice. The global assembly must include this term exactly once.

For the actual mass/log core, T_field starts at z^3 or higher with the
stated finite log degree; its even extension is C2 and has weak third
derivatives in L2. The C2,1 metric gives a C1 flow with preserved Liouville
measure. Repeated differentiation through order three uses bounded Q,Q',Q''
and the physical one-sided field jets. At a regular bolt use its even disk
jets. Thus the base lift belongs to Dom(X^3) without a smooth-double theorem.
The existing full Bessel bound gives

    0<=sum_n d_n P_Mn[T]<=4C^2 ||X^3 T||^2/(525pi^2 q).

Tonelli for its positive terms and polarization justify (4) on this core;
the fixed fourth contact is then added in the common restoration. Form
closure extends only in the already proved joint form norm. This is not a
strong-domain or metric-variation theorem for a nonlinear neighborhood.

## 4. A failed estimator on the actual mass lift, with exact replacement

Retain J'=z^3/Q^2 and the exact mass lift
n_m=-2beta QJ/z, b_m=beta z^3/Q+beta Q'J. Its actual cutoff is one through
.1 and zero by .2. BOTH local reference rows cancel before the cutoff;
the transition remains f_cut=2xz chi_m-2F n_m Dchi_m and
alpha_cut=n_m Dchi_m, with the original regular-bolt matching. This is not
a global fixed-temperature zero mode. In particular b_m=z^3+O(z^5).

On each crossing orbit, a profile of mass coefficient m therefore has

    f_j(t)=m v^3 |t|^3+O(|t|^4 log^N|t|),
    v=sqrt(beta-j^2), jump(f_j''')=J3=12m v^3.             (5)

For the exact lift the regular remainder in (5) is O(|t|^5). Hence its
fourth distributional derivative contains J3 delta_0. For m!=0 a periodic
H4 Fourier-tail bound has an INFINITE right-hand norm, for almost every
crossing j. This is a verified failure of that estimator on the actual
input, not a divergent quantum form or a construction obstruction.

Use the periodic Bernoulli polynomial B4 instead:

    s_j(t)=-J3 T(j)^3 B4(t/T(j))/24,
    (s_j)_0=0, (s_j)_k=J3/[T(j)(2pi k/T(j))^4], k!=0.     (6)

The [Bernoulli Fourier identity](https://dlmf.nist.gov/24.8.E3) applies
including the periodic endpoints for n=4. Its third-derivative jump is
exactly J3 and its regular fourth derivative is -J3/T(j). Subtracting it
makes r_j=f_j-s_j an H4 periodic profile for the retained mass/log input.
This subtraction is on the flow bundle only. It changes no physical field,
boundary condition, mass coefficient or original fast norm. Its cross terms
with r_j must be kept. Trapped trajectories have no seam jump and remain
in r_j and (4). The exact implicit slow exponent is retained throughout.
Use (6) for the high-frequency tail identity; evaluate low-frequency
coefficients from the original profile. At j=0 the unit-mass k=1 seam
coefficient is about 2.0e3, although the cut profile is much smaller.
Treating that subtraction as an independently bounded low-frequency field
would introduce a large cancellation loss. It is not a new boundary basis.

## 5. Executed complete seam tail and the remaining full-profile allowance

Set a_flow=C/(90pi^2), Omega>=10, and

    X_Omega=-F(C Omega^2/4)-4log(mu0 B).

F is the existing full-tower digamma finite part, not a truncated mass
sum. Its monotonicity and retained |d[-F(exp s)]/ds|<=4 give, for w>=Omega,
chi(w)<=X_Omega+8log(w/Omega), and X_Omega>2 here.

A coefficient-sum bound gives |Q''|<12 on the WHOLE finite interval.
Apply the periodic Poincare inequality to dot z and ddot z=Q'/2: the
signed radial period is >=2pi/sqrt(6), so every base period is at least
pi/sqrt(6). Thus N=floor(Omega T/(2pi))>=Omega T/(4pi)>=1.
For the decreasing function w^-4[X_Omega+8log(2w/Omega)], integral
comparison sums every |k|>N. Integrating every crossing j, using
int_0^sqrt(beta) (beta-j^2)^3 dj=(16/35)beta^(7/2), proves

    Qgeo,>Omega[s]/calB <= m^2 C_seam(Omega),
    C_seam=36864 a_flow beta^(7/2)/(35pi Omega^3)
                                      *[X_Omega/3+8/9].       (7)

Actual coefficient enclosures give

| Omega | Certified upper guard for C_seam |
|---:|---:|
| 100 | 790.478560 |
| 10000 | .001004145567 |
| 1000000 | 1.219301126e-9 |

These are bounds on a specified quadratic-form tail DIVIDED BY calB, per
unit squared mass coefficient. They are NOT kappa_m, physical graph errors,
or entries in the .0012069 feedback budget. The original period and action
factors must be restored in the actual assembly. No improvement by changing
normalization is claimed.

For a full profile put J4=||X^4 T_field||_reg,SM/sqrt(calB), with only the
explicit seam delta removed. Then

    ||partial_t^4 r||_SM/sqrt(calB)
       <=J4+|m| sqrt[4608 sqrt(6) beta^(7/2)/(35pi)] =: R4,
    Qgeo,>Omega[T_field]/calB
       <=[|m|sqrt(C_seam)+sqrt(a_flow X_Omega) R4/Omega^2]^2. (8)

The second term follows from A(w)<=a_flow X_Omega w^8/Omega^4 for w>=Omega;
the square retains the cross term. J4 includes all actual cutoff derivatives
and is NOT assigned a value here. At each physical z it can be integrated
without tracing individual orbits: with A=Q-j^2,

    X^4 f=A^2 f''''+3A Q' f'''+[2A Q''+3(Q')^2/4]f''
                                      +[A Q'''/2+Q'Q''/4]f'.

Angular moments are int cos^(2n)theta dtheta=2pi binom(2n,n)/4^n.
This supplies a concrete full-profile tail calculation, including the fast
remainder. Its J4 bounds and the finite-frequency part of (4) remain to be
computed for the coupled trial and residual. No finite set of modes is
reported as the complete operator.

## 6. Decision and unchanged requirements

RESEARCH DECISION: retain the construction; sharpen the global comparison
implementation using (4) and the exact seam subtraction (6), including all
trapped orbits. Reject a smooth-periodic H4 tail certificate for the unmodified
mass lift. The implementation now has certified orbit geometry and the
complete seam-tail bound (7); it is not a full H0 solver or a strong inverse.

ONE NEXT CALCULATION: evaluate the full coupled comparison return using (4),
(8), all original local/matching rows and the physical comparison remainder,
and certify its ORIGINAL strong graph residual. This remains PROPOSED, NOT
EXECUTED. It is the necessary next test of whether this route is quantitatively
attainable; another orbit or seam-tail optimization would not answer it.

No kappa_m or other post-elimination entry is measured. Conditional
kappa_m=5e-5 still spends 2.14765e-12 affine and .0024731 graph feedback,
leaving 5.34235e-12 and .0012069 for ALL remaining entries. Independent
energy feedback remains; none of the other entries is shown to fit.
The >.00906808164 original mass conversion still requires
|m|<2.20553814952e-9 if both 1e-11 graph ceilings are retained.

The admitted-family fast screen .028913+1e-11(1.0706e10)=.135973<.1599164
has only its two-fast-input scope. It controls neither independent phase
variations nor the unknown physical and nonlinear-neighborhood comparisons.
Still unproved are Ecenter+1e-8 alpha<=1e-11,
R_X alpha*.028913+Z0_other<=.25,
R_X alpha^2*1.0706e10+Z2_other<=5e7,
Y0+1e-8 Z0+5e-17 Z2<=1e-8, and Z0+1e-8 Z2<1.

The [local nonsmooth symbol theorem](https://web.math.ku.dk/~grubb/CPDE14.pdf),
Theorem 2.3, requires an identified symbol and |s|<tau; it does not supply
these graph constants. The [smooth boundary logarithm theorem](https://arxiv.org/pdf/math/0703878),
Theorem 2.4 and Proposition 2.9, is not a high-regularity theorem for this
C2,1 geodesic compression. The ordinary flow/third-jet and distributional
arguments above are used instead, only in their stated scope.

Remaining to-do list: full coupled global return and strong residual;
actual correction and original physical-norm comparisons; nonlinear
admission/contraction and exact compatibility; physical matching/causal
realization; semiclassical validity/noise; all-sector stability; charge/return,
opening, transport, preparation, separate creation and the full objective.
