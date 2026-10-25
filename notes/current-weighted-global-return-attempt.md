# Attempted weighted coupled global return: domain, ultraviolet and norm losses

2026-10-24, E003/S169. Ordinary argument before supporting verification.
Full stable semiclassical wormhole objective ACTIVE/INCOMPLETE; the closed
feasibility audit is unchanged.

OUTCOME: UNRESOLVED complete coupled return. This calculation repairs the
proposed weight's domain behavior and bounds its actual massive propagators.
It also identifies an infinite absolute-contact estimator and quantifies
the failure of a bounded-kernel similarity estimator. It DOES NOT compute
the renormalized coupled inverse or any post-elimination mass column.
Neither estimator failure is an obstruction to the actual construction.

## 1. Retained global problem and the first unsupported implication

Retain H=I_ref''-I_pot''+I_sq''+Gamma_q'', both original equations, their
original boundary completion, fixed temperature and regular bolt. Retain
the actual paired first-order domains and complementary adjoints, composed
square, mass cubic, forced quartic logarithms, correlated slow pair, common
regulator, all angular and thermal modes, contacts and finite returns.
The fixed-seed weak inverse in the positive-diagonal form norm is available.
It is not a quantified strong interior/bolt inverse or physical trace map.

The [preceding coupled calculation](current-coupled-global-return-attempt.md)
remains the reference: Phi_m=chi_m(n_m,b_m), J'=z^3/Q^2,
n_m=-2 beta Q J/z, b_m=beta z^3/Q+beta Q'J. Locally its reference return
is f=2xz, alpha=0, so BOTH reference rows vanish. The actual cut lift has
f_cut=chi_m 2xz-2F n_m Dchi_m and alpha_cut=n_m Dchi_m. The original F_l/G
prefactor gives (G/F_l)int rho dr=-xB and zero bolt shift. The cut scalar
quartic and correlated slow profiles have zero net reference mass/bolt
shift. All transition terms remain. The uncut temperature sensitivity
about -6.51517e65 is not a fixed-temperature zero mode or a kappa_m entry.

Use the complete block system, with matching on overlapping restrictions,

    [L_oo L_oI L_oc] [r_o]   [f_o]
    [L_Io L_II L_Ic] [r_I] = [f_I],
    [L_bo L_bI L_bc] [ c ]   [f_b], c=(m,c4,Delta).             (1)

Its phase rows include Ghat(3)=0, the original scalar finite part with
coefficient 20 beta Fstar=3.33333329399556..., and the original slow
residue. The known correlated B0(z^4-z^(4+e_*)) and forced log remain in
the affine source. There is no extra wall or freely chosen quantum stress.

The first unsupported implication remains a quantitative strong return
for the entire interior block and its phase traces. An outer-only norm
cannot supply it: the actual interior kernel exhibited previously remains.
This turn attempts a weighted estimate of L_II and the Schur products
L_xI L_II^-1 L_Iy in (1). The calculation below reaches actual propagator
bounds but stops at the full constrained, renormalized commutator. A gap
bound for Dsharp is not a bound for that Hessian commutator.

## 2. An admissible weight and the required TT correction

Put z=B/r, z_h=1/(1+h), h=5.057799424826589...e-34. On the actual seed

    Q=(1-z)^2 P(z)-tau z^3,
    P=beta(1+2z+3z^2)+u z^2+eps z^4,
    6 beta+u+eps=1, tau>0, eps<0.

Take s=7/5 and a fixed monotone theta: zero through z=.1, one from z=.2,
with transition I_(10z-1)(11,11). Set

    omega(z)=s int_0^z theta(t)/(h+1-t) dt, W=exp omega.       (2)

W is constant at the paired seam, C11 or better at the transition, and
smooth and even in regular bolt coordinates. It preserves the original
paired boundary traces and composed-square transmission conditions.
Multiplication by W or its inverse preserves the mass/log class because
W=1 near z=0. At fixed h both are bounded. This is an auxiliary norm
weight on the actual finite domain, not a surrogate metric or domain.

The simpler (h+1-z)^(-s) reflected evenly would have a cusp at the seam:
it preserves the H1 first-order domain but need not preserve the H2
transmission domain of the composed square. Formula (2) removes this
problem rather than applying a smooth-domain theorem to the cusp.

There is a separate constraint issue. For a TT field (b,xi),

    (r^2 F b)'=-(r^2 F'-2rF)xi.

Multiplying both components by W leaves the constraint defect W' r^2 F b.
The exact bolt-regular correction, keeping xi_W=W xi, is

    delta b(r)=-(r^2 F)^-1 int_R^r W'(t)t^2 F(t)b(t) dt.       (3)

It is regular at the bolt (delta b=O(r-R)) and has the genuine outer tail

    lim_(r->infinity) r^4 delta b
          =-lambda^-1 int_R^infinity W' r^2 F b dr.           (4)

The integral is supported away from the outer boundary but is not
identically zero. Thus a weighted TT estimate must return this quartic
phase contribution and its energy/graph cost. Equations (3)-(4) repair the
constraint algebra; they do not bound the resulting York operator. Using
the old TT coercivity on (Wb,Wxi) without (3) is invalid. In canonical
coordinates the same issue appears as noncommutation with the York split.
It precedes numerical interpretation of a weighted quantum estimate.

## 3. Executed bound for the actual massive propagators

Since P<=1+2|eps|, Q<= (1-z)^2(1+2|eps|), so in the dimensionless orbit

    |grad omega| <= L=s sqrt(1+2|eps|)<1.400000000001.        (5)

Use the actual anti-self-adjoint Dsharp, not a cap or a smooth extension.
The Leibniz identity on its paired domain gives

    W(Dsharp+M)W^-1=Dsharp+M-c(domega),
    ||c(domega)||<=L.

The original inverse has norm <=1/M. The Neumann series therefore proves
on the actual global L2 domain, for EVERY n>=1,

    ||W(Dsharp+M_n)^-1 W^-1||<=1/(M_n-L),
    M_n>18 sqrt(n), inflation <=18/(18-L)<1.084338.           (6)

The full geometry, period and all its returns are already in Dsharp.
No thermal sum or angular mode is truncated. This elementary bounded
perturbation argument needs only the retained skew-adjoint realization,
the Clifford norm and domain-preserving multipliers. It does not assume
smooth reflected mass/log coefficients. It does not apply to the massless
LLL, whose separate anomaly and global-state terms remain in H. It is a
seed result; a neighborhood metric comparison has not been proved here.

For the retained geodesic comparison, use the exact positive representation
of K_nu for real positive argument in [DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9).
Writing K1(t)^2-K0(t)^2 as a double integral has nonnegative weight
cosh v cosh w-1 and exponent -t(cosh v+cosh w). Since cosh v+cosh w>=2,

    exp(Lr)[K1(Mr)^2-K0(Mr)^2]
       <=[K1((M-L/2)r)^2-K0((M-L/2)r)^2].                  (7)

With the ORIGINAL amplitude M^4/(2pi^2) in B_M, variable substitution
in the exact sixth moment gives

    int r^7 exp(Lr)|B_M(r)| dr
       <=48/(35pi^2 M^4) (1-L/(2M))^-8.

Sum the COMPLETE tower, using sum d_n/M_n^4<4C^2/q and C/q<.0035:

    C^-1 sum d_n int r^7 exp(Lr)|B_Mn(r)| dr
       <.002671854.                                       (8)

The inflation is <1.373445. Liouville invariance and
|omega(Phi_r x)-omega(x)|<=Lr apply to every complete geodesic, including
returns and conjugate points. Formula (8) concerns the comparison kernel;
it does not identify that kernel with the actual quantum Hessian.

## 4. Why these bounds do not close the quantum/global estimate

At a common positive regulator, the complete symmetric conjugation error
is the form of

    C_W=(W H W^-1+W^-1 H W)/2-H.                           (9)

On the fixed-seed strong core the retained domain/action argument defines
its common restored limit. What is missing is a useful operator bound in
the required constrained strong and phase norms, uniformly in the intended
neighborhood. The existence of this fixed-seed form is not in question
merely because the following estimator fails.

In particular, the geodesic integrand is
g(rP)=1-r^2 P^2/2+r^4 P^4/24-cos(rP). Its polynomial contacts cannot be
bounded separately after the angular sum. The exact moments give

    sum d_n |I2|=(6pi^2)^-1 sum 2(q+2n)=infinity,
    sum d_n |I4|=(4/(15pi^2)) sum d_n/M_n^2=infinity.       (10)

The first partial sum is exactly [2qN+2N(N+1)]/(6pi^2).
For n>=q, d_n/M_n^2>=2C/n, proving the second logarithmic divergence.
For a nonconstant omega the second-order conjugation already contributes
-|Xomega|^2, so the first failed bound is not repaired by declaring that
polynomial commutator zero. Cancellation with the translated cosine and
the SAME local restoration must precede the angular limit. Equation (8)
controls a sixth remainder only AFTER the needed subtraction; by itself
it supplies neither an H2 logarithmic bound nor the constrained trace maps.
This verifies failure of the specified absolute-contact estimator, not
divergence of the physical renormalized functional.

There is also a quantified failure of a finite but excessive estimator.
The existing global row-sum-subtracted gluing form has L2 norm <=6e24 C.
The retained TT conversion is C J0/E0<5.2e-29. Blind similarity gives the
certificate 3.12e-4 cond(W) in E0 units. From (2), exactly,

    [(h+.8)/(h+1-z_h)]^s <=cond(W)<=[(h+.9)/(h+1-z_h)]^s.

These bounds are 1.141136887...e46 and 1.345709367...e46. Thus the NUMBER
produced by that certificate exceeds 3.56034e42. This is a lower bound
on the proposed upper-bound expression, NOT a lower bound on the actual
gluing operator or solution. It plainly cannot certify a .001 allowance.

A .001 allowance for this one gluing term would require a weighted row
bound <=1.923076923e25 C (round down for a sufficient decimal target),
or inflation <=3.205128205 of the existing 6e24 C. The actual-propagator
inflation in (6) is modest; even its fourth power is <1.382477. That makes
an exponentially weighted reconstruction quantitatively plausible, but
the gluing expression, its diagonal row correction and local restoration
have not been reconstructed here. Multiplying 6e24 by this fourth power
without proving the requisite factorization would be unjustified.

An explicit provisional quantum-form target is

    |C_W,quantum[T]| <= C J2[T]/(200pi)
                        +1e23 C J1*[T]+1e25 C J0*[T].      (11)

The original C/(100pi) fourth margin and conversions
C J1*/E0<1e-26, C J0*/E0<2e-28 would leave half the fourth margin and
spend at most .003 E0 on the last two terms. The old negative quantum
budget .003601000003 is below .01 by .006398999997. This is a quantified
target for ONE ingredient, not a new theorem or a global-inverse budget.
The classical weighted block, (3), the mixed scalar block, LLL, Weyl terms,
and phase traces must also be controlled. The paired mass/log estimates
currently available do not provide those constants.

[Grubb's Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf) assumes an
identified C^tau S^d_(1,0) symbol, |s|<tau, and N>dimension/2+1. It supplies
neither the missing symbol seminorms for (9), global phase traces, nor a
uniform inverse on this paired family. We do not invoke it to fill these
gaps. No smooth boundary-log theorem or finite-mode approximation is used.

## 5. Complete budgets remain conditional, in the original units

The actual Schur blocks are Lbar_xy=L_xy-L_xI L_II^-1 L_Iy, with the
analogous affine return. Every mass/scalar/slow/fast column receives this
same elimination. No such complete column has been bounded in this turn.
For actual post-elimination bounds k_j and proved c_j<=a_j+b_j E, one needs

    S_other+sum k_j a_j<=7.49e-12,
    eta_other+sum k_j b_j<=.00368, j=m,4,Delta,I,             (12)

plus the separately allocated known source <=1e-14 and total feedback
<=.25. At k_m=5e-5, the retained mass affine allowance 4.2953e-8 spends
2.14765e-12, leaving 5.34235e-12 for ALL remaining terms. The rounded
graph slope 49.462 spends .0024731, leaving .0012069. Using the earlier
unrounded slope gives .00120692478. These are conditional reserves, not
measurements that the other entries fit. All the other entries are unknown.

Differences also require the independent energy contribution
|delta m|<=390.478 E_phys+49.462 E_lifted (same forced leading data).
The mass target then adds about .0195239 times the energy increment.
The previous mass-only restricted-family target E_phys/E_lifted<=.0618179
or equal-independent-input column target k_m<=8.3647e-6 remains conditional;
other entries make both stricter. No such restricted-family comparison is
proved. The old source-dependent mass allowance is not a Lipschitz bound.

Changing the auxiliary weight does not change g_old=g_lifted+m d_chi,
||d_chi||H3>.00906808164. Preserving both original 1e-11 graph ceilings still
requires actual |m|<2.205539e-9. The conditional 4.345e-8 allowance is not
the actual mass. The ~19.70 gap is unchanged; no normalization improvement
is claimed.

The two-fast-input screen remains .028913+1e-11(1.0706e10)=.135973
<.1599164, reserve .0239434, solely within its proved admitted-family
scope. It does not bound independent boundary-phase, interior or mixed
input variations. The original physical comparisons remain unproved:

    Ecenter+1e-8 alpha<=1e-11,
    R_X alpha*.028913+Z0_other<=.25,
    R_X alpha^2*1.0706e10+Z2_other<=5e7.

Consequently Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 are still
unproved. The outer-only R_X is impossible; an augmented global return
has not been certified. Nothing here turns the seed estimate into a
neighborhood estimate.

## 6. Decision boundary

RESEARCH DECISION: retain the construction and exact cut reference lift;
sharpen the SIGNED coupled weighted estimate. Park blind similarity and
absolute-contact bounds: one loses over 42 orders of magnitude and the
other diverges. Actual massive propagation is compatible with a modest
weighted loss, but this does not decide the full feasibility question.
Prioritize the constraint/phase return and common-restored commutator,
not further cold/scalar source optimization or an isolated contact.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): evaluate and bound the
complete constrained commutator (9), including the York return (3)-(4)
and original phase rows, in the bordered system (1). Perform the common
restoration BEFORE taking absolute values or removing the regulator.
Test (11) as its quantum sub-budget and return the resulting coupled
resolvent/phase constants to (12) and the original physical comparisons.
This must be one coupled operator calculation, not a product-model test.

Executed: the admissible weight, exact TT correction/outer tail, all-mass
actual propagator bound, complete tilted comparison moment, divergence
of the specified absolute-contact estimator, and numerical failure and
replacement threshold for blind gluing similarity. Proposed/unexecuted:
the complete restored constrained commutator, actual coupled inverse and
phase coefficients. No certified k_m, actual correction, preserved nonlinear
ball, exact compatibility, physical matching, semiclassical validity/noise
or all-sector stability follows. No construction obstruction is established.

Remaining to-do list: full constrained common-restored global return and phase rows; actual correction and original physical norm comparisons; nonlinear admission/contraction and compatibility; physical matching/causal realization; semiclassical validity/noise; all-sector stability; charge/return, opening, transport, preparation, separate creation and the full stable-wormhole objective.
