# A three-factor polynomial inverse for the complete coupled return

2026-10-25 HST, E003/S169. Ordinary argument and execution contract BEFORE
supporting calculation. Full stable semiclassical wormhole objective ACTIVE
and INCOMPLETE. The closed feasibility audit is unchanged.

OUTCOME TO TEST: an exact-operator polynomial improves the conditional error
allowances without changing any norm. The ORIGINAL strong defect, numerical
output constants, signed solution and actual mass column remain UNRESOLVED.
No numerical application of the full quantum inverse is claimed.

## 1. Actual spaces, domain, and the first unsupported step

Retain the original two-row Hessian and York realization:

 H=[-A,-C;-C*,D], K=C^-1 A(C*)^-1, S=I+KD,
 (nu,T)=R(s,xi)=(s+Uxi,Vxi), R*HR=[-a,b;b*,ell],
 J(0)=0, (Q^2 J'/z^4)(z_h)=0.

These are the actual paired first-order Dirac domain and complementary
adjoint, one-sided mass cubic and finite logarithmic phases, fixed clock
and temperature, and regular rotational bolt. The reflected cubic is only
C2,1. The common regulator is removed in the retained closed forms before
the following algebra. ell is the COMPLETE positive form: exact UV logarithm,
all angular and thermal modes, coincidence and density contacts, cutoff
transitions, and finite-domain returns. No finite mode matrix replaces it.

Write B=b*a^-1b, J_F=V* K^-1, and V_ell for the ell form domain. The fixed
seed realization gives B<=q ell as forms, q=.338. The candidate P_1 and its
input/output spaces are

 P_1=V ell^-1 J_F : Y_F -> X_F,
 Y_F=J_F^-1(V_ell*), X_F=V(V_ell),
 ||y||Y_F=||J_F y||V_ell*, ||Vxi||X_F=ell(xi,xi)^(1/2).

Use completions of the common fixed-gauge core and the retained inverse
transports. These are FORM norms, not the original strong source norm Y
or physical admission norm X. The first unsupported quantitative step is
still the full original-source/output comparison. Neither a smaller form
defect nor a general smooth-manifold theorem supplies it.

In particular, Grubb Theorem 2.3 requires the relevant symbol seminorms and
|s|<tau; it does not supply our differentiated kernel constants. Moller
Theorem 1 concerns a closed smooth spin manifold with stable kernel. It is
not a theorem on the global paired mass/log realization. Both sources were
rechecked for this attempt; their limited local use remains as in the
retained ordinary argument:
https://web.math.ku.dk/~grubb/CPDE14.pdf
https://arxiv.org/html/0705.3857v3

The attempted disk-jet bootstrap did NOT produce the missing constants.
The existing interior 1e48 C estimate is H3->L2, not differentiated output.
A C5-to-point-source estimate and holomorphy in metric parameters do not
by themselves justify a spatial-derivative operator norm or its common
regulator limit. Finite translation jets with interpolation are a possible
repair, but the full actual-return estimate has not been established.
No numerical value from that proposal is used below.

## 2. Exact signed algebra, including an erratum

Let T=ell^-1 B act on V_ell with its ell inner product. It is bounded,
self-adjoint and nonnegative, with spectrum in [0,q]. This follows directly
from the Riesz-form identities and B<=q ell; no differential boundary
spectral theorem is needed. Retain its similarity on the dual space.

 J_F S V=ell+B=ell(I+T),
 E=I-SP_1=-J_F^-1 B ell^-1 J_F,
 F=(ell+B)^-1 B=(I+T)^-1 T,
 Q_T=(ell+B)^-1=(I+T)^-1 ell^-1.

Consequently the general power identity is

 P_1 E^p=(-1)^p V F^p(I-F)^(-p-1)Q_T J_F.             (1)

ERRATUM: the general-p formula in current-schur-physical-output-test.md
omits (-1)^p. Its p=2 specialization and all its absolute norm allowances
are unaffected. For p=3 the sign MUST be retained in a signed solve.
The frozen earlier artifact is preserved; this note records the correction.

## 3. Explicit degree-32 candidate, with three residual factors

Use N=32, m=N-3=29 and Chebyshev polynomials defined by
T_0(x)=1, T_1(x)=x, T_(j+1)(x)=2xT_j(x)-T_(j-1)(x).
To avoid confusion with the operator T, denote these polynomials by Ch_j.
Set

 x(t)=2t/q-1, x_minus=-1-2/q,
 d=Ch_29(x_minus),
 r(t)=-t^3 Ch_29(x(t))/d,
 pi_31(t)=(1-r(t))/(1+t),
 P_32=V pi_31(T) ell^-1 J_F : Y_F -> X_F.             (2)

r(-1)=1, so pi_31 is an EXACT polynomial of degree 31, not a rational
operator concealing S^-1. Evaluation needs applications of the full
ell^-1 and B. None of those applications has been numerically performed.
The retained realization suffices to define (2) on the actual form spaces.

The recurrence, or cos(j theta), proves |Ch_j(x)|<=1 for -1<=x<=1.
For x0=1+2/q>1, D=|d|=Ch_29(x0)>0. Thus the COMPLETE form residual obeys

 I-SP_32=J_F^-1 ell r(T) ell^-1 J_F,
 ||I-SP_32||Y_F<=q^3/D.                              (3)

This is an infinite-dimensional spectral interval proof using the full
operator; it does not sample its physical modes or discard its UV term.
The improvement over a geometric series is cancellation in an explicit
polynomial. It is not a new normalization or a new boundary basis.

The root t^3 is essential for the selected physical certificate. The
retained actual-domain test rejects two residual factors for physical
areal C5 output. Three factors give fixed-seed qualitative regularity,
with numerical graph/physical/phase constants still unknown. An ordinary
unweighted Chebyshev residual need not have this required zero at t=0;
we do not infer its form-to-physical boundedness.

## 4. Original output maps and the full coupled inequalities

Use the ORIGINAL energy E(xi)=ell(xi,xi)^(1/2)/sqrt(|W| beta), where
|W|=N_halves pi B^2/(G gamma); do not change its units. For Schur errors,
both original equations require

 R_S T=((C*)^-1 D T,T), alpha=nu+2T-r partial_r nu.

Define C_i,3=||O_i R_S V F^3||E->Z_i for EACH original output: graph,
physical areal C5, mass, scalar quartic, correlated slow pair and matching
returns. The full logarithm in (C*)^-1 D is retained. Each numerical
C_i,3 remains UNKNOWN. For a source y put u0=ell^-1 J_F y. Exact algebra
gives the error

 (S^-1-P_32)y
   =-V F^3 (I+T)^2 Ch_29(x(T)) u0/d.                 (4)

Hence, with h2=(1+q)^2/D and h3=(1+q)^3/D,

 ||O_i R_S(S^-1-P_32)y|| <= C_i,3 h2 E(u0).         (5)

The known COMPLETE affine source comparison is
E(u0_aff)<=(1+q)1.1e-10. Therefore its sufficient tail bound is

 error_i_aff <= C_i,3 (1.1e-10) h3.                 (6)

This uses the retained coupled source conversion, not an assumed unit
column. It is conservative; it does not assign the actual solution a
newly measured strong norm. For an independent source/phase column j,
I_j=E(ell^-1 J_F y_j) must instead be bounded separately. With certified
trial entry c_hat_ij and numerical application error e_ij, the complete
entry requirement is

 |c_ij| <= |c_hat_ij|+e_ij+C_i,3 I_j h2.             (7)

All columns j are retained: mass, scalar quartic, correlated slow pair,
fast remainder, and the infinite matching/return quotient. The latter
must be controlled as an operator, not three sampled phase observations.
No c_hat, e or I is assigned zero or an affine-source bound.

For the requested ORIGINAL strong residual let

 C_def,3=||J_F^-1 ell F^3||E->Y,
 I_Y=||ell^-1 J_F||Y->E.

The first expression is initially on the common core; bounded extension
is an additional required argument, not asserted here. Then the exact
sufficient screen is

 delta_Y <= C_def,3 I_Y h3 <=1/2.                    (8)

This is the missing original-source comparison. A certified (3) alone
does not close (8), original physical admission, or contraction.

For the actual post-elimination mass column the sufficient condition is

 |kappa_hat_m|+e_m+C_mass,3 I_m h2 <=5e-5.            (9)

Neither kappa_hat_m nor the actual kappa_m is computed. IF (9) holds,
the retained mass debits are 2.14765e-12 affine, .0024731 feedback and
.0195239 energy. The remaining entries must JOINTLY satisfy

 Y_other+4.084e-22 A_o+Y_phase <=5.34235e-12,
 eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase <.0012069. (10)

Use (7) for every contribution to those quantities, including signed trial
centers and all numerical errors. A full-margin conditional graph ceiling
from (6) is 5.34235e-12/[(1.1e-10)h3]. An illustrative physical ceiling
is 5e-9/[(1.1e-10)h3], ONLY IF that whole physical allowance is available.
The all-zero-other-debits feedback product ceiling is .0012069/h2; the
analogous mass product ceiling is 5e-5/h2. These are allowances, not actual
coefficients or evidence that remaining entries fit. Compare the same
graph/physical units with the previous geometric prefactor
(1+q)^5(1.1e-10)q^29/(1-q). No signed center is improved by this bound alone.

The mass lift retains J_m=int_0^z t^3/Q^2, n_m=-2beta QJ_m/z,
b_m=beta z^3/Q+beta Q'J_m. Its uncut reference cancellation is unchanged;
the cutoff is one through .1 and zero by .2. The transition terms
f_cut=2xz chi_m-2F n_m Dchi_m and alpha_cut=n_m Dchi_m, original boundary
conditions and regular bolt are kept. It is NOT a fixed-temperature zero
mode. The original fast conversion 2133<||d_chi||H3<2174 and the LINEAR
mass debit <=3.538254568e-12 remain. The scalar quartic row 20beta Fstar
and the exact correlated slow root remain in (7), not projected away.

The two-fast/admitted-family screen stays
.028913+1e-11*1.0706e10=.135973<.1599164. It supplies no independent
phase bound or nonlinear-neighborhood transfer. Original nonlinear gates:

 E_center+1e-8 alpha_embed<=1e-11,
 R_X alpha_embed*.028913+Z0_other<=.25,
 R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
 Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

All remain unproved, as do exact compatibility, physical matching,
semiclassical validity/noise, all-sector stability and the full objective.

## 5. Decision, execution boundary, and verification contract

If the supporting arithmetic confirms (2)-(8), change the residual
polynomial from the geometric series to (2); keep the construction and
three-factor physical output. This is a quantified improvement of the
COMPLETE estimator, not a passed original-space certificate. The radial
energy-only commutator shortcut remains rejected. No unknown constant
becomes known by increasing the polynomial degree.

The strongest executed result is the certified form bound and conditional
original-space error formula. It closes only the form residual target.
The first missing actual strong comparison remains unresolved, not a
demonstrably excessive bound and not an obstruction to the construction.
The failed radial shortcut and two-factor physical estimator are separate
verified failures. Fixed-seed qualitative finiteness is not quantitative
evidence that the new conditional ceilings are attained.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): certify the complete triple-
residual output operator, including its original-source defect map and
independent input/phase comparisons, against (6)-(10), on the actual paired
domain. Local Cartesian finite-jet interpolation may be used only after
its full-return/common-regulator argument is supplied. Do not apply another
isolated source or contact bound in place of that calculation.

Freeze this proof/checker/review before arithmetic. Check exact polynomial
divisibility, identities and the odd-power sign; enclose D and the whole
spectral interval bound at 512 bits; evaluate all displayed conditional
ceilings and their ratios. Retain every development failure, review the
ordinary argument, and run unchanged acceptance. Checks support, rather
than prove, the domain and spectral argument. No physical inverse is run.

Remaining to-do list: complete coupled original graph/physical/phase bound
and signed solution; original source and neighborhood comparisons;
nonlinear admission/contraction; exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; opening, transport,
preparation, separate creation and the full stable semiclassical wormhole
objective.
