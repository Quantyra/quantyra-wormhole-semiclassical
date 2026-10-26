# All-trial test of the complete P32 adjoint residual estimator

2026-10-26 HST, E003/S169. Ordinary argument frozen before numerical
execution. Previous goal turn: NO PROGRESS. Retain the stable semiclassical
wormhole objective, accepted 1ee884b, complete P32 and the closed audit.

OUTCOME CONTRACT: test a SPECIFIED positive residual estimator for every
trial. Its failure would not imply a lower bound on the actual C_b,3,
residual, correction or quantum operator error. The physical norm and
signed coupled columns remain uncomputed. This is not another contact
calculation: the lower bound below applies to the COMPLETE reported
P32 physical bolt-tail allowance, with every other debit set to zero.

## 1. Actual operators, domain and normalization

Use the unchanged physical blocks and fixed York return

 H=[-A,-C;-C*,D], R(s,xi)=(s+Uxi,Vxi),
 R*HR=[-A,b;b*,ell], B=b*A^-1b,
 N=|W| beta, e=ell/N, T=ell^-1B, 0<=T<=q=.338,
 F=(I+T)^-1T, f=q/(1+q).

The original Schur operator and polynomial are unchanged:
K=C^-1 A(C*)^-1, S=I+KD, J_F=V*K^-1,
P32=V pi31(T) ell^-1 J_F:Y_F->X_F,
Y_F=J_F^-1(V_ell*), X_F=V(V_ell),
r(t)=-t^3 Ch29(2t/q-1)/Ch29(-1-2/q), pi31(t)=(1-r(t))/(1+t).
The factor (1+t) divides exactly. These transported form spaces are not
the original strong Y and physical X. With the ORIGINAL reconstruction
R_S w=((C*)^-1 D w,w), the complete maps are
G_i=O_i R_S V F^3:E->Z_i and G_def=J_F^-1 ell F^3:E->Y.
Their outer returns, inputs I_j=E(ell^-1 J_F y_j), and numerical original
norm comparisons remain unknown. Only the affine input bound
I_aff<=(1+q)1.1e-10 and the seed form spectral interval are known.

E is the actual fixed-seed form space with inner product e. Its dual is
E*, not the original strong source Y. Let L=(ell+B)/N, viewed as E->E*,
and let L0=(ell0+B)/N be the complete geodesic comparison. Their difference
Rop=(ell-ell0)/N obeys ||Rop||_(E->E*)<=eta=.0022. The actual scalar and
mixed blocks are identical in the two operators. Thus

 ||L0||_(E->E*)<=M=1+q+eta.                         (1)

All forms have the original paired first-order matter realization,
complementary adjoint, common regulator and mass/log end data. The exact
UV logarithm, all angular/thermal modes, coincidence and density contacts,
cutoff transitions and regular-bolt returns are retained. L0 contains the
complete geodesic-flow form; it is not a flat scalar spectral surrogate.
In particular J(0)=0 and (Q^2 J'/z^4)(z_h)=0, with the original fixed
clock, temperature and scale. The reflected mass cubic remains C2,1.
The three residual equations index three inverse applications, not physical
modes. None of these fixed-seed statements supplies a nonlinear neighborhood.

For the physical bolt row define

 g(v)=ell_b[(C*)^-1 D Vv](z_h), ell_b=(1+h)/h,
 chi=V^-1 C^-1(ell_b delta_b), chi_N=N chi,
 j(v)=ell_b[U v+A^-1 b v](z_h)
     =ell_b[v(z_h)/2+(A^-1 b v)(z_h)].                (2)

All inverses and deltas in (2) have their original physical action
normalization and boundary conditions. The low-order bolt functionals
are in E* by the retained form-domain argument. The mixed Green return
is distributional at the bolt and keeps its end phase and finite moment.
It is not a smooth homogeneous adjoint. We initially use the common core
and the accepted fixed-seed closures; no fifth-derivative delta is put in E*.

In particular chi belongs to E, not necessarily to the strong operator
domain. The existing mixed Green realization has rho^2 log(rho) leading
behavior at the smooth disk center. Its local order-two logarithmic form
energy is finite (the high-frequency squared form integrand is of order
|k|^-4 log|k| in two dimensions). Its finite end/bolt returns use the
accepted paired closure. Thus v in (4) is in E, and the already proved
form comparison extends to it by core density. This does not extend that
comparison to fifth-derivative delta sources or move the mass/log seam.

The identities b=-AU-CV and
ell+B=V*(D+C*A^-1C)V give the EXACT correlated source relation

                         g=L chi_N+j.               (3)

Indeed L chi_N=V*D C^-1(ell_b delta_b)+V*C*A^-1(ell_b delta_b).
The first term is g. The second pairs with a test w as
ell_b[A^-1 C Vw](z_h)=-ell_b[Uw+A^-1 bw](z_h)=-j(w).
The N in chi_N cancels the N in L; no action normalization is dropped.

In particular j is the effective physical bolt functional of the existing
scalar Schur calculation. It is NOT g, and its lower bound is NOT a lower
bound for C_b,3. Equation (3) is the reason it can be used in this
estimator test without making that invalid transfer.

## 2. Retain the cancellation before testing the uniform estimate

For any adjoint trial p3 put v=p3-chi_N. Its third residual is

 s3=g-L p3=j-L v=s30-Rop v, s30=j-L0 v.             (4)

This already cancels the common quantum action on chi_N and p3 before
norms. Bounding their norms separately would be worse. Suppose the
chosen implementation uses the available UNIFORM form comparison to
report the positive bound

 D3=d0+eta V0,
 d0>=||s30||E*, V0>=||v||E.                         (5)

This is the estimator being tested, not a claim that the actual residual
equals its upper bound. Even granting exact evaluation of d0 and V0,
the triangle inequality and (1) imply

 ||j||E*<=d0+M V0,
 D3>= (eta/M)||j||E*, since eta<M.                  (6)

The latter follows directly from
d0+eta V0 >= (eta/M)(d0+M V0). It holds for EVERY trial, including
p3=chi_N, which eliminates Rop v but leaves s30=j. Thus changing a trial
cannot evade (6) while retaining estimator (5). Additional certified
quadrature, contact or source errors only increase its positive report.

The first missing numerical step of the DIRECT method remains a bound
on the actual correlated residual, or its transferred version below.
The uniform form theorem DOES apply to (4) on E, so (5) is well-defined.
Its numerical feasibility can therefore be tested without postulating a
new rough-boundary strong estimate or computing an inverse. This is a
different outcome from an undefined operator or a merely unknown constant.

## 3. Original-unit source lower bound from an accepted actual-domain test

Use the genuine localized TT field of
[the accepted sensitivity calculation](current-localized-bolt-adjoint-sensitivity.md).
It has J=psi/h, psi=(1-t/a)^8 on 0<=t<=a=1e-4, zero outside, with the
exact correlated K term in xi. It satisfies both York endpoint conditions,
is C6 at the join and regular at the bolt, and vanishes near the mass/log
end. Its full scalar return was bounded using the ACTUAL global scalar
inverse. Its energy bound includes all quantum terms and finite returns.
No model field or finite physical-mode subset replaces that test.

The frozen receipt and ordinary proof establish

 |j(xi_test)|>=ell_b b_low,
 ell(xi_test,xi_test)/|W|<=L_c,
 b_low=0.997806941653159481462...,
 L_c=6.325574662940601861460...e-30.                  (7)

Use the retained outward intervals, not these printed midpoints. Since
E(xi_test)^2<=L_c/beta, the original E-dual lower bound is

 Jstar=ell_b b_low sqrt(beta/L_c) <= ||j||E*.         (8)

The physical |W| cancels in this RATIO because E=ell^(1/2)/sqrt(N).
This is an explicit norm conversion, not a new physical normalization.
It agrees with multiplying the prior physical ell-dual lower bound by
sqrt(|W| beta), for either one or two folded exteriors. Both conversions
will be checked. No lower bound on the actual affine source is used.

## 4. Test the complete residual-to-output estimator

The accepted direct adjoint system has A3 diagonal I+T and upper/lower
shift as appropriate. For residuals s1,s2,s3, with Riesz_E denoted R_E,
its EXACT error transfer is

 T(p1-p1_trial)=F R_E(s1)+F^2 R_E(s2)+F^3 R_E(s3).  (9)

This follows by triangular substitution, without commuting T through an
unrelated operator. With Riesz trial a_trial and its residual tau, it gives

 C_b,3 <= ||a_trial||E+tau
              +||F R_E(s1)+F^2 R_E(s2)+F^3 R_E(s3)||E. (10)

The older .338*sigma bound is valid but loses the powers in (9). Even the
STRONGER positive uniform estimate

 U_b=||a_trial||E+tau+f D1+f^2 D2+f^3 D3            (11)

has, by (6)-(8), the all-trial lower bound on its REPORTED value

                    U_b >= f^3 eta Jstar/M.         (12)

Here D1,D2 are any nonnegative certified residual allowances. This is not
a lower bound on (10)'s signed norm, nor on C_b,3. The only asserted lower
bound is on the particular positive expression (11) using (5).

Use complete P32's unchanged p_aff=(1.1e-10)h3, with
h3=(1+q)^3/|Ch29(-1-2/q)|. The physical bolt allowance reported by (11)
is at least p_aff f^3 eta Jstar/M. Compare THIS to 1.667e-9, granting
the entire bolt margin and zero other debits. No intermediate constant
is compared directly with a physical ceiling.

Freeze the decision test: this physical lower bound on the reported
UPPER expression exceeds the entire margin. Compute its ratio, and also
the less sharp q*sigma version. If it exceeds the margin, every trial
using (11),(5) fails this sufficient P32 certificate. This does not
reject P32 with different residual estimates or the actual construction.

## 5. Quantify a sharper alternative against the same obstacle

For a new UNIFORM comparison coefficient eta_new, even the most favorable
version of (11) necessarily requires

 eta_new <= Cstar(1+q)/(f^3 Jstar-Cstar),
 Cstar=1.667e-9/p_aff,                               (13)

provided the denominator is positive. This is only a necessary threshold
for that estimator: it grants all other terms zero and uses a LOWER bound
on ||j||. It is not a sufficient estimate, not the actual eta, and not a
necessary condition for wormhole existence. Evaluate also the weaker
q*sigma threshold and the hypothetical report with eta_new=1e-10.

The specified sharper METHOD is (10): bound the COMPLETE transferred
residual before taking norms. In particular its third quantum contribution
is F^3 R_E(Rop v), not f^3 eta ||v||. For an actual trial the precise
remaining allowance is Cstar-||a_trial||E-tau minus the other certified
transferred residual debits, and must be positive. The all-trial ratio
from (12) quantifies the minimum gain a successful estimate would have
to make over the best possible positive report, even with other costs zero.

The known replacement q -> f^3 is included in this test; it is not enough
if (12) still fails. Further gains must come from actual spectral placement
and coupled cancellation, not a boundary-basis or normalization change.
No such numerical gain is certified here. F is small at high frequency,
but its uniform norm and the existing form comparison alone supply no
additional directional reduction. An isolated cold/contact optimization
cannot be credited without its contribution to (10).

ONE NEXT CALCULATION (PROPOSED): construct a complete actual-domain
adjoint/Riesz trial and enclose the JOINT transferred residual in (10),
including its commonly restored quantum term, against the actual remaining
bolt margin. No such trial or signed quantum inverse application is
claimed executed in this test. This is the decisive calculation for
whether the direct method can evade the verified estimator loss.

## 6. Scope, full coupled budget and decision boundary

The interior 10^4000 C bound and generic 10^74739.22 absorption coefficient
are not used in (12); the new test can fail even after those losses are
avoided. The former generic-expression failure remains unchanged.
For clarity its actual composition was K0=4e4003,
L0_generic=(2K0)^(4/3), M10=K0 L0_generic^(53/4), and
A_generic=2 M10 1e19 f^3. Only AFTER multiplying by the original bolt
conversion d_b>1e-177 and p_aff did its reported allowance exceed the
bolt margin by >10^74546.2101. Those losses do not represent an actual
operator lower bound. The new residual method avoids this bootstrap;
using it again to bound the dual residuals would reintroduce its losses.

Original delta<=.5, the bounded original-Y extension and I_Y remain
unproved. C_def,3 I_Y<=1.0965284986e32, C_g,3<=1.0650980044e31 and
C_X,3<=9.9684408960e33 remain conditional ceilings with their original
remaining-margin hypotheses. Every coupled entry still requires
|c_hat_ij|+e_ij+C_i,3 I_j h2. No independent input or phase entry is zero.
Actual kappa_m<=5e-5 remains untested. If it passes, other entries must
jointly satisfy
 Y_other+4.084e-22 A_o+Y_phase<=5.34235e-12,
 eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase<.0012069.
Retain the exact mass lift and cutoff transitions, original 2133--2174
conversion, scalar quartic, correlated slow pair and infinite matching
quotient. The .135973<.1599164 quantum-fast screen keeps its limited
two-fast admitted-family scope. Nothing transfers to the nonlinear ball.

The original mass test is |kappa_hat_m|+e_m+C_mass,3 I_m h2<=5e-5;
its conditional debits are 2.14765e-12 affine and .0024731 feedback,
not measured debits. The scalar-quartic row, correlated slow pair,
fast remainder and infinite-dimensional boundary/matching column all
retain unknown signed centers, residual errors, input norms and output
norms. None fits a remaining budget merely because the mass target fits.
The specified complete positive certificate fails if its bolt row fails,
even granting every other entry and error zero. That logical test does
not compute those entries or prove that their actual contributions fail.

The first missing global original-Y/strong-output argument is still the
quantitative weighted end-kernel estimate and scalar/York interface and
infinite matching-quotient returns. The accepted 1ee884b interior theorem
does not cover it. The sharper bolt calculation instead first needs the
full correlated E* residual or its transferred norm. Both gaps remain
explicit; no smooth symbol theorem supplies either on the actual seam.

DECISION CONTRACT: if (12) exceeds the bolt budget, stop estimator (11)
with the available uniform comparison. Use the signed transferred-residual
method (10), or supply a new comparison meeting at least (13) before trial
refinement under a uniform bound. Preserve the construction and P32 for
now; this is an estimator failure, not nonexistence or an actual output
lower bound. No original controlled-solution requirement is declared closed.

Verification contract: freeze this ordinary proof/checker before execution,
validate the retained source/energy receipt against its unchanged ordinary
argument, compute with outward 512-bit arithmetic, preserve all failures,
review and run unchanged acceptance. Supporting checks do not replace
the actual-domain proof behind (7), or calculate the quantum inverse.

Remaining to-do list: actual complete transferred residual and original
output/source maps; every coupled entry; nonlinear admission/contraction;
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability; opening, transport, preparation, separate creation
and the full stable semiclassical wormhole objective.
