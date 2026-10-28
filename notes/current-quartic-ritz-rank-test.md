# Full-domain residual test of a one-dimensional correction

2026-10-27 HST, E003/S169. The preceding email review supplied no new
bound toward the active goal. This calculation tests the proposed exact
quartic Ritz map before spending effort on its scalar coefficient solve.
The full stable semiclassical wormhole objective, closed audit, P32 and
original physical norms remain unchanged.

## 1. Complete operator and the estimator being tested

Use the actual paired first-order Dirac realization, complementary adjoint
and both squares, fixed clock/temperature/spin/scale, physical outer mass
cubic and weight-four logarithms, cutoff transitions and regular disk bolt.
The form space E_t is the retained TT form closure with

    e_t = ell_t/N, N=|W_action| beta,
    T_t=ell_t^-1 B, B=b* a^-1 b,
    0 <= T_t <= q_t=.338/(1-epsilon),
    F_t=(I+T_t)^-1 T_t,
    R_n=T_t-(I+T_t)F_n : E_t -> E_t.

These are the scalar/TT Schur forms of BOTH original equations in the
fixed affine canonical chart. They are not a matter-mode truncation.
The actual common-regulated quantum form, exact UV logarithm, all angular
and thermal terms, coincidence contacts and finite-domain returns remain
in ell_t. The accepted transport comparison is epsilon<1.227e-74.

The proposed exact rank-one Ritz map is

    F_1 u = psi B[u,psi]/(ell_t+B)[psi,psi], psi=z^4.

Its denominator is positive. The accepted quartic column enclosure bounds
its massive quantum application, but does not implement its coefficient
functional or bound R_1 on the entire E_t ball. We test the stronger
statement that ANY output of the form psi c(u) has an unavoidable uniform
residual. This includes inaccurate coefficient solves whose output stays
in span{psi}. It does not assume c is already computed.

## 2. Admissible two-dimensional test space

Take xi_1 to be the accepted J=z^4 field and xi_2=psi. With
Q=sum Q_k z^k, on 0<z<L=z_h the first is exactly

    J_1=z^4,
    xi_1=4[z^3 Q+z^4 sum_(k>=2) k Q_k z^(k-1)/(k-1)],
    U xi=xi/2+QJ/z, V xi=-xi-(QJ)'/2.

For xi_2, J_2=0, U xi_2=z^4/2 and V xi_2=-z^4. These solve

    -(Q^2 J'/z^4)'=Q(4xi/z-xi')/z^4,
    J(0)=0, (Q^2 J'/z^4)(L)=0.

The exact root Q(L)=0 supplies the bolt flux; an interval subtraction is
not used to prove it. The first field has nonzero admissible mass cubic,
and the second is a fast field. Both are smooth at the physical disk bolt.
Their membership in the actual form closure follows from the TT-preserving
outer cutoff/compensator construction in stationary-curved-quantum-tt.md
section 11. That construction retains the cutoff transitions and their
returns; imposing zero mass or an extra inner boundary would change it.
The scalar tests s_1=z^4,s_2=z^5 obey the same fixed-source outer and
regular-bolt conditions.

Let Mss and Mst be the complete scalar and mixed weak matrices from
current-finite-local-action-assembly.md, restricted to these tests.
That action proof includes the Weyl anomaly and boundary Hessian limits;
the remaining canonical determinant has no scalar dependence. Since
W_action<0, a/N=Mss/beta and b/N=-Mst/beta. Maximizing
2b[s,u]/N-a[s,s]/N over ONLY span{s_1,s_2} gives the actual variational
lower bound

    B[u,u]/N >= x^T B_low x,
    B_low=Mst^T Mss^-1 Mst/beta, u=x_1 xi_1+x_2 xi_2.       (1)

The full scalar inverse is not replaced by the small matrix: restriction
of the maximization gives a LOWER bound on the full form. No scalar-solve
accuracy, basis completeness or unknown inverse norm is assumed.

## 3. First additional continuum step: complete quantum energy upper bound

A local TT matrix alone cannot upper-bound the full energy in (1). Use
the accepted complete geodesic comparison ell_0, including its a2, Weyl,
LLL and shared vacuum contacts, with |ell-ell_0|<=.0022 ell and
ell_t<= (1+epsilon)ell. Its positive nonlocal part is the full flow form

    P_M[v]=int_0^infinity r[-B_M(r)] <v,g(rP)v>_SM dr,
    B_M=M^4(K0(Mr)^2-K1(Mr)^2)/(2pi^2),
    g(s)=1-s^2/2+s^4/24-cos(s), P=iX.

Here X generates the Liouville-preserving complete geodesic flow on the
actual C2,1 even double, with half-double normalization. It is not the
Laplacian multiplier on a smooth surrogate. The reference-to-actual
comparison already encloses the full thermal state and boundary returns.
Its original finite fourth contact is included exactly once via the known
a2 assembly. No angular renormalization is performed separately.

For the two polynomial base lifts v_i=V xi_i the even extension is C2,1:
its first two derivatives match at the cubic seam, its third weak derivative
is bounded, and no delta occurs before the fourth derivative. The metric
connection is C1,1. Thus Xv and X^2v have weak flow derivatives in L2,
and X^3v is L2. This places v_i in Dom(X^3): for each successive weak
derivative, integrate against a compact test along the measure-preserving
C1 flow; the distributional identity implies
U_h f-f=int_0^h U_s Xf ds in L2. Strong continuity then gives the generator
domain and its derivative. Repeat three times. Smooth disk coordinates at
the bolt remove the polar-coordinate singularity; there is no bolt trace
condition to add. This argument applies to these trials, not every rough
E_t field and not a nonlinear metric neighborhood.

Repeated integration of 0<=1-cos(s)<=s^2/2 gives
0<=g(s)<=s^6/720 for all real s. The positive-real Schwinger integral
[DLMF 10.32.10](https://dlmf.nist.gov/10.32.E10), squared and integrated
by Tonelli for mu>2|nu|, gives

    int_0^infinity x^(mu-1) K_nu(x)^2 dx
      =2^(mu-3) Gamma(mu/2)^2 Gamma(mu/2+nu)
                       Gamma(mu/2-nu)/Gamma(mu).

At mu=8 the nu=1 and nu=0 moments are 384/35 and 288/35.
The spectral theorem therefore yields

    P_M[v] <= ||X^3 v||_SM^2/(525 pi^2 M^4).                 (2)

This is a genuine upper bound on the EXACT logarithmic form on these
inputs. It does not replace its principal symbol by a sixth derivative
operator on the whole domain. The full angular sum satisfies
S2=sum_(n>=1) 2(q+2n) C^2/[n(n+q)]^2 < 4C^2/q: the summand decreases,
its integral from 1 to infinity is 2C^2/(1+q), and adding the first term
is still below 4C^2/q. All angular levels are covered.

With f=sqrt(Q), proper derivative D=-f partial_z, and angle cosine c,
the exact third geodesic derivative is

    X^3 v=-f c[c^2 A_v+D_v],
    A_v=Q v''', D_v=3Q'v''/2+Q''v'/2.

Angular integration gives pi Q[5A_v^2/8+3A_v D_v/2+D_v^2]. Integrate
this polynomial on the ACTUAL interval [0,L], not beyond the root where
Q changes sign. The physical time factor N_halves 2pi/gamma divided by
N=N_halves pi B_radius^2 beta/(G gamma) is 2G/(B_radius^2 beta).
Consequently the complete geodesic nonlocal energy on each trial is <=

    e_geo,i = 2G S2/(525 pi B_radius^2 beta)
       * int_0^L Q[5A_vi^2/8+3A_vi D_vi/2+D_vi^2] dz.       (3)

All primitives use outward intervals for the actual seed and endpoint.
The complete known canonical matrix includes the bounded shared vacuum
coefficient. If h_i is its diagonal in I/W_action units, then

    ell_0[xi_i,xi_i]/N <= E_i=-h_i/beta+e_geo,i.

ell_0 is positive by its accepted relative comparison. Its Gram matrix
on the two-dimensional test space has largest eigenvalue <= its trace.
Thus e_t[u,u]<=L_up |x|^2, where

    L_up=(1+epsilon)(E_1+E_2)/(1-.0022).

Combining this with (1) proves, on this actual two-dimensional subspace,

    <u,T_t u>_Et >= mu ||u||_Et^2,
    mu = det(B_low)/(tr(B_low) L_up)>0.                    (4)

Using det/trace is a conservative bound on the smaller eigenvalue, not
a claim that the full operator has only two eigenvalues.

## 4. Full-domain consequence, including coefficient/application errors

In the two-dimensional subspace choose a unit u perpendicular in E_t to
(I+T_t)psi. One homogeneous condition always leaves a nonzero vector.
Then for EVERY scalar c(u),

    <u,[T_t u-(I+T_t)psi c(u)]>_Et=<u,T_tu>_Et>=mu.

Hence every bounded rank-one map into span{psi} has ||R_n||>=mu.
The same argument works for any fixed one-dimensional output range.
This is a lower bound for the ACTUAL full operator, proved using admissible
test inputs; it is not an arbitrary-source counterexample or a finite
selection of quantum modes. Coefficient roundoff and quantum-entry errors
cannot remove it while the output stays in that range. If an implementation
has additional output outside the range bounded by a_leak in E_t, the
lower bound weakens by at most (1+q_t)a_leak. A useful extra direction is
a new estimator, not improved precision in the old coefficient.

Since R_n=(I+T_t)(F_t-F_n), the true approximation error also satisfies
||F_t-F_n||>=mu/(1+q_t). Compare this with the previous sufficient
1e-30 application scale; that scale is NOT a physical margin.
For the IDEAL exact Ritz map one also has

    ||R_1|| <= q_t(q_t+2)/(2sqrt(1+q_t)).                  (5)

Indeed write e=psi/||psi||, t=<e,T_t e>, and
F_1=e tensor(T_t e)/(1+t). The residual is the oblique projection
[I-(e+T_t e) tensor e/(1+t)]T_t. Its projection norm squared is
(1+2t+||T_t e||^2)/(1+t)^2, and ||T_t e||^2<=q_t t. Maximizing the
result for 0<=t<=q_t gives (5). This does not implement the coefficient.

For the original residual-to-bolt UPPER estimator, retain
Gamma(z)=W0(3z+3z^2+z^3)+W1(2z+z^2), including the unchanged radius
conversion. Its input z=rho/f_t with rho>=||R_n|| uses
||(I+T_t)^-1||<=1. Therefore Gamma(mu/f_t) is a floor on the value
this particular upper estimator can report for a rank-one map. It is
NOT a physical-output lower bound: J0/J1 may cancel on the actual error.
Compare with the entire 1.667e-9 bolt allowance only as a best-case
test; all other debits and the signed center are still unknown.

## 5. One alternative facing the same full-domain obstacle

The already defined geodesic-preconditioned ordered Neumann map has an
infinite-dimensional leading range, so the rank-one argument does not
apply. Retain its actual-domain forms, delta=1e-38 and m=12:

    L_delta=I+T_t-Q_delta, ||Q_delta||<=eta_delta,
    L0=(ell_0+B)/N, S=L0^-1 T_t,
    K=L0^-1(L_delta-L0), ||K||<=r=(e0+eta_delta)/(1-e0),
    e0=(.0022+epsilon)/(1-epsilon),
    F_delta,12=sum_(k=0)^11 (-K)^k S.

The accepted complete cutoff and truncation proof gives

    d_exact= f_t(eta_delta+r^12)/(1-eta_delta),
    ||T_t-(I+T_t)F_delta,12|| <= (1+q_t)d_exact.           (6)

This concerns EXACT applications of these complete operators. To specify
the numerical gap without replacing unknowns by targets, suppose
||S_hat-S||<=a_S and ||K_hat-K||<=a_K in E_t. Put r_hat=r+a_K<1.
For F_hat=sum_(k=0)^11(-K_hat)^k S_hat, telescoping gives

    a_app <= A12(r_hat)a_S + [q_t/(1-e0)]D12(r_hat)a_K,
    A12(x)=sum_(k=0)^11 x^k, D12(x)=sum_(k=1)^11 k x^(k-1).

Every caller/solve/quadrature/quantum/regulator error must be inside these
two OPERATOR errors (or added separately). A one-column certificate does
not establish either. The complete residual ceiling would be
(1+q_t)(d_exact+a_app). Evaluate the common threshold a_S=a_K=tau_app
needed for a sufficient 1e-30 residual screen, and the conditional example
tau_app=5e-32, without claiming either is attained. This screen is stricter
than the earlier 1e-30 map-application example and is not necessary for
wormhole existence. Also propagate into the original Gamma budget.

The exact-reference approximation is defined on E_t by the accepted
bounded-form comparison. A numerical realization with uniform errors a_S
and a_K is NOT supplied by that theorem. This is the next unsupported
implication for this alternative, and the next decisive calculation.

## 6. Executed outcome and method decision

The [first scientific execution](../evidence/current-quartic-ritz-rank-test-acceptance.json)
uses 512-bit Arb arithmetic and encloses every source-parameter interval.
The earlier [failed launch](../evidence/current-quartic-ritz-rank-test.json)
stopped on a legacy receipt-status string before the scientific calculation;
it is preserved. The final execution also checks the nonlinear application
threshold, without treating that threshold as achieved numerical accuracy.

The variational lower matrix and complete energy trace give

    B_low approximately [[ .40475970354086833, -.32324403348521980],
                          [-.32324403348521980,  .259726029048268612]],
    det(B_low) > .00063992533563025,
    det(B_low)/tr(B_low) > .00096303848863213,
    L_up < 5.908391873687.

The full geodesic nonlocal diagonal contributions in these original
energy units are <2.098e-37 and <3.016e-38. They are included, as are
the actual-to-reference .0022 comparison and retained transport error.
The known local/LLL/vacuum energies are approximately 3.628726776988005
and 2.266666634576169. This establishes

    ||T_t-(I+T_t)F_n|| > 1.62995026264e-4
       for every map into a fixed one-dimensional output range,
    ||F_t-F_n|| > 1.21819900048e-4.

For the ideal exact quartic Ritz map, (5) also gives ||R_1||<.341589.
This is a certified full-domain bracket for that exact mathematical map,
not a numerical implementation of its coefficient. Arbitrarily improving
the quartic quantum column cannot remove the missing independent range.
The true map error exceeds the earlier sufficient 1e-30 application
screen by >1.21819900048e26. This rejects the specified rank-one uniform
estimator. It does not reject a directional physical calculation, a larger
or analytic-range correction, the full coupled inverse, or the construction.

The residual-to-bolt UPPER estimator would report at least 4.95066682420e13,
over 2.96980613329e22 times the ENTIRE bolt allowance, even before other
debits. This is an estimator-report floor, not a lower bound on the actual
signed bolt output. No original physical or coupled budget closes.

For the specified geodesic/Neumann alternative, executed bound (6) is
<8.953511e-31 with EXACT applications. The two numerical applications are
UNKNOWN. If both a_S and a_K were <=5e-32, their combined error would be
<6.712269e-32, and the full residual ceiling would be <9.851613e-31.
The remaining slack in that sufficient 1e-30 screen would exceed 1.48387e-32.
For equal a_S=a_K the actual nonlinear estimator threshold is enclosed
between 5.8261192456e-32 and 5.8261192458e-32. This is a required accuracy,
not an established error or physical margin. The corresponding conditional
bolt error expression is <2.234912e-13; its signed center and other debits
remain unknown, so even this hypothetical pass would not close admission.

Decision: CHANGE THE UNIFORM ESTIMATOR from quartic-only Ritz to a method
with the complete geodesic principal range retained. Keep the successful
quartic resolvent column and all accepted results; do not optimize it or
P32 further for this decision. The dominant new limitation is the omitted
geometric range, not quartic application precision. The alternative avoids
this proved rank-one obstacle but has no certified numerical implementation.

One next calculation: a complete E_t-to-E_t application certificate for
S=L0^-1 T_t, including its infinite-dimensional complement and all source,
scalar-solve and quantum/numerical errors, tested against a_S<=5e-32.
Its multiplier in a_app is <1.002210, so this would pay <6.705e-32 of the
1.04649e-31 residual allowance left after (6). The other application K
would still require its own bound; no such reserve is claimed available
in the physical budget. This next calculation is PROPOSED, not executed.

Remaining to-do list: implemented full-domain correction and all application
errors; signed G_t and complete original coupled/admission budgets, including
actual mass and independent boundary-phase columns; nonlinear existence,
compatibility, matching, semiclassical validity/noise and all-sector stability;
the full stable semiclassical wormhole objective.
