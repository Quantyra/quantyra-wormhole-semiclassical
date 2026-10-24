# Positive-time response: retain the contacts and square the propagation defect

2026-10-23, E003/S169. Ordinary argument and protocol before execution.
Previous turn: progress. The full stable-wormhole objective is active and
incomplete; the feasibility audit stays closed. This is a method decision
inside the OPEN full-projection milestone, not its completion.

## 1. What the preceding domain result permits

Keep the actual paired Dirac domain, regular thermal bolt, period, original
scale and physical-half rank-two trace. The selected mass/log field h=b_h1
is a smearing, and v=-z^4 is the smooth metric direction. The preceding
actual-disk result justifies the common-regulator limit of this selected
response and bounds its omitted t<=.001 part by 1.5483364e-43. It does
not supply positive-time values or a twice differentiable action on a
neighborhood of rough base metrics. The latter remains a prerequisite for
nonlinear admission, even if the selected projection eventually passes.

The first remaining quantitative continuum requirement is a bound on the
ACTUAL positive-time kernel and its v variation, with the contacts kept
correlated. The short-time coefficient B_heat=1.3212e-5 has only been
proved for t<=.001. Extending that inequality to all times is unsupported.
No numerical value of the full projection is inferred from it here.

Use the ordinary proofs in current-smeared-heat-projection-domain.md,
current-disk-heat-propagation.md and current-geometric-mass-conversion.md.
The domain is never replaced by a smooth reflected metric, a wall at the
bolt, or finitely many physical angular/thermal modes. Auxiliary kernels
below must map into the original paired square domain. A matrix identity
test only checks algebra, not this continuum hypothesis.

## 2. An exact resolvent representation, including the usually lost constant

Let T(t)=Tr_red[h exp(-t L)], I_hv=int h v dz, and

    S(t)=delta_v T(t)-2 K_*(t) I_hv,
    S_sub(t)=S(t)-s0/t-s1-s2 t.

Here the s_j are the smeared local coefficients of THIS expression,
including the actual paired completion. They are not independent finite
counterterms. The selected product-subtracted mixed heat remainder is

    r_hv^prod(t)=-2t partial_t S_sub(t).                  (1)

The area contact cancels, so s0=0. The constant is generally NONZERO:

    s1=(1/(12pi)) int h(Delta v-2v) dz
       =-(1/(12pi)) int [Q h' v'+2h v] dz.             (2)

Indeed the rank-two smeared a1 density is -h R/(24pi); its density-
inclusive variation is h Delta v/(12pi). K_* has constant 1/(12pi).
The endpoint flux in (2) vanishes on the actual mass/log phase and regular
disk bolt. Although t partial_t kills s1, its resolvent subtraction is
essential. Omitting it leaves a mass-independent term before the angular
sum, and does not implement the declared common prescription. Formula (2)
is precisely the gradient/contact structure of the retained geometric-
mass conversion; it is not an additional LLL or anomaly allocation.

Put lambda_n=M_n^2=n(n+q)/C, d_n=2(q+2n). The already justified small-time
limit gives S_sub(0)=0 and, for t<=.001,

    |S_sub(t)| <= U6 t^2/4+B_heat t^3/6,
    U6<=2.3454197e-13, B_heat<=1.3212e-5.                (3)

The positive masses control infinity, even though local polynomials were
subtracted for all t. Integration by parts at each mass, followed by the
common convergent subtraction, gives the EXACT massive remainder

    P_mass=(2G/(B^2 F_lambda)) sum_n d_n lambda_n
                      int_0^infinity exp(-lambda_n t) S_sub(t) dt. (4)

This is a locally subtracted resolvent trace, not the raw mass determinant.
At a positive split t0 the corresponding expression instead is

    (2G/(B^2 F_lambda)) sum_n d_n [
       -exp(-lambda_n t0) S_sub(t0)
       +lambda_n int_t0^infinity exp(-lambda_n t) S_sub(t) dt]. (5)

The displayed finite-time boundary term cannot be discarded. If (4) is
used, its complete t<=t0 part is already present; do not add the old heat
cap twice. The fixed local action, LLL, a2 restoration, nonmagnetic Weyl,
reference/charge difference and product contacts of J_rest remain required.

## 3. Execute one complete contact, without a spectral cutoff

The retained hyperbolic product density is

    K_*(t)=1/(2pi t)+(2/pi) int_0^infinity
                            x exp(-t x^2)/(exp(2pi x)-1) dx.

This is the rank-two N=2 specialization of Camporesi--Higuchi,
[equations (3.70) and section 4.2](https://arxiv.org/pdf/gr-qc/9505009).
It is used only for the original product-gradient contact, not for the
actual finite disk. The actual thermal kernel in T(t) remains unchanged.

Subtract its coefficients through t, integrate the COMPLETE positive
spectral density, and retain the full angular tower. The contribution of
this product remainder to (4) is exactly -C_prod I_hv, where

    C_prod=8G/(pi B^2 F_lambda) sum_n d_n/M_n^2
           int_0^infinity x^5/[(M_n^2+x^2)(exp(2pi x)-1)] dx >0. (6)

There is no quadrature cut, frequency cut or dropped product contact.
To enclose (6), expand the Stieltjes denominator to a FINITE order with
its signed exact remainder. For every x>=0, the odd truncation is a lower
bound and the even truncation an upper bound. Thus J=7 and J=8 give

    sum_(j=0)^7 (-1)^j I_(5+2j) S_(j+2)
      <= pi B^2 F_lambda C_prod/(8G)
      <= sum_(j=0)^8 (-1)^j I_(5+2j) S_(j+2),        (7)
    I_(2m-1)=|B_(2m)|/(4m),
    S_p=sum_n d_n/M_n^(2p).

This inequality holds also for x>M_n; it is not an assumed convergent
power series. Each S_p is enclosed by its first N=2048 terms PLUS BOTH
integral-tail bounds. With

    f_p(x)=2 C^p(q+2x)/[x(x+q)]^p,
    int_a^infinity f_p=2 C^p/[ (p-1)(a(a+q))^(p-1) ],       (8)

f_p decreases for p>=2, so the lower tail starts at N+1 and the upper
tail at N. The actual finite q is retained. This encloses a complete
sum, not a physical finite-mode approximation. Use the prior actual
constraint bound |I_hv|<6.704e-12; it does not determine its sign.
Only C_prod, its sign, and the resulting entire contact allowance are
evaluated here. The signed actual h1 functional remains unknown.

## 4. A quadratic error certificate on the actual boundary domain

For one mass let A=L+M^2, R=A^-1. Choose any differentiable P with range
in Dom L such that F=I-A P is bounded and has square-integrable kernel
columns, including the metric derivative. For example the Laplace
transform of the preceding boundary-compatible heat parametrix is a
candidate; its full-time column estimates have NOT yet been supplied.
Its local logarithmic diagonal is permitted; F is a smoother defect.
Define

    Q_P=P+P^dagger-P^dagger A P.
    R-Q_P=F^dagger R F.                              (9)

Multiply out (9) first on the operator core, then use bounded extension.
No commutation or symmetry of P is assumed. Boundary terms are absent
because P has the ACTUAL domain, not because they have been ignored.
For the fixed-Hilbert family D_s=exp(-s v/2)D exp(-s v/2),

    delta D=-vD-c(dv)/2,
    |delta <psi,A psi>| <= kappa_M <psi,A psi>,
    kappa_M=2||v||_infinity+||dv||_infinity/(2M).      (10)

This follows from the quadratic form ||D psi||^2+M^2||psi||^2 and
2M||D psi||||psi||<=that form. Scalar v preserves the first-order
domain; its zero normal jet preserves the induced square condition.
The physical v gives ||v||<=1, ||dv||<.5. Hence kappa_M<=2+.25/M.
This uses the positive paired square, not self-adjointness of the original
unpaired D. The rough h is only a bounded smearing throughout (9)-(11).

Let f_M(z)^2 be the physical-half spin trace of F^dagger F on the diagonal,
and g_M the same column norm for delta F. The columns include the entire
original time circle, and the spin trace already includes rank two.
Differentiating (9), using ||R||<=M^-2 and (10), gives

    |delta(R-Q_P)(x,x)|_trace
           <= M^-2 [2 f_M(z) g_M(z)+kappa_M f_M(z)^2].       (11)

In particular the selected normalized remainder error is bounded by

    E_res <= (2G/(B^2 F_lambda)) sum_n d_n
              int |h| [2 f_n g_n+kappa_n f_n^2] dz.       (12)

The error is quadratic in complete-domain defects. This preserves the
outward and return propagations that a one-sided absolute gluing estimate
loses. It is an operator identity and ordinary form bound; numerical
matrix checks do not prove the domain or the column hypotheses.

## 5. Numerical targets and an exact improvement with a measurable effect

A sufficient target profile, to be VERIFIED on the actual paired disk,
is f_n<=a(M1/Mn)^4 and g_n<=b(M1/Mn)^4 uniformly in z. Then

    E_res<=K_res [2ab+kappa_1 a^2],
    K_res=2G(1.2e-8) M1^8 S4/(B^2 F_lambda).          (13)

The missing a,b are UNKNOWN actual operator estimates. The execution
calculates K_res and the threshold 5e-43/K_res. The 5e-43 reserve leaves
2.5e-43 of the existing 7.5e-43 positive-time numerical/tail allocation;
the sufficient total error allocation 1e-42 is not relaxed.

If the available P is too coarse, make the exact residual correction

    P1=P(I+F), F1=I-A P1=F^2.
    f1_n<=r a(M1/Mn)^4,
    g1_n<=(s a+r b)(M1/Mn)^4,                       (14)

provided ||F_n||<=r and ||delta F_n||<=s on the SAME domain. Noncommuting
metric derivatives give delta(F^2)=(delta F)F+F(delta F). No inverse of
I-F is asserted or needed for (14). Use Q_(P1) in (9), so all its extra
cross terms are part of the approximate projection, not dropped errors.

Execute the conditional capacity test at a=r=.01, b=s=.2, both BEFORE
and AFTER (14). These are pre-stated DESIGN TARGETS, not measured bounds.
A failed allowance in this test is only the arithmetic consequence of
those hypothetical guards; it is not even a known upper bound on the
physical response until the guards are proved. Also compute the maximum
permitted a at b/a=1,5,20 without the correction. This measures an actual
change in an error estimator, not a change in normalization.

For comparison calculate the fraction of the full positive angular
Laplace moments of t and t^2 coming from t>.001, using (8) and the exact
integer incomplete gamma formulas. These are WEIGHTS, not values of the
unknown heat response. They test whether extrapolating the short-time
certificate would affect the decision; it is not authorized by the result.

## 5a. Supply the operator estimates from the existing actual-domain parametrix

The first frozen run executed the contact and conditional capacity test.
It did not verify a,b,r,s. The following argument supplies different,
more conservative COMPLETE bounds, so those design guards need not be
assumed. It uses no positive-time extrapolation.

Let Hcal(t) be the boundary-compatible parametrix already certified in
current-disk-heat-propagation.md, and Ecal=(partial_t+L)Hcal. Its complete
kernel tables, including cutoffs, distant commutators and all thermal
images, imply the following rounded guards for 0<t<=t0=.001:

                  pointwise       both Schur directions
    Hcal          .19/t            3
    delta Hcal    .8/t             11
    Ecal          2.3 t^2          34 t^3
    delta Ecal    27 t^2           450 t^3.             (15)

Reconstruct and verify (15) from the frozen complete tables; do not read
the stored uncut table as already containing the cutoff corrections.
Pointwise bounds are fiber operator norms; a column Hilbert--Schmidt norm
therefore has the factor sqrt(2) for the physical-half rank-two trace.

For EVERY M>=M1 define, using only this certified short interval,

    P0,M=int_0^t0 exp(-M^2 t) Hcal(t) dt,
    F0,M=exp(-M^2 t0) Hcal(t0)
                         -int_0^t0 exp(-M^2 t) Ecal(t) dt. (16)

Integration of the exact defect equation gives F0=I-A P0. This proves
the actual-domain range condition: first integrate from epsilon>0,
then use the strong identity initial value and closedness of A as
epsilon decreases to zero. All kernels retain the original period and
paired boundary; at the bolt they are the regular disk kernels.
The same argument applies to the differentiable v family. The upper
time cut in (16) is an AUXILIARY approximation, not a physical heat cutoff.

Put x=34 t0^4/4, dx=450 t0^4/4. The actual heat operator H=exp(-tL)
has operator norm <=1. Its convergent defect expansion gives

    ||Hcal(t0)-H(t0)|| <= e0=3x/(1-x),
    ||delta Hcal(t0)-delta H(t0)||
                           <=e1=11x/(1-x)+3dx/(1-x)^2. (17)

To see the denominators, the time L1 norm of Ecal is <=x, and every
additional Volterra convolution contributes at most x. Differentiating
the convergent series gives the second formula. These are operator
bounds, not an attempt to continue the t^3 diagonal remainder past t0.
Consequently

    ||F0,M|| <= exp(-M1^2 t0)(1+e0)+x <.723.          (18)

For the metric derivative use (10) and the exact form-Duhamel formula.
For u>0,

    ||A^(1/2) exp(-u A)||
       <=exp(-u M^2)[M+(2e u)^(-1/2)].

The integrable endpoint factors yield, with y=M^2 t0,

    ||delta exp(-t0 A)|| <= kappa_M exp(-y)
                    [y+4 sqrt(y/(2e))+pi/(2e)]
                  <= kappa_1(3+pi/2)/e.

Thus ||delta F0,M|| <= kappa_1(3+pi/2)/e+e1+dx <4.
This derivation uses the actual positive square and form derivative;
it is uniform over the complete angular spectrum.

The endpoint columns in (16) are bounded by
sqrt(2*.19*3/t0) exp(-M^2 t0) and
sqrt(2*.8*11/t0) exp(-M^2 t0). The integrated defect columns are bounded
by sqrt(2*2.3*34) Gamma(7/2)/M^7 and
sqrt(2*27*450) Gamma(7/2)/M^7. Since M1^2 t0<2,

    sup_(M>=M1) (M/M1)^4 exp(-M^2 t0)
                              =4exp(-2)/(M1^2 t0)^2,

these give COMPLETE profile guards

    f0,M<=173(M1/M)^4, g0,M<=680(M1/M)^4,
    ||F0,M||<=.723, ||delta F0,M||<=4.                (19)

These are executed inequalities derived from available constants, unlike
the earlier hypothetical .01/.2 design guards. Higher masses improve the
integral term, so bounding it at M1 after the profile conversion is valid.

For an integer m set P_m=P0 sum_(j=0)^(m-1) F0^j. Its exact residual is
F_m=F0^m. Repeated squaring constructs m=64 using six applications of (14).
The full column bounds are

    a_m=173(.723)^(m-1),
    b_m=680(.723)^(m-1)+(m-1)4*173(.723)^(m-2),
    E_res,m<=K_res[2a_m b_m+kappa_1 a_m^2].           (20)

Execute m=16,32,64 against 5e-43. These are complete physical-domain
ERROR allowances for the exact auxiliary operators, not sampled-mode
estimates. At m=64, Q_(P64) therefore supplies an explicit finite
operator polynomial whose subtracted mixed trace differs from the actual
massive response by the certified (20). Its kernels and all their
products still contain the complete angular and thermal spectra.

This closes the previously unknown propagation-error conditions for this
SPECIFIED representation. It does not evaluate the trace of Q_(P64),
provide numerical quadrature errors, or determine the sign of h1. The
remaining numerical task must include all polynomial cross terms and
common contact subtractions before angular summation; replacing it by a
finite matrix without continuum/tail control would lose this conclusion.

## 6. Decision boundary and claim limits

This execution must report separately: the complete contact enclosure
(6)-(8); the exact-domain reduction (9)-(12), proved by the ordinary
argument; the initial CONDITIONAL targets (13)-(14); and the ACTUAL
complete propagation certificate (15)-(20). None is the complete signed
P_rest. The final execution verifies the supplied operator guards, but
does not evaluate the local/finite-geometry projection. Restored local
terms may not be assessed from the small product remainder alone.

Research decision: sharpen the actual paired method with the two-sided
resolvent certificate and its executed m=64 residual correction. Do not
extend the small-time bound or optimize cold/scalar terms. The actual
propagation error can now meet the reserved precision; the value of the
projection and nonlinear feasibility remain UNDECIDED. This removes a
specific numerical obstacle without supplying a controlled solution.
The full-projection milestone remains OPEN; no construction switch is
justified by the present evidence.

The SINGLE next calculation is a validated evaluation of the complete
subtracted Q_(P64) projection, including its actual h1 constraint pairing,
restored local terms, paired boundary, regular bolt and all angular/
thermal tails. Its complete propagation error is already certified by
(20); spend the next calculation on the signed functional and its
quadrature/tail enclosure, not another generic source bound. This TRACE
EVALUATION is PROPOSED, not executed here.

At radius 1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1
remain open. Preserve the sufficient screening target
|P_hat_rest|+E_rest<1.287e-40. A failure of a sufficient estimator would
not be a failure of wormhole existence. Rough-family action regularity,
the nonlinear physical domain/inverse and exact compatibility precede
a nonlinear solution; physical matching/causal boundary realization,
semiclassical validity/noise and all-sector stability remain independent
unresolved requirements. The closed audit and controlled components stay
unchanged. Preserve every attempt, input hash and independent verification.

Remaining to-do list: complete the actual signed quantum/contact projection
and residual certificate; controlled forward/adjoint correction, rough
nonlinear action/domain, nonlinear admission and exact compatibility;
physical matching, causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
