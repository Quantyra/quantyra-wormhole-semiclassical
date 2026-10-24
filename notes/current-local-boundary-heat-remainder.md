# A finite-time bound for the local paired boundary parametrix

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous turn: progress. Full stable-wormhole objective active/incomplete;
closed feasibility audit preserved. This calculates a complete capped
boundary-parametrix contribution on 0<t<=.001, not another heat coefficient.
It does NOT yet bound the entire actual-domain heat remainder.

## 1. Define precisely the part being estimated

Retain the actual seed Q, smooth test v=-z^4, paired boundary projector,
original thermal period and the rough h1 as a smearing. The preceding
heat-remainder-transfer audit proves that the first possible boundary
response has weight seven, but gives no numerical constant.

Let Y be inward proper distance from the seam, z'=sqrt(Q), and
f(Y)=sqrt(beta)/z'(Y). Continue the ACTUAL analytic germ to negative Y;
write f_o for that continuation and f_e(Y)=f_o(|Y|). The exact paired
image trace vanishes, so the local boundary correction is the difference
between these two diagonal spin traces. The even coefficient is only
C2,1. It is never replaced by a smooth double.

For a unit Brownian bridge beta_u, put Z_u=sqrt(2) beta_u,
X_u=Y+sqrt(t)Z_u, S=sup|Z_u|, V=sup|X_u|. Define the local term by
restricting the ALREADY DIFFERENTIATED kernel difference to crossing
paths with min X_u<0 and V<=1/2. No event indicator is differentiated.
Both auxiliary extensions can be capped away from this collar to justify
Feynman--Kac; the local contribution is independent of those caps.
It uses the exact seed on the entire collar, not a polynomial surrogate.

First bound paths inside the collar, then all excursions of explicit
auxiliary caps in section 4a. The complement contains the smooth interior
residual and the actual disk-to-parametrix propagation/gluing error,
including the difference on excursions. Those pieces remain in
R_other=R_actual-R_local. They are NOT discarded or assumed small.
In particular no cylinder estimate is transferred through the bolt.

The complete remaining target is

    D_actual=sup_(0<t<=.001) t^-3 int |R_actual(t,z)| dz
                 <=1777.7385,
    B_heat <=1.2e-8 D_actual.                       (1)

We preassign D_local<=600 as a useful local budget. This leaves at least
1177.7385 for the OTHER pieces if the local test succeeds. The execution
also reports its sharper outward bound and resulting available budget.
A local success alone does not certify (1) or the full signed projection.

## 2. Certify the full collar functions, including their analytic tails

Use Q=b^2+u z^2+c3 z^3+c4 z^4+c5 z^5+c6 z^6, b=sqrt(beta), with all
actual parameters, including epsilon and thermal tau. On the complex
disk |Y|<=3/4 use z''=Q'(z)/2, z(0)=0, z'(0)=b. The current coefficient
box gives b in [.408,.409] and, for |z|<=.4,

    |Q'(z)/2|<.225, |Q''(z)/2|<1.282,
    .409*(3/4)+(.225/2)*(3/4)^2<.4,
    (1.282/2)*(3/4)^2<1.                          (2)

Picard iteration in the analytic sup norm is therefore a contraction
of that ball. Also |z'|>=.408-(3/4)*.225>.239, so f is analytic there,
|f|<1.72 and |f'|<1.62. Set

    V0(Y)=int_0^Y v(z(x))dx,
    d(Y)=dot f(Y)=-v(Y)f(Y)-f'(Y)V0(Y).

The same disk gives |d|<.08. These are explicit Cauchy remainder bounds,
not a claim that a finite Taylor polynomial equals the geometry.

Generate the first 65 coefficients by the exact polynomial ODE and
reciprocal-series identities, with outward Arb arithmetic. At r=1/2,
R=3/4, the tail of the p-th derivative divided by |Y|^(k-p) is bounded by

    M R^-k sum_(n>N) (n)_p (r/R)^(n-k), N=64.      (3)

The sum is the p-th derivative of x^(N+1)/(1-x), times x^(p-k).
It is evaluated exactly as a positive rational number at x=2/3.
Parity and the previously proved vanishing jets give the following
pre-stated guards; every guard must be checked against the full tail:

    .95<=f<=1.06, |v(Y)|<=.04 |Y|^4,
    |f'|<=3e-11 |Y|+.6 |Y|^2,
    |f''|<=3e-11+1.2 |Y|, |f'''|<=2,
    |d|<=.03 |Y|^4, |d'|<=.13 |Y|^3, |d''|<=.4 |Y|^2,
    |f_odd^(p)|<=A_p |Y|^(3-p), (A0,A1,A2,A3)=(.15,.45,.9,1),
    |d_odd^(p)|<=B_p |Y|^(7-p), (B0,B1,B2)=(.003,.02,.11).   (4)

The first execution's absolute third-derivative allowance was about
1.660999, exceeding the proposed 1.5 guard. This is a failed helper
allowance, not a source obstruction. The replay uses 2 and propagates
that change to the shifted second derivative below; the D_local<=600
acceptance target is unchanged. All attempts and snapshots are retained.
Here odd means (f(Y)-f(-Y))/2, and similarly for d. The actual finite
geometry is retained. A failed guard stops acceptance and is preserved.

## 3. The exact path functional and both variations

The scalar Fourier potentials of the paired square are k^2 f^2 +/- k f'.
Their bounded-cap Feynman--Kac formula is covered by
[Simon, theorem 1.1](https://arxiv.org/html/math-ph/9907022), after
rescaling his generator to partial_Y^2. Integrating the continuous
frequency gives the zero-winding spin kernel. Its path function is

    H(F,G,t)=F^-1/2 exp[t G^2/(4F)],
    F=int_0^1 f(X_u)^2 du, G=int_0^1 f'(X_u)du.     (5)

The full discrete thermal kernel is restored in section 5; zero winding
is not substituted for the prescribed state.

At fixed original z, the observation proper coordinate moves by V0(Y).
Consequently define on each path

    delta g(X)=d(X)+V0(Y) f'(X),
    delta g'(X)=d'(X)+V0(Y) f''(X),
    delta F=2 int f delta g, delta G=int delta g'.   (6)

The metric/density variation per ORIGINAL canonical time, in base dY
measure, is sqrt(beta)/(2pi t) times

    M=v(Y)H+H_F delta F+H_G delta G.                (7)

The factor v(Y) includes the density and outside-warping conversion.
The known relation dz=sqrt(Q)dY cancels the outside f(Y) exactly. This
keeps the physical normalization, rather than changing the tolerance.
Applying the selected mixed identity -2t partial_t makes the integrand
sqrt(beta)/(pi t) [M-t partial_t M]. The rough h has not been differentiated.

Write Delta for even minus analytic continuation. Then exactly

    Delta M=v Delta H+H_F,e Delta(delta F)+H_G,e Delta(delta G)
                   +Delta H_F delta F_o+Delta H_G delta G_o. (8)

All terms vanish on paths that never cross the seam. Differentiate (8)
before restricting to V<=1/2. If D=t partial_t, then
D f(X)=(sqrt(t)Z/2)f'(X), and similarly for the other path functions.
On crossing, Y<=sqrt(t)S, V<=Y+sqrt(t)S, and |sqrt(t)Z|/2<=V.
Equations (4)-(6) give these positive polynomial majorants in V:

    |G|<=.31 V, |DF|<=2(1.06)(.31)V^2, |DG|<=.61 V,
    |Delta F|<=2(1.06)(.30)V^3, |Delta G|<=.90 V^2,
    |delta F|<=2(1.06)(.031)V^4, |delta G|<=.132 V^3,
    |Delta delta F|<=2[(1.06)(.0132)+(.30)(.031)]V^7,
    |Delta delta G|<=.0544 V^6,
    |D Delta F|<=2[(1.06)(.90)V^3+(.30)(.31)V^5],
    |D Delta G|<=1.8 V^2,
    |D delta F|<=2[(1.06)(.132)V^4+(.31)(.031)V^6],
    |D delta G|<=.403 V^3,
    |D Delta delta F|<=2[(1.06)(.0544)+(.30)(.132)+(.90)(.031)]V^7
                                      +2(.31)(.0132)V^9,
    |D Delta delta G|<=.236 V^6.                   (9)

For example |V0(Y)|<=.008Y^5, |delta g|<=.031V^4,
|delta g'|<=.132V^3, |delta g''|<=.403V^2. Their even-minus-analytic
differences have bounds .0132V^7, .0544V^6 and .236V^5. Thus (9)
retains the coordinate shift and derivatives of the connection vertex.

The box F in [.9025,1.1236], |G|<=.155, t<=.001 gives

    |H_F|<=.7, |H_G|<=.7t|G|,
    |H_FF|<=1.2, |H_FG|<=1.2t|G|, |H_GG|<=.7t,
    |H_FFF|<=3.5, |H_FFG|<=3.5t|G|,
    |H_FGG|<=1.2t, |H_GGG|<=1.5t^2|G|,
    |H_tF|<=1.2G^2, |H_tG|<=.7|G|,
    |H_tFF|<=3.5G^2, |H_tFG|<=1.2|G|, |H_tGG|<=.7. (10)

The executable differentiates (5), removes each displayed factor and
encloses its explicit remaining expression. No derivative of an upper
bound is used as a derivative estimate.

For any H_I (I empty,F,G), the mean-value theorem and product rule give

    |Delta H_I|<=sum_j |H_Ij| |Delta x_j|,
    |D H_I,e|<=sum_j |H_Ij| |Dx_j,e|+t|H_It|,
    |D Delta H_I|<=sum_j |H_Ij| |D Delta x_j|
       +sum_jk |H_Ijk| |Delta x_k| |Dx_j,o|
       +t sum_j |H_Itj| |Delta x_j|, x=(F,G).       (11)

Insert (9)-(11) into (8) and its full time derivative. This produces an
explicit polynomial P(V,t), with NONNEGATIVE coefficients, bounding
|Delta M|+|D Delta M| on every retained path. Every monomial V^p t^q
has p+2q>=7. This is a bound on the whole functional for the entire
time interval; it is not a truncated heat expansion.

## 4. Integrate all bridge paths and the whole observation collar

Reflection gives P(S>x)<=2exp(-x^2), hence
E S^m<=2 Gamma(1+m/2). For a monomial in P, enlarge the retained
crossing event AFTER differentiation and integrate its absolute bound:

    int_0^infty E[V^p 1_cross] dY
      <= [(2^(p+1)-1)/(p+1)] t^((p+1)/2) E S^(p+1). (12)

This follows pathwise by integrating (Y+sqrt(t)S)^p from Y=0 to
sqrt(t)S. In particular the factors for p=7 and p=5 are 1530 and 126.
The enlargement is only for integrating the positive majorant; it does
not extend the validity of (4) outside its collar. Exited paths are
bounded separately below. The in-collar absolute-diagonal term obeys

    D_local,0 <= sqrt(beta)/pi sum_(p,q) P_pq
       [(2^(p+1)-1)/(p+1)] 2 Gamma((p+3)/2)
                            (.001)^(q+(p-7)/2).    (13)

All exponents are nonnegative, so this bounds the supremum down to t=0.
There is no path sampling, spatial sample substitution or selected-mode
approximation. The numerical bound in (13) is to be executed with Arb.

## 4a. Complete the auxiliary boundary parametrix with all excursions

This extends the successful in-collar execution; its source snapshot is
retained. The acceptance target remains 600. On each side of |Y|=1/2,
continue w=f'/f by its quadratic Taylor polynomial in outward distance
r in [0,1/2], multiplied by eta(2r), eta(x)=1-10x^3+15x^4-6x^5.
Keep w=0 afterwards and integrate w to define a positive f. Reflect the
positive cap for f_e; use the actual negative germ for f_o. Both agree
on the entire positive side. They retain their original seam regularity.
This is an explicit auxiliary operator used in a parametrix, NOT a
replacement of the physical disk. Its comparison error remains unknown.

The collar guards imply |w0|<.16, |w1|<.66, |w2|<2.43 at the joins.
On a half-unit transition, |eta'|<=3.75, |eta''|<=24, hence
|w|<.8, |w'|<5, |w''|<36. It follows globally that

    .6<f<1.6, |f'/f|<2, |f''/f|<10, |f'''/f|<100.             (13a)

For v use its quadratic endpoint Taylor polynomial times the same eta,
and zero after the transition. The additional verified endpoint guards
|v|<=.0025, |v'|<.03, |v''|<.3 imply globally
|v|<1, |v'|<10, |v''|<100. The functions have the C3/C2 regularity
needed for the differentiated path identities; higher cap jets are unused.

Put ell=sqrt(t)S. The exact fixed-observation variation is
delta g=-v f-f'[V0(X)-V0(Y)]. Thus the global relative bounds are

    |delta g|/f<=1+2ell, |delta g'|/f<=14+10ell,
    |delta g''|/f<=190+100ell.

For J=f(Y)H, J<=exp(2ell+t). Direct differentiation of (5)-(7), with
D=t partial_t, gives |delta J|<=J A and
|D delta J|<=J(BA+C), where

    A=2+2ell+t(16+14ell), B=ell+t(1+7ell),
    C=10ell+11ell^2+t(16+206ell+145ell^2).

For t<=.001, A+BA+C<=33(1+ell)^2 by coefficient comparison.
There are two cap kernels and 1/f(Y)<1/.6 in the dY trace measure.
An exited crossing path has S>1/(4sqrt(t)). Integrating Y up to
sqrt(t)S, then the Brownian tail, bounds its contribution by

    D_exit <=220 sqrt(beta)/pi exp[t-c/(16t)]
       [25/(64t^4)+95/(16c t^3)+8/(c^2 t^2)], c=7/8.          (13b)

For completeness, write a=1/(4sqrt(t)). On s>=a,
2sqrt(t)s<=s^2/8. Integration by parts against P(S>s)<=2exp(-s^2)
applied to s exp(2sqrt(t)s)(1+sqrt(t)s)^2 gives
2exp(-c a^2)[25a/16+95sqrt(t)/(16c)+8t^(3/2)/c^2].
This yields (13b). Each t^-k exp[t-c/(16t)] increases up to .001,
since c/(16*.001)>4, so the endpoint is a uniform numerical bound.
All auxiliary crossing paths are now covered; no path is sampled.

## 5. Preserve every thermal winding

In the seam-normalized time the ORIGINAL period is
T=2pi sqrt(beta)/gamma>1e33. Poisson summation of the half-integer
thermal spectrum multiplies (5) by

    sum_(m in Z) (-1)^m exp[-m^2 T^2/(4tF)] cos[m T G/(2F)].  (14)

It is the exact sum of both spin potentials, not a zero-temperature
replacement. The complete caps have F<=2.56, giving the exponential
majorant exp[-m^2 x], x=T^2/(11t)>1e67. With
A_m=m^2T^2/(4tF), B_m=mTG/(2F), the relative variation bounds above give
|delta A_m|,|D A_m|<=10(1+ell)A_m,
|D delta A_m|<=100(1+ell)^2 A_m,
|delta B_m|,|D B_m|<=100(1+ell)mT and
|D delta B_m|<=1e4(1+ell)^2 mT. The differentiated prefactor is
bounded by 100(1+ell)^2 exp(2ell+t). The fewer than 16 product terms
are therefore bounded by
1e14(1+ell)^8 exp(2ell+t)(mT)^4 t^-2 exp(-m^2 x).
Integration of Y over a crossing path costs sqrt(t)S. The elementary
bound E[S(1+S)^8 exp(S/15)]<1e10 covers all t<=.001; it follows by
integrating the Gaussian tail and using s/15<=s^2/2+1/450.
Including both cap kernels and signs proves the conservative contribution

    D_local,thermal <=1e26 T^4 t^-6
                sum_(m>=1)(1+m)^8 exp[-m^2 x].      (15)

For x>=1e67, (1+m)^8 exp[-m^2 x/2]<=1 and the remaining Gaussian
series is <=2exp(-x/2). Also T^4t^-6=11^6x^6/T^8<=11^6x^6.
Using exp(x/2)>=(x/2)^100/100! makes (15) <1e-100 uniformly.
This explicit very loose bound covers both signs and all windings.

## 6. Scope, decision and next calculation

Execute (2)-(4), (9)-(13b) and (15) against the pre-stated D_local<=600.
Report the outward result and conditional remaining allowance in (1).
The ordinary path argument, not the checker, supplies the continuum
bound. Preserve all attempts, including a failed guard or budget test.

This closes the COMPLETE CAPPED boundary-parametrix contribution if it
passes, including its excursions and all thermal windings. The actual-domain
interior and differentiated propagation/comparison errors remain UNKNOWN;
the physical disk's excursions have not acquired that comparison bound.
No full B_heat, signed projection or actual
solution is thereby certified. The selected weak regulator argument is
retained, while C2 action regularity on rough base metrics remains open.
Known oversized reference/gluing allowances still do not prove an actual
source obstruction. No cold/scalar component is reoptimized.

Research decision on a pass: retain the construction and this exact
paired-boundary parametrix as a quantitatively affordable component.
The single next calculation is a certified bound on R_other in the SAME
time-reduced diagonal norm against 1777.7385-D_local, with the regular
bolt domain and all excursions/propagation errors. This calculation is
not executed by the local test. Only after it passes is the heat cap
certified; positive-time signed projection and fixed contacts still follow.

The full sufficient screen |P_hat_rest|+E_rest<1.287e-40 remains open.
Preserve the leading heat allowance <5.4970775e-48, magnetic
<1.284506e-40, complete cold <1.693010e-43 and scalar <2.331e-44.
No total-source improvement or nonlinear admission is claimed.
At radius 1e-8, Y0+1e-8Z0+5e-17Z2<=1e-8 and Z0+1e-8Z2<1 remain open.

Remaining to-do list: bound the other actual-domain heat contributions
and complete signed quantum/contact projection; controlled forward/adjoint
correction and nonlinear admission; rough-base action/domain and exact
compatibility; physical matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
