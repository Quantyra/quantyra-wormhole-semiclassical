# Actual-disk heat propagation and the selected regulator limit

2026-10-23, E003/S169. Ordinary argument before numerical acceptance.
The full stable-wormhole objective is active/incomplete. The feasibility
audit remains closed. This package bounds the actual-disk short-time heat
remainder; it does not evaluate the complete signed projection.

The retained sufficient target is D_other<=1628.7385, after the complete
capped boundary allowance 149. Unknown propagation constants are not
replaced by the old undifferentiated heat bound. The mass/log direction
h remains a smearing, |h|<1.2e-8. Only the smooth v=-z^4 is differentiated.

## 1. Interior transport on the actual metric, including the regular bolt

Write K=-Q''/2. A unit speed geodesic based at z0, with initial angle alpha,
satisfies z''=Q'(z)/2, z(0)=z0, z'(0)=sqrt(Q(z0)) cos(alpha). Its angular
momentum is ell=sqrt(Q(z0)) sin(alpha). The polar Jacobi field satisfies
J''+K(z(r))J=0, J(0)=0, J'(0)=1. These equations are regular at the bolt:
z'(0)=0 and z''(0)=-gamma there. No reciprocal cylindrical warping is
bounded through the bolt and no singular radial extension is selected.

In a radial parallel spin frame the polar metric is dr^2+J^2 d alpha^2.
For one chirality A_alpha=i(1-J_r)/2, A_r=0; the other chirality is its
complex conjugate. The positive Dirac square has radial expression

    L=-partial_rr-(J_r/J)partial_r-J^-2(partial_alpha+A_alpha)^2
       +(J_alpha/J^3)(partial_alpha+A_alpha)+K/2.

The exact heat transport amplitudes obey

    u0=(r/J)^(1/2),
    uj=-u0 r^-j int_0^r s^(j-1) u0^-1 L u_(j-1) ds.

For G=(4 pi t)^-1 exp(-r^2/(4t)), the uncut P3=G sum_(j=0)^3 t^j uj has
ENTIRE defect (partial_t+L)P3=G t^3 L u3. This equality, not an asymptotic
identification of its next coefficient, is the object being bounded.
The transport recursion is the spin version of [Vassilevich, section
4.3](https://arxiv.org/html/hep-th/0306138v3). Numerical constants below
come from the displayed actual polynomial equations, not that reference.

The smooth normal-coordinate amplitude of degree n contains angular
Fourier indices |m|<=n of the same parity as n. Hence the flat polar
Laplacian has coefficient (m^2-n^2) r^(n-2). Its n=0,1 terms vanish
exactly; for n>=2 its coefficient norm is at most n^2. Retain this
regularity cancellation BEFORE taking absolute values. This concerns
the analytic coordinate expansion, not a truncation of physical modes.

## 2. Positive majorants for complete functions and metric variations

At each closed base-z interval let F_k bound the absolute Taylor
coefficients of Q'(z0+x)/2; K_k=(k+1)F_(k+1). Set p to a factored upper
bound on sqrt(Q(z0)), and Q_abs(x)=p^2+sum_(k>=1)2F_(k-1)x^k/k.
Use nonnegative Taylor series Z,H,B,C satisfying

    Z''=F_abs(Z), Z(0)=0, Z'(0)=p,
    H''=K_abs(Z)H+2(zmax+Z)^4 F_abs(Z)
                    +4(zmax+Z)^3[Q_abs(Z)+2p^2],
    H(0)=0, H'(0)=zmax^4 p,
    kdelta=K_abs'(Z)H+2(zmax+Z)^4 K_abs(Z)
                    +8(zmax+Z)^3 F_abs(Z)+12(zmax+Z)^2 Q_abs(Z),
    B''=K_abs(Z)B, B(0)=0, B'(0)=1,
    C''=K_abs(Z)C+kdelta B, C(0)=C'(0)=0.

These follow by differentiating the actual conformal geodesic equation

    z_s''=exp(-2s v)(Q'/2-s v'Q)+2s v'ell_s^2 exp(-4s v),
    ell_s=exp(s v(z0))sqrt(Q(z0))sin(alpha).

Thus delta z''=(Q''/2)delta z-v Q'-v'Q+2v'ell^2 and
delta K=K_z delta z-2vK-(Qv')'. H bounds the actual metric derivative;
it is not obtained by differentiating an inequality for Z.

Put b=B/r-1, db=C/r, U=(1-b)^-1, W0=U^(1/2), A=(b+D b)/2, D=r partial_r.
The following positive coefficient operator bounds L on a smooth radial
spin-frame amplitude with coefficient norm W:

    Lmaj W = flat(W)+b' U W'
      +r^-2(U^2-1)D^2W+r^-2(Db)U^3 DW+2r^-2 A U^2 DW
      +r^-2(DA)U^2 W+r^-2 A^2 U^2 W+r^-2(Db)A U^3 W
      +(K_abs(Z)/2)W,
    flat(W)_n=(n+2)^2 W_(n+2).

Every factor divided by r^2 starts at least at r^2. Define Wj by replacing
L,u0 and u0^-1 in the integral recurrence by Lmaj,W0,W0, respectively.
Retain the exact coefficient integrals 1/(n+j). The same W0 bounds
u0^-1 since the absolute binomial coefficients of (1+b_actual)^(1/2)
are bounded by those of (1-b)^-1/2. Positive dual arithmetic using
independently derived db and kdelta bounds the intrinsic metric derivative.

For the whole analytic remainder, a Picard ball of radius d on |r|<=.75
must satisfy p*.75+(.75)^2 F_abs(d)/2<d and
(.75)^2 K_abs(d)/2<1. B/r is bounded by sinh(sqrt(Kmax)*.75)/(sqrt(Kmax)*.75).
Require bmax<1. A Volterra bound for H gives kdelta_max and then C/r.
Cauchy derivative estimates on radii .74,.73,.72,.71,.70 give explicit
norms for W0,W1,W2,W3,Lmaj W3 and their independently majorized variations.
These potentially large norms are used ONLY for the analytic tails.

If M is the full norm on |r|<=R=.70 and coefficients through N are retained,
then, for r<=rmax=.625,

    tail <= M (r/R)^(N+1)/(1-rmax/R),
    D^p tail <= M (r/R)^(N+1)
                 sum_(k>=0)(N+1+k)^p (rmax/R)^k.

The last sum is an exact rational geometric derivative. Integrate these
tails against the Gaussian, along with the retained coefficients. No
finite Taylor polynomial is substituted for the metric or the kernel.

## 3. Norms that propagation actually requires

For each kernel jet (1, delta, D_t, delta D_t), record pointwise and BOTH
Schur norms. For the uncut transport kernel their powers are t^-1,1;
for its entire defect their powers are t^2,t^3. Bounds are per spinor
operator norm, before rank-two trace and the selected -2D_t factor.

At fixed endpoints the conformal distance change obeys |delta r|<=r.
The normal Jacobi boundary-value equation bounds |delta alpha|<=r for
r<=.625, |grad v|<=.5 and |K|<=2. The spin parallel-transport variation
is bounded by r under the same guards. Fixed-Hilbert density factors
contribute at most 2. These guards require verification on the geodesic
neighborhoods, including the actual smooth continuation near the seam.
Consequently, with u=r^2/(4t), an amplitude W and its intrinsic bound Ws
have fixed-endpoint variation bound

    QW+2uW, QW=Ws+(1+r)D W+(2+r)W.

For a term G t^j W, its D_t bound is (|j-1|+u)W and its mixed bound is
(|j-1|+u)QW+2(|j-1|+1)uW+2u^2W. The entire defect has j=3.
The Gaussian moments are exact:

    integral G r^n u^k d^2x=(4t)^(n/2) Gamma(n/2+k+1),
    sup exp(-u) r^n u^k=(4t)^(n/2) a^a exp(-a), a=n/2+k.

The polar volume ratio uses the COMMON guard
sinh(sqrt(2)*.625)/(sqrt(2)*.625)<1.14, from |K|<=2 on every retained
geodesic. Both Schur directions use this guard even when the integration
base differs from the amplitude base. The first two development runs used
a cell-specific volume factor; that factor did not justify both Schur
directions. Their arithmetic passes are retained, and their Schur tables
are superseded by the common 1.14 guard. Analytic-tail moments are included.
These norms do not yet give an actual-domain remainder without propagation.

## 4. Boundary-compatible construction and propagation requirements

An admissible construction is P=chi_b H_cap,e psi_b+P3 psi_i, with
psi_b+psi_i=1, psi_b=1 for Y<=1 and zero for Y>=2, chi_b=1 for Y<=3
and zero for Y>=4. The cap must agree with the actual positive geometry
and smooth variation through Y=5. The geodesic kernel has a cutoff equal
to one at r<=.5 and zero at r>=.625. Thus P3 psi_i vanishes near the
physical boundary; the other summand satisfies the actual paired
projector and induced square condition. The regular bolt stays an interior
point. P has identity initial value on the original domain.

The complete defect is the transported interior defect, the distance
cutoff commutator and [L,chi_b]H_cap,e psi_b. The last has spatial
separation at least one. Its metric/time derivatives still require actual
cap derivative guards, full thermal periodization and explicit Gaussian
tail estimates. The previous C3/C2 cap bounds alone are not a bound on
every spatially differentiated mixed commutator. These steps must not be
silently inferred from the previous diagonal certificate.

On the diagonal P differs from its local heat polynomial by
psi_b(H_cap,o-P3) plus the paired boundary difference. The latter can use
the 149 budget only after checking the enlarged caps against the same
full-path argument. H_cap,o-P3 must be propagated as well. A kernel on
a smooth auxiliary surface is a parametrix component; it is not the
physical boundary realization and has its own retained comparison error.

If P has pointwise/Schur constants A,B and F has constants a,b with the
powers above, direct splitting at convolution time t/2 gives

    ||P*F||_infinity <=(A b/32+7B a/24)t^3.

At fixed endpoints D_t(P*F)=P*F+(D_t P)*F+P*(D_t F). Apply this identity
also after delta. Thus the first mixed error has SIX terms, not just
P*(delta D_t F). Higher Volterra terms require their own summable bound;
the small row parameter b t0^4/4 is useful, but its numerical derivatives
must also be included. Multiply the resulting mixed pointwise bound by
rank two, by the selected factor two, and by the time-reduced volume
integral dz<=1. This fixes the physical normalization.

The injectivity/cover argument, enlarged cap construction, all cutoff
commutators and differentiated propagation must be complete before a
numerical interior table is called D_other. Until then D_other is UNKNOWN.
This gate is deliberately explicit: checker passes of coefficient algebra
cannot fill it. The signed projection and controlled solution remain open.

For the cap propagation, suppose globally |f'/f|<=2, |f''/f|<=10,
|f'''/f|<=100 and |v|<=1, |v'|<=10, |v''|<=100. For a bridge from Y to Y',
write a=|Y-Y'|, ell=sqrt(t)S, d=a+ell, uY=a^2/(4t),
uT=(T-T')^2/(4tF). After Fourier integration the absolute kernel is

    (4 pi t)^-1 sqrt(f(Y)f(Y')) F^-1/2
      exp(-uY-uT+tG^2/(4F)).

Its phase is +/-G(T-T')/(2F), and MUST be differentiated too. The
fixed-Hilbert density derivative in this representation is
(v(Y)+v(Y'))/2. The independent derivative bounds give positive polynomials

    A=2+2d+t(16+14d)+2uY+(2+4d)uT+18sqrt(t uT)(1+d),
    B=1+ell+t(1+7ell)+uY+(1+2ell)uT+9ell sqrt(t uT),
    C=ell(10+11d)+t[(16+14d)+ell(192+145d)]+2uY
          +uT[(2+4d)+ell(24+30d)]+ell(181+150d)sqrt(t uT).

The four kernel jets are bounded by the displayed absolute kernel times
1,A,B,AB+C. For clarity, DF/F<=2ell, delta F/F<=2+4d,
D(delta F/F)<=ell(20+22d), |G|/sqrt(F)<=2,
|D(G/sqrt(F))|<=7ell, |delta(G/sqrt(F))|<=16+14d and its mixed bound is
ell(136+96d). These give C by differentiating the exact logarithm,
including its phase, not by differentiating A.

For the pointwise norm sqrt(f(Y)f(Y'))/sqrt(F)<=exp(a+2ell).
For a Schur norm, integration of the tangential Gaussian instead leaves
sqrt(f(Y)/f(Y'))<=exp(a). Put x=sqrt(uY), y=sqrt(uT), and ell=sqrt(t)S.
Young's inequality 2sqrt(t)x<=x^2/20+20t, and its S version, reduce every
positive monomial to exact Gaussian and bridge moments. Use
E S^m<=2 Gamma(1+m/2) for m>0 and E1=1. For the tilted bridge moment,
E exp(S^2/20)<=1+2/19, while for m>0 the bound is
2 Gamma(1+m/2)/(19/20)^(m/2+1). The implementation integrates the whole
positive polynomials. A paired image doubles the pointwise cap allowance;
its Schur allowance is the full-line one, since the two half-line image
integrals partition that line. Full thermal periodization still needs a
separate uniform image-tail guard on the enlarged cap.

After the actual-domain first convolution E1=-P*F is established, the
higher Volterra terms obey an elementary positive-jet bound. For m extra
defects the scalar beta factor is 6^(m+1)/(4m+3)!, at most 140^-m.
The time derivative of each added convolution supplies its integration
measure term. Thus use Fbar jets (b,b_delta,b+b_D,b_delta+b_deltaD).
With x=t0^4 sum(Fbar jets)/140<1, the ENTIRE higher-term mixed allowance
is at most sum(E1 jets) x/(1-x). This includes mixed variations of the
propagator; smallness of an undifferentiated residual alone is insufficient.

## 5. The complete product-gradient contact

The retained product subtraction includes its whole local gradient, not
only the a6 contact. For the rank-two H2 reference at unit curvature scale,

    K_*(t)=1/(2 pi t)+(2/pi)int_0^infinity
                         lambda exp(-t lambda^2)/(exp(2 pi lambda)-1)d lambda.

This is the same [Camporesi--Higuchi spectral input](https://arxiv.org/abs/gr-qc/9505009)
already retained in magnetic-product-vacuum.md. The product-gradient
contact is 4t K_*'(t) v. Subtract its retained terms through t^2.
Taylor's inequality |exp(-x)-1+x|<=x^2/2 gives the ENTIRE remainder

    |R_product(t,z)| <= |v(z)| t^3/(120 pi),
    D_product <= 1/(600 pi).

Here int lambda^7/(exp(2 pi lambda)-1)d lambda=1/480 and
int_0^zh |v|dz<=1/5. This is an all-spectral contact bound for every
t>0, not another coefficient substituted for its remainder.

Research decision pending the full propagation test: retain the actual
paired geometric route and do not optimize the controlled cold/scalar
components. No construction obstruction has been established. A positive
majorant exceeding its allocation would reject that allowance, not prove
the signed source or stable wormholes impossible.

## 6. Complete the actual-domain gate

Here are the construction and estimates supplying the hypotheses of the
preceding test. They precede the final accepted execution. Earlier runs
remain development arithmetic, including their explicitly conditional
transfer rows.

For the negative analytic cap, start at Y=-5/8 and continue the outward
logarithmic derivative w=f'/f by its quadratic endpoint Taylor polynomial
times eta(r/L), eta(x)=1-10x^3+15x^4-6x^5, with L=1/2. Set w=0 after
the transition and integrate it. For the positive cap, keep the ACTUAL
geometry through Y=5 and use the same construction with L=2. Elementary
proper-distance integral bounds give .9<z(5)<.995. At that positive join,

    .8<w<1.01, 0<w'<.2, -.2<w''<0.

Retain the sign of w'' in bounding its quadratic polynomial. Its minimum
is positive and its absolute upper bound is 1.41. On the negative join,
the full analytic-tail endpoint bounds give absolute derivative allowances
approximately (.196575,.675682,1.240259). The signed actual endpoint data,
not these rounded displays, define each cap. The metric test is extended
by integrating (v'_out+v''_out r)eta(r/L), and then kept constant.

For the positive test extension, v''>0 at the join. Its lower value is
at least v(5)+v'(5)>-1, since sqrt(Q)<=1-z and
z^4+4z^3(1-z)<=1. Its upper bound is negative. The other test bounds follow
by differentiating the explicit cutoff. The caps are C3 for f and C2
for v, with bounded piecewise fourth/third derivatives respectively.
No distributional contact arises from those next derivatives: the
preceding derivatives are continuous. The final full-cap guards are

    .6<f<1e6,
    |f^(j)/f| <= (2,10,100,10000), j=1,...,4,
    |v^(j)| <= (1,10,100,10000), j=0,...,3.

The actual finite tau is retained. For the positive join only, stable
rational formulae evaluate w',w'' at tau=0 and then enclose the change
by +/-1e-40. This is a certified perturbation enclosure, not a changed
geometry: tau<1e-60, Q>1e-5 there, and differentiating the explicit
quotients in tau bounds their first three log jets by 1e20. The stored
positive-germ checks independently use the actual finite-tau Q.

The local power-series certificates used on the negative germ extend
the actual polynomial ODE analytically to radius .75. They do not assert
that the globally capped metric equals that analytic continuation beyond
the join. The heat transport inside r<.625 agrees with the analytic germ;
its tail bound can therefore use that germ's analytic continuation.

For the chart cover, extend only the negative germ by this cap and retain
the WHOLE actual disk on the positive side, including its regular bolt.
This auxiliary complete surface of revolution has K<=10 and a flat
cylindrical end. On its actual portion, |K|<=2, |grad K|<=4, and |grad v|<.5;
the same bounds hold on all retained geodesics through the negative germ.
On the first five units of proper radius from the bolt, the exact ODE gives

    s=zh-z <= gamma[cosh(sqrt(1.01)rho)-1]/1.01,
    s(5)+1-zh<1e-30, K<-.99,
    A(rho)=sqrt(Q)/gamma >=sinh(sqrt(.99)rho)/sqrt(.99).

In particular A(3.5)>16. Outside that region A stays above 16: Q' is
negative for z>=2u, while the remaining outer strip has Q>.1. The added
negative cap has A>1e20. These are facts about the actual finite gamma,
not a zero-temperature limit.

There is no closed geodesic of length <=1.5. A geodesic confined within
rho<=5 cannot close, because A'>0 makes rho strictly convex for nonzero
angular momentum; a radial geodesic cannot close either. Any other such
short loop stays at rho>3.5, where its angular speed is <=1/16. Conserved
angular momentum fixes its sign, so it cannot complete the nonzero 2pi
winding needed to close. The conjugate radius is at least pi/sqrt(10)>.75.
The cut-locus alternative (see [Abresch--Meyer, lemma 1.8 and its following
argument](https://library.slmath.org/books/Book30/files/abresch.pdf)) then
gives injectivity radius at least .75. The minimum argument applies here:
any hypothetical injectivity radius below .75 is attained in a compact
region, since the flat cylindrical end has huge circumference. This
explicitly rules out short loops; simple connectivity alone would not.
The .125 margin beyond the support radius persists for a sufficiently
small smooth v-family. The regular spin structure at the bolt is retained.

For the geodesic cutoff use eta((r-.5)/.125). Its first three derivative
bounds are 15,384,30720. On its support, |J'/J|<3. The Jacobi variation
formula and the verified bmax<.2, dbmax<100, kdelta<100 give
|delta(J'/J)|<20000 at fixed endpoints. Differentiate the exact commutator

    [L,chi]P3=-(chi''+(J'/J)chi')P3-2chi' partial_r P3.

The distance, angle, density and transport bounds in section 3, the
displayed cutoff derivatives and t<=.001 give the common four-jet bound

    G 1e8 t^-3 sum_(j=0)^3 t^j
             [Wj+D Wj+D^2 Wj+Wj,s+D Wj,s] 1_(r>=.5).

This deliberately loose bound also dominates the metric derivative of
the cutoff in P itself. Its Gaussian tail is integrated from
u>=.5^2/(4t), using the full analytic amplitude tails. The resulting
additional defect constants are below .007726 pointwise and .111212 in
both Schur directions. For the kernel P, these same allowances gain t0^3.

For the distant cap commutator, differentiate the complete Fourier-integrated
kernel, INCLUDING its phase. With E=partial_Y at a fixed endpoint, the
normalized path jet bounds, in the order 1,delta,D,delta D,E,delta E,DE,
delta DE, are

    f/f: (1,1+2d,ell,ell(7+5d),2,16+10d,5ell,ell(100+50d)),
    f'/f: (2,14+10d,5ell,ell(95+50d),10,200+100d,
                                          50ell,ell(5950+5000d)).

They follow from E X=1-u, D X=sqrt(t)Z/2,
delta X=(1-u)V0(Y)+uV0(Y')-V0(X), and E delta X=(1-u)(v(Y)-v(X)),
together with the cap guards. For F, square the first positive jet;
for G use the second and Cauchy--Schwarz against F. Expanding the inverse,
logarithm and exponential in these THREE square-free jet variables is
just the product and chain rules. It bounds all four jets of
[L,chi_b]H_cap psi_b by

    1e12 t^-3 (1+a+ell+uY+uT)^12 times the absolute cap kernel.

The exact positive polynomial check gives maximum coefficient ratio
790300 before the 1e12 guard. The bound includes the varying commutator
coefficients; 100 dominates the cutoff/operator coefficients and their
first metric variations. Thus this check verifies an explicit ordinary
derivative argument rather than supplying a missing continuum assertion.
Spatial separation a>=1 gives pointwise/Schur defect constants below
8.788e-46 and 1.945e-46. The pointwise constant includes both paired
images; the Schur integral partitions the full line. Their time suprema
occur at t0.

The original thermal period obeys T>1e33. With f<=1e6, every nonnearest
image has uT>=m^2 T^2/(16t fmax^2), and T^2/(16t0 fmax^2)>1e55. Absorb
all differentiated polynomial factors in half its Gaussian exponent.
The explicit positive bridge/Gaussian moments above cost less than the
conservative prefactor 1e100 t0^-6. This is checked by expanding
(1+2sqrt(t0)x+sqrt(t0)S+x^2+y^2)^12(1+S), integrating the tilted bridge
moment and taking Gaussian maxima in x,y; a factor 1e14 covers the kernel,
commutator and both image/difference multiplicities. The inequality exp(x/2)>=
(x/2)^100/100! bounds the full remaining image sum by 2.367e-5194.
Schur norms can also periodize directly by partitioning the full line.
No Matsubara modes, angular sectors or cap excursions have been dropped.

These enlarged caps satisfy the previous complete boundary path argument:
the actual germ through |Y|<=.5 is identical, the global relative guards
and f>.6 are retained, and all windings obey the new tail bound. The
allowance 149 therefore holds for psi_b times THIS boundary parametrix.
The smooth cap comparison and actual-disk Volterra propagation described
in section 4 are now justified. Both are included; the auxiliary kernel
does not furnish its own transfer to the physical disk.

For continuity with the preceding checkpoint, also compare the NEW
boundary parametrix to the OLD capped boundary component. Their local
crossing-path integrands agree through |Y|<=.5. Their difference is
supported on the excursions already bounded in the two constructions;
cutting the new boundary term by psi_b changes only that same tail.
Thus the old-to-new comparison costs less than 4.041e-11 in D, including
all windings. Add it explicitly to D_other. The final rounded guard 952
therefore also bounds the ORIGINAL complement, not only a changed split.
The physical source, normalization and screening threshold are unchanged.

The actual paired positive square generates a unique heat semigroup.
The boundary-compatible P has identity initial value, and the stated
integrable defect bounds make its Volterra series converge. It therefore
equals that actual heat kernel by uniqueness. The same positive jet
majorants justify the selected mixed derivatives of this series. All
piecewise cutoff derivatives are bounded weak derivatives, so the
identities also hold by dominated smooth approximation of auxiliary
cutoffs. The PHYSICAL paired domain and its C2,1 reflection are unchanged;
no smooth double is substituted for it.

## 7. Accepted numerical target and claim boundary

The final execution must test the ENTIRE complement, including the
geodesic and distant commutators, all propagation iterations, the smooth
cap comparison, and the complete product-gradient contact:

    D_other<952<1628.7385,
    D_actual<=149+952=1101,
    B_heat<=1.2e-8*1101=1.3212e-5.

Together with U6<=2.3454197e-13 this gives a sufficient entire omitted
heat-cap allowance

    E_cap <= (2.5e-43/1.0666666e-11)
               [2.3454197e-13*.001+1.3212e-5*.001^2/2]
           <1.5483364e-43<2.5e-43.

The acceptance threshold was not relaxed. The main loss is the positive
pointwise/Schur propagation estimate, especially its mixed metric/time
jets, not the exponentially small cutoff or winding errors. No additional
cancellation is needed to close this SHORT-TIME allocation. These are
certified upper allowances, not values of the signed quantum projection.

This closes the quantitative short-time/common-regulator gate for the
SELECTED weak h,v response on the actual paired seed. It does not prove
a C2 action on rough base metrics, the full finite-time/angular/thermal
projection, an inverse on the nonlinear physical domain or a controlled
solution. At solution radius 1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and
Z0+1e-8 Z2<1 remain unproved. Exact compatibility, physical matching,
causal boundary realization, semiclassical validity/noise and all-sector
stability remain unresolved. The full objective is preserved.

Research decision: continue the actual paired geometric method. Its
previously unknown short-time propagation requirement is attainable with
the executed constants; no construction change is indicated by this test.
The SINGLE next calculation is the complete signed positive-time
quantum/contact projection on t>=.001, restoring the fixed local contacts
once and retaining the entire angular and thermal sums. Combine its
certified interval with the executed heat-cap error above and test
|P_hat_rest|+E_rest<1.287e-40. That full signed test is still PROPOSED,
NOT EXECUTED. Keep the original E_rest<=1e-42 allocation, with no more
than 7.5e-43 assigned to the remaining numerical/tail errors; the excess
short-time margin need not be spent. Passing this sufficient screen would
still not close the other requirements for an actual stable solution.

Remaining to-do list: execute the complete signed positive-time
quantum/contact projection with all local contacts and angular/thermal
tails; controlled forward/adjoint correction, rough nonlinear domain,
nonlinear admission and exact compatibility; physical matching, causal
realization, semiclassical validity/noise, all-sector stability,
charge/return, opening, transport, preparation, separate creation and
the full stable-wormhole objective.

Remaining to-do list: complete the actual-domain propagation gate and
the entire signed quantum/contact projection; controlled forward/adjoint
correction, rough nonlinear domain, nonlinear admission and exact
compatibility; physical matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
