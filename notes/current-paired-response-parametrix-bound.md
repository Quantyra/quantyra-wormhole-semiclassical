# Actual paired response: a two-power comparison and its first conormal commutator

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: NO PROGRESS on quantitative feasibility. The full stable
semiclassical wormhole objective is ACTIVE/INCOMPLETE; the audit stays closed.

OUTCOME: the complete coupled return is UNRESOLVED. This attempt bounds an
actual all-mode response comparison and its FIRST conormal commutator. Their
proper-radial L2 norms are <3.609/K and <7059/K, respectively, with K defined
below. They include the tadpole and finite-domain response. The coefficient
comparison gains two spectral powers in joint mass/frequency scaling.
Keeping the mass-to-frequency ratios additionally bounds the COMPLETE
thermal sums by 1.686 and 2923 at every angular level. The remaining
absolute angular bound diverges. No commonly renormalized strong bound,
kappa_m, signed coupled center or nonlinear solution requirement is certified.

## 1. The unsupported step being addressed

The previous [exact scalar response](current-paired-quantum-scalar-return.md)
has a uniform L2 inverse, but that bound alone does not justify commuting
output derivatives or summing the resulting angular/thermal expressions.
The complete commonly subtracted derivative bound, its actual coupled
inverse, and original phase/bolt rows remain the first unclosed numerical
continuum composition. We now exploit the coefficient equation rather than
assigning the missing derivative constants or differentiating an envelope.

Keep the actual seed, original paired realization, common regulator and
all local allocations. No finite angular/thermal set is selected. For
every n>=1 and j>=0 put

    M=sqrt(n(n+q)/C), d_n=2(q+2n), p=gamma(j+1/2),
    f=sqrt(Q), W=p/f, k=sqrt(M^2+W^2),
    K=sqrt(M^2+p^2/.17), sigma=sqrt(1.05), d_M=M-sigma/2.

Both signs of thermal frequency have the same scalar response. Use x
increasing from the paired seam toward the actual bolt, dz/dx=f. It is
the reverse of rho in the preceding note, with the same positive measure.
For its actual U,V define

    t=(1-UV)/(1+UV), a=(U+V)/(1+UV), b=(U-V)/(1+UV).

Here t,a are NOT physical scalar/phase amplitudes. The exact identities are

    t'=-2Wb, a'=2Mb, b'=2Ma-2Wt,
    t^2+a^2-b^2=1, b(0)=0,
    L=-partial_x(P partial_x)+V,
    P=1/(Wa), V=4W/a+4Mt/a^2,
    H_M=-Mt-4M^2 L^-1.                                (1)

Primes now mean x derivatives. H_M is the complete relative mass response
per thermal sign, before the common four-dimensional finite part. The
global L has the derived Neumann seam and regular positive-form bolt
branch, including the lowest thermal mode. It is not an independently
chosen square of the physical Dirac operator.

## 2. Exact coefficient cancellation with actual finite-domain data

Differentiating (1), without a frozen-potential approximation, gives

    A b=2W't, A=-partial_x^2+4k^2,
    A b'=2W''t-12WW'b, (b')'(0)=0.                   (2)

The last condition uses b(0)=W'(0)=0. It does not prescribe b'(0).
The actual values at an interior proof endpoint are retained below; that
endpoint is not a new physical wall or a replacement for the bolt.

Whole-interval Bernstein bounds on the actual degree-six Q give:

| Interval | Lower Q | Upper absolute Q' | Upper absolute Q'' | Further bound |
| --- | ---: | ---: | ---: | --- |
| 0<=z<=.4 | .13 | .2 | .75 | Absolute Q'''<4.1 |
| 0<=z<=.6 | .08 | .3 | .9 | |

Also Q<.17 on the whole actual seed, Q>.16 for z<=.2, M>18,
and Q'(0)=0 exactly. These statements retain the tiny positive temperature
and its actual root. They are not samples or an h=0 limit. In particular

    |W'/W|<=c1, |W''/W|<=c2,
    (c1,c2)=(.3,.55) on [0,.4], (.6,1.1) on [0,.6],
    |W'''/W|<2 on [0,.4].                             (3)

For example W''/W=Q'^2/(2Q)-Q''/2 and
W'''/W=5Q'Q''/(4sqrt(Q))-3Q'^3/(4Q^(3/2))-sqrt(Q)Q'''/2.
The signs of W' do not enter the bounds.

Write w=W/k^2. Direct differentiation gives
|w''/w|<=c2+3c1^2<3 on both proof intervals. Thus A w>=3W.
At an endpoint z=c, the previous GLOBAL Riccati cone implies

    |b(c)|<=B_c=min(1,W(c)/d_M),
    |b'(c)|<6.12 W(c),                                (4)

since a<=2W/d_M and 2+4M/d_M<6.12. These are bounds on the actual
solutions propagated from both bolts, not guessed endpoint data.

For a target interval ending .2 below c in z, let
Delta=.2/sqrt(.17) and r=sqrt(.17/Q_min). Then k>=K and k<=rK
on the proof interval, W(c)/W(x)<=r, and the proper separation is >=Delta.
Decompose b into the solution of (2) with zero values at both proof
endpoints plus its homogeneous actual endpoint return. The maximum
principle, proved by testing the positive part of the difference, gives

    |b_local|<=2c1 w/3,
    |b_return(x)|<=B_c exp[-2K(x_c-x)].                 (5)

For b', use its Neumann condition at zero and actual value at c.
Its local forcing in (2), after (5), is bounded by
(2c2+8c1^2)W. The remaining forcing is bounded by
12c1 W^2 B_c exp[-2K(x_c-x)]. The barrier

    [(6.12 W(c))+8c1 r^2 B_c] cosh(Kx)/cosh(Kx_c)

retains the endpoint value and dominates this last forcing, because
A cosh(Kx)>=3K^2 cosh(Kx), and the cosh ratio is between
exp[-K(x_c-x)]/2 and 2exp[-K(x_c-x)]. It has zero derivative at
the seam. This proves the comparison with the correct mixed endpoint
conditions. No statement about the physical quantum source is replaced
by this auxiliary decomposition.

After division by w, conservative return constants are

    B_tail=r^3 (18^2/d_18) exp(-36Delta),
    D_tail=2r^3[6.12+8c1r^2/d_18]18^2 exp(-18Delta).   (6)

Here d_18=18-sigma/2. The replacement K>=M>=18 is valid because
K^2 exp(-KDelta) and K^2 exp(-2KDelta) are decreasing above 18;
the remaining mass factors are decreasing there as well. Thus (6)
includes every angular level and every actual thermal frequency.

Executed upper expressions give, respectively, B_tail<7.223e-7 and
D_tail<.985575 for c=.4, and B_tail<1.497e-6 and D_tail<2.172463
for c=.6. Including the local terms proves

    |b|<.201 w, |b'|<1.7 w, 0<=z<=.2;
    |b|<.401 w, |b'|<4 w,   0<=z<=.4.                (7)

The buffers are part of the proof of the ACTUAL return. No thermal or
finite-domain contribution is discarded because its endpoint was moved.

## 3. Compare coefficients through one derivative

Put c_aux=Wa+Mt. Equations (1) give the exact cancellations

    c_aux'=W'a,
    c_aux^2=k^2(1+b^2)-(b')^2/4,
    a=W c_aux/k^2+M b'/(2k^2).                        (8)

Define the algebraic comparison coefficients

    a0=W/k, t0=M/k, P0=k/W^2, V0=4k^3/W^2.

These are coefficients, not a flat replacement of the physical geometry.
The small-p denominators below are kept correlated with W and a. From
(7)-(8), on 0<=z<=.2,

    |c_aux/k-1|<.043/k^2,
    |a/a0-1|<.9/k^2,
    |P/P0-1|<.904/k^2,
    |V/V0-1|<1.86/k^2,
    |P0/P-1|<.9/k^2, |V0/V-1|<1.85/k^2.            (9)

For instance the first numerator is bounded by
.201^2+1.7^2/(4*18^2)<.043; adding 1.7/2 gives the second.
The remaining estimates are rational consequences, with their positive
denominators retained. This is the two-power coefficient cancellation.

One derivative requires a further cancellation, not differentiation of
(9). Let b0=W't/(2k^2). Then A(b-b0)=b0'', b0(0)=0.
Using the larger interval in (7), |t'|<.802 and |t''|<8.2406 on
[0,.4]. For R=1/k^2, |R'/R|<=.6 and |R''/R|<=2.1. The product
rule in b0'' and (3) consequently gives

    |b0''| < (1/2)[2+2(.55)(.802)+.3(8.2406)
                       +2(.55+.3(.802))(.6)+.3(2.1)]w
             <3.5w.

For w4=W/k^4, |w4''/w4|<=3c2+14c1^2<3, so A w4>=3w.
At z=.4 the actual difference satisfies
|b-b0|<=W(.4)[1/d_M+.15/M^2]. The same Dirichlet comparison as (5)
therefore gives, for z<=.2,

    |b-b0| <1.21 W/k^4.                               (10)

Its endpoint-return coefficient after division by w4 is at most
r^5[18^4/d_18+.15(18^2)]exp(-36Delta)<.000309; the local
coefficient 1.2 already dominates 3.5/3. All endpoint data remain present.

Since a'/a=2Mb/a and (log a0)'=M^2(W'/W)/k^2,
substitution of b=b0+(b-b0) uses 2Ma-2Wt=b' to yield

    |(log(a/a0))'|
      <=[2(1.21)+.3(1.7)/2]/[(1-.9/324)k^2]<2.69/k^2.

Also (log(c_aux/k))'=(W'/W)[Wa/c_aux-W^2/k^2]. Therefore

    |(log(P/P0))'|<2.69/k^2,
    |(log(V/V0))'|<5.67/k^2.                          (11)

The coincidence multiplication term obeys

    |M(t-t0)|<.9/k, |[M(t-t0)]'|<2.69/k.             (12)

The second inequality follows from
(t-t0)'=-WW'(t-t0)/k^2-2W(b-b0), not from differentiating
the first bound. This preserves the tadpole in the differentiated estimate.
No claim of arbitrarily many unweighted seam derivatives is made.

## 4. A global response comparison on the same endpoint realization

Use the ORIGINAL graph cutoff chi, one through z=.05 and zero by z=.1,
with the retained beta(11,11) polynomial on the logarithmic transition.
It is not the independently retained mass-lift cutoff chi_m. Define

    P_*=P+chi(P0-P), V_*=V+chi(V0-V),
    L_*=-partial_x(P_* partial_x)+V_*,
    H_*= -M[t+chi(t0-t)]-4M^2 L_*^-1.                 (13)

L_* is an AUXILIARY response comparison. It is not asserted to be a
physical Dirac square, determinant Hessian or a newly renormalized action.
The physical H_M in (1) is unchanged. Outside the collar the coefficients
and bolt branch are exactly those of L. At the seam both natural conditions
are y'=0, with the value free. Equation (9) gives equivalent closed forms,
so the same finite-form endpoint realization is used throughout. Actual
finite-domain returns occur in both inverses and in their difference.

Writing R_L=L^-1 and S_L=L_*^-1, (9) gives the GLOBAL form estimate

    |L_*-L| <=delta L, delta=1.85/K^2<1,
    ||R_L||<=c0/K, c0=sqrt(2)/[4(1-sigma/36)],
    ||S_L-R_L||<=c0(1.85)/(1-1.85/324) K^-3
                         <.678 K^-3.                (14)

Here L>=4d_M and L>=4W imply the middle bound, using the GLOBAL
Q<.17 and max(d_M,p/sqrt(.17))>= (1-sigma/36)K/sqrt(2).
The inverse comparison follows by conjugating the form difference by
L^-1/2; no commutation of its coefficients with the inverse is assumed.
Combining the unrounded expression in (14) with (12) proves

    ||H_M-H_*||_(L2(dx)->L2(dx)) <3.609/K.             (15)

This is the complete per-mode relative response, including the tadpole
and all propagation through the actual finite geometry. It applies to
arbitrary L2 inputs; the actual bounded mass/log directions are included.
It is not an inequality for the summed renormalized physical response.

## 5. First actual conormal commutator, with cutoff costs

Put X=chi z partial_z=h(x)partial_x. It agrees with the original dilation
on z<=.05, vanishes before the bolt, and preserves the scalar Neumann
condition since h(0)=0. All commutators below start on the actual regular
core and extend by the displayed form estimates. The cutoff's C10
regularity suffices for this calculation; no analytic cutoff is assumed.

The retained Bernstein derivative reserves are
K1=21/log(2), K2=840/log(2)^2 for |Dchi|, |D^2chi|.
Using f<sqrt(.17), Q>.16 on this collar gives

    |chi'|<250, |chi''|<122000,
    |h|<=.25, |h'|<63.6, |h''|<31100.                (16)

For clarity, |(z/f)'|<1.1 and |(z/f)''|<.4, so the last two bounds
include EVERY cutoff product derivative. They do not silently change
the original graph norm or the transition.

Direct expansion gives

    [X,L]=-partial_x[(hP'-2Ph')partial_x]
                           -Ph'' partial_x+hV'.     (17)

On the support of X, (7),(9),(11) give

    |P'/P|<.704, |V'/V|<1.107, sqrt(P/V)<.502/K.

Thus the bilinear form of Q_X=[X,L], relative to L in both slots, has
norm <996. The first-order term in (17) is retained and bounded by
Cauchy--Schwarz with sqrt(P/V); it is not dropped as symmetric noise.

Let Delta P=chi(P0-P), Delta V=chi(V0-V). Their derivative bounds are

    |Delta P'|/P <(.9|chi'|+3.332)/K^2,
    |Delta V'|/V <(1.85|chi'|+7.752)/K^2.

Together with (9),(16),(17), they imply

    ||L^-1/2 [X,L_*-L] L^-1/2|| <1070/K^2.          (18)

For example a sufficient coefficient is
.25[2.75(250)+11.084]+1.8(63.6)+.0251(31100)<1070.
These are bounded form operators; neither X nor [X,L] is called bounded
on unweighted L2 by itself.

Set E=L^-1/2(L_*-L)L^-1/2 and F=(I+E)^-1. The exact resolvent
commutator identity gives, with B=L^-1/2 Q_X L^-1/2 and
C_X=L^-1/2[X,L_*-L]L^-1/2,

    [X,S_L-R_L]=L^-1/2[B-F(B+C_X)F]L^-1/2.

Hence its norm is at most

    c0[2(996)(1.85)+1070]/(1-1.85/324)^2 K^-3
          <1751 K^-3.                               (19)

The unrounded coefficient is 1750.40495722576... . Retaining (12) and
the derivative of the cutoff tadpole in (13), the complete response has

    ||[X,H_M-H_*]||_(L2(dx)->L2(dx)) <7059/K.          (20)

Its upper expression is 4 times that in (19), plus
.25[250(.9)+2.69], namely 7058.54232890304... . Thus the FIRST conormal
commutator preserves the two-power spectral gain. This is a numerical
contribution to the missing strong estimate, not a claim that its higher
commutators, common subtraction or original graph composition are closed.

## 6. Execute the complete thermal sums without dropping mass ratios

Using only 3.609/K would give a harmonic thermal upper bound and fail.
That failure is avoidable. Equations (8),(14) retain the sharper bounds

    |M(t-t0)|<.043 M^2/k^3+.85 M/k^2,
    ||H_M-H_*|| <2.752 M^2/K^3+.85 M/K^2.            (21)

In the first inequality use |b'|<1.7W/k^2 and W^2<=k^2. The second
uses the UNROUNDED resolvent coefficient in (14); its M^2 numerator is
not replaced by K^2. Similarly the exact derivative in (12) gives

    |[M(t-t0)]'|<.0129 M^2/k^3+2.675 M/k^2,
    ||[X,H_M-H_*]||<7004.311 M^2/K^3+53.794 M/K^2.   (22)

The M^2 coefficient in (22) is four times the unrounded expression in
(19), plus .25[250(.043)+.0129]. The M coefficient is
.25[250(.85)+2.675]. Thus this improvement changes an actual bound,
not its normalization or input/output units.

For any positive decreasing F, monotone integration proves
gamma sum_(j>=0)F(gamma(j+1/2)) <= integral_0^infinity F(p)dp+gamma F(0).
For the two functions in (21)-(22), the complete positive integrals are

    integral M^2/(M^2+p^2/.17)^(3/2) dp=sqrt(.17),
    integral M/(M^2+p^2/.17) dp=pi sqrt(.17)/2.

The actual gamma<1e-33 and M>18 control the discrete correction, including
j=0; no zero-temperature replacement is used. Hence, for every angular n,
the sums converge ABSOLUTELY in operator norm and

    ||gamma sum_(j>=0)(H_M-H_*)|| <1.686,
    ||gamma sum_(j>=0)[X,H_M-H_*]|| <2923.             (23)

Upward expressions using gamma<=1e-33 and M>=18 are
1.68518619776512579... and 2922.79141117585... . Both thermal signs
are supplied by the factor two in the physical action conversion below.
This is an actual infinite-sum bound for the specified relative response
comparison. It is not yet the common four-dimensional finite part.
The operator and its commutator both converge in norm at fixed n.
Closedness of the actual first-order X therefore identifies the limit
commutator on its core and domain. No exchange with the physical common
regulator or with the still divergent angular bound is inferred.

## 7. Original units, remaining angular sum, and the feasibility decision

The original action factor per half is -pi B^2/(G gamma) in dz.
For the pair of thermal signs and angular multiplicity d_n, the response
difference in (15) contributes -2G gamma d_n/(pi B^2 f) times its
proper-coordinate output. This retains the original action sign.
For near input 1_[0,.2]T the retained norm comparison is
||T_near||L2(dx)<35.534 E(u). Since f>.4 on the original outer output,
(15) gives the per-mode original canonical L2 upper expression

    [406 G gamma/(pi B^2)] (d_n/K) E(u).              (24)

The unrounded coefficient is 405.53746312451... . This conversion is
explicit; (15) or (20) is not relabelled as an original H3 or physical
C5 bound. An independent global E-to-proper-L2 input comparison is not
assumed. The previous actual far-input-to-outer estimate remains in use.

Using (23) instead of the coarser (24), the COMPLETE thermal sum costs
<190 G d_n/(pi B^2) E(u) in this original L2 source norm. Its conormal
commutator costs <329000 G d_n/(pi B^2) E(u); this also includes the
density derivative, since |X log f|<.0625. Neither is an original H3
graph estimate. The near input split commutes with X, which vanishes
before z=.2; it does not introduce a metric interface or drop a contact.

Attempting to sum these angular upper bounds gives infinity:
sum_(n=1)^N d_n=2qN+2N(N+1). The remaining bound has no angular decay.
Thus the complete common-subtracted estimate STILL DOES NOT FOLLOW.
This is failure of a specified positive upper-bound estimator, not
divergence of the physical response or a construction obstruction.
The thermal failure of the coarser 1/K estimate has been repaired by
retaining its actual mass ratios; the angular problem has not.

The missing argument is precise: match this auxiliary response comparison
to the retained short-geodesic principal part and ORIGINAL coincidence
contacts at the SAME physical regulator, then bound the differentiated
all-angular/thermal remainder and its actual coupled inverse. A heat
factor assigned directly to L_* would not be that regulator. No contact
is set to zero or newly chosen; no new state, boundary realization or
finite counterterm is introduced. The existing admitted-family weak
identification is retained in its proved scope. Equations (15),(20) do
not justify strong regulator removal or transfer to a nonlinear ball.

Both original gravitational equations and their finite-part rows remain
necessary. The [reference lift and coupled ledger](current-coupled-global-return-attempt.md)
are unchanged: J_m=int_0^z t^3/Q^2,
n_m=-2beta QJ_m/z, b_m=beta z^3/Q+beta Q'J_m. Both reference rows
cancel locally. Its ACTUAL cutoff still contributes
f_cut=2xz chi_m-2F n_m Dchi_m and alpha_cut=n_m Dchi_m, has zero
reference bolt shift and net mass -xB. It is not a global fixed-temperature
zero mode. Scalar quartic, exact correlated slow pair, fast, interior and
bolt slots stay coupled; none of their post-elimination entries is assigned
its allowance. The original scalar quartic row is not killed by D-4.

Actual m and kappa_m remain UNKNOWN. The original-fast conversion
>.00906808164|m| still requires actual |m|<2.20553814952e-9 if both
1e-11 graph ceilings are imposed. IF kappa_m<=5e-5, the debits remain
2.14765e-12 affine, .0024731 graph and .0195239 independent energy.
Including the established far guard, the residual requirements remain

    Y_remaining<=5.34234999994958e-12,
    eta_remaining+(.0195239+4.583e-13)gamma_E<.0012069.

Even with eta_remaining=0, the actual source-family energy-to-graph
constant must be <.06181654279993595; it is unknown. No finite such
comparison exists for unrestricted outer data alone. The conditional
fast screen .028913+1e-11(1.0706e10)=.135973<.1599164 retains only
its two-fast-input admitted-family scope. Independent phase variations
and original nonlinear admission are not certified by (15) or (20).

EXECUTED: exact actual coefficient cancellation, all-mode bounds through
one derivative, same-endpoint response comparison and first conormal
commutator, COMPLETE thermal sums with their actual discrete correction,
original-unit conversion, and failure of the remaining absolute angular
sum. Supporting symbolic/Arb checks verify these ordinary arguments;
they do not execute a commonly renormalized coupled operator application.

RESEARCH DECISION: retain the construction and sharpen this response
parametrix within the existing common-subtracted method. A quantitative
two-power cancellation survives the actual finite-domain return and its
first conormal derivative, and both relative thermal sums are bounded.
That supports continuing this specific method, but the common angular
subtraction and full coupled feasibility milestone are still OPEN.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete commonly
subtracted all-mode remainder, using this parametrix inside the retained
short-geodesic/contact allocation, after composition with the actual coupled
inverse and original quartic/mass/slow/bolt rows. Its higher conormal
commutators and regulator error must be included in the numerical
T_old(H_s^-1 R_s)^4 estimate and the signed all-column budget test.

Remaining to-do list: complete common-subtracted strong coupled return and
signed response; original physical-norm/neighborhood comparisons;
nonlinear admission/contraction and exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; charge/return, opening,
transport, preparation, separate creation and the full stable semiclassical
wormhole objective.
