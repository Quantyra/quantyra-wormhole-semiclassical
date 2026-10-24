# Rough paired metric variations: domain proof and endpoint test; uniform bridge open

2026-10-23, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress. The complete selected projection screen is
certified at <2.123e-42 and is not reoptimized here. The closed audit and
the full stable-wormhole objective remain unchanged.

## 1. Transfer methods, not theorems or freely chosen sources

The user selects the Jacobian--Weyl domain/adjoint reasoning first and the
Navier--Stokes explicit correction accounting next. Source inspection uses
the immutable [Jacobian--Weyl manuscript at 0010354](https://github.com/Quantyra/jacobian-weyl-quantum-phase-space/blob/001035470f8ebfa180c840e507796aec560284b8/docs/notes/A001-arxiv.tex)
and [effective-bounds manuscript at 35b102c](https://github.com/Quantyra/navier-stokes-effective-bounds/blob/35b102c2ea7b47aef091f232132f81d69e68348c/manuscript/main.tex).
Their hashes and the sections actually read are retained in the review
receipt. The former distinguishes a symmetric transport core from its
adjoint and extension domains; its withdrawn interior-indicator argument
is particularly relevant to avoiding artificial cutoff boundaries. The
latter makes a local analytic contraction's constants and remainder
inventory explicit. Neither supplies this rough-metric theorem, a quantum
source, a physical matching condition, or a wormhole existence result.
The stress remains the derivative of the specified action and state.

## 2. The same closed first-order realization throughout a rough family

Let D0 be the original paired Dirac operator on the fixed seed Hilbert
space, with its already established complementary adjoint domains. Retain
the regular disk bolt. Let b be a real identical-on-both-copies scalar
with bounded value and first weak derivative. All admitted mass/log
fields have these properties and b(0)=b'(0)=0. Multiplication by
E=exp(-b/2), and by E^-1, preserves both first-order domains: the
distributional product rule gives D0(E psi)=E D0 psi+c(dE)psi, and the
scalar trace commutes with the fixed projector. Thus

    D_b=E D0 E, Dom D_b=Dom D0,
    D_b^*=E D0^* E, Dom D_b^*=Dom D0^* .             (1)

These are closed-operator identities. To check closedness rather than
infer it from formal symmetry, put B=||b||infinity, G=||db||infinity:

    ||D_b psi|| <=e^B(||D0 psi||+G||psi||/2),
    ||D0 psi|| <=e^B||D_b psi||+G||psi||/2.         (2)

The graph norms are equivalent, and the boundary Green pairing has the
same annihilator. The exact partner J commutes with scalar E, hence
D_b^*=-J D_b J and Lbar_b=J L_b J on their COMPOSED domains.
For Dsharp_b=J D_b, (1) gives Dsharp_b^*=-Dsharp_b on the fixed original
domain. No assertion is made that the unpaired Euclidean D_b is
self-adjoint. No deficiency space or extra extension is selected by an
algebraic commutator. The paired physical projector already selects it.

Use the square's closed form with domain Dom D0. Its operator domain is
{psi in Dom D_b: D_b psi in Dom D_b^*}; it must not be frozen by fiat.
At the outer boundary b'=0 preserves K=S=0. At the bolt the original
regular disk domain is retained, not a punctured radial extension.

There is also a quantitative parameter test. For M>0,
A_b=Dsharp_0+M exp(b) has inverse norm <=e^B/M by accretivity and the
same estimate for its adjoint. For a complex bounded perturbation z,

    ||A_b^-1 M exp(b)(exp(z)-1)|| <=e^(2B)(e^r-1),
    ||z||infinity<=r.

At B<=.001 and r=1/8 the right side is <.134. The geometric resolvents
S_b,plus=exp(b/2)A_b^-1 exp(b/2), and the corresponding minus operator,
therefore have a holomorphic continuation of norm <1.32/M on that tube.
Cauchy on one or two parameter disks gives

    ||D S_b[u]|| <10.56 ||u||infinity/M,
    ||D2 S_b[u,w]|| <337.92 ||u||infinity||w||infinity/M. (3)

For R_b=(L_b+M^2)^-1=(S_plus+S_minus)/(2M), divide these bounds by M.
These are actual-domain operator estimates for every angular mass and
every thermal mode. They are NOT trace estimates: summing a bare M^-2
allowance over d_n is divergent. No ultraviolet conclusion follows from
(3) alone. The B bound is an explicit sufficient hypothesis, not an
evaluated norm of the unknown corrected geometry.

## 3. The admitted rough coefficient bounds and an honest distinction of spaces

Use proper distance y on the physical outer half, fixed at the seed,
and the original units in finite-outer-admission-and-response.md. A bounded
admissible family has, for a fixed log weight L=1+log(1/y),

    b=c_b y^3+r_b,
    |partial_y^j r_b| <=K j!100^j y^(4-j)L, j<=5,
    |c_b|<=K, 0<y<=delta=.1.                       (4)

The same linear bounds hold for directions u,w, with their own norms.
For a specified finite higher log weight L^p the following estimates have
the same proof and finite log moments. An unbounded union of log degrees
is not a bounded Banach ball; no uniform constant for that union is claimed.
Keep the existing interior/cold C5 charts and regular even disk fields.
All coefficient bounds, ellipticity margins and atlas constants must be
uniform on the family being used. In particular (4) implies a uniform
C2,1 bound on the exact reflected coefficients, but not a C4 bound.

The old X-ball admissions are NOT automatic consequences of (4). The
closed audit still requires its specified center offsets and embeddings
into the outer, middle and cold coefficient spaces. A C2 result on
bounded admitted families must not be labeled a preserved X ball.

Approximation must preserve the mass cubic. For the leading log use

    y^4 log y -> y^4 log(sqrt(y^2+epsilon^2)).        (5)

Keep c_b y^3 exactly. The physical-side approximants are smooth at zero,
preserve all first-order traces and the zero normal jet, and converge in
W3,infinity on the physical half. They have UNIFORM weighted bounds (4).
Their ordinary fourth derivatives need not be uniformly bounded.
For t=y^2/(y^2+epsilon^2), the derivatives of the regularized logarithm
are y^-j E_j(t), where E1=t and
E_(j+1)=2t(1-t)E_j'-j E_j. Leibniz gives exact finite polynomial guards
through j=5. The checker verifies these guards on the whole t interval.
No smooth even double of the cubic is introduced. General remainders can
be regularized in dyadic physical-side annuli, keeping the cubic and the
same weighted derivative bounds by fixed-scale convolution.

## 4. An exact complete-frequency calculation of the leading boundary terms

Let u=c_u |y|^3+d_u y^4 log|y| and similarly w near the seam, with a
smooth localization equal to one there. The absolute values implement
the EXACT balanced folded trace, not a smooth physical replacement.
Compute in boundary orthonormal time first. The local flat seed is the
principal tangent problem only; it does not replace the global operator.

The density-inclusive identity on the paired square is

    delta_u K(t)=-2t partial_t Tr[u exp(-tL)],

also for W2,infinity directions by (1) and the fixed positive-time form
limit. In flat Fourier variables q=p+k, the traced first square vertex
is -|p+q|^2 w_hat(k). Integrating the ENTIRE two-dimensional loop momentum
and its Duhamel parameter alpha gives the mixed trace multiplier

    H_t(k)=2/(pi t) F(t k^2),
    F(x)=I(x)-x I'(x),
    I(x)=int_0^1 [1+x(1-2alpha)^2/4]
                          exp[-alpha(1-alpha)x] d alpha.

Integration by parts in alpha also gives

    F(x)=1/2+(1/2)int_0^1 (1+a x)exp(-a x)d alpha,
    a=alpha(1-alpha),
    F(x)=1-x^2/120+x^3/840-x^4/10080+O(x^5).       (6)

This contains the second vertex through the exact time-derivative
identity; it is not a bubble-only calculation. Its t and t^2 terms
reproduce the retained a2 and a6 Hessian coefficients. The zero linear
term is the topological a1 cancellation.

The high-frequency tails are u_hat=12c_u/k^4-24pi d_u/|k|^5+lower terms.
They follow by differentiating the full-line distributional Fourier
transform of |y|^s; the cutoff changes only lower tails. Remove the
ordinary local terms through x^3 to isolate the first boundary term.
The needed integral is an ordinary convergent integral, not a finite part:

    J=int_0^infinity [F(s^2)-1+s^4/120-s^6/840] ds/s^8
      =-5 pi^(3/2)/49152.                         (7)

For example integrate the exponential Taylor remainder four times by
parts, then use the beta integral in alpha. All endpoint terms vanish;
this supplies (7), rather than assigning an analytic-continuation value
to a divergent integral. Half of the full-line trace is the physical-half
trace. Thus its first cubic-cubic boundary coefficient is

    -15/(1024 sqrt(pi)) c_u c_w t^(5/2).           (8)

The cubic-log Fourier cross term is -288pi(c_u d_w+d_u c_w)|k|^-9.
The x^4 coefficient in (6) fixes its logarithmic boundary term:

    -(c_u d_w+d_u c_w)/(70pi) t^3 log t.           (9)

The accompanying finite t^3 coefficient is NOT set to zero; it depends
on the rest of the fields and remains in the remainder. Higher finite
log powers give finite polynomials in log t at this or higher powers.
For the original reduced canonical time, multiply (8)-(9) by the actual
boundary length density sqrt(beta). If coefficients are instead expressed
in z, also convert c_y=beta^(3/2)c_z and d_y=beta^2 d_z; the log scale
change affects the finite t^3 term, not (9).

Equations (8)-(9) are EXACT for the tangent problem with the actual rough
jets and all its frequencies. Transferring these coefficients, with a
uniform error, to rough variable base metrics is NOT established here.
The power counting predicts that fixed lower boundary jets cannot change
them. That prediction does not supply a localization remainder theorem.
The original thermal images and curved geometry remain obligations in
section 5; this diagnostic does not replace the physical functional by
the tangent problem. The complete angular weight remains Theta(t)<=4C/t.
Equations (8)-(9) are integrable against dt/t^2. Their nonzero values are
not an action divergence, and no finite selection of modes was used.

There is a stronger exact test than the leading coefficients. For x>=0,

    0<=R(x):=F(x)-1+x^2/120<=min(x^3/840,x^2/120).   (10)

Indeed, with phi(z)=(1+z)exp(-z),
phi(z)-1+z^2/2=int_0^z s(1-exp(-s))ds lies between zero and
min(z^3/3,z^2/2). Integrate this identity over alpha in (6).
This bound covers ALL frequencies, not a truncated Taylor series. For
even, localized physical-half directions in H3, Parseval and Cauchy--
Schwarz give the complete tangent trace after the a0/a2 subtractions:

    |r_flat(u,w;t)|
      <=t^2/(420pi) ||u'''||_(L2,half)||w'''||_(L2,half),
    int_0^t0 |r_flat|dt/t^2
      <=t0/(420pi) ||u'''||_(L2,half)||w'''||_(L2,half). (11)

The derivatives and norms include the localization region. They cannot
be replaced by the collar-only norm without its transition terms. For
the independent mass/log jets on 0<y<delta, that collar norm is exactly

    int_0^delta |(c y^3+d y^4 log y)'''|^2 dy
     =36 delta c^2+delta^2(144 log delta+84)c d
       +delta^3[192(log delta)^2+288 log delta+388/3]d^2. (12)

It is finite, including the cross term; at delta=.1 the two diagonal
entries are 3.6 and about .484. Thus this executed complete-frequency
test does not find a mass/log regularity obstruction. It supplies an
explicit tangent bound, not a bound for the missing curved commutator.

## 5. The first unsupported continuum step remains a uniform trace remainder

The operator-domain result (1)-(3) is proved. It does NOT imply a
renormalized trace estimate. The endpoint computation (6)-(12) is also
proved for its stated tangent diagnostic, without a smooth surrogate for
the mass/log fields or a finite frequency truncation. The variable-metric
common-regulator bridge is STILL OPEN.

The pre-verification draft asserted that a finite-regularity parametrix
would give a uniform remainder and hence a C2 action. Review withdrew that
assertion before execution. Its complete text and hash are retained in
the review receipt. In particular, a Gaussian/first-gradient estimate
does not prove a fifth-jet Taylor remainder for a twice differentiated,
subtracted trace. The earlier draft did not construct the required mixed
boundary kernel, bound all differentiated Volterra terms, or establish
their uniform time derivatives. Its statement that higher terms start at
parabolic weight five was the unsupported step. A checker of beta moments
would not repair it.

To be precise, localize the exact weak mixed response on the ORIGINAL
paired square, retaining the density/second vertex and the local action
contacts. Let r_curved be its subtracted remainder after extracting the
localized tangent contribution with its cutoff terms. A sufficient target
on a bounded admitted coefficient family is

    |r_curved(g;u,w;t)|<=A(K)||u||A||w||A t^(3/2),
    0<t<=t0=.001, A(K)<infinity,                    (13)

uniformly in the physical-side regulators and in rough base metrics,
with positive-time convergence and continuity. It would give

    int_0^t0 |r_curved|dt/t^2
       <=2 sqrt(t0) A(K)||u||A||w||A
       <.063246 A(K)||u||A||w||A.                 (14)

This numerical factor is evaluated; A(K) is UNKNOWN, not a known
excessive upper bound. No regulator-uniform finite value is certified.
The tangent term (11) is already integrable. Complete angular weights,
finite-temperature images, the actual bolt, boundary contacts, the LLL
and finite local restorations must all stay in the same physical action.
They are not dropped from the missing estimate.

The dyadic weak-form method is a credible sharper method for (13): on a
collar y>=sqrt(t), the fifth coefficient jet has size L(y)/y, while the
mass smearing supplies y^3. Its candidate spatial moment is finite:

    int_0^delta y^2 L(y)dy
       =delta^3[(1-log delta)/3+1/9].             (15)

The candidate principal bridge moment has the finite scaling

    t^-1 int_0^t [s(t-s)/t]^(5/2) ds/sqrt(s(t-s))
       =B(3,3)t^(3/2).                           (16)

These exact moments support the method's plausibility; they do not prove
that every matrix, domain and commutator term admits that reduction.
In particular the y<=sqrt(t) layer must use the actual paired weak kernel
and account for the third-jet seam jump. Bounding its fourth derivatives
after smoothing would reintroduce the known nonuniform C4/H4 loss.
The missing estimate has a positive integrability margin if (13) holds,
and the executed whole-frequency tangent test removes the most direct
power-counting objection. No comparable argument yet bounds the curved
remainder or makes the nonlinear constants realistically numerical.

[Grieser's notes](https://web.math.ku.dk/grubb/notes/heat.pdf), sections
2.2--2.4 and 3, describe a smooth construction and leave minimal
regularity as an exercise. [Smith's theorem](https://arxiv.org/pdf/1809.05614)
concerns potentials on a fixed manifold. Neither is a theorem for this
rough changing metric with its paired domain. The user-selected
Jacobian--Weyl work motivates (1); the Navier--Stokes local analytic norms
do not supply (13). A regulator-uniform C2 action is not claimed here.

## 6. Decision and the full solution requirements

The selected projection screen remains certified at <2.123e-42, below
the sufficient 1.287e-40 target. Its original selected common-regulator
scope is unchanged. The complete reduced-source allowance remains
<1.3079691e-40<2.60e-40. Neither number is an operator norm or an actual
correction bound. Further cold, scalar or local-source optimization
cannot change this passing screen and is not selected.

Requirements closed here: the fixed first-order realization and its
adjoint throughout the admitted W1,infinity family, explicit positive-mass
resolvent bounds, and an exact complete-frequency tangent regularity test.
Requirements NOT closed: the variable-base uniform subtracted trace
bound (13), common-regulator metric differentiation throughout a nonlinear
neighborhood, and the map from the actual X ball into the admitted family.
The square's composed operator domain is still allowed to vary.

Keep the simultaneous conditions at rho=1e-8:

    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1,
    Y0<=5e-9, Z0<=.25, Z2<=5e7 (an allocation, not values),
    1.632e-11 Km+1.36e-9 Kc+Yrest<=5e-9,
    Km<=61.2745098, Kc<=1.4705882, Yrest<=2e-9 (allocations),
    G_(1/4)^lo E_D^2<=1/402.

The exact allocated Km threshold is 1e-9/1.632e-11=61.2745098039...;
the exact Kc threshold is 2e-9/1.36e-9=1.47058823529.... The displayed
decimals are conservative. The admission requirements are

    b45+1e-8 E45<=1e-3; b45=1.3e-4 permits E45<=87000,
    b_c+1e-8 E_c<=5e-7, b_o+1e-8 E_o<=1.5e-7,
    b_D+1e-8 E_D<=2.5e-10.

If b_D=0, E_D<=.025; at that endpoint the G threshold is 800/201.
These inverse, center-offset and embedding quantities remain UNKNOWN for
the actual correction. No first-order resolvent norm in (3) fills them.
The independent actual bolt-difference requirement is <1.667e-9.
Changing normalization changes none of these products or obligations.

Research decision: CONTINUE THE CONSTRUCTION WITH THE SHARPER DYADIC
WEAK-FORM METHOD, but do not promote the nonlinear domain bridge or start
a correction certificate that assumes it. The nonuniform smooth-norm
argument is rejected as a proof method; no construction obstruction or
wormhole nonexistence theorem has been established. The stress remains
the derivative of the specified field theory, state and boundary action.

The SINGLE NEXT CALCULATION is the full curved mixed Duhamel commutator
remainder (13), uniformly on the admitted mass/log family, including its
actual boundary contacts and common regulator limit. This calculation is
PROPOSED, NOT EXECUTED. The first-order domain proof, endpoint coefficients,
whole-frequency tangent bound and explicit moments are EXECUTED here.
If (13) succeeds, quantify the X-to-admission return and the actual
forward/adjoint correction; those are subsequent, unexecuted obligations.

Remaining to-do list: uniform curved mass/log action and source domain;
actual controlled correction and weighted return; nonlinear admission
and higher source remainder; exact compatibility; physical matching and
causal realization; semiclassical validity/noise and all-sector stability;
charge/return, opening, transport, preparation, separate creation and the
full stable-wormhole objective.
