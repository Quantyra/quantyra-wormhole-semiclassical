# Interior Dirac response with divided high metric jets

2026-10-20, E003/S154 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The fixed-chart actual-domain theorem of
[S145](general-interior-dirac-response.md) extends to divided
second-through-fifth metric jets, retaining its numerical response
constants. Small value and gradient bounds remain separate.
This removes a quantitative mismatch with the metric estimates;
it does not certify the completed endpoint or unknown solver in
the full cold family of [S153](cold-response-from-real-endpoints.md).

## 1. Norm and statement

Keep S145's seed normal charts, fixed cap, original actual paired
domain, physical clock, B,A,C, angular tower and common scale.
For a real field u define

    P5(u)=max(
        actual global and capped value/gradient norms,
        sup_chart |partial^alpha u|, |alpha|<=1,
        sup_chart |partial^alpha u|/(j!100^j), 2<=j=|alpha|<=5).

The cap/global low norms have exactly the same conventions as S145.
The high jets are those of the uncapped field in the fixed seed
chart. Take the supremum over the observation cover when using
the theorem uniformly. The actual domain and regularity conditions
are unchanged. This remains real C5 spatial regularity.

Let P5(b)<=.001, U=P5(u), V=P5(w), and use parameters

    phi=sigma0+b+s u+t w,
    |s|<=1/(200U), |t|<=1/(200V).                 (1)

Zero directions are handled by continuity. The parameter field
zeta=b+s u+t w has every low norm <=.011 and every coordinate
j-jet, j=2,...,5, <=.011 j!100^j. It is not claimed to
have raw high jets <=.011.

For the canonical scalar source T_hat with S145's action
normalization, the conclusion is

    |D T_hat[b][u]| <1e-9 U,
    |D^2 T_hat[b][u,w]| <2e-7 U V.                (2)

The nonconstant scalar source is holomorphic and bounded by
K<5e-12 on (1). The fixed flat-vacuum scalar is separated,
not assigned zero value. Lorentzian trace has the opposite
sign and the same absolute bounds.

P5 and raw C5 are equivalent norms at this finite order. The
new result is their quantitative coefficient ball, not a new
regularity class. Even at the seed, the elementary comparison
raw_C5<=1.2e12 P5 would give the much larger mixed allowance
2.88e17 U V from S145. That comparison also does not enlarge
its background ball. Equation (2) is proved directly.

## 2. Real finite-jet calculation

Write m=M exp(phi(0)) and

    h(z)=exp(phi_cap(z)-phi(0))-1.

Since all value and first-jet hypotheses are unchanged, S145's
low estimates still give

    Re m>.98M, |m|<1.02M, |exp(-2phi(0))|<1.03,
    |h|<=.06, |h(z)|<=2r, r=|z|.                 (3)

In the inner disk r<=1/64, consider any unit real ray. Its
j-th derivative of zeta is bounded by .011 j!150^j for
j>=2, because 100sqrt(2)<150. The first derivative is
bounded by .016. Also |zeta(z)-zeta(0)|<=.016/64,
so its exponential is below 1.001.

All coefficient calculations below are finite Taylor
majorants through degree five. They do not analytically
continue a C5 spatial field. Introduce a formal variable X:

    Z(X)=.016 X+.011 sum_(j=2)^5 (150X)^j,
    E(X)=exp(Z(X)), retained through degree five,
    A(X)=1.026+.026 sum_(j=1)^5 (64X)^j.           (4)

The seed analytic bound .026 on its radius-1/32 chart gives
the coefficients .026*64^j for its mass potential at any
point in the inner half chart. The finite product and
exponential derivative rules therefore bound the j-th ray
derivative of h, divided by j!, by [X^j]1.001 A(X)E(X).
For j=2,...,5 these are less than

    361, 45369, 6276299, 910342924.               (5)

Taylor's integral remainder uses the derivative at every
intermediate point on the real segment; hence (5), not only
the derivatives at the center, controls the remainders.

At the center the seed value and first jet vanish. Its
quadratic coefficient is at most 8; the cubic/quartic
coefficients are at most .026*64^3 and .026*64^4.
Multiplying this center polynomial by E gives

    |h_1|<=.02r, |h_2|<=300r^2,
    |h_3|<=5e4r^3, |h_4|<=7e6r^4.               (6)

The four unrounded coefficients are below
.016, 255.501, 43944.833, 6038268.826, with equality
allowed for the first .016 reserve.

Combining the inner remainder bounds with (3),(6) outside
r=1/64 yields the GLOBAL Taylor-division estimates

    |h-T1 h|<=400r^2,
    |h-T2 h|<=5e4r^3,
    |h-T3 h|<=7e6r^4,
    |h-T4 h|<=2e9r^5.                            (7)

For the outer parts the respective coefficients are bounded by

    2.02*64,
    2.02*64^2+300*64,
    2.02*64^3+300*64^2+5e4*64,
    2.02*64^4+300*64^3+5e4*64^2+7e6*64.

These are strictly below (7). No high derivative of the cap
transition is needed for this division. Its low norms are
already in P5. The fourth-order local subtraction remains
exactly S145's general-gradient polynomial.

## 3. Complete Born remainder and actual domain

Use S145's finite product recursion with

    p=(.02,300,5e4,7e6),
    B=(2,400,5e4,7e6,2e9), A0=2.

It gives, for products of N potentials after removal of all
total degrees through four,

    C(1,4)=2e9, C(2,4)=49140000,
    C(3,4)=446820, C(4,4)=3232.3232.             (8)

For N>=5, retain the first five linear factors 2^5 L^5
and .06 for every additional insertion. The same complex
mass kernel moments apply. Their Bessel integral input is
[DLMF 10.32.9](https://dlmf.nist.gov/10.32#E9), inspected
2026-10-20. It bounds |K_nu(mr)| by K_nu(Re(m)r);
the complete loop moment derivation is retained in S145.

With r_m=1.02, the coefficient of M^-3 is

    20 r_m^2 C(1,4)
     +sum_(N=2)^4 28*5!*10^5*5^(N-2) r_m^(N+1) C(N,4)
     +28*5!*10^5*2^5*5^3 r_m^6/(1-5r_m*.06)
       <1.84e16.                                (9)

The scalar density factor makes this less than 1.90e16,
and thus less than S145's retained allowance 1e17.
The Neumann ratio |m| ||S_m|| ||h||<1.02*.06 is unchanged.
There is no finite angular or Born truncation.

The actual-to-cap comparison uses only the low coefficient
norms, actual domain, separated cap rows and commutators.
Those assumptions are unchanged, so its scalar bound is

    7e11 M exp(-M/256).                          (10)

Higher spatial jets have not entered this operator inverse
estimate. The local remainder is what needed the new proof.

## 4. Common contact and response constants

Keep S145's exact fourth-order scalar polynomial

    Q_phi=Delta^2 phi-3(Delta phi)^2
                     -4 grad phi.grad Delta phi
                     +4|grad phi|^2 Delta phi,
    T_a2=8alpha exp(-4phi) Q_phi,
    |alpha|<100C/(960pi).

At the observation center the seed first jet vanishes and
its coordinate derivatives through four are bounded by 1e8.
The new parameter bounds therefore give

    |p_i|<=.011,
    |A|<=2(1e8+220),
    |T_i|<=2(1e8+66000),
    |B4|<=4(1e8+26400000).

Here p,A,T,B4 have S145's meanings. Direct insertion proves
|Q_phi|<2e17 and |T_a2|<1e17 C, with exp(.044)<1.05.
Thus the larger divided jets still fit the previous contact
reserve. The vacuum constant and Polyakov allocation stay
in the same common scheme, including every density term.

Using sum d_n/M_n^3<2C/3 and
sum d_n M_n exp(-M_n/256)<=24C*256^3, the complete
nonconstant normalized source is bounded by

    K=2 eta_G[(2/3)1e17+7e11*24*256^3+1e17]
      <5e-12, eta_G<8.5e-33.                    (11)

All local subtractions and full sums converge uniformly on
the parameter bidisk by (7)-(10). They define a holomorphic
bounded map into the observation supremum space. Cauchy
with radii (1) gives 200K U and 40000K UV, which prove
(2), including continuity and finite-segment Taylor bounds.
No spatial trace gradient or sixth metric derivative follows.

## 5. Consequence for the cold family and remaining scope

S153 needs an unweighted complex bound and separate REAL
cold endpoint estimates. Equation (11) supplies the same
complex bound when its family norm dominates P5 uniformly,
in place of the old raw high-jet norm. Its real endpoint,
cold direction, domain and clock hypotheses remain in force.
Its quarter-power modulus, Ward estimates and cold constraint
constants then follow unchanged. This is a quantitative
extension of S153's sufficient condition 1, not a deletion
of its unproved condition 2.

In particular the endpoint metric estimates written in
moving normal/proper charts do not automatically imply
these FIXED-chart divided-jet bounds. Proving their transfer
and preservation in a common family remains necessary.
The existing weighted outer response and correlated cubic/
Weyl allocation also remain necessary. The full inverse
must control the retained higher derivatives; these scalar
response bounds do not construct it.

## 6. Focused verification protocol

Check (4)-(11) by exact rational coefficient arithmetic and
independent finite exponential/product recursions. Verify
every inner, center, outer-division, product, loop, contact
and angular/Cauchy budget. Check that deriving (2) by the
raw-norm comparison would give the much larger 2.88e17.
Revalidate S145 and S153's current receipt/input identities.

Independently test Taylor's integral remainder on the real
fields

    zeta_n(x)=1e-6+.001x+a_n sin(nx)+|x-.005|^(11/2),
    a_n=min(.002/n,6e9/n^5), n=400,1000,4000,
    sigma0(x)=-log(1-x^2/4).

These are C5 but not C6. Their low norms and divided high
jets are checked on the full uncapped interval; their fifth
raw jets need not be small. Use h=exp(sigma0+zeta-zeta(0))-1
inside the cap and the same point-value normalization for
the capped extension. The cap cutoff is the quintic
smoothstep on 1/64<=|x|<=1/32.

Use 60 digits, 32/64 Gauss nodes per panel, real segments
ending at r=.0001,.005,1/64, and panel lengths at most 1/n.
Split at x=.005 and use a square-distance substitution on
each panel touching that join to remove the fractional endpoint
singularity from the quadrature integrand. Compare direct h-T4 against the independent
fifth-derivative integral, dividing both by 2e9 r^5.
Require discrepancy and node-refinement difference <1e-20.
Use exact differentiated formulas for the nonanalytic term,
with its continuous fifth derivative at the join. Check the
four global Taylor remainders at positive/negative radii
inside, through and outside the cap. These diagnostics test
finite-jet mathematics, not a physical Dirac determinant.

Preserve every run and raw log. Failed runs retain source
snapshots; acceptance retains exact input hashes before
and after. The ordinary argument supplies the continuum
and infinite-sum statements.

Remaining to-do list: establish a usable global metric family,
outer and other source responses, and the full inverse; solve
exterior/state, charge, return and signed source; complete
opening, stability, transport, preparation and separate creation.
