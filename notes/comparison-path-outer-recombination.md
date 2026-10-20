# Recombined canonical and Weyl outer response

2026-10-20, E003/S158 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The canonical cubic allocated in S157 combines with the local Weyl
stress and the already fixed outer vacuum allocation into an
integrable source on the actual comparison path. Its first and
second responses are bounded in the squared-log outer norm.
Together with S157, this bounds the massive-canonical plus Weyl
outer density/null contribution. The LLL, other reference terms,
independent radius equation and full nonlinear solve remain separate.

## 1. Exact local stress in proper jets

Keep the action and scheme of
[the anomaly sources](spherical-anomaly-sources.md). The underlying
anomaly framework is [Mazur--Mottola](https://arxiv.org/abs/hep-th/0106151),
rechecked 2026-10-20; the following response estimates are derived
here, not supplied by that paper.

At a real proper point write the product metric as
A^2(dY^2-a(Y)^2 dtau^2), w=log a, phi=log(r/B),
e=1/C, eta_G=GC/B^2<8.5e-33 and Qmag=q^2/C^2<160000.
Primes in this section mean proper Y derivatives. Define

    R=-2(w''+w'^2), U=phi''+w'phi', V=phi'^2,
    aA=11/360, cA=1/20, bA=1/30, dA=bA-2aA/3,
    Lagr=phi[cA(R+2e)^2/3-4aA e R-Qmag/3]
       +bA R U-3bA U^2+[bA R+(2bA-4aA)e]V
       -6dA U V-3dA V^2,
    P=phi[2cA(R+2e)/3-4aA e]+bA(U+V),
    Q=bA R-6bA U-6dA V,
    Z=bA R+(2bA-4aA)e-6dA(U+V).

The dimensionless mixed Euler components are exactly

    M_t=-Lagr/2+RP/2+P''+(Q'phi'+QU)/2,
    M_Y=-Lagr/2+RP/2+w'P'+ZV-Q'phi'/2+QU/2.       (1)

An independent radial variation checks both signs. Varying w in
int exp(w)Lagr gives the Euler density
Lagr+(D+w')(4w'P-Qphi')-2(D^2+2w'D+w''+w'^2)P=-2M_t.
Varying the radial lapse n at n=1 instead gives
Lagr-2RP-2QU-2ZV-(D+w')(2w'P-phi'Q)=-2M_Y.

For example Box P-Hess^t_t P=P'', whereas
Box P-Hess^Y_Y P=w'P'. These formulas include the variation
of Box phi, not only explicit inverse metrics.
The normalized physical hat profiles of the Weyl sector are

    W_rho=(eta_G/pi)M_t,
    W_N=(eta_G/pi)(M_t-M_Y).                       (2)

Indeed r^2 rho2,W=B^2 rho_hat,W=C^2 M_t/(2pi);
the analogous null expression gives the second relation.
All formulas use sphere-integrated sources, not 4D densities.
The full local radius source and boundary current retain the
independent allocation in the cited action.

## 2. The cubic correlation and exact recombination

Retain S157's fixed seed coordinate y, 0<y<=.08, and path
b_lambda=lambda u(y), r_lambda=B+(r0-B)exp(lambda V(y)).
Put z_lambda=B/r_lambda, d=sqrt(x/C), .4<d<.42,

    Y_lambda(y)=int_0^y exp(lambda u(v))dv,
    a_lambda(y)=a0(y)exp(lambda u(y)),
    c_lambda=c0+lambda m, |m|<1e-11.

The reconstruction's zero cubic lapse fixes the proper jets

    w=eY^2/2+c_lambda Y^3+O(Y^4 L),
    phi=-log(dY)-eY^2/6-c_lambda Y^3/4+O(Y^4 L),
    L=1+log(1/y).                                 (3)

These are finite-jet remainder statements, not differentiated
value asymptotics. Their quantitative extension is proved below.

Insert the two displayed polynomial/logarithmic parts in (1).
The negative powers, independent of the fourth-order remainder,
are

    M_t=11/(120Y^4)-11e/(180Y^2)+37c/(120Y)+O(L),
    M_t-M_Y=2c/(5Y)+O(L).                         (4)

The correlated coefficient -c/4 in phi is essential. With a
general cubic f3 instead, the 1/Y coefficients are
(12c+11f3)/30 and (23c+44f3)/30.

Define, at the ACTUAL radius and proper distance,

    Vvac=11eta_G d^4/(120pi z_lambda^4),
    Jmass=2eta_G c_lambda d/(5pi z_lambda),
    T=2eta_G c_lambda/(5pi Y_lambda),
    U_rho=W_rho-Vvac-Jmass, U_N=W_N-Jmass,
    J=Jmass-T.                                    (5)

Vvac is precisely the fixed reference vacuum source evaluated
at the same radius. It is not a new finite subtraction.
In physical units its density is
11lambda_e^2 r_lambda^2/(240pi). Jmass is the Weyl physical
-lambda_e m_geom/(5pi r_lambda) tail in normalized hat units.

Let P_cub be S157's exact cubic Ward pressure, with zero matching
value at the fixed seed cut y_*=.08. The omitted allocation in
that note is now

    Z_rho=P_cub+J+U_rho,
    Z_N=2P_cub+J+U_N.                              (6)

Adding (6) to S157's regular canonical profiles gives EXACTLY
the full massive canonical plus Weyl profiles minus the original
vacuum allocation. Both cubic pressure terms are retained.
The identity does not assign a separate absolute integral to T.

Although c_lambda is affine, T has nonzero second response:

    T''=(2eta_G/(5pi))[
       -2m Y'/Y^2-c_lambda Y''/Y^2
                          +2c_lambda(Y')^2/Y^3],
    Y'=int_0^y exp(lambda u)u,
    Y''=int_0^y exp(lambda u)u^2.                  (7)

Here and below parameter primes are identified explicitly;
they must not be confused with the proper derivatives in (1).

## 3. Uniform complex local remainder, without a moving C5 map

Allow complex |lambda|<=33. All quantities are evaluated at the
same REAL y. Proper derivatives are implemented algebraically as

    D_lambda=exp(-lambda u(y)) partial_y.

Taking four such derivatives uses at most four unknown spatial
jets. It never evaluates the unknown metric at a complex spatial
point and never differentiates a moving C5 proper-coordinate map.
The analytic reference expressions in (3) may be evaluated at
Y_lambda(y); their logarithm has the branch continued from
positive y. The actual low bounds give
|Y_lambda/y-1|<1e-6 and Y_lambda lies inside |Y|<.081.

Here are sufficient quantitative finite-jet guards, with

    N4(F)=max_(0<=j<=4)
        |y^j partial_y^j F|/[j!100^j y^4 L(y)].

The retained common metric estimates give
N4(u-m y^3)<=2e-8 and the same bound for the linear common
radius-log remainder. Since tau=1-z0>.96, the exact expression
log(r_lambda/r0)=log[1+tau(exp(lambda V)-1)] has cubic
-lambda m y^3/2 and remainder N4<1e-5 on this complex disk.
This follows from finite product/exponential/log series: their
zeroth argument is <4e-7, products through four cost at most
five and the analytic inverse-tau normalized jet sum is <1.05.
For clarity the quadratic bound includes the mass cubic.
Writing A(y)=1e-11 y^3+2e-8 y^4 L, a sufficient N4
majorant for the nonlinear part is

    1e6*33^2*[1e-22 y^2/L+4e-19 y^3+4e-16 y^4 L]<1e-7.

The normalized product and fourfold chain expansions cost
less than 2*5^6*1.05^8<1e6: the factor two bounds the
convergent exponential/log tails on the small scalar disk.
The linear part costs 33*2e-8<1e-6. Each term of the
quadratic majorant increases on this collar, so checking
the endpoint .081 gives a uniform bound. The y^2/L term
cannot be discarded merely because the cubic is small.
The analogous
q_lambda=log(Y_lambda/y) has cubic lambda m y^3/4 and
remainder N4<1e-5. Integration of y^4 L divided by y has
norm below .3 here; the remaining exponential/log terms are
quadratic in u. These are coefficient estimates, not assertions
that proper-coordinate parameter directions satisfy S146.

For completeness, the seed bounds used for the remainder in
(3) can be obtained with slack from the retained disk |Y|<=.15:

    |log a0|<.00056,
    |log[z0/(dY)]|<.00101.

The known quadratic/cubic terms are those in (3). Their fourth
and fifth analytic coefficients have modulus <1. One direct
check uses z0'^2=d^2+e z0^2-2A z0^3+Q0 z0^4+E5 z0^5+E6 z0^6,
where |A|<.334, |Q0|<.501, e<1/6000 and |E5|,|E6|<1e-20.
The coefficients follow by the finite power-series recurrence;
w0=log(z0'/d). The analytic tail starts at degree six.
For either seed remainder and j<=4, its N4 bound is at most

    [sum_(k=4)^5 binom(k,j)y^(k-4)
      +B/.15^4 sum_(k=6)^infinity
                          binom(k,j)(y/.15)^(k-4)]/[100^j L],
    B=.00101.                                    (8)

At y=.081 every expression is <.7. Positivity and the lower
bound L>=1+log(1/.081) make this a uniform bound.

Subtracting the reference jets in (3) at Y_lambda cancels the
cubic terms exactly. In the radius log the common cubic
-lambda m/2 and the log-distance cubic +lambda m/4 sum
to the required proper cubic -lambda m/4.
The reference polynomial shifts have degree at least five
and N4<1e-5, from the same small distance and finite product
bounds. Adding these to (8) and the finite-path remainders
still costs less than one.
The remaining terms have N4<1 before proper-jet conversion.
For that conversion use the finite inverse jets of F=Y_lambda:

    G1=1/F1, G2=-F2/F1^3,
    G3=3F2^2/F1^5-F3/F1^4,
    G4=-15F2^3/F1^7+10F2F3/F1^6-F4/F1^5.

In scaled jets, |F1-1|<1e-6,
|y F2|<2e-8, |y^2 F3|<3e-6, |y^3 F4|<.001.
The fourth-order inverse/composition coefficient sum is at
most 52; inverse denominators and weight ratios cost <1.01.
Use the round reserve 100. Consequently the proper remainders
in (3), denoted delta w and delta phi, satisfy

    |Y^j D_lambda^j(delta w,delta phi)|
          <=100 j!100^j |Y|^4 L, j<=4.            (9)

The actual scaled background jets separately obey
|Y^j D_lambda^j(w,phi)|<=j!100^j for j=1,...,4,
and |phi|<2L. For example the raw common scaled phi jets
are bounded by (3,4,8,20), and those of w by (1,1,1,2);
the small F2,F3,F4 bounds above preserve these generous
100-normalized guards. The same is true on the segment to
the reference jets (3).

Formula (1) is a finite polynomial in these jets, with total
derivative order at most four and at most one undifferentiated
phi. A direct count bounds the sum of absolute coefficients
times jet degree by 1000 for each of M_t and M_t-M_Y,
with e<=1 and the magnetic monomial separated. Product
factorials at derivative order four are bounded by 4!.
Thus the jet variation in (9) costs at most
1000*100*24*100^4*2 L^2; the order-zero magnetic variation
costs at most Qmag*100*.081^4 L^2. This is also below
the local-response coefficient reserve in
[the retained local proof](stationary-centered-local-response.md)
and applies to complex jet coefficients by absolute values.
These counts prove, after (9),
the conservative joint estimate

    max(|delta W_rho|,|delta W_N|)
          <4e-10 L^2.                             (10)

The explicit changes in Vvac
and Jmass are also included: their relative change is bounded
by exp(4|delta phi|)-1 and exp(|delta phi|)-1, while the
factor delta phi retains |Y|^4 L. The actual value difference
is <.001, so these exponential factors are harmless.

For the reference jets, expansion of (1) is finite in powers
of Y and log(dY). After (4) is removed, its coefficient
majorant is <1e8(1+Qmag)L^2. The exponential remainder in
Vvac and Jmass after their singular powers is bounded by
|exp(z)-1-z|<=|z|^2 exp(|z|)/2 with
z=-2eY^2/3-cY^3. With |c|<2 and |Y|<.081, the entire
reference regular remainder in normalized units is <1e-15 L^2.
This also directly checks the lower-power vacuum cancellations;
matching only the 1/Y coefficient would not suffice.

Finally |P_cub|+|J|<1e-28 L^2 on the same disk. Indeed
|T|<eta_G/y, the Ward ratio is <1.01 and its logarithmic
coefficient <3y. For J use
d/z-1/Y=[exp(phi+log(dY))-1]/Y and (3),(9).
Thus (6) is holomorphic, with the uniform weighted bound

    max(|Z_rho|,|Z_N|)<1e-9 L^2.                  (11)

The finite weighted jet algebra, exponential remainders and
bounded Ward integral justify holomorphy into this weighted
space, not only pointwise analyticity. The canonical determinant
is absent from (6); it remains in S157's regular source.

## 4. Responses and the recombined outer constraints

For any real 0<=lambda<=1 the parameter circle of radius 32
lies in |lambda_complex|<=33. Cauchy applied to (11) gives

    value: 1e-9 L^2,
    first response: 3.2e-11 L^2,
    second response: 2e-12 L^2.                  (12)

This includes the changing distance, cubic coefficient, radius
and Ward kernel. The fixed seed outer interval has L(y)<=1+log t,
t=r0/B>=40, as established in S157.

With d_g=3e-7, S152's SQUARED-log outer column applies to (12):

    m2=1.003(2e-12)+4.02d_g(3.2e-11)+5.54d_g^2(1e-9),
    a2=1.006(2e-12)+10.08d_g(3.2e-11)+26.5d_g^2(1e-9).

The allocated pair (6) has second-response bounds

    mass primitive: .00016 m2/2 <1.7e-16,
    absolute mass density /BC: .835 m2/2 <8.4e-13,
    absolute lapse density L1: .00003 a2/2 <3.1e-17.

Add these to S157's regular canonical contribution, using
its sixth-log column separately. The FULL massive canonical
plus Weyl minus the fixed vacuum allocation has outer bounds

    mass primitive <1.5e-14,
    absolute mass density /BC <7.4e-11,
    absolute lapse density L1 <1e-15.             (13)

Taylor allowances are half (13) times lambda^2. These are
density-map derivatives into absolute normed spaces. No scalar
absolute-value derivative, boundary Hamiltonian charge, full
radius-equation bound or nonlinear inverse is asserted.

## 5. Predeclared verification and remaining proof

Revalidate S157 and S152 receipts and current input identities.
Derive (1) from the radial form of the retained local stress;
check its finite jet polynomial and coefficient counts, the
singular Laurent coefficients, the independent cubic mismatch,
seed coefficients four/five, analytic tail sums, inverse jets,
complex guards and all response/constraint budgets.
Preserve every run, failed snapshot and raw output.

Use 70 digits. A manufactured radial geometry takes
w0=eY^2/2+c0Y^3, phi0=-log(dY)-eY^2/6-c0Y^3/4,
e=.001, d=.41, c0=.02, Qmag=4. In fixed y use
u=m y^3+v y^4 L with m=.001,v=.002, and
w_lambda(y)=w0(y)+lambda u(y),
phi_lambda(y)=phi0(y)-lambda m y^3/2
                         +lambda*.003 y^4 L.
This larger smooth test direction makes the changing distance
and cubic detectable; it is not the retained physical path.
Set Y_lambda=int exp(lambda u) and c_lambda=c0+lambda m.

At y=.002,.02,.06 and lambda=0,.4,1, evaluate the proper
jets with D_lambda=exp(-lambda u)partial_y. Compare the exact
recombination (6) against direct full cubic-plus-Weyl-minus-vacuum
evaluation, normalized by eta_G/pi, tolerance 1e-35. Use
adaptive and 32/64-node Gauss quadrature for Y and P_cub,
with tolerance 1e-30 in that normalization. For the distance
Gauss rule use v=y s^16 to resolve the logarithmic endpoint,
with panels [0,.5] and [.5,1], each having the stated node count;
for pressure use v=y exp(s). Evaluate Y inside the pressure
also by its independently integrated moments through power ten
in lambda. The exponential tail on the fixture disk |lambda|<=2
is bounded before use by y exp(2 max|u|)(2 max|u|)^11/11!.
Compare second
responses with 64/128-node Cauchy contours of radius .01 and
one-sided/centered real steps .001,.0005, tolerance
1e-6(1+|target|), refinement ratio .35 unless both errors
are below 1e-35. Use independent Taylor coefficient evaluation
for the target, using finite Taylor jet arithmetic and exact
parameter derivatives of the local proper jets. Require that
freezing Y_lambda in T changes
its second response and that violating phi's cubic correlation
leaves a nonzero 1/Y coefficient.

The fixture tests local tensor and parameter algebra, not a
continuum determinant or a self-consistent quantum state.
The weighted ordinary proof, not sampled geometry, establishes
the uniform actual-path estimate.

Remaining to-do list: other local/radius/LLL/reference sectors,
homogeneous allocation and E0; a general admissible family
and usable full inverse; coupled exterior/state, charge, return
and signed-source solve; opening, stability, transport,
preparation and separate creation.
