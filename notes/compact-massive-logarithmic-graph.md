# A numerical compact-chart logarithmic graph estimate

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Local component
estimate; no global TT inverse or nonlinear solution.**

The actual canonical massive Hessian has a uniform numerical
fourth/log graph estimate on the retained compact normal charts,
including charts at the bolt. It also tolerates the negative
principal fourth-order coefficient supplied by the Weyl action.
The resulting lower-order constant is extremely conservative;
it supplies no usable nonlinear contraction by itself.

## 1. Precise chart, operator and norm

Use the normalized conformal charts of
[the strong interior remainder](strong-interior-quantum-remainder.md).
Choose a center with canonical inward distance Y>=1/4 from the
seam. The chart has gbar=exp(2sigma)|dz|^2, sigma(0)=D sigma(0)=0,
radius 1/32, and |sigma|<a=.027 on the retained smaller disks.
Let V={|z|<1/256}, O={|z|<1/128}. Their images lie in Y>=1/8.
The actual smooth bolt is allowed as the chart center.

For T smooth and supported in V, extend T by zero on both the
actual orbit and the Euclidean coordinate plane as appropriate.
Use Euclidean measure for every norm in the following statement,
the unitary two-dimensional Fourier transform, and

    ||T||_4log = ||(1+|xi|^2)^2 log(2+|xi|) T_hat||_2.

H_mass retains the complete actual paired massive Hessian per
fixed gbar area, without the factor 2G/B^2. For any measurable
real b on O with |b|<=4 put

    A_b T = H_mass T
       - [C/(30pi)] b(z) exp(-4sigma(z)) Delta_0^2 T. (1)

The second term is in nondivergence form. It is an auxiliary
principal correction, not the full Weyl operator. The theorem is

    ||T||_4log <= (4000/C)||A_b T||_L2(O)
                         +1e390 ||T||_L2(R2).        (2)

Taking b=0 proves an estimate for the ACTUAL massive Hessian.
For the physical principal Weyl coefficient b=omega+1/6,
omega=log(r/B), the same estimate applies to (1). The lower
Weyl terms, LLL, classical/reference terms, TT reconstruction,
partition commutators and weighted outer region must still be
included before asserting a full TT graph or inverse estimate.

## 2. Exact flat multiplier and its scale derivative

Keep the full masses M_n, multiplicities d_n, and common finite
coefficient F_eff from
[the curved TT proof](stationary-curved-quantum-tt.md). Write
kappa4=-F_eff/4>13.9. The positive flat multiplier is

    K(p)= C kappa4 p^4/(30pi)
      +sum d M^2/(2pi) integral_0^1 R(alpha(1-alpha)p^2/M^2) d alpha,
    R(z)=log(1+z)-z+z^2/2.                          (3)

The full flat geodesic model equals this multiplier, not merely
its leading symbol. To check the normalization at every order,
for k>=3 the Bessel moment of one mode is

    -integral_0^infinity r^(2k+1) B_M(r) dr
      =2^(2k-2)(k!)^4 M^(2-2k)/
                                 [pi^2 k(2k+1)!].

Multiplication by the cosine Taylor coefficient and
integral cos^(2k)theta dtheta
=2pi(2k)!/[4^k(k!)^2] yields
(-1)^(k+1)M^(2-2k)(k!)^2/[2pi k(2k+1)!].
This is exactly the coefficient obtained from R and
integral_0^1[alpha(1-alpha)]^k d alpha.
The common fourth restoration also agrees. Analyticity in the
strip of the massive Fourier transform extends equality from a
neighborhood of zero to the real axis. The bound by p^6/M^4
on fixed frequency intervals permits the full angular sum.

These moments follow by squaring the positive Bessel integral
representation, integrating and using the beta integral, as in
the predecessor. [DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9)
and [5.9.13](https://dlmf.nist.gov/5.9.E13) were inspected
2026-10-20; the latter supports the retained digamma bound below.

Since R'=z^2/(1+z), R>=0, and

    (3R-zR')'=z^3/(1+z)^2>=0,
    (zR'-2R)'=z^2/(1+z)^2>=0,

the complete multiplier satisfies, for p>0,

    4K(p)<=p K'(p)<=6K(p),  K(p)>0.                 (4)

The retained all-frequency lower bound and benchmark give

    K(p)>=C p^4[log p+9-38/15]/(30pi),
    K(p)>=C kappa4 p^4/(30pi).

For p>=1, log(2+p)<=log p+log3. Thus Plancherel gives

    ||T||_4log <= (400/C)||K(|D|)T||_2+5||T||_2,
    |||D|^4 T||_2 <= (7/C)||K(|D|)T||_2.           (5)

No asymptotic approximation or angular truncation enters these
inequalities.

## 3. Freezing the local scale without spatial differentiation

At each point in O the flat operator for the frozen metric has
symbol K(exp(-sigma(z))|xi|). Call this point-quantized operator P.
For fixed real s, (4) gives

    ||[partial_s K(exp(-s)|D|)] K(|D|)^-1||<=6exp(6|s|).

The ratio has a finite extension at zero frequency. Integrate
this identity from 0 to sigma(z). Multiplication by the indicators
{sigma(z)>=s} or {sigma(z)<=s} has L2 norm at most one. Minkowski
therefore proves, for arbitrary bounded spatial sigma,

    ||1_O(P-K(|D|))T||_2
      <=2(exp(6a)-1)||K(|D|)T||_2
      <.36 ||K(|D|)T||_2.                           (6)

There is no spatial derivative of sigma in this argument.
Its quantization and the use of both signs of sigma are explicit.
Similarly the principal correction in (1) has norm at most

    [4exp(4a)/kappa4] ||K(|D|)T||_2
      <.33 ||K(|D|)T||_2.                           (7)

Consequently both changes together cost less than .7 of the
flat graph norm.

## 4. Geometric comparison, including all finite contacts

We now quantify Q_geo-P, where Q_geo is the actual complete
geodesic model. This comparison is on a FIXED conformal chart;
normal coordinates are not recentered separately inside (6).

For points in O and coordinate separation r<=1/1024, the analytic
seed charts give |grad sigma|<20, |sigma|<.027, and curvature
bounded by two. First-exit comparison with the straight segment
and integration of the metric along a minimizing path give

    |d(z,z')-exp(sigma(z))r|<100r^2,
    |exp(2sigma(z'))-exp(2sigma(z))|<100r.

The local exponential-map Jacobian differs from one by less
than 2d^2. All metric/distance ratios here lie between .9 and
1.2; enlarging the bounds slightly covers the transition edges.

The tower bounds from the curved proof and counting N(E)<=4CE^2
give

    sum d |B_M(r)| <=72C/r^6,
    sum d |B_M'(r)| <=1000C/r^7.

For the second bound use
B_M'<=2M^2 exp(-sqrt(2)Mr)/(pi^2 r^3) and
sum d M^2 exp(-tM)<=96C/t^4. The coordinate-density kernel
difference, with its metric Jacobians included, is consequently
bounded by

    |K_geo(z,z')-K_frozen(z,z')| <=1e12 C/r^5.       (8)

The finite-mode allocation fixes the local fourth coefficient:
the directional fourth Taylor moment has the same frozen
orthonormal angular integral, and its common restoration is
C kappa4 exp(-4sigma) Delta_0^2/(30pi).
All remaining contacts have differential order at most three.
We use this allocation as well as (8); a separated kernel bound
alone would not exclude an added fourth derivative of a delta.

Choose the small auxiliary separation R0=1e-20. Subtract the
Euclidean Taylor polynomial through degree three in the near
kernel. Taylor's integral formula, translations in the fixed
coordinate plane and Minkowski give the strong bound

    ||near remainder T||_L2(O)
      <=(2pi/24)1e12 C R0 |||D|^4 T||_2
      <1e-6 ||K(|D|)T||_2.                         (9)

There is no physical cutoff at R0; the complementary operator
and every contact are retained. The geodesic full-flow cosine
tail has row mass <=2pi(18C)/(.9R0)^4. The frozen flat tail has
both rows bounded by integrating its uniform C/r^6 envelope.
All finite restrictions of these tails have norm below 1e86 C.
This includes geodesics that leave and return, and any escape
from the coordinate patch. Coordinate and physical density
conversions fit the same bound.

To bound the remaining local coefficients, use all ten probes

    q_alpha(z')=chi(|z'-z|)(z'-z)^alpha/alpha!,
                     |alpha|<=3,                 (10)

with the beta(6,6) cutoff equal to one through 1/2048 and zero
by 1/1024. These supports lie inside |z'|<1/64. Their Euclidean
jets through order three at z form the identity matrix.
The cutoff/Leibniz bounds of the preceding strong remainder
proof, rescaled from outer radius 1/64 to 1/1024, give all
coordinate jets through six below 1e32 and supremum at most one.
The same finite geodesic-jet majorant then gives
sup|X^4q_alpha|,sup|X^6q_alpha|<1e56.
The complete model point action is below 2e56 C. The frozen
flat model point action has the same bound, including its
fourth restoration.

The near Taylor remainder on a probe costs less than 1e26 C;
the bounded tail on it costs less than 1e86 C. Evaluation of
Q_geo-P on (10), subtracting those two retained pieces, bounds
each of the ten contact coefficients by 1e87 C. This includes
the constant coefficient and coordinate-patch escape contact;
we do not assume that a zero extension preserves constants.
In particular no derivative of a contact bound is inferred.

Each coordinate derivative of order <=3 has L2 norm at most
||T||_H3, where H3 uses multiplier (1+|xi|^2)^(3/2).
Together these bounds give

    ||1_O(Q_geo-P)T||_2
      <=1e-6 ||K(|D|)T||_2+1e90 C ||T||_H3.        (11)

For supported T, the preceding actual H_mass-Q_geo estimate
applies throughout O. Cauchy bounds on sigma and its first
three derivatives, the coordinate formula for nabla^3, and
exp(|sigma|)<2 give

    sum_(j=0)^3 ||nabla^j T||_L2(gbar)
                              <1e12 ||T||_H3.

Its Euclidean output conversion is less than two. Thus replacing
Q_geo by H_mass in (11) costs less than 2e60 C ||T||_H3.
We retain 1e91 C for the total lower-order term.

## 5. The chart graph estimate and the Weyl coefficient

Since supp T is contained in V, the flat response outside O has
separation at least 1/256. Its local differential contacts vanish
there. The complete flat cosine kernel and Schur's test give

    ||1_(R2\\O) K(|D|)T||_2 <1e13 C ||T||_2.       (12)

Combining (6), (7), (11), the actual remainder and (12) leaves
more than one quarter of the flat graph norm:

    ||K(|D|)T||_2
       <=4||A_bT||_L2(O)+1e93 C ||T||_H3.

Equation (5) then gives the more useful intermediate statement

    ||T||_4log <=(1600/C)||A_bT||_L2(O)
                                      +1e96||T||_H3. (13)

For any epsilon>0 and log2>1/2, splitting the Fourier variable
at sqrt(1+|xi|^2)=1/(epsilon log2) proves

    ||T||_H3 <=epsilon ||T||_4log+(2/epsilon)^3||T||_2.

Take epsilon=1/(2e96), absorb one half, and enlarge the explicit
low-norm constant to 1e390. This proves (2).

Finally b=omega+1/6 is within the allowed interval on O.
One must not extend the earlier very-short-collar inequality to
Y=1/8 without checking it. Use instead the exact seed polynomial
f=F/(lambda r^2)=1+nu u^2-2a_seed u^3+b_seed u^4-2e u^5+e u^6,
where u=B/r, nu>=0, a_seed<2.002, b_seed>0 and |e|<1e-15.
For 0<=u<=1/4, f>.9. Since C/x<6.1,

    Y(u)=integral_0^u sqrt(C/(x f(v)))dv <3u
                                 on this interval.

At u=1/28 this is below 1/8. Monotonicity therefore gives
u(Y)>1/28 for Y>=1/8, so omega<log28<10/3 and b<7/2<4.
This checks the actual seed Weyl PRINCIPAL coefficient. The
derivatives of b in its divergence-form completion still belong
to the lower local operator and are not supplied by (7).

## 6. Verification and consequence

Freeze proof and checker before the scientific run. Verify the
all-k flat moment normalization, exact derivatives of R,
continuum scale integration, rational exponential bounds, full
contact probe jets, rescaled cutoff and geodesic reserves, kernel
and tail constants, Taylor identities, interpolation and the
actual seed bound on b. Pin the predecessor receipts and exact
inputs. Preserve every scientific run and failure. These finite
checks supplement the ordinary continuum argument.

This gives a numerical local principal graph comparison, with
an actual massive-sector case b=0. Its large L2 constant cannot
be called a useful global inverse bound or a Newton contraction.
The intermediate estimate (13), the sharper retained form
coercivity and an eventual direct residual/inverse estimate are
the appropriate starting points for further work.

Remaining to-do list: incorporate the full local TT operator,
patch commutators and higher York maps; quantify the outer fast
fourth/log estimate and full inverse/residual composition;
nonlinear signed compatibility and all exterior, opening,
stability, transport, preparation and creation requirements.
