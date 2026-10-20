# The complete outer massive response compared with its flat logarithmic term

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Massive comparison
controlled on localized fast input; full disk assembly remains open.**

The curved geodesic model can be compared strongly with the exact
flat multiplier after a radial half-density conjugation. Explicit
point probes determine its finite contacts. The conjugation is
quantified in the chosen fourth/log norm; it does not change the
physical equation or improve an inverse by an arbitrary rescaling.

## 1. Operators, norms and result

Retain the full seed, angular tower, physical time, state, paired
domain and common finite subtraction. Use the operators of
[the actual-to-geodesic estimate](outer-strong-curved-quantum-remainder.md)
and [the flat tower estimate](outer-flat-tower-normal-remainder.md).
Put

    j(Y)=a(Y)/a(0), m(Y)=sqrt(j(Y)), xi=lambda B^2,
    rho=C^2 u^4/(xi^2 Y^4)<3/2, delta=7/2, x=log Y,
    K_infty(p)=C p^4[log(p/(mu0 ell_o))-31/30]/(30pi).

K is the complete flat geodesic multiplier, with all angular
levels and its fixed fourth restoration. Work first with smooth
even stationary input U supported on 0<Y<epsilon, epsilon<=1e-15.
Its extension to the double is understood. All norms of the
weighted input use zero extension on the x-line. Then

    ||Y^-7/2 1_(Y<epsilon) rho Y^4/(16pi C)
          [m Q_geo m^-1-K]U||_L2(dx)
      <=1e44 epsilon ||Y^-7/2 U||_H4log(dx),           (1)

and, after adding the actual remainder and flat tower correction,

    ||Y^-7/2 1_(Y<epsilon) rho Y^4/(16pi C)
          [m H_mass m^-1-K_infty]U||_L2(dx)
      <=1e65 epsilon ||Y^-7/2 U||_H4log(dx).           (2)

Conjugated operators here mean: multiply the supported input by
m^-1 in the outer chart, extend by zero on the actual double,
apply the COMPLETE operator, and multiply its outer output by m.
No division by m is made at a bolt. Geodesic returns are retained.
These DIFFERENCES extend by the displayed bounds. No separate
weighted boundedness of H_mass or K_infty is asserted before
the Weyl logarithm is combined with them.

The physical equation may be written in auxiliary variables
(theta_hat,b_hat)=M_epsilon(theta,b), f_hat=M_epsilon f, with
M_epsilon=m on the input/output collar and M_epsilon=1 past
2epsilon. This is an operator similarity, not a replacement
dual variation of the physical action. Section 5 proves that
M_epsilon and its inverse have norm less than 1.01 on the
specified mixed H4/H4log space and on the output L2 space.
Thus the original physical norms remain explicitly comparable.

In these auxiliary variables the complete normal inverse gives
the component budget

    ||N^-1 lambda^2 (0,massive_difference)||
        <8.7e84 epsilon <1e-30 at epsilon=1e-120.       (3)

Equation (2) still retains the actual rho in front of K_infty.
The local Weyl/principal-density combination, other local disk
terms, localization and global matching are not included in (3).

## 2. Proper radial half-density and the separated kernel

In proper Y coordinates the radial Laplacian is
Delta=partial_Y^2+h' partial_Y, h=log j. Conjugation by m gives

    m Delta m^-1=partial_Y^2-V,
    V=h''/2+h'^2/4.                                  (4)

Set Qtilde=m Q_geo m^-1 on the outer chart. It is symmetric for
dY, as is K on the even line. The midpoint radial comparison in
section 3 of [the curved proof](stationary-curved-quantum-tt.md)
applies to their short-flow kernels. At the symmetric midpoint,
normalize the time coordinate by a_mid and put r^2=s^2+tau^2,
s=Y-Y'. The retained estimates are

    |d-r|<=17r^3,
    |sqrt(a(Y)a(Y'))/a_mid-1|<.3s^2,
    |J_exp-1|<=d^2/3.

Together with the Bessel derivative bound, the per-mode
distance/density error is below 30 M^2 exp(-Mr).
Localize the short flows with a smooth cutoff equal to one
through 1/256 and zero by 1/128, slope at most 512.
Its distance-change error is bounded by the same Bessel
envelope times 512*17r^3. Including this error and the density
reserves, 1e4 M^2 exp(-Mr) is sufficient.

The counting estimate sum d_n M_n^2 exp(-M_n r)<=96C/r^4
and integral_R (s^2+tau^2)^-2 d tau=pi/(2|s|^3) therefore
give a short-flow radial difference kernel with

    |K_rem(Y,Y')|<=1e9 C/|Y-Y'|^3.                   (5)

The near separation used below is R=1/400. Cutoff transitions
and all complementary flows have length at least R/2 when
placed in the bounded part. The geodesic tail bound
2pi(18C)/(R/2)^4, its flat counterpart and the outer density
ratios give both total-variation row bounds below 1e20 C.
This statement allows the geodesic pushforward measure at
returns; it does not assume a global injective exponential map.
Liouville invariance bounds both rows, also after restriction
to the outer chart. The finite radial cutoff is auxiliary.

## 3. The common finite contacts and exact probes

The common fourth coefficient is constant in proper half-density
coordinates. For example the restored Q4 operator, divided by
C kappa4/(30pi), is exactly

    (partial_Y^2-V)^2
      +(1/3)[Rbar partial_Y^2+Rbar' partial_Y
                           -Rbar V-Rbar' h'/2].        (6)

Its leading term is partial_Y^4, with NO third derivative.
The flat restoration has precisely the same fourth coefficient.
At each finite angular cutoff the original negative fourth
moment cancels the fourth Taylor subtraction in both models;
the remaining mode contacts have order at most two. The common
finite-part prescription preserves this allocation. Together
with (5)-(6), this leaves only contacts of order at most two
in Qtilde-K. In particular no new fourth delta derivative is
allowed. Symmetry is also consistent with the absence of a
third contact when the common fourth coefficient is constant.

Define the near Taylor operator in the flat Y measure by

    T3 U(Y)=integral_(|s|<R) K_rem(Y,Y+s)
            [U(Y+s)-U(Y)-sU'(Y)-s^2 U''(Y)/2]ds.

The bounded complementary kernel/measure operator is B.
Keep all its row sums; do not assume it annihilates constants.
For outer observation and outer supported input the exact form is

    (Qtilde-K)U=c0 U+c1 U'+c2 U''+T3 U+B U,           (7)
    ||B||_(2->2), sup row_abs(B)<=1e20 C,
    ||T3 U||_2<=1e7 C ||U'''||_2.                    (8)

For the last estimate, use the integral Taylor remainder,
translations on the Y-line and (5): 2(1e9)R/6<1e7.
No derivative of an unsigned kernel bound is used.

Fix an observation Y in (0,epsilon), and take

    phi_k(Y')=chi(Y'-Y)(Y'-Y)^k/k!, k=0,1,2,
    chi=1 for |s|<=R, chi=0 for |s|>=2R,
    chi=1-I_(|s|/R-1)(6,6) on the transition.

The jets through two at Y form the identity matrix and T3 phi_k=0
exactly. Their supports lie in |Y'|<1/64. For the full geodesic
point action use the base probes q_k=m^-1 phi_k. These are C2,1
across the seam. The beta cutoff derivative bounds and the one-sided
analytic seed bounds give, conservatively,

    sup |X^3 q_k|<1e30 on the complete flow,
    sup_(flow length<Y/2) |X^n q_k|<1e30, n<=6.        (9)

Here are explicit reserves behind (9). The cutoff derivatives
through three are at most 3960 R^-3 and its derivatives through
six at most 10644480 R^-6, with the smaller derivative orders
included by these upper bounds. Polynomial product rules give
|partial_Y^l phi_k|<1e15 for l<=3 and <1e30 for l<=6.
The one-sided analytic functions m^+/-1 have modulus below two
on |Y|<.1, using the inverse disk of the preceding proof.
Thus their first three derivatives on the probe support are
below 2 l!20^l<=1e5. Across the seam these are essential bounds
through three, not a claim of a bounded fourth derivative.
The radial third covariant derivative identity and the bounded
h',h'' give the first line of (9), including the zero exterior
of the support. For short flows the cutoff is one, and
m^-1(Y')(Y'-Y)^k/k! has analytic modulus below two on the
centered 1/20 disk. Twice the preceding geodesic Bell majorant
is below 1e30 through order six. The short flow stays on the
positive analytic side.

For r<Y/2 use the sixth-order even remainder. For r>=Y/2,
the GLOBAL third-derivative bound gives

    |q_k-cos(rP_flow)q_k+(r^2/2)X^2q_k|
                                          <=r^3(1e30)/6.

The tower envelope 72C/r^6, with measure r dr dtheta, integrates
this part to less than 1e33 C/Y. The quartic polynomial and its
restoration use the envelope 1152C r^-6 exp(-9r), giving at most
400 C(1e30)[1+log(2/Y)]+10C(1e30). The short part is below
1e30 C Y^2. Since 1+log(2/Y)<3/Y on this collar,

    |Q_geo q_k(Y)|<1e35 C/Y.                         (10)

This is a full-flow bound, including escape and return.
One may smooth the seam and then pass to this absolutely
convergent separated part to define the point action; no global
fourth/log domain for q_k is assumed. The flat point bound on
phi_k from its global sixth jet is below 1e32 C. Multiplication
at the observation by m<2 and the row bound for B now give

    |c_k(Y)|=|(Qtilde-K)phi_k-B phi_k|<1e40 C/Y,
    k=0,1,2.                                         (11)

The constant coefficient, including escape terms and half-density
contacts, is explicitly retained. It is not inferred from K(0)=0.

## 4. Weighted comparison with the flat model

At delta=7/2,

    ||Y^-delta Y^4 F||_L2(dY/Y)=||F||_L2(dY).

For U=Y^delta w(log Y), the three weighted contact derivatives
satisfy

    ||Y^-1 U^(k)||_L2(dY)
       <=40 epsilon^(3-k) ||w||_H3(dx), k=0,1,2,
    ||U'''||_L2(dY)<=40 epsilon ||w||_H3(dx),
    ||U||_L2(dY)<=epsilon^4 ||w||_L2(dx).

These follow from the exact polynomials
product_(l=0)^(k-1)(partial_x+7/2-l). Even extension costs
the same square-root factor on input and output; there is no
new half in the even model. Using rho/(16pi)<1, (7)-(11) give
less than

    [120e40 epsilon+40e7 epsilon+1e20 epsilon^4]
                                           ||w||_H3.

Finally H3<=2 H4log proves (1), with room below 1e44.
The proof controls the entire infinite end, not sampled scales.

## 5. Explicit norm conversion and the complete massive difference

For either f=m or f=m^-1, the one-sided analytic continuation
has f(0)=1, f'(0)=0 and Taylor coefficients |f_n|<=2*20^n.
For D=Y partial_Y and k<=5, the positive series gives

    |D^k(f-1)|<=2 sum_(n>=2) n^5 (20Y)^n
                 <=2(4326)(20Y)^2 <1e7 Y^2.          (12)

Indeed sum_(n>=2) n^5 (1/2)^(n-2)=4326; the ratio is increasing
in 20Y<=1/2. Define eta(x)=1 for x<=0, zero for x>=log2,
with the same beta(6,6) transition. Its derivatives through
five are below 1e8, since log2>1/2. Set

    M_epsilon(Y)=1+eta(log(Y/epsilon))(m(Y)-1).

Leibniz and (12) show that every log-coordinate derivative
through five of M_epsilon-1 is below 1e18 epsilon^2.
Multiplication therefore has H3 and H5 operator norm below
1e21 epsilon^2. The explicit orthogonal-block argument of
[the compact cover proof](compact-cover-graph-absorption.md)
costs less than 1e6 to pass from these bounds to H4log, and
less for H4. Thus, for either input component,

    ||M_epsilon-1|| <=1e27 epsilon^2<=.001.            (13)

The Neumann series for this multiplication operator gives
||M_epsilon||, ||M_epsilon^-1||<1.01. Its pointwise inverse
agrees with m^-1 on the collar. The output L2 bounds follow
directly. This also gives a positive global extension equal
to one past 2epsilon. The leading mass cubic is unchanged
because M_epsilon=1+O(Y^2). These are auxiliary variables for
the same equation, with explicit norm conversions.

For supported U, the already proved actual-to-geodesic estimate
applied to m^-1 U and multiplied by m costs less than
2e64 epsilon in the present norm. Add (1). The retained flat
comparison K-K_infty, multiplied on output by rho, costs less
than 1.5e30 epsilon in the same physical shape normalization
(we may drop its favorable factor 1/(480pi^2)). Therefore

    2e64+1e44+1.5e30 <1e65

proves (2). All actual massive contacts and state/global terms
are now inside this difference. They are not separately added
again. The normal inverse bound 8.7e19/lambda^2 gives (3).

## 6. Verification and remaining work

Freeze this proof and checker before execution. Pin the actual
remainder, flat correction, compact logarithmic transfer and
curved-model receipts. Check the half-density Laplacian and
fourth restoration, two-order kernel budget, contact-probe jets,
global/one-sided flow reserves, point integrals, weighted powers,
the exact 4326 series, logarithmic multiplication constants,
Neumann conversion and complete massive/normal composition.
Preserve every run and any failed sources. The finite checks
support the ordinary operator and finite-part arguments; they
are not formal verification or a solution of the field equations.

Next combine rho K_infty with the actual local Weyl term before
estimating its growing logarithm, and include the remaining local
disk blocks. Mass/fast localization, nonlocal cutoffs, annulus
matching, higher York maps and full inverse/residual composition
are still required. No global or nonlinear inverse is claimed.

Remaining to-do list: local Weyl/principal-density and disk
assembly; localization and annulus matching; higher reconstruction
and inverse/residual composition; nonlinear signed compatibility;
exterior/state, charge, return, opening, stability, transport,
preparation and creation.
