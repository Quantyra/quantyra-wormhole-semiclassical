# Weighted cold-cap inverse derivatives through eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
This supplies the global auxiliary inverse part of the cold eighth-order
source calculation. It does not supply the complete cold source, the
actual finite-domain inverse derivatives or any nonlinear solution.

Retain the original centered seed, temperature, paired domain and physical
subtraction. The construction below changes an auxiliary complete H2 cap
only. Use the exact uniformization and seed guards from
[the cold construction](stationary-reference-cold-bolt.md),
[the centered source](stationary-centered-exterior-source.md), and
[the fixed-seed derivative calculation](stationary-centered-canonical-response.md).
Distances and the Dirac operator in this note use unit H2.

For a cold observation x, put

    h<=s_x<=.001, a_x=s_x[1+log(.01/s_x)],
    W_x(y)=exp[-2 d_H2(x,y)], M>=18.

The explicitly defined parameter family below has two real parameters b.
Let D_H2 be anti-self-adjoint, and define the canonical inverse

    T_b=D_H2+M exp(sigma_b),
    S_b=exp(sigma_b/2) T_b^-1 exp(sigma_b/2).             (1)

The endpoint factors in (1) are part of the operator. The result is

    M ||W_x (partial_b^j S_b)|_0 W_x||
       <=a_x j!15000^j B_j, 1<=j<=8,
    B_j<2000, B_8<1987.                                (2)

Here the derivative norm is the square root of the sum of squared
operator norms over all ordered parameter indices. In the source's
usual derivative normalization, the eighth entry is consequently

    M ||W_x (partial_b^8 S_b)|_0 W_x||/(8!100^8 a_x)
       <5.093e20.                                     (3)

The zeroth comparison with the homogeneous inverse S_H2 obeys

    M ||W_x(S_0-S_H2)W_x||<1000 a_x.                   (4)

Equations (2)-(4) concern the global auxiliary inverse before the
additional local normalizing-coordinate transport. They are not bounds
on eight derivatives of the renormalized diagonal source. In particular
the earlier transported C3 estimate is not replaced by (2).

## 1. Specify the observation motions

Rotate the bolt-centered hyperboloid coordinates so
x=(C_x,S_x,0), C_x=cosh(chi_x), S_x=sinh(chi_x). Write
c_x=1/C_x, v_x=S_x/C_x. Use Lorentz boosts in the first and second
spatial directions and set

    L_b=B_1(b1) B_2(c_x b2), sigma_b(y)=sigma_cap(L_b y). (5)

The coefficients c_x and v_x are fixed at the base observation.
This is a real H2 isometry for real b. Its pullback keeps D_H2 fixed,
with the canonical spin identification, and moves the observation to
L_b x in the original capped geometry. Both parameter vectors at zero
have unit H2 length and are orthogonal. No claim of a physical normal
coordinate chart is made.

The exact H2 metric of b -> L_b x is

    g11=1+u^2, g12=-v_x u, g22=1,
    u=sinh(c_x b2)/c_x, det g=cosh(c_x b2)^2.           (6)

Thus this chart remains regular at the bolt and as c_x becomes small.
At the bolt c_x=1 and v_x=0. The c_x->0 expression is its regular
limit u=b2. Without the factor c_x, the second parameter would have
speed C_x at the observation. Using (5) makes that factor explicit.
Conversion to the physical covariant derivatives still has to retain
the connection of (6) and the seed's conformal factor.

For any real y on the hyperboloid, H_y=h y0=h cosh(chi_y). On the
complex bidisk |b1|,|b2|<=rho=1e-4, the matrix exponential series gives

    |H(L_b y)-H_y|<=H_y[exp(2rho)-1]<.000201 H_y.        (7)

Indeed |y1|,|y2|<=y0, and each boost generator has maximum row-sum
norm one; c_x<=1. This estimate is uniform in y and chi_x.

## 2. A cap with true finite jets and a weighted branch bound

Keep the analytic seed through chi1=arccosh(.01/h). For
t=(chi-chi1)/.01 in [0,1], set

    sigma_cap(chi)=eta(t) sigma_seed(chi),
    eta(t)=sum_(i=0)^10 binom(21,i)t^i(1-t)^(21-i),       (8)

and set sigma_cap=0 beyond the transition. The seed in (8) is its
analytic continuation across the seam, not a polynomial extrapolation
whose eleventh coefficient has been omitted. The step is C10 with
matching endpoint derivatives through ten. This cap is complete and
bounded, equals the original seed on the retained cold patch, and
has the higher finite regularity needed here.

Here are explicit bounds for its branch construction. The exact cold
formulas write

    exp(sigma)=A(s)[cosh(hD)+(s/h)sinh(hD)],
    H(s)=h cosh(chi)
        =[s cosh(hD)+h sinh(hD)]/[cosh(hD)+(s/h)sinh(hD)].

On the quarter s-disk |s-s0|<=(s0+h)/4, h<=s0<=.01, the earlier
rational formulas give

    |D|<=2.03 log[(.01+h)/(s0+h)]+4.5,
    |D'|<=4.5/|s+h|, |log A|<=4.1s0.

These deliberately loose bounds also cover .001<s0<=.01. Explicitly,
the divided-difference formula for q0 gives |q0-1|<=2.12s0,
|q0'|<1.42 and |J|<1.06 on these quarter disks. Hence
|J'|<1.06(2.04+1.46)<4.5. Integrating J' from h removes the factor
s-h in D', and the quarter-disk path adds less than 1.5 to the real
D bound. The log A estimate follows from |log q0|<=2.166s0 and
2|log(1+s)-log(1+h)|<=(2.5s0-1.5h)/(1-1.5s0), together with
|log q0(h)|<.501h. The centered perturbation fits the same margins.

In particular |sD|<=s0[3.045 log(.01/s0)+6.75]<=.0675.
Using |log(1+z)|<=|z|/(1-|z|), the hyperbolic-function series and
h<1e-30 gives |sigma|<=12 a_s0 on this entire quarter disk.
The negligible h terms are bounded by 1e-24 in the following derivative
estimate, rather than set to zero:

    |H'-1|<= (.135+2(.0675)+.0675^2)/(1-.0675)^2+1e-24
           <.32.                                    (9)

For (9), |s^2 D'|<=13.5s0<=.135; the leading expression is the
derivative of s/(1+sD). Differentiating the exact quotient above bounds
its omitted h terms by 1e-24, using |hD|<1e-26 and
h^2|D'|<=3h. The seed's degree-six perturbation is included by the
retained centered q0 estimates and fits these strict margins.

Equation (9) supplies an inverse on the image disk of radius at least
.68(s0+h)/4. The shift in (7) therefore changes s by less than
.000402s0. Cauchy on the remaining quarter-disk margin gives
|sigma_s|<49 a_s0/s0 along that shift. Since the retained real inner
bound is |sigma|<.021, its branch value is below
.021+49(.000402)(.01)<.022. The weighted bound 12a_s0 also remains
available from the full quarter disk.

At the transition, use the retained |sigma|<.1 on the s-disk of
radius .0075 about .01, and |sigma(.01)|<.02. The exact radial
equation ds/dchi=exp(sigma)sqrt(P) has modulus below .02 there.
For a complex chi displacement at most .0105, its solution stays
within .00021 of .01. Thus its sigma value is below

    .02+.00021*.1/(.0075-.00021)<.023.                (10)

The arcosh branch at this distant seam turns (7) into a chi shift
below .0005, hence a step-argument shift below .05. With
N=21!/(10!)^2=3879876,

    |eta'|<=N(1/4+2(.05)^2)^10<4.6,
    |eta|<=1+.05(4.6)=1.23.

Together with (10) this bounds the transition branch by .03.
For reference, the real derivative is below ten: the seed derivative
costs less than .02*.1/(.0075-.00021), and the taper term costs less
than 100(3.8)(.023). No high derivative is assigned this low bound.

At a fixed real y use the analytic branch of its actual piece of
(8), and a constant branch outside the cap. At a transition endpoint
the neighboring branches have the same jets through ten. These
branches bound true parameter jets through eight; they do not define
a globally analytic spatial cutoff or a complex geometric cap.

The cold weight a(s)=s[1+log(.01/s)] obeys |(log a)_chi|<=1 in
the inner patch. This follows from the retained s_chi/s<=1+3sL
and 3s L^2<1, keeping the exact radial correction. Extend the weight
past chi1 by .01 exp(chi-chi1). Then

    a(y)<=a_x exp[d_H2(x,y)].

The inner bound 12a(y), the transition bound .03 with a(y)>=.01,
and the zero exterior therefore imply the uniform branch bounds

    |sigma_branch(L_b y)|<.03,
    |sigma_branch(L_b y)|<=100 a_x exp[d_H2(x,y)].       (11)

The real family is C8 in the multiplier norms used below. Uniform
finite-jet remainders follow from the C10 endpoint matching and the
same local analytic branch estimates. No complex global inverse is
invoked to differentiate this real family.

## 3. Coefficient jets and real weighted inverses

Put E_b=exp(sigma_b/2), V_b=exp(sigma_b). Cauchy on (11), applied to
E_b-1 and V_b-1, gives, for 1<=j<=8,

    ||partial_b^j E_b||_infty, ||partial_b^j V_b||_infty
        <=u j! k^j,
    ||W_x(partial_b^j E_b)W_x||,
    ||W_x(partial_b^j V_b)W_x||
        <=w a_x j! k^j,
    u=1/10, w=200, k=15000.                           (12)

The factor k exceeds sqrt(2)/rho and accounts for all ordered
parameter indices. The bounds exp(.03)-1<.1 and
100 exp(.03)<200 justify the unweighted and weighted constants.
The same w a_x bound holds for W_x(E_0-1)W_x and W_x(V_0-1)W_x.
Also ||E_0||<e0=51/50.

For every real cap in this family, T_b has positive real part at
least M exp(-.03). Conjugation by exp(plus or minus 2d_H2(x,.))
changes D_H2 by a multiplication operator of norm at most two.
Therefore both weighted similarities of T_b^-1 have norm at most

    [M exp(-.03)-2]^-1<d/M, d=6/5.                  (13)

Use truncated Lipschitz distance weights and then their strong limits
to interpret the similarities with unbounded weights. The bound is
uniform in the truncation; scalar weights preserve the Dirac domain.
The same estimate holds for the homogeneous inverse. This calculation
requires no large curvature estimate for the new narrow taper.

## 4. Keep one weighted factor in every ordered term

The jth derivative of T_b^-1 is the finite sum over all ordered
compositions i1+...+il=j, with coefficient j!/(i1!...il!), of

    (-1)^l T^-1 (M V^(i1)) T^-1 ... (M V^(il)) T^-1.   (14)

There are binom(j-1,l-1) such compositions of length l. An unweighted
term is bounded by d^(l+1) u^l j! k^j/M. For its weighted version
reserve one differentiated multiplier with the w a_x bound of (12).
Move W_x and W_x^-1 through the other scalar multipliers and use
(13) on every inverse. For example, putting the weight at the first
vertex uses W_x T^-1 W_x^-1 to its left and W_x^-1 T^-1 W_x on its
right. This costs the same d powers and replaces exactly one u by
w a_x. It does not require the product of several small weights.

After division by j!k^j/M, the inverse majorants are

    r0=d,
    r_j=d^2 u(1+du)^(j-1), j>=1,

with the weighted replacement factor w a_x/u for j>=1.
For the endpoint factors let e_j=u for j>=1 and e_0=e0. The full
Leibniz sum in (1) is bounded by

    U_j=sum_(a+b+c=j) e_a r_b e_c,
    B_j=(w/u)U_j, j>=1.                            (15)

Every term for j>=1 has at least one positive derivative order, so
the reserved weighted factor is always present. It can lie in either
endpoint or in the inverse. No endpoint derivative is discarded.
The generating function is an independent way to check the sum:

    sum U_j t^j
      =d[e0+(u-e0)t]^2/[(1-t)(1-(1+du)t)].           (16)

For the zeroth difference in (4), expand the two endpoint differences
and use T_0^-1-(D_H2+M)^-1=-T_0^-1 M(V_0-1)(D_H2+M)^-1.
The resulting bound is w[d(e0+1)+d^2]a_x/M=772.8a_x/M,
which fits the stated reserve 1000.

Exact rational evaluation of (15), rounded upward, gives:

| j | B_j | B_j 150^j |
| --- | ---: | ---: |
| 1 | 790 | 1.184e5 |
| 2 | 908 | 2.043e7 |
| 3 | 1041 | 3.514e9 |
| 4 | 1190 | 6.024e11 |
| 5 | 1357 | 1.031e14 |
| 6 | 1544 | 1.758e16 |
| 7 | 1753 | 2.995e18 |
| 8 | 1987 | 5.093e20 |

The last column rounds the exact B_j, rather than the rounded middle
column. This proves (2)-(4). It supplies a uniform quantitative inverse
estimate at the bolt as well as farther along the cold patch.

## Verification and next use

Freeze the note, checker and review before the supporting run. Check
the Lorentz chart and its normalization, cutoff endpoint matching,
complex and weight budgets, the weighted similarity margin, every
ordered composition and both endpoint factors. Compare the finite
insertion formula with an independent noncommuting matrix inverse
through eight, and compare its coefficients with complex contour
integrals at two precisions. These fixtures do not evaluate the
charged field's diagonal source or formally verify continuum analysis.

The next source calculation must include the local renormalized term,
both localized rows, the normalizing map and density, the full angular
tower, the actual-domain comparison and the physical connection of
(6). It must preserve the same a_x factor in their sum. The existing
complete cold C3 allowance and all endpoint C5 contracts are unchanged.
The middle and outer source C8 results remain available. No complete
cold C8, metric order ten, finite-metric H6, transfer target, signed
compatibility, nonlinear existence or physical stability is claimed.

Remaining to-do list: complete cold weighted source derivatives through
eight; quantitative ten-derivative metric and finite-metric quantum/Schur
source estimates; positive transfer or a sharper actual transfer;
cancellation, signed compatibility, nonlinear existence; physical
validity, state/boundary, stability, matching, charge/return, opening,
transport, preparation, separate creation and audit.
