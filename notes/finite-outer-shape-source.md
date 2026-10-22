# Complete outer shape forcing and the global physical H4 source

2026-10-22, E003/S167. Ordinary argument and verification contract before
execution. Keep the original finite comparison, physical charge, HH
state, paired domain, angular tower, common restoration and fixed return.
Use the notation and exact orbit identity of the
[interior source argument](finite-interior-shape-source.md):

    O_ell=K_ell+ell Q_ell, O_0=dot O_0=0,
    Q_ell=(C/x) z_ell^4(S_M,ell+S_LLL,ell+S_W,ell),
    b_src=-O_ell/(16 pi g)-f_tau/6, g=G lambda_e.

Let O be the fixed seed region t>=40. The indicator in a restricted
norm is applied AFTER taking the actual physical tensor derivatives.
This argument supplies a finite numerical outer bound, combines it with
the interior bound, and thus quantifies the global F_b4 and F4. The
constants are deliberately conservative. They do not establish a useful
nonlinear correction, preservation of a neighborhood, or stability.

## 1. The cancellation needed for this norm

In the seed proper collar write y=Y0, L=1+log(1/y). The seed guards give
z0/y<.42, y/z0<3, z0'<.42. Thus t>=40 implies y<.075. The actual proper
coordinate X=Phi_ell(y) obeys |X/y-1|<3e-10, so X<.076.

Use the original seed and complete tangent equations BEFORE estimating:

    O_ell = K_ell-K_0-ell dot K_0 + ell(Q_ell-Q_0).       (1)

This is an identity for the original physical residual. It does not
replace that residual by an artificially subtracted source. In particular
the constant outer vacuum cancels with the affine classical allocation
in the seed/tangent equations. The radius multiplier of every quantum
term is still the actual z_ell^4.

A remainder of order y^3 L^2 is sufficient here. Physical radial
derivatives have scale y partial_y and the volume density has order
y^-4 dy. Its square therefore gives an integrable y^2 L^4 dy. This
does not supply the stronger weighted mass-integral cancellation: no
absolute integral of an individual cubic canonical or Weyl mass tail is
assigned a finite value.

## 2. A finite dilation family for the actual outer quantum source

The [inverse endpoint theorem](inverse-endpoint-ten-jets.md) supplies,
for each REAL ell, the proper logarithm b=log(f_ell/f_0), f=1/a, with
cubic b3=-ell u3 and remainder

    |b3|<ell 1e-11,
    N_(4,5)(b-b3 X^3)<ell 1e-7,
    N_(4,10)(b-b3 X^3)<ell 1e17,
    N_(m,J)(h)=max_(a<=J)|h^(a)|/[a!100^a X^(m-a)L(X)].

The small C5 and large C10 allowances remain distinct. We never
differentiate a moving proper-coordinate map in ell.

Take the degree-21 C10 step, equal to one through X=.089 and zero
from X=.099, and cap b by multiplication with it. The derivative
bound is |eta^(r)|<=4200^r through ten. Its normalized product sum
on X<=.1 is below sum_(r=0)^10 4.2^r/r!<67. The cubic cutoff terms
are supported away from zero and have an additional 1/X<12. Consequently
the capped remainder has low N_(4,5)<ell 1e-5 and high
N_(4,10)<ell 1e20. Its value/gradient extension fits 100 times these
norms. One may use the analytic mass direction .01 X^3, tapered only
outside X=.1; its coefficient is 100 b3. This has the same majorant-ten
and cap value/gradient contract as the mass direction in the retained
heat proof. The coefficient 100 b3 and the compensating exterior taper
fit the displayed low and high reserves.

Let D=X partial_X. For j<=5 and a<=5, differentiating D^j gives the
explicit coefficient bound

    K_(j,a)=1e20 sum_(k=0)^j S(j,k)
        sum_(r=0)^min(a,k) binom(a,r) k!/(k-r)!
                   (a+k-r)!/a! *100^(k-r),             (2)

where S(j,k) are Stirling numbers of the second kind. The cubic is
separated before this formula; its finite cutoff contribution fits
the factor 1e20. Formula (2) uses at most ten derivatives.

On |s|<=rho=1e-40 define P_s b=sum_(j=0)^4 s^j D^j b/j!. The exact
finite sums give

    N(P_s b)<=ell [1e-5+sum_(j=1)^4 rho^j max_a K_(j,a)/j!]
              <ell 2e-5.                              (3)

Here N includes the separated mass coefficient and cap value/gradient
contract. The fifth D derivative, supplied by C10, bounds the Taylor
remainder of the TRUE real dilation in this C5 weighted space. Thus
P_s b has its actual four dilation jets, including the five coefficient
derivatives required by the heat calculation. No complex spatial value
of b is used.

Use the complete seed cap from
[the seed outer source proof](seed-outer-eight-jets.md), but replace its
auxiliary scale family globally by exp(P_s log f_0). Its nonconstant
support is in |X|<=100 and it has the retained branch bounds through
eleven. Finite differentiation gives a change below 1e-20 in the
global first-three relative coefficients on this disk. Near the outer
boundary this seed family is spatially analytic, with the original
strict strip and heat-cone margins. The seed cap itself need not be
globally analytic. This global finite polynomial matches the true four
dilation jets on the WHOLE cap, so there is no omitted positive-seam
dilation defect. Both families use exactly the same seed polynomial.

The real heat scaling identity is

    k_f(exp(2s)t,exp(s)X;beta)
       =exp(-2s) k_(f_dil)(t,X;beta exp(-s)).           (4)

Accordingly keep Atilde(exp(2s)t), the scaled thermal circle and
exp(-4s) times the common k2 restoration. The fixed a0 scalar cancels
in the endpoint DIFFERENCE; it has not been set to zero in the source.
The cubic term is exp(-4s)c_s C^2/(5 pi X), where c_s has the true
four dilation jets, hence this term has the jets of c C^2/(5 pi exp(s)X).
Masses, physical scale and state are auxiliary-scaled only inside this
exact identity; the physical comparison data are unchanged.

Apply the complex-cap estimate of
[the finite-background response](outer-finite-background-response.md)
to the endpoint difference (3). Its small-time proof uses only C5
weighted real spatial jets, and is uniform under the seed polynomial's
1e-20 perturbation. The longer seed cap has ratio 1e14 rather than
the earlier response cap's 1e10. Enlarge its long-time row from
1e20 to 1e24 and its fixed-exit row from 1.536e18 to 1.536e22.
Together with the analytic-mass short row 2.4e24, the other retained
rows still sum to less than 5e24, below the 1e25 response reserve.
The scaled angular exponent has real part greater
than .9999; even the larger retained heat growth leaves the old
.9998 rescaling margin. The strict coefficient, inverse-cone, angular
moment and winding budgets cost less than two in total. Therefore
the normalized regular cap difference has the holomorphic bound

    |Delta R_cap(s,X)| < 1e-11 ell L(X)^2, |s|<=rho.   (5)

For example 2*(2 eta_G)*1e25*(2e-5)<1e-11. This is a bound for the
complete angular limit after the common subtraction, not a finite
mode calculation. The dominated complex continuation is proved in
the cited response argument. Simon's real Feynman--Kac theorem
(https://arxiv.org/html/math-ph/9907022, theorem 1.1, inspected
2026-10-22) supplies only the underlying real bridge identity;
it does not itself prove this complex or renormalized assertion.

## 3. Actual-domain correction with the shorter separation retained

The new endpoint cap agrees with the actual proper geometry through
X=.089. The observation is X<.076. Use a gluing cutoff supported
between .084 and .088, with first derivative below 1e4. Retain the
actual middle inverse and the constant paired boundary projector in
the two-sided resolvent identity. No complete-half-line replacement
of the physical state is made.

Set K=M+1e4. The point-evaluation ball 1/(10K) fits before the
commutator even when the original 1/(10M) ball does not. The retained
Caccioppoli/Fourier argument bounds each complex separated row by
3e5(K/M) times its exponential. The two rows, source mass, density
and cutoff give the safe per-mode value bound

    1e16 (K^2/M) exp(-M/1000).                       (6)

The exponent uses the short separation above, NOT the distance 24 of
the old long cap. The local real/complex relative coefficients obey
|f'/f|<5, |f''/f|<1e3, |f'''/f|<1e5. Tangential derivatives commute
with the paired projector. Nested tangential energy estimates through
six, followed by the first-order equation for the four normal
derivatives, cost at most j!(1e8 K)^j. Only three normal coefficient
derivatives occur in this recovery. The fourth derivative of the
outside scalar density is bounded separately by 1e10 L(X).
The endpoint Leibniz sum through j=4 fits 1e50 K^j L(X).

Since K<600M for M>18, a common bound for all four observation
derivatives is consequently

    1e85 M^5 exp(-M/1000)L(X).                       (7)

This deliberately enlarges (6), including both endpoint derivative
distributions. It never differentiates the C2 even reflection across
the boundary; all normal derivatives are on its smooth positive side.
The complete counting estimate N(E)<=4CE^2 gives, without truncation,

    sum d_n M_n^5 exp(-delta M_n)
       <=4C[7!+5*6!] delta^-7, delta=1/1000.         (8)

One obtains (8) by Stieltjes integration by parts and bounding both
terms of the absolute derivative. No assumption that the summand is
decreasing down to M=0 is needed.

To bound the endpoint DIFFERENCE, use the fixed-endpoint metric chord
from the retained response proof. Its global C5 input is at most
ell 1e-4; its relative first-order graph coefficient is at most ell.
At every real chord point a complex parameter disk of radius
1e-4/ell therefore has graph-Neumann ratio at most .001, with the same
local coefficient and row bounds. Cauchy and the real chord integral
give an additional 1e4 ell. The case ell=0 is identically zero.
Multiplying (7),(8) by 2 eta_G/C gives an allowance below
1e85 ell L^2 for the actual-minus-cap difference through four jets.

Cauchy on (5) gives j!1e-11 rho^-j ell L^2. Adding the actual-domain
correction, then using X^j partial_X^j=D(D-1)...(D-j+1), yields

    |X^j partial_X^j Delta R_M|<1e154 ell L(X)^2,
    0<=j<=4.                                      (9)

The explicit cubic difference is at most a constant times ell/X and
has its exact ordinary derivatives. Restoring it fits
1e155 ell X^-1 L^2 for the complete massive scalar difference at the
same proper X. The fixed flat-vacuum scalar cancels only in this
difference. It remains in z_ell^4 S_M,ell.

## 4. Pullback, radius multipliers and every local sector

The forward map's low relative jets through six are bounded by
3e-10 j!100^j. Their finite composition through four costs less than
1000. Pull (9) back by X=Phi_ell(y), and retain
S_M,0(Phi_ell(y))-S_M,0(y) separately. The seed source derivatives
through five are supplied by the seed outer theorem. The mean-value
identity and the forward-map remainder give a bound below
1e80 ell y^3 L^2 for this regular pullback term; its explicit cubic
pullback is also retained. A common, generous fixed-y bound is

    |y^j partial_y^j(S_M,ell-S_M,0)|
       <1e160 ell y^-1 L^2, j<=4.                 (10)

The original radius satisfies z_ell<.43y. Its low divided scaled
jets give |y^j partial_y^j z_ell^4|<j!1000^j y^4. The endpoint
radius difference has the weight-three factor ell y^3 L, with
finite derivatives supplied by the outer metric theorem. Include
both terms in

    z_ell^4 S_M,ell-z0^4 S_M,0
      =z_ell^4(S_M,ell-S_M,0)+(z_ell^4-z0^4)S_M,0.

The latter term includes the fixed vacuum. The full product bound
through four is below 1e176 ell y^3 L^2. Multiplication by C/x<6.1
and the bookkeeping ell fits 1e177 ell^2 y^3 L^2.

For completeness the following finite local ledger handles all
remaining terms in (1). It uses the original common coordinate,
with lapse exp(ell u), orbit warping log a0+ell u and
phi=phi0+log[1+(1-z0)(exp(ell V)-1)]. Keep the physical potential
and comparison potential distinct. After multiplication by z^2,
their coefficients have the bounds already proved in the interior
argument; their powers are 0,2,4,5,6.

In divided scaled jets through eight, the real finite background
has coefficient reserve H^j, H=1e6, for positive orders, and
|phi|<2L. Retain the smaller low C5 bounds when verifying this
statement; the large high allowance is not used as a value. The
endpoint logarithm differences have reserve
ell A R^j y^3 L, A=1e20, R=1000. The original path's first and
second radius derivatives have reserves 1e25 H^j y^3 L and
1e50 H^j y^6 L^2. These follow by finite products in

    S'=tau V exp(ell V)/[1+tau(exp(ell V)-1)],
    S''=tau(1-tau)V^2 exp(ell V)/[1+tau(exp(ell V)-1)]^2,
    tau=1-z0.

All scalar denominators exceed .999999; their positive-order jets
use the finite metric arrays. The second-order Taylor remainder of
K thus has divided scaled four-jet reserve below

    1e12*6!*H^6*(2e50) <1e105,

times ell^2 y^3 L^2. The explicit ell in the interpolated potential
is included; its mixed term has only one field difference.

For the Weyl term use the exact proper-jet polynomials M_t+M_Y from
the retained outer recombination, with proper derivative
exp(-ell u) partial_y. Differentiate z^4(M_t+M_Y) in y through four,
then take its endpoint difference along the scalar segment of the
finite logarithms. Each differentiated monomial has total derivative
order at most eight. Its tagged variation has one factor
ell A y^3 L; the sole undifferentiated phi supplies at most one
additional L. The magnetic term qbar phi/3, qbar<160000, remains.
The checker expands these actual polynomials and their tagged
coefficient sums; after eta_G/pi and C/x they fit 1e110.
The complete LLL scalar eta_G q R_gamma/(12 pi C) is treated by
the same finite products through six and fits the same reserve.
Neither sector has been identified with the massive determinant.

Combining (1), (10) and these local differences proves the fixed-y
ordinary derivative bound

    |y^j partial_y^j O_ell| < B ell^2 y^3 L^2,
    B=1e180, j=0,...,4.                           (11)

All cancellations in (1) precede these unsigned bounds. No physical
source or signed compatibility condition is assigned zero separately.

## 5. Full physical tensors, outer volume and global assembly

Convert (11) to canonical product tensors using the raw connection
reserves (.3,4000,2e6), now with a factor y^-j at order j. The proper
chart is normalized tangentially at the observation; stationarity
sets its ordinary time derivatives to zero, while the connection
terms retain the time slots. The raw bounds also cover the seed
proper-chart connection jets. Since y<1, each connection derivative
can be enlarged by its corresponding inverse power of y. The same
ordered-tensor recurrence as the interior proof therefore applies
to the coefficient array (B,B,B,B,B).

Multiply order j by (z0/y)^j<.42^j and divide by j!. These are the
divided D=z0 nabla_gamma arrays. Apply the actual four-dimensional
connection recurrence from the refreshed finite geometry,

    H_(n+1,k)=(k+1)H_(n,k+1)+n sum_(i=0)^k C_i H_(n,k-i),
    P_n=2(2.5)^n H_(n,0), n<=4.

Every output retains ell^2 y^3 L^2. The connection recurrence
includes the angular slots generated by the physical radius. It is
not a radial-only replacement of H4.

The actual volume ratio is Jvol<1.001. In the common seed coordinate,

    dmu_ell/gamma0 <=1.001 z0'/z0^4 dy <36 y^-4 dy,
    36 integral_0^1 y^2(1-log y)^4 dy
       =36 sum_(k=0)^4 binom(4,k) k!/3^(k+1) <8^2.   (12)

The larger integration interval only enlarges the bound. There is
no inverse temperature factor. Let C_sigma4 be the retained global
scalar H4 coefficient. With pi>3, a valid outer shape constant is

    C_O=(1/6) sum_(n=0)^4 P_n+C_sigma4/6.

The interior and outer restrictions partition the original tensors;
no derivative of their indicators is taken. By the triangle inequality
for the two restricted L2 norms,

    F_b4 <= C_b ell^2/g, C_b=C_I+C_O <1e200,
    F4 <= (1e4+1e15 C_b) ell^2/g <1e215 ell^2/g.     (13)

The exact rational receipt supplies tighter upward decimal enclosures
of C_O and C_b. This closes the numerical source input only. The
higher TT inverse C_D4 and positive transfers R2,R4 remain unknown;
the useful Schur target has not passed. No controlled nonlinear
solution or stable semiclassical wormhole is inferred from (13).

## 6. Frozen verification and limits

Freeze this note, checker and review before execution. Verify the
cutoff, D-jet formulas, cap tube, scaled angular/restoration factors,
short-separation row and full moment budgets, fixed-coordinate
composition, all local polynomial terms, and full tensor/volume
assembly. Independently test finite dilation jets on a C10 field
that is not C11, the exact heat scaling with a nonconstant exponential
coefficient, and outer residual/volume integrals at 60 and 90 digits.
Diagnostics are manufactured examples, not samples of the physical
quantum state or a machine verification of continuum operator bounds.
Preserve every run, failure, source snapshot and inherited input.

Remaining to-do list: higher TT inverse and positive transfers;
sharpen the source/inverse combination to a useful preserved radius;
nonlinear mapping, cold/outer preservation, existence and signed
compatibility; physical matching, stability, charge/return, opening,
transport, preparation and creation.
