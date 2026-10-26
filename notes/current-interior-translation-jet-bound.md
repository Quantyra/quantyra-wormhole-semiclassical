# Finite translation jets for the complete interior quantum response

2026-10-25 HST, E003/S169. Ordinary argument before supporting arithmetic.
The full stable semiclassical wormhole objective is active and incomplete;
the closed audit and physical norms are unchanged.

This supplies the differentiated INTERIOR quantum ingredient of the triple-
residual calculation. It is not a global coupled inverse, numerical signed
solve, boundary-phase estimate, or nonlinear-neighborhood theorem. The
quantitative contribution is specified below: a loss of 13/4 derivatives,
strictly less than the fourth-order inverse gain, through output order ten.
The full original-space certificate remains unresolved.

## 1. Actual operator, localization and statement

Use the current finite seed and original paired Dirac domain, complementary
adjoint, fixed period, angular multiplicities and common renormalization.
Keep its one-sided mass/log geometry and regular bolt. No smooth reflection
of the mass cubic is made. On the retained analytic conformal interior
charts, center Y>=1/4, use radius 1/32 and gbar=e^(2sigma)|dx|^2. The
holomorphic seed bound |sigma|<.027 and the zero value/first jet at the
center are those of the retained normal-chart construction. Current finite
seed guards are inside its polynomial parameter range, as checked in
current-seed-full-hessian-realization.md, section 3. In particular M_1>18,
C/q<.0035 and N(E)<=4CE^2 still hold.

H_mass is the full fixed-area massive Hessian in gbar units, WITHOUT the
factor 2G/B^2. Its density term is included. Retain the exact frozen symbol
K of compact-massive-logarithmic-graph.md, including its common fourth
coefficient and entire angular tower. Put

 P_sigma u(x)=K(e^(-sigma(x))|D|)u(x),
 A=psi (H_mass-P_sigma) chi.

The LLL and four-dimensional local Weyl/radius entries remain separate
entries of the original coupled H. They are not set to zero or included
in the lower-order difference A by calling a fourth-order term a remainder.

chi is one on radius 1/512 and zero from 1/256; psi is one on 1/256 and
zero from 1/128. Both use the increasing beta(33,33) step for their radial
transition, with the appropriate orientation. They have C32 regularity.
All norms below are Euclidean H^s(R2), with unitary Fourier transform and
weight <xi>^s. Input is extended by chi on the actual physical chart; output
is extended by psi. For identical paired fields include the matching input
on the second copy and its cross return. The same proof works for the
matrix of these compactly supported inputs; at most two copies cost a
fixed factor below sixteen in the bounds below.

The claim is

 ||A u||H^s <= 10^4000 C ||u||H^(s+13/4), 0<=s<=10, s integer. (1)

It applies to the actual fixed-seed response, including the whole finite-
domain return. Auxiliary caps only enter a resolvent identity; their
difference from the actual realization is explicitly estimated. No finite
thermal or angular sample replaces the response. The statement does NOT
cover an independent mass/log variation at the seam, nor backgrounds with
only the admitted nonlinear regularity.

## 2. A quantitative lower-order starting estimate

The retained fixed-chart proof bounds the actual H_mass-Q_geo by 10^48 C
times three covariant input derivatives. Its coordinate conversion is
below 2*10^60 C. For Q_geo-P_sigma it keeps the full contact allocation:
the kernel difference is <=10^12 C/r^5; every remaining local contact has
order <=3; the ten finite probes bound those coefficients by 10^87 C.
The separated tails are bounded, including geodesics that leave and return.

Use the SAME allocation, replacing only the near Taylor remainder bound.
For a scalar Fourier variable t,

 |exp(it)-sum_(j=0)^3 (it)^j/j!| <=4 |t|^(13/4).

For |t|<=1 use the fourth-order integral remainder; for |t|>=1 sum the
four polynomial terms and the exponential. Translation and Plancherel then
give the corresponding L2 Taylor remainder <=4|h|^(13/4)||u||H^(13/4).
At separation R0=10^-20 its kernel integral is

 2pi*4*10^12 C int_0^R0 r^(-3/4) dr
       =32pi*10^12 C R0^(1/4).

The retained tail/contact budget is below 10^91 C in H3->L2. Thus the
uncut fixed-chart difference has H^(13/4)->L2 norm below 10^100 C.
Multiplication by chi is bounded in H^(13/4) by its C4 bounds. Integer
Leibniz and the fractional difference quotient give a factor <10^80;
psi has value <=1. We retain the deliberately larger bound

                 ||A||H^(13/4)->L2 <= A0=10^300 C.          (2)

This is a strong bound for the commonly restored DIFFERENCE, not separate
H^(13/4) bounds for its fourth/log terms. No contact is differentiated by
differentiating its unsigned coefficient estimate.

## 3. Finite jets, not a holomorphic moving cutoff

Use a plane auxiliary cap c, one through 3/160 and zero by 1/40, with
the same beta(33,33) transition; phi_cap=c sigma. It agrees with the actual
seed on the entire region needed for the localization and gluing below.

The step is a degree-65 polynomial on its transition. Its normalization
65!/(32!)^2 is <=65*2^64. On a complex radial neighborhood with scaled
transition variable within .001 of [0,1], its coefficient formula bounds
the step by 10^32. The radial square-root branch is analytic away from
the center; the transition annuli stay away from that center. A complex
Cartesian radius 10^-6 works for c; radius 10^-9 works for chi and psi.
The seed Cauchy margin on the cap is at least 1/160. Consequently, for
every real unit direction v and 0<=j<=17,

 ||partial_v^j phi_cap||infinity <10^150,
 ||partial_v^j chi||infinity, ||partial_v^j psi||infinity <10^201. (3)

More specifically the cutoff Cj bounds are j!10^32(10^9)^j, and the cap
ones are j!10^32(10^6)^j. These estimates use analytic polynomial branches
only on the transition. The matching derivatives through order 32 at the
joins justify the finite global bounds; no analytic compact cutoff exists
or is presumed.

For each of phi_cap, chi and psi, replace its translated family f(x+tv)
by the finite parameter polynomial

             f_z=sum_(j=0)^12 z^j partial_v^j f/j!.

On |z|<=rho=10^-250, (3) gives a C5 difference <10^-40. The true real
translations and this polynomial have the SAME twelve parameter jets.
All input functions are evaluated at real spatial coordinates. The
complex continuation is only in this finite parameter. This does not
reuse the disproved holomorphic moving-short-travel construction.

## 4. The capped response in a sufficiently strong input norm

The complete point-source Born proof in general-interior-dirac-response.md
is holomorphic under small complex C5 metric perturbations. Here observation
points lie within radius 1/128, so the seed is not separately recentered.
Its zero value/first jet at the chart center and its analytic radius give
the following bounds for the locally matched mass potential h at EACH
such observation point, throughout the parameter tube and the source
variation circle |w|<=1/(200U):

 |h|<=.06, |h|<=10r, |h_1|<=2r, |h_2|<=900r^2,
 |h_3|<=2*10^6 r^3, |h_4|<=10^8 r^4,
 |h-T1h|<=10^4r^2, |h-T2h|<=2*10^6r^3,
 |h-T3h|<=10^10r^4, |h-T4h|<=10^14r^5.           (4)

These follow from Cauchy within radius 1/64 about the observation point,
the real C5 Taylor formula for the perturbation, and division at r=1/128
outside that disk. For example the outer remainder budgets use, successively,
12*128, 12*128^2+900*128,
12*128^3+900*128^2+2*10^6*128, and
12*128^4+900*128^3+2*10^6*128^2+10^8*128.
All fit (4). The mass has Re m>.98M and |m|<1.02M.

Every one-through-four-Born Taylor remainder is a sum of products with at
most four h factors. Inflating the corresponding retained coefficient
bounds by 100 costs at most 100^4. The five-and-higher-Born tail uses five
linear factors and then the unchanged .06 geometric bound; changing 2r to
10r costs 5^5. Thus 100^5 is a sufficient common inflation of the retained
10^17 M^-3 remainder and 10^17 C fourth-contact reserves. The full angular
sum still uses sum d/M^3<2C/3. This proves a common-restored, uniformly
holomorphic nonconstant scalar-source bound below 10^28 C on this tube.

Do NOT drop the fixed-area density contact. At the real base the retained
constant-probe estimate |H_mass 1|<2*10^26 C and the scalar first variation
|D T_mass[1]|<6*10^22 C bound the source value by <2*10^26 C. The cap/actual
source difference is below 7*10^11 sum d M exp(-M/256), and is smaller
than 10^21 C. The uniform parameter bound just proved controls the change.
With the source circle, first Cauchy and H_mass u=D T_mass[u]+2u T_mass give
the conservative point bound

                 |H_cap u(x)| <=10^31 C U.                  (5)

U is the real-coordinate C5 input norm, including global capped low norms.
The flat vacuum value is retained in the density term; it is not assigned
zero. The common a2 restoration and all its variations remain in (5).
The uniform M^-3 majorant permits the parameter derivatives and angular
sum to pass through common regulator removal, uniformly on this input ball.

H8(R2) bounds the raw C5 norm by a constant below 10, directly by Fourier
Cauchy--Schwarz. Multiplication by chi_z costs <10^120 in H8 and <10^90
in the C5 input estimate. The exact flat multiplier satisfies
|K(e^-sigma p)|<=10 C(1+p^6) in this complex sector: use the positive
Feynman remainder bounded by its cubic term and sum d/M^4<4C^2/q, as well
as the retained common fourth coefficient. Fourier Cauchy--Schwarz with
H8 bounds the corresponding point action. It is applied at real x.

Consequently the finite parameter family

           psi_z (H_cap[phi_z]-P_sigma_z) chi_z

is holomorphic H8->L2 and bounded by 10^150 C on |z|<=rho. Its twelve
jets coincide with conjugation of the fixed plane operator by translations.
The COMPLETE capped twelfth commutator is therefore bounded by

               12! rho^-12 10^150 C <10^3160 C.             (6)

This is an operator-valued Cauchy estimate justified by the uniform source
majorants, not an inference from a real pointwise inequality alone.

## 5. Retain and differentiate the actual finite-domain return

The actual inverse stays on the original paired domain. Use a fixed gluing
cutoff equal to one through 1/64 and zero by 3/160, with commutator norm
<10^5. Actual and cap coefficients agree on its support. The observation
disk of radius 1/128, even enlarged by 1/4096, has separation >1/176 from
the gluing transition. The retained geometric resolvent identity expresses
E_M=S_actual-S_cap with two cap endpoint rows and an ACTUAL middle inverse.
Its operator norm is <=1/M (or <=3/M with retained weighted factors).
No differentiation of the middle inverse or global boundary data is made.

The endpoint row derivatives are bounded by the local homogeneous Dirac
equation. In the flat conformal frame it is (D_flat+Ma)v=0, a=e^sigma.
Direct integration by parts gives

 ||(D_flat+Ma)w||^2
 =||partial w||^2+M^2||aw||^2-M<w,c(da)w>.

Here |a|>=exp(-.027), |da|<2 and M>=18, so the right side is >=.8 times
the flat gradient-plus-M^2 norm. Cauchy on the fixed analytic seed gives
A1=2, A_j=2j!100^j for j>=2, through the required finite orders. For a
target j<=12, set m=j+3, L=8192m and B0=20000m. On nested balls of initial
radius R=1/4096 the cutoff energy recurrence is

 E0<=1.3 L F0,
 E_k<=1.3[L E_(k-1)+sum_(b=1)^k binom(k,b)A_b E_(k-b)].

The coefficient ratio is below one through k=j+1. Thus derivatives through
j+2 are bounded by B0^k F0. The scaled two-dimensional H2 point estimate
on radius R/(2m) gives the following safe endpoint row constants:

 c0=10^8,
 c_j=10^14(j+3)[20000(j+3)]^j, 1<=j<=12,
 ||row_j|| <= c_j exp(-M/352).                       (7)

The original point-row bound supplies F0 uniformly on the initial ball;
its L2 volume factor cancels the inverse radius in the point estimate.
All c_j are <10^90. The same estimate applies to Hilbert-valued rows and
the other physical copy. Spin and bundle identifications are unchanged.
Leibniz, both gluing cutoffs, density/vertex factors and (7) consequently
give, for total endpoint order <=12,

 |partial_x^a partial_y^b E_M(x,y)|
       <=10^300 M^-1 exp(-M/200).                   (8)

This generous bound includes all finite returns, including the cross-copy
kernel when inputs on the two identical physical pieces are linked.

For clarity the other factor in the Hessian difference is controlled as
well. The complex flat Dirac envelope can be enlarged to
k_M(r)=.6 exp(-.6Mr)/r. In the n-link convolution, one link is at least
r/n. Splitting off half of the common exponential, integrating the other
links, and summing the n choices bounds it by

 n^2 (.6/r) exp(-.3Mr) [4pi/M]^(n-1).

The insertion norm is <=.0612M, hence the ratio is <.8. Summing n^2 times
this geometric sequence gives

            |S_cap(z;x,y)|<=200 exp(-Mr/4)/r.        (9)

The endpoint density factors fit a further factor ten. The Bessel integral
underlying the envelope is DLMF 10.32.9, with Re m>0 as required:
https://dlmf.nist.gov/10.32.E9 . No theorem on a smooth global double is
used. Finite parameter Cauchy supplies joint endpoint derivatives through
twelve of (9), with factor j!rho^-j; this uses the fixed-cap translation
jets from section 3, not derivatives of an unsigned kernel bound.

At each finite angular sum the difference H_actual-H_cap consists of the
tadpole difference -M tr E_M(x,x) times the diagonal delta, and bubble
terms with S_cap E_M, E_M S_cap and E_M E_M, with the original vertex and
density factors. Local subtraction contacts cancel in this DIFFERENCE
because the coefficients agree on the observation/input region. They
remain in H_cap. No actual boundary contact is dropped: direct seam
variations vanish for these supported inputs, and nonlocal boundary effects
remain in E_M. Both actual thermal branches are in S_actual.

Equations (8)-(9), product differentiation and Schur integration bound
each joint derivative through twelve of the bubble by an integrable
kernel. The cross term has type
10^3320 M exp(-M/200) exp(-Mr/4)/r before the localization factors.
Its row integral is <=10^3330 exp(-M/200). The square E_M term and
tadpole are smaller. Cutoff derivatives cost <10^300, and the complete
angular estimate sum d exp(-M/200)<=8C*200^2 costs <10^6 C. Rank, paired
copy and finite combinatorial factors fit 10^50. Thus the FULL actual
localized twelfth commutator, including (6), obeys

       ||ad_(partial_v)^12 A||H8->L2 <= B12=10^4000 C.       (10)

The exponential majorants justify all differentiated return sums and the
common-regulator limit. The actual paired domain has not been transported
by Cartesian differentiation. Only compact interior kernel jets were
taken; no global stationary inverse is applied to a nonstationary field.

## 6. The remote-input map, including paired boundary inputs

Restrict observation now to D={|x|<1/1024}. The support of 1-chi has
proper separation >1/2048. All local coincidence contacts vanish on this
separated input, including their derivatives. The nonlocal kernel is still
the full actual tadpole/bubble realization; the separated tadpole vanishes.

Use the FAR one-sided resolvent identity from stationary-curved-quantum-
tt.md, section 4, retaining its adjoint identity for reverse rows. Its
baseline commutator row is 9.8304*10^8 exp(-M/530); the actual inverse
adds 1/M. Shrink its evaluation neighborhoods by 256 and allow the same
factor in the gluing gradient. The retained graph/point proof costs at
most the fourth power of this factor. The new separation and a cutoff
gradient <10^5 fit this enlargement. Therefore both separated L2 rows
are bounded by

                  10^20 M^-1 exp(-M/10^6).                 (13)

The reverse row only uses the original zeroth-order estimate at its
possibly rough seam endpoint. No high derivative at that endpoint is
assumed. For derivatives in the interior observation variable, use the
homogeneous local equation on balls of radius 1/65536. Nested cutoffs
and the same .8 coercivity give L=131072m and B0=10^6m, m=j+3. The
row derivative bound

 10^35(j+3)[10^6(j+3)]^j M^-1 exp(-M/10^6), j<=10,

is adequate and is <10^120 M^-1 exp(-M/10^6). For a fixed far endpoint,
apply the local L2 derivative estimate in the observation variable to
its original reverse L2 row; this gives the same bound. Thus BOTH Schur
rows of every differentiated bubble kernel are bounded. There is no
integration of a large point bound over the whole exterior.

Indeed a jth derivative is a sum of products of two differentiated
Dirac kernels. Cauchy--Schwarz in the integrated endpoint bounds a Schur
row by the product of their L2 rows. The two mass vertices cancel the
M^-2 in (13). Density derivatives through ten, rank/copy factors and
Leibniz fit 10^50. The entire angular sum obeys

            sum d exp(-2M/10^6)<=2*10^12 C.

Schur's test and the finite H10 derivative sum prove the conservative bound

 ||H_mass(1-chi)u||H10(D) <=10^500 C ||u||L2(gbar, whole pair). (14)

This is an actual remote-input estimate, not the finite-domain return of
the localized input alone. The kernel proof extends from the physical
paired core to this L2 input norm. It does not redefine the domain of
metric variations or change the paired boundary conditions.

For T=Vxi the ORIGINAL energy converts into this global input norm without
assuming a physical C5 embedding. Retain m0=4(3-sqrt(5))>3, the actual
York disk bound ||h_disk(xi)||physical<=35||xi||physical, and
ell>=.99 E_eff>=.99 F_l lambda ||xi||physical^2/(8piG). The physical
four-volume is 4pi A_r^2 B^2 z^-4 times gbar area, and z<=1. Thus

 ||Vxi||L2(gbar)^2
 <= [2450 N_halves pi/(m0 .99 F_l gamma)] E(xi)^2,
 E(xi)^2=ell(xi,xi)/(|W| beta), |W|=N_halves pi B^2/(G gamma).

With N_halves<=2, F_l>.99 and gamma>5*10^-34 this gives

                       ||Vxi||L2(gbar)<10^19 E(xi).          (15)

Combining (14)-(15) bounds the remote-input H10 output by 10^519 C E(xi).
This does not estimate independent phase-source energies I_j, the original
source embedding I_Y, or the nonlinear graph embedding alpha_embed.

## 7. Recover a lower-order strong estimate by finite interpolation

Let U_t be ordinary translations of the ZERO-EXTENDED localized operator
on R2. This is an auxiliary conjugation A_t=U_t A U_t^-1, not a physical
moving boundary family. Translations are isometries of every H^s. Hence
(2) and (10) bound A_t and its twelfth derivative for every real t.

Let Pi_n be orthogonal Fourier annuli <xi> in [2^n,2^(n+1)), n>=0,
with the low annulus enlarged to cover <xi><2. For f in that annulus set
g(t)=A_t f. Then

 ||g||infinity <= a_n||f||2, a_n=A0 2^((13/4)(n+1)),
 ||g^(12)||infinity <= b_n||f||2, b_n=B12 2^(8(n+1)).

At the twelve points kh, k=0,...,11, Taylor through order eleven has
remainder <=(kh)^12 b_n||f||2/12!. Choose h=(a_n/b_n)^(1/12).
The matrix k^j/j! has determinant one. Bounding its inverse by cofactors,
using entries <=11^11 and eleven-factor determinants, gives a row bound
<10^140. The Taylor remainder and twelve sampled values then imply

 ||ad_(partial_v)^j A Pi_n||2->2
 <=10^150 A0^(1-j/12) B12^(j/12)
               2^([13/4+(19/48)j](n+1)), 1<=j<=10.          (11)

No sharp interpolation theorem or unverified numerical constant is used.
Polarization handles mixed coordinate directions, at cost below 10^20
through order ten. In differentiating A u, a j-fold commutator is applied
to an input derivative of order s-j. Relative to H^(s+13/4), the remaining
dyadic factor is

                       2^(-(29/48)j n).

Its square sum is <2 for j>=1. Cauchy--Schwarz therefore sums ALL input
frequencies. The j=0 term is estimated directly by (2), without a dyadic
triangle inequality. Leibniz, polarization and the H^s norm conversion
cost less than 10^50 through s=10. The largest constant exponent in (11)
is 150+300+(3700)(10/12)<3534. The finite factors fit 10^4000, proving (1).

## 8. Quantitative role and limits in the complete feasibility test

For the last triple-residual bootstrap the desired local output is H14,
and the equation source is H10. The new quantum remainder needs only
H^(53/4)=H13.25 input, leaving an actual 3/4 derivative reserve. The earlier
point-source bound alone lost too many derivatives to supply this step.
Thus (1) closes the formerly missing differentiated INTERIOR quantum
remainder bound, with its local contacts and actual finite return included.

It does not bound C_g,3, C_X,3 or C_def,3 I_Y. The original physical areal
conversion, coupled scalar/York maps, outer mass/log weights, independent
phase inputs, cutoff commutators and global matching are still required.
For observation in D, section 6 ALSO controls H_mass(1-chi)u, including
its global input-norm comparison for u=Vxi. This is distinct from the
actual finite-domain return of the LOCALIZED input, bounded in section 5.
The constant in (1) is not one of those output constants and must not be
compared to their ceilings as if their spaces agreed.

For a quantitative sensitivity of the available reserve, a local principal
inverse factor 4000/C would multiply (1) into K=4*10^4003. This is an
EXPLICIT CONTRIBUTION in such a bootstrap, not a completed coupled bound.
Frequency splitting at L=(2K)^(4/3) gives, for s<=10,

 K ||u||H^(s+13/4)
 <=(1/2)||u||H^(s+4)+K L^(s+13/4)||u||L2.          (12)

Compute L and the s=10 low-frequency coefficient in logarithmic units.
Their size documents the severe loss of this generic proof. An oversized
certificate of this kind is not a lower bound on the actual response.
No original output constant is assigned the low-frequency coefficient.

Keep the polynomial candidate of current-schur-polynomial-return-test.md:
its form defect <7.351e-35 and affine factor <5.016e-43 are unchanged.
The original-source and physical output comparisons, actual kappa_m<=5e-5,
all joint affine/feedback entries and the independent boundary phases
remain unresolved. The quantum-fast .135973<.1599164 screen remains limited
to its two-fast admitted-family scope. No nonlinear-ball transfer is made.

DECISION: use finite translation jets and interpolation for the needed
interior regularity; retain the construction, three residual factors and
Chebyshev polynomial. The generic constants are coarse. Complete the global
composition before deciding whether this proof can certify its numerical
targets; do not mistake the local remainder bound for admission.

ONE NEXT CALCULATION (PROPOSED): compose this actual interior estimate with
the original coupled mass/log and bolt returns to bound the complete
triple-residual output/source map and its independent inputs. The signed
coupled columns and original-norm residual have not been computed here.

Verification contract: freeze the ordinary proof and checker before data;
check finite cutoff/jet majorants, source inflations, kernel convolution,
endpoint recurrence, finite interpolation and absorption arithmetic. Use
512-bit interval arithmetic, preserve failures and frozen sources, review
the ordinary argument, and run unchanged acceptance. These checks do not
replace the continuum argument or certify the uncomputed global map.

Remaining to-do list: global coupled output/source and independent phase
comparisons, signed solution, nonlinear admission/contraction, exact
compatibility, physical matching, semiclassical validity/noise, all-sector
stability, opening, transport, preparation, separate creation and the full
stable semiclassical wormhole objective.
