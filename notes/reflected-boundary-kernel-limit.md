# The reflected boundary correction has an ordinary cutoff limit

2026-10-21, E003/S168. Ordinary argument before focused verification.
After the complete constant row is allocated, the nonlocal correction
caused by reflecting the actual seed at a positive cut has an ordinary
finite cutoff limit on the zero-trace H2 shape class. This includes its
remaining corner contacts. The correction cannot supply the persistent
oscillation required by the coupled trace calculation.

The comparison here is with a specified SMOOTH CONTINUATION of the same
local seed, including a symmetrized reflected kernel. That representative
is not identified with the constant-mass normal model. The smooth local
comparison, its derivative contacts and all original local shape partners
remain to be combined before a full pure-shape limit or branch exclusion.

Keep the fixed charged seed, Euclidean period, paired domain, full angular
tower, mass chart and common subtraction. Reuse the exact
[even-mass representation](moving-cut-constant-quantum-row.md), its full
constant-row limit, and [global locality](moving-cut-global-quantum-locality.md).

## 1. Isolate the reflected mass without smoothing it away

In the normalized isothermal collar at cut d, let rho_d(x) be the smooth
continuation of exp(sigma_d(x)) through x=0. Thus rho_d(0)=1 and
b_d=rho_d'(0)=O(d). Choose one fixed even cutoff and sufficiently small
collar, so the capped functions

    rho_e(x)=rho_d(|x|), rho_s(x)=rho_d(x)

are interpreted as one outside that collar and both differ from one by
at most epsilon<=1/28. They agree for x>=0. The cutoff is incorporated
in the following bound; Taylor's odd part gives, on x<0,

    |rho_e(x)-rho_s(x)|<=2|b_d||x|+K3|x|^3.           (1)

All constants are uniform for small d. There is no quadratic difference.
The even mass remains Lipschitz at the cut; rho_s is only a comparison
coefficient. Let S_e=(Dflat+M rho_e)^-1 and S_s=(Dflat+M rho_s)^-1 on
the complete plane. Changes outside the local cap are the already
controlled global contribution. The actual thermal identification is
therefore retained in that contribution rather than replaced physically.

## 2. A uniform resolvent bound with a tagged odd jet

The retained flat Dirac majorant is

    k_M(r)=(sqrt(2)/pi) exp(-Mr/sqrt(2))/r.

It follows from the positive-real Bessel representations in
[DLMF 10.32.9-10](https://dlmf.nist.gov/10.32) and the previously retained
flat kernel estimate. At exponential weight exp(Mr/4), its L1 norm
is L/M, where L=2sqrt(2)/(1/sqrt(2)-1/4)<7, and its r-weighted
supremum is at most sqrt(2)/pi<1/2.

The convolution inequality

    ||r(f*g)||_infinity
       <=||r f||_infinity||g||_1+||f||_1||r g||_infinity

and r<=r1+r2 show that the n-th Born term has weighted r-supremum
at most (n+1)(sqrt(2)/pi)(epsilon L)^n. Summing the absolute series,
with epsilon L<=1/4, gives for BOTH capped kernels

    |S_e(p,q)|, |S_s(p,q)|<=exp(-M|p-q|/4)/|p-q|.    (2)

This is a scalar envelope for matrix operator norm; spin traces only
enlarge the later constant. The small cap is independent of M and d.

The exact resolvent difference is

    S_e-S_s=-S_s M(rho_e-rho_s)S_e.                    (3)

For p=(x,0), q=(x',tau), x,x'>=0, put

    r=|p-q|, R=sqrt((x+x')^2+tau^2).

Every path through z in the negative half-plane has r1+r2>=R.
With mu=M/4 and j=1,3, the tagged positive convolution satisfies

    I_j=integral_(z_x<0) |z_x|^j exp[-mu(r1+r2)]/(r1 r2) d^2z,
    I_1<=4pi mu^-1 exp(-mu R/2),
    I_3<=32pi mu^-3 exp(-mu R/2).                    (4)

Indeed |z_x|<=min(r1,r2); cancel one denominator for j=1 and use
|z_x|^3<=r1^2 r2 for j=3. Extract exp(-mu R/2), and integrate the
remaining radial envelopes over the whole plane.

For q reflected into the negative half-plane, R=|p-q|. The same I1
bound holds because |z_x|<=r1. For I3 use
r1^2/r2<=2r2+2R^2/r2 and absorb R^2 into a quarter of the exponential.
A common bound is then 256pi mu^-3 exp(-mu R/4). Equations (1)-(4)
give, for both the direct and crossing pairs,

    |(S_e-S_s)(p,q)|
       <=K (|b_d|+K3/M^2) exp(-M R/16),              (5)

with a finite fixed K, for example 20000pi after enlarging K3 if needed.
This is a boundary correction bound, not a bound by the direct distance r.
The distinction matters when the two positive-side points approach each
other while staying away from the boundary.

## 3. The two-point correction is a lower-order corner kernel

The paired image identity makes the half-domain bubble the sum of its
direct and reflected even-double kernels. Cross terms with one reflection
involution have zero spin trace. Replace S_e by S_s in the direct kernel;
in the reflected kernel use the symmetric average of the two crossings
(p,q*) and (q,p*), with the original positive-side mass vertices. This
defines the smooth-continuation representative used here.

Expand each bubble difference into S_s(S_e-S_s), its reverse and the
square of the difference. Endpoint mass factors are uniformly bounded.
The exact angular count N(E)<=4CE^2 implies, for example,

    sum d_n M_n^j exp(-a M_n R)<=K_j C R^(-j-2), j>=0.

The positive lower mass handles negative powers. After integrating the
time separation, the resulting symmetric radial correction E_d obeys

    |E_d(x,x')|<=K C [|b_d|/(x+x')^4+1/(x+x')^3]
                     [1+|log(|x-x'|/(x+x'))|]        (6)

on the fixed small collar, up to bounded global kernels already allocated.
One source of the logarithm is the direct cross term 1/(r R^4): split
the time integral at |tau|=x+x'. It is bounded by a constant times
(x+x')^-4[1+log((x+x')/|x-x'|)]. The cubic tag and square terms are
lower order and fit (6). The crossing term has no direct-distance
singularity and satisfies the same, looser bound.

In corner coordinates s=x+x', t=x/s, dx dx'=s ds dt, the angular
factor L(t)=1+|log|2t-1|| is integrable, with integral_0^1 L(t)dt=2.
Thus the strongest correction has radial measure s^-3 ds, one order
below the leading reflected normal kernel. Its diagonal logarithm is
locally integrable, not a new interior distribution. The same cap
resolvent identity gives coefficient convergence away from s=0 as d->0.

## 4. Finite parts and the remaining contact are controlled

Allocate the complete constant row as proved previously. For a static
radial shape, the correction now pairs with (v_d(x)-v_d(x'))^2. Write
r_t=2t-1 and fix a cutoff chi(s) equal to one near zero. The finite part
can be represented as

    Q_d[v]=integral E_d(st,s(1-t)) s
        { [v(d+st)-v(d+s(1-t))]^2
                    -s^2 r_t^2 v'(d)^2 chi(s) } ds dt
                 +ell_d v'(d)^2.                    (7)

The normalization and sign of E_d can include the original -1/2 factor.
No new subtraction is imposed on the physical action. Equation (7) is
a representation of the difference in its COMMON scheme. Scaling degree
at most four in the two corner variables permits contacts of total order
at most two. A value or first derivative of the difference square is zero
at the corner; second derivatives are multiples of v'(d)^2. Thus all
remaining corner choices enter only through ell_d.

Here ell_d is bounded and converges for the fixed scheme. This is not
inferred from an off-diagonal envelope alone. Its highest third-jet
boundary normalization agrees between the two representatives because
rho_e(0)=rho_s(0)=1; it is the already fixed flat principal contact.
Each different background jet has positive derivative weight and lowers
the remaining shape-contact order to at most two, as in the mixed heat
coefficient. Evaluate the residual coefficient on one fixed smooth probe
v(x)=x chi0(x). The actual half-domain metric and probe depend smoothly
on d, with fixed paired boundary projector after pullback. They form a
compact smooth family; the finite probe has no rough-jet limit problem.

For precision about this last continuity assertion, construct the heat
parametrix on the SMOOTH HALF-domain with that fixed probe and a small
variation parameter. The Gaussian recursion through an order beyond a4
uses uniformly bounded coefficient jets and uniform ellipticity on this
compact parameter set. The Volterra remainder has a common integrable
bound after the fixed four-dimensional subtraction. A small complex
parameter circle and Cauchy's formula give the same bound for the second
variation. The massive large-time part is controlled by the retained gap
and global locality. The finite probe Hessian therefore converges.
The smooth-reference finite part converges by the same fixed-scale
Taylor subtraction. Extracting ell_d from their difference proves its
continuity. Only this finite smooth family is differentiated.

The heat expansion and common zeta subtraction are the standard ones in
[Vassilevich (2.21), (2.28)-(2.32), and (5.33)](https://arxiv.org/html/hep-th/0306138v3).
The family argument above supplies the uniformity needed here; it does
not invoke a smooth heat expansion on the Lipschitz double or a C2
theorem for arbitrary finite-amplitude energy-class metric families.

## 5. An ordinary H2 corner-limit theorem

Let v in H2(0,R0), v(0)=v'(0)=0, initially smooth at positive arguments.
Put f_d=v'(d), e(a)=||v''||_(0,a). The elementary trace bounds are

    |v'(y)|<=sqrt(y)e(y), |v(y)|<=y^(3/2)e(y)/sqrt(3).

For 0<=t<=1 and s>0, Taylor's integral remainder gives

    v(d+st)-v(d+s(1-t))=s r_t f_d+e_s,
    |e_s|<=(2/3)s^(3/2)||v''||_(d,d+s).

Consequently the braces in (7) near s=0 have absolute value at most

    (4/3)s^(5/2)|f_d| e(d+s)+(4/9)s^3 e(d+s)^2.      (8)

This is integrable against s^-3 ds. On 0<s<d, its integral is bounded
by a fixed multiple of d e(2d)^2. On d<s<a, the unsubtracted difference
square is bounded by a fixed multiple of s^3 e(3a)^2, while the
subtracted term costs at most

    f_d^2 log(a/d)<=d e(d)^2 log(a/d)<=a e(a)^2/e.

After the integrable t factor, these estimates give a common small-collar
bound K a e(3a)^2 for all d<a, with harmless adjustments at the fixed
outer cutoff. The s^-2 lower-order part is easier. Also ell_d f_d^2->0.
Away from the corner, ordinary dominated convergence applies. Therefore

    lim_(d->0) Q_d[v]
       =integral E_0(x,x')[v(x)-v(x')]^2 dx dx', finite. (9)

The right integral is absolutely convergent at the corner. This is an
ORDINARY cutoff limit for this correction, without a pointwise v'' trace
and without the extra logarithmic form-domain condition needed by the
principal normal kernel. No cutoff averaging is used in this theorem.

The physical-to-isothermal maps satisfy F_d(0)=d, uniformly bounded C2
jets and positive first derivative, and F_d(x)<=K(d+x). The chain rule
gives the same local H2 tail estimates for v(F_d(x)); the extra v' term
is lower order. Thus (9) applies in the actual collar coordinates.

For a fixed temporal frequency, the correction with cos(gamma tau)-1
has two extra time-separation powers. The same estimates lower its
radial corner degree from four to at most two. Its product pairing
v(d+x)v(d+x') needs only the constant Taylor subtraction, whose contact
is proportional to v(d)^2. The corresponding first-order remainder and
the same trace bounds give an ordinary limit. Finite sums of temporal
frequencies and polarization follow. This is a fixed Euclidean-frequency
statement, not continuation to a causal physical disturbance.

## 6. Remaining comparison

The reflected boundary correction, its corner contacts and the full
constant row cannot provide the missing persistent vector. The smooth
continuation still has the seed's curved bulk principal structure. Its
comparison with the normal representative must retain density factors,
fixed derivative contacts and all Weyl, LLL and classical shape partners.
The present result does not set those remaining terms to zero or prove
the full pure-shape cutoff limit. No physical branch exclusion, nonlinear
solution or stability is claimed.

Verification checks the positive convolution bounds and their tagged
moments, the logarithmic corner factor, finite-part kernels against
independent integral continuation, and the H2 remainder estimate including
concentrated smooth profiles. These are diagnostics for the ordinary
proof, not evaluations of the full charged two-point kernel.
The first frozen run failed in angular quadrature: finite-precision
rounding at pi/2 made cos(theta) negative and sent the incomplete-gamma
formula outside its real domain. The exact coordinate cos(theta)=1-w^2
keeps the integration domain real and removes the endpoint singularity.
No complex value is projected to the real axis and no tolerance is
relaxed. The failed sources, exception and raw run are preserved.

Remaining to-do list: smooth local quantum-to-normal comparison, density
and derivative contacts with all local shape partners; physical boundary
phase and source compatibility; continuous causal constraints, exact
nonlinear existence and all-sector stability; exterior/state, charge/return,
opening, transport, preparation, separate creation and full audit.
