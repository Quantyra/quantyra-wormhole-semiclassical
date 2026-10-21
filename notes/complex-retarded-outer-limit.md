# Complex retarded outer limit and its scaled spectral measure

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

For each fixed compact profile pair, the actual seed's retarded outer
form converges locally uniformly throughout the cut squared-frequency
plane, and so do all fixed frequency derivatives. The same holds in
matrix norm on each fixed finite profile space. Its scaled positive
spectral measure has an explicit weighted weak limit. These statements
extend the [positive-axis transfer](retarded-outer-frequency-transfer.md).
They do not give an inverse on the physical metric space or a mode.

## 1. Actual analytic family with fixed contacts

Keep the seed, HH restriction, full angular tower, neutral spherical
sector, fixed projector and common subtraction of the preceding note.
Use its normalized canonical time tau and constant C. Fix real smooth
profiles p,v supported in a compact subinterval of (0,infinity), and set
p_d(Y)=d^-1/2 p(Y/d), v_d similarly. If D_d(sigma) is their full
fixed-density retarded form, define

    F_d(z)=d^4 D_d(z/d^2)/C,
    Omega=complex plane minus (-infinity,0].              (1)

The physical frequency is omega=i sqrt(z)/d, using the square root
with positive real part. Thus Omega is exactly the square of the open
right Laplace half-plane. This is the analyticity of the prescribed
matter response plus local metric terms, not of a gravitational inverse.

The [thermal spectral derivation](thermal-shape-response-and-subtractions.md)
gives a nonnegative measure mu_d depending only on v_d and the exact
third-subtracted response. Define its scaled measure by

    nu_d(B)=d^6 mu_d(B/d^2)/C.                            (2)

Every local term remains a polynomial of degree at most two in z.
The unsubtracted integral may diverge and is not used. For any three
distinct positive anchors a,b,c, the actual response satisfies

    F_d(z)-I_(a,b,c)F_d(z)
     =(z-a)(z-b)(z-c)
       integral dnu_d(t)/[(t+z)(t+a)(t+b)(t+c)].           (3)

The fixed-probe bound integral dnu_d/(1+t)^4<infinity and local
finiteness at zero justify (3) on Omega by compact domination. Indeed
(t+1)/|t+z| is bounded when z ranges over a compact subset of Omega.
The right side is holomorphic there. It agrees with the previously
derived positive-axis response and hence with the retarded analytic
function. All seagull, chart, volume and off-shell contacts are retained
in the actual three values. No polynomial is chosen from a discontinuity.

## 2. Four positive samples control complex compact sets

Fix a1=1, a2=2, a3=3, a4=4 in the dimensionless scaled variable z.
These are actual positive-imaginary responses; they need not themselves
be thermal frequencies. The preceding theorem already proves their
limits as d->0. Let P(z)=(z-1)(z-2)(z-3), and put

    G_d=F_d[1,2,3,4]
       =(F_d(4)-3F_d(3)+3F_d(2)-F_d(1))/6
       =integral dnu_d(t)/[(t+1)(t+2)(t+3)(t+4)] >=0.    (4)

The positive local measure, rather than the sign of the gravitational
form, proves this inequality. For a compact K in Omega, set

    M_K=sup_(z in K,t>=0) (t+4)/|t+z| < infinity.

With ell_j the three Lagrange weights at 1,2,3, (3) gives

    |F_d(z)| <= sum_(j=1..3)|ell_j(z)| |F_d(j)|
                            +|P(z)| M_K G_d, z in K.    (5)

For example if delta=dist(K,(-infinity,0])>0 and R=max_K |z|,
then M_K<=1+(4+R)/delta. A sharper sector bound is

    M_K<=max(1,4/r)/sin(theta/2)                         (6)

when r<=|z|<=R and |arg z|<=pi-theta, with r,theta>0.
It follows from |t+z|>=sin(theta/2)(t+|z|). This estimate does not
assume that the radial operators at different frequencies commute.

The four F_d(j) converge by the actual positive-axis outer theorem;
therefore G_d is bounded and (5) is a bound uniform for sufficiently
small d on every K. No escaping-support Lehmann constant is assumed.
For a fixed finite real profile basis, the polarized spectral measure
is a positive matrix measure. Its version of (4) is positive
semidefinite. Entry bounds follow from polarization, or directly from
the trace of that finite matrix measure. This supplies matrix local
boundedness without an infinite-dimensional norm claim.

## 3. Locally uniform causal limit and derivatives

Every sequence d_n->0 now gives a locally bounded holomorphic family
F_(d_n). Montel compactness gives convergent subsequences on compact
sets. The existing limit at every positive real z makes any two such
limits equal by the identity theorem. If the whole family failed to
converge on a compact set, a failing sequence would contain one of
these convergent subsequences, a contradiction. Consequently

    F_d(z) -> Q(z)[p,v] locally uniformly on Omega.      (7)

This is the standard Vitali argument, as set out in
[Boas, sections 2.3 and 2.7](https://haroldpboas.gitlab.io/courses/618-2016a/notes2016.pdf).
The task-specific input is the uniform bound (5); pointwise convergence
or analyticity alone would not suffice. For a compact K and a positive
rho with its closed rho-neighborhood K_rho contained in Omega,
Cauchy's formula gives, for each fixed integer m>=0,

    sup_K |partial_z^m(F_d-Q)|
       <= m! rho^-m sup_(K_rho)|F_d-Q| ->0.              (8)

For m=0 this is just inclusion. Polarization proves the same result
entrywise and hence in matrix norm on any fixed finite profile space.
Composition with z=s^2 gives local uniform convergence and derivatives
for every compact subset of Re s>0 in the scaled Laplace variable
s=d gamma. No assertion reaches the imaginary s axis or s=0.

In particular this does not settle a fixed small unscaled gamma:
then z=d^2 gamma^2 approaches the excluded endpoint. The preceding
separate positive-axis result at fixed gamma>=2h remains valid.
No numerical convergence rate or physical annular threshold follows
from compactness. Profile dimension and radial oscillations stay fixed.

## 4. Identification of the complete analytic normal form

The positive-axis identity uniquely fixes Q(z). It is the holomorphic
extension of the FULL coupled normal form, not an extension of an
inverse or just of a favorable block. In dmu=eta^-4 d eta set

    D_z=eta^2(partial_eta^2-z)-2eta partial_eta,
    A_z=D_z-4, B_z=eta^2(partial_eta^2-z)-12,
    e=60pi/g+11/3, g=G lambda, H_z=-partial_eta^2+z.

Then

    Q(z)[p,v]=<p,(D_z+e)A_z p>_mu/(20pi)
             +<p,(D_z+e)B_z v>_mu/(30pi)+F(z)[v],       (9)

where every pairing is extended complex BILINEARLY in the frequency
dependent factors. In particular squared norms at real z are replaced
by bilinear squares, not by modulus squares at complex z. Write w_z=H_z v.
The direct shape form is

    F(z)[v]=1/(30pi){L_v(z)+integral log(eta) w_z^2
                                  +(Lbar-31/30)integral w_z^2}
           +2/(15pi)integral eta^-2(v'^2+z v^2)
           -(6/g+11/(30pi))integral eta^-4 v^2.          (10)

For the even extension v_e on the full line, use the unitary Fourier
transform and beta_v(dxi)=|vhat_e(xi)|^2 dxi/2. Its total mass is
||v||_(L2(0,infinity))^2. The nonlocal term is

    L_v(z)=integral (xi^2+z)^2 Log(xi^2+z)/2 beta_v(dxi). (11)

Log is the principal branch. Smooth compact support away from zero
makes v_e smooth compact and its Fourier transform rapidly decreasing.
On compact subsets of Omega the integrand and every fixed derivative
are dominated by an integrable polynomial-log weight times beta_v.
Thus (9)-(11) are holomorphic and agree with the established real
normal form, proving their identification in (7).

The positive real-frequency curvature-free estimate is preserved in
its previous scope. Its sign is not asserted for complex z. No complex
Schur inverse, auxiliary-tail interchange or removal of D_z+e is used.

## 5. Weighted weak limit of the actual scaled spectral measure

Only (11) has a third z derivative. Direct differentiation gives

    partial_z^3 Q(z)[p,v]
       =1/(30pi) integral beta_v(dxi)/(xi^2+z).          (12)

This identifies a positive limiting spectral density

    rho_v(t)=1/(60pi) integral (t-xi^2)_+^2 beta_v(dxi). (13)

Indeed integral_0^infinity u^2/(u+A)^4 du=1/(3A) for
Re A>0, then by analyticity for A outside the negative real axis.
Consequently 6 integral rho_v(t)/(t+z)^4 dt equals (12).
The difference of the two subtracted analytic expressions is a
quadratic polynomial, which cancels in their interpolation remainder.
This verifies the representation without taking boundary values of
the actual seed response on its cut. The sign and factor 1/(60pi)
are also consistent with the upper-cut imaginary part of (11).

Let W(t)=(t+1)(t+2)(t+3)(t+4) and lambda_d=dnu_d/W. Its total mass
is G_d. On the one-point compactification [0,infinity], any sequence
d_n->0 therefore has a weakly convergent subsequence. By (3) and (7),
every subsequential limit lambda_* has transform

    integral (t+4)/(t+z) dlambda_*(t)
          =(Q(z)-I_(1,2,3)Q(z))/P(z).                  (14)

The ratio has value 1 at infinity; removable anchor values are
understood. The measure dlambda_Q=rho_v(t)dt/W(t) is finite and has
the same transform. Finiteness follows from
rho_v(t)<=t^2||v||^2/(60pi); W is strictly positive at zero.

These transforms determine a finite measure even on the compactified
interval. Put u=1/(t+4), mapping infinity to zero. The ratio becomes
1/[1+(z-4)u]. Its power series near z=4 determines all moments in u,
including the zeroth. Polynomials are dense in C([0,1/4]); hence two
finite measures with these moments agree. Therefore lambda_*=lambda_Q,
and the entire family satisfies

    dnu_d(t)/W(t) -> rho_v(t)dt/W(t) weakly on [0,infinity]. (15)

In particular no weighted mass escapes to infinity. This also gives
vague convergence of nu_d against compactly supported continuous
tests, and weighted tail tightness as d->0 for this fixed profile.
Polarization gives the corresponding statement for a fixed finite
matrix of profiles. It is not a uniform theorem over all profiles,
pointwise density convergence, a limiting-absorption estimate, or a
finite-temperature on-cut response formula. No rate is assigned.

## 6. Verification and consequence for the next decision

The checker tests the exact scaled measure factor, third divided
difference, complex interpolation, sector denominator and polynomial
contacts. Independent continuous positive matrix-measure quadrature
checks the complex bound, remainder and third derivative. Independent
quadrature of (13) checks its subtracted transform and (12), including
complex points with negative real part. Two contour resolutions test
Cauchy derivatives. Positive scaled atomic families test complex
convergence and weighted measure convergence. All are controlled
identity fixtures, not a computation of the actual exterior spectrum.

Two controls protect the scope. The entire functions exp(-n z)
converge on the positive axis but grow at a fixed negative-real-part
point in Omega; positivity/local boundedness cannot be skipped. Also
the first compression of [[z-1,1],[1,0]] vanishes at z=1 although
the full matrix has determinant -1. A finite profile determinant zero
therefore need not be a full metric pole, even in a polynomial example.

The actual seed now has a complete analytic OUTER CORE family on the
open causal half-plane. The next transfer must choose the physical
metric boundary/initial-data realization and establish suitable full
operator estimates with variable-curvature mixing and residual
contacts. Neutral current decoupling and Ward recovery are available;
the boundary solvability test still applies to a proposed fast-branch
exclusion. This theorem does not close S168.3-S168.4, S167 or S077.

Remaining to-do list: physical charged metric boundary and quantitative
retarded operator transfer or justified exclusion; S167 signed source,
cold/higher response, nonlinear solve and compatibility; exterior/state,
charge/return, opening, full-sector stability, transport, preparation,
separate creation and final audit.
