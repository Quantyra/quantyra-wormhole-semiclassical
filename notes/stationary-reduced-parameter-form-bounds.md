# Reduced parameter forms and the combined outer logarithmic energy

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary proof before
focused verification. **PROGRESS; complete nonlinear shape estimate open.**

The actual common-renormalized a2 contact has a variable-parameter bound
of 1600 relative to the reduced TT form, in the low-jet norm specified
below. Explicit multiplication and bilinear estimates also control
endpoint-vanishing coefficients in the combined outer logarithmic
energy. Those latter estimates concern a defined scalar comparison;
the actual Dirac three-vertex and remaining contact estimates are still
required.

This accompanies the [coupled cold response](stationary-reduced-cold-response.md)
and the [covariant reduced inverse](stationary-covariant-reduced-inverse.md).
Use the same analytic seed, seed York/return chart, prepared period,
magnetic flux, physical exterior, complementary paired Dirac domain,
angular masses and common scale mu0. The radial scalar simplification
is already proved in the reduced-inverse note. No nonlinear branch,
physical error bound or original signed-source decision follows here.

## 1. The actual a2 contact and its third variation

Normalize the canonical orbit metric by gbar=ghat/ell_o^2, with
ell_o=B/sqrt(C), and write Delta for its fixed seed Laplacian. Define
J0[T]=integral_gbar |T|^2 and J2[T]=integral_gbar |Delta T|^2.
The [retained common heat allocation](stationary-seed-quantum-form.md)
has

    alpha=-C F_eff/(960pi)>0,
    F_eff=4gamma_E-2H_q+4log(mu0 B), |F_eff|<100,
    Gamma_a2(v)=alpha integral_gbar exp(-2v) P_v^2,
    P_v=Rbar-2Delta v, c2=8alpha=-C F_eff/(120pi).       (1)

This follows from the exact two-dimensional curvature transformation
under gbar_v=exp(2v)gbar. It is the fixed four-dimensional common
subtraction, with the sphere radius B retained. The integrated Delta Rbar
allocation is treated in section 3. Complex parameters mean holomorphic
coefficients at real spatial points.

Start with smooth bolt-regular parameters and compact outer support.
Use a pairing antilinear in its first entry. The sesquilinear
polarization of the real Hessian in the same affine disk shape coordinate
is

    q_v[U,T]=c2 integral exp(-2v) [conj(Delta U)Delta T
       +P_v(conj(U)Delta T+conj(Delta U)T)
       +(P_v^2/2)conj(U)T].                            (2)

Indeed the second derivative in a real T of the integrand in (1) is
exp(-2v)[8(Delta T)^2+16P_v T Delta T+4P_v^2 T^2].
Thus (2) includes its density, curvature and seagull terms. At v=0,
integration by parts gives the accepted curvature-gradient and Delta Rbar
zero contacts in the seed formula.

For ||v||infinity+||Delta v||infinity<=1/100 and |Rbar|<=2, differentiate
(2) in a complex parameter direction w. On the vector (Delta T,T), the
coefficient matrix is

    c2 exp(-2v) [ -2w                 -2wP_v-2Delta w
                  -2wP_v-2Delta w   -wP_v^2-2P_v Delta w ]. (3)

It is complex symmetric, not necessarily Hermitian. Put
a=||w||infinity and b=||Delta w||infinity. Since |P_v|<=101/50 and
|exp(-2v)|<50/49, the absolute row sums before the exponential are
bounded by (151/25)a+2b and (20301/2500)a+(151/25)b.
The column sums have the same bounds. Their maximum times 50/49 is
at most (20301/2450)(a+b)<9(a+b). The finite matrix Schur inequality
and spatial Cauchy--Schwarz prove directly

    |partial_w q_v[U,T]|
      <=9c2(a+b) sqrt(J2[U]+J0[U]) sqrt(J2[T]+J0[T]).   (4)

There is no further polarization factor. Integration on the straight
segment from 0 to v gives the same estimate for q_v-q_0 with a+b
replaced by ||v||infinity+||Delta v||infinity. This also proves
holomorphy of this contact as a bounded form in the stated fixed norms.

## 2. Relative bound on the actual reduced TT form domain

Let c be the scalar TT amplitude, (S_c,T_c) its fixed seed return,
and Dr the positive selfadjoint reduced TT operator from the
[reduced inverse](stationary-covariant-reduced-inverse.md).
Write Q_a2,v[c1,c2]=q_v[T_c1,T_c2]. No moving York map is differentiated.
The independent Weyl functional still depends on S_c as well.

The classical energy is the previously fixed physical quantity

    E0=beta/(2G) [3D/2+lambda C0],
    D=integral r^2 F|c_r+2c/r|^2 dr, C0=integral r^2 |c|^2 dr.

Section 10 of the [actual curved TT proof](stationary-curved-quantum-tt.md)
retains C/(100pi) sum_j J2[chi_j T_c] after the fourth-order comparison
and localization costs. Its total disk lower cost is below .004E0.
Reserve also the full off-shell gauge/exponential-chart cost below
1e-6 E0. The classical contribution is at least E0; the reduced Gram
subtraction costs less than 2G lambda E0. Consequently on the genuine
compact-TT core,

    Dr[c]>=(.995-2G lambda)E0[c]
                    +C/(100pi) sum_j J2[chi_j T_c].    (5)

This uses the retained fourth-order margin, in addition to Dr>=.989E0.
To recover the GLOBAL fourth norm, use the same partition
chi_o=cos(theta), chi_i=sin(theta), whose derivative is supported in
.01<=Y<=.02 and satisfies |theta'|<500. The exact identity is

    Delta T=sum_j chi_j Delta(chi_j T)+|grad theta|^2 T.

It gives C J2[T]<=2C sum_j J2[chi_j T]+2(500)^4 C J0,outer[T].
The accepted physical-weight bounds are
C J0,outer/E0<3.4e-37 and C J0/E0<5.2e-29. The extra localization
term is therefore below 4.25e-26 E0. Combining with (5) gives

    C(J2[T_c]+J0[T_c])<(200pi+1)Dr[c].                 (6)

This also proves the required domain extension. The
[critical graph-density proof](stationary-weighted-tt-inverse.md)
makes the exact compact-TT core an operator core for D. The bounded
reduced Gram subtraction preserves that operator domain and core,
which is a spectral form core for Dr. Apply (6) to differences of
core approximants. Completeness and the distributional return map
identify their limits as (Delta T_c,T_c). Thus (6) holds on the actual
Dr spectral form domain. No equality with the older norm storing
separate positive-model and absolute Weyl-log terms is asserted.

Equations (4),(6) prove the main actual parameter bound:

    ||Dr^(-1/2) (partial_w Q_a2,v) Dr^(-1/2)||
        <1600 (||w||infinity+||Delta w||infinity).      (7)

The displayed ledger gives less than 1503: 9c2/C<=15/(2pi), and
[15/(2pi)](200pi+1)<1503 using pi>3. The rounded 1600 covers this
a2 contact only. A parameter change of low-jet norm 1e-6 changes it
relatively by less than .0016 along the stated parameter ball.

## 3. Outer traces and the integrated Delta Rbar term

For the accepted smooth conormal mass/fast representatives, the seed
return has (S_c,T_c)=(3/4,-3/2)z^3+X_(7/2). The radial TT homogeneous
z^4 datum remains in the fast remainder. Since Y and z are smoothly
comparable, v=delta V and S are O(Y^3) with the stated conormal
derivatives; O(Y^2) is a weaker valid envelope. Their value and first
normal trace vanish. The original projectors and mixed-square Robin
datum therefore retain their fixed realization.

For the integrated Delta Rbar allocation, the conformal factors in
boundary measure and unit normal cancel. Its boundary integrand is
the seed-normal derivative of exp(-2v)P_v. Its second variation is
the normal derivative of

    exp(-2v)[4UT P_v+4U Delta T+4T Delta U].             (8)

Differentiate (8) once more in w for the third variation. The selected
mass/fast differentiated powers make both cutoff fluxes tend to zero.
For example U,T=O(Y^3), Delta U,Delta T=O(Y) make the mixed products
O(Y^4) before the final normal derivative. The same conclusion holds
with the weaker differentiated quadratic envelope. The smooth bolt is
an interior point. Compact-core identities extend by (4),(6), without
an H1 cutoff approximation to a nonzero bolt value.

The FIRST flux variation is retained: a mass cubic can give a nonzero
curvature-normal-derivative source. Vanishing of the second and third
fluxes does not remove that source or establish arbitrary rough traces.

## 4. Actual resolvent holomorphy and an auxiliary relative-log lemma

For the fixed-domain paired anti-selfadjoint first-order Acal_0 and
a bounded C1 conformal multiplier v, the actual mass-potential identity is

    (Acal_v+M)^(-1)
      =exp(v/2)(Acal_0+M exp(v))^(-1)exp(v/2).          (9)

If ||v||infinity<=rho<log2, the middle Neumann ratio is at most
exp(rho)-1 and its norm is at most [M(2-exp(rho))]^(-1).
Bounded C1 multipliers preserve the graph domain. This proves parameter
holomorphy with the unchanged projector and period. It does not bound
a common-subtracted angular trace by unsigned summation.

For later scalar comparisons, let A0>=c>0 be selfadjoint and define the
sectorial form A=A0^(1/2)(I+B)A0^(1/2), ||B||<=delta<1, on its common
form domain. B need not be selfadjoint. Then the principal logarithms
have equal domains and

    ||log A-log A0||<=delta/(1-delta).                 (10)

Here is a proof that does not assume a common operator core. Set
R0(t)=(A0+t)^(-1), Kt=A0(A0+t)^(-1). The form resolvent equation gives
A0^(1/2)(A+t)^(-1)=(I+Kt B)^(-1)A0^(1/2)R0(t), and similarly for
the adjoint. The spectral identity
integral_0^infinity ||A0^(1/2)R0(t)f||^2 dt=||f||^2 and the weak formula

    log A-log A0=integral_0^infinity
        (A+t)^(-1)A0^(1/2)B A0^(1/2)R0(t) dt

give (10) by Cauchy--Schwarz in t. The common form domain belongs to
both log domains: for its f and t>=1, the form resolvent remainder
(A+t)^(-1)f-f/t is O(t^(-3/2)). The operator domain of A lies in this
form domain and is a core for log A by sectorial resolvent
regularization; the same holds for A0. Approximation and closedness
in both directions give equal log domains and the bounded difference.
Taking an operator norm before the t integral would lose this proof.

The remaining logarithmic statements use the explicitly AUXILIARY
full-line scalar operator

    Q=partial_Y-h'/2, H0=Q*Q=-partial_Y^2+V, V real, |V|<=22,
    A=H0+24>=24, K[g]=<g,log(A)g>.                     (11)

The measure and half-density are fixed. For m=exp(-v), H_v=m H0 m,
A_v=H_v+24, put alpha_v=||m-1||infinity+||m'||infinity/sqrt(24).
Leibniz Q(nf)=nQf+n'f bounds the relative form perturbation by
2alpha_v+alpha_v^2, also for complex m using the conjugate left factor.
Thus ||v||<=1e-6, ||v'||<=1e-3 give alpha_v<.000206,
delta<.000413<.001 and ||log A_v-log A||<.001002.

Since V+24>=2, the scalar heat kernel is bounded by
exp(-2t)(4pi t)^(-1/2)exp(-(Y-Y')^2/(4t)). Integrating dt/t gives,
for bounded Lipschitz n,

    ||[log A,n]||<=sqrt(2)Lip(n).                      (12)

Indeed its absolute kernel is bounded by
Lip(n)exp(-sqrt(2)|Y-Y'|), whose row and column integrals give (12).
These auxiliary statements do not identify the actual Dirac bubble.

## 5. The closed combined Hardy-log form

Use fhat(xi)=(2pi)^(-1/2)integral exp(-i xi Y)f(Y)dY, D=-i partial_Y,
and first take smooth compact f away from Y=0. Define

    D_H=log2+psi(1/4)=-gamma_E-2log2-pi/2<0,
    L[f]=<f,log|D|f>+integral log|Y| |f|^2-D_H||f||^2.

Its exact nonnegative representation is

    L[f]=(1/4) integral integral
      |f(Y)/omega(Y)-f(Y')/omega(Y')|^2
             omega(Y)omega(Y')/|Y-Y'| dY dY',
    omega(Y)=|Y|^(-1/2).                              (13)

The primary input is the fractional ground-state identity in
[Frank--Lieb--Seiringer, Hardy-Lieb-Thirring inequalities for fractional
Schrodinger operators](https://arxiv.org/pdf/math/0610593), equations
(2.2), (3.1)-(3.2), and Proposition 4.1. In dimension one its constants are

    C_s=2^(2s) Gamma((1+2s)/4)^2/Gamma((1-2s)/4)^2,
    a_s=2^(2s-1)pi^(-1/2)Gamma((1+2s)/2)/|Gamma(-s)|,
    omega_s(Y)=|Y|^(-(1-2s)/2), 0<s<1/2.

Divide that identity by 2s and send s down to zero. The constants tend
to (C_s-1)/(2s)->D_H and a_s/(2s)->1/4. The limit is dominated:
on the supporting annulus the difference supplies two diagonal powers;
a free variable near zero has weight at worst |Y|^(-1/2), and a large
free variable decays at least |Y|^(-3/2). Both variables outside the
support give zero. On the Fourier side the local logarithm and Schwartz
tail are integrable. This proves (13), including its complex polarization.

The nonnegative integral is closable: for an L2-null, form-Cauchy sequence,
take an almost-everywhere-null subsequence and apply Fatou to its
differences. Let E be its form closure with ||f||_E^2=L[f]+||f||^2.
The two logarithmic terms need not separately be finite in E.
No equality with an arbitrary Sobolev completion is used. The whole
even line gives the reflected scalar subspace; no spin or physical-half
factor is included in (13).

For complex Lipschitz n supported in [-R,R], M=||n||infinity and
L_n=Lip(n), (13) gives the useful finite multiplier bound

    L[nf]<=2M^2 L[f]+50(M^2+R^2 L_n^2)||f||^2.         (14)

To check the constant, the product-rule error is at most
(1/2)integral |f(y)|^2 K_n(y)dy, where

    K_n(y)=integral sqrt(|y|/|x|)|n(x)-n(y)|^2/|x-y| dx.

For |y|<=2R, split at |x|=4R. The near part is at most
48sqrt(2)R^2L_n^2 and the far part at most 16sqrt(2)M^2.
For |y|>2R only |x|<=R contributes, bounded by 4sqrt(2)M^2.
Thus K_n<=100(M^2+R^2L_n^2), proving (14). Closure extends it
to E. For even restriction require n even; the full-line estimate
itself allows arbitrary complex n.

## 6. Endpoint-vanishing multiplication into the positive log form

Keep (11),(13), let 0<R<=1/2, and suppose n is complex Lipschitz,
supported in [-R,R], with

    M=||n||infinity, L_n=Lip(n),
    W_n=ess sup_(0<|Y|<=R) |n(Y)|^2 log(1/|Y|)<infinity.

Then multiplication maps E to Dom((log A)^(1/2)), with

    K[nf]<=4M^2 L[f]
              +[116M^2+100R^2L_n^2+2W_n]||f||^2.      (15)

The low-frequency step uses support, not a bounded comparison at xi=0.
For g supported in [-R,R], |ghat|^2<=R||g||^2/pi and
integral_R log(1+a^2/xi^2)dxi=2pi a. Hence

    0<=<g,[log(D^2+24)-2log|D|]g>
                           <=2R sqrt(24)||g||^2.

The perturbation V of Aflat=D^2+24 has relative form norm at most
11/12. Equation (10) gives ||log A-log Aflat||<=11, with equal log
domains and therefore equal log-form domains. On the initial core,

    K[g]<=2L[g]+2 integral log(1/|Y|)|g|^2
                                  +(2R sqrt(24)+11)||g||^2.

Insert g=nf, (14), 2R sqrt(24)<5 and the definition of W_n to obtain
(15). For an E-Cauchy sequence f_j, (15) makes nf_j Cauchy in the
closed positive log form and nf_j->nf in L2. This proves the full
extension without subtracting infinite logarithmic quantities.

For |n(Y)|<=kappa Y^2 and L_n<=kappa_1 R, monotonicity of
Y^4 log(1/Y) through R<=1/2 yields the concrete version

    K[nf]<=4kappa^2 R^4 L[f]
      +R^4[(116+2log(1/R))kappa^2+100kappa_1^2]||f||^2. (16)

For example take n=epsilon Y^2 chi(Y/R), with even chi equal to one
through |Y|=R/2 and zero after R, and quintic transition
1-10t^3+15t^4-6t^5, t=2|Y|/R-1. Then Lip(chi)<=15/(4R),
L_n<6|epsilon|R, and

    K[nf]<=|epsilon|^2 R^4
                    [4L[f]+(3716+2log(1/R))||f||^2].

Thus the E-to-K multiplication norm is less than .62|epsilon| at
R=.1 and .0248|epsilon| at R=.02, using log10<3 and log50<4.

If only f is supported in [-r,r], r<=1/4, and n is known on [-2r,2r]
with |n|<=kappa Y^2, |n'|<=kappa_1|Y|, multiply n by the even quintic
cutoff equal to one through r and zero after 2r. With R=2r the literal
constants for (15) are

    M<=4kappa r^2, L_n<=(2kappa_1+15kappa/2)r,
    W_n<=16kappa^2 r^4 log(1/(2r)).                    (17)

This changes only the auxiliary multiplier. It does not alter the
physical state. Differences of coefficients obey (15) with their own
M,L_n,W_n, so this is also a Lipschitz multiplication estimate in the
coefficient norm squared 116M^2+100R^2L_n^2+2W_n.

## 7. A bilinear variable-coefficient log form

For 0<R<=1/4, let n be supported in [-R,R] with
|n(Y)|<=kappa Y^2 and |n'(Y)|<=kappa_1|Y|. Choose a real even cutoff
chi equal to one through R, zero after 2R, Lip(chi)<=15/(8R), and put

    b(Y)=Y chi(Y), c(Y)=n(Y)/Y, c(0)=0, n=bc,
    ||b||<=2R, Lip(b)<5,
    ||c||<=kappa R, Lip(c)<=kappa+kappa_1.

The quotient extends Lipschitz across zero by the quadratic bound.
Equation (15), applied to these LINEARLY vanishing factors, gives
K[conj(b)g]<=B_R^2||g||_E^2 and K[cf]<=C_R^2||f||_E^2, where

    B_R=R sqrt(10464+8log(1/(2R))),
    C_R=R sqrt((116+2log(1/R))kappa^2
                                  +100(kappa+kappa_1)^2).

Here W_b<=(2R)^2log(1/(2R)) and
W_c<=kappa^2 R^2log(1/R); Y^2log(1/Y) increases through 1/2.
On the smooth core the exact complex-pairing identity is

    <g,n log(A)f>
      =<log(A)^(1/2)(conj(b)g),log(A)^(1/2)(cf)>
                               -<conj(b)g,[log(A),c]f>.

Cauchy--Schwarz and (12) therefore prove

    |<g,n log(A)f>|
      <=[B_R C_R+2sqrt(2)R(kappa+kappa_1)]||g||_E||f||_E. (18)

Both positive-log factors extend by (15); the commutator is L2-bounded.
This proves the unique continuous E-to-E-dual FORM extension, without
claiming a bounded L2 operator or its original maximal operator domain.
If V,n,f,g are even, b and c are odd. The intermediate logarithm in
(18) must act on the full line, including its odd invariant sector.
This auxiliary proof introduces no odd physical Dirac mode.

When kappa_1<=6kappa, at R=.1 or .02 one has B_R<103R and
C_R<71kappa R, by log5<2, log10<3, log25<4, log50<4. The coefficient
in (18) is less than 7313kappa R^2+20kappa R, hence less than
76kappa at .1 and 3.33kappa at .02.

Endpoint vanishing matters for this isolated form. For an L2-unit smooth
even annular f, set f_delta(Y)=delta^(-1/2)f(Y/delta). Its E norm is
constant by cancellation of the two dilation logarithms, whereas
<f_delta,log(D^2+24)f_delta>=2log(1/delta)+O(1). The Fourier logarithm
is locally integrable and the tail is Schwartz, so the remainder is
bounded. Equation (10) changes this by at most 11 for A. If n=1 near
zero, <f_delta,n log(A)f_delta> therefore diverges. This does not decide
a full physical variation whose Weyl allocation may add correlated terms.

## 8. A constant-complex-background flat majorant

One additional bound concerns the positive flat UV remainder alone, in
the rank-two orbit normalization. Let x=a(1-a)k^2/M^2 and

    f(x)=log(1+x)-x+x^2/2=x^3 integral_0^1 t^2/(1+tx)dt,
    P_M(k)=M^2/(2pi) integral_0^1 f(a(1-a)k^2/M^2)da.

Its leading term is k^6/(840pi M^4), and sum d_n/M_n^4 converges.
For a CONSTANT complex background |v|<=rho=.01, the exact expression
is P_M,v=exp(2v)M^2/(2pi)integral f(exp(-2v)x)da. Relative to the
positive real integrand its factor is

    R_v(x,t)=exp(-4v)(1+tx)/(1+exp(-2v)tx).

The denominator has positive real part. With
c=exp(-2rho)cos(2rho), its modulus is at least c(1+tx).
The two coefficients of the numerator of R_v-1 are bounded by
exp(4rho)-1; hence |R_v-1|<.05. Explicit rational guards are
exp(.04)-1<1/24 and c>(49/50)(4999/5000), whose quotient is below 1/20.
Thus

    |P_M,v(k)-P_M(k)|<=.05 P_M(k),                     (19)

uniformly in M,k. The convergent positive angular remainder inherits
this form bound with its actual d_n. The common FOURTH restoration,
curvature, variable-v commutators, global domain and W4 are excluded
from (19). Replacing constant v by v(Y) pointwise does not prove a
variable-coefficient Dirac operator estimate.

## 9. Remaining analytic task and verification scope

Equation (7) closes an actual variable-parameter common finite contact.
Equations (15),(18) give a usable combined-energy bridge for vanishing
coefficients after their actual principal normalization is derived.
The fixed physical-half, spin, time and sphere factors must be applied
at that identification, rather than inserted into the scalar lemmas.

The missing object is the derivative of the actual paired,
common-subtracted two-vertex bubble: three resolvents with a
mass-potential insertion, differentiated seagulls, and the same common
a0/a1/a2 finite parts. Its first-jet cancellation and gain-two remainder
must be quantified in fixed REAL normal charts, together with uniform
actual-minus-cap and far-boundary estimates. The principal fourth/log
piece must be matched to the growing Weyl log(r/B) at the SAME mu0.
Neither (9), real seed positivity, nor the scalar model supplies that
complete estimate. No constant for it is selected here.

The selected endpoint powers make coefficient-squared times log(1/Y)
bounded, which is exactly the mechanism used in (15),(18). This avoids
demanding separately finite positive-log and Weyl-log energies. It does
not establish arbitrary high-regularity mapping from a C5 label: on a
compact interior, n_N=N^(-5)sin(NY) has bounded C5 norm but its product
with a fixed nonzero low-frequency test has H^s norm growing like
N^(s-5) for s>5. Any chosen higher weighted level needs its own coefficient
and commutator bounds.

The proposed shared finite verification is
[stationary_reduced_cold_response.py](../checks/stationary_reduced_cold_response.py),
with [JSON receipt](../evidence/stationary-reduced-cold-response.json) and
[raw log](../evidence/stationary-reduced-cold-response.log), shared with
the companion cold-response note. No execution result is asserted here.
Finite identities and rational bounds are checkable; form closures,
fractional limits, actual domain transfer and infinite kernel estimates
remain ordinary analytic proofs. The original external helpers remain
preserved separately.

Remaining to-do list: prove the actual variable-parameter three-vertex,
common-contact and Weyl matching estimate; establish the selected higher
weighted map and the coupled product-to-actual-seam inverse; construct
the reduced nonlinear branch with physical observable-error bounds.
The original signed source condition and S139/S074-S080 gates remain open.
