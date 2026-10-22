# The actual positive inverse in the weighted second-order domain

2026-10-21, E003/S167. Ordinary argument before focused verification.

The positive auxiliary inverse used by the signed-moment certificate
preserves the required outer weight and gains two scalar derivatives.
This is a theorem about the actual full finite-comparison operator,
including its nonlocal tensor term. It supplies the first positive
solve on the actual residual. Its fourth-order certificate-domain
membership and three transferred sources remain unproved.

## 1. Operator and quantitative result

Keep the original physical Hessian, comparison g_ell, regular bolt,
physical measure dmu_ell and fixed tensor return. Use the notation of
the [moment certificate](coercive-signed-moment-certificate.md):

    T=-Box_ell/lambda, H=T+4, g=G lambda,
    E=60pi/g+11/3, e=g/(60pi), W=t_seed^2,
    A=V_T/sqrt(E), Q=A A*,
    K=H+Q-V_s J_s V_s*/E.

Here A is an auxiliary vertex, not the physical orbit radius. All
derivatives below are in units gradbar=lambda^-1/2 grad. The inverse
K^-1 is the already constructed positive form realization. For every
q with Wq in physical L2, its unique solution u=K^-1q satisfies

    ||Wu||             <=4200 ||Wq||,
    ||W gradbar u||    <=9500 ||Wq||,
    ||WTu||            <=18000 ||Wq||,
    ||W gradbar^2 u||  <=100000 ||Wq||.             (1)

Define X2_W as the closure of regular-even-bolt, smooth radial functions
compact at the outer end in the norm

    ||Wu||+||W gradbar u||+||WTu||.

Then K:X2_W -> L2_W is a bounded bijection with inverse (1).
This second-order domain is not identified with D_cert, whose norm
also contains W T(T+4)u. No physical boundary condition is changed.

For the actual density-correct Schur source p of the
[preceding starting-source theorem](finite-strong-schur-source.md), let
u0=K^-1p and N(v)=||Wv||/sqrt(gamma0), also for tensor-valued v. Then

    N(u0)<8.4e31 ell^2,
    N(gradbar u0)<1.9e32 ell^2,
    N(Tu0)<3.6e32 ell^2,
    N(gradbar^2 u0)<2e33 ell^2,                    (2)

for 0<ell<=1; u0=0 at ell=0. These are analytic upper bounds, not a
numerical evaluation or a small physical metric correction.

## 2. Strong tensor term and the exact local coefficients

The preceding theorem gives

    ||WA||<20000/sqrt(g E)<1500,
    ||Q||<1500^2, ||WQ u||<=1500||A*u||.           (3)

Only the output is weighted. W is not commuted through A, A* or Q.
The entire tensor determinant and mass/return response remain in A.

The [general-background scalar block](general-background-newton-operator.md)
and its [finite coefficient bound](finite-comparison-scalar-curvature-remainder.md)
give, on radial scalars,

    a=Zbar/36+(11/9)Jbar_rr, Zbar=(R4+12lambda)/lambda,
    alpha=a/E, b=Zbar/(3e E), s=1+alpha,
    K=-divbar(s gradbar)+4+b+Q
      =s T-gradbar(alpha).gradbar+4+b+Q.           (4)

The sign in front of the local divergence is retained. In particular
K-H is not just Q. The radial frame is parallel along the actual
proper-radius geodesics, so differentiating Jbar_rr there costs at
most the full covariant derivative of Jbar. The retained bounds give

    |a|<10, |gradbar a|<4e17, E>1.125e35,
    |alpha|<1e-32, |gradbar alpha|<4e-18,
    |b|<2e-12, |Ric|<15lambda.                    (5)

For the gradient use |gradbar Zbar|<.0017 and
|gradbar Jbar|<3e17 from the weighted coefficient theorem, since W>=1.
For b use the sharper earlier |Zbar|<3.1e-12/t and eE>1.
No derivative of |a| or of an unsigned source bound is inferred.

Equations (3)-(5) also identify the unweighted operator domain:

    Dom(K)=Dom(T).                                (6)

Indeed the symmetric perturbation of T in (4) is T-bounded. Spectral
calculus gives ||gradbar u||<= (||Tu||+||u||)/2. Its relative T-bound
is at most 1e-32+2e-18<1; 4+b+Q is bounded. The
[Kato--Rellich theorem, Teschl 6.4](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf)
therefore gives the selfadjoint operator on Dom(T). Its closed form
agrees with the original K form on the common regular core and hence
on its form closure. This identifies the existing realization rather
than choosing a new extension. Scalar Bochner and the bounded Ricci
tensor give global unweighted second derivatives on this domain.

## 3. A global weight budget on the actual metric

Write y=sqrt(lambda)rho_seed and z=sqrt(lambda)rho_ell. The original
finite metric has dz=exp(w)dy and physical radial measure density
m(z)=8pi^2 a_ell r_ell^2/sqrt(lambda). The
[finite low-jet theorem](finite-comparison-classical-tt-form.md) gives

    |w|<3e-8, |w_y|<1.03e-6, |S_y|<7.1e-7,
    A0>=1/sqrt(2), R0>=0,
    (log m)_z=exp(-w)(A0+2R0+w_y+2S_y)>2/3.        (7)

For example the lower bound
(1-3e-8)(.7071-2.45e-6)>2/3 suffices. For a regular compact scalar v,
integrate (m|v|^2)' and use Cauchy--Schwarz. Since m(0)=0 at the
regular bolt, there is no inner boundary term or zero-value condition:

    ||gradbar v||^2 >= ||v||^2/9.                  (8)

Closure gives (8) on the physical scalar form domain.

The [retained seed polynomial bound](global-compatible-scalar-response.md)
is F0<=lambda r_seed^2+1. It applies to this seed function inside W,
not to an independent finite areal-lapse polynomial. Consequently

    |gradbar log W|^2=4exp(-2w)R0^2
      <=4exp(6e-8)(1+1/x)<4.000001, x>7.9e9.      (9)

This is the sharper weight budget needed here; replacing it by the
older reserve twenty would lose the positive weighted energy. One
exact bound is exp(6e-8)<=1/(1-6e-8). The truncations W_R=min(W,R)
have the same logarithmic-gradient bound almost everywhere, with
|gradbar log W_R|<2.001. They are allowed Lipschitz form multipliers.

## 4. Retain unweighted tensor energy before applying the weight

The previous certificate proves

    K >= cH+A A*, c=1-3e-13.

For u=K^-1q, testing against u and dropping its nonnegative gradient
energy gives

    ||A*u||^2 <= ||u||||q||-4c||u||^2
              <=||q||^2/(16c),
    ||A*u||<.251||q||.                            (10)

This uses the maximum of a scalar quadratic, not commutativity of the
operators. Hence (3) gives ||WQ u||<376.5||q||<=376.5||Wq||.

One must not instead declare <W^2u,Qu> nonnegative. In the exact
radial Hilbert space L2([1,infinity),t^2dt), take W=t^2,
A c=c t^-4 and u=t^-4-2t^-6. Then Q=A A* is positive, but

    <A,u>=-3/35, <WA,Wu>=1/3,
    <W^2u,Qu>=-1/35.                              (11)

Both fields have the needed weighted norms. This control explains why
the retained unweighted energy in (10) is used before weighting.

## 5. Weighted value and first derivative

Test (4) against W_R^2u and set v_R=W_Ru. Integration by parts gives
the exact identity

    integral s|gradbar v_R|^2
      +integral [4+b-s|gradbar log W_R|^2]|v_R|^2
      =Re <v_R,W_R(q-Qu)>.                        (12)

The local scalar coefficients in (5) and (9) bound the left side below

    (1-1e-32)||gradbar v_R||^2-2e-6||v_R||^2
      >= (||gradbar v_R||^2+||v_R||^2)/11.         (13)

For the last step use (8); the remaining coefficient is
1/11-1e-32-9(2e-6)>0. The right side of (12) is at most
377.5||Wq||||v_R|| by (10). Therefore

    (||gradbar v_R||^2+||v_R||^2)^(1/2)
      <=11(377.5)||Wq||<4200||Wq||.               (14)

Monotone convergence for the values and weak derivative compactness
as R increases give v=Wu in H1 with the same bound. This is not an
assumption that W preserves the unweighted operator domain. Since
W gradbar u=gradbar(Wu)-Wu gradbar log W, the joint bound (14) and
sqrt(1+2.001^2)(4200)<9500 prove the first two inequalities of (1).

## 6. Weighted Laplacian, Hessian and core closure

Equation (6) gives Tu in ordinary L2. Using (4) distributionally,

    s Tu=q-Qu-(4+b)u+gradbar alpha.gradbar u.

Every term on the right now has its W-weighted L2 bound. Equations
(5), (10) and (14) imply

    ||WTu|| <= [377.5+(4+2e-12)4200+(4e-18)9500]
                                      ||Wq||/(1-1e-32)
              <18000||Wq||.                      (15)

A weighted scalar Bochner identity avoids derivatives of the weight
above its first order. Put L=||W_R Tu||, G=||W_R gradbar u|| and
F=||W_R gradbar^2u||. With beta=2.001, integration by parts and
|Ric|/lambda<15 give

    F^2 <= L^2+2beta L G+2beta F G+15G^2,
    F <= L+(2beta+sqrt(15))G.                     (16)

For clarity the two cross terms before absolute values are
2 integral W_R gradbar W_R.[(Boxbar u)gradbar u
                                      -Hessbar u(gradbar u)].
The scalar curvature term is -integral W_R^2 Ricbar(gradbar u,gradbar u).
Completing the square in the first line yields the second. Smooth
bounded approximations to W_R and regular core approximation justify
the integrations. The unweighted Hessian already follows from (6),
so no weighted Hessian membership is presumed in this argument.
Taking R to infinity and using (15) gives

    ||W gradbar^2u|| <=[18000+(4.002+4)9500]||Wq||
                      <100000||Wq||.

Finally use regular outer cutoffs chi_L(z), equal to one up to L and
zero after L+1. Their first two derivatives are uniformly bounded;
the radial measure logarithmic derivative is bounded on their supports
for large L by the actual outer geometry. The commutator

    T(chi_Lu)-chi_L Tu
       =-2gradbar chi_L.gradbar u-(Boxbar chi_L)u

converges to zero in L2_W by the weighted value/gradient tails. Local
H2 approximation in regular bolt Cartesian coordinates, averaged over
the static symmetries, supplies smooth even compact radial fields.
Thus u belongs to X2_W. Formula (4) is bounded X2_W->L2_W because
||WQ u||<=1500^2||u||. It is onto by this construction and injective
by the original unweighted coercivity. This proves (1) and the stated
isomorphism with the actual operator and physical measure.

## 7. Consequence and remaining obstruction

Apply (1) to the already bounded actual p, with N(p)<2e28 ell^2,
to obtain (2). This removes an unproved weighted second-order inverse
assumption from the proposed positive-solve route.

It does not prove W D0u0 in L2 for D0=T(T+4), or approximation in
the stronger D_cert norm. On suitable test inputs the exact next
source is still D0u, with the distributional identity

    D0u=T[q-(K-H)u].                              (17)

Controlling the two derivatives on the right, including those of the
actual nonlocal and local coefficients, remains necessary. The new
X2_W core closure cannot be substituted for that fourth-order closure.
The next work must establish those estimates for the actual source or
provide admissible trials with the earlier stated residual errors.
No small signed moment, signed zero, nonlinear solution or stability
claim follows from (1). No Lean run is selected.

Freeze this argument, checker and source review. Check the exact
local/weighted identities, physical coefficient budgets, independent
warped scalar and Hessian calculations, and continuum rank-two
fixtures with regular bolt values and weighted tails. Retain the
negative weighted-positive-form control (11). The fixtures are not the
actual charged kernel or an evaluation of its signed coefficient.

Remaining to-do list: fourth-order certificate-domain and weighted
transfer estimates; signed compatibility, cold/higher metric response,
nonlinear preservation and existence; physical boundary/state and
validity checks, stability, exterior matching, charge/return, opening,
transport, preparation, separate creation and full audit.
