# Annular retarded limit in the full logarithmic form norm

2026-10-20, E003/S168. Ordinary proof before focused verification.
The full research goal remains incomplete.

The actual seed's dilated retarded response converges in operator norm
on an infinite-dimensional space of perturbations supported in a fixed
annulus. The norm retains the fourth-order logarithmic shape term.
Convergence is locally uniform at complex causal frequencies and holds
for every fixed frequency derivative. This strengthens the preceding
[fixed-profile theorem](complex-retarded-outer-limit.md).

This is a localized form operator. Both test arguments are supported
in the annulus. It does not select a physical boundary realization,
control output tested outside the annulus, or prove a global inverse.

## 1. Space and precise result

Fix I=(ell,b), 0<ell<b<infinity. Let S=(-b,-ell) union (ell,b)
be its even double and L=|S|=2(b-ell). Use unitary Fourier transforms
of even extensions, and half of full-line norms. Define

    ||p||_P^2 = (1/2) integral (1+xi^2)^2 |p_e_hat|^2 dxi,
    ||v||_V^2 = (1/2) integral (1+xi^2)^2
                             [1+log(2+xi^2)] |v_e_hat|^2 dxi,
    X_I = closure of C_c^infinity(I)^2 in P x V.          (1)

This is the natural local form order: H2 for the conformal variable
and H2 with a square-root logarithmic Fourier weight for shape.
It is an auxiliary space of localized tests, not a choice of physical
Dirichlet data at ell or b. No actual exterior wall is introduced.

Let B_d(z) be the polarized actual seed form d^4 D_d(z/d^2)/C and
Q_I(z) the polarization of the full normal form in the preceding note,
restricted to X_I. The real core forms are extended sesquilinearly
in the profile arguments, with conjugation in the first slot. The
frequency dependence remains holomorphic; no modulus-square of a
frequency-dependent factor replaces its analytic bilinear formula.
Write X_I^* for the continuous anti-dual and Omega=C minus (-infinity,0].
We prove

    sup_(z in K) ||partial_z^m(B_d(z)-Q_I(z))||_(X_I -> X_I^*) ->0
       for every compact K in Omega and fixed m>=0.     (2)

The supremum includes the whole unit ball of X_I, including sequences
with increasing radial oscillation. The annulus I and the norm stay
fixed. This is not a statement uniform over annuli reaching eta=0 or
infinity, a weaker unweighted metric norm, or an exact nonlinear exterior.

## 2. What must be strengthened in the thermal argument

The established arbitrary-epsilon comparison gives, for ALL outer
thermal core tangents and each epsilon>0,

    |Qcan-Qsc| <= C[epsilon J2+A_epsilon J1+B_epsilon J0]. (3)

The constants are independent of the thermal frequency and of the
support's distance from the seam. This input is the gain-two estimate
in [the uniform comparison](uniform-logarithmic-frequency-comparison.md)
with the actual canonical contacts included as in section 3 of
[the outer-limit proof](outer-coupled-scalar-separation.md).
It is stronger than the rounded .001 inequality alone.

After dilation and with scaled squared thermal frequency a_d in a
fixed compact subset of (0,infinity), (3) gives, uniformly in v,

    |d^4(Qcan-Qsc)/C|
       <= C_I[epsilon ||v||_H2^2
               +A_epsilon d^2 ||v||_H1^2
               +B_epsilon d^4 ||v||_L2^2].             (4)

The constant C_I is independent of epsilon and small d; the lower
constants may be enlarged. Thus this difference tends to zero in
H2 form norm: take limsup as d->0 and then epsilon->0. No numerical
rate is extracted from the unevaluated A_epsilon,B_epsilon.

The fixed flat multiplier remainder obeys
|r(p)|<=epsilon p^4+M_epsilon. Its scaled spectral form likewise
has norm bounded by C_I epsilon+O(M_epsilon d^4) on H2. The -log d
term multiplies the EXACT scaled J2 and cancels its Weyl partner
before any coefficient approximation. All other local coefficients
converge, after integration by parts to form order two, in the
finite coefficient norms supplied by the actual seed jets. Since I
stays away from zero, their error is o(1) in H2 x H2 form norm.
Off-shell contacts, global returns and finite subtractions stay in
these estimates. Their omission is not the source of smallness.

The remaining issue is norm convergence of the logarithmic scalar
form. Strong resolvent convergence alone is not enough for this step.

## 3. A localized logarithm estimate with explicit distant-tail bound

For positive scaled frequencies a_d->a>0, the half-density scalar
operator on the even line is

    H_d=-partial_eta^2+d^2 V(d eta)+a_d exp(-2h(d eta)),
    H_a=-partial_eta^2+a, H_d=H_a+W_d,
    V=h''/2+h'^2/4.                                    (5)

The retained auxiliary metric satisfies |h|<=.04, |h'|<=1.2,
|h''|<=42 and |h'''|<=3270 almost everywhere. Hence |V|<=22 and
|V'|<=2000. The radial part is
(partial_eta-d h'(d eta)/2)^*(partial_eta-d h'(d eta)/2),
so H_d>=a_d exp(-.08). For a_d near a choose a fixed c>0 with
H_d,H_a>=c, c<=a. Also ||W_d||_infinity<=M for a fixed finite M.

For R>b put delta_d(R)=sup_(|eta|<=R)|W_d(eta)|. Since h(0)=h'(0)=0,
the retained jet bounds give, whenever dR<=.1,

    delta_d(R)<=22d^2+exp(.08)|a_d-a|
                         +42a exp(.08)d^2 R^2.          (6)

The first derivative on S also tends uniformly to zero:
||W_d||_(W1,infinity(S))<=C_I(|a_d-a|+d^2).
No assertion that W_d is uniformly small on the whole line is used.

Spectral calculus and the second resolvent identity give

    log H_d-log H_a
       =integral_0^infinity (H_d+t)^-1 W_d (H_a+t)^-1 dt. (7)

This integral converges in operator norm, bounded by M/c. It is
first identified on the common smooth core and then extends as a
bounded difference of the logarithms. The scalar logarithm integral
and the resolvent identity justify the order in (7); the two operators
need not commute. The background spectral-calculus framework is in
[Teschl, section 3.1 and Lemma 6.5](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf).

Let chi_S be multiplication by the indicator of S. The full-line free
resolvent kernel is exp(-sqrt(a+t)|x-y|)/(2sqrt(a+t)); its derivative
jump is -1. Its distant part obeys the Hilbert--Schmidt estimate

    ||1_(|x|>R)(H_a+t)^-1 chi_S||
       <= sqrt(L) exp(-sqrt(a+t)(R-b))/(2(a+t)^(3/4)).    (8)

This follows by integrating its squared kernel in x outside [-R,R]
and y in S. Restricting to even functions can only improve the bound.
Split W_d into |x|<=R and its complement in (7), and use
||(H_d+t)^-1||<=1/(c+t). Integration proves

    ||chi_S(log H_d-log H_a)chi_S||
       <= delta_d(R)/c
          +2M sqrt(L) exp(-sqrt(c)(R-b))/(3c^(3/4)).     (9)

For example choose R=b+2 log(1/d)/sqrt(c). For all sufficiently
small d this satisfies dR<=.1, and (6),(9) give an explicit bound
of the form C_I,a[|a_d-a|+d^2(1+|log d|^2)]. The global completion
enters the exponentially small second term and has not been dropped.
Section 7.4 of the same Teschl source supplies the free-resolvent
framework; (8)-(9) are the localized estimate proved here.

## 4. From logarithms to the critical fourth-order form

For a smooth v supported in S, write w=H_a v and e_d=W_d v. Expand

    <H_d v,log H_d H_d v>-<H_a v,log H_a H_a v>
      =<w,(log H_d-log H_a)w>
          +2 Re<e_d,log H_d w>+<e_d,log H_d e_d>.        (10)

Equation (9) controls the first term by its bound times ||v||_H2^2.
For the cross terms choose ell_c=max(0,-log c). Then
L_d=log H_d+ell_c is nonnegative. The global bound M/c in (7) and
Fourier calculus for H_a give, as forms,

    0<=L_d<=C log(2-partial_eta^2),
    <w,L_d w><=C||v||_V^2,
    <e_d,L_d e_d><=C||e_d||_H1^2.                      (11)

Here log(2+xi^2) dominates a constant and grows slower than 1+xi^2.
Cauchy--Schwarz for the positive form L_d, and the subtracted ell_c
term, bound the last two terms of (10) by
C(epsilon_d+epsilon_d^2)||v||_V^2, where
epsilon_d=||W_d||_(W1,infinity(S))->0. Therefore the logarithmic
fourth-order forms converge in V->V^* norm. The factor 1/2 for
log sqrt(H) and physical-half normalization do not affect convergence.

The original half-density profile is m_d v, m_d=exp(h(d eta)/2).
On a fixed enlargement of S, m_d-1 tends to zero in W3,infinity.
Cut it off there to obtain a global multiplier with the same property
and the same action on the profiles. It is o(1) as a multiplier on
both H2 and H3. The following elementary interpolation calculation
shows it is also o(1) on V, without adding a fractional derivative.

Let K_2(t,v)^2=inf_(v=v0+v1)(||v0||_H2^2+t^2||v1||_H3^2).
Pointwise Fourier minimization, with A=1+xi^2, yields

    K_2(t,v)^2=integral A^2 [t^2 A/(1+t^2 A)] |vhat|^2,
    ||v||_H2^2+2 integral_0^1 K_2(t,v)^2 dt/t
        =integral A^2[1+log(1+A)]|vhat|^2=||v||_V^2.   (12)

Use matching full- or half-line normalization on both sides.
A multiplier bounded by epsilon on both endpoint spaces has
K_2(t,Tv)<=epsilon K_2(t,v), so its V norm is bounded by epsilon.
The half-density change therefore contributes o(1) in the same form
norm. This proves the required scalar step at the logarithmic order,
rather than only on H3-bounded subsets of test profiles.

Combining (4)-(12) with the fixed local polynomial gives

    ||B_(d,thermal)(a_d)-Q_I(a)||_(X_I->X_I^*) ->0       (13)

whenever a_d is a nonzero scaled thermal squared frequency tending
to a>0. The estimates are uniform for a in any compact positive interval.
This is an actual seed statement with the imported common subtraction.

## 5. From thermal anchors to the analytic operator family

For any fixed a>0 the four neighboring actual nonzero thermal values
have scaled squared frequencies tending to a. The preceding real-axis
enclosure applies to every core vector. Its signed interpolation
weights have total absolute value at most 73/48. Applying (13) to
the four samples therefore gives

    ||B_d(a)-Q_I(a)||_(X_I->X_I^*) ->0.                 (14)

Hermitian form bounds give operator bounds after the fixed Riesz
identification. This step controls both sides of each anchor error;
the negative interpolation weights are not ignored.

Set Gamma_d=B_d[1,2,3,4]>=0. Equation (14) bounds this positive form
in X_I norm. Its polarized weighted spectral measure lambda_d has
total Gamma_d. Positivity gives, for every bounded scalar r(t),

    |integral r(t) dlambda_(d;u,w)(t)|
       <=||r||_infinity sqrt(Gamma_d[u] Gamma_d[w]).    (15)

For example prove it for simple r by measure Cauchy--Schwarz and
then pass to bounded limits. Apply it to r_z(t)=(t+4)/(t+z).
The same interpolation formula as before yields the DIMENSION-FREE
bound, for z in any compact K in Omega,

    ||B_d(z)|| <= sum_(j=1..3)|ell_j(z)| ||B_d(j)||
                          +|P(z)| M_K ||Gamma_d||.     (16)

There is no trace factor or finite-basis assumption. These bounds
extend the core forms uniquely to X_I->X_I^*. Locally uniform
domination of the derivatives of r_z also proves operator-norm
holomorphy. The explicit normal Fourier formula similarly defines
a norm-holomorphic bounded Q_I(z) on these spaces.

To prove norm convergence, one must NOT apply a scalar Montel theorem
directly to an infinite-dimensional operator ball. If (2) with m=0
failed, there would be d_n->0, z_n in a compact K, and unit vectors
u_n,w_n with |<u_n,(B_(d_n)(z_n)-Q_I(z_n))w_n>| bounded below.
These scalar holomorphic functions are locally uniformly bounded
by (16). They converge to zero at every positive a by (14), even
though the vectors vary with n. The scalar Vitali theorem therefore
makes them converge uniformly to zero on K, a contradiction.
The standard scalar result is in
[Boas, sections 2.3 and 2.7](https://haroldpboas.gitlab.io/courses/618-2016a/notes2016.pdf).
This varying-test argument supplies the operator upgrade.

Operator-valued Cauchy integration now gives the derivative estimate
m! rho^-m times the norm error on a compact rho-neighborhood. This
proves (2) for every fixed m. Composition with z=s^2 gives the same
statement on compact subsets of the scaled right Laplace half-plane.

## 6. Scope, verification and next transfer

This closes the infinite-profile gap on each fixed annulus in the
critical logarithmic form norm. The strong result is local: it bounds
all bilinear pairings of normalized annular tests. It does not bound
the response against tests outside I, send ell->0 or b->infinity,
choose metric traces, or show that a global inverse exists. No
physical constraint domain, admitted ultraviolet cutoff or nonlinear
background is substituted by X_I. The fast factor remains intact.

The explicit scalar tail estimate is quantitative, but the full
comparison constants A_epsilon,B_epsilon remain uncalibrated.
Consequently no numerical physical collar radius is certified.
There is no uniform claim at z=0 or on the negative spectral cut.

Focused verification checks the critical interpolation identity,
resolvent order/sign, free-kernel jump, tail integrals and seed-jet
constants. Independent noncommuting matrix logarithms test (7) and
its bound. A solvable step potential tests the distant-resolvent
mechanism against independently matched ODE solutions. Positive
matrix-measure fixtures test (15)-(16) without a trace factor.
Controls distinguish fixed-vector from operator-norm convergence
and show why the shape logarithmic weight cannot be discarded.
These fixtures verify identities and bounds, not a physical spectrum.
The continuum theorem remains the ordinary proof above; no Lean is selected.

Next establish the physical gravitational boundary/initial-data
realization and global coupling/endpoint control. Use neutral current
decoupling and Ward recovery, retaining both off-shell residual
contacts. The earlier boundary-solvability test remains necessary
for a proposed fast-branch exclusion. S168.3-S168.4 remain open.

Remaining to-do list: physical metric boundary and global retarded
transfer or justified exclusion; S167 signed source, cold/higher
response, nonlinear solve and compatibility; exterior/state,
charge/return, opening, full-sector stability, transport, preparation,
separate creation and final audit.
