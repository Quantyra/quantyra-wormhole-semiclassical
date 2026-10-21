# Outer continuous causal shape bound from four thermal anchors

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

For the retained charged seed and fixed HH state, the actual canonical
causal shape response has the required outer logarithmic compensation
at EVERY sufficiently large real Laplace frequency. The bound is uniform
as the compact radial support approaches the outer boundary. Combining
it with the unchanged local action gives the full curvature-free causal
shape coercivity bound for outer-supported shape inputs and their global
conformal reconstructions.

This advances beyond discrete thermal frequencies. It does not establish
global causal coercivity for arbitrary horizon-supported inputs, a closed
global causal shape realization, or its uniformly bounded Schur feedback.
Those are still needed by the preceding conditional resonance theorem.
No gravitational phase, nonlinear geometry or physical stability is selected.

## 1. Common normalization and the retained thermal anchors

Keep the action, state, flux, period, angular tower, paired matter domain,
common subtraction and every finite/residual contact. Use the complete
even auxiliary outer cylinder of
[the uniform logarithmic comparison](uniform-logarithmic-frequency-comparison.md).
Set tau=a0 t_E, a0=a_phys(0)>0, so its metric is

    dY^2+exp(2h(Y)) d tau^2, |h|<=.04,
    H_s=A+s B_tau,
    A=-exp(-h) partial_Y exp(h) partial_Y >=0,
    B_tau=exp(-2h), b_-=exp(-.08)<=B_tau<=b_+=exp(.08).

All radial forms use measure exp(h)dY, the even extension, and one
physical half. Equivalently conjugate by exp(h/2) to unweighted L2.
Divide the thermal forms by the tau-mode norm. Converting back to
division by the physical t-mode norm multiplies ALL forms and Jj by
the same a0; no physical clock or period is changed.

Let s=(gamma/a0)^2, omega_T=kappa_h/a0 and s_n=omega_T^2 n^2. Thus
the actual physical Laplace frequency is gamma=a0 sqrt(s), and the
thermal values are gamma_n=n kappa_h. Write

    J0=||v||^2, J1(s)=<v,H_s v>, J2(s)=||H_s v||^2.

Initially v is real smooth and compact in 0<Y<.02, extended evenly.
Let D(s)[v] be its ACTUAL renormalized canonical retarded form at i gamma,
and Qsc(s)[v]=<v,Kflat(sqrt(H_s))v>. This comparison is not a replacement
of the actual causal response. The existing two-sided thermal estimates
give finite A0,B0>=0, independent of n and the support's seam distance,

    |D(s_n)[v]-Qsc(s_n)[v]|
        <=C[J2(s_n)/500+A0 J1(s_n)+B0 J0].              (1)

The coefficient 1/500 is .001 from the actual/geodesic comparison plus
.001 from the geodesic/scalar comparison. Their lower-order constants
are retained, not evaluated numerically here. The canonical LLL and
the actual second vertex are included in D.

## 2. Two-sided interpolation of the actual causal form

The [subtracted spectral representation](thermal-shape-response-and-subtractions.md)
proves on each compact probe that the third divided difference of D
is a positive measure integral. This remains true after rescaling s,
with the positive measure rescaled accordingly. For a<b<s<c<d,

    I_(b,c,d) D(s) <= D(s) <= I_(a,b,c) D(s),            (2)

where I denotes quadratic interpolation of the ACTUAL anchor values.
Indeed D(s)-I_(a,b,c)D(s) is (s-a)(s-b)(s-c) times the positive
four-denominator integral, and the product is negative. For (b,c,d)
the product is positive. The bounds include the endpoints by continuity.
They are quadratic-form inequalities, also after polarization.

Every prescribed local polynomial of degree at most two is reproduced
exactly in (2). No subtraction coefficient is dropped or fitted. This
is a bound on the actual spectral function between its anchors, not
an assumption that the actual response equals an interpolation.

Choose a=s_(n-1), b=s_n, c=s_(n+1), d=s_(n+2), n>=20. For s in [b,c],
each interpolation has total absolute Lagrange weight at most 3/2.
For the upper one the only negative weight is the a weight, so its
maximum total is

    1+(2n+1)^2/[8n(2n-1)] <=3/2.

For the lower one it is the d weight, giving

    1+(2n+1)^2/[2(4n+4)(2n+3)] <=3/2.                 (3)

For each of these four nodes t, |t-s|/s<=4/n+4/n^2<=1/4.
Because b_+/b_-=exp(.16)<6/5 and H_s>=s b_-,

    J2(t) <= (13/10)^2 J2(s),
    J1(t) <= (13/10) J1(s).                            (4)

The first uses H_t=H_s+(t-s)B_tau and the norm triangle inequality;
the second uses B_tau<=H_s/s in quadratic forms. The negative
interpolation weights are controlled by (1), (3) and (4), not by
three lower bounds. Their total error is at most

    C[507 J2(s)/100000+(39/20)A0 J1(s)+(3/2)B0 J0].    (5)

## 3. Uniform third-derivative bound for the scalar comparison

The remaining interpolation error of Qsc can be bounded in J0, even
though Qsc itself is an unbounded fourth-order logarithmic form. Its
retained multiplier is

    Qsc(s)=-C F_eff <v,H_s^2 v>/(120pi)
       +sum_m d_m M_m^2/(2pi) integral_0^1
          <v,R(alpha(1-alpha)H_s/M_m^2)v>d alpha,
    R(z)=log(1+z)-z+z^2/2.                              (6)

Here M_m^2=m(m+q)/C and d_m=2(q+2m), with the original integer q.
The polynomial terms in (6) have zero third s derivative, including
the noncommuting expansion of H_s^2. For c_alpha=alpha(1-alpha)>0,
put m_alpha^2=M_m^2/c_alpha and R_lambda=(H_s+m_alpha^2+lambda)^-1.
The resolvent integral for log gives

    partial_s^3 log(H_s+m_alpha^2)
        =6 integral_0^infinity
                R_lambda B_tau R_lambda B_tau R_lambda B_tau R_lambda
                                                               d lambda.

The integrand is positive: it is R_lambda^(1/2) times the cube of
R_lambda^(1/2)B_tau R_lambda^(1/2) times R_lambda^(1/2).
Its norm integral is at most 2 b_+^3/(m_alpha^2)^3. This argument does
not commute H_s with B_tau. Bounded B_tau and the common nonnegative
A realization justify the resolvent derivatives on the whole auxiliary
cylinder. Compare the general resolvent derivative representation in
[Schweitzer, section 2.2](https://arxiv.org/pdf/2203.03930); the present
bounded-operator estimate follows directly from the displayed integral.

Since integral_0^1 alpha^3(1-alpha)^3 d alpha=1/140,

    0 <= partial_s^3 Qsc(s)[v] <= C K3 J0,
    K3=b_+^3 C/(35pi q).                               (7)

We used the exact telescoping identity and its bound

    d_m/M_m^4=(2C^2/q)[1/m^2-1/(m+q)^2],
    sum_m d_m/M_m^4=(2C^2/q)sum_(j=1)^q 1/j^2 <4C^2/q.

The same summable majorant permits the third derivative and finite
divided differences to pass through the entire angular tower. For fixed
compact v, (6) itself converges in the retained H3 core by its subtracted
positive remainder. No truncation of the actual tower is made.

The scalar divided-difference remainder formula now bounds either
quadratic interpolation error by

    C K3 (d-a)(c-b)^2 J0/24.                            (8)

For example |(s-a)(s-b)(s-c)|<=(d-a)(c-b)^2/4; the third derivative
contributes the factor 1/6. The other interpolant is bounded the same way.
Since J2(s)>=b_-^2 s^2 J0, n>=20 and b_-^-2<6/5, (8) is at most

    C [2K3 omega_T^2/n] J2(s).                         (9)

This follows from (d-a)=omega_T^2(6n+3), c-b=omega_T^2(2n+1) and
(6/5)(6n+3)(2n+1)^2<=48n^3 for n>=20. In particular define the
finite integer threshold

    N_I=max(20,ceil(2000K3 omega_T^2)).                 (10)

Then every interval with n>=N_I has scalar interpolation loss at most
C J2/1000. No endpoint-uniform bound for the actual spectral measure
has been assumed: only its sign on compact probes and the already
uniform thermal comparisons enter the argument.

## 4. Actual continuous outer logarithmic compensation

Combining (2), (5) and (8)-(10), for every s>=s_(N_I),

    |D(s)[v]-Qsc(s)[v]|
       <=C[607 J2(s)/100000+(39/20)A0 J1(s)+(3/2)B0 J0]. (11)

This two-sided estimate holds uniformly for all compact profiles in
the outer collar, with no dependence on the distance to Y=0.
The retained scalar logarithmic inequality holds for every REAL s>=0:

    Qsc(s)[v]-C/(30pi) integral (omega+19/60)|H_s v|^2
       >=C[33 J2(s)/(600pi)-12 J1(s)/(30pi)].            (12)

Here omega=log(r/B) is the original Weyl factor and the pairing uses
the same outer half-density convention. The reference radius B is
distinct from the multiplier B_tau.

Source-location correction: the earlier
[canonical frequency note](canonical-frequency-core-bound.md), section 4,
labels its Beckner link as Corollary 2, equation (10). That locator does
not match the linked paper. The needed logarithmic inequality follows
instead by differentiating the sharp Pitt inequality (3) on page 1 of
[Beckner's linked paper](https://arxiv.org/pdf/math/0612486v2) at alpha=0.
Its 2pi-Fourier convention gives psi(1/4)-log(pi); changing to momentum
p=2pi xi adds log(2pi), giving

    D_H=psi(1/4)+log 2=-gamma_E-2log 2-pi/2>-4.

Thus the retained normalization and numerical bound used in (12) are
unchanged. Historical frozen sources are preserved; this note corrects
their citation locator and supplies the differentiation explicitly.

The strict budget 33/(600pi)-607/100000>1/(100pi) gives

    D(s)[v]-C/(30pi) integral (omega+19/60)|H_s v|^2
       >=C[J2(s)/(100pi)-A_c J1(s)-B_c J0],             (13)
    A_c=12/(30pi)+(39/20)A0, B_c=(3/2)B0.

All constants are finite and frequency/support independent on this
outer core. Their actual numerical lower-order values remain uncalibrated.
If also

    s b_- >= max(1,400pi A_c,sqrt(400pi B_c)),            (14)

the right side of (13) is at least C J2(s)/(200pi). The proof uses
J1<=J2/(s b_-) and J0<=J2/(s b_-)^2. Equations (13)-(14) now concern
the ACTUAL causal response at every sufficiently large real gamma,
not only the thermal grid and not an arbitrary continuation.

## 5. Full curvature-free causal core for outer-supported inputs

Let v have the outer support just specified, w=r v and
p=-(1/(3r))K_gamma^-1 B_gamma w be its global curvature-free
reconstruction. These physical reconstruction symbols are those of
[the global classical theorem](global-classical-curvature-free-form.md),
which already holds at all sufficiently large real gamma. Its bounds
S_C>N_gamma/(32G), Gv<=1.001N_gamma, E_gamma<10N_gamma and the weighted
p estimates are unchanged.

The local Weyl square completion and residual budget in
[the global thermal theorem](global-thermal-curvature-free-form.md),
equations (8)-(10), are identities and inequalities for the same real
frequency polynomial. They do not require gamma/kappa_h to be an integer
once their radial integrals and integrations by parts are justified.
At the outer end p=O(y^4) with its retained differentiated jets. At the
horizon set rho proportional to exp(-kappa_h y), nu=gamma/kappa_h.
The positive reconstruction branch has

    p=rho^nu times a function smooth in rho^2.

For nu>2 its radial derivatives through fourth order are O(rho^(nu-j)).
Near the bolt, the local invariant quadratic density has at most four
derivatives and smooth background coefficients. Its worst radial
integrability power is rho^(2nu-3)d rho, and its boundary flux is
O(rho^(2nu-2)); both are harmless for nu>2. This argument does not
call a noninteger mode a smooth Euclidean harmonic. It justifies the
radial causal polynomial and its unchanged endpoint terms. The compact
v has no horizon contribution, while the full p tail is retained.

After (14), gamma>=Gamma_C and gamma>2kappa_h, the complete causal
curvature-free core form therefore satisfies, in physical normalization,

    S_full,gamma[v] >= Gv/(40G)
       +C J2_gamma[v]/(200pi)+C J2_gamma[p+v/3]/(40pi).  (15)

This follows by the same exact cancellation of the negative Weyl
logarithm, its <1e-24 N_gamma/G residual budget, and the real-frequency
classical bound. The canonical contribution is now (13)-(14), so (15)
is valid between thermal frequencies as well. The physical J2_gamma
norms have the original gbar/time normalization.

## 6. What this closes and what remains

The outer quantum compensation and full outer-supported core estimate
are now continuous-frequency results for the retained seed and state.
In particular outer annuli with gamma times their optical width tending
to a fixed positive value no longer need thermal-frequency rounding.
The constants in (13) are uniform as the annuli approach the boundary.

This does NOT give (15) for arbitrary profiles extending to the horizon.
The auxiliary outer cylinder has bounded B_tau; the full canonical disk
has B_tau unbounded at the bolt. Applying the derivative estimate there
with the same b_+ would be invalid. Global localization of the ACTUAL
causal response and its contacts, the horizon region, closure/core-to-H
mapping, and continuity of the full eliminated strong operator remain
to be established. Thus the actual continuous-family hypothesis of
[the boundary-phase theorem](coupled-boundary-phase-winding.md) is still
open. No continuous charged resonance or physical instability follows.

The constrained boundary/initial-state domain, seed residual, static
signed condition, nonlinear construction and all-sector stability remain
in scope. No new phase, finite action coefficient, source projection or
order reduction is adopted. No calibrated semiclassical frequency window
is claimed.

Focused checks verify the exact interpolation signs, negative weights,
uniform polynomial budgets, angular telescoping and Fourier normalization.
Independent positive-measure integrals test the two-sided bracket on
nonthermal points. Noncommuting positive matrix logarithms test the
resolvent third derivative against refined finite differences at two
precisions. A radial fourth-derivative density identity checks the retained
noninteger horizon power counting. These are proof diagnostics, not
charged spectral data or a formal continuum certificate. No Lean
campaign is selected.

Remaining to-do list: global continuous causal horizon/localization and
closed shape realization; physical boundary and initial/state domain;
quantitative thresholds, static signed condition, nonlinear existence
and all-sector stability; exterior/state, charge/return, opening,
transport, preparation, separate creation and full audit.
