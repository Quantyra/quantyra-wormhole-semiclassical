# A necessary cold-source budget for the original projection

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof
before frozen verification. **The actual fixed projection direction
has large cold values and derivatives. Its signed coefficient must
be estimated with cancellation before a small cold source can be
certified. This is not a nonexistence result or a metric lower bound.**

The [cutoff localization theorem](quantile-cutoff-localization-and-jets.md)
controls theta. The remaining adjoint factor cannot have a comparably
small higher norm: its established L2 proximity to a rapidly oscillating
reference already forces the following LOWER bounds. No unknown
adjoint derivative is assumed small, and no numerical oscillations
are sampled.

Keep the original physical model and seed projection. Put

    gamma0=4pi beta B^3, y=sqrt(lambda)rho0, t=r0/B,
    G0=F0/x, t_h=1+h, kbar=k/sqrt(lambda),
    P=sqrt(gamma0) p_cut, J=[1,2] in y.

The scalar disk-source component P_S satisfies, for j=0,...,5,

    sup_J |partial_y^j P_S| > .006 kbar^j/sqrt(h).    (1)

These are derivatives in the seed dimensionless proper radius. In
particular sup_J|P_S|>1.9e14. The other disk-source component is
retained in constructing p_cut and in the norm/error estimates; no
uncoupled adjoint is substituted for the full one.

On J define w_c=s[1+log(.01/s)], s=t-1, and the proper-radial
divided source norms

    X_(m,M,w)(v)=max_(0<=j<=m) sup_J |partial_y^j v|/[w j! M^j].

The same full source direction therefore has

    X_(0,100,w_c)(P)>3e43,
    X_(5,100,1)(P)>1e90,
    X_(5,100,w_c)(P)>2e119.                         (2)

The value norm in the first line does not depend on M. These norms
are explicitly defined source norms on J, not canonical-chart P45
norms and not metric norms. No C5 forcing requirement is imposed on
the nonlinear problem merely because the last two numbers are large.

## 1. Exact normalization and the cold annulus

Use the actual full seed adjoint chi0 with the common sign convention
of the [reference/error theorem](on-shell-adjoint-error-and-trace-gradient.md).
Let N=||chi0||_D, sigma0=||t^-4 chi0||_H and alpha be the original
cutoff normalization. The retained identities are

    sigma0=lambda^-1/4 N, .99<=alpha<=1,
    p_cut=theta t^-4 chi0/(alpha sigma0),
    8.7<N^2, N<4.007, 2.9<alpha N<4.007.           (3)

The physical Liouville multiplier is

    U=sqrt(4pi beta r0^2 sqrt(F0))
      =lambda^1/4 sqrt(gamma0) t G0^1/4.

For v=U chi0,S this gives the EXACT scalar source identity

    sqrt(h) P_S=theta b(y) v(y)/(alpha N),
    b=t^-5 (G0/h^2)^-1/4.                          (4)

Thus no physical volume factor is lost in comparing the unit-H
projection direction to a point value. The thermal factor h^-1/2
in (1) is a consequence of the small local volume.

The cold comparison and its proper-coordinate inversion give

    cosh(sqrt(5.98)y)<=s/h<=cosh(sqrt(6.02)y).

At y=1 the lower bound exceeds 5, and at y=2 the upper bound is
below 70. The cold comparison is valid throughout J: reaching
s=.001 would require more than y=2, since 70h<.001. The cutoff
quantile bounds give y(a)<.2 and b_quantile>1.01. Hence theta=1
on J. Here b_quantile denotes the outer quantile, distinct from
the positive multiplier b(y) in (4).

The retained bounds 5.98(s^2-h^2)<G0<6.02(s^2-h^2) imply

    1/14<b(y)<.3,  t^2<1.000001 on J.              (5)

For the upper bound use G0/h^2>5.98(25-1)>.3^-4.
For the lower bound use
(1+70h)^20 6.02(70^2-1)<14^4, with h<1e-30.
Since t and G0 increase, b is decreasing. No derivative bound on
this multiplier is required below.

Let v_ref=U phi_ref be the actual scalar reference amplitude.
The [reference comparison](cold-reference-overlap-and-adjoint-budget.md)
and full adjoint error give on J

    |v_ref-cos(kbar y-pi/4)|<1e-9,
    ||v-v_ref||_(L2(J,dy))<.007001.                (6)

For the second inequality, the D-norm error <.007 controls
integral_J t^-4 |v-v_ref|^2dy; multiplication by sup_J t^2
gives (6). The full two-component error bounds its scalar component.
In particular, with u=sqrt(h)P_S, equations (3)-(6) give

    ||u||_(L2(J))<(.3/2.9)(1.000001+.007001)<.105.  (7)

## 2. An oscillatory test detects the actual source

Put z=y-1 on J and take the nonnegative polynomial window

    zeta(y)=12012 z^6(1-z)^6.

Its integral on J is one, its maximum is below 3 and its L2 norm
is below sqrt(3)<1.733. It and its first five derivatives vanish
at the endpoints. Its total variation is below 6. With
c(y)=cos(kbar y-pi/4), monotonicity of b and (5) imply

    TV(b zeta)<(.3)(6)+(3)(.3)<3,
    integral_J b zeta c^2 >=1/28-3/(4kbar).         (8)

The oscillatory term is integrated against a bounded-variation
amplitude with zero endpoint values; its frequency is 2kbar.
This is an exact integration-by-parts estimate, not a sampled
quadrature at frequency kbar>4e17.

From (4),(6),(8), the actual source has a positive detected moment

    M=integral_J u zeta c
      > [1/28-3/(4kbar)-.3e-9-(.3)(.007001)(1.733)]/4.007
      >.007.                                      (9)

The numerator is positive, so the upper denominator in (3) may
be used for the entire expression. Thus the L2 error allowed for
the coupled adjoint cannot remove this local oscillation. This M
is a chosen test of the projection direction, not the unevaluated
signed physical residual or compatibility moment.

## 3. Derivative lower bounds without adjoint derivative estimates

For each j=1,...,5 set

    phi_j=zeta cos(kbar y-pi/4-j pi/2).

Leibniz differentiation gives

    partial_y^j phi_j=kbar^j zeta c+R_j,
    ||R_j||_2/kbar^j
      <=sum_(i=1)^j binom(j,i) ||zeta^(i)||_2/kbar^i.

The absolute sum of coefficients of zeta is 12012*2^6<1e6.
Since its degree is 12 and 0<=z<=1,

    ||zeta^(i)||_2<1e6 12^i,
    ||R_j||_2/kbar^j
      <1e6[(1+12/kbar)^j-1]<2e-10.                (10)

All boundary terms vanish through j integrations by parts. With
(7),(9),(10),

    |integral_J u^(j) phi_j|
      =|integral_J u partial_y^j phi_j|
      >kbar^j[.007-(.105)(2e-10)]>.006 kbar^j.

Since ||phi_j||_1<=integral zeta=1, this proves (1) for j>=1.
Equation (9) proves it for j=0. If a derivative were only a
distribution and had no finite L-infinity representative, the
corresponding finite-norm conclusion would already fail; the
inequality itself requires no a priori derivative upper bound.
The established smooth compact-interior adjoint makes the ordinary
derivatives applicable here.

The retained physical constants satisfy

    h<2000^-10, kbar>4e17, L_c(h)<81.

On J, w_c<=70h L_c(h)<5670h. Therefore (1) implies

    X_(m,100,w_c)(P)
      >1e-6 h^-3/2 kbar^m/(m!100^m), m=0,5.

Here .006/5670>1e-6. Using h^-3/2>2000^15 proves the cold
bounds in (2). Similarly .006 h^-1/2 kbar^5/(5!100^5)>1e90
proves the unweighted fifth-jet bound. These are lower bounds on
the actual fixed source direction, not merely large upper estimates.

## 4. A necessary budget for the signed coefficient

For the actual finite comparison, the
[range-response theorem](on-shell-finite-range-response.md) uses

    f_comp=f_ell-p_cut L_ell^0(f_ell)/sigma_ell,
    beta_ell=L_ell^0(f_ell)/(lambda^2 sqrt(gamma0) sigma_ell),
    lambda^-2(f_ell-f_comp)=beta_ell P.             (11)

The symbol beta_ell here is a dimensionless signed coefficient;
it is distinct from the physical Euclidean period beta in gamma0.
The source f_ell and functional L_ell^0 have the original FIXED
seed density. The density conversion before York transpose remains
as previously proved. No rescaling of the adjoint changes (11).

If a proposed source estimate requires either norm in (2) of the
normalized subtracted source to be at most epsilon>0, necessary
conditions are respectively

    |beta_ell|<epsilon/(3e43)       [cold value],
    |beta_ell|<epsilon/(1e90)       [unweighted proper C5],
    |beta_ell|<epsilon/(2e119)      [cold proper C5].  (12)

For epsilon=0 the necessary and sufficient coefficient condition
for the subtracted source to vanish is beta_ell=0.

If instead both lambda^-2 f_ell and lambda^-2 f_comp are bounded
in that norm by M_f and M_comp, the triangle inequality supplies
epsilon=M_f+M_comp in (12). This accounts for possible cancellation
between the original and subtracted sources. The result does not
assert that either individual norm has already been bounded.

For the existing projected Newton correction h_ell, the fixed-density
equation B_ell h_ell=-f_comp makes the same conditional statement
with M_comp a bound on lambda^-2 B_ell h_ell. Thus any claimed cold
neighborhood whose equation image has such a source bound must
also meet (12). This is a necessary test, not a sufficient inverse
estimate or a lower bound on h_ell.

The currently proved absolute estimate is only

    |beta_ell| <21.3(2.2e-13)ell^2/g
               <5.3e21 ell^2, g=G lambda>(3e-17)^2. (13)

At ell=1 this allowance exceeds the unit cold-value threshold in
(12) by a factor 1.59e65. This compares two bounds, not a measured
coefficient. It neither proves that beta_ell violates (12) nor that
the solution fails to exist. At ell=0 the residual and coefficient
vanish. The actual signed coefficient remains unevaluated.

## 5. Research decision and verification contract

Before seeking a generic high-derivative upper bound on the compact
projection source, estimate its actual signed coefficient with the
oscillatory and full coupled cancellations retained. Its large cold
norm is now proved, so the Hilbert estimate alone cannot certify a
small cold source. Match the source norm to the nonlinear operator;
do not introduce a C5 source requirement merely to reuse a C5 metric
contract. A compatible weak-source route remains possible but must
separately prove the physical cold metric response.

The prior full finite value/H1 correction bounds remain unchanged.
There is no nonlinear preserved family, projected nonlinear solve,
signed physical compatibility value or stable wormhole conclusion.

Freeze this proof and checker before execution. Reconcile the fixed
cutoff, adjoint-error, reference and finite-response input identities.
Check the Liouville/volume normalization, positive-series cold
geometry enclosures, window integral and endpoint jets, exact
finite-jet/phase derivative identities, remainder, moment and norm
budgets, and signed-coefficient scaling. Check a nonconstant test's
integration-by-parts identity independently and retain a control
detecting loss of window derivative terms. These are algebraic
fixtures, not physical adjoint or residual evaluations. Continuum
duality and interval arguments remain ordinary proofs. Preserve
every frozen run and failed source. No Lean campaign is selected.

Remaining to-do list: cancellation-sensitive signed projection
coefficient estimate; matched cold/local and higher metric response;
nonlinear mapping, preservation and solve; signed physical
compatibility; exterior/state, charge/return, opening, stability,
transport, preparation, separate creation and full-result audit.
