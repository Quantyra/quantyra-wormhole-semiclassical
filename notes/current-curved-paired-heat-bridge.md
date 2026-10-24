# Curved paired heat response: a uniform weak trace argument

2026-10-23, E003/S169. Ordinary argument before acceptance calculations.
Previous goal turn: progress. Preserve the closed audit, selected projection
screen and full stable-semiclassical-wormhole objective. This note supplies
a qualitative C2 action on bounded ADMITTED stationary coefficient families.
It does not prove that the unknown solution's X ball is admitted, compute
Y0/Z0/Z2, or establish a solution or stability.

## 1. Precisely the family and the change in method

Retain the first-order realization and complementary adjoint proved in
current-rough-paired-variation.md, the identical paired copies, fixed
thermal period/spin structure, regular disk bolt and complete Dirac tower.
Use real stationary Weyl factors in the existing weighted C5 phase:

    b=c_b x^3+r_b,
    |r_b^(j)|<=K j!100^j x^(4-j)(1+log(1/x)), j<=5,
    |c_b|<=K, 0<x<=.1.                              (1)

Keep the existing bounded C5 interior and regular-disk charts, uniform
ellipticity, fixed atlas and a bounded global W1,infinity norm. Directions
u,w satisfy the same linear bounds, with their own coefficient norms A.
All constants below depend on these stated family bounds. No uniformity
over unbounded log degrees, collapsing ellipticity, arbitrary temperatures
or an unproved X-to-A embedding is claimed. A fixed finite higher log weight
has the same vanishing powers used below.

The earlier draft's generic rough parametrix assertion remains withdrawn.
Here the seam is treated by an EXACT curved heat formula. Finite Gaussian
parametrices are used only in ordinary C5 charts away from that seam.
No smooth reflected physical metric or finite physical spectrum replaces
the original problem. Auxiliary extensions used in a local parametrix have
their entire commutator defect retained.

## 2. Exact curved collar kernel, including every thermal frequency

In proper distance y from the seam write g2=dy^2+a(y)^2 d tau^2, f=1/a.
Extend f evenly on the local double. On a fixed collar, f is positive,
bounded W3,infinity and generally not C3 across the seam. In the radial
half-density frame the actual square is

    L=-partial_y^2-f^2 partial_tau^2
                         -gamma_y gamma_tau f' partial_tau.

At frequency p its two spin sectors are the ordinary one-dimensional
operators -partial_y^2+p^2 f^2 +/- p f'. Both are defined through the
original closed form; this diagonalization introduces no boundary datum.
The paired spin trace equals the exact even-double spin trace, by the
chirality/reflection identity in stationary-massive-boundary-layer.md.
The physical integrated collar trace is one half of the even-double trace.

Choose a positive bounded extension of f outside a slightly larger collar.
It is an auxiliary operator for the parametrix, not a replacement of the
physical bolt or exterior. Its W3,infinity bounds are uniform. For a real
Brownian bridge xi on [0,1] with covariance

    E[xi_s xi_r]=2(min(s,r)-sr),

put h=sqrt(t), A_h(y)=int_0^1 f(y+h xi_s)^2 ds and
B_h(y)=int_0^1 f'(y+h xi_s) ds. Feynman--Kac in each spin sector and the
ENTIRE Gaussian frequency integral give the diagonal, per unit time,

    k_cont(t;y,y)=1/(2pi t) E[A_h^(-1/2)
                              exp(t B_h^2/(4 A_h))].           (2)

This is exact for nonconstant f. There is no expansion in curvature or
in the rough coefficient. The derivation follows directly by integrating
2 exp(-t p^2 A_h) cosh(t p B_h) over dp/(2pi). Bounded positive f and
bounded f' justify the Gaussian integration and expectation. The same
formula extends holomorphically to a small complex W3,infinity tube with
Re f^2 bounded below. Trotter's product formula gives the complex-potential
version of Feynman--Kac; its Gaussian majorant justifies this extension.

For the ACTUAL antiperiodic period T, Poisson summation gives the exact
replacement of the integrand in (2) by itself times

    sum_(k in Z) (-1)^k exp[-k^2 T^2/(4t A_h)]
                              cos[k T B_h/(2 A_h)].             (3)

Thus (2) is just the zero image of (3), not a zero-temperature physical
substitution. For real f with m<=f<=K0 and |f'|<=K1, all nonzero images
are bounded in absolute value by

    2 exp[-T^2/(4t K0^2)]/[1-exp(-3T^2/(4t K0^2))].             (4)

They retain their exponential prefactor from (2). On the complex tube
Re A_h>=m^2/2, |A_h|<=K0^2, |B_h|<=K1, Young's inequality instead gives
the bound

    2 exp(4t K0^4 K1^2/m^6)
       exp(-c/t)/[1-exp(-3c/t)], c=T^2 m^2/(16K0^4)>0.          (5)

The factor exp(t K1^2/(2m^2)) from (2) is kept as well. These bounds
include all images, and persist for metric derivatives by the argument
below. No finite set of angular or thermal modes is used.

## 3. The integrated weak-derivative estimate that replaces the lost C4 bound

Let chi be a fixed smooth even localization, equal to one near the seam,
supported strictly inside the extended collar. Work first with smooth
coefficients to justify ordinary integrations by parts, then pass to the
EXACT W3,infinity coefficient. Uniform weak derivative bounds, not a
smooth C4 norm, control the passage. The paired domain/projectors stay
fixed. Such approximations justify identities; the limiting physical
field, its cubic coefficient and its log terms are not changed.

Expand only the explicit exponential in (2) through j=2. For

    F_j(h)=int chi(y) A_h(y)^(-j-1/2) B_h(y)^(2j) dy,

Taylor-expand F_j through degree 5-2j in h. Its integral remainder uses
r=6-2j derivatives. Before integration by parts, each differentiated
term has total spatial derivative count r+2j=6, distributed among factors
f(y+h xi_s), and powers of A_h in the denominator. At most ONE factor
can carry more than three derivatives. If one does, integrate in y to
move its excess derivatives onto chi or other factors. At most three
such integrations are needed. Differentiating A_h brings another integral
of 2f f'; it preserves the total derivative count. Once the largest
factor has order three, every other factor has order at most three as
well, since the total count is at most six. No inverse power of h is
introduced: partial_y f(y+h xi_s)=f'(y+h xi_s).

This is an ordinary finite integration-by-parts algorithm for every
term, including denominator derivatives. It does not assume a formal
jet expansion across the seam. All remaining factors are bounded by
W3,infinity norms. Derivatives of chi have order at most three; their
L1 norms control the spatial integral. Path factors are bounded by
(1+Z)^6, Z=sup_s|xi_s|. The reflection estimate
P(Z>z)<=2 exp(-z^2) gives E Z^6<=12 and

    E(1+Z)^6<=32(1+12)=416.                        (6)

For explicit bookkeeping, let K bound the first three coefficient jets
on the complex tube, D=4(1+K+m^-1), and
C_chi=sum_(j=0)^3 ||chi^(j)||_1. In the algorithm there are at most nine
differentiation/transfer steps. A step has fewer than 64 branches and
coefficient magnitude at most 64: there are at most 22 coefficient
factors and the denominator power is below 12. All denominator/numerator
powers are bounded by D^64. Hence 2^120 D^64 C_chi(1+Z)^6 bounds the
finite derivative sum; the deliberately enlarged 2^124 below also covers
the three Taylor remainders and the exponential tail. Odd powers of h
vanish after expectation because xi and -xi have the same law. It follows
that the FULL localized continuous-frequency trace has

    Tr_chi H(t)=a0_chi/t+a1_chi+a2_chi t+r_chi(t),
    |r_chi(t)|<=C_out t^2,
    C_out=2^124 *416*D^64*C_chi/(2pi)
                      *exp[t0 K1^2/(2m^2)], t0=.001.           (7)

The full-line bound is also a conservative physical-half bound. Add the
complete exponential image bound (5), not a dropped thermal term.
Equation (7) supplies a finite, explicit-form majorant on bounded
admitted families. Its enormous generic constant is NOT a useful Z2
certificate. It establishes integrability rather than numerical sharpness.

For nonsmooth f the integration-by-parts statements hold distributionally.
Uniform W3,infinity approximations converge in W2,infinity and locally
H3; translation continuity and the above integrable majorant pass the
integrated formulas to the limit. A jump in f''' is allowed. No fourth
derivative is estimated separately. The formula itself is holomorphic
in f on the indicated W3,infinity tube. Cauchy gives, for tube radius
m/4 and directions eta,zeta,

    |D_f r[eta]|<=4 C_out ||eta||_(W3,infinity) t^2/m,
    |D_f^2 r[eta,zeta]|<=64 C_out ||eta|| ||zeta|| t^2/m^2.     (8)

Use tube bounds in C_out, so the constant on the contour is uniform.
The exponential image remainder obeys the corresponding derivatives.

## 4. Contacts: compute, rather than presume, the local coefficients

Wick integration of the bridge gives

    E mu2=1/3, E mu1^2=1/6,
    E mu4=2/5, E(mu3 mu1)=1/5, E mu2^2=1/5,
    E(mu2 mu1^2)=11/90, E mu1^4=1/12,
    mu_j=int_0^1 xi_s^j ds.

Write f_j=partial_y^j f. The coefficients inside the expectation in (2)
through h^4 are

    C0=1/f,
    C2=(-f f2+2 f1^2)/(6 f^3),
    C4=(-f^3 f4+7 f^2 f1 f3+4 f^2 f2^2
                                  -23 f f1^2 f2+12 f1^4)/(60 f^5). (9)

Put R=2f2/f-4f1^2/f^2 and Delta R=R''-(f'/f)R'. Direct algebra gives

    C2=-R/(12f),
    C4=-R^2/(240f)-(R'/f)'/120.                   (10)

Together with 1/(2pi), these are exactly the retained rank-two bulk
coefficients and their divergence term. Integrating the last derivative
weakly includes the distribution at the jump of f''' on the double.
On the physical half it is the actual paired boundary completion already
present in current-paired-regulator-domain-audit.md. It is not zeroed
before taking the trace. With localization, derivatives of chi are kept
and cancel with the other chart contributions in the global coefficient.
For chi=1 globally the completed a2 is -int R^2 dvol/(480pi).

Thus the weak continuation matches the ORIGINAL subtraction, rather than
defining a new action with convenient counterterms. Existing finite local
restorations, Weyl conversion, product contact and LLL remain unchanged.
The earlier t^(5/2) and t^3 log t boundary terms are consistent with (7)
and are not erased by it; the remainder contains them.

## 5. Return from the proper collar coordinate to the actual metric variables

Let x be the seed proper coordinate and f0 its fixed inverse time radius.
For a real b define

    Y_b(x)=int_0^x exp(b(s)) ds, x_b=Y_b^-1,
    F_b(x)=exp(-b(x))f0(x), f_b(y)=F_b(x_b(y)).       (11)

This is a coordinate identification of the SAME physical collar. Its
boundary, paired projector and thermal period are unchanged. At fixed y,
with I_u=int_0^x exp(b)u and I_uw=int_0^x exp(b)uw,

    xi_u=D_u x_b=-exp(-b) I_u,
    xi_uw=-exp(-b) I_uw-u xi_w-w xi_u-b' xi_u xi_w,
    D_u f_b=-u F_b+F_b' xi_u,
    D_uw f_b=uw F_b-(u'F_b+u F_b')xi_w
       -(w'F_b+w F_b')xi_u+F_b'' xi_u xi_w+F_b' xi_uw. (12)

These follow by twice differentiating Y_b(x_b)=y. Near the seam
xi_u=O(||u||A x^4), xi_uw=O(||u||A||w||A x^7), allowing fixed logs.
After three spatial derivatives, the most singular base term in (12)
is F_b^(5) xi_u xi_w. Equation (1) bounds it by O(x^7 L(x)); terms with
u^(4) have a factor xi_w=O(x^4). Every term is bounded. Away from zero,
the existing C5 charts supply ordinary bounds. Consequently (11) is C2
from the admitted weighted phase to W3,infinity on a fixed smaller proper
collar, with bounded derivatives on bounded families. Continuity follows
by splitting off a small collar, where the vanishing powers control the
norm, and using ordinary C5 continuity on its compact complement.

Equations (8) and (12), with the chain rule, give an O(t^2) mixed metric
remainder in the ORIGINAL u,w directions. No analyticity in the spatial
mass/log fields is asserted. Complex continuation was used only in the
coefficient f for (8); a complex inverse proper-distance map is not used.
This argument provides two metric derivatives. It does not automatically
provide the higher source variations required by Z2.

## 6. The actual finite geometry: regular charts and complete patching defect

Here is the separate finite-regularity construction away from the seam,
including Cartesian charts at the regular bolt. In a fixed trivialization,
write the square as L=-A^(ij)(x)partial_i partial_j+B^i(x)partial_i+C(x).
Uniform metric C5 gives principal C5, first-order C4 and zeroth-order C3
bounds. Freeze A at source y and call its Gaussian kernel P0. Expand
principal coefficients through degree four in x-y, B through degree
three, and C through degree two. Collect parabolic weight r-2 in L_r:
principal degree r, first-order degree r-1, potential degree r-2.
For j=1,...,4 define the FINITE recursion

    P_j(t)=-sum_(r=1)^j int_0^t exp[-(t-s)L0] L_r P_(j-r)(s) ds. (13)

All coefficients of L_r are evaluated at y. L acts on x, so it never
differentiates those source jets. Gaussian polynomial integration gives
P_j=t^(j/2) times a Gaussian of dimension two times a finite polynomial
in (x-y)/sqrt(t). Its parity is j, so odd j vanish on the diagonal.
The time endpoint at the first insertion is at worst s^-1/2; successive
convolutions are finite beta integrals. The finite polynomials and their
first two x derivatives have bounds depending only on the stated jets
and ellipticity. This defines the coefficients, not an unexpanded
assertion about an infinite Volterra remainder. Matrix coefficients keep
their order in (13); they are not commuted as scalar symbols.

Apply partial_t+L to P0+...+P4. The recursion cancels every weight through
four. The remaining finite products have r+j>=5 and hence Schur order
t^((r+j-2)/2)>=t^(3/2). The actual Taylor remainders are bounded by
K|x-y|^5, K|x-y|^4 and K|x-y|^3 respectively; multiplying the Gaussian's
second, first or zeroth derivative gives that SAME t^(3/2) Schur order.
Thus the defect is bounded by C(K)t^(3/2) G_c(t,x,y), uniformly in C5
physical-side approximations. This step uses no sixth metric derivative.

Patch the exact extended collar kernel of (2)-(3) and these regular-chart
kernels using cutoffs whose derivative supports are separated from their
source supports. The initial kernel is the identity because the source
cutoffs form a partition of unity. All cutoff commutators lie away from
their sources and have Gaussian exponential bounds; the full defect is
retained. The true physical kernel is obtained from this parametrix by
the exact error identity K-P=-K*(partial_t+L)P. Gaussian convolution gives

    |K-P|<=C(K)t^(5/2)G_c(t,x,y),
    |Tr(K-P)|<=C(K)t^(3/2),                       (14)

plus the retained exponentially small terms. Constants are finite for
this fixed positive temperature and bounded admitted atlas. Uniform
Gaussian/first-gradient bounds require only the lower coefficient norms:
freezing the principal coefficients gives a defect of Schur order
t^-1/2, whose Volterra powers have denominators Gamma(1+j/2). This
convergent low-order construction supplies those bounds independently of
the high-order remainder. Small complex perturbations of the coefficients
preserve ellipticity and the Gaussian bounds with enlarged constants.

For mixed metric derivatives of (14), differentiate the exact error
identity. The true kernel has the same Gaussian estimates after two
parameter derivatives: the fixed first-order domain and sectorial forms
give holomorphy on a uniform small W2,infinity parameter disk, and Cauchy
keeps the time power unchanged. The complex family is the holomorphic
continuation of the real square through its common sectorial form; it is
not defined by inserting a complex b into an antiholomorphic adjoint.
The W2 bound also controls the square's zeroth-order coefficient; the
W1 graph-domain statement alone is not used as a Gaussian kernel bound.
The admitted A norm controls this W2 norm. Regular-chart defects and their derivatives
have the same t^(3/2) Schur order by the finite recursion. The collar cutoff
defects and their two derivatives are exponentially small: (2)-(3) are
holomorphic in f, and (12) controls the real metric-coordinate derivatives.
Coordinate cutoffs may move with the proper collar; in fixed interior
coordinates they are chi(Y_b(x)), whose derivatives are bounded by (1).
All differentiated cutoff terms are included. Each differentiated
convolution therefore has the same bound in (14).

This construction includes the actual bolt through its regular disk charts;
it does not impose data at a radial buffer wall. The auxiliary collar
extension affects the parametrix defect, not the physical domain. The
finite recursion is the elementary Gaussian version of the smooth heat
construction in [Grieser's sections 2.2--2.4](https://web.math.ku.dk/grubb/notes/heat.pdf).
The explicit derivative orders and weak collar argument above supply the
finite-regularity adaptation; the reference alone is not cited as that
theorem. [Vassilevich's heat coefficients](https://arxiv.org/html/hep-th/0306138v3)
provide the independent convention/contact comparison in (10).

## 7. Common regulator, what closes, and what does not

Combining (7)-(14), with the actual completed local coefficients, gives

    |D_uw r_sub(g;t)|<=C_A ||u||A ||w||A t^(3/2),
    0<t<=.001, C_A<infinity,                       (15)

uniformly on each bounded admitted family and in its physical-side
regularizations. The same bound holds before and after the core limit;
positive-time convergence follows from the fixed closed forms, or the
exact kernel construction. For the complete angular weight Theta(t)<=4C/t,

    int_0^.001 |D_uw r_sub| dt/t^2
       <=.063245553204 C_A ||u||A ||w||A.           (16)

The fixed positive angular mass gap controls large time; for a bounded
family, graph norm equivalence retains the needed spectral bounds. All
thermal images and the LLL's original kernel/phase allocation remain in
their original prescriptions. Finite restored local terms are C2 in this
weighted phase. Dominated common-regulator removal, and Taylor integration
from the already defined seed, identify the action and its first two
metric variations on these admitted families. Continuity follows from
the same domination and positive-time form convergence.

This closes the qualitative curved C2 bridge ON THE ADMITTED STATIONARY
SCALAR FAMILY. It does NOT prove a C2 map on the unknown solution's X ball
until its admission bounds are established. The finite C_A is now supplied
by an ordinary construction; its complete numerical value, atlas losses,
X-to-A embedding and inverse remain unevaluated. The deliberately excessive
universal allowance in (7) is not evidence of a large actual correction.
No source third-variation or Z2 bound is inferred from C2 of the action.

The selected sufficient screen remains <2.123e-42<1.287e-40. The complete
reduced-source allowance remains <1.3079691e-40<2.60e-40. Neither is a
correction norm. At rho=1e-8 the unchanged requirements remain

    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1,
    Y0<=5e-9, Z0<=.25, Z2<=5e7 (unproved allocation),
    1.632e-11 Km+1.36e-9 Kc+Yrest<=5e-9,
    G_(1/4)^lo E_D^2<=1/402,
    b45+1e-8 E45<=1e-3, b_c+1e-8 E_c<=5e-7,
    b_o+1e-8 E_o<=1.5e-7, b_D+1e-8 E_D<=2.5e-10.

The existing conditional thresholds Km<=61.2745098039..., Kc<=1.47058823529...,
Yrest<=2e-9, E45<=87000 when b45=1.3e-4, and E_D<=.025 when b_D=0
are unchanged. At that last endpoint G<=800/201. Actual operator/return
constants are still UNKNOWN. The independent bolt-difference target is
<1.667e-9. No normalization change is counted as improvement.

Research decision: CONTINUE THE SAME CONSTRUCTION. Use the exact curved
collar trace and weak finite-jet argument for the domain bridge; do not
optimize the already passing cold/scalar/projection screen. The SINGLE
NEXT CALCULATION is the actual coupled Schur-reduced forward/adjoint
correction with a certified residual and return to the weighted mass/log
admission norms. Test its bolt difference and Y0/admission requirements
together. This solve is PROPOSED, NOT EXECUTED here. The curved domain and
uniform integrability argument is EXECUTED; its numerical nonlinear and
physical consequences are not presumed.

Remaining to-do list: actual controlled correction and weighted return;
X-ball admission, quantitative inverse and higher source remainder;
exact compatibility; physical matching and causal realization;
semiclassical validity/noise and all-sector stability; charge/return,
opening, transport, preparation, separate creation and the full stable-
wormhole objective.
