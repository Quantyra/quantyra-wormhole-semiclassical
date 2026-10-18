# Covariant reduced equations and their actual seed inverse

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary theorem before
focused verification. **PROGRESS; nonlinear physical error control open.**

This gives a particular first-order EFT field-redefinition prescription,
its exact reduced scalar equation and an actual coupled inverse at the
analytic centered seed. In its stated critical-energy norm the inverse
bound is 1000/989. The mass-extended fast inverse also exists, without
a numerical fast/cold constant yet. These results do not evaluate the
original signed source condition or establish a nonlinear physical branch.

Continue the [actual TT inverse](stationary-weighted-tt-inverse.md),
[complete scalar operator](stationary-scalar-outgoing-inverse.md) and
[original source condition](stationary-scalar-phase-moment.md), using
the [analytic seed](stationary-centered-covariant-seed.md),
[global disk gauge](stationary-disk-and-normal-inverse.md) and
[actual weighted map](stationary-global-weighted-map.md).
The fixed clock, positive prepared temperature, flux, retained Dirac
state/domain and subtraction remain those of that theory. No member of
the original preparation interval is chosen to tune a source moment.

## 1. The fixed field-redefinition prescription

Use the repository Euclidean convention

    I0'[2sigma g]=-(8piG)^-1 integral sigma e0 dvol,
    e0=E0_trace=R+U(r).

I0 is the accepted covariant spherical EH/Maxwell/reference-potential
theory. Maxwell has zero trace and the reference adds a function of the
areal radius only. Hence U has no metric derivatives. Since the analytic
seed g0 solves I0 at every exterior radius, its known scalar curvature
fixes this function exactly:

    U(r)=-R0(r)=12lambda-(d/B^2)(4u^5-6u^6), u=B/r.       (1)

This is a statement within the accepted covariant spherical theory, not
an invented four-dimensional extension of the reference potential.
The radius is its geometric sphere-area field.

The actual Dirac trace coefficient of Box R is epsilon/(8piG), where
epsilon=G/(60pi). In this Euclidean convention it is represented by

    Gamma_d=c_R integral R^2 dvol,
    c_R=1/(5760pi^2).

Indeed delta_sigma integral R^2=-12 integral sigma Box R (with its
actual boundary completion), so delta_sigma Gamma_d=-(480pi^2)^-1
integral sigma Box R. This agrees with the full trace normalization.

Put a=8piG c_R=G/(720pi)=epsilon/12. For bookkeeping ell define
the explicit field map

    Phi_ell(gbar)=exp[2ell a e0(gbar)] gbar.                (2)

Its derivative at ell=0 changes I0 by

    d_ell I0(Phi_ell(gbar))|0=-c_R integral e0(gbar)^2.

Therefore the coefficient through order ell of the PULLED-BACK action
I0+ell Gamma_rest is precisely

    Ired,ell=I0+ell[Gamma_rest-c_R integral e0^2].          (3)

The coefficient is dictated by the selected Box R anomaly; no new
arbitrary finite coupling is introduced. Equation (3) is the truncated
transformed action, not an assertion that the original action itself
has been changed by a permissible finite counterterm.

Its local curvature content loses the R^2 term:

    c_R[R^2-(R+U)^2]=-2c_R U R-c_R U^2.                  (4)

These are determined lower-derivative spherical dilaton couplings.
The same fixed U, B, A, lambda, flux and common scale are retained.
All other local anomaly terms and the complete nonlocal determinant
are kept. The massless four-dimensional Dirac traceless determinant
is unchanged under the conformal map; its actual Wess-Zumino change
belongs to the higher-order remainder. The paired HH domain is
transported conformally with the original time period. Maxwell is
exactly conformally invariant. No Fourier projector is being applied.

Because e0(g0)=0, the map fixes g0, and the subtracted functional has
zero first derivative there. Thus Ired,ell has the SAME seed forcing
ell Gamma_rest'(g0), the same product preparation and the same selected
vacuum to first order; no first-order observable is silently changed.
At the outer NORMAL/ASYMPTOTIC selected vacuum, the leading e0 and
Gamma_rest' both vanish, so the same vacuum remains an exact leading
stationary asymptotic datum of (3). This is not a new exact GLOBAL AdS
solution with the retained magnetic flux and thermal state. For example
pure AdS at finite r has e0=-Z_ref(r), and the magnetic tail is nonzero.
Smooth static conformal maps with asymptotic factor one preserve bolt period
and laboratory clock; radius observables must nevertheless be pulled
back through (2), rather than identified at higher order by fiat.

For clarity the observable map is explicit. Put s=ell a e0(gbar),
R=e^s r, and assume |r s_r|<1 so this areal change is invertible. Then

    F_new(R)=Fbar(r)(1+r s_r)^2,
    alpha_new(R)=alphabar(r)+s-log(1+r s_r),
    R_h=e^{s_h}r_h, kappa_new=kappa_bar.                 (4a)

The last equality follows by evaluating F_new,R at Fbar(r_h)=0:
F_new,R=e^{-s_h}(1+r_h s_r,h)Fbar_r while the lapse factor supplies
the exact inverse. For the admissible fast tail e0=O(r^-4 log^p r)
with differentiated bounds, (4a) preserves the coefficient lambda,
the constant1, the mass coefficient -2m/r, alpha_infinity, and the
absence of an r^-3 lapse term. Thus its finite geometric mass and
prepared temperature are unchanged exactly. Its magnetic two-form
and time period are not rescaled. Finite-radius lapse/radius errors
are given by (4a), not declared zero; for |r s_r|<=1/2 they satisfy
|alpha_new-alphabar|<=|s|+2|r s_r| and
|R_h/r_h-1|<=exp(|s_h|)|s_h|. These require the stated branch/jet
bounds on the actual reduced metric, not just on the seed.

## 2. Exact scalar reduction away from a stationary metric

The cancellation is an identity of the reduced action, not only a
constant-curvature symbol calculation. Write J0=c_R integral e0^2,
e0=R+U(r), H_U=2U+rU'. For delta g=2sigma g and delta r=r sigma,

    delta e0=-6Box sigma+(rU'-2R)sigma,
    delta dvol=4sigma dvol,
    delta J0=c_R integral sigma[-12Box e0+2H_U e0]dvol.

The last identity follows by integration by parts with the same boundary
completion and cutoff. In the action convention of section 1 the trace
of J0 is epsilon Box e0-(epsilon/6)H_U e0. At ell=1 the exact reduced trace is

    Ered=E_full-epsilon Box e0+(epsilon/6)H_U e0
     =R+12lambda0+epsilon[(3/2)C^2-(11/12)Euler4-20M_em
                           -Box U+(H_U/6)(R+U)].       (12)

Here Euler4=C^2-2J_Ric^2+R^2/6; J_Ric denotes trace-free Ricci.
For general bookkeeping ell, add (1-ell)(U-12lambda0) to the
classical part R+12lambda0 and multiply the displayed epsilon bracket
by ell. This restores the reference potential and gives R+U at ell=0.
The Box R terms cancel exactly at every ell. Since
Box U=U'Box r+U''|grad r|^2, the remaining scalar expression has at
most two geometric derivatives. The Weyl-invariant nonlocal Dirac
determinant has zero conformal trace identically and is retained in
the shape equation.

For unrestricted conformal fields the Euler term contains quadratic
Hessian combinations. The retained static radial sector admits the
sharper quasilinear statement below. The independent shape equation
still has its full fourth/logarithmic order. Conformal variation gives
the reduced action Hessian principal part div(C_r grad)/(8pi G), with

    C_r^{ab}=[6-(11epsilon/6)R+2epsilon(U+rU')]g^{ab}
                                        +(22epsilon/3)J_Ric^{ab},
    div C_r=2epsilon grad(U+rU').                     (13)

For example delta Euler4=-4sigma Euler4+8J_Ric:Hess sigma
-2R Box sigma, while delta C^2=-4sigma C^2. The principal contributions
from -Box U and (H_U/6)e0 are -(rU'+H_U)Box sigma.
These give (13), and div J_Ric=grad R/4 cancels all grad R terms
in its divergence. Thus there is no hidden third metric derivative
in the first-order coefficient. The fixed-volume density contact is
zero order and does not alter this conclusion.

For the actual radial problem write

    g=d rho^2+a_orb(rho)^2 dvarphi^2+r(rho)^2 dOmega^2,
    G_rhorho=2a_orb' r'/(a_orb r)+(r'^2-1)/r^2,
    C_r^{rhorho}=6+(22epsilon/3)G_rhorho
                                      +2epsilon(U+rU'). (13a)

The identity follows from J_Ric=Ric-Rg/4 and G=Ric-Rg/2.
This radial Einstein component uses only first metric jets. Equivalently,
the radial sectional curvatures obey

    K01=-a_orb''/a_orb, K02=K03=-r''/r,
    K12=K13=-a_orb' r'/(a_orb r), K23=(1-r'^2)/r^2,
    Euler4=8(K01 K23+2K02 K12).

The Euler term is linear in a_orb'' and r''. At fixed shape, Weyl
covariance removes every conformal derivative from C^2. Thus the
static radial scalar equation is quasilinear in its conformal second
derivative, with a coefficient depending on first metric jets. Passing
to the fixed disk coordinate preserves this property. This does not
remove the independent shape dependence or prove its nonlinear map.

On the concrete curvature neighborhood

    |R|,|J_Ric|<=100lambda, |U|<=13lambda, |rU'|<=lambda,

one has |C_r-6g|<1000epsilon lambda and C_r>=5g, since
epsilon lambda<2000^-10/180. This is an explicit principal ellipticity
reserve. It does not establish stability of the nonlocal TT form or a
nonlinear inverse. The local scalar expression is defined on C2 geometry;
the existing full C5 source phase can formulate the reduced equations
directly, without first applying the curvature-dependent field map.
Higher regularity is still needed for the physical pullback comparison.

## 3. Exact reduced Hessian at the seed

Let H0(k)=delta e0[k]. At e0(g0)=0,

    (c_R integral e0^2)''[k1,k2]=2c_R integral H0(k1)H0(k2). (5)

All measure, acceleration and off-shell chart contacts of this
subtraction vanish. The original action's contacts are not discarded.
For a pure conformal direction 2sigma g, with T=-Box,

    H_sigma sigma=6T sigma+(2U+rU')sigma
             =6(T+V0)sigma,
    V0=4lambda+v,
    v=(d/B^2)(2u^5-4u^6).                              (6)

The identity follows from delta R=-6Box sigma-2R sigma and
delta r=r sigma, using R0=-U. H_sigma is selfadjoint in the fixed
physical measure. For the genuine TT tensor
h=diag(-b-2c,b,c,c),

    H_c c=-J:h+(rU'/2)c.                               (7)

Although the actual disk direction is h-L_eta g, scalar covariance
gives delta e0[L_eta g]=eta.grad e0=0 at g0. Thus (7) is exact for
the returned direction too. In particular H_c has order zero through
the actual global York H map; no derivative of eta or a new density
contact has been omitted.

For the original fixed-density Hessian blocks (A,B;B*,D), equation (5)
gives the complete reduced blocks, at ell=1,

    Ar=A-2c_R H_sigma*H_sigma,
    Br=B-2c_R H_sigma*H_c,
    Dr=D-2c_R H_c*H_c.                                 (8)

This explicitly changes the TT and mixed blocks. Importing D unchanged
would be wrong. The return contact involving the full nonzero original
trace defect is still present in B. The original scalar density4E
contact is still present in A. Only contacts of the subtracted e0^2
functional vanish, by e0(g0)=0.

## 4. The actual TT gap survives

The accepted seed estimates are

    |J|<=7lambda u^4, |rU'|<=56(|d|/x)lambda u^5,
    ||h||^2<=176piG E0/lambda, D>=.99E0,
    delta0=|d|/x<1.3e-22.

Since |c|<=|h| and 28delta0<1, (7) gives
|H_c c|<=8lambda|h|. Therefore

    2c_R||H_c c||^2
    <=128c_R lambda^2||h||^2
    <=[22528/(5760pi)]G lambda E0<2G lambda E0.

Consequently

    Dr>=(.99-2G lambda)E0>.989E0.                       (9)

H:L2(c)->L2(tensor) is the accepted bounded York reconstruction.
Thus H_c*H_c is an actual bounded L2 operator, not an assumed
negative-Sobolev extension. The bounded-perturbation theorem makes Dr
selfadjoint on the SAME critical operator domain as D. Its positive
spectral form is the original form minus this bounded term. It has an
everywhere-defined positive critical L2 inverse.

On the mass-extended fast domain the perturbation is relatively compact:
J=O(z^4), rU'=O(z^5), the established York maps and their transposes
preserve the stated weights with at most their retained z^4 tail, and
local order zero is below the fourth/log TT order. A mass input z^3
is sent into the fast output. The original TT normal operator is
unchanged. Hence the fast Dr map still has index zero. Any fast/mass
kernel is critical L2 and vanishes by (9). The same critical-to-fast
bootstrap as the actual TT inverse gives its fast surjectivity and
identifies it with its critical inverse on fast sources. These are
uses of the accepted sharp analytic-seed York/parametrix framework;
they do not claim such sharp maps on an arbitrary finite C5 metric.

## 5. Negative scalar coercivity on the nonconstant seed

Write Mstar^2=1/epsilon+11lambda/3 and Z=R+12lambda.
The exact original A is

    A=(80pi^2)^-1[(T+4lambda)(T-Mstar^2)+W],
    W=div((Zg/36+11J/9)grad)-Z/(3epsilon).

Since 72c_R=1/(80pi^2), (6),(8) yield

    80pi^2 Ar=-(Mstar^2+4lambda)(T+4lambda)+W
       -[(T+4lambda)v+v(T+4lambda)+v^2].             (10)

This includes the exact lower-curvature coefficients; it does not assume
constant scalar curvature. Put delta0=|d|/x. Direct differentiation gives

    |v|<=6delta0 lambda, |grad v|<=80delta0 lambda^(3/2).

Equivalently, integrate all multiplier cross terms exactly. For radial
scalars put a_s=Z/36+11j/9. The two coefficients of -80pi^2 Ar are

    gradient: Mstar^2-4lambda+2V0+a_s,
    potential: 4lambda Mstar^2+V0^2+Z/(3epsilon)-Box V0. (10a)

The sign of the last term follows from
<sigma,[TV0+V0T]sigma>=integral[2V0|grad sigma|^2
                                      -(Box V0)sigma^2].
There has been no commutation of T with the nonconstant V0.
The explicit power coefficients give |r v_r|<=34delta0 lambda and
|r^2 v_rr|<=228delta0 lambda. With F/r^2<=5lambda and
F_r/r<=22lambda this yields |Box V0|<=2228delta0 lambda^2
<2300delta0 lambda^2. Also |a_s|<=9lambda and
|Z|<=4delta0 lambda. The gradient coefficient is at least
Mstar^2-5lambda-12delta0 lambda and the potential at least
4lambda Mstar^2-4delta0 lambda/(3epsilon)-2300delta0 lambda^2.
Use Mstar^2>=1/epsilon, delta0<1.3e-22 and
epsilon lambda<2000^-10/180. Both exceed one half of their
constant baseline (Mstar^2+4lambda) and
4lambda(Mstar^2+4lambda), respectively. Therefore

    -<sigma,Ar sigma>>=(Mstar^2+4lambda)/(160pi^2)
                 [||grad sigma||^2+4lambda||sigma||^2]. (11)

All block identities (8)-(10) are initially smooth-core identities.
The cancellation (10a) defines the combined Ar closure on scalar H1.
It is not a subtraction of two separately infinite fourth-order forms
on that larger domain. Br similarly extends by its one-integration
form identity, and Dr is the bounded perturbation already described.
The regular scalar form domain handles the smooth disk bolt; no
unweighted proper-coordinate H1 assertion at its critical endpoint is
substituted. At infinity the reduced scalar normal roots are precisely
-1 and4 (in z powers). The scalar oscillatory pair has been removed
FROM THE REDUCED EQUATIONS. Logarithms resonant at4 remain inside the
auxiliary weight delta=7/2; no log-free endpoint assertion is made.

## 6. The coupled linear inverse

Br remains a bounded pairing between the scalar H1 energy and the
actual TT form domain. For the extra subtraction, integrate H_sigma
once against H_c. Besides |H_c|<=8lambda|h|, the explicit seed bounds
give, for example,

    |grad H_c|<=8lambda|grad h|+101lambda^(3/2)|h|,
    |V0|<=5lambda.

Here differentiate J and rU' as well; |grad(rU')|<=800delta0
lambda^(3/2). The resulting extra mixed form is bounded by
12c_R integral[8lambda|grad sigma||grad h|
    +101lambda^(3/2)|grad sigma||h|+40lambda^2|sigma||h|].
The accepted tensor/return energy estimates therefore apply. This
constructs the weak pairing, without applying an undefined inverse to
a negative-Sobolev source.

Define Dr's dual-form inverse by Riesz/Lax--Milgram from its genuine
positive spectral form. It agrees with its actual critical L2 inverse
when the forcing is L2. The reduced scalar Schur form is

    Sr=Ar-Br Dr^-1 Br*.

The second term is positive. Hence -Sr is coercive by (11), without
any assumption that the mixed term is small. Lax--Milgram solves the
scalar weak equation, followed by the TT reconstruction; elliptic
regularity identifies the smooth solutions with the actual reduced
operator. This is a bounded saddle inverse in its declared energy
source/response spaces.


There is also an explicit inverse bound in a stated critical-energy
norm. Define scalar and TT lower-energy norms by

    a(sigma)^2=(Mstar^2+4lambda)/(160pi^2)
                            [||grad sigma||^2+4lambda||sigma||^2],
    d(c)^2=E0[c],
    p=||f_tau||_(a-dual), q=||f_c||_(E0-dual).

For Ar sigma+Br c=f_tau and Br* sigma+Dr c=f_c, test the equations
against (-sigma,c). The mixed pair cancels in the REAL part, so (9),(11)
give

    a^2+.989d^2<=p a+q d,
    sqrt(a^2+d^2)<=(1000/989)sqrt(p^2+q^2).           (14)

Existence follows equivalently from this signature-flipped
Lax--Milgram form on the ACTUAL closed -Ar and Dr form domains.
The real part is the sum of those positive forms; boundedness of Br
was established above. Sources bounded by p,q are continuous on these
domains. The higher positive fourth/log Dr form is not replaced by
bare E0/H1. Equation (14) is a lower-energy output estimate for the
actual inverse, not a claim that Dr acts boundedly on all of E0.
It requires no small mixed coefficient. It is not yet a numerical
fast/cold C5 bound, a bolt point-value bound or an observable error.

For the fast statement one must use the two different orders: scalar
order2, TT order4/log. The coercive scalar form first gives an actual
critical L2 inverse, and regular-disk/global elliptic regularity gives
its H2 operator domain. At the analytic end Ar has rational coefficient
errors O(z^2) and normal roots -1,4. For a fast forcing the normal
inverse and two-power improvement shift the critical solution from
weight3/2 to5/2 and then7/2, crossing no root; the compact interior is
handled by ordinary elliptic regularity. Thus the actual fast scalar
inverse is onto and unique and has index0, rather than being inferred
from a horizon/outer count alone. Dr is the mass-extended TT
isomorphism of section 4.
The Schur correction has interior order0/log (Br is order2 on each
side, Dr^-1 order-4/log), below the scalar order2, and its coefficients
have the established positive outer decay and z^4 York tails. It is
relatively compact on the scalar fast domain. Thus the reduced fast
Schur map has index0. Its kernel and adjoint kernel vanish by (11)
and the critical form realization. It has a bounded fast inverse;
the full reduced mass-extended map has the same property.

This last assertion uses the same ordinary global parametrix/finite
asymptotic bootstrap framework as the accepted fixed-seed weighted
map, with scalar order reduced to2; it is not a new numerical fast
inverse constant. The regular bolt, time period and original source
constraints remain in the domain. There is no new physical horizon
datum available to tune an original full-theory moment.

## 7. EFT scope and the finite remainder

A bounded invertible LINEAR input/output congruence preserves the
Fredholm index. The displayed original fast mass-extended seed map has
index-1, while the displayed reduced seed map has index0. These two
maps therefore cannot be related by such a linear congruence. This
does NOT rule out all nonlinear field transformations: g0 is not a
stationary point of the full action, and an exact nonlinear pullback
Hessian also contains I' times Phi'', as well as an output derivative
contact. Compactness/domain claims for those contacts have not been
proved here. The specific action (3) is not an exact pullback identity
because it retains only its coefficient through ell; its explicit
Taylor remainder below is the relevant distinction. Its solutions
are EFT predictions, not an automatically selected subset of exact
full-action solutions. The original adjoint moment remains unevaluated.

At the AdS linear-symbol level a bounded invertible LINEAR congruence
cannot erase the root T=Mstar^2 of
(T+4lambda)(T-Mstar^2). The truncated field map cancels its fourth
derivative coefficient only through EFT order. Terms formally beyond
that order are order one at Planck momentum. No analyticity or
all-orders convergence hypothesis repairs that fact automatically.

There is an exact finite-cut identity with an integral Taylor remainder,
using the SAME common cutoff and subtraction on both sides:

    I_ell(Phi_ell(gbar))=Ired,ell(gbar)+ell^2 R_ell(gbar).

Equivalently it is an identity for the action DIFFERENCE after removing
the zeroth and first Taylor terms. The fast-tail remainder displayed
below is integrable. No separate finiteness of the I0/reference volume
actions is asserted; those pieces need not be independently finite.
For the EH part, writing s=ell a e0(gbar), its integrated conformal
formula is

    I_EH(e^{2s}gbar)=-(16piG)^-1 integral
    [e^{2s}Rbar+6e^{2s}|grad s|^2+6lambda0 e^{4s}]dvolbar,

with the prescribed boundary completion. After its linear term the
absolute remainder, for |s|<=delta, is bounded by
(16piG)^-1 integral[(2e^{2delta}|Rbar|+48lambda0e^{4delta})s^2
                                      +6e^{2delta}|grad s|^2].
The reference potential has its ordinary explicit Taylor remainder.
The full Dirac remainder is its exact conformal Wess-Zumino functional;
the traceless nonlocal determinant is not approximated or projected.

This is an action-value identity, not yet an observable error theorem.
For a general C5 endpoint the curvature-dependent map (2) produces
only C3 metric coefficients. Applying the accepted C5 source response
after it needs a C7 branch or a new lower-regularity theorem. Moreover
the exact finite Wess-Zumino remainder contains (Box s)^2; with
s=a(R+U) this is not bounded as a first-variation functional in the
reduced energy dual solely by the already published C5 norm. Such
higher-order terms cannot be dropped in a claimed uniform error
bound merely because their coefficient is small.

A complete controlled-reduction theorem must quantify this remainder
in an appropriate source dual, construct a reduced nonlinear solution
in a sufficiently regular admissible phase, and use its reduced
inverse/stability estimate to bound the requested physical observables
AFTER pulling back (2). The mass/clock map has no asymptotic source
change if e0 has the required fast conormal tail, but that tail and
the moving-bolt radius correction must be checked for the actual
reduced solution. This is a concrete route with a favorable seed
linear inverse; it is not yet the required nonlinear/error theorem.

## 8. The concrete nonlinear problem that remains

The [global disk gauge](stationary-disk-and-normal-inverse.md) is

    g=e^(2Phi)[e^(2V)(d xi^2+xi^2 dvarphi^2)+B^2 dOmega^2].

One may use the fixed analytic-seed linear York map Rcal0 to define
(Phi-Phi0,V-V0)=Rcal0(tau,c), keeping its mass lift. The tau column
changes only Phi. Hence this column is a pure conformal direction
at every metric in this chart; the Weyl-invariant determinant depends
on the shape column only. No moving York projector is differentiated.
The original fixed-period paired Dirac domain and common angular
subtraction remain part of the nonlinear problem.

The reduced orders suggest scalar H_b^(s+2) and TT H_b^(s+4,log),
with the one mass lift, common weight7/2 and H_b^s source norms,
patched to genuine smooth disk norms at the bolt. These are candidate
nonlinear spaces, not an already proved C1 map for the actual
renormalized determinant family. The physical measure is z^-3 times
the b measure dz/z; these weights must not be confused with weights
in physical L2. Interior high regularity and weighted conormal outer
jets must both be tracked under the physical field map.

To construct a reduced nonlinear branch, one still needs an actual
fast/cold inverse constant M_s, a nonlinear derivative bound L_s and
a complete source defect epsilon_s in these matched norms. A
same-level Newton certificate, if applicable, requires
2M_s^2 L_s epsilon_s<1 with radius2M_s epsilon_s inside the admitted
phase. Existing source-value or integrated-constraint residual bounds
are not these three constants. The numerical energy estimate (14)
is a useful input but does not by itself preserve a cold horizon
weight or give a pointwise lapse/radius bound.

A reduced solution must then have sufficient regularity to apply
the exact field map and control the finite action-EQUATION remainder
in an observable dual. An action-value remainder or a small equation
defect alone does not supply that physical error theorem. This is
the next constructive problem; no impossibility result is inferred
from the remaining estimates.

## 9. Primary-source scope

Flanagan--Wald, gr-qc/9602052v2, IV C--D pp22--26, distinguish selecting
exact solutions from modifying equations and exhibit failure of analytic
parameter expansions even for Gaussian forcing. Their treatment of a
nonlocal equation uses a special Einstein-tensor structure and does not
license an arbitrary nonlocal reduction here:
[primary paper](https://arxiv.org/pdf/gr-qc/9602052)

Arrechea et al., 2212.09375v2 IV, evaluate a residual of the static
order-reduced AHS scalar/Boulware approximation. That comparison is not
a full-Dirac magnetic AdS HH observable-error theorem:
[primary paper](https://arxiv.org/html/2212.09375v2#S4)

## 10. Verification scope and physical limits

The [focused checker](../checks/stationary_covariant_reduced_inverse.py),
[receipt](../evidence/stationary-covariant-reduced-inverse.json) and
[raw log](../evidence/stationary-covariant-reduced-inverse.log) audit
the finite field-map, trace, block, coercivity-budget and observable
identities. Closed form domains, weighted elliptic bootstrap and
nonlinear mapping limitations remain ordinary analytic arguments.
No finite algebra check constructs the missing nonlinear branch,
evaluates the original signed moment or supplies physical error bars.

The correction to the index argument in section 7 is substantive:
the full-action seed is nonstationary, so a nonlinear pullback Hessian
has gradient contacts. A bare linear congruence argument cannot prove
an impossibility result for every nonlinear reparametrization.
The actual reduced inverse and the stated EFT truncation retain
their distinct scopes. No original physical gate is closed here.

Remaining to-do list: prove the parameter-dependent mixed map and
quantitative fast/cold nonlinear inverse, construct the reduced branch
and bound physical observable errors; separately decide the original
signed source condition; match exterior/state and full charge,
transfer S138, and complete driven stress, radius, controller,
capture, tail, spatial area and every S074-S080 gate.
