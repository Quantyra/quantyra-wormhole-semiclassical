# Combine the product contact and bound the dominant joint projection

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous turn: progress. The full stable-wormhole objective stays active
and incomplete. The feasibility audit stays closed. This follows the
changed first approximation h_tilde=h1/F_lambda, with the SAME full
finite-seed Hessian and physical normalization. It does not replace the
complete quantum response by a finite-mode or local approximation.

## 1. Eliminate the shared area contact with the exact product equation

Use the affine chart, corrected trace sign and original-G reference in
[the preceding source test](current-reference-rescaled-source-test.md).
The physical residual is r_new=J''h1/F_lambda, where I=F_lambda I0+J.
Write J in I/W units as kappa_area integral exp(2b) dz+L[phi,b].
Extract ALL terms proportional to this same area functional into
kappa_area: this includes the common angular a0 coefficient and the
constant reference/charge difference. Keep every other term in L.

At the exact homogeneous product phi=b=0, R_gamma=-2, the reference I0
is stationary and the full product equations imply the shape equation

    2 kappa_area+g_*=0,
    J=L-(g_*/2) integral exp(2b) dz,
    J''[h1,v]=L''[h1,v]-2g_* integral b_h b_v dz.       (1)

Here g_* is the normalized product shape-gradient DENSITY of L. It
includes its full quantum remainder, not just its heat polynomial.
B,C,mu0 and this product density are held fixed during the metric
variation. No derivative of a moving product is inserted. Equation (1)
is a repartition of the same action, not a finite-coefficient choice.
It cancels the a0 uncertainty exactly: Ha0-2Ga0,* times the multiplication
pairing has coefficient c0-2(c0/2)=0. The curvature-dependent product
quantum response remains. It is not permissible to set g_* to zero.

The covariant reference potential gives the nonconstant density

    -(F_lambda d/C)(4 exp(-phi)-3 exp(-2phi)) exp(2b).

The magnetic part of the existing Weyl action is
-10 S beta phi exp(2b), where beta=x/C, x=lambda B^2 and
S=Gq^2/(60pi lambda B^4). The exact product TRACE equation is
d/x=-5S+eta P_*, eta=G lambda/(120pi); P_* is the retained product
curvature polynomial. Thus separate the following EXACTLY SPECIFIED
part, adding its difference back into the remainder:

    J_M=k_M integral g(phi) exp(2b) dz,
    k_M=5 S beta,
    g(phi)=4exp(-phi)-3exp(-2phi)-1-2phi.             (2)

In particular g(0)=g'(0)=0. The difference between -F_lambda d/C
and k_M is retained in J_rest, as are I_sq, the nonmagnetic Weyl terms,
LLL, a2, the complete massive heat remainder and their product contact
from (1). We do NOT approximate d by -5Sx or omit a quantum constant.
The identity J=J_M+J_rest defines the remainder exactly.

The current selected projection uses the already controlled global
scalar inverse, v=z^4/2+A^-1 B xi, xi=z^4 and b_v=-z^4. In this
execution use its fixed fractional trial v_p and charge its full-source
pairing error at the end. No scalar inverse refinement is selected.

## 2. Keep the actual first correction's radial dependence

Let t=r/B, z=1/t, v_h=1+h and L(z)=-log(1-z)-z. On the actual seed,
the complete massive source has |Theta'(t)|<=Kq^2, K=1e15, and

    a_mass,t=-(G/B^2) J(t)/(t F(t)^2),
    J(t)=integral_v_h^t F(s)Theta'(s)/s^2 ds.

F is increasing, hence |J(t)|<=Kq^2 F(t)(t-v_h)/(v_h t).
Using F>=.99 x(t-v_h)(t-1) and integrating from the fixed outer clock,

    |a_mass(z)| <= A L(z),
    A=G Kq^2/(.99 lambda B^4 v_h)
      <=G Kq^2/(.99 lambda B^4)<3.768e-12.           (3)

The retained LLL, Weyl and added R2 lapse contributions total <4e-28.
Thus |a|<=A0 L(z)+4e-28, A0=3.768e-12. The old global 3.3e-10
bound paid the entire throat logarithm even for a bulk observation.
Equation (3) integrates the actual constraints; no quantum-source
constant or physical norm has been changed.

Use the existing complete Qrho<1.3e-10 and horizon bound
|G rho(R)/(kappa R)|<2e-12. Put q_rho=1.3e-10. Then

    |G rho(Bt)|/x <=1.4e-11 h+(q_rho/2)(t^2-v_h^2),
    f(t)=-[2kappa B v_h(s_R/B)+2G integral_v_h^t rho(Bs)ds]/t,
    fbar=f+(s_R/B)F_t, |s_R|/B<3.33e-10 h.          (4)

The factor B in the dimensional radial integral has already canceled
in this t formula. We use it only for t-1>=1e-15. On this region the
exact factored thermal seed gives

    F(t) >=(1-1e-27)x (t-1)^2(t^2+2t+3)/t^2.       (5)

Indeed the negative d term has relative magnitude <4e-29. The thermal
term has relative magnitude <=7h^2/(t-1)^2<3e-36 using h<6e-34.
For the q_rho part of (4), integration from 1 instead of v_h enlarges
the bound. Dividing by (5) gives

    |f/F|_density <=q_rho t(t+2)/[3(1-1e-27)(t^2+2t+3)]
                   <=q_rho/[3(1-1e-27)].            (6)

The horizon constant contributes <3e-29, the bolt value <1e-44, and
the translation term <5e-27. For the last bound, F''<=21lambda and
F_t(v_h)=2kappa B give F_t/F<22/(t-1) on this region. These ample
explicit reserves imply

    |fbar/F|<4.34e-11,
    |H(z)|=|a+fbar/F|<=A0 L(z)+4.34e-11+4e-28.      (7)

The exact fixed-chart displacement is w=xi_disp/B with
xi_disp=s_R-F integral_r^infinity H/F. The previously proved
F integral_r^infinity dr'/F<=2(r-B), together with monotonic L, gives

    |w(z)|<=s0+2(1-z)[A0 L(z)+4.34e-11+4e-28]/z,
    s0=3.33e-10 h.                                  (8)

Every integration point in this outward integral is in the same bulk
region. Equations (3)-(8) are bounds on the ACTUAL full-source h1.
They make no assertion that h1 equals its magnetic or local surrogate,
and do not establish a C5 inverse or a nonlinear admission bound.

## 3. The joint weak kernel removes displacement derivatives

At the seed phi=-log z write

    g=4z-3z^2-1+2log z,
    g_phi=-4z+6z^2-2=-2(1-z)(1+3z),
    g_phiphi=4z-12z^2.

The actual h1 has nu_h=z w and b_h=a/2-z w-z^2 w'/2. Differentiate
(2) in both fields before any absolute values. Its mixed density is
k_M[k_nu nu_h+k_b b_h], with

    k_nu=g_phiphi v_p-2g_phi z^4,
    k_b=2g_phi v_p-4g z^4.

Integrate w' once. The g_phiphi and g_phi terms cancel algebraically,
leaving the derivative-free weak expression on 0<z<z_c:

    J_M''[h1,(v_p,-z^4)]/W
       =k_M {integral [a K_a+w K_w] dz-[z^2 k_b w/2]_0^z_c},
    K_a=g_phi v_p-2g z^4,
    K_w=z^2 g_phi v_p'-8g z^5.                       (9)

The outer endpoint vanishes in the actual mass/log phase; the bound (8)
already suffices. The finite z_c endpoint is RETAINED, not discarded.
Use z_c=1-2^-49, so s_c=(1-z_c)/z_c<2e-15 and s_c>1e-15.
Define

    I_a=integral_0^z_c |K_a| L dz,
    I_w=integral_0^z_c (1-z)|K_w|/z dz,
    I_wL=integral_0^z_c (1-z)|K_w| L/z dz.            (10)

Equations (3),(8) then bound the bulk contribution by

    |P_M| <= (k_M/F_lambda){ A0(I_a+2I_wL)
       +2(4.34e-11+4e-28)I_w
       +4e-28 integral|K_a|+s0 integral|K_w|
       +(z_c^2/2)|k_b(z_c)| |w(z_c)|}.              (11)

Freeze the target |P_M|<1.31e-40. The exploratory use of uniform
3.3e-10 and 1.08e-9 gave 3.34e-39 for this component. That failed
UPPER allowance was not an observed source value. Equation (11) tests
the correlated weak functional with the actual radial envelopes.

For rigorous integration, write v_p=P(z)+(1-z)^nu T(z), with the
previous fixed coefficients and nu. Both v_p and (1-z)v_p' have this
form. On each cell use a degree-six Taylor enclosure for (1-z)^nu,
its explicit seventh-derivative remainder, and Bernstein bounds for
the combined polynomial before taking absolute values. Bound g using
its exact increasing endpoint values. Subdivide until each integral's
upper-minus-lower enclosure is <.0005. This avoids treating correlated
large polynomial coefficients as independent function values.

The interval 0<z<2^-20 is covered analytically by absolute coefficient
sums, |g|<=8+2|log z| and L<=z^2/[2(1-z)]. All its contributions
are included. There is no outer wall. Positivity and the old v_p
monotonicity certificate control its bolt value; no new scalar solve
or endpoint interpolation is made.

## 4. Bound the complete omitted bolt contribution

On h<=s<=s_c the already proved FULL action defect satisfies
max_i |r_old,i|<1.73e-15 a_s, a_s=s[1+log(.01/s)]. The full seed
shape row satisfies |F_b|<8.9e-11 s, and its scalar row is <1.5e-27 s.
The fixed trial has |v_p|<1.01 and |b_v|<=1. Therefore

    |P_cold| <=1/F_lambda {2.01(1.73e-15)(s_c^2/2)
                       [3/2+log(.01/s_c)]
                +epsilon_R(8.9e-11+1.01(1.5e-27))s_c^2/2}.       (12)

Here dz=ds/(1+s)^2<=ds and the positive upper envelope was extended
to zero only inside the bound. The actual temperature, domain and bolt
are not replaced. This controls BOTH source rows and the full quantum
tower in the omitted region. Freeze (12)<2e-43. It is a projection
remainder for this test, not a source or inverse norm for a whole ball.

The cutoff is an observation/integration split, not a Dirac or metric
boundary condition. The nonlocal kernel still lives on the full paired
geometry. A test multiplied by the sharp characteristic of [0,z_c] is
not automatically in the Hessian form domain. A later operator evaluation
must compute the full-domain source row before restricting its observation
integral, or use a smooth partition whose transition is confined to
s<=s_c, retain all cutoff/domain contacts and charge the entire cold
remainder (12). It must not cut the Green problem at z_c or use a
discontinuous metric variation as an admitted test.

Finally the preceding complete trace majorant and seed subtraction give
the GLOBAL scalar-trial pairing error <2.331e-44. Apply it ONCE to the
full r_new, including the cold/bulk split. Thus a sufficient working
screen is

    |P_rest,bulk|<1.287e-40,
    |P_M|<1.31e-40, |P_cold|<2e-43, E_scalar<2.331e-44,
    total <2.60e-40.                                  (13)

This allocation replaces the preceding component allocation; do not add
the earlier trace-only and principal-stress pieces again. P_rest,bulk
is the UNKNOWN signed projection of the exactly defined J_rest in (2),
with its product contact (1), on 0<z<z_c. A missing allocation is not
itself a violation of the necessary uniform-TT ceiling.

## 5. Research decision and continuum boundary of this result

Retain h_tilde and the complete common prescription. Execute (10)-(12)
now. If they pass, the single next calculation is P_rest,bulk with a
certified quantum/contact error sufficient to decide its 1.287e-40 target.
This must include the actual geometric Dirac Hessian, its product
gradient subtraction, full thermal/angular tails and all remaining local
terms. Selected mass-vertex modes are not that operator. Its regularized
trace, tadpole, Weyl conversion and common restoration cannot be replaced
by the four-dimensional trace anomaly.

For context, [Namuduri--Bilal, section 5.1](https://arxiv.org/html/2308.05802v2#S5.SS1)
separates the geometric variation into a local area term, Liouville term
and renormalized Green function. Their smooth compact-surface calculation
does not supply our paired C2,1 seam remainder or authorize a small-mass
expansion at M1>18. No quantitative bound from that paper is imported.
The actual differentiated paired domain, one-sided cubic and weight-four
logs, and uniform removal of the COMMON regulator still precede numerical
interpretation of P_rest,bulk. Any failure there takes priority over more
source optimization. The global graph inverse is not an admission bound.

The full projection and controlled solution remain UNEXECUTED here.
No total-source improvement, nonlinear neighborhood, compatibility,
matching, causal energy balance, semiclassical validity/noise or all-sector
stability follows from these component tests. The full objective is kept.
Freeze the note, checker and review before running; preserve all attempts
and exploratory snapshots. Supporting checks do not replace these ordinary
continuum arguments. No Lean is selected.

Remaining to-do list: complete the joint quantum/contact remainder and
global Schur forward/adjoint solve; physical bolt/exterior certificate,
X/admission, nonlinear closure and exact compatibility; matching and causal
realization, semiclassical validity/noise, all-sector stability, charge/return,
opening, transport, preparation, separate creation and the full stable-wormhole objective.
