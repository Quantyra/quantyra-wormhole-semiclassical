# Reconstruct a prescribed source component before choosing the next solve

2026-10-24, E003/S169. Ordinary argument before acceptance execution.
The full stable semiclassical wormhole objective remains active and incomplete;
the feasibility audit stays closed. This is a component of the actual reference
correction, not the full quantum correction or a freely chosen forcing.

## 1. What this calculation can decide

The previous localized adjoint test requires |I_g|<1.06655e-54 for the
current uniform .0022 certificate, even with the entire bolt error budget.
The complete available allowance, 1.3079691e-40, is not a measured value.
Before interpreting that gap, reconstruct the response to the exactly
specified product-matched magnetic potential J_M. Its coefficient is
k=G q^2/(12 pi C B^2), not an adjustable coupling. Use precisely the
decomposition J=J_M+J_rest of current-product-matched-source-projection.md.

Let h_M solve I0'' h_M=-J_M'(g0), with the same fixed clock, temperature,
finite geometry and reference boundary action. The actual full first
correction decomposes h1=h_M+h_rest. This linear identity does NOT identify
h_M with h1. In particular the full-action residual of h_M/F_lambda is

    f_rest + J'' h_M/F_lambda,

where f_rest=J_rest'(g0). The first-order term f_rest must not be omitted.
The selected full residual projection is

    I_g = P_MM + P_cross,
    P_MM = J_M''[h_M,(v,-z^4)]/(F_lambda W),
    P_cross = {J_M''h_rest + J_rest''h_M + J_rest''h_rest}
                                      [(v,-z^4)]/(F_lambda W).

Here v is the SAME global scalar-return test. No rest source is being
selected to cancel a residual. We execute P_MM, including the existing
scalar-trial error, and give a numerical target for the still unknown
P_cross. Passing this necessary projection screen is not a bolt certificate.

## 2. Exact reconstruction on the actual finite-temperature geometry

Put t=r/B, v_h=1+h, z=1/t, z_h=1/v_h, F(t)=C t^2 Q(1/t), and

    g(t)=4/t-3/t^2-1-2 log(t),
    A(z)=(3-2 log z)z-2z^2+z^3-2,  A'(z)=-g(1/z).

J_M is a covariant radius potential. Its complete orbit source satisfies

    rho_M=-k C g(t)/(2G t^2),  p_M=-rho_M,  a_M=0.

The angular source is the derivative of this same potential, so the
reference Ward identity is retained. Write sigma=delta R/B. The exact
fixed-temperature constraint solution is

    sigma=-k C g(v_h)/[v_h^2(v_h F_tt(v_h)+F_t(v_h))],
    f_M(t)={-2 C gamma v_h sigma
                   +k C[A(1/t)-A(z_h)]}/t,
    w_M(t)=-F(t) integral_t^infinity f_M(s)/F(s)^2 ds,
    nu_M=w_M/t,  b_M=w_M,t/2-w_M/t.                 (1)

Both horizon equations hold: f_M(v_h)+2C gamma sigma=0 and
f_M,t(v_h)+F_tt(v_h)sigma=0. Set fbar=f_M+sigma F_t and H=fbar/F.
Then t fbar and its first derivative vanish at v_h. Thus H(v_h)=0;
the equivalent formula w_M=sigma-F integral_t^infinity H/F has no
logarithmic bolt singularity. It follows that w_M(v_h)=sigma and the
fields are regular radial disk fields. No inner wall is introduced.

At infinity f_M=const/t+O(log(t)/t^2). Consequently nu_M and b_M have
their allowed cubic mass terms followed by z^4 log z terms, with all
five weighted derivatives. They have zero leading clock/scalar source.
This checks the domain of this PARTICULAR tangent; it proves neither
admission of the nonlinear X ball nor a strong inverse for the full source.
The original paired quantum realization is unchanged and is not replaced
by a smooth double or a selection of modes.

## 3. Complete, conservative bounds for this reconstructed tangent

The identities g(1)=g_t(1)=0 and
g_t=-2(t-1)(t+3)/t^3 give |g|<=4(t-1)^2 and
|(g/t^2)_t|<=16. The actual parameter intervals in (1) give
0<sigma/k<3h^2. For N=t fbar,

    |N_tt| <= k C(16+30h^2),  N(v_h)=N_t(v_h)=0.

Use the retained whole-geometry inequalities
F>=.99x(t-v_h)(t-1), |2F_tt+tF_ttt|<10C, x=lambda B^2.
They imply

    |H| <= [(8+15h^2)C/(.99x)] k <49k,
    |w_M| <=3kh^2+98k(t-1),  |nu_M|<99k.         (2)

The second estimate uses the already proved
F(t) integral_t^infinity ds/F(s)<=2(t-1). For t-v_h>=h the independent
geometric bounds F_t integral 1/F<16 and (F/t) integral 1/F<10 give
|b_M|<(1/2+10)49k+3kh^2<515k.
For 0<=y=t-v_h<=h, (1) instead gives |H|<4.1k y/h. The exact convex
lower bound F>=2C gamma y+.95x y^2 gives

    |integral_t^infinity H/F| <= (4.1+2.8*49)k/(2C gamma).

Here log(1+2C gamma/(.95xh))<2.8. Also
F_t<=2C gamma+12.1xh and F<=2C gamma h+6.05xh^2.
Substitution in b_M=H/2-sigma/t-(F_t/2-F/t)integral H/F gives
|b_M|<170k there. Use the global bound |b_M|<600k.
All constants are numerical substitutions into actual geometry bounds,
not a new source normalization. The physical bolt functional is

    |L_b h_M|=sigma/h<3kh.                         (3)

## 4. Certified signed integration, with the endpoints retained

Use the previously frozen fractional scalar trial v_p, with its original
coefficients. Let g now denote g(1/z), g_phi=-2(1-z)(1+3z), and

    K_w=z^2 g_phi v_p'-8g z^5,
    k_b=2g_phi v_p-4g z^4.

Since a_M=0, integration by parts on [0,z_c], z_c=1-2^-49, gives

    P_MM,bulk = k/F_lambda [integral_0^zc w_M K_w dz
                                  -zc^2 k_b(zc)w_M(zc)/2].   (4)

The outer endpoint vanishes in the actual mass/log phase. The finite
z_c endpoint is INCLUDED. On [z_c,z_h] use the ORIGINAL mixed density
k[(g_phiphi v_p-2g_phi z^4)nu_M+(2g_phi v_p-4g z^4)b_M], not a
discontinuous test or a new boundary condition.

The exact nested integral in z is

    w_M/k = Q(z)/z^2 integral_0^z s^3[-A(s)+A(z_h)
                                  +2 gamma v_h sigma/k]/Q(s)^2 ds. (5)

Its integrand is nonnegative. Evaluate Q with its factored polynomial
and the ORIGINAL nonzero thermal term; retain correlated gamma,h,z_h.
Monotone endpoint values bound A and g on each cell. Cumulative lower
and upper rectangle integrals bound (5), including every intermediate
point of a cell. Signed interval products then enclose (4). The scalar
trial and (1-z)v_p' use the earlier degree-six Taylor/Bernstein enclosure,
with its explicit seventh-derivative remainder. This is a whole-interval
quadrature proof, not point sampling.

On 0<z<z_o=2^-20, Q>beta/2 and -A(z)+A(z_h)+2gamma v_h sigma/k<3;
thus the initial integral lies in [0,3z_o^4/beta^2]. The existing
absolute scalar-trial coefficient sum gives the analytic contribution
bound with integral z^4 |log z|. For the last cold interval, (2) gives
|nu_M|/k<=3h^2+98(1-z_c)/z_c; use |g_phi|<=8(1-z_c),
|g|<=4(1-z_c)^2/z_c, |g_phiphi|<=8 and |v_p|<1.01.
These are full endpoint allowances, not omitted cold physics.

Finally transfer from v_p to the actual v. The retained scalar inverse
certificate bounds ||v-v_p||L2(z^-4 dz)<3.6891662e-12. The scalar row of
J_M''h_M is k(g_phiphi nu_M+2g_phi b_M). Since |g_phi|<=3,
|g_phiphi|<=8, (2) and |b_M|<600k give its dual weighted norm
<=4392 k^2/sqrt(5). Therefore add

    E_scalar,M <= (4392/sqrt(5)) k^2(3.6891662e-12)/F_lambda.  (6)

This uses the existing scalar solve; no scalar optimization is selected.

## 5. Decision rule and claims boundary

The first frozen run rejected the bolt shift: independently evaluating
g(z_h) lost its exact double zero and produced an interval wider than
h^2. Preserve that failure. Evaluate the SAME function for z>1/2 as
g=-4y^2-sum_(j=3)^48 2y^j/j, y=1-z, with tail at most
2y^49/[49(1-y)]. Likewise A=-4y^3/3-sum_(j=3)^48
2y^(j+1)/[j(j+1)], with tail at most 2y^50/[49*50*(1-y)].
Both are convergent identities with their entire tails retained. This
preserves the shared bolt parameter; it changes no source or temperature.

Freeze acceptance at a complete |P_MM|<3e-57. Independently, (2) and the
previous certified I_w bound give a conservative check <3e-57 without
the signed quadrature. Compare both with the actual localized ceiling.
If the component is far below that ceiling, the oversized original
allowance is not evidence that THIS source component forces failure.
It still supplies no bound for P_cross or the complete correction.

The decisive outstanding calculation becomes the reconstruction of
h_rest=-I0''^-1 J_rest'(g0) on the actual paired domain, retaining its
mass/log endpoint expansion and full thermal/angular/contact source.
This is the missing input to the coupled directional forward/adjoint
solve, not permission to replace that solve by h_M. First-order source
evaluation and its strong weighted return take priority over optimizing
the already passing cold/scalar/quantum projection bounds.

The known excessive upper allowance, unknown actual P_cross, and proved
failure of a uniform V* to C5 admission bound are three different facts.
None proves physical nonexistence. The qualitative C2 action on bounded
ADMITTED stationary families remains valid in its existing scope.
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 remain unclosed;
Y0<=5e-9, Z0<=.25, Z2<=5e7 remain sufficient unproved allocations.
The actual <1.667e-9 bolt difference, strong admission margins, exact
compatibility, physical matching/causal realization, semiclassical
validity/noise and all-sector stability are not established here.

Remaining to-do list: reconstruct the complete remaining prescribed source
and correction with weighted error; coupled directional forward/adjoint
certificate and nonlinear admission; exact compatibility; physical matching
and causal realization; semiclassical validity/noise and all-sector stability;
charge/return, opening, transport, preparation, separate creation and the
full stable-wormhole objective.
