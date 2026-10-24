# Actual-pair heat variation with the rough direction kept as a smearing

2026-10-23, E003/S169. Ordinary argument before execution.
Previous turn: progress, paired boundary and regulator audit. The full
stable-wormhole objective remains active and incomplete; the feasibility
audit remains closed. Keep the original theory, state, period, scale and
paired matter boundary. Do not optimize the controlled cold/scalar terms.

## 1. Resolve a selected-pair limit without asserting a smooth double

The previous audit correctly distinguishes a candidate form extension from
a uniform limit of regulated metric variations. For the SELECTED projection
there is additional structure: its quantum shape test is the globally
smooth v=-z^4. Only h=b_h1 has the mass/log regularity. Exploit this
asymmetry before requiring high derivatives of h.

The actual seed Q is a polynomial on the physical half, with a smooth
regular disk bolt. It is smooth up to the physical outer boundary; it is
its reflected cubic extension that is only C2,1. Use the actual boundary
problem on this half. Let g_s=exp(2s v)g_0. For sufficiently small s this
is a smooth, uniformly elliptic family with the same paired projector,
induced square condition, period and regular bolt. The multiplier v and
its first three outer jets vanish. In particular K_s=S_s=0 at the seam.
This is an actual allowed variation, not a replacement for h or the seed.

Let L_s=D_s^dagger D_s, Lbar_s=D_s D_s^dagger and K_s(t)=Tr exp(-tL_s),
with physical-half rank two. The exact partner intertwiner in
stationary-paired-dirac-phase.md preserves scalar smearings and gives
equal partner traces. On the fixed spinor Hilbert space the metric vertex
in h is D_h=-{h,D_s}/2. Consequently, for t>0,

    delta_h K_s(t)=-2t partial_t Tr[h exp(-t L_s)],
    delta_v delta_h K_0(t)
       =-2t partial_t partial_s Tr[h exp(-tL_s)]|_0. (1)

To derive (1), vary D^dagger D, cycle the trace-class positive-time
products, and intertwine D through its two squares. The two terms are
t Tr[h L_s exp(-tL_s)] and t Tr[h Lbar_s exp(-tLbar_s)]. This derivation
keeps the domain and both partners. It is not the variation of the raw
variable-mass characteristic without its conversion. Differentiating
the complete expression in s retains the second metric vertex; it is
not an approximation retaining only a two-insertion bubble.

For the actual h, choose physical-side core approximations h_e that
converge in W2,infinity and pointwise, with |h_e(z)|<=H z^3 on a fixed
outer collar and a common bound elsewhere. Such H is finite for the
retained finite log polynomials and their conormal remainder. The cubic
cutoffs tested in the preceding audit have these properties despite their
diverging fourth-derivative norm. Scalar multiplication preserves the
first-order domain; the vanishing first normal jet also preserves the
induced square condition where that is used. At fixed t>0 the heat trace
and its vertices pass to the h_e limit by these form/operator bounds.

## 2. The uniform short-time estimate for this ordering of variations

Use the ordinary smooth strongly elliptic mixed-boundary heat parametrix
on g_s, uniformly for s in a small compact interval. Its interior diagonal
coefficients and their s and t derivatives have a uniform expansion.
The near-boundary part has profiles in r=z/sqrt(t), with Gaussian decay;
the coefficient of order t^(j/2) depends only on boundary metric and
operator jets of weighted order at most j. A boundary profile therefore
cannot respond to v=O(z^4) before j=4. This is also seen by rescaling the
local heat equation: its first changed coefficient is of order t^2.
Keep the coordinate-density and operator vertices in that rescaling.

Thus, after s differentiation, the traced boundary diagonal has bound

    |partial_s K_boundary(t;z,z)|
        <=C t (1+r)^N exp(-r^2/c),                (2)

for some finite C,c,N, with the same bound after t partial_t. Terms beyond
the local collar are exponentially small. The constants depend on the
actual smooth family g_s and its fixed boundary problem, NOT on h_e.
The leading rank/projector cancellation is compatible with (2); no
derivative of the log smearing is needed to use it. Multiplying by
|h_e|<=H z^3 and integrating dz gives a boundary contribution O(t^3),
since t times int z^3 exp(-z^2/(ct))dz is O(t^3). The time-circle integral
is divided by its original period; it is not used as a gain in precision.

In the interior, keep the diagonal expansion through t^2 and bound the
next term and its indicated derivatives by C t^3. Smearing requires only
the common L1 bound of h_e. Formula (1), the bulk/boundary completion and
the common product contact consequently give

    r_hv^prod(t)=a_hv t^2+E_hv(t),
    |E_hv(t)|<=B_hv^(heat) t^3,  0<t<=t0,          (3)

with FINITE constants uniformly in the selected h_e approximation. Here
r is the same time-reduced geometric heat remainder after a0/a1/a2 as
the previous checkpoint. Existence of these constants is a qualitative
parametrix statement on the actual smooth HALF-domain and the smooth v
variation. No numerical value of B_hv^(heat) follows from that theorem.
The fixed positive temperature may enter its constants.

With Theta(t)<=4C/t, (3) is an integrable common-regulator majorant.
For positive times, heat smoothing gives convergence; at large times the
massive angular gap supplies the required decay, including possible
massless orbit kernels. The LLL allocation remains separate. This proves
the common heat-regulator limit of the SELECTED regulated mixed response
and its equality to the limit of the smooth-core responses. It uses the
full thermal heat kernel and complete angular counting, not selected modes.

The smooth-boundary heat construction and local jet weights are reviewed
in [Vassilevich, sections 2, 4.4 and 5](https://arxiv.org/html/hep-th/0306138v3).
They are used only on the genuine smooth family g_s above. The rough h is
handled by (1)-(2), not by applying that theorem to a reflected log metric.

Scope: (3) does NOT establish a C2 action on a two-parameter family
exp(2(t h+s v))g0, a nonlinear neighborhood, arbitrary rough test pairs,
or an operator-norm heat majorant. Those require uniform estimates when
the BASE metric itself has the log regularity. Nor does it finish the
raw thermal ODE-to-geometric representation match. This calculation
closes the selected regulated weak response limit, not those stronger
claims. The remaining projection can be defined by this geometric route;
it still requires a numerical remainder and complete evaluation.

## 3. Calculate the first omitted bulk coefficient on the actual seed

Write the integrated rank-two orbit heat coefficient at t^2 as

    A3[g]=int [|grad R|^2/840-R^3/1008] dvol
                                      +boundary completion. (4)

Reduce the general a6 formula in Vassilevich (4.29) using E=-R/4,
Ric=R g/2 and Omega_12=(R/4)gamma_1 gamma_2. After integration by parts,
the scalar metric, spin connection and E terms give gradient coefficients
-1/140, 1/480 and 1/160, whose sum is 1/840. The cubic coefficients are
1/315 and -1/240, whose sum is -1/1008. The flat sixth Bessel moment
I6=-48/(35pi^2 M^4), with int cos^6(theta)dtheta=5pi/8, independently
gives the Hessian coefficient 1/(840pi M^4).

An a6 boundary invariant has metric jet weight five. A mixed h,v term
would need at least weight 3+4=7 at this boundary. Its limiting mixed
variation is therefore zero, including completion divergences. Log jets
do not invalidate this count: every such product still vanishes as at
least z^2 times a finite log polynomial. At the bolt the disk is smooth
and has no boundary. The same reasoning gives no mixed boundary term at
a7; the boundary remainder at order t^3 is retained, not set to zero.

For b in the fixed affine chart, put Delta=(Q partial_z)', R=-Q'',
S=R-2Delta b. The bulk density in (4), per unit canonical time, is

    Q [partial_z(exp(-2b)S)]^2/840-exp(-4b)S^3/1008. (5)

Let L=Delta+R, v=-z^4. Move derivatives off the actual h. Direct variation
of (5) and integration by parts give int h O6[v] dz, with

    O6[v]=-[L Delta L v+(Delta R)Lv+Delta(v Delta R)]/105
      -[16R^3 v+24R^2 Delta v+24Delta(R^2 v)
                                      +24Delta(R Delta v)]/1008. (6)

All fluxes at the outer seam vanish by h=h'=h''=0 and the v jets above;
the regular bolt fluxes contain Q or its smooth-disk normal derivative.
At the exact product R_*=-2 the first A3 gradient density is -2/63.
Keep its product-gradient contact, so the required kernel is

    O6prod[v]=O6[v]+(4/63)v,
    a_hv=(1/(4pi)) int_0^zh h O6prod[v] dz.         (7)

This is an actual local coefficient of the retained remainder, not a new
counterterm or normalization. It does not replace the full remainder by
its leading coefficient.

## 4. Use the actual first correction and certify the coefficient bound

Use precisely the radial h1 bounds of the magnetic projection milestone:

    h=a/2-z w-z^2 w'/2,
    |a|<=A0 Lz+a_loc,  Lz=-log(1-z)-z,
    |w|<=s0+2(1-z)(A0 Lz+f0+a_loc)/z,
    A0=3.768e-12, a_loc=4e-28, f0=4.34e-11,
    s0=3.33e-10 h_bolt.                            (8)

They hold on the retained bulk 0<=z<=zc=1-2^-49. Keep the actual polynomial
Q, including C, beta, thermal tau and epsilon; do not replace it by the
product or zero-temperature geometry. Formula (7) and the constraint
identity imply on the bulk

    int h O =int [a O/2+(z^2 O'/2)w] -[z^2 O w/2]. (9)

Retain the finite-cut endpoint. For the remaining interval use the already
available global |h|<1.2e-8 and its length <=1-zc. This bounds only the
local coefficient's omitted integration interval, not a new optimization
of the full cold source.

Evaluate |O| and |O'| by outward Bernstein polynomial ranges on a complete
dyadic partition. Integrate the positive weights in (8) exactly, using

    Jk(z)=int_0^z u^k Lz(u)du
       =[(1-z^(k+1))log(1-z)+sum_(j=1)^(k+1) z^j/j]/(k+1)
                                             -z^(k+2)/(k+2). (10)

The two cell weights are A0 Delta J0/2+a_loc Delta z/2 and
A0 Delta(J1-J2)+(f0+a_loc)Delta(z^2/2-z^3/3)+s0 Delta(z^3)/6.
Start with 256 equal bulk cells; refine only if the enclosure width of
the positive integral exceeds 1e-12 after the factor 1/(4pi). Retain all
cells/endpoints and every attempted run. A direct independent quadrature
is a diagnostic of the enclosure, not a replacement for it.

Let U6 be the resulting upper bound for |a_hv|, including both end terms.
With t0=.001 and Dtarget=2.5e-43 B^2 F_lambda/(4GC), the ACTUAL requirement
on the complete heat cap is

    U6 t0+(B_hv^(heat)/2)t0^2 <=Dtarget,
    B_hv^(heat) <=2(Dtarget-U6 t0)/t0^2,            (11)

provided the right side is positive. Execute this test without assuming
it passes. The B coefficient is UNKNOWN; U6 is a known upper bound, not
the signed coefficient. An oversized U6 would fail this estimator only.
Even success of (11)'s numerical capacity test does not certify the cap
until the full remainder constant is supplied.

## 5. Decision boundary and remaining objective

This work must distinguish three results: the ordinary selected-pair
regulator-limit argument (1)-(3); the executed coefficient bound (7)-(11);
and the still unexecuted complete signed projection. Do not count checks
as the proof of the first result or use the second as the third.

If (11) leaves a positive usable capacity, continue the geometric
smearing method. The single next calculation is a VALIDATED bound on the
entire residual heat response after the A3 term, including the boundary
profiles, for 0<t<=.001. Compute B_hv^(heat) against (11), not another
local coefficient or a smooth surrogate. If the coefficient capacity
fails, first test whether its envelope rather than the actual response
is responsible before changing construction.

The full signed target remains |P_hat_rest|+E_rest<1.287e-40, sufficient
only. Thermal/angular/positive-time/local contributions still need their
complete evaluation. Preserve magnetic <1.284506e-40, complete cold
<1.693010e-43 and scalar-trial <2.331e-44. No total-source improvement,
actual corrected solution or nonlinear admission is claimed.

At rho=1e-8 the requirements Y0+1e-8Z0+5e-17Z2<=1e-8 and
Z0+1e-8Z2<1 remain open, as do a C2 action on the actual nonlinear phase,
exact compatibility, physical matching, causal boundary/energy data,
semiclassical validity/noise and all-sector stability.

Remaining to-do list: bound the full residual heat response and evaluate
the complete signed projection; controlled forward/adjoint correction
and admission; nonlinear action/domain closure and exact compatibility;
physical matching, causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
