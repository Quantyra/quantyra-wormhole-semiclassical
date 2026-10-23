# Complete separated static shape response on the finite seed

2026-10-22, E003/S169. Ordinary argument before numerical execution.
The full stable semiclassical wormhole objective remains incomplete.

The selected calculation is one off-diagonal kernel value of the actual
static canonical shape Hessian, summed over every massive angular level
and every antiperiodic thermal frequency. The points are z_A=7/8 and
z_B=1/2, with z=B/r, at the allowed preparation z_i=.265. Separation
makes this observable ultraviolet finite. It tests assembly of the
actual finite-domain response; it is neither the coincident source nor
an operator norm or metric correction. The closed feasibility audit and
its admission targets are unchanged.

## 1. Domain and normalization before optimization

Retain all parameter formulas and their interval enclosure from
[the preceding radial calculation](current-seed-radial-green-enclosure.md).
In particular G=1 is only the unit, k=2000, q=k^4, c_R=1e22,
M_n^2=n(n+q)/C, d_n=2(q+2n), p=gamma(j+1/2), j>=0, and

    Q=(1-z)^2 R-tau z^3, a=sqrt(Q)/gamma,
    R=beta(1+2z+3z^2)+u z^2+epsilon z^4,
    beta=(1-u-epsilon)/6, u=1/C, epsilon=d/C,
    d z/d rho=-sqrt(Q), period(theta)=2pi.

The exact product enclosure and original clock are used, not numerical
pressure quadrature or a product inverse. The physical tau is nonnegative
and <1e-65. The cell enclosure retains that small thermal term.

The domain step is stronger than equality of a single diagonal heat
trace. On the folded four-dimensional pair let Rf=chi P, with P the
outer reflection and chi=tau_3 gamma_n. Reflection anticommutes with
the Dirac expression, so its square preserves the two Rf sectors. The
actual complementary-adjoint square is the plus sector. The parallel
J of [the phase audit](stationary-paired-dirac-phase.md) commutes with
the Dirac expression and anticommutes with Rf. It therefore pairs the
two square sectors, also along identical scalar deformations supported
in the interior. After the same per-exterior normalization, the action
is half the doubled action along even variations. This justifies
differentiation of that identity; differentiating a diagonal identity
known only at a constant mass would not suffice.

The even double has the retained cubic seam and only C2,1 coefficients.
Use the closed first-order/square realization and one-dimensional
transmission at that seam, not a smooth-double heat expansion. Compact
interior variations leave the seam jets and the first-order projectors
fixed. This argument does not extend to unequal mouths, evolving source
angles or the unproved retarded gravitational constraint domain.

In the canonical conformal representation one massive determinant is
Gamma_M[b]=-Tr log(D+M exp(b)), with the usual local Weyl allocation;
see [the first-order vertex derivation](flat-variable-dirac-three-point.md).
At distinct radial points its static Hessian kernel is

    M^2 sum_omega tr[S_omega(A,B) S_omega(B,A)],
    S=(D+M)^-1.

The exponential-mass seagull, local Weyl allocation, local curvature
restoration, reference functional and LLL Polyakov Hessian are supported
on the diagonal for this fixed-chart shape variation. They vanish when
tested on disjoint radial supports. They cannot be dropped in the full
Hessian. Even variations and the half-action factor give a DIRECT plus
REFLECTED kernel, not the direct kernel alone. Negative thermal numbers
give the same spin trace by chirality exchange.

## 2. First-order cancellation inside the kernel

On the double choose D=[[0,partial+W],[partial-W,0]], W=p/sqrt(Q)>0.
Its square gives P_+/-=-partial^2+M^2+W^2 +/- W'. Let l_+/- be
the left regular logarithmic derivatives and r_+/- the right inward
ones. Exact intertwining, including their regular endpoint branches,
gives

    (l_+-W)(l_-+W)=M^2,
    (r_++W)(r_--W)=M^2.

Set a_L=l_-+W, b_R=r_--W and

    U=(a_L-M)/(a_L+M), V=(M-b_R)/(M+b_R).

Regular bolt behavior and scalar comparison imply 0<=U,V<=1.
For the lowest omega=1/2 the regular, nonlogarithmic branch is retained;
the endpoint is not replaced by a wall. In its propagation direction
each variable obeys the SAME stable Riccati equation

    dot U=W(1-U^2)-2M U.                              (1)

Indeed a_L'=M^2+2W a_L-a_L^2, whereas b_R evolves inward by
dot b_R=M^2-2W b_R-b_R^2. At the seam reflection gives V=U.
At an inner computational cut z=15/16, start the entire interval [0,1].
It contains every actual regular-bolt datum, even when an individual
square potential is negative. Propagate to the seam, then propagate V
back from the enclosed seam value. Neither bolt datum is guessed.

For constant W, K=sqrt(M^2+W^2), T=tanh(KL), the exact flow is

    Phi(W,M,U0,L)=[W T+(K-M T)U0]/[K+M T+W T U0].     (2)

Comparison makes it increasing in W and U0 and decreasing in M on
[0,1]. Evaluate the two extremal flows with the complete enclosed cell
length. A constant flow is monotone toward its equilibrium, so the
initial and final extremal values also enclose U throughout that cell.
This bounds integrals of M+WU without differencing large chiral terms.

Write I_AB=integral_A^B(M+WU)d rho. The diagonal and separated square
kernels obey

    G_-=(1-U)(1+V)/[2M(1+UV)],
    G_+=(1+U)(1-V)/[2M(1+UV)],
    G_-(A,B)=(1-U_A)(1+V_B) exp(-I_AB)/[2M(1+U_B V_B)].

The last identity follows by integrating
l_-=M+WU+partial log(1-U). Intertwining also gives
G_+(A,B)/G_-(A,B)=a_L(A)b_R(B)/M^2. Multiplying the full first-order
2 by 2 Green matrices before taking absolute values yields exactly

    B_direct=-4 U_A V_B exp(-2I_AB)/(1+U_B V_B)^2,
    B_image =-4 U_A U_B exp(-2I_image)/(1+U_B V_B)^2,
    I_image=integral_A^seam(M+WU)+integral_B^seam(M+WV). (3)

At mirror(B), left U is V_B and right V is U_B. That explains both
the image exponent and its numerator. In a constant-W full line the
direct formula reduces to -W^2/(M^2+W^2) exp(-2sqrt(M^2+W^2)L).
This cancellation is real: bounding the four matrix products separately
would lose the vanishing at W=0. It is not a change of normalization.

Define the dimensionless reduced sum

    R_AB=gamma sum_(n>=1) (d_n/q) M_n^2
                     sum_(j>=0) [B_direct(p_j)+B_image(p_j)]. (4)

Fourier normalization and the two signs of omega give the action
source Hessian kernel, per d rho_B at output A, as
q R_AB/(pi sqrt(Q_A)). The retained normalized disk row multiplies
this by -2G/B^2. Report both (4) and that conversion. A kernel value
per proper-distance input is not its integral against an actual
correction, nor a raw C5 estimate.

## 3. Absolute tails for every omitted mode

The exact polynomial gives R'(z)>=2beta-4|epsilon|>0, R(1)=1.
Hence Q<=(1-z)^2 on the physical interval. Also

    (1-z)^2(1+2z+3z^2)=1-4z^3+3z^4<=1,
    Q<=beta+u/16+|epsilon|<.2.

Thus W>=sqrt(5)p, L_AB>=log(4), and cut-to-A distance>=log(2)>.25.
For K*=sqrt(M^2+5p^2), the constant comparison starting at zero gives
U>=U_* tanh(K* .25), U_*=sqrt(5)p/(M+K*). Since M>=18,
tanh(18/4)>.99. After that buffer, M+WU>=.99K*. The same bound
holds along the longer reflected path. Consequently

    |B_direct+B_image|<=8 exp[-1.98 log(4) sqrt(M^2+5p^2)]. (5)

For alpha,beta_t>=0 with alpha^2+beta_t^2<=1, set
c_M=1.98 log(4) alpha, c_p=1.98 log(4) beta_t sqrt(5).
Then (5) is bounded by 8 exp(-c_M M-c_p p). The parameter box gives
q/C<400, M_n>=18sqrt(n) and

    (d_n/q) M_n^2<=800n+2400n^2/q+1600n^3/q^2.

Let S_N(c) be the sum of this polynomial times exp(-18c sqrt(n))
over n>=N. For each polynomial power r=1,2,3, monotonicity follows
from 18c sqrt(N)>2r. Bound the series by f(N)+integral_N^infinity f.
With a=18c and t=sqrt(N), the integral is exactly

    2 exp(-a t) sum_(k=0)^(2r+1)
        [(2r+1)!/(2r+1-k)!] t^(2r+1-k)/a^(k+1).    (6)

The thermal numbers remain discrete. For j>=J and p_min=gamma(J+1/2),

    gamma sum_(j>=J) exp(-c_p gamma(j+1/2))
       <=exp(-c_p p_min)(1/c_p+gamma),              (7)

using 1/(1-exp(-x))<=1+1/x. Use alpha=.99,beta_t=.14 for the
n>=3 tail and alpha=beta_t=.7 for the frequency tail at all n>=1:

    T_angular=8(1/c_p1+gamma) S_3(c_M1),
    T_frequency=8 exp(-c_p2 p_min)(1/c_p2+gamma) S_1(c_M2).

The overlap of those two tails is deliberately counted twice. Both
omitted contributions are nonpositive; append [-T_angular-T_frequency,0]
to the finite interval. Equations (5)-(7) prove absolute convergence
and justify the sum and off-diagonal variations. No common finite part
or ultraviolet restoration is inferred from this separated calculation.

## 4. Exact thermal bins and predeclared execution

Explicitly enclose n=1,2 and every j below J_max. Choose P_cut=16.
Allocate three quarters of the bins uniformly to nominal p in [0,4]
and the remaining quarter uniformly to [4,16]. If p_l is that exact
rational node, set J_l=floor[p_l/gamma_classical]. Each bin contains
the exact integer set J_l<=j<J_(l+1). Enclose the actual p interval
using gamma, and multiply its kernel interval by gamma times the exact
integer count. This is an enclosure of a discrete sum, not a
zero-temperature integral or an unbounded quadrature approximation.

Use 256-bit Arb, python-flint 0.9.0 in the existing isolated installation.
The mesh d=1-z=2^(-4+i/N_oct), i=0,...,4N_oct, contains A, B and
the seam exactly. Whole-cell Q and proper-length intervals feed (2).
Run schedules (N_oct,N_frequency)=(1024,256),(2048,512). Require finite
negative full intervals, mesh overlap and narrowing; final relative
width of -R_AB <.20 and total analytic tail <.01 of its lower magnitude.
These are precision targets for this separated matrix element only.
If they fail, retain the failure and refine; do not weaken them.

The first numerical execution used uniform nominal-frequency bins and
schedules (512,128),(1024,256). The final relative width .7579061 failed
the unchanged .20 target, although its tail ratio was <3.986e-7. Preserve
that run. The revised mesh above concentrates effort on the important
low frequencies and doubles the radial resolution. A preceding startup
import failure ran no scientific calculation; its terminal error and
source-reconstruction provenance are retained separately in the receipt.

Before interpretation check symbolic first-order matrix cancellation,
constant-coefficient transfer and semigroup identities, the exact
polynomial-tail integral, and overlap with the preceding independently
enclosed square traces in its three selected sectors. Freeze this note,
the executable and the review before every run. Preserve all numerical
bins, tails, dependency hashes, raw check rows, snapshots and failures.
The [Arb documentation](https://python-flint.readthedocs.io/en/latest/arb.html)
describes ball operations; continuum validity rests on the argument
above, not on the count of passing checks. No Lean is selected.

## 5. Research consequence and stopping boundary

A successful test closes the full mode summation and actual paired
attachment for this separated static shape kernel. It does not supply
the near-diagonal distribution, differentiated response, complete
source, projected inverse, compatibility or nonlinear contraction.
Smallness at two separated points does not bound the local operator:
the order-four logarithmic shape symbol remains.

Continue the specified sharper method: construct the commonly
subtracted near-diagonal response on this same finite geometry, then
use it in the actual differentiated H e=-r1 solve. The next decisive
calculation is the locally restored diagonal/near-diagonal action of
that response on the actual defect/correction, with a controlled
derivative remainder; a second separated pair is not a substitute.
The present run executes only the separated assembly test. It does
not execute that next calculation or certify attainability of the
full metric solution. Retain all quantitative admission and nonlinear
targets in [the closed audit](cold-defect-simultaneous-feasibility-audit.md).
Retarded boundary realization, matching, semiclassical validity and
all-sector stability remain full requirements. S169.2-S169.5 stay open.

Remaining to-do list: commonly subtracted near-diagonal response with
derivative control; actual finite-domain correction, nonlinear admission
and compatibility; physical matching, semiclassical validity, causal
boundary realization, all-sector stability and the full wormhole goal.
