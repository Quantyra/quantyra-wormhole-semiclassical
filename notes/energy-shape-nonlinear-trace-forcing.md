# Nonlinear trace forcing with the full first-order shape response

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The first-order shape in the retained full auxiliary energy space does
not change the persistent second-order trace forcing. Its additional
forcing is a sum of at most two distributional derivatives of functions
in L1+L2 in logarithmic radius. The resulting bounded second-order tail
has the same persistent trace/action combination as before. This resolves
the nonlinear trace obligation left by the
[mixed-boundary result](first-order-mixed-boundary-response.md).

It does not yet resolve the pure shape quadratic finite-cut action.
That coefficient, including the canonical SECOND variation, is now the
only unevaluated possible persistent contribution within the expansion
specified below. Its closed energy form is finite, but its literal cutoff
value has not been identified. No full fast-branch exclusion, physical
mode, nonlinear wormhole or stability theorem follows.

## 1. Scope and the useful form of the shape norm

Keep the same charged seed, fixed background cut, paired matter domain,
state, flux, physical clock, subtraction and finite coefficients. Fix
one positive Matsubara frequency above the retained auxiliary thresholds.
Write z=log(r/B), t=G lambda/(60pi)>0, k^2=1/t+17/12. At first order
the exact metric tangent is (p+s,v), where

    p=R_gamma U/r, s=Q_gamma w/r, v=w/r,
    U in D_max(Atilde), w in E_gamma.

The auxiliary shape equation and any continuous admitted shape source
produce such a w. Physical constraint and source compatibility remain
additional requirements. The field need not have pointwise fourth jets.
For a family satisfying the full equations its first-order equation is
also necessary; the calculation here does not supply that solution.

The retained estimates imply the following stronger statement in z:

    A(z)=exp(3z/2)[p+s-p_fast], V(z)=exp(3z/2)v
                         belong to H2((Z,infinity),dz),
    p_fast=exp(-3z/2)f(z), f=a cos(kz+phi).              (1)

Here a is nonzero for the branch being tested. To see the gain, first
consider either shape amplitude h=s or v. The weighted H1 bounds on
r h give H1 control of exp(3z/2)h. In the normal geometry B=lambda=1,
put y=exp(-z), h=y^(3/2)H(z). Then exactly

    integral (h/y)^2 dy/y^2 = integral H^2 dz,
    integral [(h/y)_y]^2dy = integral (H'-H/2)^2 dz,
    integral h_yy^2dy = integral (H''-2H'+3H/4)^2 dz.    (2)

The charged seed changes these identities by smooth bounded factors
with bounded inverse and coefficient corrections O(exp(-2z)). The J2
bound on v, and on s+v/3, therefore upgrades both profiles to H2.
The conformal reconstruction error has its first two differentiated
bounds from the preceding note, O(exp(-z)) after this rescaling, and
also belongs to H2. This proves (1). In particular A,V and their first
derivatives tend to zero; their second derivatives need only lie in L2.

We allow second-order logarithmic metric terms

    Phi-Phi_seed=eta exp(-3z/2)(f+A)H(tau)
                                  +eta^2 exp(-3z)h(z,tau),
    V_metric-V_seed=eta exp(-3z/2)V H(tau)
                                  +eta^2 exp(-3z)j(z,tau). (3)

H is the fixed nonzero thermal cosine. h,j have bounded radial jets
through the order used by the second-order equation, and fixed finite
thermal frequencies. Smaller remainders must have the corresponding
controlled distributional order: after normalization, their contribution
must be a sum of up to two derivatives of L1+L2 functions as in (8).
An arbitrary pointwise little-o is insufficient. No theorem that every
exact solution has (3) is asserted.
Normalize the temporal means by hbar=<h>/<H^2>, jbar=<j>/<H^2>.

## 2. The complete coupled quadratic trace in divergence form

Use the original normal local density L=L_C+L_W from the
[nonlinear trace derivation](fast-nonlinear-trace-boundary-compatibility.md).
With w=Phi-Phi_seed, v=V_metric-V_seed and D=d/dz, it is

    L_C={exp(2w)[6(1+w')^2+4(1+w')v']
                          +6(1+11t/6)exp(4w+2v)}/4,
    L_W=15t exp(-2v){(z+w)(v'+v'')^2/15
       -(1+w'+w'')(v'+v'')/15-(1+w'+w'')^2/10
       -(v'+v'')(1+w')^2/15
       -7(1+w'+w'')(1+w')^2/90-7(1+w')^4/180}.        (4)

The trace density is partial_w L-(D+3)partial_w' L
+(D+3)^2 partial_w'' L. The canonical determinant has no conformal
derivative at fixed independent shape; its omission from this trace
row is exact. It remains in the shape row and quadratic action.

For arbitrary profiles P,V define three quadratic polynomials F_i by
the eta^2 coefficient of partial_(D^i w)L, substituting

    D^j w=eta exp(-3z/2)(D-3/2)^j P,
    D^j v=eta exp(-3z/2)(D-3/2)^j V,

and removing the common exp(-3z) factor. Then the normalized complete
quadratic forcing is exactly

    N[P,V]=F_0(P,V)-D F_1(P,V)+D^2 F_2(P,V).          (5)

Each F_i is a constant-coefficient quadratic in P,V through second
derivatives. In particular no explicit z survives. The potentially
dangerous z(v'+v'')^2 term is independent of w, and its conformal
derivative comes only from the w(v'+v'')^2 term. For reference,

    F_2/t=-21P^2/8+7PP'/2-5PV/2-7(P')^2/6
                     -22P'V/3+6P''V-41V^2/6-4VV'+2VV''. (6)

F_0 contains t(V'')^2, so an L1 component really is needed. The checker
derives all three polynomials from (4) and independently verifies (5)
using the warped metric's scalar, traceless Ricci, Weyl and Box R terms.
It keeps the density factor exp(4w+2v), all cross terms and the vacuum
shift. The underlying anomaly and conformal-curvature identities are
standard [Mazur--Mottola, section IV](https://arxiv.org/html/hep-th/0106151v1);
the energy-space and tail-transfer claims here are additional arguments.

For P=f,V=0, f''=-k^2 f, (5) reduces to the previous exact forcing

    N_t=(471t/4+289/4+6/t)f^2
                         +(35t/3-5)ff'-(21t+3)(f')^2. (7)

For P=f+A with (1), define R_i=F_i(f+A,V)-F_i(f,0).
Terms linear in A,V and their first two derivatives belong to L2:
all coefficients made from f and its derivatives are bounded at this
fixed k. Quadratic terms belong to L1 by Cauchy--Schwarz. Consequently

    R_i in L1+L2, i=0,1,2,
    N[f+A,V]=N_t+R_0-D R_1+D^2 R_2.                 (8)

Their L1 norms on translated bounded intervals tend to zero. For a
fixed compact test psi the difference pairs as
integral (R_0 psi+R_1 psi'+R_2 psi''); this tends to zero on translation
toward infinity. No pointwise estimate on third or fourth jets is used.
Smooth approximation in H2 justifies this second-variation trace as a
distribution when those higher jets do not exist as functions.

The exact charged local density has the same limit in this divergence
form. Its coefficient differences and their fixed-order z derivatives
are O(exp(-2z)), from the explicit rational F and spherical curvature.
Multiplying bounded fast jets and H2 shape jets therefore adds only
L1+L2 terms to (8). Fixed-frequency time derivatives also carry the
exp(-2z) gain; a gamma^4 term has no radial derivatives. The conformal
trace is a local curvature polynomial and has no remaining explicit
logarithm. These estimates apply before differentiating R_i. They do
not incorrectly differentiate a merely L1 function pointwise. The
linear charged corrections on bounded second-order h,j have the same
decaying coefficient property. Thus the actual mean tail equation is

    L_h hbar+L_j jbar+N_t+R_0-D R_1+D^2 R_2=0,        (9)

with the aggregate R_i still in L1+L2. Both temporal source sectors
remain required in the full system; taking the mean extracts a necessary
condition and does not solve static signed compatibility.

## 3. Transfer the weak forcing through the exact tail equation

Put C(D)=3tD^2-9tD+11t+3 and

    X=3(D-4)hbar+(D-6)jbar,
    P(D)=(D+1)C(D).

The unchanged common-row identity gives

    P(D)X=3N_t+3(R_0-D R_1+D^2 R_2).                 (10)

An explicit two-sided Green function avoids an unproved pointwise
forcing estimate. With c=1/(3+23t), define

    G(x)=c exp(-x),                                    x>=0,
    G(x)=c exp(3x/2)[cos(kx)-5 sin(kx)/(2k)],           x<0. (11)

It and its first derivative are continuous at zero; its second
derivative has jump 1/(3t). Direct differentiation gives P(D)G=delta.
G,G',G'' are bounded and in L1 and L2. Write

    G'''=g_3+delta/(3t),

where g_3 is its ordinary piecewise third derivative and has the same
integrability. Extend R_i by zero below a fixed Z. A particular error
response on the whole line is

    Y=3(G*R_0-G'*R_1+G''*R_2),
    Y'=3(G'*R_0-G''*R_1+g_3*R_2)+R_2/t.              (12)

The contact R_2/t is essential. Approximate each L1+L2 input by compact
inputs in those norms. Boundedness and L2 integrability of the kernels,
followed by their exponential decay, prove Y(z)->0. The convolutions
and R_2 in the second line have vanishing translated local L1 norms.
Hence Y' tends to zero in that topology. Constants are finite for every
fixed t>0; no uniform t->0 or numerical benchmark threshold is claimed.

Let X_per be the known periodic solution of P X_per=3N_t. Every bounded
solution of (10) on the half-line satisfies

    X=X_per+Y+c_- exp(-z).                            (13)

The other homogeneous roots are 3/2 +/- i k and are excluded by
boundedness, since Y is bounded and tends to zero. Thus X-X_per tends
to zero and X'-X_per' tends to zero in translated local L1. Under the
bounded third-jet hypothesis on h,j in (3), X'' is bounded. Uniform
continuity then upgrades the latter convergence to pointwise zero:
a derivative of size at least epsilon with a fixed Lipschitz bound
would give a fixed positive L1 mass on a nearby interval. This proves
the stronger conclusion without requiring a pointwise forcing limit.

## 4. Consequence for the remaining boundary problem

The earlier unchanged second-order boundary contribution is

    M_h hbar+M_j jbar=(t/7)(7D-5)X.

Together with the pure conformal coefficient Q_t, (13) leaves the same
persistent normalized vector for unit fast amplitude:

    R_c=-(78815t^3+74391t^2+16137t+81)
                         /[168(5t+3)(11t+3)] < -3/56,
    R_s=-k(6655t^3+111t^2-3015t-567)
                              /[84(5t+3)(11t+3)].      (14)

The vector rotates with 2phi and scales with a^2. The local second-order
tails therefore cannot independently repair the persistent term when
first-order shape is included under (1),(3). Their canonical FIRST
variation still has a finite limit by the
[cubic-direction result](coupled-fast-quadratic-boundary-obstruction.md).
The first-order mixed action has a finite translated-logarithmic limit
by the preceding note. Those two finite values are retained.

Accordingly, if the total original quadratic cutoff coefficient on (3)
has an ordinary limit, the PURE SHAPE quadratic coefficient must carry
the opposite of (14), in the necessary translated-logarithmic sense,
with physical prefactor S0<H^2>a^2, S0=lambda B^3/G. This is now an
explicit requirement on the unevaluated canonical second variation and
its local shape partners. If that pure shape coefficient instead has a
finite translated-logarithmic limit, the full family is obstructed.
That premise has not been proved and is not silently adopted.

The existing canonical comparisons are on the limiting geodesic seam
and the compact core or its form closure. A cut at positive y has
different boundary jets and moves the matter boundary. Those estimates
do not automatically control its second variation uniformly in the
cutoff. No logarithmic cutoff average is prescribed as a physical action,
no boundary counterterm is changed, and no complete branch exclusion is
claimed. Establish that finite-cut quantum comparison next.

Verification includes all coupled quadratic trace terms from both the
action and independent curvature invariants, charged unexpanded trace
fixtures, concentrated H2 profiles with large pointwise forcing, and the
two-sided Green function with its distributional contact. These are
supporting calculations, not a nonlocal shape solve or physical mode.
No Lean run is selected. The goal remains incomplete.

Remaining to-do list: pure shape quadratic finite-cut limit, including
canonical second variation; physical boundary phase, source compatibility,
continuous causal constraints and all-sector stability; exact nonlinear
existence; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
