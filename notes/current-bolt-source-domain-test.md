# Actual coupled source domain: the global radial-graph shortcut fails

2026-10-26 HST, E003/S169. Ordinary argument BEFORE arithmetic. Previous
turn: PROGRESS, closing the complete transport-to-bolt error slice.
Keep the full stable-wormhole objective, closed audit, accepted interior
result and P32. The complete coupled estimate remains unresolved.

## 1. The unsupported step encountered while preparing a trial

The current backward-error proposal uses GLOBAL H3(d log z) norms for
the scalar row and differentiated shape row. It was explicitly a
stronger sufficient test, not an identification with original Y. Before
using it to recover original source comparisons, check its bolt domain.
The outer conormal derivative D=z partial_z is singular as a regular-disk
vector field at the bolt. Original regular-disk source regularity does
not automatically supply the proposed global graph norm.

This note tests the possible estimate

 Rfast_global(r) <= K [outer phase/graph observations
                         + regular-disk H3 norm of BOTH original rows]. (1)

It proves that NO FINITE K exists, on actual coupled sources r=H v.
It does not declare the physical H undefined, prove I_Y itself infinite,
or invalidate the accepted sufficient global-residual estimate. The full
original Y completion and its quantitative comparison remain to be
specified/proved. Equation (1) is one concrete proposed route to that
comparison, and fails even on the original smooth core.

Keep the original finite Q, paired mass/log domain and complementary
adjoint, common regulator and ultraviolet logarithm, complete angular
and thermal spectra, all local/coincidence contacts and finite returns.
The argument below uses a subspace of the actual domain, not a smooth
replacement of that domain or a selection of quantum modes.

## 2. A family in the actual complete coupled domain

Let L=z_h, Q(L)=0, Q'(L)=-2 gamma, and set

 w=sqrt(2(L-z)/gamma), z=L-gamma w^2/2, theta=gamma tau.

This is a regular radial bolt chart: the orbit metric is
 [gamma^2 w^2/Q] dw^2+[Q/gamma^2] dtheta^2,
with coefficient ratios analytic in w^2 and equal to one at w=0.
Its area is EXACTLY w dw dtheta; after division by the unchanged thermal
period, radial measure is dz=gamma w dw. No period or norm is rescaled.

Choose fixed smooth chi(w), equal to one for w<=1/4 and zero for w>=1/2,
and smooth psi(s), zero for s<=1 and one for s>=2. For epsilon>0 put

 v_epsilon=(nu_epsilon,T_epsilon),
 nu_epsilon=a chi(w) psi(w/epsilon) w^(17/2), T_epsilon=0. (2)

Here a is any fixed nonzero real amplitude and 0<epsilon<1/16. Every
member is smooth and zero near the actual bolt and asymptotic boundary.
It obeys the original fixed-temperature/paired conditions, has zero
mass/log/slow coefficients, and no new wall. All original .05/.1 and
.1/.2 outer transitions are retained (the fields vanish there).

The H7 norm in regular disk charts is uniformly bounded: in the inner
cutoff ring, a kth derivative is O(epsilon^(17/2-k)), and area is
O(epsilon^2). For k<=7 its L2 contribution tends to zero. The same
family is bounded in C5. It converges in H7 to a chi(w) w^(17/2).

The canonical determinant has NO nu input by the already derived exact
Ward identity. Thus H v_epsilon is the COMPLETE local scalar/mixed column,
including every Weyl, density, R2 and Maxwell term. The quantum modes
have not been dropped: their nonlocal derivative in this direction is
identically zero. Both original output equations are retained. The same
identity holds for the commonly restored retained transport operator.

## 3. Compute the full leading rows, not a principal symbol guess

Let D_R=96 pi c_R-G/(60 pi), d=D_R C/B_radius^2>0, c_R=10^22, G=1.
The complete fourth-derivative scalar/mixed column at a regular bolt is

 (F_nu,F_T)=(6d,2d) Delta_orbit^2 nu + lower radial orders. (3)

For clarity, the R2 quadratic density gives coefficients
576 pi c_R C/B_radius^2 and 192 pi c_R C/B_radius^2. The full Weyl
quadratic density adds -GC/(10 pi B_radius^2) and
-GC/(30 pi B_radius^2), respectively. All other terms have fewer
derivatives with smooth disk coefficients. The exact full local generator
will independently check (3), including the variable-Q third/second
derivative terms contributing to the SAME indicial order.

For r=z^4 F, alpha=17/2 and w outside the inner cutoff ring, (3) gives

 r_nu=a A_nu w^(9/2)+O(w^(13/2)),
 r_T =a A_T  w^(9/2)+O(w^(13/2)),
 A_nu=6d L^4 alpha^2(alpha-2)^2, A_T=A_nu/3.             (4)

The remainders have the differentiated expansions needed below, since
all fixed coefficients are analytic in w^2. This is a one-sided Taylor
calculation on the actual regular bolt, not a new physical geometry.
Both rows are uniformly bounded in regular-disk H3: their third spatial
derivatives are O(w^(3/2)), and cutoff-ring H3 contributions vanish.
The original outer graph and finite-part observations are zero.

However D=z partial_z=-(z/(gamma w)) partial_w. With sigma=4.005,

 D^3(z^-sigma r_nu)
 =-a (45/8) A_nu L^(3-sigma) gamma^-3 w^(-3/2)
                                               +O(w^(1/2)). (5)

Using |d log z|=gamma w dw/(L-gamma w^2/2),
define

 Kdiv=(45/8) A_nu L^(5/2-sigma) gamma^(-5/2)>0.

The actual nonnegative H3 norm includes the square of (5). Integrating
over 2epsilon<w<r0, then taking epsilon to zero, proves

 liminf epsilon^(1/2) Rfast_global(H v_epsilon)
                                      >= |a| Kdiv/sqrt(2)>0. (6)

Lower-order terms cannot cancel this leading power; the second source
row cannot cancel a norm of the first row. Hence (1) is unbounded even
on smooth actual-domain fields. A fixed amplitude can be made arbitrarily
small in physical C5 without changing that conclusion. This is a failure
of a SPECIFIED source-norm estimator, not of the construction.

Execute the full indicial coefficients and Kdiv. The coefficient in (6)
is an asymptotic growth coefficient, NOT a lower bound at a particular
finite cutoff without a remainder bound. No finite improvement factor
or normalization change repairs (1).

This explains the derivative mismatch: three z derivatives at the bolt
can cost six regular Cartesian derivatives; the differentiated shape
row can demand eight. The auxiliary global test therefore imposes more
than regular H3 source data. Its use on an extra-regular explicit trial
is still legitimate if that trial's FULL residual is proved to satisfy it.

## 4. A complete residual test that respects the original bolt chart

Retain the original nonnegative beta(11,11) cutoff chi_o, one through
z=.05 and zero from z=.1. For the ACTUAL full two-row residual r, after
exact prescribed outer phase cancellation, split
 r_o=chi_o r, r_i=(1-chi_o)r.
Keep all cutoff derivatives in

 f_o=z^-sigma r_o,nu,
 g_o=z^-sigma[(D-4)r_o,T-(D-6)r_o,nu/3],
 R_o=(||f_o||H3^2+||g_o||H3^2)^(1/2),
 c4=lim_(z->0) z^-4(r_T-r_nu/3),
 R_i=||(r_i,nu,r_i,T)||L2(dz).

R_i covers the ENTIRE complementary interval, including the .1/.2
transition, matching quotient and regular bolt. It is not replaced by
finite boundary observations. Its stationary L2 measure is the actual
area per thermal period, as above. The split introduces no boundary.

The accepted exact undifferentiated Volterra recovery gives
 ||r_o||L2(z^-4 dz)<=Krow R_o+|c4|/sqrt(5).
Since r_i vanishes for z<.05,
 ||r_i||L2(z^-4 dz)<=400 R_i.
Use the unchanged physical source pairing and return constant MR/beta:

 ||r||E* <= (MR/beta)[Krow R_o+|c4|/sqrt(5)+400 R_i],
 |L_b(H^-1 f-u_trial)|
     <= Knew (MR/beta)[Krow R_o+|c4|/sqrt(5)+400 R_i].     (7)

This bounds a COMPLETE source-to-bolt map. It uses both original rows
and the actual finite part, and requires no artificial z derivatives
at the bolt. Evaluate its coefficients and the separate zero-other-debit
ceilings for half the bolt allowance. They are sufficient residual
targets, NOT actual residuals, original-Y coefficients or output budgets
already attained. The same physical target and normalization are retained.

For the actual complete transport error, the previously accepted direct
E* estimate remains stronger; do not optimize its already passing tail.
Equation (7) is useful for evaluating a retained trial's residual. It
does not yet evaluate that trial or its common-restored quantum center.

## 5. Consequence for the full feasibility test

Outcome: verified failure of the proposed global radial-graph comparison,
with a finite complete mixed-chart residual alternative. This is not
an undefined physical operator or an obstruction to wormhole existence.
The earlier bolt-map bounds, degree-63 tail test and transport-error
slice stay valid. No graph/C5 output, original-source defect, actual
kappa_m, signed coupled center or nonlinear requirement closes.

In particular the original requirements Cdef3 IY<=1.0965284986e32,
the graph/physical ceilings only with their stated remaining margins,
actual kappa_m<=5e-5 and EVERY other coupled entry stay unresolved.
The two-fast quantum screen retains its limited scope. The scalar trace
remains the dominant known bolt-estimator loss, but optimizing it cannot
replace the missing original source/phase/graph maps or signed response.

Decision: retain the construction and P32; replace the attempted use of
GLOBAL H3(d log z) as an original source comparison with the original
outer graph plus regular-chart interior residual. This does not redefine
the physical admission norm. Keep the extra-regular global test only
as an optional sufficient test, with its hypotheses explicit.

One next calculation, PROPOSED and NOT EXECUTED here: construct and apply
the full commonly restored quantum operator to a phase-conforming coupled
trial, enclosing its R_o, c4 and R_i together. The signed retained center
has not been calculated. All angular/thermal sums, contacts and returns
must remain; the point of (7) is to make that complete residual test
well-defined without an unsupported bolt regularity conversion.

Verification contract: freeze ordinary proof/checker, execute the exact
full-column indicial identities and 512-bit original-unit coefficients,
retain outcomes, ordinary review, unchanged acceptance. Checks support
the continuum counterexample and integral estimates, not a quantum solve.

Remaining to-do list: full retained quantum application and coupled
residual; original source/graph/physical maps and every coupled entry;
nonlinear admission/contraction; exact compatibility; physical matching;
semiclassical validity/noise; all-sector stability; opening, transport,
preparation, separate creation and the full stable-wormhole objective.
