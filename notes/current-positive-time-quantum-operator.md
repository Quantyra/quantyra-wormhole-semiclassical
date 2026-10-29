# Complete positive-time quantum operator and first-solve debit

2026-10-28 HST, E003/S169. The preceding turn supplied a conditional
application contract, but no signed quantum bound. Preserve 5967b97,
the closed audit, P32, the original physical norms and the full stable
semiclassical wormhole objective. This tests a part of the ACTUAL quantum
operator uniformly on the whole energy space, with its effect on the
complete first-solve error. It does not evaluate the remaining short-time
action or certify the full signed response.

## Actual realization and the first extension to justify

Use e_t=ell_t/N on E_t and the original York shape T=Vxi. The accepted
comparison is ||Vxi||L2(dz)<=CLt ||xi||E_t. Its proof includes the original
mass/log domain, cutoff transitions and regular bolt; it is not an
independent boundary-phase or nonlinear-neighborhood estimate. A=iJD0 is
the self-adjoint operator on the accepted PAIRED first-order realization;
L=A^2 uses its composed domain. The complementary adjoint, physical-half
rank-two trace, temperature, spin and subtraction scale are unchanged.

The complete positive-time heat Hessian on the original admitted core is

 D2 K(t)[u,v]=sum_(j,k) q_t(a_j,a_k) u_jk v_kj,
 q_t(a,b)=-t(a+b)^2 [g_t(a^2)-g_t(b^2)]/(a^2-b^2),
 g_t(lambda)=lambda exp(-t lambda).

Use the coincident limit and zero when a+b=0. This is the differentiated
paired conformal identity from current-structured-weak-response-attempt.md,
including the second/coincidence vertex. It is not an unregulated angular
finite part or a mass-only bubble. That source proves, for all real a,b,

 |q_t(a,b)|<=32[exp(-t a^2/2)+exp(-t b^2/2)].          (1)

The sharper Cauchy--Schwarz step in current-quartic-column-energy-dual.md
does not depend on its quartic input. For real scalar u,v, each weighted
double sum in (1) is bounded by

 sqrt(Tr[u^2 exp(-tL/2)]) sqrt(Tr[v^2 exp(-tL/2)]).

Indeed sum_k |u_jk|^2=<e_j,u^2 e_j>; the same identity holds for v.
Both spectral weights therefore give the factor 64. Divide by the SAME
time period 2pi/gamma only after this inequality. The unvaried actual heat
diagonal bound M<734 at t*=.0005, from the accepted disk propagation,
is independent of the direction. Its positive spectral representation
implies the same diagonal bound for every later time. Hence, for t>=t0=.001,

 |D2 k(t)[Vxi,Veta]|
     <=64 M CLt^2 ||xi||E_t ||eta||E_t.               (2)

All thermal, radial and boundary-return modes are in this trace. The
unvaried bound is an actual-domain kernel bound, not an operator-norm
assertion about its diagonal. Absolute summability proves (2) first on
the original core; L2 convergence and the uniform bound extend it uniquely
to E_t x E_t. This addresses the first domain/regularity implication needed
here. It does not assert that every energy vector is a C5 metric path.

## Split the common regulator without dropping contacts

Let c=G/(B^2 beta), G=1 in the retained units, and

 Theta(t)=sum_(n>=1) d_n exp(-m_n t),
 m_n=n(n+q)/C, d_n=2(q+2n), q=2000^4.

The massive action uses ONE common Mellin regulator on Theta K. Its
t>=t0 part is holomorphic at regulator zero and has the unambiguous
energy-form Hessian

 Qplus[xi,eta]=c int_(t0)^infinity Theta(t)
                             D2 k(t)[Vxi,Veta] dt/t. (3)

The sign is that of the positive TT energy convention; all bounds below
are absolute and need no positivity of this Hessian. There is no angular,
thermal or late-time truncation in (3). The canonical LLL allocation is
separate and retained exactly.

It is essential that (3) is the RAW positive-time heat Hessian. The
previous commonly subtracted remainder has

 D2 r(t)[u,v]=D2 k(t)[u,v]-C0[u,v]/t-t C2[u,v],
 C0[u,v]=(2/pi) int u v dz,
 C2[u,v]=-(1/(480pi)) H4[u,v],
 H4[u,v]=int[4R^2 u v+8R(u Delta v+v Delta u)
                                      +8 Delta u Delta v]dz,
 R=-Q'', Delta=(Q partial_z)'.

The a1 variation vanishes on the fixed physical phase. The completed
a2 includes its paired boundary representative. Retain those identities,
the exact UV logarithm and every local allocation. Splitting the SAME
remainder at t0 gives the exact change to the known local energy form

 Delta_local=a0 int u v dz+a2 H4[u,v],
 a0=-2c I0/pi, a2=c I2/(480pi),
 I0=int_(t0)^infinity Theta(t)dt/t^2,
 I2=int_(t0)^infinity Theta(t)dt.                      (4)

Thus full massive energy = old local energy + Delta_local + the original
subtracted integral on (0,t0) + Qplus. Neither sign in (4) may be omitted.
This is an exact allocation within the original common regulator, not a
change of coupling, normalization, scale or physical boundary condition.
The shared physical vacuum coefficient retains its accepted enclosure;
it is not independently reset or estimated by this calculation.

## Uniform bound and complete finite-contact enclosures

The accepted all-angular estimate and actual gap m1>324 give
Theta(t)<=4C/t0 exp[-324(t-t0)] for t>=t0. Combining with (2),

 ||Qplus||_(E_t to E_t*) <=kplus,
 kplus=256 G C M CLt^2/(B^2 beta *324*t0^2).           (5)

This deliberately retains the existing simple integral bound 1/(324 t0).
No positive-time threshold or scalar component is optimized.

To make (4) implementable, enclose both coefficients with complete angular
tails and outward arithmetic. For n<=256 use

 I0_n=d_n E_2(m_n t0)/t0,
 I2_n=d_n exp(-m_n t0)/m_n,
 E_2(x)=int_1^infinity exp(-x s)ds/s^2.

These are exact positive scalar integrals, evaluated with Arb. For N=256,
n0=N+1, a=q/C and zeta=exp(-a t0), the ENTIRE omitted angular heat sum is
bounded by

 S_tail=zeta^n0[(2q+4n0)/(1-zeta)+4zeta/(1-zeta)^2],
 I2_tail<=S_tail/(a n0), I0_tail<=S_tail/(a n0 t0^2). (6)

No finite matter selection replaces the operator. All thermal modes and
finite returns remain in the actual heat Hessian; (6) sums only scalar
weights for the exact finite contacts in (4). Retain the parameter balls
and their correlations in the formulas. Implement the two resulting
coefficient enclosures by explicit rational decimal centers and charge
their full radii, including parameter, angular-tail and arithmetic errors.

The original energy comparison also gives ||Delta Vxi||2<=CDt ||xi||E_t,
with CDt obtained from the retained fourth reserve, including its original
radius/volume factors. Since |R|<=2 on the physical half,

 H4norm <=16CLt^2+32CLt CDt+8CDt^2.

If da0,da2 are the two absolute coefficient errors, the COMPLETE local
reallocation error is

 kcontact=da0 CLt^2+da2 H4norm.                       (7)

The H4 polynomial is retained as a single exact bilinear expression; a
later numerical evaluation of its functions still has to enclose its own
errors. Equation (7) certifies coefficient rounding, not arbitrary matrix
or spatial quadrature errors. The size of Delta_local is also reported,
so a small raw tail cannot conceal an omitted fourth contact.

## Original first-solve consequence and remaining budget

Let Ghat be the proved source upper from current-first-quantum-solve-attempt.md.
Subtract Qplus from ell_t and use the certified centers in (4), retaining
the SAME P32 transport difference. Let k=kplus+kcontact. The new seed
form has lower bound 1-k in the original E_t norm. Bounded form
perturbation therefore defines its exact first solution x_short and gives

 ||x_full-x_short||E_t<=k Ghat/(1-k),
 E_first,positive<=K_first k Ghat/(1-k),
 K_first approximately 3.2008398094e-44.              (8)

This includes the inverse effect, not just a residual evaluated at an
unknown trial. The scalar/mixed source is unchanged, since the canonical
massive determinant has only a shape row. The exact original inverse
ell_t^-1 is being approximated, not replaced by a later (ell_t+B)^-1 filter.
The fixed physical radius cost and two subsequent exact F_t filters are
already in K_first. E32 and the accepted source error are each charged
once; no previously certified debit is removed by renaming an allocation.

Execute (4)-(8) with 512-bit intervals and frozen inputs. Add the resulting
physical debit to the known partial subtotal; this is NOT a free remaining
margin. For a future numerical x_hat with residual rho_short in the new
form, the first-solve error is at most

 [k Ghat+||rho_short||E_t*]/(1-k).                     (9)

Report the corresponding CONDITIONAL rho_short ceiling for the existing
1e33 first-residual screen. No such short-time solve or residual has yet
been executed. This bound applies to the first inverse on this source;
it does not silently remove the positive-time operator from every later
filter or the independent correlated j_r term without their own bounds.

## Attempted complete column and decision boundary

The actual positive-time part is uniformly controlled by (5)-(8), including
all spectra, contacts and finite returns allocated above. The remaining
short-time part has its exact logarithmic principal operator. Its required
complete action on a localized Galerkin column is still UNKNOWN.

The current disk transport coefficients 953 and the boundary polynomial
149 were derived for the particular metric direction -z^4. Their unvaried
kernel table can be reused in (2); their differentiated short-time remainder
cannot be applied to a localized trial by assigning its constant to that
direction. The qualitative admitted-family C5 theorem proves neither the
needed E_t-dual constant nor the complementary Galerkin residual. The first
unsupported step for that alternative is now precise: a quantitative bound
on the full commonly subtracted short-time density functional for the
chosen localized column, uniformly over ALL original energy tests.

The full geodesic-reference method retains known form coercivity but still
requires its complete numerical inverse. A conforming Galerkin method would
avoid that primitive; it has not yet supplied the localized short-time
columns or complementary residual. The present result closes neither
alternative's full residual. Do not turn a finite matrix or this affordable
positive-time component into a signed G_t certificate.

Decision after the test: if (8) is affordable, use this certified enclosure
for the first solve and stop positive-time quadrature/optimization there.
Focus the next calculation on the complete short-time quantum column plus
its local terms and original-energy residual. Preserve the construction;
no existence obstruction or full coupled feasibility pass follows.

Remaining to-do list: complete signed first quantum residual and later
applications; Z_r/G_t and every coupled/source/graph/physical/boundary/mass
requirement; nonlinear admission/contraction and existence; exact
compatibility; physical matching/opening; semiclassical validity/noise;
all-sector stability; the full stable semiclassical wormhole objective.
