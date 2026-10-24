# A cutoff obstruction in the strong graph, and a matched principal subtraction

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The stable-wormhole objective remains ACTIVE and INCOMPLETE; audit CLOSED.

## 1. The proposed norm limit needs a correction

The previous checkpoint selected the original geometric Duhamel Hessian,
including its second vertex. That choice avoids the unproved variable
mass-chart conversion. It does not justify taking its common proper-time
cutoff to zero in operator norm after subtracting an UNREGULATED normal
operator. This note tests that precise operation before estimating its
curved remainder. It finds an actual obstruction to the operation, and
gives a subtraction that removes its identified cause.

Keep the original paired realization, canonical thermal period, massive
angular tower and mass/log phase. There is no new physical cutoff, boundary
wall, forcing or normalization. The cutoff epsilon below is only a common
computational regulator that must be removed. The weak C2 action on bounded
admitted families remains available. Uniform nonlinear X-ball admission
is not inferred. No full signed projection or controlled solution is
computed in this package.

## 2. Actual-domain wave packets disprove the raw norm limit

Use x=log z, sigma=4.005 and the full normal B_sigma from
current-outer-graph-feedback.md. Its fixed coefficient and asymptotic are

 K=GC beta^2/(15pi B^2),
 B(sigma+iN)=-K(iN)^5 log N+O(N^5), N -> infinity.        (1)

The actual all-angular shape row has local logarithmic coefficient

 -K a(z)(i xi_x)^5 log|xi_x|,  a(z)=(Q(z)/beta)^2.       (2)

The scalar/mixed rows have no such logarithmic principal part. These are
the same previously established physical symbols; no new flat-domain
identification is being used to assert (2). Local fourth contacts and the
explicit Weyl cocycle remain in the lower logarithmic order. The full
normal graph, not the insufficient H6 estimate, is used throughout.

Choose a smooth real bump psi supported inside
I=(.0019,.0021) in z, and put w_N(x)=psi(x) exp(iNx). Use the physical
direction b_N=z^sigma w_N, nu_N=0. Real and imaginary parts give real
directions. These are admissible identical-mouth variations: they vanish
in open neighborhoods of BOTH the original seam and the bolt. Their
boundary jets are zero. This tests the original domain, not a substitute
domain or a finite collection of physical modes. The existing output
cutoff equals one throughout I, since .0021<.01/e.

In the triangular input coordinates theta=nu+a_op(D)b,
a_op(s)=(s+3)/(3(s+1)), the corresponding comparison data are
(Jc_sigma a_op(sigma+partial_x)w_N, B_sigma w_N).
The factor s+1 cancels in Jc a_op, so its scalar component is order four.
Its H3 norm is O(N^7), whereas the shape component has norm
K N^8 log N ||psi||_2(1+o(1)). Thus the shape component determines the
normalized high-frequency test without prescribing an artificial source.

Let H_epsilon denote the geometric heat Hessian integrated over
t>=epsilon, with the common local subtractions/restorations retained.
For each FIXED epsilon>0 its localized quadratic form is at most order
four, with a constant allowed to depend on epsilon. To see the relevant
bound directly, split its Duhamel alpha integral at 1/2 and use the
longer heat factor. Cyclicity and Hilbert--Schmidt Cauchy--Schwarz bound
the two-vertex trace by products of

 ||L_u exp[-t(1-alpha)L/2]||_HS,

with the longer heat time at least epsilon/4. On the compact interior,
L_u has at most two derivatives of u; smooth positive-time heat kernel
jets therefore bound this by C_epsilon ||u||_H2. The complete angular
sum is finite at t>=epsilon, and its mass gap controls large t. The
second vertex obeys the same finite-order bound. The local a2 subtraction
adds an order-four form whose coefficient can grow as log epsilon.
It does not add a frequency-dependent logarithm. All thermal modes and
boundary contributions remain in the heat kernels in this argument.

Apply the original differentiated shape row and its H3 testing to these
packets. The finite-cutoff contribution is O_epsilon(N^8); the true
contribution from (2) is order K N^8 log N. Smooth coefficient derivatives
and the fixed weight z^sigma affect only the lower asymptotic terms.
Equivalently one can test the quadratic form with the adjoint row applied
to the packet: integration by parts places its one row derivative and
three Sobolev derivatives on that packet, giving the same powers.
The long-heat-factor bound just proved suffices for this packet matrix
element; a global H3 bound for H_epsilon is not assumed.

Let R_epsilon be the proposed feedback after subtraction of the full
UNREGULATED normal and already allocated logarithmic principal difference.
That removed logarithmic coefficient is exactly K a(z). Other already
allocated pieces have no N^5 log N term. Consequently the packet test gives

 ||R_epsilon||_(coupled H3 -> coupled H3)
       >= inf_(z in I) a(z),  for every epsilon>0,           (3)

where an unbounded map also satisfies the assertion. The same lower
bound holds for the norm distance between the regulated response and
its true response wherever that latter strong map exists. The true
local remainder has no logarithmic principal term after this subtraction.
For each fixed smooth packet the weak regulated forms converge by the
admitted-family argument. If R_epsilon were norm Cauchy, its norm limit
would have to agree with that weak limit. The nonzero packet lower bound
for their difference contradicts this. Thus the proposed norm-Cauchy
assertion is FALSE even after localization strictly inside the actual
collar. Fixed unregulated subtractions cannot change a Cauchy difference.

Execute the lower bound using the actual polynomial

 Q/beta=1+(u/beta)z^2+(-4-2u/beta-tau/beta)z^3
             +(3+u/beta+eps/beta)z^4-2eps z^5/beta+eps z^6/beta. (4)

Here eps in (4) is the retained geometric coefficient, not the heat
cutoff epsilon. Enclose (4) on the entire interval I, with all correlated
finite-geometry parameters. Require the conservative guard a> .9999999.
This is already incompatible with a .18 bound for the RAW regulated
feedback. It says nothing adverse about the desired unregulated matched
remainder, whose norm remains UNKNOWN. It is a method obstruction, not
wormhole nonexistence, instability or a change in semiclassical theory.

## 3. Identify the missing heat cap exactly

The previously derived density-inclusive geometric tangent multiplier is

 F(v)=1/2+(1/2)integral_0^1 (1+av)exp(-av)dalpha,
 a=alpha(1-alpha),
 R(v)=F(v)-1+v^2/120,
 0<=R(v)<=min(v^3/840,v^2/120).                       (5)

It includes the second vertex. In particular it is NOT a bubble-only
mass calculation. Its role here is to construct a principal subtraction,
not to replace the actual paired geometry by a smooth tangent domain.
The universal angular term Theta(t)=2C/t+lower terms gives the omitted
canonical heat cap

 U_epsilon(k)=-4GC k^4/(pi B^2) I(epsilon k^2),
 I(X)=integral_0^X R(v) dv/v^3.                      (6)

On freezing the ACTUAL coefficient at z_*, with k=sqrt(Q_*) xi_x/z_*,
multiplication by z_*^4 and the principal row derivative gives

 U_epsilon,*=-60K a_* (i xi_x)^5
                     I(epsilon Q_* xi_x^2/z_*^2).  (7)

This determines the normalization of the missing log before absolute
values. For a convenient exact expression set
psi(v)=(1+v)exp(-v)-1+v^2/2 and
Ein(v)=integral_0^v (1-exp(-r))dr/r. Integration by parts gives

 I(X)=1/4 integral_0^1 a^2[Ein(aX)-psi(aX)/(aX)^2]dalpha. (8)

The endpoint at aX=0 is a removable limit. The definition and
Ein(v)=gamma_E+log(v)+E1(v) are those of
[DLMF 6.2.3--6.2.4](https://dlmf.nist.gov/6.2#E3).
This source supplies the special-function identity only. It does not
supply a rough-domain or operator-norm theorem.

For all X>=1, (5) and 1/2<=F<=1 imply

 log(X)/120-1/4 <= I(X) <= log(X)/120+1/840.          (9)

For X<=1, 0<=I(X)<=X/840. Thus (7) has exactly the missing
-K a_* (i xi_x)^5 log|xi_x| at every positive cutoff. The cap tends to
zero at each fixed frequency, but not in the normal graph operator norm.
This explains why weak or pointwise convergence cannot prove the proposed
uniform strong limit. No rescaling of the graph fixes this missing term.

## 4. A quantitative full-angular test of the proposed repair

The exact frozen-coefficient cap retains Theta(t), not just 2C/t:

 U_full,epsilon(k)=-2G/(pi B^2)
                       integral_0^epsilon Theta(t)R(t k^2)dt/t^2. (10)

An elementary full-sum discrepancy bound is enough to test whether (6)
captures its problematic order. Let f_t(x)=2(q+2x)exp[-x(x+q)t/C]. Then
integral_0^infinity f_t=2C/t. Comparison on each interval [n-1,n] yields

 |Theta(t)-2C/t|<=integral_0^infinity |f_t'|dx
                    <=2q+4 sqrt(2C/t).             (11)

For the last bound f_t is unimodal. Its maximum is 2q when its stationary
point lies left of zero; otherwise it is at most 2sqrt(2C/t). Its total
variation is 2max f_t-f_t(0). This proves (11) for EVERY t>0 and retains
every angular level, without a spectral cutoff or a guessed Euler tail.

Using R(tk^2)<=t^2 k^4/120 in (10)-(6) gives

 |U_full,epsilon-U_epsilon|
 <=G k^4/(60pi B^2)[2q epsilon+8sqrt(2C epsilon)].    (12)

For the FROZEN principal row at z_*=.002, the previously certified full
normal graph constant Cq=2(17/16)^(5/2) gives

 eta_cap,model <= a_* Cq/log(2)
              [q epsilon/(2C)+2sqrt(2epsilon/C)].   (13)

This follows from K|xi|^5/|B(sigma+i xi)|<=Cq/log(2) at all frequencies.
It is a genuine uniform multiplier bound for this specified comparison,
not a finite-frequency sample. It establishes norm convergence of the
full-angular frozen cap after its universal part is restored. It does
NOT bound derivatives of the variable coefficient, the curved remainder,
thermal images, or seam/bolt transfer in the actual problem.

Execute (13) at epsilon=1e-3, 1e-4 and 1e-6. Preassign the MODEL arithmetic
target eta_cap,model<.001 and the refined guard <.00055 at 1e-6.
Also solve its sufficient cutoff threshold: with
A=q/(2C), B=2sqrt(2/C), d=.001 log(2)/(a_* Cq),

 sqrt(epsilon)<=2d/[B+sqrt(B^2+4Ad)].                (14)

A bound at 1e-3 exceeding the target is an excessive upper estimate,
not a lower bound on its actual error. This classification differs from
the actual-domain obstruction (3). No part of (13) is silently promoted
to an actual curved operator bound or deducted from its physical budget.

## 5. Corrected method and decision boundary

Use a principal comparison P_epsilon defined with the SAME heat cutoff,
metric coefficient and a0/a2 subtraction as H_epsilon. Restore its full
uncut P analytically, so the expression to control is

 (H_epsilon-P_epsilon) + P,

with the previously allocated full normal, local contacts and logarithmic
principal difference subtracted exactly once. The equivalent operation
adds P-P_epsilon, its missing principal cap, BEFORE the strong limit.
This is an exact comparison/addback prescription. It introduces no new
force or physical counterterm. The chosen quantization's lower symbols
and commutators remain in H_epsilon-P_epsilon. Equations (6)-(14) show
what this prescription must include and certify its frozen angular cap;
they do not prove its full variable-coefficient convergence.

RESEARCH DECISION: retain the construction and the geometric Duhamel
representation, but replace raw cutoff convergence by this matched
regulated-principal subtraction. The previous proposed norm-Cauchy step
was not attainable as written. The revised method has no such missing-log
obstruction; its actual curved strong constants remain UNKNOWN. This
updates the active method and preserves the closed feasibility audit.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the actual
H_epsilon-P_epsilon in the coupled H3 graph uniformly through regulator
removal on the cubic paired collar, including coefficient derivatives,
second vertex, finite diagonal contacts, all thermal images and original
global coupling. Test the complete canonical remainder against .18,
leaving .00368 for other unknown global feedback. The combined remaining
target .18368, joint affine <=7.49e-12 and known-source <=1e-14 still
have no actual certificate. Their conditional consequence remains
E<9.946377e-12. Independent Mellin Efar<=1e-10 is also open.

The prior sufficient projection enclosure E_rest<2.123e-42<1.287e-40
is unchanged and not reevaluated; zero remains its enclosure center,
not a measured signed projection. At rho=1e-8 the nonlinear requirements
remain Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with unproved
allocations Y0<=5e-9, Z0<=.25, Z2<=5e7. No controlled correction,
uniform nonlinear admission, exact compatibility, physical matching,
semiclassical validity/noise or all-sector stability has been established.

Remaining to-do list: actual matched geometric remainder and common limit;
global feedback and joint sources, mass/slow moments and independent Mellin
tail; controlled solution and nonlinear admission; compatibility, physical
matching/causal realization, validity/noise, all-sector stability,
charge/return, opening, transport, preparation, separate creation and the
full stable-wormhole objective.
