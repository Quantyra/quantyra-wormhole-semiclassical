# Finite-seed radial Green function with enclosed bolt data

2026-10-22, E003/S169. Ordinary argument before numerical execution.
The full stable semiclassical wormhole objective remains incomplete.

This implements a missing numerical input to the selected full-domain
solve: radial resolvents of the actual canonical massive square on the
current finite seed. It does not replace the full quantum Hessian by a
few modes. Selected mode kernels, not the renormalized stress, are the
output. The earlier feasibility audit remains closed and the actual
metric correction remains uncomputed.

## 1. Current parameters, retained temperature and a stable polynomial

Set G=1 solely as the length unit; k=2000, q=k^4, c_R=1e22 and
g_R=k^-4 are unchanged. For this first evaluation select the allowed
preparation z_i=.265, without claiming uniform numerical coverage of
the whole [.25,.28] preparation interval. Define

    g+11g^2/(360pi)=k^-10, lambda=g,
    F_l=1-384pi c_R g,
    C_RN=sqrt(1+12pi g k^16/F_l), B_RN^2=(C_RN-1)/(6g).

The exact product is enclosed by the retained analytic theorem:

    |eta|<5.4e-23, -1.7e-29<v<-1.57e-29,
    B^2=B_RN^2/(1+eta),
    C=[C_RN(1+v)+eta]/(1+eta), d=C_RN v/(1+eta).

No adaptive pressure quadrature is substituted for this enclosure.
The original clock remains

    C0=sqrt(1+12pi k^6), b0^2=k^10(C0-1)/6,
    kappa=2z_i q C0/(k b0^3), gamma=kappa B/C.

Let z=B/r, u=1/C, beta=lambda B^2/C, epsilon=d/C. Write

    Q(z)=(1-z)^2 R(z)-tau z^3,
    R(z)=beta(1+2z+3z^2)+u z^2+epsilon z^4,
    Q'(z)=-12beta z^2(1-z)
           +2u z(1-z)(1-2z)
           +2epsilon z^3(1-z)(2-3z)-3tau z^2.

Here Q=z^2 F/C and tau=delta/C is the actual positive thermal mass
coefficient. The proved h<=k^-10 and explicit finite seed formula imply
0<tau<1e-65. This loose enclosure keeps the thermal term; it is harmless
on the evaluated interval 0<=z<=1-2^-10. All parameter balls are allowed
to vary independently when enclosing a cell, which enlarges the result.

In proper canonical distance rho increasing from bolt to outer boundary,

    dz/d rho=-sqrt(Q), a=sqrt(Q)/gamma,
    gamma0_metric=d rho^2+a(rho)^2 dtheta^2, period(theta)=2pi.

The canonical outer boundary is z=0 and is geodesic, Q'(0)=0. The
metric has its one-sided mass cubic; no smooth even-double expansion
is assumed. Factoring (1-z)^2 avoids subtracting large almost equal
polynomials in the cold region.

## 2. Actual radial square and the paired outer boundary

After the radial half-density transformation the canonical orbit Dirac
expression is gamma_rho partial_rho+gamma_theta a^-1 partial_theta.
For thermal spin number omega=j+1/2, j integer, its positive square
splits into the two chiral radial expressions

    P_+/-=-partial_rho^2+V_+/-,
    V_+/-=M_n^2+W^2 +/- W',
    M_n^2=n(n+q)/C, W=omega/a=p/sqrt(Q), p=omega gamma,
    W'=p Q'/(2Q).

This is the square of the retained paired Dirac realization, not a
mass-dependent Robin replacement. The exact balanced image identity in
[the boundary-layer derivation](stationary-massive-boundary-layer.md)
identifies the half-domain spin trace with that on the even double.
The double is realized by its closed form, with regular bolts at both
ends. Reflection swaps V_+ and V_- since a is even and W' is odd.
The trace identity does not require smoothness through the cubic seam.

For numerical coverage choose angular n and scaled frequencies
(n,p_target)=(1,1),(1,8),(2,1). The ACTUAL discrete thermal numbers are

    j=floor(p_target/gamma_classical), omega=j+1/2,
    gamma_classical=kappa b0/C0.

The integer j is evaluated with an unambiguous Arb floor, and is retained
exactly. The actual p=omega gamma is enclosed separately; it is not set
equal to p_target. These large j sample finite proper frequencies at
the outer end, where the circle is very large. They are not low thermal
indices, and no integral replaces the discrete sum.

For the whole physical interval, 0<=z<z_h<1,

    |Q'|<=12 beta*(4/27)+u/2+6|epsilon|+3tau<.3.

The selected actual p>.99. Hence p^2 +/- p Q'/2>0 and V_+/- >0
everywhere, including the omitted inner segment. The regular bolt
solutions have positive growing Frobenius branch (omega is large);
u''=V u then keeps u and u' positive. Consequently their logarithmic
derivatives at z_cut=1-2^-10 lie in [0,+infinity]. No number for that
unknown bolt impedance is guessed.

## 3. Validated propagation, without a numerical inner boundary condition

For any regular radial solution write m=u'/u. In outward proper length

    m'=V-m^2.

If V is between positive constants v_lo,v_hi on a cell, scalar
comparison bounds m by the two constant-potential flows with its lower
and upper initial values. At positive length L the exact flow is

    T(v,m,L)=sqrt(v)[m+sqrt(v)tanh(sqrt(v)L)]
                          /[sqrt(v)+m tanh(sqrt(v)L)],
    T(v,+infinity,L)=sqrt(v)/tanh(sqrt(v)L).

This flow is increasing in v and m. Its length dependence need not be
monotone, so it is evaluated on the entire enclosed length interval.
Each outgoing lower/upper endpoint is rounded outwards. At the first
cell the lower datum is zero and the upper datum is infinity. The
result encloses the actual regular-bolt solution, rather than changing
that solution to a finite-wall problem.

The proper cell length obeys

    L=integral_(z_right)^(z_left) dz/sqrt(Q(z))
      in (z_left-z_right)/sqrt(Q_cell).

The mesh has d=1-z from 2^-10 to 1, with d_i=2^(-10+i/N_oct).
Cell endpoints and all transcendental operations are enclosed by Arb.
For each chirality first propagate from the inner cut to z=0. Reflection
sets the OTHER bolt's inward impedance at z=0 to the opposite
chirality's just-computed outward impedance. Propagate that inward
impedance back through the same cells in reverse order. The radial
Green diagonal is then

    G_+/-(rho,rho)=1/[m_left,+/-(rho)+m_right,+/-(rho)].

Both denominator terms are positive. This constructs the correct seam
matching and includes all possible regular data on both omitted bolt
segments. At the outer seam the spin trace is exactly
2/(m_left,+(0)+m_left,-(0)); the code checks this independently.

For two observed points A closer to the bolt than B, the full-double
radial coefficient also obeys

    G_+/-(A,B)=G_+/-(B,B) exp[-integral_A^B m_left,+/- d rho].

The amplitude ratio on a cell with u(start)=1 is enclosed using
cosh(sqrt(v)L)+(m/sqrt(v))sinh(sqrt(v)L). This follows by comparison
of the positive first-order system (u,u'). Accumulate its logarithm
between the observation points; the infinite starting impedance is
outside that interval. This also supplies a separated radial kernel
input, not merely a diagonal value.

The diagonal spin trace per unit canonical area of ONE thermal number
and ONE nonzero angular mode, before its angular multiplicity, is

    [G_+(rho,rho)+G_-(rho,rho)]/(2pi a(rho)).

The massive angular multiplicity is still 2(q+2n). Neither that sum nor
the thermal sum, mass weighting, heat finite part or local restoration
is performed here. In particular, these finite radial coefficients are
not renormalized stress values or action-Hessian matrix elements.

## 4. Predeclared execution and interpretation

Use python-flint 0.9.0 / Arb at 256 bits, installed in an isolated
task dependency directory. Record its binary hashes, Python/platform
identity and invocation. The [Arb wrapper documentation](https://python-flint.readthedocs.io/en/latest/arb.html)
specifies ball arithmetic and outward lower/upper endpoints. Arithmetic
validation relies on that implementation; no Lean certification is
claimed. There is no rounding of a point computation into an asserted
interval afterwards.

Run N_oct=4096 and 8192, with observations z=7/8,1/2,0 and separated
pair (7/8,1/2). Require finite positive enclosures, agreement of the
seam identity, overlap of both meshes and a smaller relative width on
refinement. Target final relative widths below .005 for the diagonal
and .05 for the separated coefficient. These are precision targets
for finite-mode inputs, not the much stronger metric-solution target.
Verify the flow against an independent constant-potential transfer
matrix, and the factored geometry/square reduction algebra. Preserve
all runs, failures, frozen sources and complete numerical output.

The first frozen run used N_oct=512 and 1024. Its enclosures were
positive and refined, but the n=1,p_target=1 separated coefficient
still had relative width .0575192 at 1024, exceeding the .05 target.
The p_target=8 coarse width was .408093. Retain that failed run and
increase the mesh resolution; do not weaken the precision targets.
The final schedule above supplies the next fully frozen execution.

A successful run closes a concrete implementation gap: it evaluates
selected actual finite-seed radial square resolvents over the accepted
parameter enclosure with the actual paired realization. It does not
prove a usable bound on the full sum by extrapolating a few modes.
Next use this evaluator in the full subtracted angular/thermal response,
including its first-order vertices and contacts, with an enclosed tail.
Only that assembled response can enter H e=-r1 and its differentiated
error/admission estimate. No cold source optimization or new generic
estimate chain is commissioned by this step.

The existing closed-audit thresholds, physical boundary/compatibility
requirements, semiclassical validity and all-sector stability remain.
S169.2-S169.5 stay unchecked and S166 remains the last completed component.

Remaining to-do list: full subtracted quantum response and tail enclosure;
actual finite-domain correction with weighted error and nonlinear
admission; exact compatibility, physical matching, semiclassical validity,
retarded stability and all remaining stable-wormhole requirements.
