# Actual paired scalar response: uniform mode return and the unsummed strong gap

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
The full stable semiclassical wormhole objective remains ACTIVE/INCOMPLETE.
The closed feasibility audit is unchanged.

OUTCOME: the complete coupled graph estimate is UNRESOLVED. This attempt
derives an exact positive scalar boundary problem for the complete mass
response of EVERY angular/thermal mode on the actual finite seed. Its
uniform inverse has a numerical bound, with no inverse-temperature loss.
This controls an auxiliary response, not the commonly subtracted quantum
Hessian or its coupled physical inverse. The lowest thermal mode also
shows why an unregulated termwise strong-norm estimator fails at the bolt.
No signed coupled solve, actual mass or kappa_m is certified.

## 1. The missing step and the hypotheses actually used

The [previous far result](current-global-far-return-energy-bound.md)
controls the complete actual input z>=.2 to original output z<=.1. The
first missing NUMERICAL continuum step is the differentiated, commonly
subtracted near/transition response, composed with the actual coupled
inverse and original matching/bolt rows. The existing fixed-seed form
estimate and admitted-family C2 regulator identification do not bound that
strong composition. The short-travel method has not supplied its higher
kernel/contact constants. No closed smooth-manifold theorem or analytic
moving-cutoff assumption supplies them on the actual paired mass/log domain.

Here the first-order paired equations themselves replace that unsupported
numerical step by an exact auxiliary problem. Use the fixed one-sided
polynomial Q, f=sqrt(Q), dz/d rho=-f, and rho increasing from the regular
bolt to the outer seam. The other half is its actual identical paired
copy, not an imposed wall. All M=M_n=sqrt(n(n+q)/C), n>=1, d_n=2(q+2n),
and omega=j+1/2, j>=0, are retained; p=gamma omega and W=p/f. The separate
massless/Weyl allocation and all common local restorations remain in H.

The [Dirac assembly](current-coincident-dirac-action-assembly.md) gives
the two actual positive Riccati solutions on this half:

    U'=W(1-U^2)-2M U,    V'=-W(1-V^2)+2M V,
    0<U,V<=1,           U(seam)=V(seam).                 (1)

U starts on the regular bolt branch; V is propagated from the OTHER
regular bolt through the paired seam. For a canonical shape direction v,
the mass family is m_s=M exp(s v), with the same paired projectors and
period. The mass/log directions at the outer seam, the cut mass lift,
scalar quartic, correlated slow pair and fast directions are allowed.
No smooth surrogate across the double or finite selection of modes is
used. This is a SEED derivative for bounded admissible directions; its
linear L2 extension below is not a nonlinear family of L2 metrics.
Independent gravitational boundary-phase directions are retained as
inputs v. Varying the fermionic transmission projector would be a
different calculation and is not asserted here.

Only elementary ODE comparison, differentiation, and a positive closed
quadratic form are used below. Their actual endpoint hypotheses are
proved in sections 2-3. The physical Dirac square remains the composed
square on the complementary adjoint domains from the
[regulator audit](current-mass-response-regulator-audit.md). The scalar
operator below is a response identity, not a replacement Dirac square.

## 2. Exact scalar response, tadpole and paired boundary

Define

    t=(1-UV)/(1+UV), a=(U+V)/(1+UV), b=(U-V)/(1+UV).

Then a>0, 0<=t<=1, a<=1, and direct differentiation of (1) gives

    t'=2W b, a'=-2M b, b'=2W t-2M a,
    t^2+a^2-b^2=1.                                    (2)

Put y=delta t and delta m=M v. The differentiated constraint gives
delta a=(b delta b-t y)/a, while y'=2W delta b. Elimination is exact:

    y''=(W'/W+a'/a)y'+(4W^2+4MWt/a)y-4MWa v,

    L y=4M v,
    L=-partial_rho[(Wa)^(-1)partial_rho]+4W/a+4Mt/a^2.  (3)

The unregulated relative source per thermal sign is F_M=-M t. Its
COMPLETE mass derivative, including the exponential-mass coincidence
tadpole, is therefore

    H_M v=-M t v-M y=-M t v-4M^2 L^(-1)v.              (4)

Both thermal signs have this same expression. An equivalent auxiliary
quadratic form is

    Q(v,y)=-integral M t v^2 +(1/4)<y,Ly>-2M<v,y>.

Its symmetric block is [[-Mt,-M],[-M,L/4]], and its Schur complement is
(4). This is a local auxiliary representation of the whole bubble plus
tadpole for each mode, including finite-domain propagation. It is not an
adjustable source. The massive determinant has zero canonical trace row;
the original trace/scalar/reference equations and their mixing remain in
the coupled gravitational system.

At the seam, both U=V AND their variations agree. Thus b=delta b=0,
and the derived scalar condition is y'=0, equivalently (Wa)^(-1)y'=0.
Crucially, y(seam) is NOT zero or prescribed. Its value contributes to
the original scalar quartic finite-part row. Taking D-4 first would lose
that row; the scalar reduction does not authorize doing so.

At the actual bolt, f=gamma rho(1+O(rho^2)), W=omega/rho+O(rho), a->1.
The leading endpoint expression is

    L=-partial_rho[(rho/omega)partial_rho]+4omega/rho
      + lower endpoint orders.                        (5)

The positive form is the closure of
integral [|y'|^2/(Wa)+(4W/a+4Mt/a^2)|y|^2] d rho, with
test functions vanishing near the bolt but free at the seam. It imposes
the finite-form branch at the bolt and the natural condition at the seam.
The homogeneous branches are rho^(2omega) and rho^(-2omega); the latter
is not L2 for omega>=1/2. The actual differentiated regular Riccati
solution belongs to the former form domain, including the inhomogeneous
terms below. This identifies its endpoint realization; it is not a
choice among arbitrary scalar boundary conditions.

For a smooth stationary input v=v0+O(rho^2), direct leading-order
substitution in (3) gives

    y=4M omega v0 rho/(4omega^2-1)+o(rho), omega>1/2;
    y=-M v0 rho log(rho)+O(rho),             omega=1/2. (6)

Indeed the leading operator maps rho to (4omega^2-1)/omega and maps
rho log(rho) to -4 when omega=1/2. Corrections to (5) are lower order
even at this resonant mode. Both particular solutions have finite form
integral, since integral rho log^2(rho) is finite. Cutoff approximation
converges in that form. The singular homogeneous branch does not.
Bounded mass/log inputs are smooth at the bolt in the admitted core and
L2 at the outer seam; density extends the response to arbitrary L2 input
after the coercivity estimate below. No derivatives of a rough input
or change of the physical domain are needed for this extension.

These statements also follow by differentiating the positive first-order
solutions, whose regular endpoint expansions select (6). Uniqueness in
the positive form identifies that derivative with (3). Consequently the
Neumann seam condition and bolt branch retain the original paired return.

## 3. A global all-mode bound, including the lowest thermal frequency

Let sigma=sqrt(1.05), d=M-sigma/2>0. On the ACTUAL seed,

    |f'/f|<=gamma/f+sigma,
    |W'/W|<=2W+sigma.                                  (7)

For z>=.2, Q'<0, Q''<=2.1, Q(z_h)=0 and Q'(z_h)=-2gamma imply
Q'^2<=4gamma^2+4.2Q by integration. Since f'/f=-Q'/(2f),
the first inequality follows. On 0<=z<=.2, Q>.16 and |Q'|<.3
give it directly. The second uses gamma/f=W/omega<=2W, so it
covers every actual discrete thermal mode, including omega=1/2.
It is not asserted for arbitrary continuous frequencies below gamma/2.
The argument works in either orientation along the paired interval.

Comparison in (1) now gives the global supersolution cone

    U,V<=min(1,W/d).                                   (8)

For W<=d write r=W/d. If F(r)=W(1-r^2)-2Mr, (7) gives
r'-F(r)>=W(r-1)^2>=0. At W=d the cap is 1>=U,V;
near the bolt the cap is 1 and F(1)=-2M. This proves comparison
on every excursion, without assuming monotonic W. The positive
solutions propagated from either bolt satisfy the same cone.

For W<=d, monotonicity of a and t in U,V gives

    a<=2r/(1+r^2), t>=(1-r^2)/(1+r^2),
    4W/a+4Mt/a^2 >= M/r^2+2d+(2d-M)r^2 >=4d.           (9)

The last difference is (1-r^2)[M/r^2-(2d-M)]>=0.
For W>=d, a<=1 and t>=0 give the same lower bound 4d.
It follows on the actual endpoint domain that

    L>=4d I,
    ||L^(-1)||_(L2(d rho)->L2(d rho)) <=1/(4d),
    ||y||_2 <=(M/d)||v||_2,
    integral |y'|^2/(Wa) <=(M^2/d)||v||_2^2.           (10)

For completeness, the form is positive and locally uniformly elliptic
away from the bolt. Its completion embeds continuously in L2 by (9).
The linear functional 4M<v,y> is continuous there. Minimization of the
strictly convex quadratic form minus this functional gives the unique
weak solution; the energy identity gives (10). Local ODE regularity
recovers (3) and the natural seam condition. This proves the inverse
claim without importing a smooth compact-manifold boundary theorem.

With the actual M_1=18.048130263729889289300... and M_n>=M_1,

    ||L_nj^(-1)|| < .014256564,
    M_n/(M_n-sigma/2) < 1.029218                       (11)

UNIFORMLY for all n,j. These are upper bounds, not measured norms.
The energy bound in (10) keeps its explicit M-dependent coefficient.
One further coefficient bound useful for a later derivative estimate is

    4W^2+4MWt/a >=1.94(M^2+W^2).                       (12)

For W<=d its lower bound is 2Md+(4-2M/d)W^2; both normalized
coefficients exceed 1.94 for M>=18. For W>=d use 4W^2 and
4d^2/(M^2+d^2)>1.94. Equation (12) is not a high-output-derivative
estimate or a coercivity claim for the physical gravitational Hessian.

## 4. Pay the original norm conversion; attempt the complete sum

The original energy E and canonical measure remain unchanged. The
[input-energy calculation](current-coupled-input-energy-bound.md) and
previous far bound give

    ||T||_L2(dz) <= C_T E(u), C_T=sqrt(500/(.99 F_l))<22.474.

For the entire near input v=1_[0,.2] T, f>.4 gives

    ||v||_L2(d rho) <35.534 E(u),
    ||y_nj||_L2(d rho) <36.573 E(u),
    ||y_nj||_L2(dz;0<=z<=.1) <23.485 E(u).             (13)

The last bound uses dz=f d rho and f<=sqrt(.17). The unrounded
upper expressions are 35.533452935244388..., 36.571643098960116...
and 23.483168979569695..., respectively. These global auxiliary
returns include the actual bolt and both paths for every mode. They
are a quantitative alternative to propagating unknown high kernel jets
or paying a spurious inverse-gamma cost at this L2 stage. They are NOT
original H3 graph or physical C5 estimates. For E(u)<1.1e-10 the
second guard gives <4.024e-9 per auxiliary mode, not a summed source.

The sharp input split is only a bounded insertion in a linear operator.
It can be approximated in L2; no discontinuous metric is introduced.
When differentiated contacts and both pieces are recombined, interface
terms cannot be omitted. The original fixed smooth cutoffs and all their
transition derivatives remain in the complete equations.

Restoring the original action per half, -pi B^2/(G gamma) in dz,
the pair of signs at angular level n contributes the canonical shape row

    2G gamma d_n M_n/(pi B^2 f) [t_nj v+y_nj].         (14)

There is no extra mouth/spin factor. The trace row of this massive part
is zero, but all other original rows remain. On the outer collar,
(13) and t<=1 give an L2 upper expression at most

    [229.795 G gamma/(pi B^2)] d_n M_n E(u)             (15)

for EVERY j. Attempting to sum (15) over the COMPLETE angular and
thermal sets gives infinity, already for one angular level and all j.
This is a verified NONSUMMABLE upper-bound estimator, not divergence
of the renormalized physical operator. No finite list of modes can repair
that estimator or certify any part of the remaining graph budget.

There is a second precise ordering failure at the actual bolt. From
(6), y/f=-(M v0/gamma)log(rho)+O(1) at omega=1/2. The raw t term
has the same leading logarithm: t=-M rho log(rho)+O(rho). Thus the
full unregulated expression (14), including its tadpole, has leading
term -4G d_n M_n^2 v0 log(rho)/(pi B^2). The estimate that sums
individual unregulated source sup/C5 norms in the bolt chart fails even
before its thermal tail is considered. This does not show that an
individual source composed with the coupled inverse is unbounded.
Admissible interior directions with v0!=0 exist (the compact interior
example in the [coupled audit](current-coupled-global-return-attempt.md)).
For the near-supported input in (13), v0=0: this particular resonance
does not disprove its outer return. It does disprove that estimator for
the COMPLETE coupled domain including independent bolt directions.

The logarithm of an individual unregulated thermal sector is not a
physical bolt singularity of the full commonly renormalized action.
The complete sum and common local restoration must precede the strong
norm and bolt limit. At a positive common heat regulator the response
uses the Duhamel vertices of the actual composed square. It is NOT
obtained by multiplying (4) by a guessed e^(-epsilon M^2) or thermal
weight. Equations (3)-(4) are exact per-mode relative response identities;
interchanging their unregulated sums, strong derivatives and endpoint
limits is precisely what has NOT been justified numerically here.
The retained admitted-family C2 action identity allows its weak use
within that scope; it supplies neither uniform C5 convergence nor a
strong neighborhood estimate. No contact or finite restoration is set
to zero in making this distinction.

## 5. Complete coupled ledger and what remains unknown

The [exact reference mass lift and coupled border](current-coupled-global-return-attempt.md)
are retained. With J_m=int_0^z t^3/Q^2,

    n_m=-2 beta Q J_m/z, b_m=beta z^3/Q+beta Q'J_m,
    f_uncut=2xz, alpha_uncut=0,
    f_cut=2xz chi_m-2F n_m D chi_m,
    alpha_cut=n_m D chi_m.                             (16)

Both local reference constraints cancel. The actual cut field vanishes
before the bolt; its reference bolt shift is zero and net mass is -xB,
including the transition. The uncut tangent changes fixed temperature
by about -6.515e65 per unit m and is not a global zero mode. Equation
(3) applies to the shape part of every such column, retaining the
transition; it does not evaluate the quantum column after elimination.

Keep the original two equations and the block Schur correction
Lbar_xy=L_xy-L_xI P_I L_Iy. The scalar quartic row still has leading
20 beta Fstar c4=3.33333329399556... c4; its other finite parts are
unknown. Keep the exact mass moment and correlated B0(z^4-z^(4+e_*))
with the exact implicit slow root, not two independently enlarged
amplitudes. Neither y(seam) nor the mass/slow/bolt rows are dropped.

| Requirement | Executed evidence / remaining status |
| --- | --- |
| All-mode auxiliary L2 return with actual endpoints | (3), (10)-(13), proved at the fixed seed |
| Complete common-subtracted near/transition strong response | UNKNOWN; (15) fails as a summable estimator and the bolt resonance excludes termwise unregulated source C5 |
| Actual strong coupled inverse and post-elimination columns | UNKNOWN; an auxiliary L inverse is not P_I or H inverse |
| Actual kappa_m<=5e-5 | UNKNOWN, not a failed bound or an assigned coefficient |
| Other scalar/slow/fast and interior feedback entries | UNKNOWN; no row has been certified to fit by this calculation |
| Original-fast basis conversion | ||d_chi||H3>.00906808164; two 1e-11 ceilings require actual |m|<2.20553814952e-9 |
| Actual mass | UNKNOWN; 4.2953e-8+49.462 E_graph is a conditional allowance, not its value |
| Complete far input to original outer graph | Previously proved <4.583e-13 E, <5.042e-23 for E<1.1e-10 |

IF kappa_m<=5e-5, its debits remain 2.14765e-12 affine,
.0024731 graph feedback and .0195239 independent energy feedback.
After the known far guard the remaining sufficient requirements are

    Y_remaining <=5.34234999994958e-12,
    eta_remaining+(.0195239+4.583e-13) gamma_E <.0012069, (17)

where E(delta u)<=gamma_E ||delta g|| is needed on the ACTUAL
source-controlled solution family. Even with eta_remaining=0,
gamma_E must be <.06181654279993595. Its actual value is unknown;
on the unrestricted outer-data-only domain no finite such comparison
exists, by the retained interior kernel. Bounds (13) do not close (17).

The conditional fast screen remains

    eta_seed+1e-11 C3 < .1599164;
    .028913+1e-11(1.0706e10)=.135973 < .1599164.         (18)

This retains ONLY the proved two-fast-input admitted-family scope. It
does not control independent phase variations, produce the actual
coupled coefficient or transfer the seed estimate to a nonlinear ball.
If alpha is the required physical-to-direction comparison, R_X the
return constant and E_center its center offset, the original requirements
still include

    E_center+1e-8 alpha <=1e-11,
    R_X alpha*.028913+Z0_other <=.25,
    R_X alpha^2*1.0706e10+Z2_other <=5e7,
    Y0+1e-8 Z0+5e-17 Z2 <=1e-8,
    Z0+1e-8 Z2 <1.                                    (19)

None of those missing comparisons or inequalities is certified here.

## 6. Decision and one next calculation

EXECUTED: exact scalar response and endpoint realization for every
actual mode; uniform inverse/energy constants and original L2 costs;
failure of the raw norm-sum estimator and its actual lowest-mode bolt
regularity test; unchanged conditional coupled ledger. These are ordinary
arguments supported by frozen symbolic and 512-bit interval arithmetic.
No physical mode sum, signed coupled center or finite-power operator
application has been executed, and no nonlinear solution gate is closed.

RESEARCH DECISION: retain the construction and sharpen the method using
the exact paired scalar resolvent, with common subtraction BEFORE strong
norms and bolt limits. Its uniform auxiliary conditioning removes a
specific unnecessary loss, but quantitative feasibility of the actual
coupled return remains undecided. A nonsummable raw upper bound or a
single-sector logarithm is not an obstruction to the physical construction.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete
common-subtracted scalar-resolvent sum and its required output
commutators AFTER composition with the actual coupled inverse, retaining
both original equations, the original quartic/mass/slow/bolt rows and the
known far entry. The output is the numerical composite
T_old(H_s^(-1)R_s)^4 on all columns, with the actual energy input norms
and common-regulator error. Test it and the signed centers against (17)
and the actual mass conversion. This is the same missing global estimate
in a more explicit representation, not an already executed solve.

Remaining to-do list: complete common-subtracted strong coupled return
and signed response; original physical-norm/neighborhood comparisons;
nonlinear admission/contraction and exact compatibility; physical matching
and causal realization; semiclassical validity/noise; all-sector stability;
charge/return, opening, transport, preparation, separate creation and the
full stable semiclassical wormhole objective.
