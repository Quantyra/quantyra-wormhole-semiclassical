# First-order shape mixing in the literal boundary action

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The complete auxiliary first-order shape response cannot supply a
persistent quadratic cutoff oscillation through its MIXED pairing with
the fast conformal lift. The mixed finite-cut coefficient converges to
its finite variational value in translated logarithmic L1. This is a
necessary-limit test for the unchanged action, not a prescription to
average or renormalize the physical action. Pointwise convergence is
not inferred from the energy norm.

The pure conformal lift retains its previous quadratic oscillation.
The pure shape coefficient, including the canonical determinant, and
the second-order trace forcing with first-order shape present remain
unresolved. Thus this result does not extend the preceding obstruction
to every coupled metric, construct a physical mode, or prove stability.

## 1. Exact split and endpoint regularity

Keep the charged seed, fixed background cut, physical time and flux,
paired matter domain, common subtraction and all finite coefficients.
Fix a positive Matsubara frequency above the retained auxiliary shape
and reconstruction thresholds. In the notation of the
[Schur calculation](thermal-curvature-shape-schur.md), write

    p=R_gamma U/r, s=Q_gamma w/r, v=w/r,
    U in D_max(Atilde), w in E_gamma.

The full metric tangent is (p+s,v). This includes the actual auxiliary
response w=w_ell-K_cross^dagger U whenever the admitted shape source is
continuous on E_gamma; no solution of the physical constraints is implied.
Use inward optical y with r_y=-F and y=0 at infinity. The retained
estimates give w,r s in H1_0, and v_yy-gamma^2 v and
(s+v/3)_yy-gamma^2(s+v/3) in L2 on a closed outer collar.
The coefficient r^2/F in J2 is bounded above and below there. Hence
v,s belong to H2 on that collar. Their first two boundary traces vanish:
their already proved o(y^(3/2)) values imply v(0)=v_y(0)=s(0)=s_y(0)=0.
For either h=v or h=s, the absolutely continuous representatives obey

    |h_y(y)| <= sqrt(y) ||h_yy||_(0,y),
    |h(y)| <= y^(3/2) ||h_yy||_(0,y)/sqrt(3).             (1)

In particular h=o(y^(3/2)), h_y=o(sqrt(y)). There is no pointwise
bound on h_yy. These facts use only the established energy estimates.
The distinction between a form domain and stronger operator regularity
is standard; compare [Teschl, section 2.3](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf).
The endpoint bounds (1) follow directly by two integrations and
Cauchy--Schwarz, rather than an unproved trace theorem for h_yy.

The [metric trace identity](coupled-fast-metric-boundary-traces.md) also
improves the pure conformal lift beyond its leading amplitude. Put
d0=E+4lambda. It gives q=r p=U/(6d0)+z with z=o(y^(3/2)),
z_y=o(sqrt(y)). Since K_gamma q=F U/6 and U has its differentiated
maximal-domain expansion, the leading fast metric p_f satisfies

    p=p_f+e, p_f=const*y^(3/2)[A_c cos(k log y)+A_s sin(k log y)],
    e=o(y^(5/2)), e_y=o(y^(3/2)), e_yy=o(sqrt(y)).        (2)

For the last assertion use q_yy=(gamma^2+V_K)q-F U/6, not a
pointwise derivative of the arbitrary L2 graph forcing. V_K=6/y^2+O(1),
F=1/(lambda y^2)+O(1), and the leading terms satisfy the normal equation.
The bounded seed corrections are of higher order. Conversion from
log y to z=log(r/B) changes only the higher-order remainder.

Every cross term between p_f and e in the pure conformal bulk is
integrable at y=0: the worst terms are p_f,yy e_yy and y^-4 p_f e,
both o(1). Its literal boundary cross terms are o(y); the e,e terms
are smaller. The retained fixed-frequency terms are lower order.
Thus the exact charged lift has the same persistent Q_t cutoff term
as the [previous conformal calculation](fast-trace-finite-cut-action.md),
plus a finite constant and o(1). This does not choose a boundary phase.

## 2. Differentiate the original action before estimating

At a fixed finite r use the independent amplitudes alpha,beta:

    g=exp(2(alpha p+beta s))
      [exp(2 beta v)(F dt^2+dr^2/F)+r^2 dOmega^2].

Let L_x and B_x be the alpha beta coefficients of the complete local
bulk density per dr and the outer boundary value. They are a Hessian
pairing, not one half of that pairing. Include Einstein, magnetic and
Weyl terms, GHY, both Einstein counterterms, loop counterterms and the
converted literal Weyl boundary polynomial. The magnetic bulk term is
conformally invariant and has zero mixed coefficient, without deleting
its shape-only contribution. The canonical determinant is independent
of alpha p at every cut, so its mixed coefficient is exactly zero.

For radial primes define

    m2=partial L_x/partial p'',
    m1=partial L_x/partial p'-(m2)',
    E_p=partial L_x/partial p-(m1)',
    E_cut=B_x+p m1+p' m2.                               (3)

Then integration on any finite annulus is exactly

    integral L_x dr + [B_x] = integral p E_p dr + [E_cut]. (4)

All surface terms must remain in (3). The general integration-by-parts
distinction is described in [Iyer--Wald, (30)-(43)](https://arxiv.org/pdf/gr-qc/9403028);
that formalism alone does not determine this action's finite-cut value.

The curvature-free reconstruction is, at zero temporal frequency,

    6F s''+6(F'+2F/r)s'+H_R s
                          +2F v''+2F'v'-C_v v=0,
    H_R=2R+rR', C_v=2F''+8F'/r+4F/r^2.                  (5)

Differentiate (5) once when reducing (3). Exact differentiation of
the original action cancels all third derivatives and all s'' terms.
There is ONE remaining second derivative:

    E_cut=A s+B s'+C v+D v'+J v'',
    J=-F p r(rF'-2F)/(60pi).                            (6)

A,B,C,D are fixed local functions of F through its fourth derivative,
p through its second derivative, G and lambda. Equations (3) and the
retained original densities specify them without a free boundary term;
the checker derives them directly. They satisfy A,C=O(r^(3/2)) and
B,D=O(r^(5/2)) for a fast p. No differentiated remainder of v is used.
For smooth interior tangents this is the literal finite-cut expression.
For a general energy-domain tangent it defines the same mixed coefficient
almost everywhere using its weak second jet; (12) concerns that extension.
This does not assert differentiability of a rough finite-amplitude metric
family or existence of the full determinant on every such family.

There is an independent check on the bulk row. With D_gamma w exactly
as in the retained Schur calculation, including the seed residual,

    E_p = r D_0 w/(2G F).                               (7)

The sign is positive in the literal areal-action convention used by
the preceding finite-cut calculations. The Euclidean Hessian convention
of the Schur note has the opposite sign. In particular the finite bulk
mixed value here is -b_gamma(U,w)=<R_gamma U,D_gamma w>/(2G), not
b_gamma. Both R_gamma U and D_gamma w are in L2. It is finite and
is not assigned zero. For a fixed interior endpoint, its contribution
and that endpoint's flux are simply fixed finite constants.

For a normal-coefficient check, set z=log(r/B), t=G lambda/(60pi),
S0=lambda B^3/G and write z derivatives as dots. The leading endpoint,
divided by S0 exp(3z), reduces to

    -(3+221t/7)p s
    +[(3+11t)p+(33t/7)pdot+3t pddot] sdot
    +[-(3+125t/7)p-(54t/7)pdot-6t pddot]v
    +[(1+11t/3)p-(3t/7)pdot+t pddot]vdot.             (8)

Its normal bulk row is zero under
3(sddot+3sdot-4s)+vddot+vdot-12v=0. Dropping the subleading
charged J term from (6) would nevertheless be unjustified pointwise.

For each fixed nonzero thermal frequency, (5) has the additional term
-gamma^2(6s+2v)/F. It changes only the lower jets in the elimination.
The coefficient J of v'' stays the same. The mixed bulk frequency
polynomial in the [optical calculation](local-coupled-frequency-polynomial.md)
has at most two radial derivatives in its gamma^2 coefficient and no
radial derivatives in its gamma^4 coefficient. Its extra integration
flux therefore uses only first jets. The full boundary anomaly has
derivative weight three: after two tangential derivatives, at most one
normal derivative remains. Periodic time integration removes tangential
divergences. The extra mixed endpoint terms are bounded on the collar by

    C_gamma(1+|log y|){|p|(|s_y|+|v_y|)
          +|p_y|(|s|+|v|)+y^-1 |p|(|s|+|v|)}=o(1).    (9)

This uses the full spherical mixed heat coefficient, including its
tangential-curvature terms, not a static boundary formula as an exact
dynamic formula. The extension and bound (9) are ordinary derivative
counting; the numerical fixtures below are static.

## 3. A necessary limit, without a new cutoff prescription

Since v_r=-v_y/F and v_rr=(v_yy+F'v_y)/F^2, (6) becomes

    E_cut(y)=o(1)+b(y)v_yy(y),
    b(y)=-p r(rF'-2F)/(60pi F)=O(y^(5/2)).             (10)

Indeed rF'-2F=-2+O(1/r), while (1),(2) make all first-jet terms
vanish. The F'v_y term in this conversion is o(y^2). For sufficiently
small d choose |b(y)|<=K y^(5/2). Cauchy--Schwarz gives

    integral_0^d |b v_yy| dy/y
                      <=(K d^2/2)||v_yy||_(0,d) ->0.   (11)

Consequently for any fixed 0<a<b<infinity,

    integral_(a d)^(b d) |E_cut(y)| dy/y ->0 as d->0.   (12)

The finite bulk term in (4) has an ordinary limit by its L2 pairing.
Thus the mixed action tends to that finite value in translated
logarithmic L1. It cannot contain a nonzero persistent sine/cosine in
log y plus a remainder tending to zero in that topology. If its ordinary
cutoff limit exists, it must be the same value. This necessary-limit
argument leaves the original regulator and physical action unchanged.

There is a real reason not to claim a pointwise limit from (11). On
d_n=2^-n, n>=2, put ell_n=d_n^8 and

    v_n(y)=d_n^(27/2) psi((y-d_n)/ell_n),
    psi(x)=x^2(1-x^2)^6/2 on [-1,1], zero elsewhere.

The supports are disjoint, psi''(0)=1, and sum v_n has H2 norm finite:
sum ||v_n,yy||^2 = ||psi''||^2 sum d_n^3. Its zero value and first
derivative traces and the weighted H1 norm of v/y also have the needed
decay. Yet d_n^(5/2)v_yy(d_n)=1, with zero values between supports.
Its logarithmic L1 contribution vanishes. This is a counterexample to
upgrading these Sobolev estimates alone; it is NOT asserted to be the
actual shape solution or a counterexample to its possible stronger
regularity. No evaluation of that full nonlocal solution is supplied.

## 4. What this resolves and what it leaves

For the exact first-order tangent (p+s,v), the unchanged quadratic
coefficient splits into pure conformal, mixed, and pure shape pieces.
Equations (2),(12) settle the persistent part of the first two: the
pure conformal coefficient is unchanged, and the mixed coefficient
cannot cancel it. This conclusion includes every shape in the retained
full energy domain, including the solved auxiliary response. It requires
no pointwise bound on its second or third normal derivative.

The remaining pure shape piece contains the canonical SECOND variation.
Finiteness of its closed form does not yet identify its literal finite-cut
limit. The preceding canonical FIRST-variation theorem for second-order
cubic directions is not a theorem about this new coefficient. Nor have
the nonlinear trace source terms involving first-order shape been shown
to preserve the earlier r^-3 tail equation. Both issues must be resolved
before asserting the full coupled obstruction or selecting physical
boundary data. Leading metric/clock data alone still do not select them.

Verification derives (3),(6),(7) for symbolic F, all finite contacts and
the reconstruction, verifies (8), and checks independent complex-amplitude
extraction of the original charged action. Normal annular integrations
and shrinking Sobolev spikes test the endpoint and necessary-limit steps.
These are supporting fixtures, not evaluations of a wormhole or the full
quantum shape inverse. No Lean run or physical completion is claimed.

Remaining to-do list: pure shape quadratic finite-cut limit and complete
second-order trace forcing; physical boundary phase, source compatibility,
continuous causal constraints and stability; exact nonlinear existence;
exterior/state, charge/return, opening, transport, preparation, separate
creation and full audit.
