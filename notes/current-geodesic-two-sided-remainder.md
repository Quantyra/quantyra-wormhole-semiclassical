# Two-sided finite-seed remainder after complete geodesic subtraction

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: progress. The full stable-wormhole goal is incomplete;
the feasibility audit remains closed and its physical thresholds fixed.

The failed heat-cutoff method omitted a full logarithmic principal term.
Retain it instead through the GLOBAL geodesic form of the existing curved
quantum argument. That argument contains absolute kernel/contact estimates,
not only its final one-sided TT inequality. This note extracts a two-sided
comparison in the actual action-form norm, with numerical constants at the
current finite seed. It does not compute a forward or adjoint solution.

## 1. Exact comparison, finite fourth contact and unchanged variables

Keep the full Dirac tower, q=2000^4, c_R=1e22, fixed mu0=1/b0, original
temperature, paired HH realization and physical mass/log phase. Use the
same frozen linear York return and affine canonical action chart as
current-finite-local-action-assembly.md. The canonical massive determinant
depends only on the shape T; hence replacing that block leaves the complete
scalar and scalar/TT mixed rows unchanged. Do not differentiate a moving
York map or change nonlinear tensor charts without their off-shell contacts.

On the physical-half unit tangent bundle of the actual C2,1 even double,
let X be the complete geodesic generator and P=iX its self-adjoint unitary
flow generator. Its Liouville measure is preserved. With the base lift
T -> T(x), retain the exact form from stationary-curved-quantum-tt.md:

    B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2),
    g(s)=1-s^2/2+s^4/24-cos s,
    P_M[T]=integral_0^infinity r[-B_M(r)] <T,g(rP)T>_SM dr,
    Qgeo[T]=sum_(n>=1) d_n P_Mn[T]+Q4geo[T],
    Q4geo=-C F_eff ||X^2 T||_SM^2/(90pi^2),
    F_eff=4 gamma_E-2H_q+4log(mu0 B).                 (1)

There is no finite angular, thermal, flow-time or frequency cutoff in
(1). The bolt is a smooth disk point, not a puncture boundary. The paired
half normalization is the same as in the actual determinant. Positivity
and the Bessel moments follow from the ordinary integral derivation in
the retained proof, using [DLMF 10.32.9-10](https://dlmf.nist.gov/10.32).
Those Bessel identities do not supply a physical domain theorem.

Define H0 to contain all the already assembled local terms, exact LLL,
the actual local canonical a2 Hessian, a shared bounded trial vacuum
coefficient c0*, and sum d_n P_Mn, pulled back by the same fixed York map.
Allow any c0* in the known interval |c0*|<=U0; a midpoint implementation
must carry that uncertainty. It never declares the physical c0 zero.

This uses the latest local assembly without double-counting Q4geo. Weak
Bochner on the physical phase gives EXACTLY

    Ha2-Q4geo = C F_eff/(72pi) integral Rbar |grad T|^2
       -C F_eff/(960pi) integral (4Rbar^2+8Delta Rbar)T^2. (2)

The fourth coefficient agrees because 1/60-1/360=1/72 after integration
by parts. The seam distribution in Delta Rbar is retained; its quadratic
pairing vanishes for T(0)=0, and the physical boundary representative uses
both T(0)=T'(0)=0. Extend by the same form-core limit, not a smooth-double
heat theorem. The finite fourth contact required by the bolt parametrix
is therefore unchanged. A mismatch there could not be bounded as the
lower-order error below.

The comparison is the geometric-flow form (1), NOT Kflat(sqrt(Laplace)).
The previously certified 2e-5 scalar spectral quadrature cannot simply
be applied to P: its flat coefficient already includes an angular average.
An implementation must evaluate (1), or prove the appropriate equivalent
flow multiplier and retain the compression by the base lift.

## 2. What the actual finite-kernel argument bounds absolutely

Read stationary-curved-quantum-tt.md sections 2-6, 8-11 and their actual
finite-seed transfer in current-seed-full-hessian-realization.md. The
following are inputs from those ORDINARY continuum arguments, not from
their supporting checker outputs:

* The paired first-order inverse has norm <=1/M and preserved grading.
  Direct/reflected spin traces and their physical rank are kept. On
  actual seam-crossing C2,1 normal caps, |v|<=8r^2, sup|v|<=.026 and
  |v-T2v|<=1e4r^3. The covariant Born comparison gives an absolute
  separated remainder <=1e8 C/d^4, or <1e9 C/|Y-Y'|^3 radially.
* The two-buffer geometric resolvent identity retains the actual global
  state and domain. Its local/far rows give a bubble Schur reserve
  3e24 C. The ROW-SUBTRACTED gluing form costs <=6e24 C J0; the full
  constant probe already includes its row sum. Neither term is erased.
* The constant probe includes tadpole, bubble, common anomaly and finite
  parts. Its actual-to-cap mass dependence is bounded holomorphically
  before Cauchy differentiation, rather than differentiating a real
  source bound. The final contact is <=1e20 q^2(1+1/Y) pointwise and
  costs <=1500e20 Gq^2/(lambda B^4) E0 on physical TT fields.
* A compact quadratic probe fixes the second contact. The local Born
  expansion, global geodesic returns and anchor change fit 1e13 C.
  The original fourth finite part and two-dimensional anomaly are
  included before this comparison. An off-diagonal kernel envelope
  alone would not fix these diagonal terms.

These statements use the actual polynomial bounds, mass gap and fixed
positive temperature of the current finite seed. They require neither
fourth seam derivatives nor a smooth family at the singular adjoint.
The actual source/domain argument is essential: a smooth closed-manifold
symbol theorem by itself proves none of the above global estimates.

For a symmetric gain-two remainder with anchor radius R0=1/400 and
K=1e13, the retained finite-difference proof gives, for epsilon>0,

    |F[T]| <= (epsilon C/3)||X^2T||_SM^2
       +K C[1/4+2log(R0/b)]||XT||_SM^2,
    b^2=min(R0^2,epsilon/(4K)).                       (3)

It uses only unitarity of transport and Taylor differences, not derivatives
of the kernel envelope. Set epsilon=1e-6. The gradient coefficient still
fits the existing rounded 1e15 reserve. With the explicitly bounded
contacts/gluing, the retained wider envelopes 1e23 C J1* and 1e25 C J0*
therefore also bound the actual comparison absolutely. They include
unneeded localization/comparison allowances, so using them here is
conservative, not an inference from the final negative-form bound.

Let Rop=H-H0 in the fixed physical action chart. The exact angular Hessian
identity and |Rbar|<=2 give
||X^2T||_SM^2 <=(3pi/4)J2+(pi/2)J1. Thus

    |Rop[xi,xi]| <= epsilon*pi*C J2/4
       +1e23 C J1*+1e25 C J0*+p0 E0+vac E0+curv E0, (4)
    p0=1500e20 Gq^2/(lambda B^4).

The epsilon*pi*J1/6 cost fits the same gradient envelope. The last two
terms conservatively cover subtracting the trial vacuum coefficient and
the explicit known curvature polynomial (2), even if already present in
the wider comparison ledger. They are not counted as improvements.

## 3. Substitute actual geometry and retain the joint positive margin

Use E_eff=F_lambda E0. The actual finite-seed proof retains MORE than its
displayed final lower bound: its positive fourth reserve can be kept while
subtracting the listed lower-order, square and chart losses. Specifically

    l[xi] >= .99 F_lambda E0[xi]+C S2[xi]/(100pi),
    S2=sum_(j=o,i) J2[chi_j T].                       (5)

The remaining E0 coefficient is at least
F_lambda(1-7.2e-8-1.37e-21-1e-15-1e-27)-.003601000003,
which exceeds .99F_lambda at the actual parameters. The C/(100pi)
fourth reserve was already left AFTER the gain-two, scalar-comparison
and localization costs in the ordinary proof; it is not counted twice.
The same canonical chart's full-gradient contacts stay in this ledger.

Let a=1/100 and dtheta=5/a=500 for the retained partition. The exact
identity sum chi_j Delta(chi_j T)=Delta T-|grad theta|^2 T gives

    C J2 <=2C S2+2 dtheta^4 C J0,outer.              (6)

All quantities are in the original physical-half action units. The actual
finite-seed return bounds give the following explicit substitutions:

    j0=1000 G C^2/(lambda B^4),  C J0/E0<=j0,
    j1=2e5 G C/B^2,              C J1/E0<=j1,
    jo=(16a^4/25)j0,             C J0,outer/E0<=jo,
    j0star=3j0,
    j1star=4j1+(3 dtheta^2+.0008)jo.                 (7)

The first two follow from ||T||^2<=500 C0 and ||T'||^2<=1e5 E in
physical r^2dr measure. The outer factor is u^4<=16a^4/25. For J1*, use
sum J1[chi_jT]=J1[T]+int theta'^2T^2 and the actual auxiliary half-density
bound J1[v]<=2J1[chi_oT]+.0008 J0,outer. Its unitary J0 equals the
outer base J0, hence J0*<=3J0. This retains the physical cutoff weights
instead of inserting a global unweighted volume bound.

Let Jbound=604q^2+(1+2|log(mu0B)|)(q^2/6-1/30). The available bounds are

    vac=Jbound*j0/(pi C^2),
    curv=|F_eff|[j1/(36pi)+816j0/(960pi)],
    Berr=1e23*j1star+1e25*j0star+p0+vac+curv
                                  +(epsilon*pi/2)dtheta^4 jo. (8)

For curv use |Delta Rbar|<100 and |Rbar|<=2 on each physical side, with
the seam delta handled above. G=1 is the current unit; lambda=g from
the actual interval constructor. The vacuum bound remains shared and
sufficient; its spectral finite part is not optimized.

Combining (4)-(8), and retaining BOTH positive pieces of (5), yields

    |Rop[xi,xi]| <= eta*l[xi],
    eta=max(50pi^2 epsilon, Berr/(.99F_lambda)).       (9)

Using a maximum is justified by the joint coercivity in (5). Adding
separate inverse bounds would lose this information. Freeze and execute
the target eta<.0022 with the current B,C,g,F_lambda and corrected fixed
scale. Record every component of Berr, including the latest j0/j1/outer
substitutions. This is not a normalization change or a frequency sample.

## 4. Full form operator and the singular adjoint

Rop has only a TT/TT block in the frozen scalar/TT congruence. It is real
symmetric on the physical core. Polarization and rescaling the two inputs
upgrade (9) to |Rop[xi,zeta]|<=eta sqrt(l[xi]l[zeta]); equivalently

    ||H-H0||_(V to V*)<=eta,
    ||(s,xi)||_V^2=a(s,s)+l(xi,xi).                  (10)

This extends uniquely to the actual form closure. The actual paired
domain and TT-preserving mass/fast approximation are those of the
finite-seed realization. The point-source bolt adjoint belongs to V,
even though it is not a smooth L2-domain field. Approximate its source
in V*, then use the form inverse and (10) to pass to the limit. Thus the
remainder acts on that adjoint without an invalid H3 or smooth-bolt
assumption. Its NORM is still not computed.

The combined H0 form is defined with the same common cut/limit; no
individually divergent outer Weyl and geodesic pieces are assigned separate
values. For J_sign(s,xi)=(-s,xi), actual coercivity and (10) give

    Re<J_sign u,H0u> >=(1-eta)||u||_V^2,
    ||H0^-1||_(V* to V)<=1/(1-eta),
    ||H0^-1(H-H0)||_(V to V)<=eta/(1-eta).           (11)

This closes a quantitative LINEAR form-comparison/inverse requirement at
the seed. It is not the admission-space Z0, K_c, K_m or a C5 return bound;
it does not transfer the nonlinear source map to a ball of finite metrics.

## 5. What now decides the physical correction

For e_app,psi_app in V and their H0 residuals d0,db0, keep the unchanged
physical bolt certificate

    |Q0|+eta E P+(Db0+eta P)(D0+eta E)<1.667e-9,     (12)

where E=||e_app||_V, P=||psi_app||_V, D0=||d0||_V*, Db0=||db0||_V*.
Those actual quantities and Q0 remain uncomputed. Even with D0=Db0=0
and half the target left after |Q0|, the sufficient uniform bound
eta_star=.0022 requires

    E P < 8.335e-10/[eta_star*(1+eta_star)].          (13)

If E0app,P0app are computed in the model norm a+l0, use
EP<=E0app P0app/(1-eta_star). Record the resulting threshold as well. These
are feasibility thresholds, not assumed sizes of the physical solution.
If this allowance fails, a direct signed pairing and directional Rop
residuals can improve it; failure of this upper bound is not nonexistence.

The dominant known losses in (8) are the global gluing/zero-order envelope,
finite-contact/gradient envelope and actual constant contact. The singular
principal heat tail has been retained exactly. Optimizing the already
tiny vacuum coefficient or source bounds is not the priority. A useful
sharper estimate must reduce the actual remainder on the chosen pair or
its certified residual product, in the original physical units.

Research decision: continue with this complete geodesic comparison and
the physical weak primal-dual solve. The single next calculation is a
global H0 forward/bolt-adjoint trial pair with certified model residuals,
E0app,P0app and Q0, using the correct full angular flow form (1), exact
local rows and contact-complete adjoint. Test (12) before commissioning
another generic remainder optimization. That solve is PROPOSED, not
executed here. The two-sided operator feasibility bound (9) IS executed.

Exterior <6.67e-10, Y0<=5e-9, the preserved nonlinear neighborhood and
exact compatibility remain open. Matching, physical causal boundary
realization, semiclassical validity/noise and all-sector stability are
still required. No actual controlled wormhole or completion is claimed.

Remaining to-do list: physical comparison solve and two-residual bolt
certificate; X/admission, nonlinear closure and exact compatibility;
matching, causal realization, semiclassical validity/noise, all-sector
stability, charge/return, opening, transport, preparation, separate creation
and the full stable-wormhole objective.
