# The full mixed return and the domain of the Schur residual powers

2026-10-25 HST, E003/S169. Ordinary argument and execution contract before
calculation. Previous goal turn: NO PROGRESS. Full stable-wormhole goal
active/incomplete; closed audit preserved. This attempts the domain step
of the COMPLETE O_graph/phase P(I-SP)^2, not an isolated quantum contact.
Its numerical norm and the signed coupled solution remain uncomputed.

The candidate and physical forms are those of current-schur-form-
preconditioner-test.md: H_Y=[-a,b;b*,ell], B=b*a^-1b, q=.338,
P=V ell^-1 J_F, E=I-SP. All operators use the actual paired mass/log
domain, common-restored full Dirac theory and fixed finite seed. ell is
the FULL quantum TT form. It has not been numerically applied here.

## 1. Exact mixed operator: include the vacuum cancellation first

Write D_z=z partial_z and use the original York auxiliary variables:

 D_z J=z(b_r+xi)/Q,
 D_z b_r=(4-zQ'/Q)b_r-zQ'xi/Q,
 U xi=xi/2+QJ/z, V xi=-(b_r+3xi+Q'J)/2.

They retain J(0)=0 and the original regular-bolt branch. With H_loc the
exact ORIGINAL dz Hessian from signed_local_core.py, the physical mixed
row multiplied by 8piG A_radius^2 is

 b_hat xi=-z^4[H_loc,00 U xi+H_loc,01 V xi].           (1)

The canonical nonlocal determinant has zero conformal/mixed derivative;
its entire anomaly and density contacts are already in H_loc. Thus (1)
uses the full mixed row. No angular or thermal mode is omitted by this
Ward reduction. The independent full quantum ell is not removed.

In the inherited G=1 units the exact vacuum identity is
lambda=C beta/B^2, lambda0=lambda+11lambda^2/(360pi).
It must be substituted in the SUM, not separately bounded in the Einstein
and Weyl terms. The seed Q is the actual polynomial, including tau and eps.
Repeated use of the two auxiliary identities reduces (1) to

 b_hat xi=sum_(j=0)^2 c_j D_z^j xi+c_b b_r+c_J J,     (2)

Retain rho_a=8piG A_radius^2: b=rho_a^-1 b_hat. In particular the
PHYSICAL adjoint is rho_a^-1 times the sum of the differential adjoints
(3-D_z)^j(c_j v) and the two integral transposes below. All compositions
B, F and all form norms in sections 4-5 use the physical b, not b_hat.

with one-sided analytic coefficients at z=0. The whole-expression test
below checks absence of derivatives three and four, and outer orders at
least three for c_j,c_b,c_J. These statements concern the COMPLETE sum.
Terms that look singular at the bolt in this coordinate expression are
kept combined in the original smooth tensor/disk expression; separate
coefficient bounds there are not inferred.

The derivative order also follows directly from the trace identity:
delta R[k]=-4 tau_R xi, so -D_R Delta(delta R) uses at most two derivatives
of xi. Curvature-anomaly variations have at most two tensor derivatives;
the density/gauge return is lower order. At infinity tau_R=O(z^4),
C4=O(z^3), J4=O(z^4), and the TOTAL background trace is O(z^4).
The cancellation of its constant term is exactly the vacuum identity.
This proves the stated outer gain; the symbolic calculation checks every
coefficient of (1), not a selected packet or a finite spectral matrix.

## 2. Keep the entire adjoint boundary moment

The physical radial measure is a constant times z^-4 dz. Hence D_z^*=3-D_z.
Let R_b xi=b_r and R_J xi=J. Their exact transposes on the original core are

 R_b* f(z)=Q'(z) int_0^z f(t)/Q(t) dt,
 R_J* f(z)=z^4 R_f(z)/Q(z)+Q'(z) int_0^z t^4 R_f(t)/Q(t)^2 dt,
 R_f(z)=int_z^z_h t^-4 f(t) dt.                    (3)

These follow by Fubini from the original regular-bolt primitive. Endpoint
terms vanish on the actual phase/core; the identities extend by the
weighted integral estimates, not by smoothing the reflected cubic.
There is an additional NECESSARY bolt cancellation: the two adjoint
terms in (3) separately can have logarithms. The complete seed coefficients
in (2) are polynomials, and the exact calculation checks

    c_J-Q'c_b=Q d_b, with d_b a polynomial.          (3a)

For f=c_b v, d=d_b v, put
I(z)=int_z^z_h Q(t)[d(t)-f'(t)+4f(t)/t]/t^4 dt.
Integration by parts before division gives exactly

 R_f_for_J=-Q f/z^4+I,
 R_b*(c_b v)+R_J*(c_J v)
   =-f+z^4 I/Q+Q' int_0^z t^4 I(t)/Q(t)^2 dt.      (3b)

At the simple actual bolt, I=O((z_h-z)^2), so both divisions in (3b)
are regular. No independent logarithmic bolt term survives. This uses
the original regular branch, not a prescribed zero moment. Smooth-core
regularity extends by the radial weighted integral estimates; one areal
derivative of v costs at most two regular-disk derivatives here. The
operator therefore has the required order at most two on that chart.

For v with scalar leading z^4 and a faster remainder, c_J=O(z^3) makes

 mu(v)=rho_a^-1 int_0^z_h t^-4 c_J(t)v(t) dt

convergent at the outer end. Interior and bolt finiteness uses (3a)-(3b)
and the original regular-disk class; the separated terms are not assigned
independent logarithmic bounds. The entire finite-domain integral is
retained. It is a finite continuous moment on the stated class, with a
numerical constant still UNKNOWN.

Equations (2)-(3) show that b* v has an analytic z^4 moment and a faster
remainder. Locally its moment profile is

 mu(v)[z^4/Q+Q' int_0^z t^4/Q^2 dt],               (4)

with the entire finite-domain contribution retained. There is NO
z^4 log z term. The coefficient is not zero and is not discarded. A
partition changes the two contributions to the moment, not their sum.

For clarity, this is not the invalid inference from a small form norm to
a boundary trace. Critical form data have two local derivatives (including
the retained log reserve) at weight 3/2. Equation (2) sends them to scalar
L2 data at weight at least 9/2. The actual positive scalar inverse crosses
the root four, retaining v4 z^4, and gives four derivatives. Its heavy root
is above 18000. Applying (3) and the local differential adjoint gives a
TT source with two derivatives, the moment (4), and a remainder of weight
at least seven. The norms of these maps are finite at this fixed seed;
none has been numerically bounded in this note.

## 3. The inverse York map preserves absence of a leading log

Let A_Y=(2/(3z^2))(partial_z-4/z) and
L_Y=-partial_z^2+6/z^2+2(Q''-4Q'/z)/(3Q), with y(0)=y(z_h)=0.
The retained physical inverse is

 y=L_Y^-1 A_Y T,
 V^-1 T=-T-z y-z^2 y'/2.                           (5)

Integration in dz gives A_Y*=-(2/(3z^2))(partial_z+2/z), so

 V^-1=-I+(3/4)z^4 A_Y* L_Y^-1 A_Y.                (6)

Consequently V^-1 is self-adjoint in the physical z^-4 dz measure.
This is an identity of the original auxiliary realization, not an
assumed symmetry of a changed boundary problem. The positive potential
L_Y>=-partial_z^2+5/z^2 and the regular-bolt endpoint justify its inverse.

For T=z^4 phi with analytic phi, A_Y T=(2/3)z^2 phi'. In particular a
pure z^4 coefficient has ZERO resonant forcing in (5). The outer roots
of L_Y are 3 and -2; a leading z^4 log T would create the resonant z
forcing, but (4) does not. Thus (V*)^-1 b*v still has a pure z^4 moment,
no leading quartic logarithm, and a faster remainder. Its moment is
generally nonzero and retains its actual finite-domain value.

## 4. Both original equations give the required homogeneous phase

Define Q_T=(ell+B)^-1, the TT diagonal of the ACTUAL coupled inverse, and
F=Q_T B. Solve H_Y(s,xi)=(0,Bw): s=a^-1 b xi, xi=F w.
Its canonical source is (0,(V*)^-1 Bw). Section 3 implies that its
original normalized second row has a constant leading coefficient and
no logarithm. The first row is zero. In the exact localization/moment
identity of current-global-boundary-moment-return.md this gives

 g0=-(c_log+2 f_nu,0/3)=0.                          (7)

This uses BOTH undifferentiated equations; the scalar quartic is still
determined by its original finite-part row. No boundary condition or
source coefficient is assigned zero merely for convenience. Equation
(7) is a property of this error source, not of the physical affine source.

The retained fixed-seed coupled graph argument therefore gives mass cubic,
scalar quartic, the EXACT slow root 4+e_*, and weight-4.005 remainder,
with no independently generated forced z^4 shape. Full angular/thermal
terms, common contacts and finite-domain returns remain in Q_T. Its
normal root/moment theorem and finite-order regularity on physical halves
are used in their fixed-seed scope; no global smooth-double theorem or
nonlinear-family extension is invoked.

For form inputs these are WEIGHTED Sobolev expansions. The leading moment
is analytic, while the faster remainder has its stated finite derivatives;
it is not asserted to be pointwise conormal to all orders. The existing
weighted kernel estimates and Mellin contour argument apply with that
regularity: the faster L2 weight makes the residue integrals continuous,
and each subsequent application supplies four more derivatives. Thus no
analyticity of an arbitrary input, beyond the fixed seed coefficients,
is being inserted into the phase argument.

Each F gains four local derivatives: b has order two, a^-1 gains four,
b* has order two, and the actual coupled inverse gains four, retaining
the log and the just identified phases. Thus F^2 maps the form space
into the original fast graph and finite phase domains with ten local
derivatives. In particular the original graph and phase outputs of
V F^2 are bounded linear maps, with finite but UNKNOWN numerical norms.
The original mass-lift conversion, including its cost below, is included
when defining those outputs. This is not a claim about arbitrary phase
inputs to the nonlinear quantum variation.

For the scalar reconstruction of a Schur error, add (C*)^-1 g to the
canonical scalar component of the coupled solution just used. The actual
mixed inverse is a fourth-order ODE with roots six and minus one (and the
heavy scalar roots). Its z^4 source is nonresonant and supplies the
retained scalar quartic; it preserves the stated derivative count.
Five physical areal jets require a separate higher composite. Their
numerical X-norm comparison is NOT established here.

No claim of bounded graph output for P alone is restored: its earlier
actual-domain counterexample remains valid.

## 5. Complete composite identity and quantitative target

In the ell form coordinates, ell^-1 B is positive with norm <=q. Spectral
calculus gives ||F||<=q/(1+q), and exact noncommuting identities give

 P E^2=V F^2(I-F)^-3 Q_T J_F,
 ||Q_T J_F||Y_F->V_ell<=1,
 ||(I-F)^-3||V_ell->V_ell<=(1+q)^3.                (8)

Section 4 proves that (8) is in the graph/phase domain, not merely on an
energy core. It supplies finiteness, not the numerical constant needed
by the requested certificate. The remaining first unsupported step is
now that NUMERICAL full strong inverse/composite estimate.

Keep the original E energy normalization, E(u)^2=||u||V^2/(|W| beta).
For the actual affine solution with E(u)<1.1e-10, the source transform is

 J_F y_aff=g_Y+b* a^-1 f_Y,
 ||J_F y_aff||V_ell*<=sqrt(1+q)||(f_Y,g_Y)||V*,
 ||(f_Y,g_Y)||V*<=sqrt(1+q)||u||V.

Hence U_aff<=(1+q)sqrt(|W| beta)1.1e-10. This is the COMPLETE coupled
affine source conversion, not a unit input assumption or a new physical
norm. With N_halves<=2 and gamma>5e-34 the actual parameter enclosure
bounds it explicitly. It does not bound independent feedback inputs.

Let C_g=||O_graph V F^2||E->H3, where E on xi means its original form
energy ell(xi)^(1/2)/sqrt(|W| beta). The original graph is unchanged.
Equations (8) and the N=32 error formula give the sufficient affine test

 C_g (1+q)^4 (1.1e-10) q^30/(1-q) <=5.34235e-12,    (9)

ONLY if signed centers leave that full conditional margin. The checker
evaluates its ceiling. This removes an unknown affine input conversion;
C_g itself remains UNKNOWN. Independent phase/feedback input norms and
all signed centers remain required. Finite C_g means eventual exact-
operator error decay, not that (9) holds at N=32 or that a trial center
fits admission. No actual inverse has been numerically applied.

## 6. Unchanged requirements and research decision

The mass lift still cancels only the uncut reference rows. Its .1/.2
transition and original regular bolt are retained. Original conversion
2133<||d_chi||H3<2174 and LINEAR debit <=3.538254568e-12 remain.
Actual post-elimination kappa_m is UNKNOWN. Conditional kappa_m<=5e-5
still spends 2.14765e-12 affine, .0024731 feedback and .0195239 energy.
All remaining entries must fit jointly:

 Y_other+4.084e-22 A_o+Y_phase<=5.34235e-12,
 eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase<.0012069.

The fast screen .135973<.1599164 retains only its two-fast/admitted-family
scope. Physical-X comparisons, E_center+1e-8 alpha<=1e-11,
R_X alpha*.028913+Z0_other<=.25,
R_X alpha^2*1.0706e10+Z2_other<=5e7,
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 remain unproved.
So do exact compatibility, physical matching, validity/noise and stability.

OUTCOME: unresolved quantitative coupled estimate. This addresses the
fixed-seed composite DOMAIN and its affine input conversion, not a passed
strong residual or actual solution. Supporting algebra must verify the
whole mixed operator and the adjoint identities; it is not a substitute
for the ordinary phase/regularity argument above.

DECISION: retain this method for the quantitative composite test, with
the boundary moment explicitly retained. Do not optimize weak residuals
or the controlled cold/scalar contributions. ONE NEXT CALCULATION:
bound the COMPLETE C_g and phase analogues in (8), including their full
common-restored quantum inverse, and the independent feedback inputs.
Compare actual products and signed centers with (9) and the joint budget.
That numerical calculation is PROPOSED, NOT EXECUTED.

Remaining to-do list: quantitative coupled composite, original source/X
comparisons and signed columns; nonlinear admission/contraction; exact
compatibility; physical matching; validity/noise; all-sector stability;
opening, transport, preparation, separate creation and stable wormholes.
