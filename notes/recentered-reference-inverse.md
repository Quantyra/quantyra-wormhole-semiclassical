# The recentered bulk reference inverse and its mass-boundary limitation

2026-10-22, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress, satellite 4ff6739 and planning d9aaf18.
The full stable semiclassical wormhole objective remains incomplete.

At the same c_R=10^22 RN reference, the actual compact-test bulk
Hessian has a closed coupled realization and an inverse on the
mass-extended weighted space. Its scalar range has no residual phase
condition of the original zero-c construction. This removes a linear
bulk-domain uncertainty; it does NOT identify that closure with the
complete physical renormalized variational problem.

The distinction is consequential. The mass coefficient is unbounded
in the critical graph and energy norms. Any surviving boundary mass
momentum must vanish exactly, or be included as a separate boundary
equation. It cannot be absorbed into the proposed small transfer defect.
The next decisive physical calculation is therefore that boundary
momentum and its variation, with the actual paired determinant and local
boundary action, before finite-geometry/Newton estimates are extended.

## 1. Actual operator and fixed spaces

Keep the reference, G, flux, prescribed positive temperature, paired
Dirac realization, fixed scale mu0 and all contacts of
[the mixed-form proof](recentered-mixed-form.md). No field equation is
imposed on this off-shell reference. With fixed physical volume dmu,
the York/outer-fixed disk return gives the compact-test action Hessian

    H = [-A_s, B; B*, D_T],
    <s,A_s s>=a(s,s), <xi,D_T xi>=l(xi,xi),
    A_s=A/(8piG), A=6(-Delta+4lambda)(-D Delta+C_v)
                              -(22epsilon/3)div(J grad),
    D=96piG c_R-epsilon, epsilon=G/(60pi),
    C_v=1+11epsilon lambda/3, F_lambda=1-384pi c_R G lambda.

The mixed form is b(s,xi)=<s,B xi>. In particular the scalar action
block is NEGATIVE; reversing only that equation reverses one mixed sign.
The nonlocal determinant remains in D_T. The density return
-div(T0 eta), every finite contact and the physical 8piG are retained.

Let z be the fixed disk end coordinate proportional to 1/r. Use exactly
the mixed spaces and physical transpose of
[the weighted TT theorem](stationary-weighted-tt-inverse.md):

    Z_delta=z^delta H_b^4, W_delta=z^delta H_b^(4,log),
    Y_delta=z^delta L_b^2, dmu=gamma(z)z^-3 dz/z,
    gamma(z)=gamma0(1+O(z^2)), gamma0>0.

At the bolt these are smooth rotational disk norms, with no artificial
inner boundary. The logarithmic multiplier in W is
<D_logz>^4 log(2+|D_logz|). The critical weight is 3/2. The scalar
sphere-component TT pairing, not the tensor Gram matrix, is used.
The genuine TT reconstruction and its outer-fixed return are unchanged
formulas with the NEW f: (r^2 f b)'=-(r^2 f'-2rf)xi. Their inverse
York map is bounded on -1<delta<4. Its homogeneous tails have weights
-1 and 4; neither is a free physical datum. Positivity of the vector
inverse uses the same d=0 RN guards already verified in the mixed proof.

## 2. Rechecking the actual parametrix hypotheses

This step uses the finite kernel construction in sections 2-4 of
[the global weighted map](stationary-global-weighted-map.md), not an
assumed compactness of a full-order coefficient perturbation. The
following checks identify its application here.

* The new RN compactified metric is analytic on each exterior side,
  has fixed constant/quadratic jets, a mass cubic, and higher powers.
  The even reflected seam can have a cubic loss of smoothness; it is
  treated by the actual half-space boundary problem, not an all-orders
  smooth double. The canonical boundary is geodesic. The finite five
  conormal jets required by the kernel proof exist with bounded values.
* The monopole connection, spin structure and paired complementary
  first-order domains are unchanged. The direct proof of
  ||(Dsharp+M_n)^-1||<=1/M_n and M_n>18 sqrt(n) is already checked
  on this new geometry in the quantum TT note. The LLL uses its
  regular HH disk allocation; no new massless gap is postulated.
* The reflected Green construction consequently has the same local
  orders: D-D_flat has derivative coefficients O(R^2) and zeroth
  coefficients O(R); its first Green correction has degree -1 versus
  -3. On comparable annuli of size h the combined finite-part Hessian
  error is O(h^2 log^N h) from H4/H4log to L2. For separated annuli
  the coordinate-density difference is O((z+y)^-3 log^N(z+y)),
  versus the leading (z+y)^-5 radial kernel. These statements include
  the inverse-square contact's O(log^N z) remainder, its differentiated
  rows, and common Weyl restoration. The changed mu0 B is retained.
* After weighting, separated kernels decay with rates 4-delta and
  1+delta. The end error is at most C h^2(1+|log h|)^N on every
  closed substrip. C is finite, not numerically estimated here. A
  sufficiently short end exists; no useful finite cutoff is claimed.
* The extra local -c_R(R-R0)^2 action has zero background gradient.
  In TT variables its Hessian is the multiplication form
  -32c_R t^2 xi^2, t=Z/r^4, with no new fourth-order TT term.
  The effective Einstein coefficient is F_lambda. Mixed terms have
  at most two derivatives and vanishing outer coefficients. The
  scalar fourth-order coefficient is -6D/(8piG), nonzero. Thus the
  exact compact-interior/disk mixed parametrix remains elliptic.

The cited finite construction extracts the direct logarithmic diagonal
before estimating the rest. Taking the absolute value of that finite
part, or estimating the boundary heat contacts by a separated kernel,
would not establish these bounds. On the present analytic reference
the needed finite differentiations and adjoint regularity have no C5
bootstrap ambiguity. State-dependent global Green coefficients are
retained in the degree-two/three expansion, not set to zero.

[Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf) supplies
the normalized-symbol Sobolev norm estimate (output order zero, positive
spatial Holder regularity, two frequency derivatives in one dimension).
[Gaarde--Grubb, Theorem 2.4 and Proposition 2.9](https://arxiv.org/pdf/math/0703878)
describe the differential boundary resolvent and the nonstandard
logarithmic reflected term. They do not by themselves prove a theorem
for this determinant Hessian; the retained finite expansion supplies
that identification. The large c_R changes finite constants, not these
orders or the strictly positive D hypothesis.

Use the actual interior parametrix and the end normal inverse below.
Invert the tail-small error by a Neumann series. Cutoff commutators
gain order and have compact support or decaying separated tails; these
remainders are compact. Full-order errors on a fixed annulus are NOT
called compact. This gives left/right parametrices and distributional
adjoint regularity on noncritical weights. The same proof applies at
3/2 because the NEW scalar normal family has no roots there.

## 3. New normal families and the complete relevant root count

Put q(s)=s(s-3), d0=1/(480pi^2), ell=lambda^-1/2, and

    F_M(s)=psi(s-1)+(pi/2)cot(pi s/2),
    kappa_M=-6/5-log(mu0 ell),
    A_T=F_lambda/(48piG lambda)+23/(8640pi^2),
    N_T(s)=d0 q(s)(q(s)+2)(F_M(s)+kappa_M)-A_T q(s),
    J_R(s)=3[C_v-Dlambda q(s)][q(s)-4]/(4piG lambda).

The scalar-coordinate TT normal family is 9lambda^2 N_T, and the
scalar one is lambda^2 J_R. Mixing has zero normal family. In disk
coordinates the reflected congruence is still
C(3-s)^T diag(J_R,N_T)C(s),
C(s)=[[1,(s+3)/(3(s+1))],[0,2]]. The transformation poles at -1,4
remain outside [-1/2,7/2]. The added square contributes no TT normal
term because J tends to zero; the effective Einstein term gives the
displayed F_lambda factor. None of the finite quantum constants is reset.

The new scalar roots are

    -1, 4, Delta_-=(3-sqrt(9+4m^2/lambda))/2,
           Delta_+=(3+sqrt(9+4m^2/lambda))/2,
    m^2/lambda=C_v/(Dlambda),
    3.39e8 < m^2/lambda < 3.40e8.

In particular Delta_-<-18000 and Delta_+>18000. There are NO scalar
roots in -1/2<=Re s<=7/2. On the critical line q=-t^2-9/4, so J_R
is strictly negative. The roots -1,4 and the fast scalar are not
discarded as unphysical: the chosen fast remainder includes decaying
scalar tails, while the growing roots are outside critical L2.

The retained fixed scale still gives 9<-log(mu0 ell)<12, hence
39/5<kappa_M<54/5. Also F_lambda>1-1.179e-8 implies

    A_T/d0=10pi F_lambda/(G lambda)+23/18>30*2000^10.

Thus the entire-strip TT root argument applies with the new A_T.
For Re s in [3/2,7/2], |Im s|<=2, its removable product obeys
|(q+2)(F_M+kappa_M)|<400. For |Im s|>=2, the digamma series and
cotangent identity give Re[(q+2)(F_M+kappa_M)]<12. Reflection
F_M(3-s)=F_M(s) treats the other half. The only TT zeros in the
strip are the simple roots 0 and 3. At 1 and 2 the product q(q+2)F_M
has limit -2; those points are not spurious roots. There are no poles
in this strip after those removals. The log symbol is elliptic at
large imaginary part. This is a continuum root exclusion, not sampling.

## 4. Actual graph closure, then the coupled inverse

First realize D_T on W_(3/2) in physical scalar L2. The parametrix
proves closedness and upgrades its maximal distributional adjoint
domain to W_(3/2). Symmetry on the compact core and the fixed measure
then identify the adjoint domain, so this is self-adjoint. These are
the actual combined Hessian and its graph domain; no separate positive
model form is substituted for their sum.

To extend l>=X^2 from the proved core, first truncate/smooth xi in
W_(3/2). A compact sphere component may leave a radial TT r^-4 tail.
Add a bump a_L psi(r/L) in [2L,3L] to cancel
int_(r+)^infinity (r^2f'-2rf)xi dr. Since r^2f'-2rf=-2r+6M-4Z/r,
a_L=O(L^-2). Its W_(3/2) norm is O(L^-1/2), including the logarithm:
after multiplication by z^-3/2 it is a fixed smooth log translate with
that amplitude. The analytic coordinate change has uniformly bounded
jets. A diagonal sequence gives a compact genuine TT graph core.
This is NOT a density assertion at weight 7/2.

The pairing <xi,D_T xi> and X^2 are graph-continuous, hence
l>=X^2 on all W_(3/2). In particular the actual D_T is positive with
L2 gap at least .07 F_lambda lambda/(8piG). Its spectral form is
closed in its OWN norm, and is the closure of the compact-test form:
an operator core is a form core for a positive self-adjoint operator.
This resolves the previous reference-TT closedness uncertainty without
identifying its completion with the stronger N space or with the full
physical boundary domain. Its inverse satisfies ||D_T^-1||_(X*->X)<=1.

The positive scalar form a is similarly the closed form of A_s on
Z_(3/2); its local elliptic estimate and scalar root exclusion give
the same graph/maximal-domain identification. Let V be the direct sum
of these two ACTUAL reference form domains, norm

    ||(s,xi)||_V^2=a(s,s)+l(xi,xi).

The mixed bound |b|^2<.336 a l extends by core density to V. Let
J_sign(s,xi)=(-s,xi). In the real pairing (or its real Hermitian part),

    Re <J_sign u,H u>=||u||_V^2.                         (1)

The two mixed terms cancel with their actual opposite signs. Thus
Lax--Milgram for J_sign H gives a unique solution for every V* source,

    ||H^-1 f||_V <=||f||_V*.                            (2)

Equivalently in the a,l Riesz coordinates H is [-1,K;K*,1],
||K||^2<.336, and H^2=diag(1+KK*,1+K*K). This argument does not
mistake Euclidean saddle coercivity for Lorentzian stability. For L2
sources the maximal-domain parametrix upgrades the form solution to
Z_(3/2) x W_(3/2), identifying it with the actual graph inverse.

This is sharper than the earlier absolute Schur estimate <2.42.
It uses a favorable sign already present in the action, not a new norm
normalization or discarded coupling. If a future transferred operator
in THESE same form spaces differs by V->V* norm at most 1/4, its
inverse is bounded by 4/3. That transfer constant remains UNKNOWN;
neither number is an admission-norm inverse bound.

## 5. Fast sources, one mass datum, and no bulk scalar obstruction

For f in Y_(7/2)^2, apply the critical inverse. Iterate the actual
minus-normal two-power improvement in steps of length one. From weight
3/2 to 5/2 no root is crossed; from 5/2 to 7/2 only the simple TT
mass root 3 is crossed. The scalar has no root on either step. Its
residue therefore gives exactly

    s in Z_(7/2), xi=m_c chi(z)z^3+v, v in W_(7/2).    (3)

The finite degree-two/three kernel/contact expansion from section 2
retains global coefficients. Corrections to this mass mode start at
degree five locally; remote York/kernel terms are O(z^4 log^N z),
which is already beyond 7/2. Thus no omitted forced term lies on the
new final weight. The cutoff mass column has output in Y_(7/2)^2
by the same local/far estimates. Its disk image is (3/4)z^3(1,-2)
modulo the fast space, with the single original correlated mass datum.

Consequently H:Z_(7/2) x [W_(7/2)+span{chi z^3}]->Y_(7/2)^2 is
a bounded bijection: onto by the critical construction/expansion, and
injective because this domain embeds in the critical domain. The
bounded inverse theorem supplies a finite weighted inverse constant,
but no numerical value for it has been proved. Fast-only index is -1;
the mass extension has index zero. There is no old scalar adjoint-phase
condition in this reference bulk realization. This is NOT permission
to drop any boundary, current or finite-geometry compatibility equation.

The equivalent contour index change counts two roots, versus four in
the zero-c construction. [Melrose's relative-index theorem, section
6.2](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/melrose.pdf) motivates
that bookkeeping for b-differential operators. The logarithmic response
here uses the finite-strip parametrix and residue construction just
stated, not an unsupported assertion of membership in that calculus.

## 6. Why the physical mass boundary calculation cannot be bypassed

Let xi_L(r)=r^-3 chi(r/L), with chi=0 below 1 and chi=1 above 2.
It is smooth at the bolt, has r^-3 coefficient exactly one (a fixed
nonzero multiple of m_c),
and belongs to the mass-extended domain. The critical W norm obeys

    ||xi_L||_W(3/2)=O(L^-3/2),
    ||xi_L||_L2=O(L^-3/2), l(xi_L,xi_L)=O(L^-3).       (4)

Indeed r^(3/2)xi_L=L^-3/2 p(log r-log L), where
p(v)=exp(-3v/2)chi(exp v). Translation preserves every indicated
Sobolev/log multiplier norm. The actual fixed-coordinate/measure
changes are bounded and tend to constants; the actual operator is
bounded W_(3/2)->L2. The York tail does not change that scaling.
Thus the coefficient m_c is NOT continuous even in the critical graph
norm. The simple model chi(r/L)=(1-L/r)^5 for r>=L (zero below)
already gives int r^2 xi_L^2 dr=1/(858 L^3); smoothing the junction
preserves the scaling and the mass coefficient. No discontinuous mass
trace is silently supplied by an energy inverse.

If the renormalized physical first variation leaves sigma_boundary
delta m_c, its norm on these unit-mass fields divided by their energy
norm grows at least as L^(3/2) when sigma_boundary is nonzero. If the
boundary Hessian includes gamma_boundary delta m_c delta m_c, its
relative form ratio grows at least as L^3 for nonzero gamma_boundary.
Therefore there is no positive numerical allowance such as 1/4 for
such a term on V. The target is EXACT cancellation/vanishing, or an
explicit augmented boundary equation in a norm controlling m_c. A
small coefficient alone does not close this requirement. General
boundary couplings to m_c likewise require their own trace analysis.

Compact tests annihilate such boundary functionals and cannot detect
them. Nor does a formal completion prove they are absent from the
physical renormalized action. The existing scalar and mixed limiting
boundary checks leave the independent shape/mass momentum unresolved.
The next calculation must differentiate the COMPLETE paired determinant,
local EH+c_R boundary action, WZ boundary terms and fixed common
counterterms along the correlated mass/fast family, at fixed leading
clock, temperature, scalar source and flux. It must compare that first
variation and Hessian with the bulk graph operator, including the
cut-removal limit. No arbitrary finite-cut Navier wall or new boundary
coefficient is introduced to enforce cancellation.

For exact nonlinear static solutions, the retained two disk equations
plus the joint current/Ward identity give
(a_orb^2 r^2 E_rho)'=0; regularity at the bolt kills the constant.
Off shell its linearization includes connection-times-residual terms.
The reference inverse is for the density-correct reduced rows and
does not separately solve Maxwell/current or remove those contacts.
Lorentzian constraint initial and incoming data remain independent.

## 7. Research decision and verification

Continue this specified c_R construction. The actual reference bulk
domain and mass-extended range now have a justified inverse, and the
old scalar phase obstruction is absent in that scope. Prioritize the
renormalized mass boundary momentum: the demonstrated discontinuity
makes it an earlier gate than numerical finite-geometry transfer.
If a nonzero physical boundary equation remains, solve its augmented
compatibility or document a construction change; do not set it to zero.

The reference bound (2) does not estimate the complete residual, cold
admission or a nonlinear correction. At rho=1e-8 the targets remain
Y0<=5e-9, Z0<=.25, Z2<=5e7, with every positive admission margin.
The original C_I=9.1591500000000000015e103 failed upper certificate
remains parked and is not a nonexistence result. Finite geometry/state,
exact compatibility, fixed-flux/temperature matching, semiclassical
omission/noise control and all-sector retarded stability remain open.

Freeze this ordinary argument, checker and source review before running.
Check the new normal polynomials, root guards, reflected congruence,
graph-core and mass-trace scaling, and the saddle identity exactly.
Independently vary the unexpanded AdS conformal trace density to check
the scalar normal symbol, and the unexpanded anisotropic local action
to check the F_lambda change of the TT second-order term. Use 70/100
digits for the actual benchmark roots and Mellin identities. Preserve
all runs, failures, hashes and raw logs. No sampled spectrum proves
the root exclusion, graph closure or physical boundary cancellation;
these continuum arguments retain their stated scope. No Lean selected.

Remaining to-do list: physical renormalized mass/shape boundary momentum
and complete current/constraint realization; finite state, residual,
admission and nonlinear solution; exact compatibility, matching,
validity, all-sector stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
