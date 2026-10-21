# Coupled thermal boundary compatibility after shape elimination

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full wormhole goal remains incomplete.

Removing BOTH fast curvature traces still leaves a source compatibility
condition after eliminating the complete auxiliary thermal shape. At each
fixed admitted Matsubara frequency, the resulting strong curvature operator
is Fredholm of index -1 on that minimal trace domain. Its cokernel has
dimension one greater than its kernel. No smallness of the feedback is
needed. This is a result for the retained seed and the specified auxiliary
realization, not a physical boundary selection or a Lorentzian mode theorem.

The distinction from the preceding [Schur form](thermal-curvature-shape-schur.md)
matters: compactness of its variational feedback alone does not imply
compactness of the forcing in the physical curvature equation. We derive
that strong equation, retaining the nonzero seed residual, before using
the [scalar range theorem](physical-curvature-source-and-boundary-range.md).

## 1. Strong trace row and density contacts

Keep the original action, HH state, magnetic charge, angular tower, clock,
subtraction and finite contacts. Fix gamma at a positive Matsubara frequency
above both retained shape and reconstruction thresholds. Use the preceding
notation H=L2(dy), U=r(delta R-r R_r p), w=r v, q=r p and

    q=R_gamma U+Q_gamma w,
    C_gamma=A1/6+(B1 partial_y+B0)R_gamma,
    Atilde=-partial_y^2+F(F_r/r-E+11Z/36-c_R).

R_gamma, Q_gamma and C_gamma are bounded on H as proved in
[physical reconstruction](physical-curvature-reconstruction.md) and
[the full shape inverse](global-thermal-shape-inverse.md). Put

    P=F(r calT_r+4calT), epsilon=G/(60pi),
    S=rF[delta calT+(4p+2v)calT].                         (1)

S is the conformal coordinate-density row with its overall action factor
removed. It is not just the gauge-invariant trace source. The exact
kinematic curvature identity and trace equation give

    S=-epsilon(Atilde+gamma^2+C_gamma)U
                          +P R_gamma U+D_gamma w.        (2)

Here D_gamma is exactly the mixed-row map of the preceding Schur proof:

    D_gamma w=epsilon g_w+P Q_gamma w+2F calT w,
    g_w=A1 U_v(v)/6-B1 partial_y(Q_gamma w)
                                  -B0 Q_gamma w-2A0 v.   (3)

To check (2) without assuming a field equation, the areal trace identity
reads rF delta calT=epsilon[-(Atilde+gamma^2)U+g_mix]
+F r calT_r q, while g_mix=-C_gamma U+g_w. Add the volume
contact F(4calT q+2calT w). The radial gauge contact and volume contact
combine into P; neither is set to zero on the seed.

For a curvature test V, the full congruent Hessian row is

    -(1/(2G))<R_gamma V,S>.

R_gamma^* is injective on H. Indeed, for compact test f,
R_gamma^* f=(1/6)F(K+gamma^2)^-1 f, with the same identity in
distributions for general f in H. F is strictly positive in the open
exterior. Thus R_gamma^* f=0 implies (K+gamma^2)^-1 f=0 and f=0.
Equivalently R_gamma has dense range. Consequently the zero congruent
curvature row is equivalent to S=0 whenever the strong S belongs to H.
This does not assert a bounded inverse of R_gamma^*, or identify arbitrary
form-dual sources with physical H sources.

## 2. Elimination in the strong source space

Let E_gamma be the full closed shape energy space, with inner product
a_gamma. Write K_cross=-(1/(2G))R_gamma^*D_gamma:E_gamma->H and let
K_cross^dagger:H->E_gamma be its form adjoint. The preceding result gives
bounded D_gamma:E_gamma->H and compact K_cross^dagger.

For an independent shape functional ell in E_gamma^*, let w_ell be its
Riesz representative. The complete congruent shape row gives

    w=w_ell-K_cross^dagger U.

Normalize the prescribed strong trace-density source by S=-epsilon f.
Substitution into (2) yields the exact reduced equation

    L_gamma U=f+epsilon^-1 D_gamma w_ell,
    L_gamma=Atilde+gamma^2+B_gamma,
    B_gamma=C_gamma-epsilon^-1 P R_gamma
                         +epsilon^-1 D_gamma K_cross^dagger. (4)

Thus all of B_gamma is bounded H->H. In particular

    ||B_gamma|| <= ||C_gamma||
           +||P||/(18epsilon lambda)
           +epsilon^-1 ||D_gamma|| ||K_cross^dagger||.   (5)

These are finite constants, not calibrated small error bounds. The strong
shape correction is D_gamma K_cross^dagger, not the variational operator
M_gamma=K_cross K_cross^dagger. It need not be self-adjoint or positive
in physical curvature H. Equation (4) is not a new positive energy claim.

Initial metric jets and independent state/boundary sources must first be
included in the original rows and transformed with the retained initial
reconstruction. They are not free curvature data. This note analyzes the
resulting fixed-frequency operator and its explicitly stated source space.

## 3. Relative compactness on the minimal curvature domain

Let D_min be the maximal scalar graph domain of Atilde with BOTH leading
outer coefficients A_c=A_s=0 and its limit-point condition at the horizon.
Use the graph norm ||U||+||Atilde U||. The preceding scalar analysis gives
the outer inverse-square potential and the exponentially integrable horizon
tail. Adding the bounded diagonal -F c_R does not change these traces.
C_c^infinity is a graph core for this minimal realization. This is a
candidate domain under examination, not an adopted gravitational boundary
condition.

The terms in B_gamma have the following compactness properties.

* (B1 partial_y+B0)R_gamma is compact H->H by the retained weighted
  derivative and local H2 estimates.
* P vanishes at both endpoints: calT=O(r^-4) and
  r calT_r+4calT=O(r^-5), so P=O(y^3) at zero; F vanishes at the
  horizon. Local H2 compactness of R_gamma and small multiplier tails
  therefore make P R_gamma compact H->H.
* D_gamma K_cross^dagger is compact H->H because the second factor is
  compact into E_gamma and the first is bounded into H.
* Multiplication by A1/6 is generally NOT compact H->H. It is compact
  D_min->H in the graph norm: A1 tends to zero at both endpoints, and
  an Atilde graph-bounded set is H2-bounded on every compact interior
  interval. Local Rellich compactness handles that interval; the small
  coefficient tails are controlled by ||U||.

Therefore B_gamma:D_min->H is compact in the graph norm. This is the
required relative compactness; no norm-smallness or unproved inverse is
used. Since B_gamma is also bounded on H, Atilde+B_gamma is a closed
operator on exactly D_min, and its graph norm is equivalent to the old one.
The core trace identity (2) extends continuously in that graph norm and
the shape energy norm. This defines the strong realization and the
associated row pairings; it does not assert that an indefinite curvature
form has a Friedrichs closure.

## 4. The missing source condition persists

At this real gamma>0, the retained scalar range result gives

    ker(Atilde|D_min+gamma^2)=0,
    Ran(Atilde|D_min+gamma^2)=ker Mtilde_gamma,
    ind(Atilde|D_min+gamma^2)=-1.                        (6)

Treat (6) as a bounded Fredholm map from the graph Hilbert space D_min
to H. Relative compactness and index invariance now imply

    L_gamma:D_min->H is Fredholm,
    ind L_gamma=-1,
    dim ker L_gamma^*=dim ker L_gamma+1 >=1.             (7)

For the standard theorem and its proof see
[Haller, Corollary I.3.3, page 12](https://www.mat.univie.ac.at/~stefan/files/ASIT/ASIT.pdf).
One can apply its parametrix proof directly to the graph-space map;
bounded perturbation leaves the domain fixed. The adjoint in (7) is the
unbounded H adjoint. Its kernel is the annihilator of the same closed
range, so its dimension is the graph-map cokernel dimension.

For any basis z_1,...,z_m of that kernel, solvability of (4) is EXACTLY

    <z_j,f+epsilon^-1 D_gamma w_ell>=0, j=1,...,m,
    m=dim ker L_gamma+1.                                (8)

These are Hermitian H pairings. The z_j are the full coupled adjoint
solutions, not the old unperturbed horizon function. Their profiles and
the actual kernel dimension have not been evaluated. In particular we
do not claim exactly one condition, uniqueness, or invertibility on the
compatible subspace without accounting for ker L_gamma.

The full pair of strong trace and congruent shape rows has the same
index -1 as a map D_min x E_gamma -> H x E_gamma^*. Eliminate shape
with the bounded triangular substitution above and subtract its bounded
source contribution; the diagonal blocks become L_gamma and the Riesz
isomorphism. This proves the statement for the complete auxiliary pair,
not just its isolated curvature row.

Since its range is a proper closed subspace, there are arbitrarily small
incompatible H trace sources for every fixed shape source. It has no
right inverse on all H x E_gamma^*. Therefore removal of both fast traces
does not by itself produce a freely forced auxiliary evolution, even
after the full thermal shape response is included. This is stronger than
the preceding fixed-shape compatibility test.

## 5. Physical decision and verification

The result does not decide which full trace/state sources the gravitational
constraints permit. Those sources could lie in (8); their compatibility
must be demonstrated, not assumed. Nor does this thermal operator theorem
give continuous causal-frequency control or a growing physical mode.
One may not infer a no-go theorem for wormholes or for every physical
boundary prescription. The background is still the off-shell seed.

The next physical step is to derive the admitted boundary traces and
source class from the retained gravitational preparation, then evaluate
the coupled compatibility or analyze a justified one-condition boundary
realization. Changing the curvature boundary domain changes this index
problem and requires an explicit physical justification. No boundary
phase, contour rotation, counterterm or order reduction is chosen here.

Focused verification checks the strong trace/contact bridge against
independent four-dimensional curvature contractions at charged metric
jets, all row/source signs and endpoint coefficient limits. Exact
finite-rank perturbations of an infinite unilateral shift check the
Schur source transformation and the distinction between index -1,
injectivity and cokernel dimension. Their unchanged shift tail proves
their infinite-dimensional index; they are controls, not approximations
to the charged Dirac kernel. The charged Fredholm theorem remains an
ordinary proof, not a consequence of a finite matrix rank test. No Lean
campaign is selected. Freeze the note, checker and review before the run;
retain every outcome and do not edit sources after successful acceptance.

Remaining to-do list: physical boundary/source domain and coupled
compatibility; continuous causal response and constrained stability;
nonlinear existence and signed compatibility; exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and full audit.
