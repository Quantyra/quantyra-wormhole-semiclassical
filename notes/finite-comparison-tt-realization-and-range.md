# Actual finite TT realization and coupled physical range

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof before
frozen verification. **The finite comparison's TT inverse and coupled
one-condition range are established. The nonlinear goal remains open.**

At each original finite comparison g_ell, the complete returned TT
operator has its actual positive selfadjoint realization. Its scalar-c
critical inverse obeys

    ||D_ell^-1||_(L2->L2)<=8piG/(.90lambda),
    sqrt(E_ell[D_ell^-1 f])<=(1/.90)||f||_(E_ell-dual). (1)

It has an isomorphism on the physical fast domain plus exactly one
correlated mass column. Applying the previously conditional mixed
estimate now gives an actual outgoing auxiliary contraction below .025.
The full fixed-density physical mass-extended Hessian has

    ker A_ell,m=0, index A_ell,m=-1, dim coker A_ell,m=1. (2)

Its range is the kernel of the actual adjoint moment, whose value for
the preparation source remains unevaluated. Neither a useful fast
fourth/log inverse constant nor a nonlinear solution is supplied.

The operator domain needs care: its logarithm belongs to the physical
disk SHAPE component. At C5 regularity one must not silently replace
that requirement by a logarithmic multiplier on the scalar TT amplitude.

## 1. Keep the physical disk domain and change its TT coordinate

Use the full action, fixed density, physical measure and mixed disk
space of the [actual weighted map](stationary-global-weighted-map.md):

    X_delta: theta=S+T/3 in Z_delta=z^delta H_b^4,
             T in W_delta=z^delta H_b^(4,log),
    Y_delta=z^delta L_b^2, dmu=gamma(z)z^-3 dz/z,
    -1<delta<4.

At the bolt use the regular two-dimensional disk norms. At delta=3/2,
Y_delta is physical scalar L2. The physical fast weight is 7/2 and its
mass vector is e_m=chi(z)z^3(1,-2). All these spaces and boundary data
are unchanged. Each finite comparison is smooth on the open exterior
and at the bolt, with the retained five conormal outer jets; uniform
higher boundary derivatives are not assumed.

At this fixed metric let A=a'/a, R=r'/r, w=A+2R and

    P Y=3(Y''+wY')/2-VY,
    V=w'/2+2A^2+4R^2,
    P Y=T'+4RT.                                      (3)

This is the actual York extraction, with its regular vector bolt and
decaying outer realization. The [finite York gap](finite-comparison-york-gap.md)
and global Green construction give its weighted inverse for -1<delta<4.
For integer 0<=j<=4 the resulting vector is H^(j+1) when T is H^j.
This uses the established ordinary C5 mapping, not its unproved
fourth-log upgrade in the scalar-c coordinate.

Define

    q=Y'-2T/3,
    K T=(q+wY)/4,
    N T=-2T/3+2K T-2RY.                              (4)

The extracted conformal amplitude is tau=theta-K T, and the sphere
component of the genuine TT tensor is c=N T. The returned TT disk
column and complete coordinate map are therefore

    R_T T=(K T-T/3,T),
    U(tau,T)=(tau+K T-T/3,T).                        (5)

The inverse to N is the actual scalar TT-to-shape return Lc=T_c.
Indeed York extraction and the outer-fixed return are inverse: their
possible disk-preserving discrepancy solves eta'=Aeta, hence is a
multiple of a~z^-1, excluded by delta>-1. This also preserves the
regular bolt and the fixed radial z^4 tail. Thus N and L are mutually
inverse bounded maps on weighted L2 and on ordinary H^j, j<=4.
Use the scalar c pairing for N*, L* and D; no tensor Gram matrix is
substituted for these physical transposes.

## 2. The cancellation that avoids a sixth derivative

Equation (3) implies the exact identity

    (K T)'=(2/3)RT+(1/3)(w'+A^2+2R^2)Y
           =(2/3)RT-(1/3)Ric_(rho rho)Y.             (6)

The derivative Y' and the dangerous wq term cancel. The intrinsic
radial form is d(KT)=(2/3)T dlog r-(1/3)Ric(Y, .); it extends
through the smooth disk bolt and involves no bound on A separately.
The Ricci coefficient needs three more derivatives for an H3 product,
hence at most five metric derivatives in total. Ordinary York bounds,
Leibniz and (6), with the weighted L2 bound fixing the constant, give

    K:z^delta H_b^3 -> Z_delta,
    N+(2/3)I:z^delta H_b^3 -> z^delta H_b^4.         (7)

For the second map use (4), with Y in H4 and dlog r in C4. Its
inverse has the analogous gain: the [finite return identity](finite-comparison-disk-return-form.md)

    (T_c+3c/2)'=2Rb-E1 eta,
    E1=a''/a-(R'+AR)

uses the same five metric derivatives. Lower integer orders and the
physical output transposes follow in the same strip by duality.
Constants are finite for the given comparisons; no numerical H4/log
constant is claimed here.

Since W_delta embeds in z^delta H_b^3, (7) proves directly

    U: Z_delta x W_delta -> X_delta is an isomorphism. (8)

This triangular map never asks K T or c to have H4log regularity.
The scalar-c TT domain is the explicitly identified pullback

    W_delta^c=N(W_delta)
      ={c in z^delta H_b^4: Lc in W_delta},
    ||c||_(W_delta^c)=||Lc||_(W_delta).              (9)

It is exactly the original physical disk domain, expressed in TT
coordinates, rather than a weaker substitute. For example the local
sequence v_n=n^-5 sin(nx) has bounded C5 norm, but v_n' has H4log
norm growing like log n. A C5 bound alone therefore cannot justify
an order-zero logarithmic multiplier involving a first metric
derivative. This is a limitation of that inference, not a statement
that the actual finite metric is nonsmooth or fails a sharper norm.
The analytic seed's simpler W_delta^c=W_delta statement is unaffected.

## 3. The TT shape operator and its actual normal family

Define, using the existing complete fixed-density action operator,

    D_T=R_T* A_disk R_T: W_delta -> Y_delta,
    D_c=L* D_T L: W_delta^c -> Y_delta.              (10)

Boundedness follows from (8), A_disk:X_delta->Y_delta^2, and the
weighted L2 return transposes. This keeps every local contact and
the actual nonlocal paired Dirac response. It does not manufacture
a different Friedrichs operator from the form lower bound.

At the fixed outer normal metric the exact York solution is

    Y(s)=-2T/[3sqrt(lambda)(s+1)],
    K(s)=-2/[3(s+1)], N(s)=-2/3,
    R_T(s)=(-(s+3)/(3(s+1)),1)^T.                    (11)

With the retained disk congruence matrix
C(s)=[[1,(s+3)/(3(s+1))],[0,2]], C(s)R_T(s)=(0,2).
Thus, with the physical reflected transpose at 3-s,

    N_(D_T)(s)=4lambda^2 N_T(s),
    N_(D_c)(s)=9lambda^2 N_T(s), N_mixed(s)=0.       (12)

The complete TT normal N_T, its only roots 0 and 3 in the relevant
strip, and its positive critical line are those of the weighted-map
theorem. The fixed-theta diagonal's critical zero is irrelevant:
R_T is the actual returned TT column, not the direction (-1/3,1).

In the compact interior and regular disk, D_T has the shape principal
symbol p^4 log|p|/(120pi^2). To justify the parametrix with finite
regularity, work in (theta,T) coordinates. Terms introduced by K
are locally compact W_delta->L2: an input-side K factors through
H3->H4 and the compact embedding H4log->H3; an output-side K*
is locally compact on L2 since K:L2->H1. Use separated cutoffs
to handle their global tails. The exact interior shape parametrix
from the weighted-map proof therefore applies. Full-order errors
in the actual shape diagonal are treated by that parametrix, not
declared compact merely because they have variable coefficients.

At infinity (6) and the York Green kernel give the same two-power
actual-minus-normal improvement in the H3->H4 K norm and L2 output
transpose. The coefficients require only the retained five conormal
jets. The global z^4 York term remains. Composition with the actual
disk kernel preserves its integrable separated tails and tail-small
O(z^2 log^N z) error. Use the inverse of the FULL scalar normal
4lambda^2 N_T on a short end; do not use the fixed-theta normal.
The Neumann and cutoff construction therefore supplies left and right
parametrices for D_T whenever that scalar normal is invertible.

The corresponding adjoint regularity is part of this scalar
construction. Locally K is order minus one through the required
integer orders; the K corrections map H3 to L2. Starting from an
L2 distributional adjoint null vector, the elliptic shape parametrix
bootstraps through H1,H2,H3 and then H4log. The weighted normal inverse
gives the same statement at the complementary weight 3-delta.
No logarithmic multiplier on the C4 coefficient in Lc has been used.

In particular D_T is Fredholm at delta=3/2. Its symmetry and adjoint
regularity make its critical index zero. This statement is independent
of the coupled operator, which has scalar phase roots on that line.

## 4. Positivity on the correct critical graph domain

The [complete disk form](finite-comparison-quantum-disk-tt.md) gives
.90 E_ell on the regular compact-outer TT core and its explicit
stronger form completion N_ell. A general critical operator-domain
field need not lie in N_ell: the positive model and negative position
logarithm may separately diverge. The extension below uses two stages.

First let T be smooth, bolt-regular and compact at the outer end.
Its York vector has its permitted z^4 tail. Consequently c=N T,
b and S have z^4 conormal tails, while the returned shape is exactly
the chosen compact T. These inputs lie in N_ell. The nonlocal shape
arguments are compact, and the remaining conformal terms are local
with integrable z^4 tails. Their N_ell form agrees directly with
the actual weighted action pairing; no ambiguous subtraction of
two infinite forms occurs for these inputs.

Here is explicit approximation by the compact-TT core. Use areal r,
which is monotone on the actual finite path, and set

    U=a^2 r^2, M=r^2 partial_r(a^2)-2ra^2,
    b=-U^-1 integral_(r_h)^r M c dr.

Fixed clock and quadratic outer data imply
M=-2r/kappa^2+O(1). Cut c off beyond r=L and add a smooth bump
a_L psi(r/L) on [2L,3L], with integral M c_L=0. Choose a nonnegative
psi with integral s psi(s)ds>0. The denominator is a nonzero constant
times L^2+O(L), so a_L=O(L^-2). The resulting genuine TT tensor and
its outer-fixed return vanish beyond 3L, and remain regular at the bolt.

The difference has E_ell=O(L^-1). On the moving outer annulus its
returned fields have size O(L^-2) on canonical width O(L^-1).
Inside it the return-vector change is a times O(L^-3); hence its
conformal component is O(L^-3/Y), and its shape component is
O(L^-3 Y). The retained coefficient jets through three suffice for
these differentiated estimates. They give

    J2[Delta T]+J2[Delta(S+T/3)]=O(L^-1),
    integral log(t)(Delta_bar Delta T)^2=O(L^-1 log L),
    ||P_ell^2 Delta T||^2=O(L^-1),
    ||P_ell^3 Delta T||^2=O(L).                     (13)

The geodesic positive tower obeys the upper bound by a constant times
C p^4[1+log(2+|p|)]. Indeed use g(z)<=C min(z^6,z^4) in its
definition, rescale the finite Bessel moments, then split the angular
sum at M=|p| using N(M)<=4CM^2. Below that split the bound is
C p^4 log(2+|p|), and above it C p^4. Spectral interpolation/Jensen
with the last two lines of (13) gives O(C L^-1 log L). Thus the
approximation converges in N_ell and establishes the lower bound
for every compact shape T. Only the already controlled form-level
logarithm is used in this approximation; H4log convergence of the
scalar-c cutoff is not asserted.

Second approximate an arbitrary T in W_(3/2) by smooth compact
shape inputs in that domain. D_T is already bounded, and
E_ell[N T] is continuous there by the ordinary H1 York/return map.
Therefore

    <T,D_T T>>=.90 E_ell[N T]                      (14)

on the entire actual critical graph domain. This is not an
identification of that domain with N_ell or with plain c H4log.

Since ||Lc||<=18||c|| from the retained return bound,

    <T,D_T T>>=.90lambda ||T||^2/(8piG*18^2).

Its critical kernel is zero; Fredholm index zero proves that D_T
is onto. The inverse is bounded, symmetric, positive and injective
on L2, with dense range. Its inverse is consequently the positive
selfadjoint D_T on the graph domain already identified above.

The scalar-c inverse is exactly

    D_c^-1=N D_T^-1 N*, Dom(D_c)=N W_(3/2).         (15)

The same bounded-inverse argument identifies its positive selfadjoint
realization. Equation (14) gives (1). Spectral truncations extend its
.90 energy coercivity to Dom(D_c^(1/2)), since the ordinary TT energy
is closed. This now licenses the spectral square root required by
the mixed estimate.

## 5. The actual mass column and fast TT isomorphism

The two-power normal remainder bootstraps a critical solution of
D_T T=f, f in Y_(7/2), first to weight 5/2 and then to weight 7/2.
Only the simple TT root 3 is crossed. The finite asymptotic construction
in the weighted-map proof gives

    T=a e_T+v, e_T=chi(z)z^3, v in W_(7/2).         (16)

Local corrections to e_T start at degree five; separated York/kernel
tails are z^4 log^N z. They lie strictly inside the fast space. Also

    R_T e_T=-(1/2)e_m+X_(7/2),
    K e_T=-(1/6)e_T+Z_(7/2).                       (17)

The physical mass mapping is therefore bounded by the established
mass theorem for A_disk and (8). This proves

    D_T: W_(7/2)+span{e_T} -> Y_(7/2)

is an isomorphism; uniqueness follows from the critical inverse.
In scalar-c language the domain is
W_(7/2)^c+span{N e_T}. Its normalized mass vector -(3/2)N e_T
has leading coefficient z^3 and returns (3/4)e_m modulo X_(7/2).
It is the SAME one physical mass parameter. Replacing this vector
by a bare scalar chi z^3 would require a separate logarithmic-domain
claim and is not needed. Equations (8),(17) identify the complete
physical mass/fast disk domain with these TT and scalar domains.

## 6. Apply the actual mixed estimate and identify the coupled range

The [weighted mixed theorem](finite-comparison-weighted-mixed-form.md)
now applies to c=D_c^-1/2 y in its proved spectral form domain.
It retains the actual scalar operator, proper length, volume,
density/gauge contact and full trace. The exact positive elimination
and limiting-absorption construction give

    ||(I+Acal Jcal)^-1 Lcal Jcal||<.025,
    ||(I+Acal Jcal-Lcal Jcal)^-1||
                      <=1.000002/(1-.025)<1.026.   (18)

This is an actual outgoing inverse in the declared weak source and
auxiliary response norms. It is not a numerical fast H4/log inverse
constant. Retain the positive elliptic subtraction at low frequency.
The actual scalar Bessel realization fixes both regular bolt conditions;
no additional bolt datum or independent radial TT tail is introduced.

The reconstruction and kernel exclusion in the
[seed outgoing argument](stationary-scalar-outgoing-inverse.md) now
use the actual finite operators. For clarity about weights, omega=t^-1/2
makes the weak scalar response an outer z^1 H_b^1 field. The scalar
equation upgrades it to z^1 H_b^3. The actual mixed row has at least
two powers of outer gain (fixed asymptotic coefficients, with its
z^4 York tails); its transpose therefore produces critical L2 TT
forcing. Only ordinary integer mapping through three is needed here.
The critical TT inverse identifies the weak reconstruction with the
operator just constructed. Fast sources receive the permitted mass
term and scalar outgoing oscillation, with all actual forced tails.

A physical homogeneous mass/fast solution has conformal component
tau in Z_(7/2) by (8),(17). Subtract its fourth-order reference
outgoing response. The difference splits exactly into regular solutions
of Tg+4lambda and Tg-Esp. The first has no growing branch and is
zero by positivity; the second is outgoing and zero by the actual
scalar reference uniqueness. Thus the homogeneous solution obeys the
auxiliary equation. Equation (18) forces tau=0, and the actual TT
inverse then forces c=0. The already proved physical index -1 gives
(2), without inferring a coupled inverse from TT positivity alone.

The [general adjoint theorem](stationary-adjoint-phase.md) now leaves
one real mass-orthogonal adjoint line chi_ell. Its conformal leading
oscillatory coefficient is nonzero; otherwise it would be a physical
mass/fast kernel field. With psi_ell its conformal component, the
physical range condition is exactly

    <chi_ell,f>_ell
      =<psi_ell,f_tau-Bmix D_c^-1 f_c>_ell=0.       (19)

This phase belongs to the actual finite metric and proper length.
Neither it nor the signed preparation moment is evaluated here.
The inverse on this closed range exists, but its useful full weighted
norm and its composition with the actual residual remain to be bounded.

## Verification and next decision

Freeze this note and checker before execution. Pin the weighted mixed,
global weighted, seed TT and actual return receipts. Reconcile every
input hash. Check (4)-(6), the inverse-coordinate/normal congruences,
mass normalization, physical transpose, finite regularity allocation,
compact-core scaling, spectral upper-bound split and numerical margins.
Use manufactured regular radial geometries to verify (6) independently
of its shortened formula: a=rho(1+rho^2/7+rho^4/101),
r=2+rho^2/5+rho^4/103, Y=rho(1+rho^2/11), and
T=r^-4[r(0)^4/13+integral_0^rho r^4 P Y]. At 90 digits
check rho=1/5,2/3,7/4 with tolerance 1e-60. These enlarged
fixtures are not the physical comparison. A C5 high-frequency control must detect the
unjustified scalar-c logarithmic upgrade. These checks support ordinary
parametrix, density and spectral arguments; they are not formal
verification, independent review or a numerical physical solution.
Preserve all scientific runs, including failures and exact snapshots.

The first numerical fixture run failed its York residual tolerance:
the uncancelled polar expression inside near-zero quadrature produced
a residual about 1e-10. Its receipt retains the failed note and checker.
For the corrected run, cancel the removable polar singularity in P Y
symbolically before quadrature, check its zero bolt limit and compare
the regularized expression with the raw differential formula at each
positive test radius. Keep the 90-digit precision and 1e-60 residual
tolerance. This changes the fixture evaluation, not the claimed identity.

The next quantitative task keeps Pi_cut fixed at the seed. In particular
the projected finite derivative also needs <chi_ell,p_cut>_ell nonzero
with a useful bound; the one-dimensional cokernel alone does not prove
that transversality. Cold/higher response, preserved nonlinear family,
useful inverse/residual composition, projected solve and signed
compatibility remain open. S167 stays active with every task unchecked.
No Euclidean estimate here proves dynamical stability or a wormhole.

Remaining to-do list: fixed-projection transversality and useful full
weighted response; cold/higher bounds and nonlinear preservation;
projected solve and signed compatibility; exterior/state, charge/return,
opening, stability, transport, preparation, separate creation and audit.
