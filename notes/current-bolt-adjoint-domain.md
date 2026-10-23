# The physical bolt functional and its distributional adjoint

2026-10-23, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress. The full stable semiclassical wormhole
objective remains incomplete. The feasibility audit stays closed.

This identifies the actual source and inner domain of the proposed
defect/adjoint-bolt calculation. It does NOT evaluate the full pairing.
The adjoint for bolt evaluation has a point-supported source and is not
a smooth homogeneous bolt field. Its conformal component has a fixed,
nonzero rho^2 log(rho) term. Imposing the smooth forward bolt domain on
this adjoint would erase the functional one intended to measure.

## 1. Physical functional, variables and measure

Keep the current action c_R=1e22, G, full Dirac tower, fixed mu0=1/b0,
prepared positive temperature, paired matter realization and stationary
mass/fast phase. Use the actual finite seed of
[the full Hessian realization](current-seed-full-hessian-realization.md).
Let R=B(1+h) be its nonextremal bolt and gamma=kappa B/C.
Here gamma is the temperature ratio, not the canonical metric.

Write a tangent after the fixed York return as

    delta g=2s g+k+Lie_eta g,
    k=diag(-b_T-2xi,b_T,xi,xi),
    (r^2 F b_T)'=-(r^2 F'-2rF)xi.

Smoothness gives b_T(0)=-xi(0). A smooth rotational vector has
eta^ell=O(ell) in physical proper distance ell and dr/dell=O(ell),
so eta^r=O(ell^2) at the bolt. It contributes no bolt-area variation.
Therefore the PHYSICAL normalized radius functional is

    L_b(s,xi)=delta R/(B h)=ell_b [s(0)+xi(0)/2],
    ell_b=(1+h)/h.                                      (1)

This also equals V(0) to first order in the common disk radius variable
r=B+(r0-B)exp(V). It is neither a lapse evaluation nor just s(0).
The factor 1/h is required by the stated <1.667e-9 cold difference
target; it must not be removed by calling a rescaled functional the same
observable. All quantities in (1) belong to the CURRENT reference.

Use canonical proper rho with angle period 2pi and gamma_metric=
d rho^2+a(rho)^2 d theta^2. Then a=rho+O(rho^3), r=R+O(rho^2),
g4=(r^2/C)gamma_metric+r^2 dOmega^2. With N identical folded exteriors,

    dmu_rad=N (8pi^2 r^4/C) a(rho)d rho,
    ell=(R/sqrt(C))rho+O(rho^3).                        (2)

Let delta_b be the radial point distribution normalized by
<delta_b,v>_dmu=v(0), including the sphere and angle integrations.
For the action Hessian H=[-A_s,B;B*,D_T], the exact adjoint equation is

    H psi=(ell_b delta_b, (ell_b/2)delta_b).              (3)

The queried radius is the common radius, not a sum over mouths. Doubling
the action/measure doubles H and halves the adjoint field; its pairing
with the doubled action defect is unchanged. We keep N explicit to
prevent a factor-two ambiguity. The previously normalized affine disk
rows require their actual York dual map before use in (3).

## 2. Which space contains the adjoint

Use the established closed form space V with
||w||_V^2=a(s,s)+l(xi,xi) and ||H^-1 f||_V<=||f||_V*.
This paragraph explicitly depends on the actual finite-kernel/closedness
arguments in the retained full-Hessian realization; a checker or a
differential-boundary theorem alone is not a replacement for them.

On a compact bolt disk, the scalar form has a positive fourth-order
coefficient. The actual TT form has a positive elliptic order-four-log
principal part. To apply the interior argument AT the bolt, use smooth
Cartesian coordinates (rho cos theta,rho sin theta); the bolt is an
interior point, not a boundary of this disk. The direct short-time kernel
has the same local symbol as a smooth closed extension. Reflected and
state-dependent terms retain their lower orders. For a nonzero orbit
covector p, the leading TT lift has orbital matrix
-2xi(I-p tensor p/|p|^2) and spherical entries xi,xi. This order-zero
projection has trace zero and squared norm 6xi^2. The positive spin-two
symbol below is consequently elliptic in xi in every orbit direction.
This checks the bolt hypothesis that a strict radial-annulus statement
alone would not establish. The local parametrix gives an H2 (respectively H2 with
the logarithmic half-weight) form estimate up to an L2 term. The proved
positive lower form bounds absorb that L2 term. Thus V embeds locally
in H2 for both components, with a finite but UNCOMPUTED constant.
Point evaluation is continuous in H2 on the TWO-dimensional disk:
in flat Fourier variables Cauchy--Schwarz uses
integral_R2 (1+|k|^2)^-2 dk=pi. Smooth localization and the actual chart
give continuity at this bolt. Hence L_b belongs to V*, and (3) has a
unique form solution psi in V. No numerical bound on ||L_b||_V* follows.

The point source is not L2. Accordingly psi is not asserted to belong
to the smooth-bolt L2 operator domain used for smooth forward forcing.
Locally delta_b belongs to H^(-1-epsilon) for every epsilon>0.
The diagonal interior parametrices of orders -4 and -4/log, and the
lower-order mixed rows, put psi locally in H^(3-epsilon). The seed is
smooth on this disk; this assertion does not use a smooth reflected
outer seam. The mixed term B psi_T has no point distribution: order
at most two maps this local regularity into H^(1-epsilon).

One concrete construction is to approximate delta_b by normalized smooth
radial densities. They converge in the dual of local H2, hence in V*.
Apply the same form inverse and take the V limit. This construction
selects the point-source adjoint with no artificial puncture boundary
condition. It also explains why smooth approximating adjoints can converge
to an adjoint whose second derivative is unbounded at the bolt.

Away from the point source, the retained actual end parametrix applies.
Its mass cubic and fast/log remainder are kept. The outer antisymmetric
pairing of two such fields is zero: the mass/mass exponents give
z^(3+3-3)=z^3, all fast terms improve that power, and the actual-minus-
normal kernels gain two powers. Separated commutator tails use the
retained end decay estimates. This is the current positive-c_R phase;
the old zero-c oscillatory scalar adjoint/projection is not imported.
The stationary boundary one-form and its linearization still vanish
within this phase. No new zero-mass or finite-wall equation is added.

## 3. Exact scalar singularity and its sign

The scalar action block is NEGATIVE. Its leading form is

    -a(s,v)=-(6D_R/(8piG)) integral Delta4 s Delta4 v dmu+lower,
    D_R=96piG c_R-G/(60pi)>0.

At rho=0, Delta4=(C/R^2)Delta_rho+lower, where
Delta_rho=partial_rho^2+rho^-1 partial_rho. Substitution of (2) gives

    -K_s integral rho Delta_rho s Delta_rho v d rho,
    K_s=6pi N D_R C/G.                                 (4)

There is no R or kappa left in K_s; the physical factors canceled
algebraically. This is not a norm rescaling. With the radial measure
rho d rho, the distributional identity is

    Delta_rho^2(rho^2 log rho)=4 delta_rad,
    <delta_rad,v>=v(0).                                 (5)

Indeed Delta_rho(rho^2 log rho)=4log rho+4 and
lim_(rho->0) rho partial_rho Delta_rho(rho^2 log rho)=4.
Integrating twice over a punctured disk gives the same factor and sign.
The corresponding full two-dimensional identity has 8pi delta_2D;
mixing these two normalizations would be an error.

Equations (3)-(5) force

    psi_s=c_s rho^2 log rho+v_s,
    c_s=-ell_b/(4K_s)=-G(1+h)/(24pi N D_R C h).          (6)

The residual after subtracting this cut-off scalar singularity is L2
locally in the scalar equation: coefficient errors are O(rho^2), the
remaining scalar differential terms are lower order, and B psi_T has
the regularity stated above. Thus v_s is locally H4, in particular C2
with any Holder exponent below one. The limiting scalar flux, in the
weak Green-identity sense, is

    -K_s lim rho partial_rho Delta_rho psi_s=ell_b.      (7)

The coefficient c_s does NOT determine the regular value v_s(0), the
TT adjoint component, or the desired pairing. The TT component must
still realize the second point source in (3) through its actual
order-four-log operator, including contacts and the global state.
It cannot be set to zero or represented by the scalar singularity.
Both components are still unknown globally.

The TT source has its own universal leading profile. The local Dirac
symbol in [the retained heat-Hessian derivation](stationary-shape-logarithmic-response.md)
has logarithmic quadratic coefficient (3T2-T1^2)/(960pi^2), where
T2=tr[(K Pi)^2], T1=tr(K Pi). The TT lift above has T2=6xi^2 and T1=0.
Thus its per-four-volume coefficient is 3/(160pi^2). Converting with
(2), including p_physical=(sqrt(C)/R)p_rho, gives

    K_T=3N C/20,
    D_T,principal=K_T |D_rho|^4 log|D_rho|

in radial measure rho d rho. The fixed scale changes an ordinary
fourth-order term, not this leading logarithm. For the radial point
source, the high-frequency inverse kernel is

    K_T^-1 integral_e^infinity J0(k rho)/(k^3 log k) dk.

Its difference from its value at zero has leading term
-rho^2 log log(1/rho)/(4K_T). Split at k=1/rho and use
J0(x)=1-x^2/4+O(x^4) below that point; both remainder regions are
O(rho^2/log(1/rho)). A smooth low-frequency cutoff changes regular terms.
Consequently the TT adjoint has leading coefficient

    c_T=-ell_b/(8K_T)=-5(1+h)/(6N C h)

multiplying rho^2 log log(1/rho). This coefficient is evaluated as well.
It is an adjoint coefficient, not a metric perturbation or an inverse
norm. The large factor 1/h must be paired with the actual defect.

Subtracting only this TT leading profile is NOT enough to declare an
L2 residual: an ordinary fourth-order symbol acts on it with a remaining
singularity. In Fourier variables a missing 1/log(k) correction to the
point-source multiplier is not L2 in two dimensions. The complete local
parametrix with its prescribed fourth-order contact and lower terms is
needed for a regular-residual adjoint implementation. The local leading
profile does not determine its finite part or global homogeneous term.
This is an ultraviolet asymptotic as rho tends to zero. No numerical
transition radius is proved for it, and it must not be extended across
the .01 source-integration collar from the preceding calculation. Lower
second- and fourth-order terms can dominate at finite frequencies.

For a smooth radial scalar, Delta_rho psi_s is finite and the left
side of (7) is zero. Thus imposing ordinary smooth homogeneous bolt
data on (3) is incompatible with the nonzero physical functional.
This is a definite obstruction to that ADJOINT implementation, not an
obstruction to a smooth forward metric or to wormhole existence.

## 4. Current finite-temperature normalization

Compute h from the actual fixed-temperature root, not the extremal limit.
With v=1+h, u=1/C, epsilon=d/C, beta=(1-u-epsilon)/6, define

    K(v)=(v+1)[u+3beta(v^2+1)]/(2v^3)
                              -epsilon(v-3)/(2v^5).
    h K(1+h)=gamma, gamma=kappa B/C,
    tau=v h^2[beta(1+2/v+3/v^2)+u/v^2+epsilon/v^4].     (8)

These follow by imposing F(v)=0 and F_t(v)=2kappa B on the exact
rational seed. K(1)=1 exactly. The retained monotonic temperature
map selects its unique nonextremal root. The checker brackets it with
gamma_lower(1-1e-28) and gamma_upper(1+1e-28), verifying BOTH endpoint
signs with the full parameter box. It verifies the retained h and tau
guards and substitutes this interval into ell_b, K_s and c_s for N=1,2.
The result is an executed physical normalization and boundary-source
calculation, not a numerical evaluation of psi or of H^-1 r1.

## 5. Correct pairing and usable error certificate

For r1 in V* and e=-H^-1 r1 in V, symmetry of this SAME form gives

    L_b e=H[psi,e]=H[e,psi]=-<r1,psi>.                  (9)

The notation is a form/dual pairing. For smooth r1 it is the corresponding
volume integral plus any retained action boundary contacts; it is not
an instruction to integrate an unsubtracted singular kernel absolutely.
The outer boundary work vanishes in the declared phase. The bolt source
is already included by (3),(7), so it must neither be dropped nor counted
again as a second boundary contribution.

For an approximate adjoint psi_app with the SAME two point sources and
physical outer phase, let d_b=L_b-H psi_app in V*. Then

    |L_b e+<r1,psi_app>|<=||r1||_V* ||d_b||_V*.          (10)

To certify the cold target eta=1.667e-9, one needs

    |<r1,psi_app>|+||r1||_V* ||d_b||_V* < eta.           (11)

The full defect norm, computed pairing and adjoint residual norm remain
unknown. A small scalar singular coefficient does not meet (11). The
unit form inverse does not control X admission or the nonlinear map.
If psi_app is incorrectly smooth, its residual retains the entire
point source; it cannot be declared zero from off-bolt residual samples.

Decision: keep the direct physical H e=-r1 solve as the next numerical
milestone, with the actual weighted derivative remainder and admission.
The singular adjoint is an equivalent route only when BOTH point sources
and the regular remainder are computed. Do not replace that work by
another generic cold-value estimate or use a homogeneous old adjoint.
The current calculation supplies the exact functional/domain contract
needed to make either a direct bolt evaluation or an adjoint check valid.
The original nonlinear, physical matching, causal, validity and stability
requirements remain unchanged.

Remaining to-do list: complete the physical finite-domain correction or
the full two-component singular adjoint pairing with certified residual,
weighted return and admission; nonlinear closure and exact compatibility;
physical matching, causal realization, semiclassical validity, all-sector
stability, charge/return, opening, transport, preparation, separate creation
and the full stable-wormhole objective.
