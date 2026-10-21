# The local coupled frequency polynomial

2026-10-20, S168. Ordinary derivation before focused verification.

The classical action and full four-dimensional Weyl restoration give
an explicit quadratic polynomial in squared imaginary frequency for
the two optical metric variables. This note computes both nonconstant
coefficients, including all conformal--shape mixing. The canonical
Dirac functional is independent of the four-dimensional conformal
variable, so these are also the COMPLETE frequency coefficients of
the conformal and mixed rows in this chart. The canonical shape row,
including its local subtraction terms, remains to be evaluated.

The local quartic coefficient is indefinite on the retained exterior.
Its direct shape entry contains a negative outer logarithm. This is
an actual term of the selected action, not the adjustable-polynomial
counterexample in the [thermal frequency note](thermal-shape-response-and-subtractions.md).
It explains what the canonical response must retain and compensate.
It does not establish an instability of a constrained physical mode.

## 1. Fixed action and optical density

Keep the selected action, common subtraction, magnetic flux, prepared
clock and domain. In Euclidean optical coordinates put

    g=e^(2Phi)[e^(2V)(dU_E^2+dx^2)+B^2 dOmega^2],
    y=B^-2, Q2=pi G q^2/g_c^2,
    A_C=B^2/(4G), A_W=B^2/(4pi),
    Delta=partial_x^2+partial_UE^2.                    (1)

Here Q2 is the classical magnetic length squared. It is not the
charge q or a spectral frequency. The selected cosmological term
has lambda0=lambda+11G lambda^2/(360pi); equivalently
lambda0=lambda+11epsilon lambda^2/6, epsilon=G/(60pi).

Use a fixed coordinate-density action and gradients. Take compact
metric variations away from the two endpoints; integrations by parts
then have no varying boundary contribution. The background itself
need not solve the equations. Differentiation in the logarithmic
variables is literal, so the action Hessian includes its chart and
volume contacts. It is not the Jacobian of a moving-volume equation
with its residual terms silently removed.

The exact [spherical Einstein--Maxwell reduction](exact-spherical-em-reduction.md),
with the selected cosmological term restored, gives

    I_C=-A_C integral dU_E dx F_C,
    F_C=e^(2Phi)[6|grad Phi|^2+4 grad Phi.grad V]
                   +e^(2V)U(Phi),
    U(Phi)=2y e^(2Phi)+6lambda0 e^(4Phi)-2Q2 y^2.    (2)

The overall Euclidean sign agrees with the retained trace convention.
For an independent unintegrated expression one may use

    F_C,raw=e^(2Phi)[-2Delta V+2y e^(2V)
                       -6Delta Phi-6|grad Phi|^2]
                 +6lambda0 e^(4Phi+2V)-2Q2 y^2 e^(2V).

The difference from (2) is a divergence. This is useful for testing
the sign without using the already integrated kinetic terms.

The [accepted spherical anomaly functional](spherical-anomaly-action.md)
has coefficients c=1/20, a=11/360, b=1/30. Write its Euclidean
four-dimensional Weyl restoration as

    I_W=-A_W integral dU_E dx F_W,
    F_W=e^(-2V){ Phi(Delta V)^2/15
                 -(Delta Phi)(Delta V)/15-(Delta Phi)^2/10
                 -(Delta V)|grad Phi|^2/15
                 -7(Delta Phi)|grad Phi|^2/90
                 -7|grad Phi|^4/180 }
             +(y/9)Phi Delta V-(y/18)|grad Phi|^2
             +y^2(1/15-q^2/3)Phi e^(2V).               (3)

To derive (3), substitute Rhat_2=-2e^(-2V)Delta V into the retained
Lpoly, use Boxhat Phi=e^(-2V)Delta Phi and
|gradhat Phi|^2=e^(-2V)|grad Phi|^2, and integrate
Phi Delta Rhat_2 twice by parts. The existing boundary Weyl action is
not replaced; it has zero variation for this compact core.
[Mazur--Mottola, section IV](https://arxiv.org/html/hep-th/0106151v1)
is the underlying established Weyl-cocycle construction. The novelty
here is its optical frequency Hessian, not a new anomaly action.

The remaining functional is Gamma_hat[V] in the same common scheme.
If it is expressed through flat optical mass Hamiltonians, its
two-dimensional measure restoration must be included there once.
Those canonical local terms have NOT been absorbed into (3) or set
to zero. No finite counterterm is changed in this decomposition.

## 2. Exact nonzero-frequency coefficients

Let Phi=P(x), V=W(x) be an arbitrary smooth static background and
take real compact amplitudes p(x),v(x). On the Euclidean circle use
delta Phi=p cos(omega U_E), delta V=v cos(omega U_E), with a nonzero
bosonic Matsubara frequency; put s=omega^2. Normalize the Hessian by
the time norm of cos, beta/2. Define the zero-frequency coefficient
using the static mode with norm beta. The resulting radial form is

    Q_loc(s)[p,v]=Q_loc(0)[p,v]+s Q_1[p,v]+s^2 Q_2[p,v]. (4)

Equivalently (4) is the Fourier-symbol quadratic form, so its local
polynomial continues to arbitrary s>0 and omega_L=i sqrt(s). The
static coefficient Q_loc(0) is the literal second variation of
(2)-(3) at time-independent fields, with the same fixed density.
Those formulas specify it fully; no unidentified static contact is
being added. Primes below mean partial_x.

The classical part is particularly simple:

    Q_1,C=-4A_C integral e^(2P)(3p^2+2pv) dx,
    Q_2,C=0.                                           (5)

The potential, including charge and lambda0, contributes only to
Q_loc(0). It affects the frequency terms through the background
P,W, not through an additional explicit time-derivative vertex.

For the Weyl restoration define the following local quadratic density:

    J_1=-24P W'' v^2+12P v v''
          +(P')^2(7p^2+14pv+12v^2)
          -14P' p p'-12P' p'v
          +P''(7p^2+36pv+12v^2)
          +W''(6p^2+24pv)
          -18p p''-6p v''-6p''v.

Then the exact answer is

    Q_1,W=A_W integral [e^(-2W)J_1/45
                             +(y/9)(p^2+2pv)] dx,
    Q_2,W=(A_W/15) integral e^(-2W)
                              [3p^2+2pv-2P v^2] dx.    (6)

Thus Q_1=Q_1,C+Q_1,W and Q_2=Q_2,W. The second spatial derivatives
in J_1 are intentional; polarization and integration by parts give
its second-order radial operator. The measure e^(-2W) must be retained
during that operation. There is no requirement that P or W be
constant and no use of background equations in (4)-(6).

One direct derivation expands the coefficient of the squared
perturbation in (3). The Laplacians of the perturbed fields give
p''-s p and v''-s v. The gradient square additionally gives
s p^2 sin^2(omega U_E); omitting that term would change Q_1.
Time integration sets the two squared trigonometric averages equal.
Expanding e^(-2V) is also necessary: it supplies several of the
P'',W'' and v^2 contacts in J_1. It does not suffice to replace Delta
by partial_x^2-s in just the explicit Laplacian factors.

Since Gamma_hat depends only on V, it adds to the V,V form only.
Consequently the p,p and p,v frequency terms in (5)-(6) are already
the complete corresponding rows of the selected action in this chart.
The full V,V row is the displayed known part plus the canonical
thermal response and its prescribed local polynomial. Independent
state sources and the initial/current constraints remain separate.

## 3. Actual sign and scale consequences

The dimensionless matrix of Q_2,W is

    M_2=[[3,1],[1,-2P]],
    det M_2=-6P-1,
    3p^2+2pv-2Pv^2=3(p+v/3)^2-(2P+1/3)v^2.             (7)

For r>=B one has P=log(r/B)>=0, so this local quartic matrix is
indefinite, including at P=0. This includes the retained exterior
comparison, whose area radius is above B. It is not a claim about
the signature of the full constrained action or physical energy.

For a purely shape amplitude, the quartic coordinate-density
coefficient equals -B^2 e^(-2W)P/(30pi). If the exterior has
N^2/r^2 -> lambda, then e^(2W)=B^2 N^2/r^2 -> B^2 lambda and

    coefficient_(v^2 s^2) ~ -log(r/B)/(30pi lambda).     (8)

This coefficient grows negatively toward the outer boundary. It is
not by itself a normalized operator estimate or an instability proof.
It rules out dropping the canonical response when proposing a
uniform positive local quartic coefficient. In the isolated quartic
matrix, adding a canonical shape term d in the same normalized units
would make it positive precisely when d>2P+1/3. This elementary Schur
test is only a diagnostic; the full canonical function contains a
frequency logarithm, and the lower-frequency/spatial terms and
constraints must still be treated. d is not a selectable counterterm.

There is an independent normalization check on the conformal entry.
With N^2=e^(2P+2W), the ratio of its quartic frequency term to the
magnitude of its classical quadratic frequency term is exactly

    [coefficient_(p^2 s^2) s^2]/
       [|coefficient_(p^2 s),classical| s]
           =G s/(60pi N^2).                            (9)

The two contact terms balance at proper imaginary frequency squared
s/N^2=60pi/G. This comparison does not discard the other coefficients
or define a mode dispersion relation. It agrees with the leading
scale E=60pi/G+11lambda/3 found in the exact AdS trace. At the retained
g=G lambda<(4e-17)^2, GE/(60pi)-1=11g/(180pi)<3.2e-35.
No reliability of the retained semiclassical action at Planck-scale
frequencies is inferred from this algebraic consistency check.

## 4. What the next estimate now contains

For any three positive-frequency anchors, (4) interpolates EXACTLY.
The previous spectral interpolation theorem therefore requires only
the remaining canonical row's frequency estimates in addition to
this now-explicit local polynomial. Its fixed common local terms,
including optical measure restoration, must accompany the nonlocal
mass-vertex response. The negative term (8) makes a constant local
positivity substitution unsuitable on the outer end; the actual
frequency/spatial dependence must be retained.

No canonical anchor value or uniform endpoint estimate has been
evaluated here. The direct shape row is still not the York TT block.
One must combine the complete two-row operator with the initial
Ward and any Maxwell/current constraints before deciding whether the
vacuum scalar obstruction transfers. Static equilibrium contacts,
off-shell background residuals and prescribed state/domain data
cannot be removed in that step. S168.3-S168.4 remain open. S167's
signed coefficient, nonlinear solve and compatibility remain required.

## Verification contract

Freeze this note, checker and source review before acceptance execution.
Reconcile the preceding thermal receipt and transitive input hashes.
Derive (5)-(6) symbolically from the local density with arbitrary
background and variation jets, retaining both time-gradient and
exponential contacts. Verify (7)-(9), the benchmark bound and exact
three-frequency interpolation of the local polynomial.

Independently evaluate the UNINTEGRATED curvature versions of the
classical action and accepted spherical Lpoly on a periodic optical
torus times S^2. The fixture fields and parameters are fixed in the
checker; this compact geometry is not an exterior or a physical state.
For Matsubara integers 1,2,3 and a static reference, recover the second
variation by complex-amplitude contours with 16/32 nodes at radius
.04 and spatial grids 32/64; use 32 time nodes. Compare with the
predicted frequency difference at relative scale 1+|answer|, tolerance
1e-9. Check real-amplitude second differences at .02/.01, tolerance
1e-3 and refinement ratio .35, allowing an absolute numerical floor
1e-9. Detect omission of time-gradient and exponential contacts.

The covariant curvature fixture uses Rhat_2, Boxhat Rhat_2 and the
unintegrated physical R4, independently of the reduced density (3).
These numerical and exact algebra checks do not prove a global
causal bound. No Lean is selected. Preserve all run outcomes and
failed exact sources; do not edit scientific sources after success.

Remaining to-do list: actual canonical positive-frequency estimates,
constraint/domain control and stability transfer or justified
exclusion; S167 signed coefficient, nonlinear solve and compatibility;
exterior/state, charge/return, opening, stability, transport,
preparation, separate creation and full-result audit.
