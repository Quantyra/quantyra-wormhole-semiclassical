# Reciprocal dipole-pair response and its local-completion requirement

S072, 2026-10-13. The physical polar normalization makes it possible to
include reciprocal dipole-pair response, rather than only forcing the
metric with a prescribed matter fluctuation. In a specified cubic
representative, the complete nonlocal pair bubble is small relative to
both polar restoring terms, including the very light branch. Its
quadratic auxiliary-field realization has a correlated stationary
vacuum and positive masses. This is a calculation of one interaction
channel, not the interacting Einstein-Maxwell/fermion state. Local
terms, other interactions and the physical boundary drive remain needed.

## Which cubic information the existing equations determine

Retain the real normalized sphere quadrupoles Q_m, free canonical
dipole triplet phi_a, and fixed global product from
[the quadratic source](dipole-quadrupole-source.md). For each m,

```
(Box2 I-K)u=(v Box2+s0)Q_m,  u=(chi,Acal),
v=G(9/(4b^2),-3E0/4)^T,
s0=(3G(C-4)/b^4,9G E0/b^2)^T,
R=s0+Kv=(3G(C-1)/b^4,9G E0/b^2)^T.                     (1)
```

W is the positive physical kinetic matrix in
[the canonical polar calculation](polar-canonical-vacuum.md), and WK
is symmetric. The source (1) was reduced using the free dipole equations.
Consequently it specifies trilinear coupling information on free
dipole data, not the complete off-shell action. A reciprocal bulk
representative that reproduces it is

```
S_rep=S_phi,2+1/2 integral sqrt(-g2) u^T W(Box2 I-K)u
      -integral sqrt(-g2) u^T W(v Box2+s0)Q_m.           (2)
```

At this order variation of the reduced action on free dipole solutions
has the physical multiplier W fixed by the quadratic action. Differences
proportional to a free dipole equation are not determined by (1). Thus
(2) selects a representative; it is not asserted to be the already
matched full four-dimensional reduced action. Its reciprocal variation
with respect to phi is essential for the bubble below.

The local change u=y+vQ_m is triangular and has unit finite-regulator
Jacobian with phi held fixed. It gives the exact algebraic identity

```
S_rep=S_phi,2+S_y,2-integral sqrt(-g2) y^T W R Q_m
      +integral sqrt(-g2) [A (D Q_m)^2/2-B Q_m^2],
A=v^T Wv=27G/(32pi b^2),
B=v^T WR-v^T WKv/2=9G(2C-5)/(16pi b^4).                (3)
```

The displayed quartic terms are generated even by this one field shift.
Equation (3) is a classical or finite-regulator identity; products Q_m^2
and (D Q_m)^2 still require their own continuum renormalization.
They are retained in the action accounting, although they do not enter
the free-dipole two-vertex bubble at its leading order. They are not the
complete physical quartic action. Further dipole-equation terms, nonlinear
constraint elimination, quadratic-polar/two-dipole vertices, tadpoles and
counterterms must be matched consistently. In particular, those omitted
local terms can contribute at the same loop order as the bubble.

This restriction follows the distinction between first-order equation-of-
motion simplification and higher-order field redefinition in
[Criado and Perez-Victoria, section 3, especially equation 10](https://arxiv.org/html/1811.09413v2#S3).
We use that caution, not their flat-space scattering conclusions for an
AdS boundary observable. The physical field is still u=y+vQ_m; its
sources, boundary data and correlators must also be transformed.

## Canonical coupling and the crucial relative scale

Let S=diag(sqrt(b^2/(3kappa^2)),1), z=S y, and let e_r be the
orthonormal eigenvectors of M=a^2 S K S^-1 for L=6. The canonical
fields are Psi_r=sqrt(C mu_r) e_r^T z. In dimensionless global AdS2
coordinates the nondifferentiated cubic term in (3) is

```
-integral d tau d sigma sec^2(sigma) (g dot Psi) Q_m,
g_r=a^2 sqrt(C mu_r) e_r^T S R,
Mcan=diag(mu_-,mu_+),
q_g=g^T Mcan^-1 g
   =3/(8pi) (G/a^2)(1+25/C+28/C^2),
g^T g=3/(4pi) (G/a^2)(1+19/C+76/C^2+66/C^3).           (4)
```

These identities include the complete electromagnetic/gravitational
mixing. For example q_g=a^2 R^T W K^-1 R. The small light mass does
not produce a large inverse-mass enhancement in this particular
coupling. Keeping W, rather than assigning unit norm to the earlier
response variables, is required for (4).

There are five independent m blocks at this quadratic order because
the ordered free Q covariance is proportional to delta_mm'. The
operator norm bound below is the same on their direct sum; no additional
factor of five or independent charge rotor is inserted.

## Complete pair bubble and an analytic bound

Use the retained global Friedrichs dipole vacuum, with
Delta=(1+sqrt(9+16/C))/2 in [2,3]. Its ordered kernel is

```
<Q_m Q_m'>=delta_mm' (16pi/15) sum_j>=0 rho_j W_Dj,
D_j=2Delta+2j, lambda_j=D_j(D_j-1),
rho_j=(4Delta+4j-1)/(4pi)
       Gamma(j+1/2) Gamma(j+Delta)^2 Gamma(j+2Delta-1/2)
       /[Gamma(j+1) Gamma(j+Delta+1/2)^2 Gamma(j+2Delta)]. (5)
```

The positive AdS product expansion and its conventions were already
checked in [the dipole vacuum note](dipole-vacuum-response.md).
The primary formula is [Cacciatori, Epstein and Moschella,
Theorem 7.2](https://arxiv.org/html/2403.13142v2#S7). Here it applies
to the ordered kernel and its retarded commutator, not to an unordered
classical noise replacement. The scalar retarded propagator supplies
the usual unit delta-function jump; with the physical volume measure
its inverse is Box2-lambda_j/a^2. This fixes the factors of a^2 in (4).

For a real dimensionless Casimir z below the pair spectrum, define

```
beta(z)=(16pi/15) sum_j>=0 rho_j/(lambda_j-z),
Sigma(z)=beta(z) g g^T,
inverse polar operator Mcan-z I-Sigma(z).                (6)
```

The minus sign is the Schur complement of a reciprocal pair coupling.
It agrees with the connected second cumulant of (2). No additional
Wick symmetry factor is needed: 16pi/15 already includes both
contractions of the normalized traceless bilinear. In this selected
nonderivative representation the sum converges without a local
subtraction. That does not set the full theory's other local terms to
zero or fix a physical renormalization prescription.

The earlier Casimir gap sharpens, for all j, to

```
lambda_j-mu_+ >= 10+4/C+14j+4j^2.                       (7)
```

Log convexity of Gamma gives Gamma(x)/Gamma(x+1/2)
<=1/sqrt(x-1/2) for x>1/2. Apply this to the three ratios in (5):

```
rho_0 <= 11/[6sqrt(3pi)],
rho_j <= (4j+11)/(4pi j^2) <=15/(4pi j), j>=1.
```

For 0<=z<=mu_+, the j>=1 contribution to beta is at most
sum_j>=1 j^-3<=3/2. The j=0 contribution is below 176/675,
using 3<pi<4. Therefore

```
0<beta(z)<2377/1350<2,
0<beta'(z)<1/5,
0 <= Mcan^-1/2 Sigma(z) Mcan^-1/2 <= epsilon I,
epsilon=2q_g=3/(4pi)(G/a^2)(1+25/C+28/C^2)
       <14 G/a^2.                                     (8)
```

The derivative bound uses (7) once more. These estimates retain every
pair representation and all of its radial descendants through the
covariant propagator. They are not finite-mode extrapolations. A sum
through j=N-1 has an explicit remainder below 1/[2(N-1)^2], N>=2,
uniform on the interval in (8). The checker uses this tail to enclose
several numerical beta values, including the actual candidate.

For k=2000, the geometry-dependent epsilon is below 1.399e-33.
Thus this channel has a small *relative* effect even on
mu_-=4.973591971e-21. The corresponding light-pole absolute shift
bound is below 6.959e-54. These are bounds in the selected quadratic
comparison, not a total physical mass correction or capture margin.

## Reciprocal quadratic realization and its correlated vacuum

Introduce auxiliary canonical scalars B_j with masses lambda_j/a^2
and couplings c_j=sqrt(16pi rho_j/15). Their free generalized field
Q_aux=sum_j c_j B_j has exactly the two-point function (5). In the
quadratic comparison the combined mass form is

```
mathcalM = [[Mcan,     g c^T],
           [c g^T, diag(lambda_j)]].                    (9)
```

The sequence c is not square summable; it must not be treated as a
bounded unweighted coupling vector. Instead sum c_j^2/lambda_j=beta(0)
is finite. On the form domain x in R^2 and sum lambda_j |b_j|^2<infinity,

```
|2(g dot x) sum c_j b_j|
 <=sqrt[q_g beta(0)] [x^T Mcan x+sum lambda_j |b_j|^2].   (10)
```

This follows from two weighted Cauchy-Schwarz inequalities and 2ab<=a^2+b^2.
If epsilon<1, (8)-(10) give a closed strictly positive form equivalent
to the uncoupled form. The self-adjoint mass operator therefore exists.
Its compact resolvent follows from lambda_j growing quadratically and
the equivalent form norm. This constructs the continuum quadratic
comparison, not just a finite-regulator matrix.

For any finite truncation the Schur complement is exactly (6) with a
partial beta sum. Below lambda_0 the two eigenvalues of
Mcan-beta(z)gg^T decrease with z, while z increases strictly. Using
(8), each has one crossing, and passage to the closed-form limit gives

```
(1-epsilon)mu_r <= mu_tilde_r <=mu_r, r=-,+.             (11)
```

The same conclusion follows directly from the infinite Schur complement;
uniform convergence holds on 0<=z<=mu_+. There are precisely two
eigenvalues below the first pair mass. All remaining eigenvalues are
positive by (10). At a pole, a polar component x has total eigenvector
norm x^T[I+beta'(z)gg^T]x, accounting for its auxiliary components.
The residue correction is bounded by
beta'(z)g^Tg <= (mu_+/10)epsilon <=6epsilon/5.

Diagonalize the positive mass operator and assign each resulting scalar
the global Friedrichs ground state with weight
h=(1+sqrt(1+4mu_tilde))/2. Pulling this covariance back constructs a
stationary, positive, canonically normalized Gaussian vacuum with
nonzero polar/auxiliary initial cross correlations. It is not the state
obtained by imposing independent zero polar data at a finite time.
No artificial gravitational switching function is used.

The representation Q_aux matches the dipole bilinear only at the level
of its two-point kernel. The real bilinear has nonzero higher cumulants.
The auxiliary fields encode that kernel; they are not additional physical
particle species. Their separate vacuum energies must not be added to
the four-dimensional stress.
Hence this Gaussian construction resums the selected quadratic bubble
kernel; it is not an exact replacement of the physical dipole theory or
an all-orders interacting-state preparation. Physical u also retains
v Q_m, as in (3); sharp equal-time composite derivatives are not assigned
finite ordinary variances by this construction.

## Verification, decision and next physical step

The checker verifies the reciprocal-action field shift and generated
quartic terms exactly, independently checks the invariant coupling
norms, and encloses the full beta sum with an analytic tail. Finite
auxiliary mass matrices test pole brackets, the Schur determinant,
weighted positivity and residues. A separate spatial Galerkin
Hamiltonian verifies stationary pure Gaussian covariances, their
nonzero initial cross blocks and agreement with mass diagonalization.
These are diagnostics for the continuum form and spectral arguments,
not substitutes for them. Only the scoped original source formulas
are imported; no general novelty verdict is claimed.

Retain the positive nonlocal pair channel. Its small relative correction
removes a specific concern that the light mass automatically makes this
vacuum bubble large. Do not infer the full interacting mass or state
from it. Derive the remaining local reduced action and a common
renormalization prescription, including polar-polar-dipole-dipole terms,
tadpoles and constraint contacts. Retain the dipole H3/H4 sector, fermion
polarization, higher modes and physical boundary stress. Their errors
must be bounded in the physical norm and compared to the actual
capture/area margins. The fixed global vacuum and reflecting ends still
need matching to the opening setup.

Remaining to-do list: complete local interaction and renormalization
matching; physical interacting initial state and conserved boundary drive;
complete stress and evolving geometry; pair preparation, sufficient
stability and finite-energy transport; distinguish opening from creation.
