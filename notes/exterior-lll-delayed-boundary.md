# Exact free LLL elimination of a static spherical exterior

2026-10-15, E003/S134.2, supporting S074. Free lowest-Landau-level
propagation from finite throat cuts through a fixed spherical exterior
to the actual paired AdS boundary gives an explicit delayed boundary
response. Both propagation legs and the source's scattering time are
retained. The cuts introduce no physical wall or new source contact.
This is a free LLL response kernel, not the complete interacting quantum
influence functional or a new opening result.

## Scope and tests fixed before execution

Use the physical paired condition of
[the coupled-domain note](coupled-fermion-domain.md), the source convention
of [boundary-vertex matching](boundary-vertex-matching.md), and fixed
static spherically symmetric magnetic RN-AdS exterior segments with
positive f from the cut to conformal infinity. The angular spinors are
normalized as in the physical source. The zero angular eigenvalue is
used exactly; every nonzero angular channel has a radial mass/potential
and requires separate elimination. No omission bound is claimed for it.

Before the checker, freeze: projector, canonical flux and spectral
identities within 3e-11; direct characteristic propagation versus the
two-leg kernel within 2e-10; independent source-timing and proper-clock
half-density tests within 2e-9; analytic Sobolev error bounds permit
relative 2e-8 rounding. Integration uses a fixed smooth compact packet
with quadrature refinement and endpoint tails controlled by compact
support. Tests use unequal delays (0.07,0.19) and (0.11,0.31), with
scaled delays 1,1/2,1/4. Incorrect source timing and omission of a
propagation leg must have a nonzero declared discrepancy >1e-4. The
instantaneous approximation is not required to pass an all-frequency
norm test; its failure is proved below. Numeric spectra supplement the
exact characteristic determinant.

## Canonical radial and clock factors

For one static exterior write

```text
ds^2=-f(r)d tau^2+dr^2/f(r)+r^2 dOmega^2,
f(r)=1+r^2/ell_AdS^2-2GM/r+r_e^2/r^2,
dx=d r/f(r), Psi=f(r)^(-1/4) psi(tau,x) tensor eta_m/r.  (1)
```

For a normalized angular zero mode eta_m the four-dimensional spinor
norm reduces exactly to integral dx psi^dagger psi, and its integrated
radial flux to the canonical chiral flux. The spherical spin connection
and conformal factor cancel in the massless radial equation. Thus the
LLL chiral amplitudes propagate without a radial potential in x. This
uses a fixed magnetic bundle and no electric potential; a prescribed
electric connection would add parallel-transport phases, and fluctuating
gauge fields are not included.

If the common dimensionless time is t=a tau with fixed a>0, then
x_t=a x and psi_t=psi_tau/sqrt(a). At a static cut r_ci its proper
time satisfies du_i=sqrt(f_ci)d tau. The canonical temporal endpoint
field is psi_ui=sqrt(d tau/du_i)psi_tau=f_ci^(-1/4)psi_tau.
Equivalently psi_ui=r_ci times the angular projection of Psi. These
half-density factors preserve flux measures and the CAR normalization;
one cannot omit them and retain the same port norm.

The one-way exterior time in t units is

```text
ell_i=a integral_(r_ci)^infinity d r/f_i(r).               (2)
```

It is finite for cuts above any horizon with f>0 throughout the segment.
For the RN-throat matching convention of the source, a=C/L_th with
C=1+6 rbar^2/ell_AdS^2. L_th is the wormhole matching length, not the
AdS curvature radius. No arbitrary return time is fitted here. Actual
cut choice, throat/exterior overlap errors and evolution of the metric
remain separate physical data.

## Scattering matrix from the actual paired projector

In the canonical collar representation let alpha=sigma_2, beta=sigma_1,
A=tau_3 tensor alpha and
chi_theta=cos(theta) tau_3 tensor beta-sin(theta) tau_1 tensor 1.
The allowed trace is the +1 eigenspace of chi_theta. For inward collar
coordinate r, A=+1 propagates from the boundary into the interior and
A=-1 toward the boundary. Put e_+=(1,i)/sqrt(2), e_-=(1,-i)/sqrt(2).
Choose incoming columns (L:e_-, R:e_+) and outgoing columns
(L:e_+, R:e_-). Projecting chi_theta between them gives

```text
B(t)=S_theta(t) A_in(t),
S_theta=-i cos(theta) I_2-sin(theta) tau_1
       =-i exp(-i theta tau_1), theta=2 atan h_b.          (3)
```

These A_in,B are amplitudes at conformal infinity; the matrix A above
is the distinct characteristic-velocity matrix. The phase -i is the
uncoupled reflection phase, not a source-dependent interaction.
S_theta is unitary and S_theta S_0^dagger=exp(-i theta tau_1).
The constant flavor rephasing diag(1,i) maps this relative rotation
to exp(J theta), J=[[0,-1],[1,0]], used in the returned-switch model.
Retaining the static propagation/reflection phases is essential to
the antiperiodic spectrum. A different incoming/outgoing basis changes
these matrices but cannot discard the total loop holonomy.

## Exact two-leg cut kernel and source timing

Let a_i(t) leave the retained throat cut toward infinity, and b_i(t)
return into that throat. Characteristics on exterior segment i give
A_in,i(s)=a_i(s-ell_i) and b_i(t)=B_i(t-ell_i). Therefore

```text
b_i(t)=sum_j S_ij(t-ell_i) a_j(t-ell_i-ell_j),
K_ij(t,t')=S_ij(t-ell_i)
                        delta(t'-t+ell_i+ell_j).        (4)
```

The source is sampled when the ray scatters at infinity, t-ell_i,
not at its return to the cut or at its earlier departure from cut j.
At finite initial time, omitted past a_j values stand for the initial
exterior field. They cannot be replaced by zero for the prepared
correlated vacuum. One may instead give the full past history; (4)
then supplies the causal boundary graph without a separate bath.

For cut clocks t_i=t_i(u), with positive derivative, set
a_ui(u)=sqrt(t_i'(u))a_i(t_i(u)), and similarly for b. Define v_ij by
t_j(v_ij(u))=t_i(u)-ell_i-ell_j. Exact transport becomes

```text
b_ui(u)=sum_j sqrt(v_ij'(u))
             S_ij(t_i(u)-ell_i) a_uj(v_ij(u)).            (5)
```

This is the required delayed clock half-density, even if the two cut
proper clocks differ. It concerns canonical fields. Quantum point-split
stress also has its geometric/Schwarzian terms; (5) alone is not their
renormalization or the nonlinear clock constraint.

## Unitarity, causal memory and boundary work

On full time-line L2 define (D_ell a)_i(t)=a_i(t-ell_i) and
(M_S a)(t)=S(t)a(t). Equation (4) is the ordered factorization
R=D_ell M_S D_ell. Both delay factors and M_S are unitary, so R is
unitary even for a time-dependent unitary matrix S. The two unequal
delays generally do not commute with the source mixing. The map is
causal with fronts ell_i+ell_j, with source support delayed by ell_i
on the outgoing leg. Unitarity on the whole time line is not an
instantaneous equality of flux at two different stored-history times.

For constant theta, R(omega)=exp(i omega L) S exp(i omega L),
L=diag(ell_i), in the exp(-i omega t) convention. For Laplace s with
Re s>0 the corresponding factors are exp(-s L), and
||R(s)||op<=exp[-2 Re(s) min ell_i]. This is the lossless/passive
static one-particle return. It does not mean a driven source does no work.

For smooth classical amplitudes let
e[A]=(i/2)(A^dagger A'-A'^dagger A). At infinity,

```text
e[B]-e[A_in]=i A_in^dagger S^dagger S' A_in,
i S_theta^dagger S_theta'=theta' tau_1.                  (6)
```

The finite exterior stores the difference between this boundary work
and the flux imbalance at the cuts. A time-dependent source can supply
or absorb energy even while the charge norm is conserved. For quantum
fields (6) requires the specified common point splitting, initial
covariance and source contacts; it is not an unsubtracted quantum
stress formula. The returned-switch correlations cannot be replaced
by independent reservoirs.

## Explicit finite derivative-expansion error

Let ell_max=max ell_i and assume S is C2 in time, unitary, with bounded
derivatives. For a in H2 of the whole time line, the zeroth-order
instantaneous approximation satisfies

```text
||R a-S a||2 <=ell_max[||S'||inf ||a||2+2||a'||2].         (7)
```

There is a controlled first-order approximation with correct ordering,

```text
T1 a=S a-L S' a-(L S+S L)a',
||R a-T1 a||2 <=ell_max^2/2 [||S''||inf ||a||2
                         +4||S'||inf||a'||2+4||a''||2]. (8)
```

Proof. Put F(s)=D_(s ell) M_S D_(s ell) a, 0<=s<=1.
The delay generator is -L partial_t. Differentiating gives F'(0)
equal to the two displayed first-order terms. Unitarity of the delays,
the product rule and ||S||op=1 bound F' by the bracket in (7).
The two same-leg second derivatives and the two mixed derivatives
give the bracket in (8), multiplied by ell_max^2. Integral Taylor
remainders prove both inequalities. Density extends the estimates
from smooth compact packets to H2. For the physical scalar angle,
||S'||op=|theta'| and ||S''||op<=|theta''|+|theta'|^2.

For constant S on a frequency band |omega|<=Omega, the degree-n
Taylor remainder is bounded by (2 ell_max Omega)^(n+1)/(n+1)! times
||a||2. The expansion retains all ordered L^p S L^q terms with
p+q<=n. No source-bandwidth assumption bounds the bandwidth of the
entire fermion vacuum. Indeed for common ell>0 choose
omega=pi/(2ell): R(omega)=-S and ||R(omega)-S||op=2.
There is no all-frequency operator-norm smallness as ell tends to zero.
Thus the finite Sobolev estimate cannot itself be assigned to the
S076 quantum stress/area budget.

## Exact spectrum, retained period and influence-functional boundary

Attach a central free interval of optical length T between the cuts.
For a stationary mode a=exp(i omega T)tau_1 b. With constant theta,

```text
det[I-exp(i omega T)R(omega)tau_1]
 =2 exp(i omega P)[cos(omega P)+sin(theta)],
P=T+ell_L+ell_R,
(1+h_b^2) cos(omega P)+2h_b=0.                           (9)
```

This reproduces the source spectrum when P=pi. Exact elimination
preserves P; if the full interval has width pi then the central
finite-cut width is pi-ell_L-ell_R. Keeping that central width at pi
and also adding the exterior legs changes the model and the exact
Fourier-zero switching condition. An actual RN-to-throat matching
must fix the common time normalization and full optical width; this
identity does not assert its corrections vanish in the opening model.
For a product-throat seam at symmetric rho_m the central width is
2 atan(rho_m); the leading one-way exterior approximation
pi/2-atan(rho_m) replaces each removed strip end. Only their total,
with the actual deformation and exterior errors retained, belongs in
P. This static normalization is not the evolving clock y(t).

Equation (4) is an explicit nonlocal response graph. In a mixed
classical polarization it admits the quadratic generating functional
G[bar b,a]=integral dt dt' bar b_i(t) K_ij(t,t')a_j(t'), with the
canonical derivatives giving b=R a and bar a=bar b R. This states
which boundary data are fixed; an overall i depends on the fermion
action convention. It is not yet the quantum influence action.
That action also requires the initial exterior covariance and its
correlations with the throat, the contour/occupation prescription,
the determinant normalization and the finite source contacts. The
retarded response by itself fixes none of those missing noise/state
terms. The prepared strip vacuum supplies some of this information
only within its declared free benchmark.

Under a common boundary-clock pullback t=F(u), a source kinetic
functional c integral dt (partial_t theta)^2 becomes
c integral du (partial_u theta)^2/F'. The two half-density factors
also pull back a classical transfer bilinear exactly. The coefficient
c remains the matched coefficient at its original boundary and in
its original subtraction. Exterior elimination does not turn a finite
throat cut into that boundary, measure c, or permit replacing the
nonlocal source timings in (4) by one instantaneous cut angle.
The distinct local product-cut contacts and exterior contacts retained
in S133 remain separate obligations.
The companion [RN finite-cut calculation](rn-ads-finite-cut-matching.md)
gives the exact RN A, B and B_J coefficients. The
[prescribed-source clock comparison](prescribed-source-clock-variation.md)
distinguishes fixed laboratory data from a fixed global-time profile;
the contact pullback above does not select one of those ensembles.

## Provenance and verification

[Bintanja--Freivogel--Rolph v2](https://arxiv.org/html/2308.00871v2),
equations (2.2)-(2.7), gives the canonical radial rescaling, angular
normalization and massless LLL characteristics; equations (2.8)-(2.14)
give the paired condition and shifted spectrum; (2.23)-(2.25) give
RN-AdS and the throat/exterior time relation. Only those physical
inputs are imported. Equations (3)-(9) above are derived directly.
The scattering/characteristic and Sobolev methods are established;
no novelty claim is made. The existing boundary-vertex and returned-
switch notes fix the free source convention and prepared covariance.

The checker compares actual-projector scattering with direct boundary
constraints, both unequal-delay legs with characteristic propagation,
source timing and proper-clock half-densities, noncommuting delay/source
ordering, exact spectral determinants and Sobolev/Taylor errors on
smooth packets. It includes a general unitary source comparator as a
test of the operator estimates; that comparator is not substituted for
the physical one-angle junction. Its UTF8/LF receipt and log record
the exact runtime, inputs, frozen tolerances and source hashes.

Remaining to-do list: match the actual optical widths and evolving
clock, initial correlated state and quantum contact/determinant terms;
eliminate or bound nonzero angular modes and interactions; propagate
the resulting stress and mixed errors to the S076 clock, energy and
area observables.
