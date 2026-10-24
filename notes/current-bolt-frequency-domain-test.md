# Bolt frequency domain test: a mass gap is not a clock-frequency gap

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The stable-wormhole objective remains ACTIVE and INCOMPLETE. The feasibility
audit remains CLOSED. No physical source, domain, contact or normalization
is changed. This is a test of a proposed thermal estimate, not a new
evaluation of the signed quantum/contact projection.

## 1. Which domain statement is already justified

Retain current-rough-paired-variation.md and its inspected, immutable
Jacobian--Weyl and Navier--Stokes sources. Their transferable methods are
domain/adjoint analysis and explicit correction accounting, respectively.
Neither supplies a wormhole theorem or permits choosing a cancelling
quantum source. This milestone uses their methodological distinction;
it imports no numerical constant or analyticity assumption.

On the actual paired first-order domain, scalar-identical mass/log
fields with bounded b and db preserve both complementary adjoint domains
under multiplication by exp(+/-b/2). For the skew representative Dsharp,
A_M=Dsharp+M exp(b) has inverse norm <=exp(B)/M when ||b||infinity<=B.
Its composed positive square has lower bound M^2 exp(-2B). Keep the
regular disk bolt and prescribed antiperiodic spin structure. The square
domain is obtained by composition, not by a separately chosen Robin wall.

The accepted common-regulator result is a weak integrated variation on
bounded admitted mass/log families. The complete common-cutoff third
derivative in current-mass-response-regulator-audit.md retains two cyclic
triangles, three pair contacts and the triple contact. The logarithmic
b(Delta b)^2 mass contact survives. These facts do not give a uniform
strong H3 output bound or admission of the original nonlinear X ball.

There is a distinct prerequisite for the tempting argument that the
large clock period T=2pi/gamma makes ALL massive thermal effects smaller
than exp(-.01 M T). One would need an appropriate off-lattice analytic
frequency family with its domain and a uniform bound in that strip.
The Euclidean mass bound above is not that result. We test the proposed
global gap on the actual finite seed before interpreting this exponent.

## 2. Actual bolt and the optical first-order operator

Keep the correlated endpoint, rather than subtracting nearly equal
rounded numbers: Q(z_h)=0 and Q'(z_h)=-2gamma, gamma>0. The seed is

 Q=beta+u z^2+(-4beta-2u-tau)z^3
                     +(3beta+u+ep)z^4-2ep z^5+ep z^6.

On 0<=z<=z_h<1 the triangle inequality gives

 |Q''|<=60 beta+26u+82|ep|+6tau<12.                 (1)

Set x=z_h-z and delta=gamma/6. Taylor's theorem with the EXACT endpoint
identities implies, for 0<x<=delta,

 gamma x <= Q(z_h-x) <=3gamma x.                   (2)

The inequalities use |Q-2gamma x|<=6x^2. Define optical distance toward
the bolt by r(x)=integral_x^delta dy/Q(z_h-y). Then r goes to infinity,

 (3gamma)^-1 log(delta/x)<=r<=gamma^-1 log(delta/x),
 x<=delta exp(-gamma r),
 sqrt(Q)<=gamma/sqrt(2) exp(-gamma r/2).             (3)

After the usual density transformation of the stationary Dirac equation,
the generator conjugate to canonical time has radial expression

 H_M=-i sigma2 partial_r + M exp(b) sqrt(Q) sigma1. (4)

Directly, the retained Euclidean radial transfer equation becomes
partial_r psi=(M exp(b) sqrt(Q) sigma3+p sigma1)psi. Multiplying by
-i sigma1 gives (-i sigma1 partial_r+M exp(b) sqrt(Q) sigma2)psi
=-i p psi; a fixed unitary rotation gives (4), up to harmless signs.
Equivalently the orbit metric is Q(-dt^2+dr^2) after stationary
continuation and its conformal density removes the spin connection.
This is a
spectral test of the stationary expression; it does NOT select a causal
state, prove Lorentzian matching or assert a stability theorem.

The full paired fold retains its prescribed transmission condition.
On each optical end the potential is bounded and tends to zero by (3).
All packets below vanish near the prescribed projector and near the bolt
itself. They therefore belong to the compact interior core before any
issue about the global clock extension arises. A whole-line H1
transmission realization of (4) is self-adjoint by bounded Hermitian
potential perturbation of the free operator. More generally the gap
counterexample holds for EVERY closed realization containing this core;
the spectral conclusion below holds for every self-adjoint such
realization. We do not infer a physical Lorentzian extension from formal
symmetry or select one to resolve the still-open causal matching. No wall
is placed at a finite optical distance.

## 3. Executed continuum counterexample to the proposed global gap

Let chi(t)=sqrt(2) sin(pi t) on [0,1], extended by zero. It is in H1,
has norm one and derivative norm pi. Its endpoints have zero trace, so
extension introduces no delta in a first derivative. Smooth compact
approximants give the same conclusion on the smooth operator core.
Choose a unit vector v with sigma2 v=v. For any real E and R,L>0 set

 psi_(R,L)(r)=L^-1/2 chi((r-R)/L) exp(i E r) v.      (5)

Each packet has compact support strictly inside the physical bolt collar.
It is zero on an open neighborhood of the regular tip in disk coordinates
and of the outer paired boundary. Pair/reflection copies, if used in the
folded representation, are normalized together and obey the same bound.
Every finite packet satisfies the actual domain. A limit packet at the
tip is NOT added to that domain.

For every angular mass M, any admitted mass field with |b|<=B obeys

 ||(H_M-E)psi_(R,L)||_2
 <=pi/L + V0 exp(-gamma R/2),
 V0=M exp(B) gamma/sqrt(2).                         (6)

For any epsilon>0 choose L=2pi/epsilon and
R>=(2/gamma) max(0,log(2V0/epsilon)). Equation (6) is <=epsilon.
Let epsilon tend to zero and choose successive R so the supports escape
and are disjoint. The packets converge weakly to zero. This directly
proves the Weyl-sequence criterion at EVERY real E and EVERY M, including
zero energy, for any self-adjoint realization containing the actual
compact core. Such a clock spectrum has no positive massive gap. The
norm-gap counterexample itself needs no self-adjointness hypothesis. This is
not a finite-mode approximation to a determinant. One sector suffices
to refute a universal gap, and the proof in fact applies to each sector
of the entire angular tower and to the stated bounded mass families.

For an explicit executed scale test use M=M1, B=.001 and
epsilon=p0/100=gamma/200, p0=gamma/2. The defining correlated parameter
enclosures give

 gamma = 5.057799424826589366111e-34 (enclosed),
 M1 = 18.048130263729889289300 (enclosed),
 |Q''| guard <10.000000000326<12,
 delta = gamma/6 <8.430e-35,
 L <2.485e36, R <3.377e34,
 gamma R <17.078, gamma L=400pi,
 ||H_M1 psi||<=gamma/200<2.529e-36.                 (7)

R and L in (7) mean their exact formulas in (6); rounded displayed
upper bounds are not substituted into those formulas. The ratio of
this residual to the proposed .01 M1 gap is <1.402e-35. A putative
lower bound ||H_M1 psi||>=.01 M1||psi|| is therefore FALSE, even at
the seed. Taking packets farther out makes the residual arbitrarily
smaller. This is an actual obstruction to THIS GLOBAL GAP METHOD.

In particular any such self-adjoint clock resolvent has
||(H_M-i p)^-1||=1/|p| for real nonzero p, since zero is in its spectrum.
At p0 the norm is >3.954e33. This is a different operator and Hilbert
representation from the accretive Euclidean mass inverse; there is no
contradiction with its <=exp(B)/M bound.

## 4. Off-lattice frequencies, local thermal response and the regulator

The Euclidean disk Fourier numbers are omega=j+1/2 and p=gamma omega.
Near its tip the radial half-density equation has leading matrix
omega sigma1/rho, with regular branch rho^|omega|. Varying p in that
indicial expression produces (log rho)/gamma, and the selected branch
changes as real p crosses zero. Non-half-integral p/gamma is not an
angular spinor with the prescribed disk spin structure. A complex shift
of p is consequently not a fixed-domain Fourier-mode operation on that
disk. An auxiliary analytic continuation would require its own argument
and endpoint terms. It is not supplied by the Euclidean mass inverse.

Numerically T=2pi/gamma is about 1.24227648814e34, and .01 M1 T exceeds
2.242e33. The apparent fantastically small exponential is NOT a certified
thermal error: its global-gap premise has just failed. We do not enter
it in the feasibility ledger, and we do not delete thermal contributions,
local contacts, angular tails, images, or endpoint normalization terms.

Neither (6) nor a large global resolvent norm establishes a large LOCAL
renormalized response. Packets escaping the bolt can have negligible
overlap with the observed outer region. A localized response might admit
an analytic continuation or a small direct thermal bound despite this
global failure. Its value and uniform strong bound are UNKNOWN; no
excessive finite upper bound for that response has been computed here.

The existing cylinder comparison is still legitimate WITH its retained
actual-minus-cylinder return. It has Q_c bounded away from zero and a
different complete-end frequency problem. A frequency-strip estimate
proved there cannot be applied to the actual bolt without that return.
The executed seed return <2.012953e-13, including its seagull and diagonal
jumps, remains valid in its existing seed scope. Its numerical value is
not a bound on its metric derivative. The historical wording 'first two
metric variations' in that transfer means source and Hessian at the
seed, not a uniform derivative of the Hessian throughout a neighborhood.

For a uniform third response, differentiate the COMPLETE same-regulator
Duhamel combination, retain the bolt-normalized sensitivities and all
coincidence terms, and prove a common strong majorant before removing the
regulator. The accepted weak remainder does not justify passing this
limit in H3. The five-dilation W3,infinity proof remains independently
obstructed on the current fast graph; no smooth replacement is selected.

## 5. Unchanged quantitative requirements and decision

The executed angular strong component remains <3.038e-9 with a 4e-9
allocation. The cold/scalar component <1.467e-13 and the closed audit
are not reopened. The remaining nonlocal requirement is still

 eta_curved_seed+1e-11 C3 < .159916396.              (8)

Both constants in (8) are UNKNOWN. C3 means the COMPLETE uniform first
metric variation in the current graph, with actual finite-domain return.
Even if eta_curved_seed were zero, this sufficient test would require
C3<1.59916396e10. More generally its ceiling is
(.159916396-eta_curved_seed)/1e-11, provided the numerator is positive.
That ceiling is a target, not an estimate or evidence it is attainable.
The separate seed contact targets M_A<=3.18e27 and M_V<=2.66e30 remain
unknown, as do the original X/graph comparisons.

A proposed allocation for the NEXT domain calculation is .01 of the
radius-times-derivative budget: C_return<1e9 at graph radius 1e-11.
If proved, the OTHER terms must then satisfy
eta_curved_seed+1e-11 C3_other<.149916396. This is conditional bookkeeping,
not a new certified allocation. The seed return's smallness and spatial
separation make a localized proof worth testing, but give no numerical
bound on C_return. There is currently no evidence adequate to call that
target attained or unattainable. A global clock-gap proof cannot attain it.

At the original rho=1e-8, the unchanged preserved-ball and contraction
requirements are

 Y0+1e-8 Z0+5e-17 Z2<=1e-8,
 Z0+1e-8 Z2<1.                                    (9)

The sufficient allocations Y0<=5e-9, Z0<=.25, Z2<=5e7 give equality
in the first bound and .75 in the second. They remain UNPROVED; no
comparison factor between (8) and Z2 is assumed to equal one.
The earlier selected screen E_rest<2.123e-42<1.287e-40 is preserved
with its original scope: zero is an enclosure center, not a measured
signed projection. No full signed projection is executed this milestone.
No new preserved-ball, projected-solution or stability requirement closes.

RESEARCH DECISION: retain the construction; reject the proposed global
Boltzmann-gap shortcut and sharpen the LOCALIZED actual-domain method.
The executed continuum test resolves a false prerequisite before numerical
thermal interpretation. It supplies no nonexistence, instability or
divergence theorem for the actual quantum action. Park controlled cold,
scalar and angular components; do not optimize them to offset an unknown
domain constant. Uniform strong common-regulator convergence precedes
further source optimization.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the first metric
variation of the actual-minus-cylinder near Hessian in the current
coupled H3 graph, uniformly on the admitted mass/log family with the
regular bolt, including the complete seagull/diagonal terms and normalized
endpoint sensitivities. Use localized propagation or the full relative
heat/resolvent identity, retain all angular and thermal modes, and test
C_return<1e9. Its proof must specify the global coefficient/admission
hypotheses and the common regulator limit; a seed value cannot replace it.
Passing that target would close one finite-domain neighborhood requirement,
while the local curved response, contact jets and original X comparison
would still be required for (8)-(9). A failure must again distinguish an
excessive estimator from an actual response obstruction.

Remaining to-do list: uniform actual-domain return, complete curved strong
response and contacts, global correction, X/graph comparison and nonlinear
admission/contraction, mass/slow matching, exact compatibility, physical
matching/causal realization, semiclassical validity/noise, all-sector
stability, charge/return, opening, transport, preparation, separate creation
and the full stable-wormhole objective.
