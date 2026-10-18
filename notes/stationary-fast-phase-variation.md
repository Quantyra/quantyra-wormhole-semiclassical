# A differentiable action on the fast stationary phase

2026-10-17, S139/S074. Ordinary argument, reviewed before focused
verification. This extends the conditional mass-channel calculation in
[the charge note](stationary-same-source-boundary-charge.md) to every
variation in an explicitly restricted stationary phase.

## Statement and boundary data

Fix G, q, B, mu0, the paired reflecting Dirac domain, theta=0 and the
bulk subtraction of [the centered seed](stationary-centered-covariant-seed.md).
Write lambda=lambda_e, with
`lambda+11G lambda^2/(360pi)=lambda0`. Fix a positive Euclidean period
beta in the coordinate U. The metrics are spherical and stationary,

```
g=exp(2alpha) F dU^2+dr^2/F+r^2 dOmega^2,
F=lambda r^2+1-2m/r+R_F,
alpha=alpha_infinity+R_alpha,
partial_r^j R_F=O(r^(-2-j) L^p),
partial_r^j R_alpha=O(r^(-4-j) L^p), 0<=j<=5,       (1)
L=1+log(r/B).
```

Choose one finite p and one corresponding weighted topology; constants
are locally uniform in that topology. The compact interior is real C5,
the bolt is smooth and nonextremal, and F>0 outside it. The period and
bolt regularity are imposed covariantly. Coordinate horizon locations
need not be fixed when the metric varies. Use regular disk charts to
compare moving bolts. No estimate uniform as temperature goes to zero
is asserted here.

Work in an open neighborhood of the accepted comparison on which the
actual paired square has no kernel. This neighborhood exists at each
fixed positive temperature: the compact conformal geometry has a
discrete spectrum, the accepted domain has no zero mode, and its
relative C1 graph perturbations preserve a local positive gap. This
does not assert that every metric satisfying (1) has no kernel.

Allow alpha_infinity on an OPEN interval about zero. Tangents have
`delta F=-2delta m/r+delta R_F` and
`delta alpha=delta alpha_infinity+delta R_alpha`, with the same weighted
bounds. Families are C1 in this fixed topology, including their
directions, after a fixed endpoint identification. At the physical
source one sets alpha_infinity=0. Variation of the asymptotic clock
can be supported away from the bolt, preserving beta and bolt
regularity. Independent mass variations also admit such outer support.

The fixed coefficient 1 in F, fixed lambda, absence of the r^-3 lapse
term and exclusion of the oscillatory scalar tails are PART OF THIS
PHASE. They are not consequences of fixing an arbitrary boundary metric.
This is neither a metric-only Dirichlet theorem at a finite boundary
nor a completion for arbitrary nonspherical or time-dependent sources.

Define the action by the covariant finite-cut Dirac determinant and
its actual boundary Weyl transgression, the Einstein GHY term, the
Einstein AdS counterterms at lambda, and the minimal loop counterterms
already derived in the charge note. No further finite extrinsic term
is added. A common finite functional of the conformal sources may set
the source energy zero.

**Restricted variational theorem.** On this phase the renormalized
stationary action is finite and C1. Its outer boundary one-form per
exterior, with the time integral understood, is

```
Theta_infinity=-exp(alpha_infinity)
                   [Z m/G+E_sources] delta alpha_infinity,
Z=1+G lambda/(5pi)
        [log(mu0/sqrt(lambda))+gamma_E+151/180].      (2)
```

There is no residual surface work under any allowed fixed-source
tangent. The mass-dependent coefficient is insensitive to additions
that preserve differentiability for every nearby stationary clock
source. At alpha_infinity=0 the stationary boundary conjugate satisfies
`Delta E_stat=Z Delta M_geo` for the SAME full reference and sources.
The accepted benchmark gives `|Z-1|<1e-33`.

The theorem constructs a stationary action on a chosen phase. It does
not select that phase as the physical gravitational boundary theory,
construct an extremum, or identify a Hamiltonian on a full dynamical
phase space. Existence, uniqueness and a coupled inverse are separate
questions. In particular (2) is a boundary presymplectic one-form:
its mass part has exterior derivative
`-exp(alpha_infinity) Z/G dm wedge d alpha_infinity`.
It is not an exact thermodynamic differential on an independently
two-dimensional mass/clock space. The flux vanishes at fixed clock.

## Einstein variation and the bolt

Use the positive-mass sign convention of the charge note. Angularly
integrated bulk and GHY densities per U are

```
R4=-F''-3alpha'F'-2F(alpha''+alpha'^2)
     -4F'/r-4F alpha'/r+2(1-F)/r^2,
L_EH=exp(alpha) r^2(R4-2Lambda0)/(4G),
B_GHY=exp(alpha)[r^2(F'+2F alpha')/(4G)+rF/G].        (3)
```

For a second-derivative radial density write its endpoint form as
`P_F delta F+Q_F delta F'+P_alpha delta alpha+Q_alpha delta alpha'`.
Direct differentiation of (3) gives

```
Q_F=-exp(alpha)r^2/(4G),
P_F=-exp(alpha)(r^2 alpha'+r)/(2G),
Q_alpha=-exp(alpha)r^2 F/(2G),
P_alpha=-exp(alpha)r^2(F'+2F alpha')/(4G).
```

Adding delta B_GHY cancels both differentiated variations exactly:

```
Theta_EH+GHY=exp(alpha)/G [rF delta alpha+(r/2)delta F].
B_EHct=-exp(alpha)sqrt(F)/G
                  [sqrt(lambda)r^2+1/(2sqrt(lambda))],
Theta_EH,ren=-exp(alpha)E_EH delta alpha
  +exp(alpha)/(2G)
     [r-(sqrt(lambda)r^2+1/(2sqrt(lambda)))/sqrt(F)]delta F,
E_EH=sqrt(F)/G
     [sqrt(lambda)r^2+1/(2sqrt(lambda))-r sqrt(F)].     (4)
```

The square bracket multiplying delta F is O(r^-2); E_EH tends to
m/G. Thus all normalizable Einstein surface terms vanish. The vacuum
centering difference between Lambda0=-3lambda0 and -3lambda is kept
in the bulk and combined with the loop vacuum term below.

The bolt is an interior point of the Euclidean disk, not an additional
boundary datum. Use smooth covariant variations there. Flux through
a small proper-radius circle vanishes with its radius. An apparent
areal-coordinate endpoint term cannot justify imposing delta r_h=0;
that chart degenerates at the bolt. The fixed regular beta and disk
charts are the appropriate contract.

## All fast normal jets in the local loop action

Retain the complete density L_W and boundary polynomial B_W of the
charge note. With sigma=log(r/B),

```
T=r(F'+F alpha'),
Rbar=-r^2F''-3r^2F'alpha'-2r^2F(alpha''+alpha'^2)+2T-2F,
L_W=exp(alpha)/(4pi r^2) {
 sigma[Rbar^2/60-Rbar/18+1/15-q^2/3]
 +Rbar T/30-T^2/10+11TF/90-11F^2/180-F/18}.           (5)
```

Here `Rbar=-2+12m/r+O(r^-2 L^p)`. Differentiating (5), including
the derivative of each second-derivative momentum, bounds its full
momenta by

```
P_F=O(r L^N), Q_F=O(r^2 L^N),
P_alpha=O(r^3 L^N), Q_alpha=O(r^4 L^N)               (6)
```

for a finite N. Their products with the respective remainder tangents
`delta R_F, delta R_F', delta R_alpha, delta R_alpha'` are all
`O(r^-1 L^N)`, hence vanish. These estimates do not set a normal
derivative of a remainder to zero; they cover every remainder in (1).

For clarity, the literal boundary transgression can also be bounded
before taking a limit. Let D_x=(r sqrt(F)/B)partial_r point outward,
N=B exp(alpha)sqrt(F)/r, a=D_x log N, s=sqrt(F)/B, Y0=B^-2. Then

```
B_W=-B^2 N Gcal/(360pi),
Gcal=sigma[12a D_x a+44a^3/7+10aY0]
      -6(D_x a)s-24a^2s/7-6a D_x s-96as^2/7
      -18s D_x s-85s^3/21.                          (7)
```

For a remainder tangent at fixed alpha_infinity and m,

```
N=O(1), a=O(r^-1), s=O(r), D_x a=O(1), D_x s=O(r^2),
delta N/N=O(r^-4 L^p), delta a=O(r^-3 L^p),
delta(D_x a)=O(r^-2 L^p), delta s=O(r^-3 L^p),
delta(D_x s)=O(r^-2 L^p).                           (8)
```

Each non-sigma monomial in delta(N Gcal) is at most O(r^-1 L^N).
The sigma terms are O(r^-3 L^N) or smaller. Varying D_x itself is
included in (8). No finite remainder surface term survives.

The loop counterterms are

```
B_ct=4pi exp(alpha)sqrt(F)(c0 r^2+2c1),
c0=-221lambda^(3/2)/(15120pi^2),
c1=29sqrt(lambda)/(20160pi^2).                      (9)
```

Their remainder variation is O(r^-1 L^N). For the one remaining
normalizable leading direction delta F=-2delta m/r, the accepted
exact endpoint calculation gives, in units
`exp(alpha_infinity)lambda delta m/(3780pi)`,

```
bulk L_W endpoint: -105,
delta B_W:         -116,
delta B_ct:        +221.                            (10)
```

Their sum is zero. Equations (6)-(10) extend the earlier mass test
to ALL allowed local stationary tangents, not to arbitrary finite-cut
metric jets. Multiplication by a constant exp(alpha_infinity) does
not change this cancellation.

## Canonical determinant: endpoint terms and C1 action

In the proper inward conformal coordinate y, write
`ghat2=dy^2+a(y)^2 dt^2`, with
`t=B sqrt(lambda)exp(alpha_infinity)U`. Its expansion is

```
a(y)=1+y^2/(2B^2)+c3 y^3+O(y^4 log^p y),
c3=-m sqrt(lambda)/B^3.                            (11)
```

At fixed clock, a general allowed tangent in a fixed background
proper chart has half-relative components
`h=-delta F/(2F)=O(y^3)` and
`j=delta alpha+delta F/(2F)=O(y^3)`. Decompose it into a conformal
variation omega and a boundary-preserving radial diffeomorphism xi:

```
omega+xi'=h, omega+(a'/a)xi=j,
xi=a integral_0^y (h-j)/a,
xi=O(y^4), omega=O(y^3).                            (12)
```

The canonical distribution consists of the cubic `A FP(1/y)` and
its delta contact plus an L1 remainder, as proved in the charge note.
Its contact vanishes on omega(0)=0, the finite-part product with
omega is locally integrable, and the finite normal pressure times
xi vanishes at the boundary. Thus there is no canonical fixed-source
endpoint term. This does not discard the nonzero bulk determinant
variation.

Pointwise existence of a source would not suffice for differentiability
of its action. Here is the uniform argument. The proper coordinate
above is used to evaluate an invariant endpoint coefficient on each
FIXED metric, not to differentiate a moving proper-coordinate family.
For parameter estimates use the fixed compactified areal coordinate
`z=B/(sqrt(lambda)r)` throughout. In that chart the orbit metric is

```
ghat2=b(z)^2 dz^2+a(z)^2 dt^2,
a=exp(alpha-alpha_infinity)sqrt(F)/(sqrt(lambda)r),
b=sqrt(lambda)r/sqrt(F),
a=1+z^2/(2B^2)-m sqrt(lambda)z^3/B^3+O(z^4 log^p z),
b=1-z^2/(2B^2)+m sqrt(lambda)z^3/B^3+O(z^4 log^p z).
```

These coefficients are algebraic functions of the areal data, so
their parameter directions have all five weighted jets directly
from (1). A constant clock variation can instead be kept in the
tangential coefficient in the fixed U chart. Neither a parameter
dependent coordinate nor a translated distance d(g) is differentiated.
Fixed-chart distance and proper distance are uniformly comparable.

Isolate the exact cubic
heat profile and its a4 meromorphic contribution FIRST, using the
finite-cut continuation of the charge note. Subtract the common a0,
a1,a2 terms from the remaining orbit heat trace. The real C5 grouped
Brownian-bridge proof in
[the canonical response](stationary-centered-canonical-response.md),
equations (6)-(9), supplies the grouped expansion used here. The
earlier proof used an analytic SEED to obtain numerical constants;
one must not transfer that assumption to an arbitrary C5 baseline.
For the qualitative assertion use the fixed-chart finite Taylor
parametrix and Duhamel expansion for BOTH coefficients a,b. Its
principal, first-order and zero-order remainders use respectively
fifth, fourth and third real jets of the metric and its direction.
The actual mixed boundary realization is retained, including the
metric dependent Robin coefficient induced by the first-order
projector in a half-density representation. Its low-weight boundary
profiles are included in the subtraction below; it is not replaced
by a newly chosen Neumann square. Uniform ellipticity and the fixed
paired projector give uniform Gaussian/derivative bounds for this
finite parametrix. No analytic spatial extension is needed.
Constants may depend on that neighborhood, fixed p and fixed
positive temperature; no global quantitative ball or zero-temperature
gap is being inferred.

For the direction after cubic subtraction, let Y denote the fixed
coordinate z, comparable to inward proper distance. On Y>=sqrt(t)
the degree-four local model has fifth-order remainder bounded by
`C Y^-1(1+|log Y|)^N`, giving
`C t^(3/2)(1+|log Y|)^N/Y` by the Duhamel counting below.
On Y<=sqrt(t), the weight-four principal, first and zero-order
remainders have sizes t^2, t^(3/2), t respectively. Their insertion
costs are respectively 1, sqrt(t), t, including the Duhamel time
integral. After the orbit normalization t^-1 each diagonal bound
is `C t(1+|log t|)^N`. The local model matches the actual boundary
jets through degree three, so their projector-induced Robin
coefficients agree; no unmatched boundary-condition remainder is
hidden in this estimate.
Integrating Y gives

```
||partial_parameter H_orbit,sub(t)||_L1
       <= C t^(3/2)(1+|log t|)^(N+1).              (13)
```

The crossing error `C Y^4/t exp[-Y^2/(8t)]` has the same integrated
power: put Y=sqrt(t)v. Fixed-chart exits are exponentially small.
On the compact middle and regular bolt use fixed regular charts,
not a differentiated moving normalizer. At an observation point
replace the metric by its degree-four Taylor polynomial, continued
to a fixed positive cap. The principal, first-order and zero-order
coefficients of the two Dirac squares differ respectively by
O(|x|^5), O(|x|^4), O(|x|^3), also after a parameter derivative.
The common a0,a1,a2 use jets through order four and therefore agree.
The fixed-chart finite parametrix for a uniformly elliptic C5
neighborhood gives Gaussian kernels and the first two derivative
bounds with costs t^(-1/2), t^-1. Duhamel's formula and bridge
moments then give, including the orbit normalization t^-1 and
time integration t, the three powers
`t^-1 t t^(5/2) t^-1`,
`t^-1 t t^2 t^(-1/2)` and `t^-1 t t^(3/2)`.
Each is t^(3/2). Split the time convolution at its midpoint and
place differentiated kernels on the longer-time factor; integration
by parts on the other factor is allowed by the stated vanishing
orders, which still exceed the derivative costs. Thus endpoint
time integrals are finite. The same finite parametrix and Duhamel
argument after one parameter derivative uses only the corresponding
C5 direction jets. The resulting constants are locally uniform.
This proves the compact/bolt counterpart of (13) without a sixth
metric derivative, and supplies the same fixed-chart local proof in
the outer collar. The expansion is differentiated before the
Gaussian/bridge bounds are applied; no moving indicator is used.

For the ACTION the angular heat trace is O(t^-1), whereas the
mass-Green trace in the earlier source argument was O(t^-2).
Multiplication by the action measure dt/t turns (13) into
`C t^-1/2(1+|log t|)^(N+1) dt`, which is integrable at zero.
The undifferentiated remainder obeys the same bound. The isolated
cubic/a4 contribution is explicitly meromorphic and C1 in c3 and
the clock; take its prescribed finite part. The nonzero angular
modes and the LLL are both included in the actual paired square.
At large t its locally uniform positive gap controls the trace and
its parameter derivative.

The action-level finite-cut limit needs more than fixed-smearing
continuity of the source distribution. At a fixed areal cutoff use
the fixed coordinate Y=d and translated chart u=Y-d. Apply the
argument to BOTH a(d+u), b(d+u) and the induced mixed-square
boundary coefficients. Constant normalizations by their positive
boundary values preserve the C1 topology. Before estimating a remainder, extract the
FULL common finite-cut heat polynomial, including the mixed-square
boundary half-integer terms, not just the cubic term surviving at
d=0. Equivalently in the reflected collar extract all boundary
profiles of weighted jet degree at most three, together with the
interior a0,a1,a2. These profiles follow by the finite Taylor/Duhamel
expansion; their coefficients are finite polynomials in the boundary
jets and in the fixed projector data. Their Mellin poles and finite
parts at fixed mu0 are C1 functions of those jets. Use one fixed
proper collar cutoff and one fixed small-time splitting scale,
neither of which depends on d. There is no inserted log(d) scale.

In this common identification both linear jets are O(d), the
quadratic jets tend to +1/(2B^2) and -1/(2B^2), and the cubic
jets tend to c3 and -c3. The projector-induced Robin data and their
jet limits follow algebraically from the same coefficients. These
convergences hold also after a parameter derivative. The explicit
low-weight profiles therefore converge with their renormalized
finite parts; factors with a vanishing linear jet disappear.
Their un-subtracted divergent time integrals are NEVER dominated
merely because their coefficients tend to zero.

After this complete subtraction the translated Taylor remainders through
degree three satisfy
`|partial_u^j remainder|<=C u^(4-j)(1+|log u|)^p`, j<=5,
uniformly in small d and in the chosen C1 family. For j<=3 use the
integral Taylor remainder; for j=4 use
`1+|log(d+u)|<=1+|log u|` in a collar with d+u<1; for j=5
also use `(d+u)^-1<=u^-1`. Lower-jet nonlinear products are
already in the extracted finite profiles or have weight at least
four. Applying the fixed-chart parametrix/Duhamel estimate to the remainder
again gives (13) uniformly in d. The low jets stay in a fixed
elliptic compact set, so the constants do not diverge when the
linear jet vanishes. Away from the collar, fixed charts and graph
perturbation give convergence and a common large-time gap. Hence
dominated differentiation is valid AFTER the explicit finite-cut
meromorphic terms have been removed. This proves the claimed C1
limit of the action, including the moving cut.

This proves local C1 dependence of the canonical action, not merely
a formal variation of a determinant on a single metric. The outer
coefficients need only weighted open-collar C5; `y^4 log y` is
allowed and no closed-boundary C5 or C6 direction is asserted.

## Finiteness of the complete action and clock response

No mass-dependent logarithmic action divergence is hidden in (5).
In its non-sigma braces, the four possible lambda m r terms are

```
Rbar T/30: +4lambda m r/5,
-T^2/10:  -4lambda m r/5,
11TF/90:  -11lambda m r/45,
-11F^2/180:+11lambda m r/45.                        (14)
```

They cancel. At alpha_infinity=0 the bulk density is therefore
`-13lambda^2 r^2/(240pi)-lambda/(60pi)+O(r^-2 L^N)`;
the sigma term is integrable. All fast remainders in (1) change
this only by integrable terms. The boundary action has just the
known R^3 and R divergences, a finite mass term and vanishing
remainders. The Weyl bulk-plus-boundary R^3 coefficient is
`653lambda^2/(15120pi)`; vacuum centering of the Einstein action
adds `231lambda^2/(15120pi)`. Their sum is
`221lambda^2/(3780pi)`, canceled by (9). The R coefficient
`67lambda/(3780pi)` is likewise canceled by the combined expansion
of both counterterms in (9). Each coefficient carries the same
exp(alpha_infinity). The Einstein part at lambda has its usual
finite limit by (4), and the canonical part is finite by (13).
The total action and its first variation thus have finite limits.

For an open clock variation, retain the complete limiting
distributions until they have been added. In the notation of the
charge note, `A=-exp(alpha_infinity)lambda m/(5pi)`. The canonical
`+A FP_ell(1/y)` and local `-A FP_ell(1/y)` cancel before a
shrinking-collar test is taken. Their delta coefficient is

```
-exp(alpha_infinity)lambda m/(5pi)
       [log(mu0/sqrt(lambda))+gamma_E+151/180].      (15)
```

The leading vacuum terms cancel using the exact lambda root before
the same limit. A separately divergent bulk term must not be
discarded by claiming it is L1. Adding (4) gives (2). The L1 state
remainder has no new endpoint charge; it remains part of the bulk
equations. Identical exteriors double the per-exterior result.

An admissible additional boundary functional must have zero
variation in every normalizable direction for EACH nearby
alpha_infinity, not just at alpha_infinity=0. Consequently it is
independent of m and the fast remainders on each connected local
fiber. Differentiating in the clock cannot change the mass slope
in (2). A term proportional to alpha_infinity m fails this test
away from the single physical clock slice. This establishes the
restricted source-offset ambiguity without assuming a completion
for all boundary metrics.

## Relation to the physical boundary and the coupled solver

The finite-boundary distinction is substantive. In metric F(R)
gravity a generalized GHY term generally leaves a curvature
variation at a finite endpoint; specifying only the induced metric
does not remove it. See [Dyer--Hinterbichler, section 5.3](https://arxiv.org/html/0809.4033v2)
and [Guarnizo--Castaneda--Tejeiro, section 3.1](https://arxiv.org/html/1002.0617v4).
Those results do not themselves establish the nonlocal Dirac
theorem here. The present proof instead uses the explicit action
and the vanishing weights of the restricted asymptotic tangents.

The resulting solver problem is now precise: seek an extremum of
this action in the fast phase at the fixed physical clock, flux,
lambda and prepared temperature. The exact trace still requires

```
C[g]=integral_(r_h)^infinity exp(alpha)r^2 phi_g S[g] dr=0,
(a_infinity+epsilon Box_g)phi_g=0, phi_g regular at the bolt,
epsilon=G/(60pi), a_infinity=1+11G lambda/(180pi).    (16)
```

See [the range theorem](stationary-trace-range-condition.md) for
the source S[g], endpoint normalization and fixed-metric inverse.
Neither differentiability nor the fast phase makes (16) automatic
on an off-shell candidate. The small complete residual proves no
distance to an extremum. The next substantive decision is whether
an allowed regular horizon datum in the FULL coupled equations
can enforce (16) while holding the physical sources and temperature
fixed, followed by a controlled inverse on the compatible
complement. Changing the preparation to solve a scalar shooting
equation would answer a different question. No transversality,
physical instability, nonexistence or physical gate closure is
claimed by this variational theorem.

## Focused verification contract

The ordinary heat/domain/limit argument and restricted scope above
precede execution. A focused checker will independently differentiate
the Einstein and local radial densities, test the exact endpoint
cancellations, derive the bulk action divergences, and check the
algebra of the boundary weights and integrated heat exponents.
Symbolic residuals must be exactly zero; rational inequalities must
be true. No numerical tolerance, sampled uniformity or new physical
run is selected. Reuse the accepted cubic distribution and charge
coefficients; do not rerun their historical checkers. Hash exact
UTF8/LF inputs before and after, retain every execution outcome and
failed input bytes, and record the actual Python -B command and
versions. The ordinary analytic proof is not certified by symbolic
checks. Evidence belongs in stationary-fast-phase-variation.json/.log.

Remaining to-do list: verify this restricted stationary theorem;
justify the physical boundary/branch choice and the compatible
coupled solver at fixed preparation; establish the exact exterior
and physical state and transfer S138; finish S074-S079 and S080.
