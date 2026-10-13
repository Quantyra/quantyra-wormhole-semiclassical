# Compact charge and rotation modes in the formation transfer

Subsequent 2026-10-13 correction to the physical transfer:
[fixed magnetic boundary-field response](angular-boundary-response.md).
The Kerr-Newman inertia used below is a family coefficient that changes
the boundary electric field. It is not the fixed-source susceptibility.
The conditional rotor calculations below are preserved, but their
single-species failure and many-species selection do not apply to the
newly specified fixed-field ensemble. Its coupled static bound permits
a single-species candidate; frequency-dependent matching remains open.

S072, 2026-10-13. **The free JT initial covariance cannot be transferred
by simply freezing the electromagnetic and angular phases of a magnetic
black hole.** A compact-mode calculation gives an exponentially small
angular correlation for one four-dimensional fermion species in the
weakly dressed, separated-LLL regime. A many-flavor AdS scaling avoids
this particular state mismatch. It remains a candidate, with its coupled
switching kernel and angular response unverified.

The preceding [capture](switching-radiation-capture.md) and
[bulk geometry](driven-bulk-geometry.md) results remain valid in their
declared free-Majorana JT model. This note tests a missing physical input
to that model. It does not revoke those ordinary theorems or promote
them to four-dimensional results.

## What is being matched

Use Einstein-Maxwell normalization `-F^2/(4g^2)` and unit electric charge
for each of `N_f` four-dimensional Dirac species. A magnetic flux integer
`q>0` gives `n=N_f q` complex LLL channels, or `N=2N_f q` real fields.
For extremal radius `b` and AdS4 radius `L`, define

```
X=b^2/L^2, C=1+6X, L2^2=b^2/C,
b^2(1+3X)=pi G q^2/g^2, phi_r=b^3/(G C),
T=N eta_f z_0/(2 pi phi_r), 1/4<=z_0<=3/4.                 (1)
```

The last relation is the temperature of the initial TFD trajectory in
the previous opening proof. The compact modes are physical boundary
charge and rotation variables, not the already excluded independent
proper-mouth clock.

[Iliesiu and Turiaci, sections 3.1--3.2](https://arxiv.org/html/2003.02860v3)
derive the compact group particle alongside the Schwarzian and distinguish
fixed-charge from fixed-potential ensembles. That derivation neglects
coupling of these gauge modes to additional light charged matter. We
therefore derive a spherical LLL screening bound below and keep the
angular dressing as an explicit uncompleted condition. We do not borrow
their electric susceptibility along an electrically charged background:
here electric charge varies at fixed **magnetic** flux.

A TFD projected to fixed individual electric charges has zero expectation
of a charged transfer operator. A TFD projected to individual angular
singlets has zero expectation of a nontrivial angular tensor. Even the
more favorable ensemble with zero chemical potentials and allowed charge
and angular fluctuations must include their quantum phase factors.
The factorization into neutral, U(1) and SU(2) factors below is the leading
compact-mode approximation. Matching the actual monopole charge-spin
lattice, neutral matrix elements and matter-dressed kernel remains required.

## Electric susceptibility, including a static spherical screening bound

Adding a small electric integer `Q_e` changes the geometric RN charge
squared by `G g^2 Q_e^2/(4 pi)`. Extremality makes the first-order radius
shift drop out of the mass variation. Thus

```
M_ext(Q_e)-M_ext(0)=E_C Q_e^2+O(Q_e^4),
E_C=g^2/(8 pi b), K_e=1/(2E_C)=4 pi b/g^2.                 (2)
```

This is also obtained by the electrostatic solution
`A_t=mu(1-b/r)` and the boundary flux. It is a susceptibility in the
electric direction, not the radial charge variation of a purely electric
extremal family. For the bare compact mode its thermal parameter is

```
x_e=E_C/T=pi q/[8 N_f eta_f z_0 C(1+3X)].                 (3)
```

Small `g` alone does not force this ratio to be small at fixed `X`.

In the static spherical massless-LLL Gaussian polarization problem, with
regular horizon potential zero and boundary potential `mu`, the positive
electrostatic functional is

```
F[A]=(1/2) integral_b^infinity
     [(4 pi r^2/g^2)(A_r)^2 + (n/(pi f(r)))A^2] dr.
```

The coefficient `n/pi` is the massless Dirac polarization; it is the
static form of the bosonized charge mode in
[MMP appendix A](https://arxiv.org/html/1807.04726v3). This functional
fixes the boundary conditions and the Gaussian screening model being
checked. It does not include nonspherical photons, gravitational mixing
or a real-time boundary interaction.

Its minimum defines `F_min=K_eff mu^2/2`. Dropping the positive polarization
term gives `K_eff>=K_e`. Using the Maxwell solution as a trial function
gives the upper bound. The exact extremal factorization is

```
f(r)=(r-b)^2[r^2+2br+L^2+3b^2]/(L^2 r^2),
D(X)=atan(sqrt(1+2X)/(2sqrt(X)))/sqrt(X(1+2X)),
epsilon_pol=(n g^2/(4 pi^2))D(X),
1<=K_eff/K_e<=1+epsilon_pol,
E_C/(1+epsilon_pol)<=E_C,eff<=E_C.                         (4)
```

To obtain `D`, integrate
`(1-b/r)^2/f=L^2/[r^2+2br+L^2+3b^2]` from `b` to infinity and divide
by `b`. The integral is finite even at the extremal horizon. For a
nonextremal outer radius `r_+`, the same upper bound with `b` replaced
by `r_+` follows because its `f` is at least the extremal `f` with that
same outer radius. No uniform flat-space bound is inferred:
`D(X)` diverges as `X->0`. At large `X`,
`D(X)~atan(1/sqrt(2))/(sqrt(2) X)`.

The proof controls the complete static spherical Gaussian polarization,
not just its first perturbative insertion. It does not determine the
time-dependent compact-phase influence on the switching stress.

## Exact electric TFD phase factor and preparation cost

For a compact rotor with single-side energies `E_C,eff m^2`, `m in Z`,
put `x=beta E_C,eff`. Its correlated two-sided TFD has amplitudes
proportional to `exp(-x m^2/2)` in the states `|m,-m>`. The expectation
of one charge transfer, relative to a phase held fixed, is

```
C_e(x) = sum_m exp[-x(m^2+(m+1)^2)/2] / sum_m exp(-x m^2)
       = exp(-x/4) theta_2(0,exp(-x))/theta_3(0,exp(-x)).    (5)
```

For `x>=1`, summing the numerator's tails gives
`C_e<=2 exp(-x/2)/(1-exp(-2x))`. Poisson summation gives, for `0<x<=1`,

```
C_e>=1-x/4-4R/(1-R), R=exp(-pi^2/x).
```

Thus a thermal electric phase can be almost coherent or almost fully
suppressed, depending on the *charging energy divided by temperature*.

An engineered phase packet has a cost too. If its transfer expectation
is real and at least `c`, the commutator of the integer charge with
`sin(alpha)` and the uncertainty inequality give
`<m^2>>=c^2/[4(1-c^2)]`. The pair's charging energy is consequently at
least `E_C,eff c^2/[2(1-c^2)]`. Compared to the earlier physical capture
margin `E_margin=(3/32)N^2 eta_f^2/phi_r`, this implies

```
E_phase/E_margin >=
 c^2/[6 N_f^2 eta_f^2 C(1+3X)(1-c^2)(1+epsilon_pol)].      (6)
```

For example, with one species, `X<=1`, `eta_f<=1/2000`,
`epsilon_pol<=0.01` and `c>=0.9`, this exceeds `10^5`.
Such a phase packet cannot be counted as a negligible preparation error
within the existing capture budget. This does not forbid a different
state-preparation protocol with a separately accounted energy supply.

## Angular phase: a stronger fixed-flavor transfer test

An LLL fermion transforms with angular spin `j=(q-1)/2`, not spin of
order one; the monopole multiplet has dimension `q`. The first higher
Landau eigenvalue, using the retained
[magnetic sphere spectrum](magnetic-product-vacuum.md), satisfies

```
mu_LL^2=m_1^2 L2^2=(q+1)/C.                               (7)
```

At fixed physical magnetic charge, expanding extremal Kerr-Newman-AdS
gives

```
M_ext(J)-M_ext(0)=J^2/(2I_0)+O(J^4),
I_0=b^3(1+2X)/[G(1+X)].                                  (8)
```

This follows independently from `J=a M_ext` and
`Omega=a(1+X)/b^2` at small rotation. The checker also expands the
extremality equation at fixed physical charge, including the `Xi` factors.
It is not sufficient to use the uncharged or fixed metric-parameter
susceptibility. The [KN conventions](https://arxiv.org/html/2003.02860v3)
are reviewed in appendix A.1 of the cited source.

Use the spin cover SU(2) to allow both integer and half-integer sectors.
For the compact rotor define `a=beta/(2I)` and

```
K_t(g)=sum_l (2l+1) exp[-t l(l+1)] chi_l(g),
Z(a)=K_a(1)=sum_l (2l+1)^2 exp[-a l(l+1)].
```

The relative-group TFD wavefunction is `K_(a/2)/sqrt(Z(a))`.
Its angular factor for a contracted spin-j transfer is exactly

```
C_j(a) = integral dg K_(a/2)(g)^2 chi_j(g)/[(2j+1)Z(a)]
       = sum_(l,l') (2l+1)(2l'+1) exp[-a(C_l+C_l')/2]
         N_(l,l')^j / [(2j+1)Z(a)],                       (9)
C_l=l(l+1).
```

Here `N` is the SU(2) fusion multiplicity (zero or one). All terms in
the second expression are nonnegative. The triangle rule implies
`C_l+C_l'>=C_j/2`. Splitting the exponential into two equal factors
and using `|chi_j|<=2j+1` therefore proves

```
0<=C_j(a)<=exp[-a j(j+1)/8] Z(a/2)/Z(a)
          <=12 exp[-a j(j+1)/8],  a>0.                   (10)
```

For completeness the heat-trace ratio bound is uniform. Write
`Z(a)=exp(a/4) sum_(n>=1)n^2 exp(-a n^2/4)`. For `a>=1`,
`Z(a/2)/Z(a)<=Z(1/2)<12`. For `a<=1`, Poisson summation gives
`Z(a)=2sqrt(pi) exp(a/4)a^(-3/2)B(a)`, where
`B(a)=1+2sum_(k>=1)(1-8pi^2 k^2/a)exp(-4pi^2 k^2/a)`.
All corrections are negative and their magnitude is below `10^-12`:
bound the tail by `144 r(1+r)/(1-r)^3`, `r=exp(-36)`.
Thus the ratio is below `3` on this interval. These loose constants
suffice for (10).

With the undressed inertia, the exponent is

```
chi_J := beta j(j+1)/(2I_0)
 = pi(q-1) mu_LL^2 (1+X)/[8 N_f q eta_f z_0(1+2X)].       (11)
```

Suppose `N_f=1`, `q>=100`, `eta_f<=1/2000`, `z_0<=3/4`,
`mu_LL^2>=1`, and matching changes the inertia by no more than one
percent upward. Equations (10)--(11) then imply
`C_j<1.66e-27`. Merely enlarging the AdS black hole cannot remove this
angular suppression while retaining the stated Landau separation.
The original phase-held-fixed cross covariance is therefore not a
small-error approximation to this thermal compact-mode state.

This is a **conditional transfer failure** in the weakly dressed compact
model. We have not bounded the actual angular inertia correction by one
percent. Strong screening, a different ensemble, neutral composite
couplings or a different preparation can change the question. In
particular, the spherical Maxwell bound (4) is not an angular inertia
bound: frame dragging mixes with the magnetic gauge field and the
fermions' large angular representation. It must be matched separately.

## A surviving many-flavor scaling, and what it does not prove

The next selected candidate keeps the free opening angle small while
increasing the number of **four-dimensional species**, rather than
identifying magnetic degeneracy alone with all required large-N limits.
For integer `k>=2000`, take, in Planck units,

```
eta_f=k^-1, g=k^-4, L/l_p=k^5, q=k^3, N_f=k^3,
z_0=1/2, X(1+3X)=pi k^4, X~sqrt(pi/3) k^2.               (12)
```

Then `N=2k^6`, `N eta_f=2k^5`, and direct substitution yields

| Quantity | Large-k behavior |
| --- | --- |
| Electric phase parameter `E_C/T` | `O(k^-3)` |
| Angular exponent `chi_J` | `O(k^-1)` |
| First higher Landau mass squared in AdS2 units | `O(k)` |
| Proven spherical screening bound `epsilon_pol` | `O(k^-4)` |
| Spherical Schwinger mass squared in AdS2 units | `O(k^-4)` |
| Angular polarization scale proxy `N_f g^2 q` | `k^-2` |
| `N G/b^2` | `O(k^-6)` |
| `N eta_f/Phi_ext`, at classical leading order | `O(k^-7)` |
| Initial `T b` | `O(k^-5)` |

The small-a heat kernel is a Gaussian near the group identity, with
mean squared class angle `3a/8+O(a^2)`. Since
`1-chi_j(alpha)/(2j+1)<=2j(j+1)alpha^2/3`, it follows that
`C_j->1` whenever `a->0` and `a j(j+1)->0`. The image corrections are
exponentially small and the normalized character is bounded by one.
Together with (5), both leading compact-mode coherence factors approach
one along (12). The angular proxy in the table is a scale diagnostic,
**not a proved dressed-inertia bound**.

The charge dependence of the extremal entropy is also distinguishable
from a constant rotor normalization. At fixed magnetic flux,
`partial_(Q_e^2) S_ext=g^2/(4C)`. Its ratio to `E_C/T` is `O(k^-7)`
along (12), so that particular correction does not spoil the displayed
leading electric limit. Other density-of-states and matter corrections
are not certified by this observation.

This is an explicit extension of the matter content. For example,
`k=2000` already means eight billion massless Dirac species. It is a
theoretical scaling test, not a proposed material or a Standard Model
construction. It proves compatibility of the listed necessary scales,
not a four-dimensional solution. Large AdS holes also have many light
angular Maxwell and gravitational modes; a Landau gap for fermions does
not remove those modes. Their response, the joint charge-spin state,
and the fully gauge-invariant switching stress are the next hard task.

## Verification and decision

[Checker](../checks/compact_phase_transfer.py),
[receipt](../evidence/compact-phase-transfer.json), and
[source scopes](../evidence/compact-phase-source-review.json).
The checker verifies RN/KN coefficients and the Landau/thermal comparison,
then tests the capacitance bound in twelve radial boundary problems with
horizon-cutoff refinement. Six electric TFD sums agree with theta functions
at 70-digit working precision. Twelve angular fusion sums agree with
independent Haar-character integrals. A warm-kernel refinement and three
exact-flux scaling examples are retained. These tests verify the stated
compact/linear-response models, not their full physical transfer.

The result rejects silently dropping compact phases in the fixed-flavor,
weak-dressing transfer. Work continues on the many-flavor candidate, with
the angular inertia and coupled phase/fermion quench as the next decision.
No cold-state universal no-go, four-dimensional stability result, novelty
acceptance, publication, outreach or Lean certificate is claimed.

Remaining to-do list: derive the joint dressed charge/angular response and
switching kernel; bound exterior, finite-cutoff and other angular sectors;
establish initial-pair preparation and finite-energy stability/transport.
