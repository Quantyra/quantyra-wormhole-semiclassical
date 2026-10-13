# Radiation-wall occupations and Euclidean continuation

S070, 2026-10-12. **The spherical collective action is a consistent
classical restriction of the ray model, but not its complete thermal
fluctuation integral.** Exchange statistics and the gravitational ensemble
prevent promotion of that action's radial index to a physical creation
rate. This completes the specified probe comparison; the full formation
saddle's physical mode count remains unresolved.

The [preceding response calculation](ads-wall-response.md) defines the
material and collisionless comparators. Here we continue the latter at
fixed angular momentum. Everything through the exchange-sector calculation
uses a fixed background, neutral noninteracting bosons and a leading
geometrical-optics action. This is a chosen comparison, not the microscopic
state of the BFR wall. Its statistics, particle number and occupations have
not been derived from the bulk fields.

## 1. A specified canonical probe problem

Use units hbar=1 and metric

    ds² = -f(r) dt² + dr²/f(r) + r² dOmega².

Choose a regulated radial interval containing the circular orbit R with
f>0 and reflecting endpoints. This makes a positive self-adjoint radial
Hamiltonian and its thermal trace meaningful after specifying its quantum
ordering. We work only with its principal classical symbol

    H_J(r,p) = sqrt(f² p² + f J²/r²),     J>0.

An exact quantum field reduction can add potentials, measures and
polarization constraints. Those are not supplied by this symbol. In a
large angular-momentum sector j, J is the leading value sqrt(j(j+1));
angular degeneracy is d=2j+1 for a scalar comparison. We do not identify
that scalar degeneracy with an exact photon or fermion spectrum.

Time slicing Tr exp(-beta H_J) gives the leading phase-space exponent

    S_E = integral [H_J(r,p) - i p r'] d tau.

Writing Q=sqrt(f+r'²/f), stationarity in p on the branch continued from
p=0 gives p*=i J r'/(r f Q), H(p*)=J f/(r Q). Substitution yields

    L_E,J = (J/r) sqrt(f+r'²/f).

This is exact elimination of the **classical symbol**, not an exact quantum
configuration-space path-integral measure. The complex momentum saddle is
why simply Wick rotating the angular advection operator is insufficient.

Let u=sqrt(f)/r, u'(R)=0 and u''(R)<0. Expanding r=R+chi gives

    H = J u(R) + ap p²/2 + br chi²/2 + O(3),
    ap = R f(R)^(3/2)/J > 0,       br = J u''(R),
    S_E^(2) = (1/2ap) integral [chi'² - gamma² chi²] d tau,
    gamma² = -ap br = f(R)²/R² - f(R) f''(R)/2.

The positive kinetic term also follows by Gaussian momentum integration.
The RN-AdS values and the source F.8 discrepancy are retained from the
preceding note. We use the independently derived gamma, not F.8.

For N rays with the same r(t), the action has J replaced by sum J_i.
Their classical radial equations agree independently of J_i, so this is
an invariant classical restriction. It excludes independent radial
motions, however. Consistency of a restriction does not make it a complete
semiclassical fluctuation space.

## 2. Boson exchange changes the radial periods

For fixed N in one j multiplet, tracing all m states, write h=h_J tensor
1_d and K=exp(-beta h). The symmetrizing projector is
P_sym=(1/N!) sum_pi U_pi. Hence

    Z_N,j = (1/N!) sum_pi product_{cycles c of pi} Tr K^c
          = (1/N!) sum_pi product_c [d Tr_rad exp(-c beta h_J)].

Proof: in a product eigenbasis, each permutation cycle identifies all
indices around that cycle, leaving one trace of K^c. Equivalently Z_N,j
is the coefficient of z^N in det(1-zK)^(-1). These are standard canonical
boson identities; see [Salmhofer, section 4.3, equations 67-70](https://link.springer.com/article/10.1007/s00220-021-04010-4).
We independently check the cycle sum against explicit Fock-state sums.

A c-cycle glues its radial paths into one path of period c beta. Its
static-orbit radial Hessian has normalized eigenvalues

    lambda_n,c = (2 pi n/(c beta))² - gamma²,   n in Z.

For real paths n=0 counts once; each positive n contributes sine and
cosine. Thus, with x_c=c beta gamma/(2 pi)>0,

    nu_c = 2 ceil(x_c) - 1.

At positive integer x_c there are also two zero modes, excluded from nu_c.
For a specified permutation saddle the radial count is sum_c nu_c.
For example at N=4 and beta gamma/(2 pi)=1/5, the identity sector has
four radial negative directions; the single four-cycle has one. At ratio
2/5 those counts are four and three. We cannot report either number as
the index of the whole quantum gas. Sectors have different path gluing,
weights and possible integration contours; their indices are not added
across the permutation sum.

If individual angular occupations N_jm are fixed, only exchanges within
each identical angular state occur: use a product of such traces with d=1.
If only N_j is fixed, sum over m occupations as above. Fixing all radial
energy occupations would specify a different state again. The boundary
conditions and the allowed state variations must therefore be declared.
Fermi statistics require permutation signs; we have not analyzed that
case. Physical radiation also need not conserve the chosen N.

These are **local leading-action radial indices**, not determinants. An
ordinary thermal trace of a self-adjoint Hamiltonian is real and positive.
A barrier-top saddle with a negative Hessian does not by itself define
an imaginary part or a decay rate. A metastable-state prescription and
its integration contour are additional information.

## 3. Why fixing boundary masses does not finish the count

Our particle Hamiltonian already demonstrates the danger of a naive
linear constraint. At the static orbit, dH=0 for radial (chi,p) variations,
but the fixed-energy equation to second order is

    ap p² + br chi² = 0.

Its level set at the saddle is not a regular hypersurface with a single
linear normal. Imposing dH=0 and then subtracting one mode is unjustified.
This observation does not determine the gravitational energy surface.

[Marolf and Santos, sections 2.1-2.2](https://arxiv.org/html/2202.12360v1)
derive a static microcanonical gravitational integral by transforming the
canonical one. The fixed quantity is an off-shell Hamiltonian containing
bulk constraint terms as well as boundary energy. They retain the
second-order constraint when comparing Hessians. Their Schwarzschild
analysis is not a radiation-wall theorem; we import the distinction,
not their stability result.

For the BFR two-boundary geometry we would need the full off-shell action,
the relation between the two time lapses, and the variables fixed by the
chosen preparation. Holding two boundary temperatures fixed differs from
fixing two separate energies or only their sum. The imposed inter-boundary
coupling also makes independent conservation assumptions nontrivial.
With a microscopic matter state selected, the metric, Maxwell and matter
constraints must be reduced together and the conformal integration
contour specified. None of these follows from the probe action.

The original [BFR sections 5 and appendix F](https://arxiv.org/html/2308.00871v2)
provide the ensemble-dependent background actions and the radial caveat.
They do not provide the microscopic occupation and full nonspherical
fluctuation problem used in this comparison. We have not refuted their
formal action difference or computed its physical rate.

## 4. Verification and research decision

[Checker](../checks/ads_wall_ensemble.py),
[receipt](../evidence/ads-wall-ensemble.json), and
[source scopes](../evidence/formation-transfer-source-review.json).
Symbolic checks verify the complex momentum elimination and kinetic
coefficient. Independent canonical sums use 24 cases with N=1..6,
degeneracies 1 and 3, and two temperatures. Thirty-six finite differences
of the nonlinear action test twelve constant/oscillating perturbations
on two backgrounds; step reduction checks convergence. Exact rational
thresholds test negative and zero mode handling.

S070 is **parked with its physical mode count unresolved**, not completed
as a proof of instability or impossibility. The unsupported promotion of
a spherical wall rate to a controlled microscopic creation process is
rejected for this research program. The material conditional obstruction
and the collisionless comparison remain useful. Reopening requires an
explicit bulk-field state/action and admissible variations, followed by
the constrained gravitational calculation.

We now execute the distinct S071 real-time formation transfer gate in
[formation-entropy-transfer.md](formation-entropy-transfer.md). This
follows the user's instruction to record an inconclusive avenue and
continue; it does not reduce the broad stable-wormhole goal.

Remaining to-do list: for S070, specify the microscopic state and solve
the full constrained fluctuation problem before claiming a rate; active
work proceeds with the S071 entropy, boundary-interaction and preparation
tests, then formation and finite-disturbance traversability.
