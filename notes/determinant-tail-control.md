# Global inner determinant bound and the gauge marginal

2026-10-11, E003/S053. The preceding [nonlinear-current estimate](nonlinear-current-determinant.md)
controlled small sources but left the determinant-weighted large-field
tail open. For the specified vector-coupled, antiperiodic free-fermion
return model, the full covariantly normalized determinant satisfies
`0<=Delta[A]<=1` at arbitrary amplitude. This gives an integrable weight
for any finite smooth Gaussian photon regulator and closes that tail
under explicit small-noise bounds. It retains the unequal fermion returns.

This is a statement about a regulated Euclidean gauge marginal in the
defined inner-source model. The physical gauge covariance, exterior
sources, nonlinear gauge/geometry action, photon-cutoff removal, local
stress and Lorentzian state reconstruction are not supplied. The general
two-dimensional determinant upper bound is established prior work; no
physical endpoint or consequential novelty is accepted.

## The determinant and its normalization

Keep q complex massless Dirac channels with spatial lengths L_i=T+D_i,
antiperiodic spatial reference and common source interval I of length T.
The unperturbed spatial Hamiltonians are h_R=-i partial_s and h_L=-h_R.
Take real Euclidean time and space components A_tau,A_s that are
Hermitian q by q matrices on I, zero outside it. Initially use smooth
sources with support away from its endpoints. Time is a circle of length
beta with antiperiodic fermions, or the zero-temperature limit with the
source confined to a fixed finite time interval. Spatial H1 profiles
that vanish at the endpoints follow by smooth approximation.

Normal order the free Fock Hamiltonian H0 and vector currents in the
neutral spatial antiperiodic sea. Thus the reference ground energy is
zero and Tr exp(-beta H0) is finite at fixed q and beta>0. Take
the thermal trace over the full Fock space at zero chemical potential;
this is not a projection onto a fixed charge or compact gauge sector.
Write
rho=psi_R^dagger psi_R+psi_L^dagger psi_L and
j=psi_R^dagger psi_R-psi_L^dagger psi_L, both normal ordered. The
covariant external-field Hamiltonian is

```text
H(A_s)=H0+integral_I Tr(A_s j) ds
            +(1/(2pi)) integral_I Tr(A_s^2) ds.             (1)
```

Tr here contracts the color indices of the source and current. The
temporal component enters the Euclidean evolution as i Q(A_tau), where
Q(A_tau)=integral_I Tr(A_tau rho) ds is self-adjoint. Define Delta_beta
as the normalized thermal trace of this evolution, or equivalently its
time-sliced limit. The trace at zero source divides it.

The last term in (1) is essential. It is the vector-current contact of
[the polarization calculation](propagating-photon-polarization.md), not
a freely chosen new interaction or a photon mass. Spatial covariant
point splitting multiplies each chiral two-point function by the short
Wilson line. Its singular term is plus or minus I/(2pi i(s-s'));
expanding the Wilson line changes the spatial current by A_s/pi after
adding the two chiralities. Thus variation of the covariant Hamiltonian
is j+A_s/pi. Integrating this variation from zero source fixes exactly
the c-number in (1). Unequal returns change the regular part of the
two-point function, not this local singular coefficient. There is no
A_tau^2 term in this Hamiltonian prescription: temporal vector charges
generate unitary gauge rotations. This prescription gives the same
vector Ward contact as the retained Euclidean polarization.

The normal-ordered current forms exist beyond the small-source domain.
For a finite fermion regulator the pair estimate in
[the driven-sea proof](driven-fermion-sea.md) bounds the pair contribution
by 2 sqrt(a(A) E), with a(A) controlled by integral Tr A^2. For the
diagonal term, Pauli exclusion and an energy split at M give
N_particles+N_holes <= E/M+N_low(M), where, for both chiralities,
N_low(M)<=2M sum_i L_i/pi+4q. Hence its expectation is bounded by
||A||_infinity [E/M+N_low(M)]. Choosing M large and completing the pair
square makes either current infinitesimally form bounded relative to
H0. These bounds are uniform in the fermion regulator at fixed q.

For smooth spatial A_tau, its positive/negative-energy off-diagonal
block is Hilbert-Schmidt: its matrix elements are Fourier transforms
at the pair energies, and the pair counting function grows at most
quadratically while the smooth transforms decay. Multiplication by
exp(-it A_tau) therefore has an implementable vector CAR action. Its
continuous unitary group, with the neutral normal-ordering phase,
defines the self-adjoint charge Q used in the trace product. The two
chiral vector charge anomalies cancel; no chiral determinant is being
integrated on its own.

Consequently (1) defines a lower-bounded self-adjoint form on the free
energy form domain. Smooth bounded time dependence gives the common
sectorial form H(A_s)+iQ(A_tau) and its Euclidean evolution. In a
finite time slicing one may use the product
exp[-dt H(A_s)] exp[-i dt Q(A_tau)]. The form product converges as the
mesh tends to zero. Free energy cutoffs can be removed first: the
infinitesimal bounds give H(A_s)>=(1-epsilon)H0-C_epsilon uniformly on
a bounded source family; the finite-volume free Fock heat trace is
finite at every positive temperature. Strong form convergence and this
heat-trace domination give the trace limits used below. This construction
does not require a gap uniform in q.

## Static holonomy bound with unequal lengths

First make A_s independent of time and set A_tau=0. The spatial
first-order equation transports a color vector through I by a unitary
matrix U=P exp(-i integral_I A_s ds). At positive imaginary momentum
omega the common scalar propagation commutes with this ordered unitary.
Completing the free exterior gives the round-trip matrix
R(omega) U, where R(omega)=diag(exp(-omega L_i)). The other chirality
gives the conjugate determinant. The antiperiodic sign is a plus sign
in I+R U.

With omega_n=(2n+1)pi/beta, n>=0, the static ratio is

```text
log[Tr exp(-beta H(A_s))/Tr exp(-beta H0)]
 =4 sum_(n>=0) [log |det(I+R(omega_n) U)|
                         -sum_i log(1+exp(-omega_n L_i))]. (2)
```

One way to obtain (2) without a choice of spectral phase is to compute
the first-order resolvent on each directed circle. The returning part
is the geometric series in -R U; its norm is less than one for omega>0.
Integrating the diagonal resolvent gives the omega derivative of
log det(I+R U). The local nonreturning term is independent of U and
cancels in the relative determinant. Fix the integration constant by
the ratio tending to one at omega=infinity. Combining the opposite
chirality and the positive/negative temporal frequencies gives the
factor four. Covariant point splitting gives the Hamiltonian (1), so
this transfer determinant and the canonical trace have the same
normalization. The sum converges exponentially at every fixed beta.
For one equal-length eigenphase alpha, its beta=infinity limit gives
the previously checked vacuum shift wrap(alpha)^2/(2pi L).

For any positive diagonal R with eigenvalues r_i and any unitary U,

```text
|det(I+R U)|
 <= sum_k |Tr exterior_power_k(R U)|
 <= sum_k Tr exterior_power_k(R)
  = product_i (1+r_i).                                     (3)
```

The middle inequality is trace-norm duality with the unitary exterior
power of U. No diagonalization of U with R occurs. Equations (2)-(3)
show that every static covariant Hamiltonian has partition function
at most the reference partition function. This need not hold after
deleting the c-number in (1), changing the antiperiodic reference, or
introducing a real chemical potential.

## Time-dependent fields and positivity

For a time partition with positive steps dt_k summing to beta, Schatten
Holder gives

```text
|Tr product_k {exp[-dt_k H(A_s,k)] exp[-i dt_k Q(A_tau,k)]}|
 <= product_k ||exp[-dt_k H(A_s,k)]||_(beta/dt_k)
  = product_k [Tr exp(-beta H(A_s,k))]^(dt_k/beta)
 <= Tr exp(-beta H0).                                      (4)
```

The unitary temporal factors have norm one. Thus arbitrary temporal
noncommutation adds no amplitude restriction. The limits just described
preserve the bound.

The vector determinant is nonnegative, not merely bounded in modulus.
In a finite spatial regulator let U_R be the one-particle Euclidean
evolution of h_R+A_s+i A_tau. Its left partner is (U_R^dagger)^-1.
After the neutral-sea charge phase is removed, the product of their
Fock traces is a positive scalar times |det(I+U_R)|^2. The real spatial
contact in (1) preserves positivity. Equivalently the Euclidean Dirac
operator is anti-Hermitian and anticommutes with chirality, pairing its
nonzero eigenvalues. Zeros give zero determinant, not a sign ambiguity.
Taking the limit in (4) proves

```text
0 <= Delta_beta[A] <= 1.                                   (5)
```

For time-compact sources the beta=infinity limit exists by the free
finite-q vacuum projection outside the support, and obeys (5) as well.
The individual spatial gap is sufficient; no q-uniform gap or joint
temperature/large-q limit is used.

## From the source estimate to a gauge marginal

Use the previous normalized electric/rotational sources
A_mu=sum_(A=0,x,y,z) phi_mu,A E_A. The vacuum one-point function vanishes.
For smooth sources the local prescription changes only the terms of
degree at most two; every higher loop is absolutely controlled by the
S3 argument. Therefore near zero source

```text
log Delta[phi]=-Q2[phi]+log det_3(I+D0^-1 V_phi).             (6)
```

Q2 is the complete covariantly normalized quadratic action of these
free unequal-return fermions. Equation (5) implies Q2>=0. At zero
temperature the existing smeared current bound, continued to imaginary
frequency and integrated by Plancherel, gives

```text
Q2 <= (B/2) integral |phi|^2,
B=8 A_T+1/pi,
A_T=4pi max(2/pi^2, (2/pi)(1/(2pi)+1/T)) (1+T^2/4).        (7)
```

The source is supported on the fixed interval I. The zero-frequency
point follows by the same positive-energy domination. Formula (7) uses
the normalized matrix trace and includes the identity vertex. It is a
zero-temperature bound; no thermal susceptibility bound is silently
inferred from it.

The factorization itself extends to arbitrary amplitude even though
the small-source logarithmic series does not. For a fixed smooth source,
both Delta[z phi] and det_3(I+z D0^-1 V_phi) are entire functions of the
complex scalar z. For the Hamiltonian trace, the infinitesimal current
form bounds give heat-trace domination uniformly on compact z sets;
the c-number is the polynomial z^2 times (1). Thus the finite-regulator
analytic traces have locally uniform limits. The Fredholm determinant
is entire by its S3 definition. The identity theorem applied to (6)
then gives Delta[z phi]=exp(-z^2 Q2[phi]) det_3(I+z D0^-1 V_phi)
everywhere. This is a determinant identity, including zeros, not a
globally chosen logarithm.

Let mu0 be any centered, finite-dimensional Gaussian distribution of
smooth real sources with the stated support. This includes a finite
smooth photon regulator after fixing its zero modes. By continuity at
zero source and (5), Z=E_mu0 Delta is positive and at most one. Thus

```text
dmu = (Delta/Z) dmu0                                       (8)
```

is an actual positive normalized gauge marginal at this regulator,
without deleting rare fields. It integrates the full inner free-fermion
determinant, not just its quadratic part.

Suppose E Q2<=kappa and E r<=rho<1, where r is the bound from the
preceding note. Put a=sqrt(rho), B_bad={r>a}, and

```text
eta=kappa+sqrt(rho) [1+1/(3(1-rho^(1/6)))].                 (9)
```

On the good event, 1-Delta<=Q2+|log det_3| by (5)-(6).
On the bad event it is at most one. Markov's inequality and the
previous all-orders bound give

```text
0<=1-Z<=eta,   total_variation(mu,mu0)<=1-Z<=eta,
mu(B_bad)<=sqrt(rho)/(1-eta) when eta<1.                    (10)
```

The total-variation claim uses the exact mixture
mu0=Z mu+(1-Z) mu_rejected, with rejected density (1-Delta)/(1-Z).
No bound on a cubic exponential is invoked.

This also resolves the particular weighted-tail test in the preceding
note. Let Z_G=E_mu0 exp(-Q2) and dmu_G=exp(-Q2)dmu0/Z_G. The globally
defined nonlinear weight is W=exp(Q2) Delta=det_3(I+D0^-1 V_phi).
Jensen gives Z_G>=exp(-kappa), and hence

```text
integral_B_bad W dmu_G
 = Z_G^-1 integral_B_bad Delta dmu0
 <=exp(kappa) sqrt(rho).                                  (11)
```

For bounded observables (10) gives expectation convergence. Unbounded
currents, stress, geometry variations and analytic continuation still
require additional estimates; (10) is not such an estimate.

## A complete regulated reference for the covariance test

On the fixed-R=1 flat-plane quadratic reference, keep the bare gauge
Gaussian rather than inserting the fermion polarization. For each
rotational axis its transverse covariance C0(p) from the preceding
note obeys C0(p)<=e_rot^2/p^2, with
e_rot^2=3 T_q(q^2+2)/(K(q+1)^2). The normalized electric singlet has
e_0^2=g_c^2 q/(4pi): integrating the microscopic Maxwell term over the
unit sphere gives coefficient 4pi/g_c^2, and A_charge=phi_0/sqrt(q).
In fact the two by two bare inverse gives exactly
C0(p)=e_common^2/p^2+(e_rot^2-e_common^2)/(p^2+6), with
e_common^2=2T_q/K and e_rot^2-e_common^2=T_q(q-2)^2/[K(q+1)^2]>=0.
The infrared cutoff below is necessary for this common rotation mode.

Retain momenta m_IR<=|p|<=Lambda, and multiply by smooth source windows
of magnitude at most one supported in I and a time interval tau. A
finite Fourier regulator on this annulus, or its smooth finite-rank
approximation, has point covariance trace bounded in the continuum
annulus reference by

```text
sigma0^2 <= (e_0^2+3e_rot^2)/(2pi) log(Lambda/m_IR).         (12)
```

For a finite quadrature use its actual covariance trace in (7)-(9),
not an unproved assertion that its discrete sum is below the integral.
The transverse projector has trace one, which fixes the coefficient.
This is a regulator/reference choice, not a derivation of the physical
mouth's infrared boundary condition.

At each fixed annulus, smooth finite-rank Gaussian approximations can
also be passed to their band-limited Gaussian limit on the compact
source region. Their samples converge in spatial and temporal smooth
norms, and (5) supplies dominated convergence of the determinant
integrals. This does not send Lambda to infinity or m_IR to zero at
fixed q.

With g_c^2=q^-2, K=pi q^4, tau=O(log^2 q), m_IR comparable to tau^-1
and Lambda any fixed power of q, (12) is O(log q/q). For the midpoint
return Lbar=O(1), Gaussian moments and (7) give

```text
kappa=O(log^3 q/q), rho=O(log^(7/2) q/q^2),
eta=O(log^3 q/q) ->0.                                     (13)
```

This closes all-amplitude fermion reweighting for this regulated inner
gauge reference with unequal fermion returns. It does not replace the
photon's physical exterior with a proven one. Compact holonomies and
zero modes were fixed, nonlinear gauge vertices were omitted, and
removing the photon cutoff at fixed q is not justified by (12). The
normalization is a finite-time source functional, not a renormalized
stationary vacuum energy. Reflection positivity and a Lorentzian
coupled state have not been proved for the chosen regulator/window.

## Comparison, verification and decision

[Fry, exact v1](https://arxiv.org/pdf/hep-th/0111244v1), section 2,
equation (3), explicitly reviews the established two-dimensional
non-Abelian upper bound and the separate lattice-to-continuum issue.
Its references lead to [Brydges, Frohlich and Seiler (1979)](https://doi.org/10.1016/0003-4916(79)90098-8)
and Seiler's constructive gauge theory work. The IHES/P/78/237 original
preprint was retrieved: Lemma 2.2 and the full Theorem 2.3 proof use
reflection positivity and iterated Schwarz inequalities on an even-sided
periodic rectangular lattice. Its remark distinguishes more general
boundary conditions. Those hypotheses do not directly describe these
unequal circles. We do not claim their
general diamagnetic bound as new or silently transplant their planar
domain to unequal circles. Equations (1)-(5) give the argument for the
present free directed-return model and its specified subtraction.

The [checker](../checks/determinant_tail_control.py) compares static
thermal traces to the ordered holonomy formula, tests noncommuting
time evolution and the trace inequality, independently integrates
a Gaussian determinant weight, and checks the global det_3 factorization
with a distinct spacetime regulator. Its [receipt](../evidence/determinant-tail-control.json)
records regulator refinement and exact inputs. Bounded numerical tests
do not supply a physical gauge covariance or quantum wormhole.
The [source review](../evidence/determinant-tail-source-review.json)
discloses the reading depth and unsuccessful retrievals.

Continue physical long-return damping as primary and matched finite-time
transfer as fallback. The specified inner Euclidean gauge marginal now
has global tail control. The next physical step is the actual constrained
gauge covariance and exterior action, followed by stress and matching.
No new physical damping, survival, transport or shortcut is accepted.

Remaining to-do list: derive the physical coupled gauge action, global
sectors and exterior response; control the quantum state, renormalized
stress and omitted modes; complete smooth matching, the physical
endpoint and novelty; independently verify, GCP-certify suitable
mathematics and prepare the reproducible manuscript.
