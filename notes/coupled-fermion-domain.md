# Coupled fermion domain, boundary work and ultraviolet subtraction

S074, 2026-10-13. The specified free four-dimensional Dirac boundary
coupling has a maximal charge-conserving domain, an explicit work term
when driven, and coupling-independent local heat coefficients for a
constant angle on identical static collars. A time-dependent angle
instead generates a boundary source-kinetic divergence. On the product
collar calculated here its unsmeared logarithmic divergence vanishes,
but the smeared heat coefficient relevant to a Weyl variation does not.
The physical evolving determinant, finite source term and interacting
state remain to be matched. No opening or stability bound follows yet.

## Domain before reducing to the lowest Landau level

Use the massless four-component Dirac field, identical boundary clocks,
round spheres and matched magnetic bundles of the selected source model.
[Bintanja--Freivogel--Rolph, sections 2.1.1--2.1.3 and appendix A](https://arxiv.org/html/2308.00871v2)
give its radial projectors and coupled condition. After the canonical
radial spinor rescaling and a constant basis rotation, write

```
H=-i alpha partial_x + beta a(x) D_q,
alpha=sigma_2, beta=sigma_1, a=exp(sigma)/R,
P_beta,+ Psi_R+h P_beta,+ Psi_L=0,
P_beta,- Psi_L+h P_beta,- Psi_R=0.                    (1)
```

D_q is the self-adjoint magnetic Dirac operator on the unit sphere.
It acts on the other two spinor components. Equation (1) retains that
operator, both nonzero eigenvalue signs and the zero modes; it is not
an LLL truncation. At conformal infinity use finite canonical boundary
data in a common conformal frame. Arbitrary finite cuts need their own
matching and are not new physical boundaries.

Pull both boundary collars back to an inward coordinate r, with flavor
order (L,R) and Pauli matrices tau_i on that flavor. Put

```
A=tau_3 tensor alpha, B=-A,
G=i tau_2 tensor beta, G^dagger=-G, G^2=-1,
theta=2 atan h, V=exp(theta G/2),
chi_0=tau_3 tensor beta,
chi_theta=cos(theta) chi_0-sin(theta) tau_1 tensor 1,
Pi=(1+chi_theta)/2=V Pi_0 V^dagger.                  (2)
```

An identity on angular spinors is understood in these matrices. The
allowed trace is Psi=Pi Psi. For beta eigenvalue + its flavor vector
is proportional to (1,-h); for beta eigenvalue - it is (-h,1).
Thus (2) is exactly (1), with the previous free-angle convention in
[boundary-vertex matching](boundary-vertex-matching.md).

The Green form is -i integral(Phi_R^dagger alpha Psi_R-
Phi_L^dagger alpha Psi_L)=-i integral Phi^dagger B Psi. Since
chi_theta^2=1, chi_theta^dagger=chi_theta and {chi_theta,B}=0,
Pi B Pi=0. Its rank is four out of the eight folded spinor components,
so it is maximal isotropic, not just zero flux for one chosen state.
The corresponding smooth compact, regular Dirac boundary problem has
the self-adjoint Hamiltonian domain; horizon ends or infrared zero modes
need their separate usual domain prescription. Arbitrary angular traces
satisfy the identity. Total charge is conserved even for theta(t).

The boundary identification includes a spin-frame and bundle map.
For a unitary map U from R to L intertwining the radial Clifford
matrices, (1) becomes P_+ Psi_R+h U^dagger P_+ Psi_L=0 and
P_- Psi_L+h U P_- Psi_R=0. Gauge changes send U to
exp(i lambda_L)U exp(-i lambda_R). In global angular coordinates the
pulled-back monopole bundles agree; opposite outward magnetic charges
do not license pairing inequivalent bundles. A prescribed boundary
coupling fixes relative global frames. It does not add an independent
link rotor to already dressed operators or solve the quantum Gauss law.
The remaining results take this map fixed and equal to the identity.

## Boundary work for a varying domain

For a classical smooth solution on a fixed background let
E=<Psi,H Psi>, and suppose there is no other boundary flux. Integration
by parts, including the fact that partial_t Psi need not obey the
instantaneous homogeneous condition, gives

```
dE/dt=<Psi,(partial_t H)Psi>+P_boundary,
P_boundary=-i integral Psi^dagger B (partial_t Pi) Psi
          =theta_dot/2 integral v^dagger C v,
v=V^dagger Psi|boundary, C=i A G, C^dagger=C, C^2=1.  (3)
```

Terms involving an arbitrary derivative within the allowed trace space
vanish by Pi B Pi=0. Hence |P_boundary| is at most
|theta_dot|/2 times the squared boundary trace norm for a smooth
one-particle wavefunction. This is not a bound on an unrenormalized
quantum trace density or on the controller's energy. A moving frame or
link has additional derivatives in Pi and in H. The general first line
retains them. Constant h on a static background has zero total boundary
work; this statement cannot be imposed during switching.

## Static local subtraction on identical collars

For identical spherical static collar data, in the radial half-density
representation the folded Hamiltonian is

```
H_fold=-i A partial_r + 1 tensor beta tensor a(r) D_q. (4)
```

The generator G commutes with both coefficients in (4). A constant V
therefore intertwines the differential expression and its boundary
domain in an entire collar. Squared Euclidean operators and their
induced domains are likewise locally unitarily equivalent. The same
argument permits an equal scalar electric potential on the two collars,
but not unequal jets, unequal clocks or arbitrary metric perturbations.

The local heat coefficients are consequently identical to those at h=0,
for fixed constant h and matched collar jets. This includes every
angular channel in the common four-dimensional regulator. The argument
uses locality of the smooth mixed-boundary heat expansion, reviewed in
[Vassilevich, sections 2.1 and 5.3](https://arxiv.org/html/hep-th/0306138v3).
Folding makes the paired condition local in the common boundary point;
it does not turn it into an unrelated reflecting physical wall.

There is no asserted global unitary equivalence. A collar rotation need
not preserve the domain where the collars join the interior. For the
massless interval of width pi, direct propagation in (1) yields

```
(1+h^2) cos(pi omega)+2h=0,
omega_k=k+1/2+(-1)^k theta/pi, k in Z, |theta|<pi/2.
K_theta(s)=2 sum_(n in Z) exp[-s(2n+1/2+theta/pi)^2]. (5)
```

These frequencies check the source spectrum. Poisson summation gives

```
K_theta(s)=sqrt(pi/s) sum_(m in Z)
 exp[-pi^2 m^2/(4s)] exp[i pi m(1/2+theta/pi)],
|K_theta-K_0| <=4 sqrt(pi/s) exp[-pi^2/(4s)]
                         /[1-exp(-3pi^2/(4s))].     (6)
```

The spectrum and finite determinant change, while their local ultraviolet
expansions agree. An infrared prescription remains necessary for a full
determinant, particularly at a zero crossing outside the selected branch.

## Driven product collar and its induced mixed condition

For the following coefficient calculation take the explicit Euclidean
product collar dr^2+d tau^2+b_c^2 dOmega^2, a=1/b_c constant, fixed
magnetic connection, no electric potential and a smooth angle theta(tau)
uniform on the sphere. The two endpoints are folded once, so the fiber
rank is eight and the boundary integral is over one common S^2 times
Euclidean time. Initial/final temporal boundaries are absent (periodic
time or a compactly supported local calculation).

Set D_E=partial_tau+H_fold. Its adjoint has the same instantaneous
allowed trace. For L=D_E^dagger D_E, membership of D_E Psi in that
adjoint domain gives

```
(1-Pi)Psi=0,
Pi(partial_r+S)Psi=0,
S=i A (partial_tau Pi)Pi=(theta'/2) C Pi.             (7)
```

Indeed (1-Pi)partial_tau Psi=(partial_tau Pi)Psi,
whereas the angular term preserves Pi. Extend Pi normally constant.
S is Hermitian and supported on Pi. These are mixed Dirichlet/Robin
conditions with scalar Laplace principal symbol; they are strongly
elliptic. Dropping S when differentiating the domain is incorrect.
The bulk squared operator is -partial_tau^2-partial_r^2+D_q^2/b_c^2;
its curvature endomorphism acts only on the angular factor.

For K_L(s;f)~sum_j s^((j-4)/2) a_j(f), evaluate the mixed coefficients
in Vassilevich (5.29)--(5.33), with the **inward** normal partial_r.
Relative to any constant angle, elementary fiber traces give

```
tr[(chi')^2]=8(theta')^2, tr S=tr S^3=0,
tr S^2=(theta')^2,
Delta a_0=Delta a_1=Delta a_2=0,
Delta a_3(f)=1/[4(4pi)^(3/2)] integral f (theta')^2,
Delta a_4(f)=4/[15(4pi)^2] integral f_;r (theta')^2.  (8)
```

The two derivative contributions to a_3 are
tr[192S^2-12(chi')^2]/384; both are essential. In a_4 all unsmeared
terms vanish on this product: K_ab=0, E_;r=R_;r=Omega_tau,r=0,
and tr(SE), tr(S^3), tr(S_:aa), tr((chi')^2 S) vanish. The surviving
normal-smearing term is tr[240S^2-18(chi')^2]/360. Thus a_4(1)
does not change, but the coefficient relevant to a general Weyl
smearing does change. Absence of an unsmeared logarithm is not absence
of a boundary Weyl contribution. Unequal/background-dependent collars
require their additional invariants and are not covered by (8).

For one complex Dirac field, the modulus effective action is
Gamma_E=-log|det D_E|=-1/2 log det L. With proper-time lower cutoff
Lambda_UV^-2, (8) gives the local divergence

```
Delta Gamma_E,div=Lambda_UV/[4(4pi)^(3/2)]
                   integral (theta')^2
                =Lambda_UV/(4pi)^(3/2)
                   integral (h')^2/(1+h^2)^2.        (9)
```

The integral includes the physical product boundary volume, not an
extra factor of two. Zeta subtraction removes this power divergence
as a convention; a finite coefficient multiplying the same source
term still needs a declared renormalization condition or physical
matching. Its metric/source variations enter the controller's stress
and work. This Euclidean positive coefficient is not a claim of
positive Lorentzian particle production. The determinant phase, full
Weyl transfer, nonlocal response and interacting state are not supplied.

## Verification and consequence for the task list

The checker verifies the complete matrix domain, Green form, work term,
collar intertwiner, induced mixed condition and heat-coefficient traces.
It also checks the interval quantization independently by propagation,
compares direct heat sums with Poisson sums, and tests differentiated
boundary traces with arbitrary angular mixing without using the
heat-coefficient formulas. The latter tests the induced mixed condition
by applying the first-order Euclidean operator to its boundary jets. Deliberate
omission of the driven Robin term fails a coefficient control. These
checks supplement the ordinary derivation and imported mixed-coefficient
theorem; they do not solve a driven four-dimensional quantum state.

Use the coupled domain rather than a guessed local fermion wall. Static
relative ultraviolet subtraction closes on matched identical collars;
the drive requires the source counterterm and its finite matching.
This advances S074.2/S074.4 and the domain/work inputs to S075 without
closing their full gates. The general local boundary action and its
smooth internal cuts remain as derived in
[spherical-anomaly-boundary.md](spherical-anomaly-boundary.md).

Remaining to-do list: complete physical curved/driven boundary matching,
finite source and bulk couplings, other S074 interactions and observable
remainder; S075 state/drive; S076 stress/geometry; S077 stability;
S078 transport; S079 preparation/creation; S080 final audit.
