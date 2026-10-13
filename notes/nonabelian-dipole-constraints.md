# Angular Gauss law with the dynamical electromagnetic dipole

S072, 2026-10-13. The round-sphere Einstein-Maxwell projection admits a
single covariant action for all three rotation generators and the
electromagnetic dipole. Its angular Gauss law includes the dipole's own
angular momentum. The one-axis Gaussian calculation loses this term
because its cross products vanish. The derived action also gives an
explicit second-order quadrupole source in the discarded Einstein
equations. Thus this projection supplies the next constrained Hamiltonian,
but is not by itself a consistent nonlinear four-dimensional truncation
or a bound on the full angular boundary correlator.

This advances the missing step after [the radial one-axis result](radial-axial-correlator.md).
The initial gravitational bounds in [the Schwarzian calculation](schwarzian-initial-state.md)
remain separate. No product of these initial correlators, three Abelian
copies, or an additional independent rotor is used here.

## Sphere geometry and an Abelian Maxwell field

Use a unit three-vector n on the sphere, a rotation connection A, and an
adjoint one-form H on the two-dimensional base. The gravitational constant
is G_N; the four-dimensional electromagnetic coupling is g. Define

```
D n = d n - A cross n,       D X = d X - A cross X,
F_A = d A - (1/2) A cross A,
(F_A)_tx = partial_t A_x - partial_x A_t - A_t cross A_x,
(D H)_tx = partial_t H_x - partial_x H_t
           - A_t cross H_x + A_x cross H_t.                 (1)
```

The cross product of vector-valued forms includes their wedge product.
For an infinitesimal rotation epsilon, delta A=D epsilon,
delta H=epsilon cross H and delta n=epsilon cross n. Use the metric and
four-dimensional Maxwell field strength

```
ds4^2 = ds2^2 + r^2 D n dot D n,
Omega_A = (1/2) epsilon_abc n^a Dn^b wedge Dn^c + n dot F_A,
F4 = f0 + P Omega_A + d(n dot H),    f0=d a0, P=q/2.        (2)
```

Here f0 is pulled back from the base. The sphere flux is 4pi P=2pi q.
In a local sphere trivialization direct expansion gives

```
Omega_A = omega_S2 + d(n dot A),
d(n dot H) = Dn dot-wedge H + n dot D H.                  (3)
```

For example, the base-base component of the first term in Omega_A is
n dot (A_t cross A_x). It cancels the commutator in n dot F_A.
The mixed components are partial_i n dot A_mu with their exterior-form
signs. Equation (3) proves closure, patchwise and hence globally; (2)
is gauge invariant and does not replace Maxwell theory by a non-Abelian
electromagnetic theory. This is the S2 global angular form, normalized
as Omega_A/(4pi). The general construction and its closure are reviewed
in [Bonetti et al., section 2.1, equations (2.18)-(2.21)](https://doi.org/10.1007/JHEP05(2023)156).
Their special consistent supergravity reductions are not assumed here.

The horizontal Maxwell curvature in the metric (2) is

```
Gcal = P F_A + D H,       F4_horizontal = f0 + n dot Gcal. (4)
```

H is an adjoint massive vector, not a second gauge connection. In
particular, if V=P A+H then

```
P F_(V/P),tx = Gcal_tx - (H_t cross H_x)/P.               (5)
```

The last term in (5) is absent from the Abelian field (2). Treating both
old axis potentials as independent non-Abelian connections would insert
a spurious interaction.

## Action and the retained matter projection

Integrating the exact substitution (2) into the classical bosonic action,
up to the spherical gravity boundary term, gives

```
S_b = S_sph + integral sqrt(-g2) [
 - kappa0 f0_mn f0^mn/4 - kappag F_A,mn dot F_A^mn/4
 - kappae Gcal_mn dot Gcal^mn/4 - mu H_m dot H^m/2 ],

kappa0=4pi r^2/g^2, kappag=r^4/(6G_N),
kappae=4pi r^2/(3g^2), mu=8pi/(3g^2),

S_sph=(1/(4G_N)) integral sqrt(-g2) [
 r^2 R2 + 2 (partial r)^2 + 2 + 6r^2/L4^2
 - 2 Qmag^2/r^2 ],       Qmag^2=4pi G_N P^2/g^2.          (6)
```

To derive the Maxwell coefficients, use integral n_a n_b=4pi delta_ab/3
and integral (delta_ab-n_a n_b)=8pi delta_ab/3. The horizontal, mixed
and vertical pieces of F4 are orthogonal. The vertical piece supplies
the magnetic term already included in S_sph. The horizontal scalar-vector
cross term integrates to zero. For gravity the connection term in R4 is
-r^2 (delta_ab-n_a n_b) F_A^a F_A^b/4. Multiplying by the sphere volume
element r^2 and 1/(16pi G_N) gives kappag. Equivalently its normalization
follows from the stationary ADM derivation in
[angular-boundary-response.md](angular-boundary-response.md), with rotation
invariance and the two-derivative connection curvature fixing its covariant
extension. Radius gradients are retained in S_sph and in the coefficients;
no constant-radius assumption has entered (6).

For the q massless Dirac lowest-Landau channels let j=(q-1)/2 and
[t_a,t_b]=i epsilon_abc t_c. The diagonal LLL projection has

```
P_LLL n_a P_LLL = alpha t_a,       alpha=2/(q+1),
B = A + alpha H,
D_mu psi = [partial_mu + i a0_mu + i B_mu dot t] psi.     (7)
```

The matrix coefficient follows by integrating the normalized monopole
wavefunctions: <j,m|cos(theta)|j,m>=m/(j+1); rotation covariance fixes
the other two components. The lifted rotation generator is t_a.
The compensating monopole piece in (2) is why the connection couples
through this generator while the residual H couples through the projected
n_a. On the Cartan axis, A=w and H=v-Pw, so
B=(w+2v)/(q+1), exactly the preceding action.

The spherical spinor rescaling removes the explicit r and conformal
factor from its massless kinetic term; the spin-j multiplet and this
rescaling are given in [Maldacena, Milekhin and Popov, appendix A.1 and C](https://arxiv.org/html/1807.04726v3).
The extra connection-curvature term of a fibred Dirac operator contains
one sphere gamma matrix and has vanishing matrix elements within the
single-chirality LLL. Virtual higher-Landau transitions are not removed
by that statement. Equation (7) is the direct LLL projection, with those
transitions and other harmonics still omitted. Fermion determinants and
their local stress require a common gauge-covariant renormalization; (6)
has not silently added a separate current-polarization mass.

## Canonical reduction and angular momentum of the dipole

In optical coordinates ds2^2=f(-dt^2+dx^2), set F=(F_A)_tx and
G=Gcal_tx. The angular boson Lagrangian density is

```
L_b,ang = kappag |F|^2/(2f) + kappae |G|^2/(2f)
          + mu (|H_t|^2-|H_x|^2)/2,
E_H=kappae G/f,       E_A=kappag F/f+P E_H.               (8)
```

Write rho=psi^dagger t psi and j_x=psi^dagger sigma3 t psi, with the
fermion Hamiltonian H_F(B_x)=H_F(0)+j_x dot B_x. Integrating the spatial
derivatives of A_t and H_t by parts in the Legendre transform yields

```
H_ang = f |E_A-P E_H|^2/(2kappag) + f |E_H|^2/(2kappae)
        + mu |H_x|^2/2 + H_F(B_x)
        + A_t dot C + H_t dot b - mu |H_t|^2/2,
C = -D_x E_A + H_x cross E_H + rho,
b = -D_x E_H + alpha rho.                               (9)
```

The boundary term is [A_t dot E_A + H_t dot E_H]_inner^outer.
It must be retained or Legendre transformed according to the selected
fixed-source/flux ensemble. Gauge transformations with compact support
give the bulk constraint; rotations at an outer boundary generate its
charge rather than imposing its vanishing. No horizon ensemble is fixed
by the local calculation.

A_t imposes C=0, while H_t=b/mu is auxiliary. Eliminating it gives

```
H_ang,reduced = f |E_A-P E_H|^2/(2kappag)
 + f |E_H|^2/(2kappae) + mu |H_x|^2/2
 + |D_x E_H-alpha rho|^2/(2mu) + H_F(A_x+alpha H_x),
D_x E_A = rho + H_x cross E_H.                          (10)
```

This proves positivity of the displayed classical boson squares. It
does not prove positivity of a normal-ordered interacting quantum
Hamiltonian or supply its counterterms. Products of quantum densities
in (10) need a regulator preserving the vector gauge constraint.

The generator integral epsilon dot C has the canonical action
delta A_x=D_x epsilon, delta H_x=epsilon cross H_x,
delta E_A=epsilon cross E_A and delta E_H=epsilon cross E_H.
The bosonic charge is the usual canonical rotation moment map; the
fermionic bilinears supply the same su(2) algebra. The vector Dirac
theory has cancelling left/right gauge anomalies, although its chiral
current algebra and quantum contact terms must still be retained.

In radial gauge A_x=0, the actual flux obeys

```
E_A(x)=E_A(inner)+integral_inner^x
       [rho(y)+H_x(y) cross E_H(y)] dy.                  (11)
```

The second term in the integrand is the angular momentum carried by
the electromagnetic dipole. It cannot be assigned to an independent
global rotor or omitted in favor of fermion spin alone. Along a single
axis it vanishes. In a static Cartan configuration E_A=-(p_w+P p_v)
and E_H=-p_v, recovering the physical rotation flux of the previous
notes. The spherical electric Gauss law remains partial_x E0=rho0.

For a quantitative diagnostic the local positive terms obey

```
|H_x cross E_H| <= r/sqrt(2f) *
 [mu |H_x|^2/2 + f |E_H|^2/(2kappae)].                  (12)
```

Weighted Cauchy-Schwarz proves this bound, with equality for perpendicular
vectors of balanced energies. The weight grows at a horizon. Therefore
these two energy terms alone do not give a uniform unweighted angular
charge bound there; the other constraints and physical state must enter.
No conclusion about a gap of the complete Hamiltonian follows from (12).

The original matter charge-spin parity is unchanged: 2J=(q-1)Q mod 2,
since adjoint dipole quanta have integer spin. A Wilson line of A and a0
to the fixed outer frame supplies a gauge-covariant dressing of a LLL
insertion before radial gauge. A line of B is another covariant choice
but differs by physical H insertions. Gauge invariance alone does not
fix this choice or the finite four-dimensional operator normalization.

## Explicit failure of closure and the next physical calculation

The exact action substitution above need not enforce omitted field
equations. Consider local data with A=0, F_A=0, H=0, f0=0 but
(D H)_tx=u !=0. Spatially constant u with rho=0 satisfies the projected
angular Gauss constraint at the initial point/patch: E_A=P E_H and
partial_x E_A=0. The Maxwell field has radial electric component n dot u.
In a horizontal orthonormal frame its excess energy density is

```
T_hat0hat0,electric = (n dot u)^2/(2g^2 f^2),
(n dot u)^2 = |u|^2/3 + (n_a n_b-delta_ab/3) u_a u_b.   (13)
```

The second term is a nonzero l=2 harmonic. With normalized sphere average
<.>_S2, its exact magnitude is

```
<[(n dot u)^2-|u|^2/3]^2>_S2 = 4 |u|^4/45,
RMS(quadrupole)/monopole = 2/sqrt(5).                    (14)
```

The four-index integral <n_a n_b n_c n_d>=(delta_ab delta_cd+
delta_ac delta_bd+delta_ad delta_bc)/15 proves (14). In the round ansatz
the horizontal Einstein tensor depends on the connection through its
curvature; its horizontal symmetric component is angle independent at
F_A=0. Spherical radius/lapse changes cannot cancel (13). With the
independent l=2 metric/Maxwell response discarded, the omitted equation
therefore has a nonzero source already at second perturbative order.
This is an off-shell closure obstruction and valid initial projected
data; it is not a global four-dimensional solution or an instability.
The same effect occurs for one-axis dipole data, so it does not contradict
the explicitly quadratic scope of the earlier result.

Consequently the next controlled comparison must retain the sourced
quadrupole (and subsequent harmonics as needed), derive its coupled
retarded Einstein-Maxwell response, and bound its renormalized contribution
in the actual driven state. Equation (14) is a ratio to the dipole's own
electric stress, not to the total magnetic background or capture margin.
It proves neither an order-one correction to the candidate nor a small
one. Rotational invariance can annul a mean quadrupole while leaving its
correlations nonzero; a quantum stress bound requires smearing and the
common subtraction. The interacting boundary state must solve (10)-(11),
not the fermion-only radial spin constraint.

## Verification and claim boundary

[The checker](../checks/nonabelian_dipole_constraints.py) independently
integrates the coordinate four-dimensional Maxwell tensor using the full
inverse fibred metric, and compares it with (6). It tests noncommuting
connection/dipole data, finite local rotations and their derivative jets,
the Legendre transform with its boundary divergence, the exact canonical
moment-map algebra, the Cartan normalization, the LLL matrix projection
and the quadrupole angular integral. These supplement the derivation;
they are not a four-dimensional evolution or a renormalized quantum test.

The evidence record reports the numerical tolerances, source identities
and any development failures. No exhaustive novelty claim is made for
the geometric identities or constrained Hamiltonian methods. The new
result for this program is the correctly coupled projected Hamiltonian
and explicit omitted-mode source that the next physical error estimate
must use. No capture, area, preparation or creation gate is closed here.

Remaining to-do list: solve or bound the non-Abelian constrained state
and sourced higher harmonics; control massive transitions and physical
operator/horizon matching; obtain driven renormalized stress and geometry;
establish preparation, quantum stability and finite-energy transport.
