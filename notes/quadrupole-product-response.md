# Sourced quadrupole response with the light Einstein-Maxwell branch retained

S072, 2026-10-13. The omitted quadrupole identified in
[the angular constraint calculation](nonabelian-dipole-constraints.md)
has a coupled gravitational/electromagnetic response, not the response of
a massive photon on a fixed metric. This note derives its gauge-invariant
source equations and retarded propagator on the exact RN-AdS product
throat. A direct four-dimensional tensor calculation checks the quadrupole
equations, including electric and magnetic effective currents. The lower
mass branch is very light in the candidate family. It must remain in the
time-dependent comparison rather than being eliminated by a heavy-mass
approximation.

This is a classical linear response to prescribed sources on a fixed
product background. The actual renormalized dipole/fermion source
correlations, its state, the evolving background and physical boundaries
remain to be supplied. No capture or creation gate closes here.

## Product background, source variables and duality

Let a be the AdS2 radius and b the sphere radius. Write

```
X=b^2/L4^2, C=b^2/a^2=1+6X, Y=1+3X=(C+1)/2,
Qgeom^2=b^2 Y, kappa^2=8pi G_N,
ds2^2=-B dT^2+d xi^2/B, B=1+xi^2/a^2,
kappa^2 E0^2=2Y/b^2.                                    (1)
```

Rescale the Maxwell field by its four-dimensional coupling g so that
its action is -Fbar^2/4, and use electromagnetic duality to make the
background electric, Fbar_T,xi=E0>0. This is a device for deriving the
response of the original magnetic throat. Under duality, an original
electric source becomes a magnetic source. It must not be discarded.
The background satisfies the full Einstein-Maxwell equations with
Lambda=-3/L4^2, for all X>=0.

For a scalar harmonic Y_lm with L=l(l+1), l>=2, choose the usual gauge
with no mixed metric components or traceless scalar sphere deformation.
Define the invariant sphere amplitude chi by
delta g_ij=2b^2 chi gamma_ij Y_lm and the invariant base metric h_ab.
Let h2=h^a_a. Define the electric master Acal by

```
delta Fbar_ab = epsilon_ab Ecal Y_lm,
Ecal = -L Acal/b^2 + (E0/2)(h2-4chi),
J_e^a = (L/b^2) epsilon^ab Jtilde_b,
(d delta Fbar)_ab i = epsilon_ab Bcal D_i Y_lm.           (2)
```

Here epsilon_T,xi=1. Bcal is the magnetic-current source in this scalar
sector, not the metric function B. With no magnetic current it vanishes.
The sign convention for J_e is nabla_nu Fbar^(mu nu)=J_e^mu.

Move all prescribed extra stress, including any nonlinear lower-order
residual, to the right of the *linear* Einstein-Maxwell equations. Use
S_ab=kappa^2 T_eff,ab for its base harmonic coefficient and
S_T=2b^2 kappa^2 tau_T/L for its scalar traceless sphere coefficient.
The traceless sphere equation fixes h2=-S_T. These sources must obey
the coupled conservation identities; an isolated choice of T_00 is not
a complete quadrupole source.

## Local source equations, including zero frequency

Eliminating the metric constraints gives

```
[Box2-L/b^2-2/a^2] chi - (L kappa^2 E0/b^2) Acal = r_F,
[Box2-L/b^2] Acal - 2E0 chi = r_A,
r_F = (1/2) S^a_a - L S_T/(4b^2),
r_A = (E0/2) S_T - D_a Jtilde^a + Bcal.                 (3)
```

The first source depends on the base stress trace and sphere anisotropy.
It follows either directly from the four-dimensional equations or by
using source conservation to cancel the inverse frequencies in the
constraint elimination. For example, with convention exp(-i omega T),
the time component of conservation is

```
-i omega S^T_T + partial_xi S^xi_T + (L/b^2) S_time
  = (L kappa^2 E0/b^2) Jtilde_T.                        (4)
```

The combination of the last three terms occurring in the frequency
constraint is therefore a local trace term. Equation (3) remains valid
at omega=0; no inverse-frequency prescription is required. The direct
tensor check does not use (4) to assume away a residual.

For the Maxwell equation, the covariant divergence in these coordinates
is

```
D_a Jtilde^a = i omega Jtilde_T/B
               + partial_xi(B Jtilde_xi).              (5)
```

The Bianchi equation supplies Bcal, giving the last line of (3).
Keeping (5) and the coefficient E0/2 matters. The covariant equation
(105) and trace constraint (110), numbered (5.20) and (5.25) in the PDF, of
[Kodama and Ishibashi](https://arxiv.org/html/hep-th/0308128v3)
give these factors. A literal use of their printed constant-radius
source expression (160), PDF (5.75), would instead give E0 S_T and omit the
B' Jtilde_xi term. We use the covariant equation and the independent
tensor check, not that literal source expression. The expression was
also visually checked on PDF page 27. The source-free spectrum is
unaffected by this discrepancy.

For reference, the remaining metric can be reconstructed without dividing
by omega. Put h^T_T=U+2chi, h^xi_xi=V+2chi, h^xi_T=Z. Then

```
U = (4b^2/L)[(omega^2/B-L/b^2-1/a^2)chi
                +(xi/a^2) partial_xi chi]
       -2kappa^2 E0 Acal -(2b^2/L)S^xi_xi,
V = -U-4chi-S_T,
Z = -i omega (4b^2/L)[-xi chi/a^2+B partial_xi chi]
       +(2b^2/L)S^xi_T.                                (6)
```

The remaining source components must satisfy conservation. For zero
sources the checker substitutes arbitrary radial master data obeying
(3) into (6) and verifies every linear Einstein and Maxwell component
exactly for l=2, arbitrary frequency and symbolic radius ratio C. For
general sources it verifies both identities (3) directly with arbitrary
metric and Maxwell perturbations, including nonzero d delta Fbar.

## Mass matrix and physical scale in the candidate family

Set s=sqrt(2b^2/(L kappa^2)) and z=(s chi,Acal)^T. In dimensionless
global coordinates tau=T/a, rho=xi/a, equations (3) become

```
[a^2 Box2-M_L] z = a^2 (s r_F,r_A)^T,
M_L = [[2+L/C, 2sqrt(LY)/C],
       [2sqrt(LY)/C, L/C]],
R_L=sqrt(C^2+4LY),
mu_plus  = (C+L+R_L)/C,
mu_minus = L(L-2)/[C(C+L+R_L)].                         (7)
```

Here mu=m^2 a^2. The rationalized expression avoids catastrophic
cancellation for the light eigenvalue. The trace is 2+2L/C and the
determinant is L(L-2)/C^2. Thus both branches are positive for l>=2.
At X=0 they reduce to l(l-1) and (l+1)(l+2). These masses also follow
independently by taking the extremal horizon limit of both polar and
axial RN-AdS potentials in
[Berti and Kokkotas, equations (13)-(17)](https://arxiv.org/html/gr-qc/0301052v2).
In that limit their W=1 and U=C+L. This does not assert equality of
the exterior spectra with different AdS boundary conditions.

The spectral projectors are

```
P_plus/minus = (1/2)[I +/- (1/R_L)
                        [[C,2sqrt(LY)],[2sqrt(LY),-C]]]. (8)
```

For fixed l and large C, mu_minus~L(L-2)/(2C^2). For the retained
candidate q=k^4, g=k^-4, L4/l_p=k^5, one has
C=sqrt(1+12pi k^6). At k=2000 the l=2 branch has

```
mu_minus = 4.9735919710e-21,
mu_plus  = 2.0000000002443013,
Delta_minus-1 = 4.9735919710e-21,
Delta=(1+sqrt(1+4mu))/2.                                (9)
```

It is not a heavy field in AdS2 units. Its metric-source projector
weight is nevertheless small:

```
(P_minus)_11 = (R_L-C)/(2R_L)
             = 2LY/[R_L(R_L+C)],                       (10)
```

equal to 6.1075314e-11 for that quadrupole. If r_A=0, this factor
appears in the metric-to-metric retarded response. It is not a bound
on the full physical source: dipole residuals and dualized fermion
currents can also drive r_A, and their normalization must be retained.

The small-mass sector is not restricted to one harmonic. Solving
mu_minus<=1 gives the exact threshold

```
L <= C+1+sqrt((C+1)^2+C^2).                             (11)
```

At k=2000 this permits 2<=l<=344361, or 118585187040 angular states
on this one polar branch. Their count divided by the q complex LLL
channels is 0.00741157419. This is a count of modes of the product
operator, not an occupation, stress, loop-size or capture-error estimate.
It rules out treating every nonspherical correction as a uniformly
heavy finite set.

## Retarded response with its boundary choice explicit

Use rho=tan(sigma), -pi/2<sigma<pi/2. Define
j=-a^2 sec^2(sigma)(s r_F,r_A)^T. The diagonalized wave equation is

```
partial_tau^2 z + H z = j,
H=-partial_sigma^2+M_L sec^2(sigma).                    (12)
```

Choose the Friedrichs realization, with the faster falloff at both ends.
For each branch, normalized eigenfunctions are

```
u_n,Delta(sigma)=N_n cos^Delta(sigma)
                         C_n^Delta(sin sigma),
omega_n=n+Delta,
N_n^-2=pi 2^(1-2Delta) Gamma(n+2Delta)
                 /[n!(n+Delta) Gamma(Delta)^2].         (13)
```

The Gegenbauer equation proves the eigenvalue and its orthogonality
integral fixes the normalization. Its polynomial completeness in the
weighted interval gives the full spectral expansion. The retarded kernel
for zero initial data is

```
G_R(t;sigma,sigma') = theta(t) sum_(branch,n)
 P_branch u_n,Delta(sigma)u_n,Delta(sigma')
                    sin[(n+Delta)t]/(n+Delta),
z(t)=integral_0^t G_R(t-v) j(v) dv.                     (14)
```

The series defines the wave propagator spectrally/distributionally,
rather than by an unjustified pointwise absolute sum. On L2 it obeys

```
||G_R(t)|| <= min(t,1/Delta_min), Delta_min>=1,
||z(t)|| <= integral_0^t min(t-v,1/Delta_min)||j(v)|| dv,
sqrt(2 E_H(t)) <= integral_0^t ||j(v)|| dv,
E_H=(||z_t||^2+<z,H z>)/2.                              (15)
```

For the last estimate, differentiate E_H and apply Cauchy-Schwarz to
<z_t,j>, then integrate from zero initial data. This is a mathematical
master energy in the normalization of (7); converting it to physical
four-dimensional energy or local area requires the reconstructed fields,
action normalization, source and boundary matching. In particular,
derivative observables need corresponding source regularity.

There is no 1/mu_minus penalty in this finite-time estimate. The global
spatial gradient and the chosen boundary domain supply its gap. For
0<=mu<3/4 both endpoint falloffs are square integrable; square
integrability alone does not select the Friedrichs domain. Its gap cannot
be assigned to a physical horizon or another boundary realization.

On the product Poincare exterior, the corresponding spatial operator is
H_P=-partial_z^2+mu/z^2 on z>0, with a Friedrichs AdS endpoint and the
horizon at z=infinity. It is nonnegative but has no positive gap. A
normalized smooth packet of width R centered at R^2 has kinetic energy
O(R^-2) and potential energy O(mu R^-4), proving the lack of a gap.
The universal bounds ||G_R(t)||<=t and the energy estimate in (15)
still hold. Thus a finite-time retarded bound is available without
pretending that the exterior has the global gap. Neither reference
settles the moving-mouth boundary problem.

## Verification, the source error found, and continuation

[The checker](../checks/quadrupole_product_response.py) derives the
linearized four-dimensional connection, Ricci tensor and Maxwell stress
in sphere coordinate c=cos(theta). It verifies the background, the two
sourced master identities and all fourteen source-free Einstein/Maxwell
components after reconstruction. It also checks the independent exterior
potential limits, normalized global modes, and retarded pulse evolution
against spectral convolution and integrated work.

The initial check failed the Maxwell source identity with residual
E0(T+R)/2, where T+R=h2. The incorrect source coefficient E0 S_T was
replaced by E0 S_T/2, as required by the covariant equation. That source
snapshot and failed terminal log are preserved. The refined electric-source
calculation passed; a further check with independent electric field
strength, rather than an exact potential, retains the magnetic source
required under duality. That final check also passes. Eight normalized
global-mode cases, three independent pulse/work comparisons and three
80-digit candidate substitutions pass. These numerical checks supplement
the ordinary spectral proof; no interval certificate or quantum-state
calculation is claimed.

The response derived here is the next operator needed after the
non-Abelian Gauss constraint, not yet a bound on the actual quantum
quadrupole. A dipole electric energy at one point supplies part of the
stress trace, but the full quadratic metric/Maxwell residual and matter
source must be projected together. Their connected correlations enter
the response variance through two copies of (14). An invariant mean can
vanish while those correlations remain nonzero. Pointwise bare source
moments cannot be assumed finite; use smearing, a common subtraction and
the actual dressed initial state before estimating them.

The next physical calculation is therefore the joint ordered source
kernel for r_F and r_A, including the dual magnetic-current term and
mixed correlations, in the state obeying the dipole-corrected angular
Gauss law. Propagate that kernel with the retained light branch and
reconstruct stress/geometry. Complete the physical horizon and radius-
dependent matching before comparing with the returned-switch margins.

Remaining to-do list: renormalized constrained source correlations and
physical boundaries; higher harmonics, massive transitions and evolving
geometry; quantum capture, pair preparation, sufficient stability and
finite-energy transport.
