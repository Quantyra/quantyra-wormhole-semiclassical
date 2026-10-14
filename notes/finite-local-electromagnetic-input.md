# Finite electromagnetic input and the axial spectrum

S074, 2026-10-13. Matching the magnetic RN-AdS reference geometry and
its two-derivative couplings does not determine the finite local action
needed for the opening argument. A concrete allowed operator leaves
that entire purely magnetic reference solution unchanged while changing
the axial spectrum. This note derives its constrained action, boundary
displacement, finite Gaussian measure factor and a uniform comparison
for every axial harmonic. It supplies an explicit conditional range for
one missing coefficient. It does not fix that coefficient or complete
the other local, matter and boundary terms.

## The finite coupling that the background does not determine

Use the Maxwell-normalized field Fbar=F/g, signature (-+++), and
kappa^2=8pi G. Consider the following direction in the local effective
action, holding the other renormalized coefficients fixed:

```
S(cJ)=S_Einstein-Maxwell-Dirac+integral sqrt(-g) Delta L,
Delta L=cJ kappa^4 J^2/16,
J=Fbar_(mu nu) dualFbar^(mu nu),
dualFbar^(mu nu)=epsilon^(mu nu rho sigma) Fbar_(rho sigma)/2.
                                                               (1)
```

Here cJ is dimensionless; kappa^4 means (kappa^2)^2. Its finite value
and the matching scale have not been supplied by the retained calculation.
The term is gauge invariant, diffeomorphism invariant and parity even.
It is part of the usual four-field electromagnetic operator space. In
four dimensions,

```
J^2=4 Fbar_(mu nu) Fbar^(nu rho) Fbar_(rho sigma) Fbar^(sigma mu)
       -2(Fbar_(mu nu) Fbar^(mu nu))^2.                 (2)
```

In the notation of
[Cheung, Liu and Remmen, equation 4](https://arxiv.org/html/1801.08546v3#S1),
this direction has c7=-cJ kappa^4/8 and c8=cJ kappa^4/4.
Their omission of Maxwell-divergence operators assumes no charged
matter. Our massless charged Dirac field prevents importing that
source-free operator list as the complete local basis. Their entropy
positivity theorem is not being applied to this driven AdS model.

For a purely magnetic field J=0, including the full static magnetic
RN-AdS exterior and product throat. Both Delta L and its first
functional variation vanish for every cJ. Thus the background solution,
its charge, and its mass/radius relation do not measure cJ. Around zero
field the operator also has no quadratic Maxwell term, so the same
zero-field gauge normalization can be retained.

The coefficient can affect an observable. Four null, transverse photon
configurations with incoming covectors k1=(1,0,0,1), k2=(1,0,0,-1),
k3=-k1, k4=-k2 and polarizations x,y,x,y give coefficient 256 for
t1 t2 t3 t4 in J(sum_i ti Fi)^2, in this unit-frequency convention.
The checker verifies their free Maxwell equations and (2). A field
redefinition can move this dependence into other operators and exchange
diagrams; it cannot remove the physical dependence while leaving those
terms and observables fixed. A finite subtraction convention alone does
not establish the physical value cJ=0.

## Constitutive relation, stress and boundary data

For the local electromagnetic action define its displacement tensor by
H^(mu nu)=-2 partial L/partial Fbar_(mu nu). Direct variation gives

```
H^(mu nu)=Fbar^(mu nu)-(cJ kappa^4/2)J dualFbar^(mu nu),
T_(mu nu)^Delta=-(cJ kappa^4/16)J^2 g_(mu nu),
delta S_EM|boundary=-integral_boundary sqrt(|h|)
                           n_mu H^(mu nu) delta Abar_nu. (3)
```

With the convention that the source-free equation is nabla_mu H^(mu nu)=0,
the Bianchi identity and
Fbar_(nu rho) dualFbar^(mu rho)=delta_nu^mu J/4 imply
nabla_mu T^Delta(mu)_nu=Fbar_(nu rho) nabla_mu H^Delta(mu rho).
The checker verifies the tensor identity for arbitrary electric and
magnetic components. This is the Ward identity for this local term;
it does not replace the full matter/current identity.

On a magnetic background of orthonormal magnitude B0, the quadratic
correction is cJ kappa^4 B0^2 (delta E_parallel)^2. The electric Hessian
has transverse entries 1 and parallel entry

```
1+theta, theta=2cJ kappa^4 B0^2.                        (4)
```

Fixed electric displacement at a boundary is therefore a condition on
H, not on the uncorrected Fbar. The associated Legendre term must use
that displacement. For the reference solution the electric displacement
is zero, but its response is not independent of cJ.

## Constrained axisymmetric action and its measure

Use the exact spacelike reduction from
[the axial action](axisymmetric-axial-action.md), now entirely in the
canonical Maxwell normalization:

```
ds4^2=gamma_MN dx^M dx^N+R^2(dphi+B)^2,
Abar4=A+sigma(dphi+B), Fh=dA+sigma dB,
K=2pi R [[R^2/(2kappa^2)+sigma^2,sigma],[sigma,1]].      (5)
```

Let v_I^M=(*3 F_I)^M, F_I=(dB,dA), and w=(sigma,1)^T.
The four-dimensional pseudoscalar is, with a common orientation,
J=(4/R) d sigma dot (v_A+sigma v_B). Therefore the curvature action is

```
L3,curv/sqrt(-gamma)
 =1/2 K_IJ gamma_MN v_I^M v_J^N
  +(2pi cJ kappa^4/R)[w_I v_I^M partial_M sigma]^2.      (6)
```

This is a rank-one update on the six curvature components. In matrix
notation put A0=K tensor gamma, u=w tensor d sigma and
r=4pi cJ kappa^4/R. The inverse and determinant are

```
eta=2cJ kappa^4 (grad sigma)^2/R^2,
A^-1=A0^-1-r(A0^-1 u)(A0^-1 u)^T/(1+eta),
det A/det A0=1+eta,
K^-1 w=(0,1/(2pi R))^T.                               (7)
```

The formulas require 1+eta nonzero. For the continuous positive
constitutive branch around the reference use 1+eta>0. A general driven
background may have timelike d sigma; this condition must be checked
there and cannot be inferred from the static product.

Add the Bianchi multipliers Phi_I dF_I, retaining the original boundary
Legendre term. Eliminating the six curvature variables in their
first-order action gives the original dual action plus

```
Delta Sdual=integral sqrt(-gamma)
 [cJ kappa^4/(2pi R^3(1+eta))]
                    (grad sigma dot grad PhiA)^2.      (8)
```

It is important to eliminate the first-order action, rather than insert
the duality equation into the original curvature term alone. The checker
verifies (7)-(8) exactly. Equation (8), with the prescribed even
background retained, generates the associated cubic and quartic
contacts. In particular it contains derivatives of the even Maxwell
background and is not captured by changing one constant mass by hand.

For example, in the earlier angular gauge with s=1-c^2,
rho=chi P2(c), sigma=-P[c+alpha2 V(c)], V=s P2',

```
eta=2cJ kappa^4 P^2 [
 (1+alpha2 V')^2/(b^4 exp(4rho))
 +V^2 h^ab partial_a alpha2 partial_b alpha2
                                /(b^2 exp(2rho) s)].   (9)
```

At a finite positive Euclidean Gaussian regulator, integrating the
curvature variables also changes the local measure by (1+eta)^(-1/2).
The corresponding formal local action term is (1/2)Tr log(1+eta).
Its continuum value and combination with other Jacobians, constraints,
ghosts and counterterms require the common regulator. A six-dimensional
Gaussian determinant is not a computation of that four-dimensional
quantum measure. The numerical checker tests the finite inverse and
determinant with both Lorentzian and positive Euclidean base matrices,
including off-diagonal components.

Equations (7)-(9) are exact algebra in the action with just this local
addition. An effective-field-theory expansion uses their first-order
cJ terms at this derivative order. The rational denominator does not
supply the unknown independent higher-order operators.

## All axial harmonics on the product

On the magnetic product B0^2=(C+1)/(kappa^2 b^2), so

```
C=b^2/a^2>=1, gamma=G/a^2,
theta=16pi cJ gamma(1+1/C).                            (10)
```

Only the horizontal Maxwell curvature Fh_tx changes at quadratic order.
In the dual action the (partial_c PhiA)^2 coefficient is divided by
1+theta, while the two-dimensional kinetic normalization remains the
one already derived. Integrating the canonical harmonic profiles gives,
for L=l(l+1), l>=2,

```
Ma_L(theta)=(1/C) [[L/(1+theta)+2(C+1), sqrt(2(C+1)(L-2))],
                   [sqrt(2(C+1)(L-2)), L-2]],
tr Ma_L(theta)=2+L(2+theta)/[(1+theta)C],
det Ma_L(theta)=L(L-2)/[(1+theta)C^2].                  (11)
```

The same expression holds for every m by spherical symmetry. The
canonical angular identity is integral (A_l')^2/wA_l=L. For l=1 the
single propagating dipole has mass squared in a^2 units
d(theta)=2(C+1)/C+2/[(1+theta)C]. Fixed monopole charge and rotational
flux are not promoted to new scalar species.

Let M0=Ma_L(0), e=(1,0)^T. The exact rank-one difference obeys

```
Ma_L(theta)-M0=-[theta/(1+theta)](L/C) e e^T,
e^T M0^-1 e=C/L.
```

Consequently (L/C) M0^-1/2 e e^T M0^-1/2 is an orthogonal rank-one
projector. For theta>-1, this proves the all-harmonic comparison

```
min(1,1/(1+theta)) M0 <= Ma_L(theta)
                           <= max(1,1/(1+theta)) M0.   (12)
```

The bound is relative even on the very light axial branch. There is no
new factor growing with C or l in this local direction. In particular,
at fixed l and large C the light mass approaches
L(L-2)/[2(1+theta)C^2]. The selected product remains positive in this
quadratic axial sector for theta>-1; this is not nonlinear stability of
the complete theory or a claim of hyperbolicity on all backgrounds.

At first order the relative perturbation norm is |theta|. The exact
selected-model remainder after its linear term is at most
theta^2/(1-|theta|) for |theta|<1. This last bound is an algebraic
remainder for (11), not a bound on the rest of the EFT.

## Comparison with the published nonlinear-electrodynamics equations

[Nomura, Yoshida and Soda](https://arxiv.org/html/2004.07560v2#S3.SS1)
derive the coupled odd-metric/even-electromagnetic master system for
general nonlinear electrodynamics, including a cosmological constant.
Their action uses -calL and invariants calF=Fbar^2/4, calG=J/4. Here
calL=calF-cJ kappa^4 calG^2, so calL_F=1, calL_FF=calL_FFF=0 and
calL_GG=-2cJ kappa^4. Their constitutive denominator is consequently
1+2cJ kappa^4 P^2/b^4=1+theta.

Their potential matrix, equations (70)-(72), multiplied by a^2 at the
product and with electromagnetic/gravitational field order exchanged,
equals (11). Equation (83) likewise gives our single dipole mass.
The checker substitutes the source expressions independently using
kappa^2 P^2/b^2=C+1 and a^2=b^2/C. The source also defines the
electromagnetic master variable with the displacement factor, consistent
with (3)-(4). This is a normalization and application check of established
perturbation theory; the magnetic-background blindness is not a novelty
claim. No exterior or driven boundary stability theorem is imported.

## Conditional coefficient range and the missing input decision

For any desired relative axial mass comparison delta>0, a sufficient
two-sided coefficient condition is

```
|theta| <= delta/(1+delta),
|cJ| <= delta/[16pi gamma(1+1/C)(1+delta)].             (13)
```

As an illustration, delta=1/100 gives |theta|<=1/101. This is an
internal spectral comparison, not an allocation of the plan's energy,
clock or area error budget. At k=2000 the retained family gives

```
C=sqrt(1+12pi k^6), gamma=6C/[k^10(C-1)],
theta/cJ=2.9452431128603523e-31,
|cJ| <=3.3616885668206477e28                           (delta=1/100).
```

The coefficient range is conditional; no such physical upper bound
has yet been established. For example, cJ=0 and the positive value
3.3953054524888543e28 give exactly the same magnetic reference solution.
The latter has theta=0.01 and changes the l=3 light axial mass by about
-0.009900990098 relative to its two-derivative value. This lies within
the one-sided positive-theta one-percent comparison and has a small
quadratic constitutive correction. It does not claim a UV completion
with that coefficient. Even cJ=1 has a nonzero axial shift, about
-2.9452431125e-31 at this k. These axial shifts are not polar self-energy
bounds and cannot be added to the earlier polar response budget.

The precise unresolved input is a finite local matching condition,
in an explicitly chosen field basis and common scheme, constraining
this direction of the action. A four-photon amplitude or a matched
parallel electric susceptibility can probe it; in the full theory these
observables constrain combinations of coefficients. Magnetic background
matching alone cannot. If a future calculation chooses cJ(mu0)=0,
that is an additional finite matching condition, and the scale-dependent
local and nonlocal terms must be combined before an observable is claimed.

The same discipline applies to curvature-squared, curvature-Maxwell,
the other electromagnetic, charged-matter and boundary operators.
Eliminating Maxwell or Einstein equations inside a source-free basis
generates charged-current and stress contacts in this model; they cannot
be discarded. The complete coefficient basis, regulator and finite
matching conditions remain unfinished. Equation (1) identifies one
independent missing direction rather than pretending to exhaust them.

## Verification and continuation

[The checker](../checks/finite_local_electromagnetic_input.py) verifies
the covariant invariant identity, background and linear variation,
electric Hessian, free photon contact, constitutive Ward identity,
six-component constrained dual inverse, harmonic mass matrix and
rank-one comparison, with an independent substitution into the published
master equations. Six finite inverse/determinant cases include three
positive Gaussian measure checks. Sixty-four independent 110-digit
mass comparisons cover four C values, four angular momenta and four
coefficient values, with three candidate substitutions. The displayed
rank-one argument establishes the continuum parameter/harmonic result;
the samples only check its implementation.

The decision is to keep S074 open and carry cJ as a physical finite
input until matched or bounded, while continuing the other action and
constraint calculations. This is not an instability, a wormhole no-go
result or an external-data impasse for all remaining research.

Remaining to-do list: fix the common local action and matching conditions,
complete nonlinear even constraints and induced contacts, include the
remaining matter/angular and boundary contributions, and control the
common measure and subtraction. Then complete S075-S080 physical state,
drive, stress/geometry, stability, transport, preparation/creation and audit.
