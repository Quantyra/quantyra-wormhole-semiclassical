# Bosonic dipole source, angular Gauss flux and correlated quadrupole response

S072, 2026-10-13. The electromagnetic dipole produces both sources of
the coupled polar quadrupole. They are correlated quadratic operators,
not a prescription for the electric energy density alone. This note
derives their coefficients from the four-dimensional equations and gives
their complete leading Gaussian ordered kernel on the product throat.
It also solves the round-projection bosonic Gauss law in radial gauge,
exhibiting the cubic and quartic interactions omitted by that Gaussian
kernel. Thus the approximation and its missing interaction are explicit.

This is the bosonic contribution around the classical product with zero
background rotation flux. Fermion polarization and its mixed correlations,
induced anomaly/Hall response, other massive modes, physical boundaries
and an interacting state remain to be added. The round projection is not
a consistent nonlinear truncation; the quadrupole calculated here is one
of its required additional fields. No capture or formation gate closes.

## Canonical dipole and the exact projected bosonic Gauss constraint

Use the action and conventions of [the all-axis constraint calculation](nonabelian-dipole-constraints.md).
On the fixed product let b be the sphere radius, a the AdS2 radius,
`C=b^2/a^2`, `Y=(C+1)/2`, and `Qmag^2=b^2 Y`. The magnetic flux is
4pi P, P=q/2; the canonical dual electric background is
`E0=P/(g b^2)`, with `kappa^2 E0^2=2Y/b^2`, `kappa^2=8pi G`.
Keep distinct the gravitational constant G and horizontal Maxwell
curvature Gcal. The retained action coefficients are

```
kappag=b^4/(6G), kappae=4pi b^2/(3g^2), mu=8pi/(3g^2).
```

In optical coordinates `ds2^2=f(-dt^2+dx^2)`, choose radial gauge A_x=0
and zero inner rotation flux. Temporarily remove fermionic sources.
The canonical transformation

```
phi=E_H/sqrt(mu),       pi=-sqrt(mu) H_x,
rho_bos=H_x cross E_H=phi cross pi,
K(x)=integral_inner^x rho_bos(y) dy,       E_A(x)=K(x)    (1)
```

solves the full bosonic Gauss constraint of the round projection. It
does not replace three noncommuting currents with three independent
charge sectors. The inner boundary choice is stated; a nonzero inner
flux adds to K. Outer rotations still generate the boundary charge.
Keep the boundary terms of the earlier canonical reduction.

Substitution gives the projected Hamiltonian

```
H = H2+H3+H4,
H2 = (1/2) integral [pi^2+(partial_x phi)^2+f m_d^2 phi^2] dx,
m_d^2 = mu(1/kappae+P^2/kappag)=(2+4Y)/b^2,
H3 = -integral [f P sqrt(mu)/kappag] phi dot K dx,
H4 = integral [f/(2kappag)] K^2 dx.                    (2)
```

Equation (2) is an exact classical rewriting of the bosonic round
projection, not a proof of a renormalized continuum quantum Hamiltonian.
At a vector-symmetry-preserving finite regulator, rho_bos is the usual
rotation generator and satisfies its su(2) algebra. The summed
commutator of phi_a with K_a is zero, fixing the displayed cubic
ordering. Products and the vacuum subtraction in H4 still require
continuum renormalization. Classical square positivity does not provide
a normal-ordered quantum error bound.

At leading order K is quadratic, so the three components of phi obey
`(Box2-m_d^2)phi_a=0`, with
`m_d^2 a^2=2+4/C`. This free triplet is the H2 term of (2). H3,H4, the
omitted quadrupole and fermion coupling are not absent from the physical
system. In particular it is not its already-dressed thermal ground state.

## Direct four-dimensional quadratic source

One axis fixes the symmetric spin-two coefficients by polarization.
Write c=cos(theta), B=1+xi^2/a^2 and use

```
ds4^2=-B dT^2+d xi^2/B
       +b^2[dc^2/(1-c^2)+(1-c^2)(dvarphi-w dT)^2],
F4=-P dc wedge dvarphi+d[c(Pw+h)dT+c v dxi],
Gcal_Txi=partial_T v-partial_xi(Pw+h),
F_A,Txi=-partial_xi w,
H^2=-h^2/B+B v^2.                                    (3)
```

Count w,h,v as first order and include the quadratic w^2 term in the
metric. The source is minus the quadratic residual of the *full*
Einstein-Maxwell equations. Use its horizontal base components, not
coordinate components that still include the fibre shift. With the
unnormalized harmonic P2(c), the quadratic effective base trace and
scalar traceless sphere coefficient are

```
(S^a_a)_l2 = -2 kappa^2 Gcal_Txi^2/(3g^2)
                           +b^2 F_A,Txi^2/3,
S_T = (2/3)[kappa^2 H^2/g^2-b^4 F_A,Txi^2/2].          (4)
```

For example the full angular base trace before projection is
`-kappa^2 c^2 Gcal_Txi^2/g^2-b^2(1-c^2)F_A,Txi^2/2`.
The second term includes the gravitational connection curvature, which
would be lost in an electromagnetic-stress-only estimate. The sphere
tensor is proportional to
`D_i c D_j c-(1-c^2)gamma_ij/2`. Its traceless Hessian relation to P2
fixes the coefficient S_T in the conventions of the
[sourced product equations](quadrupole-product-response.md).

The horizontal azimuthal Maxwell residual is
`c (partial_xi w) Gcal_Txi/g^2` in the action-normalized Maxwell equation
(or without g^(-2) in div F). The first-order t and xi Maxwell equations
vanish on the dipole solution. For a second-order original magnetic
axial potential `delta A_varphi=-a2(1-c^2)P2'(c)`, the dual polar master
is Acal=a2/g. Thus both complete sources are

```
r_F = -kappa^2 Gcal_Txi^2/(3g^2)
       -kappa^2 H^2/(b^2 g^2)+(2b^2/3)F_A,Txi^2,
r_A = E0 S_T/2 - b^2 F_A,Txi Gcal_Txi/(3g).            (5)
```

The E0 S_T/2 term is essential. The last term alone is only the extra
Maxwell residual, not the full master source. This also avoids dropping
the effective current terms that appear after electromagnetic duality.

On the first-order solution with zero background rotation flux,

```
Gcal_Txi=sqrt(mu) phi/kappae,
F_A,Txi=-P sqrt(mu) phi/kappag,
h=-B partial_xi phi/sqrt(mu),
v=-partial_T phi/[B sqrt(mu)],
H^2=-(D phi)^2/mu.                                    (6)
```

The one-axis P2 coefficients in (5) reduce to

```
r_F = (3G/b^2)(D phi)^2+[4G(4Y-1)/b^4]phi^2,
r_A = -G E0(D phi)^2+[4G E0(1-Y)/b^2]phi^2.            (7)
```

The direct checker derives the metric, inverse, connection, Ricci,
Maxwell tensor and stress through second order, retaining arbitrary
functions w,h,v before the dipole equations are imposed. It verifies
(4),(5), the first-order Maxwell equations and the linear polar Maxwell
equation including its base trace. The units a=1, b^2=C, kappa^2=2,
g=1 retain arbitrary C and its exact background charge relation. The
restoration of (7) is checked symbolically with independent b,G,g,P.

## All axes and the common composite operator

Define the traceless quadratic tensor and its normalized quadrupoles by

```
Q_ab=:phi_a phi_b:-delta_ab :phi dot phi:/3,
T_m,ab=integral dOmega Y_2m(n)(n_a n_b-delta_ab/3),
Q_m=T_m,ab :phi_a phi_b:,
Tr(T_m T_m')=(8pi/15)delta_mm'.                        (8)
```

Real orthonormal sphere harmonics are used. On one axis,
`(n dot phi)^2-phi^2/3=(2/3)phi^2 P2(c)`; hence multiplying (7) by
3/2 and polarizing gives the normalized all-axis source. Antisymmetric
products of dipole vectors are spin one and enter the Gauss sector;
they do not supply an additional symmetric spin-two coefficient.

Use `Box Q_m=2:(D phi D phi)_m:+2m_d^2 Q_m`. At this free order a
rotation-invariant Wick prescription has only delta_ab c-number
ambiguities in the quadratic field product, removed by its traceless
projection. At separated points the complete sources are

```
r_F,m = D_F Q_m,
D_F = [9G/(4b^2)]Box2 + 3G(C-4)/b^4,
r_A,m = D_A Q_m,
D_A = -(3G E0/4)Box2 + 9G E0/b^2.                    (9)
```

Their mixed correlation is fixed by the same Q_m. A zero mean in a
rotation-invariant state does not set either source variance to zero.
Local response contacts in a varied geometry and interaction-dependent
renormalization are separate from this fixed-background Wick kernel.

There is a useful exact field rearrangement at this order. Let
u=(chi,Acal), and write its mass matrix as K, so that the product
equations read `(Box2 I-K)u=r`. Define

```
v=(9G/(4b^2),-3G E0/4)^T,    u=utilde+v Q_m.
K v=(9G/b^4,0)^T,
(Box2 I-K)utilde
    =(3G(C-1)/b^4,9G E0/b^2)^T Q_m.                  (10)
```

This removes derivatives from the forcing, while retaining the local
composite contribution v Q_m in the reconstructed fields. It does not
erase that contribution or its fluctuations. Initial data transform
with it. If a switching function multiplies the interaction, derivatives
of that function introduce the corresponding commutator with Box2;
one cannot use (10) to infer zero switching radiation. The very light
quadrupole branch remains in the retarded response to the right side.

## Ordered vacuum kernel, including the Gauss charge

Specify global AdS2 with the Friedrichs ends and H2 vacuum. In optical
coordinates `ds2^2=a^2 sec^2(sigma)(-d tau^2+d sigma^2)`, let u_n be
the normalized Gegenbauer modes already used in the product response,
now with `Delta_d=(1+sqrt(9+16/C))/2`, `omega_n=n+Delta_d`. The
canonical dipole Wightman function is

```
<phi_a(X)phi_b(X')>=delta_ab W(X,X'),
W=sum_(n>=0) u_n(sigma)u_n(sigma')
          exp[-i omega_n(tau-tau'-i0)]/(2omega_n).
```

Wick contraction gives the entire two by two source matrix

```
<r_I,m(X) r_J,m'(X')>_conn
 = (16pi/15)delta_mm' D_I,X D_J,X' [W(X,X')^2],
I,J in {F,A}.                                        (11)
```

This is an ordered distribution. Its reverse ordering is not to be
substituted into (11); response uses their commutator and its local
contacts. The stress is composite and need not have Gaussian higher
moments merely because phi is Gaussian. The global vacuum pair
frequencies are omega_p+omega_q>=2Delta_d. Finite temperature adds
other spectral processes. Neither this vacuum gap nor the global
boundary domain is assigned to a physical horizon.

The Gauss density in (1) is nonzero as an operator even in this invariant
vacuum. Its ordered kernel is

```
<rho_a(X)rho_b(X')>
 =2delta_ab[W partial_tau partial_tau' W
                -(partial_tau W)(partial_tau' W)]
 =delta_ab sum_pq [(omega_p-omega_q)^2/(4omega_p omega_q)]
   u_p(sigma)u_q(sigma)u_p(sigma')u_q(sigma')
   exp[-i(omega_p+omega_q)(tau-tau'-i0)].               (12)
```

The common frequency vanishes from the numerator: equal-frequency
pairs carry no rotation charge. Integrating over the complete spatial
slice uses mode orthogonality and annihilates the vacuum global charge.
A partial-slice flux K generally fluctuates. This distinguishes a
global singlet from vanishing local gauge flux. The free Q-rho two-point
function vanishes by its symmetric/antisymmetric index contraction;
the F-A source correlation in (11) generally does not.

For a cut z put `I_pq(z)=integral_(-pi/2)^z u_p u_q d sigma`. Smear K
in time with a normalized Gaussian of standard deviation s>0. Its
per-axis variance is the positive series

```
V_K(z;s)=sum_pq [(omega_p-omega_q)^2/(4omega_p omega_q)]
              I_pq(z)^2 exp[-s^2(omega_p+omega_q)^2]
 <= sum_(m>=1) [(m+1)m^2/(4Delta_d^2)]
                         exp[-s^2(m+2Delta_d)^2].      (13)
```

Cauchy-Schwarz gives |I_pq|<=1. Group p+q=m and use
`|omega_p-omega_q|<=m`, `omega_p omega_q>=Delta_d^2` to prove the
cut-uniform bound. For modes p,q<N, the discarded part has m>=N;
when `(N-1)^2>3/(2s^2)`, it is bounded by

```
exp[-s^2 R^2] [R^2/s^2+1/s^4]/(4Delta_d^2), R=N-1.     (14)
```

The decreasing Gaussian-polynomial integral proves this bound. It
controls the spectral tail, not a separate quadrature error in evaluating
the retained overlaps. Time smearing is part of the observable; an
unsmeared point-flux or stress bound is not claimed. Equations (11)-(14)
are Gaussian reference kernels in the explicit Gauss reduction (1)-(2),
not estimates on the corrections generated by H3 and H4.

## Verification and continuation

[The checker](../checks/dipole_quadrupole_source.py) combines the direct
second-order tensor derivation with the canonical and source identities.
Independent bosonic two-particle states verify the full ordered source
matrix and Gauss kernel, including mixed F-A entries. Spatial cuts,
mode sums and their analytic tail comparison check the smeared Gauss
variance. The [receipt](../evidence/dipole-quadrupole-source.json) records
the exact source/dependency hashes and numerical scope.

The initial tensor comparison failed because its two sides applied the
background charge identity at different stages. Its residual was exactly
proportional to `C^2+C-2P^2` in the stated units. Applying the identity
to both sides made the independent tensor comparison pass, without
changing a physical coefficient or tolerance. That failed snapshot and
log are preserved. An intermediate development kernel retained only
the extra Maxwell residual; the final complete r_A includes E0 S_T/2
from the separately checked sphere/Maxwell equations. The intermediate
PASS therefore does not certify the final physical source by itself.

This is a derivation within the prior action and master conventions, not
a novelty claim or a full nonlinear reduction. The next required step
is to use the joint fermion/dipole theory, its induced response and
physical state to control the difference from (11), then propagate and
reconstruct physical observables with the correct boundaries. H3,H4
identify concrete interaction terms for that comparison. Fermionic
polarization, other harmonics and massive levels cannot be absorbed
into the bare dipole mass by assumption.

Remaining to-do list: control the joint interacting source and induced
response, including H3/H4 and mixed fermion correlations; match physical
boundaries and preparation; obtain evolving renormalized stress and
geometry, sufficient stability and finite-energy transport.
