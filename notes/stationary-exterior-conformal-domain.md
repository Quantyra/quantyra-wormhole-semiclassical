# Stationary exterior: conformal product and the actual Euclidean domain

2026-10-15, E003/S139. Ordinary argument before focused verification.
The massless four-dimensional Dirac problem on the retained magnetic RN-AdS
exterior has an exact fixed-sphere conformal description. Its physical
Euclidean boundary domain separates, and the higher angular sector has a
proper Euclidean mass gap and uniformly small curvature-to-mass ratios.
These facts remove an apparent failure of the direct-frame strong-field
expansion at infinity. They do not supply the renormalized stress on the
whole exterior, a differentiated heat remainder, the complete boundary
Weyl transgression, or a quantum-corrected exterior solution.

## Fixed classical reference and the conformal geometry

Retain the one free massless unit-charge four-component Dirac field,
the common four-dimensional subtraction and the actual paired radial
projectors of [the coupled-domain note](coupled-fermion-domain.md).
The magnetic bundle and relative spin/bundle map stay fixed. The
reference here is the identical classical RN pair used in
[S137](rn-higher-mode-boundary-spectral-control.md) and
[the complete-pulse theorem](finite-pulse-higher-mode-work.md).
No electric potential, photon/graviton loop or interaction is added.

Write m=GM, and use b for the retained classical extremal radius, a
fixed length rather than the radius of an unknown corrected horizon:

```
F(r)=1+r^2/L4^2-2m/r+r_m^2/r^2,
X=b^2/L4^2, C=1+6X, r_m^2=b^2(1+3X),
m=b(1+2X)+D/2, D=2G DeltaM>=0,
r_+=b(1+s_h), kappa0=F'(r_+)/2,
g_E=F dU_E^2+dr^2/F+r^2 dOmega_2^2.
```

The prepared range is kappa0=z_i nu, 1/4<=z_i<=7/25, with the same
laboratory source, k=2000, q=k^4 and eta=1/2000. The accepted family
also has X>=1000, C/q<7/k, s_h<1/1000 and
D<=1.001 C b s_h^2. In particular, with a_m=m/b and h_m=1+3X,

```
1+2X<=a_m<=1+2X+C/10^6.                              (1)
```

Choose the fixed conformal length b0=b; a different fixed b0 gives the
same construction with its rescaled coefficients. Set

```
u=b/r, Omega=r/b=1/u, g_E=Omega^2 ghat_E,
ghat_E=ghat_2+b^2 dOmega_2^2,
ghat_2=B(u)dU_E^2+b^2 du^2/B(u),
B(u)=u^2 F(b/u)=X+u^2-2a_m u^3+h_m u^4,
0<=u<=u_h=b/r_+<=1.                                 (2)
```

The magnetic two-form is unchanged. Its orthonormal magnitude in ghat
is Bhat_mag=q/(2b^2), constant and covariantly constant on the product.
The actual outer AdS infinity is u=0, a regular finite-distance boundary
of ghat_2, since B(0)=X>0. This is an exact massless-field change of
frame, not a replacement of the exterior by its near-horizon product.
The accepted bulk transformation and source allocation are in
[the spherical anomaly action](spherical-anomaly-action.md).

For comparison, the direct physical-frame invariants obey
|R4|/B_mag=24r^2/(qL4^2) and
|grad B_mag|^2/B_mag^3=8F/q. Both are unbounded at infinity for each
fixed q. That observation prevents a uniform direct-frame local
strong-field expansion; it is not an obstruction after (2) and the
complete conformal anomaly are retained.

## Uniform curvature and explicit geometry jets

For a two-metric as in (2), direct computation gives

```
Rhat_2=-B''(u)/b^2
      =[-2+12a_m u-12h_m u^2]/b^2.                   (3)
```

Here primes on B denote u derivatives. Its numerator P is concave.
Its minimum on [0,1] is at an endpoint, with
P(0)=-2 and P(1)=-2C+12[a_m-(1+2X)]>=-2C.
The vertex value is -2+3a_m^2/h_m. To bound it by C, observe

```
h_m(C+2)/3-(1+2X)^2=X+2X^2,
2(1+2X)(C/10^6)+(C/10^6)^2
 <=[42/10^6+49/10^12]X^2<X+2X^2,
```

using X>=1. Thus (1) implies

```
-2C<=b^2 Rhat_2<=C, |B''|<=2C.                       (4)
```

On the physical exterior B>=0, while
B<=X+1+h_m=2+4X<=C. For X>=1000, (1) also implies
a_m<2.01X and h_m<=3.01X. Differentiating the exact quartic gives

```
|B'|<=2+6a_m+4h_m<5C,
|B'''|<=12a_m+24h_m<17C,
|B''''|=24h_m<13C, B^(j)=0 for j>=5.                 (5)
```

These are coordinate jet bounds, valid through both endpoints. The
first two covariant curvature jets also have explicit regular bounds.
Let ell be proper radial distance, so the radial orthonormal derivative
is plus or minus sqrt(B) partial_u/b. For a static scalar Rhat_2,

```
|grad Rhat_2|=sqrt(B)|B'''|/b^3
                 <=17 C^(3/2)/b^3,
(Hess Rhat_2)_ell,ell
 =-[B B''''+(B'/2)B''']/b^4,
(Hess Rhat_2)_tau,tau=-B'B'''/(2b^4).
```

The signs of the displayed radial convention do not affect the bounds.
The two nonzero orthonormal Hessian components are bounded respectively
by (13+85/2)C^2/b^4 and (85/2)C^2/b^4. Hence

```
max_orthonormal_components |Hess Rhat_2|<=56 C^2/b^4. (6)
```

The formulas have regular limits at the smooth bolt. There is no
inverse kappa0 in (3)-(6). These estimates prove precisely the stated
jet orders; higher covariant jets and differentiated parametrix
constants are not silently included in (6).

At the outer boundary the outward unit normal is
nhat=-sqrt(B) partial_u/b, and the only nonzero candidate principal
curvature is the time-circle one:

```
Khat=-B'/(2b sqrt(B)), Khat(0)=0,
Rhat_2(0)=-2/b^2, Rhat_4(0)=Rhat_2(0)+2/b^2=0,
B'''(0)=-12a_m,
nhat.Rhat_4(0)=-12a_m sqrt(X)/b^3 !=0.               (7)
```

The sphere has zero extrinsic curvature. Thus the limiting boundary is
totally geodesic, but its collar is not an exact metric product in its
normal coordinate. Vanishing Khat cannot justify discarding all higher
boundary heat data.

## Physical radial projector, adjoint and the induced square

The imported microscopic projector is the eigenspace of the holographic
radial gamma matrix, rather than an independently chosen Euclidean
MIT-bag projector. See [Bintanja--Freivogel--Rolph, (2.8)-(2.12) and
appendix A.4](https://arxiv.org/html/2308.00871v2). The constant gamma
basis rotation in the existing Hamiltonian notes does not change that
domain.

Use Hermitian Euclidean gamma matrices and the anti-Hermitian bulk
Dirac expression D_E=gamma^mu nabla_mu, including the magnetic
connection. At a true outer boundary choose a normal parallel frame and
fold the identical pair. With the convention that Pi_+ is the allowed
trace,

```
chi=tau_3 gamma_n, Pi_+=(1+chi)/2, Pi_-=1-Pi_+,
Pi_- psi=0, [chi,gamma_n]=0.
```

The Green boundary form is, up to the immaterial bulk adjoint sign,
the integral of phi^dagger gamma_n psi. Since gamma_n preserves the
allowed subspace, its annihilator is the complementary trace subspace.
Thus D_E and D_E^dagger have complementary boundary domains:

```
Dom(D_E): Pi_- psi=0,
Dom(D_E^dagger): Pi_+ phi=0.                          (8)
```

In particular this D_E is not being asserted self-adjoint. The squared
operator is the positive operator D_E^dagger D_E.

The spinorial Gauss decomposition in the normal parallel frame is

```
D_E=gamma_n(nabla_n+K/2)+D_intr.
```

D_intr uses the intrinsic boundary spin connection, which keeps the
normal gamma parallel. Every tangential gamma anticommutes with chi,
so {chi,D_intr}=0. Requiring D_E psi to belong to the second domain in
(8) gives

```
Pi_- psi=0,
Pi_+(nabla_n+K/2)psi=0.                              (9)
```

There is no angular-mass Robin term: the angular Dirac expression is
part of D_intr and flips chi. In the heat convention with inward normal
n_in=-n, (9) has S=-K/2 on Pi_+, equivalently S=-L_aa/2 in the
standard mixed-boundary convention. This reconciles the outward
extrinsic-curvature sign with the inward heat formulas. At (7), S=0.

This derivation is for the static zero-angle domain. Constant paired
angles on identical exteriors are related by the retained constant
flavor rotation. A changing angle introduces derivative boundary terms;
they are not removed by (9).

## Finite-cut conformal covariance of the domain

Perform the change of frame first at a finite spherical cut, where
Omega is smooth, positive and bounded. Put sigma=log(Omega). The exact
massless Dirac covariance and spinor weight are

```
D_g=Omega^(-5/2) D_hat Omega^(3/2),
psi_g=Omega^(-3/2)psi_hat.
```

The orthonormal Clifford identification leaves Pi unchanged. With
normal-only sigma at a stationary spherical cut,

```
n_g=Omega^-1 nhat,
K_g=Omega^-1(Khat+3 nhat.sigma),
(nabla_n_g+K_g/2)psi_g
 =Omega^(-5/2)(nabla_nhat+Khat/2)psi_hat.             (10)
```

To verify the last line, the normal spin-connection Weyl correction
vanishes: its two normal sigma terms cancel. Differentiating the spinor
weight contributes -3 nhat.sigma/2, canceled by the extrinsic-curvature
term. Thus (9) transforms into the same physical mixed domain without
an extra mass Robin term.

Equation (10) concerns domains and first-order covariance. It does not
make D_g^dagger D_g and D_hat^dagger D_hat similar by a single scalar
weight, or set their determinant difference to zero. The bulk anomaly,
boundary anomaly and their finite-cut transgressions must still be
included before taking u_cut to zero. The singular limiting Omega
does not authorize discarding them.

## Exact product separation and the proper Euclidean gap

In the product frame the Clifford decomposition may be written with
the orbit gammas acting on the first factor and the sphere Dirac
operator multiplied by the orbit chirality. Its square has no cross
term. More importantly, (8)-(9) depend only on the orbit factor: the
normal gamma is gamma_n,2 tensor 1 and the sphere's extrinsic curvature
is zero. The same boundary domain therefore applies to every angular
eigenspace. Consequently, as positive squared realizations,

```
Lhat_4=Lhat_2 tensor 1+1 tensor D_S,q^2/b^2,
Lhat_2=D_hat,2^dagger D_hat,2>=0,
Lhat_4|_(n>=1)>=lambda_1^2/b^2=(q+1)/b^2.            (11)
```

The adjoint-complement condition, not an unsupported self-adjoint
first-order realization, establishes this mass-independent product
domain. The usual form closure supplies (11). The Euclidean horizon
is a smooth bolt, not an additional reflecting boundary.

The positive square and its heat kernel control the determinant modulus
and its parity-even variation. This argument does not establish that a
possible determinant phase or eta contribution vanishes in the paired
problem; no such cancellation is assumed in the global matching claim.

In the classical RN reference the regular period is beta0=2pi/kappa0.
Spinors in the rotating polar frame have the antiperiodic thermal spin
structure, with Matsubara frequencies (j+1/2)kappa0. The conformal
factor is finite and smooth at this bolt, so the period in U_E is
unchanged by (2). The gap in (11) is a proper Euclidean angular gap;
the Lorentzian Killing-frequency spectrum remains gapless at the
thermal horizon, as proved in S137.

For example the exact local heat trace factorization is

```
tr Khat_high(s;x,Omega;x,Omega)
 =tr Khat_2(s;x,x)/(4pi b^2)
       *sum_(n>=1)d_n exp[-s n(n+q)/b^2].             (12)
```

There is no separately renormalized determinant for each n. Equation
(12) is regulated before the common four-dimensional subtraction.
Its higher-angular large-time factor is explicitly bounded: with
z=exp(-q s/b^2),

```
sum_(n>=1)d_n exp[-s n(n+q)/b^2]
 <=2q z/(1-z)+4z/(1-z)^2.                            (13)
```

This follows from n(n+q)>=nq and the two geometric sums. The LLL is
excluded from (13) and must retain its nonlocal determinant. In
particular the massless infrared difficulty of an unsplit strong-field
curvature expansion is not an infrared failure of (11)-(13).

Relative to m_1=sqrt(q+1)/b, the established geometry estimates become

```
|Rhat_2|/m_1^2<=2C/(q+1)<14/k,
|grad Rhat_2|/m_1^3<=17[C/(q+1)]^(3/2),
max_components|Hess Rhat_2|/m_1^4<=56[C/(q+1)]^2.     (14)
```

They are uniform on the entire classical exterior and over the
retained thermal range. They are not error bounds for a truncated
determinant or its first metric variation.

## The nonzero boundary a4 normal flux

The general mixed Laplace-type coefficient applies directly to (9);
there is no need to replace the physical projector by a different bag
condition. Use [Vassilevich, (3.27) and (5.27)-(5.33)](https://arxiv.org/html/hep-th/0306138v3).
Here L=-(nabla^2+E), with tr E=-r_f Rhat_4/4, and r_f is the fiber
rank, four for one exterior or eight for the folded pair. The magnetic
spin term and its chi traces vanish in the traces used below.

At the strict limiting geometry (7), K=S=chi_:a=0, and Pi_+, Pi_-
have equal rank. For a smooth bounded scalar smearing f the only
surviving boundary terms of a4 are

```
f tr[(240Pi_+-120Pi_-)E_;n_in
                 +(42Pi_+-18Pi_-)Rhat_4;n_in]
 =-3 r_f f Rhat_4;n_in.
```

Indeed the two projector traces give 60r_f and 12r_f; multiplying the
first by -1/4 gives -15r_f+12r_f=-3r_f. Terms containing derivatives
of f vanish here because tr chi=tr chi E=0 and K=S=chi_:a=0. Hence

```
a4_boundary(f)=-r_f/[120(4pi)^2]
                    *integral_boundary f n_in.Rhat_4. (15)
```

This is the heat coefficient, before assigning the fermionic
effective-action/anomaly normalization. By (7) its normal curvature
factor is nonzero. The bulk derivative block in the same heat
coefficient is -r_f/[120(4pi)^2] integral f Box Rhat_4. For f=1 the
inward boundary flux in (15) exactly cancels this bulk total
derivative. For general smooth bounded f the combined derivative
block is instead

```
r_f/[120(4pi)^2] integral grad f . grad Rhat_4.        (16)
```

There is no bolt boundary flux. A further geometric check is
integral_disk Rhat_2=4pi: in smooth polar coordinates the circle radius
A(ell) has A'(bolt)=1 and A'(outer)=0, so integration of
Rhat_2=-2A''/A gives this result. Thus the orbit is not a homogeneous
negative-curvature H2 geometry even though the near-horizon limit is.

Equations (15)-(16) do not complete the boundary Weyl action. The
actual Weyl smearing sigma=-log(u) is unbounded at u=0, and the
finite-cut conformal path Omega^(2t)ghat has nonzero extrinsic
curvature for t>0. Its Euler and normal-curvature boundary terms must
be integrated along with the bulk transgression before subtraction
and the outer limit. In particular Khat(0)=0 at one endpoint cannot
be used to set the entire path's boundary anomaly to zero. The
accepted [spherical anomaly boundary note](spherical-anomaly-boundary.md)
explicitly does not yet calculate that microscopic fermion anomaly.

## A finite-coefficient direction invisible to product pressures

The selected reference prescription sets all additional finite local
coefficients to zero. That is a consistent specified input. The
larger homogeneous coefficient window in
[the finite-local reference note](finite-local-reference-window.md)
cannot by itself replace that specification on the exterior.

An explicit example in the retained, unredefined curvature basis is

```
Delta L=c(R4^2-2 Ric4^2), (c1,c2,c3)=(c,-2c,0),
other displayed c_i=0.
```

Its homogeneous combination A_g=4c1+2c2+4c3 and both magnetic
combinations K,H,c4 vanish, so its two product pressures and the
displayed seminorm N4 vanish for arbitrary c. On a direct product of
two surfaces this density equals the Euler density; on a general
four-metric it is not the Euler density. No source-free field
redefinition is being made. The standard curvature-square span is
also displayed in [Cheung--Liu--Remmen, (4)](https://arxiv.org/html/1801.08546v3);
their positivity assumptions are not imported.

The exterior effect can be derived with compact variations, avoiding
any boundary-action ambiguity. On the Lorentzian RN ansatz with
alpha=0 define

```
A=-F''/2-F'/r, B_r=(1-F-rF')/r^2,
R4=2(A+B_r), Ric4^2=2A^2+2B_r^2,
r^2(R4^2-2Ric4^2)=8A(1-F-rF').
```

Vary F with compact support. If Q=1-F-rF', its Euler derivative is

```
-8A+8(Q/r+Ar)'-4Q''
 =-4[F''+2(1-F)/r^2]
 =-16 r_m^2/r^4                                    (17)
```

on RN-AdS. Since the four-dimensional volume at alpha=0 is
4pi r^2 dr dU and
delta Gamma=4pi integral r^2(T^U_U-T^r_r)deltaF/(2F),
equation (17) implies

```
Delta T^r_r-Delta T^U_U=32c F r_m^2/r^6,
alpha'_c=4piG r(Delta T^r_r-Delta T^U_U)/F
        =128piG c r_m^2/r^5,
alpha_c(r)=-32piG c r_m^2/r^4, alpha_c(infinity)=0.   (18)
```

This is an order-reduced contribution, evaluated on the classical
background at first order in c. It is an independent ordinary
variation, not an assigned physical coefficient or a numerical
result. Thus an exterior lapse bound cannot follow from the product
seminorm alone. The current zero-extra-coefficient prescription sets
this c to zero explicitly; nothing in (18) contradicts that choice.

## What remains before a global static matching theorem

The closest primary constructions identify a viable method, not a
ready error bound for this problem. [Avramidi--Fucci](https://arxiv.org/abs/0810.4889)
construct local nonperturbative constant-field heat asymptotics.
[Fucci, section VII.2](https://arxiv.org/html/0906.2430v1) displays
infrared-singular massless spinor curvature coefficients when the LLL
is not separated. [Ribeiro--Shapiro--Zanusso, section 5](https://arxiv.org/html/1803.06948v2)
give massive two-dimensional Dirac curvature form factors through
curvature squared, not an all-orders remainder with this boundary.
[Grieser, section 3](https://web.math.ku.dk/grubb/notes/heat.pdf) constructs
boundary heat kernels by a Volterra parametrix; its existence and
asymptotic classes do not supply explicit constants uniform in this
near-extremal family and its metric variations.

The actual static-stress alternative is a common-regulator Euclidean
radial Green-function calculation, with high-mode subtraction and
explicit tails, following the kind of construction of
[Groves--Anderson--Carlson](https://arxiv.org/html/gr-qc/0207066v1).
Their uncharged-field framework is not already a magnetic RN-AdS
error theorem. The related [black-hole stress calculation](https://arxiv.org/html/gr-qc/0305045v2)
shows why its analytic approximation alone is inadequate: it can
give an incorrect horizon sign and spurious RN horizon divergences.

For this reference the unresolved quantitative object is the first
and second metric/radius variation of the *subtracted* product-frame
heat kernel, with (9), the regular HH bolt, the actual boundary layer,
and the common four-dimensional angular subtraction. A bound on a
positive heat moment, a small invariant ratio in (14), or an unsigned
determinant remainder cannot be differentiated into such a bound.
The needed estimate must also retain the regular combination
(T^r_r-T^U_U)/F that enters the lapse constraint; an absolute tensor
norm alone does not prove its horizon integrability.

The stationary solution must be centered consistently. The classical
nonextremal mass defect has D/(Cb)=O[(kappa0 b/C)^2], of order
eta^2 k^-18 in the retained family, whereas the proved stationary
product translation is of order k^-8. One cannot keep the old mass
and horizon fixed and infer HH regularity from a small ordinary
fractional correction. First determine the corrected extremal center;
then impose the prepared physical temperature
kappa_phys=exp(alpha_h)F'(r_h)/2 with alpha(infinity)=0.
The accepted homogeneous product alone fixes neither alpha_h nor
the exterior-normalized clock charge.

The outer AdS center also requires a renormalization condition. If the
full state approaches the AdS-invariant vacuum with inverse squared
radius y=L_eff^-2, the retained anomaly gives C4^2=Box R4=0 and
E4=24y^2 as the local magnetic field tends to zero. Consequently

```
T^mu_mu=-11y^2/(240pi^2),
T_mu_nu=-11y^2 g_mu_nu/(960pi^2),
y+[11G/(360pi)]y^2=L4^-2,                           (19)
```

where the last line uses Lambda_ren=-3/L4^2 and no additional constant
vacuum subtraction. The trace and AdS-invariant tensor agree with
[Ambrus--Winstanley, (22), (39) and (54)](https://arxiv.org/html/1505.04962v2).
One must either center the outer radius as in (19) or explicitly absorb
that vacuum stress into the cosmological renormalization condition.
The local decay of the magnetic field does not prove that the actual
monopole HH state has this leading invariant limit. That limit and its
subleading integrable remainder remain matching requirements; the
complete omitted functional cannot simply be declared an integrable
mass source on the unchanged bare-L4 background. This conditional
observation does not alter the declared isolated LLL theorem.

The companion [isolated Polyakov-LLL exterior](stationary-lll-quantum-exterior.md)
constructs a recentered global stationary solution for its explicitly
allocated sector, including the radius source, regular HH state and
laboratory temperature. It does not restore the remaining angular or
four-dimensional Weyl functional, so the full stress and matching
requirements above remain open.

Remaining to-do list: complete the finite-cut boundary Weyl action and
charge prescription; derive the differentiated global stress/radius
remainder and solve the corrected stationary exterior/state matching;
transfer S138, control local driven response and complete the remaining
S074-S079 requirements.
