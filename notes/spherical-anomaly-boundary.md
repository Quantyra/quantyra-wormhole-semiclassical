# Boundary variation and transparent gluing for the spherical anomaly

S074, 2026-10-13. The spherical local anomaly admits an explicit
auxiliary-field boundary action and canonical boundary momenta. The
variations and surface actions cancel across a smooth artificial
throat--exterior cut. At a true boundary, fixing only the induced metric
and radius leaves a nonzero boundary symplectic form for the unreduced
fourth-order local block. Additional matched data or perturbative order
reduction are required. This is not a computation of the fermion
boundary anomaly or a choice of the physical driven boundary state.

## Local action and invertible auxiliary representation

Use [the accepted source conventions](spherical-anomaly-sources.md):
g2=exp(2sigma)h, X=b0^2 exp(2sigma), y=b0^-2, r=R[h],
u=Box_h sigma, v=(nabla_h sigma)^2 and I=F2^2+q^2 y^2/2.
Let N=b0^2/(4pi), a=11/360, c=1/20 and d=b-2a/3. Retain arbitrary
b for the algebra and use b=1/30 only for the displayed heat convention.
The original bulk transgression is W=N integral sqrt(-h)(L+div Bstar),
where

```
L=A r^2+b r u-3b u^2+(B+bv)r-6dv u+L0,
A=c sigma/3, B=sigma(4cy/3-4ay),
L0=(2b-4a)y v-3d v^2+sigma(4cy^2/3-2I/3).             (1)
```

Bstar is the explicit retained flux in the preceding note. Introduce
independent scalar fields P and Q, with on-shell values P=L_r and Q=L_u.
Put p=P-B-bv, z=Q+6dv and Delta=b+4c sigma. The Legendre inverse is

```
r_*=(6p+z)/Delta,
u_*=(b p-2A z)/(b Delta),
H(P,Q,sigma,v,I)=(3p^2+p z-A z^2/b)/Delta-L0,
L=P r+Q u-H on P=L_r, Q=L_u.                          (2)
```

The (r,u) Hessian has determinant -b Delta. Thus (2) is valid when
b Delta is nonzero. At b=1/30,c=1/20, sigma>=-1/100 gives
Delta>=47/1500 and the inverse Hessian infinity norm is at most
10500/47. The accepted one-percent squared-radius reference box lies
inside this range, which also includes larger exterior radii. For
sigma>=0 the two inverse-Hessian row norms are at most 210 and 30;
on [-1/100,0] the first row gives the stated bound. This is an auxiliary-coordinate bound;
it is not a bound on an evolving solution or a physical finite coupling.
Other finite EFT terms must be added with their own variations. The
degenerate cases b=0 or Delta=0 require a different representation;
they are not excluded physical theories by this calculation.

The use of auxiliary fields to expose higher-derivative boundary data
is established; [Hohm--Tonni, section 2](https://arxiv.org/html/1001.3598v2)
derive the corresponding construction for curvature-squared gravity.
Here two scalar auxiliaries suffice for this spherical anomaly block.
Their bulk algebraic equations r=H_P and u=H_Q recover (1).

## Boundary action and canonical variation

Consider a fixed smooth timelike boundary, outward spacelike unit normal
n, induced one-metric gamma and K=div n. End variations vanish near
initial/final cuts, or their corner terms are retained separately.
The auxiliary Dirichlet action is

```
W_D=N [integral_M sqrt(-h)(P r-(nabla Q).(nabla sigma)-H)
       +2 integral_boundary sqrt(-gamma) P K].        (3)
```

In the higher-derivative variables, this is the original W plus

```
W_boundary,D=N integral_boundary sqrt(-gamma)
                 [2P K-Q n.sigma-n.Bstar],            (4)
```

with P,Q evaluated by their algebraic equations. The Q n.sigma term
retains the integration by parts of Q Box sigma; Bstar removes precisely
the earlier change of representative. Neither flux is silently dropped.

Write the bulk first variation plus N integral_boundary sqrt(-gamma)
[pi_ij delta gamma^ij+pi_P delta P+pi_Q delta Q+
pi_sigma delta sigma+pi_A^i delta A_i]. The boundary momenta are

```
pi_ij=-gamma_ij n.P,   pi_P=2K,   pi_Q=-n.sigma,
pi_sigma=-n.Q-2H_v n.sigma,
pi_A^i=-4H_I n_a F^ai.                               (5)
```

The term P(K_ij-K gamma_ij) vanishes identically on a one-dimensional
boundary. On eliminating the auxiliaries, H_v=-L_v=-V and
H_I=-L_I=-Z, so (5) agrees with the preceding bulk current plus (4).
Holding gamma, P, Q, sigma and the tangential gauge potential fixed
therefore makes the boundary variation vanish. A freely varied normal
gauge component does not appear in (5).

An independent Gaussian-normal calculation fixes the signs. With
h=dx^2-exp(2omega)dt^2, r=-2(omega_xx+omega_x^2) and boundaries
x=constant, (3) reduces, after canceling the explicit GHY flux, to

```
W_D=N integral exp(omega)[2P_x omega_x-Q_x sigma_x
                         +exp(-2omega)Q_t sigma_t-H]. (6)
```

The spatial momenta of (6) are respectively 2exp(omega)P_x,
2exp(omega)omega_x, -exp(omega)sigma_x and
exp(omega)(-Q_x-2H_v sigma_x) for omega,P,Q,sigma. They reproduce
(5), including delta gamma^tt=-2gamma^tt delta omega.

## What metric and radius Dirichlet data leave unresolved

Fix gamma and sigma, but leave their normal jets arbitrary. The
remaining gravitational/scalar boundary one-form is

```
theta_restricted/N=2K delta P-(n.sigma) delta Q,
Omega_restricted/N=2 delta K wedge delta P
                         -delta(n.sigma) wedge delta Q. (7)
```

For b Delta nonzero, the independent second normal derivatives can
vary P and Q while K and n.sigma are varied independently through the
first normal derivatives. On these four jet coordinates (7) has rank
four. Adding a boundary functional shifts theta by an exact variation,
whose exterior derivative vanishes. It cannot cancel this nonzero
two-form while leaving all those jets arbitrary. Constant tangential
profiles already establish the obstruction, so tangential integrations
by parts do not remove it.

This is a precise failure of a metric/radius-only Dirichlet prescription
for the unreduced local fourth-order block. Fixing P,Q is one admissible
mathematical polarization, not two newly prepared physical particles.
For a perturbative EFT, their values must instead be inherited from
the lower-order matched solution and the declared order reduction.
The independent nonlocal functional and the microscopic boundary theory
may change the complete boundary problem. Equation (7) does not prove
that that physical problem has no consistent ensemble.

## Classical action and the electric ensemble

The classical spherical GHY term is
(2G)^-1 integral_boundary sqrt(-gamma_g) X K_g. Combining it with
the spherical Einstein--Maxwell bulk action and transforming to h gives

```
I_class,D=integral sqrt(-h)[X r/(4G)+3X v/(2G)
       +exp(2sigma)V0(X)/(4G)-pi b0^2 F2^2/g_c^2]
       +(2G)^-1 integral_boundary sqrt(-gamma) X K,
V0=2-2Lambda X-2r_m^2/X.                              (8)
```

Here r_m^2=pi G q^2/g_c^2 and derivatives/F2 use h. The n.sigma term
from K_g=exp(-sigma)(K+n.sigma) cancels the bulk Weyl integration
by parts. This reproduces the retained
[classical spherical boundary reduction](collective-clock-response.md).

For (8)+(3), replace the curvature multiplier by
Fcal=X/(4G)+N P. The full displayed boundary coefficients are

```
Pi_ij=-gamma_ij n.Fcal,   Pi_P=2N K,   Pi_Q=-N n.sigma,
Pi_sigma=X K/G+3X n.sigma/G+N(-n.Q-2H_v n.sigma),
Pi_A^i=-kappa_h n_a F^ai,
kappa_h=4pi b0^2/g_c^2+8N sigma/3.                   (9)
```

The anomaly contribution to kappa_h is the same radius--electric contact
as in the accepted bulk action. Its current must not also be added a
second time after placing it in the displacement. Other finite/local
and nonlocal electric terms remain to be matched.

Electric Dirichlet data fix A_i. The boundary Legendre transform
-integral A_i Pi_hat_A^i instead fixes the density
Pi_hat_A^i=sqrt(-gamma)Pi_A^i and leaves -A_i delta Pi_hat_A^i.
Fixing that density and fixing a contravariant current component while
varying the metric are different prescriptions. In the Gaussian-normal
chart Pi_hat_A^t=-s_n kappa_h exp(-omega)F_tx, where s_n=+1 or -1
is the outward normal sign. This also retains the corrected scaled
charge convention. Neither electric ensemble is silently assigned to
the actual preparation or driven two-boundary state.

## Smooth internal cuts carry no independent local wall

Split a smooth solution domain across an artificial timelike surface.
Use continuous induced fields and their common smooth jets, and the
same local action on both sides. The outward normals are opposite.
Consequently K, every normal derivative in (4),(5),(9), and the
electric displacement reverse sign. The surface actions and their
variations cancel exactly. The two normal momenta sum to zero, which
is the transparent matching condition in this representation.

In particular, an artificial throat--exterior cut does not authorize a
new reflecting fermion wall or an independently chosen finite wall
coupling. A discontinuity of these momenta would require an actual
interface source/action. Smooth gluing does not impose P=Q=0 or
n.sigma=0: it matches the data from the same bulk fields.

For a moving cut x=c(t) in the Gaussian-normal chart, put
D=1-exp(-2omega)c'^2>0. One outward orientation has

```
n_a=(-c',1)/sqrt(D),
K=[omega_x(1-2exp(-2omega)c'^2)
       +exp(-2omega)(c''-omega_t c')]/D^(3/2).         (10)
```

Changing sides reverses both. The numerical checks split the bulk
integral along several timelike moving cuts and include their two
oriented GHY actions explicitly. The result agrees with the unsplit
action. This is a check of the local functional and smooth gluing;
it is not a determinant-gluing theorem or a computation of a physical
scattering state.

## Remaining quantum boundary input

The variational completion (4) is not, by itself, the boundary part of
the Dirac Weyl anomaly. [Herzog--Huang--Jensen, sections 3.4 and 4.1](https://arxiv.org/html/1510.00021v3)
retain genuine boundary contributions to the anomaly and its effective
action. [Faraji Astaneh--Solodukhin, section 2.4](https://arxiv.org/html/2308.10351v2)
specify the fermion projector/domain before calculating boundary heat
coefficients. Those local reflecting conditions are not automatically
the coupled two-boundary domain of this opening problem.

The selected free boundary-angle convention remains the one in
[boundary-vertex matching](boundary-vertex-matching.md). A physical
completion must specify its full four-dimensional domain and supported
drive, calculate the true boundary determinant/counterterms in the
common regulator, and match the independent state functional. The
local gluing result removes an artificial-cut ambiguity but supplies
none of those missing quantum sources or an opening-error allocation.

## Verification and decision

The checker verifies the Legendre inverse, finite invertibility range,
canonical boundary momenta, nonzero restricted symplectic form, classical
Weyl/GHY cancellation and electric Legendre transform. Independent
variations of the second-order auxiliary action with its explicit
boundary term are compared with the bulk Euler derivatives and boundary
momenta of (6), allowing nonzero spatial boundary variations. Moving-cut
quadratures check the complete local split action. Numerical tolerances
and source hashes are fixed in the accompanying receipts.

Use transparent matching for artificial cuts in this local sector.
Carry the explicit auxiliary data or a consistently derived order
reduction at an outer boundary; do not equate a variational GHY term
with the quantum boundary anomaly. S074 remains incomplete for the
physical outer-boundary, state and other interaction requirements.

Remaining to-do list: finish S074's physical boundary/domain, independent
interactions and finite matching with a controlled remainder; S075
state/drive; S076 stress/geometry; S077 stability; S078 transport;
S079 preparation/creation; S080 final audit.
