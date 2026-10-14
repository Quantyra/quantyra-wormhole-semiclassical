# Axisymmetric axial action and the mixed dipole-octupole vertex

S074, 2026-10-13. The missing l=1 to l=3 overlap in
[the local contact audit](polar-local-contact-audit.md) now has a physical
classical bulk vertex. A spacelike reduction keeps all axisymmetric axial
harmonics before dualizing their Gauss constraints. It yields the full
quadratic axial action on a prescribed even background, including its
two-polar/two-axial contact generator. At cubic order its canonically
normalized l=1/l=3 vertex has an exact inverse-polar-mass norm. The light
l=3 branch does not cause an unbounded relative coupling as C grows.

This does not yet bound the mixed loop. The nonlinear even constraints,
other angular/parity and fermion sectors, boundary ensemble, quantum
measure and common renormalization remain unmatched. S074 stays open.

## Four-dimensional origin and local dualization

Use signature (-+++) and the Einstein-Maxwell bulk action

```
S4 = integral sqrt(-g4) [(R4-2Lambda)/(2kappa^2)-F4^2/(4g^2)],
kappa^2=8pi G, C=b^2/a^2>=1,
kappa^2 P^2/(g^2 b^2)=C+1.                              (1)
```

Here P is the magnetic potential coefficient and phi has period 2pi.
On an axisymmetric patch away from the fixed points of the phi action,
introduce a spacelike fiber and three-dimensional Lorentzian base:

```
ds4^2=gamma_MN dx^M dx^N+R^2(dphi+B)^2,
A4=A+sigma(dphi+B),
F4=dA+sigma dB+d sigma wedge(dphi+B),                    (2)
R4=R3-R^2(dB)^2/4-2 R^(-1) Laplacian_gamma R.
```

The conventional reduction/dualization method is reviewed in
[Klemm, Nozawa and Rabbiosi, section 2.1, equations 3-9](https://arxiv.org/html/1506.09017v2#S2.SS1).
Their reduction is timelike and uses different Maxwell normalization;
its signs and integrability restrictions are not imported here. The
checker reconstructs (2) from four-dimensional metric jets with all
three components of dB and arbitrary first/second angular derivatives of
the even diagonal metric. It separately contracts the Maxwell tensor.

After integrating phi the local bulk action is

```
S3=2pi integral sqrt(-gamma) [
 R(R3-2Lambda)/(2kappa^2)-Laplacian_gamma R/kappa^2
 -(d sigma)^2/(2g^2 R)] + Sodd,
Sodd=-1/4 integral sqrt(-gamma) F_I K_IJ F_J,
F=(dB,dA),
K=2pi R [[R^2/(2kappa^2)+sigma^2/g^2, sigma/g^2],
          [sigma/g^2,                       1/g^2]].     (3)
```

The Laplacian term is a boundary term, not an omitted interaction in the
interior. The four-dimensional gravitational boundary term must also be
carried to a boundary comparison. The local odd sector is exactly
quadratic for prescribed gamma,R,sigma; no harmonic truncation was used.

Introduce two multipliers Phi=(PhiB,PhiA) by adding integral Phi_I dF_I.
Integrating by parts gives minus integral dPhi_I wedge F_I plus the
boundary Legendre term integral_boundary Phi_I F_I. Eliminating F gives

```
Sodd=-1/2 integral sqrt(-gamma) (K^-1)_IJ dPhiI dot dPhiJ,
K^-1=1/(2pi R) [[2kappa^2/R^2, -2kappa^2 sigma/R^2],
 [-2kappa^2 sigma/R^2, g^2+2kappa^2 sigma^2/R^2]].        (4)
```

The checker verifies the inverse and the Lorentzian first-order action
sign, rather than substituting a duality equation into the original
Maxwell action alone. This is a classical local Gauss reduction. Its
functional measure and determinant are not asserted to be unchanged.

Regularity and fixed flux select the scalar sector used below. A magnetic
patch change sigma -> sigma+q entails A -> A-q B and
PhiB -> PhiB+q PhiA; the differential combination in (4) is invariant.
Fixed magnetic/electric monopoles and the nonpropagating rotational flux
are not added as fluctuating scalar species. The scalar potentials can
have additive constants, fixed in this comparison. On the sphere axes
the regular mode functions below vanish at the required orders.

The bulk integrations by parts below hold for compactly supported test
fields. They extend to finite sums of free global Friedrichs modes when
their boundary terms vanish. This specifies the current free comparison,
not the physical driven AdS boundary ensemble. An electric-potential or
fixed-charge boundary condition must be transformed with the Legendre
term, and cannot be silently interchanged.

## Exact angular action and its contact generator

Let c=cos(theta), s=1-c^2, Y=P2(c), V=s Y'=3cs. Choose the even
angular gauge and coordinates

```
gamma_ab=h_ab(x,c), gamma_cc=b^2 exp(2rho)/s, gamma_ac=0,
R=b exp(rho) sqrt(s),
sigma=-P[c+alpha(x)V(c)], rho=chi(x)Y(c),
X=PhiB+P c PhiA, alpha=Acal/(b^2 E0), E0=P/(g b^2).     (5)
```

In expressions using Acal below set g=1, as in the retained Maxwell
normalization. The minus sign in sigma matches the magnetic convention
of [the polar canonical action](polar-canonical-vacuum.md):
delta A_phi=-Acal V. Equation (4) becomes

```
Sodd=-1/(4pi) integral d^2x dc sqrt(-h) {
 h^ab [g^2 partial_a PhiA partial_b PhiA/s
 +2kappa^2 (partial_a X+P alpha V partial_a PhiA)
             (partial_b X+P alpha V partial_b PhiA)
            /(b^2 exp(2rho) s^2)]
 +g^2 (PhiA')^2/(b^2 exp(2rho))
 +2kappa^2 (X'-P PhiA+P alpha V PhiA')^2
            /(b^4 exp(4rho) s) }.                      (6)
```

Primes in (6) mean partial_c. In particular,
dPhiB-sigma dPhiA=dX-P PhiA dc+P alpha V dPhiA.
No derivative of alpha is lost: (6) is in dual coordinates, and the
vector-to-scalar relation itself depends on the background.

Equation (6) generates every prescribed-background cubic and quartic
vertex with two axisymmetric axial fields. It is not a truncation to
l=1. For arbitrary harmonics replace the rows used below by the complete
harmonic rows. Second-order even fields determined by their own
constraints have not been substituted into (5).

To make the contacts explicit set a=1, b^2=C, g=1,
lambda=kappa^2/b^2 and P=sqrt((C+1)/lambda). Let row vectors A(c),X(c)
map a column of canonically normalized axial fields q to PhiA and X.
Define

```
D=X'-P A, Z=P V A, E=P V A',
T0=A^T A/s+2lambda X^T X/s^2,
T1=2lambda [alpha(X^T Z+Z^T X)-2chi Y X^T X]/s^2,
T2=2lambda [alpha^2 Z^T Z-2chi alpha Y(X^T Z+Z^T X)
                     +2chi^2 Y^2 X^T X]/s^2,
U0=A'^T A'/C+2lambda D^T D/(C s),
U1=-2chi Y A'^T A'/C
    +2lambda [alpha(D^T E+E^T D)-4chi Y D^T D]/(C s),
U2=2chi^2 Y^2 A'^T A'/C
    +2lambda [alpha^2 E^T E-4chi alpha Y(D^T E+E^T D)
                     +8chi^2 Y^2 D^T D]/(C s).         (7)
```

Write h_ab=g_ab+Y H_ab with tracefree H. At an orthonormal base point
H=[[p,qh],[qh,p]], J=p^2-qh^2. Raising H uses g. Through second
external order the two-dimensional kinetic and potential matrices are

```
K^ab=(1/(2pi)) integral dc {
 g^ab[T0+T1+T2+(J/2)Y^2 T0]-H^ab Y(T0+T1)},
M=(1/(2pi)) integral dc {U0+U1+U2-(J/2)Y^2 U0},
Sodd=-1/2 integral sqrt(-g) [K^ab_ij partial_a q_i
                             partial_b q_j+M_ij q_iq_j]. (8)
```

These follow from sqrt(-h)/sqrt(-g)=1-JY^2/2 and
sqrt(-h)h^ab/sqrt(-g)=g^ab-YH^ab+JY^2 g^ab/2+O(H^3).
All cross terms in chi,alpha,H are retained. The receipt prints every
matrix coefficient for q=(q1,q3A,q3X), including the quartic contacts.
Direct quadrature of the unexpanded (6) checks the expected cubic
remainder when chi, alpha and both independent components of H
are scaled together. This is a Taylor consistency check,
not a bound for a fluctuating quantum background.

## Canonical harmonics and both axial mass branches

For L=l(l+1), put

```
A_l=s P_l', X_l=s^2 P_l'',
A_l'=-L P_l, X_l'=-(L-2) A_l,
wA_l=integral A_l^2/s dc=2L/(2l+1),
wX_l=integral X_l^2/s^2 dc=wA_l(L-2),
PhiA=sum_l a_l A_l, X=sum_(l>=2) xi_l X_l,
qA_l=sqrt(g^2 wA_l/(2pi)) a_l,
qX_l=sqrt(kappa^2 wX_l/(pi b^2)) xi_l.                 (9)
```

The norm identities follow from the Legendre equation and integration
by parts; the same argument gives orthogonality. The checker evaluates
them exactly for l=1 through 8. In these coordinates the l>=2
dimensionless axial mass matrix is

```
Ma_l=(1/C) [[L+2(C+1), sqrt(2(C+1)(L-2))],
            [sqrt(2(C+1)(L-2)), L-2]],
tr Ma_l=2+2L/C, det Ma_l=L(L-2)/C^2,
mu_plus =(C+L+sqrt(C^2+2L(C+1)))/C,
mu_minus=L(L-2)/[C(C+L+sqrt(C^2+2L(C+1)))].             (10)
```

This independently reproduces the product masses in
[the sourced response note](quadrupole-product-response.md). For l=1
there is no X_1 mode: the fixed rotational flux does not supply a second
propagating scalar. The dipole mass is d=2+4/C. For l=3 the lower mass
is 120/[C(C+12+sqrt(C^2+24(C+1)))], asymptotic to 60/C^2.

For (q1,q3A,q3X), define the cubic matrices by

```
S3=-1/2 integral sqrt(-g) {
 [chi Kchi+alpha Kalpha]_ij grad q_i dot grad q_j
 +Kh_ij H^ab partial_a q_i partial_b q_j
 +[chi Vchi+alpha Valpha]_ij q_i q_j }.
```

Exact integration of (7) gives Kchi=0 in this selected submatrix and

```
Kalpha=[[0,0,6sqrt(70(C+1))/35],
        [0,0,2sqrt(5(C+1))/5],
        [6sqrt(70(C+1))/35,2sqrt(5(C+1))/5,0]],
Kh=[[1/5,-3sqrt(14)/35,0],[-3sqrt(14)/35,-1/5,0],[0,0,0]],
Vchi=[[8/5,-12sqrt(14)(2C+5)/(35C),-24sqrt(70(C+1))/(35C)],
      [-12sqrt(14)(2C+5)/(35C),-8/5-8/C,-8sqrt(5(C+1))/(5C)],
      [-24sqrt(70(C+1))/(35C),-8sqrt(5(C+1))/(5C),-8/C]],
Valpha=[[24(C+1)/(5C),-12sqrt(14)(C+1)/(35C),24sqrt(70(C+1))/(35C)],
        [-12sqrt(14)(C+1)/(35C),16(C+1)/(5C),8sqrt(5(C+1))/(5C)],
        [24sqrt(70(C+1))/(35C),8sqrt(5(C+1))/(5C),0]].    (11)
```

Thus the nonzero angular overlap does survive the local axial Gauss
reduction. It couples the dipole to both physical l=3 branches after
diagonalizing Ma_3. This result is more specific than an overlap of
unconstrained scalar harmonics.

## Free on-shell representative of the mixed cubic interaction

At cubic order only the linear even constraints are needed for the
on-shell vertex. Nonlinear elimination still changes the fourth-order
action. From the source-free quadrupole reconstruction,

```
H_ab=(2C/3)[nabla_a nabla_b chi-g_ab Box chi/2].         (12)
```

This follows from equation (6) of the sourced response note with L=6;
it is used only for free polar solutions, in a=1 units. If f,g are free
axial solutions with mass squared d,nu, respectively, integration by
parts gives

```
integral chi grad f dot grad g
 =1/2 integral [Box chi-(d+nu)chi] f g,
integral (HessTF chi)_ab grad^a f grad^b g
 =1/4 integral [(d+nu)Box chi-(d-nu)^2 chi] f g.         (13)
```

For the second identity set
T_ab=grad_(a f grad_b) g-g_ab grad f dot grad g/2.
Its divergence is [d f grad_b g+nu g grad_b f]/2; a second divergence
gives (13). No equality of the two masses is assumed. The checker
independently verifies that divergence with arbitrary radial jets in
ds2^2=-(1+x^2)dt^2+dx^2/(1+x^2), two arbitrary frequencies and two
arbitrary masses. Boundary contributions are retained as a qualification.

Let vc,va,h,k be the off-diagonal row (columns q3A,q3X) of
Vchi,Valpha,Kh,Kalpha. For a free external polar mass z and a free
axial branch nu, the mixed source row is

```
Fc(z,nu)=vc+(C/6) h [(d+nu)z-(d-nu)^2],
Fa(z,nu)=va+(1/2) k [z-d-nu].                          (14)
```

The action coefficient of q1 q3branch has no factor 1/2 because the
two internal fields are distinct. This step selects a cubic
nondifferentiated representative. Replacing equations of motion in
the action requires transformed sources, boundary terms and generated
quartic terms for any higher-order comparison; see the retained
[field-redefinition account](polar-dipole-dressing.md). The quartic
generator (7)-(8) is in the original dual coordinates, before those
additional transformations, and is not claimed to include them.

There is a useful independent regression. The q1^2 vertex from (11),
(12) and the external master equation gives the coefficients

```
chi:   -2(C^2+11C+6)/(15C),
alpha: -2(C+1)(C+8)/(5C).                              (15)
```

Converting alpha=Acal/P agrees exactly with -(8pi/15) W R from the
retained dipole source calculation. This verifies the magnetic sign,
canonical internal normalization and on-shell diagonal overlap.
Off-shell quartics need not equal those obtained by truncating the
vector action to l=1 before dualization.

## Exact mixed coupling norm, including the light branch

The external quadrupole still uses its physical kinetic matrix, not
unit-normalized response variables. Define

```
Mp=[[2+6/C,sqrt(12(C+1))/C],[sqrt(12(C+1))/C,6/C]],
Mp ep=mu_p ep, Ma_3 es=nu_s es, N=sqrt(4pi/5),
chi_p=sqrt(3lambda) ep_1/[N sqrt(C mu_p)],
alpha_p=sqrt(lambda/(C+1)) ep_2/[N sqrt(C mu_p)],
g_ps=chi_p Fc(mu_p,nu_s) es+alpha_p Fa(mu_p,nu_s) es.
                                                               (16)
```

N converts the unnormalized external P2 into the unit sphere harmonic.
The q1,q3 fields are already canonically normalized in (9). All g in
(16) are dimensionless coefficients in the a=1 action.

The inverse-mass sum can be evaluated without subtracting tiny
eigenvalues. Set tc=(sqrt(3),0)^T, ta=(0,1/sqrt(C+1))^T and Ma=Ma_3.
The algebraic matrix

```
B=tc vc+ta va-(C/6)tc h(dI-Ma)^2-(1/2)ta k(dI+Ma)
  +Mp[(C/6)tc h(dI+Ma)+(1/2)ta k]
 = [[-2sqrt(42)(2C^2+22C+27)/(35C),
                     -4(C+2)sqrt(210(C+1))/(35C)],
    [-12(C+8)sqrt(14(C+1))/(35C),
                     -6sqrt(70)(2C+1)/(35C)]]
```

satisfies g_ps=sqrt(lambda) ep^T B es/[N sqrt(C mu_p)]. Hence

```
Q0=sum_(p,s) g_ps^2/mu_p
  =lambda tr(Mp^-2 B B^T)/(N^2 C)
  =lambda (53C^2+830C+944)/(280pi C)
  =(G/a^2)(53+830/C+944/C^2)/35.                        (17)
```

Here lambda C=kappa^2/a^2=8pi G/a^2 restores the physical units.
The identity is exact for every C>=1. It bounds the operator norm of
the inverse-polar-mass weighted positive coupling Gram matrix by its
trace. In particular the small polar and axial masses did not produce
a new growing factor of C in Q0/(G/a^2).

Rotational invariance determines the other magnetic components of this
specific l=2,l=1,l=3 vertex. The scalar intertwiner occurs once in
1 tensor 3 -> 2, so its m-dependence is the Wigner 3j symbol. For a
fixed unit-normalized external m2,

```
(3j(1,2,3;0,0,0))^2=3/35,
sum_(m1,m3) (3j(1,2,3;m1,m2,m3))^2=1/5,
Qall=(7/3)Q0=(G/a^2)(53+830/C+944/C^2)/15
    <=(1827/15)(G/a^2), C>=1.                         (18)
```

This reconstructs the magnetic degeneracy of this cubic vertex only.
It does not reconstruct every non-axisymmetric quartic from (8).
There is no extra identical-field Wick factor in (18): the l=1 and
l=3 internal fields are different. Every external m has the same
norm; an additional factor of five is not an operator norm on that
direct sum. The checker verifies the 3j sum for all five external m.

For the retained k=2000 family, C=sqrt(1+12pi k^6) and
a^2/G=k^10(C-1)/(6C), (18) is 2.0703125007e-32.
This is a coupling norm, not a loop self-energy, stress error or
probability of opening. A complete mixed-pair spectral bound, the
appropriate local counterterms and the common observable are still
needed before adding it to the prior dipole-pair response.

## Verification and next obligation

[The checker](../checks/axisymmetric_axial_action.py) reconstructs the
Ricci and Maxwell contractions, verifies the Lorentzian dual action,
derives the angular contact matrices by exact polynomial integration,
checks harmonic identities and the general mass trace/determinant,
verifies the curved unequal-mass stress identity and diagonal source
regression, and proves (17) by symbolic matrix algebra. Independent
nonlinear quadrature checks all kinetic tensor components and mass
contacts through second external order. A separate 90-digit direct
eigenbasis calculation tests (16)-(17) at C=1,4,100,10^12 and the
candidate C; it uses rationalized light eigenvalues.

The receipt and source/run record retain exact inputs, terminal outcomes
and development failures. These computational checks support the displayed
ordinary identities; finite parameter samples do not prove a continuum
loop bound. No Lean, cloud spend, publication or outreach is used.

Remaining to-do list: complete the mixed spectral response and its local
subtraction; eliminate the nonlinear even constraints and carry the
resulting contacts and boundary terms; include or bound the remaining
angular/parity and fermion sectors; fix a common regulator, measure and
finite renormalization conditions. Then match the physical state and
drive, complete stress/geometry, stability, finite-energy transport and
preparation/creation before the S080 result audit.
