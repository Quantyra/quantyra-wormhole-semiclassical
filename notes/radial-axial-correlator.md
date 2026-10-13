# Radial electric and one-axis angular boundary correlation

S072, 2026-10-13. **The radial electric field and one commuting angular
generator now have a joint, normalized boundary thermal correlator in the
projected quadratic Einstein-Maxwell/massless-fermion theory.** The
angular calculation retains the dynamical electromagnetic dipole and
the full radial profile. A bound uniform in optical length and channel
count puts the candidate's initial midpoint loss below 8.334e-5 against
the same-cut free reference. Its joint fermion charge lattice is retained.

This is a one-axis result. The two other rotation generators do not
commute with it, and cannot be supplied by multiplying three copies of
this answer. Non-Abelian interactions, other angular harmonics, higher
Landau levels, nonlinear gravity, physical horizon gluing and driven
stress remain open. The global Casimir model in
[angular-reservoir-correlator.md](angular-reservoir-correlator.md) remains
a separate comparator, not the Hamiltonian derived here.

## Constrained radial model and normalization

Use q>=2 massless complex Dirac LLL channels of spin j=(q-1)/2 on
the fixed spherical RN-AdS background. Retain the spherical electric
field and the metric/Maxwell dipole along one fixed rotation axis. Keep
the Einstein-Maxwell action to quadratic order in those perturbations.
The axis couples through the diagonal generator m=-j,...,j. Massless
bosonization for this commuting current is Gaussian; this is not a
Gaussian replacement of the entire SU(2) current algebra.
The Abelian reduction and the distinct non-Abelian sector are reviewed
in [Maldacena, Milekhin and Popov, appendix A](https://arxiv.org/html/1807.04726v3#A1).

On the optical interval [-L,0], fix zero inner electric and angular flux
and the inner dipole canonical field to zero. Use reflecting fermions
with the existing antiperiodic unfolding. The outer massive dipole has
the natural zero-source boundary condition. Sum the original integer
species charges N_m; there is no independent horizon charge or rotor.
These are comparison boundary/ensemble choices, not a physical horizon
state. The spherical boson and one angular boson are orthogonal because
sum m=0. Put

```
Q=sum_m N_m,       S=sum_m m N_m,
S_J=sum_m m^2=q(q^2-1)/12,
lambda_e=q/pi,     lambda_J=S_J/pi.                     (1)
```

The [constrained axial action](causal-compact-response.md) supplies

```
P=q/2,       mu=8pi/(3g^2),
kappa_g=r^4/(6G),  kappa_e=4pi r^2/(3g^2),
K=diag(kappa_g,kappa_e),
c=(-P,1)^T,  d=(1,2)^T/(q+1),
Lmat=(sqrt(mu)c, sqrt(lambda_J)d),
V_J=f Lmat^T K^{-1} Lmat = [[A,B],[B,D]] >=0.            (2)
```

Here A,B,D include f. In particular

```
A=f[2/r^2+4 Qmag^2/r^4],  Qmag^2=pi G q^2/g^2,
D=f lambda_J/(q+1)^2 [6G/r^4+3g^2/(pi r^2)].            (3)
```

The first component is the dynamical electromagnetic dipole; there is
no additional independent gravitational dipole wave. The exceptional-mode
equation (74) in [Kodama and Ishibashi, section 4.2.2](https://arxiv.org/html/hep-th/0308128v3)
checks the pure Einstein-Maxwell potential A after duality. It does not
supply the interacting fermion correlator derived below.

Choose the flux orientation so that the canonical field
Z=-Lmat^{-1}p obeys Z_2(0)=S/sqrt(lambda_J). The conserved rotation
flux is -n^T p, n=(1,P)^T: n^T c=0 and n^T d=1. Thus the second
canonical field is precisely the bosonized enclosed spin density,
rho_J=sqrt(lambda_J) partial_x Z_2. The quadratic Hamiltonian is

```
H_J=(1/2) integral [|Pi_Z|^2+|Z_x|^2+Z^T V_J Z] dx,
Z(-L)=0,  Z_2(0)=S/sqrt(lambda_J),  Z_1,x(0)=0.         (4)
```

This follows by the constrained dualization already derived in the
causal action; retaining the compact boson sectors restores (1). The
spherical electric Hamiltonian is the previously derived scalar one,
with V_e=f lambda_e g^2/(4pi r^2). The remaining q-2 bosons are free
within this projection. No second determinant for either retained
fermion current is appended to (4).

## Mixed boundary problem and the physical ladder

Let z_s(x) be the vector minimizer of

```
K_J,L(s)=min_(z(-L)=0,z_2(0)=1)
 integral [|z_x|^2+z^T V_J z+s^2 |z|^2] dx.             (5)
```

The first outer component varies freely, so z_1,x(0)=0. Let K_e,L(s)
be the existing electric scalar Dirichlet-to-Neumann function. For
homogeneous angular oscillators, impose Dirichlet at the inner cut,
Dirichlet on the outer light component and Neumann on the outer heavy
component. If their modes and frequencies are u_a,omega_a, set
c_a=<u_a,z_0>. Completing the square in (4) gives

```
Z=S z_0/sqrt(lambda_J)+Y,
K_J,L(s)=K_J,L(0)+sum_a s^2 omega_a^2 c_a^2/(s^2+omega_a^2). (6)
```

The proof is the same spectral minimization as for the scalar electric
operator, now with the mixed boundary domain. Integration by parts
removes the cross term because z_0 solves its Euler equation and its
natural boundary condition. Parseval gives (6). A finite-element mesh
has an additional b_mesh s^2 term; the same mass matrix gives the same
b_mesh in both references, so it cancels in the subtraction below.

A creation at the outer boundary in channel m takes
(Q,S)->(Q+1,S+m) while leaving interior enclosed flux and dipole
coordinates fixed. In the centered oscillator variables it translates

```
Y_a -> Y_a-m c_a/sqrt(lambda_J),                        (7)
```

as well as the electric translations already derived. In particular it
translates the *coupled* angular normal modes. Leaving their equilibrium
profile fixed would give the wrong boundary insertion. The equal-time
Wilson dressing is fixed by the same outer frame before radial gauge
is chosen. This determines the operator in the projected model, not
its finite normalization in the full four-dimensional theory.

## One charge lattice and one oscillator subtraction

Set c0=pi beta/(2L), e=(1,...,1)^T, and let m be the column of spin
weights. The exact dimensionless charge energy is N^T A_N N, where

```
x_e=beta K_e,L(0)/(2lambda_e),
x_J=beta K_J,L(0)/(2lambda_J),
y_e=x_e-c0/q >=0,  y_J=x_J-c0/S_J >=0,
A_N=c0 I_q+y_e e e^T+y_J m m^T.                        (8)
```

Both inequalities follow from the variational lower bound K(0)>=1/L.
Equation (8) retains the allowed joint charge-spin lattice, including
half-integer S when q is even. Free bosonization gives the kinetic
energy c0 sum N_m^2. Its projections onto e and m are replaced by the
two screened energies in (8), not multiplied by separate rotor sums.
The constructive charge-sector convention is the same as in
[von Delft and Schoeller, sections 4E--4F, 6--7 and appendix B](https://arxiv.org/html/cond-mat/9805275v3).

The angular free reference contains the same independent dipole field
with potential A, and a free angular density boson:
V_J,ref=diag(A,0). Its lift has zero heavy component and
K_ref(s)=s coth(sL). Thus it carries no spurious free dipole vertex.
Use the same ultraviolet normalization for the actual and reference
boundary fermion. For s_n=2pi n/beta, define

```
D_e(t)=2/(beta lambda_e) sum_(n>=1)
 (1-cos(2pi n t)) [K_e,L(s_n)-s_n coth(s_n L)]/s_n^2,
D_J(t)=2/(beta lambda_J) sum_(n>=1)
 (1-cos(2pi n t)) [K_J,L(s_n)-s_n coth(s_n L)]/s_n^2,
Theta_A(v)=sum_(N in Z^q) exp[-(N+v)^T A(N+v)].         (9)
```

The normalized exact thermal correlator for 0<t<1 is

```
G_m(beta t)/G_m,free,same-cut(beta t)
 = exp[-D_e(t)-m^2 D_J(t)]
   [Theta_A_N(t e_m)/Theta_A_N(0)]
   /[Theta_(c0 I_q)(t e_m)/Theta_(c0 I_q)(0)].           (10)
```

Here e_m is the unit vector for the tagged species, not the all-ones
vector e. Multiplying the shifted oscillator thermal traces from (7)
and their electric counterparts gives (10): the charge-energy difference
t(1-t)(y_e+m^2 y_J) combines with the oscillator differences using (6).
The neutral oscillators cancel. Normalized traces cancel the common
charge-independent oscillator partition factors; this is not a claim
that their contribution to the gravitational stress vanishes.

## A matrix bound that does not assume scalar monotonicity

The scalar electric result used a maximum principle to prove
delta K(s)<=delta K(0). That extension is false for general positive
matrix potentials with the boundary conditions (5). For example, on
length L=40, the positive definite constant matrix

```
V=[[1,1],[1,1.0001]]
```

gives delta K(0)=0.0254138903 and delta K(0.4)=0.2289457584.
This follows independently from the exact matrix DtN function
sqrt(V+s^2 I)coth[L sqrt(V+s^2 I)] followed by its heavy-component
Schur complement. The example rules out that proposed general estimate;
it is not a physical instability or a counterexample in the RN family.

Instead, positivity of V_J and a trial with zero heavy component give

```
0 <= K_J,L(s)-s coth(sL)
  <= integral D(x) [sinh(s(x+L))/sinh(sL)]^2 dx
  <= W_J := integral D(x) dx.                          (11)
```

Use the linear trial at s=0. The lower bound follows by dropping the
entire positive potential and minimizing the remaining free derivative
and mass norm with the same boundary constraints. For the upper bound
the mixed term vanishes on the trial; its light diagonal term remains.
This proof controls all mixing strengths without a small-B expansion.
If the inner cut has radius r_c, dx=dr/f gives the exact radial integral

```
W_J/lambda_J = [2G/r_c^3+3g^2/(pi r_c)]/(q+1)^2,
W_e/lambda_e = g^2/(4pi r_c).                           (12)
```

Define b_e=beta W_e/(2lambda_e), b_J=beta W_J/(2lambda_J).
The convergent Matsubara sums then obey

```
0<=D_e(t)+m^2 D_J(t)<=(b_e+m^2 b_J)t(1-t),
0<=y_e+m^2 y_J<=b_e+m^2 b_J.                           (13)
```

These bounds also justify the common ultraviolet point normalization
in (10). The scalar static cap is unnecessary for this result.

## Finite-cut midpoint bound on the actual species lattice

The paired-reservoir proof in [electric-reservoir-midpoint.md](electric-reservoir-midpoint.md)
extends to any finite number of commuting charge components. For a
global comparison energy H0+Qvec^T Y Qvec with Y>=0, each equal-energy
particle/antiparticle pair has charges v,-v. Its free charge character
is 1-4p(1-p)sin^2(v dot alpha/2)>=0. Products remain nonnegative and
reflection symmetric. Gaussian Fourier transformation and reflection
therefore give, for f(s)=<exp[-(Qvec+s)^T Y(Qvec+s)]>_rest,

```
exp(-s^T Y s) f(0) <= f(s) <= f(0).                    (14)
```

Degenerate Y follows by continuity. Removing the tagged pair and summing
its four occupations gives exactly the old paired-reservoir formula
with lambda_m=v^T Y v. Hence its midpoint ratio satisfies

```
exp(-lambda_m/2) <= r_global
 <= exp(lambda_m/4)/cosh(lambda_m/2).                   (15)
```

Apply this to v=(1,m), Y=diag(y_e,y_J). At finite optical length,
constructive bosonization identifies the global comparison's charge
trace with (8); its unchanged oscillator factors cancel. Thus the
radial answer (10) equals

```
r_radial = r_global exp[lambda_m/4-D_e(1/2)-m^2 D_J(1/2)]. (16)
```

With b_m=b_e+m^2 b_J, equations (13)--(16) prove the uniform result

```
exp(-b_m/2) <= r_radial <= exp(b_m/2)/cosh(b_m/2).       (17)
```

For the lower bound use lambda_m<=b_m and D_total<=b_m/4.
For the upper bound drop D_total and use that
exp(lambda/2)/cosh(lambda/2) increases with lambda. These are finite-cut,
finite-q bounds, independent of an optical-volume limit. An incoherent
sum of the diagonal boundary two-point functions inherits the worst
channel bound. No statement about noncommuting spin correlators follows.

For the candidate q=k^4, g=k^-4, L_AdS/l_p=k^5, k>=2000, the extremal
radius b is a lower bound on r_c. Put Y_b=Qmag^2/b^2=1+3X_b and use
the already fixed preparation beta. Since |m|<=j,

```
b_m <= [beta g^2/(8pi b)]
       [1+j^2/(q+1)^2 (12+8Y_b/q^2)]
    < [1/(12k)] (4+2Y_b/q^2),
2Y_b/q^2 < 1e-8,
1-r_radial <= 1-exp(-b_m/2)
           < (4+1e-8)/(24k) <= 8.333334e-5.             (18)
```

The first coefficient bound uses 1/4<=z0<=7/25. The small term follows
from X_b(1+3X_b)=pi k^6, so Y_b<1+sqrt(3pi)k^3. This controls the
initial electric plus one-axis loss against its *same-cut* free state;
it is not the total correction to capture energy or area. Finite-energy
transport, driven subtraction, other angular axes and preparation are
not certified by a small thermal midpoint loss.

## Verification and next action

The checker independently compares coupled finite-element spectra with
continuum matrix Riccati integration and with primal boundary elimination
on two RN profiles. It compares canonical oscillator traces with a
Matsubara representation retaining the mesh contact and its tail, checks
the exact vector-charge Fock transition formula, sums the joint species
lattice, and evaluates the candidate bound at high precision. The
constant-matrix example records why the scalar static cap was replaced.
These are verification details for the ordinary correlator and bound,
not a numerical four-dimensional opening simulation.

All six vector-charge trace comparisons agree within 3.34e-16. Six
radial meshes have primal/spectral agreement within 9.46e-13 and thermal
spectral/Matsubara agreement within 1.12e-15; the latter also retain
explicit, larger tail bounds. Independent continuum errors decrease by
about a factor of four under mesh doubling. Five joint lattice midpoints
and six 70-digit candidate substitutions satisfy the displayed bounds.
These are supplementary numerical checks, not interval certification.
The checker passed on its first run. The matrix static-cap counterexample
is an intentional mathematical test, not a solver failure.

The next step is the noncommuting angular current theory with its actual
radial gauge/metric constraints, or a controlled bound comparing it
with this one-axis calculation. The existing current/coset split must
retain its non-Abelian Ward terms; the all-axis global Casimir result
does not supply the radial comparison. Physical horizon gluing and the
dressed four-dimensional boundary operator are still required before
constructing the driven covariance and renormalized stress.

Remaining to-do list: noncommuting angular and massive interactions;
physical state/operator matching; driven stress and evolving geometry;
preparation, quantum fluctuations, sufficient stability and finite-energy
transport.
