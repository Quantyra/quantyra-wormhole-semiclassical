# Physical polar normalization and intrinsic vacuum fluctuations

S072, 2026-10-13. This calculation supplies the quadratic symplectic
normalization missing from the retained classical polar response. It
quantizes both coupled Einstein-Maxwell branches on the fixed global
product and bounds their intrinsic, time-smeared sphere fluctuation over
every l>=2 harmonic and radial mode. It does not construct the interacting
initial state or the physical boundary drive. The existing induced-source
bounds retain their original scopes.

## Action and conventions

Use the four-dimensional action

```
S = integral sqrt(-g) [(R-2Lambda)/(2 kappa^2)-Fbar^2/4],
kappa^2=8pi G, C=b^2/a^2>=1, kappa^2 E0^2 b^2=C+1.       (1)
```

The original magnetic background is Fbar_ij=-E0 b^2 epsilon_ij.
The Maxwell field has already been divided by its gauge coupling. The
cosmological constant is fixed. Set hbar=1. Real sphere harmonics obey
integral Y_lm Y_l'm'=delta_ll' delta_mm', L=l(l+1), l>=2. Write
V_i=epsilon_i^j D_jY, so dV=L Y epsilon. In this magnetic frame
delta Abar_i=Acal V_i is the same master normalization used in the
[dual electric response](quadrupole-product-response.md). In particular,
the perturbation of the dual electric field in an orthonormal frame is
-L Acal/b^2-2E0 chi. Changing orientation changes both mixed
signs, not any eigenvalue or variance.

The canonical reduction uses the lapse, shift and electromagnetic
constraint structure of [Arnowitt, Deser and Misner, sections 3.2-3.3
and 6.1](https://arxiv.org/html/gr-qc/0405109v1). Our sign convention for
extrinsic curvature is the opposite of their displayed K definition;
the quadratic combination is unchanged. Their point-particle and
self-energy conclusions are not imported.

## Constraint reduction before quantization

Choose a proper radial coordinate y:

```
ds2^2=-N(y)^2 dT^2+dy^2, N=cosh(y/a), N''/N=1/a^2,
delta gamma_yy=D Y, delta gamma_ij=2b^2 chi Y gammaSphere_ij,
shift_y=beta Y, shift_i=alpha D_iY, u=(chi,Acal)^T.       (2)
```

The scalar spatial gauge removes the mixed radial-angular metric and
traceless sphere deformation. Retain both shifts through variation.
Dots mean T derivatives, primes mean y derivatives. The first-order
extrinsic curvature and electric field are

```
K_yy=(dot D-2beta')Y/(2N),
K_yi=-(beta+alpha')D_iY/(2N),
K_ij=[b^2 dot chi Y gammaSphere_ij-alpha D_iD_jY]/N,
Fbar_Ti-shift^j Fbar_ji=(dot Acal-E0 alpha)V_i.            (3)
```

Use integral |DY|^2=L and integral |DDY|^2=L(L-1). Sphere integration
of N sqrt(gamma)(KijK^ij-K^2)/(2kappa^2) plus Maxwell electric energy
gives the quadratic kinetic density

```
Lkin=b^2/(2kappa^2 N) [
 -2(dot chi)^2-2L alpha dot chi/b^2-L alpha^2/b^4
 -(dot D-2beta')(2dot chi+L alpha/b^2)
 +L(beta+alpha')^2/(2b^2)]
 +L(dot Acal-E0 alpha)^2/(2N).                           (4)
```

The lapse imposes the linear Hamiltonian constraint. Its kinetic
contribution begins at second order because the background has zero
extrinsic curvature and electric field. Direct spatial curvature and
magnetic-energy variations give

```
delta R3/Y=-4chi''+LD/b^2+(2L-4)chi/b^2,
delta Bmag^2/Y=-2E0 L Acal/b^2-4E0^2 chi,
delta R3=kappa^2 delta Bmag^2,
D=4b^2 chi''/L-(2+4C/L)chi-2kappa^2 E0 Acal.             (5)
```

After using the remaining time gauge to set alpha=0, the radial and
angular momentum constraints are, respectively,

```
beta=4b^2/L [dot chi'-(N'/N)dot chi],
dot D=-2dot chi+beta'+(N'/N)beta-2kappa^2 E0 dot Acal.     (6)
```

The second equation is precisely the time derivative of (5) after
substituting the first and N''/N=C/b^2. Thus this gauge choice has not
discarded an independent constraint. Maxwell's linear Gauss equation
holds because V_i is coexact and has zero angular divergence. These
l>=2 perturbations add no monopole electric or magnetic charge; no
independent charge rotor is introduced. Substitute (5)-(6) into (4) and
integrate by parts. The result is

```
Lkin,red=(dot u)^T W_L dot u/(2N),
W_L=[[2b^2/kappa^2 (1+2C/L), 2b^2 E0],
     [2b^2 E0,                L]],
det W_L=2b^2(L-2)/kappa^2>0.                            (7)
```

Writing v=dot chi, the discarded density is the exact derivative
`partial_y{4b^4/(kappa^2 L)[v v'/N-3N'v^2/(2N^2)]}`.
The checker retains this term in the pointwise identity.
For each finite-mode smooth datum its boundary value vanishes with
the Friedrichs falloff. Complete the quadratic form from that dense
domain. No unfixed boundary degree of freedom is being quantized here.

The previously independently verified source-free equations are
(Box2 I-K)u=0, with

```
K=[[(L+2C)/b^2, L kappa^2 E0/b^2],
   [2E0,        L/b^2]].                               (8)
```

W_L K is symmetric. The kinetic term just obtained and these complete
equations fix the reduced bulk quadratic action, up to fixed-source
boundary terms, as

```
S2=-1/2 integral sqrt(-g2)
   [(D_a u)^T W_L D^a u+u^T W_L K u].                  (9)
```

This identifies the potential using the verified linear field equations;
it is not a separate brute-force expansion of the full curvature
potential. The absolute kinetic coefficient comes from (1)-(7), not
from arbitrarily multiplying the classical master equations.

## Physical oscillators and correlated vacuum

Let S=diag(sqrt(2b^2/(Lkappa^2)),1), z=S u. The symmetric mass matrix
M_L=a^2 S K S^-1 and its spectral projectors P_+/- are retained from
the response note. The new physical kinetic matrix is

```
W_z=S^-T W_L S^-1=C M_L
   =[[L+2C, sqrt(2L(C+1))],[sqrt(2L(C+1)), L]],
R=sqrt(C^2+2L(C+1)),
mu_- = L(L-2)/[C(C+L+R)], mu_+=(C+L+R)/C,
Psi_r=sqrt(C mu_r) e_r^T z, r=+,-.                     (10)
```

The Psi_r have the usual positive scalar kinetic action and mass
mu_r/a^2. In global conformal coordinates their conjugate momenta are
partial_tau Psi_r. Equivalently, [z_i,partial_tau z_j]=i(W_z^-1)_ij
delta(sigma-sigma'). At C=1 the masses reduce to l(l-1) and
(l+1)(l+2), agreeing with the bosonic scalar sector of
[Michelson and Spradlin, section 3.8, equations 3.44-3.45](https://arxiv.org/html/hep-th/9906056v3).
Their supersymmetry or boundary-mode conclusions are not extended to
C>1. The l=1 zero determinant is only a gauge-sector consistency
check; formulas requiring W^-1 are restricted to l>=2.

Choose the free global Friedrichs vacuum of the Psi_r. With
h_r=(1+sqrt(1+4mu_r))/2 and canonical scalar Wightman function W_h,

```
<z z'^T> = sum_r [P_r/(C mu_r)] W_h_r(X,X'),
<chi chi'> = Lkappa^2/(2b^2)
             sum_r [(P_r)_11/(C mu_r)] W_h_r(X,X'),
sum_r (P_r)_11/(C mu_r)=(W_z^-1)_11=1/(L-2).            (11)
```

This also fixes the electromagnetic/metric cross covariance. For fixed
L and large C, (P_-)_11~L/(2C) while Cmu_-~L(L-2)/(2C).
Their ratio tends to 1/(L-2). The light branch's small response
projector therefore does not suppress its intrinsic quantum metric
variance by itself. The old classical response calculation remains valid.
As a source-normalization consistency check, its quadrupole derivative
vertex v=G(9/(4b^2),-3E0/4)^T obeys W_6 v=(3G/kappa^2,0)^T.
This does not reconstruct the complete off-shell cubic interaction or
the matter/metric dressed vacuum.

## Bound over the complete angular and radial spectrum

At the global center sigma=0, smear in dimensionless global time tau
with a normalized real Gaussian of standard deviation s=1/4. Its
squared Fourier transform is exp(-s^2 omega^2). The even canonical
radial modes have central vacuum weights

```
a_n(h)=Gamma(n+1/2)Gamma(n+h)
       /[2pi Gamma(n+1)Gamma(n+h+1/2)] <= 1/pi,
V(h)=sum_n>=0 a_n(h) exp[-s^2(2n+h)^2]
    <= B exp(-s^2 h^2), B=(1+sqrt(pi))/pi<9/10.           (12)
```

The gamma ratios decrease with their positive argument because the
digamma function increases. Their maximum is n=0,h=1. Also
(2n+h)^2>=4n^2+h^2 and the decreasing Gaussian sum is bounded by
1+integral_0^infinity exp(-x^2/4)dx. The exact normalized radial
eigenfunctions give (12); independent numerical quadrature checks
their normalization and central weights.

Define the scalar-sector fractional sphere metric
H=2 sum_l>=2,m chi_lm Y_lm at a fixed angle. It is the sphere
amplitude in the fully fixed scalar gauge, equivalently its linear
gauge-invariant harmonic completion. This is not yet a relational
measurement by a physical detector. Rotational addition and (11) give

```
Var H <= 3kappa^2 B/(4pi b^2)
         sum_l>=2 (2l+1) exp(-s^2 h_-^2).               (13)
```

Here L/(L-2)<=3/2, and h_+>=h_-. For t=L/C,

```
mu_- = 1+t-sqrt(1+2t(1+1/C)) >= 1+t-sqrt(1+4t),
(t/2+17/8)^2-(1+4t)=(4t-15)^2/64>=0,
h_-^2=mu_-+h_- >= t/2-1/8.                            (14)
```

Put alpha=1/(32C)<=1/32. The positive unimodal function
f(x)=(2x+1)exp[-alpha x(x+1)] has integral 1/alpha and maximum
sqrt(2/alpha)exp(alpha/4-1/2)<=sqrt(2/alpha). A nonnegative
unimodal function's integer sum is at most its integral plus its
maximum (compare each monotone side with the adjacent unit intervals).
Consequently

```
sum_l>=0 (2l+1)exp[-alpha l(l+1)]
 <=1/alpha+sqrt(2/alpha)<=5/(4alpha)=40C,
Var H < (27648/127) G/a^2 <225 G/a^2,
RMS H <15 sqrt(G)/a.                                   (15)
```

We used exp(1/128)<128/127. This is an infinite-sum inequality, not
a cutoff extrapolation. Both physical polar branches and all l>=2
angular multiplicities are included. The n=odd central modes vanish.
For the retained k=2000 candidate, a^2/G=k^10(C-1)/(6C), so the
bound is below 1.149e-15. The receipt supplies the unrounded value.

## Scope and continuation

The calculation is a free quadratic Einstein-Maxwell result on the
fixed global product with specified reflecting boundaries and vacuum.
It resolves one missing intrinsic fluctuation term and its physical
normalization. It is not a thermal, horizon, interacting or driven
state calculation. The observable differs from the spatial average
in the causal-envelope note. No independence between intrinsic and
induced fluctuations in a dressed state is assumed, so those bounds
are not simply added. No full metric, renormalized stress, ultraviolet
effective-theory error, trajectory-supremum probability or capture
estimate follows from (15). Monopole, dipole and physical boundary
degrees of freedom still require their own constraints and states.

The ordinary harmonic identities prove the action reduction for all
l>=2. The checker independently derives the spatial connection,
linear curvature and kinetic contractions for l=2,3,4; verifies both
shift equations and the exact boundary derivative; checks the
canonical commutator matrix and the known C=1 masses; and compares
finite spectral sums with the analytic all-mode bound. Finite sums
are diagnostics only. This is scoped primary-source comparison, not
an exhaustive novelty review or a new general quantization theorem.

The first development run passed all three direct tensor calculations
and both shift equations, then rejected an incorrect integration-by-parts
boundary term in the checker. The second run used the exact derivative
printed above and passed without changing the reduced kinetic matrix,
the fluctuation bound or any numerical tolerance. Both source snapshots
and terminal outcomes are retained in
`C:/Users/dfred/QuantyraTools/SpaceProofs/polar-canonical-vacuum-2026-10-13`.

The next physical target is the interacting constrained initial state
and conserved source from the actual fermion/dipole/boundary action.
The new normalization makes its quantum amplitudes and source vertices
well defined; it does not by itself supply those correlations.

Remaining to-do list: interacting initial dressing and induced response;
actual driven boundary matching; complete stress and evolving geometry;
pair preparation, sufficient stability and finite-energy transport;
distinguish opening a prepared bridge from creating one.
