# Initial-time obstruction and a smooth causal dipole comparison

S072, 2026-10-13. Sharp zero response at a finite initial time fails
as a continuum prescription for the retained dipole-induced polar
fluctuations. A positive subset of pair states already makes the
lowest light-master variance diverge. Smoothing a later measurement
in time does not generically remove this initial-time term.

A smooth diagnostic envelope on the complete master source gives a
different result: the spatially averaged sphere response has a uniform
finite bound during and after the ramp, with both polar branches and
every retained radial mode included. This is a response theorem for a
specified source envelope. Such an envelope is not by itself a conserved
four-dimensional preparation. The constrained dressed initial state
and actual boundary-driven source remain the physical target.

## Scope and prior argument

Retain the fixed global product, H2 dipole vacuum, Friedrichs boundaries,
complete quadratic source and its local term from
[the source derivation](dipole-quadrupole-source.md) and
[the vacuum propagation result](dipole-vacuum-response.md).
Use dimensionless global time t=tau and spatial coordinate sigma.
There are five sphere quadrupoles; all statements below concern this
leading dipole contribution, not all angular or massive field sectors.

The initial-state issue is established prior art.
[Hu, Roura and Verdaguer, Appendix D](https://arxiv.org/html/gr-qc/0402029v1#A4)
explain how a finite initial-time step can make induced metric
correlations singular, how initial matter/metric correlations can matter,
and why a finite-time gravitational interaction switch requires care
with conservation and the Bianchi identity. Appendix D was read in full;
their large-N quantum-gravity equivalence and full interacting-state
construction are not imported into this three-component dipole model.
Here the obstruction and a quantitative smooth-source bound are derived
for the actual previously retained polar equations.

For reference, set x=1/C, r=sqrt(1+12x+12x^2), and define

```
Delta=(1+sqrt(9+16x))/2,     D_j=2Delta+2j,
lambda_j=D_j(D_j-1),
mu_minus=24x^2/(1+6x+r),     mu_plus=1+6x+r,
h_r=(1+sqrt(1+4mu_r))/2,     omega_r,n=n+h_r.
```

The label r on h_r or mu_r denotes a branch; the scalar radical r
in the first line above is used only in the eigenvector formulas below.
One has `0<mu_minus<=2`, `2<=mu_plus<=12`, and omega_r,n>=1.
The endpoint C=infinity can be taken as a limit. The vacuum pair
decomposition is `W_Delta^2=sum_j rho_j W_Dj`, with rho_j>0 and
the gamma-function weights recorded in the preceding note.

## A divergent subset of the sharp initial-data prescription

Symmetrize the two masters as z=(s chi,Acal), with
`s=sqrt(b^2/(3kappa^2))` and `s/(b^2 E0)=c=1/sqrt(3(C+1))`.
Take normalized eigenvectors

```
e_minus=(-sin(theta),cos(theta)), e_plus=(cos(theta),sin(theta)),
sin(theta)^2=6x(1+x)/[r(r+1)],
cos(theta)^2=(r+1)/(2r).
```

The stationary particular light-master transfer, in units G E0, is
`f_minus(lambda)=-sin(theta)c F_chi(lambda)+cos(theta)F_A(lambda)`.
At large lambda it tends to the nonzero local coefficient

```
v_minus=-(9/4)c sin(theta)-(3/4)cos(theta)<-3/(4sqrt(2)). (1)
```

Project onto the lowest light spatial mode u_0^h, where h=h_minus.
Within each pair representation keep only its lowest mode u_0^Dj.
Their overlap has the exact positive beta integral

```
J_j=N_0^h N_0^Dj sqrt(pi)
          Gamma((h+D_j+1)/2)/Gamma((h+D_j+2)/2),
(N_0^D)^2=Gamma(D+1)/[sqrt(pi)Gamma(D+1/2)].             (2)
```

The gamma-ratio asymptotic
`Gamma(y+a)/Gamma(y+b)~y^(a-b)`
([DLMF 5.11.12](https://dlmf.nist.gov/5.11.E12)) gives, at fixed C,

```
rho_j ~1/(pi j),       J_j^2 ~2sqrt(pi)(N_0^h)^2/sqrt(D_j),
D_j~2j,               f_minus(lambda_j)->v_minus.
```

Consequently the positive contribution of this subset to the particular
initial velocity variance, divided by (G E0)^2, has individual terms

```
L_j=(16pi/15)rho_j f_minus(lambda_j)^2 (D_j/2) J_j^2
   ~ L/sqrt(j),
L=(16sqrt(2pi)/15)(N_0^h)^2 v_minus^2 >0.                (3)
```

Its partial sum grows as 2L sqrt(J). This is an ordinary asymptotic
proof of divergence, not an inference from a finite cutoff plot.
The field-position subset instead decays as j^(-5/2).

For a pair frequency Omega and master frequency omega the sharp
zero-data response multiplies the particular amplitude by

```
B_0(t)=exp(-i Omega t)-cos(omega t)
                                  +i(Omega/omega)sin(omega t). (4)
```

At t=pi/(2h), (3) makes the lowest light-mode variance diverge, with
primary terms asymptotic to `(G E0)^2 L/(h^2 sqrt(j))`.
More generally, for a real smooth test f supported at positive times
and `S_f=integral f(t)sin(h t)dt !=0`, its smeared (4) contains
`i(D_j/h)S_f`. Its Fourier term decreases while the other terms remain
bounded. The same divergence survives, multiplied by S_f^2/h^2.
There are exceptional tests with S_f=0; no assertion that every test
diverges is needed to reject the prescription.

All omitted pair contributions are nonnegative in the variance.
Deleting a positive divergent norm by a c-number subtraction would not
construct a finite-variance state. Nor is this an instability theorem
about a correctly correlated initial state. It rejects the use of
sharp independent zero master data with this vacuum source as a
continuum preparation.

## Smooth source envelope, including its derivative terms

For a diagnostic comparison set
`(Box2 I-K)u=g(t)r_source`, where r_source is the complete original
dipole master source. This is an equation for the masters, not yet a
four-dimensional matter/gravity experiment. Let g=0 before t=0 and
g=1 after T, with g'>=0, g'(0)=g'(T)=0 and
`V=integral_0^T |g''(t)|dt <infinity`.

The same source projected onto any branch mode is
`(omega^2-Omega^2) A exp(-i Omega t)`. Superpose delayed steps using
g'. For all t>=0, define
`F_t(nu)=integral_0^min(t,T) g'(u)exp(-i nu u)du`. The exact response is

```
q_g(t)=A B_g(t),
B_g(t)=g(t)exp(-i Omega t)
 -(Omega+omega)/(2omega) exp(-i omega t) F_t(Omega-omega)
 +(Omega-omega)/(2omega) exp(+i omega t) F_t(Omega+omega). (5)
```

This expression also works at Omega=omega: then F_t(0)=g(t) and
the response is zero, as required by the vanishing source overlap.
It solves the original multiplied source with zero initial data.

Equivalently, writing q_g=g q_part+delta q requires
`delta q''+omega^2 delta q=-g''q_part-2g'q_part'`.
The commutator terms are essential. Simply multiplying the stationary
response by g would not solve the stated equation.

Integration by parts, including the upper endpoint when t<T, yields

```
|F_t(nu)| <= min(1,V/|nu|).                              (6)
```

Indeed `g'(u)+integral_0^u |g''|<=V`, because g'(T)=0.
For Omega<=2omega the sum of the absolute prefactors of F_t in (5)
is max(1,Omega/omega)<=2. For Omega>2omega, using (6) bounds that
sum with its Fourier factors by 2V/omega. Since omega>=1,

```
|B_g(t)| <= K_g=1+2max(1,V),       for every t>=0.         (7)
```

This estimate is uniform in both master branches, every master radial
mode and every vacuum pair spectral descendant. It removes the
ultraviolet factor in (4), without deleting the local source term or
assuming separated frequencies.

## An explicit bound with every retained radial mode

The sphere part contributed by branch r has spectral transfer

```
chi_r=(G/b^2) A_r(lambda_j) Q_m,j,
A_r=e_r,1^2 F_chi+sqrt(3(C+1)) e_r,1 e_r,2 F_A.
```

The previous intervals `|F_chi|<=51/20`, `|F_A|<=3/4`, together with
`sqrt(3(C+1))|e_r,1 e_r,2|=3(1+x)/r<=3`, give
`|A_r|<=24/5`. This retains both fields and both branches.

For any bounded real spatial test f at equal time, the Wick covariance
and the Hilbert-Schmidt norm of the one-particle covariance give

```
Var Q_m[f] <= (4pi/15)||f||_infinity^2
                         sum_(p>=0)1/(p+Delta)^2.          (8)
```

To see this, bound |f(sigma)f(sigma')| by its supremum in the integral
of W(sigma,sigma')^2, and use orthonormality to obtain
`integral integral W^2=sum_p 1/[4(p+Delta)^2]`.
Finite spectral cutoffs justify the integrations; the square-summable
one-particle spectrum supplies the continuum limit. For each normalized
master mode u_n, one-dimensional integration and Cauchy-Schwarz give
`||u_n||_infinity^2<=2||u_n||||u_n'||<=2omega_n`.
Thus the RMS of its particular sphere coefficient is at most

```
(24/5)(G/b^2)sqrt(A0 omega_n),
A0=(8pi/15)(pi^2/6-1).                                    (9)
```

Positivity of the pair spectrum and (7) multiply this by at most K_g.
For a common real spatial test f in both operator domains D(H_r),
expand f in each master basis. Cauchy-Schwarz gives
`sum_n |f_n|sqrt(omega_n)<=||H_r f||sqrt(zeta(3,h_r))`.
Since h_r>=1, triangle inequalities in the vacuum Hilbert norm yield

```
RMS chi_m[f,t] <= (24/5)(G/b^2)K_g sqrt(A0 zeta(3))
                            (||H_minus f||+||H_plus f||). (10)
```

This controls the complete sums, including cross-correlations between
branches and modes. It is not a diagonal-in-master-mode assumption.

Choose the normalized spatial average `f=8cos^4(sigma)/(3pi)`, which
lies in both Friedrichs operator domains and has integral one. Then

```
H_r f=8/(3pi)[16cos^4(sigma)+(mu_r-12)cos^2(sigma)],
||H_r f||^2=64/(9pi)[3mu_r^2/8+mu_r+4],
||H_minus f||<=sqrt(160/(3pi)),
||H_plus f||<=sqrt(4480/(9pi)).                            (11)
```

Use the quintic monotone ramp `g=10(t/T)^3-15(t/T)^4+6(t/T)^5`
on [0,T], with constant extensions. Its variation is V=15/(4T).
At T=1/4, K_g=31. Rotational invariance gives the fractional sphere
metric `h_sphere[f]=2sum_m chi_m[f]Y_2m` RMS equal to
sqrt(5/pi) times the single-m RMS. Equations (10)-(11) imply

```
RMS h_sphere[f,t] < 4000 G/b^2,    every t>=0.              (12)
```

For example, pi>157/50, pi<22/7, zeta(3)<5/4 imply
`A0<11/10`, `sqrt(A0 zeta(3))<6/5`, `sqrt(5/pi)<13/10`,
and the sum of the two norms in (11) is less than 17.
The resulting exact rational constant is 3946.176<4000.
At the retained k=2000 geometry, (12) is below 4.772e-40.
This is an ensemble RMS bound at each time, not a probability bound on
the supremum of a trajectory. It concerns the stated spatial average
of the leading dipole response, not renormalized energy or full stress.

## Why this does not yet define a physical preparation

Multiplying a conserved stress and current by a time-dependent g gives

```
nabla_mu(g T^(mu nu))-F^nu_mu (g J^mu)
                                      =(partial_mu g)T^(mu nu),
nabla_mu(g J^mu)=(partial_mu g)J^mu.                       (13)
```

The missing switching apparatus, stress and current are physical
requirements. An envelope on the reduced master source alone cannot
be asserted to preserve the four-dimensional constraints. Therefore
(12) is not compared with the capture margin as a physical error.
The actual selected boundary interaction must supply its own conserved
source and correlated initial dressing. Intrinsic homogeneous master
fluctuations also remain to be included.

The null (3)-(4) gives an ultraviolet matching test. If the initial
position error has finite vacuum Hilbert norm, a relative initial
velocity mismatch epsilon_j on the primary subsector requires
`sum_j |epsilon_j|^2/sqrt(j)`. An asymptotic power mismatch
epsilon_j~j^(-p) needs p>1/4 for this series to converge. A finite-norm
position error cannot cancel its divergence for a later test with
nonzero S_f. This is not a sufficient all-state criterion, but a fixed
nonzero velocity mismatch fails under the stated position assumption.
A stationary particular solution only supplies formal dressing data;
it is not a constructed gauge-constrained interacting quantum state.

## Verification and continuation

[The checker](../checks/dipole_causal_preparation.py) checks the exact
commutator, ramp variation, operator norms and rational bounds; tests
the positive primary asymptotics in three geometries; and compares
(5) with independent ODE solutions during and after switching,
including an exact frequency coincidence. The
[receipt](../evidence/dipole-causal-preparation.json) records source
hashes and scopes. The [review record](../evidence/dipole-causal-review.json)
preserves the source comparison, snapshots and any development failures.
The ordinary positive-series and operator arguments, rather than finite
numerical tests, justify the null and all-mode bound.

Next, derive the initial dressing and source from the actual constrained
fermion/dipole/boundary theory, retaining H3/H4 and induced response.
The sharp independent zero-data route is rejected; the smooth diagnostic
bound identifies sufficient regularity without supplying that physical
state. No general impossibility, creation or stability result follows.

Remaining to-do list: construct constrained initial dressing and bound
interactions; match the physical driven source and boundaries; obtain
complete renormalized stress and evolving geometry; establish pair
preparation, sufficient stability and finite-energy transport.
