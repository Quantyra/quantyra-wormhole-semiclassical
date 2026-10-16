# A recentered stationary Polyakov-LLL exterior

2026-10-15, E003/S139.3, supporting S074/S075. Ordinary argument,
independently reviewed before focused verification. The explicitly
defined isolated free LLL subtheory below has a horizon-regular,
asymptotically AdS stationary
solution at the prescribed laboratory temperature. Its nonlinear lapse
obeys |alpha|<1e-20 at k=2000, below the component target 1e-8.

This is not the stationary solution of the complete four-dimensional
Dirac/EFT action. The exact four-dimensional Weyl action and the remaining
angular functional are not dropped from the physical theory or bounded
by this component theorem. The product center constructed here is not
identified with S137's full one-loop corrected product.

## Declared sector, state and subtraction allocation

Use one massless four-dimensional Dirac field with q magnetic zero modes,
so the LLL has two-dimensional central charge c=q=N/2. The canonical
angular/radial reduction and source conventions are those of
[the Dirac reduction](dirac-spherical-source-reduction.md). For a spherical
metric write

```
ds4^2=g_ab dx^a dx^b+r^2 dOmega^2
     =(r/b0)^2 [dshat2^2+b0^2 dOmega^2],
ghat_ab=(b0/r)^2 g_ab.
```

The inherited common four-dimensional prescription has the exact split

```
Gamma_D,ren=W4[g,ghat,A]+Gamma_hat[ghat2,A;state,domain].
```

W4 is the already retained
[spherical Weyl functional](spherical-anomaly-action.md), including its
boundary transgression obligations. Define Gamma_hat,LLL to be the
canonical massless two-dimensional Dirac/Polyakov functional on ghat2,
with flat-vacuum normalization, the actual zero-source reflecting
projector, and the Hartle--Hawking state. Define the remaining functional
by subtraction: Gamma_hat,rest=Gamma_hat-Gamma_hat,LLL. This is an
allocation inside the common scheme, rather than independent physical
renormalization of each angular determinant. A different finite local
allocation must be compensated in Gamma_hat,rest.

The isolated subtheory studied here is classical spherical magnetic
Einstein--Maxwell--AdS gravity plus this declared Gamma_hat,LLL. Restoring
W4 and Gamma_hat,rest is a subsequent coupled calculation. No local
counterterm is silently chosen to force an asymptotic cancellation.

The conformal LLL state is fixed by stationarity, reflection at the
outer boundary, zero electric current, and regularity at both horizons.
In null coordinates u=U-x, v=U+x its state constants are
`t_u=t_v=q kappa^2/(48pi)`. The Euclidean period is beta=2pi/kappa
in the same laboratory Killing coordinate U. The two exterior copies
have the corresponding canonical thermofield purification with the
specified bundle pairing. A constant boundary reflection phase does
not change the bulk diagonal stress. There is one horizon half-line,
not a finite optical interval with an additional Casimir length.
Stationarity here concerns each exterior's HH restriction and local
stress. The thermofield cross covariance acquires the usual phases
under forward evolution of both boundary clocks; at zero coupling
these do not change either exterior's local stress. Its relative
phases remain part of the correlated pulse preparation.

This is the established two-dimensional anomaly/conservation and
horizon-regularity construction of
[Christensen--Fulling](https://doi.org/10.1103/PhysRevD.15.2088).
The spherical zero-mode input and physical paired domain are imported
from the original RN-AdS fermion model,
[Bintanja--Freivogel--Rolph, sections 2.1.1--2.1.2](https://arxiv.org/html/2308.00871v2).
Neither source is asserted to prove the component estimates here.

The canonical optical spinor and its Killing-energy spectrum are
unchanged by the four-dimensional Weyl change: both reductions give
`phi=F^(1/4) b0^(3/2) r^(-1/2) Psi_hat` when alpha=0, with the analogous
common lapse factor generally. Thus this change of frame does not replace
the state by a different temperature or an unrelated reservoir. It does
change the allocation of the covariant local stress and radius source.

## Complete sector stress and its physical radius source

Use

```
ds2^2=-a(r)^2 F(r)dU^2+dr^2/F(r), a=exp(alpha),
dx/dr=1/(a F),
dshat2^2=exp(2psi)(-dU^2+dx^2),
exp(2psi)=a^2 F b0^2/r^2,
P=dr/dx=a F, v=dpsi/dx,
v=a[F'/2+F alpha'-F/r].
```

Here primes on F or alpha mean d/dr; a subscript x means optical
derivative. The exact bulk covariant Polyakov stress is

```
T_uu=T_vv=q[kappa^2+v_x-v^2]/(48pi),
T_uv=q v_x/(48pi),
rho2=-T_U^U=q[kappa^2+2v_x-v^2]/(24pi a^2 F),
p2=T_r^r=q[kappa^2-v^2]/(24pi a^2 F).                (1)
```

The orbit sources refer to physical g, with the sphere already
integrated. At fixed physical g, varying r also varies ghat. Hence

```
T_ab^(2)=T_hat,ab,
J_X=-2G T^(2)a_a/r^2, X=r^2,
T^(2)a_a+X J_X/(2G)=0,
nabla^a T_ab^(2)=J_X partial_b X/(4G), j^a=0.        (2)
```

This is the Weyl-invariant part of the source dictionary in
[the spherical source note](spherical-anomaly-sources.md), not an assertion
that the complete four-dimensional trace anomaly vanishes. Its full W4
contribution remains separate. Setting J_X=0 for this allocation is wrong.

On a prescribed classical RN metric a=1 let
`A=4kappa^2+2FF''-F'^2`. The formulas reduce to

```
T_uu=T_vv=q A/(192pi),
T_uv=q F[F''-2F'/r+2F/r^2]/(96pi),
rho2=q[A/F+2F''-4F'/r+4F/r^2]/(96pi),
p2=q[4kappa^2-F'^2]/(96pi F)+q(F'/r-F/r^2)/(24pi),
T^(2)a_a=q[-F''+2F'/r-2F/r^2]/(24pi).               (3)
```

Indeed the two-dimensional conformal change omega=-log(r/b0) has
`omega_uu-2sigma_u omega_u-omega_u^2=0`, where sigma=log(F)/2.
It leaves the null geometric stresses unchanged and shifts T_uv by
`-qF(F'/r-F/r^2)/(48pi)`. At infinity the constant term
`rho2_can=q/(24pi L4^2)` of the canonical physical-two-metric Polyakov
piece cancels in (3). In this product-frame allocation

```
rho2=q[2+kappa^2 L4^2]/(24pi r^2)+O(r^-3).
```

Thus the sector mass tail is integrable. The apparent linear mass tail
of an isolated differently allocated two-dimensional piece is not a
physical ADM divergence. It signals that its four-dimensional Weyl/radius
allocation has not been included consistently.

## Exact nonlinear radial equations

Let Q=r_m^2, L=L4, and introduce

```
gamma=Gq/(12pi), D(r)=r^2-3gamma,
C(r)=1+6r^2/L^2,
V(r)=2+6r^2/L^2-2Q/r^2,
U(r)=r^2 V(r)/b0^2,
H(r)=r^2+3r^4/L^2-Q-gamma C(r).
```

After integrating r^2 Rhat once, the static bulk action is

```
(1/(4G)) integral [4r r_x psi_x+6r_x^2
                    +exp(2psi)U(r)+2gamma psi_x^2]dx.
```

Its equations and the Hartle--Hawking constraint are

```
r psi_xx+3r_xx=exp(2psi)U'(r)/4,
gamma psi_xx+r r_xx+r_x^2=exp(2psi)U(r)/2,
4r r_x psi_x+6r_x^2+2gamma(psi_x^2-kappa^2)
                                      =exp(2psi)U(r). (4)
```

The kinetic determinant is D>0 in the regime below. These equations
include the metric dependence of the exact Polyakov stress. They are
not the replacement of that stress by a fixed RN value. The constraint
is conserved by the first two equations; its constant is fixed by
regularity, v_h=kappa and P_h=0.

Define

```
k(r)=H(r)/(rD(r)),
B(r,F)=[-2-3r^2/L^2+3Q/r^2+3F]/D(r),
l(r)=alpha'(r), A(r,F)=k(r)-rF/D(r),
S(r,F)=3gamma[-F(r^2+6gamma)-3Q+3gamma+r^2
                       +6gamma r^2/L^2+r^4/L^2]/[r^2 D(r)^2].
```

Eliminating the optical auxiliary variables gives the exact radial
system on F>0,

```
F'=A-F l,
l'+[r/D+2A/F]l=S/F.                                 (5)
```

The apparent l^2 term cancels. In particular the second equation is
linear in l for a given F and has the stable horizon damping. Equivalent
useful relations are `P'=a A`, `v'=a B` and
`l=2v/P-P'/P+2/r`. Formula (5) can be checked directly from (4);
it is not obtained by treating a boundary form as an ordinary operator.

## The center, regular horizon and fixed temperature

Let R denote the unique positive root H(R)=0 and C_R=1+6R^2/L^2.
The root is unique since H'(r)=2r[C(r)-6gamma/L^2]>0 for r>=R
in the stated parameter range. It obeys

```
R^2+3R^4/L^2=Q+gamma(1+6R^2/L^2),
R^2=b^2+gamma+O(gamma^2/(C_b L^2)),
R-b=Gq/(24pi b)+O(gamma^2/b^3),                     (6)
```

where b is the classical extremal radius. This center is exact in the
isolated subtheory. Equivalently its tadpole can be reorganized by adding
and subtracting the same local term `-A0/X`, with A0=q C_R/(24pi),
inside the sector action. It is not an uncompensated extension of the
reference's linear tangent to the whole exterior.

For any r_h>R, H_h>0. Put y=r-r_h. The regular horizon data for (5) are

```
F(r_h)=0, F'(r_h)=k_h>0,
l_h=S(r_h,0)/(2k_h).                                (7)
```

There is a unique local regular solution. To see the nature of the
singularity write F=y f with f_h=k_h. Its equations become

```
y f'=A-f-y f l,
y l'=[S-(r y f/D+2A)l]/f.
```

At y=0 their regular-singular Jacobian has eigenvalues -1 and -2;
the f-to-l off-diagonal derivative is zero because S_h-2A_h l_h=0.
The coefficients are analytic and f_h>0. The recursive analytic solution
has invertible matrices nI-J at every n>=1, and the ordinary convergent
power-series majorant gives a unique analytic regular branch. Equivalently
its Volterra kernels have the damping exponents one and two. The
singular independent solutions are excluded. This also gives continuous
dependence on r_h. All coefficients are finite for each r_h>R. The global
bounds below supply continuation; no uniform smallness of
(R-b)/(r_h-R) is assumed.

Normalize alpha(infinity)=0 after solving (5). The exact horizon clock
law, from a regular optical expansion, is

```
kappa(r_h)=exp(alpha_h) H_h/[2r_h D_h].               (8)
```

For example set r=r_h+B0 exp(2kappa x)+..., and
psi=kappa x+psi_h+A0 exp(2kappa x)+.... Equations (4) give
`B0/exp(2psi_h)=[2r_h U_h-gamma U'_h]/(16kappa^2 D_h)`;
using exp(2psi)=a^2 F b0^2/r^2 gives (8). Thus the period of the
regular Euclidean bolt is the stipulated laboratory period, rather
than the unchanged surface gravity at the old classical horizon.

## Global continuation and the small lapse

The benchmark parameters are unchanged:

```
k=2000, q=k^4, eta=1/2000,
xi=b^2/L^2, xi(1+3xi)=pi k^6, b^2/G=xi k^10,
C_b=1+6xi, nu=2q eta G C_b/b^3,
nu/4<=kappa_target<=7nu/25.
```

Let epsilon=gamma/R^2, h=(r_h-R)/R and
L_h=1+log(1/h). The exact root equation gives
`1<=R/b<1+1e-30`, `C_R/C_b<1+1e-30`, and
`epsilon<1/(36k^9)<6e-32`. For the shooting bracket use

```
h_low=kappa_target R/(2C_R),
h_high=2kappa_target R/C_R.
```

Both are positive, smaller than 1e-30, and
`h_low>=eta/(16k^9)`. At this benchmark `log(1/h)<80`, hence L_h<81.
These logarithmic numerical bounds are only asserted at k=2000.

The two exact integrating-factor identities underlying continuation are

```
(a F sqrt(D))'=a H/(r sqrt(D)),
(l a^2 F^2 sqrt(D))'=a^2 F sqrt(D) S.                 (9)
```

Both horizon integration constants vanish on (7). Set

```
Fbar(r)=D(r)^(-1/2) integral_(r_h)^r H(s)/(s sqrt(D(s)))ds.
```

If V_a(r)=integral_(r_h)^r |l|<=1/10, (9) compares F to Fbar
between exp(-V_a) and exp(V_a). Compare Fbar in turn to the classical
RN-AdS function F_* whose extremal center is R and whose horizon is r_h,
with effective charge Q_*=R^2+3R^4/L^2. Positivity of the kernel and

```
H=(r^2-R^2)[1+3(r^2+R^2)/L^2-6gamma/L^2]
```

give

```
exp(-V_a)(1-2epsilon) <=F/F_*<=exp(V_a)(1+6epsilon).  (10)
```

In particular F remains positive. No shift of the old horizon has been
expanded in units of the tiny thermal gap.

Here is a global scalar barrier proof with deliberately generous constants.
Use dimensionless t=r/R, xi_R=R^2/L^2>=1, Delta=t-1, and
`ell=R l`, `mathcal A=R A`, `mathcal S=R^2 S`.
Bootstrap `.99<=F/F_*<=1.01` and V_a<=1/10. These hold initially at
the regular horizon for the tiny epsilon above. Let T=F_*'/F_*,
where the derivative in this paragraph is in t. The extremal reference is

```
F_ext=(t-1)^2 P_ext(t)/t^2,
P_ext(t)=xi_R(t^2+2t+3)+1,
T_ext=2/Delta+P_ext'/P_ext-2/t.
```

The positive nonextremal subtraction is a constant divided by t, so
`T>=T_ext`. The exact identity

```
[(3/2)T_ext-1/t-1/Delta] t Delta P_ext
                 =xi_R(t^3-t+12)+4-2t >=14
```

follows from xi_R>=1 and `(t-1)^2(t+2)>=0`. Thus
`2T-1/t-1/Delta>=T/2`, and T>=1/t.

Write `a_*=F_*'+F_*/t` for the dimensionless reference forcing,
and `a_q=R k`. Exactly

```
a_q/a_*=(1-theta)/(1-3epsilon/t^2),
theta=6epsilon/[1/xi_R+3(t^2+1)]<=2epsilon/t^2.
```

Hence `1<=a_q/a_*<=1+4epsilon/t^2`. The bootstrap and
`R r/D=t/(t^2-3epsilon)` imply

```
.98T<=mathcal A/F<=1.022T.
```

For `Lbar(t)=K epsilon/[t^2(t-1)]`, K=1e6, direct differentiation and
the preceding positive identity give

```
Lbar'+[t/(t^2-3epsilon)+2mathcal A/F]Lbar
                         >=.1(mathcal A/F)Lbar.      (11)
```

The source bound needed here is uniform over the whole exterior. The
bootstrap gives `F<=1.01F_ext<=1.01(xi_R t^2+1)`; the exact center
gives `Q/R^2<=1+3xi_R`. Substitution into the numerator N of S gives
`|N/R^2|<(16.02+21.12epsilon)xi_R t^4<16.1xi_R t^4`.
Consequently `|mathcal S|<50epsilon xi_R/t^2`.
Also `mathcal A>=.97F_*'` and `F_*'>=F_ext'>=2xi_R Delta`, so

```
|mathcal S|<=5000epsilon mathcal A/[t^2(t-1)].
```

This much enlarged bound is already below the right side of (11)
after multiplying the ell equation by F. At the horizon
`ell_h=mathcal S_h/(2mathcal A_h)` is dominated by Lbar(t_h).
Scalar comparison for the exact linear ell equation therefore proves

```
|R l(r)|<=K epsilon/[t^2(t-1)], t>=1+h.              (12)
```

The barrier includes the near-horizon layer and improves to r^-3 in
the outer region. Its exact integral is

```
V_a(infinity)<=K epsilon[log((1+h)/h)-1/(1+h)]
              <=K epsilon[1+log(1/h)]<1e-20.         (13)
```

Equation (10) now improves the one-percent bootstrap strictly, and
(13) improves 1/10 strictly. The usual first-exit argument therefore
precludes either exit. Formula (10) precludes an additional zero of F
on any finite outer interval; (5) then has smooth bounded coefficients
on compact subintervals. This proves continuation to every finite r and
the existence of a finite alpha(infinity). Subtract that constant to
impose alpha(infinity)=0. It follows that `sup |alpha|<1e-20`, meeting
the predeclared component target 1e-8. The quantitative proof covers
the nonlinear system (5), not only its first derivative in a bookkeeping
coupling.

Local continuous dependence on r_h and the uniform integrable tail (12)
make alpha_h continuous throughout the bracket. Since
`H_h/[2r_hD_h]=(C_R h/R)[1+O(h+epsilon)]`, elementary bounds on this
bracket place kappa(r_h) between .9 C_R h/R and 1.1 C_R h/R.
Equation (8) is therefore below kappa_target at h_low and above it
at h_high. The intermediate value theorem supplies at least one
solution with the prescribed temperature. Uniqueness is not claimed.

## Finite asymptotic mass and a regular quantum state

The coarse relative bound (10) alone would not prove a finite ADM mass.
The asymptotic improvement is essential. There is an exact elementary
formula for the zero-lapse comparison. With u=sqrt(r^2-3gamma),
u_h=sqrt(r_h^2-3gamma), define

```
M0=u_h^3/L^2+(1+3gamma/L^2)u_h
             +(Q+gamma)/sqrt(3gamma) atan(sqrt(3gamma)/u_h).
Fbar=r^2/L^2+1-M0/u
             +(Q+gamma)/(sqrt(3gamma)u) atan(sqrt(3gamma)/u).
                                                               (14)
```

To derive it change variables in its positive integral:
`H/(r sqrt(D))dr=[3u^2/L^2+1+3gamma/L^2-(Q+gamma)/(u^2+3gamma)]du`.
Thus Fbar=r^2/L^2+1-M0/r+O(r^-2); its constant is exactly one.

From (12) and the first equation of (5), initially F=r^2/L^2+O(1).
In the explicit polynomial S the r^4/L^2 terms cancel, giving
`S/F=O(r^-6)`, while the damping coefficient is
`r/D+2A/F=5/r+O(r^-3)`. Its scalar integrating factor gives
`l=O(log(r)/r^5)`. Feed this integrable improvement into the F equation:

```
F=r^2/L^2+1+M1/r+O(log(r)/r^2).
```

The constant one is forced by its exact asymptotic coefficients. This
improves the numerator of S to O(r), so `S/F=O(r^-7)`. The same
integrating factor now gives

```
l=O(r^-5), alpha=O(r^-4),
F=r^2/L^2+1-2G M_infinity/r+O(r^-2).                 (15)
```

All estimates follow on the already constructed global positive branch;
this is not an assumed asymptotic boundary condition that replaces a
shooting result. In particular the exact spherical mass constraint
`m'=rho2` has a finite limit and a convergent exterior integral. Its
integration constant is fixed by F(r_h)=0. The sphere-integrated radius
source (2) accompanies this stress throughout.

Regularity at (7) together with (8) supplies the smooth Euclidean bolt
at beta=2pi/kappa_target. The conformal LLL covariance specified above
therefore is an actual HH/KMS state on this isolated stationary
solution, not a HH tensor at a different, unchanged classical surface
gravity. The outer static time has a(infinity)=1. Boundary-local
energy conventions and the complete four-dimensional boundary Weyl
transfer are separate from the finite geometric mass coefficient in
(15); they are not silently equated to it.

## Why classical fixed-horizon perturbation was inadequate

For comparison only, evaluate (3) on the classical thermal RN metric.
At r_0 its horizon density is
`rho_h=q[F''_h-4kappa/r_0]/(48pi)` and its lapse integral is finite,
of order `(Gq/b^2)log(1/h)`. Holding r_0 fixed gives

```
delta kappa/kappa approximately -1/(48pi z_i eta),
```

whose magnitude is between 47 and 54 at eta=1/2000. This was a
diagnostic of misplaced centering, not a physical failure. A formal
first-order fixed-temperature solution can be written completely:

```
alpha1(r)=-(Gq/(48pi)) integral_r^infinity A(s)/[s F0(s)^2]ds,
r1=[G rho_h/r_0-kappa alpha1(r_0)]/[F0''(r_0)/2+kappa/r_0],
m1(r)=kappa r_0 r1/G+integral_(r_0)^r rho0(s)ds.
```

The density in the last line is (3) and its outer integral is finite.
In the horizon-aligned radius coordinate its metric variation is
`F1bar=r1 F0'-2kappa r_0 r1/r-(2G/r)integral_(r_0)^r rho0`, with
`F1bar(r_0)=0` and `F1bar'(r_0)=-2kappa alpha1(r_0)`.
The leading center shift is r1=Gq/(24pi b). Expanding the lapse at
the old horizon can nevertheless have a nonuniform higher-order term
because this shift exceeds the thermal gap. The existence proof
(5)-(13) uses the exact center and regular horizon equations instead.
It does not infer its nonlinear conclusion from that formal Taylor
expansion.

## What this component does and does not match

The complete isolated Polyakov-LLL contribution has a consistent state,
subtraction allocation, nonzero radius source, exact nonlinear
constraints, horizon recentering, laboratory temperature, finite
geometric mass and small normalized lapse. The benchmark predeclared
component target passes by (13), before supplementary verification.

The missing W4, massive-angular and remaining local/state/boundary
functionals need not be small compared with the reference displacement.
Their action and source variations must be restored in the common scheme
and controlled on this exterior. This calculation cannot be added to the
full product shift as a second independent tadpole. It does not transfer
S138 automatically to the complete corrected operator, bound driven
pointwise stress, determine the signed finite-cut beta remainder, or
certify common-drive capture. It supplies no physical opening,
stability, transport or preparation claim and no novelty verdict.

## Focused verification contract frozen before execution

The ordinary argument above was independently reviewed before the
following supplementary checks were selected. Their execution outcomes
are recorded separately in the combined
[receipt](../evidence/stationary-quantum-exterior.json) and
[log](../evidence/stationary-quantum-exterior.log), using
[the checker](../checks/stationary_quantum_exterior.py). No numerical
example replaces the global continuation, horizon shooting or infinite
tail argument.

Require exact symbolic residuals to vanish for the action equations,
conserved HH constraint, physical radius/stress Ward identity, optical
elimination, horizon seed and clock law, center, Fbar primitive and
asymptotic cancellation, barrier polynomial and finite-coefficient
comparison. Exact rational inequalities verify the conservative
barrier/source constants, benchmark center bounds and target margin.
Record the conformal-domain and boundary checks in the same evidence
set while retaining their separate scope.

For independent ODE diagnostics set R=1, epsilon=1e-10, xi_R=3,
and h=1e-3,1e-4. Integrate (5) in log gap
`tau=log((r-r_h)/h)` through r=20, using the scaled variables
`p=l/epsilon`, `w=(F-F_*)/epsilon`, and the integral of p. Use analytic
second-order horizon seeds, with offset/h=1e-6 on the first run and
1e-8 on refinement. DOP853 tolerances are respectively
rtol=2e-9, atol=2e-11 and rtol=2e-11, atol=2e-13.

Compare the independently evaluated exact HH constraint using 80-digit
arithmetic against the ODE, with acceptance bounds 2e-4 and 2e-5 for
`abs(constraint)/(epsilon max(1,F^2,kappa^2))`. Compare the two
refinements separately for p, w and integral_p, requiring
`abs(coarse-refined)/max(1,abs(refined))<2e-5` at the frozen sample gaps
`.01h,.1h,h,10h,.01,.1,1,4,19-h`. Verify F>0 and the analytic lapse
barrier there. These synthetic checks test the actual isolated equations,
regular horizon seed and conserved state constraint; they are not a
numerical full-Dirac exterior or a discretized proof at the benchmark.

Evaluate the benchmark root and the fixed-temperature bracket at 80 and
110 digits, with scaled center-root residual below 1e-60. Preserve raw
comparison errors, refinements, input source hashes, precisions, exact
command/cwd and terminal status. Keep every failed attempt with its
exact source snapshot and do not retune the contract. Sources and outputs
are UTF-8/LF; input hashing is explicitly raw-byte SHA-256 so inherited
inputs retain their original bytes. No Lean or cloud execution is part
of this contract. The receipt and current planning audit carry the
verification outcome; the immutable note does not presume a pass.

Remaining to-do list: restore the full common four-dimensional stationary
sources and boundary charge; transfer the
state/operator and pulse theorem to the corrected evolving geometry;
complete the S074-S079 gates and S080 claim audit.
