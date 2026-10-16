# RN-AdS exterior contacts, travel time and a finite throat seam

2026-10-15, E003/S134.1, supporting S074. This is a static geometric
matching result for the retained RN-AdS reference. It corrects the
equal-radius assumption when applying S133 to this candidate and supplies
explicit exterior travel and contact bounds. It does not construct an
evolving semiclassical wormhole or match a physical controller.

## Sources, scope and acceptance fixed before execution

[Bintanja, Freivogel and Rolph v2](https://arxiv.org/html/2308.00871v2),
equations (2.18)–(2.27), specifies the unequal-radius throat, RN exterior
and clock/overlap dictionary. Equations (2.34)–(2.43) split optical length
into truncated throat and exterior pieces. Our estimates below are derived
directly from the exact extremal RN factorization, rather than adopting
their regional approximations as error bounds. The local conventions and
candidate scales are already retained in
[axisymmetric-axial-action.md](axisymmetric-axial-action.md),
[boundary-vertex-matching.md](boundary-vertex-matching.md) and
[causal-compact-response.md](causal-compact-response.md).

The checker must establish the factorization and normal/extrinsic
curvatures symbolically; verify contacts independently from the metric;
and compare two changes of variable in exterior null-time quadrature.
Before execution fix 60-decimal arithmetic, quadrature agreement tolerance
1e-30, direct-geometry tolerance 1e-40, and bound rounding allowance 1e-30.
The predeclared cases are (X,delta)=(1,1/10),(100,1/100), plus the retained
k=2000 and k=4000 candidates with delta=k^-4. Signed mass envelopes
m=1/100 are tested as synthetic comparisons; they are not physical mass
assignments. The initial classical thermal reference is checked separately.
No historical checker is run. All source and checker hashes are retained.

## Exact reference geometry and unequal radii

Distinguish the AdS4 radius L4, extremal sphere radius b, and AdS2 radius a.
For X=b^2/L4^2, C=1+6X, the extremal exterior is

```
ds^2=-f_ext(r)*dU^2+dr^2/f_ext(r)+r^2*dOmega^2,
f_ext=1+r^2/L4^2-2*b*(1+2*X)/r+b^2*(1+3*X)/r^2,
a=b/sqrt(C).
```

The extremal near-horizon product is AdS2(a) times S2(b). S133's exact
equal-radius formulas are valid for its declared product, but a=b is not
the retained RN-AdS4 candidate. With s=(r-b)/b,

```
f_ext=s^2*J(s)/(1+s)^2, J(s)=C+4*X*s+X*s^2.
```

At the timelike static cut r_m=b*(1+delta), delta>0, use the outward
normal n=sqrt(f)*partial_r and the inward heat convention
`L_ab=-<nabla_a n_in,e_b>`. The time and sphere principal curvatures are
`kt=f'/(2*sqrt(f))`, `ks=sqrt(f)/r`, so
`Lhat_U^U=(2/3)*(kt-ks)`. Sphere integration yields

```
A(r)=r^2/sqrt(f), B(r)=r^2*f'/(2*f)-r,
L_source=2*pi*zeta_s*A*theta_U^2,
L_curvature=-(8*pi/3)*c_theta*B*theta_U^2.
```

The finite scheme term -d_theta*J_theta/3 has density
`(4*pi/3)*d_theta*(B+3*r)*theta_U^2`, since the full trace is kt+2ks.
Under a time-coordinate change U=U(t), each density becomes its coefficient
divided by U'(t), multiplying theta_t^2. In particular for a fixed
reference Killing-time rescaling t_ref=kappa0*U, A_t=kappa0*A and
B_t=kappa0*B. A coordinate
rewrite alone supplies no surface/controller variation or new dynamics.

For a product cutoff T=F(u), Z=epsilon_cut*P(u), define Delta,A_jet,M
as in [S133](cutoff-clock-contact-matching.md). Its exact density becomes

```
epsilon_cut*[2*pi*zeta_s*(b^2/a)*A_jet/sqrt(Delta)
 +(-8*pi*c_theta/3+4*pi*d_theta/3)*(b^2/a^2)*A_jet*M/Delta^2].
```

Thus source prefactors carry b^2/a and curvature/scheme prefactors carry
C=b^2/a^2. A finite O(epsilon_cut) statement at fixed radii is not
uniform in flux or candidate k without these prefactors and normalization.

## A controlled static contact match in the overlap

Let lambda0=kappa0*b/C and epsilon_cut=lambda0/delta. BFR's static
wormhole dictionary is rho=L_WH*(r-b)/b^2, t=C*U/L_WH. For comparison
of scales only, choose a reference L_ref=C/kappa0 and its
rho_ref=L_ref*(r-b)/b^2; then rho_m=delta/lambda0.
The initial thermal temperature does not determine a physical wormhole
length L_WH. We have not joined the thermal black hole and the static
wormhole into one solved geometry. The product overlap requires
lambda0<<delta<<1.
It does not permit an arbitrary cutoff limit at fixed physical parameters.

For 0<delta<=1, the exact extremal contact ratios are

```
A_t/[(b^2/a)*epsilon_cut]=(1+delta)^3/sqrt(J(delta)/C),
B_t/[C*epsilon_cut]
 =1-delta^2+X*delta*(1+delta)^2*(delta+2)/J(delta).
```

Since C<=J(s)<=C*(1+s)^2 for s>=0, the first ratio lies between
(1+delta)^2 and (1+delta)^3 and differs from one by at most 7*delta.
The second differs from one by at most 3*delta, using X/C<=1/6.
The scheme ratio (B_t+3*kappa0*r_m)/(C*epsilon_cut) differs from one
by at most 9*delta. These are static radial matching bounds with known
unequal-radius prefactors; no clock directional-variation bound follows
without the actual moving-surface and field maps.

## Exact one-sided exterior null-time bounds

In the fixed static reference, the boundary-to-cut radial null time is

```
T_ext(delta)=integral_(r_m)^infinity dr/f_ext(r)
 =b*integral_delta^infinity (1+s)^2/[s^2*J(s)] ds,
b/(C*delta)<=T_ext
 <=(b/C)*[1/delta+2*log(1/delta)-delta+C/X], 0<delta<=1.
```

For the lower bound use J<=C*(1+s)^2 everywhere. On [delta,1], use
J>=C and integrate s^-2+2*s^-1+1. On [1,infinity), use
J>=X*(1+s)^2, whose remaining integral is 1/X. These explicit bounds
also establish convergence at conformal infinity. The delay in reference
clock units is kappa0*T_ext, with leading term epsilon_cut.
It is not yet the delay in the evolving t(U) of the opening experiment.

Do not add the exterior delay to an already complete pi-long throat.
The leading global-AdS2 throat truncated at rho_m has optical width
2*atan(rho_m)=pi-2*atan(epsilon_cut). A model joining only these two
leading pieces has width W=2*atan(rho_m)+2*kappa0*T_ext and satisfies

```
0<=W-pi
 <=2*[lambda0*(2*log(1/delta)-delta+C/X)+epsilon_cut^3/3].
```

Here epsilon_cut-atan(epsilon_cut)<=epsilon_cut^3/3 proves the bound.
The leading 2*epsilon_cut contributions cancel between the truncated
throat and exterior. W is a seam model, not the deformed self-consistent
wormhole metric; perturbations of its pieces and gluing remain to be
controlled. No change to the physical free-return period is inferred.

## Signed mass perturbation and exterior horizon exclusion

At the same magnetic charge take f=f_ext-2*G*DeltaM/r and set

```
m=2*G*|DeltaM|/[r_m*f_ext(r_m)]<1.
```

The denominator r*f_ext(r)=b*s^2*J(s)/(1+s) is strictly increasing for
s>0, because its logarithmic derivative is
2/s+J'/J-1/(1+s)>0. Therefore on the whole exterior r>=r_m,

```
(1-m)*f_ext<=f<=(1+m)*f_ext,
T_ext/(1+m)<=T_mass<=T_ext/(1-m).
```

In particular f>0 there: any horizon of a positive-mass perturbation
lies below the chosen cut. This is a conditional finite mass envelope,
not a global stability or interior solution. The A ratio receives at
most the multiplicative factors (1+m)^(-1/2),(1-m)^(-1/2). At the cut,
the curvature ratio relative to C*epsilon_cut changes by at most
6*m/(1-m), obtained by differentiating the exact log f and using
delta*J'/J<=1 for delta<=1. Hence a sufficient total curvature-ratio
error bound is 3*delta+6*m/(1-m).

## Retained candidate scales and reference clock comparison

The selected family has integer k>=2000,
`g=k^-4`, `L4/l_p=k^5`, `q=k^4`, `N=2*q`, `eta=1/k`,
`X*(1+3*X)=pi*k^6`, `phi_r=b^3/(G*C)`, and
`beta0=2*pi*phi_r/(N*eta*z0)`, 1/4<=z0<=7/25. Thus

```
kappa0=2*pi/beta0=N*eta*z0/phi_r,
lambda0=kappa0*b/C=2*z0/(k^7*X).
```

For all k>=2000, X>=k^3 and C/X<=6+k^-3. Choose delta=k^-4. Then
lambda0<=14/(25*k^10), epsilon_cut<=14/(25*k^6), and the explicit
delay upper bound gives

```
kappa0*T_ext<8.76e-21,
7*delta<=4.375e-13, 3*delta<=1.875e-13,
W-pi<7.5e-32.
```

For example log(k)/k^4 decreases on this range and
2*log(1/delta)+C/X<=9*log(k), proving the stated uniform numerical
bounds at k=2000. The reference-clock exterior delay is small compared
with the nominal global switch smoothing width 1/4. This comparison
does not bound the actual evolving clock's waveform, stress or geometry.
The known prefactor C is approximately 4.912e10 already at k=2000;
S133's equal-radius constants cannot be inserted as candidate error bars.

For an additional classical thermal RN reference at T0=kappa0/(2*pi),
write its horizon r_+=b*(1+s_h). The exact temperature equation is

```
2*lambda0*C*(1+s_h)^3
 =2*C*s_h+(1+18*X)*s_h^2+12*X*s_h^3+3*X*s_h^4.
```

For X>=1 its temperature increases for r_+>=b. The candidate temperature
is below that at r_+=2b, so s_h<1. Positivity of the right-hand terms
then gives s_h<=8*lambda0. Eliminating the horizon mass yields
m<=4*(s_h/delta)^2<=256*epsilon_cut^2<2e-38 throughout the candidate
range. The finite-mass delay envelope is therefore effective on this
initial classical reference. This neither identifies an interacting
quantum state nor proves that a subsequent exterior stays in that envelope.

The retained evolving definition is u=N*U/phi_r,
dt/du=eta*z(t), so dt/dU=N*eta*z(t)/phi_r. Only its initial stationary
value equals kappa0. Therefore t_ref is only an initial tangent-time
comparison, not the actual global t(U) beyond that point. It is also
not the prescribed-source note's open-loop trajectory t0(u).
Dynamic use needs the actual t(U), exterior travel
map, controller and correlated state. The geometric matching cutoff is
distinct from the ultraviolet regulator. No finite coefficient, full
anomaly, S074 pass, new opening interval or S076 error budget is inferred.

Remaining to-do list: derive the evolving exterior/throat clock and
source action map; control geometry/gluing remainders and physical
determinant response; complete S074's finite matching and total remainder.
