# Exact product-throat cutoff contacts and clock variation

2026-10-15, E003/S133, supporting S074. This resolves a boundary-model
distinction, not the exterior matching problem. We derive the finite-cut
action first and then bound its first and second clock variations uniformly.
A bounded finite coefficient of each displayed contact is suppressed on
this throat cut. A coefficient in a separately renormalized exterior
boundary action has not thereby been measured or equated to it.

## Scope and acceptance fixed before execution

Use the product metric and normal convention of
[the radiation note](source-clock-radiation.md), the source contacts and
finite scheme term of [the curved-boundary note](curved-fermion-boundary.md),
and the distinct finite-tensor boundary assumed in
[the joint-clock note](source-clock-joint.md). This calculation addresses
the unproved identification between those boundary descriptions. It
keeps the product metric fixed; it does not solve Einstein constraints,
the interacting Dirac state, the full Weyl anomaly or an exterior match.
The geometric matching cutoff epsilon is distinct from the ultraviolet
proper-time or spectral regulator. The local EFT coefficients here are
finite coefficients in the declared subtraction scheme at a fixed scale;
removing an ultraviolet regulator does not supply their dependence on
the geometric cutoff or identify them with exterior coefficients.

Before running the checker, freeze the following tests. Exact symbolic
identities must establish the induced lapse, an independent Christoffel
calculation of extrinsic curvature, all three contact densities, their
directional variations and the fixed-lapse distinction. Two prescribed
smooth polynomial/sinusoidal families on 0<=u<=1, at epsilon=1/5,1/10,1/20,
are checked by independent integration of the exact action at displaced
curves. Use 50-decimal arithmetic. Central differences at steps 2e-4 and
1e-4 must agree with the analytic integrated first derivative within 1e-8
and both pure and mixed second derivatives within 1e-6. Independent
geometric densities must agree within 1e-35. Bound checks permit 1e-30
rounding allowance. These synthetic families test the derivation; their
parameters are not assigned to a physical clock. Cutoff suppression is
proved below, not inferred from those three numerical cutoffs.
For the two test families and |s|,|t|<=1, elementary sine/cosine bounds give
|F|<=6/5, F'>=87/100, P>=22/25, and all required source/geometry jets
bounded by B=2. Thus v=p=4/5 and epsilon0=1/5 provide a common analytic
neighborhood; the finite grid is only a supplementary check of it.

## Exact finite-cut geometry and action

Write the four-dimensional product metric and a timelike curve as

```
ds^2=r_e^2*(-dT^2+dZ^2)/Z^2+r_e^2*dOmega^2,
T=F(u), Z=epsilon*P(u),
x=F', y=P', f2=F'', p2=P'',
Delta=x^2-epsilon^2*y^2>0,
N=r_e*sqrt(Delta)/(epsilon*P), gamma_uu=-N^2.
```

Assume F'>0, P>0 and epsilon>0. The outward unit normal points toward
decreasing Z. Its (T,Z) components are
`n=-(Z/(r_e*sqrt(Delta)))*(epsilon*P',F')`.
For `g_2=r_e^2 diag(-1,1)/Z^2`, the nonzero Christoffels are
`Gamma^T_TZ=Gamma^T_ZT=Gamma^Z_TT=Gamma^Z_ZZ=-1/Z`.
The checker derives these from metric differentiation rather than
inserting the proposed curvature formula. With tangent e=(x,epsilon*y),
the outward extrinsic-curvature trace is

```
K=n_mu*(e'^mu+Gamma^mu_ab*e^a*e^b)/N^2
 =[x*Delta+epsilon^2*P*(y*f2-x*p2)]/(r_e*Delta^(3/2)).
```

For a static cut it is 1/r_e. The sphere extrinsic curvature is exactly
zero because its radius is constant in both product coordinates. In the
inward heat convention `L_ab=-<nabla_a n_in,e_b>`, this gives L=K and
`Lhat_uu=(2/3)*K*gamma_uu`.

Set theta(u)=Theta(F(u)) for a fixed smooth global-time source profile.
Thus theta'=Theta'(F)*F'; the source is varied through F rather than
held fixed as a function of u. Sphere integration gives exact u-densities

```
L_s     = 2*pi*zeta_s*r_e^2*theta'^2/N,
L_c     =-(8*pi/3)*c_theta*r_e^2*K*theta'^2/N,
L_scheme= (4*pi/3)*d_theta*r_e^2*K*theta'^2/N.
```

Here `I_s=-(zeta_s/2)*integral sqrt(-gamma)*(Dtheta)^2`,
`I_theta=integral sqrt(-gamma)*Lhat^ab*theta_a*theta_b`, and the last
density is the displayed finite scheme choice `-d_theta*J_theta/3`,
`J_theta=integral sqrt(-gamma)*L*(Dtheta)^2`. The Lorentzian contact
conventions are those retained in the cited notes. This is not a
continuation or derivation of the full quantum Weyl functional.

For `r=theta'`, `A=P*r^2`, and
`M=x*Delta+epsilon^2*P*(y*f2-x*p2)`, the full displayed density is

```
L_total=epsilon*[2*pi*zeta_s*r_e*A*Delta^(-1/2)
                  +(-8*pi*c_theta/3+4*pi*d_theta/3)*A*M*Delta^(-2)].
```

This identity holds before any clock variation or cutoff expansion.
Coefficients are held fixed under the clock variation. If they depend
on the clock or controller, their derivatives are additional inputs.

## Uniform action and C2 variation bounds

Let J be a finite u-interval. Consider a compact smooth family with a
common positive lower bound `F'>=v`, `P>=p`, `0<epsilon<=epsilon0<=1`,
and `Delta>=d=v^2/2`. For example, `epsilon0*sup|P'|<=v/sqrt(2)`
ensures the latter. On the source image and a common neighborhood require

```
|P|,|F'|,|F''|,|P'|,|P''|,|Theta'|,|Theta''|,|Theta'''|<=B, B>=1.
```

Relevant variations have a common C2 bound. For affine variations
`F+s*h+t*k`, `P+s*j+t*l`, define
`H=max(||h||_C2,||j||_C2)`, `Jv=max(||k||_C2,||l||_C2)`,
where the C2 norm is the maximum of supremum norms through order two.
The source composition is included in the bounds. A C3 bound suffices
when imposing P=F', so its two spatial derivatives remain controlled.
Bounds must hold in a neighborhood of the varied curves, not just on
one unvaried path. No rapidly oscillating variation with unbounded jets
is admitted.

The following constants bound values, first derivatives divided by H,
and mixed second derivatives divided by H*Jv, respectively:

```
r0=B^2,                 r1=2*B^2,       r2=3*B^2,
A0=B^5,                 A1=5*B^5,       A2=22*B^5,
Delta1=4*B,             Delta2=4,
M0=3*B^3,               M1=11*B^2,      M2=24*B.
```

Proof. Differentiate `r=Theta'(F)*F'` twice. Its mixed second derivative
is `Theta'''(F)*h*k*F'+Theta''(F)*(h*k'+k*h')`, giving r0,r1,r2.
The product `A=P*r^2` gives A0,A1,A2. For positive Delta<=B^2,
differentiate its quadratic expression to obtain Delta1,Delta2.
Differentiate `M=x*Delta+epsilon^2*(P*y*f2-P*x*p2)`; the first term
has bounds B^3,5*B^2,12*B, and the remaining two cubic terms together
have bounds 2*B^3,6*B^2,12*B. These prove the displayed constants.

For alpha>0 define bounds on `Delta^(-alpha)` by

```
R0(alpha)=d^(-alpha),
R1(alpha)=4*alpha*B*d^(-alpha-1),
R2(alpha)=16*alpha*(alpha+1)*B^2*d^(-alpha-2)
            +4*alpha*d^(-alpha-1).
```

For a triple `X=(X0,X1,X2)`, define its product with Y by
`X*Y=(X0Y0, X1Y0+X0Y1, X2Y0+2X1Y1+X0Y2)`.
Put `T=A*R(1/2)`, `G=(A*M)*R(2)`, and

```
Cj=2*pi*|zeta_s|*r_e*Tj
       +(8*pi*|c_theta|/3+4*pi*|d_theta|/3)*Gj, j=0,1,2.
|S|<=epsilon*|J|*C0,
|DS[h,j]|<=epsilon*|J|*C1*H,
|D2S[(h,j),(k,l)]|<=epsilon*|J|*C2*H*Jv.
```

Uniform coefficient upper bounds can replace their displayed magnitudes.
These are explicit sufficient finite-cut ledger bounds on this local
action and its variations. They do not yet bound a clock solution or
area observable; that requires the matched constraint-response map.
They are not uniform in magnetic flux, r_e or the small-clock parameters
unless the displayed constants and normalization bounds are tracked in
those limits. In particular F'=y~eta*z may have a small lower bound v,
which worsens d and the Cj. The synthetic v=4/5 is not the opening
benchmark's lower bound. No S076 numerical margin follows from it.
For a nonlinear C2 path X(s), its second derivative additionally contains
`DS[X''(s)]`; the same C1 bound controls that term. No interchange of
a limit and unbounded differentiation is used: the exact differentiated
density is bounded by an epsilon factor on the declared neighborhood.

## Asymptotic fixed lapse is not exact fixed lapse

For P=F', the exact formulas reduce to

```
N=(r_e/epsilon)*sqrt(1-epsilon^2*(F''/F')^2),
K=(1-epsilon^2*F'''/F')
    /[r_e*(1-epsilon^2*(F''/F')^2)^(3/2)]
 =(1/r_e)*[1-epsilon^2*Sch(F,u)+O(epsilon^4)],
Sch(F,u)=F'''/F'-(3/2)*(F''/F')^2.
```

The expansion is uniform with bounded jets and F' bounded away from zero.
This recovers the inherited curvature expansion. But P=F' matches
`N=r_e/epsilon` only to leading order. Its exact density has

```
r_e^2*K/N=epsilon*[1-epsilon^2*
                    (F'''/F'-2*(F''/F')^2)+O(epsilon^4)].
```

Replacing that density's cubic correction by the exact-fixed-lapse one
would omit the lapse correction. Exact fixed lapse instead requires
`P^2+epsilon^2*P'^2=F'^2`; a controlled regular expansion would have to start
`P=F'-epsilon^2*F''^2/(2*F')+O(epsilon^4)`.
One must solve or control that condition and its variations before using
it as an exact gauge. The general finite-cut bounds above apply without
imposing it. We make no global existence claim for that embedding ODE.

## Matching decision

With bounded coefficients, regular source rates and uniformly controlled
geometry, all three throat-cut actions and their first/second variations
vanish as O(epsilon). This includes zeta_s, not only c_theta. A nonzero
finite limit requires violating at least one of these hypotheses: for
example coefficient growth of order 1/epsilon, source rates of order
epsilon^(-1/2), degenerating geometry/variation bounds, or an independently
matched exterior boundary action. Such scalings are different regimes
and require their own renormalization and error analysis. They are not
inferred from a finite coefficient on this cut.
Moreover the product throat is itself an approximation valid over a
matching range. Sending epsilon to zero in that formal product does not
prove the approximation remains valid, or allow the physical matching
surface to be moved outside that range. Its exterior matching error and
normalization must be bounded before this result can enter an opening
or area error budget.

Thus the finite arbitrary tensor in the S084 aggregate clock is not
identified by the S083 product-throat geometry. S084 remains a conditional
ODE window; the product-throat suppression does not assign a physical
exterior coefficient or demonstrate that it vanishes there. In particular
the zeta_s in a renormalized conformal-boundary action need not equal a
bounded throat-cut coefficient. The exact cutoff map fixes neither value.

The physical next input is a derived exterior-to-throat boundary action
and canonical clock/field map, including regulator, reference subtraction,
source/controller implementation, and a controlled limit of its clock
variations. The boundary half-density relation in
[boundary-vertex-matching.md](boundary-vertex-matching.md) is a retained
component, not that complete matching. No coefficient is fitted here.
No full S074 pass, new Lambda interval, physical opening, stability,
transport, preparation or creation follows. No Lean or paid compute is used.

Remaining to-do list: derive the exterior/throat action and clock map;
match its finite source coefficients and geometric variations; control
the physical determinant response and remaining S074 remainder.
