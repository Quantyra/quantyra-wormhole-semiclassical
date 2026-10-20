# Complete equation defect on the retained comparison path

2026-10-20, E003/S161 supporting S139/S074. Ordinary proof and
predeclared verification before acceptance.

The original comparison geometry has a complete retained-equation
defect in fixed-chart mass, lapse, angular-profile and horizon norms.
The new result is the global independent angular equation estimate
and its combination with S160. The sharper earlier mass, lapse and
horizon defect bounds are reused, not replaced by the much coarser
cold source Hessian allowance. This supplies no inverse or solution
error, and no physical completion gate.

## 1. Same action, same path, independent radius variation

Keep the exact spherical action split, reference potential, fixed
clock, prepared period, state/domain and common scale of
[the centered seed](stationary-centered-covariant-seed.md).
The retained model context is
[Bintanja--Freivogel--Rolph](https://arxiv.org/abs/2308.00871v2),
inspected 2026-10-20; that paper does not supply this residual bound.
Write

    S_lambda=S0+lambda Gamma_rest,
    sigma_lambda=sigma0+lambda u,
    t_lambda=1+(t-1)exp(lambda V), t=r0/B,
    0<=lambda<=1, C=B^2/A^2.

All observations are in the fixed seed disk. Let Delta and the
gradient be for its dimensionless product metric ghat0/A^2, and
R0 its dimensionless curvature. Define

    Q(t)=t F0(t)/C,
    P(t)=t^3 Q''(t)/2
        =3(x/C)t^4+(1+3x+d)/C-6(d/C)/t+6(d/C)/t^2.

Vary the reduced action at fixed PHYSICAL orbit metric:

    E_X=R2-2 Box2 r/r+V0'(r^2).

Since pS=J_X/(16piG), the angular geometric source is -E_X/2
in the 8piG pS normalization. The conformal relation
g2=(r/B)^2 ghat gives, in S160's normalized angular units,

    A0(lambda)=t_lambda^2 exp(-2lambda u)
      [-R0/2+lambda Delta u+2 Delta t_lambda/t_lambda
                       -|grad t_lambda|^2/t_lambda^2]
                  -P(t_lambda).                         (1)

Equivalently the local combination is
t_lambda^2[-R_lambda/2+2 Box_lambda phi_lambda
+|grad phi_lambda|^2]-P(t_lambda).
The coefficient of the gradient-square term is one. The independent
variation is necessary; its bolt value alone would not distinguish
this coefficient.

The seed has A0(0)=0. The actual tangent solves the complete
linear equations, including the radius equation by the joint Ward
identity, so

    A0'(0)=A_*,  A_*=A_rest(g0).

The complete angular residual is therefore

    R_S(lambda)=A0(lambda)-lambda A_rest(g_lambda)
      =[A0(lambda)-lambda A_*]
          -lambda[A_rest(g_lambda)-A_*].                 (2)

This is the source normalization of an equation defect, not a
linearized solution error. No source prefactor is counted twice.

## 2. Exact nonlinear identity without a varying stress gradient

Set s=t-1, tau=s/t, nu=tau V, k=exp(lambda V),
T=1+tau(k-1), a=exp(-2lambda u), and

    L=Delta+grad(log t).grad,
    B_lambda=k-1-lambda V k.

Then

    L T=lambda k L nu+B_lambda L tau
                                  +tau lambda^2 k |grad V|^2,
    grad T=lambda k grad nu+B_lambda grad tau.

Expanding (1), its tangent is

    A_*=(2nu-2u)P+t^2 Delta u+2t^2 Lnu-t P' nu.

Eliminate Lnu by this exact first-order equation. The complete
geometric angular remainder N_S=A0-lambda A_* is

    N_S=lambda(aTk-1)A_*
       +P{a[T^2+2lambda Tk(u-nu)]-1}
       +lambda aTk t P' nu-[P(tT)-P(t)]
       +lambda a tT(1-k)Delta u
       +2a t^2 T B_lambda Ltau
       +2a t^2 T tau lambda^2 k |grad V|^2
       -a t^2 |grad T|^2.                                (3)

This uses the known seed source value and local geometric jets.
No spatial derivative of the source on an unknown metric appears.

For radial scalars the seed identities are

    Delta f=t^2[(F0/C)f_tt+(F0_t/C)f_t],
    |grad f|^2=t^2(F0/C)f_t^2,
    Ltau=F0_t/C-F0/(Ct), b=s V_t,
    grad T=sqrt(F0/C)[(k-1)/t+lambda k b].                (4)

They retain the actual radius and moving-bolt geometry. In
particular the final negative gradient square in (3) must be kept.

## 3. An explicit majorant and regional bounds

In this section U=|u|, W=|V|, N=|nu|, H=|Delta u|,
and D1,D2 bound |tP'| and |t^2P''| on the segment between
t and tT. All scalar multipliers in (3) are within 1.001
of one; the 1.01 and 2.02 reserves below include them.

For 0<=lambda<=1,

    |k-1|<=1.001lambda W,
    |k-1-lambda V|<=.501lambda^2 W^2,
    |B_lambda|<=.501lambda^2 W^2,
    |aTk-1|<=1.01lambda(2U+W+1.001N).

For C_P=a[T^2+2lambda Tk(u-nu)]-1, C_P(0)=C_P'(0)=0
and differentiating twice gives

    |C_P|<=lambda^2[10(U+N)^2+5(U+N)W+2tau W^2].         (5)

For a direct bound, put n=U+N. The second derivative of
log T has magnitude at most 1.001 tau W^2. The half-second
derivative bounds for aT^2 and 2lambda aTk(u-nu) sum to
at most

    2.02n^2+1.01tau W^2
     +2.01n(2n+W)+1.01n[(2n+W)^2+1.001tau W^2].

Using n,W<1e-6 gives (5) coefficient by coefficient.
Taylor's formula for P(tT) retains both the exponential radius
acceleration and its potential Hessian. Thus |N_S|/lambda^2
is bounded by the sum

    1.01(2U+W+1.001N)|A_*|
    +|P|[10(U+N)^2+5(U+N)W+2tau W^2]
    +1.01D1 N(2U+W+1.001N)+.501D1 tau W^2+.502D2 N^2
    +1.01t W H+1.01W^2 t^2|Ltau|
    +2.02t^3 F0 b^2/(Cs)
    +1.01t^2(F0/C)(W/t+|b|)^2.                          (6)

The original invariant curvature identity gives globally
H<=2U+5e-8. On the cold patch its sharper version gives
H<=2U+5e-8 s. These use the already proved first-order
curvature variation, not raw areal second derivatives.

Here are sufficient input rows. All inequalities include the
small displacement of the P-derivative segment.

* Cold, h<=s<=.001, a_x=s L_c(s)<=.00331:
  U<=4e-6 a_x, W<=1.2e-8, N<=sW, H<=8.05e-6 a_x,
  |A_*|<=1.7e-6 a_x, |b|<=3e-7,
  |P|<=1.01, D1<=3, D2<=10,
  F0/C<=1.01s^2/t^2, t^2|Ltau|<=5s.
* Middle through t=4:
  U,W,N<=1.2e-8, H<=7.4e-8, |A_*|<=2.6e-9,
  |b|<=3e-7, F0/C<=s^2, t^2|Ltau|<=3t^3.
  Also |P|<=.51(t^4+1), D1<=2.1t^4, D2<=6.3t^4.
* t>=4:
  U<=8e-9/t^2, W,N<=2e-9/t^2,
  |b|<=1.4e-8/t^2, H<=7.4e-8,
  F0/C<=5t^2/5.9, t^2|Ltau|<=8t^3,
  with the same P bounds. Use |A_*|<=4.7e-7 L_e^2.

The curvature, low outer jets and total-variation inputs are
from [the E0 constraint proof](stationary-centered-constraint-residual.md)
and [the completed metric proof](stationary-completed-metric-response-norm.md).
The cold fixed-chart U and source rows are from S160.
The simple potential rows follow directly from P and x/C<1/5.9.
For the middle F bound, factor F_ext:
[1+x(t^2+2t+3)+d/t^2]/(Ct^2)<=1 for t>=1.
The thermal subtraction only lowers it. The derivative of this
ratio is negative, including the retained |d|/C<1e-24 margin.
The Ltau bounds follow from (4) and the seed F derivative bounds;
in the cold region F0_t/C<2.1s and t<=1.001 suffice.

In the cold row use s/a_x<=1, s<=.001 and a_x<=.00331
before summing (6), so the bound is uniform down to the bolt.
For the middle row all upper expressions are monotone up to
t=4. In the last row use s/t>=3/4 and t>=4.
These substitutions, with exact rational arithmetic, give

    |N_S| < lambda^2 times
      2e-12 a_x             (cold),
      6e-11                 (middle through t=4),
      2e-15 L_e^2           (t>=4).                      (7)

For 4<=t<=40, L_e<5, so the last bound is smaller than
the middle allowance. No maximum over sampled radii is used.

## 4. Complete residual and regularity

Integrating S160's angular first-response bounds along the
ORIGINAL path, and adding (7), proves

| Angular residual region | Bound for abs(R_S) |
| --- | ---: |
| cold | `2.7e-10 lambda^2 a_x` |
| middle, 1.001<=t<=40 | `6.3e-11 lambda^2` |
| outer, t>=40 | `5e-11 lambda^2 L_e^4` |

In particular the independent angular equation is included globally.
No division by r_chi at the bolt or by a vanishing angular source
has been used to infer it from the constraints.

For the remaining components reuse the stronger already verified
[complete candidate defect](stationary-centered-complete-residual.md):

| Complete retained-equation output | Bound |
| --- | ---: |
| signed zero-datum mass primitive / r0F0 | `2.1e-11 lambda^2` |
| absolute integrated lapse-constraint density | `2.6e-10 lambda^2` |
| each physical horizon pressure / prepared thermal force | `2e-13 lambda^2` |

These concern precisely the same finite completion, ensemble and
fixed-chart source allocation. Their predecessor receipt has 92
passed checks and sixteen matching current input identities.
Its finite-endpoint source-difference proof is compatible with
S160's actual-path C2 result. Neither requires differentiating
the original path in a moving C5 proper chart.

The stronger horizon bound takes precedence over the much coarser
bound obtained by dividing the global cold row by the thermal gap.
The earlier absolute source-change mass contribution is below
1e-13 lambda^2; the E0 mass-density remainder is absolutely
integrable, but no new numerical absolute-energy bound is assigned.
None of these geometric mass statements is a Hamiltonian-charge theorem.

The full residual is C2 in the compatible fixed density/primitive
spaces and the cold/middle/outer angular space (outer weight L_e^7).
For E0 the formulas are analytic in lambda on a complex neighborhood
of [0,1]. At the bolt use N_q(t)-N_q(v) before dividing by Q(t);
its integral formula keeps the h-weighted limit finite. The lapse
uses the total variations already proved. At infinity, the bounds
on V,b,u and their first radial derivatives give an integrable
t^-2 mass-density and t^-5 lapse-density majorant, also for the
first two parameter derivatives. Formula (3) gives a uniform
angular L_e^2 majorant on a slightly enlarged parameter disk.
Thus dominated convergence applies to the fixed maps. S160
supplies the source C2 statement, including its quarter-power
bolt modulus and weaker outer topology.

This regularity statement is not a small inverse bound.
In particular S160's cold lapse second-response allowance
3.61e-5 is much larger than the sharper finite-endpoint defect.
The two assertions use different estimates; the larger allowance
does not invalidate the smaller defect or prove solver preservation.

## 5. Next mathematical requirement

The remaining problem is a usable full inverse and a general
admissible family on which the complete nonlinear map closes,
with the actual state/domain, horizon, clock and return conditions.
The inverse must recover the spatial regularity needed by the
retained quantum operator; a low constraint integral alone does not.
The angular equation must remain in that system. The prescribed
comparison path and its residual do not certify a nearby exact branch.

## 6. Predeclared verification

Revalidate the S160, S152 and earlier complete-defect receipts and
all exact input identities. Derive the radius Euler equation from
the reduced action independently of the orbit equations. Check the
seed angular identity, exact finite decomposition (3), scalar
Taylor majorants, all regional sums and the complete allocation.

Use 100 digits for manufactured exact-potential fixtures. Choose
e=1/C=.001, d/C=.00002, x/C=(1-e-d/C)/6 and thermal
h=1e-5 or 1e-30. Build the exact seed root v=1+h.
Use u=.02(t-1)/(1+t)^3 and V=.02/(1+t)^3.
These enlarged test fields are not the retained quantum tangent.
Define their manufactured A_* by the actual linear derivative.

At s=t-1 equal to 2h,.001,.1,3,49 and lambda=0,.4,1:
compare (1) with the independent areal-coordinate angular Einstein
formula, and compare (3) with A0-lambda A_*.
Normalize by .02^2 min(s,1), and require identity tolerance
1e-30(1+abs(normalized target)). Compare symbolic second responses
with 64/128-node Cauchy contours of radius .01 and real steps
.001,.0005. Contour tolerance is 1e-30(1+abs(target)); real
tolerance is 1e-6(1+abs(target)), refinement ratio .35 or floor
1e-35. Endpoint differences use second-order one-sided stencils
inside [0,1].

Add a manufactured source with independent first/second response
to test (2) and the quadratic cancellation. Check that dropping the
gradient-square term gives a detectable error. Check the physical
angular conversion and the mass/lapse prefactors exactly once.
These fixtures test local algebra and calculus, not a continuum
determinant, physical state, full inverse or nonlinear solution.
Preserve every run, raw output, failure and superseded source snapshot.

Remaining to-do list: general admissible family and usable full
inverse; coupled exterior/state, charge, return and signed-source
solve; opening, stability, transport, preparation, separate creation
and full-result audit.
