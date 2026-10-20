# Variable Dirac three-point kernel in a flat orbit chart

2026-10-19, E003/S142, supporting S139/S074. Ordinary argument before
numerical verification. This supplies the complete flat-chart massive
shape three-point kernel at the unperturbed basepoint, with its seagulls,
an explicit local restoration and a summable angular remainder. It is
not yet the three-point function on the curved HH exterior.

## Comparison, conventions and acceptance

Use the rank-two complex orbit Dirac normalization of
[the product-contour proof](stationary-product-contour-inverse.md), with
`M_n^2=n(n+q)/C`, `d_n=2(q+2n)` and no additional spin or mouth factor.
[Carmi--Di Pietro--Komatsu, (6.4), (6.9)](https://arxiv.org/html/1810.04185#S6.SS1)
fix the fermionic minus trace-log and spin-traced bubble convention;
those source equations were inspected on 2026-10-19. Their interacting
Gross-Neveu gap equation and finite subtraction are not used here.

The comparator is the retained constant-shift flat remainder in
[parameter-form bounds, section 8](stationary-reduced-parameter-form-bounds.md),
and direct variation of its curvature-squared contact. The new calculation
keeps three generally nonzero, noncollinear conserved Euclidean momenta.
It does not replace a spatial multiplier by a constant.

The hardest step is cancellation of the two-derivative term after all
seagulls are included and restoration of the fourth-order term at the
same subtraction. The intermediate observable is the third variation
as a trilinear form. Before numerical execution, require:

* Exact symbolic agreement of the two-derivative cancellation,
  fourth-order contact, sixth-order coefficient and angular tail identity.
* Agreement of Feynman-parameter and direct two-dimensional loop
  quadratures within `3e-8+3e-7*abs(vertex)`, with successive direct-loop
  resolutions agreeing to that tolerance.
* Agreement with the independently differentiated constant-insertion
  bubble within `1e-10+1e-9*abs(vertex)`.

Numerical quadratures are diagnostics, not rigorous enclosures. The
uniform estimates below are ordinary inequalities for all real momenta.

## 1. The mass-exponential variation, including all contacts

Work on Euclidean R2 in dimensionless orbit-length units with Hermitian
Clifford matrices, `tr 1=2`, and

    D=gamma^a partial_a, S_M=(D+M)^-1,
    S_M(k)=(M-i gamma.k)/(k^2+M^2), M>0,
    Gamma_M[v]=-Tr log(D+M exp(v)).

The canonical conformal change of a two-dimensional Dirac determinant
also has its local Weyl allocation. Its flat Polyakov term is quadratic
in v and has no third variation. Vacuum potential terms and the common
four-dimensional curvature restoration are retained separately below.
For compact smooth directions u,v,w, differentiation gives

    Gamma_M'''[u,v,w]
      =-M Tr(S_M uvw)
       +M^2 Tr(S_M uv S_M w+S_M uw S_M v+S_M vw S_M u)
       -M^3 Tr(S_M u S_M v S_M w+S_M u S_M w S_M v).   (1)

The three second terms arise from the second variation of exp(v), and
the first from its third variation. A triangle alone is not (1).
These are regulated identities before subtracting local terms.

Let `p1+p2+p3=0`, `Qi=|pi|^2` and `K=max |pi|`. Removing only the
zero-momentum value of the kernel eliminates the local vacuum vertex.
Call the resulting kernel `V_M(p1,p2,p3)`. No physical vacuum term has
been set to zero: this defines the part to which the estimates apply.

## 2. Exact parameter integrals at arbitrary external momenta

The bubble difference follows by shifting its two-denominator integral:

    B_M(p)-B_M(0)
      =(1/(2pi)) integral_0^1 log(1+a(1-a)Q/M^2) da.   (2)

Indeed the spin numerator is `2(M^2-k.(k+p))`; after the Feynman
shift it is `2(Delta-l^2)`, `Delta=M^2+a(1-a)Q`.
The divergent constant cancels in (2).

For the triangle use propagator shifts `r1=0, r2=p1, r3=-p3`.
On the simplex `x,y,z>=0, x+y+z=1`, with measure `dx dy` and area 1/2,
define

    Dq=xy Q1+yz Q2+zx Q3,
    Jq=(1/2-z)Q1+(1/2-x)Q2+(1/2-y)Q3.                (3)

The trace numerator is `2M[M^2-k1.k2-k1.k3-k2.k3]`.
After shifting by `rbar=x r1+y r2+z r3`, put `ai=ri-rbar`.
The identity

    sum_(i<j) ai.aj + 3Dq = Jq

and the integrals

    integral d^2l/(2pi)^2 (l^2+Delta)^-3=1/(8pi Delta^2),
    integral d^2l/(2pi)^2 l^2(l^2+Delta)^-3=1/(8pi Delta)

give

    C_M(p1,p2,p3)
      =-M/(2pi) integral_simplex (2M^2+Jq)/(M^2+Dq)^2,
    C_M(0,0,0)=-1/(2pi M).                           (4)

Permutation of the simplex variables makes (4) symmetric in the three
Qi. Thus the two oriented triangles in (1) integrate to the same value.
Equations (1)-(4) give the exact subtracted vertex

    V_M = M^2/(2pi) {
       sum_i integral_0^1 log(1+a(1-a)Qi/M^2) da
       +2 integral_simplex [(2+j)/(1+d)^2-2] },
    d=Dq/M^2, j=Jq/M^2.                              (5)

This integral depends on all three momenta; it includes both seagull
and triangle effects. Its denominators are positive for real Euclidean
momenta, with no low-momentum expansion needed for (5).

## 3. Local cancellation and common fourth-order restoration

Scale all momenta by sqrt(t). Then d,j scale by t and

    (2+t j)/(1+t d)^2-2
      =t(j-4d)+t^2(6d^2-2jd)+t^3(-8d^3+3jd^2)+...

Use the exact simplex moment

    integral_simplex x^a y^b z^c
      =a! b! c!/(a+b+c+2)!.

The t coefficient of the braces in (5) is identically zero. This is a
seagull-dependent cancellation, not an assumption about the background.
The fourth-order kernel is

    V_M^(4)=(Q1 Q2+Q1 Q3+Q2 Q3)/(60pi M^2).          (6)

As an independent action check, the single-mode fourth-order quadratic
term of the retained bubble is `-Q^2/(120pi M^2)`. Its covariant local
action on a flat conformal chart is

    -1/(960pi M^2) integral sqrt(g) R^2
      =-1/(240pi M^2) integral exp(-2v)(Delta v)^2.

Its third variation is precisely (6). A sign or factor error in the
triangle would violate this check.

The sum of (6) over angular levels diverges logarithmically. It must
not be summed as an unsigned remainder. Define `R_M=V_M-V_M^(4)`
and restore the common finite curvature allocation separately:

    Gamma_a2[v]=alpha integral exp(-2v)(Rbar-2Delta v)^2,
    alpha=-C F_eff/(960pi), c2=8alpha,
    F_eff=4gamma_E-2H_q+4log(mu0 B).

At `Rbar=0` its third-variation kernel is
`-16alpha (Q1Q2+Q1Q3+Q2Q3)`. This is the same local functional
and finite scale as the existing parameter-form bound, not a new
per-level subtraction prescription. The zero-momentum vacuum
restoration also remains separate. Curvature and boundary completions
of these contacts are not inferred from their flat values.

For an additional low-frequency diagnostic the sixth-order kernel is

    -[2 sum_(i!=j) Qi^2 Qj +5Q1Q2Q3]/(840pi M^4).     (7)

Here `sum_(i!=j)` has six terms. Equation (7) is the leading part of
R_M; no definite sign is claimed for the full remainder at arbitrary K.

## 4. A global sixth-order bound with a summable angular tail

The positive log remainder obeys, for a>=0,

    f(a)=log(1+a)-a+a^2/2
        =a^3 integral_0^1 t^2/(1+ta) dt <= a^3/3.

The triangle remainder after its degree-two polynomial in d,j is exactly

    r(d,j)=d^2 [j(3+2d)-2d(4+3d)]/(1+d)^2.           (8)

Since d>=0,

    |r(d,j)|<=3|j|d^2+8d^3,
    d<=K^2/(3M^2), |j|<=3K^2/(2M^2).

The simplex has area 1/2. The bubble contribution is at most
`K^6/(280pi M^4)` and the triangle contribution at most
`43K^6/(108pi M^4)`. Consequently

    |R_M| <= a6 K^6/M^4, a6=3037/(7560pi).            (9)

No restriction `K<M` is imposed here. For integer magnetic charge q,

    d_n/M_n^4=(2C^2/q)[1/n^2-1/(n+q)^2],
    sum_(n>=1) d_n/M_n^4=(2C^2/q) H_q^(2)<4C^2/q,
    sum_(n>N) d_n/M_n^4
      =(2C^2/q) sum_(n=N+1)^(N+q) 1/n^2
      <=(2C^2/q) min(1/N, q/(N+1)^2), N>=1.         (10)

Thus the complete massive remainder `R=sum_n d_n R_Mn` converges
absolutely and locally uniformly in the external momenta. Equation (10)
gives a rigorous truncation error for any proposed numerical angular sum.
At the retained bounds `C/q<.00352`, `eta_G=GC/B^2<8.5e-33`,
the physical normalization `Rhat=(2G/B^2)R` obeys

    |Rhat|<3.3e-35 K^6.                              (11)

This is a flat kernel coefficient bound, not an error for an observable.

## 5. The full tower has fourth-order logarithmic growth

One need not lose six derivatives at high momentum. From (5),
`log(1+a)<=a` and

    |(2+j)/(1+d)^2-2|<=|j|+4d

give `|V_M|<=5K^2/(3pi)`. Equation (6) then yields, when K>=M,

    |R_M|<=a4 K^4/M^2, a4=103/(60pi).                (12)

The actual counting function satisfies

    N(E)=sum_(M_n<=E) d_n<=4C E^2,

because its finite sum through n is `2n(q+n+1)<=4n(q+n)`.
Stieltjes integration by parts therefore gives, for K>=M1,

    sum_(M_n<=K) d_n/M_n^2 <=4C+8C log(K/M1),
    sum_(M_n>K) d_n/M_n^4 <=8C/K^2.                  (13)

Use (12) in the first sum and (9) in the second. For K<M1 use (9)
and `sum d_n/M_n^4<=8C/M1^2`. Since M1>18 and pi>3, both cases give

    |Rhat| <= eta_G K^4 [6.8+9.2 log_+(K/18)]
             <20 eta_G K^4 log(2+K).                (14)

In particular the infinite angular tail does not turn this flat
three-point remainder into an arbitrary sixth-order operator.

There is a direct trilinear form consequence. Use Fourier transform
`fhat(p)=integral exp(-ip.x)f(x) dx` and measure `dp/(2pi)^2`. Set

    ||f||_E^2=integral (1+|p|^2)^2 log(2+|p|)
                             |fhat(p)|^2 dp/(2pi)^2,
    A0(f)=integral |fhat(p)| dp/(2pi)^2.

For three conserved momenta, the second-largest length is at least
K/2. Also `log(2+K/2)>=log(2+K)/2`. With
`Wi=(1+Qi)sqrt(log(2+sqrt(Qi)))`, it follows that

    K^4 log(2+K) <=6 sum_(i<j) Wi Wj.

Young's inequality and Cauchy--Schwarz bound each paired term by
`||u||_E ||v||_E A0(w)`, with the corresponding permutation. Moreover

    A0(w)<.36 ||w||_E,

since `integral (1+|p|^2)^-2 dp/(2pi)^2=1/(4pi)` and
`log 2>2/3`. Thus the normalized remainder has a continuous trilinear
extension from compact smooth directions with

    |Rhat[u,v,w]| <130 eta_G ||u||_E ||v||_E ||w||_E
                 <1.11e-30 ||u||_E ||v||_E ||w||_E.  (15)

This is a basepoint shape form on flat R2. E is an explicitly defined
inhomogeneous Fourier energy, not the exterior's combined Weyl-log
energy or its cold cosh norm.

## 6. Independent checks and limits of transfer

If p3=0 and p2=-p1, differentiating the constant-background two-point
kernel gives independently

    V_M(p,-p,0)=M^2/pi integral_0^1
       [log(1+x)-x/(1+x)] da, x=a(1-a)|p|^2/M^2.     (16)

The two-derivative term vanishes; the fourth coefficient is
`|p|^4/(60pi M^2)` and the sixth `-|p|^6/(210pi M^4)`,
consistent with (6)-(7). This special case verifies the full formula
without substituting it for a variable insertion.

The numerical checker also integrates the original rational spin
traces over the complete loop plane and compares with (5).
It uses neither a fitted contact nor the parameter integral to compute
the direct-loop value. Its run outcome and input hashes are retained
separately. The analytic bounds, infinite summation and form extension
are ordinary proofs, not consequences of sampled numerical tests.

The previous goal turn made progress by correcting the benchmark
dependency. This result advances the missing variable response by
providing an all-momentum, all-massive-level flat three-point form.
It is not a finite-background Lipschitz theorem: (15) is evaluated at
v=0, and higher insertions are needed to make it uniform along a path.

Remaining to-do list: derive the actual curved/direct-plus-reflected
three-point remainder with the common contacts and thermal winding;
match the principal logarithmic form to its Weyl partner in the
exterior norm; control finite-background and higher variations; solve
the coupled exterior and its signed source condition. S139.1-S139.2
and S074.2-S074.6 remain open. Opening, stability, transport, preparation
and creation are unproved, and the area budget remains unevaluated.
