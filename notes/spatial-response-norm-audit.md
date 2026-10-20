# The spatial-gradient endpoint and a sufficient flat response norm

2026-10-19, E003/S149 supporting S139/S074. Ordinary argument and
diagnostic protocol before execution.

The retained flat canonical linear response is not a bounded map
from the ordinary uniform C5 norm to the uniform source-gradient
norm. Its angular tower has a positive fourth-order logarithmic
tail. A logarithmically weighted Fourier norm controls both the
first and mixed scalar-gradient responses at the flat basepoint.

This is a regularity result for the exact flat comparison operator.
It does not identify the actual curved operator with that comparison,
supply S148's missing finite-background J, or disprove a physical
wormhole. In particular a bounded C5 source-value theorem is not
contradicted by failure of the stronger spatial-gradient endpoint.

## 1. Fixed convention and the complete linear multiplier

Retain S142's rank-two complex convention, the whole massive tower
and the common scale/contact:

    M_n^2=n(n+q)/C, d_n=2(q+2n), q=2000^4,
    eta_G=GC/B^2<8.5e-33, M_1>18, |F_eff|<100.

Let K be the flat action Hessian after the constant-vacuum and
two-derivative allocations, with the common fourth contact restored:

    K(k)=-C F_eff k^4/(120pi)
       +sum_n d_n M_n^2/(2pi) integral_0^1 R(a k^2/M_n^2) dz,
    a=z(1-z), R(v)=log(1+v)-v+v^2/2.                (1)

This is the multiplier already derived in the curved-form note,
section 7. No per-level unsigned inverse-square sum replaces its
common finite part. The normalized multiplier is
Khat=(2G/B^2)K. The Lorentzian scalar trace has the opposite action
sign; all estimates here concern absolute values.

At a flat basepoint the vacuum scalar is a spatial constant and
has zero scalar variation. Its coordinate-density variations cancel
when converting the action derivatives to scalar derivatives. No
physical vacuum value is assigned zero.

The fermionic minus trace-log and bubble conventions agree with
[Carmi--Di Pietro--Komatsu, section 6.1](https://arxiv.org/html/1810.04185#S6.SS1),
re-inspected on 2026-10-19. Their interacting model is not used
as a theorem about this exterior. The proof below uses elementary
positive parameter integrals and exact angular sums; it does not
need to assume a high-frequency asymptotic or a curved symbol theorem.

## 2. Upper bound with the entire tower

For v>=0,

    0<=R(v)<=min(v^2/2,v^3/3),
    R(v)=v^3 integral_0^1 t^2/(1+tv)dt.

Since integral a^2=1/30 and integral a^3=1/140, a single
remainder term L_M(k) obeys

    0<=L_M(k)<=min[k^4/(120pi M^2),k^6/(840pi M^4)].

The exact counting function N(E)<=4CE^2 gives

    sum_(M<=k) d/M^2 <=4C+8C log(k/M1),
    sum_(M>k) d/M^4 <=8C/k^2,                       (2)

for k>=M1. Split there. The normalized remainder is at most

    (2eta_G/pi) k^4 [3/70+log(k/M1)/15].

For k<M1 use sum d/M^4<=8C/M1^2 and the sixth-order bound.
Together with the common contact, pi>3 and M1>18, both cases give

    |Khat(k)|<=eta_G k^4[.6+.05 log_+(k/18)]
              <=eta_G k^4 log(2+k).                (3)

The contact alone costs less than (5/9)eta_G k^4.
The second inequality follows from log 2>.6 for k<=18,
and log k>=log 18>2 for k>=18. These are continuum estimates;
no angular truncation or finite grid establishes them.

## 3. Lower bound and the failed C5 gradient endpoint

On z in [1/4,3/4], a>=3/16. If k>=4M then
v=a k^2/M^2>=3, so R(v)>=v^2/2-v>=v^2/6.
Thus

    L_M(k)>=3 k^4/(2048pi M^2), k>=4M.             (4)

For any integer N>=1, put k_N=4M_N. Exact partial fractions give

    d_n/M_n^2=2C[1/n+1/(n+q)],
    sum_(n<=N) d_n/M_n^2
       =2C[H_N+H_(N+q)-H_q]>=2C H_N.

All unretained remainder terms are positive. Therefore

    Khat(k_N)/(eta_G k_N^4)
       >=-|F_eff|/(60pi)+(3/(512pi))H_N
       >=-|F_eff|/(60pi)+(3/(512pi))log(N+1).       (5)

The right side grows without bound at fixed positive eta_G.
A finite local fourth-order counterterm changes the constant in
(5), not this growth.

Consider the smooth bounded functions on the flat plane

    u_N(x1,x2)=k_N^-5 cos(k_N x1).

For k_N>=1, max_(|beta|<=5)||partial^beta u_N||inf=1.
As a Fourier multiplier on periodic bounded functions,

    ||grad Khat(D)u_N||inf=|Khat(k_N)|/k_N^4,

which is unbounded by (5). The use of bounded periodic fields
avoids a false compact-support claim. This is an exact flat
comparison, not a direction asserted to satisfy the curved paired
domain. No localization remainder or equivalence of curved norms
is inferred here.

By contrast (3) gives
||Khat(D)u_N||inf<=eta_G log(2+k_N)/k_N, tending to zero.
This distinguishes the source-value endpoint from its gradient.
For any exponent 5+delta, delta>0, the single-wave gradient is
bounded, but a single-wave test alone is not a general function-space
mapping theorem.

The conclusion concerns the retained operator at arbitrarily high
frequency. It is not a physical ultraviolet impossibility statement.
An explicitly specified cutoff would replace the divergent factor
by a cutoff-dependent logarithm; introducing one would be an
additional model decision, not a proof in the present norm.

## 4. A sufficient norm for first and mixed scalar gradients

Use Fourier measure dp/(2pi)^2. For integrable Fourier transforms set

    A0(u)=integral |uhat(p)| dp/(2pi)^2,
    A5log(u)=integral (1+|p|)^5 log(2+|p|)|uhat(p)| dp/(2pi)^2.

The identical inequalities hold for discrete Fourier sums.
Equation (3) gives the actual linear mapping estimate

    ||grad Khat(D)u||inf<=eta_G A5log(u).            (6)

This Fourier norm is a sufficient norm, not another name for
ordinary C5. Also log(2+r)<=1+r, so A5log<=A6, where
A6 uses weight (1+|p|)^6. This last statement concerns Fourier
absolute integrability, not an unproved C6-to-A6 embedding.

The mixed source must include the scalar-area conversion. Let
Vhat(p,q,r) be S142's complete third action variation with
p+q+r=0, vacuum removed and the common fourth contact restored.
The scalar mixed kernel at the flat basepoint is

    Bhat(q,r)=Vhat(-q-r,q,r)-2Khat(|q|)-2Khat(|r|).  (7)

It follows by differentiating exp(-2b) delta Gamma/delta b.
The zero-momentum vacuum terms cancel exactly. Omitting the two
terms in (7) would give a coordinate-density response instead.

Put Kmax=max(|q+r|,|q|,|r|). S142 bounds the third-action remainder
by 20 eta_G Kmax^4 log(2+Kmax). Its restored contact costs
at most (10/3)eta_G Kmax^4, hence less than
5 eta_G Kmax^4 log(2+Kmax). Equation (3) consequently gives

    |Bhat(q,r)|<30 eta_G Kmax^4 log(2+Kmax).        (8)

Write m=max(|q|,|r|). Then |q+r|<=2m, Kmax<=2m and
log(2+2m)<=2log(2+m). Multiplying (8) by the output
gradient costs at most 1920 eta_G m^5 log(2+m).
Fubini and the elementary maximum bound prove

    ||grad D^2 S_flat[ u,v ]||inf
       <2000 eta_G [A5log(u) A0(v)+A0(u) A5log(v)]. (9)

The action/trace sign does not change (9). This is a full
bilinear Fourier mapping estimate at the flat basepoint, including
both seagulls, both triangle orientations, the common contact and
the density terms. It is not a finite-background curved theorem.

There is a useful independent constant-direction check. If the
second direction is spatially constant, covariance of the flat
scalar under a constant scale gives

    Bhat(0,p)=-k dKhat(k)/dk.

For the action vertex the corresponding formula is
Vhat(0,p,-p)=2Khat(k)-kKhat'(k). Both formulas retain the
fourth contact, and differ by the required density term.

## 5. Consequence for S148 and the cold problem

S148 correctly leaves the moving-cut Hessian conditional on a
spatial trace bound. The ordinary C5 source-value estimate alone
does not supply the missing spatial-response norm: even the exact
flat comparison has the endpoint failure (5). Equations (6),(9)
provide a sufficient flat norm in which to seek the actual
finite-background estimates. They cannot be inserted as the
actual curved J without that proof and a norm-transfer argument.

The prior centered-canonical-response note already proves a
cold oscillation estimate for its particular endpoint contract,
proportional to epsilon h L_h sqrt(chi). The complete-residual
note already carries it through the lapse integral. Those results
are retained, not counted as a new cold theorem here.
They concern the specified endpoint comparison and do not supply
an arbitrary finite-background mixed response on S145's raw C5 ball.

The distinction matters quantitatively: S148's H_bolt scales as
1/h. A spatial regularity repair alone does not preserve the h
weight that cancels this amplification. The next coupled estimate
must retain both the logarithmic spatial control and the physical
cold weight, or use a fixed-interface formulation that proves the
needed integrated response directly. No new radius, temperature,
counterterm or ultraviolet cutoff is selected by this audit.

## 6. Verification protocol

Before execution freeze exact algebraic checks of the mode and
counting bounds, common contact, scalar/density conversion and
constant-direction identity. Check the high-frequency sequence
by its proven harmonic lower bound, not a finite cutoff extrapolation.

Independent parameter quadratures use 48 and 96 Gauss nodes,
with a 1e-10 relative-to-(1+|value|) comparison to adaptive
quadrature. A manufactured q=7,C=2 tower diagnostic sums 512
and 1024 levels and retains the exact inverse-fourth tail bound;
these are not replacement physical parameters. Two finite
differences at steps 1e-4 and 5e-5 must agree with the constant-
direction identity to 2e-7 relative-to-(1+|target|), with error
ratio below .35 unless both errors are below 1e-35. Fourier
mode-pair diagnostics must detect omission of the density terms.
All failed acceptance runs retain their source snapshots.

Remaining to-do list: prove the curved finite-background spatial
mapping and compatible norm transfer; retain cold weighted response;
complete other source maps and the global exterior/state solve;
charge and signed source; opening, stability, transport, preparation
and creation.
