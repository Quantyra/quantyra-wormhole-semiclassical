# A quantitative finite mass integral on the stationary RN reference

2026-10-15, E003/S139. Ordinary fixed-reference argument before focused
verification. No numerical run is reported here. This note reuses the
frozen complete source bound and derives a quantitative absolutely
convergent vacuum-subtracted exterior mass integral. It does not change
any previous artifact, the subtraction scheme or the prescribed state.

For the specified classical RN/HH reference, one exterior of the identical
pair obeys

```
rho2,ren = rho2,mass + rho2,LLL + rho2,W - rho2,AdS,
rho2,AdS(r)=11 r^2/(240pi L4^4),
M_abs = integral_(r_h)^infinity |rho2,ren(r)| dr,
b M_abs < 1e16 q^2.                                    (1)
```

The signed exterior source integral DeltaM_ext=integral rho2,ren dr
therefore exists and has absolute value at most M_abs. At k=2000,

```
G M_abs/(b C) < 1e16/(6 k^8) < 6.511e-12,
G M_abs/b     < 1e16/k^5 = 5/16.                       (2)
```

These are bounds on the complete fermion bulk source evaluated on the
same classical geometry. In the first Einstein-constraint response,
delta mu(infinity)=delta mu(r_h)+G DeltaM_ext, where mu is the
vacuum-subtracted geometric mass parameter defined below. The horizon
integration datum is not fixed by (1). Neither
this signed geometric mass coefficient nor the source integral is
identified with the Hamiltonian of the full renormalized action.

## Fixed inputs and allocations

Use precisely the field, monopole bundle, identical paired projectors,
HH state, laboratory clock and common four-dimensional prescription of
[the complete reference source theorem](stationary-reference-higher-mode-clock.md).
The common scale is mu0 b=1 and all additional finite local coefficients
remain zero. The split is the retained bulk Weyl functional W4, the
canonical product-frame massless LLL functional, and its complementary
common-renormalized higher angular functional. No angular determinant is
renormalized separately.

The symbols are

```
t=r/b, v=r_h/b=1+h, u=1/t,
xi=b^2/L4^2, C=1+6xi, ell=b/sqrt(C),
q=k^4, k=2000, eta=1/2000, z_i in [1/4,7/25],
xi(1+3xi)=pi k^6, b^2/G=xi k^10,
F(t)=1+xi t^2-2a_m/t+Qhat/t^2,
a_m=m/b, Qhat=Q/b^2=1+3xi,
B(u)=u^2 F(1/u), or B(t)=F(t)/t^2,
Theta(t)=b^2 tau_hat,mass(bt),
Y=y/ell=sqrt(C) integral_0^u B(w)^(-1/2) dw.
```

The notation B(t) means the same scalar composed with u=1/t; a prime
on B below means the explicitly displayed argument derivative. Frozen
parameter bounds include

```
xi>=k^3>=1000, C/q<r0=.00352, 0<h<1e-30,
a_m<2.01xi, Qhat<3.01xi,
a_m<C<q, Qhat<C<q, xi<C<q,
M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18.
```

The d_n count both nonzero angular signs; the orbit heat kernel has
rank two. These are proper Euclidean angular masses, not Killing gaps.
The already proved complete reference bounds are

```
|Theta'(t)|<=K q^2, K=1e15,
|Theta(2)|<=K0 q^2, K0=1e11,
|rho2,mass|, |p2,mass| <= S q^2/b^2,
S=K+K0/2=1.00005e15.                                 (3)
```

The HH integration constant is fixed by smoothness at the actual
classical bolt. All constants here are uniform over the declared
preparation interval at the fixed k. No response estimate for an
unknown corrected metric is imported.

The vacuum term in (1) is the already established asymptotic AdS Dirac
vacuum, in orbit normalization T2=4pi r^2 T4. Subtracting it defines
the vacuum-subtracted geometric mass on this reference. It is an
explicit relative source quantity; it is not a new finite adjustment
of one angular sector or a solution of the shifted outer-curvature
equation. If the vacuum is instead retained in the geometric equation,
its cosmological response must be restored separately.

## An explicit angular discrepancy bound

Let tau denote heat proper time in ell units and define

```
Atilde(tau)=sum_(n>=1) d_n M_n^2 exp(-M_n^2 tau),
E(x)=x(x+q)/C,
g(x)=E'(x)E(x)exp[-tau E(x)].
```

Then Atilde=2C sum g(n). On using E as the variable,

```
g=E sqrt(q^2+4CE) exp(-tau E)/C,
d_E log g=1/E+2C/(q^2+4CE)-tau.
```

The last expression strictly decreases from positive infinity to -tau.
Thus g is unimodal and vanishes at both endpoints. The rectangle-rule
error is bounded directly, without an asymptotic expansion:

```
|sum_(n>=1)g(n)-integral_0^infinity g(x)dx|
 <= integral_0^infinity |g'(x)|dx =2 max g.
```

For example apply the fundamental theorem of calculus on each interval
[n-1,n], integrate its absolute derivative bound and sum. The continuous
integral is 1/tau^2. Since sqrt(q^2+4CE)<=q+2sqrt(CE),
max E exp(-tau E)<=tau^-1 and
max E^(3/2)exp(-tau E)<tau^-3/2, it follows that

```
|Atilde(tau)-2C/tau^2|
 <=4q/tau+8sqrt(C)/tau^(3/2),  tau>0.                  (4)
```

This weaker-than-optimal UV estimate suffices because it gains a half
power relative to the universal leading term. Its q/tau part is used
below only on 0<tau<=.01. Integrating that discrepancy bound to infinite
heat time would not be valid for the estimate used here.

## Subtracting the full universal cubic trace term

The frozen [non-bolt proof](stationary-reference-nonbolt-stress.md)
uses the complete even RN cap, a smooth analytic continuation at its
outer seam, and actual-domain resolvent localization. Its rank-two
cubic heat profile K3 in ell units has coefficient

```
c=c3 ell^3=-a_m sqrt(xi)/C^(3/2), |c|<1,
|partial_Y K3(tau,Y)|<=exp[-Y^2/(2tau)].
```

The exactly integrated universal angular profile is

```
Theta_sing(Y)
 = C integral_0^infinity (2C/tau^2)K3(tau,Y)dtau
 = c C^2/(5pi Y).                                    (5)
```

This is the exact cubic coefficient proved in
[the outer boundary-layer calculation](stationary-massive-boundary-layer.md).
Both the profile and its derivative have absolutely convergent model
integrals at every Y>0. The comparison below uses the full angular
Atilde, not a replacement of its higher terms by zero.

To evaluate the derivative after (5), split heat time at tau0=.01.
On the short interval insert

```
actual - even cap,
even cap - smooth continuation - K3,
smooth continuation - common interior coefficients,
(Atilde-2C/tau^2) K3.
```

For long time keep the even cap, the actual-to-cap resolvent difference
and the common polynomial restorations. Subtract the long-time part of
the universal model (5) explicitly. This prevents any use of (4) at
infinite heat time. The common constant a0 finite restoration has zero
spatial derivative; its value is already retained in the anchor (3).
The common k1 cancellation and finite k2 term are unchanged.

Here is a ledger directly from the frozen non-bolt inequalities, in
units q^2 for Theta_Y. It is uniform for 0<Y<=1:

| Term | Upper bound |
| --- | ---: |
| Short smooth subtracted kernel, Cs=1e15 | 8Cs(C/q)^3 <3.49e8 |
| Long even cap, R_even<=1e10 | 1400 R_even C/q <4.928e10 |
| Common local long subtractions and finite k2 restoration | less than 1 |
| Actual positive-square resolvent versus the cap | less than 1 |
| Removed universal-model tail tau>=.01 | 200(C/q)^2 <1 |

The sum is less than 5e10. For the local-restoration row, one may use
the explicit upper bound
r0(3/8+5/600)+500r0^2<1 after dropping the favorable 4pi.
For the actual-cap row the frozen bound is at most .2(C/q)<1.
The long model row follows from |K3_Y|<=1 and
integral_.01^infinity 2C/tau^2 dtau=200C.
Thermal windings are already included in the frozen Cs, Cb and cap
estimates; the thermal state has not been replaced by a vacuum kernel.

The reflected-kernel remainder obeys, in proper gradient units,

```
|partial_y(k_even-k_smooth-K3)|
 <= Cb ell^-3 sqrt(tau)exp[-Y^2/(4tau)], Cb=1e16.
```

Together with Atilde<=24C/tau^2 and
integral_0^infinity tau^-3/2 exp[-Y^2/(4tau)]dtau<4/Y,
its contribution to |Theta_Y| is at most

```
96 Cb C^2/Y.                                         (6)
```

The short angular discrepancy (4) contributes at most

```
C [4q E1(50Y^2)+24sqrt(C)/Y]
 <=C(12q+24sqrt(C))/Y < q^2/Y.                        (7)
```

Here E1(z)=integral_z^infinity e^-s ds/s. If z>=1 it is below one;
if z<1 split at one to obtain
E1(50Y^2)<=1+2log(1/Y)<=3/Y for 0<Y<=1.
The other integral is
integral_0^infinity tau^-3/2 e^[-Y^2/(2tau)]dtau
=sqrt(2pi)/Y<3/Y. For the last strict bound in (7),
12C/q<.04224 and 24C^(3/2)/q^2<24/sqrt(q)=6e-6.

Since 96Cb r0^2=1.1894784e13, the ledger, (6) and (7) prove

```
|partial_Y[Theta-Theta_sing]|
 <= D q^2/Y, D=2e13, 0<Y<=1.                         (8)
```

In particular the logarithmic residual now has an explicit uniform
constant. No qualitative big-O has been differentiated.

## From the proper-distance profile to a trace tail in r

The frozen geometry gives

```
B<=C/5, B>=C(1-u)^2/8,
sqrt(5)u<=Y<=sqrt(8)log[1/(1-u)].                    (9)
```

These inequalities hold on the entire outer region used here. At
u=1/2 the lower bound exceeds one, while at u=1/4 the upper bound
is less than sqrt(8)/3<1. Therefore the unique Y=1 point has
t1 in (2,4). Monotonicity of Y supplies the uniqueness. The frozen
anchor and global gradient give

```
|Theta(t1)| <= (K0+2K)q^2.
```

Also |Theta_sing(1)|<=C^2/(5pi)<q^2. Integrating (8) from Y=1
therefore yields

```
|Theta(Y)-Theta_sing(Y)|
 <=q^2[3e15+2e13 log(1/Y)], 0<Y<=1.                  (10)
```

The relation between 1/Y and the exact linear-r trace also needs an
estimate. Put alpha=sqrt(C/xi). When Y<=1, (9) implies
u<=1/sqrt(5). For 0<=w<=u the exact quartic gives

```
|B(w)/xi-1|
 <=w^2[.001+4.02w+3.01w^2]
 <(5/2)w^2<=1/2.
```

The reciprocal square root on this interval obeys
|(B/xi)^-1/2-1|<=2|B/xi-1|. Hence

```
|Y-alpha u| <=(5alpha/3)u^3,
|1/Y-1/(alpha u)| <=sqrt(5)u/3<1.                   (11)
```

The second step uses Y>=sqrt(5)u, so no endpoint division is
uncontrolled. Define

```
L=-a_m xi/(5pi).
```

Then c C^2/(5pi alpha u)=Lt exactly. Equations (5),(11) show
|Theta_sing-Lt|<q^2. Moreover Y>=sqrt(5)/t implies
log(1/Y)<=log t. Since Y(1/4)<1, all t>=4 lie in this region.
Thus the explicit trace remainder

```
R(t)=Theta(t)-Lt
```

satisfies

```
|R(t)|<=T q^2(1+log t), T=4e15, t>=4.                (12)
```

The sign and coefficient L are fixed by the actual balanced domain;
they have not been selected to cancel another sector after the fact.

## Exact conserved mass source and its finite pressure integral

Use the frozen HH reconstruction, writing a dimensionless integral

```
I(t)=integral_v^t B'(s)Theta(s)ds,
b^2 p2,mass=I/(2F),
b^2 rho2,mass=I/(2F)-Theta/t^2.                       (13)
```

Here B(t)=F(t)/t^2 and

```
B'(t)=-2/t^3+6a_m/t^4-4Qhat/t^5.
```

HH regularity fixes I(v)=0. On [v,4], B'>=0. To verify this,
t^5B'/2=-t^2+3a_m t-2Qhat is concave; at v it is positive
because F'(v)>0, and at four it is at least -6+18xi>0.
Furthermore

```
B(4)=xi+1/16-a_m/32+Qhat/256
 <=xi+(9-13xi)/256 <xi.
```

The anchor and gradient imply |Theta|<=(K0+2K)q^2 on this compact
interval. Consequently

```
|I(4)|<=xi(K0+2K)q^2.                               (14)
```

For t>=4, elementary reference inequalities give

```
|B'(t)| <= [2+(3/2)a_m+Qhat/4]/t^3
         <4xi/t^3,
F(t)>=xi(t-v)(t-1)>=xi t^2/2.                        (15)
```

The first uses xi>=1000, a_m<2.01xi and Qhat<3.01xi. The second
uses v<1.001; already (1-1.001/4)(1-1/4)>1/2. The underlying
exact horizon factorization and monotonicity were established in
the reference clock theorem.

No derivative of R in (12) is needed. Reconstruct directly as

```
I(t)=I(4)+L integral_4^t s B'(s)ds
                +integral_4^t B'(s)R(s)ds.           (16)
```

For an independent exact check the first integral has primitive
2/t-3a_m/t^2+4Qhat/(3t^3). Bounds (12),(15) give

```
|I(t)| <= xi q^2[K0+2K+1/15+(3/8)T]
        <4e15 xi q^2, t>=4.                          (17)
```

Indeed |L|<=q^2/15 since a_m,xi<=q and pi>3; its contribution
is at most xi|L|. The remainder contribution uses the exact integral

```
4 integral_4^infinity (1+log s)/s^3 ds
 =(3+2log4)/16 <3/8.
```

Thus I has an absolutely convergent outer limit and the pressure
term has the explicit integrable bound

```
|I/(2F)| <=4e15 q^2/t^2,
integral_4^infinity |I/(2F)|dt <=1e15 q^2.            (18)
```

The frozen individual-source bound (3) handles the lower endpoint
without estimating the separate singular-looking 1/F factors:

```
integral_v^4 |b^2 rho2,mass|dt <=3S q^2.              (19)
```

In particular no inverse temperature has entered the mass norm.

## The exact Weyl cancellation inside the absolute integral

The RN compact-variation source is fixed in
[the radial Weyl-source note](stationary-weyl-radial-sources.md). In
the present variables it can be written exactly as

```
b^2 rho2,W(t)=V(t)+L/t+RW(t),
V(t)=11xi^2t^2/(240pi).
```

The remainder, not an asymptotic symbol, is

```
RW=-(Alog log t+Aconst)/(720pi t^6),

Alog=432Qhat^2+288Qhat xi t^4-864Qhat a_m t
       +432Qhat t^2-60q^2t^4,

Aconst=-593Qhat^2-338Qhat xi t^4+1408Qhat a_m t
       -504Qhat t^2-588a_m^2t^2+288a_m t^3-11t^4.
```

The omitted two terms of the original Euler numerator are precisely
-33xi^2t^8 and +144xi a_m t^5; they give V and L/t respectively.
This checks the sign against (5) rather than relying on a proposed
subtraction of a logarithmic mass divergence.

Since a_m,Qhat,xi<=q, q>=1 and t>=1, the absolute coefficients
of Alog are bounded by 2076q^2t^4, and those of Aconst by
3730q^2t^4. Since 720pi>2160,

```
|RW(t)|<=2q^2(1+log t)/t^2,
integral_v^infinity |RW(t)|dt<=4q^2.                  (20)
```

Now combine the sectors before taking the absolute value. By (12)-(13),

```
b^2[rho2,mass+rho2,W-rho2,AdS]
 = I/(2F)-R(t)/t^2+RW(t), t>=4.                      (21)
```

The separate +/-L/t terms have canceled exactly. Neither separate
mass integral is being asserted absolutely convergent. Equations
(12),(18),(20) bound the three remaining outer terms, because

```
integral_4^infinity |R(t)|/t^2 dt
 <=Tq^2(2+log4)/4 <4e15q^2.                          (22)
```

On [v,4] use the regular mass bound (19), the finite term
integral_v^4 |L|dt/t<q^2, and the applicable part of (20).
The same single RW norm (20) can cover both intervals without
double counting.

## Same-reference canonical LLL contribution

The LLL piece here is the exact canonical Polyakov source on this
same classical geometry, not the lapse or mass of the separately
solved nonlinear LLL geometry. Its trace is

```
Theta_LLL(t)=q[-2+12a_m/t-12Qhat/t^2]/(24pi),
|Theta_LLL|<=qC, |Theta_LLL'|<=qC/4.
```

The second inequality is the already retained canonical derivative
bound. The first follows immediately from a_m,Qhat<C and pi>3.
Use the same HH stress reconstruction as (13). The global
anchor-plus-gradient estimate gives |rho2,LLL|<=qC/b^2; hence
the interior integral is at most 3qC. At t=4, positivity of B'
and B(4)<xi give |I_LLL(4)|<=xi qC. Equations (15) then give

```
|I_LLL(t)| <=(9/8)xi qC,
|b^2 rho2,LLL(t)| <=(17/8)qC/t^2 <3qC/t^2,
integral_v^infinity |b^2 rho2,LLL(t)|dt <4qC<q^2.     (23)
```

The factor 9/8 uses integral_4^infinity 4s^-3ds=1/8.
This is an actual finite same-reference contribution in the declared
allocation, with no new Casimir length or reservoir state.

## Absolute mass ledger and its meaning

Use (21) only on the outer interval, the regular source bound on
the inner interval, and (20),(23) globally. The absolute integrated
source is bounded, in units q^2, by

| Contribution | Dimensionless integral bound |
| --- | ---: |
| Massive source on [v,4] | 3S |
| Weyl linear term on [v,4] | less than 1 |
| Outer pressure I/(2F) | 1e15 |
| Outer trace remainder R/t^2 | 4e15 |
| Entire Weyl remainder RW | 4 |
| Entire canonical LLL source | less than 1 |

Their sum is less than 9e15, so the larger 1e16 in (1) is a strict
bound. Here the cancellation has occurred inside the combined source
before its absolute integration. It is stronger than a bound obtained
by assigning two divergent signed integrals opposite values.

Finally b M_abs=integral_v^infinity |b^2 rho2,ren|dt, and

```
Gq^2/(b^2 C)=1/(xi C k^2)<=1/(6k^8),
Gq^2/b^2=1/(xi k^2)<=1/k^5.
```

This proves (2) with the retained exact parameter family. The
finite source integral has now been bounded quantitatively, rather
than inferred from the qualitative outer O(log r)/r^2 theorem.

The signed geometric coefficient still contains a horizon datum. To
track the vacuum allocation explicitly put A_v=11/(240pi L4^4), so
rho2,AdS=A_v r^2, and write

```
delta mu(r)=delta m_geo(r)-G A_v r^3/3,
delta mu_h=delta m_geo,h-G A_v r_h^3/3,
delta mu(r)=delta mu_h+G integral_(r_h)^r rho2,ren(R)dR,
delta mu_infinity=delta mu_h+G DeltaM_ext.
```

Here delta m_geo is the geometric mass response before removing the
known vacuum volume term. The prepared-temperature condition determines
the horizon datum only after
horizon recentering and the required metric/state response are
controlled. A small ratio in (2) to C is not smallness relative to
the thermal gap h, which is of order k^-9 here. In particular it
does not justify expanding a fixed classical horizon through a
source displacement larger than that gap. No nonlinear continuation,
uniqueness or physical state transfer follows from this first-order
reference integral.

The full-action Hamiltonian has additional variational and boundary
obligations identified by [the charge audit](stationary-exterior-charge-audit.md).
This integral does not discharge them. The higher-order solver
obstruction also remains unchanged: fixed-reference source amplitudes
are not renormalized metric-response norms.

Remaining to-do list: control horizon recentering and the complete
horizon-weighted metric response; justify the selected higher-order
or controlled reduced solver; match the corrected exterior/state and
full boundary charge; transfer the driven reference results before
closing the physical S074-S079 gates.
