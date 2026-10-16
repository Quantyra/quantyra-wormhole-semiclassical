# The massive fixed-RN trace through the cold bolt

2026-10-15, S139/S074. Ordinary proof independently reviewed before
focused verification.
For the retained classical RN Hartle--Hawking reference, including the
complete nonzero angular tower in the common four-dimensional scheme,
the canonical product-frame trace satisfies

```
|d[b^2 tau_hat,mass]/d(r/b)| < 3e14 q^2,
             r_+ <= r <= 1.001 b.                         (1)
```

The bound includes the cold bolt and the matching overlap. It is a
fixed-reference result. It is not a Frechet response estimate, a
nonlinear corrected-exterior theorem, a boundary-charge calculation,
or a bound on driven stress. The selected finite coefficients, source,
preparation rectangle and scheme are unchanged. In particular the
[higher-order solver obstruction](stationary-exterior-solver-obstruction.md)
is not evaded by (1). No computation or new receipt is reported here.

## Reference, normalization and the exact radial polynomial

Use the identical paired magnetic exteriors and actual projector of
[the conformal-domain note](stationary-exterior-conformal-domain.md),
with the global bundle pairing in
[the phase audit](stationary-paired-dirac-phase.md). The fixed sphere
in the conformal product is the classical extremal radius b. Write

```
k=2000, q=k^4, eta=1/2000,
C=sqrt(1+12 pi k^6), xi=(C-1)/6, b^2/G=xi k^10,
nu=2q eta G C/b^3, kappa=z_i nu, 1/4<=z_i<=7/25,
ell=b/sqrt(C), M_n^2=n(n+q)/C,
d_n=2(q+2n), n>=1.                                       (2)
```

The state is the regular static HH restriction on each exterior, with
the canonical paired purification. Euclidean angle phi=kappa U_E has
period 2pi. Distances and curvatures in the proof are in ell units.
The proper angular gap is M_n>=18 sqrt(n); this is not an optical
Hamiltonian gap. Indeed C<(123/20)k^3 and k=2000 give M_1>18.
We also use C/q<.00352 and xi>=1000.

Set s=r/b-1 and h=(r_+-b)/b. With

```
a=4xi/C<=2/3, c=xi/C<=1/6,
delta=h^2(1+ah+ch^2)/(1+h),
P(s)=s^2+a s^3+c s^4-delta(1+s),
q0(s)=1+a[s+h^2/(s+h)]+c(s^2+h^2)-delta/(s+h),
P=(s-h)(s+h)q0(s),
lambda=h q0(h)/(1+h)^2=kappa b/C.                         (3)
```

Here q0 is a polynomial-factor function, not the magnetic charge q.
The orbit metric divided by ell^2 is

```
d rho^2+f(rho)^2 dphi^2,
ds/d rho=sqrt(P), f=sqrt(P)/[lambda(1+s)^2].               (4)
```

The preparation gives eta/(8k^9)<=h<=eta/k^9<1e-30 and
log(1/h)<80. To see the lower bound without an extremal limit,
lambda=2z_i eta/(xi k^6), k^3<xi<2k^3, and
q0(h)/(1+h)^2 is between .998 and 1.001 on the accepted small-h
branch. These also give the displayed upper bound. All estimates are
at this fixed positive temperature; no uniform zero-temperature
determinant or new extremal state is asserted.

On h<=s<=s1=.01, direct differentiation of (3) gives

```
1<=q0<=1.01004, 0<=q0'<.921,
0<=q0(h)-1<=.501h.                                       (5)
```

For example q0'=a[1-h^2/(s+h)^2]+2cs+delta/(s+h)^2.
At h, positivity follows by multiplying (q0(h)-1)/h by 1+h:
the result is (3a/2-1/2)+(a+2c)h+(3c/2)h^2. The same
formula and a<=2/3 give the upper bound.

## An exact small H2 cap, with the thermal scaling retained

Let v=log tanh(chi/2), where chi is the radial coordinate of unit H2.
Choose its free scale by v(s1)=-atanh(h/s1), and solve

```
v_s=lambda(1+s)^2/P,
J(s)=q0(h)(1+s)^2/[q0(s)(1+h)^2],
v=-atanh(h/s)-h D(s),
D(s)=integral_s^s1 (J(t)-1)/[(t-h)(t+h)] dt.               (6)
```

The factor J-1 vanishes at h; its apparent pole there is removable.
Since 0< (log J)'=2/(1+s)-q0'/q0 <=2, one has

```
0<=hD<=2.03h log[(s1+h)/(s+h)].                           (7)
```

The actual metric is exactly exp(2sigma)g_H2 on this patch. A formula
regular even at the bolt is

```
A(s)=sqrt(q0(s))(1+h)^2/[q0(h)(1+s)^2],
exp(sigma)=A(s)[cosh(hD)+(s/h)sinh(hD)].                   (8)
```

This follows from exp(sigma)=f/sinh(chi) and
sinh(chi)=1/sinh[atanh(h/s)+hD]. Thus it does not compare
different temperatures or discard the long throat. Equations (5)-(8)
give |sigma|<.038 on the complete patch. Conservative derivative bounds
are

```
|sigma_s|<=20L, |sigma_ss|<=100/(s+h),
L=1+log[(s1+h)/(s+h)].                                   (9)
```

One way to check (9) is to use the rational derivative
D'=-(J-1)/[(s-h)(s+h)], then differentiate (8); no derivative of a
big-O remainder is used. For s0<=.001, on the complex half disk
|s-s0|<=(s0+h)/2, q0 stays within .004 of 1,

```
|D|<=2.03 log[(s1+h)/(s0+h)]+4.5,
|D'|<=4.5/|s+h|, |D''|<=20/|s+h|^2.                     (10)
```

These follow by integrating from the real s0 and then along a straight
segment; q0'/q0 and its derivative are rational with their nearest
remaining pole at -h. The hD terms in (8) are bounded directly, using
h<1e-30. On the quarter disk, |q0-1|<=2.113s0,
|log q0|<=2.118s0, and

```
2|log(1+s)-log(1+h)|
 <=(2.5s0-1.5h)/(1-1.5s0).
```

Consequently |log A|<.006 and |sigma|<.020 for s0<=.001.
These complex guards will be used for actual analytic chart estimates.

Extend this metric to a complete small conformal cap on H2. At chi1,
where s=s1, (6) gives D=0 and direct substitution gives
|sigma1|<.02, |sigma1'|<.05, |sigma1''|<.1. For
t=chi-chi1 in [0,1] set

```
sigma_cap=eta5(t)[sigma1+sigma1' t+sigma1'' t^2/2],
eta5=1-10t^3+15t^4-6t^5; sigma_cap=0 for t>=1.            (11)
```

Primes here mean chi derivatives. The inequalities
eta5 t<=.3, eta5 t^2<=.14 show |sigma_cap|<.042<.05.
It matches through two derivatives and has |sigma_cap'|<.4,
|sigma_cap''|<1.5. C2 suffices for this massive resolvent construction;
the whole retained inner patch remains analytic.
For the seam bounds one may use |D'(s1)|<135,
|D''(s1)|<1.4e4, |(log A)'|<2 and |(log A)''|<3.
Equation (8) then gives |sigma_s|<3.35, |sigma_ss|<416;
P(s1)<.00010067 and P'(s1)<.020202 yield the stated chi jets.

For a unit hyperbolic transvection, chi_t=cos(phi) and
chi_tt=coth(chi)sin(phi)^2. Combining these exact formulas with (9)
removes the apparent bolt singularity. In particular
sqrt(P)coth(chi)<=1.01s, and the source variations obey

```
|sigma_t(y)|<=200 a_x exp(d(x,y)),
|sigma_tt(y)|<=4000 a_x exp(d(x,y)), a_x=s_x L_x,
sup|sigma_t|<.4, sup|sigma_tt|<4.                         (12)
```

The taper (11) obeys the same generous constants. For the spatial
weight, do not assume ds/dchi<=s. Instead (8) gives
s_chi/s<=1+3sL, whereas
(log(sL))_chi=(s_chi/s)[1-s/((s+h)L)]<=1, since
6sL^2<1. The last inequality follows by maximizing sL^2 on s<=s1;
its maximum is below .015. This proves the exp(d) weight in (12).
Beyond the seam extend this comparison weight by exp(chi-chi1).

## The paired first-order realization and exact locality identity

For the retained Euclidean D2, J is the parallel involution of the
phase audit: it commutes with the differential expression and exchanges
the complementary adjoint domains. Thus

```
Dsharp=-J D2, Dsharp^dagger=-Dsharp,
-Dsharp^2=L2=D2^dagger D2.                               (13)
```

This is an operator statement, not just a formal square. The normal
Clifford matrix of Dsharp anticommutes with the actual projector
involution chi_b, while its tangential matrix commutes with chi_b.
The projected realization is therefore selfadjoint after multiplication
by i and elliptic. Squaring imposes the complementary D2 condition,
precisely the accepted Dirichlet/Robin square without a mass Robin
term. The regular bolt is an interior point. In particular
G_M=(Dsharp+M)^-1 has norm<=1/M on the actual paired exterior.
This argument is essential: the original D2 on its original domain
was not assumed anti-selfadjoint or accretive.

Per exterior, with the paired normalization already divided by two,
the massive canonical trace is M tr G_M=M^2 tr(L2+M^2)^-1,
plus the specified local anomaly. The rank is two. The symmetry T
which anticommutes with Dsharp makes the odd local trace vanish.

Suppose two such operators agree on a patch, and chi_c is one near
the observation point and zero before the edge of the common patch.
Put C_c=[Dsharp,chi_c]. The exact smooth diagonal identity is

```
(G_actual-G_cap)(x,x)
 =-[G_cap C_c G_actual C_c G_cap](x,x).                   (14)
```

It follows by applying the resolvent identity on both sides. The
remaining first-order term is
G_cap chi_c C_c G_cap=G_cap[Dsharp,chi_c^2]G_cap/2;
its kernel vanishes near the diagonal where chi_c=1. Thus (14)
does not cycle divergent traces, and the actual operator enters only
through its rigorous massive inverse bound. All insertions in (14)
are separated from x.

## Kernel bounds and the two different caps

For the complete H2 spin Green kernel write its invariant coefficients
as positive A,B. The equations, with the flat short-distance
normalization B~1/(2pi d), are

```
A'+(tanh(d/2)/2)A=-M B,
B'+(coth(d/2)/2)B=-M A.                                 (15)
```

They follow from the invariant Green construction of
[Muck, equations (14)-(15)](https://arxiv.org/html/hep-th/9912059).
Its decaying hypergeometric solution has A/B tending to 1.
The Riccati equation shows A<=B: a ratio above 1 has positive
derivative and cannot return to that limit. The flux equation gives
B<=1/[4pi sinh(d/2)]. Applying the equation for A+B from d/2
to d gives

```
k_M(d)=||G_H2,M||<=2 exp(-Md/2)/(pi d).                   (16)
```

For the sharper L1 bound, compare A/B with K0(Md)/K1(Md).
Their Riccati forcing terms are 1/sinh(d) and 1/d; backward
comparison from infinity gives A/B at least the flat ratio. Therefore
B<=M K1(Md)/(2pi). Integrating the A equation then gives
A<=M^2 K0((M-.5)d)/[2pi(M-.5)]. For M>=18 this implies

```
k_M<=1.03 M[K0((M-.5)d)+K1((M-.5)d)]/(2pi).
```

Using dV<=2pi d exp(d)dd and the integral representation of K_nu,
one obtains

```
||exp(Md/4) k_M||_1<=8/M,
||exp(Md/4) k_M||_(4/3)<=7/sqrt(M).                      (17)
```

For example the L1 estimate is at most
1.03M(1+pi/2)/(.75M-1.5)^2<8/M. The L4/3 estimate follows
from (16) and the elementary Gamma(2/3) integral. Young's inequalities
are legitimate on H2: lift radial majorants to bi-K-invariant
functions on the unimodular group SO_0(1,2), with K Haar mass one,
and descend the convolution bound to the quotient.

In the global cap, the conformal mass perturbation has
|exp(sigma_cap)-1|<=exp(.05)-1. Its weighted Neumann ratio is
less than .411 by (17). Terms with at least one insertion have
weighted row L2 norm below
49(exp(.05)-1)/(1-8(exp(.05)-1))<4.3.
The free row restricted to distances at least one has norm below .011.
Including conformal endpoint factors, use 10 for either offdiagonal
row. Interior Dirac estimates give first two point derivatives bounded
by 1e8 M^j, j=1,2. For an explicit check, rescale X=Mx in an analytic
normal chart. Fourier's identity ||grad u||_2=||D_flat u||_2 and
four nested cutoffs from B2 to B1 give
||partial^alpha u||_2<=10^j j!||u||_2, j<=4. Apply the whole-plane
bound |u(0)|<=(2sqrt(pi))^-1||(1-Delta)u||_2 to
(1-|X|^2)^3_+ partial^alpha u; its gradient is at most 6 and its
Laplacian at most 12. The resulting coefficient is below 2e5 for
j<=2. The nearby weighted row, conformal frame and area factors fit
inside the stated 1e8. These are explicit local estimates, not an
assumed derivative bound on the physical stress.

The second cap below is different: it is a point-centered normal
metric capped to FLAT outside radius 1/4. It is used for the common
UV subtraction. Its finite gluing difference from the global H2 cap
must be estimated; the two caps are not identified.

## Analytic normal charts and the complete local UV remainder

Put u=b/r and Bcal=B(u)/C, with B from the conformal-domain note.
The exact normalized identities are

```
|grad u|^2=Bcal(u), Hess u=Bcal'(u)g/2, R=-Bcal''(u).
```

In the holomorphic conformal normal chart g=exp(2sigma_N) dz dw,
sigma_N(z,0)=sigma_N(0,w)=0, these become the Goursat problem

```
u(z,0)=u0+pz, u(0,w)=u0+pbar w,
(sigma_N)_zw=Bcal''(u)exp(2sigma_N)/8,
u_zw=Bcal'(u)exp(2sigma_N)/4.                             (18)
```

For independent complex |u0-1|<=.01, |p|,|pbar|<=.01 and
|z|,|w|<=R0=1/4, the integral equations close
|u-1|<=.03, |sigma_N|<=.025. Indeed |Bcal'''|<10,
|Bcal''|<2.31 and |Bcal'|<.07 there; the two Picard row sums
are below .13. This proves analytic dependence rather than assuming
finite jet bounds. For physical parameters
p pbar=Bcal(u0)/4 and |p|<=.501s.

Write v_N=exp(sigma_N)-1, so |v_N|<delta_N=.026 and
v_N(0)=grad v_N(0)=0. Cap v_N by any fixed radial cutoff which
is one through 3/16 and zero past 1/4. Let v_c denote the result.
No cutoff derivatives are needed in the following global bounds:

```
|v_c(y)|<=c2 |y|^2, c2=4delta_N/R0^2<1.664,
|v_c-T5 v_N|<=10011 |y|^6,
|v_c(y)v_c(z)-T5[v_N(y)v_N(z)]|
 <=200(|y|+|z|)^6.                                     (19)
```

Here T5 is total degree at most five. For |y|<R0/2, Cauchy's
geometric tail gives 2delta_N R0^-6 |y|^6. Outside, bounded
v_c and the degree 2,3,4,5 polynomials give
(64+16+8+4+2)delta_N R0^-6<10011. For the product, the same
argument gives (64+4+4)delta_N^2 R0^-6<200. Odd total degrees
integrate to zero by simultaneous reflection. In particular degree
five is removed only after its Taylor coefficient is isolated.

The flat rank-two Green kernel is
M[K0(Mr)+gamma dot x K1(Mr)/r]/(2pi). Its weighted L1 norm
at weight M/4 is below 5/M and its weighted L3/2 norm cubed is
below 100/M. For a closed N-edge loop, writing L for total edge
length, Young and Holder give

```
integral L^6 product(k_M)
 <=720*4^6*100*5^(N-3)/M^(N+4), N>=3.                   (20)
```

This follows from L^6<=6!(4/M)^6 exp(ML/4). The spin trace
costs rank two. Three quadratic vertices supply the six distances;
the others are bounded by delta_N. Thus all Born orders at least
three, including their complete Neumann sum, contribute at most

```
2,949,120,000 c2^3/[1-5delta_N] * M^-4.
```

The Born-one degree-six remainder is below 612*10011 M^-4.
The Born-two degree-six remainder is below
589,824,000*200 M^-4. These numbers follow respectively by the
two-edge flat integral and (20) with three edges. Altogether the
local remainder is less than 1.4e11 M^-4, uniformly on the full
complex parameter polydisk. It includes every Born order.

For completeness the retained low terms can be checked directly.
At the normal center set A=Delta sigma_N, H=Hess sigma_N,
B4=Delta^2 sigma_N. The moments

```
integral r^2 tr[G(-x)G(x)] dx=-1/(3pi M^2),
integral r^4 tr[G(-x)G(x)] dx=-8/(15pi M^4)
```

give Born one A/(12pi)+(B4+A^2+2trH^2)/(120pi M^2).
The quadratic-quadratic Born-two term is
-(2A^2+trH^2)/(60pi M^2). One direct Fourier check uses
z=|p|^2, D=1+z, T=trH, S=trH^2: before the final factor 1/4
the angular trace is

```
16z^2(T^2+2S)(3-z)/D^7-64zT^2(1-z)/D^6
 +8(1-3z)T^2/D^5-16zS/D^5.
```

Its beta integrals give -(2T^2+S)/(15pi). The sum is therefore

```
-R/(24pi) -(R^2+2 Delta R)/(480pi M^2).                  (21)
```

The first term cancels the canonical per-mode anomaly. The flat
vacuum term is spatially constant at fixed M. This identifies the
common a0,a1,a2 subtraction; the nonzero modes have not been
independently renormalized with unrelated finite constants.

The scalar remainder is rotationally invariant in the normal chart,
so its analytic dependence on p,pbar is through beta=p pbar. Its
beta disk has radius 1e-4, while physical beta<2.6e-7.
At fixed h, |beta_s|<.00125 and |u0_s|<=1. Cauchy gives a
radial derivative below 125 times the remainder bound. Finally

```
sum_n d_n/M_n^4 <4C^2/q,                               (22)
```

by telescoping the identity
(q+2n)/[n^2(n+q)^2]=(1/n^2-1/(n+q)^2)/q.
Multiplication by C converts the normalized trace to b^2 tau.
The local remainder's derivative is consequently below
7e13(C/q)^3 q^2<3.1e6 q^2.

## The common local finite part

The M^-2 coefficient in (21) uses the common angular zeta function
Z(p)=sum d_n[n(n+q)/C]^-p. Subtracting
4C^p zeta(2p-1) leaves a convergent function at p=1 and yields

```
Res Z(1)=2C,
FP Z(1)=4C gamma_E-2C H_q+2C log C.                      (23)
```

The common proper-time Gamma/scale factors add the fixed elementary
constant/log multiple of 2C. Here mu0 b=1, as declared in the
retained reference; no new scale is chosen. All these factors are
bounded by 1000C at k=2000. The local derivative has the exact
quartic simplification

```
Delta R=-Bcal Bcal''''-Bcal' Bcal''',
d_u(R^2+2Delta R)=-4 Bcal' Bcal''''.                     (24)
```

On the bolt patch |Bcal'|<.005, |Bcal''''|<13 and
|u_s|<=1. Equations (21),(23),(24) bound this contribution by
100 q^2 with wide room: the displayed guards in fact give less
than 3e-6 q^2. Spatially constant vacuum terms do not
contribute to (1); they remain in the absolute stress/centering
accounting. The selected finite local prescription is retained.

## Distant-seam error: actual RN versus the global H2 cap

Use (14) with a cutoff which changes over chi width 1/4 before
chi1. Its physical commutator norm is below 9. For s<=.001,
the separation from its support is

```
D_x>=.946 log(.01/s)-.25.                               (25)
```

This follows from dchi/ds=exp(-sigma)/sqrt(P) and (5),(8).
The global cap row estimates above give, with the trace's M factor
included, a per-mode gradient error at most
1e13 M exp(-MD_x/2), and Hessian error at most
1e20 M^2 exp(-MD_x/2). The j=0 row is 10, not the much
larger generic point-evaluation constant. The latter is only used
for j=1,2. These estimates follow by differentiating the separated
rows of (14); the actual inverse remains bounded by 1/M.

For chi>=1, ds/dchi>.5s; (25) converts the gradient bound to
2e16 M exp(-.964M). The elementary angular estimate

```
sum d_n M_n exp(-.964M_n)<3.2e-6 q                     (26)
```

follows from M_n>=18sqrt(n), the decreasing exponential, and an
integral bound after n=2. Thus its normalized radial contribution
is below 2.3e8 q^2. For chi<=1, rotational invariance makes the
first derivative vanish at the bolt. Integrate the Hessian bound
from 0 to chi. The exact chi1=arcosh(.01/h) gives the factor
(h/.01)^(M/2)/h times at most exp(.625M) and the displayed
polynomial M^2; with h<1e-30 and M>=18 it is negligible.
This is an actual state/domain localization bound, not replacement
of an unknown actual kernel by the homogeneous kernel.

## Finite gluing: global H2 cap versus the local flat cap

This step is needed in addition to (25): the UV cap and the global
comparison cap are different metrics. Work in common H2-centered
Poincare coordinates and transport the local normal cap explicitly.
If v is the total conformal factor on the holomorphic axis, the
normalizing map satisfies

```
F(0)=0, F'(w)=exp[v(0)-2v(F(w),0)].                      (27)
```

For H2 itself v(z,0)=0. The complex guards (10) and the exact
identity

```
h cosh(chi(s))=
[s cosh(hD)+h sinh(hD)]/[cosh(hD)+(s/h)sinh(hD)]
```

give inverse s derivative below 2 on the required quarter disk.
Along the H2 holomorphic axis,
h cosh(chi_y)=h cosh(chi_x)+(h sinh(chi_x)/2)z.
Thus |z|<=.25 remains inside that disk; its first two center
derivatives give |v_t|<=200sL, |v_tt|<=4000sL. Since
|v(0)|<.008 and |v(z,0)|<.020, (27) gives
|F'-1|<.05 and |F-id|<.05|w| on |w|<=7/32.
The image lies inside .232<.25, so this bootstrap is strict.

Use a gluing cutoff which is one to 1/8 and linearly zero by 5/32.
It is Lipschitz with slope 32; this suffices for (14). The flat metric
cap is unchanged to 3/16, so the whole cutoff lies in the common
metric. Extend F-id by multiplying with the quintic radial step
which is one to 5/32 and zero past 7/32. Its radial derivative
operator bounds through order three are 30,1500,3e5. The extension G
therefore satisfies

```
||DG-I||<.439, ||DG^-1||<=2, ||D^2G||<22.                 (28)
```

It is a global diffeomorphism equal to the identity outside this
disk. First-center variations obey, with a_x=s_xL_x,

```
|G_t|<=260a_x, |DG_t|<=9200a_x, |D^2G_t|<=6e5a_x.
```

On chi_x<=1, the second variation has
|G_tt|<=7000a_x, |DG_tt|<=3.1e5a_x. Quadratic terms in
a_x are harmless there because s<=1.62h and a_x<1e-27.
These estimates follow by differentiating (27), using Cauchy on
its image margin .018, and then differentiating the explicit quintic
extension. For example |v_z|<.02/.018, and Gronwall in the
first-center equation gives |F_t|<=260a_x,
|F_t'|<=1400a_x, |F_t''|<=1.2e5a_x. Multiplication by the
quintic gives 1400+30*260<9200 and
1.2e5+60*1400+1500*260<6e5. At the baseline,
|F''|<2.5 and therefore |D^2G|<2.5+60*.05+
1500*.05*(7/32)<22. Second-center differentiation gives the
displayed 7000 and 3.1e5 bounds; its quadratic terms are
negligible in the explicitly restricted near-bolt range. These are
bounds on the actual normalizer, not freely specified coordinate jets.

For clarity about representation, in fixed H2 L2 use

```
A_sigma=exp(-sigma/2) D_H2 exp(-sigma/2),
S=(A_sigma+M)^-1,
delta S=(S f+f S)/2-M S f S, f=delta sigma.               (29)
```

The last identity is exact. Its first two terms are the endpoint and
density factors accompanying the conformal mass-potential insertion;
they are not dropped or counted as another anomaly. In the compact
normal-chart pullback, use the canonical unitary
U_G psi=Jac(G)^(1/2) psi composed with G. The transported differential
expression is exactly

```
D_G=1/2 {B^a,partial_a},
B^a=exp(-nu) gamma^i(DG)^-1_i^a,
nu=(v_H2+sigma) composed with G.                         (30)
```

The Piola identity partial_a[Jac(G)(DG)^-1_i^a]=0 gives the
zeroth-order divergence term. When differentiating transport, freeze
the explicit global sigma variation: nu_t=grad(v) dot G_t. Handle
the explicit sigma_t once, by (29). This avoids mixing canonical
and conformal representations or double counting the source.

On the image, Cauchy gives |grad v|<=4, |Hess v|<=400 and
|grad sigma_t|<=2e4a_x. Differentiating (30), with (28), gives

```
|delta B|_Cliff<6.1e4a_x,
|grad delta B|<2e7a_x,
|delta^2 B|_Cliff<2.1e6a_x       (chi_x<=1).               (31)
```

The divergence form makes its bilinear bound especially direct:
one half of |delta B| multiplies each of the two first-derivative
energy terms. The physical covariant resolvent obeys weighted
gradient bound 4 by Lichnerowicz and DS=I-MS; (28), including the
half-density connection, enlarges the ordinary transported energy
bound to 20. Explicitly canonical flat-volume conversion has
|v_H2|<.02, |dv_H2|<.13 and gives energy below 4.2. In the
pullback, DG<1.439 and
|grad log Jac(G)|/2<=||DG^-1|| |D^2G|<=44. Consequently
the transported energy is less than
1.439*4.2+44*1.2/18<10<20.
For centered weights W=exp(-2d_H2(x, .)), the resulting
two-sided inverse bounds are

```
M||W S_t W||<=2e6a_x,
M||W S_tt W||<=1e8a_x       (chi_x<=1).                  (32)
```

These are weighted bounds. No unweighted global O(a_x) norm is
asserted. A detailed constant ledger is as follows. Conjugated
resolvents with one or two distance weights have bound 1.2/M.
The compact first insertion costs at most
1.2*20*6.1e4 a_x=1.464e6a_x. Equation (29), using
|f|<=200a_x exp(d), costs another 528a_x. The endpoint and
weight reserves fit below 2e6a_x. The direct compact second
insertion costs 1.2*20*2.1e6a_x=5.04e7a_x. The global
linear second source costs 10560a_x. For double insertions use
the weighted O(a_x) bound for one vertex and the uniform global
first-vertex bound 3 for the other; this costs at most
2*(2e6a_x)*3=1.2e7a_x. Compact quadratic terms are negligible
in the stated near-bolt range. Remaining endpoint terms fit within
1e6a_x, leaving total below 6.4e7a_x<1e8a_x.

Here the coefficient bounds in (31) have an elementary matrix
check. Put E=DG and K=E^-1. Then K_t=-K E_t K and
B_t=exp(-nu) gamma[-nu_t K-K E_t K]. The two terms have
sizes at most 2080a_x and 36800a_x before the Clifford/Frobenius
factor. Multiplying by sqrt(2)exp(.07) gives less than 6.1e4a_x.
For its gradient use |K_j|<=88 and
|nu_t,j|<=400*1.439*260a_x+4*9200a_x<1.87e5a_x.
The product-rule sum before the two-dimensional factors is below
6.4e6a_x, hence below 2e7a_x afterwards. In B_tt the linear
terms are bounded by 56000a_x+1.24e6a_x; multiplying by
sqrt(2)exp(.07) is below 2.1e6a_x. The remaining terms are
O(a_x^2) with coefficients below 1e12, negligible for a_x<1e-27.
This also specifies the norm factors used in the transport ledger.

For the double-insertion estimate with the exact weights, factor
W S A_t S A_t S W as
(W S A_t S W)(W^-1 A_t S W). For the global conformal piece,
A_t S=-f+M f S-[A,f]S/2. Equations (9),(11) imply |df|<3
(even 10 would suffice), so its second factor is below 3; the
compact part is O(a_x) and negligible in this near-bolt use.
Thus the second-product ledger does not assume an unweighted
global perturbation is O(a_x).

The local flat-cap rows also reserve both weights in (32). At weight
M/4+2 their flat L1 bound is below 7/M and L4/3 norm squared
below 56/M. Hence their nonfree row part is below
56delta_N/(1-7delta_N)<1.8. The free far row at r>=1/8 is
below 1.5, by (8/pi)E1(.625)<2.181 for its squared norm.
Endpoint factors fit below 4 in total. After the commutator slope
33 and the small distance-conversion factors, use

```
||W^-1 r_M||_2<=R_M=400 exp(-M/32).                      (33)
```

In this conversion d_H2(F(w))=2atanh(|F(w)|/2)<=1.06|w|
on the row annulus. The extra factor exp[2(d_H2-|w|)] is
below 1.02. In particular the Euclidean normal radius is not
silently identified with the hyperbolic radius.

Intrinsic normal-chart row derivatives do not carry the global log:
Cauchy in (18), with margins .0089 and physical parameter derivative
sums at most 2.3s and 3s, gives conservatively
|r_t|<=1e3s R_M and |r_tt|<=1e4s R_M. These are L2 bounds
with the same reserved weights. They include the cutoff commutator;
the cutoff is fixed in its own normal chart.

Insert (32),(33) in (14), keeping the external trace-M and rank two.
The resulting finite gluing scalar has first derivative below
1e12 a_x exp(-M/16) and near-bolt second derivative below
4e13 a_x exp(-M/16). These envelopes include the own-chart
rows and endpoint normalizations. For the exact large charge in (2),

```
sum_n d_n exp(-M_n/16)<4q.                              (34)
```

Indeed sum exp[-(18/16)sqrt(n)]<1.42 by the decreasing
integral test. The 4sum n exp[-(18/16)sqrt(n)] term is less
than 100 and is negligible compared with q=2000^4.

For chi>=1 divide the first derivative by ds/dchi>.5s.
For chi<=1, use rotational invariance of the scalar gluing difference:
its first derivative vanishes at zero. Integrate its second bound,
use s<=1.62h and ds/dchi>=exp(-.05)h sinh(chi), and use
L<81. With the rank and endpoint reserves already specified this
bounds the physical b^2-normalized gluing derivative by

```
K_glue q^2, K_glue<2.2e14.                              (35)
```

For example the direct near-bolt arithmetic is below
4e13*1.62*81*exp(.05)*4*(C/q)<7.9e13;
(35) retains further room rather than using that sharp endpoint.
The first-derivative region is much smaller. This is the finite
transport estimate missing from an identification of the two caps.

## Consequence and remaining scope

The local complete remainder, common local restoration, distant seam
and finite gluing bounds sum to less than 3e14q^2, proving (1).
There is no assumption that the actual state equals the homogeneous
H2 state: the actual inverse is retained in (14), its domain is (13),
and the complete cap comparison and common subtraction are bounded.

The independently reviewed
[non-bolt continuation](stationary-reference-nonbolt-stress.md)
matches (1) at s=.001. The
[complete reference result](stationary-reference-higher-mode-clock.md)
combines the two regions, the common trace-value anchor and exact HH
conservation into individual source bounds and the integrated
first-order reference clock bound. The
[outer-boundary analysis](stationary-massive-boundary-layer.md)
retains the separate boundary layer and Weyl cancellation. The
four-dimensional Weyl sources and radius allocation are retained
once, as in [the radial-source note](stationary-weyl-radial-sources.md).
The actual higher-order inverse or controlled spatial order reduction,
metric response, corrected exterior/state and full-action charge
remain separate open obligations.

Remaining to-do list: complete focused verification of the combined
ordinary reference result; quantify the finite mass and horizon
recentering; establish the corrected-exterior solver, its metric
response bounds and full boundary charge/state; retain the original
driven preparation, capture, tail and renormalized spatial-area tasks.
