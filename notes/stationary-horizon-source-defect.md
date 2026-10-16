# Quantitative horizon source defects on two explicit references

2026-10-15, S139/S074. Ordinary argument independently reviewed before
focused verification. On the retained classical RN exterior, the complete
Dirac HH source at the bolt differs from the homogeneous product
source, at the same physical sphere and orbit-curvature radii, by
less than 5.1e-12 of the prepared thermal horizon forcing. Both
orbit and sphere pressures obey this bound. The comparison includes
the complete massive tower, the LLL, the retained Weyl functional,
and their common subtraction at fixed mu0=1/b0.

The final section extends the comparison to the exact quantum-centered
covariant seed in
[the comparison note](stationary-centered-covariant-seed.md), with its
independently determined two product radii. The horizon source in the
selected HH state on each explicit reference is controlled. Neither reference
is thereby a corrected exterior. The response to the subsequent metric
correction and a justified inverse remain separate obligations.

## Data and the statement

Use one exterior's normalization, as in the frozen
[cold-bolt proof](stationary-reference-cold-bolt.md),
[common subtraction proof](stationary-reference-nonbolt-stress.md),
and [horizon equations](stationary-horizon-product-matching.md).
The actual Euclidean domain and its paired first-order realization
are the ones in those notes. The bolt is smooth, the exterior state
is the selected static HH state, and the paired phase allocation is
the retained one. No reflecting throat cut is substituted for the
physical outer boundary.

```
k=2000, q=k^4, eta=1/2000, z_i in [1/4,7/25],
xi=b^2/L4^2, C=1+6xi=sqrt(1+12pi k^6),
b^2/G=xi k^10, ell=b/sqrt(C), mu0=1/b,
F(r)=1+r^2/L4^2-2m/r+Q/r^2,
Q/b^2=1+3xi, m/b=1+2xi+delta/2,
R=r_h=b(1+h), kappa=F'(R)/2=z_i nu.
```

The accepted candidate inequalities give

```
eta/(8k^9)<=h<=eta/k^9<1e-30,
xi>=k^3, C/q<.00352, M_n^2=n(n+q)/C, M_1>18,
d_n=2(q+2n), n>=1.
```

Let a_h^-2=F''(R)/2. Denote by pA_RN,pS_RN the actual
four-dimensional horizon pressures, and by pA_prod,pS_prod the
complete homogeneous Dirac pressures on AdS2(a_h) x S2(R), in the
same theory, fixed flux and scale. Then

```
|pA_RN-pA_prod|, |pS_RN-pS_prod|
       < 1e14 h q^2/R^4,                               (1)

max(|8piG delta pA|,|8piG delta pS|)/(2kappa/R)
       < 5.1e-12 < .01.                                (2)
```

The numerical inequalities are uniform in the stated z_i interval.
All source comparisons below are at the same R. The homogeneous
comparison is not obtained by resetting mu0 to 1/R.

## A real metric homotopy with a small weighted amplitude

Recall the exact cold-cap variables s=r/b-1, s1=.01 and

```
P=(s-h)(s+h)q0(s),
lambda=kappa b/C=h q0(h)/(1+h)^2,
v=log tanh(chi/2)=-atanh(h/s)-hD(s),
D(s)=integral_s^s1 (J(t)-1)/[(t-h)(t+h)] dt,
J=q0(h)(1+s)^2/[q0(s)(1+h)^2],
exp(sigma)=A(s)[cosh(hD)+(s/h)sinh(hD)],
A(s)=sqrt(q0(s))(1+h)^2/[q0(h)(1+s)^2].                 (3)
```

The global H2 cap is exp(2sigma_cap)gH2 in units ell. It equals the
actual orbit metric through s1 and becomes exactly H2 after the
accepted unit-width quintic extension. Its real bounds are

```
|sigma_cap|<.05, |grad sigma_cap|<.4,
|Hess sigma_cap|<1.5.
```

Set

```
L0=1+log[(s1+h)/(2h)]<81,  A0=h L0<1e-27.             (4)
```

The notation A0 here is an amplitude, not a radius or heat
coefficient. At the bolt, (3) gives exactly

```
sigma_h=-log(q0(h))/2+hD(h), |sigma_h|<3A0.             (5)
```

For completeness, the weighted envelope used below follows directly
from (3), rather than from a spatial derivative estimate. The real
bounds 0<=D<=2.03 log[(s1+h)/(s+h)],
q0(h)-1<=.501h and q0'<=.921 give
|log A(s)|<=5s. If H=cosh(hD)+(s/h)sinh(hD), then
log H<=(s+h)D, so |sigma(s)|<=10sL(s), where
L(s)=1+log[(s1+h)/(s+h)]<=L0. Also

```
h cosh chi=[s cosh(hD)+h sinh(hD)]/H,
sD<=2.03(s1+h)/e<.008.
```

Thus s<=1.008h cosh chi. This proves the following envelope in
the actual patch. In the cap extension it follows instead from
|sigma_cap|<.05 and h exp(chi1)>=h cosh(chi1)=s1.

```
|sigma_cap(y)| <=20 A0 exp(d_H2(0,y)).                  (6)
```

Consider the genuine metric homotopy sigma_zeta=zeta sigma_cap,
0<=zeta<=1. This is a real family of complete capped metrics. Its
variation is the amplitude sigma_cap, not a Killing derivative of
sigma_cap and not a derivative with respect to temperature of an
unknown exterior solution.

## The local complex amplitude disk

Use the H2 holomorphic normal coordinates with

```
gH2=dz dw/(1-zw/4)^2,
cosh chi=(1+zw/4)/(1-zw/4).
```

On |z|,|w|<=.26, |cosh chi-1|<.035. The accepted inverse-s
bound of 2 on the quarter s-disk therefore gives |s-h|<.07h.
This remains well inside that disk's radius h/2. In this smaller
disk the explicit formulas (3) give |sigma_s|<100L0. To see the
guard, use |D|<=7L0, |D'|<=4.5/|s+h|,
|q0'/q0|<3 and |1+s|>.99. Then |H-1|<1e-24,
|H'/H|<10L0, and |(log A)'|<4. The stated 100L0
bound has ample room. Equations (4)-(5) imply

```
|sigma(z,w)|<10 A0, |z|,|w|<=.26.                     (7)
```

Radial analyticity gives sigma(z,0)=sigma_h exactly. The
holomorphic normalizer of the amplitude family is consequently
especially simple:

```
F_zeta(z)=exp(-zeta sigma_h) z.                         (8)
```

In the new normal coordinates its conformal factor is

```
v_zeta(z,w)=vH2(F_zeta z,F_zeta w)
 +zeta[sigma(F_zeta z,F_zeta w)-sigma_h].               (9)
```

For every real center zeta in [0,1], extend zeta to its complex disk
of radius

```
rho=1e-4/A0.                                          (10)
```

This extension is used ONLY for the local normal metric and its
flat cap. It is not an extension of the global metric by that
enormous amplitude. On (10), |zeta sigma_h|<.000301,
so the image of |z|<=.25 lies inside .251<.26. Equations
(7)-(9) give |v_zeta|<.018<.025 on the normal quarter disk.
In particular exp(v_zeta)-1 has the same .026 bound, quadratic
vanishing and convergent flat Born expansion as the frozen local
proof.

Let E_M(zeta) be its complete local canonical trace remainder after
the flat vacuum and the two displayed local terms in (21) of the
cold-bolt proof. That proof, valid for these complex normal metrics,
gives |E_M|<1.4e11 M^-4. Cauchy's inequality on (10), integrated
along the real interval, therefore gives conservatively

```
|E_M(1)-E_M(0)|<3e15 A0 M^-4.                         (11)
```

This includes every Born order. The same local complex disk applies
to the separated flat-cap rows r_left,r_right in their own normal
chart. With the frozen row bound R_M=400 exp(-M/32), it gives

```
|partial_zeta r_left|_2, |partial_zeta r_right|_2
            <=1e5 A0 R_M.                             (12)
```

Left and right rows are separately holomorphic kernel families.
No antiholomorphic adjoint is differentiated as if it were
holomorphic. Their real-family adjoint relation may be used only
after these estimates are obtained.

## Global inverse variation and finite gluing

The middle, global inverse is differentiated ONLY along the real
homotopy. In fixed-H2 L2 put

```
A_sigma=exp(-sigma/2) D_H2 exp(-sigma/2),
S=(A_sigma+M)^-1, W=exp(-2d_H2(0,.)),
S_zeta=(S f+f S)/2-M S f S, f=sigma_cap.                (13)
```

The weighted similarities have norm at most 1.2/M, uniformly for
this real family. If |f|<=a exp(d), then |W f W|<=a and (13)
implies

```
M ||W S_zeta W|| <=(1.2+1.44)a=2.64a.                 (14)
```

The endpoint terms in (13) are retained. In particular (14) is an
actual metric response of the inverse, not a relabeling of the
frozen spatial-source derivative.

To compare with the local flat cap, extend the dilation (8) by the
same radial quintic used in the cold proof: one through 5/32 and
zero after 7/32, with derivative bounds 30 and 1500. Writing
G_zeta(y)=y+chi(|y|)[exp(-zeta sigma_h)-1]y gives

```
|G_zeta'|<A0, |DG_zeta'|<24A0,
|D^2G_zeta'|<1300A0, ||DG_zeta-I||<31A0.               (15)
```

Primes in (15) denote amplitude variation. For example
|partial_zeta exp(-zeta sigma_h)|<3.1A0; multiplying by
1+(7/32)30 and by 60+(7/32)1500 gives the latter two
constants. The map is a global real diffeomorphism for 0<=zeta<=1.

Use the canonical unitary pullback by G_zeta. In its flat-volume
representation the differential expression is
D_G={B^a,partial_a}/2, with
B^a=exp(-v composed G) gamma^i(DG)^-1_i^a. Split the explicit
sigma amplitude, already handled in (13), from the coordinate
transport of v. The frozen guards |dv|<=4, |D^2v|<=400 and
(15) give

```
|delta B|_Cliff<45A0, |grad delta B|<1e4A0.             (16)
```

Here is a direct ledger. With K=(DG)^-1, |K|<1.01,
|K_zeta|<25A0, |nu_zeta|<4A0 and
|grad nu_zeta|<500A0. Even using the coarser frozen
|grad K|<=44, differentiation of K_zeta=-K DG_zeta K
gives |grad K_zeta|<3500A0. In the product rule for delta B,
the sum before the Clifford/exponential factor is below 4300A0;
that factor is below 1.52. This proves (16).

The divergence form makes the first-order insertion a bilinear
energy form; the weighted graph-energy bound is 20. Its inverse
variation costs at most 1.2*20*45A0=1080A0 in M units.
Equation (14) costs less than 56A0, including the tiny G change
of the distance weight in (6). Endpoint and transport reserves give

```
M ||W S_transported,zeta W||<1e4 A0.                   (17)
```

Keep the point-density conversion explicit. In fixed-H2 L2 the
physical diagonal scalar is exp(-2zeta sigma_h) M tr S(0,0).
The canonical pullback has Jac(G)(0)=exp(-2zeta sigma_h), so
its normal-coordinate kernel already contains exactly this factor.
Equivalently, each normal row contains the endpoint
exp(-zeta sigma_h). The factor has norm below 1.01 and its
derivative is below 6.1A0. Its contribution is included in (12)
and may also be bounded separately by 6.1A0 times the unvaried
finite gluing scalar. It is not a second anomaly, a change of mass,
or a resetting of the renormalization scale.

The exact separated identity (14) of the cold proof now applies
between the transported global cap and the local flat cap. Its
scalar trace has one left row, one middle inverse and one right
row. Equations (12),(17), rank two and the external M give

```
|partial_zeta Phi_M|
 <=2(2*1e5*1.2+1e4) A0 R_M^2
    + the stated density reserve
 <1e12 A0 exp(-M/16).                                 (18)
```

For reference the first term is 8e10 A0 exp(-M/16), and the
separate density reserve is below 3e6 A0 exp(-M/16).
Thus (18) has a strict margin. At zeta=0 the global cap is exactly
H2, rather than an arbitrary comparison state.

## Common subtraction and the massive value ledger

The common angular zeta and scale remain those of the frozen
nonbolt proof. In dimensionful notation

```
zeta_h(s)=sum_n d_n[n(n+q)]^-s,
zeta_h(-1)=1/30-q^2/6,
Fang=FP zeta_h(1)=4 gamma_E-2 H_q, |Fang|<100,

tau_hat,mass= integral A(t)[k2D-(k0/t+k1+t k2)/(4pi)]dt
 +(1/(4pi)){(k0/b^2)[zeta_h'(-1)+zeta_h(-1)]
                                      +b^2 Fang k2},
k0=2, k1=-Rhat/6,
k2=-Rhat^2/120-Boxhat Rhat/60.                         (19)
```

The a0 restoration in (19) is not set to zero: it is the same
constant for both orbit metrics and therefore cancels in this
particular value difference. The canonical per-mode anomaly and
the common k1 restoration cancel as in the frozen proof. For k2,
the finite +2 from Gamma(s+1)/Gamma(s) and Res zeta_h(1)=2
is canceled by the common a4 subtraction. Thus Fang, not a
separately chosen modewise finite part, occurs in (19).

On the local complex disks above, the normal-factor Cauchy bound
gives |Rnorm^2+2Delta Rnorm|<1e4. Indeed
|Rnorm|<4 and |Delta Rnorm|<1000 follow from four derivatives
of |v|<.025 on radius 1/4. Applying amplitude Cauchy and
|Fang|<100 bounds the complete retained k2 value difference by
1e12 A0 C^2 in Theta=b^2 tau_hat units. This deliberately
overbounds the direct factor 1/(480pi).

The exact angular sums needed for (11),(18) are

```
sum d_n/M_n^4=2C^2 H_q^(2)/q<4C^2/q,
sum d_n exp(-M_n/16)<4q.                              (20)
```

Multiplication by b^2/ell^2=C converts the normalized orbit
trace to Theta units. Hence the local complete remainder, finite
cap gluing and local restoration contribute respectively

```
1.2e16 A0 C^3/q,   4e12 A0 Cq,   1e12 A0 C^2.          (21)
```

The actual exterior has not been replaced by its cap. At the bolt,
the physical boundary/attachment is separated by
D_h>=.946 log(.01/h)-.25>=.9 log(.01/h). The exact separated
identity, inverse norm 1/M and the offdiagonal row constant 10
give its value error at most 2e4 exp(-M D_h/2) per mode.
This includes both commutators, each of norm below 9:
2M(10*9)^2/M=16200<2e4.
For h<1e-30 and M>=18,

```
exp(-M D_h/2)<1e-180 h exp(-M/4),
sum d_n exp(-M_n/4)<4q.
```

This is below hq^2 after the C conversion. It is a bound on the
actual-domain inverse, not a thermal-state identification.

Dividing (21) by hq^2 and using L0<81,C/q<.00352 gives

```
1.2e16*81*.00352^3 <4.24e10,
4e12*81*.00352      <1.141e12,
1e12*81*.00352^2    <1.004e9.
```

Their sum, with the actual-domain error, is less than 1e14.
Consequently the complete massive comparison satisfies

```
|Theta_RN,h-Theta_H2(ell)|<1e14 h q^2.                 (22)
```

The homogeneous side of (22) has fixed product sphere b in the hat
frame. After the same constant physical Weyl rescaling at the bolt,
its sphere is R and its orbit radius is R/sqrt(C).
This intermediate radius is not identified with a_h: the actual
hat curvature is Rhat_h=-2(w-j)/b^2, whereas the product obtained
by physically matching a_h has Rhat_product=-2w/b^2. Their
difference 2j/b^2=4kappa R/b^2 must be retained.

## LLL, local Weyl source and physical product-radius transfer

Static HH regularity gives t_null(R)=0. Therefore the orbit pressure
from a product-frame trace is

```
pA_hat,h=b^2 tau_hat,h/(8pi R^4)=Theta_h/(8pi R^4).      (23)
```

For the massless LLL its trace is q Rhat/(24pi); it is state
independent, and (23) applies in the selected regular state. The
curvature difference in the cap comparison is bounded directly
from the same local complex disk, giving conservatively
|delta Theta_LLL|<1e5 hq^2. A finite Casimir integration constant
is not introduced at the bolt: HH regularity fixes t_null(R)=0.

For the local Weyl contribution use the exact compact-variation
sources in [the Weyl radial note](stationary-weyl-radial-sources.md).
Define

```
x=R^2/L4^2, y=Q/R^2, j=R F'(R)=1+3x-y=2kappa R,
sigma_h=log(R/b).
```

The sigma_h in this paragraph is the physical four-dimensional
Weyl factor, distinct from the cap factor in (3)-(18). Substituting
m/R=(1+x+y)/2 into the frozen E_alpha gives

```
R^2 E_alpha,RN
 =sigma_h[-4xy/5-q^2/3]+j(72x+22-36j)/180.             (24)
```

For the constant-Weyl product with sphere R and orbit R/sqrt(C),
the orbit Weyl pressure is

```
pA_W,prod=sigma_h[(1-C^2)/15-q^2/3]/(16pi^2 R^4).
```

This follows by varying the constant-Weyl product density
sigma_h[b^2 Rhat^2/60-Rhat/18+(1/15-q^2/3)/b^2]/(4pi).
The linear curvature term is topological in this compact variation;
the curvature-squared term contributes with the opposite sign to
the constant potential in the orbit pressure.

Since xy=xi(1+3xi)=(C^2-1)/12, the logarithmic and gauge
terms in (24) cancel exactly in the comparison. Thus

```
delta pA_W=j(72x+22-36j)/(2880pi^2 R^4).                (25)
```

For the stated candidate this has magnitude less than hq^2/R^4.
It is a compact-source result; no finite-cut boundary charge is
inferred from it.

Finally transfer the complete homogeneous pressure from orbit
radius R/sqrt(C) to the physically matched a_h. On RN

```
R^2/a_h^2=-1+2y,
|log[a_h^2/(R^2/C)]|<4h.                              (26)
```

The complete homogeneous theorem, with flux, sphere R and mu0 fixed,
gives |D A|<=100(1+|L|), where D is the orbit-curvature logarithmic
derivative and L=log(mu0^2 R^2/q). Here |L|<32. Thus the full
product transfer is below 2e4 hq^2/R^4. It includes every angular
level exactly once and does not require a separate estimate of
each product allocation.

Equations (22)-(26) imply the orbit inequality in (1), with a large
margin: its massive term has the additional factor 1/(8pi).

## Sphere source, derivative anomaly and the thermal target

The full retained trace fixes pS=A4/2-pA at the regular bolt.
On this particular RN reference R4=-12/L4^2 is constant. The
homogeneous product also has constant R4; hence Box R4 vanishes
on both. This is not a justification for omitting the retained
derivative anomaly on other metrics.

Set w=R^2/a_h^2=2y-1. The two complete trace anomalies give

```
16pi^2 R^4 A4_RN
 =7(1+x-y)^2/30-11x^2/15+11y^2/45-q^2/3,
16pi^2 R^4 A4_prod
 =(w^2+1)/15+w/9-q^2/3,
16pi^2 R^4 delta A4=j(19j-72x)/90.                    (27)
```

The magnetic term cancels. The scalar curvatures themselves differ
by -4j/R^2; constancy, rather than equality, makes their Box terms
zero. The anomaly difference in (27), together with the already
much sharper margin in the orbit estimate, proves the sphere
inequality in (1). No independent nonlocal angular-source estimate
is required for this fixed-reference comparison.

Using the exact thermal law from (3), q0(h)>=1 and R=b(1+h),

```
h/(kappa R^3)=1/[C b^2 q0(h)(1+h)]<=1/(C b^2),
Gq^2/(Cb^2)=1/(C xi k^2)<=k^-8.
```

Consequently each ratio in (2) is bounded by
4pi*1e14/k^8<5.1e-12 at k=2000. This compares the complete
physical source defect to the prepared nonzero thermal gap; it
does not expand a quantum center translation in that smaller gap.

## Extension to the covariant quantum-centered seed

Here the reference itself changes; it is not inferred from a spatial
derivative of the RN source. Keep mu0=1/b0, the original flux and the
selected finite coefficients fixed. Use the exact product radii and
outer curvature of the covariant-seed note:

```
B=b_star, A=a_star, x=lambda_e B^2,
C_star=B^2/A^2=1+6x+d, e=d/C_star,
|d|<1, |e|<1e-24,
F_ext(t)=(t-1)^2/t^2[1+x(t^2+2t+3)+d/t^2],
F_delta=F_ext-delta/t, t=r/B,
R=B(1+h), kappa=F_delta'(R)/2.
```

The prime in the last line is a physical r derivative. In this section
replace b,C,ell by B,C_star,A in the preceding cap estimates. The
conclusion is the SAME bound (1)-(2), with those replacements and
the actual HH source on this explicit seed. The product comparison
has sphere R and orbit radius a_h^-2=F_delta''(R)/2.

The seed proof supplies global exterior positivity and its unique
thermal root at the retained laboratory kappa. Its proper angular
gap and paired domain have the same construction as before. The
relative product-radius changes are below 1e-21; hence the loose
preparation bounds eta/(8k^9)<=h<=eta/k^9, L0<81,
C_star/q<.00352 and M_1>18 all survive. In detail
gamma=kappa A^2/B=(kappa b0/C0)u_star/sqrt(z_star), and
h=gamma[1+O(gamma)] with the explicit bound in the seed note.
No fixed-horizon expansion across the center displacement is used.

The exact normalized cold polynomial is now

```
P(s)=s^2 f(s)-delta_tilde(1+s),
f(s)=1+a s+c s^2+e g(s), g(s)=(1+s)^-2-1,
a=4x/C_star, c=x/C_star,
delta_tilde=h^2 f(h)/(1+h),
P=(s^2-h^2)q0(s),
q0=f(s)+h^2/(s+h)[f^[1](s,h)-f(h)/(1+h)].             (28)
```

Because |d|<1, 0<a<2/3 and 0<c<1/6. Let q_RN use the same
x,h but d=0, and set a_RN=4x/(1+6x), c_RN=x/(1+6x).
Then a=(1-e)a_RN,c=(1-e)c_RN, so q0=q_RN+e r.
The exact extra divided-difference term from g is

```
q_e=g(s)+h^2/(s+h){
 -(2+s+h)/[(1+s)^2(1+h)^2]+h(h+2)/(1+h)^3}.
```

Together with the a,c rescaling this gives directly

```
|r(s)|<4.2s, |r'(s)|<4, |r(h)|<4.2h,
                 h<=s<=.01.                           (29)
```

For example |g|<=2s, the divided difference is bounded by its
displayed positive denominators, and the polynomial a,c terms
contribute less than 1.1s. On the accepted complex half s-disk,
the same rational formula gives |r|<10s0. Cauchy on its quarter
disk then gives
|r^(j)|<=j! 10s0[4/(s0+h)]^j. These bounds make the dependence
on h explicit; an O(e) estimate alone would not suffice here.

Applying (29) to J,D and (3), with the same seam condition, yields

```
|J_seed-J_RN|<=5|e|(s-h),
|D_seed-D_RN|<=5|e|log[(s1+h)/(s+h)],
|sigma_seed,h-sigma_RN,h|<=8|e|hL0,
|sigma_seed(s)-sigma_RN(s)|<=20|e|sL(s).               (30)
```

The zero in J_seed-J_RN at s=h is exact; differentiating its
explicit quotient and using (29) proves the first inequality. Divide
it by (s-h)(s+h) and integrate for the second; insert these in (3)
for the latter two. The seam chi1=arccosh(.01/h) is EXACTLY
unchanged because D(s1)=0. Thus (30) preserves the h-weighted
envelope even though |e| is much larger than h.

Here are explicit jet and global-domain guards. On the quarter
s-disk, |delta sigma|<=100|e|s0L(s0). At s1, Cauchy with
margin s1/4 gives |delta sigma_s|<=400|e|L and
|delta sigma_ss|<=3200|e|L/s1. Since s_chi<=1.01s and
|s_chichi|<=2s, the seam chi-jets through order two, including
the change of s_chi, are below 1e8|e|. The same linear quintic
extension has C2 norm below 100, so its total global C2 change
is below 1e10|e|<1e-14. We may use

```
|sigma_cap|<.051, |grad sigma_cap|<.401,
|Hess sigma_cap|<1.501.                                (31)
```

The nonfree weighted H2 row correction is then
49delta_N/(1-8delta_N)<4.5, where here
delta_N=exp(.051)-1. The free distant row and endpoint factors
still fit inside 10. The weighted similarity bound 1.2 follows
from M/(M-2exp(.051))<1.2 for M>=18. The graph-energy
reserve 20 and local |dv|<=4,|D2v|<=400 have larger margins.
The normalized seed conformal polynomial on |u|<=1.25 is

```
Bcal(u)=[x+u^2-(2+4x+delta)u^3
                +(1+3x+d)u^4-2d u^5+d u^6]/C_star.
```

Its coefficient and derivatives through order six differ from the
normalized RN polynomial by less than 1e6|e| on that disk. This
follows by summing the at most degree-six monomials and the
coefficient change C_star/(1+6x). The previous strict complex
chart guards therefore persist. Equations (30)-(31) give
|sigma_h|<3.001A0, the local complex bound can be enlarged to
10.1A0, and |v_zeta|<.019<.025 in the local amplitude disks.
The derivative 3.1A0 in (15), row constants, restoration bounds
and the density derivative 6.1A0 are unchanged. The distant-seam
bound may be weakened to
D_h>=.94log(.01/h)-.3>=.9log(.01/h), which suffices for
the same negligible actual-domain tail estimate. Thus the actual positive-square
inverse and all three terms in (21) apply to the seed itself.

The scale needs explicit bookkeeping. For hat sphere B at the
retained mu0, (19) has the replacements

```
zeta_h'(-1)+zeta_h(-1)
 -> zeta_h'(-1)+[1+2log(mu0 B)]zeta_h(-1),
Fang -> Fang+4log(mu0 B).                              (32)
```

They follow respectively from Gamma(s)/Gamma(s+1)=1/s and
from (s+1)zeta_h(1+s)(mu0B)^(2s), with residue 2. The
finite +2 is still canceled by the common anomaly. Since
|log(mu0B)|<1e-21, the effective Fang is still below 100.
The a0 term, including its scale dependence, still cancels between
the two orbit metrics with this same sphere B. There is no new
per-mode subtraction or a reset of mu0. Thus (22) holds with
B,C_star and the seed's actual inverse.

For clarity the remaining local source transfer can be done without
an RN identity. At any alpha=0 smooth bolt write

```
j=R F'_h, w=R^2 F''_h/2, l3=R^3 F'''_h,
c_h=w-j, sigma4=log(R/B).
```

Varying the same W_bulk before evaluation gives exactly

```
R^2 E_alpha,W
 =sigma4[(1-c_h^2)/15-q^2/3+j l3/30]
                         +j(c_h/15+1/18).             (33)
```

One direct Euler check is to write H for the braces in the Weyl
Lagrangian. At the bolt
R^2 E_alpha=H-j(r H_Rbar')+(2w-2j)H_Rbar-jH_T;
rH_Rbar'=f'(Rbar)-sigma4 l3/30+(j+2w)/30,
where f(Rbar)=Rbar^2/60-Rbar/18+1/15-q^2/3.
This reduces to (33), including its non-logarithmic term.

On the seed put v=1+h and Qhat=1+3x. Eliminating delta using
the horizon equation gives

```
w=-1+2(Qhat+d)/v^2-10d/v^3+9d/v^4,
l3=6+6xv^2-18(Qhat+d)/v^2+108d/v^3-114d/v^4.
```

Consequently |w-C_star|<3C_star h,
|log(C_star/w)|<4h, |l3|<20C_star and j<3C_star h.
The LLL trace transfer remains within its earlier bound. Equation
(33), compared to the constant-Weyl product with curvature C_star,
gives a Weyl orbit-pressure difference below 2hq^2/R^4:
use |sigma4|<=h and
|C_star^2-(w-j)^2|<=18C_star^2 h. The complete
homogeneous orbit-radius transfer is again below 2e4 hq^2/R^4
at fixed mu0. Thus the orbit inequality (1) holds for the seed.

The angular transfer must retain the derivative anomaly. The exact
generic horizon identities are

```
Box R4= j(-l3-8w+2j-4)/R^4,
16pi^2 R^4(A4_h-A4_product(a_h,R))
                      =j(j-36w-3l3)/90.               (34)
```

For an independent algebra check, the orthonormal curvature sections
give Ric4^2 R^4=2(w+j)^2+2(1-j)^2,
C4^2 R^4=4(w-j-1)^2/3 and E4 R^4=4j^2-8w.
The radial derivative of R4 at a smooth alpha=0 bolt gives the
first line of (34). Inserting these in the retained Dirac anomaly,
including its fixed Box coefficient, gives the second. With the
displayed seed jet bounds it is below hq^2/R^4 in absolute trace
units. Therefore pS=A4/2-pA obeys (1) as well.

Finally, the exact thermal identity still gives
h/(kappa R^3)<=1/(C_star B^2) up to a factor less than
1+1e-21; alternatively q0(h)>=1-1e-21 suffices. The product
bounds imply Gq^2/(C_star B^2)<k^-8 with a margin greater
than a factor five. Hence (2) holds for this seed too. This is a
quantitative source/state input on an explicit quantum-centered
geometry. It does not say that Gamma_rest vanishes there or that
the seed solves the original retained Dirac equations.

## The actual Gamma_rest horizon forcing on the centered seed

The preceding product comparison can be converted to the source of
the exact action split Gamma_rest=Gamma_Dirac-Gamma_ref. This
additional step uses the finite-temperature reference potential; it
does not identify the finite-temperature product with Gamma_ref.

Let p0A,p0S be the complete homogeneous pressures at the exact
center B,A. Define

```
p_v=3(lambda_e-lambda0)/(8piG),
D_P=d/(8piG B^2)=p0A+p0S-2p_v,
t=r/B.
```

The product Einstein equations and the fixed outer-vacuum equation
give exactly

```
q_e-Q=3(lambda_e-lambda0)B^4-8piG B^4 p0A,
d=8piG B^2(p0A+p0S)-6(lambda_e-lambda0)B^2.
```

Insert these in the explicit DeltaV of the covariant seed. Both
reference pressures reduce to

```
pA_ref(t)=p0A t^-4+p_v(1-t^-4)
                       +D_P(-t^-4+4t^-5-3t^-6),
pS_ref(t)=-p0A t^-4+p_v(1+t^-4)
                       +D_P(t^-4-6t^-5+6t^-6).         (35)
```

In particular pA_ref(1)=p0A and pS_ref(1)=p0S exactly. The
large classical terms have canceled before taking a bound. The
second formula also follows from the exact reference Ward identity
pS_ref=pA_ref+(t/2)partial_t pA_ref.

The homogeneous bound at the retained scale gives
|p0A|<3300 q^2/B^4. Its joint trace gives
|p0A+p0S|<=q^2/(2B^4). The selected vacuum root gives
|p_v|=11lambda_e^2/(960pi^2)<q^2/B^4, so
|D_P|<3q^2/B^4. On 1<=t<=1.001, termwise differentiation
of (35) therefore bounds each reference-pressure derivative by

```
|partial_t pA_ref|,|partial_t pS_ref|
 < [4*3300+4+70*3]q^2/B^4
 <1.4e4 q^2/B^4.                                     (36)
```

For the finite-temperature homogeneous product on this same seed,
use r=w(t)/q and L=log(mu0^2 B^2t^2/q). The explicit formula for
w above gives |w'/w|<4. On this interval the complete homogeneous
coefficients and their r-logarithmic derivatives are bounded by
3400 in both pressure channels; their L derivatives are bounded
by 1. Since L'=2/t and the dilution factor is t^-4,

```
|partial_t pA_prod|,|partial_t pS_prod|
 <(4*3400+4*3400+2)q^2/B^4
 <3e4 q^2/B^4.                                       (37)
```

Both pairs agree at t=1. Integrating (36)-(37) therefore proves

```
|pA_prod-pA_ref|,|pS_prod-pS_ref|
 <1e5 h q^2/B^4 <2e6 h q^2/R^4,                       (38)
```

where (R/B)^4=(1+h)^4<1.005. Combining (38) with (1) gives
the complete actual residual source on the explicit centered seed:

```
|pA_rest,h|,|pS_rest,h|
 <(1e14+2e6)h q^2/R^4
 <1.001e14 h q^2/R^4,

max(|8piG pA_rest,h|,|8piG pS_rest,h|)/(2kappa/R)
 <5.1e-12<.01.                                       (39)
```

The tiny rounding increase still fits inside the strict margin in
(2). Equation (39) is a source VALUE on the specified regular HH
seed. It is not the derivative of that source with respect to a
metric correction, an estimate of the full corrected residual, or
a bounded inverse for the remaining equations.

## What remains

The proof controls the state and source on the fixed RN reference
and on the explicit centered covariant seed, against their respective
products with their actual physical bolt radii. Its new amplitude
response is only each explicit cap-to-H2 homotopy used above. A
subsequent unknown corrected exterior must have its own response
estimate; the existing spatial Theta' theorem cannot serve as that
parameter derivative.

In particular the covariant centered seed's independent a_star is
not set equal to its RN value by this calculation. The response of
its source to the subsequent metric correction has not been bounded
here. The cubic-jet derivative
anomaly witness, the same-regularity solver obstruction, the moving
state/domain problem and the full finite boundary charge retain
their stated force.

Focused verification outcomes are recorded separately in its receipt
and run log; the ordinary theorem does not claim an execution outcome.

Remaining to-do list: control the centered seed's exterior source away
from the bolt and the h-weighted response to its metric correction;
justify the spatial approximation or inverse; complete the corrected
exterior/state and boundary charge.
