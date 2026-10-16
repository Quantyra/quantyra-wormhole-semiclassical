# A quantitative massive trace off the RN bolt

2026-10-15, E003/S139. Ordinary fixed-reference argument, written for
review before focused verification. No numerical run is reported here.
All quantities below are per exterior of the identical pair. The two
angular signs and the rank-two orbit spinor are counted separately as
specified below; an identical second exterior does not double the source
in the first exterior's constraint.

On the retained classical RN reference at k=2000 the entire nonzero
angular tower, in the common four-dimensional subtraction, satisfies

```
|d_t [b^2 tau_hat,mass(bt)]| <= 3e14 q^2,  t=r/b>=1.001,
|b^2 tau_hat,mass(2b)|       <= 1e11 q^2.                 (1)
```

Here tau_hat,mass is the orbit trace of the complementary massive
functional after the canonical product-frame LLL allocation. It is not
the joint four-dimensional trace anomaly. The proof includes the actual
reflected boundary layer rather than imposing a bounded raw trace at
conformal infinity. Constants are uniform over the retained preparation
range on this fixed classical family. Neither (1) nor the auxiliary caps
is a statement about a corrected geometry, a metric-response norm, a
semiclassical inverse or the boundary Hamiltonian.

## Reference and exact scalar reduction

Use the common field, state, domain and finite coefficients of
[the conformal-domain note](stationary-exterior-conformal-domain.md),
[the paired modulus/phase note](stationary-paired-dirac-phase.md), and
[the boundary-layer note](stationary-massive-boundary-layer.md). In
particular the common scale is **mu0 b=1**, with zero additional finite
local coefficients in the retained representative. Write

```
xi=b^2/L4^2, C=1+6xi, q=k^4, k=2000, eta=1/2000,
xi(1+3xi)=pi k^6, b^2/G=xi k^10,
m/b=a_m=1+2xi+d/2, Q/b^2=Qhat=1+3xi,
d=D/b, r_h=b(1+h), d<=1.001 C h^2, 0<h<1e-30,
F(r)=1+r^2/L4^2-2m/r+Q/r^2,
u=b/r, B(u)=u^2 F(b/u)=xi+u^2-2a_m u^3+Qhat u^4,
ell=b/sqrt(C).
```

The actual thermal family, not merely the coarse inequality C/q<7/k,
also gives

```
xi>=k^3, C<(123/20)k^3, C/q<7/2000=.0035,
M_n^2=n(n+q)/C, M_1>18, q/C>300,
d_n=2(q+2n), n>=1.                                  (2)
```

For explicit rational guards, pi>25/8 implies xi>=k^3 here, while
pi<22/7 implies sqrt(pi/3)<128/125. Thus
C<1+(768/125)k^3<(123/20)k^3. The stronger C bound is also
retained in [the RN spectral theorem](rn-higher-mode-boundary-spectral-control.md).
M_n is the angular mass in proper orbit units; it is not a Killing
spectral gap. The Euclidean positive square has the added mass
m_n=M_n/ell=sqrt(n(n+q))/b. The d_n count both nonzero sphere-Dirac
signs. The orbit heat trace still has spin rank two.

Let y be inward proper orbit distance from infinity and set Y=y/ell.
After rescaling the outer time by sqrt(xi),

```
ghat_2=dy^2+a(y)^2 dt_E^2, a=sqrt(B/xi), f=1/a,
y=b integral_0^u B(v)^(-1/2) dv.
```

The HH period is unchanged by all comparisons. After the radial
half-density transformation, the differential expression is

```
Dtilde=gamma_y partial_y+gamma_t f partial_t,
Ltilde=-partial_y^2-f^2 partial_t^2
       -gamma_y gamma_t f' partial_t.                (3)
```

The actual complementary-adjoint Dirac square induces constant
Dirichlet/Neumann projectors at the geodesic outer boundary. Evenly
reflect a. The mixed heat trace equals the full even-double diagonal
spin trace: its image term vanishes by chirality, as proved in the
boundary-layer note. This equality has no extra factor one half.
The even coefficient is C2, with a bounded piecewise third derivative;
it is not being treated as a smooth metric through the seam.

In units ell=1, diagonalizing gamma_y gamma_t and Fourier transforming
the thermal circle gives scalar potentials k^2 f^2 +/- k f'. They are
uniformly bounded below by -L^2/4 when |f'/f|<=L. The scalar
Feynman--Kac formula consequently gives, for zero time winding,

```
k_f(tau;Y,Y)=(2pi tau)^(-1) E J,
J=f(Y) Fbar^(-1/2) exp[tau Gbar^2/(4Fbar)],
Fbar=integral_0^1 f(Y+sqrt(tau) Z_s)^2 ds,
Gbar=integral_0^1 f'(Y+sqrt(tau) Z_s) ds,
Z_s=sqrt(2) beta_s.                                  (4)
```

beta is the unit standard Brownian bridge. The outside f(Y) is the
covariant density factor and is included in every derivative bound.
For S=sup|Z_s|,

```
P(S>x)<=2exp(-x^2),
P(min_s[Y+sqrt(tau)Z_s]<0)=exp(-Y^2/tau).             (5)
```

The scalar semibounded Feynman--Kac input is
[Simon, theorem 1.1](https://arxiv.org/html/math-ph/9907022), with its
generator rescaled to Delta. All constants specific to this problem
are derived below. Thermal windings are controlled explicitly; the
zero-winding formula is not substituted for the thermal state.

## Two auxiliary caps with different ratios

All unlabelled derivatives in the cap estimates are with respect to Y.
Put z=1-u, s0=1/1000, z0=1/1001 and zc=z0 exp(-16). The extremal
polynomial factors as

```
Bext=C z^2 P(u),
P(u)=[xi(1+2u+3u^2)+u^2]/C,  1/7<=P<=1,
B=Bext-d u^3.
```

On z>=zc/2 the d term is negligible in the following explicit sense:
zc>1e-11 and d/(Cz^2)<5e-38. Thus

```
C z^2/8 <= B <= C z^2,
B<=xi+1<=C/5.                                        (6)
```

To verify the numerical exponential bound without an unevaluated
constant, e<30/11, (30/11)^4<56 and 56^4<1e7 give e^16<1e7.
Every observation r>=1.001b has z>=z0, so its proper distance to
z=zc is at least ell log(z0/zc)=16ell.

The reference curvature bounds, including the next derivative obtained
directly from the quartic, are

```
|Rhat_2|<=2ell^-2, |Rhat_2,y|<=17ell^-3,
|Rhat_2,yy|<=56ell^-4, |Rhat_2,yyy|<=115ell^-5.        (7)
```

For example the last numerator is
sqrt(B)[(3/2)B'B''''+(1/2)B''B''']; the coordinate bounds
|B'|<=5C, |B''|<=2C, |B'''|<=17C, |B''''|<=13C
give 114.5 C^(5/2)/b^5. Here primes on B mean u derivatives.

On the extremal interval Bext'^2<=4C Bext. One direct proof
differentiates Bext'^2-4C Bext and uses Bext''<=2C and the single
interior maximum of Bext. Including d changes its consequence on
z>=zc/2 by less than the margin in

```
w=f'/f: |w|<=1.01, |w'|<=2.021, |w''|<=13.           (8)
```

Indeed w'=w^2+ell^2 Rhat_2/2 and
w''=2ww'+ell^3 Rhat_2,y/2. The bound on w follows by using
|B'|<=2sqrt(C Bext)+3d in |w|=|B'|/(2sqrt(CB)).

At zc continue the actual coefficient for a distance ell/2 while
tapering w to zero. Use eta(v)=1-10v^3+15v^4-6v^5 on 0<=v<=1,
with v the transition coordinate. It has 0<=eta<=1 and, for a
transition of length ell/2, |eta'|<=4 and |eta''|<=24 in Y units.
Beyond the transition keep w=0. Reflect this construction evenly at
Y=0. It defines the **even RN cap**, on a complete cylinder, with

```
|f'/f|<=2, |f''/f|<=8, |f'''/f|<=100,
R_even=max(f)/min(f)<=1e10.                          (9)
```

For the jets, |wcap'|<=4(1.01)+2.021=6.061 and
|wcap''|<=24(1.01)+8(2.021)+13=53.408. Applying
f''/f=w'+w^2 and f'''/f=w''+3ww'+w^3 proves (9).
The transition stays at z>zc/2 because |d log z/dY|<=1 before
the cap. At its start sqrt(B)>=.99sqrt(C)zc. The taper loses less
than a factor two in a: integral|wcap|<=.505 and
e^.6<1+.6+.18/(1-.2)=1.825<2. Consequently

```
a_min >= .495 sqrt(C) zc/sqrt(xi),
a_max <= sqrt(C/5)/sqrt(xi),
R_even < .91/zc < .91(1001)e^16 < 1e10.              (10)
```

The period in these units has minimum circumference at least zc/h,
using kappa b/C<=2h. In particular it is greater than 1e19.

A second, **smooth outer continuation**, agrees with the even cap
everywhere on Y>=0, but extends the analytic RN polynomial through
Y=0 rather than reflecting it. It can be continued through Y=-1
before flattening. The following direct bounds establish this rather
than assuming a global continuation constant. On -1/2<=u<=0,

```
a_m/C<.334, Qhat/C<.501, C>6000,
B/C<.424^2 for |u|<=1/4,
B/C<.531^2 for |u|<=1/2.
```

The distance needed to reach u=-1/2 is therefore at least
(1/4)/.424+(1/4)/.531>1. On that interval the quartic gives
|w|<=1, |w'|<=3 and |w''|<=9. For instance
|B'|/C<.752, B/C>=1/7, |B''|/C<3.51 and
sqrt(B/C)|B'''|/C<6. Taper w on [-1,-1/2] with the same quintic.
Then |wcap'|<=7, |wcap''|<=57, so again the relative f jets in
(9) hold, with the sharper last estimate 79<100. Its global ratio
need only satisfy

```
R_smooth<=4e10.                                     (11)
```

The negative continuation changes the outer warping by less than e<3;
this proves (11). The thermal minimum circumference remains at least
1e18. This smooth continuation is used only in the short-time local
subtraction. The long-time bounds use R_even, not R_smooth.

The square (3) on either complete cylinder is positive, by its Dirac
factorization. Piecewise C3 caps and the C2 reflection suffice for the
form and Feynman--Kac constructions. If needed, smooth even
approximations preserve the strict margins in (9)-(11); dominated
convergence passes the displayed estimates to these coefficients.
No cap is asserted to be another physical HH exterior.

## Actual analytic strips and the short smooth remainder

The smooth continuation agrees with analytic RN throughout the real
.4-neighborhood of each non-bolt observation. Its inner cap is at
least 16 units away and its negative outer flattening begins at -1/2.
For the required uniform complex strip use the exact normalized ODE

```
u''=B'(u)/(2C).
```

On |u|<=1.25 the retained quartic coefficients give
|B'/(2C)|<3.54 and |B''/(2C)|<8. Starting on the physical real
patch, |u'|<=1. Picard iteration on a complex disk of radius .1
stays within |delta u|<=.1+3.54(.1)^2/2<.118<.25;
its Lipschitz factor is at most 8(.1)^2/2=.04. Thus the disk
is genuinely analytic. The real warping is positive and
|a'/a|<=19/8 on this patch, already following from (8) and the
negative-continuation bounds. The equation a''=-K a, with |K|<8,
then gives on the complex disk of radius .05

```
|a/a_real-1| < cosh(sqrt(8)/20)-1
                  +(19/(8sqrt(8)))sinh(sqrt(8)/20) < .13,
|a'|/a_real < sqrt(8)sinh(sqrt(8)/20)
                  +(19/8)cosh(sqrt(8)/20) < 2.81.
```

Consequently |f/f_real-1|<.15 and |f'|/f_real<4.
Cauchy's inequality on the nested disk of radius .025 gives
|f''|/f_real<=160 there. The same strips give the fourth-order
Taylor bound needed below; they are not a postulated bounded-geometry
constant.

Put zeta=sqrt(tau)<=.1. For a fixed bridge with supremum S, center a
complex zeta disk at any real zeta0 in [0,.1] with radius
rho=[40(1+S)]^-1. On safe paths zeta S<=.4 all real bridge
points lie in the unchanged analytic patch, and the disk moves them
by at most .025. With Fbar,Gbar as in (4), their complex values obey

```
Re Fbar>=.6775 Fbar_real,
|Gbar|^2/|Fbar|<=16/.6775<24,
|zeta_complex|<=.125,
|Fbar_y/Fbar|<14,
|partial_y[zeta^2 Gbar^2/(4Fbar)]|<9.                (12)
```

For example |Fbar_y|<=9.2 Fbar_real. The last inequality follows
from the upper bound
(.125^2/4)[2(4)(160)/.6775+16(9.2)/.6775^2]<9.
The prefactor in J has absolute value at most 2e^(.4S).
Including the outside density derivative gives

```
|J|<=2e^(.4S), |J_y|<=40e^(.4S).                    (13)
```

In particular the inner derivative uses at most
2(14/2+9)e^(.4S)=32e^(.4S); the density derivative adds less
than (19/8)2e^(.4S), leaving a strict margin below 40.

Here is a finite bridge-moment estimate. For
I_j=integral_0^infinity s^j exp(-s^2+.4s)ds,
I0<2, I1<.9 and I_j=.2 I_(j-1)+(j-1)I_(j-2)/2.
This gives the successive upper bounds

```
I0..I6 < (2,.9,1.18,1.136,1.9972,2.67144,5.527288).
```

Integrating the tail (5) against the derivative of
(1+S)^6e^(.4S) then gives expectation less than 600, and in
particular

```
E[(1+S)^6e^(.4S)]<3000.                            (14)
```

Cauchy's sixth-derivative estimate and the real Taylor remainder
through degree five give, on safe paths, an expected remainder for
J_y at most

```
40(40^6)(3000) zeta^6 = 4.9152e14 tau^3.             (15)
```

The value J has a smaller bound and may use the same constant.
Expected odd coefficients vanish by bridge symmetry. The even
coefficients of degrees zero, two and four are exactly the ordinary
local Dirac heat coefficients

```
k0=2, k1=-Rhat_2/6,
k2=-Rhat_2^2/120-Box_hat Rhat_2/60,
k(t;y,y)=(4pi)^-1[k0/t+k1+t k2]+remainder.           (16)
```

These are interior coefficients on the actual positive half, not a
smooth-double expansion through the reflected seam. They also follow
from the standard Dirac local coefficients in
[Vassilevich, (4.27)-(4.28)](https://arxiv.org/html/hep-th/0306138v3).

For completeness, unsafe paths do not hide an extension constant.
By (9)-(11) the full zero-winding functional satisfies

```
|J_y|<=R_smooth e^tau(4+12tau)<1e12, tau<=.01.        (17)
```

This follows from |Fbar_y/Fbar|<=4,
|Gbar|/sqrt(Fbar)<=2 and |Gbar_y|/sqrt(Fbar)<=8, with the
outside density term included. At zeta=0 the Taylor coefficients
are bounded by 40[40(1+S)]^j without an e^(.4S) factor. If
a=.4/zeta>=4 and j<=5, tail integration gives
E[(1+S)^j 1_(S>a)]<=3(1+a)^j e^-a^2. One can check this
by bounding the tail integral with logarithmic decay
2a-(j-1)/(1+a); its endpoint correction is less than one
additional endpoint term. The entire Taylor polynomial on unsafe
paths is therefore bounded by

```
120 sum_(j=0)^5 (40zeta+16)^j e^(-.16/tau)
 <4.1e8 e^(-.16/tau).
```

The actual functional contributes at most 2e12e^(-.16/tau).
Since sup_(tau<=.01) tau^-3e^(-.16/tau)<.113, the additional
expected error is below 2.31e11tau^3. Dividing by 2pi tau and
enlarging the strict bound proves, for the smooth continuation,

```
|k_smooth - (k0/t+k1+t k2)/(4pi)| <= Cs ell^-2 tau^2,
|partial_y of that difference|  <= Cs ell^-3 tau^2,
Cs=1e15, 0<tau=t/ell^2<=.01.                        (18)
```

This proves both value and gradient estimates; the value statement
has not been inferred from the gradient statement.

## Thermal windings

Poisson summation adds terms with
Z_w=beta^2 w^2/(4tau Fbar) and phase beta w Gbar/(2Fbar).
The antiperiodic sign has absolute value one. Minimum circumference
Bmin=beta/max(f)>=1e18 gives Z_w>=Bmin^2 w^2/(4tau).
Using (9), the derivatives of the exponential and phase are bounded
by 8Z_w and 72sqrt(tau Z_w), respectively. The full density/phase
derivative is therefore bounded by 20(1+Z_w) times the undifferentiated
positive functional for tau<=.01. Summing the two winding signs yields

```
|nonzero-winding heat gradient|
 <=320 R_smooth ell^-3 tau^-1 exp[-Bmin^2/(8tau)].    (19)
```

The value has a smaller bound with ell^-2 instead. For example
e^-x<=6!/x^6 and Bmin>=1e18 make (19) less than
ell^-3tau^2 on this time range. On the non-bolt observation patch
Y<sqrt(8)log(1001)<21. Thus (19) also fits the localized boundary
estimate below: its exponent exceeds Y^2/(4tau) by more than
1e34/tau. These strict margins absorb all thermal terms in the
constants in (18) and (25); no zero-temperature replacement occurs.

## The differentiated reflected boundary remainder

Let k_even be the even RN cap and k_smooth its smooth outer
continuation. Their coefficients agree on the entire positive half.
At Y=0,

```
a=1+Y^2/(2C)+c Y^3+O(Y^4), |c|<=1,
f_smooth=1+a2 Y^2+b3 Y^3+R4,
f_even  =1+a2 Y^2+b3 |Y|^3+R4,even,
|a2|<=1/2, |b3|<=1, b3=-c.
```

The actual coefficient is c=-a_m sqrt(xi)/C^(3/2). On |Y|<=.4,
the preceding analytic strips imply |f^(4)|<=24e6, hence

```
|R4|<=1e6 |Y|^4, |R4'|<=4e6 |Y|^3,
|R4''|<=12e6 Y^2,
1/3<=f<=3, |f'|<=24|Y|, |f''|<=24.                 (20)
```

Even reflection obeys the same one-sided remainder bounds. On a
bridge contained in this collar put V=sup|Y_s| and

```
DeltaF=Fbar_even-Fbar_smooth, DeltaG=Gbar_even-Gbar_smooth,
DeltaF3=2b3 integral_0^1 (|Y_s|^3-Y_s^3)ds.
```

Elementary multiplication of (20), with no heat-calculus constant,
gives

```
|DeltaF-DeltaF3|<=2e7 V^4,
|DeltaF_y-DeltaF3_y|<=1e8 V^3,
|DeltaF|<=1e7 V^3, |DeltaF_y|<=5e7 V^2,
|DeltaG|<=4e6 V^2, |DeltaG_y|<=1e7 V,
|Fbar_y|<=144V, |Gbar_y|<=24.                        (21)
```

For H(F,G)=F^-1/2 exp(tau G^2/(4F)), interpolate F and G linearly
between the two paths. The global relative bound |f'|<=2f gives
G^2<=4F on each path and, by convexity, throughout this interpolation.
Along it F>=1/9 and tau<=.01. Direct differentiation gives

```
|H|<=6, |H_F|<=30, |H_G|<=18tau,
|H_FF|<=400, |H_FG|<=300tau, |H_GG|<=30tau,
|H_F+1/2|<=30000(V^2+tau).                           (22)
```

The last bound uses |F-1|<=48V^2 and |G|<=24V.
For example H<=3e^.01<6; the exponential is bounded using G^2<=4F,
not by the weaker separate local estimates in (20).

Define R=H_even-H_smooth+DeltaF3/2. Its derivative consists of
(H_F+1/2)DeltaF_y, the error
-(DeltaF_y-DeltaF3_y)/2, and the five remaining product derivatives
of H_F DeltaF+H_G DeltaG. Substitution of (21)-(22) gives

```
|R_y|<=1e12(V^3+tau V),
|R|<=1.3e11(V^4+tau V^2).
```

For instance the two largest derivative contributions before using
V<=.4 are 1.5e12(V^4+tau V^2) and 5.76e11V^4; they fit the
displayed bound after that substitution. Multiplying by the common
outside f(Y), and including f'(Y) and f(Y)-1, gives the complete
covariant remainder bound

```
|partial_Y [f(Y)R-(f(Y)-1)DeltaF3/2]|
 <=1e13(V^3+tau V).                                (23)
```

The removed term is exactly the cubic profile already established
in the boundary-layer note, now in ell units:

```
K3(tau,Y)=c/(2pi tau^2) integral_0^tau
  E[|Y+X|^3-(Y+X)^3] ds,
Var(X)=2s(tau-s)/tau.
```

No event indicator is differentiated. First differentiate the FK
functional under the uniform bounds. Its difference vanishes
pathwise whenever the bridge stays nonnegative, because f and its
first two derivatives agree there. Only afterwards use the crossing
event to bound it. On crossing V<=2sqrt(tau)S. Formula (5) and
E S6<=12, E S4<=4, E S2<=2 imply
E(8S^3+2S)^2<=904<1552. Cauchy--Schwarz therefore bounds the
expected local expression (23) by
4e14tau^(3/2)e^[-Y^2/(2tau)].

For 0<=Y<=.1 a collar exit requires S>=.3/sqrt(tau). Use (17)
on the exited paths and the Gaussian polynomial bound for the cubic
profile. Explicitly, with a=.3/sqrt(tau),
E[S^2 1_(S>a)]<=2(a^2+1)e^-a^2 and
V<=.1+sqrt(tau)S give E[V^2 1_exit]<=.44e^-a^2.
The cubic-profile derivative on these paths is therefore at most
[3/(pi tau)](.44)e^-a^2<.5tau^-1e^-a^2. The two full
functional derivatives contribute at most 4e12e^-a^2 before
division by 2pi tau. Their combined heat-gradient contribution is
less than 1e12tau^-1e^(-.09/tau). Since

```
sup_(tau<=.01) tau^(-3/2)e^(-.085/tau)<.21,
```

it also fits the constant below. This is an estimate after a pathwise
derivative, not a differentiation of the exit event. For Y>=.1,
the global crossing estimate suffices:

```
|partial_Y(k_even-k_smooth)|
 <=1e12tau^-1 e^(-Y^2/tau),
|partial_Y K3|<=e^[-Y^2/(2tau)].                     (24)
```

For the latter inequality use the fourth Gaussian moment and
|c|<=1; its prefactor is at most 3sqrt(3)/(2pi)<1.
The ratios to sqrt(tau)e^[-Y^2/(4tau)] are bounded by 631 and 10,
respectively, because the first maximum is
2^(3/2)e^(-3/2)/Y^3 and the second is
sqrt(2)e^(-1/2)/Y. Combining the regimes and restoring units yields

```
|partial_y[k_even-k_smooth-K3]|
 <=Cb ell^-3 sqrt(tau) exp[-Y^2/(4tau)],
Cb=1e16, 0<tau<=.01.                                (25)
```

The much smaller thermal bound (19) is included. This is a direct
differentiated estimate, not differentiation of a qualitative
O(log Y) trace remainder.

## Actual geometry versus the positive cap square

This comparison avoids an unspecified global diagonal heat bound on
the actual long disk. Work on the even double. The actual compact
positive square and the complete even cap agree before both cap
transitions, with the same time circle and spin structure. A scalar
cutoff chi is one on the observation patch and vanishes before the
cap changes the coefficient. Its transition can have length ell,
with

```
|chi'|<=2/ell, |chi''|<=6/ell^2,
dist(observation,supp dchi)>=15ell.                 (26)
```

Let Rf=(Lf+m_n^2)^-1, Rc=(Lc+m_n^2)^-1. The former has norm at
most m_n^-2 by positivity. The cutoff identifies the two Hilbert
spaces only where their coefficients coincide. With Cchi=[L,chi],
the exact finite diagonal difference is

```
(Rf-Rc)(x,x)=-[Rc |chi'|^2 Rc](x,x)
             +[Rc Cchi Rf Cchi^dagger Rc](x,x).      (27)
```

Indeed chi Rf chi=Rc chi^2+Rc Cchi chi Rc+
Rc Cchi Rf Cchi^dagger Rc, while
Cchi chi=[L,chi^2]/2-|chi'|^2. The remaining commutator has zero
diagonal. This derivation is valid between the cutoff-identified
spaces and does not assert equality of the two global domains.
It defines the smooth difference even though individual two-dimensional
resolvent diagonals need a UV subtraction.

Here are explicit off-diagonal row constants for (27). Let r=j+k,
j,k in {0,1}, and D be radial separation in ell units. After Fourier
transformation, a bridge between endpoints has mode exponent
V=-k_time^2 tau Fbar +/- k_time tau Gbar. For qk=k_time sqrt(tau Fbar),

```
|V_x|<=4qk^2+8sqrt(tau)|qk|,
|V_xz|<=24qk^2+100sqrt(tau)|qk|.
```

Radial Gaussian differentiation bounds each mode row by

```
512 p_tau tau^(-r/2)(1+tau)^2(1+U)^2
 E[(1+|qk|)^4 exp(V)], U=|x-z|/sqrt(tau).
```

The physical evaluation factor a(x)^-1/2 and its derivative require
at most doubling 512. Square using Jensen and use Parseval in the
time circle. Completing the square shifts its Gaussian by at most
sqrt(tau). The polynomial after the shift is bounded by
2^11(1+tau)^4(1+|X|)^8. Its Gaussian integral and total variation
are each less than 2^18. For example
(1+|X|)^8<=2^7(1+|X|^8), and direct Gaussian moments bound the
integral by 2^11 and its total variation by 2^14. The shifted lattice
sum is at most the integral divided by its mesh plus that total
variation. Its mesh satisfies 1+mesh<=5e^tau, since the minimum
thermal circumference is at least one and
9sqrt(tau)<=4e^tau.

The radial factors obey
(1+U)^4e^(-U^2/2)<=256e^(-U^2/4), with exterior radial integral
at most 6e^[-D^2/(8tau)]. The physical Parseval normalization gives
only one ratio, f(x)/sqrt(Fbar)<=R_even. Combining these factors
gives a row-square constant less than 2^67 R_even. Thus a deliberately
larger valid row bound, including the spin sum, is

```
||1_far partial_x^j partial_z^k Kc(t;x,.)||_2
 <=1e12 sqrt(R_even) ell^(-r-1)
   tau^[-(r+1)/2](1+tau)^4
   exp[2tau-D^2/(16tau)].                           (28)
```

Laplace integration for M=M_n>=18 and D>=15 uses

```
(M^2-2)tau+D^2/(16tau)
 >=MD/5+161tau+7/tau.
```

The residual integral of tau^[-(r+1)/2](1+tau)^4e^(-161tau-7/tau)
is less than one: on (0,1) bound it by
16tau^-3/2e^(-7/tau), and on (1,infinity) by
16tau^4e^(-161tau). Consequently (28) for a resolvent row becomes

```
1e12 sqrt(R_even) ell^(1-r) exp(-MD/5).              (29)
```

The Cchi rows can conservatively be bounded by 14 times the
corresponding rows divided by ell, using (26). Differentiate (27)
at its diagonal. The direct term costs at most 8, and the middle
term costs at most 2(14)^2/M^2<1.21. We retain the bound

```
|partial_y tr(Rf-Rc)(x,x)|
 <=20(1e24) R_even ell^-1 exp(-6M),
|tr(Rf-Rc)(x,x)|
 <=20(1e24) R_even exp(-6M).                         (30)
```

Every use of the actual geometry in this comparison is through its
positive resolvent; neither a full-disk heat constant nor a temperature
dependent evaluation constant has been introduced.

The angular tail is explicitly summable. Since d_n<=6qn and
M_n>=M_1sqrt(n), monotonicity of M^2e^-6M and an elementary
integral comparison give

```
sum d_n M_n^2 e^-6M_n
 <=7q M_1^2e^-6M_1 <q 1e-36.                       (31)
```

For the integral comparison use
sum n^2e^-A sqrt(n)<=1.03e^-A for A>=108; its integral tail is
2e^-A[1/A+5/A^2+20/A^3+60/A^4+120/A^5+120/A^6].
Finally 7(18)^2e^-108<1e-36 follows already from e^5>100.
Equations (30)-(31) are sufficient for both statements in (1).

## Common angular zeta subtraction, including its value constant

Write

```
H(t)=sum d_n e^(-m_n^2 t), A(t)=sum d_n m_n^2 e^(-m_n^2 t),
zeta_h(s)=sum d_n [n(n+q)]^-s.
```

This is one angular zeta function inside the common four-dimensional
regulator, not a selection of separately renormalized angular actions.
Its expansion is

```
H(t)=2b^2/t+H0+H1 t+...,
H0=-q-1/3=zeta_h(0),
H1=(q^2/6-1/30)/b^2=-zeta_h(-1)/b^2,
Res_(s=1) zeta_h(s)=2,
Fang=FP_(s=1) zeta_h(s)=4gamma_E-2H_q.               (32)
```

The last identity follows by subtracting 4zeta_R(2s-1): their
difference at s=1 is -2H_q. The resulting finite value is at the
specified dimensionless scale mu0 b=1; it is not a new subtraction
freedom. In particular |Fang|<100 for q=2000^4.

The Lorentzian trace has the fermionic anomaly sign
-a4,orbit+FP(m^2 Green). Its orbit local anomaly is
(k0 H1+k1 H0+2b^2k2)/(4pi). In the common mass-weighted Green
regulator, the three added-back terms are exactly

```
(4pi)^-1{
 b^-2 k0 zeta_h(s-1)/s,
 k1 zeta_h(s),
 (s+1)b^2 zeta_h(s+1) k2 }.
```

The first equality uses Gamma(s)/Gamma(s+1)=1/s. The third finite
part contains **+2** in addition to Fang, from its pole; this +2
cancels the anomaly's 2. The k1 term cancels H0 k1. The full
common-scheme interior value is therefore

```
tau_hat,mass(y)
 =integral_0^infinity A(t){k(t;y,y)
              -(k0/t+k1+t k2)/(4pi)}dt
  +(4pi)^-1{(k0/b^2)[zeta_h'(-1)+zeta_h(-1)]
                        +b^2 Fang k2}.             (33)
```

Its radial derivative is

```
partial_y tau_hat,mass
 =integral_0^infinity A(t){partial_y k
                    -(k1_y+t k2_y)/(4pi)}dt
   +b^2 Fang k2_y/(4pi).                            (34)
```

The constant a0 term must be kept in (33), although it disappears
in (34). At another scale its finite part would acquire the
corresponding 2log(mu b)zeta_h(-1) term. No such scale change is
made here. Compact bulk formulas (33)-(34) do not compute boundary
contact distributions or a boundary charge.

Here is a direct bound on the constant in (33). It does not expand
q/n outside its convergence range. For Re s>1 put

```
g_s(x)=2(2x+q)[x(x+q)]^-s,
I(s)=integral_1^infinity g_s(x)dx
    =2(q+1)^(1-s)/(s-1).
```

Euler--Maclaurin through B6 gives
zeta_h=I+g_s(1)/2-sum_(j=1)^3 B_(2j)g_s^(2j-1)(1)/(2j)!
plus its periodic-B6 integral remainder. That remainder is analytic
for Re s>-2, since its integrand decays as x^(-5-2Re s).
At s=-1 set

```
P=g_-1=4x^3+6qx^2+2q^2x, L=log[x(x+q)],
partial_s g_s|_-1=-PL.
```

For q>=1 and x>=1,

```
|P^(a)|<=12[3!/(3-a)!]q^2 x^(3-a), 0<=a<=3,
|L^(j)|<=2(j-1)!x^-j, j>=1,
L(1)<32 at q=2000^4.
```

For the last numerical guard e^4>50 gives e^32>50^8>q+1.
Leibniz differentiation now yields

```
|I'(-1)|<=130q^2, |PL(1)|/2<=192q^2,
|(PL)'(1)|<=1176q^2,
|(PL)'''(1)|<=3000q^2,
|(PL)^(5)(1)|<=7056q^2,
|(PL)^(6)(x)|<=31968q^2 x^-3.                       (35)
```

For example the four sixth-derivative coefficients are
2880,10368,12960,5760. The periodic Bernoulli polynomial
B6(x)=x^6-3x^5+(5/2)x^4-x^2/2+1/42 has absolute value less
than eight. Thus the derivative remainder is at most
8(15984)/720 q^2<178q^2. The three endpoint coefficients are
1/12,1/720,1/30240, so (35) gives

```
|zeta_h'(-1)|<(130+192+98+5+1+178)q^2
             =604q^2<1000q^2.
```

At s=-1, P^(6)=0, so the same formula gives exactly
zeta_h(-1)=1/30-q^2/6. Its absolute value is at most q^2/6.
Consequently the dimensionless a0 restoration in (33) is bounded by

```
|k0[zeta_h'(-1)+zeta_h(-1)]/(4pi)|<167q^2.           (36)
```

This includes its finite common-regulator contribution; it has not
been set to zero by a separate mode normalization.

## Angular majorants and the non-bolt normalization

For E(x)=x(x+q)/b^2, H=2b^2 sum E'(n)e^-tE(n). The summand is
decreasing or has a single maximum. Its integral from zero is 1/t.
If the maximum is interior, it is at most sqrt(2)/(b sqrt(t));
the interior-maximum condition implies this is at most 2/(qt).
The sum is bounded by its integral plus its maximum. Hence, at all t,

```
H(t)<=6b^2/t, A(t)<=(2/t)H(t/2)<=24b^2/t^2.          (37)
```

The exact rational telescoping identity gives a second useful bound:

```
sum d_n/M_n^4=(2C^2/q)H_q^(2)<4C^2/q,
integral_0^infinity Atilde(tau)tau^2 dtau<8C^2/q,
Atilde=sum d_n M_n^2e^-M_n^2tau.                    (38)
```

At tau0=.01, q/C>300 gives e^(-.01q/C)<1/20. Direct geometric
sums then imply H(.01ell^2)<=q/8, and

```
integral_(.01ell^2)^infinity t A(t)dt<=ell^2 q/600.  (39)
```

The coefficients (16) obey, on the actual non-bolt patch,

```
|k1|<=ell^-2/3, |k2|<=2ell^-4,
|k1_y|<3ell^-3, |k2_y|<5ell^-5.                     (40)
```

For example |Box R|<74ell^-4 and
|(Box R)_y|<=115+2.021(17)+1.01(56)<206 in ell units,
which proves the last bounds.

The conversion to the required areal derivative is essential:

```
Jgeom=b^3 ell^-3 dy/dr=C^(3/2)u^2/sqrt(B).
Jgeom<=4C/s0 everywhere r>=1.001b,
Jgeom<=2CY^2 when Y<=1.                              (41)
```

For the second inequality, B<=C/5 gives u<=Y/sqrt(5); combine
this with B>=C(1-u)^2/8. These inequalities retain the u^2 factor
from dy/dr=u^2/sqrt(B). Dropping it would incorrectly magnify the
outer boundary layer.

The cap long-time heat value and gradient have the direct bounds

```
|k_even|<=100R_even ell^-2 tau^-1(1+tau)e^(3tau),
|partial_y k_even|
         <=100R_even ell^-3 tau^-1(1+tau)e^(3tau).    (42)
```

For the value, complete the Fourier Gaussian square and use
sum e^-X_j^2<=sqrt(pi)/mesh+2. For the gradient its bracket is
at most 22(1+tau)(1+X^2), whose integral and variation are at
most three and two. The period estimate used in (28) applies again.
These calculations include the outside density factor. In particular
neither line is inferred from the other.

Integrating (42) for tau>=.01 gives, respectively,
1400R_even q ell^-2 and 1400R_even q ell^-3. Explicitly
(1+tau)/tau<=101 and

```
integral_.01^infinity Atilde(tau)e^(3tau)dtau
 <=(300/297)e^.03 H(.01ell^2)<(300/297)e^.03 q/8,
100(101)(300/297)e^.03/8<1315<1400.                  (43)
```

This particular 1400 bound uses the shifted exponential integral,
not a cruder estimate by integral tA.

We can now sum (34), in units of q^2 for the left side of (1).

* Short smooth term: (18),(38),(41) give
  32 Cs(C/q)^3/s0<1.4e12.
* Short reflected remainder: (25),(37),(41) give, for Y<=1,
  192Cb(C/q)^2Y<=2.352e13, since
  integral_0^infinity tau^-3/2e^[-Y^2/(4tau)]dtau<4/Y.
  For Y>=1 extract e^(-12.5Y^2) before integration, obtaining
  576Cb(C/q)^2e^-12.5/s0<7.1e11. These are disjoint regimes;
  the global contribution is their maximum, bounded by 3e13.
* Whole angular cubic profile: (24),(37) give 48C/Y^2 before
  multiplying (41). The resulting normalized contribution is
  at most 96(C/q)^2 for Y<=1 and
  192(C/q)^2/(s0Y^2) for Y>=1, hence less than three.
  The full A(t) is used, not only its 2b^2/t^2 asymptotic.
* Long even cap: (41),(43) give
  5600R_even C/(s0q)<=1.96e14.
* Actual-minus-cap comparison: (30),(31),(41) give less than three.
* Long polynomial subtractions and the finite k2 term in (34):
  (32),(39),(40),(41) give less than 100.

The total is less than 2.28e14. Rounding upwards proves the first
inequality in (1) with K=3e14.

## The value anchor at r=2b

The inverse radius here is u*=1/2. Equations (6) give
Y*=sqrt(C)integral_0^(1/2)du/sqrt(B)>=sqrt(5)/2>1.
For short times decompose the kernel as actual-minus-even,
even-minus-smooth, and smooth; for long times use actual-minus-even
and even only. Thus the long term has R_even<=1e10, whereas the
short Cauchy and crossing terms may use R_smooth<=4e10.

Apply (33). The short smooth term (18),(38), with the value units
ell^-2 rather than gradient units, contributes at most

```
8Cs(C/q)^3 q^2 <3.43e8 q^2.                        (44)
```

The direct short boundary comparison uses values and needs no
differentiated cubic model. Each covariant FK functional has
|J|<=R_smooth e^tau<5e10. Their difference vanishes on paths
remaining positive, so

```
|k_even-k_smooth|
 <=1e11 ell^-2 tau^-1 exp(-1/tau), tau<=.01.
```

The thermal terms fit this bound by (19). Equations (37) and
b^2=C ell^2 give its normalized trace contribution at most

```
24(1e11)(C/q)^2 integral_0^.01 tau^-3e^(-1/tau)dtau
 =24(1e11)(C/q)^2 (101)e^-100 <1.                   (45)
```

Here e^100>1e40 already follows from e^5>100. The positive-square
value comparison (30)-(31) also contributes less than q^2.

The long even heat term (43) contributes

```
1400R_even(C/q)q^2<=4.9e10q^2.                     (46)
```

This bound would be different if the smooth continuation were used
in the long term; it is not used there. The long polynomial
subtractions in (33) contribute less than q^2: their combined bound
is 13(C/q)q^2, using (39)-(40) and
integral_(.01ell^2)^infinity A(t)/t dt<=100ell^-2 q/8.
The finite k2 restoration is also less than q^2, since its normalized
bound is 200(C/q)^2/(4pi). Finally the finite a0 term is bounded
by 167q^2 in (36). All terms in (44)-(46) and these restorations
sum to less than 5e10q^2. The larger stated anchor 1e11q^2 in
(1) therefore holds in the common scheme.

## Consequence and remaining scope

The anchor supplies information genuinely absent from a trace-gradient
norm alone. Adding a constant c/b^2 to tau changes the conserved
physical sources by delta p2=c/(2r^2), delta rho2=-c/(2r^2) and
delta J_X=-2Gc/r^4, without changing tau' or the lapse source.
This illustrates the logical blind direction of the derivative bound;
it is not permission to add a finite coefficient to the selected model.
Equation (33) and the anchor fix and bound that finite value in the
specified common prescription.

With the HH integration constant fixed, the stationary sources follow
from the trace by

```
I(r)=integral_(r_h)^r B'(R)tau(R)dR,
p2=I/(2F), rho2=p2-b^2tau/r^2,
J_X,mass=-2Gb^2tau/r^4.
```

The full radius source also contains the known four-dimensional
anomaly/Weyl allocation. The previously proved exact leading
massive/Weyl boundary cancellation and qualitative integrable mass
tail retain their scope. A numerical full geometric mass allowance
requires the corresponding weighted remainder integral; (1) alone
does not substitute for that calculation or the full boundary charge.

This spatial proof stops at r>=1.001b. The complementary region is
controlled in [the cold-bolt argument](stationary-reference-cold-bolt.md).
[The complete reference theorem](stationary-reference-higher-mode-clock.md)
combines the two regions and this anchor, and derives the fixed-reference
first-order clock and pointwise-source bounds. First and second
renormalized metric/radius variations, the appropriate
higher-order or controlled reduced gravitational inverse, recentered
full exterior/state and Hamiltonian remain separate problems. The
same-regularity solver rejected in
[the solver audit](stationary-exterior-solver-obstruction.md) is not
restored by the small fixed-reference source bounds here.

Remaining to-do list: quantify the weighted full mass allowance and
horizon recentering; establish the required metric responses and
solver; match the corrected exterior, state and boundary charge.
