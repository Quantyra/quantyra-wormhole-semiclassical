# Actual metric norms for the centered completion

2026-10-16, S139/S074. Ordinary proof before focused verification.

The finite comparison geometry in
[the reconstruction proof](stationary-centered-first-order-geometry.md)
has the real five-derivative metric bounds needed for the source-response
calculation. This is a bound on an actual specified endpoint, not a
regularity assumption on an unknown solution. It uses the fixed-seed
third-derivative trace estimates, the invariant curvature and radius
equations, and a fixed endpoint gauge. No sixth metric derivative or
derivative of a moving proper-coordinate map is assumed.

The metric and radius contracts are deliberately separate:

| Region | Orbit metric direction | Radius-log direction |
| --- | ---: | ---: |
| physical cold s<=.001, with weight a_x=s L_c(s) | 1e-6 | 1e-5 |
| compact non-bolt s>=.001, Y>=.05 | 1e-4 | 1e-4 |
| outer, after the correlated mass cubic | 1e-7 | 1e-7 |

Each bound is multiplied by the original endpoint parameter ell,
0<=ell<=1. The local derivative convention is j!100^j through order
five; the explicit local Weyl source needs only four radius derivatives.
All low coefficient bounds needed by the actual-domain resolvent are
kept separately from the larger middle high-derivative allowance.

## Objects and supporting source estimates

Use the potential S0, quantum radii B,A, C=B^2/A^2, prepared h,
physical q=k^4, k=2000, eta=1/2000 and fixed common scale mu0=1/b0
of [the centered seed](stationary-centered-covariant-seed.md).
Write t=r0/B, z=1/t and x=lambda_e B^2. Then

```
5.9<C/x<6.1,
gamma=Gq^2/(CB^2)<6.8e-28,
gamma_x=Gq^2/(xB^2)<4e-27,
eta/(8k^9)<=h<eta/k^9, log(1/h)<80.
```

The full trace, canonical state and paired domain are those of
[the fixed-seed source proof](stationary-centered-exterior-source.md).
The new differentiated fixed-seed estimates are proved, including
their C3 cap and actual-domain rows, in
[the canonical source-response proof](stationary-centered-canonical-response.md).
Their supporting-source statements, rather than its conditional
metric-response contract, are the inputs here:

* In the physical cold patch s<=.001 the trace difference from the frozen homogeneous
  product obeys spatial jets through three bounded by
  1e20 q^2 a_x j!100^j, a_x=s L_c(s),
  L_c(s)=1+log(.01/s).
* On the compact non-bolt cover s>=.001, Y>=.05 the corresponding unweighted
  fixed-seed trace jets are bounded by K_m q^2 j!100^j,
  K_m=1e26. This is the larger actual finite row estimate; it is
  not replaced by the outer constant 1e20.
* In the outer collar, with proper distance Y in A units,
  D=Y partial_Y and L(Y)=1+log(1/Y),
  |D^j[Theta-cC^2/(5piY)]|<=1e20 q^2 j!100^j L(Y),
  j<=3 and 0<Y<=.1. The coefficient c and the common local
  restoration are retained together. Ordinary derivatives through
  three cost a factor below two.

For clarity, the cold supporting proof uses a NEW width-.01 C3 cap,
not the older width-one C2 cap. Its seam third jet is below three;
the ordinary global conformal-center jets have reserves 5,6e4,5e6,
and centered weighted reserves 500a_x exp(d),7e6a_x exp(d),
1e9a_x exp(d). On the WHOLE physical cold patch, retain the
nonlinear normalizer terms with a_x<=.00331: the global second
map arrays are (2e5,7e6,6e8)a_x and the third arrays
(3e9,2e11,1e13)a_x. The uniform first/second inverse insertions
are below 2e4 and 1e9, respectively. Ordered products with one
weighted factor give M||W S_3 W||<=1e18a_x. No O(a_x^2)
term is discarded by a bound valid only near the bolt. The three-times
differentiated gluing trace is below 1e26a_x exp(-M/16).
With sum d_n exp(-M_n/16)<4q and division by 3!100^3,
its physical normalized contribution is below 2.4e17q^2a_x.
The analytic Goursat/local-UV Cauchy term and common restoration
fit inside 1e20. The C3 cap is not silently assigned the old C2
cap's energy-twenty bound.

The middle supporting proof uses normal charts R0=1/32, a gluing
cutoff one through R0/2 and zero by 3R0/4, slope below 256,
and includes BOTH commutators. Its separated point-row estimate is
3e4 j!(2000M)^j exp(-MD/2), D>=1/65 in proper distance. The small evaluation and
nested energy balls occupy less than .19/M<=.01056, so they fit
even at M=18; the old 2/M evaluation ball would not fit.
The third Leibniz sum, two commutators and rank two cost below
2.3e25 M^3 before normalization. Division by 3!100^3 and
the angular estimate sum d_n M_n^3 exp(-M_n/65)
<=768 C 65^5<8.92e11C give less than 4.4e25q^2.
Local UV and finite terms fit the reserve K_m=1e26. Thus this
supporting fact is a complete fixed-reference spatial estimate,
not an inference from a first derivative or from chart existence.

The sharper VALUE data remain important. They are

```
|Theta-L0 t|<=4e15q^2(1+log t), L0=-a_seed x/(5pi),
|I|<=4e15xq^2,
|B^2 rho_rest|<=1e16q^2(1+log t)/t^2,
|B^2(p-rho)_rest|<=1e16q^2(1+log t)/t^2, t>=4.
```

These are complete massive-plus-Weyl/reference cancellations, not
separate divergent mass integrals. They will control low outer jets
without replacing them by the much coarser 1e20 estimate.

## Sharper real low jets

Let f,a,s_R be the frozen Eulerian linear correction, and put
fbar=f+s_R F0', H=a+fbar/F0 and J=integral_r^infinity H/F0.
The actual common-disk conformal and radius tangents are

```
xi=s_R-F0 J, u=a+fbar/(2F0)-s_R/r-(F0'/2-F0/r)J,
V=xi/(r0-B), nu=xi/r0, b=(r0-B)V_r.
```

The frozen bounds include |H|<1.08e-9, |a|<3.3e-10,
|fbar/F0|<7.5e-10, |V|<1.2e-8 and |b|<3e-7.
There is a useful refinement, requiring no new source derivative.
For y=r-R>=Bh put c=kappa B/(xh), Z=y/(Bh). Convexity gives

```
F0' integral_r^infinity dr/F0
 <=[1+12.1Z/(2c)]log[1+2c/(.95Z)]<16.
```

Indeed the first logarithm is below log(1+14.74) and the
Z log term costs at most 12.1/.95. Also
(F0/r) integral dr/F0<10. Since the absolute DIFFERENCE
|F0'/2-F0/r| times that integral is bounded by the maximum
of the two positive terms, not their sum, the displayed u formula
gives |u|<1.2e-8. The regular inner J bound gives |u|<5e-9
on y<=Bh.

The common physical disk has the two-sided conformal guard

```
.39<exp(sigma0)<7.31, exp(2sigma0)<60.                 (1)
```

In the following conformal comparison reset y=t-v to be dimensionless;
it is not the earlier physical displacement r-R.
The upper bound follows by the convex comparison in
[the E0 residual proof](stationary-centered-constraint-residual.md).
For the lower bound put gamma0=kappa B/C and
w=h+C gamma0 integral_t^infinity dt/F0. The upper quadratic
bound on F0 gives
C gamma0 integral dt/F0>=.5log[1+2C gamma0/(6.05xy)].
In exp(sigma0)=sqrt(F0/C)sinh(w)/(gamma0 t), use
sinh(h+z)>=sinh z+h and the lower quadratic bound for F0.
The two contributions are respectively at least
sqrt(.95C/x)/(6.05t)>.394/t and
(h/gamma0)sqrt(.95x/C)y/t>.394y/t.
Their sum exceeds .39 because 1+y=t-h.

The invariant curvature calculation gives the sharper global
|A^2 delta Rhat|<5.2e-8 and |A^2 Rhat0|<=2. Therefore

```
(sinh chi u_chi)'=sinh chi exp(2sigma0)
                   [-A^2 Rhat0 u-A^2 delta Rhat/2],
|u_chi|<3e-6, |du|proper<7.7e-6.                    (2)
```

The proper derivative here is in A units; (1) performs the
conversion. These actual small low jets are distinct from the
larger high C5 constants used in the middle.

For the scalar radius,

```
nu=tau V, tau=(r0-B)/r0,
nu_rho=sqrt(F0/C)[b+(1-tau)V].
```

On Y>=.05, sqrt(F0/C)<20.1, giving |nu_rho|<6.3e-6.
At Y=.05 the strict outer guard gives
sqrt(F0/C)=z_Y/z<1.002/Y=20.04; monotonic F0 makes it
smaller farther inward. This uses the proper Y relation, rather
than the looser separate upper bound on t.
Also |nu|<1.2e-8. Local normal-coordinate changes below preserve
the independent low relative coefficient bounds C0<1e-6,
C1<1e-4 and C2<.01. The response's resolvent domain control
uses these small numbers, not the middle high-jet allowance.

## Invariant two-derivative transfer

The product-frame curvature and physical scalar trace are the
appropriate variables. They avoid raw f''' or a'' at the cold bolt.
In sphere-integrated source conventions,

```
A^2 delta Rhat_E=[-8piG r^2 Arest-6f+6G(p-rho)]/C,
A^2 delta Rhat_common=A^2 delta Rhat_E
                                  +nu r partial_r Rbar0/C,
delta R4_E=-8piG Arest,
delta R4_common=-8piG Arest+nu r partial_r R40.         (3)
```

Here Rbar0=B^2 Rhat0, and Arest is the exact local anomaly minus
the fixed reference trace. It is an explicitly rational seed
function; it is not another nonlocal unknown source.
For phi0=log(r0/B) the exact scalar identity gives

```
Delta_hat nu+2 grad(phi0).grad(nu)
 =(delta Rhat_common-exp(2phi0)delta R4_common)/6
   +2u[Delta_hat phi0+|grad phi0|^2]
   -(nu/3)exp(2phi0)R40.                              (4)
```

In (4) both curvature variations are in the same common physical
coordinate. The Lie term in delta R4 is retained: for the centered
seed B^2R40=-12x+4d/t^5-6d/t^6 its derivative is nonzero when
d is nonzero. Its explicit tiny d-dependent contribution fits the
existing local-force reserve. In a proper-coordinate use, replace
nu by the proper-coordinate radius displacement in both Lie terms.
The conformal equation is
Delta_hat u=-Rhat0 u-delta Rhat/2.
Both are second-order equations in the invariant variables.
Their radial versions are regular Volterra equations at the bolt.
For example u_rho=N^-1 integral_0^rho N F, where N is the
orbit warping factor and F is its curvature forcing.
No division of unrelated raw horizon jets is made.

Three source derivatives suffice for five metric derivatives.
Away from the bolt differentiate the radial equations three times.
At the bolt use the integral equation before differentiating.
For the leading N=rho term, the potentially fourth source
derivative integrates away exactly:

```
u^(5)(rho)=F'''(rho)-integral_0^1 t^4 F'''(t rho)dt.
```

The analytic N/rho factors add only the known seed-coefficient
products through three derivatives. This argument also applies to
the coupled radius equation. It is not an application of a generic
C3-to-C5 elliptic sup estimate, which would require care about
the modulus of continuity. Radial symmetry supplies the explicit
one-dimensional inverse here.

## Cold normal charts

The physical cold contract is s=r/B-1<=.001. The larger .01
seam below is an AUXILIARY barrier/cap seam, not an extension of
that response region. On the interval up to .01 the frozen barrier
can be strengthened directly.
Its seam data obey |u|<3e-9 and W_seam>=.01, whereas its
forcing requires only .56e-7/2.7<2.1e-8 times W. Thus

```
|u|<=3.1e-7 W,
W=h[(Dc+1)cosh chi-chi sinh chi],
Dc=arccosh(.01/h).                                   (5)
```

The bound retains the fixed baseline transition from the cold
auxiliary chart to the full physical disk. It does not identify
the auxiliary cap with the physical exterior.
To control low radial derivatives integrate the same equation.
The exact moment is

```
integral_0^chi sinh y W(y)dy/[h sinh chi]
 =.5[(Dc+1.5)sinh chi-chi cosh chi]
                           +chi^2/[4sinh chi].
```

With |V_pot|<3 this bounds |du| and Hess u by generous
1e-6 s(L_c+2) and 3e-6 s(L_c+2), respectively. Proper-frame
conversion costs less than 1.2 on this inner patch.
The source jets in (3) through three have the same a_x ideal:
subtract the exact product value first, use the differentiated
canonical Ward inverse for rho,p, and the local anomaly for Arest.
The fixed-seed source estimate 1e20 couples with gamma; all
coefficient products in (3)-(4) fit below 1e-4 a_x in normalized
force jets. Two integrations then divide the new high jets by
(j+1)(j+2)100^2. These high metric jets are below 1e-8 a_x;
they are smaller than the low barrier contribution in (5).

Normalize each orbit chart by a real quadratic conformal map
removing its constant and linear metric jets. On the cold normal
disk of radius 1/4 the value Taylor remainder costs at most
(1/32) times the Hessian bound. The baseline analytic normal
metric has coefficient reserve .026 times 4^j; therefore its
Lie/normalization terms, after division by 100^j, are smaller
than the remaining margin. The normal metric direction satisfies

```
|partial^j delta sigma_N|<=1e-6 a_x j!100^j, j<=5.     (6)
```

The one-center mixed jets through four use total derivative order
at most five and spend one additional factor 100, as in the
response contract; no sixth metric derivative is introduced.
The radius direction starts at |nu|<=1.2e-8s. Coordinate
normalization multiplies the baseline |dphi0|=O(s) by the
small chart displacement, so its value gains another s factor.
Equation (4) controls its high jets with the same two-order gain.
A reserve 1e-5 a_x j!100^j is valid through five.

## Compact non-bolt high jets

Use the fixed normal charts R0=1/32 on r>=1.001B and Y>=.05.
Their full physical disks fit before the boundary; no interval
between Y=.05 and the cold patch is omitted. The supporting
source proof uses this same cover.

It is useful to track old low jets separately from K_m=1e26.
Normalize the canonical quantities explicitly:

```
theta=Theta/q^2,
P=I/(2calB q^2)=B^2 p_mass/(q^2 z^2),
R=B^2 rho_mass/q^2, R/z^2=P-theta.
```

With proper radial coordinate rho,

```
P_rho+(log calB)_rho P=.5(log calB)_rho theta.          (7)
```

The compact polynomial calB is unimodal, and its total variation
from the bolt to the observation point is at most three times
its current value. Thus |P|<=1.5 sup|theta|. The normalized
proper jets of the logarithmic coefficient in (7) are bounded
by j!100^(j+1) on this cover. Differentiating (7) gives
successive reserves 2K_m,4K_m,8K_m for P derivatives one,
two and three. Consequently a reserve 16K_m bounds the
weighted R/z^2 jets. Explicit local terms are negligible
relative to K_m, and are still added, not independently discarded.

For the high curvature derivatives a crucial weight prevents an
artificial large-t loss:

```
f_rho/C=-L_phi f/C-2gamma L_phi z^2(R/z^2),
L_phi=sqrt(F0/C), L_phi z^2<.5.
```

The genuinely new third f derivative, sourced by the second
R derivative, therefore has normalized coefficient below

```
[2gamma(.5)(2!100^2)16K_m]/[3!100^3]
                         <.06gamma K_m.             (8)
```

The other product terms contain only the older, much sharper
rho value/first derivative. They are not inflated to K_m.
The direct trace term in (3) costs less than 6.6gamma K_m;
the new f term less than .36gamma K_m. The Lie term through
three uses nu through three, hence only old source value/first
jets by (4). The local anomaly is explicit. All of them fit in

```
high curvature force coefficient <10gamma K_m<.68.    (9)
```

For reference, the old low-jet arithmetic uses t<50.1,
L_phi<20.1, |nu|<1.2e-8, |nu_rho|<6.3e-6 and (2).
The background identities

```
Delta phi0+|dphi0|^2=t F0_t/C,
|t F0_t/C|<3t^2,
|partial_rho(t F0_t/C)|<10t^3,
|A^2 exp(2phi0)R40|<3t^2
```

bound the first three metric jets without using K_m. After
normalization the u jets two and three are below 1e-7 and
the radius jets below 1e-6. For example the drift coefficient
in the u equation can be bounded by 100 and its derivative
by 100^2: u''<8e-4 and u'''<.2 already suffice after
division by 2!100^2 and 3!100^3. Equation (4) similarly
gives nu''<1e-3, nu'''<1; these very loose numbers suffice.

For the new fourth and fifth u derivatives the exact ratios are

```
u4_norm <=.68/(24*100^4)*(2!*100^2)/2 <2.84e-6,
u5_norm <=.68/(120*100^5)*(3!*100^3)/2 <1.7e-6.        (10)
```

Lower coefficient products use the preceding old low jets;
the fifth-order drift feedback is below .1 and is absorbed.
The radial-to-normal-coordinate finite jet conversion has
factor below two: its leading orthonormal coefficient is one,
the baseline curvature/normal-coordinate derivatives are the
bounded seed jets, and their normalized lower products are
dominated by sum_(j<=5)(10/100)^j<1.12. A factor two also
reserves the quadratic normalizer. Thus all middle orbit
normal jets fit 1e-4 j!100^j, with more than a factor ten
over the high-jet values in (10).

For the radius, (4) divides the new curvature forcing by six.
Its remaining terms involve only u,nu through three at this
stage. Their old low bounds and the same finite normal-coordinate
conversion fit 1e-4 j!100^j. The larger radius contract is kept
separate from the orbit metric contract in the local-response
calculation. It is not an extra canonical Dirac response.

## Outer proper gauge and its correlated cubic

Define the ACTUAL baseline proper coordinate by

```
calB(z)=z^2 F0(1/z),
Y=sqrt(C) integral_0^z calB(w)^(-1/2)dw.
```

On |Y|<=.15 its analytic initial-value equation z_Y=sqrt(calB/C)
has strict guards

```
|calB/x-1|<.0011,
|z/Y-sqrt(x/C)|<.001sqrt(x/C),
|z/Y|<.412, |Y/z|<3, |(calB_Y/calB)/Y|<1.
```

For example bootstrap |z|<=.063 in the rational compact
polynomial; its nonconstant x-normalized coefficient sum there
is below .0011, and integrate z_Y-sqrt(x/C). On every
|Y-Y0|<=Y0/2 Cauchy bounds a bounded analytic multiplier's
normalized jet sum by sum(.02)^j<1.021. All products below
retain this factor for each coefficient; it is not treated as
one after multiplying several coefficients.

Let j=f/F0 and J_p(Y)=integral_0^Y j. Linear proper-gauge
directions of inverse orbit warping and radius log are EXACTLY

```
w=-z_Y J_p/(2z),
v=-a-j/2-calB_Y J_p/(4calB).                         (11)
```

They follow from deltaY=-J_p/2 and fixed-Y delta z=z_YJ_p/2.
The invariant equation for the first direction is

```
v''+(calB_Y/calB)v'=.5 A^2 delta Rhat_proper,
A^2 delta Rhat_proper
 =[-8piG r^2 Arest-6f+6G(p-rho)]/C
                                  -z w Rbar0_z/C.   (12)
```

Its signs include the actual proper-coordinate displacement.
For a mass cubic j=j3Y^3, (11) gives

```
w3=-j3/8, v3=-j3/2, v3=4w3.                        (13)
```

These are proper-gauge coefficients; they are not confused with
the different common-disk u3,V3 coefficients.

Here is a finite source-to-metric ledger. Use
N_(m,J)(F)=max_(j<=J) sup |F^(j)|/[j!100^j Y^(m-j)L].
The outer source theorem and explicit local allocations give

```
N_(2,3)(B^2rho/q^2)<2K,
N_(2,3)(B^2(p-rho)/q^2)<2K, K=1e20.                 (14)
```

For an explicit count, replace the c/Y principal term by L0/z;
the difference is analytic and below 1e8 in this norm. Thus
the remaining trace has N_(0,3)<3K. Its z^2 multiple costs
less than .63K, including the power/product jet factors.
The constant I0, bounded by 4e15xq^2, contributes below
4.2e14. For I_R=integral calB_Y Rtrace,
|calB_Y/(CY)|<1 gives N_(2,3)(I_R)<=2.5CKq^2;
integration uses integral_0^Y yL dy=(Y^2/2)(L+.5).
Multiplication by z^2/(2calB) gives N_(4,3)<1.6K,
hence N_(2,3)<.016K on Y<=.1. The principal I and explicit
Weyl/reference/LLL log polynomials contribute below 1e8.
This proves (14) with substantial slack.

Split the exact analytic mass mode BEFORE subtracting its cubic.
Write f=-2xz[m_infinity+mu(Y)], |m_infinity|<4e-11.
Then

```
mu_Y=-gamma_x z_Y (B^2rho/q^2)/z^2,
|mu^(j+1)|<=8gamma_x K j!100^j Y^-j L, j<=3,
|mu|<7e-11YL.                                      (15)
```

The last sharper value uses the sharper source VALUE data, not
K. Equation (11), or (12) for v, now gives the following ledger:
f_tail/C has N_(2,3)<.013gamma_xK; the proper curvature tail
has N_(2,3)<3gamma_xK<1.2e-6; |zRbar0_z/C|<3Y makes
the radius-Lie term lower order. The local trace polynomial
fits by many orders. The high v jets two through five are
below 6e-11 after their two-order factorial/100 gain and
the harmless first-order coefficient feedback. The low v
jets zero and one are below 1e-9 using the sharp value in
(15). The integral in w provides its second integration;
its fifth derivative uses mu^(4), hence only source derivative
three. Its high and low remainders are below 1e-9 as well.

The exact constant-mass mode has analytic higher coefficients.
Its remainder after its cubic is below 1e-9 in N_(4,5),
because calB/x-1 has no linear term and z/Y is analytic with
the strict guards above. Combining these bounds yields

```
N_(4,5)(v-v3Y^3)<1e-8,
N_(4,5)(w-w3Y^3)<1e-8.                              (16)
```

The 100^-2 gain is asserted only for the HIGH metric jets.
The low remainder coefficients in (16) use the sharper value
tails. A uniform 100^-2 gain for value norms would be false.

## Finite endpoints and fixed response paths

Equation (16) initially describes the tangent. The original finite
completion is sigma_ell=sigma0+ell u and
r_ell-B=(r0-B)exp(ell V). We need a finite endpoint estimate,
not a hidden sixth derivative from a moving gauge.

At the boundary let U(Y)=integral_0^Y u. The proper-to-common
linear identities give
U'=-v+(calB_Y/(2calB))U. The integral operator divided by Y
has N_(4,5) norm below .24; its higher recurrence is
T_j=[u_(j-1)-jT_(j-1)]/Y. The coefficient feedback is below
.0024, so the common conformal remainder is below
1.003 times (16); reserve two. The common radius-log remainder
is below 1.25 times (16); again reserve two.

For EACH fixed ell define the endpoint proper map

```
Phi_ell(Y)=integral_0^Y exp(ell u(Z))dZ.
```

It is weighted conormal C6 on the OPEN collar because u has
weighted conormal C5 jets. It is not asserted to be classical
C6 at Y=0: the allowed Y^4 logY term in u integrates to
Y^5 logY in Phi. Phi has a C4 boundary extension, with
Phi(0)=0 and Phi'(0)=1. No ell derivative of Phi is taken
in a spatial C6 norm. Its relative jets mean explicitly

```
max_(j<=6) sup |Y^j partial_Y^j[(Phi_ell-Y)/Y]|
                                      /(j!100^j)<1e-10.
```

The mass cubic is below
1e-11 and 2e-8 max(Y^4L)<6.7e-12 for the remaining value.
Finite inverse/composition formulas through order six have
coefficient sum below (6!)^2 2^6<3.4e7. Thus their
nonidentity terms enlarge coefficients by less than .01;
reserve two. The baseline radius pole is handled as
-log[Phi_ell(Y)/Y], not by a Taylor series at a singular logY.
The original log-gap radius completion adds only quadratic
terms of order Y^6 or higher, below .01ell times (16).
Consequently the ACTUAL endpoint proper chords satisfy

```
N_(4,5)(v_ell-v3,ell Y^3)<=ell 1e-7,
N_(4,5)(w_ell-w3,ell Y^3)<=ell 1e-7,
v3,ell=4w3,ell.                                    (17)
```

This proof does not differentiate the original path in moving
proper gauge. Doing that would produce xi u^(6) in its fifth
spatial derivative, which is not supplied by the present inputs.

One fixed global identification can be made explicit. Let eta=1
for Y<=.05 and zero for Y>=.1. On the intervening interval
use P(s)=1-I_s(7,7), so P'=-12012s^6(1-s)^6. On a
complex quarter-neighborhood |s(1-s)|<=.3125, |P'|<12
and |P|<4. Hence |P^(j)|<=4j!4^j; the width-.05 cutoff
has |eta^(j)|<=4j!80^j through six. Put

```
Psi_ell(Y)=Y+eta(Y)[Phi_ell(Y)-Y].
```

Its weighted relative jets are below 5e-10 and Psi'>0. It is
identity in the cold region. Away from Y=0 its C6 regularity
suffices for C5 metric coefficients. At Y=0 use the displayed
conormal calculus and its C4 extension; the resulting metric
preserves the paired boundary domain's C2 data, not a fictitious
bounded C6 closed-collar ball. The outer logarithms and exact
endpoint identification are retained. Pull the endpoint back once, then interpolate the
static diagonal metric logarithms and log radius in this ONE
fixed manifold. In the outer collar use the proper log chords
(17); in the cold region the orbit conformal path is the original
one. The radius log interpolation need not equal the original
log-gap path at intermediate parameters. It has the same
endpoints and retains the positive gap and monotone radius.

Both endpoints share the bolt period, outer clock, geodesic
boundary, pole residue and correlated cubic. These data are
preserved by this fixed interpolation. The transition lies away
from the bolt; its coefficient C1 change is below 1e-6.
The actual low coefficient guards stated after (2) therefore
remain valid. The complete quantum source is a functional of
the endpoints, so this auxiliary response interpolation changes
neither the desired source difference nor the E0 residual of
the ORIGINAL completed endpoint.

## Conclusion and limits

The completed endpoint has the actual cold, middle and outer
norms in the table, with independent small low coefficients.
For every original ell the fixed auxiliary chord has norm
bounded by ell times those constants. Therefore a first-response
bound applied on that chord gives an ell source difference;
the explicit ell multiplying Gamma_rest in an order-reduction
ledger supplies a second ell. This uniformity does not require
C5 differentiability of the original moving-gauge path.

These are the metric/state-map INPUT bounds. Combining the
reviewed response theorem, local Weyl/reference variation,
coordinate/radius prefactors and the E0 residual still requires
the finite integration ledger. No contraction of the full
higher-derivative retained equations, no exact stationary
solution and no complete Hamiltonian charge is inferred.

Remaining to-do list: independently review the complete ordinary
metric-transfer and supporting source-jet proofs; combine source
response and E0 residual in the physical inverse norms; perform
focused verification after the ordinary combination is accepted;
justify the approximation or a higher-order inverse and complete
the physical exterior/state/charge match.
