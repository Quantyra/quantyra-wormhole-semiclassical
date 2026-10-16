# Quantum reference centering at the clock and seam scales

2026-10-15, E003/S137.1. The small fractional shift of the stationary
product must be compared with the much smaller opening scale. In the
specified one-loop prescription of [the uniform reference](uniform-dirac-reference.md),
the k=2000 sphere shift exceeds 1e5 in units of the dynamical area scale
N*eta. This does not destroy the reference or the opening. It requires
consistent centering before applying an absolute remainder ledger.

The S136 identity remains exact about its explicitly classical reference.
Its illustrative weighted-jet window was conditional. The following
comparison does not turn its failed enclosure into a physical failure,
or determine the omitted finite couplings of the full theory.

## A pressure lower bound in the retained scheme

Keep exactly the previous one-free-Dirac homogeneous prescription:
q=k^4, g_R=k^-4, L4=k^5*l_p, mu0=1/b0, and zero additional finite local
coefficients. The source determinant is the complete summed angular heat
trace, not an independently renormalized LLL determinant. Write

```text
C0=sqrt(1+12*pi*k^6), xi0=(C0-1)/6,
r=C0/q, L=-log(q), epsilon=4*k^-8*(1+1/C0),
mu_open=2G*N*eta/b0^2=4/(xi0*k^7), eta=1/2000, k=2000.
```

The accepted exact reference theorem gives

```text
z_star=b_star^2/b0^2,
|z_star-1+epsilon*A(r,L)|<5e-43.
```

Its earlier pressure quadratures were numerical illustrations. Here the
predeclared test is 0.028<A<0.030 using rational inequalities and explicit
infinite tails. No tolerance for an adaptive quadrature is promoted to a
certificate. The pressure kernel and subtractions are

```text
S_q(t)=[1+2 sum_(n>=1)(1+2n/q)exp(-tn(1+n/q))]/(4pi),
hP(t,r)=1/(2pi*t)+(2r^2*t/pi) int_0^infinity
                         lambda^3 exp(-rt lambda^2)/(exp(2pi lambda)-1) d lambda,
p0=1/(4pi^2), p1=-1/(24pi^2*q),
p2=(15-3/q^2+3r^2)/(720pi^2),
I=int_0^1 [hP*S_q-p0/t^2-p1/t-p2]dt/t
  +int_1^infinity hP*S_q dt/t,
A=-(1/2)[I-p0/2-p1+p2*(gamma_E-log(q))].
```

All of these definitions are unchanged from the accepted reference.
The proof below compares this pressure with q=infinity,r=0 while keeping
the same finite scale logarithm L. This is a comparison kernel, not a
change of physical charge or renormalization prescription.

For that comparison, S_infinity=coth(t/2)/(4pi), hP_0=1/(2pi*t).
The subtracted integral is

```text
I_flat=(J_uv+J_ir)/(8pi^2),
J_uv=int_0^1 [coth(t/2)/t-2/t^2-1/6]dt/t,
J_ir=int_0^1 coth(1/(2v))dv,
-1/720 <= J_uv <= -1/720+1/60480.
```

The last inequality follows from the alternating expansion of each
positive denominator in the partial-fraction series for coth; the first
two retained terms suffice on 0<t<=1. At v=0 define the integrand as 1.
It increases on [0,1], so 64 equal left and right Riemann sums enclose
J_ir. Each coth value is 1+2/(exp(1/v)-1). Finite positive Taylor sums
for exp and a geometric bound on their remaining terms give rational
enclosures. Machin alternating sums enclose pi. Exponential comparisons
verify 30.403<log(q)<30.404; the same method and the usual one-term
harmonic-number remainder verify 0.577<gamma_E<0.578. For the latter,
H_n-log(n)-1/(2n)<gamma_E<H_n-log(n)-1/(2n)+1/(12n^2), n=100,
follows from the positive-real digamma remainder in
[DLMF 5.11(ii)](https://dlmf.nist.gov/5.11#ii).

Here is a finite bound for the difference from the actual kernel. Set
tau=1/4, r_max=31/10000 and use 0<r<r_max, q=16000000000000.
On 0<t<=tau use the third Euler--Maclaurin sphere polynomial S3 from
the accepted reference and

```text
HP3=1/(2pi*t)+r^2*t/(120pi)-r^3*t^2/(252pi),
|hP-HP3|<=r^4*t^3/(480pi), |S3|<=1/t, hP<=1/t.
```

The thermal remainder follows from |exp(-x)-1+x|<=x^2/2 and the exact
Bose moments. The sphere remainder is bounded by

```text
R_S(t)=1/69984 sum_(j=0)^3 sum_(v=0)^(3-j)
 [7! binomial(3-j,v)4^v v!/(j!(7-2j)!)] q^(-j-v)t^(5-j-v).
```

This uses pi>3 and zeta(6)<2 in the previously proved uniform remainder.
For q=infinity retain only j=v=0. Integrate the absolute coefficient
difference of HP3*S3(q) and (1/(2pi*t))*S3(infinity), after removing
their t^-2,t^-1,t^0 coefficients, divided by t. Add
int_0^tau [R_S(q)+R_S(infinity)]dt/t^2
and r_max^4*tau^2/(2880). These bound the complete UV difference without
differentiating a remainder or using a finite angular cutoff.

For t>=tau, direct comparison of each exponential gives

```text
|S_q-S_infinity| <= (1/(2pi*q))sum_(n>=1)(2n+t*n^2)exp(-nt)
                 +(1/(pi*q^2))sum_(n>=1)t*n^3 exp(-nt).
```

Since exp(-tau)<4/5, the sums of exp(-n*tau), n*exp(-n*tau), and
n^2*exp(-n*tau) are at most 4,20,180. Therefore

```text
int_tau^infinity hP*|S_q-S_infinity|dt/t
 <= [1/(2pi*tau^2)+r_max^2/(120pi)]
                       *[17/(2pi*q)+65/(pi*q^2)].
int_tau^infinity (hP-hP_0)*S_infinity dt/t
 <= 3*r_max/(16pi^2).
```

For the second inequality, S_infinity<=9/(4pi) and the positive Bose
integral gives int_tau^infinity(hP-hP_0)dt/t<=r/(12pi).
The remaining subtraction differences on [tau,1] are bounded by
3*|delta p1|+2*|delta p2|, since log(4)<2. The finite terms in A add
|delta p1|+31*|delta p2| before multiplication by 1/2.
Here |delta p1|<=1/(24pi^2*q) and
|delta p2|<=(r_max^2+q^-2)/(240pi^2).
Together these give |A-A_flat|<1e-4, using only rational coefficient sums.

## Consequence at the fixed opening scale

The pressure bracket and the exact first-order remainder imply

```text
(z_star-1)/mu_open < -1e5,
abs((b_star-b0)/b0)>2e-28.
```

Indeed epsilon/mu_open=2pi*k^5/C0. The elementary bounds
3.14<pi<22/7 and 6k^3<C0<(123/20)k^3 give the first strict inequality
with substantial margin; the 5e-43 reference remainder contributes
less than 1.4e-10 in these normalized area units. The radius difference
uses |sqrt(z_star)-1|=|z_star-1|/(sqrt(z_star)+1), with z_star<1.
The S136 classical seam delta<=5.46875e-29 is smaller than this radius
translation. Its coordinate and charge dictionary cannot simply be
identified with the same seam measured from b_star.

This refutes applying S136's illustrative bound |sin(sigma)Phi_A/B|<=1e5
to the raw homogeneous corrected product about b0: at the strip center
the constant shift alone exceeds it. It does not refute a centered
dynamical jet bound, change the source, or prove the true residual large.
The physical radius remains close to b0 in ordinary relative units.
Unknown finite coefficients or a different reference state require a
different reference solution; the stated pressure bracket does not fix them.

## Exact centering with both reference sources retained

Let b_star,a_star be any stationary product solving the specified
spherical theory, with fixed physical parameters and magnetic potential

```text
V(X)=2+6X/L4^2-2r_m^2/X, X_star=b_star^2,
T_star,ab=-V(X_star)*g_star,ab/(4G),
J_star=2/a_star^2-V'(X_star).
```

These equations follow from the orbit and angular Einstein equations;
they are not an additional assumption of zero angular source. Define a
fixed bookkeeping functional on arbitrary nearby spherical fields,

```text
Gamma_ref=int sqrt(-g) [-V(X_star)+J_star*(X-X_star)]/(4G),
Gamma_rest=Gamma-Gamma_ref,
V_eff(X)=V(X)-V(X_star)+J_star*(X-X_star).
```

Then V_eff(X_star)=0 and V_eff'(X_star)=2/a_star^2.
The reference stress is its displayed scalar density times g_ab and
its radius source is J_star. It obeys
nabla^a T_ref,ab=J_star*partial_b X/(4G), so subtraction preserves the
joint Ward identity. This is an exact reorganization of one fixed action,
not a new finite coupling or a claim that the true vacuum response is
local. The full derivative and nonlocal response remains in Gamma_rest.
At the exact stationary product its stress and radius sources vanish.

For Phi_star=(X-X_star)/(2G), g=a_star^2 exp(2omega_star)g0, the exact
S136 form now reads

```text
L0 Phi_star=T_rest+G_star,
U_star=V_eff/(4G)-Phi_star/a_star^2
      =-2G r_m^2 Phi_star^2/(b_star^4 X).
```

G_star is precisely S136's displayed kinetic/conformal/potential tensor
with b_star,a_star,Phi_star and U_star. No relation
b_star^2/a_star^2=1+6b_star^2/L4^2 is imposed on a quantum product.
The geometric estimate instead uses
mu_star=2G B/b_star^2 and vartheta_star=a_star^2 r_m^2/b_star^4.
For the stationary centered solution Phi_star=omega_star=0, both sides
vanish without subtracting two enormous absolute norms.

Retain the known 2D anomaly contact in the same allocation as the chosen
Gamma/Gamma_ref; never add an independent LLL anomaly to the already
complete four-dimensional determinant. The remaining free benchmark
stress and its dilaton shift must be transferred together as in S136.
Centering the full stationary product is more than subtracting that
single known 2D constant contact.

## The exterior normalization is still physical data

A corrected product does not determine a corrected RN exterior. For a
matched stationary metric -exp(2alpha)F dU^2+dR^2/F with alpha(infinity)=0,
near-horizon F~(R-b_star)^2/a_star^2, the leading kinematic normalization is

```text
phi_r,geom=b_star*a_star^2*exp(-alpha_star)/G,
Phi_star/rho=phi_r,geom*dt/dU.
```

Alpha_star is the limiting exterior lapse, not fixed by a product
calculation. Complete effective charges and nonlocal stress can supply
additional response, retained in the boundary defect. Substituting b_star
into the classical formula b^3/[G(1+6xi)] silently imposes both a classical
radius relation and a lapse normalization. Neither follows from the
accepted homogeneous determinant. The matched exterior stress, initial
state and moving-cut charge must still be calculated.

## Verification fixed before execution

The checker must verify the pressure bracket by exact rational arithmetic:
64 monotone cells, positive exponential Taylor sums with their geometric
tails, Machin pi bounds, explicit harmonic/log bounds and all UV/IR tails.
No floating quadrature is used for acceptance. Decimal output is a display
of rational bounds. Symbolic source, Ward, potential and normalization
identities must vanish exactly. Three independent 60-digit numerical
product/cut comparisons must agree within 1e-35; these test the identities,
not an unknown physical exterior. Preserve any failed assertion in the log.

Independent review superseded the first checker receipt: specializing the
sphere polynomial to inverse-q=1 before taking absolute coefficients can
cancel terms and does not establish a uniform majorant. The corrected
checker takes absolute multivariate coefficients first and then bounds
each inverse-q power by one. It also explicitly proves C>6k^3+1 before
using xi>k^3. The original source, receipt and log are retained in
`evidence/quantum-reference-clock-centering-superseded.json`. The targets,
64 cells and numerical tolerances are unchanged; only a fresh execution
of the corrected checker can supply the accepted receipt.

Remaining to-do list: match the centered physical exterior and state,
bound Gamma_rest's stress/radius/boundary response, prove common-drive
capture and complete the physical tail, area and S074-S079 gates.
