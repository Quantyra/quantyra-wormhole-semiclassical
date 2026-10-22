# Direct cold bound for the product-centered conformal row

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The conformal row on the actual product-centered finite seed has a direct
source-specific inverse bound. If s solves that row with the complete
trace residual and zero independent shape input, then

    |s(t)|<2.4e-28 w(t), w(t)=[(t-1+h)/t]^(2/5),
    sup_cold max(|s|,|s_chi|/8)/a_s <1.2e-9,
    a_s=(t-1)[1+log(.01/(t-1))], h<=t-1<=.001.        (1)

This controls two actual cold admission components without converting
the small energy norm through a large generic bolt embedding. It is NOT
the coupled response: the shape input to the trace row remains required.
Raw C5, the intersected real endpoint contract, global moments and the
nonlinear ball are not certified by (1).

## 1. Fixed geometry and complete trace

Use the exact product and finite seed of
[the preceding audit](recentered-product-and-cold-residual.md). Write
b=b_star, x=lambda b^2, C=1+6x+d, t=r/b, V=1+h. Keep all theory data,
original prescribed kappa, actual paired domain and subtraction scale.

    f_ext=(t-1)^2/t^2[1+x(t^2+2t+3)+d/t^2],
    f=f_ext-delta/t, delta=V f_ext(V),
    x>8e9, |d|<8.6e-19, |d|/C<1.71e-29,
    1/(8*2000^10)<h<2000^-10, log(1/h)<80,
    C/q<.0035, q=2000^4.

Let epsilon=G/(60pi), D=96piG c_R-epsilon,
sigma=epsilon/b^2, Dbar=D/b^2, C_v=1+11epsilon lambda/3.
The full four-dimensional trace, including the retained derivative
anomaly and additional R^2 coefficient, is

    T=R+12lambda0+epsilon[7 C4^2/12+11 J^2/6
                         -20 I-11 R^2/72]-D Delta R.   (2)

For this magnetic field I=q^2/(2r^4). The shape-dependent, traceless
quantum stress is not determined by (2). No LLL truncation is made.
Put Rbar=b^2R and define dimensionless curvature invariants. Directly,

    Rbar=-12x+d r(t), r(t)=4t^-5-6t^-6,
    Cbar4^2=[t^2 f''-2t f'+2f-2]^2/(3t^4),
    JA=-f''/2-f'/t, JS=(1-f-t f')/t^2,
    Jbar^2=(JA-JS)^2,
    P(t)=7 Cbar4^2/12+11 Jbar^2/6
                           -11(Rbar^2-144x^2)/72.

Primes here are t derivatives. The exact outer vacuum relation cancels
the constant part, giving

    b^2 T=d r(t)+sigma[P(t)-10q^2/t^4]
                                  -Dbar Delta_bar Rbar.            (3)

At the extremal product center, the exact product trace gives

    d=-5sigma q^2+sigma P_*/2,
    P_*=2(6x+d)^2+(22/3)(6x+d)+22/3+22x^2,
    0<P_*<100x^2, |d|/(sigma q^2)<5.0001.              (4)

Thus the leading magnetic and radius terms cancel as FUNCTIONS:

    b^2 T=-10sigma q^2 (t-1)(t+3)/t^6
        +sigma[P(t)+P_* r(t)/2]-Dbar Delta_bar Rbar.   (5)

This is the crucial source-specific improvement. Bounding the two
leading terms separately would lose their zero at the cold center.
The actual finite-temperature mass delta is retained in P and Delta R.

## 2. A global radial barrier and source bound

Set p=2/5, a=1-h and w=(1-a/t)^p. It is smooth at the positive-temperature
bolt, positive and increasing, and tends to one at infinity. The exact
radial Laplacian satisfies

    b^2 Delta w/w
      =p a [f(p-1)a+t(t-a)f'] / [t^2(t-a)^2]
      <=3.37 x.                                      (6)

Here is a finite algebraic proof, rather than a sampled inequality.
Substitute the exact delta from the root t=V, clear the positive
denominator t^6(t-a)^2 V^3, and set t=1+h+y. The numerator of
3.37x-b^2 Delta w/w is a polynomial in y,h,x,d. Set x=8e9+X,
X>=0. For each y,X coefficient factor its lowest h power; bound
|d|<=1e-12 and every negative higher h monomial using h<=1e-30.
Every remaining coefficient is nonnegative. This proves (6) for
the larger stated parameter box; the checker records the exact
coefficient lower bounds. No limit h=0 is substituted at the bolt.

The following two bounds follow from the displayed rational seed:

    |P(t)+P_* r(t)/2|<=1e5 x^2 w(t),
    |Delta_bar Rbar|<=1e6 |d| x w(t).                  (7)

For completeness their finite coefficient proof is as follows. At
delta=0 both expressions vanish at t=1. Expand in z=1/t and write
each sum as sum c_n(z^n-z^4), using sum c_n=0. Since
|z^n-z^4|<=|n-4|(1-z)<=|n-4|w, absolute coefficient sums suffice.
For P divide the sums by x^2, using x>=1 and |d|<=1; for Delta R
divide first by d and then by x. The thermal terms are polynomials
in delta, with delta<=7x h^2. Their absolute coefficient sums after
this substitution have h^2 or higher powers, bounded by w. These
operations prove (7) on the entire exterior, including the bolt.

Also (t-1)(t+3)/t^6<=4w. The inherited guards give

    x^2/q^2<(.0035/6)^2,
    Dlambda<3e-9, Gq^2/(lambda b^4)<4e-27,
    40+1e5(.0035/6)^2+1e6(3e-9)(5.0001)<41.

Consequently the complete trace residual obeys

    |T(t)|<=41 sigma q^2 w(t)/b^2
                              <9e-28 lambda w(t).     (8)

This is a full-exterior bound on the exact local trace of the complete
Dirac theory evaluated on the seed. It is more information than a
horizon value and is not an estimate of its traceless stress.

## 3. Actual conformal Hessian with nonconstant curvature

For g_s=exp(2s)g, differentiate exp(4s)T[g_s], not just T[g_s].
The density identities are

    delta(sqrt(g)R)/sqrt(g)=2Rs-6Delta s,
    delta(sqrt(g)J^2)/sqrt(g)=-4J:Hess s,
    delta(sqrt(g)R^2)/sqrt(g)=-12R Delta s,
    delta(sqrt(g)Delta R)/sqrt(g)
                           =-2div(R grad s)-6Delta^2 s.

The C4^2 and magnetic density variations vanish. Using
div J=grad R/4, the exact scalar row is

    A s=6D Delta^2 s-div(K grad s)+(2R+48lambda0)s,
    K=[6-(2D+11epsilon/6)R]g+(22epsilon/3)J.            (9)

This includes the nonconstant R and its derivatives, as well as the
off-shell density contact. At the seed let

    t_c=(1+3x+d)/(b^2 t^4), R_c=R+12lambda=d r(t)/b^2,
    J_r=-t_c+J_c, J_c=d(5t^-5-(9/2)t^-6)/b^2,
    k_c=-(2D+11epsilon/6)R_c+(22epsilon/3)J_c.

For radial functions, (9) factors up to explicit lower terms:

    A s=6(-D Delta+C_v)(-Delta+4lambda)s+Q s,
    Q s=(22epsilon/3)div(t_c grad s)
                             -div(k_c grad s)+2R_c s.               (10)

The coefficient controls needed below are

    |R_c|<=10|d|/b^2, |r R_c'|<=56|d|/b^2,
    |J_c|<10|d|/b^2, |r J_c'|<=52|d|/b^2.

No quantum-source spatial or parameter differentiability estimate is
used in deriving this local trace row. Those estimates remain necessary
for the independent shape equation and the full nonlinear map.

## 4. Direct inverse in the weighted value space

Write L4=-Delta+4lambda, LD=-DDelta+C_v, u=L4 s,
and ||u||_w=sup |u|/w. By (6),

    L4 w>=.63lambda w,
    LD w>m w, m=1-1.1e-8.

The maximum principle on finite regular-bolt exhaustions, followed by
the positive Green limit, gives

    ||L4^-1 u||_w<=||u||_w/(.63lambda),
    ||LD^-1 v||_w<=||v||_w/m.                         (11)

Use the bounded regular solution and the zero-source end; for the actual
decaying residual one can work in the subspace whose ratio to w tends
to zero. The perturbation coefficients decay, so that subspace is
preserved. The homogeneous positive-mass maximum principle ensures
uniqueness. In the actual asymptotically AdS radial equation the source
starts at t^-4, giving the allowed t^-4 log t trace tail and the faster
massive branch. The nonconstant seed coefficients do not change the
indicial roots. The regular, decaying solution agrees with the positive
scalar form realization of (9); no extra scalar boundary phase is added.

A gradient estimate with no inverse temperature follows from the radial
equation and regularity at the bolt:

    r^2 f s'=integral_(r_+)^r y^2(4lambda s-u)dy,
    |f s'/r|<H ||u||_w w/3, H=1+4/.63.                (12)

The increasing weight w was used inside the integral. Hence

    ||div(t_c grad s)||_w< (7H/3) sup(t_c) ||u||_w,
    ||div(k_c grad s)||_w
             <=H[sup|k_c|+sup|r k_c'|/3]||u||_w,
    ||2R_c s||_w<=2sup|R_c| ||u||_w/(.63lambda).

Substitution of x>8e9, |d|<8.6e-19, Dlambda<3e-9,
epsilon lambda<2000^-10/(60pi) and the coefficient bounds gives

    ||(1/6)LD^-1 Q L4^-1||_(w->w)<1e-26.             (13)

This is an actual operator bound for this SCALAR row on the finite
seed, not an assumption that a nearby RN inverse transfers. Its
Neumann series, (8) and (11), yield

    ||u||_w<9e-28lambda/[6m(1-1e-26)],
    ||s||_w<2.4e-28.                                 (14)

Changing normalization did not remove a norm factor: (8), (11) and
(13) use the same physical weight and retain lambda explicitly.

## 5. Two cold admission components and the remaining coupling

Put s0=t-1, using s only for the conformal response elsewhere. On
h<=s0<=.001 the seed cap has |sigma_cap|<.051,
q0>=.99, q0<1.011, t_chi<1.06 sqrt(s0^2-h^2).
The integral in (12) is at most
H ||u||_w w t^2(t-V) in dimensionless radius. Together with
f>=.99 C(s0^2-h^2)/t^2, t^2<1.003 and C/b^2>5.99lambda,
this gives |s_chi|<2.4e-28 w. It uses a smooth disk coordinate,
not the singular areal-radius derivative at a moving horizon.

The ratio w/a_s decreases on this cold interval. Indeed
the logarithmic derivative is bounded by
[p-1+1/L]/s0<0, with L=1+log(.01/s0)>3.
At the bolt L_h>64 and h>1/(8*2000^10), so

    ||s||_cold,01
      :=sup max(|s|,|s_chi|/8)/a_s
      <=(2.4e-28) 2^(2/5)/[64 h_min^(3/5)]<1.2e-9.   (15)

The SAME norm and proof give the useful general estimate

    ||A^-1 F||_cold,01 < (4/3)e18 ||F||_w/lambda.     (16)

For the actual coupled equations A s+B xi=-T, (15) alone cannot
bound s: B xi is an additional signed source. If no cancellation with
T is used, the sufficient threshold for this row's 5e-9 allocation is

    ||B xi||_w/lambda <=2.8e-27,
    (4/3)e18(9e-28+2.8e-27)<5e-9.                    (17)

That number is a TARGET, not a measured or certified mixed-source
bound. The earlier relative energy-form estimate <.336 does not give
(17). The mixed row includes local charged-geometry terms; calling it
only a small quantum perturbation would be incorrect. An auxiliary
Einstein-frame elimination or a direct joint scalar/shape estimate
may preserve cancellations lost by separating these sources. Any such
change must return to the original metric and its admission norm.

This is a genuine direct response result for one necessary part of
the construction, not a completion of S169.2 or S169.3. In particular
the raw C5 and general real-endpoint admission hypotheses do not follow
from (15), and neither does a bound on the whole quantum determinant
in a neighborhood. Retain the global mass and cancelling outer moments.

## 6. Research decision and verification scope

Keep the product-centered route. The cold energy embedding failure is
not an obstacle for this actual scalar-row source: the explicit source
zero and positive radial factors close its value/first-disk-derivative
budget. Next execute the joint shape response or an exact local
auxiliary elimination, retaining all residual and Hessian contacts;
test (17) or replace it by a sharper combined bound.

The full contraction targets remain Y0<=5e-9, Z0<=.25 and Z2<=5e7
at rho=1e-8, plus every positive admission margin and exact compatibility.
No full X-norm inverse, nonlinear neighborhood, physical validity band
or all-sector retarded stability is established. The product remains
nontraversable without the separate opening construction.

The checker independently differentiates the full curvature invariants,
verifies the nonconstant-curvature Hessian, proves the polynomial bounds
by rational coefficient inequalities, and substitutes the benchmark at
90/140 digits. It does not replace the maximum-principle/domain argument,
sample the finite nonlocal stress, or certify a coupled solution.
The standard local constant-curvature f(R)-Maxwell solution described
by [Moon--Myung--Son](https://arxiv.org/abs/1101.1153) is a comparator;
it does not supply this nonconstant-seed scalar inverse or its quantum
shape coupling. The retained full trace scheme remains unchanged.

Remaining to-do list: joint scalar/shape cold response and full admission
norm, raw-jet and domain transfer, global mass/kernel moments; nonlinear
solution and exact matching; semiclassical validity, all-sector stability,
charge/return, opening, transport, preparation, separate creation and
full-result audit.
