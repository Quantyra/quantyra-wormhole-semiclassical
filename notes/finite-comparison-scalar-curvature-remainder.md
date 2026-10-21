# Scalar curvature cancellation and local Schur absorption at finite comparisons

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary argument
before frozen verification. **Finite-background scalar-block progress;
the mixed estimate, coupled inverse and full goal remain open.**

On every original finite comparison, with t=r0/B and lambda=lambda_e,

    |R4(g_ell)+12lambda| < 3.1e-12 lambda/t.             (1)

This is a bound on the actual scalar curvature, not on the complete
semiclassical trace defect. The latter also contains epsilon Box R4.
The small value in (1) does not bound that derivative.

For the local scalar block of the exact low-frequency elimination,
the corresponding auxiliary norm is less than 2.6e-13. This replaces
an uncontrolled use of the seed's much smaller Z at a finite endpoint.
It does not assume or construct the candidate TT operator inverse.

## 1. Fixed geometry and the linear curvature equation

Keep the full theory, physical flux, prepared period, clock, state,
paired domain and finite coefficients. Use the original path of the
[classical TT proof](finite-comparison-classical-tt-form.md):

    g_ell=e^(2w)(d rho^2+a0^2 dtheta^2)
                            +r0^2 e^(2S)dOmega^2,
    y=sqrt(lambda)rho, t=r0/B, tau=(t-1)/t,
    k=exp(ell V), T=1+tau(k-1), S=log T,
    nu=tau V, w=ell u+S, 0<=ell<=1.

Let D2=partial_y^2+A0 partial_y, A0=(log a0)_y,
R0=(log r0)_y, and L4=D2+2R0 partial_y. These are fixed-seed
operators. A bare R below denotes R40/lambda, not R0.
The [invariant radius equation](finite-comparison-classical-disk-contact.md)
and its input proofs give

    |u|,|V|,|nu|<1.2e-8, |nu_y|<7e-7,
    |D2 u|<3.05e-7/t^2,
    |nu_yy+(A0+R0)nu_y|<6.3e-7,
    |L4 nu|<2.2e-6,                                (2)

where the last line uses R0<2.24 and
6.3e-7+2.24(7e-7)<2.2e-6. No derivative of a varying source is used.
Writing dotR for the derivative of R4/lambda at ell=0, the same
fixed-coordinate trace equation is

    dotR=-8piG Arest/lambda+nu t partial_t R,
    |8piG r0^2 Arest/C|<1e-20,
    R=-12+(d/x)(4t^-5-6t^-6),
    x>7.9e9, |d|<1e-12, C/x<6.1.

Consequently, with the sharper actual d/x rather than |d|<=1,

    |R+12|<1.3e-21/t^5,
    |dotR|<6.2e-20/t^2, |R|<12.01.                 (3)

The coordinate Lie term nu t partial_t R is retained: its modulus
is below 1.2e-8(56)(1e-12/7.9e9)t^-5. The trace here is the
complete seed anomaly-minus-reference allocation, already proved in
[the constraint estimate](stationary-centered-constraint-residual.md).

## 2. An exact nonlinear remainder

Warped-product scalar curvature, including the physical radial lapse,
gives the following identity without a background field equation:

    R_ell/lambda=e^(-2w)[R-2ell D2u-6L4S-6S_y^2
                                  +2(e^(2ell u)-1)/(xt^2)]. (4)

It follows either from R4=R2+2/r^2-4 Box2 r/r-2|dr|^2/r^2,
or directly from the four sectional curvatures. In particular the
term with S_y^2 has coefficient six. Its linearization is

    dotR=-2(u+nu)R-2D2u-6L4nu+4u/(xt^2).           (5)

Set n=u+nu, B=k-1-ell V k. The exact exponential radius identities are

    T_y=ell k nu_y+B tau_y,
    L4T=ell k L4nu+B L4tau+ell^2 k tau V_y^2,
    L4S+S_y^2=(L4T)/T,
    k/T-1=(k-1)/(tT), L4tau=2A0R0/t.              (6)

Thus N=R_ell/lambda-R-ell dotR is exactly

    N=R[e^(-2w)(1+2ell n)-1]
       +ell(e^(-2w)-1)dotR
       -6e^(-2w)[ell(k/T-1)L4nu+(B/T)L4tau
                                      +ell^2(k/T)tau V_y^2]
       +2e^(-2w)[e^(2ell u)-1-2ell u]/(xt^2).       (7)

The S_y^2 term has cancelled against the logarithm derivative in
(6). Keeping these terms together is essential. No division by the
thermal gap, or by r0_y at the bolt, occurs in (7).

## 3. Uniform weighted remainder budget

Use the retained bounds |w|<3e-8, 0<=tau<1. In (7), each small
exponential and inverse T is within 1.001 of one. In particular

    |B|<.501 ell^2 V^2,
    |S-ell nu|<.502 ell^2 tau(1-tau)V^2,
    |e^(-2w)(1+2ell n)-1|
          <2.01 ell^2[n^2+tau(1-tau)V^2/2],       (8)

where using .502 in place of 1/2 in the last bracket is a safe
common majorant. The second bound follows by differentiating
log[1-tau+tau exp(ell V)] twice. The third follows by writing
w=ell n+(S-ell nu), using the double zero of
exp(-2ell n)(1+2ell n)-1, then the first bound on S-ell nu.
We use that safe .502 majorant in every numerical sum below.

Put b=(t-1)V_t. The exact seed polynomial and thermal subtraction
F0<=F_ext imply

    t tau V_y^2=F0 b^2/[x(t-1)]
       <=(t-1)[1+2/t+(3+1/x)/t^2+d/(xt^4)]b^2.

For 1<=t<=4 this coefficient is below 5.064. Indeed
(t-1)(1+2/t+3/t^2) has derivative 1-1/t^2+6/t^3>0
and value 81/16 at four. The 1/x and d/x terms fit below .0015.
For t>=4 the coefficient is below 1.689t.
This is sharper than using the unweighted global reserve seven.

On the first region use |u|,|nu|,|V|<1.2e-8 and |b|<3e-7.
On t>=4 use the independently retained low outer bounds

    |u|<8e-9/t^2, |nu|,|V|<2e-9/t^2,
    |b|<1.4e-8/t^2.                                (9)

They are from the [complete defect proof](comparison-path-complete-equation-residual.md),
not extrapolations of the small-Y jet theorem. Keep the global
|L4nu|<2.2e-6 in both regions; its further sharpening is unnecessary.
Also t|L4tau|<=22 by (6).

The following bounds for t|N|/ell^2 are obtained term by term from
(7), in its displayed order. All six rows are added, not optimized
separately.

| Term | 1<=t<=4 | t>=4 |
|---|---:|---:|
| R times exponential/value remainder | 5.75e-14 | 5e-17 |
| dotR times exponential change | 4e-27 | 4e-27 |
| L4nu | 1.60e-13 | 1.66e-15 |
| L4tau | 9.57e-15 | 1.1e-18 |
| tau V_y^2 | 2.75e-12 | 3.12e-17 |
| final u exponential | 3e-25 | 3e-25 |
| Total, rounded upward | 3e-12 | 1.8e-15 |

All replacements are monotone rational majorants on their indicated
regions. For example the third row is bounded by
6.02|V||L4nu|, and the fifth by 6.02 times the two polynomial
b^2 estimates. There is no sampling over physical radii.
Equations (3) and (7) therefore give the stronger split estimate

    t|R_ell/lambda+12|
       <1.3e-21+6.2e-20 ell+3e-12 ell^2<3.1e-12.  (10)

At ell=0 retain the first term; the full physical trace does not
vanish there merely because the comparison uses a classical seed.
Regular parity extends the formulas through the positive-temperature
bolt. All quantities refer to the original finite path, independently
of any interpolation used to estimate the quantum source.

## 4. The local low-frequency block is now small

Let J=Ric-R4g/4. To control its value without a new high-jet estimate,
compare the sectionals in the physical orthonormal frame. Put

    K=a''/a, A=a'/a, Rr=r'/r, E2=Rr'+ARr.

The four sectional values are -K, -(E2-ARr+Rr^2), -ARr,
and r^-2-Rr^2. The earlier actual pullbacks and the sharp low
bounds in the classical disk proof give changes from the seed,
after division by lambda, below respectively

    3.3e-6, 1.5e-5, 7e-6, 3.6e-6.                 (11)

For example |Delta E2|/lambda<3.7e-6,
|Delta ARr|/lambda<7e-6 and |Delta Rr^2|/lambda<3.6e-6.
Each Ricci eigenvalue is a sum of three sectionals; hence its
Frobenius change is below 9e-5 lambda. The seed scalar change
is at most the sum of (3) and (10), so |Delta J|<.001lambda.
The seed |J|<=7lambda/t^4 therefore implies |J_ell|<8lambda.
Only its radial eigenvalue enters the static spherical scalar form.

Use the exact scalar block from
[the general-background theorem](general-background-newton-operator.md):

    80pi^2 L=(Tg+4lambda)(Tg-Esp)+W,
    W=div[(Zg/36+11J/9)grad]-Z/(3epsilon),
    Z=R4+12lambda, epsilon=G/(60pi),
    Esp=1/epsilon+11lambda/3, Tg=-Box_g>=0.

For radial scalars define a_s=Z/36+11J_rr/9,
U_s v=(sqrt(|a_s|)grad v,sqrt(|Z|/(3epsilon))v), with the
same sign involution as in the
[seed outgoing construction](stationary-scalar-outgoing-inverse.md).
We have |a_s|<10lambda. The actual regular scalar realization gives
||(Tg+4lambda)^-1/2||<=1/(2sqrt(lambda)) and
||grad(Tg+4lambda)^-1/2||<=1. Thus

    A22=U_s(Tg+4lambda)^-1 U_s*/(Esp+4lambda)>=0,
    ||A22|| <=[10lambda+|Z|_infinity/(12epsilon lambda)]
                                            /(Esp+4lambda)
             <10epsilon lambda+3.1e-12/12<2.6e-13.  (12)

If the actual positive TT block and its mixed form have subsequently
been realized, the positive-block algebra of that construction gives

    ||(I+Acal Jcal)^-1||
       <=(1+sqrt(b0)/2)^2/(1-b0)<1.000002,
    b0=2.6e-13.                                   (13)

This is a conditional algebraic consequence with a newly proved
local-block input. It is not an assertion that the full Acal or the
finite TT spectral square root has already been constructed.

## 5. Exact mixed row and the remaining estimate

The TT trace-variation formulas (14)-(17) of
[the TT inverse note](stationary-weighted-tt-inverse.md) were derived
without a seed field equation. They therefore apply at each finite
comparison to its genuine TT lift h and outer-fixed return Y:

    Bmix c=-(Ttrace h-div(Etr Y))/(8piG),
    Etr=(1+11epsilon lambda/3)Z+epsilon Box Z
          -11epsilon Z^2/72+(7epsilon/12)C4^2
                              +(11epsilon/6)J^2-20epsilon M_em.

Here Ttrace is the exact local TT variation of the full trace.
Since delta_h R4=-J:h, its highest input derivatives are second
order. Weyl variation is second order too. The nonlocal determinant
has no additional compact conformal trace. The divergence term
retains the density and gauge contact, and can be paired weakly
without differentiating Etr. Its transpose must use the scalar-c
return transpose and the actual volume.

Equation (10) alone does not establish a small Etr: epsilon Box Z
and the other local terms still need their actual weighted bounds.
The retained C5 geometry makes these coefficients defined, but a
numerical estimate must be proved before importing the seed's mixed
contraction constant. The same applies to differentiated Weyl/Ricci
coefficients after integration by parts. The preceding .90 TT form
bound still needs its full logarithmic operator-domain realization.

## Verification and next decision

Freeze this note and checker before execution. Pin the physical-TT,
complete-defect and seed-outgoing receipts and reconcile their exact
inputs. Check the warped scalar formula independently using sectional
curvatures with a radial lapse, linearization and the remainder algebra,
the exponential-radius identities and exact rational regional sums.
At 90 digits use the independent manufactured warped fields
a=rho(1+rho^2/7+rho^4/101), r=2+rho^2/5+rho^4/103,
u=.02rho^2/(1+rho^2), V=.03/(1+rho^2), B=lambda=1.
Compare at rho=1/5,2/3,7/4 and ell=.2,.7,1 with absolute
tolerance 1e-70. These enlarged fields test identities only; they
are not approximations to the retained physical comparison.
Check the seed trace conversion, sectional-to-Ricci budget and local
positive-block constant. Omission controls must detect the missing
gradient square and the missing radius acceleration. Finite algebra
checks do not certify the continuum source proofs, a numerical physical
solution, or an operator-domain theorem. Preserve every run and failure.

Next establish the weighted anomaly/curvature derivative and mixed
energy-dual estimates, then realize the finite TT operator in its
physical logarithmic domains and complete the coupled range argument.
Pi_cut stays fixed at the seed. Every S167 task remains unchecked.
No dynamical stability or physical wormhole claim follows.

Remaining to-do list: complete mixed coefficients and operator domains/
range; cold/higher response and nonlinear preservation; projected solve
and signed compatibility; exterior/state, charge/return, opening,
stability, transport, preparation, separate creation and final audit.
