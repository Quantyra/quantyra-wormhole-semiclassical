# Cold local renormalized source derivatives through eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original centered seed, full angular tower, paired normalization
and common physical subtraction scale. This note bounds the local term
in its moving physical normal chart, denoted U_mov. It does not include
the finite gluing term or the comparison with the actual outer domain.

For h<=s_x<=.001 and a_x=s_x[1+log(.01/s_x)], the result is

    |nabla^j(U_mov-U_H2)|/(j!100^j q^2 a_x)<1e8,
    0<=j<=8.                                           (1)

These are full physical covariant tensors, including their nonsymmetric
higher components. U_H2 is the local homogeneous term defined with the
same sphere, mass tower, cutoff and subtraction scale. It is a constant
in the observation variable. The complete cold source is not identified
with U_mov. The global auxiliary inverse estimate from
[the preceding calculation](seed-cold-weighted-inverse-jets.md) is a
separate input to the remaining gluing calculation.

## 1. Reduce the local invariant parameters before differentiating

The [cold Born calculation](stationary-reference-cold-bolt.md) proves
that the complete per-mode local remainder E_M is holomorphic in
(u0,p,pbar) on |u0-1|<=.01, |p|,|pbar|<=.01 and satisfies

    |E_M|<1.4e11 M^-4.                                 (2)

It includes all Born orders after the same local a0,a1,a2 subtraction.
The [centered transfer](stationary-centered-exterior-source.md) retains
the rational-seed coefficient guards and the physical scale mu0 B.

Fix the radial flat-cap cutoff once. For definiteness use

    chi(r)=eta(16r-3),
    eta(t)=sum_(i=0)^10 binom(21,i)t^i(1-t)^(21-i),

with eta equal to one to the left and zero to the right. As in the
original Born construction, cap exp(sigma_N)-1 by multiplication by
chi. The cutoff equals one through 3/16 and zero from 1/4. The proof
of (2) needs only its fixed radial support and values in [0,1]; its
spatial nonanalyticity does not obstruct parameter holomorphy.

Rotation acts on (p,pbar) with opposite weights. The scalar trace,
radial cutoff and spin trace are invariant. In its convergent power
series only equal powers of p and pbar remain. Hence E_M is a
holomorphic function of u0 and beta=p pbar, with beta radius 1e-4.
There is no need to differentiate a choice of square root of beta.

Write z=u-1 and delta_bar=delta/C for the thermal coefficient of the
normalized seed polynomial. Its exact form is

    Bcal(u)=z^2+b z^3+c z^4+4e z^5+e z^6-delta_bar u^3,
    b=(8x+2+4d)/C<1.334, c=(3x+1+6d)/C<.501,
    |e|<1e-24, 0<=delta_bar<=1.001h^2.                 (3)

Here x is the retained physical model parameter, not the observation
point used in a_x. Positive b,c also satisfy the displayed absolute
bounds. Continue delta_bar as an auxiliary local parameter on

    |u-1|<=r_u=.005, |delta_bar|<=r_delta=1e-4.          (4)

On the larger Goursat spatial range |z|<=.03, direct coefficient sums
give |Bcal'|<.07, |Bcal''|<2.31 and |Bcal'''|<10, including the
complex thermal variation in (4). Thus the original Picard and Born
estimates continue to hold. On (4) the sharper bounds are

    |Bcal|<.00013, |Bcal'|<.011, |Bcal''|<2.05,
    |Bcal'''|<8.1, |Bcal''''|<13.                      (5)

The physical invariant relation is beta=Bcal(u)/4. Equation (5)
keeps it strictly inside the beta disk. One can evaluate the analytic
function without a square root by setting p=.008 and pbar=beta/.008;
both arguments stay inside their .01 disks. Invariance gives the same
scalar value. Denote this one-variable composition by F_M,delta(u).
The bound (2) holds throughout (4).

At delta_bar=0,u=1,beta=0, the Goursat solution has u identically one
and Bcal''=2. Its normal metric is exactly H2. This comparison concerns
the auxiliary local term; it does not construct or vary a global HH
state at zero temperature.

## 2. Sum the complete local term in the common scheme

The exact full angular estimate is

    C sum_(n>=1) d_n/M_n^4 <4 C^3/q,
    C/q<.00352.

It follows from the telescoping angular identity in the cold proof.
It makes the series of F_M,delta locally uniformly convergent on (4).
No angular truncation or mode-by-mode choice of finite constants is
introduced.

The surviving local finite term is a fixed coefficient times

    P(u)=Bcal''(u)^2-2Bcal(u)Bcal''''(u)
                         -2Bcal'(u)Bcal'''(u).         (6)

This is R^2+2 Delta R. The coefficient has absolute value at most
1000 C^2/(480pi) in the common subtraction, including the retained
log(mu0 B). Equation (5) gives |P|<4.4. The spatially constant a0
term cancels in U_mov-U_H2. The canonical per-mode anomaly cancels
the same local R term as in the original derivation. Weyl/reference
terms belonging to the separate physical action split are not added
again here.

After removing that common constant, let F_delta(u) be the angular
sum plus (6). Equations (2),(5)-(6) imply the simple joint bound

    |F_delta(u)|/q^2
      <4(.00352)^3(1.4e11)
          +[1000/(480pi)] * (.00352)^2 * 4.4
      <25000                                             (7)

on the whole polydisk (4). This is the bound used below, rather than
a differentiated unsigned source estimate. U_mov=F_delta(u) plus the
common constant, and U_H2=F_0(1) plus that same constant.

At a real cold observation |u-1|=s/(1+s)<=s<=.001. Cauchy in u and
delta_bar, with delta_bar<=2h^2 and h<1e-30, now gives

    |U_mov-U_H2|/q^2
      <=25000[s/(.005-s)+2h^2/(.0001-2h^2)]
      <1e7 a_x.                                      (8)

The thermal difference is retained; it contributes less than
1e-20 a_x to this loose bound. The physical temperature and scale
are restored to their original fixed values after this local estimate.

## 3. Observation-parameter derivatives with the cold weight

Use exactly the observation motions of the weighted inverse note:

    L_b=B_1(b1) B_2(c_x b2), c_x=sech(chi_x),
    H_x=h cosh(chi_x), H(b)=h(L_b x)_0.

Only the analytic seed at the observation is continued here. No
global cutoff or complete physical-domain inverse is continued on
the larger complex bidisk |b1|,|b2|<=rho=.01.

The Lorentz matrix estimate gives

    |H(b)-H_x|<H_x[exp(.02)-1]<.020409s_x.

The previously proved uniformizing inverse has derivative below two
and an image disk of radius at least .17(s_x+h). Consequently

    |s(b)-s_x|<.040818s_x,
    |u(b)-u_x|<.05s_x, |u(b)-1|<.0011.               (9)

All arguments remain inside the same quarter s-disk and inside (4).
Applying the u derivative bound from (7) along the segment in (9),

    |F_delta(u(b))-F_delta(u_x)|/q^2
       <25000(.05s_x)/(.005-.0011)<4e5 a_x.           (10)

Cauchy in the b bidisk, including the Hilbert norm over all ordered
parameter indices, therefore gives

    |partial_b^j(U_mov-U_H2)|/(j!150^j q^2 a_x)<4e5,
    1<=j<=8,                                        (11)

since sqrt(2)/rho<150. The zeroth coefficient is (8). These
estimates already include every Born order, the full angular sum and
the finite local restoration. They do not include gluing.

## 4. Keep the nonzero connection of the observation chart

The H2 metric in b coordinates is

    g11=1+w^2, g12=-v_x w, g22=1,
    w=sinh(c_x b2)/c_x, v_x=tanh(chi_x),
    det g=cosh(c_x b2)^2.

The physical seed metric is exp(2sigma(b)) times this metric.
For |b1|,|b2|<=.01 the quarter-disk uniformization gives |sigma|<.02.
For example its upper estimate on s_x<=.001 is

    .0041+t/(1-t)+1e-24<.02,
    t=.001[3.045 log(10)+6.75].                       (12)

Both terms increase with s_x in the retained interval, so this is a
uniform bound, not a sample at one radius.

On the smaller bidisk |b1|,|b2|<=.005, |w|<.00501 and
|cosh(c_x b2)-1|<.000013. Thus ||g-I||<.006 and
||g^-1||<1.007. The nonzero H2 Christoffel entries, with symmetric
lower indices understood, are

    Gamma^1_11=-v_x w^2/z, Gamma^2_11=-(1+w^2)w/z,
    Gamma^1_12=w/z,        Gamma^2_12=v_x w^2/z,
    Gamma^1_22=-v_x/z,     Gamma^2_22=-v_x^2 w/z,
    z=cosh(c_x b2).

Their full Hilbert norm is below two on this bidisk. In particular
Gamma^1_22(0)=-v_x; this is not a normal chart.

The conformal correction to the connection is
delta^i_j sigma_k+delta^i_k sigma_j-g_jk g^il sigma_l.
At g=I its Hilbert norm is exactly 2|partial sigma|. The bounds on
g and g^-1 enlarge the coefficient to less than 2.1. Cauchy from
(12) gives |partial sigma|<=sqrt(2)(.02)/.005. Hence the full
physical connection has norm below twenty. A further Cauchy estimate
on the radius-.005 bidisk yields

    |partial^a Gamma|<=gamma_a=20 a!300^a, a>=0.        (13)

Put B_(0,0)=1e7 and B_(0,k)=4e5 k!150^k for 1<=k<=8. For
r+k+1<=8 define

    B_(r+1,k)=B_(r,k+1)
       +r sum_(a=0)^k binom(k,a) gamma_a B_(r,k-a).    (14)

The a=0 term must be retained. The r connection actions account for
all covariant tensor indices, with the full product rule. Equation
(14) bounds |partial^k nabla^r(U_mov-U_H2)|/(q^2 a_x) at the base
point and does not symmetrize the higher covariant tensors.

At that point g_physical=exp(2sigma_x)I. For j<=8 the conversion
from coordinate Hilbert norm to physical tensor norm costs less than
exp(8*.02)<6/5. Thus define

    J_j=(6/5) B_(j,0)/(j!100^j).

Exact rational evaluation gives the following upward-rounded bounds:

| j | J_j |
| --- | ---: |
| 0 | 1.2e7 |
| 1 | 7.2e5 |
| 2 | 1.152e6 |
| 3 | 1.918e6 |
| 4 | 3.307e6 |
| 5 | 5.909e6 |
| 6 | 1.098e7 |
| 7 | 2.127e7 |
| 8 | 4.310e7 |

Every entry is below 1e8, proving (1).

## Verification and scope of the remaining source calculation

Freeze the note, checker and review before execution. Verify the joint
polynomial and beta domains, all local and angular budgets, cold
observation shifts, explicit nonnormal connection and full tensor
recurrence. Compare independent invariant compositions and complex
contour coefficients at two precisions. A scalar b1 has a nonzero
covariant 22 Hessian at the base point; omitting Gamma(0) must fail
that control. These fixtures are not evaluations of the actual signed
charged source and do not machine-verify the continuum Born argument.

Use this bound only for U_mov with the specified moving normal cap.
The cold global-to-local gluing term must retain that same family,
including its normalizer, density and spin identification. If a later
calculation instead fixes its local cap while moving the observation,
include the fixed/moving correction; the local counterterms alone do
not remove it. The separate actual-domain error must also be bounded
through eight. No complete cold C8 source bound is asserted here.
The existing complete cold C3 and endpoint C5 contracts remain unchanged.
No ten-derivative metric, finite-metric H6 source, transfer target,
signed zero, nonlinear solution or physical stability is supplied.

Remaining to-do list: cold gluing, coordinate/density/frame transport
and actual-domain estimates through eight; quantitative ten-derivative
metric and finite-metric quantum/Schur source estimates; positive
transfer or a sharper actual transfer; cancellation, signed
compatibility, nonlinear existence; physical validity, state/boundary,
stability, matching, charge/return, opening, transport, preparation,
separate creation and audit.
