# Signed local curvature feedback at the finite seed

2026-10-22, E003/S169. Ordinary argument before focused verification.
The full stable semiclassical wormhole objective remains incomplete.

An executed part of the selected correction calculation gives an exact
cancellation in the local R2 contribution to the REFERENCE-CONSTRAINT
bolt response. Its normalized magnitude is <3.927e-18. This is not
the full action-Hessian defect pairing or the actual correction error.
The nonlocal Dirac shape response and the conversion to the fixed-disk
action rows remain required. The distinction is explicit below.

## 1. Defined response and fixed physical inputs

Use the current finite seed and its complete first reference correction
from fixed-clock-seed-bolt-feasibility.md and current-seed-full-trace-defect.md.
The Lorentzian areal metric is

    ds^2=-exp(2alpha) F dt^2+dr^2/F+r^2 dOmega^2,
    alpha0=0, R=B(1+h), F(R)=0, F'(R)=2kappa,
    D_b=F''(R)/2+kappa/R,
    R4=-12lambda+w, w=d(4t^-5-6t^-6)/B^2, t=r/B.

R is the bolt radius; R4 is scalar curvature. The symbol D_b here is
the bolt denominator, not the fourth-order trace coefficient D_R.
The action, q, mu0, c_R=1e22, prescribed kappa and outer clock are fixed.
The first correction (f,a,s) obeys

    (r f)'=-2G rho, a'=G(rho+p)/(rF), a(infinity)=0,
    f(R)=-2kappa s, f'(R)=-F''(R)s-2kappa a(R),
    sigma=G rho(R)/(kappa R), s=(kappa/D_b)[sigma-a(R)]. (1)

Here rho,p are the COMPLETE evaluated remainder source. Exact reference
constraints give delta R4=-T0 and delta Ric^t_t=-2G rho/r^2-T0/2.
No approximation to the Dirac shape derivative follows from that trace
identity. The fixed-radius reference potential has zero variation of
its mixed energy density and zero rho+p; its action-density contacts
are a separate matter, addressed in section 4.

For +c_R int R4^2 the source-side stress is -2c_R H_ab, where

    H_ab=2R4 Ric_ab-g_ab R4^2/2+2(g_ab Box-Hess_ab)R4.

This agrees with the metric f(R) equation of
[De Felice and Tsujikawa, equation (2.4)](https://arxiv.org/html/1002.4928v2#S2.SS1).
Their general equation does not supply this boundary functional or estimate.
Sphere-integrated local energy/pressure satisfy
rho_L=8pi c_R r^2 H^t_t, p_L=-8pi c_R r^2 H^r_r.

Define the local input to the next REFERENCE-CONSTRAINT response by

    rho_2,L=delta rho_L[h1],
    a_2,L'=G/r delta[(rho_L+p_L)/F][h1], a_2,L(infinity)=0,
    s_2,L=(kappa/D_b)[G rho_2,L(R)/(kappa R)-a_2,L(R)].   (2)

The variation is at fixed areal r. The denominator variation in (2)
is included; replacing it by delta(rho+p)/F alone is incorrect.
Equation (2) defines a scalar diagnostic on the actual finite seed.
It does not apply a second-order inverse to arbitrary quantum data,
and does not assert convergence of reference iteration.

## 2. Compute the local column before estimating it

Direct curvature algebra for arbitrary F,alpha gives

    H^t_t=2R4 Ric^t_t-R4^2/2
                   +2[F R4''+(F'/2+2F/r)R4'],
    (rho_L+p_L)/F=16pi c_R r^2
                          [R4''-(2R4/r+R4')alpha'].       (3)

The lapse factors cancel in the derivative term of the first line.
Differentiating the second line, using delta R4=-T0, gives

    a_2,L'=-16piG c_R[r T0''+(2R4+rR4')a'].

The known physical mass/log tail makes the following boundary terms
at infinity vanish. Integration by parts therefore gives EXACTLY

    a_2,L(R)/(16piG c_R)
      =-R T0'(R)+T0(R)-(2R4(R)+R R4'(R))a(R)-I,
    I=int_R^infinity [3R4'+rR4'']a dr.                    (4)

Only the already proved value bound for a is needed to bound I; no
total-variation bound for the complete quantum lapse is assumed.
The horizon energy row from (3), divided by 16piG c_R, is

    [G rho_2,L(R)/(kappa R)]/(16piG c_R)
      =-R Ric^t_t(R)T0(R)/kappa-2R4(R)sigma
       -2R s R4''(R)
       -[R F''(R)s/(2kappa)+R a(R)+4s]R4'(R)-R T0'(R).
                                                               (5)

In (5) use the actual moving-bolt tangent constraints (1), not f(R)=0.
Since Ric^t_t(R)=-D_b-kappa/R, subtracting (4) from (5) cancels
T0'(R) and the a(R)R4'(R) terms. The resulting signed functional is

    s_2,L/(16piG c_R)
      =R T0(R)-2R4(R)s-(2kappa R/D_b)s R4''(R)
                   -(R+3kappa/D_b)s R4'(R)+(kappa/D_b)I. (6)

This is the useful cancellation: the large separate horizon-derivative
terms never need absolute estimates. It includes the nonconstant
curvature and the actual lapse extending to the outer end.
For w=0, (6) reduces to

    s_2,L=(1-F_l)s+16piG c_R R T0(R),
    F_l=1-384piG c_R lambda.                              (7)

Equation (7) is only a specialization of (6); w is NOT set zero for
the numerical bound on the current seed. The large c_R does not cause
a large reference-bolt feedback in this particular direction.

## 3. Substitute current constants without dividing a loose weight by h

The available current bounds are

    |s|/(Bh)<3.33e-10, ||a||infinity<3.3e-10+9.182e-37,
    kappa/(D_b B h)<1.001, 1-F_l<1.179e-8,
    |d|<8.6e-19, x=lambda B^2>8e9,
    h_min=2000^-10/8<=h<=h_max=2000^-10,
    16piG c_R lambda<4.91e-10.                          (8)

The trace polynomial in current-seed-full-trace-defect.md also supplies
a sharper type of bound than its weighted supremum at the bolt:

    |T0(R)|/(lambda h)<8e-27.                           (9)

To prove (9), retain its exact product zero at t=1 and zero temperature.
Its previously evaluated coefficient sum for t partial_t is <7e-27
on t>=1. Integrating from 1 to 1+h costs at most h. The finite-temperature
polynomial has beta<=7h_max^2; its remaining constant at t=1, divided
by h_min, is <1e-60 in lambda units. The exact rational calculation
checks this last inequality with all u=1/x and d terms retained.
Thus (9) follows directly, without the spurious h^(-3/5) from dividing
the weighted trace bound by h. This is not further optimization of the
already closed scalar-response allocation; it is the endpoint quantity
needed in the signed shape/energy functional (6).

The curvature polynomial gives

    |w|<=10|d|/B^2, |r w'|<=56|d|/B^2,
    |r^2 w''|<=372|d|/B^2,
    int_R^infinity |3w'+r w''|dr<=36|d|/B^2.

Using the small positive h in (8), the three s-weighted w terms in
(6) have combined coefficient <78|d|/B^2. Hence

    |s_2,L|/(Bh)
      <(1.179e-8)(3.33e-10)
       +(4.91e-10){1.001(8e-27)
          +(8.6e-19/8e9)[78(3.33e-10)
                         +36(1.001)(3.3e-10+9.182e-37)]}
      <3.927e-18.                                      (10)

This is more than 4e8 below the outstanding 1.667e-9 cold-difference
target as a numerical comparison. It is NOT a bound on that difference.

## 4. Exact connection to, and limits of, the full defect calculation

Let E(g) denote original constraint equations in a fixed chart and let
the action gradient in the physical solver be F_eq(g)=M(g)E(g), where
M includes the equation-to-variation pairing and density. For the
off-shell seed, its actual action Hessian satisfies

    F_eq+H h1=M(E+DE h1)+(DM[h1])E.                    (11)

A coordinate return to the fixed disk adds the corresponding pullback
and chart contacts. They are not zero merely because the reference
action is stationary. Equation (2) belongs to the first constraint term
of (11); it has NOT evaluated the second term or the entire local
fixed-disk action contribution. The direct continuum inverse also acts
on the FULL defect, not on (2) alone. These limitations prevent promoting
(10) into an actual physical correction or an admission certificate.

The complete Dirac directional shape response is still nonlocal. In
fixed canonical conformal coordinates it involves the differentiated
resolvent, tadpole and bubble, common four-dimensional subtraction,
LLL/Weyl terms and actual paired-domain boundary contribution. The
first-order matter graph bounds cannot be traced over the infinite
tower to obtain it: that unrenormalized sum diverges. The trace anomaly
and the local calculation above do not determine it.

Decision: continue the SAME selected full finite-domain residual/adjoint
calculation. This turn executes (2)-(10), a signed local contribution,
and identifies the necessary conversion (11). It does not reopen the
closed feasibility audit or substitute a reference inverse for H^-1.
Next evaluate the full quantum shape column along h1, carry (11) to
the fixed-disk physical rows, and pair the complete residual with the
actual adjoint or solve the equivalent full BVP. That coupled calculation
remains unexecuted. Preserve <1.667e-9 cold difference, <=1.25014e-27
shape-to-trace source, Y0<=5e-9, Z0<=.25, Z2<=5e7 and all admission
and physical requirements. S169.2-S169.5 remain unchecked.

Remaining to-do list: complete the actual quantum shape column and
fixed-disk contacts, full signed correction with weighted admission;
nonlinear compatibility, matching, causal boundary realization, validity,
all-sector stability and every remaining stable-wormhole requirement.
