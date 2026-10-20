# Strong outer bounds for the global term and constant contact

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Two actual quantum
components controlled; the full curved outer inverse remains open.**

The bounded actual-domain gluing term and the complete constant
contact have strong weighted bounds at the selected fast weight.
The same bounds include the mass cubic and input from the rest of
the disk. They require no first-contact estimate and do not infer
strong regularity from the existing quadratic-form inequality.

## 1. Exact components and conventions

Use the seed, common finite subtraction, complete angular tower,
physical period, paired domain and canonical orbit of
[the curved quantum proof](stationary-curved-quantum-tt.md).
Write

    gbar=dY^2+a(Y)^2 dt^2, h=log(a/a(0)),
    u=B/r, xi=lambda B^2, delta=7/2, x=log Y.

Here t is the retained physical time coordinate and its period is
beta. Thus a(0)=sqrt(C lambda) and

    A0=beta a(0), ||v||_0^2=A0^-1 integral_gbar |v|^2
                             =integral e^h |v(Y)|^2 dY.

Only stationary scalar fields are considered. On |Y|<=1/10 the
accepted Riccati bounds give |h|<=Y^2; hence on Y<=1/16,

    1/2 < e^h < 2.                                      (1)

Let G be exactly the sum of the near actual-to-cap bubble difference
and the far actual bubble in section 4 of the curved proof. Each
absolute row, in the actual canonical measure, is below 3e24 C.
This includes the retained global state/domain and thermal effects
in that decomposition. It is not a replacement of the actual
Green function by its cap. Keep the row-sum-subtracted operator

    G0 v=Gv-(G1)v.

Let p0=H_mass 1=Theta_p0/C be the complete pointwise constant
contact of sections 4-5 there. It includes the gluing row sum,
vacuum, tadpole and density contributions:

    |p0(Y)| <=1e20 (q^2/C)(1+1/Y), Y>0.                 (2)

In the exact allocation

    H_mass=Q_geo+L_sing+p0 Id+G0,                       (3)

L_sing is the remaining local singular comparison, with its
finite derivative contacts, chosen to annihilate constants.
This notation does not estimate L_sing or its first derivative
contact. The common fourth restoration is unchanged.
The original seam-supported contact annihilates the selected
compact core and mass cubic. Its extension is the original
completion, not a new boundary condition.

Define the actual dimensionless physical shape component

    Z_epsilon v=1_(0<Y<epsilon) C u^4/(16pi xi^2)
                                      [p0 v+G0 v],       (4)

and its output norm ||Z_epsilon v||_out=
||Y^-delta Z_epsilon v||_L2(dx).
For b=2T, (4) is precisely this contribution to f_b/lambda^2,
where f_b=(f_T-f_S/3)/2 is the dual output for
theta=S+T/3. The massive action depends on T alone.
This physical factor is not replaced by its normal limit.

## 2. Strong component estimate

For 0<epsilon<=1/16 and v supported in 0<Y<epsilon, put
w=Y^-delta v. Then

    ||Z_epsilon v||_out <=1e26 epsilon^3 ||w||_L2(dx).
                                                               (5)

For a mass piece v_m=m chi(Y)Y^3, with |chi|<=1 and support in
0<Y<epsilon,

    ||Z_epsilon v_m||_out <=1e26 |m| epsilon^(5/2).      (6)

For input v_i vanishing on 0<Y<epsilon,

    ||Z_epsilon v_i||_out <=1e25 ||v_i||_0.             (7)

All three statements extend as bounded maps in their displayed
norms. Hard input restrictions are legitimate for this bounded-
kernel/multiplication component. They are not inserted into the
domain of the full fourth-order Hessian.

Here is a direct proof. Cauchy--Schwarz against |G(X,X')| followed
by integration and the reverse row bound gives

    ||G||_(0->0)<=3e24 C, ||G1||_infinity<=3e24 C,
    ||G0||_(0->0)<=6e24 C.                             (8)

The normalization by A0 cancels from this estimate. Restriction
to stationary fields cannot increase its norm. The even-double
and physical-half normalization is already included in the
predecessor's row bound; no further rank or image factor is added.

The seed inequality Y>=sqrt(5)u and C/xi<6.1 imply

    rho(Y)=C^2 u^4/(xi^2 Y^4) <6.1^2/25<3/2,
    C u^4/(16pi xi^2)=rho(Y)Y^4/(16pi C).              (9)

At the selected weight there is the exact measure identity

    ||Y^-7/2 Y^4 F||_L2(dY/Y)=||F||_L2(dY).            (10)

Consequently (1), (8), (9) and rho/(16pi)<1 give

    ||1_(Y<epsilon) Y^-delta
              C u^4 G0 v/(16pi xi^2)||_L2(dx)
        <=(sqrt(2)/C)||G0 v||_0 <=9e24 ||v||_0.        (11)

For a fast input, ||v||_0<=sqrt(2)epsilon^4||w||_2,
so the sharper product of the two square-root factors in (11)
is 12e24 epsilon^4||w||_2.
For the mass input it is 12e24 |m|epsilon^(7/2)/sqrt(7).
The outside input has no multiplication contribution on the
output interval, proving (7).

For the contact use q/C<400. Equations (2) and (9) give the
pointwise multiplier bound, in the conjugated log norm,

    rho Y^4 |p0|/(16pi C)
         <=1.6e25 (Y^4+Y^3).                         (12)

Thus its fast norm is at most
1.6e25(1+epsilon)epsilon^3||w||_2.
Its mass norm is at most

    1.6e25 |m| [epsilon^(7/2)/sqrt(7)
                             +epsilon^(5/2)/sqrt(5)].

These integrals use dY/Y; no finite lower cutoff is present.
Adding the corresponding gluing bounds proves (5)-(6), with
ample reserve below 1e26. In particular the 1/Y contact is
not set to zero or treated as a bounded unweighted potential.

## 3. Normal composition and what it does not include

The H4log multiplier is (1+p^2)^2 log(2+|p|), so
||w||_2<=2||w||_H4log. The complete normal inverse from
[S163](weighted-normal-inverse-scale.md) satisfies
||N^-1||<8.7e19/lambda^2 in these exact theta,b and output norms.
For inputs supported in the collar, extend the output by zero
on the x-line. The physical component is lambda^2 Z_epsilon.
Therefore

    ||N^-1 lambda^2 (0,Z_epsilon)|| <=8.7e46 epsilon^3,
    <1e-130 at epsilon=1e-60.                         (13)

We used 2e26<1e27 in this display. This is a bound for the two
components in (4), not for the complete actual-minus-normal map.

The earlier [flat tower correction](outer-flat-tower-normal-remainder.md)
can also retain the actual density factor rho. Multiplication on
the output costs less than 3/2, with no differentiated rho.
Thus that correction together with (4) has normal-inverse budget
below 2e-10 at the same epsilon. The term (rho-1)K_infty,
its necessary Weyl cancellation, the local curved singular
comparison and derivative contacts are NOT part of this sum.
No Neumann inversion of the actual operator is asserted.

## 4. Mass extension and the global value tail

For an actual compatible coupled solution, the leading disk mass
vector is (3/4)m_c u^3(1,-2). Since
u=sqrt(xi/C)Y+O(Y^3), the leading coefficient of b=2T in Y is

    m_Y=-3(xi/C)^(3/2)m_c, |m_Y|<|m_c|/4.             (14)

The inequality uses C/xi>5.9; no higher York norm is inferred.
The remaining fast coefficient is the actual remainder.
Split b algebraically at Y=epsilon into
1_(Y<epsilon)(b-m_Y Y^3), 1_(Y<epsilon)m_Y Y^3 and
1_(Y>=epsilon)b. This is only a decomposition for (4).

The physical four-dimensional measure is
dmu=4pi r^4/C dvol_gbar. With gamma0=4pi beta B^3 and r>=B,

    ||b||_0 <=(C/xi)^(1/4) gamma0^(-1/2)||b||_phys
             <2 gamma0^(-1/2)||b||_phys.               (15)

This follows also by substituting A0=beta sqrt(C lambda)
before taking a square root. The accepted solution value and
[actual mass bound](tt-mass-numerical-bound.md) are

    ||b||_phys <=1.116e10 g lambda^-2 Fdisk,
    |m_c| <=5.40144e14 g lambda^-2 gamma0^(-1/2)Fdisk,
    g=G lambda.

Combining (5)-(7) and (14)-(15) gives

    ||Z_epsilon b||_out
      <=1e26 epsilon^3
          ||1_(Y<epsilon)Y^-7/2(b-m_Y Y^3)||_L2(dx)
       +(2e40 epsilon^(5/2)+3e35)
                       g lambda^-2 gamma0^(-1/2)Fdisk. (16)

The final 3e35 term is a global value/source bound, not a small
collar perturbation. It keeps the nonlocal input from the rest
of the disk. Equations (6), (14) do not put the mass cubic in
the fast input space: its output is integrable after (4).

## 5. Frozen verification and next step

Freeze this proof and its checker before the scientific run.
Pin the curved-response, flat-correction, physical-source and
mass receipts and their exact inputs. Verify the density/dual
normalization, row-subtraction and constant preservation,
stationary and physical measure conversions, exact weight
powers and infinite-end integrals, the separate fast/mass/value
budgets, mass-coordinate conversion and complete normal
composition factors. Retain every run and failed source snapshot.

The checker supports the arithmetic and algebra in this ordinary
operator argument. It is not formal verification of the continuum
theorem or a numerical solution of the actual field equations.

The remaining strong outer task is now specifically the local
curved singular comparison, its first/second derivative contacts,
the actual density/principal and Weyl combination, and cutoff
matching. The old interior H3 point-probe estimate is not uniform
at the seam and cannot be used there without a new argument.

Remaining to-do list: local curved outer response and derivative
contacts; annulus matching; higher York maps and full inverse/
residual composition; nonlinear signed compatibility; exterior/
state, charge, return, opening, stability, transport, preparation
and creation.
