# Global thermal curvature-free shape bound

2026-10-21, E003/S168. Ordinary continuum argument before focused checks.
The full research goal remains incomplete.

For the retained charged seed, fixed HH state and common subtraction,
the COMPLETE shape Hessian restricted by Psi=delta R-r R_r p=0 obeys

    S_full,gamma[v] >= (1/(40G)) integral r^2(v_y^2+gamma^2 v^2)dy
                     +C J2[v]/(200pi)+C J2[p+v/3]/(40pi)           (1)

at every sufficiently large positive Matsubara frequency gamma=2pi n/beta.
Here v is any smooth compact shape profile in the entire open exterior,
and p is its global conformal reconstruction, including both tails.
The three components are the classical action, fixed four-dimensional
Weyl action and complete canonical Dirac response. No angular tower,
second vertex, seed residual or finite counterterm is removed.

This is global coercivity on the stated core, with a continuous extension
to the stronger completion specified below. It does not establish a
closed full form in its own energy norm or a full variational inverse.
The constants defining the frequency threshold are finite but not
numerically calibrated. Discrete thermal positivity does not determine
continuous causal frequencies, the curvature cross block, a physical
gravitational boundary domain or stability of a nonlinear solution.

## 1. Coordinates and classical input

Use physical optical time t and r_y=-F, with y=0 outer and y=infinity
at the horizon. Set w=r v, q_c=r p and

    N_gamma=integral[w_y^2+(gamma^2+2/y^2)w^2]dy,
    Gv=integral r^2(v_y^2+gamma^2v^2)dy,
    Gp=integral r^2(p_y^2+gamma^2p^2)dy.

The [global classical theorem](global-classical-curvature-free-form.md)
gives, above its finite Gamma_C,

    S_C>N_gamma/(32G), Gv<=1.001N_gamma, E_gamma<10N_gamma,
    integral F q_c^2<=10N_gamma/(27lambda),
    integral F w^2<=N_gamma/lambda.                              (2)

Keep its exact physical charge and residual contacts. Its positive
reconstruction is q_c=-K_gamma^-1 B_gamma w/3. Since
H_R=2R+r R_r=-24lambda+(d/B^2)(-12z^5+24z^6)<0,
z=B/r, integration by parts also gives

    Gp=<q_c,K_gamma q_c>+integral F H_R q_c^2/6
      <=E_gamma/9<10N_gamma/9.                                  (3)

The seed guards x=lambda B^2>1 and |d|<1e-12 suffice for this sign.
They are much weaker than the retained benchmark guards.

For the canonical orbit metric and its constant rescaling,

    ghat=B^2 F/r^2 (dt_E^2+dy^2), gbar=C ghat/B^2,
    dY/dy=a_phys=sqrt(CF)/r, gbar=dY^2+a_phys^2 dt_E^2,
    a_phys,*^2=sup a_phys^2<=5C lambda.

J0,J1,J2 denote the gbar squared field, gradient and Laplacian norms.
All forms here are divided by the actual time-mode norm beta/2. Thus

    C J2[T]=B^2 integral_hat(Delta_hat T)^2,
    J1[T]=integral(T_y^2+gamma^2 T^2)dy,
    J0[T]=C integral F T^2/r^2 dy.                               (4)

There is no conversion of gamma to an arbitrarily rescaled clock.

## 2. Global canonical logarithmic compensation

We first prove, for all compact smooth thermal Fourier sums T, that
finite nonnegative A,B independent of frequency and support satisfy

    Qcan[T]-C/(30pi) integral_bar (omega+19/60)(Delta_bar T)^2
       >=C[J2[T]/(100pi)-A J1[T]-B J0[T]], omega=log(r/B).         (5)

The [uniform outer comparison](uniform-logarithmic-frequency-comparison.md)
already proves the required two-dimensional scalar comparison on
0<Y<.02. The [curved quantum form](stationary-curved-quantum-tt.md),
sections 4-7, gives the global actual/geodesic comparison

    |Qcan-Qgeo|<=C[.001J2+3e24J1+3e26J0],

including contacts and global returns. Those two-dimensional operator
statements precede the static TT restriction in that note.

Use its partition chi_o=cos(theta), chi_i=sin(theta), where theta is
zero below a=.01, pi/2 above 2a, and (pi/2)I_s(6,6) between, with
s=(Y-a)/a. This is C5 with bounded piecewise sixth derivative, constant
near seam and bolt. The needed geodesic localization estimate is

    |Qgeo[T]-sum_j Qgeo[chi_j T]|
       <=C[epsilon J2+K_epsilon J1+L_epsilon J0]                 (6)

for arbitrary two-dimensional T. Here is the extra contact argument
needed beyond the earlier radial calculation. The separated commutator
kernel contains sum_j(chi_j(x)-chi_j(x'))^2, gaining two powers against
the fourth-order r^-6 kernel. It is bounded by KC r^-4 locally.
The common regulated finite part fixes its contact of order at most two;
its fourth principal contact has canceled in the double commutator.
To bound ALL second contacts, use compact normal-coordinate probes
1, z1, z2, z1^2, z1 z2, z2^2 times a fixed cutoff, translated through
the transition neighborhood. Both tangent directions are included.
This neighborhood is a compact smooth region away from seam and bolt,
so the geodesic derivatives of these probes and their products with chi
through order six have finite uniform bounds. The exact point estimates

    |P_M u(x)|<=2 sup|X^6u|/(525pi M^4),
    sum d_M |P_Mu(x)|<2e-5 C sup|X^6u|,
    |Q4 u(x)|<=C sup|X^4u|

from the same geodesic representation bound those contact coefficients.
When a probe's support misses the transition its local commutator
vanishes. The zero contact still includes -sum chi_j Qgeo(chi_j);
the global derivative bounds control it, including its nonlocal row
tail. Finite long-flow row mass follows from the retained Bessel tail
moments and Liouville invariance, including returns and conjugate points.
No bound on only a separated kernel is used to choose a new contact.

For the anchored symmetric KC r^-4 remainder in two dimensions, Taylor
subtraction through its fixed quadratic moment gives the same estimate
as section 6 of the curved quantum note: with X=||grad T|| and
Z=||Hess T||, split at b to bound it by
KC[bXZ+b^2Z^2/8+2log(R/b)X^2]. Choose b small for epsilon and use
the integrated Bochner identity and bounded curvature to replace Z^2
by J2 plus a fixed multiple of J1. The bounded contacts and tail add
finite J1,J0 constants. This proves (6) without a frequency-dependent
constant or derivatives of the kernel envelope.

There is also a time term in the fourth-order localization identity.
Write h=log a_phys and W=omega+19/60. After radial and periodic time
integration by parts, the difference
sum integral W(Delta(chi_j T))^2-integral W(Delta T)^2 equals the
integral of

    6W theta'^2 T_Y^2 + 2W theta'^2 a_phys^-2 T_t^2
    +[8W theta'theta''+4W h'theta'^2+2W'theta'^2]T T_Y
    +W[theta'^4+(theta''+h'theta')^2]T^2.                        (7)

All coefficient derivatives are supported in the fixed transition.
Consequently (7) is bounded by a finite constant times J1+J0. The
time term is essential for nonzero modes. The unweighted version also
gives sum J2[chi_j T]>=J2[T]-K(J1+J0). The unbounded logarithm is
never differentiated at the outer seam.

On the outer piece the scalar logarithmic bound is

    Qsc-C/(30pi) integral W(Delta T)^2
       >=C[33J2/(600pi)-12J1/(30pi)].

Its comparison with Qgeo costs .001J2 plus finite lower-order terms.
On the interior omega<6 and kappa4>13.9. The positive geodesic model
there gives a fourth-order reserve >7.5C J2/(30pi), again with a
finite J1 cost. Take epsilon=.001 in (6), apply (7), and then use the
global actual/geodesic comparison. The total fourth-order loss is
at most .003, and 33/(600pi)-.003>1/(100pi). All other localization
costs are finite multiples of J1,J0, proving (5). Enlarge A,B to be
nonnegative. Their finiteness, rather than a numerical calibration,
is the assertion.

## 3. Exact Weyl completion and a curvature-free residual budget

Use physical ghat derivatives in this section. Set S=p,T=v and

    Rhat=R[ghat], U=Delta_hat omega, V=|grad_hat omega|^2,
    yB=B^-2, k=B^-4(1/15-qmag^2/3), qmag=2000^4,
    dT=Delta_hat T, f=Delta_hat(S+T/3),
    womega=grad_hat omega.grad_hat S, jS=|grad_hat S|^2.

The unrestricted Hessian identity in section 8 of the curved quantum
note holds before any static TT estimates are made:

    W4''=B^2/(20pi) integral f^2
      -B^2/(30pi) integral(omega+1/6)dT^2
      -B^2/(2pi) integral[dT D+f L+J jS+Wc T womega
                                             -7womega^2/45+Z],
    D=S(-Rhat/15+yB/9)
      +T[(6omega+1)Rhat/45+11V/135]-11womega/135,
    L=T[-Rhat/15+2U/5+7V/45]-7womega/45,
    J=Rhat/30-7(U+V)/90-yB/18,
    Wc=-2Rhat/15+14(U+V)/45,
    Z=T^2[omega Rhat^2/30+Rhat(U+V)/15-U^2/5
                    -7UV/45-7V^2/90+2k omega]
      +ST[-Rhat^2/30+2k].                                      (8)

This follows directly by twice differentiating the fixed bulk density
and its volume factor, with no reference-potential subtraction. The
full decomposition used here is S_C+W4+Gamma_hat. Young's inequalities
|dT D|<=dT^2/100+25D^2 and f^2/10-|fL|>=f^2/20-5L^2 give

    W4''>=-C/(30pi) integral_bar(omega+19/60)(Delta_bar v)^2
       +C J2[p+v/3]/(40pi)
       -B^2/(2pi) integral_hat[25D^2+5L^2+|J|jS+|Wc v womega|+|Z|]. (9)

The favorable womega^2 is dropped only after completing the Hessian.
In particular 19/60=1/6+3/20; the constant contact remains.

Put t_r=r/B>=1. Retained global geometry bounds are
|Rhat|<=13lambda, |U|<=22lambda t_r^2, V<=5lambda t_r^2,
yB<=lambda, omega<=t_r^2, F/r^2<=5lambda. Equations (2)-(4),
and womega=-r p_y/B^2, imply the following integral bounds in units
lambda N_gamma/B^2:

    integral lambda^2 t_r^4 p^2 <=10/27,
    integral lambda^2 t_r^4 v^2 <=1,
    integral womega^2 <=50/9,
    integral lambda t_r^2 jS <=10/9,
    integral lambda t_r^2 |v womega| <=sqrt(50/9)<3.

The coefficient of S in D has absolute value <lambda; its T
coefficient is <3lambda t_r^2. The T coefficient of L is
<11lambda t_r^2, |J|<3lambda t_r^2 and |Wc|<11lambda t_r^2.
Three-term and two-term square inequalities therefore yield

| Integral | Bound in lambda N_gamma/B^2 units |
| --- | ---: |
| D^2 | 30 |
| L^2 | 250 |
| abs(J) jS | 4 |
| abs(Wc v womega) | 33 |
| abs(Z without k) | 160 |

For the last row the pure v^2 coefficient is bounded by
169/30+351/15+484/5+770/45+175/90<145; the mixed pv
coefficient adds at most (169/30)sqrt(10/27). The k terms cost
at most 4qmag^2 N_gamma/(3lambda B^6) before multiplying by B^2/(2pi).
Thus the entire last line of (9) is bounded by

    [2197lambda/(2pi)+2qmag^2/(3pi lambda B^4)]N_gamma
       <[400lambda+qmag^2/(lambda B^4)]N_gamma
       <1e-24 N_gamma/G.                                      (10)

The last inequality uses G lambda<(4e-17)^2 and
G qmag^2/(lambda B^4)<4e-27. No fourth-order term is treated as
an H1 remainder: its logarithmic negative part and positive square
were separated explicitly in (9).

## 4. Thermal threshold, endpoints and precise domain scope

In the n-th thermal sector, the positive scalar Laplacian on gbar has
L_gamma>=z=gamma^2/a_phys,*^2. Hence J1<=J2/z and J0<=J2/z^2.
Define

    Z_Q=max(1,400pi A,sqrt(400pi B)),
    Gamma_T^2=max(Gamma_C^2,5C lambda Z_Q).                     (11)

For every integer n>=1 with gamma_n^2>=Gamma_T^2, (5) is at least
C J2[v]/(200pi). Combining it with (2), (9) and (10), then using
Gv<=1.001N_gamma, proves (1). Thermal Fourier sums supported above
this threshold obey the summed inequality. No interpolation to the
gaps between the thermal frequencies is asserted.

For compact v, the reconstructed p has outer behavior O(y^4), with
the corresponding differentiated estimates. Near the bolt it is a
smooth n-th harmonic, proportional to rho^n times a smooth function
of rho^2. This follows from the homogeneous reconstruction ODE away
from the source, its decaying horizon branch and integer gamma/kappa_h.
The Weyl integrals and their boundary fluxes converge. At the outer
end its jets make every boundary flux vanish; the bolt is an ordinary
interior point in smooth polar coordinates. Do not truncate p
independently or impose a punctured-bolt fourth-order cutoff.

An optional explicit stronger completion is the completion of this
compact-v core in

    Nplus=N_gamma+sum d_M P_M[v]
           +C{J2[v]+J2[p+v/3]+integral_bar omega(Delta_bar v)^2}.

Here P_M is the retained positive geodesic spectral form. Its closed
positive summands, weak Laplacians, weighted Laplacian, and the fixed
continuous reconstruction map show that this completion embeds in the
N_gamma space: if a Cauchy sequence tends to zero there, local
distributional convergence and closedness force all its graph limits
to zero. Equations (2)-(10), their absolute versions, and polarization
make S_full continuous in Nplus, so (1) passes to this completion.
This does not show that the full form norm is equivalent to Nplus.
In particular, the two unbounded logarithmic pieces can cancel at
relative bound one. No closed-form representation theorem or inverse
for the complete shape operator is inferred from positivity alone.
This auxiliary completion also does not select the physical boundary
or initial-state data of the charged metric problem.

## 5. Sources, verification and consequence

[Beckner, Corollary 2](https://arxiv.org/pdf/math/0612486v2) supplies
the established logarithmic uncertainty comparison behind the retained
outer bound. [Mazur--Mottola, equations (54)-(59)](https://arxiv.org/html/hep-th/0106151v1)
separate the fixed local curvature representative, anomaly action and
independent Weyl-invariant functional. Neither source proves this
charged-seed theorem. The contribution here is the global localization
with all time contacts, the curvature-free Weyl budget and their
combination with the retained classical estimate. The closest retained
comparators are the outer-only quantum bound and global classical form.

The checker verifies the unrestricted Hessian, full weighted localization
identity, physical rescalings, exact budget arithmetic and threshold
implication. Independent warped-cylinder quadrature detects omission
of the time term. Charged-seed reconstruction fixtures compare the
completed Weyl Hessian with symmetric finite differences of its full
nonlinear density, refine the ODE and quadrature, and retain conformal
tails. They do not evaluate the canonical Dirac kernel or certify the
benchmark threshold. Continuum estimates remain ordinary proofs;
no Lean run or formal continuum certification is claimed.

Remaining to-do list: close the full form in its own energy norm or
construct its required inverse; control continuous causal frequencies
and curvature coupling; impose physical boundary/initial data and
transfer to a nonlinear solution; complete signed compatibility,
exterior/state, charge/return, opening, all-sector stability, transport,
preparation, separate creation and full audit.
