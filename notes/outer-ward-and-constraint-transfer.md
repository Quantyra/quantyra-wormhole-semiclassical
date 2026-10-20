# Outer scalar-to-constraint transfer at finite geometry

2026-10-19, E003/S147 supporting S139/S074. Ordinary argument and
verification protocol before execution.

This carries the regular canonical scalar estimate of S146 into
the outer physical geometric-mass and lapse integrals. Geometry
varies inside both the Ward kernel and the constraint prefactors.
The pressure inherited from the interior is kept as an explicit
datum and is not assigned a value by this calculation.

## 1. Exact outer proper gauge and physical variables

Use the retained B,A,C=B^2/A^2,x=lambda_e B^2 and fixed physical
time t. Set d=sqrt(x/C), so .4<d<.42. Define tau=sqrt(x)t/A,
whose period is the original physical period times sqrt(x)/A.
This constant change does not select or rescale the physical clock.
On the outer collar write

    g4=u(Y)^(-2){A^2[dY^2+a(Y)^2 d tau^2]+B^2 dOmega^2},
    r=B/u, u(0)=0, u'(0)=d, a(0)=1.                 (1)

Y increases inward and u'>0. In the physical areal gauge,

    F=C u'^2/u^2, exp(alpha)=d a/u'.                 (2)

Indeed dr/dY=-B u'/u^2, and the time coefficient is x a^2/u^2.
These are exact identities, not first-order coordinate formulae.
We use a and u as independent functions in this fixed collar.
There is no differentiation of a moving C5 proper-coordinate map.
An extension through the interior, bolt and fixed paired domain
must still satisfy the S146 contract and the eventual global equations.

Let F_ref(r) be the fixed covariant reference-potential solution
used by the centered reconstruction. Its mass parameter drops out
of F_ref+r F_ref'. With the repository's sphere-integrated
rho,p convention, the exact mass and lapse constraints become

    2C u u''-3C u'^2
       +u^2[F_ref(B/u)+(B/u)F_ref'(B/u)]
          =2G rho u^2,
    C[u''-(a'/a)u']=G u(rho+p).                     (3)

For example the first follows by differentiating
m_rest=r(F_ref-F)/2 and using dm_rest/dr=G rho.
The second follows from alpha_r=G(rho+p)/(rF).
Equivalently their radial-pressure combination is

    3C u'^2-2C u(a'/a)u'
       -u^2[F_ref+rF_ref']=2G p u^2.

We do not replace (3) by a scalar trace equation. The separate
angular equation follows only when the complete joint Ward
identity and source prescription hold.

## 2. The physical mass-cubic correlation

Write the fixed outer jets as

    a=1+a2 Y^2+c Y^3+O(Y^4 log^k Y),
    u=dY+e3 Y^3+e4 Y^4+O(Y^5 log^k Y).

The retained constant term in F and the absence of quadratic and
cubic areal-lapse terms imply, directly from (2),

    a2=1/(2C), e3=d/(6C), e4=d c/4,
    m=-BC c/d.                                      (4)

Here m is the geometric mass coefficient in
F=1+lambda_e r^2-2m/r+O(r^-2 log^k r).
To check the last relation, C u'^2-x has quadratic coefficient
6C d e3 and cubic coefficient 8C d e4. Thus the constant and
linear-in-Y terms of F-x/u^2 are 6C e3/d and 8C e4/d.
The latter is -2m d Y/B. Meanwhile

    alpha=(a2-3e3/d)Y^2+(c-4e4/d)Y^3+O(Y^4 log^k Y).

In S146 the analytic mass direction has gamma=d^3 exactly.
It follows from its explicit proper-distance formula, whose
leading term is (x/C)^(3/2)Y^3. Therefore the parameter
p=delta m/(Bx) changes the finite jets by

    delta c=-d^3 p,
    delta e4=-d^4 p/4,
    delta(log u)=-d^3 p Y^3/4+O(Y^4 log^k Y).        (5)

The radius quartic is consequently not a second free mass datum.
This supplies an exact local identification between S146's mass
parameter and the physical mass coefficient, provided the stated
asymptotic lapse condition is retained. Without that condition,
the cubic lapse coefficient is c-4e4/d and the identification
cannot be inferred from a alone.

## 3. Ward reconstruction, anomaly and matching memory

Put

    rho=(u^2/B^2) H_rho, p=(u^2/B^2) P,
    T=P-H_rho,
    Abar=(B^2/u^2)(-rho+p+2pS),

where pS is the sphere-integrated angular pressure. The four-
dimensional Ward identity in (1) is exactly

    P'+2(a'/a)P=(a'/a)T-(u'/u)Abar.                 (6)

For example the original identity is
p'+(rho+p)(a'/a-u'/u)+2pS u'/u=0. Substitution gives (6);
the radius terms cancel only when the full trace allocation
is retained. The canonical conformally invariant sector has
Abar=0; Weyl and reference contributions must supply their
own Abar and stress.

Let y0=1/20 and P_m=P(y0). At fixed Y in (0,y0],

    P(Y)=H_b(Y) P_m+L_b[T](Y)+J_(b,u)[Abar](Y),
    H_b(Y)=a(y0)^2/a(Y)^2,
    L_b[T]=-a(Y)^(-2) integral_Y^y0 a(s)a'(s)T(s) ds,
    J_(b,u)[Abar]=a(Y)^(-2)
                   integral_Y^y0 a(s)^2 u'(s)Abar(s)/u(s) ds. (7)

In the full solution P_m is fixed by the actual interior state
and regular bolt. It is not a new adjustable boundary source.
Equation (7) separates the local outer map from that unresolved
matching data. Its derivatives must also include derivatives of
P_m when P_m depends on the geometry/state.

For the canonical calculation below, write T=R+c C^2/(5pi Y)
as in S146. By linearity, decompose its Ward pressure into the
two corresponding pressures and their matching values. The
R-sector is the regular component used below. The explicit cubic
sector stays allocated to its correlated Weyl partner before
physical absolute integration. This decomposition neither
changes the renormalization scheme nor deletes a physical term.

## 4. Finite-background Ward operator bounds

Take a=a0 exp(-b) with the S146 real background N(b)<=1e-7.
The centered seed has |Rhat|<=2, a0(0)=1,a0'(0)=0. Hence
cos(Y)<=a0(Y)<=cosh(Y), and |a0'/a0|<=tan(Y)<2Y on
this collar. The S146 cap bound gives |b|<=1e-5.
For its analytic mass direction, Cauchy on the half-width .005
strip gives |psi'''|<=4.8e8. Since psi'(0)=psi''(0)=0,
|psi'|<=2.4e8Y^2. The weighted remainder has
|v'|<=100 N(b)Y^3 L(Y). It follows that

    .99<a<1.01, |a'/a|<4Y.                          (8)

Assume the independent radius field satisfies

    dY/2<=u<=2dY, d/2<=u'<=2d.                     (9)

These are explicit low-jet radius hypotheses; they are not supplied
by the canonical determinant. They hold on an open neighborhood
of the seed but need preservation by the future coupled solve.

Let H(Y)=1+log^2 Y and ||T||H=sup |T|/H.
For shape directions h,k set U>=N(h), W>=N(k); then
|h'|<=100U, |k'|<=100W on the cap. In particular
|h(Y)-h(s)|<=100U|Y-s|. Define

    I_j=integral_0^y0 Y^j H(Y)dY
       =y0^(j+1)/(j+1)
          [1+log^2(y0)-2log(y0)/(j+1)+2/(j+1)^2].

The useful enclosures are

    I0<.899, I1<.01684, I3<.00001813,
    J3=integral_0^y0 Y^3 dY=1/640000.                (10)

In (7), K_b(Y,s)=a(s)^2/a(Y)^2<1.05 and

    D_h K_b=2[h(Y)-h(s)]K_b,
    D_h D_k K_b=4[h(Y)-h(s)][k(Y)-k(s)]K_b.

Also D_h(a'/a)=-h'. Thus the operator norms, for all real
backgrounds in the ball, satisfy

| Operator | Value | First shape derivative / U | Mixed derivative / UW |
| --- | ---: | ---: | ---: |
| H_b: scalar to bounded pressure | 1.05 | 11 | 110 |
| L_b: weighted trace to bounded pressure | .071 | 100 | 2000 |

Here the second row follows from

    1.05*4 I1<.071,
    1.05*(10*4 I1+100 I0)<100,
    1.05*(100*4 I1+2000 I0)<2000.                   (11)

No derivative of T is used. The factor h(Y)-h(s) is retained;
replacing it by two independent value bounds would lose these
constants. The derivative of a'/a is also retained.

For completeness, a full-trace input
|Abar(Y)|<=Z Y^2 H(Y) has |J[Abar]|<.071 Z by
(8)-(10) and u'/u<=4/Y. This only proves its value map.
No finite-background derivative bounds on the independent
Weyl/reference Abar are claimed here.

## 5. Canonical source jets and composite pressure

The centered source theorem gives
|R[0]|<=1e20 q^2 L(Y). Since q/C<400 and L<=2H,
its weighted norm is below 3.2e25 C^2. S146's difference
bound and actual-domain correction then give, for N(b)<=1e-7,

    ||R[b]||H <=K0 C^2, K0=4e25,
    ||D_h R[b]||H <=K1 C^2 U, K1=2e25,
    ||D_h D_k R[b]||H <=K2 C^2 UW, K2=1e30.         (12)

The first derivative is not inferred from the second derivative
alone. Apply Cauchy to R_cap[b+z h]-R_cap[0] at radius
1/(40000U): its bound is
1e25(1e-7+1/40000)*40000 C^2 U=1.004e25 C^2 U.
The actual-cap first derivative is below .48 C^2 U.
Both fit K1. The fixed common a0 value stays inside K0.

Define P_reg=L_b[R], the trace-forcing part of (7).
This is one summand of the pressure, not a decision that P_m=0.
The preceding operator product rule gives

    ||P_reg||inf <=p0 C^2, p0=.071 K0=2.84e24,
    ||D_h P_reg||inf <=p1 C^2 U,
        p1=.071 K1+100 K0=4.00142e27,
    ||D_h D_k P_reg||inf <=p2 C^2 UW,
        p2=.071 K2+200 K1+2000 K0=1.55e29.          (13)

These include the changing Ward kernel and its first two derivatives.
The matching contribution and its derivatives remain H_b P_m,
with the first row of the operator table and the ordinary product rule.

## 6. Outer mass and lapse maps, including radius variations

For any source sector with Abar=0 and pressure P, the outer
constraint integrals from (2)-(3) are

    M_out=[m_rest(infinity)-m_rest(r(y0))]/(Bx)
         =G/(B^2 x) integral_0^y0 u'(P-T)dY,
    A_out=alpha(y0)-alpha(0)
         =-G/(C B^2) integral_0^y0 (u^3/u')(2P-T)dY. (14)

The first is a geometric mass difference. It is not a Hamiltonian
charge, and it does not determine the total mass without the interior.
Absolute values of the integrands are integrable for T=R because H
is integrable. The unrecombined cubic trace would produce a logarithm
in the mass integral and must not be bounded separately.

In addition to the shape norm, require the radius direction v to obey
|v/u|,|v'/u'|<=U, and similarly w with W. Use affine variations
of u and logarithmic variations of a at fixed Y. The combined
direction norms are the maxima of these bounds and the shape norm.
In particular u' has zero second variation. For Q_u=u^3/u',

    |D_v Q_u|<=4U Q_u,
    |D_v D_w Q_u|<=14UW Q_u,
    0<Q_u<=16d^2 Y^3.                              (15)

The mixed factor is
6(v/u)(w/u)-3[(v/u)(w'/u')+(w/u)(v'/u')]
+2(v'/u')(w'/u'), which proves the coefficient 14.
This includes variations of the physical source Jacobian; it is
not an estimate on the scalar trace alone.

Insert P_reg,R into (14), retaining the memory sector separately.
The fixed parameter bounds GC^2/(B^2 x)<5e-32 and
GC/B^2<8.5e-33 give

    |D^2 M_reg[h,k]| <.039 U W,
    |D^2 A_reg[h,k]| <4.6e-7 U W.                   (16)

The explicit unrounded upper bounds are

    (2*.42)(5e-32){y0(p2+2p1)+I0(K2+2K1)} <.039,
    16*.42^2*(8.5e-33)
      {J3(2p2+16p1+28p0)+I3(K2+8K1+14K0)} <4.6e-7.

Their values are below .038073 and 4.491e-7. Along an affine/log
segment that preserves the hypotheses, the Taylor remainders are
therefore below .0195 epsilon^2 and 2.3e-7 epsilon^2.
At epsilon=1e-7 these are 1.95e-16 and 2.3e-21 for this component.
They do not bound the uncomputed matching-memory change.

The memory contribution has the explicit form

    M_mem=K_M(b,u) P_m,
    K_M=G/(B^2 x) integral_0^y0 u' H_b dY,
    A_mem=K_A(b,u) P_m,
    K_A=-2G/(C B^2) integral_0^y0 (u^3/u')H_b dY.

Useful coefficient bounds through order two are

| Coefficient, with physical prefactor removed | Value | First / U | Mixed / UW |
| --- | ---: | ---: | ---: |
| K_M / [G/(B^2 x)] | .05 | .51 | 5.6 |
| K_A / [G/(C B^2)] in absolute value | 1e-5 | 1.4e-4 | .002 |

For example the mass entries follow from
2d y0 times (1.05,1.05+11,110+22); the lapse entries follow
from 32d^2 J3 times (1.05,4*1.05+11,14*1.05+88+110).
When P_m varies, retain K P_m,hw+K_h P_m,w+K_w P_m,h+K_hw P_m.
This names the exact matching information still required.

## 7. Verification and consequence for the research plan

The checker is to verify (2)-(6) independently from the areal
constraints and four-dimensional Ward identity; the correlated
outer series and mass coefficient; the differentiated kernels;
the exact radial factors and all integral budgets. Numerical
tests use a manufactured weighted logarithmic trace and
nonzero matching pressure, direct quadrature of the Ward solution,
and independently differenced composite mass/lapse maps at two
steps and two quadrature resolutions. Omitting the Ward-kernel
variation or the matching pressure must be detected.

Those checks are algebraic and numerical diagnostics. The actual
Dirac source estimates remain the ordinary arguments in the cited
notes. No physical state or solution is inferred from a manufactured
trace. The source/frame allocation here agrees with the retained
Einstein-Maxwell/charged-fermion framework; the original
[Maldacena--Milekhin--Popov construction](https://arxiv.org/abs/1807.04726)
does not supply these new finite-background norm constants.

This proves the outer mass/lapse part of the nonlinear transfer
for the regular canonical sector and gives its precise dependence
on matching memory. It also fixes the correlated radius jet for
the physical mass. The interior/cold reconstruction and actual
inverse norm, the complete Weyl/LLL/reference/radius response,
global charge and state matching remain unresolved.

Remaining to-do list: determine and control the interior matching
memory and bolt response; complete the other source mappings and
the actual global inverse composition; solve the coupled exterior,
state and signed source; establish opening, stability, transport,
preparation and creation.
