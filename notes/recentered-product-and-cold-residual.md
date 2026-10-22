# Quantum product mismatch and the actual cold horizon residual

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The full quantum product differs from the CURRENT c_R=10^22 RN reference
by less than 5.5e-23 in squared sphere radius. Its squared orbit-radius
increment is positive and between 1.56e-29 and 1.71e-29. These replace
neither a global response estimate nor the nonlinear admission test.
The old 9.43e-8 allowance used a different center and included the local
R^2 shift. Subtracting that shift exactly is substantive cancellation,
not improvement by renormalizing a norm.

Centering an explicit finite comparison on this exact product gives a
complete horizon residual below 2e-12 of the prescribed thermal forcing.
This includes the local R^2 derivative terms and the actual paired HH
Dirac source, by the retained fixed-seed comparison with hypotheses
checked below. It is a SOURCE VALUE on one explicit geometry. The
inverse and its nonlinear neighborhood have not been transferred to it.

## 1. Fixed theory and the two product equations

Keep G, q=2000^4, g_R=2000^-4, lambda0=2000^-10/G, c_R=10^22,
mu0=1/b0, the prescribed temperature and original paired matter domain.
All other additional finite coefficients remain zero. Write

    g+11g^2/(360pi)=2000^-10, lambda=g/G,
    F=1-384pi c_R G lambda, Z=piGq^2/(g_R^2 F),
    Z=B^2+3lambda B^4, C=1+6lambda B^2,
    Y0=B^-2, X0=C/B^2=6lambda+Y0,
    p_v=-3(lambda0-lambda)/(8piG)=-11lambda^2/(960pi^2).

For the product set X=a^-2,Y=b^-2. The full determinant supplies
pA=q^2 Y^2 A(r,L), pS=q^2 Y^2 S(r,L), with

    r=X/(qY), L=log(mu0^2/(qY)), A+S=-a2/2,
    a2=(15-3/q^2-5r/q-3r^2)/(720pi^2).

The two equations of the UNCHANGED action are

    0=-Y-3lambda0+(piGq^2/g_R^2)Y^2
                           +32piG c_R(X^2-Y^2)-8piG pA,
    0= X-3lambda0-(piGq^2/g_R^2)Y^2
                           -32piG c_R(X^2-Y^2)-8piG pS.       (1)

At (X0,Y0) the local left sides before the quantum terms equal 8piG p_v.
This follows from FZ=piGq^2/g_R^2 and X0-Y0=6lambda. In particular,
using lambda in place of lambda0 requires subtracting p_v as well.

Use dimensionless variables that retain this cancellation:

    eta=Y/Y0-1, v=(X-Y-6lambda)/X0, chi=alpha eta+v,
    alpha=1/C, e=8piGq^2/(CB^2), beta=64piG c_R X0,
    a_v=p_v/(q^2 Y0^2), F=1-beta(1-alpha).

Here v is a product variable, not the finite horizon radius used below.
Then X=X0(1+chi), Y=Y0(1+eta),
r=(C/q)(1+chi)/(1+eta), L=L0-log(1+eta), and (1) is exactly

    v=e[-(1+eta)^2 a2/2-2a_v],
    F eta+beta(1+alpha eta)v+(beta/2)v^2
           +F(1+alpha)eta^2/2=e[(1+eta)^2 A-a_v].       (2)

The R^2 contribution cancels from the trace equation; its mixing beta v
remains in the other equation. Do not omit it or treat c_R as tunable.

## 2. Quantitative enclosure about the current reference

The retained full-angular theorem gives |A|,|D A|<=100(1+|L|),
|partial_L A|<=1, D=r partial_r. It is an analytic subtracted-integral
bound, not a fit to its floating-point quadrature. On
|eta|,|v|<=1e-20, the benchmark and its existing guards give

    alpha<2.1e-11, 1.5e-26<e<1.6e-26,
    0<beta<1.18e-8, F>1-1.179e-8,
    r<.0036, |A|,|D A|<3334, |a_v|<4e-10,
    .0021<a2<.00212, |D a2|<2e-8.                     (3)

The logarithm uses |log(mu0 B)|<1e-8. In (2) solve the second row
for eta and use the first row for v. This map sends the box into
itself. Its maximum derivative row sum is less than 1.2e-8: the
largest entry is beta/F, while pressure derivatives cost at most
e(5*3334+3)/F and the trace row costs less than 1e-27. The denominators
1+eta and 1+chi stay positive. Banach contraction proves the unique
full product in this box. It is the branch of the earlier product
theorem, since the boxes lie in its uniqueness region.

At the fixed point, the trace row and then the other row sharpen to

    -1.7e-29<v<-1.57e-29, |eta|<5.4e-23,
    |b_star^2/B^2-1|<5.5e-23,
    1.56e-29<a_star^2/(B^2/C)-1<1.71e-29.              (4)

The last lower bound follows from |alpha eta|<1.134e-33. It is a
product-radius lower bound, NOT a lower bound on the actual finite
metric correction. For reproducible numerical illustrations define

    v1=e(-a20/2-2a_v),
    eta1=[e(A0-a_v)-beta v1]/F, chi1=alpha eta1+v1.

Using the analytic pressure Lipschitz bound in (2) gives

    |eta-eta1|<2e-44, |v-v1|<2e-51,
    |b_star^2/B^2-1+eta1|<3e-44,
    |a_star^2/(B^2/C)-1+chi1|<3e-51.                  (5)

These remainders refer to the EXACT pressure integral A0. They do not
certify its adaptive quadrature error. Two proper-time UV splits are
used in the supporting illustration; (3)-(4) need neither result.

The trace shift alone has |v|/(2h)>2500 using the retained
h<3*2000^-10. Numerically this ratio is about 1.5e4--1.8e4 in the
prepared interval. Dividing an uncentered product defect by the tiny
thermal forcing is therefore the wrong perturbation problem. This is
an earlier matching issue, not evidence against existence.

## 3. A finite comparison retaining both product radii

Now write b=b_star, a=a_star, x=lambda b^2, C_s=b^2/a^2,
d=C_s-1-6x. Equations (2)-(4) imply exactly

    d=C v/(1+eta), d/C_s=v/(1+chi),
    |d|<8.6e-19, |d|/C_s<1.71e-29.                    (6)

Use the geometric comparison from the retained covariant-seed note:

    f_ext(t)=(t-1)^2/t^2 [1+x(t^2+2t+3)+d/t^2],
    f(t)=f_ext(t)-delta/t, t=r/b,
    r_+=b V, V=1+h, delta=V f_ext(V),
    kappa b=[(t f_ext)'(V)]/(2V).                     (7)

This fixes the SAME laboratory kappa and outer clock; the small-h
equation has the unique positive root and exterior positivity already
proved for |d|<1 and x>1000. Here x>8e9. The original b0 in kappa and
mu0 is not replaced. These data still give the loose positive-temperature
interval, log(1/h)<80, C_s/q<.0035, M_n>18 sqrt(n), and
|log(mu0 b)|<1e-8. They do not assert an extremal-state limit.

The r^-2 coefficient b^2(1+3x+d) in (7) is a metric coefficient, not a
change in magnetic flux or g_R. We evaluate the ORIGINAL full action
on (7). No extra reference potential is adopted as a physical source.
This comparison has a nonconstant scalar curvature

    b^2 R=-12x+d(4t^-5-6t^-6),                        (8)

so the R^2 derivative terms must be included. It has the prescribed
outer curvature datum and fast curvature tail. Fixed-clock mass
matching and the scalar/shape phase at infinity still require the
complete inverse and return, not just (8).

## 4. Exact local cancellation at the finite bolt

All quantities in this paragraph are dimensionless in b units. Put
Qc=1+3x+d, Delta_x=(lambda0-lambda)b^2, w_c=16piG c_R/b^2.
At t=V define

    j=V f'(V)=1+3xV^2-Qc/V^2+4d/V^3-3d/V^4,
    w=V^2 f''(V)/2=-1+2Qc/V^2-10d/V^3+9d/V^4,
    Rbar=-12x+d(4V^-5-6V^-6),
    Rbar'=d(-20V^-6+36V^-7),
    RicA=-(w+j)/V^2, RicS=(1-j)/V^2,
    HA=2 Rbar RicA-Rbar^2/2+j Rbar'/V,
    HS=2 Rbar RicS-Rbar^2/2+2j Rbar'/V.               (9)

These are the actual R^2 Euler tensor components at the smooth bolt.
The factors one and two on j Rbar'/V follow from
Box R=f'R' and nabla_tau nabla^tau R=f'R'/2 there.
In particular the local derivative terms are not set to zero by using
the homogeneous product equations.

Let LA,LS be b^2 times the EH+cosmological+Maxwell+c_R field equations
before the Dirac subtraction. Eliminating the physical Maxwell
coefficient between its bolt value and its exact product value gives

    LA(V)-LA(1)/V^4
      =d(4V^-5-3V^-6-V^-4)-3Delta_x(1-V^-4)
                                  +w_c[HA(V)-HA(1)/V^4],
    LS(V)-LS(1)/V^4
      =d(-6V^-5+6V^-6)-3Delta_x(1-V^-4)
                                  +w_c[HS(V)-HS(1)/V^4].           (10)

LA(1),LS(1) are 8piG b^2 times the EXACT full-product pressures.
The EH terms of order xh have canceled. For d=0, HA=24xQc/V^4,
HS=-24xQc/V^4, so those large R^2 terms cancel as well.

For an explicit bound expand HA-24xQc/V^4 and HS+24xQc/V^4 as d
times polynomials in z=1/V, x and d. Bound each monomial difference
z^n-z^4 by (n+4)h on 1<=V<=1.001. The coefficient sums are below
1e5(1+x+|d|); this is an elementary finite polynomial bound, not a
sampled derivative claim. The EH derivative sums are at most 66.
Since w_c(1+x+|d|)<1e-9 and |1-V^-4|<=4h,

    |LA,S(V)-LA,S(1)/V^4|<67|d|h+12Delta_x h.

Moreover j>=C_s h, V^2<1.003, Delta_x/C_s<2e-36. Hence

    |local difference|/(2kappa/r_+)
      <68|d|/C_s+13Delta_x/C_s<1.2e-27.              (11)

Here the numerator is restored to its physical b^-2 units. The
denominator is j/(b^2 V^2), not a freely rescaled forcing.

## 5. The quantum transfer and what its proof actually uses

The retained [horizon source proof](stationary-horizon-source-defect.md)
compares the ACTUAL HH source of (7) to the full invariant product with
sphere r_+ and orbit inverse radius f''(r_+)/2. Its geometry extension
uses |d|<1, |d|/C_s<1e-24; (6) is stronger. Its real and complex cap
bounds, angular masses and fixed positive-period bounds therefore hold.
The changes of x and C from the original classical reference are about
1e-8, but these proofs use the displayed strict guards, not equality
with that reference. The scale restoration changes Fang to
4gamma_E-2H_q+4log(mu0 b), still of magnitude below 100. Its a0 term
cancels in the fixed-radius comparison; it is not discarded globally.

The main massive ledger, in h q^2 units, is

    1.2e16*81*.00352^3 <4.24e10,
    4e12*81*.00352 <1.141e12,
    1e12*81*.00352^2 <1.004e9.

These retain the full angular sum after common subtraction, the cap
amplitude homotopy, and separated actual-domain inverse. The massive
inverse <=1/M_n relies on the paired first-order realization, not on
formally squaring an operator with the wrong domain. The distance to
the physical boundary supplies the exponentially suppressed attachment
term. There is no analogous massive-gap argument for the LLL: its
bolt value instead uses its trace and exact HH regularity t_null=0.
The Weyl terms and the nonconstant Box R anomaly of (8) are retained.
These are essential distinctions between the fixed-geometry result
and a hypothetical uniform nonlinear response theorem.

Thus each physical pressure differs from its matched product pressure
by <1e14 h q^2/r_+^4. From (9), |log(w/C_s)|<4h; also
L(r_+)-L(b)=2log V. The same full-angular pressure bound in (3) gives
the further product change <2e4 h q^2/r_+^4. Consequently

    |pA,S[seed](r_+)-V^-4 pA,S[exact product]|
                       <1.001e14 h q^2/r_+^4,
    max |full field residual_A,S(r_+)|/(2kappa/r_+)
                       <1.2e-27+1.004*1.001e14*1.6e-26
                       <2e-12.                       (12)

This closes the horizon VALUE test, below the old .01 target by over
5e9. It does not bound spatial derivatives, metric derivatives or the
global coupled response. No supporting checker can supply those missing
continuum arguments.

The external review retains [Vassilevich, sections 3.3 and 5.3](https://arxiv.org/html/hep-th/0306138v3#S5.SS3)
for complementary projectors and local heat terms. That source does not
prove the nonlocal estimates above. [Camporesi--Higuchi](https://arxiv.org/abs/gr-qc/9505009)
derives homogeneous spinor heat kernels; it does not identify the HH
state of a finite exterior. The separated-resolvent, cap and regularity
arguments in the retained note are required for that use.

## 6. Feasibility decision and numerical requirements

Continue the c_R=10^22 route using (7) for the next direct cold-response
calculation. Use the exact full product, its fixed-scale source and
the SAME temperature; preserve the outer combined cancellation and
mass/kernel moments. This changes the geometric center, not the action
or the norm. Do not transfer the RN inverse by a declaration that the
two metrics are close: the pulled-back domain, mass column, curvature
datum, Hessian contacts and quantum determinant must be compared in
the actual solver/admission spaces.

For a complete physical space X, preconditioner B_p and any REQUIRED
projection Pi, use

    Y0=||B_p Pi F(0)||_X, Z0=||I-B_p Pi F'(0)||,
    ||B_p Pi(F'(u)-F'(v))||<=Z2||u-v||_X,
    Y0+Z0 rho+(Z2/2)rho^2<=rho, Z0+Z2 rho<1.           (13)

At rho=1e-8 the sufficient targets remain Y0<=5e-9, Z0<=.25,
Z2<=5e7. B_p must be injective on the projected source range and the
unprojected compatibility rows must vanish exactly. No old signed
projection is imported into this new realization. A ball of radius rho
must also lie inside EVERY positive geometric, state and source margin.

| Requirement | Quantitative threshold | Status |
| --- | --- | --- |
| Product mismatch from current RN | <5e-9 in radius component | (4) passes that component; no full X statement |
| Product-centered horizon residual | <.01 of thermal forcing | (12) <2e-12 passes on this explicit seed |
| RN energy-to-bolt conversion | <=625/18 for the available 1.44e-10 allowance | Known generic lower bound >2.6e15 fails that certificate; C5 embedding is not finite |
| Direct cold response, if using the scalar allowance 2e-12 alone | response factor <=2500 to allocate 5e-9 | UNKNOWN, and a horizon value is not a complete source norm; first prove a valid full-source response inequality |
| Inverse/domain transfer to (7) | ||B_p Pi(H_seed-H_RN)||<=.25 in the SAME admitted spaces | UNKNOWN; the <.336 mixed form bound and energy inverse <=1 do not supply it |
| Combined outer response | contribution <=5e-9 with its assigned weight and mass correlations | UNKNOWN actual global moments; finite endpoint regularity does not quantify them |
| General cold theorem admission | N_c(seed)+E_c rho<=5e-7 inside its real 1e-6 contract | UNKNOWN transferred background and E_c; E_c=1 would require N_c(seed)<=4.9e-7 |
| General outer theorem admission | N_o(seed)+E_o rho<=1e-7 | UNKNOWN transferred background/E_o; E_o=1 would require <=9e-8; old 1.01e-7 allowance fails to certify slack |
| Interior source family | P45(seed)+E45 rho<=.001 | Old P5<1.3e-4 would allow E45<=8.7e4 only AFTER transfer to this seed |
| Complete source differentiability | full B_p Pi F'' norm <=5e7 on the ball | UNKNOWN; spatial Theta' and fixed-seed pressure bounds are not this norm |

These are simultaneous conditions. The component targets cannot each
consume 5e-9 in an additive norm; either choose a max norm with its
actual operator constants or allocate a sum budget explicitly.

The original high-derivative chain remains parked: C_I=9.1591500000000000015e103,
C_F,I=4.141944868545865e118 and C_F,global=1.069626977892751e201,
with Cstart=3.123143830205975e28 and m_hi=2.577600186401461e12.
Its known certificate failure is not a lower bound on the actual solution.
Nothing here makes those old bounds applicable to the changed theory.
The useful improvements are the exact product trace and local tensor
cancellations, the actual horizon comparison, and ultimately a direct
bound on the correlated correction. Additional generic derivative
chains and unit changes do not meet (13).

Static current/constraint recovery and the mass boundary cancellation
remain restricted results on the retained phase. Their nonlinear
domain hypotheses must survive the new seed and ball. The extra R^2
scalar is part of this specified theory, not a measured EFT coefficient.
Omitted interactions, stress fluctuations and the retarded dynamics
still need independent bounds. Euclidean positivity is not all-sector
Lorentzian stability; the invariant product has no opening null stress.
No new physical obstacle is proved here that warrants abandoning the
construction. The mandatory order is now product-centered full-source
and direct response, then the preserved nonlinear ball and exact matching,
then validity/stability and the full opening/transport/creation tests.

Remaining to-do list: complete direct cold and combined outer admission
response, including global mass/kernel moments and inverse/domain
transfer; nonlinear solution and exact compatibility/matching;
semiclassical validity, all-sector stability, charge/return, opening,
transport, preparation, separate creation and full-result audit.
