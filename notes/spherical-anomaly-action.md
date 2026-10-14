# Nonlinear spherical anomaly action and the joint source identity

S074, 2026-10-13. The four-dimensional Dirac anomaly fixes a local
nonlinear functional for a varying sphere radius. Its spherical form
contains explicit quadratic, cubic and quartic derivative terms and a
radius--electric contact. Its variation reproduces the full bulk trace
anomaly, including the derivative curvature term. This supplies a
previously unspecified part of the local spherical action and its
stress/dilaton/current matching. It does not determine the independent
state-dependent functional or solve the gravitational constraints.

## Common regulator and conformal frame

Use one massless unit-charge four-component Dirac field, mostly-plus
signature and the common four-dimensional zeta prescription of the
[uniform reference](uniform-dirac-reference.md). The background is

```
g4=g_ab dx^a dx^b+X dOmega_2^2,
X=b0^2 exp(2sigma), g_ab=exp(2sigma) h_ab,
g4=exp(2sigma) h4, h4=h_ab dx^a dx^b+b0^2 dOmega_2^2,
F=F2+(q/2) sin(theta) dtheta wedge dphi,
y=b0^-2, r=R[h], H=r+2y, v=(nabla sigma)^2, u=Box sigma. (1)
```

b0 is a fixed reference length; contractions and derivatives without a
subscript below use h. Keep the same gauge potential and magnetic flux
under the Weyl transformation. Assume smooth positive X. Bulk identities
hold locally; for integrated equalities take compactly supported Weyl
variations, or retain the stated surface terms. A finite-mouth boundary
condition and its boundary anomaly are not chosen by these assumptions.

The earlier [mode reduction](dirac-spherical-source-reduction.md) uses
Psi=X^-1/2 chi. For the covariant spinor norm this is an isometry:
the X in the four-dimensional volume cancels X^-1 from the two spinors.
An orthonormal angular expansion and this isometry conjugate the full
regulated operator. They do not justify multiplying separately
renormalized two-dimensional determinants. Keep their summed heat trace
and one four-dimensional subtraction before removing the regulator.

The distinct Weyl transformation of a massless Dirac field is
Psi_g=exp(-3sigma/2)Psi_h. Together with (1), it gives
chi_g=exp(-sigma/2)chi_h and converts each angular mass lambda_n/sqrt(X)
into lambda_n/b0 in the h frame. Its nontrivial regulated determinant
variation is the anomaly functional below. In this organization it
must be included once, not supplemented by an independently guessed
radius Jacobian. Other field redefinitions, boundary measures and
gravitational constraint determinants are not covered by this statement.

## Bulk anomaly and its finite curvature convention

Write the determinant's trace as

```
16pi^2 A4=c C4^2-a E4+b Box4 R4-(2/3)F4^2,
c=1/20, a=11/360, b=1/30.                             (2)
```

The c,a and gauge coefficients are the retained
[four-dimensional trace result](throat-trace-anomaly.md). For the local
derivative convention used here, the universal heat coefficient gives
tr[60 Box E+12 Box R I]=(-60+48)Box R=-12 Box R, since tr E=-R.
Multiplying by the fermion anomaly sign -1 and dividing by 360 fixes
b=1/30. This follows directly from
[Vassilevich (3.27), (4.28) and section 7.1](https://arxiv.org/html/hep-th/0306138v3).
It specifies the bulk local-heat representative, not a physical
determination of finite curvature-squared couplings.

An added Lorentzian term cR integral sqrt(-g4) R4^2 shifts
b by -192pi^2 cR. In the preceding unredefined bosonic span,
cR=c1+(c2+c3)/3, after separating Weyl-square and Euler terms.
The derivation below retains arbitrary b until the displayed numerical
coefficients are substituted, so this finite ambiguity stays explicit.
Other finite local EFT terms have their own trace variations and must
be added separately to the determinant's identity (2).

## Spherical reduction of the four-dimensional transgression

Define Q4=E4-(2/3)Box4 R4. The standard Weyl identity and operator are

```
sqrt(-g4) Q4[g4]=sqrt(-h4)[Q4[h4]+4 P4[h4]sigma],
P4=Box4^2+2 Ric^(mu nu) nabla_mu nabla_nu
                  -(2/3)R4 Box4+(1/3)(nabla R4).nabla.
```

[Mazur--Mottola (40), (44), (46) and (56)](https://arxiv.org/html/hep-th/0106151v1)
give these established bulk identities and the Wess--Zumino functional.
Our reduction and source normalization below use their bulk construction;
their discarded surface terms do not fix this problem's boundary ensemble.
Set d=b-2a/3. On h4,

```
C4^2=H^2/3, E4=4yr, Q4=4yr-(2/3)Box r,
F4^2=F2^2+q^2 y^2/2=:Ih,
P4 sigma=Box^2 sigma+nabla_a[(r-4y)nabla^a sigma/3]=:P sigma.
```

The full spherical bulk transgression, before extra integrations by
parts, is

```
W[sigma;h,A]=b0^2/(4pi) integral sqrt(-h) Lraw,
Lraw=sigma[c H^2/3-a(4yr-(2/3)Box r)-(2/3)Ih]
     -2a sigma P sigma
     -d[(H-6u-6v)^2-H^2]/12.                         (3)
```

This is the anomaly-controlled difference of the Dirac effective
actions at g4 and h4, at fixed scale and with conformally related state
and domain data where that comparison is defined. In an in-in problem
the same local functional occurs on the two contour branches; the
initial-state and boundary terms still have to be matched. The anomaly
does not construct those data.

Integrating by parts gives the useful vertex generator

```
Ahat=c H^2/3-4ayr+b Box r-(2/3)Ih,
Lpoly=sigma Ahat-3b u^2+[br+(2b-4a)y]v-6d u v-3d v^2,
Lraw=Lpoly+nabla_a B^a,
B^a=d[H nabla^a sigma-sigma nabla^a H]
    -2a[sigma nabla^a u-u nabla^a sigma
                         +(r-4y)sigma nabla^a sigma/3]. (4)
```

For the local heat representative of (2),

```
Lpoly=sigma[H^2/60-11yr/90+Box r/30-(2/3)Ih]
      -u^2/10+(r/30-y/18)v-7u v/90-7v^2/180.         (5)
```

Thus the anomaly sector is exactly quartic in sigma at fixed h and A;
there is no fifth-order sigma term in this functional. This does not
make the full physical action quartic in delta X: sigma=log(X/b0^2)/2,
and the remaining determinant depends on h=X^-1 b0^2 g. Neither its
nonlocal vertices nor the even gravitational constraints have been
eliminated. Equation (5) is a definite component of that calculation.

## Direct variation and the trace check

For arbitrary b, put K=br+(2b-4a)y and

```
U=-6b u-6d v, V=K-6d u-6d v,
E_sigma=Ahat+Box U-2nabla_a(V nabla^a sigma),
Theta_sigma^a=U nabla^a(delta sigma)
                 -nabla^a U delta sigma+2V nabla^a sigma delta sigma.
delta Lpoly=E_sigma delta sigma+nabla_a Theta_sigma^a. (6)
```

The checker compares this Euler derivative with the four-dimensional
anomaly independently of the Q4 transgression formula. Under (1),
let T=H-6u-6v. The physical curvature satisfies

```
exp(2sigma) R4[g4]=T,
Ric_ab[g4]=(r/2-u-2v)h_ab-2nabla_a nabla_b sigma
                                      +2sigma_a sigma_b,
Ric_ij[g4]=(y-u-2v)h_ij,
exp(4sigma) Box4 R4[g4]=Box T-2nabla sigma.nabla T-2uT.
```

Use E4=C4^2-2Ric^2+(2/3)R4^2 and the two-dimensional Bochner identity
Box v=2|Hess sigma|^2+2nabla sigma.nabla u+r v. Exact algebra then gives

```
E_sigma=16pi^2 exp(4sigma) A4[g4].                    (7)
```

The identity holds for arbitrary local two-dimensional curvature and
sigma jets, not only for static or homogeneous backgrounds. Direct
four-dimensional metric-jet contractions separately check the curvature
relations used here. At constant sigma and r=-2x, (7) reduces to the
already accepted magnetic-product trace, including its negative gauge
term. This is a regression against the full-tower reference, not an LLL
trace substituted for a four-dimensional anomaly.

## Current, stress/dilaton combination and boundary obligations

The radius--electric contact in (3)-(5) is

```
W_gauge=-b0^2/(6pi) integral sqrt(-h) sigma F2^2,
j_W^b=2b0^2/(3pi) nabla_a(sigma F2^(ab)).             (8)
```

Its current is identically conserved. The associated h-frame stress
and sigma source obey the diffeomorphism identity including the Lorentz
force. In particular, omitting this contact while retaining the varying
radius would omit part of the required current/source matching. Its
constant-sigma shift of the inverse Maxwell coefficient is
sigma/(6pi^2), consistent with the scale coefficient in the
[common electric matching](curved-electric-matching.md). It is not an
additional copy of that loop; both come from the same determinant.

For the physical spherical sources defined in the earlier reduction,

```
T_ab^(2)=-2/sqrt(-g2) delta Gamma_D/delta g2^(ab),
J_X=4G/sqrt(-g2) delta Gamma_D/delta X,
j^a=1/sqrt(-g2) delta Gamma_D/delta A_a,
T^(2)a_a+X J_X/(2G)=4pi X A4[g4],
nabla^a T_ab^(2)=J_X partial_b X/(4G)+F_ba j^a,
nabla_a j^a=0.                                      (9)
```

The trace identity follows by varying sigma at fixed h and A:
delta g2_ab=2delta sigma g2_ab and delta X=2delta sigma X. Equations
(3)-(7) supply its local anomalous part. The full source still includes
variation of Gamma_D[h4,A] with respect to h and the matched state.
In particular, (9) fixes a combination of stress and dilaton source;
it does not determine each separately or justify zeroing J_X. Gauge and
diffeomorphism covariance supply the other two identities. Extra finite
local terms require their corresponding sources on both sides.

For a nonzero boundary variation, (4) retains the integration-by-parts
flux, and (6) displays the sigma boundary variation of Lpoly. Fixing
both delta sigma and its normal derivative is one way to kill (6), but
is not imposed as the physical mouth ensemble. A driven boundary needs
the appropriate boundary action/anomaly, joint metric and Maxwell
variations, and actual matching data. These bulk formulae cannot select
them. No new propagating conformal mode or stability claim is inferred
from the fourth derivatives in this component of the effective action.

## Verification and consequence

[The checker](../checks/spherical_anomaly_action.py) verifies the heat
coefficient, finite R^2 shift, reduced polynomial and its retained
divergence, the arbitrary-jet trace identity, direct four-dimensional
curvature contractions, gauge/source Ward identity, and constant-product
regression. Independent compact-profile action variation checks (6)-(7)
without a surface contribution. Numerical checks supplement the exact
identities; they do not bound the unknown nonlocal functional.

This gives explicit nonlinear local terms and a common source identity
needed by S074.3-S074.5. It organizes the radius-rescaling anomaly
without declaring the complete measure or interacting action finished.
The boundary ensemble, independent matter operators, nonlinear lapse,
shift and Gauss reconstruction, general state-dependent response and
the opening error ledger remain unfinished. Closed predecessors remain
closed.

Remaining to-do list: complete the other local/matter and boundary terms,
nonlinear constraints and state-dependent source matching; then S075-S080
state/drive, stress/geometry, stability, transport, preparation/creation
and final audit.
