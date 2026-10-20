# Interior massive response without zero value or first-jet restrictions

2026-10-19, E003/S145 supporting S139/S074. Ordinary argument and
verification protocol before execution.

The interior canonical massive SCALAR source has a uniform mixed
second-variation bound for general real-C5 conformal directions.
Neither direction is required to vanish at the observation point.
The proof keeps the coordinates fixed and resums the point value
into the reference mass. Linear-gradient insertions supply the
additional second-, third- and fourth-Born local terms.

## 1. Fixed setting, normalization and statement

Use the interior seed charts, cap, paired domain, prepared period,
angular multiplicities and common renormalization of
[S144](coincident-normal-dirac-response.md). In particular

    M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18,
    eta_G=GC/B^2<8.5e-33.

The seed chart has sigma_0(0)=grad sigma_0(0)=0 and radius 1/32.
The analytic seed bounds and the fixed radial cap, one through
1/64 and zero by 1/32, are unchanged. Directions and backgrounds
preserve the original identical paired domain and seam condition.
This is a conformal-shape estimate on interior charts, not a theorem
for arbitrary radius or nonspherical metric variations.

Use exactly S144's norm: U is the maximum of the actual global
value/gradient norms, capped value/gradient norms, and coordinate
derivatives through order five in the uncapped chart. Define W for w.
Let the real background b have norm at most 1e-3. REMOVE the
zero-value and zero-first-jet conditions on all three fields.

Let T_hat[b](0) be the canonical massive action source per deformed
metric area, normalized by 2G/B^2. Equivalently, if F_hat is the
source per fixed seed area,

    T_hat[b](0)=exp(-2b(0)) F_hat[b](0).                (1)

The result is

    |D^2 T_hat[b](0)[u,w]| <2e-7 U W.                 (2)

It includes the full massive angular tower, actual exterior-domain
comparison and the common a2 contact. Its fixed-scheme vacuum scalar
source is constant and has zero variation. The independent 4D Weyl
functional, LLL and radius response remain separate.

Integrating (2) gives the same Lipschitz bound for D T_hat along
any real segment staying in the background ball. Equation (2) is
not a bound for F_hat: differentiating (1) backwards retains
4u(0)w(0)T_hat and the two first-source/density cross terms.
It is also not yet the norm required by the global exterior inverse.

## 2. Fixed coordinates, variable reference mass and coefficient bounds

Take complex parameters only:

    phi=sigma_0+b+s u+t w,
    |s|<=1/(200U), |t|<=1/(200W),
    phi_c=phi(0), m=M exp(phi_c),
    h(z)=exp(phi_cap(z)-phi_c)-1.                     (3)

The cap equals the actual field near zero. Now h(0)=0, but Dh(0)
can be nonzero. Outside the cap h is a small constant; it need not
be compactly supported. Its value remains bounded, and the loop
estimates below apply on the whole plane.

The parameter field zeta=b+s u+t w has every local coordinate
derivative through order five at most .011. Its actual and cap
value/gradient norms obey the same bound. Thus

    Re m>.98M, |m|<1.02M, |m|/Re m<1.02,
    |exp(-2phi_c)|<1.03.                             (4)

Use the following global Taylor bounds, where h_j is the homogeneous
Taylor part of degree j at zero and r=|z|:

    |h|<=.06, |h|<=2r,
    |h_1|<=.02r, |h_2|<=9r^2,
    |h_3|<=2e4r^3, |h_4|<=1e6r^4,
    |h-T1h|<=100r^2,
    |h-T2h|<=2e4r^3,
    |h-T3h|<=1e8r^4,
    |h-T4h|<=1e12r^5.                               (5)

Here is a direct verification of these generous reserves. In the
inner disk zeta-zeta(0) has modulus below .016r, ray derivatives
through order five below .063, and zero point value. The seed's
mass potential has the bounds 8r^2 and 1e4r^3. Multiplying it by
exp(zeta-zeta(0)) gives the displayed Taylor coefficients and inner
remainders by the same finite product/Bell calculation as S144.
For the quadratic remainder one may use 10r^2 in the inner disk.
The fifth seed remainder is below .052*32^5 r^5; the degree-five
product reserve 1e8 remains sufficient.

On the transition, write
phi_cap-phi_c=chi[sigma_0+zeta-zeta(0)]-(1-chi)zeta(0).
The seed logarithm has modulus at most 8r^2/(1-.026).
For r>=1/64 this bounds |h|/r by

    exp(.05)[8/(32*(1-.026))+.016+.011*64]<2.

Also |phi_cap-phi_c|<.05, so |h|<.06. For the quadratic
remainder on r>=1/64 the sharper bound
|h|<=exp(.05)[8r^2/(1-.026)+.016r+.011]
gives |h-h_1|/r^2<100. The remaining global remainders follow
by division at r=1/64:

    (2+.02)*64^2+9*64 <2e4,
    (2+.02)*64^3+9*64^2+2e4*64 <1e8,
    (2+.02)*64^4+9*64^3+2e4*64^2+1e6*64 <1e12.

Only real spatial C5 jets enter. No parameter-dependent spatial
coordinate is evaluated at a complex point.

## 3. All local terms through four derivatives

Put S_m=(Dflat+m)^-1. The mass-source difference from its locally
matched flat mass is

    F_cap-F_flat(m)
      =sum_(N>=1)(-1)^(N+1)m^(N+1)
         tr[S_m h S_m ... h S_m](0,0).               (6)

For h(0)=0 define the real-coordinate jet contractions

    l=grad h, H=Hess h, a=Delta h,
    t=grad Delta h, b4=Delta^2 h, all at zero.

The two-derivative local terms in the first and second Born orders
are -a/(12pi) and |l|^2/(12pi). There are no others at that order.
The complete fourth-order terms, with denominator 120pi m^2, are:

| Born order | Numerator |
| --- | --- |
| 1 | `-b4` |
| 2 | `4a^2+2 tr(H^2)+8 l.t` |
| 3 | `-18a |l|^2-16 l.H.l` |
| 4 | `21 |l|^4` |

Odd total Taylor degree integrates to zero. These identities concern
the action source -m tr R; no Majorana or extra mouth factor appears.

A reproducible derivation uses h(i partial_k) on the right-hand
resolvent in each ordered product. At m=1 set

    S(k)=(1-i gamma.k)/(1+|k|^2).

Apply the polynomial derivatives recursively and integrate the
resulting rational spin trace. For each even monomial,

    integral d^2k/(2pi)^2
      kx^(2i) ky^(2j)/(1+|k|^2)^p
      =Gamma(i+1/2)Gamma(j+1/2)Gamma(p-i-j-1)
          /[4pi^2 Gamma(p)].                         (7)

Odd monomials vanish. The checker evaluates all ordered compositions
needed through four insertions. Rotation/reflection invariance leaves
only l.t at linear/cubic order, a|l|^2 and l.H.l at linear/linear/
quadratic order, and |l|^4 at fourth-linear order. Independent
Cartesian cases determine these coefficients, including a mixed
xy^2 cubic. The normal-gradient limit reproduces S144's full Hessian
tensor. For complex m in (4), holomorphic continuation in Re m>0
gives the same coefficients and powers of m.

To identify the covariant answer, put

    A=Delta phi, p=grad phi, T=grad Delta phi,
    B4=Delta^2 phi,
    Q_phi=B4-3A^2-4p.T+4|p|^2 A.                    (8)

Since phi-phi_c=log(1+h),

    A=a-|l|^2,
    Q_phi=b4-4a^2-2tr(H^2)-8l.t
                  +18a|l|^2+16l.H.l-21|l|^4.

The table consequently sums to

    F_local=-A/(12pi)-Q_phi/(120pi m^2).              (9)

Direct variation of
-integral exp(-2phi)(Delta phi)^2/(240pi M^2)
gives exactly the fourth term of (9). Thus the new gradient terms
restore the same curvature action; no finite coupling is retuned.

## 4. A uniform, absolutely convergent remainder

The complex-mass kernel has the positive envelope

    ||S_m(z)|| <=1.02 (sqrt(2)/pi)
                       exp(-.98Mr/sqrt(2))/r.        (10)

Indeed |K_nu(mr)|<=K_nu((Re m)r) follows from
[DLMF 10.32.9](https://dlmf.nist.gov/10.32#E9), valid in the right
half-plane. The usual K0<=K1 bound then gives (10).
This primary representation was inspected on 2026-10-19.

At exponential weight M/10 the envelope's L1 norm is below 5/M
and its cubed L^(3/2) norm below 14/M. Hence for N>=2,

    integral L^j product_(N+1 edges) k
      <=14 j!10^j 5^(N-2)/M^(N+j-1).                 (11)

Here L is total loop length and |tr(product)| costs at most two.
The two-edge fifth moment is below 20/M^5.

Subtract every total Taylor degree at most four in N=1,...,4.
The following recursion bounds the product remainder of degree d+1:
let A0=2, p=(.02,9,2e4,1e6) and
B=(2,100,2e4,1e8,1e12), indexed from one. Define

    C(1,d)=B_(d+1), C(N,N-1)=A0^N,
    C(N,d)=B_(d-N+2) A0^(N-1)
             +sum_(a=1)^(d-N+1) p_a C(N-1,d-a).

This is repeated Taylor division of the first factor, with every
vertex distance bounded by L. It gives

    C(1,4)=1e12, C(2,4)=204180000,
    C(3,4)=82644, C(4,4)=808.0808.                   (12)

For N>=5 use the first five insertions' linear bounds 2^5 L^5
and .06 for every remaining insertion. With r_m=1.02, the
complete coefficient of M^-3 is bounded by

    20 r_m^2 C(1,4)
      +sum_(N=2)^4 [28*5!*10^5*5^(N-2) r_m^(N+1) C(N,4)]
      +28*5!*10^5*2^5*5^3*r_m^6/(1-5*r_m*.06)
      <8e16.                                       (13)

The denominator is positive. The underlying Neumann series also
converges as an L2 operator since |m| ||S_m|| ||h||<1.02*.06.
The moment bounds justify coincident evaluation, local subtraction
and locally uniform holomorphy. There is no discarded large-mass
tail or unestimated fifth-and-higher Born contribution.

After the scalar-area factor in (1), enlarge (13) to

    |E_scalar,cap|<1e17 M^-3.

Two-parameter Cauchy on (3) gives

    |D^2 E_scalar,cap[u,w]|<4e21 U W M^-3.            (14)

## 5. Actual exterior comparison and the full angular sum

The actual operator and the original-coordinate cap agree through
radius 1/64. Use S144's geometric-resolvent identity with cutoff
gradient below 256 and separation above 1/256. The fixed-coordinate
complex actual and cap resolvents retain the weighted bound 3/M and
the point-evaluation row 3.2e4 exp(-M distance/2).

For the SCALAR source, the external factor is M exp(-phi_c).
Its modulus is below 1.02M. Thus, with the additional density
factor included rather than frozen at one,

    |E_scalar,actual-E_scalar,cap|
       <7e11 M exp(-M/256).

Cauchy now costs 200^2, so its mixed derivative is bounded by

    2.8e16 U W M exp(-M/256).                        (15)

This uses the actual global domain throughout. No moving normalizer,
boundary image truncation or replacement of the physical temperature
is made.

The counting bound N(E)<=4CE^2 and M_1>18 yield

    sum d_n/M_n^3<2C/3,
    sum d_n M_n exp(-delta M_n)<=24C/delta^3.

As in S144, the latter bound discards the negative part of the
integration-by-parts weight before applying the counting estimate.
Equations (14),(15) therefore give the normalized remainder budget

    2 eta_G [(2/3)*4e21+2.8e16*24*256^3] U W.        (16)

Absolute uniform convergence permits both variations and real
finite-path integration to pass through the full angular sum.

## 6. Vacuum, anomaly and the common curvature contact

The output is the scalar source, not a coordinate-density source.
This matters for the locally matched mass m. In point splitting,
the spin trace of the flat kernel is m K0(mr)/pi. Multiplication
by the scalar-source factor -M exp(-phi_c) gives

    -M^2 K0(M exp(phi_c)r)/pi.

The proper separation is exp(phi_c)r+O(r^2). Subtracting the
covariant logarithmic singularity therefore leaves the same
flat-vacuum constant at fixed proper M and common scale, independent
of phi_c. The angular common finite part preserves this statement.
It does NOT set that vacuum constant to zero. Its scalar-source
variations vanish; coordinate-density variations do not.

The per-mode Polyakov allocation cancels the first term of (9)
after the scalar-area factor. This is the rank-two complex convention
already fixed in S142. The primary conformal Green-function and
proper-distance relations in
[Namuduri--Bilal, sections 4.5 and 5.2.1](https://arxiv.org/html/2308.05802#S5.SS2.SSS1)
were inspected on 2026-10-19. Their Majorana normalization and
small-mass expansion are not substituted for the present complex
Dirac normalization or the convergent large-gap construction.

Restore the common action, at the unchanged mu0,

    Gamma_a2=4alpha integral exp(-2phi)(Delta phi)^2,
    alpha=-C F_eff/(960pi), |F_eff|<100.

Its SCALAR source is

    T_a2=8alpha exp(-4phi) Q_phi.                    (17)

The formula retains all gradient, density and curvature terms.
For a simple bound on the same parameter bidisk, the seed's real
coordinate derivatives through order four at zero are below 1e8
by its analytic bidisk bound. Since the seed first jet vanishes,
|p_i|<=.011, |A|<2(1e8+1), |T_i|<2(1e8+1),
and |B4|<4(1e8+1). Equations (8),(17) give

    |Q_phi|<2e17, |T_a2|<1e17 C.

Cauchy and the physical normalization bound its mixed second
variation by 8e21 eta_G U W. Adding this to (16) gives a strict
bound below 2e-7 U W and proves (2). The deliberately loose
contact reserve is finite and uses the actual common coefficient;
it is not an unsigned sum of divergent inverse masses.

## 7. Verification protocol and remaining scope

Before execution require exact momentum-integral agreement with
all local coefficients through four derivatives, the logarithmic
jet conversion, direct covariant curvature-action variation,
Taylor-product recursion and complete rational budget.

Independently integrate the new second-, third- and fourth-Born
rational loop expressions by radial/angular quadrature at two
resolutions. Require error below 1e-10+1e-8*abs(target).
Also check complex-mass kernel majorants and the scalar/density
variation distinction using a graded finite-matrix example with
nonzero observation-point values. Require central mixed differences
at steps 1e-4 and 5e-5 to have error below 1e-8 and a refinement
ratio below .35, unless both errors are below 1e-45. These finite
diagnostics do not certify the continuum geometry.

The prior restriction on conformal value and first jet is removed
in the interior, without asking a C5 field for a sixth derivative.
The result remains a scalar-source C5 estimate. It does not give
the outer weighted response, full shape/radius coupling, Weyl/LLL
control or the complete nonlinear mapping in the exterior inverse
norm. The signed source and physical stability remain unresolved.

Remaining to-do list: outer finite-background weighted response and
global norm transfer; coupled exterior/state, charge and signed source;
driven geometry, opening, stability, transport, preparation and creation.
