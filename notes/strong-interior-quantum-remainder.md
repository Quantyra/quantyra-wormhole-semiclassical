# A strong interior estimate for the actual quantum remainder

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof and
verification contract before execution. **Component estimate; the
global fast TT inverse remains open.**

The difference between the actual canonical massive Hessian and its
positive geodesic comparison has a numerical strong third-derivative
bound on the compact region, including the regular bolt:

    ||1_Omega (H_mass-Q_geo)T||_2
      <= 1e48 C sum_(j=0)^3 ||nabla^j T||_2.          (1)

Here Omega={Y>=1/8}, with Y the inward canonical proper distance
from the outer seam. Norms and covariant derivatives use the
dimensionless canonical orbit gbar, on one physical half, and include
the time integral. C=B^2/ell_o^2 is the retained sphere/orbit ratio.
H_mass is the Hessian per fixed gbar area of the COMPLETE canonical
massive action, without the normalization 2G/B^2. Q_geo is exactly
the operator of the quadratic form (2) in
[the curved TT proof](stationary-curved-quantum-tt.md).
No finite contact, mode, thermal contribution or return is omitted.

Initially T is smooth, in the original paired domain, and compact
away from the seam. The difference in (1) extends continuously to
the completion in the displayed H3 norm. This does not separately
extend H_mass or Q_geo to H3. The independent Weyl action, LLL,
reference potential, TT reconstruction and outer weights are not
part of (1). Its deliberately coarse constant is not a Newton
contraction estimate.

## 1. Extracting a first variation from the retained point proof

Use the actual interior normal conformal charts of
[S145](general-interior-dirac-response.md), including the cold bolt
charts. They have radius 1/32, sigma(0)=D sigma(0)=0 and the
complex seed bound |exp(sigma)-1|<.026. The restriction Y>=1/8
keeps these charts inside the physical interior.

For a real C5 conformal direction v, let U be S145's norm: the
maximum of the actual and capped global value/gradient norms and
the coordinate derivatives through five in the chart. At background
zero, the SAME one-parameter complex disk |s|<=1/(200U) gives:

    |E_cap| < 1e17 M^-3,
    |E_actual-E_cap| < 7e11 M exp(-M/256),
    |T_a2| < 1e17 C.                                (2)

These are the holomorphic scalar-source bounds proved in sections
4-6 there, before applying its second-derivative Cauchy estimate.
Using first-derivative Cauchy instead, and the retained full tower
bounds sum d/M^3<2C/3 and sum d M exp(-M/256)<=24C 256^3, gives

    |D T_mass[0](0)[v]| < 6e22 C U.                 (3)

Indeed the coefficient is

    (2/3)2e19 + (1.4e14)24(256^3) + 2e19 < 6e22.

This is the UNNORMALIZED scalar source in gbar units. Its flat
vacuum value has zero scalar variation; it has not been subtracted
from the action. The fixed-area Hessian is instead

    H_mass v = D T_mass[v] + 2v T_mass.

For the linear and quadratic probes below v(0)=0, so (3) bounds
the actual Hessian on those probes. We do not apply (3) to a
constant probe or omit the density term for a general input.
The constant probe is supplied independently in section 4.

The complex Bessel envelope underlying (2) uses
[DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9), inspected
2026-10-20. This primary representation supports that envelope;
it does not supply the geometric or operator estimates here.

## 2. Fixed probes and their complete geodesic action

Write P(s)=I_s(6,6)=462s^6-1980s^7+3465s^8-3080s^9
+1386s^10-252s^11. In each normalized chart choose chi=1 for
|z|<=1/128, chi=0 for |z|>=1/64, and
chi=1-P(128|z|-1) between. Use the five probes

    chi z1, chi z2, chi z1^2/2, chi z2^2/2, chi z1 z2. (4)

Their value vanishes at the center. Their first and second
covariant jets there are the coordinate polynomial jets, since
the Christoffel symbols vanish at that point. In particular they
extract all two first-order and three symmetric second-order
coefficients, with no unknown change of frame.

Here are reproducible generous derivative reserves. The derivatives
1..6 of P have Bernstein bounds
(11,220,3960,63360,887040,10644480).
For the Euclidean radius use the Cartesian tensor component bounds
(1,2,6,36,300,4000) times |z|^(1-j). If d_k denotes the first
array and r_j the second, the cutoff's j-th Cartesian derivative
is at most 64^j C_j, where

    C_0=1,
    C_j=2^j sum_(k=1)^j d_k Bell_(j,k)(r_1,...).

Leibniz gives for the linear and quadratic probes respectively

    64^(j-1)(C_j+j C_(j-1)),
    64^(j-2)(C_j+2j C_(j-1)+j(j-1)C_(j-2)).

These component bounds hold globally, including the inner polynomial
piece, by taking their maximum with 2. Thus every probe has U<1e18
and coordinate derivatives through six below 1e22. The probes are
C5 with bounded piecewise sixth derivatives; smooth approximations
preserve these strict bounds. The cap used by S145 is identically
one on their support.

For completeness we quantify their derivatives along unit geodesics.
On |z|<=1/64, the seed holomorphic bound, log(1+v) and Cauchy on
the remaining radius give a majorant for EACH Christoffel component

    sum_(|alpha|=j) |D^alpha Gamma|/alpha! <=
      100 (j+1) 128^j, j<=4.

This follows also directly from
Gamma^k_ij=delta^k_j sigma_i+delta^k_i sigma_j-delta_ij sigma_k:
|sigma|<.027 and the real-coordinate conversion costs at most 2
per derivative. The spare factor 100 covers that conversion and
all three summands. Unit initial coordinate velocities are <2.

The positive scalar formal solution

    y(0)=0, y'(0)=2,
    y''=400(1-128y)^-2 (y')^2                         (5)

therefore majorizes every coordinate derivative of a geodesic
through order six. The factor four sums its two velocity indices.
Coefficient induction in the geodesic equation proves this finite
jet statement; no spatial analyticity of the probe is assumed.
For a probe q, multivariate chain rule then gives

    |X^n q| <= 1e22 sum_(k=1)^n 2^k Bell_(n,k)
                                  (y'(0),...,y^(n-k+1)(0)).

The rational coefficients of (5) give, conservatively,

    sup|X^3 q| < 1e32,   sup|X^4 q|,sup|X^6 q| < 1e46. (6)

Outside its support the probe and its derivatives vanish.
These are full-flow bounds even for geodesics leaving and returning
to the support. No artificial compact-cap geodesic is substituted.

Section 9 of the curved TT proof, using the exact sixth Bessel
moment, gives the full point action

    |Q_geo q| <= C[2e-5 sup|X^6 q|+sup|X^4 q|]
                 < 2e46 C.                         (7)

Its ordinary fourth restoration is included. Equation (3) gives
|H_mass q|<6e40 C. We will retain the convenient combined bound
|(H_mass-Q_geo)q|<3e46 C.

## 3. Strong bound for the singular remainder

Put E=H_mass-Q_geo and R=1/400. The curved TT proof's local
normal-cap/geodesic comparison bounds the separated error by
1e8 C/d(x,y)^4 in two dimensions. Its common finite-mode
fourth restoration agrees; the remaining local distribution has
order at most two, including its explicit curvature contact.
This fact is essential: a separated kernel estimate alone would
not identify or bound those contacts.

Retain the near kernel for d<R. In exponential polar coordinates
the Jacobian is at most 2r. Remove the covariant Taylor polynomial
through degree two and define the absolutely convergent operator

    T_3 u(x)=integral_(d<R) K(x,y) [
      u(y)-u(x)-r nabla_v u(x)-(r^2/2)nabla^2_vv u(x)]dy. (8)

The integral Taylor remainder is
r^3/2 integral_0^1 (1-t)^2 X^3u(Phi_(tr)(x,v))dt.
The complete geodesic flow preserves Liouville measure. Minkowski,
the angular measure 2pi, and the radius integral therefore give

    ||1_Omega T_3u||_2 <= 2e6 C ||nabla^3u||_2.       (9)

The reserve includes the even-double/physical-half conversion.
In particular this is an L2 strong estimate; it does not
differentiate the unsigned kernel envelope or infer a strong
estimate from a quadratic form. The pointwise version and (6) give
|T_3 q(x)|<2e38 C for every probe.

Keep the actual gluing operator as a separate bounded kernel.
Its two row bounds are 3e24 C. Extending the near comparison
from 1/128 down to R costs less than 1e15 C in both rows,
by integrating 1e8 C/d^4. The full geodesic cosine tail has row
mass at most 2pi(18C/R^4)<1e13 C, including returns.
Consequently the remaining bounded kernel B has both rows
below 4e24 C. Replace it by B_0=B-(B1)Id:

    B_0 1=0, ||B_0||_(2->2)<=8e24 C,
    |B_0 q(x)|<=8e24 C sup|q|.                       (10)

The subtraction transfers its row sum to the constant contact;
it does not remove a physical contribution. All polynomial
geodesic restorations remain in the local distribution.

## 4. Fixing every local coefficient and proving (1)

The actual constant-probe theorem gives
|Theta_p0|<=1e20 q^2(1+1/Y).
Conversion to the present area is p0=Theta_p0/C.
Since q/C<400 and Y>=1/8,

    |p0| < 2e26 C.                                  (11)

Q_geo 1=0. Equations (8) and (10) also annihilate constants.
Thus this is the constant coefficient of E itself, including
the actual tadpole, vacuum, density and row-sum contacts.

At each point, the exact local representation is

    E u=p0 u + a1 u_;1+a2 u_;2
                 +a11 u_;11+a22 u_;22+a12 u_;12
                 +T_3u+B_0u.                       (12)

It follows from the retained degree-two extension just described;
the Hessian components are in an orthonormal center frame.
Evaluating on (4) fixes each coefficient uniquely:
it is E q-T_3q-B_0q at the observation point.
Equations (7)-(10) imply |a_i|,|a_ij|<4e46 C.
In particular the first-order contacts are bounded by their
actual point action. Symmetry has not been used to differentiate
an unproved bound on the second-order coefficients.

Summing the five terms in (12), (9)-(11), and the global input
norms proves (1), with ample reserve in 1e48. The bolt is an
ordinary interior point throughout this argument: neither the
singular polar radial frame nor a new bolt boundary condition
was used.

## 5. Verification and use

Freeze this proof and its checker before execution. Check the
complete first-Cauchy tower budget; the fixed-area density
conversion on zero-value probes; every Bernstein cutoff and
Cartesian radial derivative bound through six; the exact
geodesic majorant recurrence and chain-rule reserves; independent
Taylor-integral identities on polynomials; the radial integrals,
Schur subtraction and final coefficient aggregation. Pin all
predecessor receipts and source hashes. Preserve every run and
failure. Arithmetic checks support, rather than replace, the
ordinary operator and domain arguments above.

This closes a strong interior remainder estimate for the actual
canonical massive response. It does not prove ellipticity of
Q_geo in a numerical Fourier norm, control the weighted outer
fourth/log remainder, or complete the TT inverse.
The sharper accepted form bounds remain the bounds to use for
low frequencies; (1) is an unoptimized strong remainder reserve.

Remaining to-do list: quantify the outer fast TT estimate and
the principal graph comparison, higher York reconstruction and
full inverse/residual composition; nonlinear signed compatibility
and the exterior, opening, stability, transport, preparation and
creation requirements.
