# Local sources in the regular centered response

2026-10-16, S139/S074. Ordinary supporting argument before focused
verification. This note specifies the local part of
the full-response estimate. Its metric-jet hypotheses must be proved
for the actual completed family before applying the estimate.

Use the fixed action split, scale, product radii B,A, and original
prepared temperature of
[the centered seed](stationary-centered-covariant-seed.md).
Put C=B^2/A^2. The product frame is the sphere of radius B, with
the same fixed common scale mu0. Thus the Weyl field below is
phi=log(r/B); changing this reference does not reset mu0.

## Local variational formula and derivative counting

The exact formulas (1)-(3) of
[the anomaly sources](spherical-anomaly-sources.md) apply with their
fixed reference length equal to B. In their notation write
R=R[hat g], U=Box phi, V=(nabla phi)^2, y=B^-2, and no orbit
electric field. The coefficients are c=1/20, a=11/360, b=1/30,
d=b-2a/3. In particular

```
P=phi[2c(R+2y)/3-4ay]+b(U+V),
Q=bR-6bU-6dV,
Z=bR+(2b-4a)y-6d(U+V),
M_ab=-hat g_ab L/2+(R/2)hat g_ab P
     +(hat g_ab Box-nabla_a nabla_b)P
     +Z phi_a phi_b-phi_(a nabla_b)Q
     +hat g_ab[(nabla Q).(nabla phi)+QU]/2,
E_phi=L_phi+Box Q-2nabla_a(Z nabla^a phi).
```

The local stress is -B^2 M_ab/(2pi). Its independent radius
source follows from E_phi and the joint trace identity. These are
compact bulk variations. The finite-cut surface terms and boundary
Hamiltonian remain separate.

Here is a sufficient quantitative jet contract. At an observation
choose a real conformal chart for the product metric and normalize
its value to the Euclidean metric. Let omega be its conformal factor.
For total orders j<=4, background jets of omega and differentiated
phi are bounded by j!100^j; allow |phi|<=L0. A parameter direction
has jets bounded by epsilon W j!100^j for both fields. Use the
same bounds along the entire real parameter interval, including
the chart transport. In a rescaled chart of radius d, these are
jets with respect to the rescaled coordinates: an order-j physical
derivative has an additional factor d^-j. Assume L0>=1, d<=1.

In such a chart

```
A^2 R=-2 exp(-2omega)(omega_11+omega_22),
A^2 U=exp(-2omega)(phi_11+phi_22),
A^2 V=exp(-2omega)(phi_1^2+phi_2^2).
```

The displayed sources have total derivative order at most four.
There is at most one undifferentiated phi in each term. Curvature
is replaced by the displayed expression before counting derivatives.
The exponential inverse-metric factors have powers at most four;
their value at the observation is one. Their parameter derivatives
are included, rather than treating mixed stress components as fixed
coordinate components.

For an explicit conservative coefficient count, the expanded L has
coefficient norm below 60; P and each of Q,Z below 10. Each of two
coordinate differentiations of a P monomial increases its coefficient
norm by at most six, including its exponential factor. Tracing over
two coordinates and adding the connection terms gives a bound below
2000 for its complete Hessian contribution. The R P term, L term,
Z-gradient term and Q terms together cost below 1000. The analogous
E_phi expansion is below 3000. Enlarging both to 1e5 covers all real
components and the orthonormal stress conversion. Differentiating
one parameter, including inverse-metric factors and frame conversion,
has coefficient norm below 1e7.

A product of jets whose orders sum to at most four is bounded by
4!100^4=2.4e9. Consequently the following deliberately enlarged
bound applies to the dimensionless mixed stress and E_phi, with
the algebraic magnetic term separated:

```
|delta(B^2 rho_hat,W)|, |delta(B^2 p_hat,W)|,
|delta(B^4 exp(4phi) A4,W)|
 <=1e20 epsilon C^2 W L0 d^-4
   +10 epsilon q^2 W L0.                           (1)
```

The magnetic term is algebraic in phi and contains no metric
derivatives; it does not acquire d^-4. Terms carrying y or y^2
obey (1) because C>1. Equation (1) also bounds the independent
angular combination after its stated physical-radius factors are
included. The trace of the full Dirac action, rather than the Weyl
piece alone, is used when reconstructing the full angular pressure.

## Cold patch and reference action

For the cold oscillation estimate assume the same direction and
background jet contract through order five. One further spatial
derivative has coefficient norm below 1e8 in the preceding expanded
parameter variation. Its jet product is at most 5!100^5=1.2e12.
The enlarged constant 1e22 therefore bounds the first spatial
derivative of the local dimensionless sources. This uses a fifth
metric/radius jet, not a sixth one.

In the actual cold patch use s=r0/B-1<=.001 and
L_c(s)=1+log(.01/s). A completed direction with the cold contract
W=s L_c(s), d=1 has |phi|<1, so L0=1 in (1). Thus no second
cold logarithm is introduced by the local functional. At the bolt
s=h and the result retains the thermal-gap factor h. Normalized
regular radial jets are meant here; an unweighted areal-coordinate
derivative is not substituted for them.
The value .01 in L_c is the auxiliary cap seam; the compact
non-bolt estimate supplies physical coverage starting at s=.001.

For 0<=chi<=1 the weight at any point on the connecting radial
segment is at most 1.62 h L_c(h), and proper length is below
1.1 chi in the original cold chart. The local response consequently
satisfies

```
|delta H_local(chi)-delta H_local(0)|
 <=2e22 epsilon C^2 h L_c(h) sqrt(chi)
   +1e12 epsilon q^2 h L_c(h) sqrt(chi).             (1a)
```

Here H_local denotes either normalized mixed stress component, their
sum rho+p, or the normalized anomaly combination in (1), with its smooth physical
radius factors. Regular rotational symmetry sets rho+p to zero at
the bolt. The stronger Lipschitz factor chi was enlarged to sqrt(chi)
only to use the same integrable Ward budget as the canonical source.

The reference mixed pressures are the fixed rational functions in
[the complete seed source](stationary-centered-exterior-source.md),
with |p0A|<3300q^2/B^4, |pv|<q^2/B^4 and |DP|<3q^2/B^4.
Direct differentiation of their powers t^-2 through t^-4, and of
the corresponding pS profile, costs less than 1e8 through order
four. Therefore their parameter response, including the physical
radius factors, fits 1e10 epsilon q^2 W in the dimensionless
source norm. At the bolt this also has the factor h.

The position-independent common a0 restoration has zero metric
response at fixed B, masses and mu0. A frozen product trace may also
be subtracted inside the Ward reconstruction; its constant kernel
has zero response. This does not set the curvature response of the
full homogeneous-product determinant to zero: that response remains
inside the canonical estimate. The physical radius-factor response
of the frozen reference must be combined with the reference action
at the same moving radius. Since
delta log(r/B)=O(epsilon s L_c), this operation cannot create a
constant product error. The fixed B,A and reference potential are
never varied as independent couplings.

The canonical LLL has trace proportional to q R[hat g]. Its
parameter response is a local second-order curvature expression;
its HH stress is reconstructed with the same Ward integral as the
massive tower. Its trace and first spatial response fit
1e10 epsilon q C W in the above chart contract. With q/C<400,
this is below 4e12 epsilon C^2 W and can be retained separately
inside the complete canonical response budget.

## Outer cubic cancellation and integrable remainder

Let Y be proper product-frame distance to the outer endpoint.
For the outer estimate the chart is explicitly rescaled by Y.
The background phi behaves as -log Y; its ordinary derivatives
are not uniformly bounded. The contract in that rescaled chart
has d=Y, L0=1+|log Y|. After subtracting the correlated mass cubic,
the actual directional remainder must have

```
|partial_Y^j delta(omega,phi)_rem|
 <=epsilon j!100^j Y^(4-j)(1+|log Y|), j<=4.         (2)
```

Thus W=Y^4(1+|log Y|) in the rescaled chart. The Y^-4 in (1)
cancels this fourth-order weight. Multiplication by the physical
factor (B/r)^2 gives a density/pressure remainder bounded by a
constant times

```
epsilon C^2 (1+log^2 t)/t^2, t=r/B.                 (3)
```

The conversion between Y and 1/t and its finite derivatives must
be included in the actual metric-transfer lemma. An upper bound
for a remainder is not differentiated to obtain (2).

The removed cubic is fixed by the actual geometric mass m_lambda.
The regular completion preserves its outer clock and has no
r^-3 lapse coefficient, as proved in
[the finite completion](stationary-centered-first-order-geometry.md).
Consequently its outer jets have

```
F_lambda=lambda_e r^2+1-2m_lambda/r+O(r^-2 log r),
alpha_lambda=O(r^-4 log r),
c_lambda=-m_lambda sqrt(lambda_e)/C^(3/2).
```

These statements include the derivative control (2), rather than
only value asymptotics. Substitution in the exact radial Euler
formula of [the Weyl sources](stationary-weyl-radial-sources.md)
gives

```
rho2,W=11lambda_e^2 r^2/(240pi)
       -lambda_e m_lambda/(5pi r)
       +O(epsilon C^2(1+log^2 t)/(B^2t^2))
```

for the parameter difference after retaining the unchanged seed
remainder. Terms involving alpha start at the indicated integrable
order because its cubic coefficient is zero. The matched canonical
boundary layer has the opposite +lambda_e m_lambda/(5pi r).
The universal singular trace is therefore removed only after these
two sources are added. The reference vacuum term is the same
11lambda_e^2 r^2/(240pi) for the whole family and cancels exactly.

For a finite coefficient bound on the mass-cubic part, subtract the
explicit RN mass variation with the same m_lambda before using
(2). Its radial Weyl formulas contain only a quadratic polynomial
in m_lambda and one logarithm; the coefficient counts 2076 and
3730 in the accepted seed-source proof bound its integrable
remainder. The residual local jets then obey (1)-(3). This procedure
keeps the singular coefficient correlated and assigns neither
component a separate divergent mass integral.

## Allocation in the complete response

The source is a function of the actual metric and its HH domain.
Its fixed-coordinate response also includes the physical r^-2
factor, the mass-density Jacobian r_chi, and the lapse prefactor
r_chi/(rF). Each is included once in the complete source response;
none is charged again to the separately defined E0 residual.
Angular pressure follows from the complete anomaly and canonical
trace with this same allocation.

This local lemma does not prove the actual finite-family jet
contract, the massive determinant response, an inverse of the full
higher-derivative action, or a boundary Hamiltonian. Its role is to
make the finite local terms and outer cancellation explicit in
the complete residual argument.

Remaining to-do list: audit the coefficient and chart contracts;
prove them for the finite completion and combine the full canonical
response, source prefactors and E0 residual before verification.
