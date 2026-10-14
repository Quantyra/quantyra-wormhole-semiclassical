# A finite local coefficient window for the magnetic reference

S074, 2026-10-13. The nearby homogeneous product theorem survives an
explicit nonzero range of the displayed bosonic four-derivative
coefficients. The range is sufficient and conditional; it does not
determine the coefficients of the physical theory. The proof retains
the full Dirac determinant in its accepted scheme and gives first-order
radius shifts with a quadratic remainder. It supplies a reference
geometry, not opening, stability or a global wormhole.

## Action, scheme and scope

Add the unredefined density in the [matter-aware local basis](matter-aware-local-basis.md)
to the [uniform Dirac reference](uniform-dirac-reference.md):

```
L4=c1 R^2+c2 Ric^2+c3 Riem^2
   +kappa^2(c4 R I+c5 Ric.X+c6 Riem.FF+c9 D.D)
   +kappa^4(c7 I^2+c8 I4),
I=Fbar^2, X_mu_nu=Fbar_mu_rho Fbar_nu^rho,
I4=Fbar_mu_nu Fbar^(nu rho) Fbar_rho_sigma Fbar^(sigma mu),
D^nu=nabla_mu Fbar^(mu nu), kappa^2=8pi G.              (1)
```

Fbar=F/g_R is canonical; its magnetic flux is 2pi q/g_R. The c_i are
additional finite coefficients at the fixed subtraction scale mu0=1/b0.
The determinant, g_R, Lambda and Newton coupling retain the same zeta
prescription as the preceding reference. The determinant is included
once, with its existing local subtraction; (1) adds independent finite
inputs, not another copy of the fermion ultraviolet terms.

[Cheung, Liu and Remmen, equation (4)](https://arxiv.org/html/1801.08546v3#S2)
provide the standard eight nonderivative terms, in conventions with
kappa^2 and kappa^4 absorbed into the corresponding coefficients. Their
source-free field redefinitions are not used here. We retain the original
metric, fixed magnetic flux and the ninth divergence term. The present
calculation derives its homogeneous pressures and a sufficient coefficient
region using the already accepted Dirac pressure bound; it makes no claim
to a new operator basis or to that source's positivity theorem.

The model being solved is S0+integral L4+Gamma_Dirac with the specified
invariant H2 spectral state. Independent Dirac operators, photon/graviton
loops, higher local operators, boundary matching and nonhomogeneous
response remain outside this truncation. In particular, the remainder
below is the nonlinear radius remainder of this specified model, not an
error bound for omitted EFT terms or for an order-changing field
redefinition. No physical coefficient values are assumed known.

## Homogeneous pressures from the complete displayed span

For AdS2(a) x S2(b) put x=a^-2, y=b^-2 and I=q^2/(2g_R^2 b^4).
The covariantly constant magnetic field gives D=0 and Fbar.dualFbar=0.
Direct constant-curvature contractions give

```
R=2(y-x), Ric^2=2(x^2+y^2), Riem^2=4(x^2+y^2),
Ric.X=y I, Riem.FF=2y I, I4=I^2/2.
Ag=4c1+2c2+4c3, K=2c4+c5+2c6, H=c7+c8/2,
L4mag=Ag(x^2+y^2)-8c1 xy
       +kappa^2 I(-2c4 x+K y)+kappa^4 H I^2.           (2)
```

For a local Lorentzian density L(x,y,I), fixed-flux radius variation
of a^2 b^2 L gives the invariant pressure eigenvalues

```
pA=L-x partial_x L,
pS=L-y partial_y L-2I partial_I L.
```

Consequently

```
pA_loc=Ag(y^2-x^2)+kappa^2 K y I+kappa^4 H I^2,
pS_loc=Ag(x^2-y^2)+2kappa^2 c4 x I
                         -2kappa^2 K y I-3kappa^4 H I^2. (3)
```

The first variations of D.D and J^2 vanish on this background. An
invariant current cannot have a nonzero vector value, and the constant
magnetic constitutive tensor is divergence-free. Thus the Maxwell
equation is also satisfied. The stress in (3) is blockwise constant,
so its covariant divergence vanishes. Symmetry leaves exactly the two
metric equations varied above. This reasoning concerns local bulk
equations; it is not a variational prescription for an exterior boundary.

The term xy has zero homogeneous pressure. The Euler direction
(c1,c2,c3)=(t,-4t,t) has Ag=0, while the isolated electromagnetic
J^2 direction (c7,c8)=(-cJ/8,cJ/4) has H=0. These regressions agree with
the prior audits. The latter direction still changes electric/axial
response, and c9 still generates a current contact under a matter-aware
field redefinition. Their absence from (3) is not a physical bound on
them. Boundary and topology terms are not discarded.

## The extended stationary equations

Keep the family and scales from the accepted reference:

```
q=k^4, g_R=k^-4, L4/l_p=k^5, Lambda=-3/L4^2,
C0=sqrt(1+12pi k^6), b0^2=L4^2(C0-1)/6, a0^2=b0^2/C0,
alpha=1/C0, z=b^2/b0^2, u=a^2/a0^2,
epsilon=4g_R^2(1+alpha), lambda=kappa^2/a0^2,
r=(C0/q)z/u, L=-log q+log z, k integer >=2000.          (4)
```

The coupling varies with k; this is the same family of theories as
before. Coefficients are fixed while varying radii in any one member.
For coefficients depending on k the condition below must hold for each
member. Let A(r,L), S(r,L) be the retained scaled Dirac pressures, with
A+S=-a2(r)/2 and
a2=(15-3/q^2-5r/q-3r^2)/(720pi^2). Define

```
g(z,u)=z^2 pA_loc/(a0^-4)
 =Ag(alpha^2-z^2/u^2)+2alpha(1+alpha)K/z+4(1+alpha)^2 H/z^2,
sigma(z,u)=(pA_loc+pS_loc)/(a0^-4)
 =4(1+alpha)c4/(u z^2)-2alpha(1+alpha)K/z^3
                                      -8(1+alpha)^2 H/z^4. (5)
```

The Einstein equations, in the preceding normalization, are

```
0=-alpha z-(1-alpha)z^2/2+(1+alpha)/2-epsilon A-lambda g,
0=z^2/u-(1-alpha)z^2/2-(1+alpha)/2
                    -epsilon S-lambda(z^2 sigma-g).   (6)
```

They are equivalent to

```
Tz=1-epsilon A-lambda g-(1-alpha)(z-1)^2/2,
Tu=[1-alpha+alpha/z-epsilon a2/(2z^2)+lambda sigma]^-1. (7)
```

The positive sign of lambda sigma in the denominator follows by adding
the two equations in (6); the checker verifies both identities directly.

## A coefficient window uniform in the candidate family

Define the nonnegative weighted seminorm

```
N4=8|Ag|+32|c4|+4|K|+128|H|,
lambda N4 <= 1/10000.                                 (8)
```

On the box |z-1|,|u-1|<=rho=1/100 and for 0<=alpha<=1/4, each of
|g|, |sigma|, |partial_z g|+|partial_u g| and
|partial_z sigma|+|partial_u sigma| is at most N4. Bound every positive
monomial using z,u in [99/100,101/100]. The exact rational coefficient
majorants for the four combinations are checked separately, including
the derivative sums; sampling the box is not the proof.

The accepted Dirac theorem supplies B=101(1+log q), |A|<=B,
|partial_z A|,|partial_u A|<=2B, |a2|,|r partial_r a2|<=1,
alpha<1e-10 and epsilon B<7e-23 for all k>=2000. Hence the denominator
in (7) is at least 0.99, and

```
|Tz-1| <=epsilon B+lambda N4+rho^2/2<rho,
|Tu-1| <=[alpha rho/(1-rho)
             +epsilon/(2(1-rho)^2)+lambda N4]/0.99<rho.
```

The two absolute derivative row sums are at most

```
rho+4epsilon B+lambda N4,
(2alpha+3epsilon+lambda N4)/0.99^2,
```

both below 0.02. Banach contraction proves a unique stationary product
in this box in the declared homogeneous truncation. No statement about
other branches or time-dependent higher-derivative modes follows.

Put D0=epsilon B+lambda N4. At the fixed point,

```
|z-1|, |u-1| <=2D0<0.0002001.                         (9)
```

For z, absorb (rho/2)|z-1| in (7). For u, bound its denominator change
by 2alpha|z-1|+epsilon+lambda N4 <=(1+4alpha)D0 and divide by 0.99.
These are squared-radius bounds about the classical reference, not an
allocation of the plan's opening error budget.

With g0=g(1,1), sigma0=sigma(1,1) and A0=A(C0/q,-log q), the joint
first-order shifts are

```
delta z1=-epsilon A0-lambda g0,
delta u1=epsilon(a20/2-alpha A0)-lambda(sigma0+alpha g0),
|z-1-delta z1|, |u-1-delta u1| <=100 D0^2.             (10)
```

Here a20=a2(C0/q). Set R=2D0. The z remainder is bounded by
(4epsilon B+lambda N4)R+R^2/2<=10D0^2. Write the denominator of Tu
as 1+h, expand alpha/z at z=1, and use
|a2(r)/z^2-a20|<=8R and |sigma-sigma0|<=N4 R. The u remainder is
at most

```
alpha(10D0^2)+alpha R^2/0.99+4epsilon R
                    +lambda N4 R+(1+4alpha)^2 D0^2/0.99,
```

which is below 100D0^2 since B>=100. This bound is relative to the
exact Dirac pressure integrals. It does not upgrade their earlier
adaptive quadrature errors to interval certificates.

The coefficient range is explicit even without computing the pressure:
lambda(k)=48pi C0/[k^10(C0-1)] decreases with k. Thus one fixed set
of coefficients satisfying N4<=1/(10000 lambda(2000)) works for the
whole family. If |c_i|<=cmax for i=1,...,8, then N4<=324 cmax.
For example, cmax=2e24 is a sufficient common bound for (8). This very
large dimensionless bound reflects the small reference curvature in
Planck units; it is neither an inferred physical coefficient bound nor
a guarantee of small nonspherical response. c9 and the J^2 direction
are unrestricted by this homogeneous test and require their own matching.

## Verification and consequence

[The checker](../checks/finite_local_reference_window.py) independently
contracts the product curvature and magnetic tensor, varies the
fixed-flux action, checks the two Einstein equations and first-order
identities, and verifies every coefficient and contraction majorant
with rational arithmetic. High-precision local-only roots for separate
coefficient directions and mixed signs test the nonlinear shifts and
quadratic remainder. They supplement the analytic theorem; the full
Dirac input is reused from its accepted, hashed receipt and proof.

This removes the need to set the displayed finite bosonic inputs to zero
to obtain the conditional homogeneous reference. Their physical values,
the independent matter/boundary action, regulator matching, nonlinear
constraints and driven state still have to be supplied or controlled.
The invariant product stress still has zero AdS2 null contraction.
Closed predecessor workstreams remain closed.

Remaining to-do list: complete finite physical and boundary matching,
independent Dirac terms, nonlinear constraints and remaining interactions;
then S075-S080 state, stress/geometry, stability, transport,
preparation/creation and final audit.
