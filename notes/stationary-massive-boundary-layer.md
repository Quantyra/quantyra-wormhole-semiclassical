# The stationary massive tower's outer boundary layer

2026-10-15, E003/S139. Ordinary fixed-reference argument reviewed before
focused verification. This note determines an outer asymptotic term of the actual
balanced Dirac domain. It does not bound the global massive stress or solve
the full semiclassical exterior. No numerical verification is reported here.

On the prescribed nonextremal RN reference, the massive angular tower has
an orbit-trace boundary term proportional to the inverse conformal distance
to infinity. Its coefficient cancels the logarithmic mass integral and the
leading lapse term of the retained bulk Weyl action. After subtraction of
this term the local remainder is at most logarithmic. The constants in that
remainder are finite for each fixed reference, but are not evaluated or
proved uniform in its very small temperature. This is an outer-endpoint
result, not a quantitative solution of S139.2.

## Reference, domain and subtraction

Use the same single four-component massless charged Dirac field, identical
paired RN exteriors, magnetic bundle, laboratory clock and zero-coupling
stationary domain as in [the conformal-domain note](stationary-exterior-conformal-domain.md).
All formulas below are per exterior. Write

```
m=GM, Q=r_m^2, lambda=L4^-2,
F=1+lambda r^2-2m/r+Q/r^2,
xi=b^2/L4^2, a_m=m/b, h_m=Q/b^2,
u=b/r, B=u^2 F=xi+u^2-2a_m u^3+h_m u^4,
ghat_2=B dU_E^2+b^2 du^2/B.
```

The outer boundary is u=0. It is totally geodesic, with B'(0)=0. The HH
Euclidean period is the prescribed finite period of the smooth RN bolt;
fermions have the corresponding spin structure. The higher angular masses
and signed multiplicities are

```
m_n^2=n(n+q)/b^2, d_n=2(q+2n), n>=1.                (1)
```

The factor two in d_n includes the two nonzero signs of the sphere Dirac
eigenvalue. The remaining two-dimensional spinor has rank two. The folded
mouth space duplicates both: division by two gives the per-exterior
normalization in this note. The LLL has zero angular mass and does not
contribute to the mass-weighted sum used below.

Let L2=D2^dagger D2 be the positive two-dimensional square on the actual
complementary-adjoint domain. At the geodesic boundary its half-density
form has constant complementary Dirichlet/Neumann projectors. This is not
a replacement by a self-adjoint Euclidean Dirac or by mass-dependent Robin
data. The exact separated square is L2+m_n^2. The first-order modulus
variation initially averages the two partner squares; the exact identical
pair intertwiner in [the modulus/phase audit](stationary-paired-dirac-phase.md)
makes their scalar local traces equal and justifies using either square
here. That audit supplies the stationary phase allocation separately.

Keep the summed heat kernel and one common four-dimensional subtraction.
For an interior metric Weyl variation at fixed sphere radius, the
nonlocal part of the Lorentzian orbit trace has the sign

```
tau_hat,mass(y) = FP integral_0^infty A(t) k2(t;y,y) dt
                  + smooth local terms,
A(t)=sum_(n>=1) d_n m_n^2 exp(-m_n^2 t),             (2)
```

where k2 is the spin trace per unit two-volume. FP denotes the retained
bulk heat subtraction at that interior point. Possible local boundary
terms are supported at the boundary, not at y>0. Their charges are not
being determined by (2). The sign in (2) is important: with the accepted
Lorentzian convention the local fermion heat anomaly has the opposite
sign to the explicit +m_n^2 Green trace. For example tr a1=-R2/6 gives
the massless anomaly +R2/(24pi); the large-mass curvature term of the
Green trace cancels it. Using the Euclidean effective-action variation
without the Lorentzian stress convention would reverse this sign.

## Exact balanced image cancellation

Let y be inward proper distance from the outer boundary, and rescale the
time coordinate there by sqrt(xi). In a normal spin frame write

```
ghat_2=dy^2+a(y)^2 dt_E^2,
a(y)=sqrt(B(u(y))/xi),
y=b integral_0^u B(v)^(-1/2) dv.
```

After conjugation by the radial half-density sqrt(a), the differential
expressions are

```
Dtilde=gamma_y partial_y+gamma_t a^-1 partial_t,
Ltilde=-partial_y^2-a^-2 partial_t^2
       -gamma_y gamma_t (a^-1)' partial_t.          (3)
```

Extend a evenly. This extension is generally only C2; no smooth-double
theorem is applied to it. Reflection P:y->-y combined with
chi=tau3 gamma_y commutes with (3). Its even sector has the actual
Dirichlet/Neumann conditions. The method of images is therefore exact
for its heat kernel, including a nonsmooth even extension defined by its
closed quadratic form. Its reflected contribution to the spin trace is
zero. Indeed Ltilde and its heat kernel commute with the two-dimensional
chirality i gamma_y gamma_t, whereas chi anticommutes with chirality.
Conjugating chi K by chirality reverses its sign and leaves its trace
unchanged. Thus

```
tr K_half(t;y,y)=tr K_even-double(t;y,y).           (4)
```

There is no extra factor of one half in this diagonal identity. A smooth
even metric, including B0=xi+u^2, has no reflected spin-trace layer.
Equation (4) does not erase the RN boundary contribution: the RN even
extension contains a |y|^3 term and is not smooth enough to invoke the
closed smooth-manifold expansion across y=0.

## The first non-even jet and its coefficient

Normal-coordinate expansion gives

```
a(y)=1+y^2/(2b^2)+c3 y^3+O(y^4),
c3=-a_m sqrt(xi)/b^3,                              (5)
```

on the physical half. Its even continuation has c3 |y|^3. The even
quadratic jet can first affect the cubic calculation at higher parabolic
weight, so the leading difference from a smooth continuation is obtained
by one flat Duhamel insertion of this cubic jet. In (3) that insertion is

```
delta L=2 c3 |y|^3 partial_t^2
        +c3 gamma_y gamma_t (|y|^3)' partial_t.
```

The connection term is Clifford-even but traceless at this order. For
the scalar term the tangential Gaussian momentum has mean k^2=1/(2t).
The density conversion at the observation point cancels between the even
and smooth continuations, since both agree for y>=0. Hence the boundary
profile is

```
K3(t,y)=c3/(2pi t^2) integral_0^t
           D(y,sqrt(2s(t-s)/t)) ds,
D(y,sigma)=E[|y+X|^3-(y+X)^3], X~Normal(0,sigma^2). (6)
```

This is a spin trace of rank two. D is nonnegative and rapidly decreasing
when y/sigma grows. Formula (6) is the whole parabolic-weight-three
profile, not just its boundary value.

The angular heat expansion on the fixed smooth magnetic sphere gives

```
A(t)=2b^2/t^2+O(1) as t->0.                        (7)
```

One can also obtain (7) by differentiating the sphere heat trace
2b^2/t+constant+O(t). The constant term disappears under differentiation;
there is no t^-1 term in A. Magnetic and finite-q corrections do not
change the leading coefficient.

The coefficient in (2) from (6)-(7) can be integrated without asymptotic
interchanges of signed functions: D>=0 permits Tonelli, followed by the
finite scalar multiplier c3. For v>0,

```
integral_0^infty t^-3 D(y,sqrt(vt)) dt
 =6v^2 integral_y^infty (x-y)^3/x^5 dx
 =3v^2/(2y),
integral_0^1 [2z(1-z)]^2 dz=2/15.
```

Using s=tz in (6) therefore yields exactly

```
integral_0^infty (2b^2/t^2) K3(t,y) dt
 =c3 b^2/(5pi y).                                 (8)
```

Extending this model integral to all t only changes a bounded remainder
when comparing it with a fixed small-time cutoff. It is not an assertion
that the cubic model is the global geometry or the complete heat kernel.

## Why the remaining local trace is at most logarithmic

Here the needed result is a qualitative local remainder with finite
fixed-reference constants, not an explicit uniform heat estimate. For
this specific stationary domain it follows directly from a scalar
Brownian bridge formula. Choose bounded positive continuations f_s,f_e
of f=1/a, with bounded first derivative: f_s is smooth across zero and
f_e is even and C2. They agree with each other and the actual geometry
on a fixed positive collar. They can be frozen outside a larger collar;
no behavior at the RN bolt is inferred from this auxiliary extension.
Local heat-kernel uniqueness makes changes beyond that collar rapidly
decreasing at small t. Choose k_int to be the smooth continuation's
kernel, and k_bd to be the actual traced kernel minus k_int.

In (3), Fourier transforming time and diagonalizing gamma_y gamma_t
gives two scalar operators with potentials k^2 f^2 +/- k f'. These
potentials are bounded below uniformly in k, by completing the square.
The Brownian-bridge Feynman--Kac formula therefore applies. Its standard
semibounded version is [Simon, theorem 1.1](https://arxiv.org/html/math-ph/9907022),
with time rescaled from his generator Delta/2 to our generator Delta.
Summing the two scalar components and performing the Gaussian k integral
gives the exact zero-winding expression

```
k_f(t;y,y)=1/[2pi t a(y)] E[ (A/t)^(-1/2) exp(Cb^2/(4A)) ],
A=integral_0^t f(Y_s)^2 ds,
Cb=integral_0^t f'(Y_s) ds,
Y_s=y+sqrt(2t) beta_(s/t),                          (9a)
```

where beta is the standard unit Brownian bridge. This formula is valid
for the C2 even extension as well as for the smooth extension. On the
fixed thermal circle the momentum sum is obtained by Poisson summation.
All nonzero antiperiodic windings are exponentially small as t->0,
uniformly in this collar, since A is between positive constants times t.
The bounded factor exp(Cb^2/(4A)) does not affect that conclusion.

To prove the required remainder put y=sqrt(t)Y0. In the collar,

```
f_s(z)=1-z^2/(2b^2)-c3 z^3+O(z^4),
f_e(z)=1-z^2/(2b^2)-c3 |z|^3+O(z^4),
f_s'(z),f_e'(z)=-z/b^2+O(z^2).
```

Taylor expansion of the smooth scalar function in (9a) through degree
three in sqrt(t) leaves precisely
c3 t^-1 integral[|Y_s|^3-Y_s^3]ds between the two continuations.
The Cb^2/A contribution starts at degree four, since f'(0)=0.
The remaining difference is bounded by C t^2 times a polynomial in
Y0 and the supremum of the unit bridge. It vanishes on every path
staying on the positive half, where f_s=f_e. The crossing probability
for a generator-Delta bridge from y to y in time t is exp(-y^2/t).
Cauchy--Schwarz and the finite Gaussian moments of its supremum bound
the expected residual by C t^2 times a polynomial in Y0 times
exp(-Y0^2/2). Polynomial factors can be absorbed into an arbitrarily
large inverse power of 1+Y0. Paths exiting the fixed collar have the
usual exponentially small bound; bounded positive continuations make
the integrand in (9a) uniformly bounded by a constant times exp(Ct).
The common density a(y)^-1 equals 1+O(y^2), changing the cubic profile
only at higher degree. Equations (6) and (9a) consequently prove that,
for some y0,t0>0 and each fixed N, there is a finite C_N such that

```
|k_bd(t,y)-K3(t,y)|
 <=C_N t (1+y/sqrt(t))^-N,
0<t<=t0, 0<=y<=y0.                                (9)
```

The ordinary boundary heat calculus gives the same profile order.
[Grieser, theorem 3.1, proposition 3.11 and theorem 3.12](https://web.math.ku.dk/grubb/notes/heat.pdf)
give the smooth interior part, scaled boundary profiles and local
Volterra correction. His displayed theorem is for Dirichlet data; it
is not being quoted as an already proved matrix mixed theorem here.
The constant-projector reduction (3)-(4) and scalar formula (9a) supply
the domain-specific proof instead. The mixed local coefficients are
given by [Vassilevich, section 5.3](https://arxiv.org/html/hep-th/0306138v3).
None of these assertions evaluates a temperature-uniform variation
constant. As an independent normalization check, the order-t coefficient
of (9a) after division by a is a''/(6a)=-R2/12, giving the rank-two
interior heat term -R2/(24pi) with the expected sign.

Subtract the smooth interior heat coefficients through a2 before the
four-dimensional regulator is removed. This can equivalently be done by
adding and subtracting their exact mass-weighted angular sum with a
small-time cutoff and restoring its common-scheme finite part; those
restored terms are smooth functions of y. It does not select separately
renormalized angular determinants. Their remainder in k2 is O(t^2),
uniformly in this fixed collar. Multiplying it by (7) is integrable at
t=0. The O(1) correction in (7) times K3 is also integrable and bounded
as y->0. For the residual boundary part choose N>0 in (9). The only
possibly unbounded integral is bounded by

```
C integral_0^t0 t^-1 (1+y/sqrt(t))^-N dt
 <=C' [1+|log(y/y0)|].                             (10)
```

Split at t=y^2 to prove the last inequality. Higher masses make the
large-time part smooth at the boundary: L2>=0 and m1>0, while the fixed
smooth compact Euclidean reference has a bounded heat kernel at any
positive time. Local finite bulk subtractions are smooth to y=0.
Boundary counterterms have no interior support. Thus (8)-(10) prove

```
tau_hat,mass(y)=c3 b^2/(5pi y)+O(1+|log(y/y0)|).
                                                               (11)
```

No value for C_N, C' or the corresponding metric derivatives is claimed.
In particular (11) cannot be used as a uniform norm bound through a
varying bolt or as a numerical clock allowance at k=2000.

## Conserved stress and the outer mass/lapse cancellation

The HH Euclidean state is smooth at the fixed bolt. Write
tau=trace_hat T for the massive orbit trace and t=T_uu=T_vv in static
optical coordinates. With B=e^(2psi)=b^2 F/r^2, exact two-dimensional
conservation and HH regularity give

```
T_uv=-B tau/4,
dt/dr=-B tau'/4, t(r_+)=0,
t(r)=-B(r)tau(r)/4+(1/4) integral_(r_+)^r B'(R)tau(R)dR.
                                                               (12)
```

The last expression is essential: it reconstructs the stress from (11)
without differentiating a big-O remainder. The corresponding physical
orbit sources at fixed sphere are

```
p2(r)=I(r)/(2F(r)),
rho2(r)=I(r)/(2F(r))-b^2 tau(r)/r^2,
I(r)=integral_(r_+)^r B'(R)tau(R)dR,
J_X=-2G trace_g T/r^2.                             (13)
```

These are the same stress/radius source dictionary as
[the isolated LLL theorem](stationary-lll-quantum-exterior.md); setting
the radius source to zero is not allowed. Since

```
y=b^2/(sqrt(xi) r)+O(r^-3),
tau(r)=-m lambda r/(5pi b^2)+O(1+log r),
B'(r)=O(r^-3),
```

the integral I has a finite limit I_infinity. Equations (12)-(13) imply

```
rho2_hat,mass=+m lambda/(5pi r)+O((1+log r)/r^2),
p2_hat,mass=O(r^-2),
delta alpha_hat,mass=-Gm/(15pi r^3)
                     +O(G(1+log r)/r^4).          (14)
```

Here delta alpha is the first-order constraint contribution on the fixed
RN background, normalized to zero at infinity; alpha'=G(rho2+p2)/(rF).
It is not the nonlinear lapse of an independently corrected geometry.

For comparison, compactly supported metric variation of the retained
bulk [Weyl polynomial](spherical-anomaly-action.md) gives the following
RN result. Set sigma=log(r/b), T=r(F'+F alpha') and

```
Rbar=-r^2 F''-3r^2 alpha' F'-2r^2 F(alpha''+alpha'^2)
     +2T-2F.
```

After integrating sigma Box Rhat by parts, its static density is

```
4pi W4=integral dr e^alpha/r^2 {
 sigma(Rbar^2/60-Rbar/18+1/15-q^2/3)
 +Rbar T/30-T^2/10+11TF/90-11F^2/180-F/18 }.
```

The compact-support Euler derivative with respect to F at alpha=0 is
2(6Q sigma-8Q+3mr)/(15r^4). The stress and lapse consequences are

```
rho2_W=11lambda^2 r^2/(240pi)-m lambda/(5pi r)
       +O((1+log r)/r^2),
delta alpha_W=G/(15pi) [m/r^3
                +Q((3/2)sigma-13/8)/r^4].        (15)
```

The first term in rho2_W is the asymptotic AdS vacuum term. The massive
trace (11) and the LLL trace have no r^2 physical-orbit contribution, so
the full fixed-reference source has this leading AdS-invariant vacuum
automatically. Whether it is retained as vacuum stress or absorbed into
the measured cosmological constant is part of the common outer vacuum
prescription, not an independent subtraction permitted for just one
sector. The full compact-variation sources and radius Ward check are in
[the radial Weyl-source note](stationary-weyl-radial-sources.md).

Combining (14) and (15), the coefficient of the logarithmically divergent
mass integral cancels, as does the leading r^-3 lapse:

```
rho2_hat,mass+rho2_W-rho2_AdS=O((1+log r)/r^2),
delta alpha_hat,mass+delta alpha_W
 =O(G(1+log r)/r^4).                               (16)
```

More explicitly p2_W=-rho2_AdS+O((1+log r)/r^2), while the massive
pressure is O(r^-2) and the LLL orbit sources have integrable O(r^-2)
tails. Dividing by 4pi r^2, and using the full trace/radius Ward identity
for the angular component, proves the fixed-reference asymptotic vacuum
limit in orthonormal components:

```
T4_(mu)(nu) -> -[11lambda^2/(960pi^2)] g_(mu)(nu).  (17)
```

This is a result for the specified state on the classical reference.
If the same asymptotic vacuum prescription is imposed on a corrected
geometry, its leading curvature equation would instead be
lambda_eff+11G lambda_eff^2/(360pi)=L4^-2. Neither the existence of that
global corrected geometry nor a uniform transfer of (11) to it follows
from this local fixed-reference theorem.

Consequently the outer endpoint of the corresponding fixed-reference
first-order geometric mass integral converges after the common AdS
vacuum has been handled consistently. Its finite value and the boundary
charge have not been calculated. The LLL contribution has an integrable
O(r^-2) mass tail and does not alter this cancellation. A raw product-frame
trace bounded in C1 up to u=0 would miss the term responsible for (16).

## Remaining scientific work

The local asymptotic theorem supplies neither an evaluated global stress
bound nor its first and second metric/radius responses. The near-bolt
estimate must retain the HH cancellation uniformly in the prepared small
temperature; a generic covariant derivative estimate can lose inverse
surface gravity. The boundary estimate must also be stable under the
allowed metric and boundary-domain variations. Both are required before
a full semiclassical fixed point or a quantitative lapse claim can be
made. Boundary Weyl transgression, the complete finite charge, the outer
vacuum prescription and the modulus/phase allocation retain their own
obligations. S139.1-S139.2 and the physical state transfer remain open.

Remaining to-do list: obtain a quantified differentiated global massive
remainder with the bolt and boundary weights; compute the complete
stationary exterior and common state/charge matching; transfer the
driven reference estimates to that physical state.
