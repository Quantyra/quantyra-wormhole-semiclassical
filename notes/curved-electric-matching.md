# Common electric subtraction and the massive AdS coefficient

S074, 2026-10-13. The retained spectral Maxwell subtraction agrees
exactly with the homogeneous determinant's zeta prescription in their
common flat product limit. On AdS2 the complete massive tower adds a
finite, negative curvature correction, with an explicit all-level bound.
The constant-field electric coefficient can therefore be expressed in
the reference scheme, including the displayed finite bosonic terms.
This resolves one common bulk matching calculation; it does not supply
the boundary charge susceptibility, full driven kernel or physical state.

## Model and observable

Use one unit-charge massless four-dimensional Dirac field on
AdS2(a) x S2(b), fixed integer magnetic flux q>=1 and the complete
Euclidean H2 spectral prescription of the
[uniform reference](uniform-dirac-reference.md). The massive sphere
channels have m_n^2=n(n+q)/b^2 and d_n=2(q+2n), n>=1. Remove the q
massless LLL channels from the following coefficient: their nonlocal
current response is already retained in the light sector.

The observable is the coefficient c of f^2/2 per longitudinal physical
volume in the weak, covariantly constant Euclidean field f, with both
radii held fixed. Analytic continuation gives the quadratic Lorentzian
electric density with the same coefficient. This is a local invariant
background derivative, not a finite-energy global drive or a boundary
capacitance. Flux, subtraction scale mu and renormalized couplings stay
fixed during the derivative. No external Dirac operators or additional
loops are introduced.

The source spectrum is standard:
[Castellano, Montella and Zatti, (2.63)-(2.68), (2.74), and (B.43)](https://arxiv.org/html/2603.15730v1)
give the hyperbolic magnetic levels, continuation and planar heat trace.
We use the massive minimal Dirac operator, without their supersymmetric
specialization. The response and bounds below are derived from these
inputs; they are not a claim to a first determinant on this geometry.

## Finite Maxwell matching in the shared flat limit

On R2 x S2, a weak longitudinal field has plane heat trace
f coth(fs)/(2pi)=1/(2pi s)+f^2 s/(6pi)+O(f^4).
Define the sphere spectral function with the LLL omitted,

```
Z_H(w)=sum_(n>=1) d_n (m_n^2)^(-w).
```

Its continuation at w=1 is fixed by an absolutely convergent difference:

```
Z_H(1+t)=b^(2+2t) [4 zeta_R(1+2t)
 +sum_(n>=1){2(q+2n)/[n(n+q)]^(1+t)-4/n^(1+2t)}],
Z_H(1+t)=2b^2/t+2b^2[2log b+2gamma_E-H_q]+O(t).       (1)
```

At t=0 the series in braces sums to -2H_q. Uniform convergence near
t=0 follows from its O(q n^-2-2Re(t)) large-n tail at fixed q.
The f^2 term in the full zeta function, after integrating the sphere,
is f^2 t Z_H(1+t)/(6pi), since Gamma(t+1)/Gamma(t)=t.
Use the same action convention (zeta'(0)+log(mu^2)zeta(0))/2 as the
retained determinant. Equation (1) gives

```
c_H,flat^zeta(0;mu)=b^2/(3pi)[2log(b mu)+2gamma_E-H_q]. (2)
```

This is exactly the [spectral bubble's coefficient](massive-bulk-polarization.md):
its finite-cutoff sum is b^2[H_N+H_(N+q)-H_q]/(3pi), and it subtracts
2b^2 log(sqrt(N(N+q))/(b mu))/(3pi). No finite F^2 conversion remains
between these two explicitly defined prescriptions at the same mu.
This equality identifies their renormalized Maxwell parameter, not the
flat effective electric coefficient with its curved analogue.

The common scale cancellation is
d(1/g_R^2)/dlog mu=-1/(6pi^2). Independently, multiplying the two plane
and sphere heat expansions gives the f^2 heat coefficient 1/(12pi^2)
per four-volume and the same running. Finite physical EFT coefficients
are still independent inputs; matching a subtraction does not fix them.

## Exact massive constant-field coefficient on AdS2

Let M=ma>0, beta=f a^2 and nu=sqrt(M^2+beta^2). For 0<beta<1 the
dimensionless hyperbolic heat trace contains the zero-level density
beta/(2pi) and continuum density

```
rho_beta(l)=l sinh(2pi l)/[pi(cosh(2pi l)-cos(2pi beta))],
eigenvalue=l^2+beta^2.
```

The discrete term is essential, even for arbitrarily small beta.
Subtract the zero-field resolvent before integrating. Its exact
difference is

```
Delta R_beta(M)=beta/(2pi M^2)
 +integral_0^infinity [rho_beta(l)/(l^2+M^2+beta^2)
                          -rho_0(l)/(l^2+M^2)] dl
 =-[psi(nu+beta)+psi(nu-beta)-2psi(M)
                              +(nu-M)/M^2]/(2pi).    (3)
```

To derive it, expand the density using
sinh(2pi l)/(cosh(2pi l)-cos(2pi beta))
=(l/pi) sum_(j in Z)[l^2+(j+beta)^2]^-1. Integrating each rational
resolvent difference uses
integral_0^infinity l^2/[(l^2+A^2)(l^2+nu^2)] dl
=pi/[2(A+nu)]. The resulting paired harmonic differences converge;
their digamma sum plus the discrete term gives (3). This also shows
why differentiating the continuum density alone at beta=0 is invalid.

The dimensionless action w obeys partial_(M^2) w=-R/2. Expanding the
finite expression (3), then imposing decoupling as M tends to infinity,
gives the one-channel coefficient

```
c_m,AdS=a^2/pi [M psi_1(M)-1-1/(2M)].                 (4)
```

Indeed the beta^2 resolvent coefficient is
-[psi_1(M)/M+psi_2(M)+1/(2M^3)]/(2pi). The beta^2 action coefficient
is [M psi_1(M)-1-1/(2M)]/(2pi); its mass derivative is minus half that
resolvent coefficient and it vanishes at infinite mass. This fixes the
finite coefficient without a separate subtraction for each sphere mode.

The positive-mass expansion has a neighborhood |e|<M after beta=-ie;
only its zero-field quadratic coefficient is used. Formula (4) diverges
as M approaches zero. It must not be applied to the LLL or be promoted
to a local replacement for its massless response.

## A signed curvature correction with a complete tower bound

There is a positive Laplace representation

```
c_m,AdS=(a^2/pi) integral_0^infinity exp(-M t) F(t) dt,
F(t)=d/dt[(t/2)coth(t/2)]
    =4t sum_(j>=1) A_j^2/(t^2+A_j^2)^2, A_j=2pi j.   (5)
```

The partial fractions follow from [DLMF 4.36.3](https://dlmf.nist.gov/4.36.E3).
For x>=0, 1-2x<=(1+x)^-2<=1 and
0<=(1+x)^-2-(1-2x)<=3x^2. Positivity and the known even zeta sums
therefore give the global inequalities

```
0<F(t)<=t/6,
0<=t/6-F(t)<=t^3/180,
0<=F(t)-t/6+t^3/180<=t^5/5040.
```

Laplace integration proves, for every M>0,

```
0<c_m,AdS<=c_m,flat=1/(6pi m^2),
0<=c_m,flat-c_m,AdS<=1/(30pi a^2 m^4),
0<=c_m,AdS-c_m,flat+1/(30pi a^2 m^4)
                                      <=1/(42pi a^4 m^6). (6)
```

Hence Delta_curv=sum_(n>=1)d_n(c_mn,AdS-c_mn,flat) converges absolutely.
The retained telescoping identity sum d_n/m_n^4=2b^4 H_q^(2)/q gives

```
c_H,AdS^ren=c_H,flat^zeta+Delta_curv,
-B_curv<=Delta_curv<=0,
B_curv=b^4 H_q^(2)/(15pi a^2 q),
0<=Delta_curv+B_curv<=b^6 H_q^(2)/[21pi a^4 q(q+1)].   (7)
```

The last estimate uses m_n^2>=m_1^2. It bounds the remainder of the first
curvature term by 5C/[7(q+1)] times that term, where C=b^2/a^2.
The uncomputed tail after sphere level N obeys

```
|Delta_curv-Delta_curv,N| <= b^4/(15pi a^2 q)
                     [psi_1(N+1)-psi_1(N+q+1)].       (8)
```

Here Delta_curv,N denotes the partial sum through N. Equation
(8) is an analytic angular tail, not a quadrature error estimate.
No finite curvature term is absorbed silently into the Maxwell coupling.
The difference is UV finite, so the scale running and the finite matching
in (2) also apply to the curved coefficient in this fixed prescription.
At large n the bound on each weighted difference is O(n^-3). Thus the
weak-field coefficient of the four-dimensional subtracted determinant
can be compared with its flat counterpart before taking the sphere sum:
the difference has no remaining local subtraction or finite summation
ambiguity at this order. This does not assert such a result for the
unsubtracted determinant or for every other composite observable.

## Combining the loop with the displayed finite local terms

Use canonical Fbar=F/g_R, kappa^2=8pi G and the unredefined coefficients
of the [finite local reference window](finite-local-reference-window.md).
At fixed metric let Ebar be the electric field parallel to the magnetic
field Bbar=q/(2g_R b^2). Direct contraction gives

```
I=2(Bbar^2-Ebar^2), I4=2(Bbar^4+Ebar^4),
Ric.X=2x Ebar^2+2y Bbar^2,
Riem.FF=4x Ebar^2+4y Bbar^2, x=a^-2, y=b^-2.
theta_loc=4kappa^2(Kx-2c4 y)-16kappa^4 c7 Bbar^2,
K=2c4+c5+2c6.                                       (9)
```

The Lorentzian local electric density is (1+theta_loc)Ebar^2/2.
Curvature squares have no direct fixed-metric electric derivative;
D.D vanishes for this covariantly constant field. Their other responses
are not discarded. In particular, c8 affects the background but has no
direct Ebar^2 term here. The isolated J^2 direction reproduces exactly
theta_loc=2kappa^4 cJ Bbar^2 from the preceding electromagnetic audit.

In the unit-charge field normalization, the heavy/local coefficient is

```
kappa_el=4pi b^2/g_R^2 (1+theta_loc)+c_H,AdS^ren.       (10)
```

This does not include the nonlocal LLL response. Additional finite
operators would have to be added in the same scheme; (10) is complete
only for the specified one-Dirac-loop plus displayed bosonic span.

For the candidate family and the one-percent reference box, write
lambda=kappa^2/a0^2, alpha=1/C0, u=a^2/a0^2, z=b^2/b0^2. Then

```
theta_loc=4lambda(K/u-2alpha c4/z)
                         -16lambda(1+alpha)c7/z^2,
|theta_loc|<=lambda NE,
NE=5|K|+3|c4|+21|c7|.                               (11)
```

The bound uses alpha<=1/4 and u,z in [0.99,1.01]. Since mu0=1/b0,
the flat loop fraction is
g_R^2[log z+2gamma_E-H_q]/(12pi^2). Its magnitude plus the curvature
bound is below 2e-27 for every k>=2000. For example, H_q<=1+log q,
|log z|<1, gamma_E<1, pi>3 and C/q<8/k reduce this to decreasing
multiples of k^-8(1+log k) and k^-9. At k=2000, the curvature correction
is within 0.22 percent of its leading term; the checker records its
explicit analytic interval and the flat-loop fraction.

Thus the previous reference condition lambda N4<=1e-4, supplemented by
lambda NE<=1/200, guarantees kappa_el>0.99(4pi b^2/g_R^2) throughout
the candidate box. The previous simple cube |c_i|<=2e24, i=1,...,8,
satisfies both conditions. This is a conditional bulk coefficient range,
not a measurement of c_i, a boundary inertia or a stability theorem.
Unlike N4, NE detects the previously invisible J^2 direction.

## Verification and consequence

[The checker](../checks/curved_electric_matching.py) verifies the finite
zeta term, spectral-subtraction limit, resolvent/action differentiation,
positive partial-fraction remainder bounds, tower identity, running and
local invariant contractions. Direct continuum-plus-discrete resolvent
quadrature is compared with (3), and independent Laplace quadrature
with (4). Finite tower sums retain (8). Numerical quadratures supplement
the analytic proof and are not interval certificates.

The Maxwell parameter used by the homogeneous determinant and the flat
bubble is now matched explicitly. The massive constant-field AdS
correction is controlled in that same scheme. Neither result transfers
the flat momentum-dependent bubble or its light pole to a curved,
boundary-driven problem. The global absorptive work estimate, its state
assumptions and the physical charge/clock ensemble retain their separate
matching obligations. Full S074 and the creation objective remain open.

Remaining to-do list: complete independent matter/boundary terms,
nonlinear constraints, physical state and general response matching;
then S075-S080 stress/geometry, stability, transport,
preparation/creation and final audit.
