# Shifted RN boundary moments for the complete source pulse

2026-10-15, E003/S138. Ordinary proof, independently reviewed before
focused verification.
This note supplies the static trace estimates used by the complete-pulse
argument. It treats the actual gapless classical magnetic RN-AdS exterior,
the paired Dirac domain and every nonzero angular mode. It proves neither
transfer to the quantum-corrected exterior nor a spatial stress estimate.

## Reference, domain and normalization

Retain exactly the geometry, angular multiplicities and thermal reference
of [the static RN comparison](rn-higher-mode-boundary-spectral-control.md).
In particular, eta=1/2000 is fixed, k=2000 for the proposed work target,
and the uniform family estimates below also hold for k>=2000. Write

```
q=k^4, X(1+3X)=pi k^6, C=1+6X,
lambda_n=sqrt[n(n+q)], M_n=lambda_n/sqrt(C),
d_n=2(q+2n), n>=1, d=b/sqrt(CX), delta=nu d,
M_n>=18, delta<=2 eta sqrt(7)/k^9<1e-32,
kappa0<=7nu/25.
```

Both nonzero angular signs are included in d_n. Zero modes are excluded.
Use the canonical radial half-density and physical Killing time U;
h_n=H_n/nu is the Hamiltonian in laboratory time nu U. The dimensionless
boundary trace is Gamma_lab=Gamma_U/sqrt(nu). All boundary traces below
include both exterior allowed components, of total rank two per angular
state. There is no additional factor for a third or fourth exterior.

In tortoise coordinate x<=0 the one-exterior expression is
H=-i alpha partial_x+beta W, alpha=sigma_2, beta=sigma_1,
W=lambda_n sqrt(f)/r. Folding gives radial velocity A=tau_3 tensor alpha
and allowed projector Pi0=(1+tau_3 tensor beta)/2. Radial velocity
anticommutes with the projector involution; the angular potential
commutes with it. Requiring both psi and H psi to have allowed traces
therefore gives, for the square,

```
(1-Pi0)psi(0)=0, Pi0 psi'(0)=0.
```

These are the physical mixed Dirichlet/Neumann conditions. There is no
Robin W term. The potential of the square is W^2+Y W', where Y is a
Hermitian involution; it need not commute with Pi0. The horizon is the
limit-point end, with its usual self-adjoint prescription. A driven
Euclidean squared domain has other terms and is not being substituted
for this static Hamiltonian square.

The fixed domain-frame generator is K=-iG/2, where
G=i tau_2 tensor beta. Thus K is Hermitian and has norm 1/2. Its bulk
commutator with the identical-exterior differential expression vanishes;
it does not preserve the reflecting domain. The Green identity instead
defines the weak boundary commutator. With the canonical trace, its
boundary matrix is a phase times C_boundary/2, where
C_boundary=i A G is Hermitian of norm one. Consequently a separated
rational divided-difference term (1/2)r(h) Gamma_lab^dagger C_boundary
Gamma_lab s(h), up to a unit phase, has trace norm at most

```
(1/2) ||Gamma_lab r(h)^dagger||_HS
      *||Gamma_lab s(h)||_HS.                         (1)
```

Equation (1), not a bounded-operator assertion about [h,K], fixes the
normalization in the calculations below. It follows first by the Green
identity on smooth resolvent images, then by completion: all the trace
maps on those images have finite Hilbert--Schmidt norms by the moments
proved below. The resolvent commutator identity in this weak form gives
the rational divided differences. No product of a boundary delta
function with an arbitrary Hilbert-space vector is used.

## Arbitrary positive shift and the scalar comparison

The retained RN calculation proves, uniformly in the preparation range,

```
W^2-|W'| >= W^2/2-kappa0^2,
W(-y)>=M_n/(3d), 0<=y<=d.
```

For any dimensionless shift s>=1 put

```
a_s=sqrt(nu^2 s^2-kappa0^2)>=(24/25)nu s,
mu_n=M_n/(3sqrt(2)d),
B_s=-partial_y^2+a_s^2+mu_n^2 1_[0,d], y>=0,
```

with scalar Neumann boundary at y=0. Denote b0=24/25,
Lambda_n=mu_n/nu=M_n/(3sqrt(2)delta),
Z_n=exp(-M_n/sqrt(6)) and
k_s=sqrt(Lambda_n^2+b0^2 s^2). In particular Lambda_n>1.

For clarity, the heat comparison used here is a kernel comparison, not
operator monotonicity of an inverse power. The free mixed kernel is
K_N Pi0+K_D(1-Pi0), with 0<=K_D<=K_N. Its matrix norm is at most K_N.
Every matrix potential exponential has norm bounded by the exponential
of its scalar lower bound. Trotter products and kernel composition
therefore dominate the matrix heat kernel by the scalar Neumann heat
kernel of B_s. For fixed n the potentials are bounded; the usual bounded
potential expansion proves the kernel limit. Scalar positive-potential
approximation permits the step. The boundary diagonal is positive and
supported on rank-two Pi0. Integrating against
t^(r-1)/Gamma(r), for every r>0 for which it is finite, gives

```
tr_boundary Gamma_U(H_n^2+nu^2 s^2)^(-r)Gamma_U^dagger
 <= 2 [B_s^(-r)](0,0).                              (2)
```

Here tr_boundary indicates the finite fiber trace of the boundary
kernel. Equation (2) is only for these positive moments; it does not
assert domination of the actual RN spectral measure.

## Scalar moment bound, including the horizon tail

Write mu=mu_n and a=a_s temporarily. The scalar boundary resolvent is

```
g(a)=[cosh(kd)+(a/k)sinh(kd)]/[k sinh(kd)+a cosh(kd)],
k=sqrt(mu^2+a^2).
```

Hence g(a)<=coth(kd)/k<=2/k, because mu d>=3sqrt(2)>2.
This bound includes the entire scattering continuum, including its
above-barrier resonances.

The scalar continuum is parameterized by energy a^2+p^2, p>=0. There
are no bound states for the nonnegative step. Below the barrier its
Neumann boundary spectral density in dp is

```
w(p)=(2/pi)/[cosh^2(g)+(kappa/p)^2 sinh^2(g)],
kappa=sqrt(mu^2-p^2), g=kappa d, 0<p<mu.             (3)
```

To see the normalization, take the Neumann solution cosh(kappa y) on
the step. At infinity its cosine and sine amplitudes are cosh(g) and
(kappa/p)sinh(g). Dividing their squared sum into the free Neumann
boundary density 2/pi gives (3); there is no extra 2p Jacobian when
integrating in p. Completeness is the ordinary half-line scattering
resolution for this finite step.

On 0<p<=mu/2 the denominator in (3) is exactly
1+(mu^2/p^2)sinh^2(g). Also g>2 and
sinh(g)^(-2)<= (16/3)exp(-2g): use exp(-4)<1/8 and
(7/8)^2>3/4 in 4 exp(-2g)/(1-exp(-2g))^2. Therefore

```
w(p) <= [32/(3pi mu^2)] p^2 exp(-sqrt(3)mu d).       (4)
```

For the complementary spectral interval p>=mu/2 and r>=1, positivity
of the scalar spectral measure bounds its r-th moment by
(a^2+mu^2/4)^(1-r) g(a). For the low interval use (4) and extend its
positive integral to infinity. Since

```
integral_0^infinity p^2/(p^2+a^2)^r dp
 =sqrt(pi) Gamma(r-3/2) a^(3-2r)/(4 Gamma(r)), r>3/2,
```

the full scalar moment satisfies

```
[B_s^(-r)](0,0)
 <= 2*4^(r-1)(mu^2+a^2)^(1/2-r)
    +c_r exp(-sqrt(3)mu d) mu^(-2) a^(3-2r),
c_r=8 Gamma(r-3/2)/(3sqrt(pi)Gamma(r)), r>3/2.       (5)
```

For r=1 simply use g(a)<=2/sqrt(mu^2+a^2).
The exponential term in (5) is retained. This thermal horizon is
gapless: replacing (5) by a hard spectral gap or importing an extremal
zero-energy inverse moment is invalid.

## Dimensionless moments used by the source

Define m_r(s)=tr Gamma_lab(h_n^2+s^2)^(-r)Gamma_lab^dagger.
The normalization is m_r=nu^(2r-1) times the corresponding physical
moment. Equations (2)-(5) give

```
m_1(s)<=4/k_s,
m_r(s)<=4^r k_s^(1-2r)
       +2c_r Z_n Lambda_n^(-2)(b0 s)^(3-2r), r>3/2.
```

In particular, writing Lambda=Lambda_n and Z=Z_n,

```
m_2 <=16/k_s^3 +(16/3)Z Lambda^-2(b0 s)^-1,
m_3 <=64/k_s^5 +(4/3) Z Lambda^-2(b0 s)^-3,
m_4 <=256/k_s^7+(2/3) Z Lambda^-2(b0 s)^-5.          (6)
```

These are arbitrary-shift bounds, including shifts above the angular
barrier scale. The powers of k_s must be kept when integrating over s.

## Exact rational divided differences

Let g_s(x)=x/(x^2+s^2)^2 and b_s(x)=1/(x^2+s^2)^2.
Their divided differences have the common denominator
(x^2+s^2)^2(y^2+s^2)^2 and exact numerators

```
g_s: s^4-xy(x^2+xy+y^2+2s^2),
b_s: -(x+y)(x^2+y^2+2s^2).                          (7)
```

The formula remains valid on x=y by continuity. Put
ell_j=tr Gamma_lab h^(2j)(h^2+s^2)^(-4)Gamma_lab^dagger.
Spectral positivity gives ell_j<=m_(4-j), j=0,1,2,3.
Factor each monomial in (7) and apply (1). This proves

```
||[g_s(h),K]||_1
 <=(1/2)[s^4 m_4+2s^2 m_3+2sqrt(m_1 m_3)+m_2],
||[b_s(h),K]||_1
 <=sqrt(m_1 m_4)+sqrt(m_2 m_3)+2s^2 sqrt(m_3 m_4).   (8)
```

The four moments in (8) refer to (6), not to the ell_j. In particular
the second line is an independently bounded commutator; no universal
comparison to the first line is asserted.

For an explicit elementary envelope, (6), b0=24/25, k_s>=Lambda,
k_s>=b0 s and sqrt(u+v)<=sqrt(u)+sqrt(v) give

```
||[g_s(h),K]||_1
 <=250 k_s^-3+5 Z Lambda^-2 s^-1
                    +3sqrt(Z)Lambda^-3/2 s^-3/2,
||[b_s(h),K]||_1
 <=350 k_s^-4+6 Z Lambda^-2 s^-2
                    +100sqrt(Z)Lambda^-3/2 s^-3/2.  (9)
```

Some detail makes the mixed terms in the second line inspectable.
The pure local coefficient is
32+32+256 b0^-2<350. The additional mixed and pure-tail terms in its
three summands are bounded respectively by

```
sqrt(m_1 m_4):
  2sqrt(Z)Lambda^-3/2 s^-5/2,
sqrt(m_2 m_3):
 25sqrt(Z)Lambda^-5/2 s^-3/2+3Z Lambda^-2 s^-2,
2s^2 sqrt(m_3 m_4):
 59sqrt(Z)Lambda^-5/2 s^-3/2+3Z Lambda^-2 s^-2.
```

Use Lambda>=1 and s>=1 to obtain (9). For example the potentially
troublesome mixed terms use
k_s^-5/2<=Lambda^-3/2(b0 s)^-1 and
k_s^-7/2<=Lambda^-3/2(b0 s)^-2 before integration.
Discarding this s dependence prematurely would give an artificial
logarithmic divergence.

## Fractional functions and a common trace-class constant

For the functions used in the endpoint work theorem,

```
f(x)=x(1+x^2)^(-7/4), a(x)=(1+x^2)^(-7/4),
f(h)=c_* integral_0^infinity t^-3/4 g_s(h)dt,
a(h)=c_* integral_0^infinity t^-3/4 b_s(h)dt,
s=sqrt(1+t), c_*=1/Beta(1/4,7/4)<=1.                (10)
```

The beta integral giving c_* is exact. For the stated inequality its
integral on 0<t<1 alone is greater than
(1/4)integral_0^1 t^-3/4 dt=1. The scalar functional-calculus
integrals converge in operator norm. Bounds (9) also give absolute
trace-norm convergence for their commutators.

For completeness the needed estimates can all be obtained by splitting
at t=1. For p>1/4,

```
integral_0^infinity t^-3/4(1+t)^(-p)dt
 <=4+1/(p-1/4).
```

Thus the integrals with p=1/2,3/4,1 are respectively at most 8,6,16/3.
For the local terms use k_s^2>=Lambda^2+b0^2 t, rescale
t=Lambda^2 v/b0^2 and use the same split:

```
integral t^-3/4 k_s^-3 dt <10 Lambda^-5/2,
integral t^-3/4 k_s^-4 dt <10 Lambda^-7/2.
```

Indeed the remaining integrals are Beta(1/4,5/4)<24/5 and
Beta(1/4,7/4)<32/7, multiplied by b0^-1/2<2.
Combining (9)-(10), with Lambda>=1, proves BOTH estimates

```
||[f(h_n),K]||_1 <= C_n,
||[a(h_n),K]||_1 <= C_n,
C_n=10000 Lambda_n^-5/2+100 Z_n Lambda_n^-2
                         +1000sqrt(Z_n)Lambda_n^-3/2. (11)
```

For example the second line of (9) integrates to less than
3500 Lambda^-7/2+32 Z Lambda^-2
+600sqrt(Z)Lambda^-3/2, within (11). The large constants in (11)
are fixed analytic slack, not fitted numerical errors.

## Exact full angular sum

The retained exact sphere counting inequality is

```
D(Y)=sum_(M_n^2<=Y)d_n<=K_ang Y,
K_ang=2C(1+1/q), K_ang/q<15/k.
```

For p>2, integration of this positive counting measure gives
sum d_n M_n^-p <=K_ang p M_1^(2-p)/(p-2). For p>=1 and c>0,

```
sum d_n M_n^-p exp(-cM_n)
 <=K_ang exp(-cM_1) M_1^(2-p)[1+2/(cM_1)].          (12)
```

To verify (12), integrate -D(M^2) times the derivative of
M^-p exp(-cM), then integrate M^2 times that derivative explicitly
once. The remaining integral of M^(1-p)exp(-cM) is bounded by
M_1^(1-p)exp(-cM_1)/c. No angular truncation occurs.

At k>=2000 and M_1>=18 these estimates imply

```
(1/q)sum d_n M_n^-5/2 <1/100,
(1/q)sum d_n Z_n M_n^-2 <1/50,
(1/q)sum d_n sqrt(Z_n)M_n^-3/2 <1/400.              (13)
```

For the last inequality use sqrt(18)<5,
1+4sqrt(6)/18<2 and exp[-18/(2sqrt(6))]<1/32.
The last exponential inequality follows from sqrt(6)<5/2 and
sum_(j=0)^7 (18/5)^j/j!>32. The right side of (12) decreases for
M_1>=18 in these cases. For the middle inequality it is enough to
discard its exponential and bound its bracket by two.

Finally 3sqrt(2)<5, delta<1e-32 and (11)-(13) give the complete tower

```
(1/q)sum_(n>=1)d_n C_n
 <100(5delta)^(5/2)+50delta^2+30delta^(3/2)
 <1e-45.                                           (14)
```

The last line is deliberately loose: its largest displayed term is
less than 3e-47. This is a bound on the two specified static trace
commutators in the physical classical RN reference, uniform over the
original preparation range. The full-pulse endpoint argument must still
use them; smallness of the bounded domain-frame perturbation alone
would not prove finite quantum work. Real-axis resonance-free transfer,
a spectral gap and a zero-energy inverse moment have not been assumed.

Remaining to-do list: transfer the result to the quantum-corrected
evolving exterior and prepared
state, and match controller contacts and physical stress.
