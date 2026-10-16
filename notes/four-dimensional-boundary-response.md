# Four-dimensional paired-domain source response and the metric obstruction

2026-10-15, E003/S133, supporting S074. The actual matched paired
fermion boundary angle can be transferred to a fixed domain, with an
explicit bulk insertion. In four dimensions the appropriate generally
available integer regularization is det_5. This gives a conditional
response bound for the realized angle on a frozen geometry. A uniform
lapse counterexample shows why this is still not a full clock-stress bound.

## Predeclared scope and tests

Retain the minimal massless Dirac field, paired spin/gauge frames and
identical pulled-back collars of [the domain note](coupled-fermion-domain.md)
and [the curved-boundary note](curved-fermion-boundary.md). No LLL
truncation is made in the collar identity. The compact, invertible
Euclidean reference specified below is an explicit additional hypothesis;
no physical exterior, infrared state, heat constant or spectral gap is
silently assigned to the wormhole.

Before checker execution fix: matrix/domain/Clifford and lapse identities
must agree within 2e-11; independent finite differences of the differential
transfer within 2e-7; det_5 first derivatives within 2e-6 and mixed Hessians
within 2e-5 at steps 2e-4 and 1e-4; trace-series errors must lie below
their analytic bounds plus 5e-11 rounding allowance. Inequalities permit
relative rounding margin 5e-10. Omitting the normal extension derivative
must fail the transfer by more than 1e-4 in the declared test. Fourier
cutoffs are 2,4,8,16; they supplement, rather than establish, the infinite
dyadic argument. The lapse example uses N=3/2 and temporal momenta
1,2,4,8,16, with a fixed unitary identification of L2 measures.

## Exact domain and codomain transfer

Write D=Gamma^i nabla_i in the convention without an overall i. On the
folded collar use the rank-eight Clifford representation and parallel
anti-Hermitian G from the curved-boundary note:
G^2=-1, {G,Gamma^i}=0. Fix Pi0=(1+Gamma_n)/2. For boundary angle
theta(x), Pi_theta=U Pi0 U^dagger at the boundary, where

```text
alpha(r,x)=f(r) theta(x), U=exp(alpha G/2),
f=1 near r=0, f=0 before the paired collars meet their interiors. (1)
```

Choose a smooth fixed f, and set U=1 off those collars. This is a
global unitary on the original Hilbert space: it mixes paired collar
points only where the spin/gauge map has been specified, and is the
identity near their joins. It does not create an internal boundary.
The bundle generator need only be parallel on the support of alpha.
The extension is normally constant near the physical boundary, as in
the existing boundary calculation, but it cannot be normally constant
through the entire interior.

Let D_theta denote the same bulk differential expression with allowed
trace Pi_theta. Right multiplication by U maps Dom(D_0) onto
Dom(D_theta). The adjoint allowed trace is

```text
Pi_adj,theta=1-Gamma_n Pi_theta Gamma_n
            =U^dagger(1-Pi0)U.                            (2)
```

Indeed Gamma_n U Gamma_n=U^dagger. Thus U^dagger maps the reference
adjoint domain onto the physical adjoint domain. The fixed-domain
operator is a left/right transfer, not unitary conjugation:

```text
T_theta=U D_theta U =D_0+W_alpha,
W_alpha=(1/2) Gamma^i G partial_i alpha,
T_theta^dagger=U^dagger D_theta^dagger U^dagger.           (3)
```

To prove (3), use U Gamma^i=Gamma^i U^dagger and
U^dagger nabla_i U=G partial_i alpha/2. This also fixes the matrix
ordering and all signs. An overall i in a different Dirac convention
multiplies both D and W by that same i.

The domain of T_theta is exactly Dom(D_0), and its adjoint domain is
the reference complementary domain; W is bounded for smooth theta.
On the corresponding induced squared domains,
T_theta^dagger T_theta=U^dagger D_theta^dagger D_theta U.
The singular spectrum, and hence a consistently defined zeta modulus,
is therefore unchanged by this transfer. Determinant phase and fermionic
measure/Jacobian are not inferred from singular values.

The insertion contains

```text
W_alpha=(1/2)[ f Gamma^a G partial_a theta
                         +f' theta Gamma_n G ].          (4)
```

The normal f' theta term is essential, even for constant boundary
theta. Dropping it would falsely turn the known local constant-angle
equivalence into global spectral equivalence; the interval spectrum
in the domain note explicitly changes with a constant boundary angle.
Different smooth extensions f give equivalent full singular problems.
An isolated regularized remainder need not be extension independent:
its lower-order matching terms must be transformed with it.

## Actual realized-angle insertion on the frozen geometry

The retained radiation/joint clock has
rho(t)=r0(t)+kappa z(t), theta_z(t)=theta_initial+integral_0^t rho(s)ds.
Use that actual source dependence in (1), while freezing the Euclidean
background, spin frame, f and time coordinate. For a real test variation
h of z, define P h(t)=integral_0^t h(s)ds. Then

```text
delta W[h]=(kappa/2)[ f Gamma^t G h(t)
                             +f' Gamma_n G P h(t) ],
delta^2 W[h,j]=0.                                        (5)
```

On a product collar with orthonormal time and 0<=t<=H,

```text
||delta W[h]||op <= (|kappa|/2)
       (||f||infinity+H ||f'||infinity) ||h||infinity.    (6)
```

This follows from the triangle inequality and ||Gamma^i G||op=1.
Equation (6) is a window bound. For a compact Euclidean reference,
choose a smooth temporal completion explicitly. For example, extend
h in C_infinity,c(0,H) by zero, continue P h constantly after H, and
multiply P h by a fixed smooth eta with 0<=eta<=1, which is one on [0,H] and zero
before the ends of a larger periodic time chart. The global bound is

```text
||delta W[h]||op <= C_completion ||h||infinity,
C_completion=(|kappa|/2)[||f||inf(1+H||eta'||inf)
                                          +H||f'||inf]. (6a)
```

The additional eta' P h insertion outside the window cannot be
discarded. The simple global version of (6) applies instead when h
has zero integral and compact support in (0,H), so P h already returns
to zero and no compensating completion is needed. Complete the base
theta smoothly as well; its late-time constant need not be periodic
with its initial value. Completion is reference data, not a derived
Euclidean state for the returned Lorentzian switch.

For a general fixed collar, replace the first coefficient by the
supremum norm of Clifford multiplication by dt and retain the fixed
metric weight in the radial coefficient. Formula (5) supplies the
previously absent domain/source map, including memory through P h.
It does not identify a Lorentzian clock history with a real Euclidean
source: continuation and the physical state remain separate tasks.

## Four-dimensional ideal and explicit infrared hypothesis

For this subsection assume a smooth compact Euclidean realization of
D_0, the common first-order boundary domain above, and a bounded
two-sided inverse. Put A=D_0^dagger D_0 and assume

```text
A>=mu^2>0,
Tr exp(-t A)<=C t^(-2), 0<t<=1.                          (7)
```

The numerical value of C and the gap mu must come from the chosen
reference problem; neither is supplied by the local collar algebra.
For t>=1 the trace is <=C exp[-mu^2(t-1)]. Mellin integration yields

```text
||D_0^(-1)||5^5=Tr A^(-5/2)
 <= B5^5 := C/Gamma(5/2) [2+exp(mu^2) mu^(-5)
                                  Gamma(5/2,mu^2)].     (8)
```

More generally (7) gives S_p for p>4. For W_alpha and its source
variation, K=D_0^(-1)W_alpha and J_h=D_0^(-1)delta W[h] obey

```text
||K||5<=B5 ||W_alpha||op,
||J_h||5<=B5 C_completion ||h||inf.                       (9)
```

Invertibility, temporal compactification and (7) have not been proved
for the physical wormhole. A local heat coefficient does not establish
its infrared gap. Projecting away zero modes is not automatically
legitimate when the insertion mixes them with other modes.

The order-minus-one four-dimensional inverse does not generically
belong to S3 or S4 after a nonzero bounded source insertion. A direct
dimensional counterexample makes the issue decisive without importing
a full boundary spectral theorem. On a flat four-torus, omit k=0 and
let D(k)=i Gamma.k for k in Z^4\{0}. For constant W=c Gamma^t G, c>0,
each mode of D^(-1)W has eight singular values c/|k|.
For R=2^j, the disjoint positive blocks R<=k_i<2R contain R^4
points and have 2R<=|k|<4R. Their p-power sums are at least
8 c^p 4^(-p) R^(4-p). Therefore the S3 and S4 sums diverge.
For S5, the max-norm shell R<=||k||inf<2R has at most 256R^4
points and |k|>=R; its sum is at most 2048c^5/R. Summing the
dyadic shells from R onward gives an explicit tail <=4096c^5/R.
This example is a source-ideal test, not the wormhole or its boundary.
No LLL det_3 estimate is promoted to a full four-dimensional estimate.

## det_5 response and cutoff removal

For K,J,H in S5, write s=||K||5, a=||K||op<1 and R=(1+K)^(-1).
The analytic logarithm based at zero is

```text
F5(K)=log det_5(1+K)=sum_(n>=5)(-1)^(n+1)Tr(K^n)/n,
DF5[H]=Tr(K^4 R H),
|DF5[H]|<=s^4 ||H||5/(1-a),
D2F5[J,H]=Tr([(sum_(l=0)^3 K^l J K^(3-l))R
                                     -K^4 R J R]H),
|D2F5[J,H]|<=s^3(4-3a)||J||5||H||5/(1-a)^2.             (10)
```

Proof. Differentiate the degree-n trace, keeping the order of all
noncommuting factors. The first derivative has one direction and n-1
copies of K; place that direction and
four copies of K in S5 and bound remaining copies in operator norm.
Its bound is s^4||H||5 a^(n-5). A second derivative has n-1 ordered
terms containing J,H and n-2 copies of K. Place J,H and three copies
of K in S5. The sum of their bounds is
s^3||J||5||H||5 sum_(n>=5)(n-1)a^(n-5), yielding (10).
Local uniform convergence proves Frechet differentiability and the
symmetry of the Hessian. Differentiation of K^4 R gives its displayed
ordered formula. This is the same established determinant mechanism
as S132, with the exponent required by the actual dimension.

For a C2 path add DF5[K''] to D2F5[K',K']. For the linear source map
(5), K''=0 under a straight source variation. Equations (9)-(10)
therefore give an explicit conditional bound for the actual frozen-
geometry angle contribution to the clock variation.

For truncation at degree N>=5 the value, first and mixed second tails
are respectively bounded by

```text
s^5 a^(N-4)/[(N+1)(1-a)],
s^4 ||H||5 a^(N-4)/(1-a),
s^3 ||J||5||H||5 a^(N-4)[N/(1-a)+a/(1-a)^2].            (11)
```

Increasing finite-rank orthogonal compressions converge in S5 for K
and the directions. They contract operator norm. Every finite ordered
trace then converges by Schatten Holder; the uniform tails (11) give
convergence of F5,DF5,D2F5. On compact C2 source paths, finite-net
approximation makes this convergence uniform, as in S132. This removes
the fermion cutoff under (7); it removes no physical infrared limit.

In finite dimensions,
log det(1+K)=F5+Tr K-Tr K^2/2+Tr K^3/3-Tr K^4/4.
In four dimensions those first four terms require common regulated,
renormalized matching, including boundary terms and finite conditions.
They are not set to zero. det_5 begins at fifth order; it cannot by
itself supply the retained cubic and quartic physical responses or
resolve the existing source counterterms. Their regulator and possible
transfer anomaly/phase must be handled consistently with the heat
calculation and with extension independence.
The first four Born terms contain finite nonlocal pieces as well as
local ultraviolet counterterms. A local heat subtraction alone does
not reconstruct them.

## Why a full clock/metric change is a different problem

On the same flat torus take g_N=N^2 dt^2+dx^2, with N positive and
constant. The unitary half-density identification of spinor L2 norms
gives D_N=i(Gamma^t k_t/N+sum_a Gamma^a k_a) on each Fourier mode.
The coordinate time period is held fixed, so varying N changes its
proper circumference. This is a genuine metric-family comparison;
it is not a coordinate relabeling with the period transformed as well.
For K_N=D_1^(-1)(D_N-D_1) and k=(n,0,0,0), n!=0,

```text
K_N(k)=(N^(-1)-1) I_8.                                  (12)
```

Thus an infinite orthonormal temporal sequence has image norm
|N^(-1)-1|. When N!=1 the relative operator is not compact and
belongs to no finite Schatten ideal. This is an exact obstruction to
the naive fixed-reference det_p construction for metric variation;
changing det_3 to det_5 does not cure it. Source relative operators
D_N^(-1)W remain order minus one at each fixed geometry, but varying
the full determinant also varies the reference geometry/measure and
its subtraction. One needs principal-symbol cancellation under a
justified identification or an appropriate relative elliptic
determinant/heat construction. A conformal or diffeomorphism identity
must be proved for the actual family, not presumed from this source map.

There is a useful conditional interface for genuine geometry paths.
Assume the actual transport has produced one fixed graph domain X,
D_z:X->L2 is C2 with a bounded inverse, W_z is C2 in operator norm,
and one (hence locally all) inverse belongs to S5. Set

```text
R_z=D_z^(-1), Q=||R_z||5, b1=||D'_z R_z||op,
b2=||D''_z R_z||op, w_j=||W_z^(j)||op, K_z=R_z W_z.
||K'_z||5<=Q(b1 w_0+w_1),
||K''_z||5<=Q[(2b1^2+b2)w_0+2b1w_1+w_2].               (13)
```

The inverse identities R'=-R D' R and
R''=2R D' R D' R-R D'' R prove (13) by the ideal inequality.
Local inverse continuity follows in the fixed graph norm; factoring
R_z through one reference inverse then proves local C2 dependence
in S5 as well. Equations (10) and (13) consequently control the source
remainder along a geometry path if all these hypotheses and a uniform
operator-norm radius less than one are established. The noncompact
full-reference ratio (12) does not make geometry response impossible.
It separates that ratio from the source-relative path. The full
background determinant Gamma[D_z], the first four matched orders,
metric measure, mixed counterterms and phase still need their own
variations. No actual wormhole graph transport or bounds b1,b2 have
been supplied by this conditional interface.

## Source provenance and decision

[Bintanja--Freivogel--Rolph v2](https://arxiv.org/html/2308.00871v2),
section 2.1.2, equations (2.8)-(2.12), gives the physical paired boundary
condition. The source does not supply the driven transfer proved here.
[Vassilevich v3](https://arxiv.org/html/hep-th/0306138v3), equation
(2.21) and section 5.3, supplies the compact mixed heat expansion
framework, with the induced domain already derived in the local notes.
The Mellin bound (8) is proved here under the explicit gap hypothesis.
[Britz et al. v3](https://arxiv.org/html/2007.12834v3), equation (1.3),
defines det_5 with A=-K. [Ponge's published theorem 6.1, p.748](https://ems.press/content/serial-article-files/31026)
gives the order-minus-m singular-value law with p=n/m on closed
manifolds. Its scope supports the dimensional comparison, not a
spectral theorem for this physical exterior; the explicit torus argument
above suffices for the stated counterexample. No novelty claim is made.

The new checker independently differentiates local spinor jets, tests
both trace-domain maps and their adjoints, checks singular-value
equivalence under left/right unitaries, verifies the memory insertion,
compares det_5 ordered responses against finite determinants, and
checks the Fourier and lapse identities. Its receipt and raw log retain
the numerical inputs, exact command, tolerances and normalized source
hashes. No Lean, new physical state, opening or new clock interval is
claimed.

Remaining to-do list: obtain the physical Euclidean/retarded reference,
gap and heat constants; match orders one through four and transfer
measure/phase terms; derive the complete metric/domain variation and
unequal-collar response; control the interacting state and propagate
the complete remainder to the S076 clock, energy and area observables.
