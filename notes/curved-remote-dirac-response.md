# Remote massive Dirac response on the actual curved exterior

2026-10-19, E003/S143 supporting S139/S074. Ordinary argument and
diagnostic protocol fixed before execution.

Two shape variations supported near the throat have exponentially small
massive Hessian dependence on a distant conformal variation. This holds
uniformly on a specified finite-background ball on the actual paired
exterior. It controls this separated configuration only. In particular,
it does not estimate the coincident curved three-point vertex or prove
contraction of the exterior equations.

## 1. Fixed operator, family and comparison

Use the retained centered seed, prepared positive temperature, magnetic
bundle and paired domain from
[the curved form](stationary-curved-quantum-tt.md), sections 1-3, and
[the canonical response](stationary-centered-canonical-response.md).
The normalization and proper angular gap are checked in
[the recovered benchmark](retained-benchmark-parameter-recovery.md):

    gbar = ghat_2/A^2, A=B/sqrt(C),
    M_n^2=n(n+q)/C, d_n=2(q+2n), n>=1,
    M_1>18, eta_G=GC/B^2<8.5e-33.

All distances, volumes and gradients below use this FIXED gbar. No
optical spectral gap or freely chosen temperature is substituted.
The actual even double has the inherited C2,1 seam and smooth bolts.
Its anti-self-adjoint orbit realization D0 has a preserved unitary
grading J with JD0J=-D0. Scalar Lipschitz weights preserve its domain.
These are properties of the paired representative, not of the original
unpaired Euclidean boundary operator. See
[the phase allocation](stationary-paired-dirac-phase.md) and
[the boundary-layer identity](stationary-massive-boundary-layer.md).

Use real, identical, static conformal paths on the two physical copies,
even at the outer seam, preserving its geodesic condition and the
original projector. Take piecewise smooth paths with the inherited
seam regularity, smooth at the bolt. On the fixed Hilbert space,

    D_v = exp(-v/2) D0 exp(-v/2),
    R_v = (D0 + m)^-1, m=M exp(v),
    Gamma_M[v] = -Tr log(D0+M exp(v)) + local terms.       (1)

The determinant equality here means equality of the nonlocal
variations. The Weyl Jacobian, common counterterms and finite local
restorations are retained. The paired massive phase has no varying
nonlocal contribution on this family. The physical functional is half
the evenly varied double functional; the rank-two orbit trace and d_n
already account for the sphere signs. We use double traces to prove
bounds, then divide by two. X,Y below denote the even extensions of
physical supports; Vol(X) in the final bound is the physical-half area.

The determinant convention agrees with the flat variable calculation
in [S142](flat-variable-dirac-three-point.md). The boundary projectors
are not replaced by new reflecting conditions; the original source
[Bintanja--Freivogel--Rolph, section 2.1.2](https://arxiv.org/html/2308.00871v2#S2.SS1.SSS2)
specifies the coupled boundary problem. That primary section was
re-inspected on 2026-10-19. It does not supply the bound proved here.

The preceding canonical response bounds a different source-parameter
derivative. S142 treats arbitrary momenta at a flat basepoint. Neither
already proves the following finite-background, separated form bound.
No improvement for a coincident or global response is claimed.

Fix rho=1/100 and a background v0 with

    ||v0||_infinity <= rho, ||grad v0||_infinity <= 1/100. (2)

Let u,w have support in a finite-area interior region X, and v have
support in Y. Their even extensions have distance d>=24 on the actual
double. Local coordinate disks of radius 1/(10M) around X must remain
inside its smooth interior. The cold-bolt localization below satisfies
this. Directions are smooth and bounded and preserve the stated
domain/path conditions. They need not have a small derivative for the
trilinear estimate; a finite integrated path must satisfy (2).

## 2. Resolvent and local evaluation bounds

For every real v0 in (2), accretivity and the same bound for the
adjoint give

    ||R_v0|| <= exp(rho)/M.

This is a bounded perturbation of the anti-self-adjoint realization:
the positive real part is at least M exp(-rho), giving injectivity,
closed range and, by the adjoint, surjectivity. For a bounded real
Lipschitz weight phi with |grad phi|<=M/2, scalar conjugation changes
D0 by Clifford multiplication of norm at most M/2. Thus

    ||exp(phi) R_v0 exp(-phi)|| <=
      1/[M(exp(-rho)-1/2)] < 3/M.                       (3)

Truncation of distance weights and a limit give the same off-diagonal
bound. This argument uses the actual global operator throughout;
no uncomputed boundary image remainder is set to zero.

Here is an explicit, deliberately loose local evaluation constant.
The seed's conformal normal charts give, on |z|<=1/(10M),

    gbar=exp(2 sigma)|dz|^2, sigma(0)=0,
    |sigma|<=.05, |grad_flat sigma|<=10.                 (4)

For non-bolt charts this follows from Hess sigma<5 on radius 1/16
in section 3 of the curved-form note. In the cold and compact analytic
charts the complex extension of exp(sigma)-1 is bounded by .026 on
the bidisk of radius at least 1/32. On radius 1/180, Cauchy in its
two variables bounds the real gradient of sigma by
2 sqrt(2)(.026)/[(1/32-1/180)(1-.026)]<3. The real value is below
.05. This uses smooth bolt coordinates, not the singular polar frame.
The radius 1/(10M)<=1/180 is inside these charts.

If (D0+m)psi=0 in that disk, conformal covariance gives
chi=exp(sigma/2)psi and

    (Dflat+m_flat)chi=0, m_flat=M exp(sigma+v0).

The gradient in (2) converts to a flat gradient by at most exp(.05).
Consequently |m_flat|<1.07M, |grad m_flat|<11M, and

    (-Delta+V)chi=0,
    V=m_flat^2-gamma.grad m_flat, ||V||<3M^2, M>=18.     (5)

Take R=1/(20M). A cutoff equal to one on B_R and supported in B_2R
has gradient at most 40M. Integration by parts and Cauchy--Schwarz
give the Caccioppoli estimate

    ||grad chi||_(B_R)^2 <= (4*40^2+6) M^2
                            ||chi||_(B_2R)^2
                         <81^2 M^2 ||chi||_(B_2R)^2.

For a second cutoff zeta, one on B_(R/2), zero outside B_R, use
|grad zeta|<=80M and |Delta zeta|<=16000M^2. For example the
radial quintic H(t)=1-10t^3+15t^4-6t^5 on each transition has
|H'|<=15/8 and |H''|<=10/sqrt(3)<6, which gives both cutoffs.
Equation (5) then implies

    ||(-Delta+M^2)(zeta chi)||_2
       < (4+2*80*81+16000) M^2 ||chi||_(B_2R)
       < 30000 M^2 ||chi||_(B_2R).

The whole-plane inverse has L2-to-point norm 1/(2 sqrt(pi) M),
by its Fourier integral. Spin components obey the same vector norm.
The conformal density conversion costs at most exp(.025) at a
normalized center (using exp(.05) is also sufficient). We may retain
the much larger reserve

    |psi(x)| <= 10^6 M ||psi||_(B_2R,gbar).             (6)

The disk B_2R has gbar radius below exp(.05)/(10M)<1/(9M).
For a source f supported in Y, R_v0 f is homogeneous in this disk.
Apply (3) with the distance to Y and (6). The row operator from
L2(Y) into the rank-two fiber at x has norm at most

    3*10^6 exp(1/18) exp(-M d/2)
      <3.2*10^6 exp(-M d/2).

The row's squared Hilbert--Schmidt norm is at most twice its squared
operator norm. Integration over X yields

    H_M^2 := ||1_X R_v0 1_Y||_HS^2
       <=2.1*10^13 Vol_double(X) exp(-M d).             (7)

No area of Y occurs. Since R_v0^*=J R_v0 J, scalar support projections
commute with J; the reverse cross block has the same HS norm. The
estimate also proves the required cross blocks are Hilbert--Schmidt.

## 3. The surviving mixed third variation

Differentiate (1) on the fixed Hilbert space. With R=R_v0,

    Gamma_M'''[u,v,w] =
      -Tr(R m uvw)
      +Tr(R m uv R m w + R m uw R m v + R m vw R m u)
      -Tr(R m u R m v R m w + R m u R m w R m v).

Since X and Y are separated, uv=vw=uvw=0. ONE seagull and BOTH
triangle orientations survive:

    Tr(R m uw R m v)
      -Tr(R m u R m v R m w + R m u R m w R m v).     (8)

Every term has two cross blocks. Their product is trace class by
(7); the third resolvent in each triangle is bounded by exp(rho)/M.
Cyclicity is legitimate for these products. Thus the double trace
is bounded by

    [exp(2rho)+2 exp(4rho)] M^2 H_M^2
         ||u||_infinity ||v||_infinity ||w||_infinity
    <6.72*10^13 Vol_double(X) M^2 exp(-M d)
         ||u||_infinity ||v||_infinity ||w||_infinity.   (9)

Here exp(2rho)+2 exp(4rho)<3.2. The physical-half division cancels
the factor two in Vol_double(X).

Every mixed third variation of a LOCAL bulk or boundary functional
is supported where all three directions' jets meet. The positive
distance makes that set empty; u,w also vanish near the outer
boundary. This includes the vacuum restoration, the common a2 term,
and the local Weyl allocation. They contribute exactly zero to (8),
at every background in (2), and have not been removed from the action.
One may first differentiate with a regulator and then remove it:
the separated products are trace class and the angular bound below
is summable uniformly in the whole ball. There is no local contact
ambiguity in this particular mixed variation.

## 4. Complete tower and finite Hessian change

For N(E)=sum_(M_n<=E) d_n, direct summation to integer N gives
2N(q+N+1). Since M_N^2=N(N+q)/C,

    N(E)<=4 C E^2.

For d>=24, f(E)=E^2 exp(-d E) is decreasing on E>=18.
There are no modes at or below 18. Stieltjes integration by parts,
using the polynomial counting bound and the vanishing upper limit,
gives

    sum_n d_n M_n^2 exp(-d M_n)
      <=4 C integral_18^infinity (d E^4-2E^3) exp(-dE) dE
      =4 C exp(-18d) P(d),

    P(d)=18^4+2*18^3/d+6*18^2/d^2+12*18/d^3+12/d^4.   (10)

Normalize Gamma_hat=(2G/B^2) sum_n d_n Gamma_Mn on the physical
half. Combining (9)-(10), uniformly in (2),

    |Gamma_hat'''[u,v,w]| <=
       5.376*10^14 eta_G Vol(X) exp(-18d) P(d)
         ||u||_infinity ||v||_infinity ||w||_infinity
      <2*10^-200 Vol(X)
         ||u||_infinity ||v||_infinity ||w||_infinity, d>=24. (11)

The tiny last number is specific to this support separation and
proper massive gap. The checker encloses it using rational arithmetic:
e>sum_(j=0)^10 1/j!, hence exp(-432) is bounded above by that
rational number to power -432. It is not a floating-point underflow
argument. The symbolic integral establishes the infinite tower,
not a numerical cutoff extrapolation.

If v0+t v, 0<=t<=1, stays in (2) and in the domain-preserving
family, integrate (11):

    |Gamma_hat''_(v0+v)[u,w]-Gamma_hat''_v0[u,w]|
      <2*10^-200 Vol(X) ||u||_infinity ||v||_infinity ||w||_infinity. (12)

This is a finite remote-background result. Only the path is small;
u,w are arbitrary smooth bounded directions with the stated supports.
The bound is in these fixed support/supremum norms, not in the global
cosh-weighted inverse norm.

The actual cold chart supplies a concrete nonempty support class.
As recorded in [the cold response](stationary-reduced-cold-response.md),
section 3, gbar=exp(2 sigma_cold)g_H2, |sigma_cold|<.05, and the
seam s=.01 occurs at 60<L_c<84. If X lies in chi<=1 and Y lies
beyond that seam, every path, including one crossing the doubled
outer boundary, first crosses this collar. Its length is at least
exp(-.05)(L_c-1)>56>24. The smooth-bolt cap contains all the small
evaluation disks required above. Vol(X) stays explicit.

## 5. Verification protocol and unresolved step

Run checks/curved_remote_dirac_response.py once after this argument
and the checker are fixed. It records input hashes, raw check rows and
any failed source snapshot. Exact checks cover cutoff derivatives,
constant reserves, the counting inequality, Stieltjes polynomial and
the rational exponential enclosure.

An independent finite-dimensional diagnostic uses a spatially coupled
anti-Hermitian graded matrix with a nonconstant exponential mass.
Two insertions share one spatial region; the third is disjoint.
At 70-digit precision, compare the ordered traces in (8) with a
central mixed third difference of -log det, at step 10^-4 and 5*10^-5.
The absolute error must be below 10^-7 and the second error below
.35 of the first (unless both are below 10^-45). Also compare a
finite Hessian change with a quadrature of (8), tolerance 10^-30,
and check grading, the resolvent norm and cross-block identities.
These are diagnostics of algebra and finite-background conventions;
they do not certify the PDE, its local charts or the continuum domain.

The ordinary proof controls remote massive conformal dependence with
the actual boundary present. It does not establish the local curved
three-point cancellation at finite backgrounds, the full Weyl-log
form mapping, or a general mixed shape/radius metric response. It
does not bound the LLL or a source that overlaps the throat. A far
perturbation of a solved exterior may alter constraints, state or
temperature; (12) alone does not solve those coupled effects.

Remaining to-do list: coincident curved finite-background response in
the actual exterior norm; coupled exterior/state, charge and signed
source; driven geometry, stability, finite-energy transport,
preparation and separate creation requirements.
