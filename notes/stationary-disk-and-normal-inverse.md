# Global disk equations and an invertible coupled outer normal family

2026-10-17 HST, E003/S139 supporting S074/S075. **PROGRESS.**

There is a global fixed-period disk gauge in which two scalar Euler
equations recover the full static spherical tensor equation at a smooth
bolt. The [mixed interior parametrix](stationary-shape-logarithmic-response.md)
extends through that bolt. At infinity, the actual reflecting Dirac
normal problem has an explicit Mellin symbol, including its image and
Weyl contact. The complete coupled normal family is invertible on weight
7/2 in the retained benchmark. These are endpoint and normal statements;
bounded actual-minus-normal remainders and a global inverse remain unproved.

Keep the identical stationary copies, fixed flux, positive prepared
temperature, boundary clock, mu0=1/b_classical and zero additional finite
local coefficients of the [retained model](uniform-dirac-reference.md).
Set lambda=lambda_e, ell=lambda^(-1/2), epsilon=G/(60pi),
a=1+11epsilon lambda/3. The centered vacuum obeys
lambda0=lambda+11epsilon lambda^2/6. No physical gate is closed.

## 1. A global smooth disk, including a moving horizon

Write the Euclidean comparison in proper radial distance as

```
g=d rho^2+a_orb(rho)^2 d varphi^2+r(rho)^2 dOmega^2,
varphi=kappa tau, period(varphi)=2pi,
a_orb=exp(alpha)sqrt(F)/kappa,
a_orb=rho+O(rho^3), r=r_h+O(rho^2).
```

Define a globally normalized isothermal radius

```
xi=exp[-integral_rho^infinity d rho'/a_orb(rho')]
  =exp[-kappa integral_r^infinity exp(-alpha)/F dr'],
Phi=log(r/B), V=log[B a_orb/(r xi)].                  (1)
```

The outer integral converges because a_orb grows exponentially with
proper AdS distance. At the bolt, subtracting 1/rho from 1/a_orb leaves
a smooth odd function. Thus xi=c rho+O(rho^3), c>0. It increases from
0 to 1 and gives exactly

```
g=exp(2Phi)[exp(2V)(d xi^2+xi^2 d varphi^2)+B^2 dOmega^2]. (2)
```

Both Phi and V are smooth functions of xi^2 at the origin. The endpoint
xi=1 fixes radial dilation; no artificial value of Phi or V is fixed at
the bolt. Smooth parameter families with fixed kappa pull back to this
one disk. The regular part of the integral in (1) varies smoothly on
finite regions, so a moving physical horizon is included in that pullback.

Conversely, where Phi_xi>0, recover the areal data by

```
r=B exp(Phi), F=B^2 exp(-2V) Phi_xi^2,
alpha=log(kappa xi/B)+Phi+2V-log(Phi_xi).             (3)
```

For Phi=Phi_h+Phi_2 xi^2+..., Phi_2>0, and smooth V, these give a simple
areal horizon and exp(alpha_h)F'_h/2=kappa. The apparent logarithms in
alpha cancel at xi=0. These facts concern the nonextremal comparison
neighborhood, not an extremal limit or geometries with a second turning
point of r.

This also explains the pole in fixed-u shape variations of the previous
Schwarzschild representative. At a moving bolt,
delta f_h=-2kappa delta u_h and delta f'_h=-f''_h delta u_h. In y=u-u_h,
the translated pullback has delta fhat=O(y^2), hence delta fhat/f=O(y),
smooth on the disk since y=O(rho^2). The fixed-u term
-delta u_h/(u-u_h) was a coordinate singularity. Physical horizon data
still obey delta log r_h=delta sigma_h+delta u_h/u_h; this grants no free
parameter of full solutions.

For x=-log xi, the fast phase gives

```
x=kappa[1/(lambda r)-1/(3lambda^2 r^3)
                         +m/(2lambda^2 r^4)+O(r^-5 log^p r)],
Phi=log[kappa/(B lambda x)]-lambda x^2/(3kappa^2)
                         +m lambda^2 x^3/(2kappa^3)+O(x^4 log^p x),
V=log(B sqrt(lambda)/kappa)+x+lambda x^2/(2kappa^2)
                         -m lambda^2 x^3/kappa^3+O(x^4 log^p x). (4)
```

These expansions follow by integrating exp(-alpha)/F and reverting the
series. They preserve clock and period. In the product orbit cylinder,
psi=V-x has psi'(0)=0 and fixed psi(0), so the limiting boundary remains
geodesic with its actual mixed Dirac data. There is exactly one mass
cubic: delta V_3+2 delta Phi_3=0. It must not be replaced by two unrelated
asymptotic parameters.

## 2. Exact reduced equations and the bolt parametrix

Let E_rho,E_varphi,E_theta be mixed components of the full tensor residual.
Tangents S=delta Phi,T=delta V have orthonormal metric tangent
diag(2S+2T,2S+2T,2S,2S). Per physical volume their Euler gradients are

```
F_S=(E_rho+E_varphi+2E_theta)/(8pi G),
F_T=(E_rho+E_varphi)/(8pi G).                         (5)
```

Use the complete action and its actual state/domain variation. At an
exact zero of both scalar equations, E_theta=0 and E_varphi=-E_rho.
The Ward identity is

```
E_rho'+(a_orb'/a_orb)(E_rho-E_varphi)
       +2(r'/r)(E_rho-E_theta)=0,
(a_orb^2 r^2 E_rho)'=0.                              (6)
```

Smoothness at a_orb=0 forces the integration constant to vanish. Thus
(5)=0 is equivalent to the full tensor equation in this regular gauge.
The [earlier anisotropy reduction](stationary-weighted-conformal-range.md)
still requires its separate C/r^4 constraint: orbit trace and anisotropy
are different lower outputs. This result does not discard that constraint
while retaining the old output. At a nonzero-defect comparison, derivatives
of (5) include density, pullback and defect terms; a candidate Newton
operator cannot be replaced by an exact-solution derivative.

In Cartesian disk coordinates, for every orbit-plane covector the
accepted rank-four symbol has b=2T and theta=S+T/3. Its principal form is

```
p^4 theta^2/(80pi^2)
 + p^4 T^2[log(p/reference)-6/5]/(120pi^2),            (7)
```

with the retained ordinary shape coefficient understood. Smooth disk
coefficients and the same lower-contact analysis give a disk parametrix
and the mixed Garding estimate through xi=0. Average the parametrix over
rotations to preserve invariant fields. Norms here are smooth two-dimensional
disk Sobolev/Holder norms, with radial measure xi dxi, not unweighted
one-dimensional u norms. The bolt is not an added reflecting boundary.

## 3. The actual reflecting flat normal problem

Flatten the outer product collar and dilate its normal and three tangential
coordinates together. The limiting physical metric is
g0=(ell/z)^2(dz^2+sum_i dx_i^2), z>0. Curvature of the sphere, magnetic
curvature, thermal windings and the O(z^2) metric jets have lower dilation
degree. The limiting projector is chi=+/-gamma_z on each identical copy,
with complementary first-order adjoint domains as previously established.

For a tangential transverse-traceless (TT) radial family take

```
g=exp(2Omega) gbar, Omega=log(ell/z),
gbar=dz^2+sum_i exp(2beta_i(z)) dx_i^2,
(beta_1,beta_2,beta_3)=(b/3,-b/6,-b/6), sum beta_i=0.   (8)
```

This is TT at first order with respect to AdS, including covariant
divergence. Its determinant is independent of b. In the fixed spinor
density the exact Dirac expression and square are

```
D=gamma_z partial_z+sum_i gamma_i exp(-beta_i)partial_i,
Q=-D^2=-partial_z^2-sum_i exp(-2beta_i)partial_i^2
           -sum_i gamma_z gamma_i(exp(-beta_i))'partial_i. (9)
```

The spin connection's sum beta_i'/2 vanishes. No extra R/4 potential
should be appended to this square. K=sum beta_i'=0, so the actual square
domain is constant complementary Dirichlet/Neumann. Extend beta_i evenly.
Reflection times chi commutes with Q. Four-dimensional chirality commutes
with Q and its heat kernel and anticommutes with chi. Consequently the
reflected spin trace vanishes throughout this metric family, not merely
at the flat background. Thus Tr_half exp(-tQ)=Tr_double exp(-tQ)/2.
Differentiating this identity retains all contact vertices. Restricting
the doubled Hessian to positive z gives its direct PLUS image kernel,
with no extra factor 1/2 in the radial response.

This agrees with the separated-point free-fermion image construction in
[Herzog--Huang, section 5.2](https://arxiv.org/html/1707.06224).
The determinant argument above fixes the realization and its contacts;
the separated-point correlator alone would not do so. Initially use
smooth admissible tests. Even power families below are obtained by
Mellin continuation, not an assertion that all even doubles are smooth.

The flat TT determinant symbol, per Dirac copy, follows from the
[previously derived zeta Hessian](stationary-shape-logarithmic-response.md):

```
d0 partial_z^2[log(|D_z|/mu0)-6/5]partial_z^2,
d0=1/(480pi^2),                                     (10)
```

acting on the even extension. The ordinary constant is -6/5 in the TT
channel; -31/30 belonged to the non-TT fixed-sigma shape channel.

## 4. Weyl contacts and the full TT Mellin symbol

The Euclidean anomaly is -c C^2+a_E E4-d Box R, where
c=1/(320pi^2), a_E=11/(5760pi^2), d=1/(480pi^2).
Put U=sum beta_i'^2=b'^2/6. For compact TT variations the retained
Weyl and Einstein quadratic actions are determined as follows:

```
W_C, Hessian = -d0 partial_z^2(Omega partial_z^2),
(sqrt(g_t) E4[g_t])_(beta^2)=-4t(Omega' U)',
W_Euler, quadratic = 2a_E integral Omega'^2 U dz,
W_BoxR, quadratic = 2d integral U/z^2 dz,
I_Einstein, quadratic = ell^2/(16pi G) integral U/z^2 dz. (11)
```

Here g_t=exp(2tOmega)gbar. The Euler identity follows from the warped
diagonal curvature contractions; integration of Omega times that
identity in t and z gives the displayed sign. The Box R primitive is
+(d/12)[integral sqrt(g)R^2-integral sqrt(gbar)Rbar^2]; its quadratic
term uses Omega''+Omega'^2=2/z^2. The cosmological volume has no TT
variation in (8). The fixed magnetic field vanishes in the normal limit.
Boundary functionals are retained in the actual problem; these compact
normal variations have no extra boundary support to replace (11).

For an even power, differentiating
Fourier[|z|^u]=-2 Gamma(u+1)sin(pi u/2)|k|^(-u-1)
under analytic continuation gives

```
log|D_z| |z|^u
 =[digamma(u+1)+(pi/2)cot(pi u/2)-log z]|z|^u.        (12)
```

Take u=s-2 and combine (10) with the first line of (11) BEFORE applying
the outer two derivatives. The log z terms cancel exactly. Define

```
q(s)=s(s-3), P4(s)=q(s)[q(s)+2],
F_M(s)=digamma(s-1)+(pi/2)cot(pi s/2),
kappa_M=-6/5-log(mu0 ell),
A_T=ell^2/(48pi G)+23/(8640pi^2).
N_T(s)=d0 P4(s)[F_M(s)+kappa_M]-A_T q(s).            (13)
```

The coordinate-density Hessian maps z^s to N_T(s)z^(s-4). Its physical
volume gradient has coefficient lambda^2 N_T(s) on z^s. The last three
lines of (11) give -A_T q, including their Hessian factor two.

The reflection identity for digamma proves F_M(3-s)=F_M(s), hence the
correct weighted formal symmetry. F_M(0)=F_M(3)=1-gamma_E by removable
limits, so s=0 and s=3 are roots. At s=1,2 one must take the product:
P4 F_M tends to -2, not zero. The mass power s=3 is especially instructive:
the flat even logarithm supplies -6/z and the Weyl contact +6/z. Their
cancellation recovers the accepted mass boundary-layer cancellation.
There is a genuine Mellin pole at s=4; it is not another permissible
root obtained by separately evaluating factors.

For the normal response to a cut-off z^s, the homogeneous radial kernel
has degree -5. Its integral at the remote end converges for Re s<4 and
at the boundary for Re s>-1, after the coincident finite part. Removing
the far tail contributes physical O(z^4), so it cannot change the
coefficient (13) in this strip. The mixed Dirac Green parametrix has
the direct-plus-reflected flat kernel as its leading dilation term.
Two stress insertions have degree -8; tangential integration gives
degree -5. Lower geometric jets, connection, magnetic and winding terms
gain powers under dilation. Local contact terms of this degree are fixed
by the same common Mellin prescription. This identifies (13) as the
outer normal family of the retained problem. It does NOT yet bound the
full remainder as an operator on the proposed weighted spaces. In
particular state-dependent lower contacts are not all declared smoothing.

## 5. A proved zero-free weighted contour

Use delta=7/2, s=7/2+it. Digamma recurrence and reflection give exactly

```
Im F_M(s)=-t S_t,
S_t=1/(t^2+1/4)+1/(t^2+9/4) <=40/9,
q+2=15/4-t^2+4it.                                  (14)
```

The following elementary bounds avoid numerical special-function tests:

```
Re F_M(7/2+it)>-2,
Re F_M(7/2+it)<4+log(1+2|t|).                       (15)
```

For the lower bound, Re digamma(5/2+it)>=digamma(5/2)
=-gamma_E-2log2+8/3>1/6, and the real cot term is
-(pi/2)sech(pi t)>-2. For the upper bound use
Re digamma(1/2+it)=digamma(1/2)+sum_{n>=0} f_t(n+1/2),
where f_t(y)=t^2/[y(y^2+t^2)] is nonnegative and decreasing. The sum
lies between its integral, log(1+4t^2)/2, and that integral plus 2.
Adding the two recurrence terms, whose real parts sum to at most 8/3,
and dropping the negative cot term proves the stated loose upper bound.

The [retained benchmark](stationary-same-source-boundary-charge.md) gives
G lambda<k^-10, k=2000, lambda b_classical^2>k^3/4 and
|log(mu0 ell)|<12. Since e<3 and 3^18<2000^3/4,

```
9<-log(mu0 ell)<12,
39/5<kappa_M<54/5,
A_T/d0=10pi/(G lambda)+23/18>30k^10>144.             (16)
```

Set H=F_M+kappa_M. Equations (14)-(16) imply Re H>29/5 and, for t!=0,

```
Im[(q+2)H]/t
 =4 Re H-(15/4-t^2)S_t
 >116/5-50/3=98/15.                                (17)
```

If the bracketed factor is negative the subtracted term is positive;
otherwise the displayed bound follows from S_t<=40/9. Thus
-A_T/d0+(q+2)H cannot vanish for t!=0. For |t|<=1, (15) and log3<2
give Re H<84/5, and

```
Re[(q+2)H]=(15/4-t^2)Re H+4t^2 S_t<63+8=71<72.
```

Consequently the TT bracket has modulus greater than A_T/d0-72 on
|t|<=1 and greater than 98|t|/15 on |t|>=1. Also |q|>=7/4.
This proves a genuine nonzero normal gap on the whole vertical contour,
not merely on its real intersection. No exclusion of all complex roots
in an entire strip is asserted.

## 6. The coupled normal inverse and the mass extension

Decompose the disk tangent (S,T) into a full conformal tangent tau,
the TT tangent b of (8), and an infinitesimal radial diffeomorphism
eta partial_z. Comparing metric components gives

```
b=2T, eta'=b/3, eta(z)=(1/3)integral_0^z b,
tau=S+b/6+eta/z,
[tau,b]^T=C(s)[S,T]^T,
C(s)=[[1,(s+3)/(3(s+1))],[0,2]].                    (18)
```

At the exact centered AdS vacuum, diffeomorphism directions annihilate
the Hessian, and the conformal/TT cross term vanishes by the trace
identity and delta R_TT=0. The scalar coordinate-density normal symbol is

```
J_c(s)=3[a+epsilon lambda q(s)][q(s)-4]/(4pi G lambda).
H_disk(s)=C(3-s)^T diag(J_c(s),N_T(s)) C(s).          (19)
```

The adjoint reflection is 3-s because the raw density response shifts
degree by four. Using C(s)^T in its place would be incorrect. The
transformation poles -1 and 4 both avoid delta=7/2; det C=2.
Scalar roots are -1,4 and 3/2+/-i omega, with
omega^2=a/(epsilon lambda)-9/4>0. They also avoid this contour.
In fact the four root distances give |J_c|>=9/(80pi^2) there.
Together with (17), det H_disk=4 J_c N_T never vanishes.

At large |t|, J_c~t^4/(80pi^2) and
N_T~t^4 log|t|/(480pi^2), with the differentiated symbol bounds from
the digamma expansion. The bounded transformations (18) tend to the
constant variables theta=S+T/3 and b=2T, with order-minus-one errors.
Fourier transformation in log z therefore supplies a bounded two-sided
NORMAL multiplier inverse: ordinary four-derivative regularity for
theta, and four derivatives plus one logarithm for shape, into a common
weighted Sobolev output. The gap handles finite frequencies and the
symbol asymptotics handle large frequencies. This is an inverse of the
translation-invariant normal operator, not of the actual exterior.

The auxiliary delta=7/2 domain excludes the physical mass cubic at s=3.
One must split off its correlated datum from (4). It is pure TT:
S_3+(3+3)T_3/[3(3+1)]=S_3+T_3/2=0. Its cut-off extension must be included
in the actual global residual map and cannot be treated as an already
solved parameter; its required weighted bound remains to be proved.
The physical s=4 logarithmic tail also requires matching; neither the
normal pole nor the conformal resonance permits a fixed-log-power
Newton closure by declaration. The earlier nonzero anisotropy at the
conformal trace root remains consistent with the singular change of
normal decomposition at s=4.

The result now removes a spurious bolt boundary, retains the correct
number of physical asymptotic data, and proves both normal blocks have
an inverse on one common contour. To obtain an actual global Fredholm
map still requires bounded/compact actual-minus-normal remainders,
compatible global domains and mass extension, and a kernel/range audit.
Quantitative inversion on the preparation interval and nonlinear/state
matching are additional obligations. No residual is converted into a
distance to an exact solution, and no stability or nonexistence follows.

## 7. Focused verification contract

The ordinary argument precedes scientific execution. A focused checker
may verify the coordinate/temperature and Ward identities, series
reversion, disk symbol, even-reflection Clifford identities, Euler
curvature identity, exact normal polynomials, reflected congruence and
rational benchmark guards. Meromorphic Fourier continuation, heat
localization, disk estimates, digamma-series bounds and multiplier
mapping remain ordinary proofs. No numerical root shooting, fitted
kernel or global inverse test is selected. Freeze reviewed source bytes,
retain all runs and failed input snapshots, and hash raw UTF8/LF bytes.

Remaining to-do list: control actual-minus-normal weighted mappings and
the correlated mass extension; establish the global kernel/range and a
quantitative coupled inverse, physical boundary branch and logarithmic
phase, or justify controlled spatial reduction; match exterior/state
and charge, transfer S138, and complete driven stress, radius, controller,
capture, tail, spatial area and S074-S080.
