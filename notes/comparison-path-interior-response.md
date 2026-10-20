# Canonical interior response on the retained comparison path

2026-10-20, E003/S155 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The retained finite first-order comparison geometry has a controlled
canonical scalar/HH pressure Hessian on a fixed interior cover.
Its middle-region canonical constraint remainder is consequently
bounded. This is an application to a specified metric path, not a
proof that this path solves the full semiclassical equations.

## 1. Retained geometry and observation cover

Use the original common-disk completion of
[the reconstruction](stationary-centered-first-order-geometry.md):

    sigma_lambda=sigma0+lambda u,
    r_lambda=B+(r0-B)exp(lambda V), 0<=lambda<=1.   (1)

The physical clock, actual paired domain, bolt and asymptotic
normalization are fixed as proved there. The retained
[metric estimates](stationary-completed-metric-response-norm.md) give

    |u|<1.2e-8, |grad_gbar u|<7.7e-6,
    |u_chi|<3e-6, |V|<1.2e-8, |s V_t|<3e-7.      (2)

Let I be the fixed seed cover with inward proper distance
Y0>=.051. It includes the bolt and all of S152's middle
region 1.001<=t=r0/B<=40. Indeed the seed outer guard gives

    Y0(40)>=sqrt(5.9)/(40*1.001)>.0606.

For every point of I, its endpoint proper distance is at least
exp(-1.2e-8)Y0>.05. Thus the cold/middle endpoint chart
bounds apply there, using the larger common reserve

    |partial^alpha delta sigma_N|<=epsilon j!100^j,
    epsilon=1e-4, j=|alpha|<=5.                   (3)

The small shift in the cold/middle boundary causes no missing
region: their union supplies (3). We make no coverage claim
for the remaining strip Y0<.051 in this interior theorem.
The future global argument must cover it with the outer map.

The normal charts in (3) are endpoint charts. They are not
identified with fixed seed charts without the transfer below.

## 2. Invert the finite quadratic normalizer

At each observation use the seed normal chart, with
sigma0(0)=grad sigma0(0)=0. Put b=u for the endpoint lambda=1.
The real endpoint is normalized by the holomorphic quadratic map

    F(z)=a z+c z^2,
    a=exp(-b(0)), c=-a^2 b_z(0).

Only a coordinate map is holomorphic here; b is real C5.
Equation (2) and the unit seed metric at the center give

    |a-1|<2e-8, |c|<4e-6.

The inverse branch G with G(0)=0 is single-valued on |x|<=1/40
and lies inside |G|<1/32. For example on |z|=1/32,
|F(z)|> (1-2e-8)/32-4e-6/32^2 >1/40; injectivity also follows
from F(z1)-F(z2)=(z1-z2)[a+c(z1+z2)].

Exactly, in real coordinates,

    b(x)=delta sigma_N(G(x))+sigma0(G(x))-sigma0(x)
                           -log|F'(G(x))|.        (4)

No parameter derivative of F or G is used. The first term
needs at most five derivatives of the unknown endpoint.
The shift of the analytic seed may use higher SEED jets.

Here is a quantitative two-variable coefficient bound for (4).
Set A_min=.999999. Around any point x in the smaller disk,
the holomorphic inverse Taylor coefficients are

    g_j=(-c)^(j-1) Cat_(j-1)/(a+2cG(x))^(2j-1).

These follow directly by inserting a power series into
F(G(x+z))=x+z. The Catalan generating function is also recorded
in [DLMF 26.5.2](https://dlmf.nist.gov/26.5#E2), inspected
2026-10-20. Only coefficients through degree five are used.
In particular |a+2cG|>A_min, |G(x)-x|<5e-9 and
the real coordinate-column norm of DG-I is below 1e-6.

For T=|delta x_1|+|delta x_2| use the positive formal polynomials

    H(T)=(1+1e-6)T
          +sum_(j=2)^5 [2 Cat_(j-1)(4e-6)^(j-1)
                                      /A_min^(2j-1)] T^j,
    D(T)=5e-9+H(T)-T, d=1/32-1/40=1/160.

The seed mass potential has complex bidisk modulus <.026,
so |sigma0|<.027 there. At a real point in the smaller
disk its two-variable Cauchy series is majorized by
.027/(1-T/d)^2. The divided j-jet bounds for the three
terms on the right side of (4) are therefore respectively

    epsilon [T^j]sum_(k=0)^5(100H)^k /100^j,
    [T^j].027{[1-(T+D)/d]^-2-[1-T/d]^-2}/100^j,
    [T^j]sum_(k=1)^5 2(8e-6 H/A_min)^k/k /100^j. (5)

For the middle expression the constant shift in D is retained
before expansion. Its infinite seed series is convergent,
and the rational expression sums it exactly. Grouping total
degrees in two variables accounts for the factor k+1.
The last expression bounds the real logarithmic Jacobian;
its constant value is irrelevant for j>=2.

For j=2,3,4,5 the sum in (5) is below, respectively,

    .000101086, .000103131, .000107893, .000118290.

All are below 1.3e-4. The actual low norms come from (2),
not from the weaker zeroth/first estimates in (3).
In the seed chart |grad_flat u|<8e-6. Choose a fixed smooth
cap equal to one through 1/80 and zero by 1/40, with slope
at most 160. Its direction gradient is at most

    exp(.027)[8e-6+160(1.2e-8)]<1.1e-5.

Consequently, with P5_40 denoting S154's divided-jet norm
on this smaller chart and cap, uniformly over I,

    P5_40(u)<1.3e-4,
    P5_40(lambda u)<1.3e-4 lambda.               (6)

The path is affine in the fixed conformal chart. The bound
on its direction follows from one finite endpoint transfer;
it does not require differentiating an endpoint normalizer.

## 3. Account for the smaller source chart

Apply [S154](divided-jet-dirac-response.md) with outer cap
radius 1/40 and inner radius 1/80. Its inner Taylor bounds
are unchanged, since the analytic seed still has radius
1/32 and 1/80<1/64. The outer Taylor division now uses 80:

    2.02*80 <400,
    2.02*80^2+300*80 <5e4,
    2.02*80^3+300*80^2+5e4*80 <7e6,
    2.02*80^4+300*80^3+5e4*80^2+7e6*80 <2e9.

The global linear mass-potential bound is still 2r:
exp(.05)[8/(40*.974)+.016+.011*80]<2. The value,
mass and local contact bounds are unchanged. Thus the
complete local Born remainder retains its 1e17 allowance.

Choose the actual-cap gluing cutoff one through 1/160
and zero by 1/80. A smooth cutoff with Euclidean slope
<=310 has seed-metric slope <320. The point-evaluation
disk 1/(10M)<=1/180 fits before its support. Its proper
separation exceeds exp(-.05)/160>1/320.

The same two actual-domain resolvents and point rows
therefore give the scalar comparison bound

    9e11 M exp(-M/320).                           (7)

For example its row/commutator/scalar-factor product is
2(3.2e4)^2(320)(1.02)<9e11. The actual boundary and
HH state are retained; this is not a cap replacement.
Using the full angular sum gives the new complex bound

    K40=2 eta_G[(2/3)1e17+9e11*24*320^3+1e17]
        <1.25e-11.

On P5_40 backgrounds <=.001, Cauchy therefore yields

    |D S[u]|<2.5e-9 P5_40(u),
    |D^2 S[u,w]|<5e-7 P5_40(u)P5_40(w).          (8)

Here S is the normalized Lorentzian scalar trace, with the
same absolute constants as the opposite action source.
The cap change has been charged explicitly; the original
S154 constants are not silently reused.

Since (6) is below .001 throughout the actual path (1),

    |S_lambda|<5e-10,
    |partial_lambda S_lambda|<3.25e-13,
    |partial_lambda^2 S_lambda|<9e-15             (9)

uniformly on I. The seed value is below 4.352e-10 there;
integrating the first bound in (8) keeps the displayed
value allowance. These are actual canonical-source
response bounds on the retained comparison path.

## 4. HH pressure and the middle constraint remainder

The fixed-domain source is continuous and twice
differentiable into the interior supremum space.
[S148](interior-pressure-and-bolt-primitive.md) gives

    P_lambda(c)=F_lambda[S_lambda](c),
    F_lambda[f]=(2W_lambda(c))^-1 integral_0^c W_lambda' f,
    W_lambda=exp(2sigma0+2lambda u)sinh^2 chi.

For c in I the full interval [0,c] remains in I.
The HH integration constant is fixed by bolt regularity.
S148's bounds J<2.4, integral W(t)/W(c)<84 and ||F||<1.2
apply because the actual low geometry bounds are smaller
than its hypotheses. Using the sharper low direction
bounds a=1.2e-8 and g=3e-6 from (2), its exact derivatives give

    ||partial_lambda F||<=2(2.4)a+84g<.000253,
    ||partial_lambda^2 F||<=8(2.4)a^2+8(84)ag<2.5e-11.

Keeping every source/kernel cross term in (9) proves

    |P_lambda|<6e-10,
    |partial_lambda P_lambda|<5.2e-13,
    |partial_lambda^2 P_lambda|<1.1e-14.          (10)

The scalar and pressure Taylor remainders about lambda=0
are thus at most 4.5e-15 lambda^2 and 5.5e-15 lambda^2.
These are unweighted interior bounds, not cold lapse bounds.

For the full canonical normalized profiles rho=P-S,
n=2P-S, the joint value/first/second bounds are

    (N0,N1,N2)=(1.7e-9,1.4e-12,3.1e-14).

The fixed homogeneous canonical contribution remains in
N0. It is not subtracted from the changing geometry factors.
Use [S152](fixed-chart-global-constraint-transfer.md) on its
fixed middle interval 1.001<=t<=40, which lies in I.
For path (1) the geometry direction size is d_g=3e-7.
Set

    m2=1.003 N2+4.02 d_g N1+5.54 d_g^2 N0,
    a2=1.006 N2+10.08 d_g N1+26.5 d_g^2 N0.

The normalization cancels one 2 eta_G factor as in S153.
The middle canonical constraint SECOND responses satisfy

    mass primitive: 279 m2/2 <4.4e-12,
    absolute mass density in BC units: m2/2 <1.6e-14,
    absolute lapse density L1 norm: 1300 a2/2 <2.1e-11. (11)

Their Taylor allowances are half these constants times
lambda^2. The absolute norms apply to density maps; no
scalar absolute-value derivative is taken. The radius and
Jacobian cross terms are included exactly once.

This closes the middle canonical response contribution
on the specified path. It does not close the cold weighted
or outer source contributions, other action sectors, or
the E0 residual. The comparison path remains a comparison
geometry, not a solution certified by its small Hessian.

## 5. Verification protocol and remaining proof

Before acceptance: verify the inverse-quadratic Taylor
coefficients, finite metric identity, two-variable majorants,
chart margins, cap/evaluation compatibility, complete
source constants, changing HH kernel and middle budgets.
Revalidate S154's receipt and its exact current inputs.

Use independent real inversion and coordinate reconstruction
on smooth and C5/non-C6 manufactured metric fields at three
frequencies. Test nonzero point value and both gradient
components; do not restrict normalizers to real quadratic
coefficients. Compare all coordinate derivatives through
five at interior sample points, normalized by j!100^j,
with tolerance 1e-20 at 70 digits. The fractional fixture
contains |x-.005|^(11/2); its fifth derivative is continuous
and its sixth unbounded.

Independently integrate a manufactured HH family and
differentiate its pressure at background parameters .3,.8,
cuts .1,1,4, and real steps .001,.0005. Use panels of
length at most .25 with 32/64 Gauss nodes per panel,
and independent adaptive quadrature for the analytic
second-response formula. Normalize second differences by
1e-14, use tolerance 1e-7(1+|target|), and require refinement
ratio <.35 unless both errors <1e-40. Quadrature/node
agreement must be below 1e-25 in the same normalization.
These tests check transfer and integration algebra, not a
numerically evaluated continuum quantum determinant.

Preserve all runs, failed source snapshots, raw output and
input hashes before/after. Historical source proofs and
receipts remain unchanged.

Remaining to-do list: establish a general admissible family
and cold/outer/other source responses; prove the full
higher-order inverse and coupled exterior/state, charge,
return and signed-source solve; complete opening, stability,
transport, preparation and separate creation.
