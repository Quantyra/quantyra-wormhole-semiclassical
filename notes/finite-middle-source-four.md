# Four spatial derivatives of the finite middle massive source

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, paired domain, HH state, complete
massive angular tower, common subtraction and fixed physical parameters.
This proves a finite numerical middle-region estimate. It does not yet
evaluate the complete physical shape-source H4 norm.

Let gamma_ell be the dimensionless canonical orbit metric, and let

    S_M,ell=(2G/(B^2 C)) Theta_M,ell,
    eta_G=GC/B^2<8.5e-33.

Theta_M includes all massive modes with their original common fourth
restoration. Let S_vac be its fixed flat-vacuum scalar, not zero. Uniformly
for 0<=ell<=1 and 1.001<=t_seed<=40,

    max_(0<=j<=4) |nabla_gamma_ell^j(S_M,ell-S_vac)| <1e110. (1)

This very conservative bound is for the full canonical orbit tensor.
The physical four-dimensional tensor conversion, cold and outer regions,
LLL, Weyl/classical assembly and quadratic residual allocation remain.
In particular (1) is not a value for F_b4 in the
[physical York source split](finite-higher-york-source.md).

## 1. One fixed cap and the real finite metric

Use the [finite interior theorem](comparison-path-interior-response.md)
and the [seed normal chart](seed-middle-cap-transfer.md). The observation
region has Y0>.0606, so the chart and its actual paired realization fit.
In its fixed seed conformal chart, phi=sigma0+ell u, and

    |sigma0|<.002, |partial sigma0|<.067,
    |u|<1.2e-8, |partial u|<8e-6,
    |partial^alpha u|/(j!100^j)<1.3e-4, 2<=j<=5.

The last bounds are component bounds on the radius-1/40 chart, not an
assertion that the unknown finite field is spatially analytic. The seed
has complex modulus .027 on its radius-1/32 chart. Cauchy at distance
from that boundary gives the ordered-tensor reserves used below.

Choose the C11 step eta23 from the
[ninth-source proof](seed-nine-metric-eleven.md), and set

    chi(x)=eta23(160|x|-3), phi_cap=chi phi.

It is one through 3/160 and zero by 1/40. Its slope is below 640.
Consequently |partial phi_cap|<1.35 and |phi_cap|<.003. This auxiliary
choice leaves the actual source unchanged. Fix this cap and the gluing
cutoff throughout observation differentiation. Do not identify it with
a succession of independently recentered caps.

The full metric derivative ledger gives, for j<=9 and D=z0 nabla_gamma,

    |D^j u| <=4e24 j!100^j.

The large reserve deliberately includes every lower-order entry. On this
real chart t<200: a path of chart length at most exp(.027)/40 leaves
Y>.018. For z<=.005 the exact seed polynomial b obeys b>1/7, so
Y(.005)<sqrt(7)/200<.018. Monotonicity then proves the assertion.
Converting orthonormal to coordinate slots costs at most 2^j. Thus the
coordinate components of nabla_gamma^j u are below
4e24 j!40000^j. Seed connection derivatives obey

    |partial^i Gamma_seed| <=1e4 i!1000^i, i<=8.

For example Cauchy and the conformal connection give the smaller
8(.027)(i+1)!640^(i+1). To convert all covariant derivatives, use

    B(n,0)=4e24 n!40000^n,
    B(n,k+1)=B(n+1,k)
             +n sum_(i=0)^k binom(k,i) Gamma_i B(n,k-i).   (2)

This retains every tensor slot. Induction through n+k<=9 gives
B(n,k)<=4e24(n+k)!40000^n(1e6)^k; the recurrence costs less
than .04+.011 of the right-hand side. In particular all ordinary
derivatives of u through nine are below 1e85.

Here are adequate finite bounds on the cap itself. On a complex
Cartesian bidisk of radius 1e-5 centered in its transition, the local
square-root branch changes 160|x| by less than .005. Hence the step
polynomial has modulus below 2^23(1.005)^23<1e8. Cauchy gives

    |partial^j chi| <=1e8 j!(2e5)^j <1e38, j<=5.

Only the matching finite jets of these polynomial branches are used at
the two joins. The uncapped phi derivatives through five are below
1e16 on the cap support. Thus the product derivatives through five
are below 32(1e38)(1e16). No analytic cutoff is presumed.

## 2. Translation jets and a tame fourth-derivative bound

Define a Banach norm X on real plane fields as the maximum of global
value/gradient and the component jets of orders 2..5, divided by
j!100^j, on |x|<=1/80. For any real unit vector v set

    g_v(s)(x)=phi_cap(x+s v), |s|<=.001.

The inner disk stays inside the region where chi=1. The seed local
coefficient bound .027(2)^j, the finite low bounds and the cap gradient
give ||g_v(s)||_X<2. The preceding real C9 and cap bounds give

    ||g_v^(4)(s)||_X <=M4=1e88.                         (3)

These are genuine translation derivatives of the entire cap, including
its transition. Put h=M4^(-1/4)=1e-22. Taylor's formula at 0 at the four
points s=0,h,2h,3h has remainders at most k^4/24. Inverting the fixed
matrix (k^j)_(k,j=0..3), with right-hand-side bounds 2+k^4/24, gives

    ||g_v^(j)(0)||_X <=1000 h^-j, j=1,2,3;
    ||g_v^(4)(0)||_X <=h^-4.                           (4)

This elementary interpolation is what avoids the fourth power of a
coarse C6 allowance. It is linear in the high C9 allowance at final
order four. It does not differentiate a source-value inequality.

Introduce ONLY the finite parameter polynomial

    p_v(z)=g_v(0)+sum_(j=1)^4 z^j g_v^(j)(0)/j!.

For |z|<=rho=1e-28, its X-distance from the real base is below
1000 sum_(j=1)^4 (1e-6)^j/j! <.002. The true translated field and
this polynomial have the same four parameter jets at zero. The true
field is not evaluated at any complex spatial point.

The [divided Born proof](divided-jet-dirac-response.md) applies directly
to this plane family. Its proof requires the low global h bounds and
the inner finite jets, not a factorization of each perturbation by a
fixed cutoff. The background plus the polynomial change has divided
perturbation size below .011. At the center the seed jets retain their
old values. Globally |exp(phi_cap-phi_cap(0))-1|<.06 and its linear
bound is 2|x|, also on this complex parameter disk. The inner Taylor
argument works through 1/80, and the four outer-division inequalities
in S155 retain the same 1e17 per-mode M^-3 reserve.

The common fourth restoration has the unchanged 1e17 C reserve. The
scalar flat-vacuum value is separated before taking a bound, not set
to zero. Therefore the complete nonconstant capped source obeys

    Kcap=2 eta_G[(2/3)1e17+1e17] <3e-15.             (5)

This retains the full angular sum, local subtraction and common fourth
coefficient. Uniform Born convergence makes this an X-holomorphic
functional near the real base, so composition respects the four jets.

Translation covariance on the plane identifies its derivative along
g_v with the observation derivative of the SAME fixed-cap scalar.
The underlying Dirac naturality is described by
[Dabrowski--Dossena](https://arxiv.org/abs/1209.2021), inspected
2026-10-22; here only ordinary plane translations are needed and the
identity follows directly by conjugating the kernel. Common proper
point splitting and the local restoration have the same covariance.
Cauchy on p_v gives j! Kcap rho^-j. Real polarization and the two
coordinate slots convert this to a full Cartesian tensor bound

    |partial^j(S_cap-S_vac)|
       <=2^(ceil(j/2)) j^j Kcap rho^-j <1e101,
       1<=j<=4.                                      (6)

At j=0 use (5). Every cap derivative is included in this calculation.

## 3. Differentiate the actual-domain gluing term

Use the fixed S155 cutoff, one through 1/160 and zero by 1/80, with
commutator norm below 320. The actual and capped coefficients coincide
on its support. The actual resolvent and the cap resolvent retain the
point rows and inverse bounds from that theorem.

Let R=1/4096. On |x|<=R the proper separation is greater than
.95(1/160-1/4096)>1/176 and also greater than .1/18. Thus the original
evaluation disks fit for every M>=18. After the harmless endpoint
identification the initial flat row is bounded by
33000 exp(-M/352), throughout this ball.

For a=e^phi the homogeneous local equation is (D_flat+Ma)v=0.
Conformal covariance is [Walpuski, Proposition 15.5](https://walpu.ski/Teaching/SpinGeometry.pdf),
inspected 2026-10-22. Direct integration by parts in the skew-adjoint
flat convention gives

    ||(D_flat+Ma)w||^2
      =||partial w||^2+M^2||aw||^2-M<w,c(da)w>.

Here a>=exp(-.003), |da|<.07, so the right side is at least
.98(||partial w||^2+M^2||w||^2). This is a compact interior
estimate and imposes no new global boundary condition.

On this small ball the seed Cauchy and finite divided jets give
the raw derivative bounds phi_1<=.07 and

    phi_j <=j![.027(128)^j+.00013(150)^j], 2<=j<=5.

Finite exponential coefficients bound a and its endpoint factor
exp(-phi/2). For a use A1=.07 and A_j=.1 j!300^j for j>=2;
for the endpoint use zeroth factor 1.01 and .1 j!300^j above zero.
All these finite comparisons are checked coefficientwise.

For target j<=4 put m=j+3, L=8192m and B=10000m. With nested
balls as in the [coercive row proof](seed-middle-coercive-jets.md),

    E0<=1.1 L F0,
    E_k<=1.1[L E_(k-1)+sum_(b=1)^k binom(k,b) A_b E_(k-b)].

The coefficient ratio is below one through k=j+1. Hence
F_k<=B^k F0 through j+2. The same two-dimensional Fourier point
estimate, with rho_point=R/(2m), bounds the physical endpoint rows by

    c0=34000,
    c_j=600000(j+3)[10000(j+3)]^j, 1<=j<=4,
    ||row_j||<=c_j exp(-M/352).                        (7)

Both rows may be Hilbert-valued. Keep the rank-two spin trace. The
source factor M exp(-phi) is allocated as M with an exp(-phi/2)
endpoint on each row. The single fixed commutator costs 320. Thus
all derivatives of the actual-minus-cap source have bounds

    G_j=2 eta_G *2*320*24*176^3
                     sum_(k=0)^j binom(j,k)c_k c_(j-k) <1e20. (8)

Indeed N(E)<=4CE^2 gives sum d_n M_n exp(-M_n/176)<=24C176^3.
All angular modes, both observation rows and the original actual
paired domain remain. These estimates justify termwise observation
differentiation; no derivative of an unsigned value bound is used.

## 4. Actual canonical tensor conversion and consequence

Combine (6) and (8) before changing derivatives. The coordinate
derivatives through four are below 2e101. The conformal connection
and its first two derivatives have reserves

    Gamma0=.3, Gamma1=4000, Gamma2=2e6.

They follow from the raw phi jets above; Gamma is linear in its first
derivatives. The tensor recurrence, starting C(0,k)=1, is

    C(r+1,k)=C(r,k+1)
          +r sum_(i=0)^k binom(k,i) Gamma_i C(r,k-i), r+k<=3.

All C(j,0), j<=4, are below 1e8. The actual center metric changes
the norm by less than two. Therefore 4e101 max C(j,0)<1e110,
proving (1). This uses the actual canonical metric, not the seed
connection as if it were unchanged.

The bound supplies the previously missing finite middle massive
spatial estimate. It is deliberately coarse and asserts no useful
Schur-target pass. Next carry the finite-jet method to the cold and
outer patches, where endpoint weights and the correlated cubic must
be retained, and assemble the complete physical shape residual.
The fixed vacuum scalar must participate in that assembly.

## 5. Verification contract

Freeze this note, checker and review before execution. Check the
ordinary-to-covariant conversion, finite cap bounds, Vandermonde
interpolation, complete Cauchy/row/angular budgets and final tensor
reserve by exact rational arithmetic. Independently verify the Dirac
energy identity with two Gaussian spinors and nonconstant positive
mass potentials. Check the fixed-cap translation identity by direct
kernel derivatives and integration by parts at 60/90 digits and
48/96 Gauss nodes; require errors below 1e-25. Dropping cap derivatives
must fail. Test the finite-translation polynomial on a real C9 field
which is not C10, preserving the fourth source jet without assigning
it a holomorphic spatial extension. These finite diagnostics support
the ordinary argument; they do not evaluate the charged determinant.
Preserve all runs, source snapshots and inherited files.

Remaining to-do list: cold and outer finite massive spatial bounds;
complete LLL/Weyl/classical shape assembly and physical H4 norm;
higher TT inverse and positive transfers; nonlinear mapping,
preservation, existence and signed compatibility; physical matching,
stability, charge/return, opening, transport, preparation and creation.
