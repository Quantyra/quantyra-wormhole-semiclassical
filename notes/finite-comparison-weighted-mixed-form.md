# Weighted mixed form on the actual finite comparison

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof before
frozen verification. **The actual mixed energy-dual estimate is proved.
The physical TT operator realization and coupled range remain open.**

Use the full physical action at each original finite comparison g_ell,
0<=ell<=1. In its physical volume, put t=r0/B and omega=t^(-1/2).
The complete local trace, including its anomaly, obeys

    |Etr|<3.2e-12 lambda/t.                              (1)

For a genuine spherical TT lift h=H_ell c and its outer-fixed return
Y=eta partial_rho, the exact mixed row has a weak decomposition

    sqrt(80pi^2) Bmix c=div f[c]+g[c],
    ||omega^-1 f[c]||<4e6 sqrt(lambda) sqrt(.90 E_ell[c]),
    ||omega^-1 g[c]||<1000 sqrt(lambda/G) sqrt(.90 E_ell[c]). (2)

The estimates also accommodate the local scalar auxiliary channels as
specified below. They retain the off-shell trace contact and full TT
return; no derivative of a nonlocal stress functional is needed.

If the intended TT operator D_ell is subsequently realized as a positive
selfadjoint operator with spectral form >=.90 E_ell, these estimates
give an outgoing auxiliary contraction below .025. This conditional
operator consequence must not be called a constructed coupled inverse.

## 1. Inputs and value bounds

Keep the prepared temperature, physical flux, clock, paired domain and
finite coefficients. The [preceding curvature theorem](finite-comparison-scalar-curvature-remainder.md)
gives Z=R4+12lambda with |Z|<3.1e-12 lambda/t. We use the actual
[C5 metric contract](stationary-completed-metric-response-norm.md), its
fixed endpoint charts and outer proper chords, not a sixth derivative
or differentiability of the moving proper chart in ell.

First sharpen the curvature VALUE weights with the already available
low jets. On t<=40, the sectional differences from the seed have
maximum below 1.5e-5 lambda. The corresponding Frobenius differences
of Rm, J and C4 are below .001lambda. The seed has
|J0|<=7lambda/t^4 and |C0|<=40lambda/t^3.

On t>=40 the [scalar comparison](finite-comparison-scalar-resolvent.md)
gives |w|,|S|,|w_y|,|S_y|,|w_yy|,|S_yy|<10/t^3.
Here y is seed proper distance in lambda units. The seed polynomial
gives A0<2 there, R0<2.24, |K0|<10.5, |E20|<22, A0R0<=11.
The exact pullbacks in the classical disk proof therefore give

    |Delta K|/lambda<250/t^3,
    |Delta E2|/lambda<500/t^3,
    |Delta(AR)|/lambda<300/t^3,
    |Delta(Rr^2)|/lambda<150/t^3.

The radial-sphere sectional change is below 950/t^3 and the
sphere-sphere one below 151/t^3. Thus |Delta Rm|<6000lambda/t^3
and |Delta J|<12000lambda/t^3. Weyl projection is orthogonal in the
pointwise curvature norm, so |Delta C4|<=|Delta Rm|. Comparing the
orthonormal frames is an isometry, not a change of physical metric.
Combining the two regions gives the convenient global bounds

    |J|<8lambda/t, |C4|<50lambda/t,
    |Rm|<55lambda, |Ric|<15lambda,
    |P_schouten|<6lambda, |R4|<13lambda.              (3)

For J in the outer region use 12000/40^2+7/40^3<8. On the inner
region use 7/t^3+.001t<=7.04. The unweighted Rm, Ric and Schouten
bounds use the global sharp low sectional differences, not the coarser
outer estimates. For P=J/2+R4g/24 use the orthogonal decomposition.

## 2. Four metric derivatives suffice for the remaining coefficients

Here are explicit deliberately loose bounds, in the actual covariant
tensor norms:

    |nabla J|, |nabla C4|<1e15 lambda^(3/2) omega,
    |Hess R4|<1e20 lambda^2/t.                       (4)

We give the local-coordinate count, including the bolt and outer end,
so that (4) does not follow from an unquantified C5 assertion.

In a point-normalized coordinate system with g_ij=delta_ij at its
center, let Mj bound each component of the jth metric derivative.
Set

    P3=M3+M1 M2+M1^3,
    P4=M4+M1 M3+M2^2+M1^2 M2+M1^4.

Direct differentiation of the inverse metric and Christoffel formula
gives, for each component (coordinate derivatives in these lines),

    N1=M1, N2=M2+8M1^2,
    N3=M3+24M1 M2+96M1^3,
    G0=3M1/2,
    G1=3M2/2+6M1^2,
    G2=3M3/2+18M1 M2+48M1^3,
    G3=3M4/2+24M1 M3+18M2^2+288M1^2 M2+576M1^4,
    Q0=2G1+8G0^2,
    Q1=2G2+16G0 G1,
    Q2=2G3+16G0 G2+16G1^2.

N bounds inverse-metric derivatives, G bounds connection derivatives,
and Q bounds derivatives of R^a_bcd. Lowering the first curvature
index and adding every connection term gives

    |nabla Rm|<=32[Q1+(4M1+16G0)Q0],
    |Hess R4|<=4[16Q2+128M1 Q1+64N2 Q0
                                   +4G0(16Q1+64M1 Q0)].

The last polynomial has coefficients 192,8064,5376,165888,506880
in the order defining P4. Ricci contraction costs at most two and
Weyl/tracefree projections at most one. Consequently

    |nabla J|,|nabla C4|<=1e6 P3,
    |Hess R4|<=1e6 P4.                              (5)

These are component counts in four dimensions, not elliptic estimates.

### The bolt and compact cover

Take the fixed seed orbit normal charts in canonical A units,
and sphere normal coordinates in radius sqrt(C) units. At a point
rescale all four coordinate distances by t. The physical derivative
conversion is less than 2.5 sqrt(lambda), since C/x<6.1 and the
finite radius/metric value change is below 3e-8. Constant frame
normalization factors are reserved below.

For completeness, the necessary seed coordinate jets are elementary.
With f=F0/C and D=t partial_t, the canonical orbit derivatives of
phi0=log t have p=sqrt(f), Hessian eigenvalues a=Df/2 and b=a-f.
We have f<t^2, |a|<2t^2, |b|<3t^2,
|Da|<2t^2, |Db|<5t^2, |D^2a|<4t^2, |D^2b|<8t^2.
The third covariant derivative components are p Da, p Db, pb, pb.
Their radial derivatives and angular rotations give bounds

    |nabla^j phi0|< (1,5,13,100)_j t^j, j=1,...,4.

The apparent polar connection is harmless: A_hat p=b. The canonical
orbit curvature and its first two derivatives obey bounds 2,40t,
252t^2, obtained from Rhat=-D^2f+3Df-2f. Converting scalar
covariant jets to normal-coordinate partial jets permits the reserves
(1,5,40,6000)t^j. The normal metric jets through four are bounded
by (0,10,1e4,1e7) after the point rescaling. These bounds also
cover the constant-curvature sphere. No polar chart is used at the bolt.

The finite metric contract bounds the two local log changes through
four by 1e-4 j!100^j; the cold row is smaller. Enlarge to
.01 j!100^j to include the fixed normal chart and constant frame
conversion. In the point-rescaled coordinates this bounds the log
multipliers 2phi and 2phi+2u through four by

    L=(10,1000,3e5,1e8).

Their exponential Bell polynomials, multiplied by the seed normal
metric jets, are strictly below the following common component ledger:

    (M1,M2,M3,M4)=(30,3000,1e6,4e8).                (6)

For example the fourth exponential polynomial is
L4+4L1L3+3L2^2+6L1^2L2+L1^4=115610000. Adding the
normal-metric product terms and constant frame factor still leaves
less than 4e8. The core t<=40 is covered by the cold and middle
charts: when needed, the outer guard Y>1/(.412t) puts this region
inside Y>=.05. The argument evaluates each endpoint chart once.

In (5), P3=1117000 and P4=442510000. After the physical conversion,

    sqrt(40) 2.5^3 1e6 P3<1e15,
    40 2.5^4 1e6 P4<1e20.                          (7)

This proves (4) on the compact cover with its stated weights.

### The outer end

Use the fixed endpoint proper coordinate Y, with 1/t<Y<3/t up to
a factor 1.001 reserved here. The proper-map estimate in the metric
proof compares its seed t with the original t by the same factor.
The endpoint log changes v_ell,w_ell have cubic coefficients below
1e-9 and remainders bounded through four by

    1e-7 j!100^j Y^(4-j)(1+log(1/Y)).

They have the retained correlated cubic, although the estimate below
does not need its cancellation. Since (1+log t)/t<.12 for t>=40,
their point-rescaled partial jets are bounded by
1.2e-6 j!100^j/t^3. Work in radial proper coordinate and local time
arc, with sphere normal coordinates; this chart is away from the bolt.
The seed log warping and radius pole have normalized jet sums below
two through four: use the analytic calB guards, and differentiate
-log Y after the point rescaling. Leibniz, exponential and constant
frame products through four have coefficient sum below 1000.
Thus the metric-jet differences from the seed satisfy, with M0=1,

    |Delta partial^j g_ab|<.001 Mj/t^2, j=0,...,4.   (8)

Indeed .0012 j!100^j/t^3<.001 Mj/t^2 for every j<=4 and
t>=40. This uses the actual conormal remainder, not a bounded fourth
derivative at Y=0.

The explicit curvature/connection formulas behind (5) are finite
products with at most sixteen factors, including inverse-metric and
constant frame factors. Along the short metric segment their inverse
norms are below two. The telescoping product estimate is bounded by
16*2^15<1e6 times the metric-jet change. Together with (8), this
allows the conservative outer derivative bounds

    |Delta nabla J|,|Delta nabla C4|
                           <=1000(2.5^3)(1e6 P3)lambda^(3/2)/t^2,
    |Delta Hess R4|<=1000(2.5^4)(1e6 P4)lambda^2/t^2.

The comparisons include connection and frame changes, not just radial
derivatives of curvature components. The seed bounds 1000/t^3 and
4/t^5 from the outgoing proof may be added. Multiplying by sqrt(t)
or t and maximizing these decreasing powers at t=40 proves (4).
Only four finite metric derivatives have been used.

## 3. Full trace and exact weak mixed row

The full physical trace is local; its artificial reference cancels:

    Etr=(1+11epsilon lambda/3)Z+epsilon Box Z
       -11epsilon Z^2/72+(7epsilon/12)C4^2
                           +(11epsilon/6)J^2-20epsilon M_em,
    epsilon=G/(60pi), M_em=q^2/(2r_ell^4).

Use |Box R4|<=2|Hess R4|, (3)-(4),
.999*2000^-10<G lambda<2000^-10, and
Gq^2/(lambda B^4)<4e-27. The radius change contributes less than
1.001 to r_ell^-4/r0^-4. The coefficient of lambda/t is at most

    (1+11g/(180pi))(3.1e-12)
     +(g/(60pi))[2e20+11(3.1e-12)^2/72+7(50)^2/12+11(8)^2/6]
     +10.01(4e-27)/(60pi)<3.2e-12,
    g=G lambda.

This proves (1), including the previously unbounded anomaly derivative.
It is a value bound on the full trace; no derivative of Etr is taken.

Use the exact TT trace variation in
[the seed inverse proof](stationary-weighted-tt-inverse.md), which did
not require a background field equation. Its integration by parts is

    Ttrace h-div(Etr Y)=-div g_h+f_h,
    g_h=epsilon grad(J:h)-(7epsilon/3)C^{abcd}K_adb
                    +(11epsilon/6)J^{ab}grad h_ab+Etr Y,
    K^a_bc=(nabla_b h_c^a+nabla_c h_b^a-nabla^a h_bc)/2.

The free Weyl-vector index is c. The algebraic Weyl terms are bounded
by [2|C||Rm|+8|C||P|+4|C|^2]|h|. The algebraic Lichnerowicz,
R4 and J terms are bounded in the same actual frame. Equations
(3)-(4) and |K|<=3|nabla h|/2 consequently give

    |g_h|<=epsilon[200lambda omega |nabla h|
                         +1e15 lambda^(3/2)omega |h|]
                                      +3.2e-12 lambda omega |Y|,
    |f_h|<=6e15 epsilon lambda^(3/2)omega |nabla h|
          +[8lambda omega+2e20 epsilon lambda^2 omega
                         +20.02epsilon q^2 B^-4 omega]|h|. (9)

For the derivative budgets use (17/6)8+(7/2)50<200 and
(7/2+11/6)1e15<6e15. The 2e20 allowance includes the 1e20
Hessian of R4 and every curvature algebra term, which totals less
than 1e6lambda^2 omega. All t^-1 weights may be weakened to omega.
The off-shell divergence is retained, including its physical transpose.

## 4. Actual TT energy bounds and quantitative consequence

Put e=8piG E_ell[c]. The [finite return proof](finite-comparison-classical-disk-contact.md)
gives ||h||<4.8sqrt(e/lambda), ||Y||<4.83sqrt(e)/lambda,
||h'_rad||<75sqrt(e), ||A(b+c)||<9.6sqrt(e), and
||b-c||<3.002sqrt(e/lambda). The full covariant norm is

    ||nabla h||^2=||h'_rad||^2+8||A(b+c)||^2+4||Rr(b-c)||^2,
    ||nabla h||<82sqrt(e), Rr<2.25sqrt(lambda).      (10)

In particular the angular connection terms at the bolt have not been
dropped. Converting sqrt(8pi) and using pi<22/7 gives the convenient
coefficients 25,25,420 for h,Y,nabla h in units sqrt(G E_ell).

Since sqrt(80pi^2)/(8piG)=sqrt(5)/(2G), set
f=sqrt(5)g_h/(2G), g=-sqrt(5)f_h/(2G). Divide by sqrt(.90 E_ell).
The factor sqrt(5)/(2sqrt(.90)) is below 1.2. Equations (9)-(10),
3e-17<sqrt(g)<4e-17 and pi>3 give

    A_f<=1.2[(sqrt(g)/180)(200*420+1e15*25)
                                    +25(3.2e-12)/sqrt(g)]<3.3e6,
    A_g<=1.2[8*25+(g/180)(6e15*420+2e20*25)
                                +(20.02*25/180)(4e-27)]<241.

These prove (2) with slack. The local scalar auxiliary channels
U_s from the preceding curvature theorem satisfy
||omega^-1 f_s||<4sqrt(lambda) and
||omega^-1 g_s||<2e-5 sqrt(lambda/G), using |a_s|<10lambda/t
and |Z|<3.1e-12lambda/t. Direct-sum Cauchy therefore preserves
the rounded coefficients 4e6 and 1000 for the combined channels.

The estimate is first proved on the regular compact TT core, then
extends by energy continuity, and in particular to the preceding
stronger form completion. This does not identify that completion with
the physical logarithmic operator domain or prove mass/fast membership.

## 5. Conditional outgoing bound with the actual scalar reference

The [finite scalar theorem](finite-comparison-scalar-resolvent.md) gives
||q_ell||_1<72000, kbar>4e17 and |D_rho-D_Bessel|<5sqrt(lambda).
Differentiating its same Bessel resolvent series as in the seed proof
therefore gives, in physical Liouville variables,

    |G|<1.001/k, |D_rho G|<3,
    D_rho R+ D_rho*=I+K, |K|<3k.

For the last line retain the identity distribution. Before rounding
the kernel bound is 2k+[4(72000sqrt(lambda))+20sqrt(lambda)
+25lambda/k]/(1-2e-13). There is no new bolt boundary term. The
norm limits follow by the same weighted dominated-convergence argument.
The actual weight has ||omega||_infinity<=1 and
integral omega^2 d rho_ell<200/sqrt(lambda): use z_40<178 and
dz<=1.5dt/t on the outer tail, whose contribution is below 1.5/40.

Now ASSUME that the intended D_ell has its positive selfadjoint
realization with spectral form >=.90 E_ell. Then c=D_ell^-1/2 y
is a legitimate use of (2); no such square root is assumed earlier.
Combine its mixed channel with U_s. The actual scalar reference and
the coefficients above bound the outgoing auxiliary Lcal by

    ||Lcal|| < 201(1000)^2/(180 kbar)
       +(1200/180)(1000)(4e6)sqrt(g)
       +601(4e6)^2/kbar <.02405.

The previously proved ||A22||<2.6e-13 gives the positive
low-frequency elimination norm <1.000002. Hence

    ||(I+Acal Jcal)^-1 Lcal Jcal||<.025<1.           (11)

The positive elliptic subtraction in the reference fourth-order
resolvent is retained; Q^-1 is not replaced by its high-energy value
at low frequencies. This supplies a usable uniform numerical condition
for the pending operator construction, not that construction itself.

## Verification and remaining scope

Freeze this note and checker before execution. Pin the preceding
curvature, scalar and quantum-disk receipts and reconcile their input
union. Check the metric/curvature derivative component count, normal
exponential jets, outer weights, physical conversions, full trace,
weak mixed budgets, covariant TT gradient and conditional contraction
with exact rational arithmetic. Check the physical divergence pairing
and signs on independent smooth compact radial fixtures, retaining
the density/gauge divergence. These finite checks support ordinary
continuum arguments; they are not formal verification or a physical
numerical solution. Preserve all runs and failed source snapshots.

Next realize the actual TT operator and York/return maps in the full
logarithmic domains, prove physical mass/fast membership and coupled
range/adjoint with the actual phase. The mixed contraction estimate
can then be applied to that operator. Pi_cut stays fixed at the seed;
all S167 acceptance tasks and physical completion requirements remain
open. Restricted Euclidean estimates do not imply dynamical stability.

Remaining to-do list: TT/logarithmic operator realization and candidate
coupled range; cold/higher response and nonlinear preservation;
projected solve and signed compatibility; exterior/state, charge/return,
opening, stability, transport, preparation, separate creation and audit.
