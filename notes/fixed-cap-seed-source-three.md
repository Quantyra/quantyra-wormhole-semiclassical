# A fixed-cap middle source bound for the product-centered seed

2026-10-22, E003/S169. Ordinary argument before frozen verification.
The full stable semiclassical wormhole objective remains incomplete.

On the ACTUAL product-centered seed, in the middle cover s>=.001,
Y>=.05, the complete canonical scalar trace obeys

    max_(0<=j<=3) |nabla_gamma^j Theta|/(j!100^j q^2)<6e15. (1)

Here gamma is the canonical orbit metric in A units, Theta=B^2 tau_hat,
and Y is its inward proper distance from the actual outer seam. This
includes the massive tower, common restoration and LLL. The independent
four-dimensional Weyl/radius source and complete auxiliary metric
equations are NOT replaced by this canonical scalar. This is a genuine
source derivative improvement in the original units, not a metric solve.

The proof replaces a large point-evaluation row constant by a direct
massive-resolvent estimate and differentiates the UV term with ONE fixed
cap. The earlier complete 1e22 claim was withdrawn in
[the splitting audit](seed-middle-splitting-audit.md); its fixed-cap
gluing estimates survived. The repaired complete allowance was 1e26.
We do not use or reinstate the withdrawn combination.

## 1. Actual seed and the fixed cap

Keep the product, positive temperature, flux, paired domain and fixed
mu0 of [the product-centered construction](recentered-product-and-cold-residual.md).
The source is evaluated on its original metric argument, not g_E.
The actual normalized canonical polynomial is

    Bcal(z)=x/C+z^2/C-(2+4x+delta)z^3/C
             +(1+3x+d)z^4/C-2d z^5/C+d z^6/C,
    C=1+6x+d, x>8e9, |d|<8.6e-19, delta/C<1.001h^2.

The retained strict guards hold here:
1/C<1e-10, x/C<1/5.9, |d|/C<1e-24, C/q<.00352,
q/C<400, M_n^2=n(n+q)/C, M_1>18, d_n=2(q+2n).
Also |log(mu0 B)|<1e-8 and |4gamma_E-2H_q+4log(mu0 B)|<100.
The latter bound, not the obsolete 1e-21 scale displacement, is used
in the common restoration. No finite coefficient is changed.

The middle normal-chart proof in
[the cap transfer](seed-middle-cap-transfer.md) uses precisely these
polynomial guards. On the real interval |Bcal''|<2.001; its complex
extension of radius .02 has |Bcal'''|<17, so |Bcal''|<8 there.
The coefficient sums give |Bcal'|<7.08 and |Bcal''''|<13. The same
Goursat iteration, with row sums <.02, gives on radius R*=1/30

    |sigma|<.002, |partial sigma|<.067, |partial^2 sigma|<3,
    |z_z|,|z_w|<.284, |sigma_z|,|sigma_w|<.034.

Use the SAME fixed degree-nine cap as that proof:
chi(r)=eta(128r-3), one through 3/128 and zero by 1/32,
eta=sum_(i=0)^4 binom(9,i)t^i(1-t)^(9-i) on [0,1].
It is C4. Put phi=chi sigma on the plane. Its global value is <.002;
the uncapped derivative norms of chi through three are <=315,7e5,3e9.
The physical gluing cutoffs and both commutators remain fixed; their
norms are <=256 and their supports are separated from the center by
at least 1/65 in proper units. No actual endpoint is replaced by the cap.

Two additional elementary seed-jet guards will be useful. Differentiating
sigma_zw=Bcal'' exp(2sigma)/8 and z_zw=Bcal' exp(2sigma)/4 gives

    |z_zw|<1.78, |z_zz|<.03, |sigma_zw|<1.005,
    |sigma_zzw|<.68, |sigma_zz|<.023,
    |sigma_zzzw|<.34, |sigma_zzz|<.012,
    |sigma_zzww|<6.1.

For example the numerator for sigma_zzzw is bounded by
13(.284)^2+17(.03)+4*17(.284)(.034)+16(.023)+32(.034)^2,
times 1.005/8. For sigma_zzww replace 17(.03) by 17(1.78)
and 16(.023) by 16(1.005). Pure fourth derivatives on radius 1/32
are <480(.012)<6.1 by one-variable Cauchy from R*. Symmetry supplies
the remaining mixed derivatives. Cartesian conversion costs at most
2^j per component and 2^(j/2) for the ordered tensor. Consequently

    |partial^3 sigma|<100, |partial^4 sigma|<1000        (2)

on the real cap support. These are finite analytic seed estimates;
no analyticity is assigned to a finite corrected metric or to the cutoff.

## 2. Replace the crude cap row by a direct resolvent row

In the flat cap representation let R_phi=(D_flat+M exp(phi))^-1.
At an observation x with |x|<=1/256 put m=M exp(phi(x)),
h_x(y)=exp(phi(y)-phi(x))-1 and r=|y-x|. Then

    m>17.9, |h_x|<.0041,
    |h_x(y)|<=.07 r+100r^2.                            (3)

For r<=1/128, both points are in the uncapped chart and Taylor's
formula with the gradient/Hessian guards proves (3). For r>=1/128,
100r^2>.0041 suffices. Set theta=M/2. Conjugate by exp(theta r), first
with truncated Lipschitz r and then pass to the limit. The conjugated
flat operator differs from D_flat+M by a bounded multiplication of norm
at most M(exp(.002)-1)+theta<M. Its Neumann inverse exists and

    ||exp(-theta r) R_phi exp(theta r)||
       <=1/[M(exp(-.002)-1/2)],
    m/[M(exp(-.002)-1/2)]<2.1.                        (4)

Accretivity supplies the displayed sharper inverse bound; the adjoint
has the same coercivity, so the range is all L2. This is a PLANE cap
estimate. The actual paired-domain inverse in the middle of the gluing
identity remains the retained <=1/M_n bound.

For S_m=(D_flat+m)^-1, its rank-two kernel has squared Hilbert--Schmidt
norm m^2[K0(mr)^2+K1(mr)^2]/(2pi^2). The positive integral formulas
[DLMF 10.32.8-9](https://dlmf.nist.gov/10.32) give, for t>0,

    K0(t)<=sqrt(pi/(2t)) exp(-t),
    K1(t)<=sqrt(pi/(2t)) exp(-t)[1+3/(8t)].

The first follows from cosh v>=1+v^2/2. For the second set s=1+v
in the nu=1 integral and use sqrt(1+v/2)<=1+v/4. Since
2-M/m>.99, elementary gamma integrals give

    ||r S_m(x,.) exp(theta r)||_HS <1.6/m,
    ||r^2 S_m(x,.) exp(theta r)||_HS <5.3/m^2.

Indeed the squared constants are bounded by

    A_j^2=(1/2)[2(2j)!/.99^(2j+1)
             +(3/4)(2j-1)!/.99^(2j)
             +(9/64)(2j-2)!/.99^(2j-1)], j=1,2.      (5)

On a remote set with r>=D and mD>=.19 the weighted free row is <3;
its squared bound is at most
[2+3/(4*.19)+9/(64*.19^2)]/(2*.99)<9.

Use the exact identity R_phi=S_m-S_m m h_x R_phi, restricting the
free row to that remote set before taking its L2 norm. Equations
(3)-(5) yield a row constant at most

    3+2.1[(.07/17.9)*1.6+(100/17.9^2)*5.3]<7.

The endpoint conformal factors fit within 1.01; reserve 10. On the
observation ball the remaining separation is >1/88, hence

    ||row_0||<=10 exp(-M/176).                        (6)

This controls an actual cap row, not an assumed pointwise bound on the
actual global Green kernel. Both source and adjoint rows satisfy it.
The rank-two trace will still be retained in the gluing budget, even
though the Hilbert--Schmidt argument above already overestimates a
single operator row. The apparent free diagonal singularity is never
integrated unweighted: its remote part is separated and h_x vanishes
linearly at x in the remaining term.

## 3. Fixed observation derivatives and full angular gluing

Reuse the local coercive recurrence in
[the retained row proof](seed-middle-coercive-jets.md), sections 2-4;
these fixed-cap results were NOT withdrawn. The local positive mass
estimate, nested balls R=1/256, L=512(j+3), B_j=600(j+3), and endpoint
jets are unchanged. Only its initial uniform row is replaced by (6).
The Fourier point estimate and endpoint product therefore give

    c0=10,
    c_j=180(j+3)[600(j+3)]^j, 1<=j<=3,
    ||row_j||<=c_j exp(-M/176).                       (7)

The coefficient 180 is larger than
10*1.01*[9+4(75/64)+2(75/64)^2]. No derivative of a varying cap or
of the distant actual inverse is hidden in (7).

The unchanged angular argument, d_n=2(q+2n), M_n^2=n(n+q)/C, gives
sum d_n exp(-M_n/88)<=4C88^2. Thus the actual gluing contribution has

    G_j=8*256^2*88^2*(31/2500000)
           sum_(k=0)^j binom(j,k)c_k c_(j-k)/(j!100^j) (8)

as an upper bound for |partial^j G_fix|/(j!100^j q^2). This retains
both commutators, the actual paired inverse, spin multiplicity and the
full angular sum. The constant 31/2500000 bounds (C/q)^2. It is not
a truncated-mode calculation or an infinite-temperature approximation.

## 4. Differentiate the SAME UV cap, with no moving-cap substitution

Let U_fix(x) be the renormalized canonical trace of this cap, including
the common restoration. Exactly Theta(x)=U_fix(x)+G_fix(x). We will not
replace U_fix(x) by a cap recentered independently at x.

For any real unit direction e form g_e(s)(y)=phi(y+s e), and the finite
parameter polynomial p_e(z)=sum_(j=0)^3 z^j g_e^(j)(0)/j!.
Use |z|<=rho=1/500. This polynomial matches the first three jets of
the REAL translation; no complex spatial evaluation of chi is made.

From (2) and the fixed cutoff bounds the global translation tensors
of phi obey

    T1=.697, T2=1445.21, T3=6143635.

Hence |p_e|<=.002+rho T1+rho^2 T2/2+rho^3 T3/6<.0145.
At y=0 its point value is <7e-6, using sigma(0)=grad sigma(0)=0.
For |y|<=1/80 every required translation jet is uncapped. There the
coefficient perturbation has gradient at most
3rho+100rho^2/2+1000rho^3/6<.0063.

The Cartesian complex bidisk radius 1/96 fits inside R*=1/30 around
this inner real disk. For k=2,...,5, Cauchy's finite jet formula gives

    |partial^k(p_e-phi)|/(k!100^k)
      <=.002(.96)^k sum_(j=1)^3 binom(k+j,j)(96sqrt(2)rho)^j
      <.011.                                         (9)

At the center the matched complex mass has Re m>.98M and |m|<1.02M.
Globally |exp(p_e-p_e(0))-1|<.015. For |y|>=1/80 this is <2|y|;
inside, the local seed Taylor estimate and the <.0063 gradient change
give the same linear bound. Thus all global mass-potential bounds and
inner fifth-jet bounds in [the divided Born proof](divided-jet-dirac-response.md)
hold for this PLANE family. The 1/80 inner disk uses the same explicit
outer-division reserves as section 2 of
[the fixed-cap finite-source proof](finite-middle-source-four.md).

The cap's global GRADIENT need not be .011 in this argument. This is
not an application of the entire actual-domain Frechet theorem to a
larger ball. In the plane Born proof the vertices are multiplication
by h, its norm <.06 sums the Born series, its linear/global Taylor
remainders and the inner five jets give the M^-3 remainder. No remote
gradient occurs in those loop estimates. The actual-domain gluing is
estimated separately by (8) at the fixed real cap. This is the precise
scope in which the larger translation disk is justified.

The unchanged per-mode remainder 1e17 M^-3, sum d_n/M_n^3<=2C/3,
and common restoration reserve 1e17 C consequently give

    |U_fix,nonconstant[p_e]|/q^2
        <=(5/3)1e17(C/q)^2<2.1e12.                  (10)

The fixed vacuum value is separated, not set to zero; it has no
observation derivative. Local LLL curvature derivatives may be added
with the generous reserve 1e12 in the normalized jet units. They are
not counted as part of the massive gap. Translation covariance on the
plane identifies the parameter jets with the observation jets of
U_fix. Common local counterterms have the same covariance. Uniform
Born/angular convergence justifies the derivatives. Real polarization
then gives for j=1,2,3

    L_j=2^(ceil(j/2)) j^j (2.1e12)/(j!100^j rho^j).   (11)

In particular L3=4.725e15. Unlike the withdrawn proof, both (8) and
(11) differentiate one and the same cap. No unestimated cap-center
term has been discarded.

## 5. Actual source budget and consequence for the physical solve

At the normal center first and second scalar covariant derivatives
equal coordinate derivatives. At order three the connection correction
is bounded by the first normalized derivative divided by 200, as in
the retained proof (the unnormalized coefficient is <300). Therefore

    max(1e14, G1+L1+1e12, G2+L2+1e12,
        G3+L3+1e12+(G1+L1+1e12)/200)<6e15.           (12)

The zero-order input is the retained complete point bound, valid under
the same geometry/scale guards. This proves (1) in canonical tensor
norms on the stated MIDDLE cover. It does not extend that cover to the
cold bolt or to Y<.05, and it is not the full physical shape covector.

In the actual source normalization S_can=2G Theta/(B^2 C),
eta_G=GC/B^2<8.5e-33 and q/C<400 give

    max_(j<=3) |nabla_gamma^j S_can|/(j!100^j)
       <2eta_G(q/C)^2(6e15)<=1.632e-11.               (13)

This is an absolute seed source bound, not a residual after solving the
field equation. If this component alone were assigned the entire
5e-9 correction budget, its complete physical return operator would
need norm <=306.372549... in these source units; half that correction
budget requires <=153.186274.... Those operator bounds remain unknown.
Other components and regions must be included before computing Y0.

For comparison only, in the ORIGINAL first-constraint ledger the
source-dependent terms .06 gamma K_m and 6.6 gamma K_m, gamma<6.8e-28,
would now be <2.448e-13 and <2.6928e-11. They are not a transferred
inverse for the present c_R theory and are not added as if they were.
The improvement changes the actual source estimate, not its normalization.

Continue the direct physical (b,V) graph solve selected by
[the return audit](auxiliary-physical-return-feasibility.md). Use this
fixed-cap method for source derivatives; do not repeat moving-cap
Cauchy without its correction. The unresolved dominant steps are the
cold/outer differentiated sources and the full nonlocal metric inverse
with actual boundary/domain transfer. This source theorem closes a
quantitative input for that solve; it does not satisfy S169.2-S169.5.
Full nonlinear Y0/Z0/Z2, positive admission margins, exact matching,
semiclassical validity and retarded stability remain required.

Remaining to-do list: cold/outer and complete physical source estimates;
actual physical metric/graph solve and domain transfer; nonlinear solution,
matching, validity, all-sector stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
