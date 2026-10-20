# A seam-safe annulus and bolt graph with a collar overlap

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof and
verification contract before execution. **A global interior graph,
with an explicit collar interface; the full inverse remains open.**

The earlier compact graph starts at Y=1/4. It cannot be applied to
the intervening annulus by merely changing that endpoint: its charts
would cross the seam, and its Weyl allowance |b|<=4 would fail.
Here both restrictions are replaced quantitatively. The resulting
interface is wholly inside the region where the accepted outer
cutoff equals one. No new seam or bolt boundary condition is used.

## 1. Fixed cover and statement

Keep the physical action, actual massive tower, state, period and
paired domain of [the complete disk graph](complete-compact-disk-graph.md).
Fix, rather than vary, the following dimensionless canonical scales:

    epsilon=10^-160, r=epsilon/32, t=10^-170,
    K={Y>=r/4}, d=t/4096.

Take a maximal d-separated set p_i in K, normalized conformal
charts z_i of radius t/32, and cutoffs

    chi_i=c_(t/1024)(z_i), eta_i=c_(t/128)(z_i).

The function c_a is the same beta(6,6) radial cutoff, one up to a
and zero at 2a. Use precisely the partition of
[the compact cover](compact-cover-graph-absorption.md):

    s=sum chi_i^2, a(s)=0 for s<=1/4,
    a(s)=I_(4s-1)(6,6) for 1/4<s<1/2, a(s)=1 for s>=1/2,
    f(s)=a(s)/s, f(0)=0, zeta_0=1-a(s), zeta_i=f(s)chi_i^2.

For w=(theta,T), theta=S+T/3, use the SAME unscaled Euclidean
Fourier norms in each z_i chart: X_i is the sum of theta H4 and
T H4log, with logarithm log(2+|xi|). Define

    P(w)=(sum X_i(chi_i w)^2)^(1/2),
    B(w)=(sum X_i(eta_i zeta_0 w)^2)^(1/2),
    N(w)=||theta||_L2(gbar)+||T||_L2(gbar),
    Omega=union supp chi_i,
    F(w)=||L_theta w||_L2(gbar;Omega)+||L_T w||_L2(gbar;Omega).

All gbar norms include the time integral on the physical half.
L is the complete canonical fixed-area Hessian, without 2G/B^2.
Then

    P(w) <= (10^8/C) F(w) + 10^-18000 B(w)
                                      +10^381000 N(w).       (1)

Every summand in B is supported in r/8<Y<r/4. The entire outer
transition support r/4<=Y<=16r is in K. In particular the B
support lies where the preceding outer cutoff chi=c_r is one.
The regular bolt belongs to the cover. This closes the previously
uncovered annulus in the canonical patch graph.

Equation (1) is initially an a priori estimate for smooth fields
in the retained domain, and extends by graph closure with its
displayed norms finite. A mass-extended field can be used: all
charts and probes stay a positive distance from the seam. This
does not assign its cubic mass term a fast weighted norm.

The constants are deliberately very large. Matching B and the
outer transition norm requires quantitative one-dimensional to
two-dimensional chart transfers, including the logarithmic norm
and half-density transformation. The support inclusions alone
do NOT complete that step or establish a global inverse.

## 2. Small charts without crossing the physical seam

For centers Y>=1/16, restrict the retained cold/middle analytic
normal charts. Their radius 1/32 closures stay inside the half;
|exp(sigma)-1|<.026 on the complex bidisk. Bolt charts are regular.

For 0<Y<=1/16, use the one-sided analytic seed, not its C2,1
double, to construct an auxiliary analytic chart. The exact seed
argument in [the outer strong remainder](outer-strong-curved-quantum-remainder.md)
gives analytic j on |Y|<.1, |j-1|<.001. At Y0 put

    x(Y)=integral_(Y0)^Y j(Y0)/j(v) dv,
    sigma(x)=log[j(Y(x))/j(Y0)].

On |Y-Y0|<.02, |x'-1|<.003, so integration on chords gives
injectivity. On its boundary |x|>.019, and the argument principle
gives an analytic inverse for |x|<.01. There |sigma|<.003.
Cauchy on the remaining .01 x-margin bounds |sigma_x(0)|<1.
The real coefficient b=sigma_x(0) is removed by H(z)=z-bz^2/2:

    sigma_N(z,w)=sigma((H(z)+H(w))/2)
                +[log(1-bz)+log(1-bw)]/2.

Here z,w are independent complex variables; on the real slice
w=conjugate(z). This has zero value and first jet at the center.
On |z|,|w|<t/32 its exponential differs from one by less than
.026. All its ACTUAL uses have Y>0: the chart image is within
1.1t/32 of a center Y>=r/4. No analytic continuation of the
physical seam metric is claimed or used as the actual resolvent.

Both constructions give |sigma|<.027, length ratios .9..1.1,
and Cauchy jet reserves j! (1000/t)^j through the needed orders.
The retained real quadratic normal bounds also hold on these
restrictions. Restriction and the explicit quadratic normalization
give a uniformly bounded seed/cap gradient; the coarse bound 10
suffices. All cap gradients of directions are included below.

The seed inequality Y(u)<3u for 0<=u<=1/4, with u=B/r_phys,
implies u>r/24 on all observation/probe charts, since Y>r/8.
For u>=1/4 the same lower bound is immediate. Thus

    u>r/24>t, 0<=omega=-log u<400, |omega+1/6|<400.  (2)

The last strict bound follows from log(24/r)+1/6<400; it does
not use the old |b|<=4 hypothesis. Below b denotes omega+1/6,
independently of the temporary chart-normalizing coefficient.

## 3. Actual first Cauchy estimate at the smaller radius

This step prevents an unsupported use of the fixed-radius
interior source theorem near the seam. Repeat the fixed-coordinate
calculation in [S145](general-interior-dirac-response.md), with
cap one through t/64 and zero at t/32. Let U be its maximum of
actual and capped value/gradient norms and local coordinate jets
through five, all in the physical, unscaled coordinates. Use the
parameter disk |s|<=1/(200U), with background zero.

The reference mass m=M exp(phi(0)) and |h|<.06 are unchanged.
For the global Taylor division, multiply S145's coefficient of
degree k by t^-k and its degree-j remainder coefficient by t^-j.
The proof remains valid on the smaller inner disk; outside it
the division point is now t/64. Its product recursion of total
degree five is homogeneous of degree five in these scale factors.
Consequently

    |E_scalar,cap|<10^17 t^-5 M^-3,
    |T_a2|<10^17 C t^-4.                             (3)

The second bound uses total jet order four, including every
gradient term of Q_phi, rather than four separate worst-jet bounds.
The common vacuum scalar variation is still zero, not its value.

For the actual/cap comparison, set Q=M+256/t. The local evaluation
disk has radius 1/(10Q), hence fits before the cutoff commutator
whose distance is at least t/256. The proof of the point estimate
in [S143](curved-remote-dirac-response.md), with Q in its cutoffs
and in (-Delta+Q^2)^-1, gives

    |psi(0)|<=10^6 Q ||psi||_L2,
    ||1_0 R 1_separated||<=4e6(Q/M) exp(-M distance/2)
                          <10^8 t^-1 exp(-M distance/2).

Indeed the potential bound is <4M^2<=4Q^2, the Caccioppoli
and second-cutoff budgets stay below 30000Q^2, Q/M<16/t,
and the exponential evaluation-disk loss is <exp(1/18).
The global weighted resolvent remains the ACTUAL one with norm
<3/M; t does not replace the physical mass or its spectral gap.
Two evaluation rows, rank two, gradient 256/t and the scalar
source factor <1.02M now give

    |E_scalar,actual-E_scalar,cap|
                         <10^20 t^-3 M exp(-Mt/256).       (4)

The full angular bounds sum d/M^3<2C/3 and
sum d M exp(-delta M)<=24C delta^-3 apply at delta=t/256,
including the modes with Mt<1. First-derivative Cauchy in (3),(4)
therefore proves the deliberately enlarged bound

    |D T_mass[0][v](0)|<10^40 C t^-10 U.             (5)

The explicit budget is
200[(2/3)10^17 t^-5+10^20 t^-3*24(256/t)^3
     +10^17 t^-4], below (5) for 0<t<=1.
H_mass v=D T_mass[v]+2v T_mass. Use (5) only on probes whose
center value vanishes; keep the constant contact separately.

## 4. Complete massive graph and the larger Weyl allowance

First bound E=H_mass-Q_geo on the new chart outputs. Use exactly
the fixed degree-two allocation of
[the strong interior proof](strong-interior-quantum-remainder.md),
now split at R=t/400. The near kernel is <=10^8 C/d^4;
Taylor through degree two costs at most 2e6 C t ||nabla^3 v||.
Both full tail rows, including actual gluing and geodesic returns,
fit 10^30 C t^-4 after row subtraction.

The same five linear/quadratic probes with radii scaled by t have
U<10^18 t^-6. Their full geodesic jets obey
sup|X^3 q|<10^32 t^-6 and sup|X^4 q|,sup|X^6 q|<10^46 t^-12.
These bounds follow from the retained finite Bell/Leibniz ledger:
at derivative order n<=6 there are at most n coordinate derivatives
and n powers from the geodesic-jet recurrence; 0<t<=1 permits the
common t^-12 bound. The cap is one on the probe support.
Equation (5), the full Q_geo point action, Taylor and tail bounds
fix each of the five first/second contacts below 10^2800 C.
The constant contact obeys
|p0|<=10^20(q^2/C)(1+1/Y)<10^200 C.
No first-contact probe is divided by a derivative that can vanish.

Sum the contacts and convert covariant H3 to the chart H3.
The conservative conversion 10^12 t^-6 includes the measure and
connection jets. This yields

    ||1_O(H_mass-Q_geo)v||_2 <10^4000 C ||v||_H3,     (6)

for v supported in V={|z|<t/256}, O={|z|<t/128}.

Let K(|D|) be the complete positive flat massive multiplier in
[the logarithmic graph proof](compact-massive-logarithmic-graph.md).
Keep its exact bounds 4K<=pK'<=6K, kappa4>13.9, and
K>=Cp^4(log p+9-38/15)/(30pi). Freezing the bounded conformal
scale still costs <.36||Kv||. For the larger |b|<400, split the
INPUT frequency at p=exp(2000). The high-frequency part of
C b exp(-4sigma)Delta^2/(30pi) costs <.23||Kv||, because
400 exp(.108)/(2000+9-38/15)<.23. Its low-frequency part is
bounded by 10^4010 C ||v||_2. Multiplication by variable b occurs
after the Fourier split; no commutation with b is assumed.

The comparison Q_geo-P, P=K(exp(-sigma(z))|D|), uses its full
fourth allocation. With |Dsigma|<20/t its near separated kernel
is bounded by 10^12 C t^-1 |z-z'|^-5. Split at R0=10^-1000,
which is smaller than t/1024. Taylor through three has coefficient
<10^-6 ||Kv||. Both full tails fit 10^30 C R0^-4=10^4030 C.
The ten identity-jet probes now have jets <=10^32 t^-6 and full
geodesic fourth/sixth jets <=10^56 t^-12. Evaluation of the two
complete model actions fixes every contact below 10^4040 C.
Thus the total Q_geo-P lower term is <10^4050 C ||v||_H3.
The flat output outside O has gap t/256 and costs <10^13 C t^-4
||v||_2. Together with (6) this leaves more than one quarter of
the flat graph norm. No angular mode or exterior tail is dropped.

The complete local disk decomposition remains the exact one in
the predecessor: fourth diagonal C exp(-4sigma)Delta^2/(20pi)
for theta and H_mass-Cb exp(-4sigma)Delta^2/(30pi) for T.
All other terms have input order <=3. Assign background jets the
weights j!(10^6/t)^j, E=exp(-2sigma)<2, Z=exp(2omega)<t^-2,
|omega|<400 and the unchanged benchmark constants. The full
two-dimensional Euler ledger, including the radius and density
contacts, bounds the remaining coefficient sum by 10^2000 C.
The Weyl terms have total background order <=4 and at most one
undifferentiated omega; the classical terms have order <=2 and
Z degree <=2. The checker also evaluates the reweighted complete
operators directly. This verifies the physical allocation rather
than selected principal entries.

Using the flat H4log graph constant 400/C and the scalar H4
estimate therefore gives the supported complete graph

    X_i(v)<=(2000/C)||Lv||_L2(O),sum +10^4120 U3_i(v), (7)

where U3 is the sum of the two H3 norms. The numerical constants
are uniform over every new center, including the bolt.

## 5. Commutator, cover and absorption

Repeat [the actual commutator proof](compact-disk-cutoff-commutator.md)
with the scaled cutoffs and probes. The Lipschitz constant of chi
is <10^6/t, so the near commutator kernel is <10^18 C t^-1/d^5.
Choose the auxiliary split R=10^-20000. Taylor through three gives
10^20 t^-1 R C<10^-19800 C times X_i(eta_i w). The full tail
has both rows <10^30 C R^-4=10^80030 C.

To extract its ten contacts use v_x=(chi-chi(x))q_alpha, so the
center value and the density contact vanish exactly. Products of
the scaled probe/cutoff jets and recentering by the analytic chart
transition give U(v_x)<10^100 t^-12. Equation (5) bounds the
actual commutator point action by 10^140 t^-22 C<10^4000 C.
The near action is smaller; the full tail dominates. Each contact
is <10^80040 C, and the total H3 coefficient is <10^80050 C.
The local differential commutators and the cutoff multiplication
in (7) fit the same reserve. Consequently

    X_i(chi_i w) <= (2000/C)F_i+10^-19700 X_i(eta_i w)
                   +10^80060 U3_i(eta_i w)+10^80040 N(w).

Use U3<=h X+(2/h)^3 U0 with h=10^-100000, and U0_i(eta_i w)
<=2N(w). This gives the simple one-patch estimate

    X_i(chi_i w)<=(2000/C)F_i+10^-19000 X_i(eta_i w)
                                            +10^380100 N(w). (8)

These auxiliary radii change the decomposition, never the operator.

The total orbit area remains <10^36. The same packing proof with
d=t/4096 gives N_centers<10^44 t^-2=10^384, and each forward
and reverse neighbor count remains <10^6. Each point of K is
where a chi is one. The eta supports lie above r/8. Wherever
zeta_0 is nonzero they are below r/4, proving the interface claim.

For clarity retain explicit chart-conversion costs. Native cutoff
jets are <=j!(10^7/t)^j, and transition jets <=j!(10^4/t)^j.
The compact-cover Bell ledger consequently gives composed cutoff
constant D=10^12/t^2, partition constant H=2e18/t^2, and localized
coefficient constant E=10^28/t^2. Its finite product/chain budget
on H3 and H5 is <10^160 t^-10. The retained dyadic argument then
bounds H4log and H4 by <10^170 t^-10=10^1870. Thus

    (sum X_i(eta_i w)^2)^(1/2)<=B(w)+10^1876 P(w),
    (sum F_i^2)^(1/2)<=2000 F(w).

The 10^6 factor here is the two-sided adjacency norm, not the
total number of patches. Sum (8), absorb
10^-19000*10^1876<10^-17000, and use sqrt(N_centers)<10^192.
Enlarging the forcing, interface and value coefficients proves (1).

The original physical source conversion is now, using u>r/24,

    F(w)<10^330 B_phys^2 C^-1/2 Fdisk,
    N(w)<=sqrt(C) B_phys^-2 ||(S,T)||_physical.        (9)

B_phys is the seed radius called B elsewhere; it is not B(w).
The measure/dual conversion costs less than 10 u^-2, so (9)
follows directly from 10(24/r)^2<10^330. In particular (1)
has forcing coefficient 10^338 B_phys^2 C^-3/2. The accepted
compatible-solution value estimate gives the further value term
10^381010 sqrt(C) B_phys^-2 g lambda^-2 Fdisk. B(w) remains.

## 6. Frozen verification and next obligation

Freeze this proof and its checker before scientific execution.
Pin all predecessor inputs and receipts. Verify the seam margin,
complex chart inequalities, scale-homogeneous Taylor recursion,
small evaluation-disk budget, full angular Cauchy sum, probe and
contact reserves, high/low Fourier split, reweighted complete
Euler operators, logarithmic chart transfer, packing, interpolation
and absorption. Preserve every run, failure and source snapshot.
These checks support the ordinary continuum argument; they are
not a formal verification of an elliptic or nonlinear theorem.

The remaining interface conversion must use the actual normalized
stationary time measure, radial/logarithmic coordinate change,
half-density multiplier and specified mass column. It must be
proved with constants before absorption with the outer graph.
The existence of overlapping supports is insufficient by itself.
Higher reconstruction, inverse/residual composition, nonlinear
signed compatibility and the physical goal remain unresolved.

Remaining to-do list: quantitative collar/patch norm transfers and
global graph absorption; higher reconstruction and inverse/residual
composition; nonlinear family and signed compatibility; exterior,
opening, stability, transport, preparation and creation requirements.
