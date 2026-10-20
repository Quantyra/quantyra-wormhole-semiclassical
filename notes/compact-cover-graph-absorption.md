# A compact fourth/log graph with an explicit outer interface

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Compact patching;
the global weighted inverse and physical solution remain open.**

The [single-patch estimate](compact-disk-cutoff-commutator.md)
leaves a larger-patch derivative term. Here it is absorbed over
the complete compact region, including the regular bolt. The only
derivative term left on the right is supported in the explicit
outer annulus 7/32<Y<1/4. The norm, overlap, chart conversions and
finite cover size are specified. No end estimate is assumed.

## 1. Cover, partition and statement

Keep the same analytic seed, complete fixed-area disk Hessian L,
canonical metric gbar, period, state and domain. Let

    K={Y>=1/4}, d=1/4096.

Fix any maximal d-separated set of centers p_i in K, using canonical
geodesic distance. At each center use the retained normalized analytic
conformal chart z_i of radius 1/32 and its cutoffs

    chi_i=c_(1/1024)(z_i), eta_i=c_(1/128)(z_i),

where c_a is the beta(6,6) cutoff of the preceding note. Each is
extended by zero outside its chart. The closed thermal bolt is an
ordinary chart center. There are N<1e44 centers, and all overlaps
needed below have at most M=1e6 neighbors, as proved in section 2.
Each point of K lies where at least one chi_i is one.

Define scalar functions on the physical orbit

    s=sum_i chi_i^2,
    a(s)=0 (s<=1/4), P(4s-1) (1/4<s<1/2), 1 (s>=1/2),
    f(s)=a(s)/s for s>0, f(0)=0,
    zeta_0=1-a(s), zeta_i=f(s)chi_i^2.

Thus zeta_0+sum_i zeta_i=1 exactly. The functions are C5 with
bounded piecewise sixth derivatives; smooth approximation within
the strict reserves below has the same limits as in the preceding
note. On K, s>=1 and zeta_0=0. It is important to keep zeta_0;
the compact cover does not cover the whole outer end.

For w=(theta,T), theta=S+T/3, let X4_i be the sum of the
Euclidean theta H4 and T H4log norms in chart i, with exactly
the Fourier multipliers in the preceding note. Define

    P4(w)=(sum_i X4_i(chi_i w)^2)^(1/2),
    B4(w)=(sum_i X4_i(eta_i zeta_0 w)^2)^(1/2),
    N0(w)=||theta||_L2(gbar)+||T||_L2(gbar),
    Omega=union_i supp chi_i,
    F=||L_theta w||_L2(gbar;Omega)+||L_T w||_L2(gbar;Omega).

N0 uses the whole physical half, with the time integral. P4 is
the fixed compact patch norm for this result; it is not a claimed
numerical equivalence with every possible atlas norm. Then

    P4(w) <= (1e8/C) F + 1e-179 B4(w) + 1e4270 N0(w).       (1)

Every summand of B4 is supported in 7/32<Y<1/4.
For inputs supported in K it vanishes, giving a complete compact
source/value graph estimate. For a global input it is an interface
term to be controlled by the outer problem, not a solved boundary
condition. The constants are not a useful nonlinear contraction.

The statement is first an a priori estimate for smooth global inputs
in the retained domain, compact away from the seam, and extends by
graph closure with the displayed norms finite. No new boundary
condition is imposed at the bolt or the seam.

## 2. A finite cover with quantitative overlap

On all the chart regions used here the metric length ratio is between
.9 and 1.1. In particular a canonical d-ball about its center lies
where |z_i|<1.1d<1/1024, so maximality gives the asserted cover.
The disjoint canonical balls of radius d/2 have area at least
pi(d/2)^2/2. They stay inside the physical half.

The TOTAL canonical orbit area is finite, although its proper
length is large at the prescribed cold temperature. If R_h>B is
the thermal horizon radius, the exact area is

    A_orb=beta integral_(R_h)^infinity C r^-2 dr=beta C/R_h.

This is only a covering measure, not the physical spatial area
observable or a wormhole opening claim. Use the retained benchmark:
K0=2000, q=K0^4, eta=1/K0, nu=2q eta G C0/b0^3,
beta<=8pi/nu, b0^2=G K0^10 xi0, xi0=(C0-1)/6<2K0^3,
C/C0=z_star/u_star, B/b0=sqrt(z_star), and
sqrt(z_star)/u_star<2. Consequently

    A_orb <16pi K0^10 <1e36,
    N <=2 A_orb/[pi(d/2)^2] <1e8 A_orb <1e44.          (2)

The thermal seed and benchmark formulas are the ones already fixed
in [the seed](stationary-centered-covariant-seed.md) and
[the recovered parameters](retained-benchmark-parameter-recovery.md).
There is no choice of a new temperature in (2).

If supp eta_i meets supp chi_j, the centers are separated by at
most D=1.1(1/64+1/512). Their disjoint d/2-balls lie in a ball
of radius D+d/2. This lies strictly inside a normalized chart
by first-exit length comparison. Its area is at most
4pi(D+d/2)^2. Dividing by the small-ball lower bound gives

    number of neighbors <=8[(D+d/2)/(d/2)]^2 <1e6=M.    (3)

The same argument bounds reverse neighbors and the number of
chi supports through any point. Thus no estimate below uses N
where a uniform overlap bound suffices.

Finally a point of supp eta_i has
Y>=1/4-1.1/64>7/32. If zeta_0 is nonzero, that point is outside
K, since s>=1 on K. This proves the interface support assertion.

## 3. Quantitative multiplication and chart changes

We need bounds in the FOURTH/LOG norm, not just at integer order.
First bound the local coordinate maps on H3 and H5.
The cutoff jet ledger gives

    |D^j chi_i|,|D^j eta_i| <=j! A^j, A=1e7, 1<=j<=5,

in their own charts, with zeroth norms at most one.
On a relevant overlap, the analytic conformal transition F_ij
has |F_ij'|<2 and

    |D^j F_ij|<=j! b^j, b=1e4, j<=5.

Both charts have margin at least 1/512 around the support on which
the localized operator acts, so these bounds follow from the same
Cauchy argument as in the preceding note. The Jacobian change in
L2 costs less than two. For a composed cutoff, the two-coordinate
chain rule and

    Bell_(j,k)(1!b,...)=j! binom(j-1,k-1)b^j/k!

give a bound j![3Ab]^j, hence at most j! D^j with D=1e12.
All chi_j expressed in chart i obey that bound.

From (3), for j>=1,

    |D^j s| <=M(j+1)!D^j <=j! H^j, H=2MD=2e18.       (4)

On 1/4<s<1/2, use the explicit polynomial P(4s-1) and
differentiate 1/s. Together with the constant regions this gives

    |f^(k)(s)| <=4 k! Q^k, Q=1e6, 0<=k<=5.            (5)

Indeed a sufficient bound there is
k!4^(k+1) sum_(j=0)^k d_j/j!, where d_0=1 and
(d_1,...,d_5)=(11,220,3960,63360,887040).
For s>=1/2 the bound follows directly from k!s^(-k-1).
The composed derivative bound from (4)-(5) is
4j!(2QH)^j. Leibniz then bounds the coefficient

    a_ij=eta_i f(s)chi_j

by |a_ij|<=4 and |D^j a_ij|<=4j! E^j, E=1e28, j<=5.
There are three factors; their normalized derivative sum is at most
4j!(2D+2QH)^j, less than this reserve.

Define the localized coordinate operator

    T_ij v(z_i)=a_ij(z_i) v(F_ij(z_i)),

extended by zero. Its coefficient has support strictly inside both
chart domains. The vanishing cutoff jets make this extension valid
on H5; a global extension of the coordinate transition is unnecessary.
Product rule, the chain rule and the L2 Jacobian bound yield

    ||T_ij||_(H3->H3), ||T_ij||_(H5->H5) <1e160.       (6)

For a reproducible finite reserve, let

    b_0=2,
    b_k=2 sum_(l=1)^k 2^l Bell_(k,l)(1!b,...), k>=1.

Each coordinate derivative of order j<=5 of T_ij v is bounded in
L2 by

    sum_(k=0)^j binom(j,k) 4(j-k)! E^(j-k) b_k ||v||_Hs,
                                                          j<=s.

The coefficient sum of (1+xi1^2+xi2^2)^s is 3^s.
These explicit finite expressions prove (6), with ample margin.

Here is the logarithmic transfer, including its constant. Let P_j
be the orthogonal Fourier projection to
2^j<=sqrt(1+|xi|^2)<2^(j+1), j=0,1,... . If an operator T is
bounded by A0 at H3 and H5, then

    ||P_j T P_k||_(2->2)
       <=A0 min[2^(3(k+1-j)),2^(5(k+1-j))].

On these blocks log(2+|xi|) is between (1+j)/4 and 2(1+j).
Thus the H4log norm is between 1/4 and 32 times the square-sum
norm with weights 2^(4j)(1+j). Its weighted block matrix is bounded
by

    32 A0 2^(-|j-k|)(1+|j-k|).

Both Schur sums are at most
32 A0 [1+2 sum_(n>=1)(n+1)2^-n]=224 A0. Consequently

    ||T||_(H4log->H4log) <=28672 A0 <1e6 A0.           (7)

The same argument without the logarithmic ratio bounds its H4 norm
by 1536 A0. No H4-to-H4log bound is asserted. With (6), both input
components therefore satisfy

    X4_i(T_ij v) <=1e170 X4_j(v).                     (8)

The general interpolation/chart framework is described by
[Mikhailets--Murach, Theorem 5.4](https://arxiv.org/pdf/1206.6041#page=18),
inspected 2026-10-20. Equations (4)-(8) give the constants used here;
an unspecified equivalence of atlas norms would not suffice.

## 4. A sufficiently small one-patch remainder

The preceding commutator proof keeps its auxiliary near radius
explicit until its final substitution. Use the SAME cutoffs and
probes, now split the kernel at R=1e-220 instead of 1e-30.
There is no change in the physical operator. The preserved ledgers
give, successively,

    near coefficient <=1e20 R C=1e-200 C,
    full bounded-tail coefficient <=1e30 R^-4 C=1e910 C,
    each of ten contact coefficients <1e912 C,
    complete H3 commutator coefficient <1e915 C.

For the contacts the actual probe action remains <6e122 C;
its near action is at most 1e20 R 1e32 C. Thus the tail dominates
the same ten-probe extraction. No new coefficient or derivative
bound is inferred from a separated-kernel bound.

The intermediate complete disk graph, with the unchanged cutoff
H3 multiplication bound, now gives

    X4_i(chi_i w)
       <=(2000/C)F_i+2e-197 X4_i(eta_i w)
                      +1e920 U3_i(eta_i w)+2e913 N0(w).

Here F_i is the sum of the Euclidean L2 norms of chi_i Lw.
Interpolate with epsilon=1e-1105 and
U0_i(eta_i w)<=2N0(w). The result is

    X4_i(chi_i w)
       <=(2000/C)F_i+1e-180 X4_i(eta_i w)+1e4240 N0(w). (9)

This tightening is essential: the earlier coefficient 1e-6 was
not proved sufficient after the chart and partition factors in (8).
All earlier successful sources and receipts remain unchanged.

## 5. Absorption, source conversion and the actual linear solution

The exact partition identity gives in chart i

    eta_i w=eta_i zeta_0 w+sum_j T_ij(chi_j w).

The adjacency matrix in this sum has at most M nonzero entries
in each row AND column. Its l2 norm is at most M. Using (8) yields

    (sum_i X4_i(eta_i w)^2)^(1/2) <=B4(w)+1e176 P4(w). (10)

Also the overlap and measure conversion give

    (sum_i F_i^2)^(1/2) <=2 sqrt(M) F=2000 F.

Take the l2 norm of (9) across all centers, using Minkowski:

    P4 <=(4e6/C)F+1e-180 B4+1e-4 P4+1e4240 sqrt(N) N0.

Absorb 1e-4 P4, use sqrt(N)<1e22 and enlarge the constants.
This proves (1). No compact-region derivative remains on its right.

The source conversion from the complete disk graph applies
throughout Omega, whose points have Y>1/8 and B<r<28B:

    F <=1e4 B^2 C^-1/2 Fdisk,
    N0 <=sqrt(C) B^-2 ||(S,T)||_L2(dmu_4).

Consequently

    P4 <=1e12 B^2 C^-3/2 Fdisk+1e-179 B4
                         +1e4270 sqrt(C) B^-2 ||(S,T)||_physical. (11)

For the compatible coupled linear solution, the accepted
[source conversion](weighted-york-source-conversion.md) gives
||(S,T)||_physical<=5.58e9 g lambda^-2 Fdisk, g=G lambda.
Its substitution in (11) bounds the last term by
1e4280 sqrt(C) B^-2 g lambda^-2 Fdisk. The annulus norm B4
remains. No bound on the full fast mass-extended norm, higher
York transfer or residual/inverse composition is concluded.

## 6. Frozen verification and remaining work

Freeze this ordinary proof and its checker before execution.
Pin the cutoff, benchmark and source-conversion receipts and inputs.
Verify the area/count/overlap inequalities, partition identity and
interface support, all native and composed cutoff jets, rational
partition derivatives, Bell identities, finite H3/H5 budgets,
infinite dyadic Schur sums and logarithmic weight bounds, the
retained-probe budgets at the smaller auxiliary radius, and the
final absorption and physical-source arithmetic.
Preserve every run and failed source snapshot.

These finite identities and bounds support the ordinary cover,
operator and domain arguments; they do not formally certify
the continuum theorem. The prescribed low-frequency form bounds
remain preferable to the very large value constant here.

The next mathematical obligation is the quantitative outer
estimate and its interface with B4, followed by higher York maps
and complete inverse/residual composition. The nonlinear map,
signed compatibility and every physical completion requirement
remain unproved.

Remaining to-do list: outer annulus and weighted fast end; higher
York reconstruction and inverse/residual composition; nonlinear
signed compatibility; exterior/state, opening, stability, transport,
preparation, separate creation and full-result audit.
