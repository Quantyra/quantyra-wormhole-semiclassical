# Cold auxiliary gluing derivatives through eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Use the exact seed, angular tower, temperature and subtraction of the
[local source calculation](seed-cold-local-eight-jets.md). The complete
H2 cap is exactly the C10 cap in the
[weighted inverse calculation](seed-cold-weighted-inverse-jets.md).
The local flat cap is exactly the moving normal cap of the local result.

For h<=s_x<=.001 write a=a_x=s_x[1+log(.01/s_x)]<.00331. If G is the
finite global-H2-cap minus local-flat-cap source and G_H2 is its homogeneous
counterpart, the result is

    |nabla^j(G-G_H2)|/(j!100^j q^2 a)<1e44, 0<=j<=8.       (1)

The derivatives are full physical covariant tensors. Adding the already
bounded local difference U_mov-U_H2, with reserve 1e8, still gives 1e44
for the COMPLETE AUXILIARY CAP difference. The eighth gluing entry is
below 2.37e43. The comparison with the actual outer domain is not included.
Thus (1) is not the complete actual cold-source C8 theorem.

## 1. One moving normal chart throughout the locality identity

Keep the real boost family L_b=B1(b1)B2(c_x b2), c_x=sech(chi_x), and
sigma_b(y)=sigma_cap(L_b y). The cap multiplies the analytic seed
continuation by the degree-21 C10 step across a radial interval of width
.01 at chi1=arcosh(.01/h). It has |sigma_cap|<.03 and real gradient <10.
The local flat cap multiplies exp(sigma_N)-1 by eta21(16r-3), equal to one
through 3/16 and zero from 1/4. These two caps remain distinct.

The gluing cutoff kappa is one to 1/8 and falls linearly to zero by 5/32.
Its Lipschitz constant is 32, and its physical commutator norm is <33.
Extend the holomorphic normalizing coordinate by

    F_b(0)=0, F_b'(z)=exp[sigma_b(x)-2sigma_b(F_b(z),0)],
    G_b(w)=w+zeta(|w|)(F_b(w)-w).                         (2)

The H2 conformal factor vanishes on each holomorphic axis. The second
branch is the independent antiholomorphic equation; complexifying b does
not conjugate b. The fixed radial quintic zeta equals one to 5/32 and
zero from 7/32. Its first and second real directional derivative bounds
are 30 and 1500. Its C2 regularity suffices here: pullback coefficients
use DG_b at w, not composition with the spatial inverse G_b^-1. Eight
parameter derivatives therefore do not require eight spatial derivatives
of this fixed cutoff. Outside its support, the global cap is C10.

Let Q_b align the normal spin frame on the inner disk and extend its
rotation angle by zeta. On canonical flat-volume spinors the transported
global operator is

    A_b=P_b^a partial_a+(partial_a P_b^a)/2,
    P_b^a=Q_b exp(-nu_b) gamma^i K_b[a,i] Q_b^-1,
    nu_b=(v_H2+sigma_b) composed with G_b,
    K_b=(DG_b)^-1, S_b=(A_b+M)^-1.                       (3)

The gamma matrices are anti-Hermitian Clifford matrices. The canonical
unitary coordinate pullback is sqrt(det DG_b) times composition with G_b.
Piola's identity div(det DG_b K_b)=0 proves the half-divergence in (3).
For Spin(2), Q_b^-1 dQ_b is a multiple of gamma1 gamma2, which
anticommutes with both gamma matrices. The extra anticommutator in the
transported zeroth term vanishes; derivatives of Q_b within div P_b do
not vanish. Both density and spin terms in (3) are retained. The
dimension-two conformal law is also recorded in
[Walpuski, Proposition 15.5, printed page 79](https://walpu.ski/Teaching/SpinGeometry.pdf).
The coordinate formula and all constants here are derived directly.

At w=0, exp(nu_b)sqrt(det DG_b)=1 by (2). Thus the canonical diagonal
has the required physical endpoint normalization. No extra conformal
endpoint factors are added to S_b: they are already included in A_b.

The first-order locality identity of the
[cold construction](stationary-reference-cold-bolt.md), equation (14),
now gives, per angular mode,

    G_M(b)=-M tr[rL_b S_b rR_b],
    rL_b=Rlocal_b(0,.)[D,kappa],
    rR_b=[D,kappa]Rlocal_b(.,0).                         (4)

Both operators agree on the support of kappa in these coordinates.
The single-commutator remainder is
Rlocal[D,kappa^2]Rlocal/2 and has zero kernel near the diagonal where
kappa=1. The two rows in (4) are varied independently as analytic rows;
one must not differentiate a complex adjoint. In particular there is no
fixed-cap/moving-cap mismatch of the kind exposed in the
[middle correction](seed-middle-splitting-audit.md).

## 2. Uniform and cold-weighted coefficient bounds

Use a large local analytic parameter bidisk R=.005 and a smaller jet
bidisk rho=.0001. Put I=(1-rho/R)^-2-1<.041233. A holomorphic function
bounded by B on the large bidisk differs from its value at zero by at
most BI on the smaller one, by its two-variable Taylor series.

On an H2 holomorphic spatial axis |z|<=.25, the height H=h cosh(chi)
changes by at most .125 H_x. The boosts add at most
1.125[exp(.01)-1]H_x, leaving total <.137 H_x. This is inside the
quarter s-disk inverse-image radius .17(s_x+h). The retained rational
seed formulas therefore give |sigma_b|<.02 and |sigma_b|<=12a there.
Equation (2) then gives |F_b'|<1.062 and |F_b-id|<9a on |z|<=.23.
The image stays inside .245<.25, closing the bootstrap.

At real b=0 the retained normalizer estimates are

    |F0'-1|<.05, |F0''|<2.5,
    |DG0-I|<.439, |DG0|<1.439, |K0|<=2, |D^2G0|<22.

Here D^2G is a bilinear real-direction operator norm, not a Frobenius
norm. On |w|<=7/32 the small-b image remains below
1.05(7/32)+9aI<.232. The spatial analytic margin is consequently .018.
Cauchy and differentiation of the ODE give

    |F_b-F0|<=9aI,
    |F_b'-F0'|<60aI, |F_b''-F0''|<4400aI,
    |G_b-G0|<=9aI,
    |DG_b-DG0|<330aI, |D^2G_b-D^2G0|<21500aI.           (5)

For the first derivative, exp(.06)[36+2(.02/.018)9]<60.
For the second, the change of sigma_z(F) is bounded by
[12/.018+(.04/.018^2)9]aI; apply F''=-2sigma_z(F)(F')^2.
The cutoff adds 30*9 to the first-map bound and 60*60+1500*9 to
the second. Thus, on the small bidisk,

    |DG_b-I|<.485, |DG_b|<1.49, |K_b|<2, |D^2G_b|<25,
    |F_b'-1|<.059, |F_b''|<3.11.                        (6)

Absolute cold-weighted bounds follow on the same disk:

    |F_b'-1|<40a, |F_b''|<1600a, |G_b-id|<=9a,
    |DG_b-I|<310a, |D^2G_b|<17500a, |K_b-I|<620a.       (7)

For the scalar and frame factors, use the full independent spatial
bidisk |z|,|zbar|<=.25. Its centered hyperboloid coefficients satisfy
|A-1|<=.031747 and |B|,|C|<=.253969. Including the boosts, the height
change is <.30 H_x. The half s-disk has

    |sD|<=s_x[4.06 log(.01/s_x)+9]<.0184,
    |s^2 D'|<=36s_x<=.036,
    |H'-1|<.081, |log A_seed|<5s_x.                    (8)

These follow from the exact rational D, D' and q0 expressions of the
cold construction. In detail |q0-1|<3.001s_x and |log q0|<3.011s_x;
the log(1+s) contribution is at most (3s_x-h)/(1-2s_x), and the
bolt q0(h) contribution at most .501h. In differentiating the exact
hyperbolic quotient for H, the small-h remainder stays below 1e-24.
The half-disk image radius is >.459(s_x+h), covering .30 H_x.
Equation (8) gives |sigma_b|<.025 and |sigma_b|<=20a on this full
spatial bidisk. Also |v_H2|<.016 there. With the .018 image margin,
spatial Cauchy gives the following convenient reserves:

    |nu_b|<.07, |grad nu_b|<6.5,
    |nu_b-v_H2|<22a,
    |grad nu_b-grad v_H2|<3500a.                       (9)

For example |grad sigma_b|<=2(.025)/.018, or <=40a/.018 in the
weighted estimate. Combine these with |DG_b|<1.49, (7),
|grad v_H2|<.13 and |Hess v_H2|<2.

The complexified normal-frame angle has |theta|<.062 and |theta|<=36a.
For Q_b=exp[(gamma1 gamma2)zeta theta/2], (6)-(7) imply

    |Q_b|, |Q_b^-1|<1.032, |Q_b-I|<20a,
    |Q_b^-1 dQ_b|<3, |Q_b^-1 dQ_b|<1500a.              (10)

Use half of 3.11/.941+30*.062 for the uniform last bound, and half
of 1600/.941+30*36 for its weighted version. The angle uses independent
logs of the two holomorphic derivatives, preserving analytic b dependence.

Here are bounds for the complete differential coefficients of
A_b-D_H2 on the compact patch:

| coefficient | uniform | cold-weighted |
| --- | ---: | ---: |
| principal contraction | 3 | 2000a |
| divergence | 520 | 400000a |

For clarity about complex norms, the Piola contribution satisfies
|Cliff(div K)|<=4 |K|^2 |D^2G|. This includes the conversion from
complex directional bounds to a two-component Euclidean gradient;
using the smaller real-vector factor here would not be justified.
Since exp(.07)(1.032)^2<1.15, the uniform divergence is below

    1.15[4*4*25+1.415*2*6.5+2*3*2*1.415]+.3<520.

For its weighted version the leading term is 4*4*17500a. The scalar
gradient difference contributes below 11000a before the 1.15 factor,
the differentiated rotation below 8500a, and the remaining reference
coefficient differences below 1000a. Their total is below 400000a.
The principal bounds follow from (6)-(10); its weighted product sum
is below 2.3(620+44+80)a<2000a.

Outside the coordinate extension, use the exact conformal expression
exp(-sigma)D_H2-exp(-sigma)c(d sigma)/2. The true local cap branches
have |sigma|<.03, |grad sigma|<12, and weighted bounds
100a exp(d), 2000a exp(d). In the inner region use |sigma_s|<49a/s,
|ds/dH|<2 and |dH|<=2H; in the taper use its explicit derivative
bound <4.6 and width .01. These give the stated larger gradient
allowances. Outside contributions to the graph bounds below are <.5
uniformly and <300a exp(d) with the cold weight.

## 3. Real graph inverses and all ordered derivatives

Freeze W=exp[-2d_H2(x,G0(w))] in the transported coordinates. This is
the pullback of the original real weight, not a varying b weight.
On the compact patch |d log W|<3 since 2(1.439)exp(.04)<3; outside
it has its usual bound 2. Truncate the distance weight first if needed.
For the real inverse S0 and the homogeneous inverse SH,

    ||W^(+/-1) S0 W^(-/+1)||<=d/M,
    ||W^(+/-1) SH W^(-/+1)||<=d/M, d=6/5.              (11)

For S0 this is the unitary real G0,Q0 transport of the physical cap
bound. For SH use 1/(M-3)<1.2/M, since M>18.

Use E_M(v)=||grad_H2 v||+M||v|| on the fixed reference Hilbert space.
The graph inverse bound for both conjugations and both operators is

    ||W^(+/-1) S0 W^(-/+1)||_(L2->E_M)<12=G.           (12)

An explicit estimate follows by conformal factorization, D S=I-MS
and the H2 Lichnerowicz identity. Before coordinate transport the
weighted gradient is bounded by
exp(.03)[1+exp(.03)d+10d/M]<4; the 10 includes half the real sigma
gradient, the distance-weight commutator and the curvature allowance.
Flat-volume conversion gives <4.2. Coordinate and spin transport cost
at most 1.439*4.2+(44+3)d/18. Allow another 3d/18 for the weight and
d for the M L2 term. The sum is <12. The homogeneous estimate is
smaller. These are real base inverse estimates, not a complex inverse.

Flat half-density derivative conversion on the compact patch costs
at most 1.021 on the gradient and .2 on the value. The coefficient
table then bounds each conjugated difference vertex by

    3(1.021)+3(.2)/18+520/(2*18)+3*3/18<20=u.

The double-weighted vertex has bound <25000a=w a, using 2000 and
400000 in the same expression. W^2 exp(d)<=1 retains the cold weight
in the noncompact part. Applying coefficient Cauchy with
k=15000>sqrt(2)/rho gives, for the true j-th parameter tensors,

    ||W^(+/-1) A^(j) W^(-/+1)||_(E_M->L2)<=u j!k^j,
    ||W A^(j) W||_(E_M->L2)<=w a j!k^j, 1<=j<=8.      (13)

The same bounds at order zero hold for A0-D_H2. The local parameter
branches and C10 global endpoint matching give a C8 real graph family;
they do not assert a globally holomorphic cutoff. The common graph
domain and (12) justify differentiating its inverse.

In each ordered inverse derivative keep one weighted vertex and bound
the others by their conjugated uniform norms. For j>=1 the sum over
ordered positive compositions of j yields

    M||W S^(j) W||<=a j!k^j D_j,
    D_j=d G w (1+G u)^(j-1)=360000*241^(j-1).          (14)

There are binom(j-1,m-1) compositions of length m. The factorials of
their vertices cancel those in the differentiation coefficient. The
zeroth resolvent identity gives the same bound for S0-SH with
D0=d G w=360000. Only one cold factor is used in each product; no
unweighted global O(a) claim is made.

## 4. Both separated rows and the angular sum

The retained local flat Green bounds give, for each row in (4),

    ||W^-1 r||<=R_M=400 exp(-M/32).                    (15)

They hold throughout the original analytic normal-parameter polydisk
|u0-1|,|p|,|pbar|<=.01, including the auxiliary thermal disk used in
the local source calculation. The weight is frozen at G0; on the row
annulus d_H2(G0w)<=1.06|w| and the extra weight factor is <1.02.
No large-mass derivative is taken inside a free point singularity.

On the larger local observation bidisk |b_i|<=.01, the preceding
calculation gives |u0(b)-1|<1.1s_x. The centered H2 frame has
|h e0^T L_b e1|<1.02s_x and |h e0^T L_b e2|<.011s_x, whence
|partial_z H|<.5155s_x. The exact chain
p=exp(-sigma_b(x))(du/ds)(ds/dH)partial_z H and the inverse-height
bound give |p|,|pbar|<1.1s_x. All three normal parameters stay within
.0011, leaving margin .0089 to their original polydisk.

The parameter distance from the homogeneous data (1,0,0) is <3.3s_x.
Integrate parameter derivatives of each row along that straight path.
The row difference is <(3.3/.0089)s_x R_M, plus the thermal coefficient
comparison (<1e-20 a R_M). Reserve 1000a R_M. Observation Cauchy gives

    ||W^-1 r^(j)||<=j!150^j R_M,
    ||W^-1(r-rH)^(j)||<=1000a j!150^j R_M, 0<=j<=8.   (16)

The homogeneous rows are constant in b. The thermal comparison is an
auxiliary local calculation only; it does not vary the global HH state.

Let epsilon=150/15000=.01. Expanding both rows and the central inverse
in (4), retain the weighted inverse whenever its derivative is positive.
When that derivative is zero and j>0, at least one differentiated row
supplies the cold factor. At j=0 subtract the homogeneous product in
three terms. The exact majorants are

    T0=D0+2*1000*d=362400,
    Tj=sum_(b=1)^j (j-b+1)epsilon^(j-b)D_b
                         +1000*d*(j+1)epsilon^j, j>=1.

Thus |partial_b^j(G_M-GH_M)| is at most
2*400^2 a j!15000^j Tj exp(-M/16). The factor two is the spin rank;
the external M in (4) cancels the inverse M. The full retained angular
bound sum d_n exp(-M_n/16)<4q and C/q<.00352 give

    |partial_b^j(G-GH)|/(j!100^j q^2 a)<=H_j,
    H_j=4505.6*150^j*Tj.                              (17)

This sums every angular sector. The inverse bound alone would not
have supplied the exponential needed in this step.

## 5. Full physical derivatives and independent locality control

Use the exact observation connection from the local-source calculation,
including its nonzero zeroth Christoffel term. With gamma_a=20a!300^a,
define B(0,k)=k!100^k H_k and

    B(r+1,k)=B(r,k+1)+r sum_(a=0)^k binom(k,a)gamma_a B(r,k-a),
    J_j=(6/5)B(j,0)/(j!100^j).                        (18)

The sum starts at a=0. It controls all tensor components rather than
only a symmetrized radial jet. The physical frame factor 6/5 is the
retained one for j<=8. Exact rational evaluation gives J8<2.37e43
and max J_j<1e44. Adding the local allowance 1e8 leaves this reserve
strict, proving (1) and the stated auxiliary-source consequence.

An independent first-order locality control uses a circle of length 4,
D0=partial_t and Dalpha=partial_t+i alpha 1_[2,3], with M>0. Choose
x=1 and a cutoff rising with slope 4 on [.5,.75] and falling with
slope -4 on [1.25,1.5]. Both operators agree on its support. Put
r0=exp(-4M), ra=exp(-4M-i alpha). On the two row intervals the Green
kernel is exp[-M(z-w)]/(1-ra) for z>=w and
ra exp[-M(z-w)]/(1-ra) for z<w. The diagonal difference is
(ra-r0)/[(1-ra)(1-r0)]. The four double-commutator blocks, with common
denominator (1-r0)^2(1-ra), are r0(1+ra)/2 twice, -ra and -r0^2.
Their sum is the negative diagonal difference, as required by (4).
Direct kernel quadrature at 60/90 digits is compared with that identity.
This tests a true local first-order operator; a generic matrix graph
does not have [D,kappa] commuting with kappa and is not a substitute.

The supporting checker also tests the Piola/density formula on an
explicit nonlinear diffeomorphism, the nonzero omitted density and
spin terms, coefficient margins, every ordered inverse composition,
and the full row/angular/covariant ledger. These checks do not prove
the continuum estimates by themselves or sample the actual charged
quantum source. No Lean certificate is asserted.

The original complete cold C3 and endpoint C5 contracts are unchanged.
The middle/outer complete C8 bounds remain available. Actual-domain
cold C8, quantitative metric order ten, the finite-metric H6 source and
the sufficient transfer target remain open. No signed compatibility,
nonlinear geometry or physical stability conclusion follows yet.

Remaining to-do list: actual-domain cold source comparison through eight;
metric order ten and finite-metric quantum/Schur estimates; transfer,
cancellation, signed compatibility and nonlinear existence; physical
state/boundary, stability, matching, charge/return, opening, transport,
preparation, separate creation and audit.
