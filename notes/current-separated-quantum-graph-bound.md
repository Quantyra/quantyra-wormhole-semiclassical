# Differentiated cross-collar quantum response on the actual paired domain

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
Full stable-wormhole objective ACTIVE and INCOMPLETE; feasibility audit CLOSED.

The remaining .187 quantum/global graph budget is not yet certified. The
existing weak-form and H3-to-L2 remainder estimates do not supply the four
observation derivatives needed here. Differentiating their envelopes would
be invalid. Instead this note derives those derivatives for the COMPLETE
separated kernel from its actual first-order realization, then bounds the
cross-collar OPERATOR. This is not another separated point value, a finite
physical-mode approximation or the full near-diagonal quantum estimate.

## 1. Operator, domain, and precisely which contribution is tested

Keep z, Q, B, C, gamma, q, the common regulator and the paired regular-bolt
realization of current-seed-separated-quantum-response.md. Set

    a=.01, b=.02, sigma=4.005, x=log z, D=z partial_z,
    Y(z)=integral_0^z Q(t)^(-1/2)dt, L_bolt=Y(z_h).

Here b denotes the INPUT cut, not the shape field; use v for that field.
Output is on 0<z<=a and input is on b<=z'<z_h. The actual Q obeys
0<Q<.2 on the open physical interval, so their proper separation is at
least d=sqrt(5)(b-a)>0. The actual bolt is included by its regular limit.

Let Hm be the normalized canonical massive shape Hessian (before z4
multiplication), and Hm0 its logarithmic normal operator. Define

    Rfar v=(Hm-Hm0)[1_(b<=z'<z_h) v],
    Afar v=chi_a z^-sigma (D-4)[z4 Rfar v]
           =chi_a z^(5-sigma) partial_z Rfar v.       (1)

Use the SAME W3,infinity OUTPUT cutoff chi_a as in the previous two notes;
its derivative guards are (1,3,10,60). No derivative is applied to the
input indicator. Splitting an already separated bounded integral is not
evaluation of the full action at a discontinuous metric or a new boundary
condition. Start on physical smooth directions, then use the bound below
and dominated convergence. The original paired projectors, mass/log
completion and bolt data are unchanged.

All local seagulls, Weyl, LLL and finite restoration contacts have zero
cross-support action here. They are retained in the COMPLEMENTARY
near/diagonal operator, not deleted from the theory. No independent
regulator is assigned to a mode or a contact. The scalar row of (1) is
zero: the canonical determinant depends only on shape. Its shape input
under the normal coupled inverse is v=z^sigma B_sigma^-1 g, restricted
to the physical interval. Here B_sigma is the previous FULL normal
multiplier, including all its finite contacts.

This does not bound auxiliary continuation beyond z_h, original-row
moments, mass/slow coefficients or their affine returns. Those remain in
the actual global assembly. In particular (1) is NOT the independent
Efar=||z^-3.5 G_R||_(L2,z>1) target from the Mellin moment calculation.

## 2. Actual full-domain kernel and uniform decay without a guessed bolt datum

The previously derived first-order formulas use

 M_n^2=n(n+q)/C, d_n=2(q+2n), p_j=gamma(j+1/2), W=p/sqrt(Q),
 Udot=W(1-U2)-2M U, with 0<=U,V<=1 on the actual paired double.

V evolves by the same equation in the opposite propagation direction.
These invariant intervals follow from the original regular endpoint
branches and transmission; they are not chosen boundary data. In the
orientation A inner, B outer the spin traces are

 Bdir=-4 U_A V_B exp(-2I_AB)/(1+U_B V_B)^2,
 Bimg=-4 U_A U_B exp(-2I_img)/(1+U_B V_B)^2.          (2)

Both thermal signs and the physical-half image normalization are already
included in the strong normalized kernel per input dY':

 Hm(z,Y')=-2G/[pi B^2 sqrt(Q(z))]
             gamma sum_(n>=1,j>=0) d_n M_n^2 [Bdir+Bimg].          (3)

Use (2), not separate absolute values of chiral matrix products.
Since W>=sqrt(5)p, put k=sqrt(M2+5p2) and
U_*=sqrt(5)p/(M+k). Comparison from the least allowed value zero gives
U(t)>=U_* tanh(kt). Consequently, on ANY path segment of length ell,

 integral(M+WU) >= M ell+(k-M)log(cosh(k ell))/k
                         >= k ell-log 2.            (4)

This is valid even without a pre-propagated buffer. For the image use
(4) on its two legs, retaining a loss 2log 2. Its length is Y+Y'>=Y'-Y.
Thus exp(-2I_AB)<=4 exp[-2k(Y'-Y)] and
exp(-2I_img)<=16 exp[-2k(Y'-Y)]. Actual state and bolt dependence remain
in U,V; the estimate holds uniformly for all their allowed actual values.

## 3. Four observation derivatives from analytic ODE germs

Only the FIXED SEED coefficients on the physical side are analytically
continued for this estimate. No smooth double replaces the C2,1 paired
domain and no analyticity of arbitrary mass/log perturbations is assumed.
At each real 0<=z<=a, the actual solution germs U,V satisfy their original
ODE with the one-sided polynomial Q. Continue these germs to

    |zeta-z|<=R=1/[1000(M+p)].                       (5)

M>=18 puts the disk inside |zeta|<.02 and away from all input points and
their images. The actual polynomial coefficients give |Q/beta-1|<.001
on |zeta|<=.02, so take the consistent square root and
|Q^-1/2|<4, |Q^-1|<16.

In |U-U_real|,|V-V_real|<=.1, the z-ODE right sides have magnitude
<=36(M+p) and Lipschitz constant <=36(M+p). Picard's contraction on
(5) therefore gives |U-U_real|,|V-V_real|<=.036<.1 throughout the disk.
In particular |1+UV|>=.79. The exponent derivatives satisfy
|partial_z I|<=18(M+p), with either orientation. The analytic continuation
of each expression in (2) is bounded by

    4(1.1)/(.79)^2 exp(.036) exp(-2I_real)
                                      <8 exp(-2I_real).

Including |Q^-1/2|<4, (4), and both paths, Cauchy's formula gives

 |partial_z^j{(Bdir+Bimg)/sqrt(Q)}|
 <=640 j! [1000(M+p)]^j exp[-2k(Y'-Y)], j=0,...,4.   (6)

Equation (6) comes from the actual ODE, not a differentiated real kernel
envelope. Continuation through the coordinate value zero is only a local
analytic germ used to bound one-sided derivatives. The physical double
still has its cubic seam and original transmission conditions.

## 4. Entire angular and discrete thermal sums

Integrating the right side of (6) over the whole actual input interval
gives at most exp(-2kd)/(2k)<=exp(-2kd)/(2M). Take

    c_M=sqrt(2)d, c_p=sqrt(10)d,
    2kd>=c_M M+c_p p.

The counting function N(m)=sum_(M_n<=m)d_n obeys N(m)<=4C m2.
Integration by parts, dropping only its negative contribution, gives
for every integer r>=1

    sum d_n M_n^r exp(-c_M M_n)
                         <=4C (r+2)!/c_M^(r+2).     (7)

For each integer k>=0, exact midpoint-cell integration of the unimodal
function p^k exp(-c_p p) gives

 gamma sum_(j>=0) p_j^k exp(-c_p p_j)
             <=k! c_p^(-k-1)(1+gamma c_p).           (8)

The error is <=gamma times the maximum. The inequality
(k/e)^k<=k! bounds that maximum; k=0 is immediate. Thus (8) bounds the
DISCRETE sum at the original positive temperature, not a zero-temperature
replacement. No mode is omitted from (7) or (8).

Expand (M+p)^j in (6). Equations (3), (6)--(8) give, for j=1,...,4,

 sup_(z<=a)|partial_z^j Hm[1_far v](z)| <=A_j ||v||_infinity,

 A_j=2560 G C/(pi B^2) j! 1000^j (1+gamma c_p)
       sum_(k=0)^j binom(j,k)(j-k+3)! k!
                         /[c_M^(j-k+3)c_p^(k+1)].   (9)

The positive exponential majorants justify all four derivatives, integration
over the actual bolt limit and common regulator removal by dominated
convergence. This closes that continuum step for this separated operator.
It supplies no coincident finite-part estimate.

## 5. Matched normal kernel and H3 graph return

Away from the diagonal, Fourier inversion of p4 log|p| gives -12/|Y-Y'|5.
The canonical normal coefficient C/(30pi), the physical factor -2G/B2,
Y=z/sqrt(beta) and the even image therefore give

 Hm0(z,z')=12K[(z'-z)^-5+(z'+z)^-5],
 K=G C beta2/(15pi B2),                             (10)

per input dz'. All finite fourth contacts are local, hence absent only
from this separated expression. Integrating differentiated (10) yields

 sup|partial_z^j Hm0[1_far v]|<=N_j||v||_infinity,
    N_j=K(j+3)!/(b-a)^(j+4).                        (11)

Using the larger auxiliary range z'>=b in (11) merely bounds the same
physical restricted input extended by zero; it does not change its domain.
Let F_j=A_j+N_j and alpha=5-sigma=.995. Define C_00=1 and
C_(m+1),j=(alpha+j)C_mj+C_m,j-1, with missing entries zero. For k<=3,

 L_k=sum_(l=0)^k sum_(j=0)^(k-l) binom(k,l) h_l C_(k-l),j
              a^(alpha+j) F_(j+1)/sqrt(2(alpha+j)),
    (h_0,h_1,h_2,h_3)=(1,3,10,60).

This is the product rule in (1) and the exact integral of
z^(2(alpha+j)) dz/z. Thus

 ||Afar v||_H3 <= Lfar ||v||_infinity,
 Lfar=[sum_(k=0)^3 binom(3,k)L_k^2]^(1/2).           (12)

For the physical restriction of v=z^sigma B_sigma^-1 g, z<z_h<1 gives
||v||_infinity<=C3 ||g||_H3/sqrt(2), with the retained C3<1801.
The factor 1/sqrt(2) is the H1(R) point-evaluation constant, and H6
control follows from the complete normal multiplier. No graph norm is
renamed to make its inverse one. Therefore the full separated quantum
feedback block has maximum row sum

    eta_far <=C3 Lfar/sqrt(2).                       (13)

The frozen acceptance target is eta_far<2e-6 on the unchanged .01 output
collar with input from the entire physical region z'>=.02. Also record
Lfar, all F_j and the remaining allocation .187-eta_far.

## 6. Boundary coefficients, research decision, and limits

The same (12) bounds a prescribed finite-dimensional shape contribution
v_F by Lfar||v_F||_infinity. Its sufficient entire-affine-budget threshold
is ||v_F||_infinity<=7.49e-12/Lfar; do not infer its actual value.
In particular, estimating the known B0 z4 and unknown slow component
separately would destroy their needed correlation. For the retained pair

 v_pair=z4[B0(1-z^e_*)+Delta z^e_*],
 ||v_pair||_infinity<=B0 e_*/(4e)+|Delta| on 0<z<1.  (14)

Use 1-z^e_*<=e_* log(1/z) and max z4 log(1/z)=1/(4e). Execute (14) with
the existing enclosed B0,e_* and the CONDITIONAL prior Delta cap
4.500000054e-10. Compare with the bound that incorrectly estimates the
two large coefficients separately. The actual Delta, mass coefficient and
global moment are still unknown; this is a sensitivity calculation, not
evaluation of the physical correction.

RESEARCH DECISION: retain the construction. If (13) passes, remote paired
boundary/bolt response is quantitatively compatible with the contraction
allocation in this graph. The remaining unknown is not licensed by that
pass: the near/diagonal quantum response, finite contacts and transition
input must still be controlled. The old lower-order strong bounds and .0022
form bound do not provide their H3 output estimate.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): estimate the commonly
subtracted near/diagonal Weyl-plus-paired determinant response through FOUR
observation derivatives, including every diagonal finite contact and the
transition input 0<z'<.02, in this coupled graph. Use the remaining
feedback target <=.18699 if (13) leaves that reserve. Preserve the explicit
outer logarithmic cancellation before bounding its terms. A point kernel,
finite angular/thermal selection or supporting checker cannot replace it.

The complete coupled inequality remains conditional on that bound and all
other original-domain contributions, S_coupled<=7.49e-12 and
k_coupled<=1e-14. The independent Mellin Efar<=1e-10 target stays open.
The sufficient signed-projection screen E_rest<2.123e-42<1.287e-40 retains
its scope and unmeasured sign. No actual full correction is computed.

At rho=1e-8 the unproved conditions are still Y0+1e-8 Z0+5e-17 Z2<=1e-8
and Z0+1e-8 Z2<1, with Y0<=5e-9, Z0<=.25, Z2<=5e7 only allocations.
The C2 action on bounded ADMITTED stationary mass/log families does not
prove nonlinear X-ball admission. Exact compatibility, physical matching/
causal realization, semiclassical validity/noise and all-sector stability
remain unresolved. No construction switch or nonexistence theorem follows.

Remaining to-do list: near/diagonal quantum graph estimate and remaining
global assembly; joint source and independent Mellin far-tail bounds;
actual correction and nonlinear admission; compatibility, physical matching,
validity/noise and all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
