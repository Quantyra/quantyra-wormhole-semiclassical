# Ninth seed-source and eleventh original metric derivatives

2026-10-21, E003/S167. Ordinary argument before supporting verification.
This extends the existing regional source and original finite-comparison
metric estimates by one order. The physical seed, angular tower, HH state,
paired domain, mass, common subtraction and finite path are unchanged.
The need for this order is the off-shell return in the
[Schur derivative audit](schur-source-derivative-budget.md).

Normalize source derivatives by j!100^j q^2. Through j=9 the complete
middle source has reserve 1e48, the cold product-subtracted source has
reserve 1e48 a(s), and the outer remainder after S=c C^2/(5pi Y) has
Euler and weighted ordinary reserves 1e48 L(Y). Here
a(s)=s[1+log(.01/s)] and L(Y)=1+log(1/Y). All original lower-order
estimates, including source C8 reserve 1e44, remain available.

Normalize original metric logarithm derivatives by j!100^j and retain
the finite path parameter ell. Through j=11 the physical covariant
reserves are middle 1e20 ell, cold 1e27 ell a(s), and outer
1e17 ell Y^(4-j)L(Y), the last after the correlated mass cubic.
These are large upper derivative allowances, not metric amplitudes.
Original C10 and smaller low-order contracts are retained separately.
At the outer endpoint this is conormal control on the open collar;
it does not assert ordinary closed-boundary C11 regularity.

## 1. An explicit smoother auxiliary step

Replace the auxiliary degree-21 step used in the three source proofs by

    eta(t)=sum_(i=0)^11 binom(23,i)t^i(1-t)^(23-i), 0<=t<=1,
    eta=1 for t<=0, eta=0 for t>=1,
    eta'=-c t^11(1-t)^11, c=23!/(11!)^2=16224936.

The piecewise step is C11; both endpoint jets through eleven match.
No globally analytic cutoff is assumed. Analytic polynomial branches
at each real point bound only the finite jets that match the true step.
The finite-jet argument of the preceding source notes therefore applies
at seams as well as in the interiors. The new auxiliary caps are used
consistently in every term of each exact locality identity. This is not
a change of physical source or an inference from the old cap's regularity.

## 2. Middle source

Use the exact moving-cap, half-density and spin formulas of
[the middle C8 proof](seed-middle-eight-jets.md). The physical-center
bidisk still has radius rho=1e-5. For the step argument displacement
delta=.003 put B=1/4+2delta^2. The first derivative is below 4, and

    |eta''|<=22c(20B/21)^10/8+44c delta B^10<40.

This follows by maximizing v(B-v^2)^10 and using sqrt(B/21)<1/8.
The same complex radial square-root estimates give the new branch table:

| Cutoff | value | gradient | Hessian |
| --- | ---: | ---: | ---: |
| chi, scale 128 | 1.011 | 520 | 720000 |
| zeta, scale 64 | 1.006 | 260 | 210000 |

The normalizer itself has the old displacement 1.5rho, first change
.2rho and second change 10rho. Hence DG-I<400rho, D2G<316000rho,
the inverse norm is below 1.005, and the inverse displacement below
1.6rho. The base scalar gradient is below 1.12 and Hessian below 1600.
Before transport its parameter gradient change is below 170rho; after
transport the scalar change is below 2.2rho and gradient below 3300rho.
The spin gradient change stays below 65rho. In the exact operator

    A_b=P_b^i partial_i+(partial_i P_b^i)/2

these give principal change below 900rho and divergence change below
1.2e6rho. Thus 900rho+1.2e6rho/(2*18)<A=7/20. The old graph inverse
4, weighted inverse 2.01/M, weighted gradient 4.1 and incoming
divergence bound 5.2 remain valid with the larger base gradient.

Nine center derivatives plus the two spatial map derivatives use at
most eleven cutoff derivatives. Cauchy therefore gives true ninth
vertices K_j=A j!k^j, k=150000. The ordered inverse majorant is still

    V_r=A r!k^r(1+4A)^(r-1).

The [interior row recurrence](seed-middle-coercive-jets.md) already
covers observation order nine and its two extra energy orders. Its
incoming-gradient version uses the same Hilbert-valued recurrence after
the separated duality estimate in [the cap proof](seed-middle-cap-transfer.md),
not an additional observation derivative. Thus set, now through nine,

    c0=30000, c_j=520000(j+3)[600(j+3)]^j for j>=1,
    P_j=sum_(a=0)^j binom(j,a)c_a c_(j-a),
    kappa_j=82*25*88^4*(31/2500000)/(j!100^j)
                *sum_(r=1)^j binom(j,r)V_r P_(j-r),
    G_j=8*256^2*88^2*(31/2500000) P_j/(j!100^j),
    H0=1e14, H_j=G_j+3e13*1500^j+kappa_j for j>=1.

All angular sums and local restorations are exactly the previous ones.
For the physical tensor conversion use Gamma0=0,
Gamma_a=.134 a!100^a for a>=1, and the recurrence in section 5 below.
The resulting maximum through nine is below 1e48. Lower coefficients
through eight are identical to the old ledger, despite the smoother cap,
because its new low coefficients fit the same operator reserves.

## 3. Outer source

Use both caps and the exact scaled heat/local restoration identity of
[the outer C8 proof](seed-outer-eight-jets.md). Taper w=f'/f with eta23.
On a width-half interval the step derivatives are below 8,160,60000;
on the width-.1 negative blend the third derivative is below 8e6.
The real seed w,w',w'' guards still give f'/f<2, f''/f<20, f'''/f<400.
On a complex spatial disk of radius .01, |eta'|<4.1, |eta|<1.082,
|w|<4.328, and |f/f_real|<exp(.04328)<2. Therefore

    |f^(m)|/f<=2m!100^m, m<=12.

The cap f is now C12. Replace the negative finite dilation family by
sum_(j=0)^9 lambda^j D^j f/j!, with D=Y partial_Y. Its ninth parameter
jets and three spatial derivatives are the true ones. The explicit
Stirling/Leibniz sum gives the same relative perturbation bounds
.003,2.5,150,50000 through spatial order three. The negative blend
remainder is below 1e6*.001^10/(1-.001), and its product cost below
1e8. The positive analytic displacement and differences remain
.00005 and (.01,.4,50,7000); the new cutoff derivatives still give
global complex guards 5,1000,1e5 and phase below .05.

The angular weight is scaled together with proper time; the local k2
restoration carries exp(-4lambda), the cubic carries exp(3lambda), and
the thermal circle is scaled as before. In particular the complete
singular term remains S(exp(lambda)Y). The full auxiliary remainder
is bounded by 4e19 q^2 L, yielding A_j=4e19*1000^j through nine.

The actual/cap row estimate needs energy orders through eleven. Nested
half-balls now have cutoff cost at most 520M. Its coefficient convolution

    2 sum_(r=1)^11 11!/(11-r)! [40/(1e5*18)]^r<.001

and cutoff contribution .0052 still fit the absorption budget. The
one-sided PDE therefore gives the same per-mode row estimate through
nine, without differentiating the even reflection nine times. Extend
the far angular moments to p=11 and the source formulas to j=9:

    H_p=6*18^p*1e-74[1+2/(172.8-p-3)],
    G_j=.00352*1e74/(j!100^j)
              *sum_(m=0)^j S2(j,m)(m+1)m!(1e4)^m H_(m+2),
    Q0=0, Q_j=1e16*(1e5)^j*.00352*2.1e-27 for j>=1,
    E_j=A_j+G_j+Q_j,
    O_j=sum_(k=0)^j |s1(j,k)|k!100^k E_k/(j!100^j).

S2 and s1 are the second and first Stirling numbers. The positive
dilation defect still has distance 3.9, not 24. Its telescoped inverse
bound 4j!1e5^j 21^(j-2)(20j+1)<j!1e7^j also holds at nine.
All moment denominators remain positive. Both max E_j and max O_j
are below 1e48, with their old entries through eight unchanged.

## 4. Cold source

Use the exact local, gluing and actual-domain decomposition of the
[completed cold C8 proof](seed-cold-actual-domain-eight-jets.md), including
its corrected Hermitian positive Clifford convention. Upgrade both the
global and local flat auxiliary steps to eta23. The local normal family,
Born and separated row bounds use the fixed cutoff's amplitude and
rotational symmetry; their parameter disks and bounds are unchanged.
The global taper's complex step derivative at delta=.05 is below 4.9,
so its value is below 1.245. With the old seed value .023 its scalar
cap remains below .03. Its real/complex gradients remain below 10/12:

    .02*.1/(.0075-.00021)+100*4*.023<10,
    1.245*.02*.1/(.0075-.00021)+100*4.9*.023<12.

The weighted outside bounds 100a exp(d),2000a exp(d) remain valid:
inside use the old seed height estimate, and on the taper a>=.01.
The local normalizer extension zeta remains the fixed C2 quintic.
Parameter derivatives act on F_b at fixed w; its spatial inverse is
not composed with zeta. The graph coefficients need at most ten cap
derivatives for nine parameter derivatives, available from C11.

The same graph bounds k=15000,u=20,w=25000,d=6/5,G=12 give

    D0=360000, D_j=360000*241^(j-1),
    T0=362400,
    T_j=sum_(b=1)^j (j-b+1)*.01^(j-b)D_b
                         +1200(j+1)*.01^j,
    H_j=4505.6*150^j*T_j.

Extend the physical connection ledger with Gamma_a=20a!300^a,
including a=0. The frame factor 6/5 is still sufficient because
exp(9*.02)=exp(.18)<1.2. For the local contribution use H0=1e7,
H_j=4e5*1.5^j with the same connection and frame factor.

For the actual-domain rows extend R_n=2-n/12 to n=11; R11=13/12>1.
The same induction U_n<=n!50^n U0 follows from its coefficient
[26+1.1 sum_(r=1)^(n-1)(2/75)^r]/(50n)<1. At observation order nine,
the point evaluation factor is 13+600*10+5000*10*11=556013<6e5.
Since 6e5*200*1.1/(1-2/75)<2e8, the previous row bound remains
2e8 j!(100M)^j exp(-MD/4). No global cap analyticity is used.
The cold extraction and entire angular tail therefore remain

    L_j=6*18^j*6e-8[1+2/(16.776-j-3)],
    H_j=1e19(j+1)*.00352*(1000/3)*L_j.

At j=9 the denominator is 4.776>0. Convert this term using Gamma0=0,
Gamma_a=2a!25^a. Add the three physical tensor bounds. Their maximum
through nine is below 1e48 a. This retains the bolt weight a(h), all
angular sectors, and the original product subtraction.

## 5. Tensor recurrence and original metric reconstruction

For every scalar ledger H, use the full ordered tensor recurrence

    B(0,k)=k!100^k H_k,
    B(r+1,k)=B(r,k+1)+r sum_(i=0)^k binom(k,i)Gamma_i B(r,k-i),
    J_j=B(j,0)/(j!100^j).

It retains nonsymmetric covariant derivatives and all connection
actions. A stated physical frame factor multiplies J only where needed.

Append only the new highest-order entries to the prior metric ledgers:
[middle](middle-metric-ten-jets.md), [cold](cold-metric-ten-jets.md),
and [outer](outer-metric-ten-jets.md). In their normalization, the
new source entry is 2e48 for middle/cold and 1e48 for outer; entries
through eight keep 2e44 and 1e44. All analytic seed coefficient bounds
hold at arbitrary order, so the same reconstruction formulas apply.
No additional source derivative is introduced by the regular bolt
integral identity. It gives v^(11)=G^(9)-integral_0^1 t^10 G^(9)(t rho)dt.

For the middle append P9,F9 using n=8 in its pressure/mass formulas,
then append U11,N11 using n=9 in its second-order formulas. The old
arrays through ten remain exact. The combined radial tangent reserve
is H=1e17. Through degree eleven the ideal of terms of degree at least
six still squares to zero. The tau inverse coefficient sum stays below
1500. Thus V_high<9001H, X_high<.0015+28.5015H and the finite logarithm
reserve is below 40H. Use radial reserve 4e18; Gamma_k=6k!200^k gives
physical reserve 1e20 through eleven.

For the cold append F9 using n=8, the existing pressure bounds through
eight, and the old source entries. Append U11,N11 using n=9 and the
new ninth source entry. The radial tangent reserve is H=1e17. The
regular bolt integral uses no tenth source derivative. In the finite
algebra V_high<985H since a/tau<82. With the same low .00072a,
X_high/a<.00072+(7.2+.00072*985)H<8H, giving finite radial reserve
3e18 a. Extend the radial-to-Cartesian coefficients

    C(1,1)=1, C(m+1,j)=C(m,j-1)+(4m-j)C(m,j),
    A0=A1=1, A_m=sum_(j=1)^m C(m,j)/(m-j)!

through eleven. Taylor subtraction at the center uses the even degree
ten polynomial for m=11, whose eleventh Cartesian derivative vanishes.
Then use Gamma_k=100k!20^k near the bolt and 15k!200^k away from it.
The resulting reserve is below 1e27 a, including the bolt by continuity.

For the outer use the previous coefficient convolutions and weighted
primitive I4. Append mu10=8*(4e-27)*1e48/(100*10), then compute F
through nine, the small j through ten, its primitive and proper W
through eleven. Keep the old proper V entries through ten and append
V11 with n=9 in its second-order equation, source 1e48. Transform once
back to the original common coordinate using the old exact integral
identities. The combined ordinary remainder fits H=1e16. The low
finite parameter is still bounded by epsilon(Y)=2.1e-7Y^3; high terms
start at degree six and their squares vanish through eleven. A bound
10*(2.1e-7)^2/100+200*epsilon(.1)*H<1e-6 H covers the nonlinear
logarithm remainder. Use finite reserve 3H and the unchanged weighted
connection Gamma_k=(9/400)k!2^k. This gives physical reserve 1e17.
The correlated cubic and its analytic remainder are retained.

## 6. Verification and consequence

The checker freezes all inherited inputs, checks the new cutoff jets,
operator margins, complete angular/source ledgers, exact appended metric
entries and finite/tensor conversions. Independent noncommuting matrix
inverse derivatives test ninth-order insertion counting. Raw Leibniz
products test the new ODE entries. Manufactured finite logarithm contours
at 60/90 digits test eleventh-order nonlinear coefficient extraction;
they do not sample the actual charged source. Prior exact heat, locality,
curvature and coordinate identities are retained, not replaced by a
numerical model. These are supporting checks for an ordinary argument,
not a Lean or machine-verified continuum proof.

The original metric C11 input for the direct off-shell row is now
quantified. The existing endpoint inverse/pullback result is still C10
for the metric, with its stated C11 auxiliary map. This note assigns
neither C12 regularity to that map nor an eleventh-order pullback to it.
The actual finite-metric quantum H4 bound, local H6 source, mixed/lift
and higher TT inverse constants and positive transfer are not proved.
No nonlinear solution, signed compatibility or stability follows yet.

Remaining to-do list: finite-metric quantum H4, local/Schur H6 and TT
transfer estimates; cancellation, signed compatibility and nonlinear
existence; physical validity, state/boundary, stability, matching,
charge/return, opening, transport, preparation, separate creation and audit.
