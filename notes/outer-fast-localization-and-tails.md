# Outer fast localization with the actual nonlocal tails

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Global-input collar
estimate; the mass column and global matching remain open.**

The complete supported-input comparison now yields a quantitative
collar graph estimate for a global FAST field. Nonlocal input from
the rest of the disk enters through a global value norm. The only
remaining derivative term is small and supported in a specified
transition annulus. A nonzero mass cubic is not put in this fast space.

## 1. Fixed cutoffs and statement

Retain the physical action, state, clock, paired domain and all
contacts of the [complete outer comparison](outer-complete-disk-normal-comparison.md).
Use its quantified similarity vhat=M_epsilon v, Ahat=M_epsilon A M_epsilon^-1,
with theta=S+T/3, b=2T and the dual physical output. Put
fhat=lambda^-2 Ahat vhat. The input X and output Z are the same
Hilbert H4/H4log and L2 pairs at weight delta=7/2 in x=log Y.

Choose 0<epsilon<=1e-160 and r=epsilon/32. Let

    P(t)=462t^6-1980t^7+3465t^8-3080t^9+1386t^10-252t^11,
    c_a(Y)=1 (Y<=a), 1-P(Y/a-1) (a<Y<2a), 0 (Y>=2a),
    chi=c_r, eta=c_(4r), sigma=(1-c_(r/4))c_(8r).

Thus sigma is one on [r/2,8r] and supported in [r/4,16r].
Extend these functions through the actual disk as specified by
their constant/zero regions. Define the global canonical value norm

    V0(vhat)^2=A0^-1 integral_gbar (abs(theta_hat)^2+abs(b_hat)^2),
    A0=beta a(0), a(0)=sqrt(C lambda).

For stationary fields this is integral j(Y)abs(vhat(Y))^2 dY,
j=a/a(0), on the WHOLE physical half, including the regular bolt.
It is not the weighted outer value norm. The theorem is

    ||chi vhat||_X <=2e20 ||chi fhat||_Z
                    +1e-40 ||sigma vhat||_X
                    +1e2000 r^-4 V0(vhat).                    (1)

Initially take smooth fast fields in the retained domain. The
estimate extends to their graph closure when the displayed norms
and the local source converge. It is an a priori estimate for a
global input, not a collar inverse or a new boundary condition.
The C5 cutoffs have bounded piecewise sixth derivatives; smooth
approximation preserves the strict reserves used below.

The physical mass-extended solution is vhat=vhat_fast+m ehat_m.
Applying (1) to vhat_fast requires its actual source
lambda^-2 Ahat vhat_fast=fhat-m lambda^-2 Ahat ehat_m.
No bound for this full mass-column source is supplied here.
The previously estimated global/constant part alone is insufficient.

## 2. The exact normal kernel and its two kinds of tail

Let N0 be the dimensionless complete normal operator. Its only
nonlocal entry is the b,b entry. From the flat massive continuum
kernel in the [flat tower proof](outer-flat-tower-normal-remainder.md),
integration of the time difference gives, away from the radial diagonal,

    K_N(Y,Z)=-Y^4/(40pi^2)
                      [abs(Y-Z)^-5+(Y+Z)^-5], Y,Z>0.         (2)

Indeed integral_R (s^2+t^2)^-3 dt=3pi/(8abs(s)^5), and
the continuum kernel is -16C/(15pi^2 d^6), with physical b
factor Y^4/(16pi C). All normal Weyl and finite terms are local.
The even image has no extra half. Equation (2) is used only at
separated points, with the retained diagonal prescription unchanged.

After weight conjugation and Z=Y exp(t), the absolute log kernel is

    k_delta(t)=exp((1+delta)t)/(40pi^2)
                    [abs(1-exp(t))^-5+(1+exp(t))^-5].        (3)

For abs(t)<1/10 it is below 1e6 abs(t)^-5. For abs(t)>=1,
it is below 1e3 exp(-t/2) on the positive side and
1e3 exp(9t/2) on the negative side. On the remaining compact
interval use 1e6 abs(t)^-5. Consequently, for 0<R<1/10,

    integral_(abs(t)>=R) k_delta(t)dt <1e10 R^-4.           (4)

There is also a different, unweighted tail estimate. Suppose
chi(Z)-chi(Y) is nonzero and abs(log(Z/Y))>=R. At least one
of Y,Z exceeds r, and their separation is at least
r(1-exp(-R))>=rR/2. Their sum is at least r. Direct integration
of the symmetric ABSOLUTE bare kernel in (2) after removing
Y^4 gives both dY endpoint rows below 1e10(rR)^-4.
At delta=7/2 there is the exact output identity

    ||Y^-delta Y^4 F||_L2(dY/Y)=||F||_L2(dY).              (5)

Hence this tail maps ordinary input L2(dY) to the weighted
output norm. On the supported outer input j>1/2, so its input
norm is below sqrt(2) V0. This is why the tail in (1) does not
require an unknown weighted value norm extending to infinity.

If the input is chi vhat and Y>=epsilon, all local normal
terms vanish. The separation from support Y<=2r and the same
row argument give

    ||1_(Y>=epsilon) N0(chi vhat)||_Z
                                  <=1e10 r^-4 V0(vhat).    (6)

This controls the normal output on the rest of the full log line.
It is not discarded when applying the normal inverse.

## 3. A complete normal commutator with localized derivatives

Write W=Y^-delta eta vhat, and conjugate N0 to the x-line.
For the cutoff chi(e^x), its first derivative is below 22 and
its derivatives through six are below 1e14. Chain rule gives
the latter from the beta derivative bounds and the finite
Stirling sums sum_j Stirling(k,j)2^j d_j, k<=6, where
d=(1,11,220,3960,63360,887040,10644480).

The separated commutator kernel is exactly
(chi(x+t)-chi(x)) times (3). On abs(t)<R=1e-100 it is below
1e8 abs(t)^-4. Subtract the input Taylor polynomial through
degree three. Integral Taylor remainder and Minkowski give

    ||T_R W||_2 <=(2e8 R/24)||W''''||_2.

The output and every Taylor segment of a nonzero near term lie
where sigma=1. Thus W may be replaced here by Y^-delta sigma vhat,
and the near contribution is below 1e-90||sigma vhat||_X.

Keep the entire complementary kernel as the raw integral B_R;
do not subtract its row sum. All local commutator contacts have
input order at most three. To see this with the fixed subtraction,
multiplication by chi(x+t)-chi(x) removes the fourth input jet
from the original fourth-order finite part. Its ordinary fourth
restoration also commutes to order three. If input jets through
three vanish, the original distribution is the convergent integral.
This fixes the contact order without differentiating a kernel bound.

Bound every contact by the four local probes

    q_j(x+t)=c_(1/8)(abs(t)) t^j/j!, j=0,1,2,3,

in each of the two input components. Their jets form the identity
and T_R q_j=0. The q_j derivatives through six are below 1e16.
For g_j(t)=(chi(x+t)-chi(x))q_j(x+t), the H6 norm is below 1e35.
The accepted forward normal norm is below 1e38 from X to Z.
Commuting one x derivative and using H5log<=2H6, then the
one-dimensional H1 point bound, gives

    abs([N0,chi]q_j(x))=abs(N0 g_j(x))<2e73.

The raw tail on the probe is at most 1e10 R^-4=1e410 by (4).
Every contact coefficient is therefore below 1e420. Outside
the R-neighborhood of the transition, the near commutator and
all contacts vanish: chi is constant there and the entire action
is the raw separated tail. Thus contact derivatives also use sigma.
Summing the sixteen possible row/input/jet coefficients and using
the unweighted tail estimate from section 2 proves

    ||[N0,chi](eta vhat)||_Z
       <=1e-90 ||sigma vhat||_X
          +1e423 ||Y^-delta sigma vhat||_(H3 pair)
          +1e430 r^-4 V0(vhat).                         (7)

All normal finite contacts are included in this point-probe argument.
The commutator is not treated as an absolute diagonal kernel.

## 4. The actual response from the rest of the disk

On supp chi, an input (1-eta)vhat is separated by at least 2r.
The local classical, Weyl, LLL and massive derivative contacts
have zero action there. For the remaining actual massive kernel,
retain the geodesic envelope 72C/d^6, the local comparison
1e8C/d^4 and the actual gluing/far endpoint rows 3e24 C from
[the curved response](stationary-curved-quantum-tt.md).

All geodesic flows reaching that input have length at least 2r.
Integration of the full flow, including returns, gives the bound
2pi(18C)/(2r)^4. The local comparison contributes at most
pi(1e8C)/(2r)^2, with harmless enlargement for the second row.
The actual global term stays separate. The combined rows are
below 1e29 C r^-4. Liouville invariance supplies the reverse
flow row; the actual paired kernel supplies its reverse row.
No global injective exponential map is assumed.

Use the physical massive factor rho Y^4/(16pi C), (5), the
bounded pointwise M_epsilon and inverse, and j>1/2 on output.
The normalized full-disk value norm cancels the fixed A0. Then

    ||chi lambda^-2 Ahat[(1-eta)vhat]||_Z
                                <=1e31 r^-4 V0(vhat).     (8)

This includes input near the regular bolt and the thermal/domain
response. It is a large value term, not a small collar error.

## 5. Normal inversion and absorption

On Y<epsilon the exact identity is

    N0(chi vhat)=chi fhat
       -chi lambda^-2 Ahat[(1-eta)vhat]
       +[N0,chi](eta vhat)
       +chi(N0-lambda^-2 Ahat)(eta vhat).                 (9)

The supported-input comparison applies to eta vhat, whose support
is Y<=8r=epsilon/4. It bounds the last term by
1e67 epsilon ||eta vhat||_X. Let zeta=eta-chi. Its support is
where sigma=1. Its C5 log jets are below 2e14; product estimates
give H3/H5 multiplication norm below 2e17. The explicit
logarithmic norm transfer costs less than 1e6, hence

    ||eta vhat||_X<=||chi vhat||_X+1e24||sigma vhat||_X.   (10)

For the H3 pair in (7), the retained Fourier interpolation gives

    U3 <=tau X+(2/tau)^3 U0, tau=1e-500.

Here U0 is only on supp sigma, where Y>=r/4 and j>1/2.
Thus U0<=400 r^-4 V0. This uses the actual inner support r/4,
not the narrower region where sigma is one.

Apply the normal inverse bound 8.7e19 to (9), also adding (6).
The self coefficient from (10) is below 8.7e-74, so it can
be absorbed with factor below two. The transition coefficient
after this absorption is bounded by

    2(8.7e19)[1e67(1e-160)(1e24)+1e-90+1e423(1e-500)]
                                                          <1e-40.

The value coefficient is below

    2(8.7e19)[1e31+1e10+1e430
                         +400(1e423)(2e500)^3] <1e2000.

The forcing coefficient is below 2e20. This proves (1).
Nothing in the absorption claims the transition annulus is
already controlled, or that N0^-1 preserves collar support.

## 6. Verification contract and remaining work

Freeze the proof and checker before execution. Pin the complete
outer comparison, forward/inverse normal, compact logarithmic
transfer and actual global-tail receipts. Check the radial kernel
constant, two weighted decays, weight/measure identity, cutoff
jets and support margins, exact probe jets, H6 point reserves,
finite contact count, both tail types, multiplication/interpolation
budgets and final absorption. Preserve every scientific run and
failed source snapshot. Finite checks support the ordinary kernel
and distribution argument; they are not formal verification.

The next matching step must control the specified annulus and
the full actual mass-column source before combining this with
the compact/bolt graph. A qualitative mass-domain theorem does
not by itself provide that numerical source norm. The global
inverse/residual product and nonlinear requirements remain open.

Remaining to-do list: quantitative full mass-column image;
transition/intervening annulus and compact/bolt matching; higher
York reconstruction and inverse/residual composition; nonlinear
family and signed compatibility; exterior/state, charge, return,
opening, stability, transport, preparation and creation.
