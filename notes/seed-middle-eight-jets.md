# Complete middle seed-source derivatives through order eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original seed, state, full angular tower, physical action and
common subtraction. In the middle cover s>=.001, Y>=.05, the result is

    |nabla^j Theta_seed|/(j!100^j q^2)<1e44,
    0<=j<=8.                                      (1)

The individual bounds below are much smaller at the lower orders.
Retain the previously repaired C3 allowance 1e26; the withdrawn 1e22
claim is not reinstated. This gives the quantitative middle seed input
identified by [the source-to-metric review](repaired-source-metric-transfer.md).
Cold and outer derivatives through eight, the ten-derivative metric
ledger and the actual finite-metric H6 Schur source remain open.
No cancellation, nonlinear solution or stability follows from (1).

## 1. Upgrade the auxiliary cutoff without changing the source

The [common-cap construction](seed-middle-cap-transfer.md) supplied
true cap derivatives through three. Its C4 cutoff is not sufficient
for eight cap derivatives. Replace only that auxiliary cutoff by

    eta(t)=sum_(i=0)^10 binom(21,i)t^i(1-t)^(21-i), 0<=t<=1,
    eta=1 to the left, eta=0 to the right,
    chi(r)=eta(128r-3), zeta(r)=eta(64r-1).         (2)

The degree-21 step is C10, takes values in [0,1], and has zero
derivatives through ten at both endpoints. Its derivative is

    eta'(t)=-c t^10(1-t)^10, c=21!/(10!)^2=3879876.

It preserves the previous support radii, observation neighborhood and
local Born argument. It does not replace the physical exterior or
the quantum state. At the base center the fixed-cap gluing estimate
uses this same new cap. Its row proof depends on the local seed and
the low cap coefficients checked below, not on the old polynomial.

Use the existing analytic normal chart R*=1/30, |sigma|<.002,
|partial sigma|<.067 and |partial^2 sigma|<3. In the complex physical
center bidisk |b1|,|b2|<=rho=1e-5, keep the exact normalizer

    F_b(0)=beta, beta=b1+i b2, betabar=b1-i b2,
    F_b'=exp[sigma(beta,betabar)-2sigma(F_b,betabar)].

Its earlier bounds remain

    |F_b-id|<1.5rho, |DF_b-I|<.2rho, |D^2F_b|<10rho,
    |sigma_b-sigma|<.3rho, |D(sigma_b-sigma)|<10rho. (3)

All spatial arguments displaced by at most 2rho still lie in R*.
The common physical family is again

    G_b(w)=w+zeta(|w|)(F_b(w)-w),
    u_b(w)=chi(|w|)sigma_b(w),
    cap_b=(G_b)_*[exp(2u_b)|dw|^2].                (4)

Use the same extended spin rotation Q_b. For real small b this is
a global coordinate identification, equal to the identity outside
R0=1/32. The cap equals the seed in a common neighborhood of the
observation. All its b derivatives are supported outside r=1/64.

## 2. Uniform analytic branches bound exact finite jets

No artificial cutoff is globally analytic. At a fixed real spatial
point x in a transition, use the analytic polynomial branch in (2),
with the local square-root branch r(w)=(w1^2+w2^2)^(1/2). Outside
the transition use the constant branch. At an endpoint either
branch has the same jets through ten. For an interior transition
point the analytic branch matches the true cutoff on some real
neighborhood, even if that neighborhood is smaller than 2rho.
Its continuation is only a device to bound the jets at that point.

For |w-x|<=d=2rho and r_x>=r_min>=1/64, set
q_r=2d/r_min+(d/r_min)^2. Since sqrt(1-q_r)>.9987,

    |r(w)-r_x|<1.002d, |Dr(w)|<1.003,
    |D^2r(w)|<2.5/r_min.                          (5)

These follow from the square-root identity, Dr=w/r and
D^2r=I/r-ww^T/r^3, in complex Hilbert tensor norms. In particular
the step argument changes by less than delta=.003. If its real
base point is t0 in [0,1], then, with B=1/4+2delta^2,

    |t(1-t)|<=B,
    |eta'|<=c B^10<3.8.

For the second derivative write a=|t0-1/2| and complete the square
B-(a-delta)^2. Maximizing v(B-v^2)^9 for v>=0 gives

    |eta''|<=20c(18B/19)^9 sqrt(B/19)+40c delta B^9<40.          (6)

For a<delta the second term alone suffices. The rational enclosure
sqrt(B/19)<1/8 is enough to verify the displayed bound.
The chain rule and (5)-(6) therefore give the branch budgets

| Cutoff | value | first tensor | second tensor |
| --- | ---: | ---: | ---: |
| chi | 1.010 | 490 | 720000 |
| zeta | 1.005 | 245 | 210000 |

The value estimate integrates eta' from the real value in [0,1].
All higher finite jets are bounded by Cauchy on these branches;
large global higher derivatives of the piecewise cutoff are not
silently assigned the low-order bounds in this table.

For fixed x, construct the local inverse G_b^-1(x) on the 2rho
ball using these branches. Equations (3)-(4) give

    |DG_b-I|<370rho, |D^2G_b|<316000rho,
    |DG_b^-1|<1.004, |G_b^-1(x)-x|<1.6rho.         (7)

The map displacement is at most 1.005(1.5)rho and its contraction
constant is below .0037. Thus this local analytic inverse exists.
There is still no claim of a global complex cap inverse.

The branch coefficient jets through b order eight, including the
one spatial derivative in the operator, are exactly the true jets.
The scalar multiplier uses cutoff jets through nine and the
twice-spatially-differentiated map uses jets through ten. C10
therefore gives a C8 common-operator family in the first-order graph
norm, including at transition endpoints. The piecewise-polynomial
remainders start at order eleven; after the required differentiations
they cannot affect these jets. The same observation applies to the
extended spin frame. No special differentiability only at b=0 is needed.

## 3. The full half-density vertex bound through eight

Keep the exact common-space operator, including its density term:

    A_b=P_b^i partial_i+(partial_i P_b^i)/2,
    P_b^i=Q_b exp(-u_b(w)) gamma^a DG_b^i_a(w) Q_b^-1,
    w=G_b^-1(x).                                   (8)

This is the coordinate/spin transport derived in the common-cap
proof. The Spin(2) anticommutator vanishes; its absence does not
remove the half-divergence term in (8).

The branch table gives |Du_0|<1.1 and |D^2u_0|<1600, and before
composition |u_b-u_0|<.303rho, |D(u_b-u_0)|<160rho.
After inverse composition the scalar multiplier difference is
below 2.1rho and its gradient difference below 3200rho.
The spin factor satisfies |Q_b-I|<.21rho and |D_xQ_b|<65rho.
Using all inverse, density and rotation products in (8) gives

    |P_b-P_0|<800rho, |div P_b-div P_0|<1.2e6rho.  (9)

For the largest divergence term one may use
2(1.004)^5(316000)rho, allowing the scalar, inverse and frame
factors. The scalar-gradient and remaining products each fit
10000rho. These sums are below the reserve in (9).

In the graph norm ||v||_EM=||partial v||_2+M||v||_2, M>=18,
the analytic branch operator difference consequently has bound

    800rho+1.2e6rho/(2*18)<A=7/20.

Cauchy in the two center variables, including their Hilbert tensor
factor, gives true cap vertices

    ||partial_b^j A_b|_0 v||_2<=K_j||v||_EM,
    K_j=A j! k^j, k=150000, 1<=j<=8.              (10)

Here sqrt(2)/rho<k. Keeping A in every vertex is useful; it is
not raised only once in a product that contains several vertices.

The real base cap has |u_0|<.002, |Du_0|<1.1. The same weighted
inverse calculation gives ||S||<=1/M, ||S||_(L2->EM)<4, weighted
inverse <2.01/M, weighted flat gradient <4.1, and weighted inverse
followed by divergence <5.2. The slightly larger cap gradient
fits these old reserves. The interior row recurrence from
[the coercive estimate](seed-middle-coercive-jets.md) is valid through
nine. Thus for all j<=8 the separated plain/input-gradient rows
have bounds c_j exp(-M/176) and 40M c_j exp(-M/176), where

    c0=30000, c_j=520000(j+3)[600(j+3)]^j, j>=1.    (11)

The cap support and observation ball leave the same separation.
No derivative of a remote cutoff is assigned to a pointwise row.

## 4. All inverse insertions and the complete angular sum

The rth inverse derivative is the sum over every ordered composition
i1+...+il=r, ia>=1, with coefficient r!/(i1!...il!) and sign
(-1)^l. It contains l vertices and l-1 intermediate graph inverses.
Using (10) and the graph constant four gives exactly the majorant

    V_r=r! k^r sum_(l=1)^r binom(r-1,l-1) A^l 4^(l-1)
       =A r! k^r (1+4A)^(r-1).                   (12)

The final vertex acting on an endpoint row costs 41M; the trace
contributes M and rank two contributes two. Hence the factor is
82M^2, independently of the number of intermediate insertions.
The already proved complete angular moment remains

    sum d_n M_n^2 exp(-M_n/88)<25 C 88^4.          (13)

Set P_j=sum_(a=0)^j binom(j,a)c_a c_(j-a). With the exact split
Theta=U_mov+G_fix-Ccap and Ccap(x)=U(x,x)-U(x,0), define

    kappa_j=82*25*88^4*(31/2500000)/(j!100^j)
             *sum_(r=1)^j binom(j,r)V_r P_(j-r), j=1,...,8.       (14)

This bounds every mixed cap/observation term in the missing
correction. Local counterterms have zero cap derivative on the
observation neighborhood. The exponential angular majorant
justifies these differentiations without truncating the tower.

For the fixed-cap gluing derivatives retain

    G_j=8*256^2*88^2*(31/2500000) P_j/(j!100^j).

For U_mov, the normal factor in its own chart is analytic in the
physical center bidisk, below .003 with quadratic reserve eight.
The common-cap Born and local restoration bound is unchanged:
L_j=3e13*1500^j after j!100^j normalization, for j>=1.
The cutoff is fixed in that chart, so its spatial nonanalyticity
does not obstruct center holomorphy. The a0 term is constant.
Consequently the complete Cartesian bounds are

    H_0=1e14, H_j=G_j+L_j+kappa_j, 1<=j<=8.        (15)

No old moving/fixed split is reused without its correction.

## 5. Conversion to full covariant derivatives

At each base normal center sigma=partial sigma=0 and g=I. On the
complex Cartesian bidisk |x1|,|x2|<=1/64, the corresponding z,w
remain within R*. The conformal Christoffel tensor has Hilbert norm
2|partial sigma|<.134. Cauchy therefore gives at the center

    Gamma_0=0, |partial^a Gamma|<=gamma_a,
    gamma_a=(67/500)a!100^a, a>=1,                (16)

since 64sqrt(2)<100. This uses the analytic seed, not a derivative
of an unsigned real bound. It controls all connection indices.

To retain nonsymmetric higher covariant tensors as well, define
B_(0,k)=k!100^k H_k for 0<=k<=8. Recursively, for r+k+1<=8, use

    B_(r+1,k)=B_(r,k+1)
       +r sum_(a=1)^k binom(k,a) gamma_a B_(r,k-a).              (17)

The r terms differentiate the r covariant tensor indices; the
zeroth connection term vanishes at the center. Tensor contraction
is bounded by the Hilbert product norm. Induction bounds
|partial^k nabla^r Theta|/q^2 by B_(r,k).
Thus the physical normalized entries are

    J_j=B_(j,0)/(j!100^j), j=0,...,8.             (18)

This is a full tensor estimate, not merely the diagonal geodesic
derivative of its symmetrization. Exact rational aggregation of
(12)-(18) gives max J_j<1e44, proving (1). Proper radial jets are
contractions of these tensors, so the reconstruction can use them.

## Verification and limits

Freeze this proof, checker and review before execution. Verify the
cutoff endpoint jets, complex radial/step and coefficient budgets,
all ordered insertion counts, full angular ledger, and complete
coordinate/covariant aggregation. Check the connection recurrence
against explicit independent low-order covariant derivatives.
The local branch and continuum row arguments remain ordinary proofs.
No actual signed charged kernel or finite-metric residual is evaluated.

This is one regional higher-source input, not a global ten-derivative
metric bound or the actual N6(p). It preserves all prior low-order
contracts. No Lean run is selected and no physical gate is completed.

Remaining to-do list: cold and outer seed jets through eight; the
quantitative ten-derivative metric and finite-metric quantum/Schur
source estimates; positive transfer constants or a sharper actual
transfer; cancellation, signed compatibility, nonlinear existence;
physical validity, state/boundary, stability, matching, charge/return,
opening, transport, preparation, separate creation and audit.
