# Common-regulator strong-limit test in the actual coupled graph

2026-10-25 HST, E003/S169. Ordinary argument frozen before computation.
The full stable semiclassical wormhole objective remains active and incomplete.
The closed feasibility audit and its original targets remain unchanged.

OUTCOME: the complete coupled estimate is UNRESOLVED. This attempt derives
the entire positive-regulator shape Hessian, including its second vertex,
and tests a specified way of taking its strong limit. A fixed positive
proper-time cutoff with only local heat counterterms cannot approximate
the physical Hessian uniformly in the original normal graph: the error
norm is greater than .9999997 for EVERY positive cutoff. This is a
failure of that approximation, not an excessive bound on the physical
response, not the actual mass coefficient, and not a construction obstacle.

The qualitative warning against raw-cutoff strong limits already appears
in current-original-cutoff-coupled-estimator.md. The new results here are
the combined whole-spectrum formula and its quantitative original-norm
failure test. The angular-resummed common-zeta identity in
current-angular-resummed-quantum-operator.md already retains the missing
nonlocal ultraviolet term. That identity is NOT refuted. This test does
not complete the user's requested coupled numerical feasibility milestone.

## 1. First unsupported step and exact estimator being tested

The existing admitted-family weak C2 result concerns a scalar action with
two admitted inputs. The missing implication is a bound for its variational
gradient in the ORIGINAL two-row strong near/transition source norm and
phase/interface border, uniformly on the requisite nonlinear neighborhood.
One cannot obtain that implication just by sending a scalar proper-time
cutoff to zero under a strong operator norm.

Precisely, let H_Q be the canonical massive quantum shape Hessian, divided
by the retained physical action factor W=-N pi B^2/(G gamma). Let
H_Q,e,ren be the Hessian of the SAME full angular/thermal proper-time
integral over t>=e, e>0, plus its local divergent counterterms and the
finite local terms which match the selected common zeta scheme. All
other allocations of the full two-row operator are identical on both
sides of this comparison. The test asks whether H_Q,e,ren converges to
H_Q in the original normal-preconditioned fast graph operator norm.

This is a candidate evaluation route, not a statement that an earlier
accepted implementation used it. In particular we do not replace the
zeta_A'(-1) term in the exact angular regrouping by zero. No independent
renormalization of angular levels or thermal frequencies is permitted.

## 2. Actual paired realization and complete positive-cutoff formula

Use the fixed original paired first-order realization. As proved in
current-rough-paired-variation.md, Dsharp=J D has Dsharp*=-Dsharp.
Thus K=i Dsharp is self-adjoint and L=K^2 is the actual nonnegative
orbit square. For a real identical-on-both-copies shape b,

    K_b=exp(-b/2) K exp(-b/2), Dom K_b=Dom K.

Bounded W1,infinity multipliers preserve the first-order trace projector
and graph norm. The square is interpreted through its closed form on
Dom K; its composed operator domain is NOT asserted fixed. The regular
disk bolt, fixed period and spin structure remain those of the original
problem. The admitted mass cubic plus fixed-finite-log weighted C5 phase
is unchanged. This argument does not smooth the reflected seam.

Put q=2000^4, C=1/u, M_n^2=n(n+q)/C, d_n=2(q+2n), and

    Theta(t)=sum_(n>=1) d_n exp(-t M_n^2),
    Gamma_e(b)=1/2 int_e^infinity Theta(t) Tr_ph exp(-t K_b^2) dt/t,
    f_e(a)=1/2 int_e^infinity Theta(t) exp(-t a^2) dt/t,
    R_e(lambda)=int_e^infinity Theta(t) exp(-t lambda) dt,
    F_e(lambda)=lambda R_e(lambda).

Tr_ph is the original physical-half trace. Equivalently take the full
paired trace and its fixed physical multiplicity; scalar variations
preserve the pairing. Keep that factor through to W. All thermal labels
and finite-domain returns are in K. At e>0 every displayed trace/integral
converges, since M_1^2>324 and the orbit has the established heat bound.

For commuting real multiplication directions v,w the exact vertices are

    K_v=-(vK+Kv)/2,
    K_vw=(vwK+Kvw+vKw+wKv)/4.

Duhamel differentiation on the common first-order form, followed by
spectral calculus, gives the following weight for matrix elements v_ab
and w_ba in the COMPLETE K spectrum (a,b are spectral values):

    h_e(a,b)=(a+b)^2 [f_e'(a)-f_e'(b)]/[4(a-b)]
                       +(a+b)[f_e'(a)+f_e'(b)]/4
             =-(a+b)^2 [F_e(a^2)-F_e(b^2)]/[2(a^2-b^2)].       (1)

The first line includes the bubble AND the second metric vertex. The
diagonal is -2a^2 F_e'(a^2); the value for a=-b is zero, including
the continuous value at a=b=0. These are spectral identities, not a
selection or numerical truncation of physical modes. One can derive
them first on smooth compact directions and then use the bound below;
the admitted-family positive-time form argument retains the original
domain for rough directions. No spectral differentiability of individual
multiple eigenvectors is needed: divided differences include degeneracy.

An integration by parts exposes another contact that must be retained:

    F_e(lambda)=Theta(e)exp(-e lambda)
                         +int_e^infinity Theta'(t)exp(-t lambda)dt,
    d_t(lambda,mu)=int_0^1 exp[-t((1-a)mu+a lambda)]da,
    [F_e(lambda)-F_e(mu)]/(lambda-mu)
      =-e Theta(e)d_e(lambda,mu)
                         +int_e^infinity t[-Theta'(t)]d_t(lambda,mu)dt. (2)

Dropping the regulator endpoint in (2) changes (1). It is not a harmless
coincidence convention. The angular sum in Theta is complete in both
terms; (2) is taken before regulator removal.

## 3. What the positive-cutoff operator can control

Define the nonnegative function

    G_e(lambda)=lambda R_e(lambda)
                   +lambda^2 int_e^infinity t Theta(t)exp(-t lambda)dt.

For ALL real a,b, writing lambda=a^2 and mu=b^2,

    |h_e(a,b)|<=32[G_e(lambda)+G_e(mu)].                    (3)

Here is an elementary proof, including both spectral regimes. If
lambda>=4mu, (sqrt(lambda)+sqrt(mu))^2/(lambda-mu)<=3, so (1) is at
most (3/2)(F_e(lambda)+F_e(mu)); interchange the variables if needed.
If mu<=lambda<=4mu, the mean value integral for F_e and

    |F_e'(nu)|<=int_e^infinity Theta(t)(1+t nu)exp(-t nu)dt

give |h_e|<=8mu R_e(mu)+32mu^2 int t Theta(t)exp(-t mu)dt.
The case mu=0 follows from the separated regime or continuity. This
proves (3) without a spectral cutoff or an unknown comparison constant.

At each FIXED e>0, G_e(L) is positive trace class and has a bounded
diagonal on every compact physical interval I away from the ends.
Indeed its spectral decay is exponential times a polynomial; local
positive-time elliptic smoothing applies on I, and the original heat
trace controls the full spectrum. Thus for real v supported in I,

    |Gamma_e''[v,v]|<=64 Tr_ph(G_e(L)v^2)
      <=64 sup_I tr G_e(L)(x,x) ||v||_L2(I,dvol)^2.          (4)

Polarization supplies a localized L2 operator. Its constant may diverge
as e decreases; (4) is NOT a uniform regulator estimate. In particular
the contact cancellation makes the fixed-cutoff shape Hessian order zero
on these localized inputs. No unknown numerical value of that constant
is needed for the high-frequency test, since e is fixed first.

## 4. The actual ultraviolet term and its original normalization

The retained interior transfer argument in
stationary-shape-logarithmic-response.md uses the local Dirac symbol from
[Moller, Theorem 1 and Lemma 6.1](https://arxiv.org/html/0705.3857v3).
That theorem assumes a CLOSED SMOOTH spin manifold and stable kernel.
It is not a global theorem for our rough paired seam. On a compact
physical interval the actual seed is smooth and positive; nested
cutoffs identify its simultaneous short-time singularity with that local
symbol. The fixed twisting is lower order. Heat paths reaching the
boundary, thermal images and regions with only one short heat time are
retained in the lower-order term; they cannot change the fourth-order
logarithm. Fixed spin structure and the massive gap prevent a kernel
crossing in this use. No global closed-manifold factorization or stability
claim is imported. The transfer requires no smooth reflected seam.

In the actual affine chart

    g=exp(2phi)[(B^2/C)exp(2T)gamma0+B^2 dOmega^2],
    phi=-log z+nu,

the metric tangent is diag(2nu+2T,2nu+2T,2nu,2nu).
For radial covectors its transverse traces satisfy

    T2=12nu^2+8nu T+4T^2, T1=6nu+2T,
    3T2-T1^2=8T^2.

The rank-four Dirac logarithmic Hessian, including the fermionic 1/2,
is p_phys^4 log(p_phys) T^2/(120pi^2) per physical volume. Integrate
over time and sphere, divide by the ORIGINAL negative W, and use
p_phys^2=z^2 C Q xi_z^2/B^2. The strong dz shape density has symbol

    H_Q,TT(z,xi_z)=-G C Q(z)^2 xi_z^4 log|xi_z|/(15pi B^2)
                          + ordinary order four and lower log terms. (5)

Let x=log z, D=partial_x, sigma=4.005. The original normalized shape
row is z^-sigma[(D-4)z^4 H_T-(D-6)z^4 H_nu/3], with the output cutoff
one throughout I. Canonical Gamma has no nu variation. On T=z^sigma w,
its leading multiplier is therefore

    -i K0 (Q/beta)^2 k^5 log|k|,
    K0=G C beta^2/(15pi B^2)=G lambda beta/(15pi)>0.        (6)

The complete normal shape multiplier B(sigma+ik) has leading term
-i K0 k^5 log|k|. This is precisely its retained normalization.

For fixed e the local heat counterterms and scheme restorations have
differential Hessians of order at most four on I. Their coefficients may
depend on e and include log e. They contain no log|k|. This follows also
from the local four-dimensional divergent coefficients (see
[Vassilevich, equations 1.20-1.21](https://arxiv.org/html/hep-th/0306138v3));
the actual paired heat proof, not a smooth-boundary citation, supplies
our common regulator. Boundary counterterms vanish on the compact test
directions. All local Weyl coefficients, including log z, are bounded
on I and are retained once. The separate LLL term is lower order. The
other physical local allocations are identical in the compared operators.

## 5. Executed necessary test on the actual original fast graph

Choose I=(.002,.003) in z and any nonzero real smooth psi compactly
supported in its x image. Consider w_k=psi exp(ikx), T_k=z^sigma w_k,
nu_k=-T_k/3. These directions are identical on the two physical copies
and vanish near BOTH the seam and bolt. Every asymptotic phase coefficient
is zero. They lie IN the original paired mass/log domain; the domain is
not replaced by their subspace. A bound for the full operator must hold
on this subspace. Real sine/cosine packets give the same lower bound
(or use the isometric complexification of this real Hilbert operator).

Define their exact full-line normal data by

    g_k=B_sigma w_k,
    f_k=Jc_sigma [2/(3(D+sigma+1))]w_k.                    (7)

Because a(s)=(s+3)/(3(s+1)), these are exactly the original inverse
coordinates theta=nu+a(D)T. The denominator in f cancels the factor
s+1 of Jc, leaving order three. Thus

    max(||f_k||H3,||g_k||H3)
                 =K0 |k|^8 log|k| ||psi||L2 (1+o(1)).     (8)

The data are allowed full comparison-line inputs; the resulting physical
fields have the compact support just specified. Normalize (7) by (8).
Their physical C5 size tends to zero as O(|k|^-3/log|k|), with fixed
finite seed constants. Hence the test persists in arbitrarily small
admitted families. It does not assume the unknown physical X embedding.

Let R_e be the ORIGINAL two-row graph of H_Q-H_Q,e,ren composed with
the normal inverse. Both other original equations/allocations are kept
the same, so their difference is zero; the scalar equation is not removed
from the full physical problem. Test the third x derivative of R_e's
shape output against the normalized matching packet. Equations (5)-(8)
give the leading pairing

    K0 |k|^8 log|k| int (Q/beta)^2 psi^2 dx.

An order-four local Hessian contributes O_e(|k|^8). By (4), the
fixed-cutoff nonlocal Hessian contributes at most O_e(|k|^4) to this
pairing: move the four output derivatives onto the test function.
Lower logarithmic symbols are o(|k|^8 log|k|). This distributional
test is valid even before a full strong output theorem is known: if
that output is not H3, a finite operator-norm certificate already fails.
It follows that, whenever the strong map exists,

    ||R_e||_(H3xH3 -> H3xH3)
       >= int (Q/beta)^2 psi^2 dx / int psi^2 dx
       >= inf_I (Q/beta)^2 > .9999997, for EVERY e>0.     (9)

This is an operator-norm nonconvergence statement. It does NOT assert
failure of regulator removal on each fixed smooth field, nor contradict
the existing admitted-family weak C2 identity. No finite k witness or
numerical evaluation of the orbit kernel is claimed. Taking k to infinity
after fixing e is essential; a finite set of probe modes cannot certify
the operator norm or refute (9).

The supporting interval calculation uses the actual beta,u,epsilon,tau,
lambda balls from current-global-mass-conversion.json, not nominal seed
values. With b=.003, epsilon<0 and 3beta+u+epsilon>0, the exact polynomial

    Q=beta+u z^2-(4beta+2u+tau)z^3
                         +(3beta+u+epsilon)z^4-2epsilon z^5+epsilon z^6

obeys throughout I

    Q/beta>=q_lower
      :=1-(4+2u/beta+tau/beta)b^3-(|epsilon|/beta)b^6.

The frozen target is q_lower^2>.9999997. The exact polynomial and the
discarded terms' signs are checked; no unrecorded tau guard is used.

This necessary test belongs to the complete graph estimate because it
rules out a claimed uniformly small approximation error even with every
other column solved exactly. It does NOT compute a Schur-eliminated mass
column. Scalar elimination subtracts ordinary order four, so it cannot
remove the shape logarithm. Finite boundary/bolt rows cannot change this
interior principal coefficient. For a proposed full parametrix P with
the TRUE inverse principal symbol, the same argument gives a localized
principal defect of one if H_Q is replaced by H_Q,e,ren; therefore a
raw-cutoff certificate ||I-S P||<=1/2 cannot follow from that replacement.
Existence or a bound for such a P is not assumed proved here. The
original interior Schur return, all phase maps and their norms remain
required. No undefined norm comparison is replaced by (9).

## 6. Quantified sharper alternative and complete remaining obligations

Retain the full analytic ultraviolet operator before attempting a norm
limit. Schematically H_e,matched=H_e,ren+H_UV-H_UV,e, where the paired
finite-geometry realization and original common regulator define BOTH
UV terms. Matching (5) exactly removes the order-one principal defect
in (9); a relative principal mismatch rho would leave the necessary
error floor |rho| in true principal inverse coordinates. Removing this
floor is not a bound on the remaining operator. The actual analytic
route already available is to keep all derivatives of the zeta_A'(-1),
zeta_A'(0), pole/residue and convergent remainder terms of the exact
angular identity, with their shared contacts. A flat symbol alone is
insufficient for the actual boundary realization.

For the complete coupled quotient and original border S, a conforming P
must still certify delta=||I-S P||<1, with sufficient target delta<=1/2,
and the ORIGINAL output maps M_i=||T_i P||. For a column right side b,
put e_b=(I-S P)b. Then the original row error is at most

    M_i ||e_b||/(1-delta) <=2M_i ||e_b|| if delta<=1/2.   (10)

Neither delta, the M_i, nor e_b for the mass, scalar quartic, exact
correlated slow pair and fast columns was evaluated in this attempt.
Equation (10) is the quantitative acceptance target, not a measured bound.
The alternative is therefore specified but its success remains UNKNOWN.

The exact mass primitive and lift remain

    J_m=int_0^z t^3/Q^2, n_m=-2beta Q J_m/z,
    b_m=beta z^3/Q+beta Q' J_m.

Both uncut REFERENCE rows cancel. The ORIGINAL .1/.2 mass cutoff retains
f_cut=2xz chi_m-2F n_m Dchi_m and alpha_cut=n_m Dchi_m, all transition
derivatives, original boundary data and zero reference bolt shift. Its
net reference mass is -xB; it is not a global fixed-temperature zero mode.
The separate original .05/.1 graph cutoff and conversion
2133<||d_chi||H3<2174 still cost <=3.538254568e-12 using the known LINEAR
|m|<1.627532e-15. No change of basis lowers the admission target.
The scalar quartic and B0(z^4-z^(4+e_*))+Delta z^(4+e_*) remain independent
coupled components; e_* is the exact implicit root, never rounded to zero.

Actual post-elimination kappa_m is UNKNOWN. If kappa_m<=5e-5, its
conditional debits are 2.14765e-12 affine, .0024731 feedback and .0195239
energy. All other entries must still jointly satisfy

    Y_other+4.084e-22 A_o+Y_phase<=5.34235e-12,
    eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase<.0012069.

With every other debit zero, the respective individual ceilings are
A_o<=1.30811704211557e10, A_o J_E<3.25134698275862e8 and
gamma_E<.061816542801387. These are conditional allowances, not actual
operator estimates or simultaneously available reserves. The quantum-fast
screen eta_seed+1e-11 C3<.1599164 has the retained admitted-family
two-fast bound .028913+.10706=.135973; reserve .0239434. It does not
control independent phase inputs, the .05/.1 transition or an unknown
nonlinear physical X ball. Equation (9) concerns approximation error,
not eta_seed or C3, and does not revoke their proved scope.

Keep the original preserved-ball and contraction obligations:

    E_center+1e-8 alpha_embed<=1e-11,
    R_X alpha_embed*.028913+Z0_other<=.25,
    R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8,
    Z0+1e-8 Z2<1.

The actual X comparison, center, complete phase/return maps and nonlinear
constants are unknown. The available LINEAR upper bound E<1.1e-10 does
not certify admission at 1e-11; it does not prove that the actual center
is too large. Exact compatibility, physical matching, validity/noise and
all-sector stability remain unproved. No actual controlled solution is
established and no physical feasibility budget closes from this test.

## 7. Decision boundary

Executed: the combined all-spectrum positive-cutoff identity, an ordinary
actual-domain necessary strong-limit test, and its original-normalization
interval bound. NOT executed: a complete signed coupled return, a strong
bound for the matched common-regulator remainder, or actual kappa_m.

DECISION: retain the construction and the exact angular-resummed analytic
method. Exclude a raw positive-cutoff plus local-counterterm approximation
as a uniform strong certificate. Sharpen the analytic UV-matched application,
not the already controlled cold/scalar estimates. The new test strengthens
an existing warning; it does not close the requested feasibility milestone.

ONE DECISIVE NEXT CALCULATION (PROPOSED, NOT EXECUTED): apply the full
angular-resummed common-regulator Hessian with its analytic UV term to
the actual coupled quotient/border parametrix, and enclose the single
strong residual block I-S P together with its original output maps.
Use (10) on mass, quartic, correlated-slow and fast right sides jointly;
test the unchanged residual budgets, including every independent phase
column. A bound for another isolated local allocation cannot decide this.

Remaining to-do list: the complete coupled strong residual and original
output maps; physical-X/neighborhood comparisons; nonlinear admission and
contraction; exact compatibility; physical matching; semiclassical
validity/noise; all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
