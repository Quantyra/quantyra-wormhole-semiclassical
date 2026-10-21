# Fast traces in the retained finite-cut action

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The retained boundary action does not automatically extend to the simplest
fast conformal paths. A static fast direction has a divergent linear
boundary term. For a nonzero thermal harmonic its time average vanishes,
but the quadratic relative action still has a nonconstant oscillatory
cutoff term, for every nonzero fixed radial phase. Thus the preceding
isotropic trace relation alone is insufficient to complete this action.

A faster, second-order metric tail can change the quadratic cutoff term.
We derive its required leading cancellation equation, without selecting
it as a physical solution or modifying any counterterm. The result calls
for a nonlinear boundary completion; it does not exclude every completed
fast metric, decide stability, or establish a wormhole no-go theorem.

## 1. Finite-cut contract and conformal test paths

Keep the same charged seed, action, physical flux, subtraction and paired
Dirac prescription as the [boundary pairing](fast-curvature-boundary-symplectic.md).
Use a fixed background optical/areal chart, and a cut at background r=R.
The embedding of the regulating surface is held fixed when varying the
metric. Its physical area changes. No equivalence with a field-dependent
fixed-physical-area cutoff is asserted without matching its variation.

Write g_eta=exp(2eta p)g_seed, with the independent optical shape fixed.
Choose p smooth on the open exterior, zero near the horizon, and on an
outer collar put

    z=log(r/B), p=exp(-3z/2) f(z) H(tau),
    f=A_c cos(kz)+A_s sin(kz),
    k^2=1/t+17/12, t=epsilon lambda, epsilon=G/(60pi).

Here B is the unchanged reference radius and H is either 1 or a fixed
nonzero thermal cosine. The sign of the sine amplitude relative to an
inward-coordinate trace convention is immaterial, and is not silently
used to change a boundary phase. Take relative actions I_R[g_eta]-I_R[g_0].
At the same finite cut the canonical conformal metric, domain and period
are unchanged, so the independent canonical determinant cancels exactly
from this difference. Its local Weyl cocycle is retained in full.

Use the stationary charge action sign throughout. The Euclidean action
has the opposite overall sign; convergence conclusions are unchanged.
Coefficients below are Taylor coefficients in eta, not Hessians with an
extra factor of two. Time averages are divided by the fixed period.

The boundary terms are exactly the Einstein GHY and AdS counterterms,
the actual paired-domain Weyl transgression with its conversion fluxes,
and the required loop terms c0,c1 from the
[same-source charge calculation](stationary-same-source-boundary-charge.md):

    c0=-221lambda^(3/2)/(15120pi^2),
    c1=29sqrt(lambda)/(20160pi^2).

No new finite extrinsic functional or independent Robin matter condition
is introduced. The earlier stationary phase excluded these oscillatory
directions; its differentiability theorem is not contradicted.

## 2. Complete leading densities, with boundary terms retained

For the static calculation use the fixed canonical proper coordinate x,
increasing outwards, with D_x=(r sqrt(F)/B)partial_r. Let

    N=B sqrt(F)/r, a=D_x log N, sigma=log(r/B)+eta p,
    s=D_x sigma, u=D_x s+a s, v=s^2, Y=B^-2,
    Rhat=-2(D_x a+a^2).

The converted Weyl bulk representative accompanying the retained
boundary polynomial is

    L_W=B^4/(4pi r^2) {
      sigma[(Rhat+2Y)^2/60-11Y Rhat/90-q_magnetic^2 Y^2/3]
      +Rhat u/30-u^2/10+(Rhat/30-Y/18)v-7uv/90-7v^2/180 }.

This uses the two integrations converting sigma Box Rhat into Rhat
Box sigma, with their boundary flux included. On the unperturbed areal
chart, B^2 u=r F_r-F; omitting the -F would change the result.
The boundary contribution is

    B_W=-B^2 N Gcal/(360pi),
    Gcal=sigma[12a D_x a+44a^3/7+10aY]
      -6(D_x a)s-24a^2s/7-6a D_x s-96as^2/7
      -18s D_x s-85s^3/21.

These are the [actual finite-cut bulk/boundary pair](stationary-exterior-boundary-weyl.md),
not the literal boundary heat path added to an already converted bulk.

For independent normalization, the conformal Einstein bulk density is

    L_EH=r^2/(4G){exp(2eta p)[R_seed-6eta Box_seed p
                       -6eta^2 F p_r^2]+6lambda0 exp(4eta p)},
    lambda0=lambda+11epsilon lambda^2/6.

The conformally invariant Maxwell contribution cancels from the relative
action. At a fixed cut the GHY density is
r^2 sqrt(F)exp(2eta p)[K_seed+3eta sqrt(F)p_r]/(2G).
The leading Einstein counterterm is -sqrt(lambda)r^2 sqrt(F)exp(3eta p)/G;
its intrinsic-curvature partner and c1 are retained but subleading below.

Divide the leading z-bulk density and boundary density by
S0=lambda B^3/G, and factor out exp(3z). For a generic static p(z),
write primes for z derivatives. The first two Taylor coefficients are

    L1=11t p-(9+25t)p''/6-(27+53t)p'/6,
    B1=3p'/2+(t/6)[18p''+(337/7)p'-(442/7)p],          (1)

    L2=(6+22t)p^2-9pp'-3pp''
               -(3+17t)(p')^2/2-(16t/3)p'p''-(3t/2)(p'')^2,
    B2=3pp'-3p^2/2+(t/6)[(211/7)(p')^2+18p'p''-(663/7)p^2]. (2)

The lambda0 shift and covariant volume counterterm both enter these
coefficients. Neither can be dropped before the limit.

## 3. Linear static failure and quadratic thermal failure

The exact leading primitive of exp(3z)L1 is

    exp(3z)[-(9+25t)p'/6+(11t/3)p].

Adding exp(3z)B1 leaves exp(3z)(3t/7)[7p''+9p'-16p].
Thus the linear static relative-action coefficient has outer behavior

    constant + lambda^2 B^3/(140pi) exp(3z/2)
                [-(7k^2+55/4)f-12f']+o(1).             (3)

The linear map on (A_c,A_s) has determinant
(7k^2+55/4)^2+144k^2>0. Every nonzero static fast direction therefore
has an oscillating unbounded linear term. The interior contributes a
finite constant, not a cancelling outer power.

For a nonzero thermal cosine, all first-order terms integrate to zero
because the seed is static. The quadratic term remains. Substituting
p=exp(-3z/2)f and f''=-k^2 f into exp(3z)L2 gives

    l2=(17t/2+65/8+3/(2t))f^2+(41t/9+5/6)ff'
                                               -(6t+3/2)(f')^2.

Its exact primitive is

    P2=-(6t+3/2)ff'+(41t/18+5/12)f^2.

The zero mean of l2 over a radial cycle uses k^2=1/t+17/12; it is not
a prescription to discard the remaining endpoint term. Adding the full
boundary coefficient gives

    Q_t(f,f')=-(3001t+546)f^2/504
                 -(71t+21)ff'/14-(167t/42)(f')^2.       (4)

Consequently, for any fixed nonzero thermal harmonic,

    [eta^2](I_R[g_eta]-I_R[g_0])
       =constant + S0 <H^2> Q_t(f(z_R),f'(z_R))+o(1).   (5)

Here <H^2>=1/2 for a unit cosine. The mixed coefficient in (4) is
strictly nonzero for t>0. For f=a cos(kz+phi), its sine-of-double-phase
coefficient is a^2 k(71t+21)/28. Thus no nonzero phase makes Q_t
constant. Two sequences of arbitrarily large cutoffs have distinct
quadratic limits; averaging the oscillation is not the retained action.

The same conclusion holds for the finite relative action along this
simple logarithmic path at sufficiently small fixed nonzero eta. Terms
of cubic and higher order have integrable bulk tails and vanishing outer
boundary tails, so they cannot remove the quadratic oscillation.

## 4. Endpoint estimates and fixed nonzero frequencies

For the exact charged seed F=lambda r^2[1+O(r^-2)], all differentiated
corrections have the same improved radial powers. In canonical proper
coordinates N tends to B sqrt(lambda), a=O(r^-1), and Rhat is bounded.
Each normal derivative raises an outer power by at most one. At linear
order the omitted bulk terms are O(r^-3/2 log^m r) in dr and boundary
terms O(r^-1/2 log^m r); hence (3)'s remainder tends to zero. At quadratic
order the corresponding bounds are O(r^-3 log^m r) and
O(r^-2 log^m r), respectively. The leading interior constants in (3),(5)
exist without any cutoff averaging.

For a fixed thermal frequency, time derivatives do not raise the radial
power in the canonical metric. The parity-even scalar density pairs
tangential time derivatives. Replacing two normal derivatives by two
time derivatives therefore gains two radial powers. The bulk a4 terms
have derivative weight four and the boundary terms weight three,
including derivatives of the smearing function. This follows from the
full mixed heat coefficient in [Vassilevich, (5.33)](https://arxiv.org/html/hep-th/0306138v3).
Its tangential-divergence and smearing terms are included in this count,
not set to zero individually. Thus frequency terms change the finite
constant in (5), but not its nonconvergent leading coefficient.
This is not an estimate uniform as frequency grows with the cutoff.

For the nonlinear remainder, the worst local bulk bound is
O(r^2 |p|^3)dr=O(r^-5/2)dr, and the worst boundary bound is
O(r^3 |p|^3)=O(r^-3/2), with fixed-frequency derivative constants.
The Weyl action is polynomial in conformal derivatives and the Einstein
exponentials have uniform Taylor bounds for small eta. These bounds also
control amplitude derivatives through order two. The same canonical
determinant cancels at each finite cut, so no separate nonlocal cutoff
interchange is used for this conformal family.

## 5. The faster nonlinear tail is part of the missing completion

Failure of the simple path does not determine the whole nonlinear
boundary phase. Consider only as a diagnostic

    sigma=log(r/B)+eta exp(-3z/2)f(z)H(tau)
                            +eta^2 exp(-3z)h(z)H(tau)^2.

The added tail is smaller and does not change the linear fast amplitudes.
Its linear boundary term contributes at quadratic order

    S0 <H^2> (3t/7)[7h''-33h'+20h].                   (6)

Thus the oscillatory action value in (5) can be cancelled at leading
order if

    (3t/7)[7h''-33h'+20h]=-(Q_t-mean Q_t).             (7)

There is a unique zero-mean periodic solution on the doubled radial
frequency. On its cosine/sine coefficients, the matrix is

    [[20-28k^2,-66k],[66k,20-28k^2]],

whose determinant is (20-28k^2)^2+4356k^2>0. A constant h changes a
finite term instead of the oscillation. This calculation identifies a
concrete second-order boundary matching obligation. It does not prove
that the nonlinear field equations, the shape response, constraints,
matter preparation or a fixed-source variational principle admit that
tail. It is not adopted as a counterterm, solution or physical phase.

The preceding nonzero boundary symplectic pairing still requires an
appropriate trace relation. Cancellation of an action value is not a
replacement for its boundary variation. A nonlinear boundary condition
must match both, including any coupled shape or exterior data.

## 6. Verification and next decision

The checker recovers (1)-(4) by differentiation of the retained nonlinear
densities, checks the literal boundary-path conversion, and integrates
the exact charged finite-cut coefficients on refined outer annuli at two
precisions. Independent complex-amplitude extraction checks the Taylor
normalization against the unexpanded finite-cut action. Phase tests
exhibit distinct cutoff limits and verify (7) algebraically. Fixtures
are not a calibrated physical benchmark or a causal mode calculation.
The continuum estimates remain ordinary proofs; no Lean run is selected.

This closes the direct promotion of the simple fast conformal paths into
the existing stationary action. It leaves a concrete nonlinear boundary
completion task, rather than selecting a new phase or removing the fast
branch. No physical boundary exclusion, stability, exact nonlinear
solution, or wormhole impossibility result follows. All physical gates
remain unchecked and the full research goal remains active.

Remaining to-do list: nonlinear gravitational boundary completion and
source compatibility; continuous causal response and constrained
stability; nonlinear existence and signed compatibility; exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and full audit.
