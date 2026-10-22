# Finite-centered cold admission and the global transport threshold

2026-10-22, E003/S167. Quantitative audit and frozen verification contract.

The [construction audit](construction-feasibility-audit.md) leaves cold
admission open. This note tests a specific way to close it: center the
family at the actual finite comparison, use the divided-jet complex source
theorem, and reapply real endpoint normalization separately at each real
endpoint. The local coefficient calculation succeeds. The complete
global, center-dependent transport estimate is NOT yet proved. Consequently
the numerical neighborhood and response targets below are conditional,
not a newly certified cold ball. No solver acceptance gate is changed.

This distinction matters. An upper bound for the finite comparison is not
an upper bound for all nearby endpoints, and a supporting polynomial check
cannot supply the missing transported-domain argument.

## 1. What the current construction actually supplies

Keep the actual common disk, paired domain, period, B,A,C,q,mu0 and scheme.
The [finite metric estimate](stationary-completed-metric-response-norm.md)
and [fixed-chart transfer](comparison-path-interior-response.md) give

    sigma_ell=sigma0+ell u, 0<=ell<=1,
    |u|<1.2e-8, |grad u|<7.7e-6,
    P5_40(u)<1.3e-4.

P5_40 retains low value/gradient bounds and divides high jets by j!100^j.
On the cold endpoint normal charts the separately normalized difference
e_ell has bounds epsilon0 a_x j!100^j, epsilon0=1e-6 ell, and the
one-center mixed jets through four cost one further factor 100. Here

    a_x=s[1+log(.01/s)], h<=s<=.001,
    exp(-80)<h<1e-30, L_h<81, a_x<.00331,
    a_x<=2h L_h on chi<=1.

The earlier [general cold theorem](cold-response-from-real-endpoints.md)
requires a different intersected norm: raw fixed-chart C5, actual/cap low
norms, every real endpoint of its norm ball admitted, and cold direction
value/gradient bounds. Neither P5_40(u)<1.3e-4 nor the endpoint normal
bound proves that contract. Even the elementary raw-jet comparison gives

    raw_C5 <= (5!100^5) P5_40 <1.56e8,

whereas the old background threshold is 5e-7. This is a failed upper
allowance, not a lower bound on raw_C5(u). Replacing raw jets by divided
jets only in the statement of that theorem would be invalid.

There is a usable replacement for its COMPLEX source input. The actual
fixed-chart divided-jet proof gives K40<1.25e-11 throughout the complex
coefficient cone P5_40<=.011. Thus a perturbation of norm at most 5e-10
about the finite comparison stays far inside this cone:

    1.3e-4+5e-10 <.001 <.011.

The nonconstant source difference and its two-point oscillation are bounded
by B_c=5e-11 on such parameter disks. This statement retains the actual
paired inverse and the uniformly subtracted full angular sum. It does not
make their real cold bound automatic. The fixed vacuum value cancels in
differences and is not assigned zero value.

## 2. Execute the local endpoint-coordinate calculation

Fix ell and one real cold observation. Freeze its actual endpoint normal
chart once. Write the actual factor as sigma_N,seed+e_ell, with
e_ell(0)=grad e_ell(0)=0. Perturb by a real C5 field v in that fixed chart.
Suppose its local divided jets are bounded by d a_x j!100^j through five.
Put t=d a_x. For d<=5e-10, t<=tmax=1.655e-12.

The NEW real endpoint is normalized by the quadratic map

    H(z)=a z+c z^2,
    a=exp[-v(0)], c=-a^2 v_z(0).

Only this finite coordinate polynomial is holomorphic in space. No spatial
analyticity of v or e_ell is assumed. The bounds |a-1|<=2t, |c|<=100t
hold. On |z|<=7/32 they imply |H-z|<6t, |H'-1|<100t, and H stays
inside the original radius-1/4 chart. H is injective there. Its real
coordinate-column increment is majorized, for T=|dx|+|dy|, by

    J_t(T)=(1+100t)T+200t T^2,
    D_t(T)=6t+100t T+200t T^2.

The new normalized seed-to-endpoint difference is exactly

    e_ell(H)+v(H)+sigma_N,seed(H)-sigma_N,seed
                                                +log|H'|.             (1)

The retained cold seed is analytic on the radius-1/4 bidisk, with
|sigma_N,seed|<.025; use .027. The distance from 7/32 to 1/4 is 1/32.
At each real point the following finite coefficient majorants give the
divided derivatives of (1), for j<=5:

    e_ell(H): epsilon0 a_x [T^j]sum_(k=0)^5(100 J_t(T))^k /100^j,
    v(H):    d a_x [T^j]sum_(k=0)^5(100 J_t(T))^k /100^j,
    seed:    [T^j].027{[1-32(T+D_t(T))]^-2
                                      -[1-32T]^-2}/100^j,
    log:     100t + 2 sum_(k=1)^5 [200t T/(1-100t)]^k/k.

For the log expression divide its j coefficient by 100^j too. Its
constant term deliberately uses 100t. The real part of the complex log
has no larger bound, since H' remains near one. Coefficients of the seed
expression divided by t are increasing in t; evaluate at tmax. This
sums the analytic seed tail instead of truncating it without an error.

Exact rational coefficient evaluation gives the uniform allowances

    composition multiplier <1.000001,
    seed shift /t <20,
    logarithmic Jacobian /t <=100,
    local new endpoint size
        <=1.000001 epsilon0+122 d.                 (2)

This is a local finite endpoint bound in the unchanged derivative units.
It is not a Lipschitz estimate for e_ell(H)-e_ell. The latter would need
additional regularity: for e(x)=|x-x0|^(11/2), translations change its
fifth derivative by order sqrt(|shift|), not order |shift|. The calculation
instead bounds the whole e_ell(H) using its known C5 norm and spends an
explicit fixed background margin. It never differentiates H in a family
parameter while taking five spatial derivatives of e_ell.

The same approach cannot simply be declared valid for the center-dependent
global arrays. In that problem H depends on the observation, the extension
cutoff and its baseline chart, and the density/spin coefficients change.
Their bounds are the missing calculation identified next.

## 3. Exact remaining admission target

Define the intended direction space on the fixed actual common manifold,
with real C5 cold fields, regular even bolt, unchanged outer projector and
clock. A direction norm D must bound all of:

1. P5_40 of the field, using the original fixed charts and caps;
2. its local cold divided jets in the FROZEN finite endpoint charts,
   divided by a_x, and one-center mixed jets through four, with the
   additional factor 100;
3. |v|/a_x and |v_chi|/(8a_x) for the Ward kernel;
4. the actual and H2-cap low coefficient fields with the original
   exp(distance) envelopes, including their fixed cutoff extensions and
   the required one-center transported coefficients.

Items 1-3 are explicit coefficient quantities. Item 4 must be written as
linear fixed-chart data and bounded under finite endpoint normalization;
it must not be defined as 'whatever makes the source theorem hold'. No
inverse, determinant or desired quantum-source bound belongs in D.

Let K_tr be an UPPER constant for that complete finite transport, in the
same units, so that every real perturbation D(v)<=d has the original
cold endpoint contract with

    epsilon_end <=1.000001 epsilon0+K_tr d.        (3)

The local calculation proves the local contribution with 122. It does
NOT prove K_tr<=122, or even K_tr<=1000, for the global contract.

For a proposed real endpoint tube d<=5e-10 and epsilon_end<=2e-6,
the exact threshold at ell=1 is

    K_tr <= (2e-6-1.000001e-6)/(5e-10)
          =1999.998.                              (4)

Thus K_tr<=1000 would leave a strict reserve and admit a background ball
D(v)<=rho=2.5e-10 with another rho for parameter variation. This is a
target, not an assigned value for an unknown constant. The local 122
allowance leaves more than a factor sixteen to the threshold, so the
remaining transport is a concrete question rather than a loss of hundreds
of orders from the global inverse chain.

A second required check is the REAL quantum endpoint theorem on
epsilon<=2e-6, twice its published contract. Its local potential cone has
room: exp(.025+2e-6*.00331)-1<.026, and the local curvature perturbation
2!100^2(2e-6)(.00331)<.000133. These checks do not alone extend the
theorem. The H2-cap weighted graph estimate, actual-to-cap separation,
one-center mixed gluing and finite restoration must be revalidated on
that entire tube. Keep this as a separate hypothesis until done.

## 4. Consequences IF the two remaining hypotheses are proved

Assume (3) with K_tr<=1000 and the real endpoint bounds on epsilon<=2e-6
with their retained coefficients 1e23 and 1e24. For every real endpoint
in the proposed tube they would give

    |S[v]-S_seed| <=3.4e-15 a_x,
    |delta S(chi)-delta S(0)| <=3.4e-14 h L_h sqrt(chi). (5)

The actual fixed-chart complex bound B_c=5e-11 is already available.
For the first line log(B_c/(3.4e-15 a_x))<91. For the second use
91+.5log(1/chi). Apply the half-disk maximum-principle argument in
[the general cold proof](cold-response-from-real-endpoints.md), with
parameter radius rho/D(u), rho=2.5e-10. Its scalar derivative constants
are 2e a L/rho and 8e a L^2/rho^2. Polarization costs a further factor
two for the mixed derivative. Consequently the conditional bounds are

    |D S[u]| <.007 D(u) a_x,
    |D^2 S[u,w]| <2e10 D(u)D(w) a_x,
    |D_u S(chi)-D_u S(0)| <.07 D(u)hL_h chi^(1/4),
    |D^2_(u,w) S(chi)-D^2_(u,w) S(0)|
                                 <2e11 D(u)D(w)hL_h chi^(1/4).        (6)

The logarithmic loss and quarter-power modulus are necessary in this
argument. The log branch is recorded in [DLMF 4.37](https://dlmf.nist.gov/4.37).
Use exp(-t/4)(91+t/2)^j<=91^j, j=1,2, t>=0. The sharper tail tends
uniformly to zero after division by chi^(1/4), so the conditional theorem
would be C2 into that weighted space, not just pointwise differentiable.

Using the original canonical DEVIATION Ward kernel and its bounds
1.2,700D(u)a_c,5500D(u)D(w)a_c^2, the corresponding conditional
density/null allowances are

| Weight | Value | First /D(u) | Mixed /D(u)D(w) |
| --- | ---: | ---: | ---: |
| a_x | 1.1e-6 | .024 | 6.9e10 |
| hL_h chi^(1/4), null only | .0017 | .24 | 6.9e11 |

The homogeneous scalar is retained separately. This calculation uses
the massive canonical source; independent Weyl, LLL and radius responses
are still required for the full operator. With the separate low radius
hypotheses of [the constraint transfer](fixed-chart-global-constraint-transfer.md),
the cold canonical deviation mixed allowances would be <1.6e12 for the
zero-datum mass primitive, <7e4 for the absolute mass density and <8e14
for the absolute lapse-density L1 norm. These are conditional upper bounds.

At a correction radius 1e-10 the last Taylor allowance is 4e-6. It is
not a preconditioned Newton bound. This method could establish a
neighborhood, but its coarse response allowance still does not certify
a small full solution. A norm change D'=cD rescales rho and the first/
mixed coefficients together; the products L1*r and L2*r^2 do not improve.

## 5. Execute the global coefficient budget, without assuming its inputs

There is a useful second decisive calculation. The transported half-density
operator is D_G={B^a,partial_a}/2, with B=exp(-nu) gamma K, K=(DG)^-1.
For a direction let A=delta DG, k=delta K=-K A K. Then exactly

    delta B=exp(-nu)(k-delta nu K),
    grad k=-(grad K)AK-K(grad A)K-KA(grad K),
    grad delta B=exp(-nu)[grad k-(grad delta nu)K
                  -delta nu grad K-grad nu(k-delta nu K)].             (7)

The first-center formula is the same product rule with grad replaced by
the center derivative. These identities retain the Jacobian/connection
terms through the symmetric divergence expression; deleting grad k or
the two endpoint factors would not describe the actual realization.

Use the following SUFFICIENT input arrays, with weighted direction size
t=d a_x exp(distance), and mixed arrays needed near the bolt only:

    ||DG||<=1.5, ||K||<=3, ||D2G||<=23,
    exp(-nu)<=1.1, |grad nu|<=8, |grad K|<=210,
    |delta nu|<=100t, |grad delta nu|<=21000t,
    |A|<=1000t, |grad A|<=2e6t,
    |delta nu_x|<=10000t, |A_x|<=1e5t.

For the center cross products additionally require |nu_x|<=10 and
|K_x|<=1 on chi<=1. The seed arrays have much stronger a_x factors
there. These two bounds must still be verified on the transported finite
family, including its cap, rather than inferred just from seed smallness.

Substitution into (7) yields

    |k|/t<=9000,
    |grad k|/t<=2(210)(1000)(3)+9(2e6)=19260000,
    |delta B|/t<=1.1(9000+100*3)=10230,
    |grad delta B|/t
        <=1.1[19260000+21000*3+100*210+8(9000+100*3)]
         =21360240.

For the center term, |k_x|/t<=2(1)(1000)(3)+9(1e5)=906000.
The corresponding allowance is

    |delta B_x|/t
        <=1.1[906000+10000*3+100*1+10(9000+100*3)]
         =1132010.                                (8)

Relative to the retained coefficient contract 1e4,2e7,1e6 these costs
are all below two. Consequently the operator algebra itself would fit
comfortably inside K_tr=1000. It is not the numerical source of a
large loss. What remains is to derive these input arrays uniformly for
the actual center-dependent maps and the global cutoff extensions from
the explicit linear direction norm. Seed arrays and scalar endpoint
jets are not that derivation. In particular a candidate bound for
|D2G| does not establish |A_x|, and a local graph bound does not establish
the distant comparison's weighted envelope.

The real doubled-tube graph test also has numerical room. If the actual
H2 cap preserves |sigma_cap|<.05 and |sigma_cap'|<.5, then for its
conformal mass operator T=D_H2+M exp(sigma_cap), a distance weight of
slope two gives

    M||T_weight^-1|| <=1/(exp(-.05)-2/18)<1.191,
    M||S_weight|| <=exp(.05)/(exp(-.05)-2/18)<1.3.

The D_H2 T^-1 identity and its first-order connection term give gradient
allowance <4.2. With the above finite-map bounds, transport costs at most
1.5(4.2)+(3*23)(1.3)/18<12, below twenty. Thus the inverse threshold
is not algebraically inconsistent. Verifying those cap value/gradient
hypotheses AND the actual far-seam agreement is still required.

## 6. Decision, verification and next calculation

Continue the finite-centered admission route only through its decisive
transport calculation: derive the actual/cap one-center coefficient arrays
from fixed finite charts and prove or refute the numerical target
K_tr<=1999.998, while rechecking the doubled real endpoint tube. The
local normalizer and global operator-budget calculations have been
executed here. The local allowance is 122; the operator algebra costs
less than two times its retained coefficient reserves, PROVIDED its
listed input arrays hold. The full K_tr remains UNKNOWN. Do not record cold admission
as complete or silently carry this conditional result into the solver.

The alternative, if this transfer does not fit, is a direct weighted
fixed-chart resolvent response around the finite metric, avoiding repeated
normalization and the endpoint-to-derivative loss. A large transport upper
bound would call for that sharper method; it would not disprove wormholes.
The direct actual correction/moment route and physical boundary/validity
work remain the controlling program. No expansion of the parked generic
source-times-inverse chain is selected.

Freeze this note, review and checker before execution. Check rational
majorants, operator identities and admission thresholds; finite real endpoint normalization on
a manufactured C5 non-C6 field; the failure of a C5 composition Lipschitz
claim; analytic smallness and Ward cross terms at 60/90 digits. Normalize
tiny cold model values before comparisons. These are supporting mechanisms,
not samples of the actual quantum state or proofs of the missing transport.
Preserve every run, source snapshot and input hash. No Lean is invoked.

Remaining to-do list: complete global/center-dependent cold transport and
doubled endpoint tube; actual correction and signed moment; full nonlinear
map, projected solution and exact compatibility; physical boundary/state,
validity, matching, stability, charge/return, opening, transport,
preparation, separate creation and full-result audit.
