# Proper mouth clock, static anomaly and the moving cutoff

2026-10-11, E003/S051-S053. The fixed proper-time mouth condition removes
the length-dependent `q/h` static term found in the canonical lapse-only
diagnostic. For the regulated AdS2 core and the explicitly controlled
static joins below, the covariant free-channel energy reproduces MMP's
`+q/(24 ell)` anomaly contribution, with force and stiffness errors
`O(q epsilon/ell^2)`. This is a reproduction and a controlled static
completion within a specified two-dimensional metric class. It is not a
four-dimensional matched solution, a physical damping theorem or accepted
novelty.

The moving-cutoff calculation also evaluates the *local Weyl cocycle of
the isolated throat patch*. Its two-derivative kinetic coefficient is
`O(q epsilon)`. That patch functional is not the transparent CFT influence
action: boundary determinants and the exterior gluing contribution still
have to be combined. In particular this calculation does not yet determine
the physical clock's inertia. The distinction resolves the static part of
the [collar diagnostic](collar-conservative-response.md), while retaining
its central warning about inferring stability from absorption alone.

## Static energy from one covariant action

Use central charge `c` for one nonchiral channel; the charged complex
fermion has `c=1`. On a smooth static spatial circle let

`g=-N(s)^2 dt^2+B(s)^2 ds^2`, `v=N/B`, `J=int ds/v`,

with positive `N,B`, a fixed antiperiodic spin structure/ground-state
holonomy, and the conformal renormalization used in the preceding note.
No new finite dimensionful metric or boundary counterterm is added here.
The state's conformally mapped cylinder vacuum has energy

`E[g]=-pi c/(6J)-c/(24pi) int v (partial_s log N)^2 ds`.       (1)

This is the total conserved Killing energy, including the anomaly. To
derive it rather than discard the canonical contact term, write
`g=B^2 g_A`, `g_A=-v^2dt^2+ds^2`, `rho=log B`. The previous exact
Hamiltonian result gives

`E_A=-pi c/(6J)-c/(24pi) int v_s^2/v ds`.

Since `R_A=-2v_ss/v`, the static Weyl cocycle per unit time is

`(W[g]-W[g_A])/time`
` =c/(24pi) int [-2rho v_ss+v rho_s^2]ds`
` =c/(24pi) int [2v_s rho_s+v rho_s^2]ds`.

Energy is minus the static effective Lagrangian. Completing the square
gives (1). Thus the canonical `v_s^2/v` term is neither set to zero nor
subtracted by hand. It combines with the spatial metric's anomaly. In
optical distance `dy=ds/v`, (1) is

`E[g]=-pi c/(6J)-c/(24pi) int (partial_y log N)^2dy`.          (2)

The integrated second derivative in the local stress tensor vanishes on
the complete smooth circle. Applying that cancellation separately to an
open throat without its joining terms would be incorrect. Equation (1)
also exhibits spatial-coordinate invariance; it uses the physical lapse
`N`, rather than the canonical velocity `v` alone.

## The proper-time cutoff fixes the static answer

The leading constant-curvature throat is

`g=r_e^2 sec(sigma)^2[-x^2 dt^2+d sigma^2]`, `x=1/ell`.

Fix the induced mouth metric to `-r_e^2 dt^2/epsilon^2`, where
`r_e << epsilon << ell`. Both cutoffs are then at `sigma=+-a`, with

`a=acos(epsilon x)`, `b=sqrt(1-epsilon^2 x^2)`,
`T_epsilon=2a/x`.                                           (3)

The optical-coordinate lapse is `N=r_e x sec(x y)`. Direct integration,
including both halves, gives

`int_core (partial_y log N)^2dy=2[b/epsilon-x a]`,
`E_an,core=-qc/(12pi)[b/epsilon-x a]`
` =-qc/(12pi epsilon)+qc x/24-qc epsilon x^2/(24pi)`
`   -qc epsilon^3 x^4/(288pi)+O(qc epsilon^5 x^6)`.            (4)

The leading cutoff term is independent of `x`. Omitting it from the
length potential does not alter any force or stiffness. It is not the
profile-dependent `q/h` subtraction prohibited by the earlier diagnostic.
For `D_phi=x partial_x`, exact differentiation gives

`D_phi E_an,core=qc x a/(12pi)`,
`D_phi^2 E_an,core=qc/(12pi)[x a-epsilon x^2/b]`.              (5)

Both differ from `qc x/24` by `O(qc epsilon x^2)` uniformly when
`epsilon x <= 1/4`. For example the second error is at most
`qc epsilon x^2/(6pi sqrt(15/16))`: use
`asin z <= z/sqrt(1-z^2)` for `0<=z<=1/4`.
This reproduces the finite MMP anomaly term, not a new support mechanism.

For exterior optical lengths `D_j>=0`, the free-circle Casimir energy is
`-pi c/6 sum_j (T_epsilon+D_j)^-1`. Replacing `T_epsilon` by `pi/x`
changes this energy and its first two `D_phi` derivatives by
`O(qc epsilon x^2)`, uniformly over all the `D_j`. To see uniformity,
write each summand as `x/[2acos z+D_j x]`, compare with
`x/[pi+D_j x]`, and differentiate twice. The denominators are bounded
below by `2acos(1/4)>2`; ratios `D_j x/(2acos z+D_j x)` are bounded by
one. The difference `pi-2acos z` and its first two `D_z` derivatives are
`O(z)`. No finite mean or upper bound on the long returns is needed.

After combining core, fixed exterior and gravitational terms, the
specified static metric class therefore has potential

`E_static(x)=A x^2+qc x/24-pi c/6 sum_j (pi/x+D_j)^-1`
`             +constant+O(qc epsilon x^2)`,                   (6)

where the remainder has the same bound through two clock derivatives.
The gravitational coefficient `A=r_e^3/G` is derived in the
[collective-clock note](collective-clock-response.md). We have not
rederived its higher-order matching corrections in (6).

## A smooth join and a derivative bound

A derivative jump at the mouth is not a matched Einstein solution. The
following construction shows precisely what static two-dimensional
smoothing costs; it does not assert the Einstein constraints.

At the right mouth set `w=y-a/x`, `u=w/epsilon`, `z=epsilon x`. In the
last optical distance `epsilon` of the core, `-1<=u<=0`,

`N_in/N_b=[cos(zu)-sqrt(1-z^2) sin(zu)/z]^-1`,
`N_b=r_e/epsilon`, `N_0/N_b=(1-u)^-1`.                        (7)

Use the continuous value at `z=0`. The fixed exterior can continue
`N_0` near `u=0` and then connect to any specified positive smooth static
return geometry independently of `x`. Define

`r(z,u)=log(1-u)-log[cos(zu)-sqrt(1-z^2) sin(zu)/z]`,
`log N_blend=log N_0+chi(u) r(z,u)`,

where `chi=1` at `u=-1`, `chi=0` at `u=0`, with vanishing first two
endpoint derivatives. One explicit C2 choice is
`chi=1-10v^3+15v^4-6v^5`, `v=u+1`; a smooth step can be used instead.
The two joins are treated symmetrically.

On `0<=z<=1/4`, `-1<=u<=0`, the bracket in (7) is at least
`cos(1/4)>0`. The function `r` is even analytic in `z`, and `r(0,u)=0`.
Consequently

`sup_u |partial_u^j D_z^k r(z,u)| <= C z^2`,
`j<=2`, `k<=2`,                                             (8)

for a finite constant independent of `x,epsilon,D_j`. This is an ordinary
compactness/Taylor bound, not an inference from sampling: extend
`r/z^2` to zero using `int_0^1(1-t) r_zz(tz,u)dt` and take the finite
suprema of its listed derivatives on the compact rectangle.

The unblended `partial_u log N_in` is uniformly bounded there. Expanding
the square in (2), the difference in the core anomaly energy and each of
its first two clock derivatives is bounded by

`C qc z^2/epsilon = C qc epsilon x^2`.                       (9)

The support has optical length `epsilon`, and `D_phi=D_z` in this
coordinate integral. This proves the asserted static error for this
class of joins. The relative curvature error in the collar is
`O(z^2)` compared with `r_e^-2`; that estimate is not a solution of the
four-dimensional constraint equations. Other finite local terms from
the four-dimensional reduction and changes in the exterior remain to
be bounded in their physical scheme.

The error in (6) is smaller than the long-return scale `q delta/ell`
if `epsilon/ell=o(delta)`, where `delta=d/ell`. Together with the
previous nominal classical error `r_e/epsilon=o(delta)`, this requires
`r_e ell/d << epsilon << d`, already identified in the earlier note.
For the joint weak-coupling example `delta=q^-1/16`, `g^2=q^-2`, one
may use `epsilon/ell=q^-1/4`. Then `epsilon/d=q^-3/16` and
`r_e/epsilon=O(q^-11/4)`. These are compatibility checks, not bounds
on all omitted physics.

## What the time-dependent cutoff does and does not fix

Let `tau'=x(t)>0` and take symmetric boundaries `sigma=+-s(t)` in the
same global throat metric. Their *induced* metric, including the motion,
imposes exactly

`cos^2 s=epsilon^2[x^2-sdot^2]`,
`beta=sdot/x`, `|beta|<1`.                                  (10)

Keeping only `cos s=epsilon x` in a dynamical action misses a term at
the same two-derivative order as the boundary anomaly. For slowly
varying histories, use derivative counting at fixed `epsilon x<1`:

`s=acos(epsilon x)+epsilon^3 xdot^2/(2x b^3)+O(partial_t^4)`.
                                                                  (11)

The reference metric `g_hat=-x(t)^2dt^2+d sigma^2` is flat, while
`rho=log r_e-log cos sigma`. Use the boundary Weyl cocycle

`W_loc=qc/(24pi)[int sqrt(-g_hat)(partial rho)^2`
`                         +2 int_boundary rho K_hat ds_hat]`.      (12)

This fixes the boundary scheme of this *patch comparison*. For either
outward timelike boundary, `K_hat ds_hat=d(atanh beta)`. For periodic
histories, or histories with the endpoint term retained separately,
spatial integration and one time integration by parts give exactly

`W_loc=qc/(12pi) int x[tan s-s-2 beta tan s atanh beta]dt`,           (13)

up to the endpoint `qc/(6pi)[rho_b atanh beta]`. In deriving (13),
`rhodot_b=tan s sdot`; both boundaries have the same outward curvature.
Equations (10)-(13) involve no field equation for the clock.

Substituting (11), and writing `phi=log x`, yields

`W_loc=qc/(12pi) int [b/epsilon-x acos(epsilon x)`
`                     -3epsilon phidot^2/(2b)]dt`
`          +O(partial_t^4)`.                                (14)

The bulk domain displacement contributes `+epsilon phidot^2/(2b)`;
the boundary term contributes `-2epsilon phidot^2/b`. Their sum fixes
the coefficient and sign. The fourth-derivative notation means a
derivative expansion, not a bound for arbitrary fast motion. For smooth
profiles with `x`, `phidot` and successive derivatives of order
`ell^-1`, `ell^-1`, etc., separated from `epsilon x=1`, the next
local terms have size `O(qc epsilon^3/ell^4)` per unit time.

In this patch functional only, the mass correction is
`Delta M_loc=-qc epsilon/(4pi b)`. Compared with `M_clock=q ell/8`
at `c=1`, its absolute relative size is
`2epsilon/(pi ell b)`, which vanishes in the overlap limit.
This coefficient must not be inserted as the full physical mass.

For an actual transparent mouth, the throat boundary is an artificial
cut. The CFT continues into the exterior. Boundary terms cancel between
pieces only when the reference metrics, induced data, normals and the
state-dependent determinant are glued consistently. Neither a reflecting
strip determinant nor a disk determinant is that gluing functional.
Equation (14) controls one necessary local contribution; it neither
proves cancellation of all other terms nor rules out a leading real
response from the full geometry. The earlier positive-frequency
absorption and sharp subtracted kernel remain separate results.

## Source review and verification

[MMP v3, section 5.2](https://arxiv.org/abs/1807.04726v3) already obtains
the finite anomaly contribution by comparison with the reflecting AdS2
vacuum. Equations (4)-(6) reproduce that contribution with a proper-time
cutoff and keep derivatives of the regulated answer. The leading
constant-curvature and boundary-clock framework is the previously
reviewed MSY/Iliesiu--Turiaci reduction; it is not a new construction.

[Yang v3, appendix C](https://arxiv.org/html/1809.08647v3#A3), read in
full, gives the boundary Weyl functional and the induced Schwarzian
correction for a Euclidean disk. Its boundary parameter called `q` is
not our magnetic flux. The last step to its equation C.118 is stated as
a direct evaluation rather than expanded into intermediate algebra; we
have not independently reconstructed that disk determinant. Its scope
does not establish transparent two-mouth gluing.
[Blommaert--Mertens--Verschelde, section 2 through equation 2.6 and its
footnotes](https://doi.org/10.1007/JHEP09(2019)060), read for this issue,
explicitly limits the subleading anomaly statement to its setup and
cites Yang. The [Nakayama--Suzuki PTEP 2017 discussion](https://doi.org/10.1093/ptep/ptx118)
of a Liouville field was screened in sections 2-3 through equation 23;
its gravity/ghost central-charge
prescription is not imported into the matter-only coefficient here.

The [verification script](../checks/proper_clock_anomaly.py) independently
checks the static energy identity by optical-coordinate and canonical-
plus-Weyl quadratures, the core integral and its derivatives, the smooth
join and its first two clock derivatives, the return-uniform estimates,
the exact induced lapse, both forms of the moving-boundary action, and
the predicted order of the two-derivative remainder. The numerical
samples support the analytic calculations; the uniform claims follow
from the bounds above. Exact source and dependency hashes are in the
[receipt](../evidence/proper-clock-anomaly.json); scoped source provenance
is in the [source supplement](../evidence/proper-clock-source-review.json).

The decision is to continue the existing physical primary/fallback with
the static correction fixed, and next derive the *whole* Lorentzian
transparent influence action. Do not promote (14) into a complete
inertia or damping calculation. No Lean/Lake, GCP or external review
invocation occurred in this milestone. Full novelty remains undecided.

Remaining to-do list: complete time-dependent transparent metric/state
gluing and its conservative response, bound interactions and channel
matching, finish primary/fallback novelty and the consequential result,
independently verify and GCP-certify suitable endpoints, and prepare the
reproducible manuscript.
