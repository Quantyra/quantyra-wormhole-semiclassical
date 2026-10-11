# Static load normalization and the dynamical research boundary

2026-10-10. S051/S052/S053 baseline and source-proof audit. This note does
not establish a new stability result. It resolves two normalization
problems in a decisive static-load comparator and excludes its leading
load potential as a new research target. The finite-time quantum-state
question in [the preceding decision](rotating-dirac-and-state.md) remains
open; a stationary load is not a substitute for that state construction.

## Model and source editions

Use the leading, long-throat MMP approximation with Einstein-Maxwell
gravity and `N_f` free massless charged Dirac flavors, effective lowest
Landau-level central charge `c=N_f Q`, and fixed magnetic charge. Write
`G=G_N`, `r=r_e`, and the dimensionless throat coordinates as `(tau,rho)`:

\[
 ds^2=r^2[-f\,d\tau^2+d\rho^2/f+d\Omega_2^2],\quad
 f=1+\rho^2,\quad t=\ell\tau,\quad
 \alpha={Gc\over4\pi r^2},\quad\beta={Gm\over r}.
\]

The source is static, electrically and magnetically neutral for the
calculation here, angularly smeared, with fixed total *proper rest mass*
`m`. An off-center source needs a radial supporting string. The MMP
ground-state Casimir prescription and reference subtraction are imported
at their stated leading order; no evolving four-dimensional state is
derived. Require `alpha,beta<<1`, `d<<ell`, a throat/mouth overlap, and
small source size relative to `r`. Smearing describes the far field of a
localized object, not its near-zone geometry. Formal `rho -> infinity`
below extracts overlap coefficients; the throat approximation cannot
itself extend to the asymptotically flat exterior.

The decisive source is Emparan, Grado-White, Marolf and Tomasevic,
[*Multi-mouth traversable wormholes*](https://doi.org/10.1007/JHEP05(2021)032),
JHEP 05 (2021) 032, here called EGMM. Read arXiv v2 sections 2-3 and
Appendix B; compare the centered calculation in v1 and the publisher
PDF, particularly printed pages 10, 12 and 13. The publisher PDF is
available through [INSPIRE/SCOAP3](https://inspirehep.net/files/bcd01c1c0d090e073c0ad34a5f460a9b),
SHA-256 `53cac622417261aac1fc7101096c5c8a2af0dd2176acd4fe4126068bf1c92d4d`.
The equations discussed below agree between v2 and the publisher PDF;
v1 has the centered calculation but no lowering subsection.
This is a selected proof/edition comparison, not acceptance of the full
multi-mouth matching, quantum-state or causality arguments.

## 1. Invariant normalization of the held mass

Let `rho=a>=0`, `f_a=1+a^2`, and use the coordinate delta function
`delta(rho-a)`, normalized with `d rho`. The lapse, proper volume element
and unit static velocity are

\[
 N=r\sqrt f,\quad dV={r^3\over\sqrt f}\,d\rho\,d\Omega,
 \qquad u^\tau={1\over r\sqrt f}.
\]

The scalar rest-energy density and covariant stress therefore are

\[
 \epsilon={m\sqrt{f_a}\over4\pi r^3}\delta(\rho-a),\qquad
 T_{\tau\tau}^{\rm body}={m f_a^{3/2}\over4\pi r}\delta(\rho-a)
 ={\beta\over4\pi G}f_a^{3/2}\delta(\rho-a).                 \tag{1}
\]

Indeed `integral epsilon dV=m`. Independently, variation of
`S_body=-m integral ds=-m r sqrt(f_a) integral d tau` with respect
to the lapse gives the same density. The conserved Killing energy is

\[
 E_\tau=m r\sqrt{f_a},\qquad E_t={m r\sqrt{f_a}\over\ell}. \tag{2}
\]

EGMM (3.15) instead prints `beta f_a delta/(4 pi G)` while describing
`m` as the same fixed locally measured mass as at the center. With the
coordinate delta used in its ODE (3.17), that tensor integrates to
proper mass `m/sqrt(f_a)` and Killing energy `m r`. It is internally
consistent as a family with that varying mass, but not with the stated
fixed-rest-mass interpretation. A proper-distance delta convention
cannot rescue both printed equations simultaneously: changing the delta
also changes its ODE jump coefficient.

For a string of physical tension `mu` extending toward increasing rho,
its angular average has

\[
 T^\tau{}_\tau=T^\rho{}_\rho
 =-{\mu\over4\pi r^2}\Theta(\rho-a).
\]

Distributional radial conservation gives

\[
 \nabla_\nu T^\nu{}_\rho=
 {1\over4\pi r^2}\left[-\mu+{m a\over r\sqrt{f_a}}\right]
 \delta(\rho-a)=0,
 \qquad \mu={m a\over r\sqrt{f_a}}.                       \tag{3}
\]

This independently equals `m` times the proper acceleration
`a/(r sqrt(f_a))` of a static worldline. In EGMM's dimensionless tension
notation `mathcal T=8 pi G mu`, the result is
`mathcal T=8 pi beta a/sqrt(f_a)`, instead of its (3.18)
`8 pi beta a/f_a`. The fixed-mass tension approaches `m/r` deep in the
large-rho throat overlap; it does not tend to zero. This does not specify
the force measured at a winch outside the mouth or its exterior support.

## 2. Corrected off-center solution and its actual consequence

In exactly the linear throat equations used by EGMM, the replacement
`beta -> beta sqrt(f_a)` in the held-body terms gives

\[
 f\phi''+\rho\phi'-\phi
 ={\alpha\over f}-2\beta\sqrt{f_a}\delta(\rho-a)
 -{2\beta a\over\sqrt{f_a}}\Theta(\rho-a),                 \tag{4}
\]

and the radial constraint is
`rho phi' - phi = -alpha/f - 2 beta a Theta/sqrt(f_a)`.
The complete scalar solution is

\[
 \phi=\alpha(1+\rho\arctan\rho)
 -\kappa[(\rho-a)\Theta(\rho-a)-k\rho],\qquad
 \kappa={2\beta\over\sqrt{f_a}}.                          \tag{5}
\]

It is continuous, has derivative jump `-kappa`, and solves both
equations. The angular equation is
`gamma''+f phi''+2 rho phi'+4 phi=0`; its solution follows from
EGMM Appendix B by the same replacement:

\[
\begin{split}
 \gamma={}&-\alpha[(\rho^2+3)\rho\arctan\rho+\rho^2-\log f]\
 &+{\beta\over\sqrt{f_a}}\bigl[|\rho-a|f
 +2a(\rho-a)^2\Theta(a-\rho)-3a\rho^2+(1-2k)\rho^3\bigr].  \tag{6}
\end{split}
\]

The two outward sphere-growth coefficients are

\[
 A_R={\pi\alpha\over2}-\kappa(1-k),\qquad
 A_L={\pi\alpha\over2}-\kappa k.                          \tag{7}
\]

Both must be positive for the proposed two-sided throat/mouth matching.
If one retains EGMM's *assumption* `0<=k<=1/2`, this reduces to
`beta < pi alpha sqrt(f_a)/[4(1-k)]`. A height-independent sufficient
bound is still `beta < pi alpha/4`. At the symmetric center `k=1/2`,
the sharp leading condition is still `beta < pi alpha/2`, or

\[
 m<m_{\rm crit}={c\over8r}.                              \tag{8}
\]

Thus the normalization correction changes the height dependence and
local support force but does **not** invalidate the paper's centered
mass bound or its conservative factor-of-two envelope. Determining
`k(a)` requires global boundary and support data; neither the source nor
this local calculation fixes it. Positive overlap slopes are a matching
criterion, not a proof of nonlinear stability or signal arrival.

## 3. Centered energy: three different quantities

EGMM (3.8) is (8). Multiplication by its stated `r/ell` redshift gives
`E_body<c/(8 ell)`. The next printed equation (3.9) is instead
`E_bh<c/(16 ell)`. There is a factor-of-two discrepancy in that displayed
inference. It is present in v1, v2 and the publisher PDF. One should not
repair it by calling all notions of added energy the same quantity.

For a self-adjusting centered equilibrium, let `ell_0` be the unloaded
length. The leading MMP energy plus the rest-source term is

\[
 E(x;m)=A x^2+(r m-B)x,\quad x=1/\ell>0,\quad
 A={r^3\over G},\quad B={c\over8}.                         \tag{9}
\]

Here the state prescription and short exterior-loop approximation are
the same as in the unperturbed energy; the source is static and neutral.
Minimization, or independent matching of (5) at `a=0,k=1/2`, gives

\[
 \lambda={rm\over B},\quad x_0={B\over2A},\quad
 x_m=x_0(1-\lambda),\quad
 E_m=-K(1-\lambda)^2,\quad K={B^2\over4A}={c\over16\ell_0},
 \qquad 0\le\lambda<1.                                  \tag{10}
\]

For the independent matching, `phi ~ (pi alpha/2-beta)|rho|`,
`R_areal-r ~ r phi/2`, and the RN overlap relation is
`R_areal-r ~ r^2 |rho|/ell`. Thus `r/ell=(pi alpha/2-beta)/2`,
exactly as in (10). This does not compute the full exterior geometry.

The three energy quantities now differ:

| Quantity | Leading value |
| --- | --- |
| Total equilibrium excess above the unloaded state | `Delta E=K(2 lambda-lambda^2)` |
| Body energy in the loaded geometry | `E_body=m r/ell_m=2K lambda(1-lambda)` |
| Body energy evaluated at the unloaded redshift | `m r/ell_0=2K lambda` |

In particular `Delta E<K`, while `E_body<=K/2`, with the latter maximum
at `lambda=1/2`. The equilibrium energy satisfies
`dE_m/dm=r/ell_m`, consistent with reversible work. The total excess
includes geometry and vacuum readjustment; it is not the instantaneous
body term. The statement that the total equilibrium excess stays below
the original binding energy is consistent in this leading model. It
does not follow by the factor-of-two multiplication printed in EGMM.
For example, `lambda=3/4` still has positive leading overlap slopes:
`E_body=3K/8`, `c/(16 ell_m)=K/4`, and `Delta E=15K/16`.
The printed bound on the body's energy fails in this leading family,
while the bound on the total equilibrium excess remains satisfied.

Keep a fixed margin `lambda<=1-delta`, `delta>0`, when using this
truncated approximation. As `lambda->1`, the leading outward slope
vanishes and neglected terms can set the threshold. The formal limit
does not certify the actual collapse energy. A released pulse, a held
source, and a sequence of equilibrium sources require different dynamics.

## 4. Novelty and the next decisive problem

Maldacena and Qi's [v3](https://arxiv.org/html/1804.00491v3), section 4.1,
equations (4.37)-(4.43), already reduce the constrained symmetric dynamics
to `V=e^(2 varphi)-eta e^(2 Delta varphi)` and add the bulk load term
`E_bulk e^varphi/N`. At `Delta=1/2` it has exactly the quadratic-in-redshift
structure of (9), after rescaling coefficients. Appendix B's charge
formula and its symmetric reduction reproduce that equation. Sections
4.1-4.2 and Appendix B were read for this comparison; the full microscopic
derivation, other appendices, edition chain and dynamical competitors
are not all accepted. The two-boundary interaction and its global causal
setup differ from a local four-dimensional exterior return path.

Consequently (9)-(10), an elementary load-induced shift of the minimum,
or the existence of a finite static capacity are **already-known
mechanisms**, not new primary or fallback results. Equations (1)-(7)
provide a checked correction under explicit mass/delta conventions, but
no standalone new-physics paper is selected on that basis. No general
wormhole or finite-disturbance no-go has been proved.

Continue the original finite-time support question only with a specified
admissible initial state and controlled evolved stress over a transit.
If a load/transport fallback is retained, it must introduce a consequential
finite-duration effect with a conserved source and closed quantum response,
beyond the static EGMM load and the MQ potential. The actual retarded
exterior connection, its initial occupations and the common-observer
observable must be included or bounded. Final primary/fallback novelty
acceptance is still open. Extra stationary determinant or load-potential
lemmas do not close this gap.

The [independent check](../checks/static_load_normalization.py) derives
the mass, Killing energy and support force, tests every piece/jump of
the corrected throat solution, and compares the energy minimum with
the independent overlap relation and numerical minimization. The
[receipt](../evidence/static-load-normalization.json) records exact
source/script hashes and results. These checks certify the displayed
reduced-model algebra, not four-dimensional evolution. No Lean was run.

Remaining to-do list: finish the decisive dynamical proof/edition review
and primary/fallback selection; construct the finite-time support state
or a distinct controlled transport result; close physical matching and
observer comparison; then independently verify, GCP-certify suitable
complete endpoints and prepare the manuscript.
