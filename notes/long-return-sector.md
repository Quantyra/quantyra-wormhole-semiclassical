# Propagation in the shrinking long-return sector

2026-10-10, E003/S051-S053. The flux labels responsible for the proposed
absorptive clock response shrink as `u=O(d/ell)`. Fixed-trim estimates do
not apply there. Rescaling the actual point-dipole operator shows that
the relevant large parameter is **`Q=q u`**, rather than `q` alone.
The corrected packet construction gives a uniform local propagation
error `O(Q^(-1/2))` on a fixed portion of the rescaled long curve.

This supplies a previously missing local estimate for that sector. It
does not identify the entire physical return, its vacuum, or a damping
rate in four dimensions. The endpoint portions of the curve still need
matching. In particular, a local residual cannot be evolved repeatedly
around an unconstructed global loop.

## Exact rescaling and a regular limiting curve

Use the flat, axisymmetric two-component Weyl reduction already derived
in [the exterior operator note](exterior-dirac-window.md). In units of
mouth separation `d`, write

`U(r,z)=[(z+1/2)/sqrt(r^2+(z+1/2)^2)`
`        -(z-1/2)/sqrt(r^2+(z-1/2)^2)]/2`,
`H=-i sigma_x partial_r-i sigma_z partial_z+q sigma_y(U-u)/r`.

The azimuthal sector satisfies `q(1-u)` half-integral in this convention.
Restore dimensions by multiplying this displayed operator by `1/d`.
For `R=u r`, `Z=u z`, the unitary dilation of the two-dimensional flat
measure gives the exact identity

`(d/u) H_phys = -i sigma_x partial_R-i sigma_z partial_Z`
`                  +Q sigma_y M_u(R,Z)`,
`Q=q u`, `M_u=(V_u-1)/R`, `V_u=U(R/u,Z/u)/u`.                 (1)

There is a nonsingular representation, including `u=0`:

`V_u(R,Z)=(1/2) int_(-1/2)^(1/2)`
`               R^2/[R^2+(Z+u t)^2]^(3/2) dt`,
`V_0=R^2/[2(R^2+Z^2)^(3/2)]`.                               (2)

On every compact set away from the axis and origin, (2) proves
`V_u=V_0+O(u^2)` in every fixed `C^k` norm. Differentiate under the
integral; the term odd in `t` vanishes and the Taylor remainder is
uniformly bounded. The estimate is for the *rescaled* coefficients.
Replacing `M_u` by `M_0` inside (1) would introduce `O(Q u^2)`;
we keep `M_u` exactly and do not assume that error is small.

An exact parametrization of `V_u=1`, for `-1<eta<1`, is

`a=1-eta^2`, `b=(a+sqrt(a^2+4u^2 eta^2))/2`,
`R_u=(1/2)sqrt((b^2-u^2)(1-eta^2))`, `Z_u=b eta/2`.

For any fixed `eta_*<1`, restrict to `|eta|<=eta_*` and take `u`
sufficiently small. The curves converge in every fixed `C^k` norm to

`R_0=(1-eta^2)^(3/2)/2`, `Z_0=eta(1-eta^2)/2`.              (3)

Their speed never vanishes on this compact interval. The normal mass
slope `b_u=|grad V_u|/R_u` is bounded above and below by positive
constants there. At `u=0` this follows directly from (2): the radial
and vertical derivatives cannot both vanish at `V_0=1`. Compactness
and smooth convergence give the uniform result. A slightly larger
compact curve supplies a fixed embedded normal tube, with uniform
curvature and mass-derivative bounds. No claim extends this tube to
`eta=+-1`, where (3) reaches the excluded origin.

## A uniform local packet theorem

Let `s` be rescaled arclength on that larger curve, `y` its signed normal
coordinate, `h=1-k(s)y`, and choose the normal so `b_u=M_{u,y}(s,0)>0`.
The exact spin and measure transformations give

`H_F=i sigma_x partial_y-i sigma_z h^(-1/2) partial_s h^(-1/2)`
`                                           +Q sigma_y M_u(s,y)`.

Define `c_u=M_{u,yy}(s,0)`, `a_u=(partial_s b_u)/(4b_u)`,

`g=(Q b_u/pi)^(1/4) exp(-Q b_u y^2/2)`,
`P=1-Q c_u y^3/6`, `J=a_u y`,
`Psi_app(tau,s,y)=F(s-tau) chi(y) g(s,y) (P,J)`.             (4)

Here `F` is a fixed smooth compactly supported envelope, `chi` is even, one
near zero and vanishes inside the fixed tube boundary, and the moving
support remains inside the larger curve during `0<=tau<=tau_*`.
All these choices are independent of small `u` and large `Q`.
The exact Weyl field is obtained by undoing the spin/measure maps.

**Local theorem.** There are constants `C,c>0`, independent of small
`u` and `Q>=1`, such that

`||(i partial_tau-H_F)Psi_app|| <= C Q^(-1/2)+C exp(-c Q)`,
`||Psi_app||^2=||F||^2+O(Q^(-1))`.                           (5)

The constants depend on the compact curve, envelope and time interval.
For a self-adjoint Dirac evolution agreeing with this operator on the
causal neighborhood of the packet, the exact solution with initial
data (4) satisfies

`||Psi_exact(tau)-Psi_app(tau)||`
`       <= C tau [Q^(-1/2)+exp(-c Q)]`.                      (6)

One may realize this local comparison by extending the smooth real
mass in the original Cartesian patch to a bounded smooth function on
the plane. The free planar Dirac operator plus that bounded Hermitian
potential is self-adjoint on its usual first Sobolev domain. Its speed
bound is one, so extensions agreeing on the causal neighborhood yield
the same packet evolution there. This construction is a local PDE
comparison, not a replacement of the wormhole exterior by a boundary
condition or a quantum-state prescription.

**Proof.** Use the exact residual identities (5)-(8) in
[the interior trial argument](exterior-trial-residual.md), with `q`
replaced by `Q` and the uniformly bounded coefficients just constructed.
For completeness, the width derivative of the leading Gaussian is
`a_u(1-2Qb_u y^2)g`. The lower component `a_u y g` cancels it under
`i sigma_x partial_y+sigma_y Qb_u y`. The upper cubic correction
cancels the quadratic normal mass source. The remaining terms have
Gaussian norms `O(Q^(-1/2))`. For example,
`||Q^p y^j g||=O(Q^(p-j/2))`, and each remaining monomial has
`j>=2p+1`. Curvature denominators obey `h>=1/2` on the tube.

For a nonconstant envelope, `partial_tau F=-partial_s F` cancels the
leading upper longitudinal term. The extra upper remainder contains
`F'(h^(-1)-1)P`, and the lower one contains
`F'(h^(-1)+1)J`; both are `O(Q^(-1/2))`. Thus cutting off a moving
envelope does not incur the order-one error of a stationary longitudinal
cutoff. The transverse cutoff contributes the exponential term as in
the earlier proof. In the norm, the cubic cross term is odd in `y`;
the squared cubic and lower corrections are `O(Q^(-1))`, with only
exponentially small symmetric-cutoff tails. This proves (5).
Variation of constants and unitarity integrate its residual and give
(6). The estimate is uniform in `u` because every geometric constant
was obtained on the same compact rescaled family.

Restoring units, `tau=t_phys u/d`. On a passage of duration `O(d/u)`,
the local relative norm error is `O((q u)^(-1/2))`. This proves a
local propagated-packet estimate; it is stronger than a localization
criterion but does not prove an operator-norm return estimate.

## Resonant labels, omitted harmonics and the time budget

Let `delta=d/ell`, `T=pi ell`, and test frequency `omega=1/ell`.
The inverse-tail response has roots

`u_n=kappa delta/[pi(2n-1)]`, `n>=2`,
`kappa=1+asinh(sqrt(3))/(2sqrt(3))`.

For any fixed `2<=n<=N`, the actual dipole roots have the same leading
scaling, by the derivative result in
[the collective-clock note](collective-clock-response.md). Thus
`Q_n` is of order `q delta/n`, and the rescaled frequency `omega d/u_n`
is of order `n`. Equations (5)-(6) apply uniformly to every *fixed*
finite collection of these harmonics as `q delta->infinity`.

In the independent-channel model their positive absorption weights are
proportional to

`s_n=(n^2-1)/[n^2(2n-1)^2]`,
`S=sum_(n>=2)s_n=8log(2)-13pi^2/24`.

For `N>=2`, the integral comparison gives the exact bound

`0<sum_(n>N)s_n <= int_N^infinity dx/(2x-1)^2`
`                         =1/[2(2N-1)]`.                    (7)

It justifies retaining a finite number of harmonics at any fixed
relative accuracy *in that model*. It does not bound the response of
discarded modes of the four-dimensional Dirac problem. The required
order of limits is to choose `N` from (7), then make the local estimates
small at that `N`; the constants need not be uniform as `N->infinity`.

The model width is `Gamma=C_G delta/ell`,
`C_G=4kappa S/(3pi^2)>0`. A claimed physical damping law would need
control for times `t_G=ell/(C_G delta)`, much longer than one transit.
The following are scale requirements, with fixed positive numerical
constants suppressed; they are **not** additional proved error bounds.

| Issue | Ratio that must vanish | Necessary/sufficient meaning |
| --- | --- | --- |
| Nearby finite-channel lines | `t_G/(q d) ~ 1/(q delta^2)` | Necessary continuum time resolution at the low harmonics, not a recurrence theorem. |
| Accumulating a hypothetical global residual of the order in (5) | `1/[delta sqrt(q delta)]` | Sufficient by unitarity *if* matching supplies the same global residual throughout that time. |
| Frozen orbital orientation over `t_G` | `Omega ell/delta` | A useful sufficient frozen-geometry condition; not required by every possible rotating treatment. |
| Both nominal throat-overlap errors below `delta` | `sqrt(r_e/ell)/delta` | Choose optical cutoff `epsilon=sqrt(r_e ell)`; actual matching error estimates remain required. |

Using the established MMP leading scales
`ell/r_e=16pi q/g^2`, `Omega=2 sqrt(r_e/d^3)`, and
`d/r_e=q^nu`, at fixed small `g`, the four ratios have powers

`q^(1-2nu)`, `q^(1-3nu/2)`, `q^(2-5nu/2)`, `q^(1/2-nu)`.

They all tend to zero when `4/5<nu<1`. For the previous candidate
`nu=15/16`, these powers are respectively
`-7/8, -13/32, -11/32, -7/16`. Large constants, running couplings,
and interaction errors are not controlled by these power comparisons.
The radiation lifetime scale `d^3/r_e^2` also exceeds `t_G` in this
range, with ratio `t_G/t_life` of order `q^(2-4nu)`.

This rules out a simple incompatibility of those necessary scales.
It does not establish a finite admissible parameter point. Fixed small
`g` is a leading-model comparison, not an unbounded controlled limit
of a specified interacting four-dimensional theory.

## Predecessors, checks and decision

[MMP v3](https://arxiv.org/abs/1807.04726v3), Appendix E (E.22)-(E.27),
already gives the long-line localization criterion `q nu_MMP >>1`;
our label is `u=nu_MMP/2`. The effective parameter `q u` is therefore
not a new physical discovery. The source's estimate temporarily omits
spin; (1)-(6) retain the Weyl spin and geometry and give the local error
needed by this particular response calculation. Sections 6.1.1-6.1.3
also supply the orbital and radiation limitations. Those passages were
read again in full. Their estimates are not a proof of the proposed
long-return damping law.

The curved-interface packet method remains established prior work,
including [Bal et al. v1](https://arxiv.org/abs/2106.00729v1). No claim
of first adiabatic propagation or first dissipative Schwarzian model
is available. Eight supplemental search queries on collective damping,
magnetic return lines and Schwarzian baths did not settle novelty;
they are discovery evidence only, with no completed two-pass review exit.

Run `python checks/long_return_sector.py --output
evidence/long-return-sector.json`. It compares two flux representations,
the exact and limiting geometry, direct Pauli residuals and their
expanded cancellations at three physical sector labels, normal
quadrature and longitudinal derivative refinements, the absorption
tail, and the exact rational scaling exponents. These checks supplement
the ordinary local theorem. They do not simulate a full return or test
a Hadamard state. The earlier collective-clock script was also rerun.
The first global Hermite quadrature did not resolve the transverse
cutoff transition; its failed input is preserved outside Git. Splitting
the integration at the cutoff breakpoints resolves the refinement test
without loosening its tolerance. Low-flux cases with a large cutoff
residual are retained, followed by four flux refinements. The data
therefore do not assert a small error throughout the tested range.

Decision: continue the physical response test, with one primary question
(does long-return absorption damp the symmetric throat clock in a
controlled MMP regime?) and one fallback investigation (finite-duration
band-resolved quantum energy transfer for a matched physical pulse).
Both remain provisional pending model closure and final novelty review.
The local estimate removes the fixed-trim obstacle for bounded rescaled
arcs; the next hard step is matching these arcs through the endpoint
regions to the throat, while controlling the state and total vertex.

Remaining to-do list: close shrinking-sector endpoint and throat matching,
construct the physical finite-time state and bound interaction/local terms,
complete primary/fallback novelty selection and the consequential result,
then independently verify, certify suitable endpoints on GCP and prepare
the reproducible manuscript.
