# Long magnetic returns and the absorptive throat response

2026-10-10, E003/S051-S053. **Result in a specified channel model:** the
continuum of point-dipole return lengths has a nonzero low-frequency
absorptive response of order `q d/T^2`, although every fixed short loop
has a gap and its nonconstant lapse matrix elements vanish quadratically
as `d/T -> 0`. A fixed positive flux-label cutoff therefore removes the
leading absorption. The proof below includes the point-dipole tail and
weak spectral limit; it does not assume differentiability of its inverse
length distribution.

This is a candidate mechanism for a physical stability result, not an
accepted four-dimensional damping theorem. In particular, the initial
state, channel reduction in the shrinking tail, and constraint-compatible
collective mode still need quantitative control. A prescribed lapse is
not automatically an Einstein-equation solution. Generic Virasoro
response and driven-CFT evolution are established theory, not our novelty.

## Closed auxiliary problem and established baseline

Use units with light speed and hbar one. Each of `q` independent channels
is a massless complex Dirac field with both chiralities (`c=1` in the usual
two-dimensional convention) on an antiperiodic circle of optical length
`L=T+D`. The initial state is its NS vacuum, with zero charge/holonomy and
vacuum energy `-pi/(6L)`. `T=pi ell` is the throat optical length and `D`
is the exterior return length. The continuum convention is `q int_0^1 du`;
it is a large-channel approximation, not an additional degeneracy factor.

Perturb the relative throat lapse by

`H(t)=H0+epsilon a(t) O_h`,
`O_h=int_0^L f_h(s) T00(s) ds`.

Here `a` is real and smooth, `|epsilon| sup|a|<1`, and `f_h` is the periodic Gaussian convolution
of `1_[0,T]`, with strictly positive optical smoothing width `h`. Thus

`f_n=(1-exp(-i k_n T))/(i k_n L) exp(-h^2 k_n^2/2)`,
`k_n=2pi n/L` for `n!=0`.

The collar is part of the specified perturbation, not a regulator silently
dropped from a static susceptibility. In particular, a sharp step has
ultraviolet problems as a Hamiltonian perturbation. Fixed finite-frequency
absorption has the displayed controlled `h -> 0` limit. At positive `h`
the spatial operator is smoothly smeared. A finite-duration smooth `a`
defines the usual perturbative transition problem. Statements below concern
its second derivative at `epsilon=0`; no uniform finite-amplitude remainder
in all limits is claimed.

For one chirality, an occupied half-integer mode `r` is promoted to `r+n`.
For `n>=1`, the occupied-to-empty values are
`r=-n+1/2,...,-1/2`. The symmetrized fermion energy-density vertex is
`k_1 (r+n/2) f_n`. Directly summing squared matrix elements gives

`sum_r (r+n/2)^2 = n(n^2-1)/12`.

The two chiralities yield the positive-frequency spectral weight

`W_n=(2pi/L)^2 n(n^2-1)|f_n|^2/6`
`   =2(n^2-1) sin^2(pi n T/L) exp(-h^2 omega_n^2)/(3 L^2 n)`,
`omega_n=2pi n/L`, `n>=2`.                                      (1)

For `c` identical Dirac species multiply by `c`; hereafter formulas display
that factor, while the numerical checks use `c=1`. The `n=1` weight is
zero, as required by the vacuum's global conformal symmetry. A uniform
lapse over the whole circle has no excited-state matrix elements. These
are useful normalization and gauge checks, not new CFT identities.

With `chi_R(t)=-i theta(t)<[O_h(t),O_h(0)]>` and Fourier transform using
`exp(i omega t)`, the convention is

`chi_R(omega)=sum_n W_n[(omega-omega_n+i0)^-1-(omega+omega_n+i0)^-1]`,
`rho(omega)=-Im chi_R(omega)=pi sum_n W_n delta(omega-omega_n)>0`
for positive frequency. For a pulse returning to zero, the coefficient
of `epsilon^2` in its final excitation energy is

`E^(2)=sum_n omega_n W_n |a_hat(omega_n)|^2`
`     =pi^-1 int_0^infinity omega rho(omega)|a_hat(omega)|^2 d omega`. (2)

This is energy supplied by a prescribed drive. It is not yet an energy
loss theorem for a self-consistent gravitational oscillator.

## Point-dipole tail, including the endpoint

Use the flux label `u=nu/2` of the established MMP point-dipole exterior.
The [earlier delay note](dipole-delay-bound.md) defines the same geometry.
With prolate coordinate `eta in [-1,1]`, set

`a=1-eta^2`, `Q=sqrt(a^2+4u^2 eta^2)`, `mu=(a+Q)/(2u)`.

The meridional curve has
`z/d=mu eta/2`, `r/d=sqrt((mu^2-1)(1-eta^2))/2`. Differentiation gives

`4 (ds/deta)^2/d^2=mu/u+4 mu eta^2(mu-u)/Q^2`.

Thus `D=d f(u)` and

`u f(u)=int_0^1 sqrt[u mu+4u^2 mu eta^2(mu-u)/Q^2] d eta`.       (3)

For `0<u<1`, `mu>=1`, `u mu<=1` and `u mu<=Q`. The first upper bound
follows by evaluating `u mu^2-a mu-u eta^2` at `mu=1/u`; the polynomial
is increasing above its positive root. The second follows from `Q>=a`.
Consequently the integrand in (3) is at most `sqrt(1+4eta^2)<=sqrt(5)`.
For each interior `eta`, it tends to `sqrt(1+3eta^2)`. Dominated
convergence, including the endpoint region, proves

`lim_(u->0) u f(u)=kappa=1+asinh(sqrt(3))/(2sqrt(3))`,
`0<u f(u)<=sqrt(5)`.                                         (4)

This proof does not require a numerical endpoint fit, or the unproved
derivative assertion `f'(u)~-kappa/u^2`. An independent arclength
quadrature differentiates the Cartesian cylindrical coordinates instead
of using (3).

## Exact inverse-tail response and the point-dipole limit

First consider the exactly defined inverse-tail model `D=kappa d/u` on
`0<u<1`. At a fixed positive frequency with
`omega(T+kappa d)<4pi`, every `n>=2` has a root

`L_n=2pi n/omega`, `D_n=L_n-T`, `u_n=kappa d/D_n`,
`|du/domega_n|=kappa d L_n^2/(2pi n D_n^2)`.

Substituting (1), including the factor `pi` in `rho`, proves

`rho_inv(omega)=(q c kappa d/3) sin^2(omega T/2) exp(-h^2 omega^2)`
`                 * sum_(n=2)^infinity (1-1/n^2)/(2pi n/omega-T)^2`. (5)

The series converges absolutely. It is positive except at the zeros of
the sine. On compact frequency intervals within the stated band it is
exactly linear in `d`. It has no positive spectral gap.

For the actual point-dipole `D=d f(u)`, equation (5) is the leading
**weak spectral asymptotic**, rather than an asserted exact density.
Precisely, let `g` be continuous with compact support in `(0,4pi/T)`.
Keep `T,h>0` fixed and let `d -> 0`. Then

`lim (q d)^-1 int g(omega) rho_d(d omega)`
` = int g(omega) [rho_inv(omega)/(q d)] d omega`.               (6)

Here `rho_d(d omega)` denotes a spectral measure, to avoid assuming a
differentiable length inverse. A proof which also handles the infinitely
long loops is short. Write
`F(D)=pi sum_(n>=2) W_n(T+D) g(2pi n/(T+D))`.
Let `omega_max<4pi/T` contain the support and
`Delta=4pi/omega_max-T>0`. Then `F(D)=0` for `D<Delta`. The bound
`W_n<=2c n/(3L^2)` and the number of terms `n<=omega_max L/(2pi)`
bound `F` uniformly in `L>=T`. It is continuous for finite `L`, since
`g` extended by zero is continuous and only finitely many terms occur
locally. In the integral substitute `u=d v`. By (4), its integrand
vanishes when `v>sqrt(5)/Delta`, while
`d f(d v) -> kappa/v` for each `v>0`. Dominated convergence on that
fixed finite interval proves (6). Applying the exact inverse-tail change
of variables then gives (5). No interchange of an uncontrolled length
derivative or ultraviolet series is involved.

The theorem includes nonnegative windows whose limiting integral is
strictly positive. It therefore disproves a uniform quadratic bound on
this band-integrated absorption. With a fixed trim `u>=u0>0`, `D`
is uniformly small as `d->0`, so the same below-gap window has zero
response for sufficiently small `d`. That trim discards the entire
leading effect. A mode-by-mode small-`D` expansion cannot be averaged
without controlling the shrinking flux sector `u=O(d/T)`.

Two useful consequences of (5), with the limit taken in the order just
specified, are

`rho_inv(omega)~q c kappa d T^2 (15-pi^2) omega^4/4320`
as `omega T->0`, and

`rho_inv(pi/T)=q c kappa d exp(-h^2 pi^2/T^2) S/(3 T^2)`,
`S=sum_(n=2)^infinity (n^2-1)/(n^2(2n-1)^2)`
` =8 log(2)-13 pi^2/24 =0.1991417272...`.                       (7)

For the last identity, the summand is
`8/(2n-1)-3/(2n-1)^2-4/n-1/n^2`; pair the harmonic differences
before summing. Neither (6) nor the numerics alone proves pointwise
density convergence for the actual dipole, or interchange with the
`omega->0` limit.

## Metric vertex and the remaining gravitational question

There is a genuine metric realization of the **prescribed** drive. In
`ds_2^2=Omega(s,t)^2[-v(s,t)^2 dt^2+ds^2]`, canonical Weyl rescaling
of a massless Dirac field gives the two one-particle generators
`h_+ = {v,-i partial_s}/2` and `h_-=-h_+`. For
`v=1+epsilon a f_h`, the first variation is exactly `O_h` above.
The `v'/2` term explains the mean momentum `r+n/2` in the independent
fermion calculation. Renormalization and the conformal anomaly contribute
state-independent terms; those do not change excited-state transition
matrix elements. They must still be restored for stress one-point
functions, conservative forces, and Einstein constraints.

In the leading constant-radius throat ansatz
`ds_2^2=r_e^2[-(1+rho^2)e^(2phi(t))dt^2+drho^2/(1+rho^2)]`,
take `phi0=-log ell` and the fixed background optical coordinate
`s=ell(arctan rho+pi/2)`. Then `v=e^(phi-phi0)` in the throat.
Holding the exterior laboratory lapse fixed gives the relative vertex
above, with a chosen matching collar. This verifies that vertex for the
specified ansatz. It does **not** show that this ansatz and the chosen
collar solve the four-dimensional time-dependent constraints, nor that
additional mouth/gauge variations have negligible transition matrix
elements. A time reparametrization of the whole loop has zero absorption,
as (1) correctly checks.

One can state the prospective damping test without promoting it to a
physical conclusion. **If** the physical collective action has the
leading quadratic form `A(delta_phi_dot^2-ell^-2 delta_phi^2)`, with
`A=q c ell/16`, and its matter vertex is this unit relative lapse, then
the absorptive pole shift would give amplitude rate
`Gamma=rho(1/ell)/(4A/ell)`. In the inverse-tail continuum, (7) would give

`Gamma/omega0=[4 kappa S/(3pi^2)] (d/ell) exp(-h^2/ell^2)`,
`omega0=1/ell`.                                               (8)

The coefficient is evaluated in the receipt. Equation (8) is a
conditional estimate, not a verified mode frequency or stability theorem.
The known static MMP/MQ potential does not by itself fix the kinetic term.
Static finite-return corrections can also shift the real frequency.
The required action and source variation must be derived, with their
contact terms and constraints, before using a damping pole.

At finite channel number the spectrum is discrete. Around resonant
`u=O(d/ell)`, neighboring channel lines have spacing of order `1/(q d)`.
Any continuum decay interpretation consequently needs times much shorter
than `q d`; indefinite irreversible decay of a closed finite system does
not follow. A prospective decay time `ell^2/d` fits this window only if
`q(d/ell)^2>>1`. These are necessary scale comparisons, not an error bound.

For magnetic charges with `d/r_e=q^nu` and `ell/r_e` proportional to `q`,
resonant flux labels satisfy `q u` of order `q^nu`, even though `u->0`.
Thus the earlier fixed-trim approximation is inapplicable and a uniform
shrinking-sector reduction is the next hard calculation. Rotation, the
finite-time support state, outer returns and four-dimensional mode leakage
must be compared with the order-`d/ell` observable, not merely shown small
in absolute terms. No causal-shortcut conclusion is involved.

## Relevant predecessors and reproducibility

- [Maldacena--Milekhin--Popov, arXiv:1807.04726v3](https://arxiv.org/abs/1807.04726v3):
  the magnetic return geometry and fermion Casimir baseline are established.
  The earlier source and edition audits remain in force. We do not claim
  the long field lines or their guiding-center description are new.
- [Moosavi, arXiv:1912.04821v3](https://arxiv.org/abs/1912.04821v3):
  inhomogeneous-CFT Hamiltonians and diffeomorphism evolution are established;
  section 2.1's algebra/transformation was read. Full dependency acceptance
  is not inferred from that section.
- [Driven inhomogeneous CFT, arXiv:2508.18350v2](https://arxiv.org/abs/2508.18350v2):
  section 3.1 was read in full, including the Hamiltonian pictures,
  scheme-dependent c-numbers and the Virasoro representation argument.
  This prevents identifying a chirally split operator expectation with
  a covariantly renormalized gravitational stress without its local terms.
  Version 2 corrects signs and stress expressions from version 1; the
  latter was retrieved but is not used for the response formula. The full
  Ward-identity appendices and publisher-edition comparison remain open.
- [Floquet SYK wormholes, arXiv:2404.08394v3](https://arxiv.org/abs/2404.08394v3):
  its Schwarzian forced oscillator and resonance/heating discussion are
  competing established dynamics. Section 4 was read fully in v2; the
  current v3 section has identical extracted text and equations. Sections
  2-3 were inspected in part in v2 and are not accepted as a full current
  proof review. Generic resonance is not our candidate novelty.
  Its coupled-boundary drive is not a distribution of exterior path lengths.
- [How traversable is a traversable wormhole?, arXiv:2606.12528v2](https://arxiv.org/abs/2606.12528v2):
  section 5 was read fully in addition to the earlier fermion-transmission
  review. It explicitly identifies collective-mode response and Schwarzian
  corrections near resonances as further work. Its one-sided scattering
  experiment excludes participation of the other mouth; it does not close
  the vacuum-return problem studied here. Its expected gapped collective
  phase is not a theorem excluding this matter continuum. This is a close
  active research direction and requires a refreshed citation search before
  accepting novelty or manuscript scope.

The possible contribution is the nonuniform long-return contribution to
the physical quantum response, if the remaining model and novelty checks
survive. Supplemental searches have not established absence of prior work.
The planning audit records discovery and reading separately.

Run `python checks/long_return_response.py`. The
[receipt](../evidence/long-return-response.json) checks exact particle-hole
normalization, Fourier quadratures, two independent geometric arclength
formulas, the Jacobian and closed coefficients, direct loop integration
against frequency integration, a point-dipole weak-limit sequence, and
finite-channel convergence. The frequency sum has an explicit omitted-mode
bound; quadrature refinement is numerical evidence, not formal interval
certification. These checks support the ordinary proof, not the unresolved
four-dimensional interpretation. No Lean or cloud process was invoked.

Remaining to-do list: derive the constraint-compatible collective action
and total vertex, control the shrinking long-return sector and its initial
state, finish the closest-prior-work review, then decide whether this
mechanism supplies the primary physical stability result or a sharply
scoped finite-duration transport fallback.
