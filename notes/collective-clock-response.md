# Constrained throat clock and the complete channel response

2026-10-10, E003/S051-S053. The leading classical throat reduction fixes
the symmetric clock's kinetic coefficient. Its zero-charge constraint
couples that clock to the fermion energy, with the same unit relative-lapse
vertex used in the [long-return calculation](long-return-response.md).
The subtracted fermion response has a real dispersive part as well as
absorption. Dropping it would give an incomplete mode calculation.

Within the explicitly closed inverse-tail JT/channel model below, the
retarded resonance has both a calculable logarithmic frequency shift and
the previously conditional order-`d/ell` width. The actual point-dipole
absorption also has the same leading coefficient: a derivative estimate
here upgrades the previous weak spectral limit to pointwise convergence.
These results do not bound the four-dimensional state, matching or
projection errors. They are not yet a physical MMP damping theorem or an
accepted novelty claim. Classical Schwarzian dynamics and CFT correlators
are established predecessors.

## Gravitational normalization from the action

Start with the repository's fixed magnetic Einstein--Maxwell convention:
`S4=int sqrt(-g4)[R4/(16pi G)-F^2/(4g^2)]`,
`F_theta,phi=(q/2) sin(theta)`, and `r_e^2=pi G q^2/g^2`.
For `ds4^2=g_ab dx^a dx^b+R(x)^2 dOmega2^2`, spherical integration,
one integration by parts, and the Gibbons--Hawking term give

`S2=(4G)^-1 int sqrt(-g)[R^2 R2+2(partial R)^2+2-2r_e^2/R^2]`
`    +(2G)^-1 int_boundary sqrt(-gamma) R^2 K`.                (1)

The integrated normal derivative of `R` cancels its counterpart in the
four-dimensional boundary term. This reduction keeps fixed magnetic flux;
it does not import the electric-charge ensemble without its boundary term.
With `Phi=(R^2-r_e^2)/(2G)`, the linear non-topological action is

`S_JT=(1/2) int sqrt(-g) Phi(R2+2/r_e^2)`
`      +int_boundary sqrt(-gamma) Phi(K-1/r_e)`.               (2)

Terms quadratic in the dilaton fluctuation remain outside this leading
action. In the extremal overlap use `z=r_e^2/(R-r_e)` and the exterior
time `t`. Then `ds2^2=r_e^2(-dt^2+dz^2)/z^2` and
`Phi=A/z+...`, with **`A=r_e^3/G`**.

For a timelike cutoff curve `(F(t),z(t))`, fixed induced metric
`gamma_tt=-r_e^2/epsilon^2` requires

`z=epsilon F'-(epsilon^3/2) F''^2/F'+O(epsilon^5)`.

Its outward extrinsic curvature is

`K=[F'(F'^2-z'^2-z z'')+z z' F'']`
`   /[r_e(F'^2-z'^2)^(3/2)]`
` =r_e^-1[1-epsilon^2 Sch(F,t)+O(epsilon^4)]`.                 (3)

The minus sign is Lorentzian. Substitution in (2) gives one boundary
action `-A int Sch(F,t)dt`. There are two ends. Write their common global
AdS time as `tau(t)`, with `F=tan(tau/2)` and `x=tau'=exp(phi)>0`.
Using `Sch(tan(tau/2),t)=phi''-phi'^2/2+x^2/2`, the symmetric action is

`S_g=int A[phi'^2-exp(2phi)]dt`, up to an endpoint term.       (4)

As an independent normalization check, the exact RN relations
`r_+ r_-=r_e^2`, `M=(r_++r_-)/(2G)` and
`Theta=(r_+-r_-)/(4pi r_+^2)` give
`M-M_ext=2pi^2 A Theta^2+o(Theta^2)` for one hole. This fixes the
same coefficient without using the off-shell wormhole energy potential.

Equations (1)-(4) reproduce the established leading near-extremal theory.
[MSY, 1606.01857v2, sections 3 and 5](https://arxiv.org/abs/1606.01857v2)
provide the boundary and charge framework.
[Iliesiu--Turiaci, 2003.02860v3, sections 2.2, 2.4 and 3.1](https://arxiv.org/abs/2003.02860v3)
give the matching normalization and its approximation limits. Those
sections were read completely. Their mass-gap/partition-function result
is not used as a wormhole state construction. Version 3 explicitly
corrects earlier statements about the gluing procedure.
Their section 2.2 gauge-mode expansion assumes zero magnetic flux
(footnote 9). We use it as a normalization comparator for the RN
Schwarzian coefficient, not as the magnetic reduction: (1) derives
that reduction at fixed magnetic flux directly.

## The charge constraint must survive the change of variable

`phi=log tau'` is not an unconstrained replacement of the original clock.
For matter with canonical action `L_m=i psi^dagger psi'-H(phi,psi)`,
introduce a multiplier `P` and use

`L=A phi'^2-A exp(2phi)+L_m+P[tau'-exp(phi)]`.                (5)

The `tau` variation gives `P'=0`; the remaining total AdS time-translation
charge must be set to zero. Reflection symmetry removes the other two
charges in this symmetric sector. In the original higher-derivative
variables the same charge is

`P=-exp(-phi)[2A(phi''+exp(2phi))+H_phi]`.                    (6)

Thus `P=0` gives the positive-kinetic-energy reduced system

`2A(phi''+exp(2phi))+<H_phi>=0`,
`i psi'=H psi`.                                             (7)

It conserves
`E=A phi'^2+A exp(2phi)+<H>`.
Indeed the field equation implies `<H>'=<H_phi> phi'`; substituting (7)
sets `E'` to zero. The original clock Hamiltonian differs by `P exp(phi)`
and agrees when the constraint is imposed. Dropping the charge before
varying the clock would miss this condition. In the zero-exterior-length
vacuum, `H_vac=-B exp(phi)`, `B=q c/8`. Equations (7) give
`x0=B/(2A)` and `omega0=x0`, with quadratic coefficient `A`.

This is the dimension-one-half MQ baseline, not a new oscillator.
[Maldacena--Qi v3, Appendix B](https://arxiv.org/abs/1804.00491v3) was
re-read completely: its equal-clock sector has precisely the remaining
charge constraint. The opposite signs of the two ends' nonconstant
SL(2) transformations also matter. A common oscillating clock is not that
gauge transformation. For a finite exterior the matter contribution to
(6) must be the full dynamical one, not its static vacuum replacement.

## An observable check of the vertex and its collar sensitivity

In the leading throat, a radial null ray satisfies
`tau(t_out)-tau(t_in)=pi`. For
`tau'=x0[1+a(t)]`, its first-order arrival-time change, measured in the
same exterior time, is

`delta t_out=-int_(t_in)^(t_in+T) a(t)dt`, `T=pi/x0`.        (8)

A Fourier component has transfer factor
`-(1-exp(-i omega T))/(i omega)`. This independently recovers the
integrated-lapse form factor in the fermion Hamiltonian. A simultaneous
coordinate change of the whole circle changes no transition probability;
the exterior clock is held fixed in (8).

At a finite overlap cutoff the two ends satisfy
`sigma_R=pi/2-epsilon tau'+...`,
`sigma_L=-pi/2+epsilon tau'+...` in global conformal coordinates.
The missing optical pieces and their first variations are of order
`epsilon`. They cannot simply be omitted from a precision coefficient.
There is, however, a useful explicit robustness statement. Suppose the
two chiral lapse/shift profiles differ from the ideal interval profile by
functions whose `L1` norms are at most `C epsilon` at the frequency being
tested. Then their normalized Fourier errors are at most `C epsilon/L`.
For the ideal coefficient
`|f_n|=2|sin(omega T/2)|/(|omega|L)`, the relative error is bounded by

`r=C epsilon |omega|/[2|sin(omega T/2)|]`.

Each chiral spectral weight, and their sum, consequently lies between
`(1-r)^2` and `(1+r)^2` times its ideal value when `r<1`. This bound is
uniform in the return length and mode number at fixed frequency. Near
`omega T=pi` it tends to zero with `epsilon/T`. It rules out cancellation
by corrections within this specified collar class. It does not prove
that the actual four-dimensional matched solution belongs to that class.

The overlap conditions are `r_e << epsilon << ell`, with relative
classical throat error of order `r_e/epsilon`. To make both nominal
matching scales smaller than `d/ell` requires
`r_e ell/d << epsilon << d`, hence `d^2 >> r_e ell`.
This is a compatibility condition, not a bound on every omitted term.
Covariant anomaly/contact terms and the corrected reference geometry must
still be included; their conservative contribution cannot be inferred
from excited-state matrix elements.

## The full subtracted free-channel kernel

For one NS Dirac circle set `L=T+D`, `alpha=T/L` and
`z=omega L/(2pi)`. Use the weights and retarded convention of the previous
note. The difference `Delta chi=chi_R(omega)-chi_R(0)` has a convergent
sharp-profile limit even though the two unsubtracted static expressions
are ultraviolet-sensitive. Direct summation gives

`Delta chi=-c omega^2 L/(6pi^3)`
`  *sum_(n=1)^infinity [(n^2-1)/(n^2(n^2-z^2))] sin^2(pi n alpha)`
` =-c pi T D/(3L^3)`
`  -c[omega-4pi^2/(omega L^2)]/(6pi)`
`    *sin(omega T/2) sin(omega D/2)/sin(omega L/2)`.            (9)

Initially take `Im omega>0`; boundary values are retarded. At the apparent
`n=1` pole the prefactor cancels it. Real poles with `n>=2` have the
positive spectral weights already obtained from fermions. To check (9),
split its rational summand into
`1/(z^2 n^2)+(1-1/z^2)/(n^2-z^2)` and use the circle Green-function
sum, or solve the periodic Helmholtz Green equation directly. The two
needed sums are
`sum sin^2(pi n alpha)/n^2=pi^2 alpha(1-alpha)/2` and
`sum sin^2(pi n alpha)/(n^2-z^2)`
` =pi sin(pi z alpha)sin(pi z(1-alpha))/(2z sin(pi z))`.

In a compact complex frequency neighborhood of `omega T=pi`, expansion at
`D=0` gives the particularly simple result

`Delta chi(D,omega)=-c omega^2 D/(12pi)+O(D^2)`.              (10)

The cancellation of the apparent `D/T^2` pieces in (9) is essential.
For the inverse-tail distribution `D=kappa d/u`, averaging is exactly
`q kappa d int_(kappa d)^infinity Delta chi(D,omega)/D^2 dD`.
Splitting this integral at a fixed small positive `D0` proves

`Delta Chi=-q c kappa d omega^2/(12pi) log[T/(kappa d)]`
`           +O(q c d/T^2)`.                                (11)

For complex frequency the remainder is analytic on the upper half-plane.
Near `omega T=pi` it has the usual local retarded continuation across the
real cut. To justify the real principal value, write
`sin A sin B/sin(A+B)=sin A cos A-sin^2 A cot(A+B)`.
Subtract the smooth coefficient at each cotangent pole. The odd constant
term integrates to zero over a complete period; the remaining cell
integrals are bounded by a summable derivative bound. The positive
frequency density is analytic locally, from the absolutely convergent
series in the previous note, so the Cauchy-transform continuation exists.
No pole is discarded by taking a real-valued integral.

The imaginary part in (11)'s remainder is of order `d`, rather than
`d log d`, and equals minus the earlier spectral density. The logarithm
is a real inertia correction. A finite smooth collar has its own local
matching terms; the sharp, subtracted expression defines this reduced
kernel and is not a prescription for an unsmoothed full Hamiltonian.

## Closed inverse-tail model: the complex mode asymptotic

Specify the quadratic closure explicitly: the static energy is the MMP
length energy for independent inverse-tail loops, and the dynamic kernel
is (9) with its static part subtracted. Possible additional local kinetic
terms are excluded from this reduced model, not proved absent in 4D.
With `x=exp(phi)`,

`E_s(x)=A x^2-q c x/8`
`       +q c kappa d x^2/(6pi) log[1+pi/(kappa d x)]`.        (12)

This follows by integrating `1/(pi/x+kappa d/u)` over `u`, including
the `+q c x/24` throat anomaly term. Let `x_*` minimize (12),
`x0=q c/(16A)`, `delta=d x0`, and take `delta->0`. Differentiation gives

`x_*/x0=1-[8kappa/(3pi)] delta log(1/delta)+O(delta)`.

The retarded linear kernel at that equilibrium is

`K(omega)=-2A omega^2+x_*^2 E_s''(x_*)+Delta Chi(omega)`,
`T=pi/x_*`.                                                 (13)

The subtraction avoids counting the static susceptibility twice. The
coefficient in (11) contributes positive inertia. Combining it with
(12), the resonance on the retarded continuation has

`Re omega_p/x0=1-[5kappa/(3pi)] delta log(1/delta)+O(delta)`,
`-Im omega_p/x0=[4kappa S/(3pi^2)] delta+o(delta)`,
`S=8log(2)-13pi^2/24`.                                      (14)

Relative to the loaded clock, the first formula is instead
`Re omega_p/x_*=1+(kappa/pi)delta log(1/delta)+O(delta)`.
For the width, insert the positive spectral density at `omega=x0` in
`Gamma=rho/(4A x0)`; frequency and residue corrections are relatively
`O(delta log(1/delta))`. The local analytic continuation just described
and the simple zero of the unperturbed kernel justify perturbing that
zero. Equivalently apply the implicit-function argument with perturbation
size `delta log(1/delta)`. The imaginary sign places it below the real
axis. Equation (14) concerns a continuum resonance. It does not establish
indefinite exponential decay, a finite-channel decay window, nonlinear
stability, or the existence of the four-dimensional background state.

## Upgrading the actual dipole density

In the notation of the previous note, put `g(u)=u f(u)`,
`a=1-eta^2`, `Q=sqrt(a^2+4u^2 eta^2)`, `r=(a+Q)/2`, and
`F=r+4eta^2 r(r-u^2)/Q^2`. Then `g=int_0^1 sqrt(F)d eta`.
The field-line relation gives

`F-1=(1-r)(4r^2-Q^2)/Q^2>=0`.

Set `zeta=a/Q`, `b0=4eta^2 r^2/Q^2`,
`b1=4eta^2 r u^2/Q^2`. Since `r<=1`, `b0<=4`, `b1<=1`,

`u F_u=r(1-zeta)-2b0 zeta(1-zeta)-b1(1-zeta+2zeta^2)`

is bounded in absolute value by five. Since `F>=1`, the logarithmic
derivative of `sqrt(F)` is bounded by `5/2`, uniformly at the endpoint.
It tends to zero for each interior `eta`. Differentiation under the
integral followed by dominated convergence proves

`u g'(u)->0`, hence `-u^2 f'(u)=g(u)-u g'(u)->kappa`.         (15)

In particular `f` is strictly decreasing sufficiently near zero. For any
compact frequency band inside `(0,4pi/T)`, all roots are in this tail for
small enough `d`; the bound `f(u)<=sqrt(5)/u` excludes roots outside it.
At `D_n=2pi n/omega-T`, the exact Jacobian factor is

`1/[d |f'(u_n)|]=d g(u_n)^2/[D_n^2(g(u_n)-u_n g'(u_n))]`.

The fraction multiplying `d/D_n^2` tends uniformly to `kappa`, since all
`u_n<=sqrt(5)d/min D_n`. The summable `D_n^-2` bound therefore transfers
the complete leading positive-frequency density, pointwise and locally
uniformly, from the inverse-tail model to the actual dipole. At sine
zeros use absolute convergence; away from them the relative error tends
to zero. This closes the former length-derivative gap. It does not on its
own transfer the entire complex resonance theorem or control the physical
Dirac channel approximation.

## Verification and next physical test

Run `python checks/collective_clock_response.py`; see the
[receipt](../evidence/collective-clock-response.json). It checks the
reduction, Lorentzian boundary sign, RN normalization and charge/energy
identities; compares (9) with direct spectral sums with analytic tail
bounds; compares the geometric derivative with independent finite
differences; checks the actual dipole spectral Jacobians; and evaluates
the principal value by explicit pole subtraction and a bounded infinite
tail. A separate coupled clock/two-level evolution checks energy exchange
at two tolerances. That finite toy system tests the constraint and energy
identity, not continuum damping. The logarithmic frequency coefficient
also uses the ordinary analytic argument; a complex-pole numerical solve
has not been substituted for that argument.

Continuation verification reran the complete script successfully and
checked the retarded pole residues separately, including the vanishing
first harmonic. The [shrinking-sector argument](long-return-sector.md)
now controls local packets on bounded rescaled arcs with effective
parameter `q u`. It does not complete a loop or justify its vacuum.

The next decisive physical task is to bound the total matched vertex and
local kinetic terms, construct the finite-time state, and control the
shrinking `u=O(d/ell)` channel reduction. Near-horizon classical matching,
one-Dirac-loop reference corrections, exterior motion and finite-channel
recurrences must share one parameter hierarchy. The new pointwise density
does not discharge these requirements. A refreshed citation review must
also distinguish this candidate from existing dissipative JT/bath models.
No Lean or cloud execution occurred; no new scope or support mechanism is
selected.

Remaining to-do list: control the physical matching, state and projection
errors for this mode; complete novelty and primary/fallback acceptance;
obtain and independently verify the consequential physical result; certify
suitable endpoints on GCP and prepare the reproducible manuscript.
