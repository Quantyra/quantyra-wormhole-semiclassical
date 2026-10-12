# Finite spherical mouth modes in the point-dipole exterior

2026-10-10. S052/S054. This calculation identifies exact intrinsic angular
zero modes on a finite spherical port in the **leading flat, static,
axisymmetric point-dipole model**. Their profiles include the other mouth's
field. In the selected sectors, the angular complement has a quantitative
gap and the radial source has an angular inverse of norm at most
`27/sqrt(160) (rho/d)^2`. This is an ordinary analytic estimate, with
supplementary independent numerical checks. It is not a bound on the full
mouth return or the four-dimensional energy.

The result changes the matching prescription: use the deformed modes, or
pay an explicit profile error. At fixed nonzero `rho/d`, replacing them by
the isolated-monopole modes is not uniformly accurate as flux increases.
The compact-sphere gap does not contradict the gapless full-plane result
in [the transverse audit](dipole-transverse-modes.md).

## Geometry, spin frame and flux normalization

Use the conventions of [the exact exterior reduction](exterior-dirac-window.md).
Around the positive monopole at `z=-a`, put `d=2a` and

```text
r=rho sin(theta), z=-a+rho cos(theta), epsilon=rho/d,
x=sin^2(theta/2), c=cos(theta)=1-2x,
S=sqrt(1+epsilon^2-2 epsilon c).
```

Assume integer `q>=3` and `0<epsilon<=1/10` for the quantitative bounds.
The exact angular-mode formulas themselves hold for integer `q>=1` and
`0<epsilon<1`. The port lies in the intended overlap region only if the
mouth's gravitational radius is also much smaller than `rho`; that
additional approximation has not been bounded here.

Restricting the actual point-dipole potential to this sphere gives

```text
A_phi = q x + delta A,
delta A = q/2 [1-(1-epsilon c)/S]
        = 2q epsilon^2 x(1-x)/[S(S+1-epsilon c)].          (1)
```

The correction vanishes at both poles and adds no net flux. The normal
magnetic field obeys

```text
2 rho^2 B_normal = (d A_phi/dtheta)/[sin(theta)/2]
 = q [1+epsilon^2(c-epsilon)/S^3]
 >= q [1-(epsilon/(1-epsilon))^2] >= 80q/81.             (2)
```

The first lower bound follows independently by bounding the remote
monopole's whole field by `q/[2(d-rho)^2]` before taking its normal
component. Hence `A_phi` increases strictly from zero to `q`.

For half-integer `ell`, the unitary transformation is

```text
Psi_ell = exp(i ell phi) exp(-i sigma_z phi/2)
          exp(-i sigma_y theta/2) eta_ell(rho,theta)
          /[rho sqrt(2 pi sin(theta))],
H_tilde = -i sigma_z partial_rho + L_ell(rho)/rho,
L_ell = -i sigma_x partial_theta + sigma_y s_ell,
s_ell = [ell-A_phi]/sin(theta).                         (3)
```

The norm is `integral |eta|^2 dtheta drho`. The integrated outward radial
flux is `integral eta^dagger sigma_z eta dtheta`. Rotation of the previous
meridional operator produces `-i sigma_z/(2rho)`; multiplication by
`rho^(-1/2)` cancels it. The previous cylindrical measure contributes the
remaining factor in (3). There is no radial gauge component in this gauge.
Usual gauge patches and smooth spinors at the poles specify the domain;
these are not adjustable reflecting boundary conditions.

## Exact zero modes and the nonzero angular gap

For `ell=m+1/2`, `m=0,...,q-1`, there is one normalized upper-spin mode
`(g,0)`, where `g` is real and positive. It solves `g'=s_ell g`. Write

```text
g_0(theta) proportional to
  sin(theta/2)^ell cos(theta/2)^(q-ell),
g_epsilon = g_0 exp(-W)/sqrt(E_0 exp(-2W)),
W(x) = (q/2) log[(S+1-epsilon)/((1-epsilon)(S+1+epsilon))],
W(0)=0, W(1)=-(q/2)log(1-epsilon^2),
W_x = q epsilon^2/[S(S+1-epsilon c)].                   (4)
```

Here `E_0` uses the probability measure

```text
g_0(theta)^2 dtheta = p_0(x) dx,
p_0(x)=x^m (1-x)^(q-m-1)/Beta(m+1,q-m).                 (5)
```

Since `W_theta=delta A/sin(theta)`, (4) solves the full angular equation.
The finite multiplier in (4) leaves the pole powers unchanged. Regularity
of the physical spinor selects exactly these `q` half-integer sectors;
the lower-spin reciprocal solution is singular at a pole. Thus the
angular kernel has dimension `q` and unit outward flux for each normalized
mode. Charge/orientation reversal supplies the other mouth's convention.
These are intrinsic angular kernels at fixed radius, not exact radial
solutions of the full stationary equation.

The following gap applies in each of these selected sectors and hence
their direct sum. The diagonal squares are

```text
(L^2)_upper=-partial_theta^2+s^2+s',
(L^2)_lower=-partial_theta^2+s^2-s'.
```

Set `t=s+cot(theta)`. Direct factorization gives

```text
(L^2)_lower = (partial_theta+t)(-partial_theta+t)
                         +1+2 A_phi'/sin(theta).       (6)
```

The regular lower component has pole powers at least `3/2` at each end in
the selected sectors. Integration by parts has no boundary term, and
the first term in (6) is nonnegative. Smooth sections are a core for the
compact angular Dirac operator; the bound extends to its square-form
domain. Its nonzero eigenvalues pair between the two chiralities, while
the upper zero mode is removed. Equations (2) and (6) give

```text
gap(L/rho on kernel complement)
 >= sqrt(1+q[1-(epsilon/(1-epsilon))^2])/rho
 >= sqrt(1+80q/81)/rho.                                (7)
```

At zero deformation this recovers the known first gap `sqrt(q+1)/rho`.
No gap of the noncompact exterior operator follows from (7).

The maximum of `g` occurs at `A_phi=ell`, exactly the mass-zero contour of
the meridional operator. At that maximum,
`partial_theta^2 log(g)=-rho^2 B_normal`. Its leading Gaussian transverse
width in physical length is `1/sqrt(B_normal)`. This identifies the local
Gaussian profile used in a possible field-aligned construction.

## Uniform estimates for profiles and the radial source

Put `w=W/q`, `D=S(S+1-epsilon c)`. On the stated interval,

```text
D>=81/50, |w_x|<=(50/81)epsilon^2,
w_xx<0, |w_xx|<=4epsilon^3,
|partial_epsilon w_x|<=(3/2)epsilon.                   (8)
```

For explicit constants, `S` lies in `[9/10,11/10]` and
`D_x=epsilon(3S^2+4S+1-epsilon^2)/S`. Thus
`|w_xx|/epsilon^3<=225750/59049<4`.
At fixed `x`, `|S_epsilon|<=11/9` and
`|(S+1-epsilon c)_epsilon|<=20/9`, so `|D_epsilon|<=77/15`.
Differentiating `w_x=epsilon^2/D` then gives the last bound since
`2/(81/50)+(77/150)/(81/50)^2<3/2`.

Interpolate `g_t=g_0 exp(-tW)/sqrt(E_0 exp(-2tW))`, `0<=t<=1`.
The negative logarithm `V_t` of its `x` density satisfies

```text
V_t''=m/x^2+(q-m-1)/(1-x)^2+2t W_xx
 >=q-1-8q epsilon^3 >=q/2.                            (9)
```

Consequently `Var_t(f)<=(2/q)E_t(f_x^2)`. One elementary proof first
works on a truncated interval with normalized density proportional to
`exp(-V_t)` and solves `h''-V_t'h'=f-E_t f` with Neumann endpoints.
Integration by parts gives
`E[(h''-V_t'h')^2]=E[h''^2]+E[V_t'' h'^2]` and
`Var(f)=-E[f'h']`. Cauchy-Schwarz proves the bound. Letting the truncated
endpoints tend to zero and one proves it for the smooth bounded functions
used below; their expectations converge by dominated convergence.

Since `partial_t g_t=-(W-E_t W)g_t`, (8)-(9) and integration in `t` yield

```text
||g_epsilon-g_0|| <= (50 sqrt(2)/81) sqrt(q) epsilon^2. (10)
```

Similarly `partial_rho g=-(W_rho-E W_rho)g`, with `d` fixed, so

```text
||partial_rho g|| <= (3/sqrt(2)) sqrt(q) rho/d^2.       (11)
```

The derivative is orthogonal to the angular kernel: normalization gives
orthogonality in its own sector and azimuthal labels give it in the other
sectors. This step uses the unitary measure in (3), rather than
differentiating in a radius-dependent physical measure.

The factor `sqrt(q)` in (10) cannot simply be dropped for isolated-mode
matching at a fixed nonzero `epsilon`. For any sequence `m/q -> p` with
`0<p<1`, the normalized undeformed density concentrates at `x=p`, whereas
the deformed density concentrates at the unique minimizer of

```text
-p log(x)-(1-p)log(1-x)+2w(x).
```

Its strict convexity follows from (9) divided by `q`; its derivative at
`p` is `2w_x(p)>0`, so its minimizer is strictly below `p`. Laplace
concentration on disjoint neighborhoods of the two minima and
Cauchy-Schwarz on each neighborhood imply

```text
<g_epsilon,g_0> -> 0,
||g_epsilon-g_0|| -> sqrt(2), as q -> infinity.         (12)
```

This is a failure of the identity approximation between two chosen
bases. It does not disprove adiabatic transport in the deformed basis or
the MMP wormhole. Condition `sqrt(q)epsilon^2 << 1` is sufficient for
(10) to justify the isolated basis; we do not assert it is a sharp
necessary condition in every angular sector.

## Source-specific inverse and the unresolved propagation residual

Let `eta=(0,h)` be the unique kernel-orthogonal solution of

```text
(L/rho)eta = i(partial_rho g,0).
```

Equations (7) and (11) imply the flux-independent bound

```text
||eta|| <= (27/sqrt(160))epsilon^2.                    (13)
```

At a fixed port, the corrected vector `(g,0)+eta` has outward flux
`1-||eta||^2 >= 1-(729/160)epsilon^4 > 0`. Dividing by the square root
of that flux gives unit flux. This pointwise normalization does not
control the radial derivative of the normalizing factor.

An explicit formula, useful for checking regularity, is

```text
h(theta)=-rho/g(theta) integral_0^theta g(t)partial_rho g(t)dt
        =-rho/[2g(theta)] partial_rho F_rho(theta),
F_rho(theta)=integral_0^theta g(t)^2 dt.                (14)
```

The numerator vanishes at both endpoints, since `g` is normalized. Its
local powers make `h` a regular lower component. Differentiating (14)
gives `h'+s h=-rho partial_rho g`, the required angular inverse. For fixed
`q,m`, expansion at small `epsilon` gives

```text
h/g= -[q/(q+1)]epsilon^2 sqrt(x(1-x))+O(epsilon^3).     (15)
```

For example, (15) follows from `w=epsilon^2 x/2+O(epsilon^3)` and
`integral_0^x (t-(m+1)/(q+1))p_0(t)dt=-x(1-x)p_0(x)/(q+1)`.
The leading sign and angular shape agree with the spin tilt into the
distant monopole's field. No uniform-in-`q` remainder in (15) is asserted.

There is a direct one-dimensional target for the missing derivative
estimate. In the original angular measure, write the normalized density as
`r_epsilon(theta)=g_epsilon(theta)^2`, its CDF as
`F_epsilon(theta)=integral_0^theta r_epsilon(t)dt`, and set
`u=partial_epsilon log r_epsilon` and
`v=partial_epsilon^2 log r_epsilon`. The exact corrector formula (14) then
gives

```text
h = -epsilon F_epsilon/(2 sqrt(r_epsilon)),
partial_rho h = (1/d) partial_epsilon h,
partial_epsilon h = -[partial_epsilon F_epsilon
  + epsilon partial_epsilon^2 F_epsilon
  - (epsilon/2) u partial_epsilon F_epsilon]/(2 sqrt(r_epsilon)),
partial_epsilon F_epsilon = integral_0^theta r_epsilon u,
partial_epsilon^2 F_epsilon = integral_0^theta r_epsilon (u^2+v).
```

Thus the norm needed in (16) is an explicit weighted cumulative-score
integral in `dtheta`, with weight `1/r_epsilon`; a bound on the angular
inverse norm does not bound it. The existing convexity estimate controls
variances of bounded observables, but `partial_epsilon h` has the
endpoint-sensitive factor `1/sqrt(r_epsilon)`. The next subsection rules out
a flux-independent derivative estimate at fixed nonzero deformation; the
physically relevant joint-scale estimate remains open. This reduction does
not supply the throat/port scattering estimate.

### Large-flux obstruction to a flux-independent derivative bound

The weighted integral has a concrete large-`q` scaling obstruction. Fix
`m/q -> p` with `0<p<1` and a small fixed `epsilon>0`. Laplace's method
places `p_epsilon` near the unique minimizer `x_epsilon` of

```text
I_epsilon(x)=-p log(x)-(1-p)log(1-x)+2w(x,epsilon),
w=W/q,
```

with curvature `kappa=partial_x^2 I_epsilon(x_epsilon)>0`. In the central
window `x=x_epsilon+z/sqrt(q kappa)`, the normalized density is Gaussian to
leading order with respect to `dx`. Set `J(x)=sqrt(x(1-x))`; since
`dx/dtheta=J(x)`, the angular density in (14) is
`r_epsilon(theta)=J(x) p_epsilon(x)`. Differentiating the probability CDF
at fixed `x` gives
`partial_epsilon F_epsilon/p_epsilon -> -partial_epsilon x_epsilon`
uniformly for bounded `z`. The differentiated Laplace expansion applies on
these compact central windows because the minimizer stays interior and
`kappa` is positive. Hence there

```text
h = (epsilon/2J(x_epsilon))(partial_epsilon x_epsilon) sqrt(r_epsilon)
    + O(q^(-1/2) sqrt(r_epsilon)),
partial_epsilon h = (epsilon/4)(partial_epsilon x_epsilon)^2
    sqrt(q kappa) z sqrt(r_epsilon)/J(x_epsilon)
    + O(sqrt(r_epsilon))
```

on each fixed central window. The leading derivative has nonzero Gaussian
`L2` mass there, so its norm grows at least as a positive constant times
`sqrt(q)` for fixed nonzero `epsilon` whenever
`partial_epsilon x_epsilon != 0`. The latter condition holds for all
sufficiently small fixed positive `epsilon`: since
`w_x=epsilon^2/2+O(epsilon^3)`, stationarity gives
`x_epsilon=p-epsilon^2 p(1-p)+O(epsilon^3)` and
`partial_epsilon x_epsilon=-2 epsilon p(1-p)+O(epsilon^2)`.
For example, the leading coefficient in the central-window derivative
norm at small `epsilon` scales as
`epsilon^3 p(1-p) sqrt(q)`.

This rules out a flux-independent bound on `||partial_epsilon h||` at fixed
deformation in the large-flux limit; it does not rule out a useful bound
under a joint scaling such as `sqrt(q)/d` bounded or tending to zero. In
(16), the derivative contribution therefore scales at least as
`sqrt(q)/d` in this regime and must be compared with `xi ||h||` using the
physical scaling of `d`, `xi`, and `q`. This is a local Laplace asymptotic,
not yet a uniform error theorem or a counterexample to the complete MMP
scattering problem; the joint physical parameter regime remains to be
checked.

For imaginary frequency `i xi`, the trial radial state
`exp[-xi(rho-rho_0)] [(g,0)+eta]` cancels the leading radial source in
`(H_tilde-i xi)`. Its remaining lower-component residual is exactly

```text
i exp[-xi(rho-rho_0)] [partial_rho eta-2xi eta].        (16)
```

A small corrector norm does not bound its derivative or the global
boundary-value solution. Equation (16), transfer from the throat to the
port, propagation beyond the spherical neighborhood, escaping continuum,
and flux-normalized extraction at the opposite port still need control.
In particular, (13) is not the previously targeted
`0.001 exp(-xi D)` return-error estimate.

## Prior work and verification

[Maldacena, Milekhin and Popov, Appendix C](https://arxiv.org/html/1807.04726v3#A3)
already derive the spherical Dirac zero modes and reduced massless
channels. Their charge/gauge and chirality conventions differ; the mode
count and isolated profiles agree after changing the patch and spin frame.
[Dolan and Hunter-McCabe, sections II.2, III and IV](https://arxiv.org/html/2001.02208v3)
give the monopole zero modes, pole regularity and higher spectrum,
including the isolated gap used here. Those sections were read; their
many-body quantum Hall claims are not used.
[Hasebe, section 6.3](https://arxiv.org/html/2410.16193v2#S6.SS3)
discusses zero modes under ellipsoidal deformation. This is an adjacent
precedent, not a proof of the point-dipole matching estimates.

[Bal et al., Theorem 2 and sections 3.1-3.6](https://arxiv.org/html/2106.00729v1#S3)
construct propagation for prepared Schwartz packets with Gaussian
transverse profile. That proof chain has now been read in the preprint;
the publisher/dependency comparison remains open. The constants depend
on preparation seminorms. Their Fourier-Hermite complement block is
`[[0,sqrt(2n+2)],[sqrt(2n+2),2k]]`: its inverse grows with `|k|`, and its
small eigenvalue tends to zero. Continuity on Schwartz space therefore
does not imply a uniform `L2` inverse or an arbitrary-input port estimate.
The recorded check reconstructs that inverse directly; the preprint's
displayed inverse has the opposite sign. This does not affect the
existence or seminorm estimates. The proposed all-input dispersion
statement is presented there as a conjecture, not the theorem imported
by this program. No global dipole scattering estimate is inferred.

The zero-mode mechanism, compact angular gap and prepared-packet method
are established mathematics. The explicit dipole deformation and constants
here serve the physical matching audit. No standalone novelty acceptance
or new physical obstruction follows from this milestone.

Run

```powershell
python checks/mouth_mode_matching.py --output evidence/mouth-mode-matching.json
```

The [receipt](../evidence/mouth-mode-matching.json) pins the script and
library versions. It checks exact identities and rational bounds,
independently normalized profiles and radial derivatives, finite-difference
angular spectra at three resolutions, and the full Cartesian Weyl
operator after the spherical transform. Large-flux overlap calculations
supplement (12). Numerical quadrature and finite differences are not
interval certificates or a substitute for the ordinary proofs above.
No Lean invocation occurred; any later Lean checks must run on GCP.

Remaining to-do list: control (16) and actual throat/port matching, prove
the exterior return error on the selected frequency window, bound total
energy and common-observer corrections, complete the decisive source and
novelty comparison, then certify suitable endpoints and prepare the manuscript.
