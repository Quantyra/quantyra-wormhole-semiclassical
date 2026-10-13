# Rotation with the AdS magnetic boundary field fixed

S072, 2026-10-13. **The Kerr-Newman moment of inertia used in the previous
compact-mode comparison is not the fixed-boundary-field susceptibility.**
Magnetic Kerr-Newman-AdS changes the electric field on a nonrotating AdS
boundary as its angular momentum changes. Keeping that source fixed gives
a different coupled metric/Maxwell problem. Its static response is bounded
below and above here, including the projected Gaussian LLL chemical
polarization. This supplies a single-species candidate for the next
frequency-dependent test. It does not yet supply its compact quantum
kernel, a four-dimensional opening solution or a stability theorem.

The previous [compact-phase comparison](compact-phase-transfer.md) remains
a conditional calculation with the stated Kerr-Newman inertia. Its
single-species conclusion must not be transferred to the boundary ensemble
specified here. The earlier many-species family is retained as a historical
candidate, but is no longer necessary for the present static scale test.

## Boundary condition before susceptibility

Use `-F^2/(4g^2)`, magnetic potential `A_phi=-P cos(theta)`, `P=q/2>0`,
and `Q^2=4 pi G P^2/g^2`. Patch constants in the monopole potential are
understood; the expressions below use its central patch. Write the RN-AdS
exterior as

```
f(r)=1-2m/r+Q^2/r^2+r^2/L^2,
R=r_+, X=R^2/L^2, Y=Q^2/R^2,
tau=1+3X-Y=4 pi R T >= 0,  Y>0.
```

The ensemble fixes the static round boundary metric, magnetic flux and
boundary electric dipole source to zero. It allows angular momentum to
fluctuate; `Omega` is its infinitesimal chemical potential. In nonrotating
coordinates the stationary axisymmetric perturbation is

```
ds^2=-f dt^2+dr^2/f+r^2[dtheta^2+sin^2(theta)(dphi-omega(r)dt)^2],
A=-P cos(theta)dphi+V(r)cos(theta)dt,
omega(R)=Omega, V(R)=P Omega, omega(infinity)=V(infinity)=0.       (1)
```

Only terms linear in the fields in the equations, and quadratic in the
response functional, are retained. The square in the metric specifies a
convenient off-shell extension, not a nonlinear solution. The horizon
condition makes `A_t+Omega A_phi` angle independent and regular in the
corotating patch. In corotating coordinates it is the familiar boundary
rotation source, together with the transformed monopole potential.

For comparison, the magnetic Kerr-Newman perturbation with parameter `a`
in these **nonrotating** coordinates is

```
omega_KN=a(2m/r^3-Q^2/r^4),
V_KN=a P(1/r^2+1/L^2),
Omega_KN=a(1+X)/R^2,
V_KN(infinity)/(P Omega_KN)=X/(1+X).                            (2)
```

Thus its boundary `F_(theta t)=-aP sin(theta)/L^2` is nonzero. It cannot
be removed by a gauge transformation preserving the field strength.
[Caldarelli, Dias and Klemm, section 3.2, equations 3.69--3.72](https://arxiv.org/html/0812.0801v2)
explicitly exhibit the boundary coordinate transformation and the induced
electric field. [Chen and Kim, section 2](https://link.springer.com/article/10.1140/epjc/s10052-023-11380-7)
provide an independent dyonic metric/potential convention; their magnetic
sign is reversed relative to ours. Neither reference supplies the new
fixed-source response bound below. No novelty acceptance is claimed.

At extremality `R=b`, `Y=1+3X`, the previously used coefficient is
`I_KN=b^3(1+2X)/[G(1+X)]`. This is a correct Kerr-Newman family coefficient.
It is not the coefficient of (1) at `X>0`. The distinction concerns a
physical boundary field, not just a choice of angular coordinates.

## Coupled positive static functional

The Einstein-Maxwell quadratic functional in this axial sector is

```
F_0[omega,V] = (1/2) integral_R^infinity [
  r^4 omega'^2/(6G) + 4 pi r^2 V'^2/(3g^2)
  + 8 pi (V-P omega)^2/(3g^2 f) ] dr.                         (3)
```

Here `F_0` denotes the positive response functional, equal on shell to
`I Omega^2/2`; the grand potential changes by its negative. It is not the
unconstrained Euclidean action for every gravitational perturbation.
For a direct derivation, the stationary ADM extrinsic curvature obeys
`K_ij K^ij=r^2 sin^2(theta) omega'^2/2`. Integrating its Einstein term over
the sphere gives `r^4 omega'^2/(12G)`. The radial electric field is
`V' cos(theta)`, while its transverse component in the normal frame is
`-(V-P omega)sin(theta)`. Their spherical integrals are respectively
`4pi/3` and `8pi/3`, giving the other two terms. The RN background
constraints remove variations of second-order lapse/radius fields at
this quadratic, fixed-source order.

Define `p_omega=r^4 omega'/(6G)`, `p_V=4pi r^2 V'/(3g^2)`. Variation gives

```
p_omega'=-8 pi P(V-P omega)/(3g^2 f),
p_V'=8 pi(V-P omega)/(3g^2 f),
(p_omega+P p_V)'=0.                                         (4)
```

The Kerr-Newman profiles (2) satisfy both equations exactly, with their
own nonzero electric boundary value. Their quadratic functional reproduces
`I_KN` at extremality. This is a check of normalization and mixing, not a
license to use their boundary data in (1).

The conserved angular charge is
`J=-[p_omega+P p_V]_infinity`. The second term is the electromagnetic
contribution associated with the lifted rotation of the monopole bundle;
retaining just the metric flux gives the wrong conjugate response in
general. Varying the horizon data in (1), while keeping the infinity
sources fixed, gives `dF_0,min/dOmega=J`. Homogeneity then gives
`J=I_0,fix Omega`, where `F_0,min=I_0,fix Omega^2/2`.

All coefficients in (3) are positive outside the horizon. The affine
space of finite-energy fields with (1) has a nonempty trial space and
a strictly positive derivative norm on variations with zero boundary
values. Completing it in the quadratic norm gives a unique weak
minimizer by the Hilbert-space projection theorem. This proves the
linear static boundary problem; it does not prove a nonlinear black-hole
family. At positive temperature the regular horizon branch is used.
At zero temperature only the finite-energy variational limit is needed
below; no new smoothness theorem for a matter-dressed extremal horizon
is assumed.

## Fermions see a combined rotation and electric dipole

For each of `N_f` massless four-dimensional Dirac species, the LLL has
`j=(q-1)/2`, `m=-j,...,j`. With our sign of magnetic flux, its angular
density is proportional to

```
(1-c)^(j+m) (1+c)^(j-m),  c=cos(theta),
<c>_m=-m/(j+1)=-2m/(q+1),
S_J=N_f sum_m m^2=N_f q(q^2-1)/12.                          (5)
```

This follows by a beta integral from the monopole zero modes in
[MMP appendix C](https://arxiv.org/html/1807.04726v3), with flux sign
reversed. In the coframe of (1), the temporal covariant derivative
contains `partial_t+omega partial_phi-i(V-P omega)cos(theta)`.
The extra connection term proportional to `omega'` contains one sphere
gamma matrix and has zero matrix element within the chiral LLL. It can
mix higher levels, which are omitted from the projected determinant.
The projected one-particle shift is consequently

```
m omega-(V-P omega)<c>_m = m [omega+2V]/(q+1).
```

In the ensemble `H-Omega J`, the LLL chemical field is `m B(r)`, where

```
B=Omega-[omega+2V]/(q+1).
```

It vanishes at the horizon and approaches `Omega` at the boundary. The
static massless-Dirac chemical polarization gives the positive response
term

```
F_LLL=(S_J/(2pi)) integral_R^infinity B^2/f dr.              (6)
```

The local coefficient is the Gaussian bosonized response used in the
previous spherical calculation, now with the correctly projected
angular generator. This is the full minimization of (3)+(6), not just
one insertion on a chosen metric profile. It includes both axial
Einstein-Maxwell fields in this projected static sector. It excludes
higher-Landau mixing, the remaining photon/graviton determinant, thermal
and vacuum background-stress changes and the real-time influence kernel.

The charge normalization continues to hold with (6). The field equations
give `(p_omega+P p_V)'=-S_J B/(pi f)`. Direct differentiation of the
minimum includes the explicit `Omega` derivative of (6) and yields

```
dF_min/dOmega=-[p_omega+P p_V]_R + integral S_J B/(pi f) dr
            =-[p_omega+P p_V]_infinity = I_fix Omega.        (7)
```

Spherical electric and angular chemical fields have no mixed quadratic
term, since `sum_m m=0`. Rotational symmetry makes the three infinitesimal
angular coefficients equal. Neither fact determines the global compact
charge-spin lattice or the fermion transfer matrix elements.

## A uniform static interval, including positive temperature

Put `u=R/r`, `H_X(u)=X+2Xu+(1+3X)u^2`. The exact RN identity is

```
f=(1-u)^2 H_X/u^2 + tau u(1-u) >= (1-u)^2 H_X/u^2.           (8)
```

Define dimensionless `d=G I_fix/R^3` and

```
d_lower=1/2+Y/3,
gamma=G S_J/(pi R^2)
     =N_f g^2 (q^2-1)Y/(12 pi^2 q),
D(X)=integral_0^1 du/H_X
    =atan(sqrt(1+2X)/(2sqrt(X)))/sqrt(X(1+2X)),
B_X=integral_0^1 u^2(1+u)^2/H_X du,
P_Xq=integral_0^1 [1+(u+u^2)/(q+1)]^2/H_X du.
```

Then, for `X>0`, `tau>=0`, `Y>0`, in the specified static model,

```
d_lower <= d <= d_lower + (2Y/3) B_X + gamma P_Xq
             <= d_lower + 2Y/(9X)
                + gamma [1+2/(q+1)]^2 D(X).                (9)
```

Proof: drop the nonnegative mixing and matter terms. Cauchy-Schwarz
with the endpoint differences in (1) bounds the two derivative terms
by `R^3/(2G)` and `4pi R P^2/(3g^2)=Y R^3/(3G)`, respectively, in units
of the inertia. This proves the lower bound. For the upper bound use
`omega=Omega u^3`, `V=P Omega u`. These fields separately attain the
two derivative lower bounds. Their residual mixing and chemical fields
are `P Omega u(1-u)(1+u)` and
`Omega(1-u)[1+(u+u^2)/(q+1)]`. Apply (8) to their remaining integrals.
Finally `H_X>=X(1+u)^2` gives `B_X<=1/(3X)`, and the bracket in `P_Xq`
is at most `1+2/(q+1)`. This completes the ordinary interval proof.

At extremality, `d_lower=X+5/6`. Thus

```
I_fix = (b^3/G)[X+O(1)+O(N_f g^2 q)]                 (X large),
```

with explicit one-sided errors in (9). The purely classical coefficient
is already of order `X b^3/G`; this enhancement is not caused by an
assumed large fermion correction to `I_KN`. There is no uniform flat-space
matter bound from this argument because `D(X)` diverges as `X->0`.
The finite-temperature inequality avoids promoting a possibly nonuniform
extremal matter profile to a regular physical horizon.

## Revised candidate and the remaining transfer gate

If this susceptibility is the inertia of the relevant low-frequency
compact mode, its thermal exponent obeys, on the extremal reference
background and with the previous initial temperature,

```
chi_J=beta j(j+1)/(2 I_fix)
 <= pi(q^2-1)/[8 N_f q eta_f z_0 (1+6X)(X+5/6)].             (10)
```

It has an additional factor `1/X` relative to the earlier Kerr-Newman
comparison. That *conditional* compact-mode use still requires frequency
and state matching. It cannot be inferred merely from static positivity.

An explicit next family, with just one four-dimensional Dirac species, is

```
k integer >=2000,
eta_f=k^-1, g=k^-4, L/l_p=k^5, q=k^4, N_f=1, z_0=1/2,
X(1+3X)=pi k^6, X~sqrt(pi/3) k^3.                          (11)
```

The electric thermal exponent tends to `1/(24k)`, the angular upper bound
in (10) tends to `1/(8k)`, and `(q+1)/(1+6X)` grows as
`k/[6sqrt(pi/3)]`. Spherical screening and the angular LLL trial
increment divided by the inertia lower bound are `O(k^-7)`; the total
relative static interval width is `O(k^-3)`. Moreover
`NG/b^2=O(k^-9)`, `N eta_f/Phi_ext=O(k^-10)`, `Tb=O(k^-7)`, and
`N eta_f=2k^3` grows. Here `N=2q`, `Phi_ext=b^2/(2G)` and Newton's
constant and couplings are the chosen renormalized reference parameters.

At `k=2000`, the necessary-scale example has

| Quantity | Value or bound |
| --- | --- |
| Four-dimensional species | 1 |
| Magnetic flux integer | `1.6e13` |
| `X` | `8.1866136634e9` |
| `G I_fix/b^3` interval | `[8186613664.2386, 8186613664.7306]` |
| Electric thermal exponent upper bound | `2.08334e-5` |
| Conditional angular thermal exponent upper bound | `6.25001e-5` |
| First higher-LLL mass squared in AdS2 units | `325.7350` |
| Angular LLL trial increment / inertia lower bound | `<8.417e-26` |

The last row is the matter contribution evaluated on the trial fields,
not a separate bound on `I_fix-I_0,fix`. That difference is bounded by the
entire interval width in (9); its relative upper bound in this example
is `8.144e-11`.

These are extremal reference scales. The prepared pair is slightly
nonextremal; (9) applies at its actual `R,T,Y`, with no need to set `T=0`.
At fixed magnetic flux the radius displacement is
`(R-b)/b=O(Tb/(1+6X))=O(k^-10)`, so it does not change these necessary
asymptotic comparisons. Background quantum and thermal stress errors have
not been bounded by this observation. The enormous magnetic flux and weak
coupling are theoretical choices, not a material or laboratory proposal.

The next decisive calculation is the **frequency-dependent coupled
metric/Maxwell/fermion response with exactly (1)'s boundary source**,
including the joint compact charge-spin state. In particular, a static
inertia is not automatically the coefficient appropriate to a switch on
the initial thermal time scale. The other light angular sectors at large
`X` remain present. Errors must be bounded against the previous capture
and total-area margins, rather than just shown to be small in isolation.
Initial-pair preparation, nonspherical stability and finite-energy signal
backreaction remain separate unresolved physical requirements.

## Verification and scope

[Checker](../checks/angular_boundary_response.py),
[receipt](../evidence/angular-boundary-response.json), and
[source review](../evidence/angular-boundary-source-review.json).
Symbolic variation reproduces the two coupled Einstein-Maxwell equations
on magnetic Kerr-Newman and its action coefficient under its own boundary
conditions. Twelve independent angular-density integrals verify the LLL
projection. Nested Ritz minimization checks 24 fixed-source problems,
including positive-temperature cases and substantial projected matter
polarization; four Kerr-Newman comparisons recover the analytic answer.
The rigorously sufficient bounds are the ordinary inequalities (9), not
the numerical minimizers or their refinement differences. Three exact-flux
scaling examples verify the displayed necessary comparisons.

The previous turn made progress by exposing missing compact phases. This
turn makes progress by matching the axial static fields and identifying
the boundary-ensemble mismatch that changes the selected candidate. No
full physical compact-mode matching, novelty acceptance, publication,
outreach, paid computation, Lean certificate or stable creation is claimed.

Remaining to-do list: derive the retarded joint response and switching
kernel at fixed boundary field; bound remaining four-dimensional errors;
establish initial-pair preparation, sufficient stability and finite-energy
transport.
