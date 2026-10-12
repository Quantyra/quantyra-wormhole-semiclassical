# Stationary electric dipole, rotation charge and two-mouth matching

S053, 2026-10-11. The [binary drive](binary-mouth-drive.md) includes an
electric field along the orbital normal. This note solves its classical
electric/axial RN response, retaining the rotation integration constant,
and solves a specified reflection-symmetric RN/product-throat junction.
It does not select the physical quantum state or construct a smooth
four-dimensional binary wormhole. These classical ingredients are prior
work; no consequential novelty is accepted.

## Four-dimensional equations and the rotation datum

First use dimensionless RN units `P=1`, with `r` the areal radius in units
of the positive magnetic charge magnitude and `T=t/P`. At the positive
mouth, take

```text
f=(r-1)^2/r^2,                         c=cos(theta),
ds^2=-f dT^2+dr^2/f+r^2[dtheta^2+sin^2(theta)(dphi-omega dT)^2],
A=-c dphi+V(r)c dT,                    Phi=V-omega.                (1)
```

Keep only first order in `V,omega`. A constant rotation changes both by
the same constant, leaving `Phi` unchanged. Direct linearization of all
Einstein-Maxwell components gives precisely

```text
(r^2 V')'=2Phi/f,
(r^4 omega')'=-4Phi/f,
r^4 omega'+2r^2 V'=-6 J.                                         (2)
```

Here `J` is a dimensionless conserved rotation parameter. Its
normalization is fixed by the asymptotically flat slow Kerr-Newman mode:
`J=J_KN/P^2`. For a uniform external field there is no asymptotically
flat end at infinity; this calibration alone does not define an ADM spin
there. We retain the same rotational lift and electromagnetic gauge
convention when comparing mouths. The first integral, including its
Maxwell contribution, must not be replaced by a metric-only spin.

Put `u=r^2 V'`, the radial electric flux coefficient. Then

```text
(f u')'-(2/r^2+4/r^4)u=12J/r^4,
Phi=f u'/2,      V'=u/r^2,      omega'=(-6J-2u)/r^4.              (3)
```

This reconstruction retains the stationary rotation mode, which can be
lost by dividing a dynamical reconstruction by frequency. Once (3)
holds, one common additive rotation constant completes `V,omega`.

The regular uniform-field family, with dimensionless electric drive
`e=P E_external`, is

```text
V=e(r-2+3/r)+(J-e)/r^2+beta_RN,
omega=2e/r+(J-e)(2r-1)/r^4+beta_RN,
u=e u0-2J/r,                u0=(r-1)^2(r+2)/r.                   (4)
```

The terms with `J-e` are the slow Kerr-Newman solution. In particular,
the uniform-field Ernst branch has `J=e`, even though its seed rotation
parameter vanishes. Zero *conserved rotation charge*, `J=0`, requires a
compensating Kerr-Newman term. The nonspinning point-body orbit by itself
does not choose between these notions of spin.

Restoring general mass `M` and geometric magnetic charge `P`, a direct
check gives the established regular solution

```text
V=E(r-2M+3P^2/r)+P(J_phys-E P^3)/(M r^2)+P omega0,
omega=2PE/r+(J_phys-E P^3)(2Mr-P^2)/(M r^4)+omega0.                (5)
```

In (5) `r,t` have length units and `J_phys` has length squared. The
charge is `J_phys=E P^3` on the uniform-field seed branch. Its magnitude
and sign agree, after electric-magnetic duality and a constant frame
shift, with the linear term of [Brenner, Giribet and Montecchio, section 4
and (30)](https://arxiv.org/pdf/2103.06983v2). Their discussion also makes
the charge's gauge convention explicit. We do not import their nonlinear
asymptotics, mass prescription or thermodynamics into the binary model.

The homogeneous potential in (3) is the static version of the potential
displayed in [Burko (1999), (26)](https://arxiv.org/pdf/gr-qc/9902002v1).
The printed physical/constraint assignment in the paragraph following
(27) is not used: in the neutral limit, (26) has the Maxwell potential
`2f/r^2`, whereas (27) has the nonradiative gravitational dipole potential
`f(2/r^2-6M/r^3)`. The source's own (19)-(21) also distinguish these
limits. Our stationary reconstruction is instead checked against every
four-dimensional equation. The prior stationary coupled RN multipoles
and uniform-field solution are also explicitly identified in
[Bicak and Dvorak (1980)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.22.2933).
Their full article was not retrieved; its proof remains a review gap.

## Regularity and the additional finite-annulus solution

For `J=0`, (4) has

```text
u=3e z^2+O(e z^3),       Phi=3e z^3+O(e z^4),
V=beta_RN+e+e z^3+O(e z^4),
omega=beta_RN+e-2e z^3+O(e z^4),       z=r-1.                    (6)
```

Thus the common constant `beta_RN+e` is a rotating frame, not a local
electric field at the end of the throat. For general `J`, the horizon
limit is `omega_H=beta_RN+e+J` and `u_H=-2J`. A free spin mode has a
different near-throat order and cannot be absorbed into the error of (6).

On a finite RN annulus we must retain the second homogeneous solution,

```text
u1=(2r-1)(5r^2-5r+2)/[30r(r-1)^3],
u=e u0+b u1-2J/r,
f(u0 u1'-u0' u1)=-1.                                           (7)
```

It follows either by substitution or from
`u1=u0 integral_r^infinity [1/(f u0^2)] ds`. It is positive and strictly
decreasing for `r>1`, with limits `u1~1/(15z^3)` near the horizon and
`u1~1/(3r)` at large radius. The singular horizon behavior is no reason
to discard it on a truncated annulus. Its coefficient is fixed by the
inner load. The coefficient `e` in (7) specifies the growing outer datum;
at any finite outer radius it is an explicit Robin functional of `u,u'`
using the Wronskian. It does not impose an actual uniform field to
infinite distance, where this linear approximation would fail.

## Product throat and a conditional reflection-symmetric junction

In the product region use
`ds^2=P^2[-(1+rho^2)d tau^2+d rho^2/(1+rho^2)+dOmega^2]` and the
same infinitesimal rotational/electric ansatz. Here `tau=t/ell` and
`lambda=P/ell`; potentials and angular velocities in the common time
coordinate are `V_th=V_RN/lambda`, `omega_th=omega_RN/lambda`.
Direct four-dimensional calculation now gives

```text
V''=2Phi/(1+rho^2),          omega''=-4Phi/(1+rho^2),
omega'+2V'=C,               Phi=V-omega,
(1+rho^2)Phi''-6Phi=0.                                        (8)
```

All derivatives in (8) are with respect to `rho`. The full stationary
solution can be reconstructed as
`V=Phi/3+C rho/3+beta`, `omega=-2Phi/3+C rho/3+beta`.

Choose antipodal angular gluing and the leading binary electric fields
`E_+=-E_-`, with `P_+=-P_-`. In common throat angular coordinates their
electric scalar potentials and axial rotation coefficients have equal
outer data. If one additionally imposes reflection symmetry through the
core and no extra angular-current source, `V,omega` are even. Their
constant `C` must then vanish. Matching the conserved rotation flux
with opposite outward radial orientations gives `J_++J_-=0`; reflection
gives equal charges, hence both vanish. This is a conditional sector
selection. A quantum angular source or different gluing/state can change
it; it is not a derivation of the actual binary's spin preparation.

In this sector define

```text
H(rho)=rho(1+rho^2)atan(rho)+rho^2+2/3,
g(rho)=H'=(1+3rho^2)atan(rho)+3rho,
V=k H+beta,      omega=-2k H+beta,      u_th=V'=k g,
Phi=3k H.                                                       (9)
```

The flux `u_th` is odd, as required by equal electric data measured in
the two outward radial coordinates. `H` is even. Its other independent
solution is odd, `rho(1+rho^2)`, and is excluded by this reflection choice.

For a precise auxiliary matching problem take `|rho|<=R`, two RN
annuli `s<=r-1<=L`, and `s=lambda R`. Match electric flux and `Phi`:

```text
u_RN(1+s)=u_th(R),
f_RN u_RN'(1+s)=lambda(1+R^2)u_th'(R),                          (10)
```

with outward signs reversed at the lower end. These match the radial
Maxwell flux and invariant electric potential in the common time.
They do not prove matching of the entire induced metric, its extrinsic
curvature or the semiclassical source across a smooth junction.

Equations (7), (9) and (10), at `J=0`, give the exact solution

```text
D=lambda(1+R^2)g'(R)u1(1+s)-g(R)f(1+s)u1'(1+s),
k=e/D,                  b=[k g(R)-e u0(1+s)]/u1(1+s).           (11)
```

For every `lambda,R>0`, `D>0`: `g,g',u1` are positive and `u1'` is
negative. Thus this specified static boundary problem has a unique
solution, with no static resonance in this sector. As `s->0`, `R->infinity`,

```text
k=(2e lambda^2/pi)[1+O(s+R^-2)],
b=O(e s^5(s+R^-2)).                                             (12)
```

To verify the order, use `u0=3s^2[1+O(s)]`,
`u1=(15s^3)^-1[1+O(s)]`, `f u1'=-(5s^2)^-1[1+O(s)]`,
`g=(3pi/2)R^2[1+O(R^-2)]`, and `g'=3pi R[1+O(R^-2)]` in (11).
The two positive leading terms of `D` are `pi/(5lambda^2)` and
`3pi/(10lambda^2)`. The rational and arctangent remainders are uniform
on any sufficiently small fixed `s` range with sufficiently large `R`.
The estimate for `b` follows by cancellation of the leading matched values.

## Matched frame and the residual spin condition

Fix the RN angular frame by `omega(r_out)=0`, `r_out=1+L`. Equation (3)
and continuity of both potentials determine the previously free constant:

```text
beta=(2/lambda) integral_(1+s)^r_out u(r)/r^4 dr +2k H(R).        (13)
```

For the zero-charge regular mode define
`w0(r)=2/r-(2r-1)/r^4`. Since `w0(1)=1` and `w0'=-2u0/r^4`,
(12)-(13) imply

```text
beta=(e/lambda)[1-w0(r_out)+O(s+R^-2)],                          (14)
```

and `w0(r_out)=O(1/L)` for large `L`. One elementary remainder check is
`integral_(1+s)^r_out u1/r^4=O(s^-2)`: multiplying by the bound on `b`
controls its contribution to (13). The omitted regular near-horizon
integral and the `k H(R)` term have size `O(e s^3/lambda)`.
This gives (14) without invoking a horizon boundary condition for the
finite-annulus problem.

For the leading binary drive, `e=P E_+=2v(P/d)^2=2v chi^2`, so

```text
beta~e/lambda=chi Omega ell,
Omega_frame~chi Omega,
[(q-1)/2]|beta|=O(q^-7/2 (log q)^5)                             (15)
```

in the retained family `lambda=(16pi q^3)^-1`,
`d/ell=(log q)^-2`, `chi=lambda/(d/ell)`.
This is the electric-induced contribution relative to the chosen outer
frame. Companion gravitomagnetism and tetrad transport enter at the same
PN order; no complete binary precession coefficient is established.
The last expression is the norm of the *constant rigid-rotation*
generator in the LLL spin representation, not the full Dirac response.
The radial parts of (9), their spin-connection contributions, higher
angular modes and the quantum state still need treatment.

The spin sector needs a separate matching calculation; the regular RN
horizon frequency cannot simply be imposed throughout the open throat.
Take `J_+=J`, `J_-=-J`, no additional growing electric source, and fix
the center's rotating frame by `beta=0`. The odd-potential solution is

```text
V_spin=kS rho(1+rho^2)-2J rho,
omega_spin=-2kS rho(1+rho^2)-2J rho,
u_RN=bS u1-2J/r.                                                (16a)
```

The same flux/potential matching fixes both coefficients at `r_s=1+s`:

```text
bS u1-kS(1+3R^2)=2J(1/r_s-1),
bS f u1'-6lambda kS R(1+R^2)=-2J f/r_s^2.                       (16b)
```

Its determinant is nonzero by the same sign argument as in (11).
Expanding gives `kS~(8/15)J lambda/R`, `bS~-6J s^4`.
Reconstruction then determines the upper outer frame coefficient

```text
omega_out=lambda omega_spin(R)
           +integral_(1+s)^r_out [-6J-2u_RN(r)]/r^4 dr
         =-J[1-wKN(r_out)+O(s+R^-2)],
wKN(r)=(2r-1)/r^4.                                             (16c)
```

The lower coefficient has the opposite sign. In this specified model,
spin is conjugate to a *relative boundary rotation*. Equal fixed outer
frames set `J=0` in the asymptotic regime. A nonzero spin requires altered
boundary rotation data or a dynamical response; it is not an arbitrary
static addition with both frames held fixed. Equations (16a)-(16c) are
verified independently by inhomogeneous ODE integration as well as the
closed basis. They still do not construct the actual one-exterior binary.

For large `r_out`, the upper-port versus center frequency difference in
throat time is therefore of size `|J|/lambda`. Smallness of its LLL
rigid-rotation comparison requires `q|J|/lambda<<1`, much stronger than
`|J|<<1`. This is a condition for this particular perturbative comparison,
not an instability threshold or a universal traversability bound. It is
not a claim that the central metric has a constant shift `J/lambda`.

There is a sharper diagnostic when angular momentum is measured in
units of `hbar=1`. With physical angular momentum `S`, the calibrated
dimensionless parameter is `J=G S/P^2`. The retained leading relation
`ell=16P^3/(Gq)` gives

```text
q J/lambda=16 S,
[(q-1)/2]|J/lambda|=8(1-1/q)|S|.                                 (17)
```

Consequently a one-unit increment of physical angular momentum does not
become small in this operator-norm diagnostic as `q` grows. Large flux
alone cannot justify uniform perturbative control over such spin sectors.
This is a failure of that sufficient small-rotation estimate, not a
proof of collapse, nontraversability or inevitable spin fluctuations.

As an independent classical scale check, extremal Kerr-Newman satisfies
`M_geo^2=(P^2+sqrt(P^4+4G^2 S^2))/2`. Thus its physical energy increment
is `G S^2/(2P^3)+O(G^3 S^4/P^7)`, or `ell delta E=8S^2/q+...`.
Soft rotational sectors and the distinction between fixed spin and fixed
boundary angular potential are established ingredients of
[Iliesiu and Turiaci, sections 2.3, 3.4 and appendix A](https://arxiv.org/html/2003.02860v3).
Their isolated-black-hole ensemble and neglect of additional charged
matter are not the interacting wormhole state. We do not import their
partition function, spectrum, or uplift as a proved wormhole response.
Equation (17) makes the coupled angular-current sector a concrete next
feasibility question; it is not a standalone novelty claim.

[The checker](../checks/stationary_electric_dipole.py) verifies the full
background and linear Einstein-Maxwell tensors, regular and decaying
reconstructions, the exact bases and Wronskian, and the general-mass
reference branch. Independent ODE integrations compare the closed matching
formula on fifteen annuli. [The receipt](../evidence/stationary-electric-dipole.json)
pins the inputs; [the source review](../evidence/stationary-electric-dipole-source-review.json)
records exact editions and review limits. No new stability, quantum-source,
transport or shortcut result follows from these classical checks.

Remaining to-do list: derive the quantum rotational sector coupled to the
LLL angular current and determine the admissible spin sector; obtain
transported frames and holonomy; complete smooth metric/source matching
and the full Dirac response to the stationary fields; control the evolving
quantum state and renormalized exterior response; finish the consequential
endpoint and novelty review, independently verify, GCP-certify suitable
mathematics and prepare the manuscript.
