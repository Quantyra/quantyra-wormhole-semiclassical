# Prepared-input scattering through a finite mouth collar

2026-10-10. S052/S054. In the leading flat point-dipole Weyl model, the
deformed angular modes now give a **uniform operator bound for prepared
inputs across a finite spherical annulus**. The proof retains all angular
excitations in the selected azimuthal sectors. It bounds the radial
corrector derivative and uses the current identity to control the actual
annular boundary-value solution, rather than identifying a small trial
correction with a scattering error.

This closes a local propagation step. The annulus has independently
specified incoming data at both surfaces. Attaching the outer surface to
the rest of the exterior changes those data and introduces feedback. The
result is not yet the two-mouth exterior return, curved throat matching,
or an energy/shortcut theorem. The classical boundary and dissipativity
methods used here are established; no novelty acceptance is claimed.

## Model and statement

Keep all conventions and ordinary estimates of
[finite spherical mouth modes](mouth-mode-matching.md). The monopole
separation is `d`, the annulus is `rho_0<=rho<=rho_1`, and

```text
0<epsilon_0=rho_0/d<epsilon_1=rho_1/d<=1/10,
q>=3 integer, Xi=xi d>0,
H=-i sigma_z partial_rho+A(rho), A=L(rho)/rho.
```

Write the radial spinor as `(u,v)`. Inward data for the **annular domain**
are `(u(rho_0),v(rho_1))`; outward data are
`(v(rho_0),u(rho_1))`. Norms are the angular Hilbert norms of the unitary
spin/measure frame, so squared boundary norms are fluxes. Smooth spinor
regularity in gauge patches at the poles is retained. No wall is placed
at an angular cutoff.

For each `ell=m+1/2`, `m=0,...,q-1`, let `g_m(rho)` be the exact normalized
upper angular zero mode and `eta_m=(0,h_m)` its real lower-spin corrector,

```text
A eta_m=i(partial_rho g_m,0).
```

Let `S_xi` map incoming to outgoing boundary data for `(H-i xi)psi=0`.
Let `J_0: C^q -> incoming space` insert `sum alpha_m g_m(rho_0)` at the
inner surface and zero lower-spin input at the outer surface. Define

```text
C_xi alpha = sum_m alpha_m
  (h_m(rho_0), exp[-xi(rho_1-rho_0)]g_m(rho_1)).        (1)
```

Then the ordinary bound is

```text
||S_xi J_0-C_xi|| <= E(q,epsilon_0,epsilon_1,Xi),

E^2 = (729/160)epsilon_1^4 exp[-2Xi(epsilon_1-epsilon_0)]
 + (1/(2Xi)) integral_epsilon_0^epsilon_1
       exp[-2Xi(t-epsilon_0)]
       [(27/5)t+12 sqrt(q)t^3+(22/5)Xi t^2]^2 dt.      (2)
```

There is no extra factor `sqrt(q)` from superposing the inputs: the
azimuthal sectors and their residuals are orthogonal. The explicit
`sqrt(q)` already present in (2) is a radial-derivative cost. Equation
(2) controls the full outgoing angular state, including reflected and
transmitted higher modes, for these prepared inputs.

A conservative exact consequence is

```text
3<=q<=10^12, epsilon_1<=1/12500, 0.05<=Xi<=5
    => E<8*10^(-6), uniformly in 0<epsilon_0<epsilon_1. (3)
```

For equal spherical ports with `D=d-2rho_1`, the previous window
`0.05<=xi D<=4` lies in the `Xi` interval used in (3). This is a numerical
budget for the finite collar only. It is not a bound on composition with
the exterior, and the required gravitational overlap
`r_e << rho_0 < rho_1 << d` must still be justified in a selected physical
parameter family.

## The radial derivative estimate

Use `x=sin^2(theta/2)`, `epsilon=rho/d`, `W=q w`, and
`D_s=S(S+1-epsilon c)` from the preceding note. The subscript distinguishes
this algebraic denominator from an exterior front time. Differentiating
at fixed `x` gives

```text
S_epsilon=(epsilon-c)/S,
S_epsilon_epsilon=(1-c^2)/S^3,
|D_s,epsilon_epsilon| <=2420/243,
|w_x,epsilon_epsilon| <=9/4.                          (4)
```

Indeed, the previous bounds `D_s>=81/50`, `|D_s,epsilon|<=77/15`,
`|S_epsilon|<=11/9`, `|S_epsilon-c|<=20/9` imply
`|D_s,epsilon_epsilon| <= (33/10)(1000/729)+440/81`.
Twice differentiating `w_x=epsilon^2/D_s` bounds the last expression
in (4) by

```text
2/(81/50)+4(1/10)(77/15)/(81/50)^2
 +(1/100)(2420/243)/(81/50)^2
 +2(1/100)(77/15)^2/(81/50)^3
 =10421500/4782969 <9/4.                              (5)
```

All expectations below use `g(rho)^2 dtheta`, or equivalently its `x`
density. It has Poincare constant at most `C_P=2/q`. Put
`a=W_rho-E W_rho` and `v=E a^2`. Normalization differentiation yields

```text
g_rho=-a g,
g_rho_rho=[a^2-2v-(W_rho_rho-E W_rho_rho)]g.          (6)
```

If `|a_x|<=L_a`, Poincare applied first to `a` and then to `a^2` gives
`v<=C_P L_a^2` and `E a^4<=5 C_P^2 L_a^4`. Also
`E[(a^2-2v)^2]=E a^4`. With the previous
`L_a<=(3/2)q rho/d^2` and (4), this proves

```text
||g_rho_rho|| <= (9 sqrt(2)/4) sqrt(q)/d^2
                   +(9 sqrt(5)/2)q rho^2/d^4.         (7)
```

The angular coefficient has the particularly simple derivative

```text
s_rho=-q rho sin(theta)/(2 d^2 S^3),
||s_rho||_infinity <=(500/729)q rho/d^2.               (8)
```

Differentiate `A eta=i(g_rho,0)`, using
`A_rho=-A/rho+sigma_y s_rho/rho`:

```text
A eta_rho = i(g_rho_rho+g_rho/rho,0)
                         -sigma_y s_rho eta/rho.      (9)
```

The derivative `eta_rho` is lower-spin and hence orthogonal to the kernel.
Smooth parameter dependence follows also from its explicit cumulative
integral formula. Apply the compact selected-sector gap from the previous
note, as well as its `||g_rho||` and `||eta||` bounds. This yields

```text
||eta_rho|| <= [135/(8 sqrt(10))]rho/d^2
                +[81/8+25/(12 sqrt(2))]sqrt(q)rho^3/d^4
             <= (27/5)rho/d^2+12 sqrt(q)rho^3/d^4.     (10)
```

This is a derivative estimate for the actual angular inverse, with no
derivative of a numerically truncated inverse used in its proof.

For the trial spinor

```text
psi_app(rho)=exp[-xi(rho-rho_0)] [(g(rho),0)+eta(rho)],
f_app=(H-i xi)psi_app
     =i exp[-xi(rho-rho_0)] (eta_rho-2xi eta),          (11)
```

equations (10) and `||eta||<=(11/5)(rho/d)^2` bound its residual by the
polynomial appearing in (2). The integral in (2) is the squared bulk
`L2(d rho)` residual bound divided by `2xi`, after rescaling by `d`.

## Boundary-value existence and the current estimate

The argument below fixes the annular domain and justifies the relevant
solution and trace maps. It does not assume a spectral gap for the full
noncompact exterior.

First impose zero incoming data. For the isolated monopole, the compact
angular Dirac operator anticommutes with `sigma_z`. Its spectral
decomposition gives the zero modes and paired nonzero angular modes.
Each nonzero pair reduces the radial operator on this finite annulus to a
two-component ordinary differential system with continuous coefficient
`lambda/rho`. A zero mode gives one scalar radial equation. In every
finite block the two incoming endpoint conditions, or the single
condition for a zero mode, determine a unique solution at `i xi`:
uniqueness follows from the current estimate below, and invertibility of
the finite-dimensional boundary linear system supplies existence.

The same estimate gives resolvent norm at most `1/xi`, uniformly over
the angular blocks. Their direct sum therefore defines a closed isolated
operator `H_0` with bounded inverse `(H_0-i xi)^(-1)` on bulk `L2`.
Outgoing traces are continuous in its graph norm by the flux estimate.

The full point-dipole operator on the annulus differs by the bounded
self-adjoint multiplication operator
`B=sigma_y(s-s_0)/rho`: the remote field is smooth, the gauge difference
vanishes at the poles, and `rho_0>0`. For `xi>||B||`, the Neumann series
against the isolated resolvent gives the inverse of `H=H_0+B` on the same
zero-incoming domain. The current estimate still gives norm `1/xi`.
The resolvent identity continues this inverse to smaller positive `xi`
in steps of less than the current `xi` (halving suffices); at each step
the current estimate restores the `1/xi` bound. This proves existence for
every positive `xi`, without importing a general Dirac boundary theorem.

Smooth incoming data can be lifted to a smooth radial interpolation of
the endpoint spinors. Subtracting the zero-incoming resolvent solution of
its residual gives the desired homogeneous solution. The flux estimate
extends the resulting boundary map by density to angular `L2` incoming
data. Boundary traces in this completion are the limits of smooth traces,
not arbitrary point evaluations of a bulk `L2` representative.

For a smooth solution of `(H-i xi)psi=f`, self-adjointness of `A(rho)`
implies, with the inner product conjugate-linear in the first entry,

```text
j(rho)=<psi,sigma_z psi>,
j'=-2xi ||psi||^2-2 Im<psi,f>,
||out||^2+2xi ||psi||_bulk^2
       =||in||^2-2 Im<psi,f>_bulk.                    (12)
```

The identity extends to the graph/trace completion just constructed.
For `f=0`, `S_xi` is a contraction. With nonzero `f`, Cauchy-Schwarz and
maximizing `-2xi z^2+2z||f||` over `z>=0` give

```text
||out||^2 <= ||in||^2+||f||_bulk^2/(2xi).              (13)
```

Compare the actual solution with inner input `g(rho_0)` and zero outer
input to (11). Their incoming difference is only
`-exp[-xi(rho_1-rho_0)]h(rho_1)` at the outer surface. Their residual
difference is `-f_app`. Applying (13), using
`||h(rho_1)||^2 <=(729/160)epsilon_1^4`, proves (2) in one sector.
Orthogonality and the same constants in every selected sector prove the
operator bound on all of `C^q`.

To verify (3) without uncertified quadrature, drop the exponential factors,
enlarge the integral to `[0,epsilon_1]`, use `1/(2Xi)<=10`, and replace
the polynomial by
`(27/5)t+22t^2+12000000t^3`. Its squared polynomial integral plus
`(729/160)epsilon_1^4`, evaluated at `epsilon_1=1/12500`, is an exact
rational number smaller than `(1/125000)^2`. The receipt retains that
rational certificate. All bounds are monotone in the enlarged parameters.

## What this does and does not close

The outgoing upper component differs from the corrected-mode transmission
in (1) by at most `E`. The reflected component differs from `h(rho_0)`
by at most `E`, so its norm is at most
`(27/sqrt(160))epsilon_0^2+E`. These statements concern unit incoming
flux in a known finite-dimensional port subspace, not just a prepared
wave packet with an unspecified seminorm constant.

Setting the outer incoming component to zero selects a column of the
annular scattering map. The surface is an interface for defining incoming
and outgoing data; no absorbing detector or wall is added to the physical
wormhole. When the exterior is attached, its return can feed that input.
A bound for one column and contractivity of the full map do not by
themselves bound the resulting feedback inverse. The remaining model
closure needs the exterior load, continuum escape, opposite-mouth
extraction, and curvature/state corrections, including their dependence
on throat length in the energy derivative. Traversability and a
fixed-observer travel-time advantage remain separate questions.

## Literature and independent checks

[Tumulka, sections 2.1-2.2](https://arxiv.org/html/1601.04571v1#S2)
uses the positive outward-normal Dirac projector and the corresponding
current balance for detection boundaries. Those sections were read;
their stated existence result refers to a separate dependency. We use
the current/projector formalism, not the proposed detector interpretation.
[Frolov, Teufel and Tumulka v5](https://arxiv.org/html/1912.12057v5)
has a Schrödinger existence theorem (section 2.1); its introduction defers
the analogous Dirac investigation. The older search abstract describing
Dirac existence is therefore insufficient support for that claim in the
latest edition. Section 2.1, the relevant introduction passages and
section 3's semigroup theorem statement were read. Their complete
Schrödinger proof and the earlier-edition changes were not imported.
Our annular existence argument above supplies the particular construction
needed here. General dissipative boundary methods are prior art.

Run

```powershell
python checks/mouth_collar_scattering.py --output evidence/mouth-collar-scattering.json
```

The [receipt](../evidence/mouth-collar-scattering.json) records symbolic
identities, exact rational constants, and independent current checks.
It computes the corrector from its cumulative integral and checks its
radial derivative at two step sizes. A separate Jacobi angular expansion
solves the coupled annular boundary problem at three angular orders,
checking current balance, passivity, convergence and the error bound.
The numerical systems retain increasing sets of angular excitations;
their agreement supplements the untruncated ordinary proof and is not
an interval-certified infinite-dimensional simulation. No Lean was run.

Remaining to-do list: derive the exterior load and control composition
with both mouth collars, justify curved throat and state matching, bound
total energy and common-observer corrections, complete decisive source
and novelty selection, then certify suitable endpoints on GCP and prepare
the reproducible manuscript.
