# Cauchy solution of the traceless retarded JT constraints

2026-10-12. This solves the leading `l=0` JT constraints for any smooth,
conserved, traceless two-dimensional source on fixed global `AdS2`, in terms
of initial data and three embedding-mode constants. It applies directly to
the fixed-throat retarded LLL source in
`finite-cut-lll-boundary-map.md`. It does not fix the constants from the
curved mouths or include the traceful absolute anomaly and massive tower.

## Null constraints and their initial-data reduction

Use unit-radius global conformal coordinates

```text
ds2^2=sec^2(x)(-d tau^2+dx^2),
U=tau-x,  V=tau+x,  -pi/2<x<pi/2.
```

For the traceless state-dependent LLL perturbation let
`S_UU=F(U)`, `S_VV=G(V)`, `S_UV=0`. The source is conserved because each
null component depends only on its own null coordinate. The JT equation
`L_ab[Phi]=S_ab`, with
`L_ab=g_ab Box Phi-nabla_a nabla_b Phi-g_ab Phi`, has components

```text
L_UU=-(Phi_UU+tan(x) Phi_U)=F(U),
L_VV=-(Phi_VV-tan(x) Phi_V)=G(V),
L_UV=Phi_UV+(1/2)sec^2(x)Phi=0.                            (1)
```

The last equation is `Box Phi=2 Phi`. On the initial slice `tau=0`, set
`phi_0(x)=Phi(0,x)`, `pi_0(x)=partial_tau Phi(0,x)`,
`H(x)=F(-x)+G(x)`, and `J(x)=F(-x)-G(x)`. Combining the two null
constraints with the trace equation gives the complete initial-data system

```text
phi_0''-tan(x)phi_0'-sec^2(x)phi_0 = -H(x),
pi_0'-tan(x)pi_0 = J(x).                                  (2)
```

Let `P(x)=int_0^x H(s) ds` and `Q(x)=int_0^x sin(s)H(s) ds`. Variation
of parameters, using the homogeneous solutions `tan(x)` and `sec(x)`, gives

```text
phi_0(x)=C_1 tan(x)+C_2 sec(x)-tan(x)P(x)+sec(x)Q(x),
pi_0(x)=sec(x)[C_3+int_0^x cos(s)J(s) ds].                 (3)
```

Direct differentiation verifies (2). The constants `C_1,C_2,C_3` give
the full source-free family
`Phi_h=C_1 tan(x)+C_2 sec(x)cos(tau)+C_3 sec(x)sin(tau)`; they are the
three constrained embedding modes already found in the local audit.

## Evolution and constraint propagation

Evolve (3) with the trace equation

```text
partial_tau^2 Phi-partial_x^2 Phi+2 sec^2(x)Phi=0.          (4)
```

For smooth `F,G` and finite-cut endpoints `|x|<=a<pi/2`, this is an
ordinary linear hyperbolic initial-boundary problem. Equation (2) makes
both null constraints hold on the initial slice. Their propagation follows
from `nabla^a L_ab=0` on unit `AdS2`, `nabla^a S_ab=0`, and (4). Thus the
retarded source determines the interior dilaton response once the three
initial embedding data and the boundary conditions at the two cuts are
fixed. The latter conditions must come from the matched mouth mass/time
data; regularity at ideal `AdS2` boundaries alone is not a substitute for
the asymptotically flat matching.

For the finite-cut channel source, evaluate `F,G` from (10) of
`finite-cut-lll-boundary-map.md`, summing over the actual return channels.
Equation (3) is then a direct source-to-dilaton map. It closes the local
constraint solve for prescribed characteristic histories but does not
close the feedback problem: the characteristic maps themselves depend on
the moving clock, while the clock is determined by the matched JT
boundary equation and the exterior mass/charge balance.

## Verification and limits

`checks/jt_retarded_source_cauchy.ps1` tests (2)-(3) on deterministic
analytic source profiles using independently differentiated expressions.
The test covers the sourced spatial equation, the initial momentum
constraint, and the three homogeneous modes. It verifies the local
ordinary-algebra formulas, not evolution, global matching, or the physical
LLL state. No new stability result or novelty claim is made.

Remaining to-do list: impose both physical mouth boundary conditions on
(4); solve the coupled clock/source feedback with the RN mass/charge
balances; include the complete renormalized Dirac response; then determine
whether the selected stability target survives.
