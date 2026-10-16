# Finite input identifiability on the retained clock

S132, E003, 2026-10-15. This is a finite algebraic closure decision for the
derivative-only S104–S130 continuation. S131 is a preserved, unaccepted
twenty-seventh-input draft. This note supersedes the instruction to continue
that derivative sequence. It supplies no new physical measurement or novel
wormhole theorem and does not complete S074.2–S074.6.

## Scope and test fixed before execution

The comparison is the retained beta function, prescribed-angle clock sum,
parallel electric susceptibility, and absence of exterior channel lengths
from that clock. The hardest obligation is distinguishing dependence of a
formula on a parameter from independently supplied data identifying it.
An admissibility interval is a constraint, not a measurement selecting a
coefficient. We prove the structural ambiguity once for every derivative
order, rather than testing one more derivative as a candidate coupling.

Sources are [the running-coupling note](source-clock-running-inputs.md),
[the controller reduction](source-clock-controller.md),
[the joint clock checker](../checks/source_clock_joint.py),
[the finite electromagnetic action](finite-local-electromagnetic-input.md),
and [the finite mouth determinant](finite-mouth-determinant.md).
The checker hashes these inputs, this note and itself. It does not modify
or run the historical clock checkers, integrate a new trajectory or fit a
coefficient. The two uses of kappa are distinguished below.

Acceptance is exact symbolic identities plus independent rational arithmetic
examples. A 60-decimal numerical central difference, step 1e-8 and absolute
tolerance 1e-14, checks the first mu and mu0 partial derivatives and the susceptibility
variation independently of symbolic differentiation. The all-order result
is justified by induction, not by the finite numerical tests. Rational
examples are synthetic demonstrations, not physical parameter assignments.

## The scale derivatives cannot identify a finite constant

At fixed named scale mu0>0 let C=c_theta(mu0), b=1/(20 pi^2), and mu>0.
The retained equation mu*d_mu c_theta=-b has the full solution

```
c_theta(mu)=C-b*log(mu/mu0).
d_mu^n c_theta(mu)=(-1)^n b (n-1)! / mu^n, n>=1.
```

The base case is -b/mu. Differentiating the displayed nth expression
gives the (n+1)st expression, proving the formula for all positive integers.
Every such derivative is independent of C. Replacing C by C+a changes
the coupling by a and preserves the beta function and all positive-order
scale derivatives. Thus these derivatives alone cannot select C, even
when their numerical values are known exactly. This does not identify C
with b, a heat coefficient or zero. A finite renormalization condition or
an independently matched physical observable at a specified scale and
scheme is required. Choosing a subtraction convention alone does not
establish a physical value while other matched terms are held fixed.

The S105–S130 loop instead differentiates the logarithmic endpoint in mu0
while holding C as an independent symbol. Its all-order formal partial is

```
partial_mu0^n [C-b*log(mu/mu0)]
    =b*(-1)^(n-1)*(n-1)!/mu0^n, n>=1, with C held fixed.
```

The base is b/mu0; differentiating the nth expression gives the (n+1)st.
This also eliminates C at every positive order and cannot identify it.
These partial derivatives are not physical total derivatives under a
change of reference scale. For a fixed physical coupling the constant is
itself a running reference value, dC(mu0)/dmu0=-b/mu0. Its contribution
cancels the logarithm's b/mu0, so the total reference-scale derivative
of the same coupling vanishes. Neither calculation supplies a finite
matching condition. The checker verifies both inductions and this
cancellation, with an independent numerical mu0 partial check.

## What one aggregate clock coefficient identifies

Use the controller note's homogeneous, clock-independent Lhat^{uu}=L
reduction at a common matching scale, with A=pi*ell^2/N>0. Then

```
Lambda=A*(zeta_s+2*C*L),
d_(zeta_s,C) Lambda=A*[1,2*L].
```

This row has rank one. Its kernel contains (-2*L,1), and the exact shift
(zeta_s,C)->(zeta_s-2*L*a,C+a) leaves Lambda unchanged. If L=0, C is
entirely absent. If L is unknown, the missing geometric input adds to the
ambiguity. The realized-drive equation

```
z_tt=g/2-(1+Lambda*r0^2)*z
     -4*Lambda*kappa_clock*r0*z^2-3*Lambda*kappa_clock^2*z^3
```

therefore cannot distinguish this family at fixed other inputs. This is
a statement about the retained aggregate model, not the fully varied
finite boundary action. Lapse-dependent L or unmatched geometry needs a
separate action variation; treating L as fixed does not prove that scope.

For two independently matched geometries at the same scale and scheme,
assume common couplings zeta_s,C and known A_i!=0,L_i. If measured or
matched Lambda_i are supplied, set q_i=Lambda_i/A_i. Then

```
det [[A1,2*A1*L1],[A2,2*A2*L2]]=2*A1*A2*(L2-L1),
C=(q2-q1)/(2*(L2-L1)),
zeta_s=(L2*q1-L1*q2)/(L2-L1).
```

Identification holds exactly when L1!=L2 under these assumptions.
Repeated measurements with the same L retain the kernel. Near-equal
weights amplify input errors: with exact weights and |delta q_i|<=e_i,
|delta C|<=(e1+e2)/(2*|L2-L1|). Uncertain geometry requires its own bound.
The retained clock windows do not supply these two independent matched
observations. We do not assert that such physical geometries or data exist.

## The electromagnetic coefficient needs independent response data

Write kappa_grav^2=8 pi G; it is not kappa_clock. The retained local term
is Delta L=cJ*kappa_grav^4*J^2/16. A purely magnetic reference has J=0;
both this term and its first variation vanish for every cJ. Those
reference data alone cannot identify cJ. Its parallel electric Hessian is

```
chi=1+2*cJ*kappa_grav^4*B0^2,
cJ=(chi-1)/(2*kappa_grav^4*B0^2).
```

The inversion requires independently specified chi at known nonzero B0
and kappa_grav, with other operators, normalization, scale and scheme
matched. At B0=0 all cJ give chi=1. Substituting cJ=2 and differentiating
the resulting polynomial cannot infer cJ: the substitution removed it.
A measured response derivative could identify cJ when its coefficient is
nonzero; differentiating a formula does not supply that measurement.
A positive quadratic response is not a proof of nonlinear stability.

## Channel lengths absent from the clock remain unassigned

The retained acceleration above contains no D_j. Any two positive channel
length assignments therefore give the same retained clock expression.
Its D_j derivative is zero for that reason; this is not a statement that
the physical nonlocal determinant is D_j-independent. For example, the
mouth note's eta*coth(k*(T+D_j)/2)/(2*k) changes with D_j. Its nonzero
derivative demonstrates dependence but does not assign D_j. A matched
exterior return map, geometry and time normalization must supply those
lengths and controlled response bounds. Differentiating a ceiling in
a,h,epsilon cannot replace this missing match.

## Decision and remaining obligations

The four algebraic audits close derivative-only coefficient identification
as a continuation route. They identify the missing matching inputs without
determining their physical values. Do not open a twenty-eighth derivative
story. Preserve S104–S131 artifacts and failed assumptions.

S074 still needs its declared-order action, ensemble, constraints, complete
local matching and observable remainder. Unequal driven collars and the
nonlocal interacting determinant remain proof obligations. A separate
response-control result must bound the derivatives used in the stress and
geometry equations, rather than only a determinant value. S075–S079 remain
behind their existing gates; the area budget is unevaluated and the latest
S084 intervals are unchanged. No opening, stability, transport, preparation
or creation result is claimed. No Lean or paid compute is used.

Remaining to-do list: supply or match finite renormalization and response
inputs; derive unequal-collar and interacting-state control; complete
S074's full pass gate before the dependent physical-result tasks.
