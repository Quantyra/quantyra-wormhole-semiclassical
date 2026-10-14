# Source-contact allocation into the opening clock

S074, E003, 2026-10-13. **Scoped obstruction.** The finite source
coefficient is not bounded in the observables that feed S076. Two
extensions of the retained contact are consistent with the completed
identities, and they disagree on the switch clock by more than the
frozen budget. This is an unresolved physical input. It is not an
impossibility theorem for wormholes, and it is not a physical opening.

## Comparison, hardest step, observable, and test

These were fixed before the allocation below.

The closest comparison has three pieces. The returned-switch benchmark
solves `y_tt+y=-E(t)/(2N)`, or `z_tt+z=g(t)/2` with `z=y/eta`, on
`0<eta<=1/2000`, `1/4<=z(0)<=7/25` and `|z_t(0)|<=1/2000`. It keeps
`z>1/100` on the switch, `z>1/25` afterward and `e_f<-1/25` at
`T_star=pi+1/4`. The completion plan asks later corrections to obey
`|delta e_f|<1/50`, switch `sup|delta z|<1/200` and post-switch
`sup|delta z|<1/50`. The curved-collar contact (5) of
[curved-fermion-boundary.md](curved-fermion-boundary.md) is proved for a
prescribed angle on a fixed boundary metric. Its stored energy returns
to the initial value when the angle derivative vanishes at both
endpoints, so that contact does not change endpoint matter excitation.
The charge dictionary is the one already used for returned radiation:
the global-time charge is the boundary energy divided by `y=dt/du`.

The hardest step is whether that contact's stress enters
`Q_0,grav+q_0=0` once the boundary metric is built from the clock.
The observable is the resulting `delta z` and `delta e_f` on the
benchmark window. The test is direct variation of the contact
Lagrangian, an independent variation-of-constants evaluation, and a
second ODE method. A pass would be an admissible range of the finite
coefficient on which every budget holds. A vanishing endpoint contact
energy, a positive quadratic mass or a small bulk loop is not that range.

## Two extensions that the completed identities allow

Write the renormalized boundary metric as `gamma=-du^2+ell^2 dOmega^2`
and the prescribed global-time angle as `Theta(t)`, with `y=dt/du`.
The accepted lapse form of the contact, at lapse `1/y` and sphere
radius `ell`, is the u-time Lagrangian

```
L=A y^2 (dTheta/dt)^2,    A=2 pi zeta_s ell^2.          (1)
```

Its u-energy and global-time charge are

```
H_u=y dL/dy-L=A y^2 (dTheta/dt)^2,
E_t=H_u/y=A y (dTheta/dt)^2.                            (2)
```

The factor `1/y` is the same redshift that takes a boundary power to
the global-time charge in the returned-radiation identity. It is not a
new normalization of `N`.

Rule P leaves the boundary metric fixed and does not insert `E_t` into
the clock constraint. The matter endpoint excitation is unchanged when
`dTheta/dt` vanishes at the endpoints. The clock correction from this
contact is then zero. Choosing rule P sets `zeta_s` out of the clock
equation by hand. The heat polynomials do not imply that choice. Rule P
also does not cover a controller for which `Theta` depends on the clock.

Rule D inserts `E_matter+E_t` into the retained constraint
`y_tt+y=-E/(2N)`. With `E_matter=-N eta g`, `y=eta z` and
`lambda=A/(2N)=pi zeta_s ell^2/N`,

```
z_tt+(1+lambda (dTheta/dt)^2) z=g/2.                    (3)
```

Both rules reproduce the fixed-metric Ward identity and the curved-collar
derivative heat polynomials. Those polynomials are functions of the
second fundamental form and the angle gradient only. `zeta_s` is not
one of their variables. No completed axial, electric or reference-radius
bound supplies `|lambda|`.

## The endpoint cancellation is not the clock bound

On the shaped switch `Theta_t=4 eta[S(t/tau)-S((t-pi)/tau)]`, with the
quintic `S` and `tau=1/4`, both endpoint rates vanish. Therefore `H_u`
vanishes at `t=0` and at `T_star` for every `zeta_s`. The final matter
excitation is not shifted by the stored contact energy. The clock is a
driven oscillator, so an intermediate change of its frequency still
changes `z(T_star)` and `z_t(T_star)`. Equation (3) keeps that memory
after `H_u` has returned to zero. Reporting the endpoint cancellation
as `delta z=0` substitutes a different observable for the budget.

## Budget test on the benchmark window

The free drive `g=a-r` is the returned-switch drive. Six trajectories
at `eta=1/2000`, the initial-window corners `z(0) in {1/4,7/25}` and
`z_t(0) in {-1/2000,0,1/2000}`, reproduce `z>1/100` and `e_f<-1/25`.
Their final reduced energies run from about `-0.06635` to `-0.04745`,
and the smallest switch value is about `0.02535`, in agreement with the
retained benchmark sample. Differentiating (3) at `lambda=0` gives

```
z1_tt+z1=-(dTheta/dt)^2 z_0.                            (4)
```

A Dormand-Prince integration, an independent Runge-Kutta 4(5)
integration and the variation-of-constants formula at `T_star` agree
to better than `1.09e-10` on the free solution and `3.47e-11` on `z1`.
The work identity `e_t=-g_t z` for rule P, and its rule-D extension
with the extra power `2 lambda (dTheta/dt)^2 z z_t`, agree with an
independent quadrature.

On these six trajectories the switch response is the most restrictive.
`sup|z1|` lies between `9.512e-7` and `1.038e-6`, so the switch budget
`1/200` fails for `lambda` above a threshold between `4816` and `5257`.
The post-switch amplitude reaches the budget `1/50` only for `lambda`
above a threshold between `21974` and `23588`. The final-energy
derivative reaches `1/50` only above a threshold between `26215` and
`33576`. The heat calculation fixes none of these three numbers.

The most sensitive corner in the sample is `z(0)=7/25`,
`z_t(0)=1/2000`, with switch threshold `4815.9793639883765`. At twice
that value, direct integration of (3) produces

```
sup|delta z|=0.009762475724334141 > 1/200
```

on the switch. At this same finite value the post-switch amplitude is
`0.008472084538506725<1/50` and `|delta e_f|=0.00731626457742085<1/50`.
Those two budgets are not saved: each fails at its own larger finite
threshold, and neither threshold is an output of the local heat
calculation. The area budget is not reached. `Phi` depends on `y` and
on `E`, and both remain unallocated while `lambda` is free.

Rule P and rule D therefore disagree on a budgeted observable for an
open set of `zeta_s=lambda N/(pi ell^2)`. No retained identity selects
the rule or the range.

## Decision

Keep S074.2 through S074.6 open. The precise unresolved input is the
insertion rule for this finite source stress, together with an
admissible range of `zeta_s` after that rule is chosen. Setting
`lambda=0`, or stopping at the endpoint cancellation, does not close
the gate. The same obstruction blocks S075 through S078: there is no
accepted local model whose residuals fit the opening budgets, and a
conditional estimate is not a physical opening. S077 and S078 are not
given a geometry. Preparation and creation stay separate and unachieved;
the bridge remains an input, and the fixed-manifold equations do not
supply a topology change.

Other unmatched inputs stay in force, including unequal collars, the
nonlocal interacting determinant, `cJ`, the finite Weyl contact and a
dynamical controller whose angle is not a prescribed c-number. This
note does not convert those gaps into a no-go.

[Checker](../checks/source_clock_allocation.py),
[receipt](../evidence/source-clock-allocation.json) and
[source/run review](../evidence/source-clock-allocation-review.json).
No Lean, paid compute, publication or outreach was used.

Remaining to-do list: choose the source-stress insertion rule and bound
or match `zeta_s` in the opening observables; unmatched boundary,
nonlocal and controller terms remain. S075 state and drive, S076 stress
and geometry, S077 stability, S078 transport, and S079 preparation and
creation stay blocked by that open S074 gate. S080 records this endpoint.
