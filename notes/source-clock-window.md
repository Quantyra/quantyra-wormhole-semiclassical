# Admissible range for the varied-clock source contact

S081, E003, 2026-10-14. The opening benchmark varies the clock `y`
while the switch angle remains the prescribed function in that
benchmark. The contact charge therefore enters the clock constraint.
This is rule D. The fixed-metric matter statement stays in the domain
where the metric is not varied. Setting `lambda=0` is that fixed-metric
limit. It is not selected as the opening value.

## Charge identity

The contact Lagrangian `L=A y^2 (dTheta/dt)^2`, with
`A=2 pi zeta_s ell^2`, has u-energy `H_u=L` and global-time charge
`E_t=A y (dTheta/dt)^2`. The clock coupling is
`lambda=A/(2N)=pi zeta_s ell^2/N`. The linear source carries a factor
`eta^2`. Both acceptance runs check these identities.

## Benchmark window

The sample is `eta=1/2000` and the six retained initial corners
`z(0) in {1/4, 7/25}`, `z_t(0) in {-1/2000, 0, 1/2000}`. At each corner
the checker finds the largest positive and negative `lambda` whose
switch, post-switch, and final-energy budgets all hold. The switch
budget is the first to meet its threshold. At half that eta, every
sampled linear threshold is larger.

The intersection of the six admissible endpoints is

```
-4757.61079788208 <= lambda <= 4875.4563331604.
```

The most restrictive corner is `z(0)=7/25`, `z_t(0)=1/2000`. Its linear
threshold is `4815.979518870512`, its nonlinear positive endpoint is
`4875.4563331604`, and its negative endpoint is `-4757.61079788208`.
At the positive endpoint the switch deviation is
`0.0049999996964726935`, the post-switch amplitude is
`0.004360059005706085`, and `|delta e_f|` is
`0.0037124239313687524`. An independent RK45 integration of that
endpoint agrees with DOP853 to `1.4938819625776034e-10`. On this same
corner the switch deviation increases through one quarter, one half,
three quarters, and the endpoint, on both signs. The other corners were
bracketed, and their interiors were not probed one by one.

`zeta_s=lambda N/(pi ell^2)` on this interval. No retained identity
selects a point inside it. The area budget is not evaluated. `Phi`
still depends on `y` and `E`.

## What remains open

S074.2–S074.6 stay open. This interval is a range for one finite
coefficient on the benchmark sample. It does not match unequal collars,
the nonlocal interacting determinant, `cJ`, the finite Weyl contact
`c_theta`, or a controller whose angle depends on the clock. S075
through S079 stay blocked. This is not an impossibility theorem and not
a physical opening.

[Checker](../checks/source_clock_window.py),
[receipt](../evidence/source-clock-window.json), and
[source/run review](../evidence/source-clock-window-review.json).
No Lean, paid compute, publication, or outreach was used.

Remaining to-do list: carry `c_theta` and a non-prescribed controller
angle onto this clock. S074.2–S074.6 stay open. S075 through S079 stay
blocked.
