# Dipole field line delays and the conditional shortcut margin

2026-10-10. S051/S052/S054. The leading Gaussian length-energy comparison
can retain a distribution of exterior delays even when channels mix. Under
the exterior-power bound below, the unmixed comparison maximizes its positive length
derivative. In the point-dipole continuum this gives a stationary transit
ratio strictly greater than a defined constant `t_* = 2.1489119899...`.
It also gives a correction budget proportional to `q/d^2`, which remains
useful when the previous shortest-return budget `2A/d^3` becomes small.
These are conditional model results. The physical projected propagation law,
four-dimensional energy and actual omitted derivatives remain unproved.

MMP already obtains the stronger numerical ratio approximately `2.35` for
its unmixed dipole model. We reproduce `2.3471178600` independently. The
candidate addition is an analytic, mixing-tolerant comparison and error
budget, not discovery of the original no-shortcut conclusion. The matrix
comparison uses the established finite-dimensional spectral-majorization
tools of Weyl and Pólya; the delay inequalities themselves must still be
proved for the propagation model.

## Model and delay hypothesis

Retain the two-chirality Gaussian energy and its subtraction convention
from [the passive-return argument](passive-return-bound.md):

```text
E(T) = A/T^2 + q pi/(24T) + E_loop(T),    A>0,
E_loop = -(1/pi) sum_sigma integral_0^infinity
                   Re log det[I+exp(-xi T) R_sigma(xi)] dxi.
```

Specify positive delays `d_j`, independent of `T`, and require each return
to admit the physical half-delay factorization

```text
R_sigma(s) = H(s) G_sigma(s) H(s),
H(s) = diag(exp(-s d_j/2)),
||G_sigma(s)|| <= 1, Re s>0.                            (1)
```

The cores are analytic, passive and independent of `T`; they may mix
channels, depend on frequency and lose amplitude within the completed
Gaussian model. Constant input/output unitaries and additional passive
factors do not affect the exterior-power estimates if placed outside the
two delay banks. Equivalently, the following weaker condition suffices
on the imaginary-frequency axis, with delays in increasing order:

```text
||wedge^k R_sigma(xi)|| <= exp(-xi sum_{j=1}^k d_j),
                                       k=1,...,q.       (2)
```

Equation (1) is realized by a prescribed network with two mirrored banks
of free propagation and an intermediate passive core. It is a substantive
geometric assumption. A common front `d` and passivity alone imply only
the earlier equal-delay bound. For example `R(s)=exp(-sd)I` satisfies both
of those conditions but violates (2) when some assigned `d_j>d`.
Arbitrary mixing distributed throughout an exterior need not factor this
way. A local chiral propagation model implying (2) without factorization
is derived below. No same-mouth reflections are added to the model.

## Matrix comparison with all delays retained

Set `X=exp(-xi T)R_sigma(xi)`. For eigenvalues in descending absolute value,
the spectral radius of `wedge^k X` is their first `k` moduli multiplied.
It is at most its operator norm. Equation (2) therefore implies

```text
product_{j=1}^k |lambda_j(X)|
    <= exp[-xi(kT+sum_{j=1}^k d_j)].                    (3)
```

Weyl's 1949 proof obtains eigenvalue-product bounds by applying the ordinary
eigenvalue/singular-value inequality to exterior powers of a finite matrix
(pp. 409-410). Pólya's 1950 remark proves the weak-partial-sum lemma that
transfers those product bounds to sums of convex increasing functions (p. 50);
unlike equality-at-the-final-index majorization, it needs only the ordered
partial-sum bounds used here. See [Weyl (1949)](https://doi.org/10.1073/pnas.35.7.408)
and [Pólya (1950)](https://doi.org/10.1073/pnas.36.1.49). The following
argument also proves the particular convex step directly, so the comparison
does not rely on an unstated majorization convention.
Write `a_j=log|lambda_j|` and `b_j=-xi(T+d_j)`. All are negative, and (3)
is the inequality between their ordered partial sums. For any real `v`,

```text
sum_j (a_j-v)_+ = max_k [sum_{j=1}^k a_j-kv]
                <= max_k [sum_{j=1}^k b_j-kv]
                 = sum_j (b_j-v)_+.
```

For `h(x)=exp(x)/(1+exp(x))`,
`h''(x)=exp(x)(1-exp(x))/(1+exp(x))^3 >= 0` on `x<=0`, and
`h(x)=integral_{-infinity}^0 (x-v)_+ h''(v) dv` there. Integration proves
`sum h(a_j)<=sum h(b_j)`. Zero eigenvalues are handled by a limiting
negative logarithm and contribute zero. The scalar disk bound in the
previous note and the trace identity, including Jordan multiplicities,
give

```text
Re Tr[X(I+X)^(-1)]
  <= sum_j |lambda_j|/(1+|lambda_j|)
  <= sum_j 1/(exp[xi(T+d_j)]+1).
```

Differentiation at fixed exterior and the same integrable domination as
in the previous note yield the sharp bound

```text
E_loop'(T) <= (pi/6) sum_j 1/(T+d_j)^2.                 (4)
```

The identical unmixed returns `G_sigma=I` attain equality. This proves
that mixing cannot improve the length derivative within (1)-(2); it does
not prove a universal force sign or gravitational stability.

## Point dipole geometry

Take MMP's leading flat exterior with point monopoles separated by `d`.
Label a connecting line by `u=nu/2` in `(0,1]`. At its crossing of the
midplane the field-line equation gives

```text
rho_mid^2 + d^2/4 = d^2/(4u^2).
```

The two straight segments joining the monopoles to this crossing have
combined length `d/u`. Each actual half-line is at least as long as its
straight chord. Thus the exterior line length satisfies

```text
d_line(u) >= d/u.                                      (5)
```

This follows from the geometry in
[MMP section 3.1 and Appendix D](https://arxiv.org/html/1807.04726v3#A4),
without evaluating its arc-length integral.

## Distributed mixing in a specified chiral propagation model

There is also a local form of (5). In dimensionless prolate coordinates,
write the distances to the monopoles as `r_1,r_2`, with
`mu=(r_1+r_2)/d` and `eta=(r_1-r_2)/d`. Up to an orientation sign,

```text
z/d = mu eta/2,
rho/d = sqrt[(mu^2-1)(1-eta^2)]/2,
u(mu^2-eta^2)=mu(1-eta^2),   -1<=eta<=1.
```

For `0<u<1`, the positive solution obeys `1<=mu<=1/u`. Set
`Q=2u mu-1+eta^2`; implicit differentiation gives
`mu'=2 eta(u-mu)/Q`. The Euclidean arc-length rate `v=ds/deta` satisfies

```text
4v^2/d^2 = mu/u + 4 mu eta^2(mu-u)/Q^2,
4u^2 v^2/d^2 - 1
 = (1-u mu)(mu^2-1)(3mu^2-4u mu+1)/(mu^2-2u mu+1)^2 >=0.
```

The last factor is nonnegative because
`3mu^2-4u mu+1 >= (mu-1)(3mu-1)`. The axial `u=1` line is the continuous
straight-line limit. Thus `v(u,eta)>=d/(2u)` pointwise, not just after
integration.

Consider the specified free chiral action with time coefficient
`B(eta)=diag(v(u_j,eta))` and an arbitrary Hermitian mixing potential
`K(eta)`. Its equation is `i B partial_t psi+i partial_eta psi-K psi=0`.
The time norm uses the positive weight `B`, and `psi^dagger psi` is the
spatial probability flux. At frequency `omega=i xi`, propagation obeys

```text
partial_eta psi = [-xi B(eta)-i K(eta)] psi.             (5a)
```

The Hermitian part of the generator on `wedge^k` is `-xi` times the
additive exterior action of `B`. Its least eigenvalue before the minus
sign is the sum of the `k` smallest eigenvalues of `B`, at least
`(d/2) sum_{j=1}^k 1/u_j` with `u_j` in descending order. Differentiating
the squared norm and integrating from `eta=-1` to `eta=1` gives

```text
||wedge^k R(xi)|| <= exp[-xi d sum_{j=1}^k 1/u_j].       (5b)
```

Noncommuting mixing at different positions is allowed. A local dissipative
term with positive Hermitian part can only strengthen the norm bound.
The same proof works for the opposite chiral direction with its coordinate
reversed. Equations (5a)-(5b) establish the needed delay bound within this
local chiral model; no midpoint factorization is needed. They do not show
that a finite-flux four-dimensional Dirac projection has this action,
that nonlocal or chirality-changing terms are negligible, or that its
complete quantum energy has the assumed determinant and reference term.

For the factored network, (4) may retain the full assigned arc lengths.
For distributed mixing, (5b) establishes only the chord-delay comparison
`d_j=d/u_j`. It does not prove that arbitrary mixing preserves MMP's
stronger unmixed numerical ratio `2.3471`.

## Stationary ratio in the continuum comparison

In the continuum model the
channel measure is `q du`. Combining (4)-(5), or taking empirical delay
distributions converging to this measure, gives

```text
E'(T) <= -2A/T^3 - q pi/[24T^2] [1-4 J(T/d)],
J(t) = integral_0^1 [tu/(1+tu)]^2 du
     = 1+1/(1+t)-2 log(1+t)/t.                          (6)
```

The integrand increases strictly in `t>0`; consequently `J` rises from
zero to one. There is a unique positive `t_*` with `4J(t_*)=1`.
In particular `4J(2)<1`, so every stationary point of this model obeys
`T/d>t_*>2`. The displayed decimal `2.1489119899...` is numerical;
the ordinary theorem uses the uniquely defined root and the exact
inequality `4 log(3)>13/3`. That last inequality follows, for example,
from the positive series `log(3)=2 sum_{n>=0} (1/2)^(2n+1)/(2n+1)`:
its first three terms already exceed `13/12`.

MMP's full arc-length distribution gives the larger numerical limit
`2.3471178600` as `A/(qd)->0`. Our constant is a conservative chord
comparison, not a correction to their number or a sharp dipole constant.

## Correction and finite-distribution budgets

For `0<T<=d`, `1-4J(T/d)>=kappa`, where

```text
kappa = 8 log(2)-5 = 0.5451774445... >0.
E'(T) <= -2A/d^3 - q pi kappa/(24d^2).                  (7)
```

Positivity already follows from the first two terms of the positive
`log(2)=2 sum_{n>=0} (1/3)^(2n+1)/(2n+1)` series. An omitted contribution
preserves exclusion of stationary `T<=d` if its positive derivative is
uniformly smaller than the positive magnitude in (7). Unlike the
shortest-return comparison, the new term does not vanish relative to
`q/d^2` when `d` is large. No actual four-dimensional correction has yet
been bounded by it.

For discrete labels `u_j`, let the empirical CDF differ from the uniform
CDF by at most `delta`. Integration by parts for the increasing function
`g_t(u)=[tu/(1+tu)]^2` gives

```text
|q^(-1) sum_j g_t(u_j)-J(t)| <= delta [t/(1+t)]^2.       (8)
```

Thus (7) holds with `kappa-delta` when this is positive. This also
quantifies an error in the delay distribution; it is not a proof that
finite-flux quantum modes are localized on these lines. For the prescribed
midpoint discretization `u_j=(j+1/2)/q`, `delta=1/(2q)`.

A weaker but simple factor-two result has an explicit budget on `T<=2d`:

```text
kappa_2 = 4 log(3)-13/3 = 0.06111582134...,
E'(T) <= -A/(4d^3) - q pi[kappa_2-16 delta/9]/(96d^2). (9)
```

This requires the bracket to be positive. It holds for midpoint labels
with `q>=15`; this is a sufficient conservative condition, not a sharp
finite-channel threshold. The first five terms of the positive series for
`log(3)` give an exact rational certificate at `q=15`. Positive omitted
derivatives smaller than the magnitude in (9) preserve the factor-two
exclusion in the defined model.

## Physical limits and contribution decision

The flat point-dipole coordinate separation `d` approximates the exterior
observer comparison only in the leading MMP matching regime. Finite mouth
surfaces, lapse corrections and the throat/exterior partition must be
accounted for before (7) or (9) becomes a physical arrival-time theorem.
The Gaussian energy, the local chiral propagation law or sufficient delay
bound, the empirical spectral distribution and additional quantum sectors
are separate requirements.
MMP Appendix A already explains that gauging removes a fixed number of
massless modes; this does not supply a bound on all interaction derivatives
or justify replacing the complete state by a free determinant.

The new budget resolves the loss of margin caused by replacing every
dipole field line with the shortest path. It does not resolve the actual
four-dimensional error. Retain (4)-(9) as a candidate quantitative extension
for contribution review. MMP owns the dipole no-shortcut conclusion and
its stronger unmixed numerical ratio; Weyl and passive-scattering theory
own the comparison methods. A separate paper or large Lean campaign is
not selected on this milestone alone.

## Verification

The [check](../checks/dipole_delay_bound.py) and
[receipt](../evidence/dipole-delay-bound.json) verify eight exact identities,
75 pointwise matrix comparisons, 375 exterior-power product comparisons,
15 integrated energies/derivatives and 15 independent energy variations.
Signed permutation cycles also agree with the independent spectrum of
their unfolded circles. Nonnormal, lossy and resonant returns are included.
The matrix bound is sharp for the unmixed return.

A separate prolate-coordinate arc-length calculation reproduces the MMP
limit at three quadrature orders and checks the chord inequality. Eighteen
finite-quantile cases check (8). Numerical convergence is a check, not a
certified quadrature error; the comparison proof uses (5) and (8).
Six direct ODE propagators with noncommuting distributed mixing also obey
the exterior-power estimate; 168 local-speed cases check the geometric rate.
No Lean invocation, physical-stability result or accepted novelty claim.

Remaining to-do list: derive or bound the physical delay operator and total
energy corrections, complete the contribution verdict and select primary
plus fallback before certification and manuscript preparation.
