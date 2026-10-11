# Dipole transverse modes and the finite frequency closure test

2026-10-10. S052/S054 feasibility calculation. In the flat point-dipole
midplane, the exact transverse Dirac kernel has `q-1` normalizable modes,
not the `q` modes of a compact monopole sphere. There is no positive
transverse spectral gap. The mode count and gapless phenomenon are known
Aharonov-Casher results; the explicit profile below checks their consequence
for our proposed reduction. They prevent a globally gapped justification
of the channel model, without refuting MMP's leading approximation.

A separate conditional energy estimate shows that the
[dipole-delay bound](dipole-delay-bound.md) only needs stronger delay control
over a finite imaginary-frequency band. It can also reserve poorly
controlled modes at the common-front bound. This makes the next physical
test more precise; neither a full spatial reduction nor the required
quantum-energy error estimate has been proved.

## Exact magnetic midplane

Take point monopoles at `z=+-a`, separation `d=2a`, positive integer flux
`q`, and charge convention `e=1`. Choose the orientation giving positive
midplane field. In polar coordinates on the full plane, set

```text
R = sqrt(r^2+a^2),                u = a/R,
B(r) = q a/(r^2+a^2)^(3/2),
A = A_phi dtheta,               A_phi = q(1-a/R),
Phi = q log[(R+a)/(2a)],         Delta Phi = B.
```

The physical azimuthal vector potential is `A_phi/r`; its apparent origin
singularity is removable. Its components are smooth and bounded for fixed
`a>0`. The total plane flux is `2 pi q`. The self-adjoint massless Dirac
operator `D_a = sigma . (-i grad-A)` has common domain `H^1(R^2,C^2)`.
Its square is the Pauli operator. A bounded magnetic perturbation of the
free Dirac operator suffices here; a positive gap is not assumed.

Writing `w=x+iy`, zero modes of the supported spin are `h(w) exp(-Phi)`
with entire holomorphic `h`. At infinity the weight behaves as `r^(-q)`;
angular orthogonality and radial integrability restrict `h` to a polynomial
with degrees `0<=m<q-1`. The opposite-spin factor `exp(+Phi)` gives no
nonzero square-integrable entire solution. Thus, for `q>=2`, an orthonormal
basis is

```text
phi_m = C_m r^m exp(i m theta) [(R+a)/(2a)]^(-q) (1,0),
m=0,...,q-2;                   j=q-m-1=1,...,q-1.       (1)
```

At `q=1` the kernel is empty. The candidate `m=q-1` has a logarithmically
divergent norm. The missing plane mode is not a theorem about loss of a
global wormhole channel: the noncompact cross-section and its continuum
differ from the compact sphere. [Aharonov and Casher (1979)](https://doi.org/10.1103/PhysRevA.19.2461)
give the factorization; [Kenneth and Avron, equation (1.3) and section 2](https://arxiv.org/html/1403.2216v2)
state the strict integer-flux convention and normalizability argument.

## Localization and the weak tail

Set `t=(R-a)/(R+a)=(1-u)/(1+u)`. The angularly integrated probability is

```text
p_j(t) = t^(q-j-1) (1-t)^(2j-1) (1+t)/Z_j,
Z_j = Beta(q-j,2j)+Beta(q-j+1,2j)
    = [2q/(q+j)] Beta(q-j,2j),                         (2)
C_m^(-2) = 4 pi a^2 (2a)^(2m) Z_j.
```

This is a mixture of `Beta(q-j,2j)` and `Beta(q-j+1,2j)` with weights
`(q+j)/(2q)` and `(q-j)/(2q)`. Beta recurrences give

```text
E_j[u] = j/q,
E_j[1/u] = (2q^2-j)/(q(2j-1)) >= q/j,
E_j[(u-j/q)^2/u] = j(q^2-j^2)/(q^3(2j-1)) <= 1/q.     (3)
```

Consequently `Var_j(u)<=1/q`. This is absolute concentration, not a uniform
relative estimate for small `j`. The positive-radius moments obey

```text
E_j[r^p] < infinity iff p<2j, for p>=0,
E_j[r^2] = 4 a^2(q^2-j^2)/[(2j-1)(2j-2)], j>=2.      (4)
```

The outermost normalizable mode `j=1` has a logarithmically divergent
second moment. The `j=2` fourth moment likewise diverges. A finite-box
diagonalization can conceal these tails.

At the midplane the classical prolate arc-length rate is `R=a/u`.
Its quadratic-form compression to (1) is diagonal and at least `a q/j`
on mode `j`. For `j=1` use the form, since that mode is outside the operator
domain of multiplication by `R`. This is a single-slice comparison. It
does not identify the four-dimensional longitudinal generator with this
compression, nor prove propagation along a particular flux line.

## Absence of a positive spectral gap

Let `f(x)=(x-1)^2(2-x)^2` on `[1,2]` and zero elsewhere. Normalize
`psi_L = c L^(-1) f(r/L) exp(-i theta) (1,0)`. Its angular index `-1`
makes it orthogonal to every zero mode in (1). It lies in the Dirac form
domain, and the radial first-order factor is

```text
Q psi_L = c L^(-2) [f'(x)+(1+A_phi(Lx)) f(x)/x],
||D_a psi_L||^2 <= [24+2(q+1)^2]/L^2.                 (5)
```

Here `integral x f^2 dx=1/420`, `integral x f'^2 dx=1/35`, and
`0<=A_phi<=q`. The unit norm states have arbitrarily small Pauli energy
on the kernel's orthogonal complement. Therefore neither `D_a^2` nor
`|D_a|` has a positive gap above the kernel. Increasing fixed `q` cannot
remove this full-plane obstruction; choosing `L` sufficiently large always
defeats a proposed positive gap.

This excludes a uniform bounded reduced inverse and the standard gapped
adiabatic route. It does not exclude estimates for particular sources or
localized spectral windows. Gapless adiabatic control is established prior
work, including [Avron and Elgart](https://arxiv.org/html/math-ph/9805022v4).
Their section 5.1 extension assumes an operator bounded below, so it must
not be imported unchanged for a Dirac operator unbounded in both directions.

## A source-specific inverse diagnostic

For normalized modes, define the dilation source and opposite-spin function

```text
chi_j = (r partial_r+1) phi_j = (q u-j) phi_j,
eta_j = (0, (i/2) r exp(i theta) phi_j,upper).
```

The radial Dirac equations give `D_a eta_j=chi_j` pointwise. Also
`<phi_j,chi_j>=0`. For `j>=2`, this is a legitimate reduced solution with

```text
||eta_j||^2 = a^2(q^2-j^2)/[(2j-1)(2j-2)].            (6)
```

It lies in `H^1` by its origin and infinity behavior. Thus the lack of a
global gap does not make every reduced-source problem divergent.
For `j=1`, this particular solution is not in `L^2`. Cut it off by
`c_L(r)=1` below `L`, zero above `2L`, with `|c_L'|<=C_c/L`. Then

```text
Prob_1(r>L) <= 2 a^2(q^2-1)/L^2,
||c_L eta_1||^2 <= a^2(q^2-1) log(1+L/a),
||D_a(c_L eta_1)-chi_1||
    <= sqrt(2) a sqrt(q^2-1) (q+C_c)/L.              (7)
```

For the first inequality use (2) and `1-t_L<=2a/L`. For the norm use
`t^(q-1)<=1`, `1+t<=2` in its truncated integral and
`log[(sqrt(1+4(L/a)^2)+1)/2]<=log(1+L/a)`.
The last estimate follows from the commutator with the cutoff and
`|q u-1|<=q`; on the cutoff annulus `|eta_1|<=L|phi_1|`.
The cosine-squared transition used in the check has `C_c=pi/2`.

These identities supply a test case for a gapless method. A slow *temporal*
dilation of a planar Hamiltonian would be a different model from static
spatial transport through the dipole. Developing its adiabatic theorem
alone would not close the wormhole problem, so that extension is not the
selected next task. A physical field-aligned projection must account for
changing field direction, spin connection, continuum coupling and ports.

The threshold-resolvent route also needs an assumption check.
[Kovarik's Theorem 6.5](https://arxiv.org/html/2103.08262v3) and
[Baur's Assumption 1.1 and Theorem 2.4](https://arxiv.org/html/2409.17787v2)
give integer-flux expansions under magnetic decay faster than `r^(-7)`.
Our field decays as `q a r^(-3)`, so those displayed theorems do not apply
directly. Their general spectral setup has weaker hypotheses; that does
not waive the stronger expansion hypothesis. Both logarithmic threshold
terms and nonnormalizable resonant states must be considered. In (1),
formally taking `m=q-1` and `m=q` gives bounded non-L2 solutions behaving
as `r^(-1)` and a nonzero constant, respectively. These are additional
threshold diagnostics, not two extra propagating channels. Extending a
weighted-resolvent estimate to this exact tail is a possible method, not
an established ingredient. The cited expansion proofs and edition details
remain under review; their estimates are not imported here.

## Finite band and reserved mode energy bound

Retain the leading energy and its assumptions from the
[passive-return argument](passive-return-bound.md): `q` complex channels,
two chiralities, fixed-exterior variation, and

```text
E(T) = A/T^2 + q pi/(24T) + E_loop(T),  A>0.
```

For each chirality assume a causal passive return with front `d`. Suppose
the stronger exterior-power delay estimate of the dipole-delay note holds
only for `0<xi<=Lambda`, for a prescribed list `d_j>=d`. Arbitrary mixing
is allowed when it satisfies that estimate; designating a few modes as
uncontrolled does not by itself prove the estimate on the remaining ones.

The low-band trace comparison, together with the common-front bound on
the high band, gives

```text
E_loop'(T) <= (pi/6) sum_j 1/(T+d_j)^2 + U_Lambda(T),
U_Lambda(T) = (2q/pi) exp[-Lambda(T+d)]
                 * [Lambda/(T+d)+1/(T+d)^2].         (8)
```

Proof: integrate the sharper comparison up to `Lambda`, and use
`Re Tr[X(I+X)^(-1)] <= q/(exp[xi(T+d)]+1)` above it. Replace the latter
by `q exp[-xi(T+d)]` and integrate. Extending the positive low-band
comparison to infinity gives the displayed conservative bound. No
approximation to the high-frequency channel dynamics is required beyond
the assumed common-front passivity.

If `d_j=d/u_j` and the empirical CDF of `u_j` differs from uniform by
at most `delta`, the existing distribution estimate implies on `0<T<=d`

```text
E'(T) <= -2A/d^3 - (q/d^2) M(delta,x),
M(delta,x) = pi(8 log(2)-5-delta)/24
                 - (2/pi) exp(-x)(x+1), x=Lambda d. (9)
```

At `x=4`, `M(0,4)=0.01306307130...>0`. As a concrete target for the
mode audit, assign `u_j=j/q` for `j=2,...,q-1`, and two further entries
`u=1`. This reserves both the critical `j=1` mode and the missing plane
mode at the weaker common-front delay. Its CDF discrepancy is `2/q`.
Then `M(2/q,4)>0` for all integer `q>=21`; the check gives an exact
positive rational certificate using lower series bounds for `log(2)` and
`exp(4)`, with the standard rational lower bound `pi>3.14159`.
At `q=21` the margin is approximately `0.000596434 q/d^2`, in addition
to `2A/d^3`. More generally replacing `b` entries of the upper-endpoint
grid by `u=1` has discrepancy at most `(b+1)/q`.

The reservation is a proposed comparison spectrum, not a derivation of a
`q`-channel physical return from a `q-1` plane kernel. Coupling the reserved
sector to the others can invalidate a putative strong delay condition;
the exterior-power bound must hold for the full return. Equation (9)
gives a quantitative target for that test and for any positive omitted
energy derivative. It does not certify the physical corrections.

## Contribution and next decision

Known zero-mode counting, gapless adiabatic theory and passive-force methods
own the general mechanisms. The beta moments, explicit reduced source and
finite-band budget are checked calculations for our model-closure audit.
A final novelty verdict remains open. None independently establishes a
new wormhole construction, instability, or four-dimensional no-shortcut
theorem. The proposed physical obstruction remains a candidate under S051,
not an accepted primary/fallback selection.

The next decisive calculation is the static exterior Dirac return between
specified finite mouth ports, restricted initially to `xi d<=4`, with
field-aligned spin transport and continuum coupling. Establish or refute
the exterior-power comparison with the reserved-mode budget, and then
bound the determinant/reference and interacting-state corrections against
(9). Do not enlarge a purely planar adiabatic or formalization campaign
unless it supplies that missing estimate. A negative result must distinguish
failure of this comparison from failure of the wormhole approximation.

## Verification

[The check](../checks/dipole_transverse_modes.py) and
[receipt](../evidence/dipole-transverse-modes.json) verify 14 exact identities,
independent radial-coordinate integrals for 15 mode cases, nine critical
tail/cutoff cases, nine finite-band envelopes and the positive rational
`q=21` certificate. Radial quadratures do not use the beta-density formulas;
the critical-mode checks separately use the beta variable. The band
envelopes are pointwise bounds, not purported analytic scattering matrices.
Numerical quadrature is supplementary and has no certified full error
bound. No Lean invocation occurred.

Remaining to-do list: derive the physical finite-band return and total
energy correction bound, complete the decisive literature comparison and
target selection, then establish and verify the consequential result.
