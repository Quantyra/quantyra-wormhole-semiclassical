# Angular transfer matching and its finite-order limit

S072, 2026-10-13. **The angular transfer weights and a separated-insertion
transport bound are now matched to the LLL representation. The bound is
small at fixed order but does not control the complete physical switch.**
This reuses the exact electric/angular/coset stress decomposition already
established in [the closed S053 work](angular-current-coset.md), retained
in commit `fe68067`. That local algebra is not a new result of S072 and
its reuse does not reopen S053. The sections below restate its conventions
before deriving the transfer weights and transport bound.

This is an algebraic matching result and a decision about the next
calculation. It does **not** supply the fully projected neutral transfer
correlators, nonlinear gauged boundary action, compact state weights or
renormalized interacting stress. The [electric neutral kernels](neutral-transfer-current.md)
remain valid within their stated scope; they still contain the angular
sector identified here.

## Current algebra, including its normalization

Use the same unfolded free chiral reference, now on the complex plane
for local operator products. There are `q>=2` complex fermions in spin
`j=(q-1)/2`. Normalize

```
psi_m(z) psi_n^dagger(w) ~ delta_mn/(z-w),
[t^a,t^b]=i epsilon_abc t^c,   sum_a t^a t^a=C_j I,
C_j=j(j+1),    tr(t^a t^b)=S_J delta_ab,
S_J=q(q^2-1)/12.

J_c=:psi^dagger psi:,    J^a=:psi^dagger t^a psi:.
```

Wick contraction gives

```
J_c(z) J_c(w) ~ q/(z-w)^2,
J_c(z) J^a(w) ~ 0,
J^a(z) J^b(w) ~ (kappa/2) delta_ab/(z-w)^2
                         + i epsilon_abc J^c(w)/(z-w),
kappa=2 S_J=q(q^2-1)/6.                              (1)
```

For `N_f` identical spin multiplets, replace `q` in the electric level
by `n=N_f q` and `kappa` by `N_f q(q^2-1)/6`. All quantitative bounds
below select `N_f=1`. The angular coefficient in
[causal-compact-response.md](causal-compact-response.md) is therefore
exactly `lambda_J=S_J/pi=kappa/(2pi)` in that action's conventions.
Its coefficient is the current anomaly level, not `kappa+2`.

The non-Abelian term is already visible without any continuum fit:

```
<0|J_1^x J_0^y J_-1^z|0>=i S_J.                      (2)
```

This follows from `[J_1^x,J_0^y]=i J_1^z`, the vacuum conditions and
`<J_1^z J_-1^z>=S_J`. Three centered independent Gaussian currents
would give zero. Linearization in each axis can reproduce the two-point
function and the previously derived positive Gaussian action; it cannot
be promoted to the complete angular current theory on that basis.

The underlying affine embedding and stress subtraction are established
methods; see [Delmastro, Gomis and Yu, section 4.2 and appendix A.1](https://arxiv.org/html/2108.02202v2).
Equation (1) is rederived here in complex-fermion and standard spin
normalizations to avoid a real-fermion or generator factor of two.
No claim about the infrared gap of the curved RN problem is imported
from flat two-dimensional QCD.

## Exact stress decomposition

Let `T_F` be the free complex-fermion stress, with central charge `q`.
Define the normal-ordered Sugawara tensors

```
T_c = :J_c J_c:/(2q),
T_J = sum_a :J^a J^a:/(kappa+2),
T_C = T_F-T_c-T_J.                                    (3)
```

The current algebra gives a Virasoro tensor `T_J` of central charge
`c_J=3 kappa/(kappa+2)`. Both `T_F` and `T_J` make `J^a` a primary of
weight one; their singular products with `J^a` cancel in the
difference. The same statement holds for `J_c` and `T_c`, and the two
current algebras commute. Thus `T_C` commutes with both current
algebras and is a Virasoro tensor with

```
c_C = q-1-3 kappa/(kappa+2)
    = q-4+6/(kappa+2).                                (4)
```

The three Virasoro algebras commute. The coset is the commutant of
`SU(2)_kappa` inside the electrically neutral `SU(q)_1` sector. The
construction uses the positive inner product of the free fermion
Hilbert space; its branching multiplicity spaces carry the unitary
coset representations. It is not an assertion of an unrestricted
tensor product of three independent physical Hilbert spaces.

One useful direct diagnostic is the level-two vacuum Gram matrix.
With `|T_A>=L_-2^A|0>`,

```
<T_A|T_B>=delta_AB (c_A/2),
(c_c,c_J,c_C)=(1,3 kappa/(kappa+2),q-1-c_J).          (5)
```

The inherited S053 checker built these states from canonical fermion
creation and annihilation operators, independently of Sugawara
commutator formulas:

```
|T_c>=J_c,-1^2|0>/(2q),
|T_J>=sum_a J_-1^a J_-1^a|0>/(kappa+2),
|T_F>=sum_(p,m) (p+1)c_(p+2,m)^dagger c_(p,m)|0>.
```

Only single-particle modes `p=+-1/2,+-3/2` can contribute through
excitation energy two. The preserved prior receipt verifies (5) and
annihilation of `|T_C>` by both current algebras' modes `n=0,1,2`.
The present checker reuses those CAR operators to verify (2); it does
not replace or reissue the historical receipt.

For `q=2,3`, `c_C=0`; the coset stress state vanishes. At `q=4`,
`c_C=1/2`. These small cases sharply distinguish the exact result from
subtracting three unit-central-charge Gaussian bosons, which would
give a negative central charge at `q=2,3`. We do not infer the entire
coset operator algebra or its sewing data from this table alone.

For the candidate `q=k^4`, `k>=2000`,

```
0 < (c_c+c_J)/q < 4/q <= 2.5e-13.                    (6)
```

This is the fraction of the *reference stress central charge* assigned
to the retained current sectors. It is not a bound on their energy in
an interacting state, and does not authorize replacing `N` by `2c_C`
in the earlier JT switching/capture theorem.

## Transfer weights and branching constraints

A fermion is an affine primary of charge one and spin `j`: its OPE
with either current has only the expected simple pole. Acting on its
lowest state with (3) gives

```
h_c=1/(2q),
h_J=C_j/(kappa+2),
h_C=1/2-1/(2q)-C_j/(kappa+2).                         (7)
```

The weights add to `1/2`, as required. Each endpoint transfer has two
fermion insertions, so the charge and angular parts together carry

```
d_retained=1/q+2 C_j/(kappa+2) < 4/q,
d_coset=1-d_retained > 1-4/q.                         (8)
```

For this irreducible single-species embedding, `h_C=c_C/(2q)` and
`c_C=(q-3)(q-2)(q-1)(q+2)/(q^3-q+12)`. This also verifies the
nonnegative coset weight and the two small-`q` zeros directly.

The strict inequality follows from
`2 C_j/(kappa+2)=3(q^2-1)/(q(q^2-1)+12)<3/q`.
At `q=2,3`, `h_C=0`; at `q=4`, `h_C=1/16`. These are local reference
weights. They do not fix the renormalized normal-ordering constants,
the physical boundary coupling or the RN-to-strip clock map.

The representation decomposition has the form

```
H_F = direct_sum_(Q,ell) H_(U(1),Q)
                         tensor H_(SU(2)_kappa,ell)
                         tensor M_(Q,ell),            (9)
```

with only the branches present in the free-fermion character and its
chosen spin structure. `M_(Q,ell)` is the coset multiplicity space.
The necessary center selection rule is
`2ell=(q-1)Q mod 2`, consistent with
[quantum-junction-state.md](quantum-junction-state.md); it is not a
complete list of branches. In particular, zero coset central charge
does not erase all representation/sector information.

Correlation functions are sewn sums of angular and coset conformal
blocks in these branches. One cannot divide a summed fermion
four-point function by a single averaged SU(2) character and assume
the quotient is a physical neutral correlator. Nor is the WZW spin
primary an unrenormalized unitary rotor matrix. The previous electric
normal-ordering issue has a non-Abelian counterpart.
[Delmastro and Gomis, introduction and section 2](https://arxiv.org/html/2211.09036v1)
describe the current/coset operator and branching method used for
this distinction. Their proposed infrared theory and dynamical
condensates are not premises of this calculation.

## A bound that is genuinely uniform in spin, but only at fixed order

Consider `p` angular primary insertions, each in spin `j` or its dual,
on a smooth path of distinct complex positions `z_i(s)`. Use a fixed
orthonormal tensor basis and a continuous analytic branch. The
Knizhnik-Zamolodchikov connection in the conventions (1) is

```
dF/ds = A(s) F,
A(s) = [2/(kappa+2)] sum_(i<l)
                  (t_i dot t_l) d/ds log(z_i-z_l).     (10)
```

For example, in the two-point singlet, `t_1 dot t_2=-C_j`, so (10)
gives the power `(z_1-z_2)^(-2h_J)`. This checks the factor two against
(7). [Bouwknegt, Ludwig and Schoutens, sections 2-3](https://arxiv.org/html/hep-th/9412199v1)
give the affine null-state derivation and a four-point example; their
generator contraction equals twice the usual Cartesian spin dot
product used here.

The eigenvalues of `t_i dot t_l` are
`[ell(ell+1)-2C_j]/2`, `ell=0,...,2j`. Hence its operator norm is
exactly `C_j`, also for the dual representation. Define

```
L_path = sum_(i<l) integral |d log(z_i-z_l)|,
b_q = 2 C_j/(kappa+2) < 3/q.
```

If `U(s_1,s_0)` is parallel transport with initial value the identity,
the convergent path-ordered series and the operator norm inequality
give

```
||U|| <= exp(b_q L_path),
||U-I|| <= exp(b_q L_path)-1.                         (11)
```

This is valid on the ambient invariant tensor space and therefore
bounds transport of any actual physical block vector. It does not
choose the physical fusion channels, initial OPE coefficients or
coset sewing. Normalized ratios of correlators close to a zero need
not satisfy the same relative-error bound. Coincident-point limits
are explicitly excluded: `L_path` can diverge there.

For four positions `(0,x,1,infinity)`, after the usual definition of
the insertion at infinity, use `x:1/5 -> 4/5`. Then

```
dF/dx = [2/(kappa+2)] [t_1 dot t_2/x
                          +t_2 dot t_3/(x-1)] F,
L_path=log(16).
```

At `q=2000^4`, (11) is less than `5.199e-13`. The checker integrates
this matrix equation in a recoupled singlet basis and independently
verifies its pair matrices with four-spin Clebsch-Gordan tensors.
The two pair matrices fail to commute, and the resulting transport
has nonzero off-diagonal entries. A scalar coherence factor therefore
does not describe even this fixed-order angular transport.

## Why this does not yet justify the physical switch

A term with `r` transfer vertices has `p=2r` endpoint primaries. If
each pair has log-path length at most `L_*`, (11) only gives

```
b_q L_path <= (3/q) r(2r-1) L_*.                     (12)
```

At fixed order and fixed separation this tends to zero as `q` grows.
Uniform control by this estimate would instead require
`r^2 L_*/q << 1`. Moreover, summing a majorant of the form
`A^r exp(const*r^2/q)/r!` over all `r` diverges for any `A>0` and
positive constant. The finite-order estimate alone cannot be summed
into an all-order error theorem.

One can cut the series and use a separate bounded-Hamiltonian tail
estimate at a finite regulator. That still requires keeping enough
orders relative to `A=integral ||H_int(t)|| dt`. For a channel sum
with a fixed per-channel vertex bound, the elementary norm estimate
has `A=O(q integral |h| dt)`. Even the favorable illustrative scaling
`integral |h| dt=O(1/k)` gives `A=O(k^3)`, whereas `sqrt(q)=k^2`.
For the diagnostic choice `r=q/k` and unit pair log length, the exact
exponent `b_q r(2r-1)` is approximately `2.4e7` at `k=2000`.

This is a failure of that sufficient *absolute-value estimate*, not
evidence that the physical correction is large. It is not a measured
event count; the physical `h` normalization is still unmatched.
Disconnected contributions, connected cumulants, unitarity and Ward
identities can produce cancellations that this estimate discards.
Likewise, the bound on (6) does not bound energy released into
current sectors or the mixed bath-current term in
[phase-influence-action.md](phase-influence-action.md).

The resulting decision is to retain the exact current/coset branching
and use a connected influence or constrained Hamiltonian resummation
for the switch. Three independent angular Gaussian noises appended
to the electrically neutral fermion kernels are not an accepted
completion. The quadratic RN response remains a valid comparison
sector; its nonlinear error and transfer vertices require matching.

## Reproducible evidence and scope

Run `checks/angular_transfer_matching.py`. Its receipt is
`evidence/angular-transfer-matching.json`; source versions, scopes and
hashes are in `evidence/angular-transfer-source-review.json`.

Six direct CAR cases `q=2,3,4,5,6,8` verify the non-Gaussian Ward
coefficient using the preserved implementation. Three independent
tensor-product cases verify angular recoupling. Five integrated KZ
transports `q=3,4,5,8,12` obey (11); a tighter solve checks numerical
error. Exact rational arithmetic checks the weights and small cases,
and high precision evaluates the selected large-`q` comparison.
The receipt pins the inherited checker, its geometry dependency and
the historical coset receipt. These computations test the ordinary
derivation, not its physical four-dimensional transfer. No Lean
certification or novelty claim.

Remaining to-do list: construct the correctly branched coupled influence
and resum the switch; match its clock, vertices and continuum stress;
control evolving geometry and omitted modes; establish pair preparation,
stability and finite-energy transport.
