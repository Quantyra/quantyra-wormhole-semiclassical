# Identical stationary pair: Dirac modulus, partner squares and phase

S139/S074, 2026-10-15. Ordinary argument before focused verification.
This note completes the partner-square and phase audit for the retained
free, vectorlike Dirac field on two identical stationary magnetic
exteriors. It does not compute the massive stress remainder, identify
the interacting physical state or close the stationary matching task.
The already frozen exterior notes and receipts remain unchanged.

## Problem, domains and the orientation of the two copies

Use the fixed pulled-back spin and monopole bundle identification of
[the coupled domain](coupled-fermion-domain.md), with no electric field,
identical clocks and identical metric jets on the two copies. Work first
at a finite smooth spherical outer cut, at zero constant source angle.
The Euclidean bolt is smooth with its regular thermal spin structure;
it is not another boundary. The reference sphere radius is b0.
Hermitian Euclidean Clifford matrices have square plus one. The bulk
expression D=gamma^mu nabla_mu is formally anti-Hermitian.

Let tau_i act on the folded pair, sigma_i on the orbit spinor and
gamma_A,S on the sphere spinor. Identity factors are implicit. A useful
orientation-aware representation in a normal orthonormal frame is

```
gamma_tau = -tau_3 tensor sigma_3 tensor 1,
gamma_n   = 1 tensor sigma_1 tensor 1,
gamma_A   = -tau_3 tensor sigma_2 tensor gamma_A,S.
```

Here n is the common pulled-back outer normal. These matrices satisfy
the four-dimensional positive Clifford algebra. With the convention
D_S,E=i D_S,herm, multiplying the covariant Dirac expression by
gamma_tau gives the same radial and angular coefficients as the
retained canonical folded Hamiltonian:

```
gamma_tau gamma_n = -i tau_3 tensor sigma_2,
gamma_tau[-tau_3 tensor sigma_2 D_S,E/b0]
 =1 tensor sigma_1 D_S,herm/b0.
```

The opposite orbit orientations are therefore included; no common
gamma matrix has been guessed independently of the Hamiltonian.
Canonical half-density rescalings add only the corresponding scalar
connections and do not alter the algebra below.

Set

```
chi=tau_3 gamma_n, Pi_+=(1+chi)/2, Pi_-=1-Pi_+.
Dom D:       Pi_- psi=0,
Dom Ddagger: Pi_+ phi=0.                            (1)
```

These complementary adjoint domains follow from the Green form
integral phi^dagger gamma_n psi: gamma_n preserves each chi eigenspace,
so the annihilator of the allowed trace is its orthogonal complement.
The first-order D on Dom D is not self-adjoint. The positive square is

```
L_+=Ddagger D:
 Pi_- psi=0, Pi_+(nabla_n+K/2)psi=0.                (2)
```

The partner L_-=D Ddagger exchanges the two projectors and has the same
K/2 Robin expression on its allowed trace. See
[the conformal-domain derivation](stationary-exterior-conformal-domain.md)
for the Gauss formula, inward heat convention and finite-cut Weyl map.
In particular, no angular-mass Robin term is introduced.

## Exact global partner intertwiner

Define two parallel Hermitian involutions

```
J=tau_1 tensor sigma_1 tensor 1,  T=tau_3 tensor 1 tensor 1.
```

Their elementary identities are

```
[J,gamma_mu]=[T,gamma_mu]=0,
J chi J=-chi, [T,chi]=0, {J,T}=0.                  (3)
```

For gamma_tau and gamma_A the two anticommuting Pauli factors in J
give a commuting product. Both operators consequently commute with
the spin connection. They also commute with the identical pulled-back
magnetic connection. Thus they commute with the differential expression
D globally, including the regular bolt. Only J exchanges the two
domains in (1). As closed operators,

```
Ddagger=-J D J,
L_-=J L_+ J.                                      (4)
```

This is more than equality of local heat coefficients. At each finite
cut the partner heat kernels are conjugate. For any scalar smearing f
which is identical on the copies and for positive heat time,

```
Tr[f exp(-t L_-)]=Tr[f exp(-t L_+)].               (5)
```

The same conclusion holds for local traced diagonal kernels under the
identified bundles. It remains valid along arbitrary identical scalar
Weyl paths, since the rescaling is scalar and J,T remain parallel in
the identified spinor bundle. It is not necessary to replace the actual
projector by a different bag condition to obtain this result.

## Why a modulus Weyl variation initially averages both squares

For a smooth finite-cut Weyl path g_t=exp(2t sigma)g_0, identify all
spinor L2 spaces unitarily. Four-dimensional Dirac covariance gives

```
D_t=exp(-t sigma/2) D_0 exp(-t sigma/2).
```

The first-order projector domains are fixed in this representation;
the induced square domains vary as prescribed by composition. Direct
differentiation gives

```
delta L_+=-(sigma L_++L_+ sigma)/2-Ddagger sigma D.
```

Initially for Re z large and with zero eigenspaces removed,
intertwining through D and cyclicity of trace-class products give

```
delta zeta_(L_+)(z)
 =z[Tr(sigma L_+^(-z))+Tr(sigma L_-^(-z))].         (6)
```

The traces in (6) mean positive spectral subspaces. Meromorphic
continuation supplies the regularized identity. For an invertible
problem, the Euclidean modulus action Gamma_E,mod=-logdet(L_+)/2
therefore has Weyl derivative equal to the average of the two smeared
a4 coefficients. In the presence of kernels the corresponding smeared
kernel projections must be subtracted from the a4 terms, with the
normalization of the zero-mode measure stated separately. A local
ultraviolet coefficient cannot encode that finite spectral datum.

Equation (5) makes the two coefficients, and also their matched kernel
terms, equal in the present problem. Thus either square suffices after
the partner audit. Using only one square before establishing this
identity would omit a step in the variation of the first-order modulus.
The explicit finite-cut local cocycle is derived separately in
[the boundary Weyl note](stationary-exterior-boundary-weyl.md).

## Self-adjoint representative, eta symmetry and the phase convention

On Dom D define Q=i J D. Equations (1)-(4) imply

```
Qdagger=Q, Q^2=L_+, T Q T=-Q.                      (7)
```

Self-adjointness here is a property of Q, not of D on its original
domain: J maps the complementary adjoint domain to Dom D. The unitary
T preserves Dom D and pairs every nonzero eigenvalue lambda of Q with
-lambda. Consequently its eta function vanishes where its defining
trace converges and by continuation:

```
eta_Q(z)=0.                                       (8)
```

In the exact product frame, the common domain separates and every
nonzero angular sector obeys

```
Q^2=L_+ >= n(n+q)/b0^2, n>=1.                     (9)
```

There are no zero crossings in those sectors along a smooth identical
finite-cut Weyl deformation; bounded positive Weyl factors preserve
invertibility as well. Hence there is no massive spectral-flow or
nonlocal eta-phase ambiguity in this paired family. Possible LLL
zero modes and their integration measures remain separate. Smooth
finite-cut conclusions alone do not prove convergence of an infinite-cut
determinant or its metric derivatives.

There is a necessary regularization qualification. Eta=0 does not, by
itself, make every complex spectral-cut definition of det Q positive.
A one-sided cut for a self-adjoint operator can retain a local phase
proportional to zeta_(Q^2)(0), in addition to its eta term. Likewise
replacing D by iJD involves a fixed unitary identification of its
codomain. These local branch/measure conventions must be fixed with
the common four-dimensional determinant prescription.

The vectorlike paired prescription is the positive singular-value
determinant, with phase fixed continuously from one reference pair.
Equivalently, the two flavor blocks are adjoint-conjugate up to a
fixed sign: J exchanges the flavor blocks and their allowed/adjoint
domains. In a finite paired regulator their determinant product is a
fixed regulator sign times the squared modulus of one block. Fixing
that sign at the reference gives a positive product. Equations (7)-(9)
then show that no varying massive eta phase must be added to this
prescription. If a complex zeta cut for Q is used instead, its local
cut phase must be included or translated to this convention explicitly;
it is not an additional unknown state-dependent stress.

An independent static anchor is the canonical antiperiodic thermal
determinant det(partial_tau+H), where H is the self-adjoint retained
Hamiltonian with its maximal isotropic domain. In a finite paired
regulator each real energy E contributes a positive factor
2cosh(beta E/2). Real vacuum subtraction preserves that phase choice.
The gamma_tau conversion to the covariant expression uses the same
local measure convention as the declared parity-even four-dimensional
prescription. An independently chosen topological or local phase is
not fixed by eta cancellation and must be declared separately; none is
silently identified with the computed modulus.

## Scope of the cancellation

A fixed constant source angle is transported by the constant flavor
rotation V only when it extends globally and preserves the complete
regularity and end domains. This additional hypothesis holds for the
identical paired HH disks with the bundle and frame identification used
here: its generator i tau_2 tensor sigma_1 commutes with every displayed
gamma matrix, is parallel, and preserves bolt regularity and the thermal
spin structure. Conjugating J,T and the complete domains by that V
transports the proof. No extra angle-dependent phase is generated in
this specified constant-angle family. The coupled-domain note's local
collar equivalence alone would not establish global equivalence on a
different interior join.

An evolving angle gives derivatives of V. Unequal metric or clock data,
unequal electric potentials or a nonintertwining bundle identification
can also invalidate (3)-(4). This note supplies no eta cancellation or
complete determinant statement in those cases. It does not determine
finite local counterterms, the Weyl-invariant massive stress, the outer
vacuum subtraction, the boundary charge or the physical interacting
state. Those obligations persist even though the identical stationary
pair no longer has an unexamined massive eta contribution.

Remaining to-do list: verify this ordinary algebra within the focused
S139 package; control the common-subtraction massive boundary remainder
and its horizon-weighted response; combine the full Weyl and boundary
variation with the exterior state and charge prescription; complete
stationary matching and the remaining S074-S079 gates.
