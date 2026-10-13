# Boundary wave packets and thermal spectral matching

S072, 2026-10-13. **The compact CAR first-moment bound does not by itself
transfer to the physical point boundary field.** Its radial term already
diverges in the free thermal reference while the boundary midpoint
correlator is finite. Subtracting the reference first moment does not
repair that inference: an explicit same-moment counterexample has a
smaller midpoint correlator, even with a shared ultraviolet tail.

The [fermionic compact theorem](fermionic-compact-transfer.md) remains
valid for its normalized canonical modes and stated Hamiltonian. This
calculation closes the proposed use of that bound alone as a continuum
boundary certificate. It fixes a spatially smeared operator, gives the
exact remaining terms, and identifies the thermal spectral quantity
needed instead. It does not reject the wormhole candidate.

## A canonical boundary packet with its dressing retained

At a finite spatial regulator in the projected gauge/fermion description,
transport each fermion to a common boundary frame by an equal-time
Wilson matrix `W(b,x)`. Let

```
c_f,m = integral dx f(x)^* [W(b,x) psi(x)]_m,
integral |f(x)|^2 dx=1.
```

The link coordinates commute with the undressed fermions, and their
matrix entries are commuting functions of the canonical link
coordinates. Using the equal-time field anticommutator and
`W(b,x) W(b,x)^dagger=1` gives
`{c_f,m,c_f,n^dagger}=delta_mn` and `{c_f,m,c_f,n}=0`.
This is a spatial, equal-time construction. Smearing an interacting
Heisenberg field in time does not automatically give the same CAR.
No independent extra link Hilbert space is appended to an already
dressed physical operator. Gravitational dressing and the moving
boundary still require their own matching.

There are `q=2j+1` packet components of charge one and spin `j`. Define
`Q_f=sum_m c_f,m^dagger c_f,m-q/2`, `J_f=c_f^dagger t c_f`.
The total asymptotic generators satisfy the same charge/spin Ward
commutators on this packet. Hence `Q_e=Q_tot-Q_f` and
`J_e=J_tot-J_f` commute with its CAR algebra. Applying the packet
Majorana channel `E_f`, normalized by `2q`, gives exactly

```
E_f(K)-K = lambda 1-(4/q)[x Q_tot Q_f+a J_tot dot J_f],
K=x Q_tot^2+a J_tot^2,       lambda=x+a j(j+1).           (1)
```

The centered packet charge need not separately have the full physical
charge lattice; its commuting complement carries the compensating
offset. Only the total charge and total spin are physical generators.
Equation (1) therefore does not replace the total charge by a selected
packet charge or discard the other radial modes.

For a frozen connection, let `R=beta dGamma(h)` be a quadratic radial
Hamiltonian, up to an irrelevant normal-ordering constant. Here `h`
is its one-particle operator, `P=|f><f| tensor 1_q`, and `dGamma` is
fermionic second quantization. The exact channel variation is

```
E_f(R)-R = (beta/q)[Tr(P h) 1-dGamma(P h+h P)].           (2)
```

To prove it, choose a one-particle basis starting with the packet.
Conjugation attenuates a centered packet bilinear by `1-2/q`, a
packet/complement bilinear by `1-1/q`, and leaves a complement bilinear
unchanged. The scalar is the packet trace. Thus (2) retains the
off-diagonal radial coherences even when the packet is not an energy
eigenmode. A dynamical Wilson dressing also has a commutator with the
gauge-field Hamiltonian; (2) is not that full dressed variation.

For example, with one canonical fermion and `b=exp(i g X)c`, an
oscillator Hamiltonian `H_B=(P_X^2+Omega^2 X^2)/2` gives

```
E_b(H_B)-H_B = g(2c^dagger c-1)P_X+g^2/2.               (3)
```

The mixed term and recoil constant follow by conjugating `P_X`.
This is the canonical-channel form of the recoil already retained in
[phase-influence-action.md](phase-influence-action.md). It explains why
including the Wilson line only in the transfer operator, while using
the undressed radial-energy variation, would be inconsistent.

## The exact thermal kernel

For a stationary Gibbs state `rho=exp(-beta H)/Z` and a canonical
fermion `c`, its positive anticommutator spectral measure is

```
dA_c(w)=sum_mn (p_n+p_m)|<m|c|n>|^2
                        delta(w-(E_n-E_m)) dw,
integral dA_c=1.
```

For the packet, average this measure over its `q` components. Its
normalized midpoint affinity and the channel first moment are

```
C_f = integral sech(beta w/2) dA_f(w),
S_f = 2 beta integral w tanh(beta w/2) dA_f(w)
    = 2 <E_f(beta H)-beta H>,
1-C_f <= S_f/8.                                         (4)
```

Indeed `sqrt(p_n p_m)=(p_n+p_m)/(2 cosh(beta w/2))`, and the two
transfer directions in the Majorana average give (4). These statements
are first finite-regulator identities; their continuum use requires the
indicated integrals. The established fermionic Lehmann kernel in
[Kaye, Chen and Parcollet, equations (1)--(2)](https://arxiv.org/html/2107.13094v2)
is `exp(-w tau)/(1+exp(-beta w))`; setting `tau=beta/2` gives the same
midpoint weight. No spectral reconstruction algorithm or energy cutoff
assumption from that paper is required here.

## Free localization test with the physical normalization

Use a free massless chiral field on the line at temperature `1/beta`,
with single-particle spectral density `1/(2pi)`. At equal time smear
the spatial Fourier modes by `exp(-epsilon |k|/2)`, `epsilon>0`.
For this free chiral Hamiltonian only, `w=k`, so

```
psi_epsilon(0)=integral dk/(2pi) exp(-epsilon |k|/2) a_k,
N_epsilon={psi_epsilon,psi_epsilon^dagger}=1/(pi epsilon),
c_epsilon=psi_epsilon/sqrt(N_epsilon),
dA_epsilon=(epsilon/2) exp(-epsilon |w|) dw.              (5)
```

The field is obtained from finite volume and bandwidth followed by
the local KMS limit; no trace over an unregulated infinite Dirac sea
is invoked. Put `d=2epsilon/beta` and

```
I(d)=integral_0^infinity exp(-d v) sech(v) dv
    = [digamma((d+3)/4)-digamma((d+1)/4)]/2.
```

Equation (4) now gives

```
C_epsilon=d I(d) ~ pi epsilon/beta ->0,
S_epsilon=(4/d) integral_0^infinity v exp(-v) tanh(v/d) dv
          ~2 beta/epsilon ->infinity.                  (6)
```

The physical opposite-side midpoint contraction, with the usual TFD
phase absorbed, is instead

```
b_epsilon=(N_epsilon/2) C_epsilon=I(d)/(pi beta)
          ->1/(2beta).                                 (7)
```

The constant initial clock `t=kappa U`, `kappa=2pi/beta`, therefore
gives `b_U=kappa/(4pi)`, agreeing with the earlier strip value
`b_t=1/(4pi)`. This is a normalization check at the preparation slice,
not a stationary-state prescription during the switch.

Thus the CAR deficit approaching one in (6) is not physical destruction
of the finite cross-boundary correlation (7). The free radial term in
(2) gives the large positive first moment exactly. A cutoff-independent
canonical-mode error `|Delta C_epsilon|<=delta` would only imply
`|Delta b_epsilon|<=delta/(2pi epsilon)`. To control the point correlator
this way, the error must instead scale at least as `epsilon` with a
controlled coefficient. The compact coefficient alone does not do so.

## Why subtracting first-moment bounds fails

The two positive inequalities in (4) cannot be subtracted to obtain
`C-C_ref >= -(S-S_ref)/8`. Here is an exact finite counterexample.
In units `beta=1`, set

```
F(E)=E tanh(E/2),     G(E)=sech(E/2),
r=F(1)/F(3),
A_ref=(1-r) delta_0+r delta_3,    A=delta_1.
```

Both measures have mass one and positive weights. Both give exactly
`S=2F(1)`. They can be realized by the same canonical packet
`c=sqrt(1-r)c_1+sqrt(r)c_2`, with free one-body Hamiltonians
`h_ref=diag(0,3)` and `h=diag(1,1)`, respectively. Nevertheless, putting
`s=sech(1/2)` and using the triple-angle identities gives

```
C_ref-C = 4(1-s)^2(s+2)/[3(4-s^2)] >0.                  (8)
```

The same construction works at energies `E,3E` for every `E>0`.
Symmetrizing either spectral measure in frequency changes neither `C`
nor `S`. This is a failure of the proposed inference, not a physical
RN spectral model or a wormhole no-go theorem.

The failure persists with a common ultraviolet tail. Mix each measure
with the same free packet measure (5), assigning the discrete part
weight `p_epsilon=nu epsilon<1`. The normalized measures still have
exactly equal first moments and the same high-frequency component.
After multiplication by `N_epsilon/2`, their midpoint difference is
the nonzero finite number `nu(C_ref-C)/(2pi)` for every epsilon.
Consequently even an exact cancellation of the subtracted first moment
does not bound the finite physical midpoint difference. Positivity
belongs to each full spectral measure, not to their signed difference.

## The usable spectral target and its limits

For the actual spatially regulated boundary field in a common operator
normalization, let `dnu_epsilon=N_epsilon dA_f`. This definition uses
its full interacting Lehmann measure. It does not assume that spatial
smearing multiplies the interacting energy spectrum by
`exp(-epsilon |w|)`. For the free comparator that factor is known from
(5). At equilibrium the exact reference-subtracted criterion is

```
Delta b_epsilon = integral [dnu_epsilon-dnu_ref,epsilon]
                                      /(2 cosh(beta w/2)),
|Delta b_epsilon| <= integral |dnu_epsilon-dnu_ref,epsilon|
                                      /(2 cosh(beta w/2)). (9)
```

A uniform bound on the right that survives removing the shared
regulator controls the finite midpoint contraction. For example, an
upper bound `delta/(2beta)` and existence of the matched limit imply
`b>=(1-delta)/(2beta)` in this reference convention. Low-frequency
spectral weight must be controlled; the thermal kernel suppresses high
frequencies exponentially, subject to a uniform tail bound. The prior
Gaussian charged normalization can only be inserted after proving the
required factorization and common finite operator matching.

Equation (9) concerns the prepared equilibrium state. Once switching
starts, an equilibrium spectral measure alone does not specify the
two-time lesser covariance or stress. The continuation must determine
that covariance with the same gauge constraints, initial dressing and
subtraction, then take the point-split stress limit. The exact mixed
current terms derived previously remain present.

Decision: retain the finite-mode compact theorem and packet identities,
but park the raw or reference-subtracted CAR first moment as a standalone
continuum capture certificate. Use the weighted midpoint spectral
quantity for the initial boundary matching and the driven two-time
covariance for stress. This is a narrower failure of an estimate, not
failure of the interacting candidate or the earlier free-JT opening.

Evidence: [checker](../checks/boundary_packet_spectral_matching.py),
[receipt](../evidence/boundary-packet-spectral-matching.json), and
[source record](../evidence/boundary-packet-source-review.json).
The checks compare direct CAR matrices, thermal spectral sums, radial
one-body formulas, digamma/quadrature values and the exact counterexample.
No novelty acceptance or Lean certification is claimed.

Remaining to-do list: bound the actual prepared boundary spectral
difference and derive the driven continuum covariance/stress on the
evolving geometry; include angular, massive and fluctuation effects;
establish preparation, sufficient general stability and finite-energy
transport.
