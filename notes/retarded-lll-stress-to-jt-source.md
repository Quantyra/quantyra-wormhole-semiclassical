# Fixed-throat retarded LLL stress as an (l=0) JT source

2026-10-12. This converts a prescribed conformal state map into the
conserved chiral stress that sources the leading spherical JT constraint.
It closes the local fixed-throat source map for free massless LLL modes.
It does not determine the map from the physical curved mouths, solve the
four-dimensional state problem, or establish stability.

## Stress from a characteristic state map

Use global throat conformal null coordinates `u=tau-x`, `v=tau+x`, with
`x=atan(rho)` and `ds2^2=r_e^2 sec^2(x)(-d tau^2+dx^2)`. For `q` complex
massless LLL channels the fixed-domain state difference has central charge
`c=q`; its direct area variation vanishes after the canonical spinor
rescaling. Let `t_+^0(u)` and `t_-^0(v)` be the chiral stress functions of
the chosen reference state, including its Casimir constant in the MMP
normalization. For the state difference in `lll-state-difference-source.md`,
the null-coordinate reference values are `t_+^0=t_-^0=-q/(16 pi)` in the
displayed unit-throat coordinates. An orientation-preserving characteristic state map
`R_+(u)` or `R_-(v)` pulls back the chiral functions as

```text
t_+(u) = R_+'(u)^2 t_+^0(R_+(u)) - q/(24 pi) Sch(R_+,u),
t_-(v) = R_-'(v)^2 t_-^0(R_-(v)) - q/(24 pi) Sch(R_-,v).       (1)
```

The coefficient and sign use the same point-splitting convention as the
transparent-clock stress formula. For maps equal to the identity in the
static past, (1) preserves that past reference exactly. The state-dependent
change on the fixed throat is
`delta T_uu=t_+-t_+^0(u)`, `delta T_vv=t_--t_-^0(v)`, and
`delta T_uv=0`. It is traceless and conserved on the fixed conformal metric:
the two nonzero components depend only on their respective null coordinate.
This is precisely the source condition required by the linearized Bianchi
identity for the leading JT operator.

At first order, write `R_+(u)=u+eta_+(u)` and
`R_-(v)=v+eta_-(v)`. Then

```text
delta T_uu = eta_+ (t_+^0)' + 2 eta_+' t_+^0 - q/(24 pi) eta_+''',
delta T_vv = eta_- (t_-^0)' + 2 eta_-' t_-^0 - q/(24 pi) eta_-'''. (2)
```

For constant reference Casimir functions, the first term vanishes. Equation
(2) is a retarded *source map* once the causal characteristic maps are
known; it is not yet a retarded kernel in physical mouth time. In particular,
the finite-return relation in `transparent-clock-flux.md` determines
boundary histories in that reduced model, but its map has not been shown to
equal the characteristic map induced by the actual curved MMP collars and
state.

## Coupling to the monopole constraint

In the leading near-product normalization of
`l0-throat-constraint-audit.md`, feed the state difference into
`L_ab[delta Phi]=delta S_ab` with

```text
delta S_ab = (2G/r_e^2) delta T_ab,
L_ab[Phi]=g_ab Box Phi-nabla_a nabla_b Phi-g_ab Phi.          (3)
```

The two chiral source functions from (1) determine the null projections of
the constraint. Their conservation is exact, so (3) is locally compatible.
The trace equation remains `Box delta Phi-2 delta Phi=0` for this traceless
state difference, while the null constraints fix the allowed solution from
the chiral source and integration data. This separates an arbitrary clock
history from the physical question: it must first induce `R_+,R_-` through
the two-mouth propagation and then solve all components of (3). The three
source-free embedding modes remain integration data; they are not a
frequency-dependent radiative monopole.

The static reference is the conserved MMP LLL source recorded in
`lll-state-difference-source.md`; equation (1) reduces to it when both maps
are identities. This is a reproduction of the conformal stress-tensor
transformation law applied to the known state, not a new quantum response
result. The absolute anomaly/Jacobian, massive angular tower, electric
current, moving-boundary contact and corrected background are not included.

Remaining to-do list: derive the actual two-mouth characteristic maps and
their boundary/state variations; solve all JT components with the resulting
retarded source and matched asymptotic data; then include the renormalized
massive tower and absolute anomaly before making a stability claim.
