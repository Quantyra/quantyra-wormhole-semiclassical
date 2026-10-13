# Real-time formation: entropy and throat hierarchy

S071, 2026-10-12. **A direct transfer of the SYK entropy-removal mechanism
does not follow from the presence of many magnetic Landau channels.**
The conditional scaling test below makes the obstacle quantitative.
The actual interacting boundary response remains to be derived.

[Maldacena and Milekhin](https://arxiv.org/html/1912.03276v2) simulate cooling
two coupled SYK systems into their wormhole-like ground state. Their
gravity comparison relates the ease of formation to removal of the
extremal entropy by a boundary-condition flow. They explicitly leave
four-dimensional magnetic-wormhole formation unanswered. We use this
as a comparator, not a new result or a laboratory spacetime construction.

## Conditional ordinary argument

Adopt the unit-charge Dirac-field conventions of
[MMP, equations 2.3-2.6 and section 7](https://arxiv.org/html/1807.04726v3).
For Nf such four-dimensional fields, the light two-dimensional Dirac
channel count is c=Nf q. For one extremal mouth,

    r_e² = pi q² G/g²,
    S0 = pi r_e²/G = pi² q²/g².

In the small-separation approximation underlying the section 7 length,

    ell = 16 r_e³/(G q Nf),
    ell/r_e = 16 pi q/(Nf g²).

Here ell is the throat redshift/time scale in MMP's metric conventions;
it is not the spatial proper length of the throat. The near-extremal
long-throat matching uses r_e/ell << 1. The channel condition q Nf >> 1
alone does not imply this separate hierarchy. Exterior distance and
other constraints also remain necessary.

Now **assume**, only to test the proposed transfer, that changing the
boundary condition produces an additive entropy change delta_s for each
real Majorana field, uniformly bounded as the hierarchy is taken. There
are 2c Majoranas, so for one disk/mouth the hypothesized change is

    DeltaS = 2 c delta_s,
    DeltaS/S0 = 2 delta_s Nf g²/(pi² q)
              = (32 delta_s/pi) r_e/ell.

The ratio is unchanged if both mouths are counted consistently. The
[symbolic checker](../checks/formation_entropy_transfer.py) and
[receipt](../evidence/formation-entropy-transfer.json) verify this identity.

For bounded delta_s, **parametrically improving the long-throat hierarchy
drives this entropy fraction to zero**. Increasing Nf cannot change the
fraction independently of that hierarchy. Requiring an order-one entropy
fraction with fixed nonzero delta_s instead keeps ell/r_e finite. This
does not quantify the accuracy at any particular finite value; no absolute
finite-parameter exclusion is inferred. We have neither evaluated delta_s
in the MMP theory nor assumed that all ways of shedding entropy obey this
boundary-condition accounting.

This is an algebraic specialization of the known few-fields versus
extremal-entropy distinction. It is not a new general obstruction, a rate
bound, or a statement that radiation cannot carry away entropy over time.
One must distinguish ordinary bath entropy flux from a change of the
topological entropy coefficient used in the proposed gravity analogy.

## Next physical test

The inspected real-time source's gravity model uses a relevant fermion
bilinear with individual dimension Delta<1/2. MMP's lowest Landau channels
are massless two-dimensional fermions. Before inserting them into that
flow, derive their boundary dimensions and actual coupling through the
four-dimensional exterior, including the charged-field and gauge
constraints. Determine whether the proposed interaction is relevant,
marginal or irrelevant and whether any infrared enhancement invalidates
the bounded, additive delta_s premise. A field-count comparison cannot
answer that question.

If the actual theory fails this gate, record a null for this transfer
and proceed to a distinct state-preparation route. If it survives, specify
an initial density matrix, bath coupling, energy and entropy flux, and
the real-time equations across the transition before running a large
simulation. The SYK result alone is insufficient motivation to reproduce
its numerics when the four-dimensional transfer is the hard question.

Source reading scopes and hashes are recorded in the
[supplement](../evidence/formation-transfer-source-review.json).
No manuscript, publication, paid computation or Lean certification is
selected. The stable-wormhole objective remains active and incomplete.

Remaining to-do list: derive the actual boundary interaction and entropy
response, assess whether a controlled transfer survives, then establish
formation dynamics and finite-disturbance traversability for a survivor.
