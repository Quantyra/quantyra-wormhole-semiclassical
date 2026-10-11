# Quantyra semiclassical wormhole research

Technical satellite for the [Quantyra Space research program](https://github.com/Quantyra/quantyra-space-planning/blob/main/docs/wormhole-research-program.md).

The current [finite-flux response audit](notes/finite-flux-response-audit.md)
derives the free lowest-Landau-level charge form factor and an axial metric
vertex on a warped spherical background. Including the reciprocal metric
source preserves the [conditional positive-energy theorem](notes/global-drift-energy.md)
for the prescribed drift model. Finite-flux corrections can precede the
retained backreaction order, so the physical approximation needs a joint
expansion. The known Berezin form factor and planar magnetic response are
credited. Full Casimir-state response, higher-level/contact matching, global
MMP completion and novelty remain open; no physical stability or shortcut
result is claimed.

The baseline is Sadhukhan, [arXiv:2610.09847v1](https://arxiv.org/abs/2610.09847v1).
The [mouth factorization](notes/mouth-energy-factorization.md) and
[baseline/throat audit](notes/axial-baseline-and-throat-energy.md) are preserved.
Classical mouth and corrected throat tensors and their master reductions
now reproduce independently, conditional on the stated quantum response.

Install the pinned Python dependencies and run:

```text
python -m pip install -r requirements.txt
python checks/mouth_energy.py --output evidence/mouth-energy.json
python checks/mouth_reduction.py --output evidence/mouth-reduction.json
python checks/mouth_geometry.py --output evidence/mouth-geometry.json
python checks/throat_energy.py --output evidence/throat-energy.json
python checks/throat_geometry.py --output evidence/throat-geometry.json
python checks/throat_reduction.py --output evidence/throat-reduction.json
python checks/general_drift_energy.py --output evidence/general-drift-energy.json
python checks/lll_projection.py --output evidence/lll-projection.json
python checks/lll_axial_vertex.py --output evidence/lll-axial-vertex.json
python checks/finite_flux_energy.py --output evidence/finite-flux-energy.json
```

The checks combine symbolic tensor/reduction identities with independent
quadrature of the mouth energy and its endpoint term. The original throat
squares agree only through first order. The later exact theorem uses new
variables and a defined response class; it does not bound quantum terms
outside that class. The checks use no Lean compiler. All Lean and Lake
development and acceptance must run on GCP, as specified by the planning
repository. No Lean certification exists for this satellite yet.

Public source papers are referenced, not redistributed. Original code and
notes use Apache-2.0. Review is open source and informal; specialist sign-off
is not a workflow gate. The research checkpoint is [checkpoint.json](checkpoint.json).

Remaining to-do list: finish the source and novelty audit,
match the Casimir-state and summed higher-level/contact response in a joint
finite-flux/backreaction expansion, resolve the physical domain,
then pursue and verify the selected primary and fallback questions.
