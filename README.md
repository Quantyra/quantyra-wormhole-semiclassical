# Quantyra semiclassical wormhole research

Technical satellite for the [Quantyra Space research program](https://github.com/Quantyra/quantyra-space-planning/blob/main/docs/wormhole-research-program.md).

The first calculation is an exact energy factorization of the three-channel
mouth operator in Sadhukhan, [arXiv:2610.09847v1](https://arxiv.org/abs/2610.09847v1),
equations (5.20)-(5.22). The [derivation](notes/mouth-energy-factorization.md)
gives a conditional mathematical positivity result, including the endpoint
term. The [baseline and throat audit](notes/axial-baseline-and-throat-energy.md)
now reproduces the classical mouth equations independently and identifies
a first-order energy structure for the corrected throat. Quantum closure,
global matching, remainder control and novelty remain under review. It does not
establish stability of the full wormhole or tolerance of finite disturbances.

Install the pinned Python dependencies and run:

```text
python -m pip install -r requirements.txt
python checks/mouth_energy.py --output evidence/mouth-energy.json
python checks/mouth_reduction.py --output evidence/mouth-reduction.json
python checks/mouth_geometry.py --output evidence/mouth-geometry.json
python checks/throat_energy.py --output evidence/throat-energy.json
```

The checks combine symbolic tensor/reduction identities with independent
quadrature of the mouth energy and its endpoint term. The throat squares
agree only through first order; their second-order cost is retained as an
unresolved stability issue. The checks use no Lean compiler. All Lean and Lake
development and acceptance must run on GCP, as specified by the planning
repository. No Lean certification exists for this satellite yet.

Public source papers are referenced, not redistributed. Original code and
notes use Apache-2.0. Review is open source and informal; specialist sign-off
is not a workflow gate. The research checkpoint is [checkpoint.json](checkpoint.json).

Remaining to-do list: finish the source and novelty audit, reproduce and audit
the corrected throat, resolve global matching and quantum-response limits,
then pursue and verify the selected primary and fallback questions.
