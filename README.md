# Quantyra semiclassical wormhole research

Technical satellite for the [Quantyra Space research program](https://github.com/Quantyra/quantyra-space-planning/blob/main/docs/wormhole-research-program.md).

The current [LLL return-phase diagnostic](notes/lll-return-holonomy.md)
shows that a general unitary return need not preserve the angularly
averaged negative support in the free equal-length channel model. A
common electromagnetic phase cannot cancel relative eigenphases. A
cyclic return has best support coefficient `-1/3+4/(3q^2)`, while a
sufficiently narrow phase spread preserves support. MMP already discusses
the qualitative long-cycle mechanism. The actual exterior return map,
interacting vacuum and consequential new physical result remain open.

The preserved [warped-circle winding bound](notes/warped-winding-bound.md)
controls both the energy and local radial null stress of every higher
angular fermion mode in a prescribed static spherical circle. It uses an
integrated throat barrier and permits a small angular mass elsewhere.
The antiperiodic winding null stress is non-positive in this model.
Angular mixing in the real MMP exterior and comparison of the winding
reference with the physical support state remain open. This conditional
result has cross-checks; novelty and physical MMP acceptance are unresolved.

The preserved [reference/Casimir matching](notes/casimir-reference-matching.md)
reproduces the complete first-order throat for arbitrary invariant vacuum
pressures and a specified leading LLL Casimir source. The invariant part
shifts the radii; the Casimir part supplies the opening. Direct contractions
verify all Einstein-Maxwell components, conservation and known limiting
cases. MMP and the later near-horizon EFT analysis already establish the
underlying separation, so this is retained as a baseline rather than a
standalone new-physics candidate. Corrected global transit time and total
quantum response remain unresolved.

The preserved [homogeneous vacuum calculation](notes/magnetic-product-vacuum.md)
derives both renormalized reference pressures from the full magnetic Dirac
spectrum. It reproduces the anomaly and supplies a nearby stationary
product branch for sufficiently small loop strength at each fixed charge.
Numerical radius variations independently check the pressure integrals.
This is a specified Euclidean reference state and subtraction scheme;
inhomogeneous Casimir matching, global completion and stability remain open.
The determinant construction has direct prior literature and is credited.

The preserved [trace-matching audit](notes/throat-trace-anomaly.md) identifies
a necessary correction to the quantum reference background. The full
massless charged Dirac anomaly is nonzero on the equal-radius magnetic
throat, while the displayed LLL lift is traceless. The required leading
radius-ratio correction is of order `g_c^2`, preceding the retained
`g_c^2/q` Casimir backreaction at large flux. This is an application of the
known anomaly, not a wormhole no-go or a complete corrected solution.

The preserved [Ward and sign audit](notes/casimir-cylinder-ward-completion.md)
restores the angular partners of the selected axial-metric/polar-gauge
sector, including exact background stress tadpoles. It proves that two
metric families with the same linear perturbation have opposite static
Casimir curvatures: the affine family is negative and the completed-shift
family positive, with all massive levels bounded. This isolated off-shell
matter sign is therefore insufficient to establish a physical instability.
The [Euclidean axial response](notes/casimir-cylinder-dynamical-response.md)
and its subthreshold analytic domain are preserved. Curved total response,
other perturbation sectors and global gravitational stability remain open.

The [finite-flux audit](notes/finite-flux-response-audit.md) and
[conditional positive-energy theorem](notes/global-drift-energy.md) are
preserved. The new response reproduces the leading spatial drift but
supplies a finite-flux correction outside that assumed constitutive class.
Known sphere spectra and geometric Casimir-response precedents are credited.
Novelty, physical MMP stability, global completion and shortcuts remain open.

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
python checks/casimir_cylinder_twist.py --output evidence/casimir-cylinder-twist.json
python checks/casimir_cylinder_response.py --output evidence/casimir-cylinder-response.json
python checks/casimir_cylinder_ward.py --output evidence/casimir-cylinder-ward.json
python checks/throat_trace_anomaly.py --output evidence/throat-trace-anomaly.json
python checks/magnetic_product_vacuum.py --output evidence/magnetic-product-vacuum.json
python checks/casimir_reference_matching.py --output evidence/casimir-reference-matching.json
python checks/warped_winding_bound.py --output evidence/warped-winding-bound.json
python checks/lll_return_holonomy.py --output evidence/lll-return-holonomy.json
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
derive or bound the physical exterior return phases, delays and leakage,
close the reference/interacting-state comparison,
derive the curved total response,
then pursue and verify the selected primary and fallback questions.
