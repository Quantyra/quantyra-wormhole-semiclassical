# Quantyra semiclassical wormhole research

Current continuation, 2026-10-11: the [return-dependent interaction stress](notes/winding-interaction-stress.md)
controls the first photon vacuum diagram's return contribution and its
first two static clock derivatives. An exact endpoint formula retains
unequal returns and angular mixing. The renormalized reference,
higher-order and dynamical stress, and physical matching remain open.
[Planning audit](https://github.com/Quantyra/quantyra-space-planning/blob/main/docs/S053-winding-interaction-stress-audit.md).

Remaining to-do list: finish reference and higher-order static/dynamical
stress, match the coupled geometry and state, finish novelty and the
consequential result, then verify, GCP-certify and prepare the manuscript.

Technical satellite for the [Quantyra Space research program](https://github.com/Quantyra/quantyra-space-planning/blob/main/docs/wormhole-research-program.md).

The [global trial and spherical traces](notes/global-mouth-trial.md)
extend the corrected dipole field through both mouth boundaries. End
cutoff derivatives lie inside the removed balls, so the exterior weighted
residual tends to zero at fixed geometry. Exact spin/measure matching
reveals a persistent wrong-spin component of order `(rho/d)^2` and a
smaller profile mismatch. Including both gives a conditional trimmed
large-flux return bound on the selected frequency window. A finite flux
threshold, the boundary realization, exceptional-channel energy and
curved-throat/state matching remain unaccepted. No novelty claim follows.

The [interior calculation](notes/exterior-trial-residual.md) supplies the
residual cancellation and weighted current criterion used here. Its
historical global-matching gap is addressed by the new trial, while the
physical energy and complete return problem remain open.

The [finite mouth-collar calculation](notes/mouth-collar-scattering.md)
controls the radial corrector derivative and converts its residual into a
prepared-input scattering bound for the actual flat dipole annulus. The
ordinary proof retains every angular excitation in the selected sectors.
For `3<=q<=10^12`, `rho_1/d<=1/12500` and `0.05<=xi*d<=5`, the collar
error is below `8e-6`. Coupling this local map to the actual exterior load
still requires higher-mode feedback and curved throat/state control;
the full two-mouth return and energy result remain unresolved.

The [finite spherical mouth calculation](notes/mouth-mode-matching.md)
identifies exact deformed angular zero modes in the leading flat dipole
field. A compact selected-sector gap bounds their radial source correction
by `27/sqrt(160) (rho/d)^2` for `q>=3` and `rho/d<=0.1`.
Using isolated-monopole profiles at fixed `rho/d` is not uniformly accurate
as flux grows. Use the deformed modes for the next throat/port matching
step. The [collar result](notes/mouth-collar-scattering.md) controls the correction
derivative; full physical return and energy remain unresolved; no novel physical theorem or Lean certificate follows.

The [exterior Dirac audit](notes/exterior-dirac-window.md) derives the
exact two-dimensional operator in each angular sector of the leading flat
dipole exterior, including spin/measure and curvature terms. Its mass-zero
contours have the midpoint flux labels. A conditional energy bound now
tolerates a return error of `0.001 exp(-xi*D)` on `0.05<=xi*D<=4`, with
delay-CDF error at most 0.05. The remaining no-shortcut margin is about
`0.00357637 q/D^2`, plus the classical contribution. The actual port return
and its error remain unproved; known curved-edge propagation methods are
credited and under comparison.

Here `D` is the common exterior front time between the selected ports.

The preserved [transverse-mode audit](notes/dipole-transverse-modes.md) finds the
known gapless Aharonov-Casher structure in the exact dipole midplane:
`q-1` normalizable modes and a weak outer tail. A globally gapped channel
reduction is unavailable. A conditional finite-frequency estimate retains
the no-shortcut margin with stronger delay control only for `xi*d<=4`,
reserving two modes at the common-front bound for `q>=21`. The physical
exterior Dirac return and total energy corrections remain unresolved;
that finite-band return is the next test. No new physical no-go or novelty
acceptance is claimed.

The preserved [dipole-delay comparison](notes/dipole-delay-bound.md) gives a
conditional stationary transit ratio above 2.1489 in the point-dipole
continuum, and a correction budget
of order `q/d^2`. It allows distributed Hermitian mixing in a specified
local chiral propagation model. Matrix, ODE and spectral checks agree; an
independent arc-length calculation reproduces MMP's known unmixed ratio
2.3471. The physical four-dimensional projection, energy corrections and
novelty remain open. The earlier [passive-return bound](notes/passive-return-bound.md)
and its weaker common-front assumption are preserved.

The preserved [LLL return-phase diagnostic](notes/lll-return-holonomy.md)
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
python checks/passive_return_bound.py --output evidence/passive-return-bound.json
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
justify the physical exterior return and Gaussian energy,
bound omitted derivatives against the passive-return margin,
derive the curved total response,
then pursue and verify the selected primary and fallback questions.
