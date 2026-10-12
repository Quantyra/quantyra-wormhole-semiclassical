# Claim to evidence

Endpoint: the symmetric MMP proper-mouth clock is not an independent boundary mode once both asymptotic clocks are normalized. Verdict: already known. No shortcut, formation, or engineering claim.

| Claim | Evidence |
| --- | --- |
| Retarded stress (1), RN reduction, trace anomaly | `checks/retarded_dirac_spherical_response.py` conservation and RN blocks; `evidence/retarded-dirac-spherical-response.json` |
| Dilaton source and the \(j=0\) Ward identity | same script, `massive_ward` and zeta blocks |
| Lorentz force in the Ward identity | same script, `check_lorentz_ward`: Christoffel divergence minus \(j^a F_{ba}\) |
| Current conservation, Gauss law, Schwinger mass | same script, `schwinger_mass`: \(\kappa^2\) solved from the Maxwell residual at \(\epsilon_{tx}=+\sqrt{-g}\), matched to arXiv:2610.09847v1 (1.1) |
| Both-mouth mass and charge integrals | same script, `mouths` block; residuals in the receipt |
| Clock coefficients fixed | same script, `modes` block; `delta C = -2 beta a`, Killing time-component limits |
| Not a new theorem | `docs/S051-clock-mode-obstruction-2026-10-12.md` prior-art section; pinned arXiv:1807.04726v3, 2609.33511v1, 2610.09847v1, 1804.00491v3 |
| No Lean acceptance | receipt field `lean_invocations: 0`; formalization boundary in the note |
