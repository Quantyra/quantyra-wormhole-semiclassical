"""Ordinary interpretation and conditional ledger, frozen before arithmetic.

The 440f2e0 certificate has an ACTUAL paired-norm lower bound R for its
specified averaged REFERENCE vector, after the whole quartic subtraction.
Its already bounded first-contribution approximation costs at most K/pref
in energy; using the entire known subtotal K is conservative. Therefore
the actual residual Zr satisfies

    Zr >= R-K/pref-Dj.

This is not a lower bound on the original output operator, since Zr enters
the specified adjoint-trial UPPER estimator for that operator. A complete
Dj<=1.49e34 would make pref*Zr exceed 1.667e-9 and reject that particular
signed-trial certificate. It would not certify admission or obstruct the
physical construction. Dj and Zr themselves remain uncomputed here.

The present full-to-retained discrepancy account is

    Ebridge+conversion*(Eshape+E_local_scalar+E_fast+E_rest+E_app)<=1.49e34.

Only Ebridge, conversion and the entire massive/LLL shape block Eshape
have bounds in this record. Their remaining conditional account must NOT
be assigned to unknown terms. In the complete scalar/local identities,
the isolated observed k_m term has coefficient at most q*K_g*alpha_radius;
so, granting all other unknown terms zero, |A_r k_m| must be bounded by
remaining/(q*K_g*alpha_radius) to make this additive estimator sufficient.
This is a diagnostic requirement, not an attained value or physical margin.

All formulae use unchanged original normalization, actual paired domain,
complete action and prescribed averaged source. No new PDE theorem or
numerical operator application is asserted by this arithmetic ledger.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_mass_pairing_derivative_return as m
ROOT=m.ROOT;arb,ctx,dec,enc=m.arb,m.ctx,m.dec,m.enc
digest=lambda p:sha256(p.read_bytes()).hexdigest()
stem='current-mass-pairing-diagnostic-ledger'
paths=['checks/current_mass_pairing_diagnostic_ledger.py',
    'evidence/current-mass-pairing-derivative-return.json',
    'evidence/current-mass-pairing-derivative-return-verification.json',
    'evidence/current-correlated-signed-trial-test-attempt3.json',
    'notes/current-correlated-signed-trial-test.md',
    'notes/current-filtered-j-regularity-test.md',
    'evidence/current-filtered-mass-adjoint-bound-attempt2.json',
    'evidence/current-filtered-mixed-green-residual.json']
before={p:digest(ROOT/p) for p in paths}
out=ROOT/'evidence'/(stem+'.json');assert not out.exists()
ctx.prec=768
c=m.m.constants('current-mass-pairing-derivative-return')
ref=m.m.constants('current-correlated-signed-trial-test-attempt3')
adj=m.m.constants('current-filtered-mass-adjoint-bound-attempt2')
green=m.m.constants('current-filtered-mixed-green-residual')
target=arb('1.49e34');pref=green['pref_t']
R=ref['complete_trial_energy_lower'];K=ref['known_error_subtotal_NOT_margin']
bolt=ref['whole_bolt_allowance']
assert bolt.overlaps(arb('1.667e-9'))
known_gap=R-K/pref-bolt/pref
assert known_gap.overlaps(ref['decisive_remaining_j_comparison_threshold'])
conditional_Z_lower=R-K/pref-target
conditional_physical_lower=pref*conditional_Z_lower
assert conditional_physical_lower>bolt
remaining=c['conditional_full_E_account_after_bridge_and_shape']
coefficient=arb('.338')*adj['actual_filtered_mass_functional_norm_upper']*green['alpha_original_radius_cost']
scalar_ceiling=remaining/coefficient
assert remaining>0 and scalar_ceiling>0
after={p:digest(ROOT/p) for p in paths};assert before==after
protected=json.loads(m.BASELINE.read_bytes())
assert all(digest(Path(p))==h for p,h in protected.items())
rec=dict(status='diagnostic_target_and_remaining_unknowns_explicit',
    executed_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,input_sha256_after=after,
    ordinary_argument_and_source_utf8=Path(__file__).read_text(encoding='utf-8'),
    protected_paths_unchanged=len(protected),
    constants={k:enc(v) for k,v in dict(
        accepted_reference_center_energy_lower=R,unchanged_known_subtotal_NOT_margin=K,
        accepted_diagnostic_gap=known_gap,conditional_Dj_target=target,
        conditional_Zr_energy_lower_IF_complete_Dj_target_passes=conditional_Z_lower,
        conditional_original_bolt_estimator_lower_IF_complete_Dj_target_passes=conditional_physical_lower,
        conditional_estimator_overrun_IF_complete_Dj_target_passes=conditional_physical_lower/bolt,
        attained_massive_LLL_discrepancy_upper=c['complete_massive_and_LLL_mass_pairing_upper'],
        remaining_unknown_discrepancy_account_CONDITIONAL=remaining,
        actual_scalar_observation_coefficient=coefficient,
        scalar_observation_ceiling_ALL_other_unknown_errors_zero=scalar_ceiling).items()},
    actual_Dj=None,actual_Zr=None,actual_scalar_observation=None,
    complete_signed_trial_screen_rejected=False,complete_physical_output_certified=False,
    original_admission=False,nonlinear_existence=False,exact_compatibility=False,
    physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False,
    construction_obstruction=False,full_goal_complete=False)
out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
print(json.dumps(dict(status=rec['status'],constants={k:v['ball'] for k,v in rec['constants'].items()}),indent=2))
