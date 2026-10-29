"""Author review and quantitative sensitivity of the coupled mass estimator."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_filtered_mass_adjoint_bound as m
arb, ctx, dec, enc = m.arb, m.ctx, m.dec, m.enc
ROOT = m.ROOT
digest = lambda p: sha256(p.read_bytes()).hexdigest()
path = ROOT/'evidence'/(m.STEM+'-attempt2.json')
raw = json.loads(path.read_bytes())
assert raw['status'] == 'filtered_mass_functional_bounded_additive_gluing_estimator_rejected'
assert raw['input_sha256_before'] == raw['input_sha256_after']
assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items())
assert raw['ordinary_argument_utf8'] == (ROOT/'notes'/(m.STEM+'.md')).read_text(encoding='utf-8')
assert raw['arithmetic_source_utf8'] == (ROOT/'checks/current_filtered_mass_adjoint_bound.py').read_text(encoding='utf-8')
protected = json.loads(m.BASELINE.read_bytes())
assert all(digest(Path(p))==h for p,h in protected.items())
assert digest(Path(sys.executable)) == m.runtime['python_sha256']
for pkg in m.runtime['packages'].values():
    assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())

ctx.prec = 768
r = raw['result']
c = {k:dec(v) for k,v in r['constants'].items()}
q, K, R = arb('.338'), arb('1.9999999764'), arb('2.1e-5')
fm = 1-arb('1.179e-8')
J = 8*arb('1.000001')*arb('12.001')
s0 = J/(fm*arb('23.99'))
ar2 = J*s0
d0 = (arb('24.001')*s0+J/fm)/6
sref = ((arb('24.001')/(fm*arb('23.99')))*ar2+6*arb('2.946e-9')*d0*d0).sqrt()
assert sref.overlaps(c['actual_scalar_reference_energy_upper'])
kf = 2*q.sqrt()*(sref+R)/(3*K)
kg = q/(1+q)**2*kf
eg = q/(1+q)**2*2*q.sqrt()*R/(3*K)
assert kg.overlaps(c['actual_filtered_mass_functional_norm_upper'])
assert eg.overlaps(c['filtered_mass_functional_reference_error_upper'])
assert kg < 1 and eg < arb('1e-6')

target = arb('1.49e34')
factor = (c['quantum_action_conversion_alpha_q']*arb('6e24')
          *c['original_normalized_mass_shape_L2_lower']
          *c['complete_structured_profile_combination'])
floor = kg*factor
assert floor.overlaps(c['separated_gluing_estimator_expression_lower'])
assert (floor/target).overlaps(c['current_gluing_coefficient_over_optimistic_ceiling'])
assert floor/target > arb('1e5')

# Sensitivity follows directly from ordinary equations (4) and (6):
# even granting the reference scalar component zero leaves the retained
# adjoint-residual debit eg. This is still a LOWER TEST OF AN UPPER
# ESTIMATOR EXPRESSION, not an actual error or nonzero physical lower bound.
residual_only_floor = eg*factor
assert residual_only_floor > target
residual_ceiling = R*target/residual_only_floor
assert residual_ceiling < R

for key in ['signed_r0_mass_pairing','signed_r1_mass_pairing',
            'actual_signed_mass_vector_norm','actual_Dj','actual_signed_Zr']:
    assert r[key] is None
for key in ['full_signed_mass_term_certified','actual_signed_response_certified',
            'retained_boundary_trace_transferred','available_physical_margin',
            'original_P32_defect','actual_kappa_m','complete_original_admission',
            'nonlinear_existence','exact_compatibility','physical_matching',
            'semiclassical_validity_noise','all_sector_stability',
            'construction_obstruction','full_goal_complete']:
    assert r['scope'][key] is False

ordinary = [
    'The FULL coupled source is (-b w,0). Both rows give xi=-F w. The canonical-to-TT mass factor is -3/2, so the original flux gives m_xi(Fw)=-2 b(s_Psi,w)/(3K_m).',
    'The reference auxiliary has the same root-zero trace and a regular bolt. The affine York solution has y=-2/(3z)+O(z log^k z); weighted Hardy and the z^4 test have vanishing endpoint terms and give ||xi_ref||L2<=8||t_ref||L2.',
    'Both reference rows vanish. Exact full reference-density differentiation gives the derivative-free mixed row and scalar equation (1); no full physical stationarity is asserted.',
    'The scalar coordinate, unlike the whole root-zero adjoint, is energy-admissible. Its positive reference scalar equation controls its Hmu norm and actual Laplacian, so the actual fourth-order scalar energy follows without a strong output comparison.',
    'The already certified COMPLETE commonly regulated adjoint residual and coupled inverse control s_Psi-s_ref in the actual scalar norm. All angular/thermal modes, contact and boundary terms remain in that residual.',
    'Extend the mass identity from original core inputs through the continuous scalar mixed pairing and the accepted full graph realization. This is not a trace on unfiltered E and is not a retained Ft boundary theorem.',
    'The spectral multiplier t/(1+t)^2 is increasing on [0,.338]. Its bounds give the actual filtered mass functional and the full-residual reference error; neither is an applied signed center.',
    'The fixed inner part of the original mass lift gives H_min. Its transition is not removed, and no derivative or global norm of an unspecified cutoff is assigned a numerical value.',
    'The global row-subtracted gluing coefficient is converted using 2G C/(B^2 beta), exactly as in the accepted weak residual. The complete same-source spectral profile bounds preserve original units and averaging radius.',
    'The positive estimator entry exceeds the entire conditional target even using H_min and granting every other entry zero. This rejects the specified additive global-L2 estimate only; it does not lower-bound actual residuals.',
    'Sensitivity: equation (4) leaves the positive certified residual term e_g even if S_ref is granted zero. Multiplying e_g by the same factor in (6) still exceeds the target. Refining S_ref alone cannot rescue this estimator.',
    'A separation-weighted kernel bound is not supplied by the fixed-buffer e^-M/270 and e^-M/530 estimates. Its joint local row and both signed pairings remain to be evaluated. No available physical margin is inferred.'
]
failed_path = ROOT/'evidence'/(m.STEM+'.json')
failed = json.loads(failed_path.read_bytes())
assert failed['status']=='failed' and 'joint_supremum' in failed['traceback']
assert 'assert dec(oldmass' in failed['arithmetic_source_utf8']
assert failed['ordinary_argument_utf8'] == raw['ordinary_argument_utf8']
out = ROOT/'evidence'/(m.STEM+'-attempt2-verification.json')
assert not out.exists(), 'Preserve immutable verification receipts'
rec = dict(status='coupled_mass_transfer_and_estimator_failure_reviewed',
    reviewed_utc=datetime.now(timezone.utc).isoformat(), receipt_sha256=digest(path),
    review_source_sha256=digest(Path(__file__)),
    frozen_inputs_checked=len(raw['input_sha256_after']), protected_paths_checked=len(protected),
    ordinary_review=ordinary, review_type='author ordinary argument and symbolic/interval arithmetic',
    sensitivity_constants=dict(
        unchanged_adjoint_residual_only_gluing_expression_lower=enc(residual_only_floor),
        unchanged_residual_only_expression_over_complete_target=enc(residual_only_floor/target),
        optimistic_adjoint_residual_ceiling_reference_scalar_and_other_errors_ZERO=enc(residual_ceiling)),
    mass_functional_norm_bounded=True, complete_signed_mass_term_certified=False,
    specified_global_L2_gluing_estimator_rejected=True,
    reference_scalar_refinement_alone_cannot_rescue_specified_estimator=True,
    independent_specialist_review=False, lean_verification=False,
    preserved_failed_attempt=dict(path=str(failed_path.relative_to(ROOT)),sha256=digest(failed_path),
        cause='Receipt re-encoding interval compared <= the exact declared B0 guard; use exact frozen guard and overlap consistency.'),
    original_admission=False, construction_obstruction=False, full_goal_complete=False,
    review_source_utf8=Path(__file__).read_text(encoding='utf-8'))
out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
print(json.dumps({k:v for k,v in rec.items() if k not in ['ordinary_review','review_source_utf8']},indent=2))
