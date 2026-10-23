"""Exact requirement ledger and provenance audit; NOT a physical solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-controlled-solution-audit-closeout'
PREVIOUS = 'current-coincident-dirac-action-assembly'
OWN = ['notes/' + STEM + '.md', 'checks/current_controlled_solution_audit_closeout.py',
       'evidence/' + STEM + '-review.json']
ROWS = []
NUMERICS = []


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def check(name, condition):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name, passed=bool(condition)))
    if not condition:
        raise AssertionError(name)


def interval(encoded):
    a = encoded['decimal_enclosure']
    # Exact rational scale: no rounded printed Arb ball is used.
    scale = F(10) ** int(a['exponent'])
    return (F(a['mid']) - F(a['rad'])) * scale, (F(a['mid']) + F(a['rad'])) * scale


def receipt(stem):
    run = json.loads((ROOT / ('evidence/' + stem + '.json')).read_bytes())['runs'][-1]
    check(stem + '_passed', run['status'] == 'passed')
    return run


def threshold(name, value, status='unknown_operator_or_embedding'):
    NUMERICS.append(dict(name=name, exact_threshold=str(value), status=status))


def arithmetic():
    rho, y, z0, z2 = map(F, ['1e-8', '5e-9', '.25', '5e7'])
    check('preservation_at_selected_ceilings', y + rho*z0 + rho*rho*z2/2 == rho)
    check('strict_contraction_at_selected_ceilings', z0 + rho*z2 == F('.75') < 1)
    for name, value in [
        ('K_m', F('1e-9')/F('1.632e-11')),
        ('K_c', F('2e-9')/F('1.36e-9')),
        ('G_half_low', F('1e-9')/F('1.73e-15')),
        ('E45_given_baseline', (F('1e-3')-F('1.3e-4'))/rho),
        ('cold_baseline_given_unit_embedding', F('5e-7')-rho),
        ('outer_baseline_given_unit_embedding', F('1.5e-7')-rho),
        ('E_D_given_zero_baseline', F('2.5e-10')/rho),
        ('quarter_gain_embedding_product', z2/F('2.01e10')),
        ('quarter_gain_given_E_D_025', z2/F('2.01e10')/F('.025')**2),
        ('half_bolt_residual_product', F('1.667e-9')/2)
    ]:
        threshold(name, value)
    check('quarter_product_exact', z2/F('2.01e10') == F(1, 402))
    check('unit_direction_gain_ceiling', z2/F('2.01e10')/F('.025')**2 == F(800, 201))
    scalar = F('1.1e-31')*F(4, 3)*10**18
    check('scalar_only_return', scalar < F('1.466666667e-13'))
    required = F('1.667e-9')/(F(4, 3)*10**18) - F('1.1e-31')
    # The trace input is STRICTLY below its allowance. Equality in the
    # remaining shape allocation therefore still gives a strict output.
    check('shape_to_trace_target_sufficient', F('1.25014e-27') == required)
    threshold('shape_to_trace_allowance_from_scalar_budget', required)
    threshold('scalar_only_return_upper', scalar, 'established_scoped_upper_bound')
    check('R2_refinement_materiality', F('1.06e-18')/F('1.73e-15') < F('.00062'))
    for gain in [F('469.12'), F('1e-9')/F('1.73e-15')]:
        allowance = F('2.01e10') * gain * F('.025')**2
        NUMERICS.append(dict(name='conditional_quarter_product', hypothetical_gain=str(gain),
                            upper_allowance=str(allowance), ratio_to_budget=str(allowance/z2),
                            E_D_squared_threshold=str(F(1, 402)/gain),
                            status='hypothetical_gain_not_actual_operator_constant'))
        check('hypothetical_gain_fails_' + str(gain), allowance > z2)
    check('conditional_469_ratio', F('2.01e10')*F('469.12')*F('.025')**2/z2 == F('117.8664'))
    for name, bound, lower in [
        ('old_interior', '4.141944868545865e118', '1.7374e143'),
        ('old_global', '1.069626977892751e201', '4.4868e225')
    ]:
        ratio = F(bound)/F('2.38393e-25')
        check(name + '_known_failed_certificate', ratio > F(lower))
        NUMERICS.append(dict(name=name, ratio_to_old_ceiling=str(ratio),
                            status='failed_upper_certificate_not_solution_lower_bound'))
    check('old_low_term_known_failed_certificate', F('.97') > F('3.90858e-8'))


def scoped_evidence(old):
    regularity = receipt('current-reference-regularity-feasibility')
    row = next(r for r in regularity['numerics'] if r['kind'] == 'current_direction_simultaneous_integration_budget' and r['delta'] == '.01')
    for key, printed in [
        ('all_mode_value_remainder_upper', '3.28756e-38'),
        ('principal_Lipschitz_derivative_remainder_upper', '7.54498e-37')
    ]:
        upper = interval(row[key])[1]
        check(key + '_fits_display_and_budget', upper < F(printed) < F('1e-28'))
        NUMERICS.append(dict(name=key, source=row[key], display_upper=printed,
                            status='established_actual_direction_component_only'))
    contact = receipt('current-vacuum-contact-feasibility-closeout')
    row = next(r for r in contact['numerics'] if r['kind'] == 'available_vacuum_contact_budget')
    for key, printed in [('interior_value_allowance', '2.86945e-33'),
                         ('interior_derivative_allowance', '1.84123e-30')]:
        check(key + '_fits_display_and_budget', interval(row[key])[1] < F(printed) < F('1e-28'))
    row = next(r for r in contact['numerics'] if r['kind'] == 'physical_bolt_metric_combination_tolerance')
    lo, hi = interval(row['metric_combination_tolerance'])
    check('physical_tolerance_display_encloses_threshold', F('8.43135e-43') < lo < hi < F('8.43136e-43'))
    NUMERICS.append(dict(name='metric_combination_tolerance', source=row['metric_combination_tolerance'],
                        display_interval=['8.43135e-43', '8.43136e-43'], actual_correction_computed=False))
    diagonal = receipt('current-seed-diagonal-restoration-contract')
    row = next(r for r in diagonal['numerics'] if r['kind'] == 'contact_geometry_enclosure' and r['cells'] == 128)
    check('a2_value_uses_its_own_target', interval(row['contact_value_upper'])[1] < F('1e-27'))
    check('a2_derivative_uses_its_own_target', interval(row['contact_first_derivative_upper'])[1] < F('1e-24'))
    check('raw_mode_scope_not_promoted', old['coincident_finite_mode_action_assembly_executed']
          and old['lowest_thermal_physical_bolt_matrix_enclosed']
          and not old['complete_renormalized_operator']
          and not old['coincident_subtraction_executed']
          and not old['actual_physical_correction_computed'])
    final = [r for r in old['numerics'] if
             (r['kind'] == 'continuous_finite_mode_action_matrix' and r['subdivisions'] == 1024) or
             (r['kind'] == 'lowest_thermal_bolt_action_matrix' and r['subdivisions'] == 256)]
    check('finite_mode_development_width_only', len(final) == 5 and all(
        interval(w)[1] < F('.03252') for r in final for w in r['relative_widths']))


def main():
    review = json.loads((ROOT / OWN[2]).read_bytes())
    old = json.loads((ROOT / ('evidence/' + PREVIOUS + '.json')).read_bytes())['runs'][-1]
    paths = sorted(set(old['input_sha256_after']) | set(OWN) | set(review['retained_inputs']) |
                   {'evidence/' + PREVIOUS + '.json', 'evidence/' + PREVIOUS + '.log'})
    before = {p: digest(p) for p in paths}
    run = dict(started_utc=datetime.now(timezone.utc).isoformat(), input_sha256_before=before,
               frozen_sources_utf8={p: (ROOT / p).read_text(encoding='utf-8') for p in OWN},
               protocol=review['protocol'], invocation=sys.argv,
               python_version=sys.version, executable_sha256=sha256(Path(sys.executable).read_bytes()).hexdigest())
    error = None
    try:
        check('previous_execution_passed', old['status'] == 'passed')
        check('all_inherited_inputs_unchanged', all(digest(p) == h for p, h in old['input_sha256_after'].items()))
        arithmetic()
        scoped_evidence(old)
    except Exception:
        error = traceback.format_exc()
    after = {p: digest(p) for p in paths}
    status = 'passed' if error is None and before == after else 'failed'
    run.update(status=status, exception=error, input_sha256_after=after, unchanged_inputs=before == after,
               checks=ROWS, numerics=NUMERICS, finished_utc=datetime.now(timezone.utc).isoformat(),
               audit_closed=True, weak_primal_dual_method_selected=True,
               full_solution_feasibility_test_executed=False, actual_physical_correction_computed=False,
               complete_renormalized_operator=False, full_forward_residual_norm_bounded=False,
               full_adjoint_residual_norm_bounded=False, raw_C5_admission_certified=False,
               nonlinear_solution_constructed=False, physical_stability_proved=False,
               semiclassical_validity_proved=False, component_complete=False, goal_complete=False,
               formal_verification=False)
    out = ROOT / ('evidence/' + STEM + '.json')
    result = json.loads(out.read_bytes()) if out.exists() else {'runs': []}
    result['runs'].append(run)
    out.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    with (ROOT / ('evidence/' + STEM + '.log')).open('a', encoding='utf-8', newline='\n') as log:
        for row in [dict(event='start', run=len(result['runs']), source_sha256={p: before[p] for p in OWN}),
                    *ROWS, dict(event='finish', status=status, exception=error)]:
            log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status, checks=len(ROWS), inputs=len(paths), exception=error,
                         full_solution_feasibility_test_executed=False)), flush=True)
    if status != 'passed':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
