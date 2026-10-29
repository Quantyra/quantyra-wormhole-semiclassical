"""Conditional first-quantum application contract; DOES NOT apply the inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-first-quantum-solve-attempt'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = ['notes/'+STEM+'.md', 'checks/current_first_quantum_solve_contract.py',
          'checks/current_quartic_ritz_rank_test.py',
          'evidence/current-spatial-inverse-output-test.json',
          'evidence/current-first-inverse-source-screen.json',
          'evidence/current-scalar-green-source-residual.json',
          'evidence/current-scalar-green-source-residual-verification.json',
          'notes/current-geodesic-two-sided-remainder.md',
          'notes/current-global-geodesic-return.md',
          'notes/current-local-reference-range-test.md',
          'operators/radial_geodesic_return.py',
          'operators/angular_resummed_response.py',
          'checks/current_subtracted_thermal_action.py',
          'evidence/current-coupled-mass-adjoint-runtime.json']


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def result(name):
    d = json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
    run = d.get('runs',[d])[-1]
    assert run['status'] != 'failed'
    return run['result']


def execute():
    ctx.prec = 512
    spatial = result('current-spatial-inverse-output-test')
    c = {k:dec(v) for k,v in spatial['constants'].items()}
    first = {k:dec(v) for k,v in result('current-first-inverse-source-screen')['constants'].items()}
    source = {k:dec(v) for k,v in result('current-scalar-green-source-residual')['constants'].items()}
    eps,eta,e0 = c['epsilon'],c['alternative_eta'],c['geodesic_Et_comparison_error']
    assert c['alternative_delta'].contains(arb('1e-38'))
    assert spatial['internal_Neumann_terms'] == 12
    r = (e0+eta)/(1-e0)
    assert r.overlaps(c['Neumann_contraction_upper']) and r < 1
    q = arb('.338')
    source_error = (q/(1-eps)).sqrt()*source['full_scalar_residual_Es_dual_upper']
    Ghat = first['complete_source_dual_upper']+source_error
    analytic_relative = (eta+r**12)/(1-eta)
    analytic = analytic_relative*Ghat
    screen = arb('1e33')
    stage_sum = (1-r**12)/(1-r)
    conditional_stage = (screen-analytic)/stage_sum
    assert analytic < screen and conditional_stage > 0
    Kfirst = first['first_residual_to_original_bolt_multiplier']
    constants = dict(
        proved_source_dual_upper=Ghat,
        scalar_source_dual_replacement_error=source_error,
        retained_spatial_cutoff_eta=eta,complete_geodesic_form_error=e0,
        first_inverse_contraction_upper=r,
        exact_application_analytic_relative_error=analytic_relative,
        exact_application_first_residual_error=analytic,
        exact_application_analytic_physical_debit=Kfirst*analytic,
        proposed_residual_screen=screen,
        conditional_common_stage_output_error_ceiling=conditional_stage,
        conditional_stage_ceiling_over_source_upper=conditional_stage/Ghat,
        stage_error_amplification_sum=stage_sum,
        physical_debit_if_complete_screen_attained=Kfirst*screen,
        partial_subtotal_if_screen_attained_NOT_margin=
            source['known_partial_subtotal_not_free_margin']+Kfirst*screen)
    return dict(precision_bits=ctx.prec,internal_iterations_unchanged=12,
                constants={k:enc(v) for k,v in constants.items()},
                actual_first_geodesic_application=None,
                actual_quantum_remainder_application=None,
                actual_stage_errors=None,actual_signed_residual=None,
                complete_signed_response=None,
                scope=dict(seed_form_first_inverse_defined=True,
                           correct_reference_is_ell0_not_ell0_plus_B=True,
                           conditional_application_contract_only=True,
                           actual_quantum_application_executed=False,
                           proposed_screen_certified=False,
                           new_signed_feasibility_bound=False,
                           full_bolt_or_coupled_admission=False,
                           nonlinear_existence_matching_validity_stability=False,
                           construction_obstruction=False,full_goal_complete=False))


if __name__ == '__main__':
    path = ROOT/'evidence'/(STEM+'.json')
    if path.exists():
        raise RuntimeError('Preserve immutable attempt receipt.')
    before = {p:digest(ROOT/p) for p in INPUTS}
    protected = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p)) == h for p,h in protected.items())
    runtime = json.loads((ROOT/INPUTS[-1]).read_bytes())
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    pkg = runtime['packages']['flint']
    assert all(digest(Path(pkg['root'])/p) == h for p,h in pkg['files'].items())
    rec = dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
               ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
               arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
               python_sha256=digest(Path(sys.executable)))
    try:
        rec['result'] = execute()
        after = {p:digest(ROOT/p) for p in INPUTS}
        assert after == before
        assert all(digest(Path(p)) == h for p,h in protected.items())
        rec.update(status='first_quantum_application_unresolved_contract_only',
                   input_sha256_after=after,protected_paths_unchanged=len(protected),
                   pinned_flint_files_unchanged=len(pkg['files']))
    except Exception:
        import traceback
        rec.update(status='failed',traceback=traceback.format_exc())
    rec['finished_utc'] = datetime.now(timezone.utc).isoformat()
    path.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    cs = rec.get('result',{}).get('constants',{})
    print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
                         constants={k:v['ball'] for k,v in cs.items()}),indent=2))
    if rec['status']=='failed':
        raise SystemExit(1)
