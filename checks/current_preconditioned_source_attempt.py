"""Audit reference/error accounting; does not evaluate the complete source map."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-preconditioned-source-attempt'
INPUTS = ['notes/'+STEM+'.md', 'checks/current_preconditioned_source_attempt.py',
          'notes/current-geodesic-flow-multiplier.md',
          'notes/current-global-geodesic-return.md',
          'checks/current_geodesic_flow_multiplier.py',
          'checks/current_global_geodesic_return.py',
          'evidence/current-geodesic-flow-multiplier.json',
          'evidence/current-global-geodesic-return.json',
          'evidence/current-quartic-ritz-rank-test-final.json',
          'evidence/current-coupled-mass-adjoint-runtime.json']
runtime = json.loads((ROOT/INPUTS[-1]).read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def dec(x):
    a=x['decimal_enclosure']
    return arb(int(a['mid']), int(a['rad']))*arb(10)**int(a['exponent'])


def enc(x):
    m,r,e=x.mid_rad_10exp(60)
    return {'ball':x.str(40),'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}


def main():
    ctx.prec=512
    out=ROOT/'evidence'/(STEM+'.json')
    assert not out.exists(), 'Preserve immutable execution evidence.'
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for package in runtime['packages'].values():
        for name,h in package['files'].items():
            assert digest(Path(package['root'])/name)==h
    protected=json.loads(Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json').read_bytes())
    assert all(digest(Path(p))==h for p,h in protected.items())
    before={p:digest(ROOT/p) for p in INPUTS}
    started=datetime.now(timezone.utc).isoformat()
    flow=json.loads((ROOT/'evidence/current-geodesic-flow-multiplier.json').read_bytes())['runs'][-1]
    rank=json.loads((ROOT/'evidence/current-quartic-ritz-rank-test-final.json').read_bytes())
    assert flow['status']=='passed' and rank['status']=='passed'
    assert not flow['geodesic_flow_matrix_evaluated']
    f=next(x for x in flow['numerics'] if x['kind']=='actual_flow_quadrature_form_budget')
    c=rank['result']['constants']
    eps,qt,ft,ed=[dec(c[k]) for k in ['epsilon','q_t','f_t','alternative_cutoff_eta']]
    eq=(dec(f['previous_eta_upper'])+dec(f['eta_quad'])+eps)/(1-eps)
    equad=dec(f['eta_quad'])/(1-eps)
    e0=(arb('.0022')+eps)/(1-eps)
    external=equad*ft/((1-e0)*(1-e0-equad))
    contraction=(eq+ed)/(1-eq)
    residual=(1+qt)*ft*(ed+contraction**12)/(1-ed)
    target=arb('5e-32')
    rho=(1-eq)*target
    assert eq<arb('.0022') and rho>arb('4.989e-32')
    assert external/target>arb('1e26')
    values=dict(combined_Et_reference_error=eq,
                slack_in_existing_point_0022_guard=arb('.0022')-eq,
                external_fixed_reference_upper_error=external,
                external_upper_report_over_application_target=external/target,
                sufficient_source_residual_threshold_NOT_attained=rho,
                safe_rounded_source_residual_target_NOT_attained=arb('4.989e-32'),
                exact_application_residual_with_recorded_comparison=residual,
                unchanged_rounded_guard_exact_application_residual=dec(c['alternative_exact_application_residual_upper']))
    after={p:digest(ROOT/p) for p in INPUTS}
    assert before==after and all(digest(Path(p))==h for p,h in protected.items())
    data={'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat(),
          'status':'accounting_checked_complete_source_application_unresolved',
          'input_sha256_before':before,'input_sha256_after':after,
          'python_sha256':runtime['python_sha256'],'protected_paths_unchanged':len(protected),
          'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
          'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8'),
          'constants':{k:enc(v) for k,v in values.items()},
          'actual_source_residual':None,'actual_source_application_error':None,
          'actual_K_application_error':None,'complete_numerical_residual':None,
          'new_consequential_complete_bound':False,'physical_margin_closed':False,
          'goal_complete':False,'formal_verification':False}
    out.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    for k,v in values.items():print(k+': '+v.str(30))
    print('Complete S application remains unresolved; no new complete bound.')


if __name__=='__main__':
    main()
