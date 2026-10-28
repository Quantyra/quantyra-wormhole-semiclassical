"""Test a positive estimator's extension; do not claim a signed evaluation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-heat-window-extension-test'
runtime = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = ['notes/'+STEM+'.md', 'checks/current_heat_window_extension_test.py',
    'checks/current_disk_heat_propagation.py', 'notes/current-disk-heat-propagation.md',
    'notes/current-quartic-positive-time-window.md',
    'evidence/current-quartic-column-energy-dual-assembly.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def dec(x):
    d=x['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def enc(x):
    m,r,e=x.mid_rad_10exp(60)
    return {'ball':x.str(40),'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def execute():
    ctx.prec=512
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for package in runtime['packages'].values():
        for name,h in package['files'].items():
            assert digest(Path(package['root'])/name)==h
    def cut(t):
        return arb('1.14e8')/t**6*(-1/(16*t)).exp()
    def floor(t):
        return 6*cut(t)*t**4/140
    lo,hi=arb('.002'),arb('.003')
    assert floor(lo)<1 and floor(hi)>1 and hi<arb(1)/32
    for _ in range(100):
        mid=(lo+hi)/2
        if floor(mid)<1:
            lo=mid
        elif floor(mid)>1:
            hi=mid
        else:
            raise AssertionError('Ambiguous interval bisection')
    old=json.loads((ROOT/'evidence/current-quartic-column-energy-dual-assembly.json').read_bytes())
    assert old['status']=='column_bound_arithmetic_passed_full_residual_unresolved'
    c=old['result']['constants']
    prior=json.loads((ROOT/'evidence/current-coupled-truncation-energy-error.json').read_bytes())['runs'][-1]
    assert prior['status']=='passed'
    p=prior['result']['parameters']
    C,B2,beta=[dec(p[k]) for k in ['C','B2','beta']]
    CL,M,P4=[dec(c[k]) for k in ['CLt','heat_diagonal_bound_M','P4_L2']]
    CH,C0,C2=64*M*CL/3,2*CL/(3*arb.pi()),CL*P4/(480*arb.pi())
    a=(arb(2000)**4+1)/C
    alpha=dec(c['sharper_exact_one_dimensional_coefficient_envelope'])
    T=arb('.05')
    tail=4*C/arb('.001')*(-a*(hi-arb('.001'))).exp()/a*(CH/hi+C0/hi**2+C2)/(B2*beta)
    y=(840*floor(T)).root(4)
    phi=3*(y.sinh()-y.sin())/y**3
    rho=4*T**3*cut(T)/(128*arb.pi())*(phi-1)
    log10=lambda x:x.log()/arb(10).log()
    radius=(4*T*(arb(6)*arb('1.14e8')/140/T**2).log()).sqrt()
    assert floor(T)>arb('5.59912e8')
    assert alpha*tail/arb('1e-30')>arb('1.0012e7')
    assert log10(phi)>351 and log10(rho/arb('.002'))>363
    assert radius>arb('2.06849')
    values={'critical_time_lower':lo,'critical_time_upper':hi,
        'ratio_floor_at_T':floor(T),'required_cutoff_distance_at_T':radius,
        'tail_upper_report_at_critical_time':tail,
        'one_dimensional_tail_report':alpha*tail,
        'one_dimensional_tail_report_over_screen':alpha*tail/arb('1e-30'),
        'log10_entire_series_factor_floor':log10(phi),
        'log10_uniform_kernel_error_report_floor':log10(rho),
        'log10_report_over_conditional_kernel_allocation':log10(rho/arb('.002'))}
    return {'constants':{k:enc(x) for k,x in values.items()},
        'ratio_floors':{t:enc(floor(arb(t))) for t in ['.001','.002','.0025','.01','.05']},
        'scope':'Necessary-precondition failure of the specified larger-time positive estimator recipe; no larger-time operator bound asserted.',
        'actual_signed_window':None,'actual_application_error':None,
        'full_residual_bound':None,'coupled_budget_closed':False,
        'construction_obstruction_proved':False}

if __name__=='__main__':
    path=ROOT/'evidence'/(STEM+'.json')
    assert not path.exists(), 'Preserve immutable attempts.'
    before={p:digest(ROOT/p) for p in INPUTS}
    receipt={'started_utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256_before':before,'python':sys.executable,
        'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result']=execute()
        receipt['status']='specified_estimator_rejected_signed_window_and_full_residual_unresolved'
    except BaseException as exc:
        receipt['status']='failed';receipt['failure']=repr(exc)
        raise
    finally:
        receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert before==receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'],'constants':{
        k:x['ball'] for k,x in receipt['result']['constants'].items()}},indent=2))
