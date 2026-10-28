"""Certify full omitted tails; the retained signed heat window is not evaluated."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-quartic-positive-time-window'
runtime = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = ['notes/'+STEM+'.md', 'checks/current_quartic_positive_time_window.py',
    'notes/current-quartic-column-energy-dual.md',
    'evidence/current-quartic-column-energy-dual-assembly.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
    'checks/current_subtracted_thermal_action.py',
    'notes/current-subtracted-thermal-action.md',
    'checks/current_disk_heat_propagation.py',
    'notes/current-disk-heat-propagation.md',
    'notes/current-canonical-heat-restoration-assembly.md']


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
    old=json.loads((ROOT/'evidence/current-quartic-column-energy-dual-assembly.json').read_bytes())
    assert old['status']=='column_bound_arithmetic_passed_full_residual_unresolved'
    c=old['result']['constants']
    er=json.loads((ROOT/'evidence/current-coupled-truncation-energy-error.json').read_bytes())['runs'][-1]
    assert er['status']=='passed'
    p=er['result']['parameters']
    C,B2,beta=[dec(p[k]) for k in ['C','B2','beta']]
    CL,M,P4=[dec(c[k]) for k in ['CLt','heat_diagonal_bound_M','P4_L2']]
    CH,C0,C2=64*M*CL/3,2*CL/(3*arb.pi()),CL*P4/(480*arb.pi())
    t0=arb('.001');pref=1/(B2*beta);q=arb(2000)**4;a=(q+1)/C
    assert a>324
    alpha=dec(c['sharper_exact_one_dimensional_coefficient_envelope'])
    cut=dec(c['sharper_one_dimensional_cut_residual_debit'])
    target=arb('1e-30');need=(target-cut)/alpha

    def time_tail(T):
        return pref*4*C/t0*(-a*(T-t0)).exp()/a*(CH/T+C0/T**2+C2)

    def angular_tail(N):
        x=(-a*t0).exp()
        A=(2*q+4*(N+1))/(1-x)+4*x/(1-x)**2
        return pref*A*x**(N+1)/(a*(N+1))*(CH/t0+C0/t0**2+C2)

    time_trials={s:time_tail(arb(s)) for s in ['.04','.045','.05','.055']}
    angular_trials={str(n):angular_tail(n) for n in [40,45,50,55,60]}
    T,N=arb('.05'),60
    dt,da=time_tail(T),angular_tail(N)
    W=arb(0)
    cells=[]
    edges=[t0*(T/t0)**(arb(i)/256) for i in range(257)]
    for lo,hi in zip(edges,edges[1:]):
        theta=sum((2*(q+2*n)*(-arb(n)*(n+q)*lo/C).exp()
                   for n in range(1,N+1)),arb(0))
        contribution=(hi-lo)*theta/lo
        assert contribution>0
        W+=contribution
        cells.append({'lo':enc(lo),'hi':enc(hi),'upper_sum_contribution':enc(contribution)})
    remaining=(need-dt-da)/CL
    rho=arb('.002');enum=arb('9e-32')
    prospective=cut+alpha*(dt+da+CL*(pref*W*rho+enum))
    assert dt<arb('3.552e-31') and da<arb('3.468e-33')
    assert prospective<target
    assert remaining>0
    constants={'actual_angular_gap_lower':a,'full_time_tail':dt,
        'full_angular_tail':da,'truncated_time_weight_upper_expression':W,
        'cap_transport_time_and_angular_debit':cut+alpha*(dt+da),
        'remaining_column_error_threshold':need-dt-da,
        'remaining_normalized_L2_density_error_threshold':remaining,
        'uniform_kernel_L2_threshold_if_all_middle_error_spent_there':remaining/(pref*W),
        'prospective_uniform_kernel_error_guard':rho,
        'prospective_extra_weighted_L2_numerical_error_guard':enum,
        'conditional_one_column_residual_upper':prospective,
        'conditional_fraction_of_one_column_target':prospective/target,
        'whole_positive_time_upper_still_current':dec(c['canonical_positive_time_bound'])}
    return {'t0':'.001','T':'.05','angular_N':N,
        'constants':{k:enc(x) for k,x in constants.items()},
        'time_trials':{k:enc(x) for k,x in time_trials.items()},
        'angular_trials':{k:enc(x) for k,x in angular_trials.items()},
        'weight_cells':cells,
        'actual_signed_window':None,'actual_uniform_kernel_error':None,
        'actual_weighted_numerical_error':None,'full_residual_bound':None,
        'scope':{'complete_omitted_time_and_angular_tails':True,
            'full_thermal_kernel_still_required':True,
            'retained_window_evaluated':False,'actual_application_error_certified':False,
            'coupled_budget_closed':False,'physical_obstruction_proved':False,
            'nonlinear_existence_or_stability_proved':False}}


if __name__=='__main__':
    path=ROOT/'evidence'/(STEM+'.json')
    assert not path.exists(), 'Preserve every immutable attempt.'
    before={p:digest(ROOT/p) for p in INPUTS}
    receipt={'started_utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256_before':before,'python':sys.executable,
        'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result']=execute()
        receipt['status']='tail_contract_passed_signed_window_and_full_residual_unresolved'
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
