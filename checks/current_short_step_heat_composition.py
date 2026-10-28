"""Certify analytic composition bias; numerical kernel errors remain unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
STEM='current-short-step-heat-composition'
runtime=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb,ctx

INPUTS=['notes/'+STEM+'.md','checks/current_short_step_heat_composition.py',
    'notes/current-rough-paired-variation.md','notes/current-disk-heat-propagation.md',
    'checks/current_disk_heat_propagation.py','evidence/current-disk-heat-propagation.json',
    'evidence/current-quartic-positive-time-window.json',
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
    old=json.loads((ROOT/'evidence/current-disk-heat-propagation.json').read_bytes())['runs'][-1]
    assert old['status']=='passed' and old['complete_actual_domain_remainder_bounded']
    rows={a['kind']:a for a in old['numerics'] if 'kind' in a}
    assert next(a for a in old['checks'] if a['name']=='real_geodesic_metric_gradient_guard')['passed']
    keys=['base','metric','time','mixed']
    r,sigma,g,theta=arb('.1'),arb('.25'),arb('.5'),arb('.1')
    weights=[arb(1),r,sigma,r*sigma]
    uncut=rows['full_uncut_transport_norms']['tables']
    cap=rows['conditional_complete_cap_kernel_table']['table']
    budget=rows['complete_actual_domain_propagation_budget']
    comm=rows['complete_commutator_and_thermal_tails']
    P,F={},{}
    for k in keys:
        P[k],F[k]={},{}
        for typ in ['inf','row']:
            cut=dec(budget['cutoff_jet_norm'][typ])
            co=dec(comm['commutator'][typ])
            thermal=dec(comm['full_thermal_image_guard'])
            F[k][typ]=dec(uncut['F'][k][typ])+cut+co+thermal
            cp=dec(cap[k][typ])+thermal
            interior=dec(uncut['P'][k][typ])+cut*arb('.001')**3
            P[k][typ]=(arb(max((2*cp).upper(),interior.upper()))
                        if typ=='inf' else cp+interior)
    b,mt,tm,mx=keys
    bbar=2*F[b]['row']+2*F[mt]['row']+F[tm]['row']+F[mx]['row']
    ratio=arb('.001')**4*bbar/140
    assert ratio<arb('1e-9')
    cv=lambda i,j:P[i]['row']*F[j]['row']/4
    row_first=[cv(b,b),cv(mt,b)+cv(b,mt),
        cv(b,b)+cv(tm,b)+cv(b,tm),
        cv(mt,b)+cv(b,mt)+cv(mx,b)+cv(tm,mt)+cv(mt,tm)+cv(b,mx)]
    row_higher=sum(row_first,arb(0))*ratio/(1-ratio)
    row=[x+row_higher for x in row_first]
    point=[dec(budget['first_actual'][k])+dec(budget['higher_Volterra']) for k in keys]
    B0=sum((w*(P[k]['inf']*P[k]['row']).sqrt() for w,k in zip(weights,keys)),arb(0))
    d0=sum((w*(a*b).sqrt() for w,a,b in zip(weights,point,row)),arb(0))
    e0=sum((w*a for w,a in zip(weights,row)),arb(0))
    f0=sum((w*a for w,a in zip(weights,point)),arb(0))
    difference=(2*r).exp()-1
    cross=(1+sigma)*r*(2*r).exp()*g
    omega=cross**2/(4*theta)+(1+sigma)*r**2*(2*r).exp()*g**2/4
    coercivity=1-sigma-(1+sigma)*difference-theta
    T=arb('.05');C=4*(omega*T).exp()
    assert coercivity>arb('.37324') and omega<arb('.015523')
    assert e0<arb('131.451') and d0<arb('37.115') and B0<arb('1.426561')

    def bias(N):
        h=T/N
        B=B0/h.sqrt()+d0*h**arb('3.5')
        e=e0*h**4;d=d0*h**arb('3.5');z=1+C*e
        U=B**2*C*(z**N-1)+2*d*B*C*z**(N-1)+d*d*C*z**(N-2)
        return 4*U/(r*sigma)

    trials={str(N):bias(N) for N in [50,100,200,500]}
    N=500;rho=bias(N)
    assert rho<arb('.000343016') and trials['200']>arb('.002')
    hmin=arb('.001')/N;hmax=T/N
    Z=1+C*e0*hmax**4;CP=C*Z**N
    Bmax=B0/hmin.sqrt()+d0*hmax**arb('3.5')
    DP=C*(Bmax*Z**N+d0*hmax**arb('3.5')*Z**(N-1))
    eta_inf,eta_row=arb('1e-11'),arb('1e-16')
    eta_d=(eta_inf*eta_row).sqrt();Znum=1+CP*eta_row
    Unum=DP**2/CP*(Znum**N-1)+2*eta_d*DP*Znum**(N-1)+eta_d**2*CP*Znum**(N-2)
    rho_num=4*Unum/(r*sigma)
    window=json.loads((ROOT/'evidence/current-quartic-positive-time-window.json').read_bytes())
    assert window['status']=='tail_contract_passed_signed_window_and_full_residual_unresolved'
    win=window['result']['constants']
    assembly=json.loads((ROOT/'evidence/current-quartic-column-energy-dual-assembly.json').read_bytes())
    assert assembly['status']=='column_bound_arithmetic_passed_full_residual_unresolved'
    col=assembly['result']['constants']
    previous=json.loads((ROOT/'evidence/current-coupled-truncation-energy-error.json').read_bytes())['runs'][-1]
    assert previous['status']=='passed'
    pars=previous['result']['parameters']
    pref=1/(dec(pars['B2'])*dec(pars['beta']))
    alpha=dec(col['sharper_exact_one_dimensional_coefficient_envelope'])
    CL=dec(col['CLt']);W=dec(win['truncated_time_weight_upper_expression'])
    tails=dec(win['cap_transport_time_and_angular_debit'])
    weighted_other=arb('9e-32');target=arb('1e-30')
    bias_charge=alpha*CL*pref*W*rho
    conditional=tails+alpha*CL*(pref*W*(rho+rho_num)+weighted_other)
    assert conditional<target and rho+rho_num<arb('.002')
    constants={'coercivity':coercivity,'complex_semigroup_growth':omega,
        'exact_evolution_weighted_jet_bound':C,'weighted_short_step_smoothing_constant':B0,
        'weighted_short_step_error_smoothing_constant':d0,
        'weighted_short_step_error_Schur_constant':e0,
        'weighted_short_step_error_pointwise_constant':f0,
        'uniform_analytic_density_bias':rho,
        'analytic_bias_original_one_column_residual_charge':bias_charge,
        'tails_plus_analytic_bias_one_column_debit':tails+bias_charge,
        'conditional_step_pointwise_error':eta_inf,'conditional_step_Schur_error':eta_row,
        'conditional_propagated_numerical_density_error':rho_num,
        'conditional_total_density_error':rho+rho_num,
        'conditional_other_weighted_error':weighted_other,
        'conditional_one_column_residual_upper':conditional,
        'conditional_one_column_screen_slack':target-conditional}
    return {'N':N,'h_min':enc(hmin),'h_max':enc(hmax),
        'constants':{k:enc(v) for k,v in constants.items()},
        'bias_trials':{k:enc(v) for k,v in trials.items()},
        'short_step_error_pointwise_constants':{k:enc(v) for k,v in zip(keys,point)},
        'short_step_error_Schur_constants':{k:enc(v) for k,v in zip(keys,row)},
        'full_thermal_and_angular_scope_preserved':True,
        'analytic_composition_bias_certified':True,
        'actual_signed_kernel':None,'actual_step_pointwise_error':None,
        'actual_step_Schur_error':None,'actual_other_weighted_error':None,
        'actual_signed_window':None,'full_residual_bound':None,
        'coupled_admission_or_existence_certified':False}

if __name__=='__main__':
    path=ROOT/'evidence'/(STEM+'.json')
    assert not path.exists(),'Preserve immutable attempts.'
    before={p:digest(ROOT/p) for p in INPUTS}
    receipt={'started_utc':datetime.now(timezone.utc).isoformat(),
        'input_sha256_before':before,'python':sys.executable,
        'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'source_utf8':Path(__file__).read_text(encoding='utf-8'),
        'primary_theorem_source':'https://arxiv.org/html/1707.00923',
        'theorem_hypotheses_checked':'Theorems 0.1 and 3.1; dense common domain, closed sectorial forms, holomorphic coefficients, uniform semigroup growth.'}
    try:
        receipt['result']=execute()
        receipt['status']='analytic_composition_bias_passed_numerical_application_and_full_residual_unresolved'
    except BaseException as exc:
        receipt['status']='failed';receipt['failure']=repr(exc)
        raise
    finally:
        receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert before==receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'],'constants':{
        k:v['ball'] for k,v in receipt['result']['constants'].items()}},indent=2))
