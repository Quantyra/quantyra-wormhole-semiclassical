"""Arithmetic support for the ordinary scalar/geodesic comparison.
No inverse or complete quantum application is implemented by this checker.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json

ROOT=Path(__file__).resolve().parents[1]
STEM='current-global-scalar-log-comparison'
INPUTS=['notes/'+STEM+'.md','checks/current_global_scalar_log_comparison.py']
INPUTS+=['notes/'+n+'.md' for n in [
 'stationary-curved-quantum-tt','uniform-logarithmic-frequency-comparison',
 'stationary-paired-dirac-phase','current-seed-full-hessian-realization',
 'current-geodesic-two-sided-remainder','current-geodesic-flow-multiplier',
 'current-outer-input-source-test','current-quartic-ritz-rank-test']]
INPUTS+=['evidence/'+n+'.json' for n in [
 'current-geodesic-two-sided-remainder','current-geodesic-flow-multiplier',
 'current-structured-source-graph-test','current-quartic-ritz-rank-test-final',
 'current-outer-input-source-test-final','current-coupled-mass-adjoint-runtime']]
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
def digest(p): return sha256(p.read_bytes()).hexdigest()

def execute():
    from pathlib import Path
    from hashlib import sha256
    from fractions import Fraction as Q
    from math import factorial as fac
    from datetime import datetime,timezone
    import json,sys
    ROOT=Path(__file__).resolve().parents[1]
    def read(n): return json.loads((ROOT/('evidence/'+n+'.json')).read_bytes())
    runtime=read('current-coupled-mass-adjoint-runtime')
    sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
    from flint import arb,ctx
    ctx.prec=512
    def digest(p): return sha256(Path(p).read_bytes()).hexdigest()
    assert digest(sys.executable)==runtime['python_sha256']
    for pkg in runtime['packages'].values():
     for n,want in pkg['files'].items(): assert digest(Path(pkg['root'])/n)==want
    protected=json.loads(Path(r'C:\Users\dfred\AppData\Local\Temp\s169-weak-response-execution-baseline.json').read_bytes())
    assert all(digest(p)==h for p,h in protected.items())
    def dec(v):
     d=v['decimal_enclosure']
     return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
    def enc(v):
     m,r,e=v.mid_rad_10exp(60)
     return {'ball':v.str(40),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
    def ball(v): return arb(Q(v).numerator)/Q(v).denominator
    geo=read('current-geodesic-two-sided-remainder')['runs'][-1]
    flow=read('current-geodesic-flow-multiplier')['runs'][-1]
    source=read('current-structured-source-graph-test')
    rank=read('current-quartic-ritz-rank-test-final')
    outer=read('current-outer-input-source-test-final')
    assert geo['status']==flow['status']==rank['status']=='passed'
    gp=next(x for x in geo['numerics'] if x['kind']=='actual_finite_seed_two_sided_parameters')
    fq=next(x for x in flow['numerics'] if x['kind']=='actual_flow_quadrature_form_budget')
    rc,oc=rank['result']['constants'],outer['result']['constants']
    eps,qt=dec(rc['epsilon']),dec(rc['q_t'])
    eq=dec(source['constants']['actual_reference_comparison'])
    assert dec(gp['F_eff']).abs_upper()<100
    R=arb(1)/400
    t0=R**2/8
    v=[arb('.026')*fac(k)*128**k for k in range(4)]
    sig=[v[1]/(1-v[0]),v[2]/(1-v[0])+v[1]**2/(1-v[0])**2,
     v[3]/(1-v[0])+3*v[1]*v[2]/(1-v[0])**2+2*v[1]**3/(1-v[0])**3]
    assert all(x<b for x,b in zip(sig,[10,10000,10000000]))
    assert 4*(4*t0+R**2/2).exp()<5 and R**2<2
    near=6*arb('.2')/arb.pi()**2*(16+128)
    rows=36*arb.pi()+arb('6.4')*64**2/arb.pi()+arb('409.6')/arb.pi()+arb('51.2')/arb.pi()
    assert near<100 and rows<10000 and rows/R**4<arb('1e16')
    coeff={6:462,7:-1980,8:3465,9:-3080,10:1386,11:-252}
    assert sum(coeff.values())==1
    for k in range(1,6): assert sum(c*fac(n)//fac(n-k) for n,c in coeff.items())==0
    B=[Q(1)]+[sum((abs(c)*Q(fac(n),fac(n-k)) for n,c in coeff.items()),Q(0)) for k in range(1,6)]
    H=[Q(1)]
    for k in range(1,6):
     H.append(sum((Q(fac(k),fac(k-2*j)*fac(j))*B[k-j]*Q(4,3)**(k-2*j)*Q(1,3)**j for j in range(k//2+1)),Q(0)))
     assert Q(k**k,fac(k))<=2**k
    q=[ball(2**k*(4*H[k]+(4*k*H[k-1] if k else 0)+(k*(k-1)*H[k-2] if k>=2 else 0)))*R**(2-k) for k in range(6)]
    a0,a1,a2,a3=map(arb,[2,40,40800,42416000])
    m2=4*a0**2*q[4]+8*a0*a1*q[3]+4*a0*a2*q[2]
    lip2=2*(4*a0**2*q[5]+16*a0*a1*q[4]+(8*a1**2+12*a0*a2)*q[3]+(4*a1*a2+4*a0*a3)*q[2])
    h1=3*(lip2+2*m2/R)
    assert arb(5).sqrt()<3
    assert Q(4,15)*Q(1,2)*Q(3,4)*Q(5,1024)==Q(1,2048)
    scalar_probe=h1*arb.pi().sqrt()/3072
    geo_far=136/R**3+10/R
    assert geo_far<arb('1e16')
    geo_probe=arb('1e16')*2*arb.pi()*225/4096*arb(2)/3
    anchor=400*arb.pi()*(5*R**2/2+arb('1e16')*R**3/3)
    assert anchor<arb('1e13')
    contact=4*(scalar_probe+geo_probe+arb('1e16')*4*R**2+anchor+100000)
    assert contact<arb('1e19')
    split_eps,K=arb('4e-6')/arb.pi(),arb('1e8')
    b2=split_eps/(4*K)
    assert b2<R**2
    gradient=K*arb.pi()*(arb(1)/4+2*(R/b2.sqrt()).log())+split_eps*arb.pi()/6
    assert gradient+arb('1e19')+100000<arb('1e20')
    assert 2*arb('1e16')<arb('1e18')
    fourth=200*arb.pi()*arb('1e-6')
    lower=(arb('1e20')*dec(gp['j1'])+arb('1e18')*dec(gp['j0'])+2*arb('1e-6')*500**4*dec(gp['j0_outer']))/(arb('.99')*dec(gp['F_lambda']))
    assert fourth>lower
    esc=fourth/(1-eps)
    esq=esc+dec(fq['eta_quad'])/(1-eps)
    total=eq+esq
    inv=1/(1-total)
    d=esq*inv
    assert esc<arb('.001') and total<arb('.004') and d<arb('.001')
    w=qt+dec(oc['complete_outer_source_Et_to_dual_bound'])
    ceiling=dec(oc['remaining_inner_application_allowance_NOT_attained'])
    weak=(1-eq)*ceiling
    m=1
    while not d**m*w<weak/2:
     m+=1
     assert m<1000
    assert d**(m-1)*w>weak/2
    exact_residual=d**m*w
    exact_application=exact_residual/(1-eq)
    target=arb('1e-32')
    dhat=d+target
    am=sum((dhat**j for j in range(m)),arb(0))
    dm=sum((j*dhat**(j-1) for j in range(1,m)),arb(0))
    am0=sum((d**j for j in range(m)),arb(0))
    cm,dmul,bm=am*(w+target),inv*dm*(w+target),inv*am0
    conditional_s=exact_application+(cm+dmul+bm)*target+dec(oc['attained_outer_S_omission_error'])
    assert conditional_s<arb('5e-32')
    conditional_full=dec(rc['alternative_exact_application_residual_upper'])+(1+qt)*(dec(rc['source_application_error_multiplier'])*conditional_s+dec(rc['K_application_error_multiplier'])*arb('5e-32'))
    assert conditional_full<arb('1e-30')
    # Actual-domain obstruction to a fixed finite-output UNIFORM inverse.
    # This does not lower-bound its restriction to the structured B range.
    finite_inverse_error_lower=1/(1+qt+total)
    assert finite_inverse_error_lower>arb('.74')
    finite_inverse_ratio=finite_inverse_error_lower/target
    assert finite_inverse_ratio>arb('7.4e31')
    kb=dec(source['constants']['complete_B_Et_to_H_and_S_operator_graph_constant'])
    spectral_tail_prefactor=kb/(1-eq).sqrt()
    spectral_cutoff_sufficient=(spectral_tail_prefactor/arb('5e-32'))**2
    vals=dict(near_kernel_coefficient=near,return_row_coefficient=rows/R**4,
     scalar_probe_M2=m2,scalar_probe_Lipschitz=lip2,scalar_probe_H1=h1,
     complete_scalar_probe_coefficient=scalar_probe,full_geo_far_remainder=geo_far,
     complete_geo_probe_coefficient=geo_probe,anchor_coefficient=anchor,
     second_contact_coefficient=contact,near_gradient_coefficient=gradient,
     fourth_energy_ratio=fourth,lower_energy_ratio=lower,scalar_geo_form_error=esc,
     scalar_quantized_geo_form_error=esq,actual_scalar_reference_comparison=total,
     exact_reference_inverse_norm=inv,exact_complement_correction_norm=d,
     complete_inner_source_norm=w,exact_reference_one_step_residual_upper=d*w,
     one_step_upper_over_weak_screen=d*w/weak,remaining_inner_application_ceiling_NOT_attained=ceiling,
     weak_screen_NOT_attained=weak,exact_application_only_m_step_residual=exact_residual,
     exact_application_only_m_step_map_error=exact_application,
     proposed_common_C_D_B_error_NOT_attained=target,C_error_multiplier=cm,
     D_error_multiplier=dmul,B_error_multiplier=bm,conditional_S_error_NOT_attained=conditional_s,
     conditional_complete_residual_NOT_attained=conditional_full,
     finite_output_uniform_inverse_error_lower=finite_inverse_error_lower,
     finite_output_uniform_inverse_error_over_target=finite_inverse_ratio,
     structured_spectral_tail_prefactor=spectral_tail_prefactor,
     sufficient_spectral_cutoff_NOT_implemented=spectral_cutoff_sufficient)
    inputs=['notes/current-global-scalar-log-comparison.md']
    inputs+=['evidence/'+n+'.json' for n in ['current-geodesic-two-sided-remainder','current-geodesic-flow-multiplier','current-structured-source-graph-test','current-quartic-ritz-rank-test-final','current-outer-input-source-test-final','current-coupled-mass-adjoint-runtime']]
    return dict(status='arithmetic_passed_ordinary_argument_required',precision_bits=512,
     timestamp_utc=datetime.now(timezone.utc).isoformat(),input_sha256={n:digest(ROOT/n) for n in inputs},
     constants={k:enc(v) for k,v in vals.items()},minimum_exact_application_correction_count=m,
     protected_paths_unchanged=len(protected),numerical_inverse_implemented=False,
     actual_full_residual_computed=False,physical_margin_closed=False,goal_complete=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='evidence/'+STEM+'-final.json')
    args=parser.parse_args()
    out=ROOT/args.output
    if out.exists():
        raise SystemExit('Refusing to replace an immutable receipt.')
    before={p:digest(ROOT/p) for p in INPUTS}
    protected=json.loads(BASELINE.read_bytes())
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),
        precision_bits=512,input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=(ROOT/INPUTS[1]).read_text(encoding='utf-8'))
    try:
        assert all(digest(Path(p))==h for p,h in protected.items())
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        assert before==run['input_sha256_after']
        assert all(digest(Path(p))==h for p,h in protected.items())
        run['protected_paths_unchanged']=len(protected)
        run['status']='form_comparison_checked_numerical_application_unresolved'
    except Exception as exc:
        run['status']='failed'
        run['exception']=type(exc).__name__+': '+str(exc)
    run['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_text(json.dumps(run,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'status':run['status'],'receipt':str(out),
        'exception':run.get('exception'),
        'constants':run.get('result',{}).get('constants',{})},indent=2))
    if run['status']=='failed': raise SystemExit(1)
