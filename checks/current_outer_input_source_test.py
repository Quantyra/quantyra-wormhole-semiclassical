"""Complete outer-input source omission; no global inverse implemented."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import argparse

ROOT=Path(__file__).resolve().parents[1]
STEM='current-outer-input-source-test'
RUNTIME=ROOT/'evidence/current-coupled-mass-adjoint-runtime.json'
runtime=json.loads(RUNTIME.read_bytes())
sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb,ctx

INPUTS=['notes/'+STEM+'.md','checks/current_outer_input_source_test.py',
        'notes/current-structured-source-graph-test.md',
        'notes/current-seed-full-hessian-realization.md',
        'notes/current-seed-boundary-phase-and-trace-defect.md',
        'notes/strong-effective-tt-source.md',
        'notes/current-weak-complement-inverse-attempt.md',
        'evidence/current-geodesic-two-sided-remainder.json',
        'evidence/current-structured-source-graph-test.json',
        'evidence/current-quartic-ritz-rank-test-final.json',
        'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')

def digest(p):return sha256(p.read_bytes()).hexdigest()
def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))
def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items():out[k]=out.get(k,arb(0))+v
    return {k:v for k,v in out.items() if not v.is_zero()}
def scale(p,a):return {k:a*v for k,v in p.items()}
def diff(p):return {k-1:k*v for k,v in p.items() if k}
def shift(p,n):return {k+n:v for k,v in p.items()}
def mul(p,q):
    out={}
    for i,a in p.items():
        for j,b in q.items():out[i+j]=out.get(i+j,arb(0))+a*b
    return out
def coeffsum(p):return sum((v.abs_upper() for v in p.values()),arb(0))

def execute():
    ctx.prec=512
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for p,h in pkg['files'].items():assert digest(Path(pkg['root'])/p)==h
    geo=json.loads((ROOT/INPUTS[7]).read_bytes())['runs'][-1]
    source=json.loads((ROOT/INPUTS[8]).read_bytes())
    rank=json.loads((ROOT/INPUTS[9]).read_bytes())
    assert geo['status']==rank['status']=='passed'
    assert source['status']=='source_graph_bound_checked_complete_application_unresolved'
    sc=source['constants'];rc=rank['result']['constants']
    ep=next(x for x in geo['numerics'] if x['kind']=='actual_seed_correlated_endpoint')
    gp=next(x for x in geo['numerics'] if x['kind']=='actual_finite_seed_two_sided_parameters')
    Q={k:dec(v) for k,v in zip([0,2,3,4,5,6],ep['Q_coefficients'])}
    beta=Q[0];C=dec(gp['C']);B2=dec(gp['B2']);lam=dec(gp['lambda_value'])
    assert beta.overlaps(lam*B2/C)
    # The endpoint lies below 1 by the accepted correlated h>0 identity,
    # not by subtracting two independently rounded near-unit endpoints.
    assert dec(ep['h'])>0 and ep['exact_Q_bolt']=='0'
    x=lam*B2;ds=C*Q[6];q=dec(gp['q'])
    e=lam/(60*arb.pi());d=dec(sc['d']);Cv=dec(sc['C_v'])
    assert x>arb('8e9') and ds.abs_upper()<arb('8.6e-19')
    alpha=arb('1e-9')
    qmin=beta-sum((v.abs_upper()*alpha**k for k,v in Q.items() if k),arb(0))
    iq=sum((v.abs_upper()/(k+3) for k,v in Q.items()),arb(0))
    assert qmin>beta/2
    keta=3*(beta*iq/(5*qmin*qmin)).sqrt()

    w={5:4*ds/x,6:-6*ds/x}
    tau={4:(1+3*x+ds)/x,5:-5*ds/x,6:arb(9)/2*ds/x}
    wc={3:-C*Q[3]/(2*x),4:-C*Q[4]/x,5:arb(10)/3*ds/x,6:-arb(5)/2*ds/x}
    guards={}
    for name,p,limits in [('j',tau,['3.6','20','108']),('Weyl_amplitude',wc,['5.1','40','153'])]:
        tests=[coeffsum(p),coeffsum(shift(diff(p),1)),coeffsum(shift(diff(diff(p)),2))]
        assert all(t<arb(limit) for t,limit in zip(tests,limits))
        guards[name]=tests
    assert coeffsum(tau)<arb('3.01')
    assert coeffsum(w)<arb('1.1e-27') and 4000*ds.abs_upper()/x<arb('4.3e-25')
    drift=add(diff(Q),scale(shift(Q,-1),-2))
    lapw=scale(shift(add(mul(Q,diff(diff(w))),mul(drift,diff(w))),2),1/beta)
    v=q*q/(60*arb.pi()*lam*B2**2)
    # Exact original vacuum identity removes the constant; no numerical reset.
    trace=add(scale(w,Cv),scale(mul(w,w),-11*e/72),
              scale(mul(wc,wc),28*e),scale(mul(tau,tau),arb(22)/3*e),
              {4:-10*v},scale(lapw,-d))
    assert min(trace)>=4
    kt=coeffsum(trace)
    assert kt<arb('2e-27')
    h1=1/arb('5.99').sqrt()
    pieces=dict(main=dec(sc['mixed_main'])*alpha**4,
                curvature=dec(sc['mixed_curvature'])*alpha**3,
                w=dec(sc['mixed_w'])*alpha**5,
                Hess_w=dec(sc['mixed_Hess_w'])*alpha**5,
                magnetic=dec(sc['mixed_magnetic'])*alpha**4,
                GLOBAL_return=keta*kt*h1*alpha**2*alpha.sqrt())
    kalpha=sum(pieces.values(),arb(0))
    kb=dec(sc['complete_mixed_A_to_H_constant']);ch=dec(sc['c_H'])
    eq=dec(sc['actual_reference_comparison']);qt=dec(rc['q_t'])
    a12=dec(rc['source_application_error_multiplier'])
    rho=kb*kalpha/ch
    aouter=rho/(1-eq)
    debit=(1+qt)*a12*aouter
    target=arb('5e-32')
    remaining=target-aouter
    epsilon=dec(rc['epsilon'])
    e0=(arb('.0022')+epsilon)/(1-epsilon)
    inner_norm_bound=qt/(1-eq)+aouter
    previous_K_caller_bound=qt/(1-e0)
    assert inner_norm_bound<previous_K_caller_bound
    assert 0<aouter<target and debit<arb('1e-35')
    values=dict(alpha=alpha,original_beta=beta,outer_Q_lower=qmin,
                full_weighted_Q_integral_upper=iq,weighted_York_return_constant=keta,
                complete_trace_z4_coefficient_upper=kt,
                complete_outer_mixed_constant=kalpha,
                complete_outer_source_Et_to_dual_bound=rho,
                attained_outer_S_omission_error=aouter,
                attained_outer_FULL_residual_debit=debit,
                fraction_of_sufficient_S_application_allowance=aouter/target,
                remaining_inner_application_allowance_NOT_attained=remaining,
                inner_operator_norm_upper=inner_norm_bound,
                existing_K_application_caller_norm_upper=previous_K_caller_bound,
                caller_norm_margin=previous_K_caller_bound-inner_norm_bound,
                entire_residual_screen_fraction=debit/arb('1e-30'))
    return dict(constants={k:enc(v) for k,v in values.items()},
                complete_mixed_allocations={k:enc(v) for k,v in pieces.items()},
                curvature_weight_guards={k:[enc(v) for v in vs] for k,vs in guards.items()},
                full_trace_polynomial=[dict(power=k,coefficient=enc(v)) for k,v in sorted(trace.items())],
                outer_input_omission_certified=True,global_return_retained=True,
                full_output_domain_unchanged=True,numerical_inner_application_error=None,
                actual_K_application_error=None,actual_complete_residual=None,
                physical_admission_requirement_closed=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='evidence/'+STEM+'.json')
    args=parser.parse_args()
    out=ROOT/args.output
    assert not out.exists(),'Preserve immutable receipts.'
    protected=json.loads(BASELINE.read_bytes())
    verify=lambda:all(digest(Path(p))==h for p,h in protected.items())
    assert verify()
    before={p:digest(ROOT/p) for p in INPUTS}
    data=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
              ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
              arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
              previous_goal_turn='no progress: conditional complement accounting only')
    try:
        data['result']=execute()
        data['status']='complete_outer_input_omission_bounded_full_application_unresolved'
    except Exception as error:
        data.update(status='failed',error=repr(error))
        raise
    finally:
        data['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        data['inputs_unchanged']=before==data['input_sha256_after']
        data['protected_paths_unchanged']=verify()
        data['protected_path_count']=len(protected)
        data['finished_utc']=datetime.now(timezone.utc).isoformat()
        out.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert data['inputs_unchanged'] and data['protected_paths_unchanged']
    for k,v in data['result']['constants'].items():print(k+': '+v['ball'])
    for k,v in data['result']['complete_mixed_allocations'].items():print(k+': '+v['ball'])

if __name__=='__main__':main()
