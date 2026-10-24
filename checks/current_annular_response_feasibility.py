"""Arithmetic for the ordinary all-order estimator test; no full strong certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
STEM='current-annular-response-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_annular_response_feasibility.py','evidence/'+STEM+'-review.json']
checks=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(n,b):
    checks.append({'name':n,'passed':bool(b)})
    if not b:raise AssertionError(n)
def calc():
    pars,_,_,_=old.scalar.prev.york.actual_parameters();beta=pars['beta']
    q=[arb(1),arb(0),pars['u']/beta,(-4*beta-2*pars['u']-pars['tau'])/beta,
       (3*beta+pars['u']+pars['eps'])/beta,-2*pars['eps']/beta,pars['eps']/beta]
    H=arb('.005006');c=arb('4.004');a=arb('.01');R=arb('.000251');sig=arb('4.005')
    B=arb('.001')*(c*H)**3/6;ee=B.exp();jm=arb('.999').sqrt();jp=arb('1.001').sqrt()
    d2=abs(q[2])*c**2;d3=sum((abs(q[k])*c**k*H**(k-3) for k in range(3,7)),arb(0))
    e2=2*abs(q[2])*c;e3=sum((k*abs(q[k])*c**(k-1)*H**(k-3) for k in range(3,7)),arb(0))
    f2=ee*d2/(1+jm);f3=ee*(d3/(1+jm)+arb('.001')*c**3/6)
    fp2=ee*arb('1.001')*e2/(2*jm)
    fp3=ee*arb('1.001')*(e3/(2*jm)+jp*arb('.001')*c**2/2)
    check('actual_annular_mass_guards',f2<arb('1e-8') and f3<140 and fp2<arb('1e-8') and fp3<110)
    check('all_output_and_input_supports',2*H*arb('.999')>a and arb('3.5')*H*arb('1.001')<arb('.02') and c*H<arb('.04'))
    check('cutoff_derivative_guard',arb(factorial(17))*4/(factorial(8)**2*2**16)<14)
    # Coefficientwise comparison, valid at every h>=0.
    check('whole_annular_Fourier_radius_polynomial',
        arb('7.5')*140*2070<2000**2 and
        arb('7.5')*(140*arb('1.5e-7')+2070*arb('1e-8'))<arb('.004') and
        arb('7.5')*arb('1e-8')*arb('1.5e-7')<arb('1e-12'))
    r0=2000*H**3+arb('1e-6')*H**2
    check('uniform_all_annuli_series_radius',r0<R)
    def Ak(k):
        if k==2:return ((arb(79)/3)/arb(2).log()+10)/arb.pi()
        bb=arb(comb(k+3,4))/(k+2);cc=arb(comb(k+5,6))/(k+4)
        return ((arb(8)/(k-2)+2*(k+1)+4*bb+8*cc)/arb(2).log()+8*bb)/arb.pi()
    def Dn(n,degree):
        out=arb(0)
        for mm in range(1,degree+1):
            for j in range(mm+1):
                k=n+mm-j
                if k<max(2,mm):continue
                if degree==3:coef=[0,1,3*(k-1),(k-1)*(k-2)][mm]
                else:coef=[0,1,k-1][mm]
                out+=coef*comb(mm,j)*Ak(k)
        return out
    N=16;ratio15=R*(arb(N+10)/(N+9))**15;ratio16=R*(arb(N+10)/(N+9))**16
    tail=3*arb(N+9)**15*R**(N+1)/(1-ratio15)
    dtail=3*arb(N+9)**16*R**N/(1-ratio16)
    check('infinite_order_tail',ratio16<arb('.001') and tail<arb('1e-38') and dtail<arb('1e-33'))
    series={}
    for degree in [2,3]:
        for ell in range(4,9):
            val=arb(0);der=arb(0)
            for n in range(1,N+1):
                slots=n+degree-1
                coeff=Dn(n,degree)*arb(slots)**(ell-1)*(1+arb(slots).log()/arb(2).log())/16
                val+=coeff*R**n;der+=n*coeff*R**(n-1)
            series[f'{degree}_{ell}']={'upper':old.encoded(val+tail),'derivative_upper':old.encoded(der+dtail)}
    check('all_order_third_response_guards',old.decode(series['3_8']['upper'])<arb('85.321') and old.decode(series['3_8']['derivative_upper'])<343246)
    check('entire_chain_rule_first_coefficient',Dn(1,3).overlaps(7*Ak(2)+12*Ak(3)+6*Ak(4)) and Dn(1,2).overlaps(3*Ak(2)+2*Ak(3)))
    qcomp=4+48*(8*arb(3).log()*arb(2).log()+4*arb(3).log()**2)
    check('actual_log_coordinate_bound',48*(1+arb(3).log()/arb(2).log())<125 and qcomp<529 and 19*529*8<512**2)
    check('actual_top_coefficient_floor',arb('.999')**8*arb('1.001')**(sig-8)>arb('.98'))
    K=pars['C']*beta**2/(15*arb.pi()*pars['B2']);p0=2*pars['C']*beta**2/pars['B2'];Cq=arb('2.328')
    Eabs=2*2*100*arb(2)**(8-sig)*512*arb('.98')*(30*arb.pi()*Cq)*16*Dn(1,2)*r0
    check('specified_positive_estimator_fails_before_other_terms',Eabs>arb('3.2e8'))
    check('finite_guard_refinement_alone_insufficient_in_this_formula',Eabs/(512*100*2)>3000)
    check('same_physical_normalization', (p0*Cq/K).overlaps(30*arb.pi()*Cq))
    poly=[sum((q[i]*q[k-i] for i in range(max(0,k-6),min(6,k)+1)),arb(0)) for k in range(13)]
    rr=[sum((abs(poly[k])*k**j*a**k for k in range(1,13)),arb(0)) for j in range(4)]
    mult=sum((w*v for w,v in zip([100,51,12,1],rr)),arb(0));eta=Cq*mult
    dd=[sum((abs(q[k])*k**j*a**k for k in range(2,7)),arb(0)) for j in range(4)]
    ri=[2*dd[j]+sum((comb(j,k)*dd[k]*dd[j-k] for k in range(j+1)),arb(0)) for j in range(4)]
    independent=Cq*sum((w*v for w,v in zip([100,51,12,1],ri)),arb(0))
    check('actual_signed_principal_guard',eta<arb('.00731') and independent<arb('.00731'))
    rest=arb('.1599164')-arb('4.001e-9')-arb('1e-11')*359734306-arb('.00731')
    check('original_target_repartition',rest.overlaps(arb('.149009052939')))
    return dict(outcome='specified_annular_estimator_failure_with_quantified_principal_alternative',
        complete_strong_estimate_status='unresolved',actual_complete_bound_executed=False,
        coefficient_guards={k:old.encoded(v) for k,v in [('f2',f2),('f3',f3),('fp2',fp2),('fp3',fp3)]},
        all_annuli_Wiener_radius=old.encoded(r0),series_radius_guard='.000251',series=series,
        all_orders_included=True,angular_spectrum_truncated=False,
        order_tail_upper=old.encoded(tail),order_derivative_tail_upper=old.encoded(dtail),
        specified_estimator_positive_allowance_lower=old.encoded(Eabs),
        idealized_guard_sensitivity=old.encoded(Eabs/(512*100*2)),
        actual_response_lower_bound_claimed=False,construction_obstruction=False,
        principal_multiplier=old.encoded(mult),principal_norm_upper=old.encoded(eta),
        independent_principal_upper=old.encoded(independent),principal_guard='.00731',
        remaining_budget=old.encoded(rest),equal_split_seed_remaining=old.encoded(rest/2),
        equal_split_C_rem=old.encoded(rest/arb('2e-11')),best_case_C_rem=old.encoded(rest/arb('1e-11')),
        full_principal_subprincipal_remainder_known=False,seed_remaining_known=False,C_rem_known=False,
        common_strong_regulator_limit_established=False,original_norm_comparisons_established=False,
        local_contact_targets_established=False,actual_controlled_solution=False,
        nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,
        all_sector_stability=False,physical_domain_changed=False,normalization_changed=False,
        selected_projection_screen_retained=True,projection_reevaluated=False,
        audit_status='closed',goal_status='active',next_symmetric_remainder_calculation_executed=False)
def main():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
    prev=json.loads((ROOT/'evidence/current-complete-response-series-attempt.json').read_bytes())['runs'][-1]
    paths.update(prev['input_sha256_after']);before={p:digest(p) for p in sorted(paths)}
    out=ROOT/'evidence'/f'{STEM}.json';archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running','input_sha256_before':before,
         'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},'precision_bits':ctx.prec,
         'python':sys.version,'invocation':sys.argv};archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        check('previous_sources_preserved',all(digest(p)==h for p,h in prev['input_sha256_after'].items()))
        run['result']=calc();run['status']='passed'
    except Exception:run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=checks;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    result=run.get('result',{})
    print(json.dumps({'status':run['status'],'checks':len(checks),'inputs':len(before),'exception':run.get('exception'),
        'estimator_term':result.get('specified_estimator_positive_allowance_lower'),
        'principal':result.get('principal_norm_upper'),'remaining':result.get('remaining_budget')},indent=2))
    return 0 if run['status']=='passed' else 1
if __name__=='__main__':raise SystemExit(main())
