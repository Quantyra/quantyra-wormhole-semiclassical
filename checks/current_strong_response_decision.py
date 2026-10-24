"""Arithmetic support for the complete strong-response attempt; not a kernel proof."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
STEM='current-strong-response-decision'
OWN=['notes/'+STEM+'.md','checks/current_strong_response_decision.py','evidence/'+STEM+'-review.json']
checks=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(n,b):
    checks.append({'name':n,'passed':bool(b)})
    if not b:raise AssertionError(n)
def calc():
    prev=json.loads((ROOT/'evidence/current-uniform-principal-coefficient.json').read_bytes())['runs'][-1]
    flat=json.loads((ROOT/'evidence/current-weighted-third-comparison.json').read_bytes())['runs'][-1]
    check('prior_component_runs_passed',prev['status']=='passed' and flat['status']=='passed')
    check('prior_component_inputs_unchanged',all(digest(f)==h for f,h in prev['input_sha256_after'].items()))
    pars,_,_,_=old.scalar.prev.york.actual_parameters()
    beta=pars['beta'];q=[arb(1),arb(0),pars['u']/beta,(-4*beta-2*pars['u']-pars['tau'])/beta,(3*beta+pars['u']+pars['eps'])/beta,-2*pars['eps']/beta,pars['eps']/beta]
    a=arb('.01');A=arb('.001')
    poly=[sum((q[i]*q[k-i] for i in range(max(0,k-6),min(6,k)+1)),arb(0)) for k in range(13)]
    R=[sum((abs(poly[k])*k**j*a**k for k in range(1,13)),arb(0)) for j in range(4)]
    B=[A*3**j*a**3/6 for j in range(4)];e=(2*B[0]).exp()
    E=[e-1,2*e*B[1],e*(2*B[2]+4*B[1]**2),e*(2*B[3]+12*B[1]*B[2]+8*B[1]**3)]
    F=[e]+E[1:]
    AA=[sum((comb(j,k)*R[k]*F[j-k] for k in range(j+1)),arb(0))+E[j] for j in range(4)]
    weights=[100,51,12,1]
    mult=sum((weights[j]*AA[j] for j in range(4)),arb(0))
    Cdelta=arb('1.295e10')*mult/100
    check('replacement_complete_comparison_guard',Cdelta>0 and Cdelta<407000)
    z=arb('.002','.0001');d=sum((q[k]*z**k for k in range(1,7)),arb(0));rho=2*d+d*d
    check('whole_actual_patch_nonzero_negative_curvature_coefficient',rho<0 and -2*rho>arb('1.07e-7') and -2*rho<arb('1.49e-7'))
    # Independent coefficient assembly: product rule for 2d+d^2, never the
    # squared-Q polynomial. Formal exponential jets evaluated by convolution.
    dcap=[sum((abs(q[k])*k**j*a**k for k in range(1,7)),arb(0)) for j in range(4)]
    Rind=[2*dcap[j]+sum((comb(j,k)*dcap[k]*dcap[j-k] for k in range(j+1)),arb(0)) for j in range(4)]
    jet=[arb(0),2*B[1],B[2],B[3]/3]
    expjet=[arb(1),arb(0),arb(0),arb(0)];power=expjet.copy()
    factorial=1
    for n in range(1,4):
        factorial*=n
        power=[sum((power[k]*jet[j-k] for k in range(j+1)),arb(0)) for j in range(4)]
        expjet=[expjet[j]+power[j]/factorial for j in range(4)]
    Fi=[e*expjet[j]*[1,1,2,6][j] for j in range(4)]
    Ei=[Fi[0]-1]+Fi[1:]
    Ai=[sum((comb(j,k)*Rind[k]*Fi[j-k] for k in range(j+1)),arb(0))+Ei[j] for j in range(4)]
    Cind=arb('1.295e10')*sum((weights[j]*Ai[j] for j in range(4)),arb(0))/100
    check('independent_full_exponential_jet_guard',Cind>0 and Cind<407000)
    check('independent_bound_not_smaller',Cind.overlaps(Cdelta) or Cind>Cdelta)
    cells=[]
    for k in range(64):
        left=arb('.0019')+arb('.0002')*k/64;right=arb('.0019')+arb('.0002')*(k+1)/64
        zz=arb((left+right)/2,(right-left)/2)
        dd=sum((q[j]*zz**j for j in range(1,7)),arb(0))
        rr=(1+dd)**2-1
        if not (-2*rr>arb('1.07e-7') and -2*rr<arb('1.49e-7')):raise AssertionError('packet_subdivision')
        cells.append(old.encoded(-2*rr))
    check('independent_complete_patch_subdivision',len(cells)==64)
    reserve=arb('.1599164')-arb('.159916395999')
    rest=arb('.1599164')-reserve-arb('1e-11')*(arb('1.295e10')+407000)
    check('full_target_ledger_exact',rest.overlaps(arb('.030412325999')))
    check('no_previous_principal_debit_counted_twice',rest>arb('.03041232599899'))
    cap=rest/arb('1e-11');ab=arb('5e7')/(arb('1.295e10')+407000)
    check('positive_but_unknown_remainder_budget',cap>arb('3.0412325998e9') and cap<arb('3.0412326000e9'))
    check('original_comparison_product_threshold',ab>arb('.00386') and ab<arb('.003862'))
    check('original_preservation_allocation',arb('5e-9')+arb('1e-8')*arb('.25')+arb('5e-17')*arb('5e7')==arb('1e-8') or (arb('5e-9')+arb('1e-8')*arb('.25')+arb('5e-17')*arb('5e7')).overlaps(arb('1e-8')))
    return dict(outcome='complete_strong_estimate_unresolved',complete_actual_bound_executed=False,
                ordinary_actual_domain_cutoff_inference_failure=True,packet_uniform_coefficient=old.encoded(-2*rho),
                packet_bilinear_lower_bound_requires_fixed_low_datum=True,finite_packet_onset_certified=False,
                multiplier=old.encoded(mult),complete_comparison_defect_upper=old.encoded(Cdelta),
                independent_comparison_defect_upper=old.encoded(Cind),independent_patch_cells=cells,
                new_comparison_guard='407000',comparison_guard='1.295e10',reserve=old.encoded(reserve),
                replacement_feedback_cost='0.00000407',remaining_budget=old.encoded(rest),
                best_case_C_rem_ceiling=old.encoded(cap),equal_split_eta=old.encoded(rest/2),equal_split_C_rem=old.encoded(cap/2),
                beta_alpha_squared_ceiling_comparison_only=old.encoded(ab),beta_ceiling_if_alpha_point001=old.encoded(ab/arb('1e-6')),
                eta_seed_known=False,C_rem_known=False,original_norm_comparisons_established=False,
                strong_common_limit_established=False,global_admission=False,actual_controlled_solution=False,
                exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
                audit_status='closed',goal_status='active',physical_domain_changed=False,normalization_changed=False)
def main():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    prior=json.loads((ROOT/'evidence/current-uniform-principal-coefficient.json').read_bytes())['runs'][-1]
    paths.update(prior['input_sha256_after'])
    before={p:digest(p) for p in sorted(paths)}
    out=ROOT/'evidence'/f'{STEM}.json'
    archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running','input_sha256_before':before,
         'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},'precision_bits':ctx.prec,'python':sys.version,'invocation':sys.argv}
    archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        run['result']=calc();run['status']='passed'
    except Exception:run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=checks;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    r=run.get('result',{})
    print(json.dumps({k:v for k,v in {'status':run['status'],'checks':len(checks),'inputs':len(before),'exception':run.get('exception'),'comparison_defect':r.get('complete_comparison_defect_upper'),'independent':r.get('independent_comparison_defect_upper'),'remaining':r.get('remaining_budget'),'full_estimate_closed':False}.items()},indent=2))
    if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
