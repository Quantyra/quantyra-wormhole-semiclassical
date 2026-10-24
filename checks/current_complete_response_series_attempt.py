"""Arithmetic support; does not certify the missing combined strong estimate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
STEM='current-complete-response-series-attempt'
OWN=['notes/'+STEM+'.md','checks/current_complete_response_series_attempt.py',
     'evidence/'+STEM+'-review.json']
checks=[]
def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name, passed):
    checks.append({'name':name,'passed':bool(passed)})
    if not passed: raise AssertionError(name)
def calc():
    pars,_,_,_=old.scalar.prev.york.actual_parameters()
    beta=pars['beta']; b2=beta**2
    flat=arb('1.295e10')*b2; delta=arb(407000)*b2
    check('physical_beta_not_rescaled',beta>0 and beta<arb(1)/6)
    check('complete_comparison_guards',flat<359723000 and delta<11306)
    total=arb(359734306)
    remainder=arb('.1599164')-arb('4.001e-9')-arb('1e-11')*total
    check('new_remaining_budget',remainder.overlaps(arb('.156319052939')))
    ratio=remainder/arb('.030412325999')
    check('quantitative_change_material',ratio>5 and ratio<arb('5.15'))
    radius=arb('1.001').sqrt()*arb('.001').exp()-1
    neg=1-arb('.999').sqrt()*(-arb('.001')).exp()
    check('both_mass_perturbation_endpoints',neg<radius and radius<arb('.001500876'))
    check('uniform_neumann_ratio',radius<arb('.0016'))
    cubic=(-4*beta-2*pars['u']-pars['tau'])/(2*beta)
    lo=cubic-arb('.001')/6; hi=cubic+arb('.001')/6
    check('actual_nonzero_cubic_on_entire_admitted_family',lo>arb('-2.000167') and hi<arb('-1.999833'))
    k=4
    constant=Fraction(8,k-2)+2*(k+1)+Fraction(4*comb(k+3,4),k+2)+Fraction(8*comb(k+5,6),k+4)
    logcoef=Fraction(8*comb(k+3,4),k+2)
    check('four_vertex_full_angular_majorant_coefficients',constant==Fraction(364,3) and logcoef==Fraction(140,3))
    check('cubic_fourth_distributional_derivative',6-(-6)==12)
    return dict(outcome='unresolved_estimate_with_executed_attempt_evidence',
        complete_actual_bound_executed=False,eta_seed_known=False,C_rem_known=False,
        global_Sobolev_Born_majorant='infinite on actual paired cubic seam; ordinary proof, not numerical infinity',
        global_Sobolev_estimator_rejected=True,construction_obstruction=False,
        actual_paired_first_order_conjugation_proved=True,uniform_resolvent_series=True,
        trace_series_identified_on_actual_graph=False,strong_common_regulator_limit_proved=False,
        beta=old.encoded(beta),beta_squared=old.encoded(b2),
        scaled_flat_upper=old.encoded(flat),scaled_delta_upper=old.encoded(delta),
        flat_guard='359723000',delta_guard='11306',total_comparison_guard='359734306',
        comparison_feedback_guard='.00359734306',remaining_budget=old.encoded(remainder),
        best_case_C_rem_ceiling=old.encoded(remainder/arb('1e-11')),
        equal_split_eta=old.encoded(remainder/2),equal_split_C_rem=old.encoded(remainder/arb('2e-11')),
        remaining_budget_gain=old.encoded(ratio),resolvent_ratio=old.encoded(radius),
        cubic_lower=old.encoded(lo),cubic_upper=old.encoded(hi),
        R_X_alpha_squared_comparison_only=old.encoded(arb('5e7')/total),
        R_X_comparison_only_if_alpha_point001=old.encoded(arb('5e7')/total/arb('1e-6')),
        original_norm_comparisons_established=False,actual_controlled_solution=False,
        nonlinear_admission=False,exact_compatibility=False,physical_matching=False,
        validity_noise=False,all_sector_stability=False,
        selected_full_projection_screen_retained=True,projection_newly_evaluated=False,
        audit_status='closed',goal_status='active',physical_domain_changed=False,
        normalization_changed=False,next_complete_kernel_estimate_executed=False)
def main():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    prior=json.loads((ROOT/'evidence/current-strong-response-decision.json').read_bytes())['runs'][-1]
    paths.update(prior['input_sha256_after'])
    before={p:digest(p) for p in sorted(paths)}
    out=ROOT/'evidence'/f'{STEM}.json'
    archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running',
         'input_sha256_before':before,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
         'precision_bits':ctx.prec,'python':sys.version,'invocation':sys.argv}
    archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        check('prior_sources_unchanged',all(digest(p)==h for p,h in prior['input_sha256_after'].items()))
        run['result']=calc();run['status']='passed'
    except Exception:
        run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=checks;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps({'status':run['status'],'checks':len(checks),'inputs':len(before),
                     'exception':run.get('exception'),'result':run.get('result')},indent=2))
    return 0 if run['status']=='passed' else 1
if __name__=='__main__':raise SystemExit(main())
