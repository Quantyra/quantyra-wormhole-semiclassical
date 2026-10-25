"""Arithmetic supporting an unresolved complete-response attempt, not a PDE proof."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
DEVELOPMENT = '--development' in sys.argv
if DEVELOPMENT:
    sys.argv.remove('--development')  # retained dependency parses argv on import
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb, ctx
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'checks'))
import current_product_matched_source_projection as old
STEM = 'current-symmetric-strong-response-attempt'
OWN = ['notes/'+STEM+'.md', 'checks/current_symmetric_strong_response_attempt.py',
       'evidence/'+STEM+'-review.json']
checks = []

def digest(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()

def check(name, value):
    checks.append({'name': name, 'passed': bool(value)})
    if not value:
        raise AssertionError(name)

def calc():
    h = arb('.005006'); delta = arb('.0001'); floor = arb('.999')
    f = 140*h**3+arb('1e-8')*h**2
    check('actual_all_annuli_cap_radius', 1-f>floor and 2*f/(1-f)<delta)
    def sn(n):
        return arb((n+2)**5-(n+2))/30
    check('arm_sum_closed_formula', all(
        sn(n).overlaps(arb(sum(a*a*(n+2-a)**2 for a in range(1,n+2))))
        for n in range(1,30)))
    nmax = 17
    partial = sum((arb(8)**n*delta**(n-2)*sn(n) for n in range(2,nmax+1)), arb(0))
    ratio = 8*delta*(arb(nmax+4)/(nmax+3))**6
    tail = 8**(nmax+1)*delta**(nmax-1)*sn(nmax+1)/(1-ratio)
    check('entire_background_order_tail', ratio<arb('.002') and tail<arb('1e-40'))
    c2 = 2/arb.pi()**2
    d0 = 8192*c2*(1+delta)**2/floor**2
    ep = 12*c2/floor**2
    aa = 32*d0+ep*(1+delta)
    bb = d0*(partial+tail)+ep
    check('full_angular_midpoint_envelopes', aa<54000 and bb<arb('3.7e6'))
    check('whole_time_integral_guards', 60000*arb.pi()/2<arb('1e5') and arb('4e6')*arb.pi()/2<arb('7e6'))
    schur = 8*arb(106).sqrt()/arb.pi()
    theorem = 27+27*arb(12).sqrt()+324/(arb(2).sqrt()-1)
    check('explicit_rough_symbol_operator_constant', schur<27 and theorem<1000)
    check('positive_lowpass_half_moment', 192/(80*arb.pi())<1)
    rest = arb('.1599164')-arb('.00731')-arb('4.001e-9')-arb('1e-11')*359734306
    check('unchanged_complete_response_budget', rest.overlaps(arb('.149009052939')))
    theorem_guard = 1000*arb(8).sqrt()
    target = rest/theorem_guard
    check('same_H3_units', (theorem_guard*target).overlaps(rest))
    return {
        'outcome': 'unresolved_complete_strong_estimate',
        'actual_complete_bound_executed': False,
        'midpoint_linear_gradient_cancellation_proved_in_ordinary_argument': True,
        'cap_zero_winding_value_bound_only': True,
        'all_background_orders_in_this_bound_included': True,
        'physical_modes_truncated': False,
        'cap_radius': old.encoded(2*f/(1-f)),
        'order_sum_upper': old.encoded(partial+tail),
        'order_tail_upper': old.encoded(tail),
        'A2_coefficient_upper': old.encoded(aa),
        'A1_squared_coefficient_upper': old.encoded(bb),
        'rough_symbol_theorem_constant_upper': old.encoded(theorem),
        'rough_symbol_guard': '1000',
        'remaining_budget': old.encoded(rest),
        'complete_symbol_target': old.encoded(target),
        'equal_split_M_seed': old.encoded(target/2),
        'equal_split_M_var': old.encoded(target/arb('2e-11')),
        'best_case_M_var': old.encoded(target/arb('1e-11')),
        'M_seed_known': False, 'M_var_known': False,
        'eta_seed_remaining_known': False, 'C_rem_known': False,
        'strong_regulator_bound_established': False,
        'local_contact_targets_established': False,
        'original_X_comparisons_established': False,
        'construction_obstruction': False,
        'actual_controlled_solution': False,
        'nonlinear_admission': False, 'exact_compatibility': False,
        'physical_matching': False, 'validity_noise': False, 'all_sector_stability': False,
        'physical_domain_changed': False, 'normalization_changed': False,
        'selected_projection_screen_retained': True, 'projection_reevaluated': False,
        'audit_status': 'closed', 'goal_status': 'active',
        'next_full_conormal_commutator_estimate_executed': False}

def main():
    ctx.prec = 512
    if DEVELOPMENT:
        result = calc()
        (ROOT/'evidence'/f'{STEM}-development.json').write_bytes((json.dumps(
            {'scope': 'arithmetic development only; no complete strong estimate',
             'result': result, 'checks': checks}, indent=2)+'\n').encode())
        print(json.dumps({'development_checks':len(checks), 'target':result['complete_symbol_target']}))
        return 0
    review = json.loads((ROOT/OWN[2]).read_bytes())
    prev = json.loads((ROOT/'evidence/current-annular-response-feasibility.json').read_bytes())['runs'][-1]
    paths = set(OWN+review['retained_inputs']) | set(prev['input_sha256_after'])
    before = {p:digest(p) for p in sorted(paths)}
    out = ROOT/'evidence'/f'{STEM}.json'
    archive = json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run = {'started_utc':datetime.now(timezone.utc).isoformat(), 'status':'running',
           'input_sha256_before':before, 'precision_bits':ctx.prec,
           'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
           'python':sys.version, 'invocation':sys.argv}
    archive['runs'].append(run)
    def save():
        out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen', digest(OWN[0])==review['ordinary_argument_sha256'])
        check('previous_scientific_sources_preserved', all(digest(p)==h for p,h in prev['input_sha256_after'].items()))
        run['result']=calc(); run['status']='passed'
    except Exception:
        run['status']='failed'; run['exception']=traceback.format_exc()
    finally:
        run['checks']=checks; run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:
            run['status']='failed'; run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat(); save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:
        f.write((json.dumps(run)+'\n').encode())
    result=run.get('result',{})
    print(json.dumps({'status':run['status'],'checks':len(checks),'inputs':len(before),
          'exception':run.get('exception'),'target':result.get('complete_symbol_target'),
          'complete_strong_bound_executed':False},indent=2))
    return 0 if run['status']=='passed' else 1

if __name__=='__main__':
    raise SystemExit(main())
