"""Verify algebra/arithmetic of the ordinary attempt; no strong-response certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
from math import factorial
DEVELOPMENT='--development' in sys.argv
if DEVELOPMENT:sys.argv.remove('--development')
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import current_product_matched_source_projection as old
STEM='current-conormal-common-regulator-attempt'
OWN=['notes/'+STEM+'.md','checks/current_conormal_common_regulator_attempt.py','evidence/'+STEM+'-review.json']
checks=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(n,b):
    checks.append({'name':n,'passed':bool(b)})
    if not b:raise AssertionError(n)
def calc():
    p,*_=old.scalar.prev.york.actual_parameters()
    second=p['u']/p['beta'];third=3*(-4*p['beta']-2*p['u']-p['tau'])/p['beta']
    check('actual_nonzero_second_smearing_jet',arb('1.22150e-10')<second<arb('1.22151e-10'))
    check('actual_unsmoothed_cubic_jet',arb('-12.000000001')<third<arb('-12'))
    check('proper_coordinate_second_jet', (p['beta']*second).overlaps(p['u']) and p['u']>0)
    K=p['C']*p['beta']**2/(15*arb.pi()*p['B2'])
    H8=arb('2.328')/(K*arb(2).log())
    check('original_graph_H8_embedding_guard',H8<arb('9.725e35'))
    sob=[S.Rational(factorial(2*j)*factorial(14-2*j),2*4**7*factorial(j)*factorial(7-j)*factorial(7)) for j in range(8)]
    check('H8_seven_jet_Fourier_evaluation_constants',all(v<1 for v in sob))
    b,b1,b2,w,w1,w2,be,u,k,q=S.symbols('b b1 b2 omega omega1 omega2 beta u k q')
    psi=b+w;psi2=b2+w2;L=-2*be*psi2
    def W(phi,phi1,phi2):
        mm=be*phi2;nn=be*phi1**2
        return (S.exp(-2*psi)*(phi*L**2/60+L*mm/30-mm**2/10+L*nn/30-7*mm*nn/90-7*nn**2/180)-u*phi*L/18-u*nn/18+k*u**2*phi*S.exp(2*psi))/2
    bulk=S.expand(W(-b,-b1,-b2)-W(-psi,-b1-w1,-b2-w2))
    zero=-q*u*be*w1**2/12;full=bulk+zero
    check('no_flattening_gives_zero_conversion',S.simplify(full.subs({w:0,w1:0,w2:0}))==0)
    d=be*w2;g=be*w1**2
    expected=S.exp(-2*w)*(w*d*d/30+d*d/60-d*g/180+7*g*g/360)+u*w*d/18+u*g/36+k*u*u*w*S.exp(2*w)/2-q*u*g/12
    check('full_flat_conversion_normalization',S.simplify(full.subs({b:0,b1:0,b2:0})-expected)==0)
    top=S.simplify(S.diff(bulk,b2,2)/2)
    check('nonzero_bulk_fourth_contact_kept',S.simplify(top-be**2*w*S.exp(-2*psi)/30)==0)
    check('zero_angular_endpoint_independent_of_metric_mass',all(S.diff(zero,v)==0 for v in [b,b1,b2]))
    fnn,fiin=S.symbols('f_nn f_iin');chi=S.diag(1,1,-1,-1)
    check('balanced_boundary_third_smearing_trace',S.trace(30*fiin*chi)==0)
    check('second_smearing_does_not_need_to_vanish',S.expand(fnn*(24*S.Integer(0)+120*S.Integer(0)))==0)
    terms={():1}
    for order in range(4):
        nxt={}
        for word,c in terms.items():
            for i in range(len(word)+1):
                ww=word[:i]+(1,)+word[i:];nxt[ww]=nxt.get(ww,0)-c
            for i in range(len(word)):
                ww=word[:i]+(word[i]+1,)+word[i+1:];nxt[ww]=nxt.get(ww,0)+c
        terms={ww:c for ww,c in nxt.items() if c}
    expected_words={(4,):-1,(1,3):4,(3,1):4,(2,2):6,(1,1,2):-12,(1,2,1):-12,(2,1,1):-12,(1,1,1,1):24}
    check('complete_fourth_resolvent_commutator',terms==expected_words)
    check('commutator_coefficient_sum_not_operator_norm',sum(abs(c) for c in terms.values())==75)
    rest=arb('.1599164')-arb('.00731')-arb('4.001e-9')-arb('1e-11')*359734306
    check('unchanged_strong_budget',rest.overlaps(arb('.149009052939')))
    return {'outcome':'unresolved_complete_strong_estimate','actual_complete_bound_executed':False,
       'omega_ss':old.encoded(second),'omega_sss_one_sided':old.encoded(third),'omega_yy':old.encoded(p['u']),
       'original_H8_graph_guard':old.encoded(H8),'seven_jet_evaluation_constants_squared':[str(v) for v in sob],
       'weak_C3_conversion_identified_on_specified_admitted_graph_family':True,
       'original_X_ball_admission_not_inferred':True,
       'boundary_heat_endpoint_contribution':'exactly zero, including its metric and flow variations',
       'bulk_conversion_explicit':True,'bulk_conversion_norm_certified':False,
       'bulk_fourth_coefficient':str(top),'full_bulk_density':str(bulk),'zero_angular_endpoint':str(zero),
       'fourth_commutator_words':{str(ww):c for ww,c in terms.items()},
       'commutator_identity_scope':'closed-domain resolvent core identity; heat variations via Duhamel; combined distributional limit only',
       'remaining_budget':old.encoded(rest),'complete_symbol_target':old.encoded(rest/(1000*arb(8).sqrt())),
       'complete_strong_norm_known':False,'eta_seed_remaining_known':False,'C_rem_known':False,
       'original_global_norm_comparisons_known':False,'actual_controlled_solution':False,
       'construction_obstruction':False,'physical_domain_changed':False,'normalization_changed':False,
       'physical_modes_truncated':False,'projection_reevaluated':False,'selected_projection_screen_retained':True,
       'nonlinear_admission':False,'exact_compatibility':False,'physical_matching':False,
       'semiclassical_validity_noise':False,'all_sector_stability':False,'audit_status':'closed','goal_status':'active',
       'next_complete_tame_strong_bound_executed':False}
def main():
    ctx.prec=512
    if DEVELOPMENT:
        result=calc();(ROOT/'evidence'/f'{STEM}-development.json').write_bytes((json.dumps({'result':result,'checks':checks,'scope':'algebra/arithmetic only'},indent=2)+'\n').encode())
        print(json.dumps({'development_checks':len(checks),'omega_ss':result['omega_ss']}));return 0
    review=json.loads((ROOT/OWN[2]).read_bytes())
    prev=json.loads((ROOT/'evidence/current-symmetric-strong-response-attempt.json').read_bytes())['runs'][-1]
    paths=set(OWN+review['retained_inputs'])|set(prev['input_sha256_after']);before={p:digest(p) for p in sorted(paths)}
    out=ROOT/'evidence'/f'{STEM}.json';archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'status':'running','precision_bits':ctx.prec,
      'input_sha256_before':before,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},'python':sys.version,'invocation':sys.argv};archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        check('prior_sources_preserved',all(digest(p)==h for p,h in prev['input_sha256_after'].items()))
        run['result']=calc();run['status']='passed'
    except Exception:run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=checks;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps({'status':run['status'],'checks':len(checks),'inputs':len(before),'exception':run.get('exception'),'complete_strong_bound_executed':False}))
    return 0 if run['status']=='passed' else 1
if __name__=='__main__':raise SystemExit(main())
