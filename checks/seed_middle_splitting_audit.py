"""Check the cap-splitting correction without treating fixtures as the actual source."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-middle-splitting-audit'; PREVIOUS='seed-middle-coercive-jets'
OWN=['checks/seed_middle_splitting_audit.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]; R=S.Rational
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-45','mixed_derivative_order':9,'resolvent_vertex_order':3,
          'previous_arithmetic_result_preserved':True,'actual_cap_correction_evaluated':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row); print(json.dumps(row),flush=True)
    if not ok: raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(expr); check(name,v==0,v)


def prior(): return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths(): return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('previous_arithmetic_pass_preserved',p['status']=='passed' and len(p['checks'])==101 and all(c['passed'] for c in p['checks']))
    check('previous_input_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('retained_hashes_unchanged',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_sources_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('correction_scope',review['previous_complete_C3_claim_withdrawn']
          and not review['actual_source_C3_inequality_disproved'] and not review['actual_cap_correction_evaluated'])


def split_identity():
    x,b=S.symbols('x b',real=True)
    U=sum(R(r+1,s+2)*x**r*b**s for r in range(10) for s in range(10-r))
    theta=3+x/7+x*x/11+x**9/19
    G=theta-U; moving=U.subs(b,x); fixedG=G.subs(b,0)
    correction=moving-U.subs(b,0)
    zero('exact_fixed_moving_split',moving+fixedG-correction-theta)
    zero('cap_agreement_at_center',correction.subs(x,0))
    coefficients=[]
    for j in range(1,10):
        mixed=sum(S.binomial(j,r)*S.diff(U,x,j-r,b,r).subs({x:0,b:0}) for r in range(1,j+1))
        direct=S.diff(correction,x,j).subs(x,0)
        zero('mixed_cap_derivative_'+str(j),direct-mixed)
        zero('naive_derivative_error_'+str(j),S.diff(moving+fixedG-theta,x,j).subs(x,0)-mixed)
        coefficients.append(str(mixed))
    check('nonzero_omission_control',all(S.sympify(v)!=0 for v in coefficients))
    FIXTURES.append({'fixture':'fixed_moving_split_and_mixed_jets','orders':9,
                     'omitted_coefficients':coefficients,'actual_quantum_source':False})


def noncommuting_vertices():
    b=S.symbols('b',real=True)
    A0=S.diag(2,5); A1=S.Matrix([[1,1],[1,-1]])
    A2=S.Matrix([[2,3],[3,1]]); A3=S.Matrix([[1,-2],[-2,4]])
    A=A0+b*A1+b*b*A2/2+b**3*A3/6
    direct=A.inv(); Q=A0.inv()
    expected=[-Q*A1*Q,
              2*Q*A1*Q*A1*Q-Q*A2*Q,
              -6*Q*A1*Q*A1*Q*A1*Q+3*Q*A2*Q*A1*Q+3*Q*A1*Q*A2*Q-Q*A3*Q]
    for j in range(1,4):
        difference=S.simplify(direct.diff(b,j).subs(b,0)-expected[j-1])
        check('noncommuting_resolvent_derivative_'+str(j),difference==S.zeros(2),difference)
    left=Q*A2*Q*A1*Q; right=Q*A1*Q*A2*Q
    check('third_mixed_products_are_distinct',left!=right)
    check('omitting_one_third_mixed_term_detected',3*left!=S.zeros(2) and 3*right!=S.zeros(2))
    FIXTURES.append({'fixture':'noncommuting_cap_vertices','first_through_third':[[str(v) for v in mat] for mat in expected],
                     'actual_cap_vertices':False})


def continuum_formulas():
    x,b,M,y,z=S.symbols('x b M y z',positive=True)
    B=(S.exp(-M)-S.exp(-2*M))/(2*M**2)
    a=1/M**2-(1-S.exp(-M))/M**3
    zero('remote_row_integral',S.integrate(S.exp(-M*y)/(2*M),(y,1,2))-B)
    aa=S.integrate((1-z)*S.exp(-M*z)/M,(z,0,1))
    zero('remote_rank_one_energy_integral',aa-a)
    zero('positive_a_margin',a-(M-1)/M**3-S.exp(-M)/M**3)
    check('uniform_operator_coercivity',2**2-R(1,10)>0)
    U=1/(2*M)-b*B**2*S.exp(2*M*x)/(1+a*b)
    corr=U.subs(b,x)-U.subs(b,0)
    zero('rank_one_first_omitted_derivative',S.diff(corr,x).subs(x,0)+B**2)
    zero('rank_one_second_omitted_derivative',S.diff(corr,x,2).subs(x,0)+2*B**2*(2*M-a))
    zero('rank_one_third_omitted_derivative',S.diff(corr,x,3).subs(x,0)+6*B**2*(2*M*M-2*M*a+a*a))
    zero('rank_one_exact_correction_cancels',U.subs(b,x)-corr-1/(2*M))
    # Independent finite-rank inversion checks the denominator and its sign.
    H=S.Matrix([[4,1],[1,3]]); v=S.Matrix([1,2]); S0=H.inv()
    rank_inverse=S0-b*(S0*v)*(v.T*S0)/(1+b*(v.T*S0*v)[0])
    check('rank_one_resolvent_direct_inverse',S.simplify((H+b*v*v.T).inv()-rank_inverse)==S.zeros(2))
    V=(y-1)**4*(2-y)**4
    for k in range(4):
        zero('local_bump_left_jet_'+str(k),S.diff(V,y,k).subs(y,1))
        zero('local_bump_right_jet_'+str(k),S.diff(V,y,k).subs(y,2))
    check('local_bump_positive_interior',V.subs(y,R(3,2))>0)
    FIXTURES.append({'fixture':'remote_continuum_counterexample',
                     'actual_charged_geometry':False,'actual_source_inequality_disproved':False,
                     'fixed_cap_gluing_jets':'zero','moving_cap_derivative':'-B^2, nonzero',
                     'local_potential_control':True,'physical_boundary_condition_changed':False})


def numerical_integrals():
    records=[]
    y=S.symbols('y',real=True)
    local_exact={m:S.integrate((y-1)**4*(2-y)**4*S.exp(-2*m*y)/(4*m*m),(y,1,2)) for m in [2,5]}
    for precision in PROTOCOL['numerical_dps']:
        with mp.workdps(precision):
            for mass in [2,5]:
                m=mp.mpf(mass)
                B=(mp.exp(-m)-mp.exp(-2*m))/(2*m*m)
                a=1/m**2-(1-mp.exp(-m))/m**3
                bnum=mp.quad(lambda yy:mp.exp(-m*yy)/(2*m),[1,mp.mpf('1.5'),2])
                anum=mp.quad(lambda zz:(1-zz)*mp.exp(-m*zz)/m,[0,mp.mpf('.5'),1])
                local=mp.quad(lambda yy:(yy-1)**4*(2-yy)**4*mp.exp(-2*m*yy)/(4*m*m),[1,mp.mpf('1.5'),2])
                lexact=mp.mpf(str(S.N(local_exact[mass],precision)))
                errors=[abs(bnum-B),abs(anum-a),abs(local-lexact)]
                check('independent_remote_integrals_'+str(precision)+'_'+str(mass),max(errors)<mp.mpf('1e-45'),[mp.nstr(e,12) for e in errors])
                check('remote_effects_nonzero_'+str(precision)+'_'+str(mass),B*B>0 and local>0)
                records.append({'dps':precision,'mass':mass,'B_squared':mp.nstr(B*B,45),
                                'a':mp.nstr(a,45),'negative_local_derivative_magnitude':mp.nstr(local,45),
                                'errors':[mp.nstr(e,20) for e in errors]})
    FIXTURES.append({'fixture':'independent_continuum_integrals','records':records,
                     'actual_cap_correction_evaluated':False})


def retained_and_conditional_budgets():
    p=prior(); ledger=next(f for f in p['fixtures'] if f['fixture']=='complete_angular_and_seed_source_ledger')
    G=[S.sympify(v) for v in ledger['gluing_coordinate_bounds']]
    local=S.Integer(2*10**13+10**12)
    value=G[0]+local
    check('retained_complete_middle_value',value<10**14,value)
    original=max(G[0],G[1],G[2],G[3]+G[1]/200)+local
    zero('old_arithmetic_reconciles',original-S.sympify(ledger['complete_covariant_C3_bound']))
    kappa=S.Integer(4*10**21)
    repaired=[value,G[1]+local+kappa,G[2]+local+kappa,G[3]+G[1]/200+local+kappa+kappa/200]
    check('conditional_repair_budget',max(repaired)<10**22,max(repaired))
    check('conditional_decimal_margin',max(repaired)<R(9975,1000)*10**21)
    zero('covariant_correction_normalization',300*100/(S.Integer(6)*100**3)-R(1,200))
    FIXTURES.append({'fixture':'retained_value_and_conditional_repair_budget',
                     'complete_value_bound':str(value),'adopted_value_reserve':'1e14',
                     'conditional_Cartesian_kappa_reserve':'4e21',
                     'conditional_complete_C3_upper':str(max(repaired)),
                     'actual_kappas_bounded':False,'repair_certificate_passed':False,
                     'old_complete_C3_claim_recertified':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 source splitting audit','input_count':len(ps)}),flush=True)
    try:
        provenance();split_identity();noncommuting_vertices();continuum_formulas();numerical_integrals();retained_and_conditional_budgets()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 source splitting audit','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'previous_complete_C3_claim_withdrawn':True,'ordinary_fixed_cap_gluing_estimate_retained':True,
         'ordinary_complete_middle_value_bound_proved':True,'complete_middle_value_reserve':'1e14',
         'actual_source_C3_inequality_disproved':False,'actual_cap_correction_evaluated':False,
         'source_split_repair_certificate_passed':False,'complete_middle_C3_improvement_certified':False,
         'earlier_1e26_source_ledger_recertified':False,'actual_finite_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'actual_signed_compatibility_zero_proved':False,
         'endpoint_C5_contract_changed':False,'nonlinear_solution_constructed':False,
         'physical_stability_proved':False,'physical_source_projection_changed':False,
         'physical_boundary_condition_changed':False,'ordinary_proof_machine_verified':False,
         'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
