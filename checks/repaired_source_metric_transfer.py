"""Supporting identities and budgets for the restored endpoint contract."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='repaired-source-metric-transfer';PREVIOUS='seed-middle-cap-transfer'
OWN=['checks/repaired_source_metric_transfer.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':3,'endpoint_metric_order':5,
          'ordinary_continuum_review_separate':True,'actual_signed_source_evaluated':False}


def check(name,ok,detail=''):
    assert name not in {row['name'] for row in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    value=S.simplify(expr);check(name,value==0,value)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_cap_repair_passed',p['status']=='passed' and len(p['checks'])==91 and all(row['passed'] for row in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('review_inputs_in_snapshot',set(review()['retained_inputs'])<=set(paths()))
    check('restored_input_is_old_allowance',p['complete_middle_C3_reserve']=='1e26' and not p['withdrawn_1e22_claim_restored'])
    check('review_scope',review()['endpoint_C5_source_dependency_resolved'] and not review()['all_descendant_estimates_reverified']
          and not review()['actual_finite_metric_source_H6_bounded'] and not review()['goal_complete'])


def regional_budgets():
    cold=R(10)**26*4*R(352,100000)/(S.factorial(3)*100**3)
    check('cold_gluing_normalized',cold<R(24,10)*10**17,cold)
    check('cold_complete_three_allowance',cold+10**8<10**20)
    check('cold_global_weight_guard',R(1,1000)*(1+R(231,100))==R(331,100000))
    check('outer_complete_three_allowance',2*2*10**19+4*10**11+10**8+R(1,10**100)<10**20)
    outer_conversion=[1,1+R(1,200),1+R(1,100)+R(1,30000)]
    check('outer_Euler_to_ordinary_factor',max(outer_conversion)<2)
    # Exact identities, independent of the geometric/kernel estimates reviewed in the note.
    y=S.symbols('y',positive=True);f=S.Function('f')(y);D=lambda v:y*S.diff(v,y)
    zero('Euler_second_derivative',D(D(f))-D(f)-y*y*S.diff(f,y,2))
    zero('Euler_third_derivative',D(D(D(f)))-3*D(D(f))+2*D(f)-y**3*S.diff(f,y,3))
    FIXTURES.append({'fixture':'distinct_regional_source_budgets','cold_gluing':str(cold),
                     'cold_and_outer_cap_derivatives_reviewed_in_ordinary_proof':True,
                     'continuum_kernel_arguments_machine_verified':False})


def ward_and_middle():
    x=S.symbols('x');b=S.Function('b')(x);I=S.Function('I')(x);theta=S.Function('theta')(x)
    P=I/(2*b);a=S.diff(b,x)/b
    zero('canonical_Ward_derivative',(S.diff(P,x)+a*P-a*theta/2).subs(S.diff(I,x),S.diff(b,x)*theta))
    ps=[R(3,2)]
    for j in range(3):
        ps.append(S.simplify((sum(ps)+R(j+1,2))/(j+1)))
    check('normalized_pressure_recurrence',ps==[R(3,2),S.Integer(2),R(9,4),R(29,12)],ps)
    check('original_pressure_reserves_hold',all(p<=bound for p,bound in zip(ps[1:],[2,4,8])))
    check('density_reserve_holds',max(p+1 for p in ps)<16)
    check('full_Ward_interval_low_bound',10**11+49*10**15<10**26)
    cap=next(f for f in previous()['fixtures'] if f['fixture']=='actual_cap_and_complete_source_budget')
    B=[S.sympify(v) for v in cap['complete_bounds']]
    check('repaired_complete_source_fits_original_input',max(B)<10**26)
    check('withdrawn_improvement_not_used',B[3]>10**22)
    f3=2*R(1,2)*S.factorial(2)*100**2*16/(S.factorial(3)*100**3)
    check('weighted_third_metric_source_coefficient',f3<R(6,100),f3)
    gamma=R(68,10**29);Km=S.Integer(10)**26;force=10*gamma*Km
    check('original_curvature_force_reserve',force==R(68,100),force)
    u4=force*S.factorial(2)*100**2/(2*S.factorial(4)*100**4)
    u5=force*S.factorial(3)*100**3/(2*S.factorial(5)*100**5)
    check('fourth_metric_high_jet',u4<R(284,10**8),u4)
    check('fifth_metric_high_jet',u5<=R(17,10**7),u5)
    check('middle_feedback_and_conversion_margin',2*max(u4,u5)/(1-R(1,10))<R(1,10000))
    FIXTURES.append({'fixture':'repaired_source_Ward_and_middle_transfer','pressure_bounds':[str(p) for p in ps],
                     'curvature_force_reserve':str(force),'high_metric_jets':[str(u4),str(u5)],
                     'original_endpoint_constants_changed':False})


def regular_radial_and_endpoint():
    r,t=S.symbols('r t',positive=True)
    for power in [0,2,4,6,8]:
        force=r**power
        u=r**(power+2)/(power+2)**2
        zero('regular_radial_equation_'+str(power),S.diff(u,r,2)+S.diff(u,r)/r-force)
        zero('regular_fifth_derivative_'+str(power),S.diff(u,r,5)-S.diff(force,r,3)
             +S.integrate(t**4*S.diff(force,r,3).subs(r,t*r),(t,0,1)))
    ell=S.symbols('ell');u=S.Function('u')(r)
    sixth=S.diff(S.exp(ell*u),r,5)
    orders=[sum(count for _,count in d.variable_count) for d in sixth.atoms(S.Derivative)]
    check('endpoint_Phi_sixth_uses_u_through_five',max(orders)==5)
    A=S.Function('A')(r)
    relative=S.diff(A/r,r,6)
    check('relative_endpoint_sixth_uses_Phi_through_six',max(sum(count for _,count in d.variable_count)
          for d in relative.atoms(S.Derivative))==6)
    # The permitted logarithm remains conormal rather than being assigned a closed-collar C6 norm.
    logterm=r**4*S.log(r)
    zero('log_endpoint_integral',S.diff(r**5*(S.log(r)/5-R(1,25)),r)-logterm)
    check('allowed_log_is_not_classical_C5',S.limit(S.diff(logterm,r,5),r,0,dir='+')==S.oo)
    FIXTURES.append({'fixture':'independent_regular_radial_and_endpoint_order',
                     'even_polynomial_forcing_orders':[0,2,4,6,8],
                     'endpoint_requires_sixth_metric_derivative':False,
                     'outer_logarithms_retained':True})


def higher_source_order():
    x=S.symbols('x');v=S.Function('v')(x);eps=S.symbols('eps',positive=True)
    scalar=S.exp(-2*v)*(-6*S.diff(v,x,2)-6*S.diff(v,x)**2)
    lap=S.exp(-2*v)*(S.diff(scalar,x,2)+2*S.diff(v,x)*S.diff(scalar,x))
    zero('fourth_trace_principal_coefficient',S.expand(eps*lap).coeff(S.diff(v,x,4))+6*eps*S.exp(-4*v))
    tenth=S.diff(eps*lap,x,6).expand().coeff(S.diff(v,x,10))
    zero('sixth_source_tenth_metric_coefficient',tenth+6*eps*S.exp(-4*v))
    check('principal_coefficient_nonzero',(-6*eps*S.exp(-4*v))!=0)
    records=[]
    for n in [8,16,32]:
        amplitude=R(1,n**5)
        check('oscillatory_C5_bound_'+str(n),max(amplitude*n**j for j in range(6))==1)
        records.append({'frequency':n,'normalized_source_H6_amplitude':str(6*amplitude*n**10)})
    check('source_H6_grows_by_fifth_power',S.sympify(records[1]['normalized_source_H6_amplitude'])/
          S.sympify(records[0]['normalized_source_H6_amplitude'])==32)
    check('direct_reconstruction_source_order',10-2==8)
    FIXTURES.append({'fixture':'nonzero_fourth_order_trace_and_norm_gap','records':records,
                     'actual_comparison_family_evaluated':False,'physical_instability_established':False,
                     'direct_metric_order_required':10,'direct_seed_source_order_required':8,
                     'these_orders_are_sufficient_for_full_Schur_bound':False})


def independent_radial_integrals():
    records=[]
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            F3=lambda x:(12*x-8*x**3)*mp.exp(-x*x)+8*mp.sin(2*x)
            for rr in ['0.3','0.7','1.1']:
                r=mp.mpf(rr)
                # Solve the regular radial equation coefficient by coefficient for exp(-r^2)+cos(2r).
                direct=mp.fsum(((-1)**n/mp.factorial(n)+(-1)**n*2**(2*n)/mp.factorial(2*n))
                     *mp.factorial(2*n+2)/mp.factorial(2*n-3)*r**(2*n-3)/(2*n+2)**2 for n in range(2,100))
                integral=F3(r)-mp.quad(lambda t:t**4*F3(t*r),[0,1])
                error=abs(direct-integral)
                check('independent_radial_integral_'+str(dps)+'_'+rr,error<mp.mpf('1e-42'),mp.nstr(error,8))
                records.append({'dps':dps,'radius':rr,'series_value':mp.nstr(direct,45),
                                'quadrature_value':mp.nstr(integral,45),'absolute_difference':mp.nstr(error,12)})
    FIXTURES.append({'fixture':'independent_nonpolynomial_radial_integrals','records':records,
                     'forcing':'exp(-r^2)+cos(2r)','actual_charged_source_evaluated':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 repaired source metric transfer','input_count':len(ps)}),flush=True)
    try:
        provenance();regional_budgets();ward_and_middle();regular_radial_and_endpoint();higher_source_order();independent_radial_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 repaired source metric transfer','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_endpoint_C5_source_dependency_resolved':True,'endpoint_C5_contract_changed':False,
         'withdrawn_1e22_claim_restored':False,'all_descendant_estimates_reverified':False,
         'actual_finite_metric_source_H6_bounded':False,'actual_transfer_certificate_passed':False,
         'actual_signed_compatibility_zero_proved':False,'nonlinear_solution_constructed':False,
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
