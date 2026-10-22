"""Supporting checks for the weighted auxiliary cold inverse through eight."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-cold-weighted-inverse-jets';PREVIOUS='seed-outer-eight-jets'
OWN=['checks/seed_cold_weighted_inverse_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','inverse_order':8,'parameter_radius':'1e-4',
          'cap_regularity':10,'ordinary_continuum_proof_separate':True,
          'actual_charged_source_evaluated':False,'additional_normalizer_transport_included':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
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
    check('prior_outer_eight_passed',p['status']=='passed' and len(p['checks'])==76 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('inverse_scope_only',review()['auxiliary_cold_inverse_order']==8
          and review()['auxiliary_cold_inverse_weight_preserved']
          and not review()['complete_cold_source_C8_bounded']
          and not review()['additional_normalizer_transport_bounded'] and not review()['goal_complete'])


def lorentz_chart():
    A,B,P,Q,C,Z=S.symbols('A B P Q C Z')
    eta=S.diag(-1,1,1)
    basis=S.groebner([A*A-B*B-1,P*P-Q*Q-1,C*C-Z*Z-1],A,B,P,Q,C,Z)
    reduce=lambda e:basis.reduce(S.expand(e))[1]
    point=S.Matrix([A*P*C+B*Z,B*P*C+A*Z,Q*C])
    dx=S.Matrix([B*P*C+A*Z,A*P*C+B*Z,0])
    dy=S.Matrix([A*Q,B*Q,P])
    zero('hyperboloid_preserved',reduce((point.T*eta*point)[0]+1))
    zero('chart_g11',reduce((dx.T*eta*dx)[0]-(1+C*C*Q*Q)))
    zero('chart_g12',reduce((dx.T*eta*dy)[0]+Z*Q))
    zero('chart_g22',reduce((dy.T*eta*dy)[0]-1))
    zero('chart_determinant',reduce(1+C*C*Q*Q-Z*Z*Q*Q-P*P))
    zero('unscaled_transverse_speed',reduce((C*dy).dot(eta*(C*dy))-C*C))
    c,t=S.symbols('c t',real=True)
    zero('regular_long_throat_limit',S.limit(S.sinh(c*t)/c,c,0)-t)
    check('unit_parameter_basis',S.Matrix([[1+C*C*Q*Q,-Z*Q],[-Z*Q,1]]).subs(Q,0)==S.eye(2))
    FIXTURES.append({'fixture':'explicit_Lorentz_observation_chart','metric':[['1+u^2','-v*u'],['-v*u','1']],
                     'unit_H2_basis':True,'physical_normal_chart':False,'unscaled_transverse_squared_speed':'C_x^2'})


def branch_and_weight_budgets():
    t=S.symbols('t');n=S.factorial(21)/S.factorial(10)**2
    step=S.expand(sum(S.binomial(21,i)*t**i*(1-t)**(21-i) for i in range(11)))
    check('C10_step_endpoints',step.subs(t,0)==1 and step.subs(t,1)==0
          and all(S.diff(step,t,j).subs(t,v)==0 for j in range(1,11) for v in [0,1]))
    zero('beta_step_derivative',S.diff(step,t)+n*t**10*(1-t)**10)
    check('real_step_first',n/R(4)**10<R(38,10))
    rho=R(1,10000);shift=2*rho/(1-2*rho)
    check('boost_relative_H_shift',shift<R(201,10**6),shift)
    check('quarter_q0_bound',R(19,9)+R(13,24)*R(1,100)<R(212,100))
    check('quarter_q0_derivative',R(26,27)+R(5,1000)+R(4,9)<R(142,100))
    check('quarter_J_bound',R(1015,1000)**2/R(9788,10000)<R(106,100))
    check('quarter_J_derivative',R(106,100)*(R(204,100)+R(146,100))<R(45,10))
    check('quarter_log_q0',R(212,100)/R(9788,10000)<R(2166,1000))
    check('quarter_log_A',R(2166,2000)+R(25,10)/R(985,1000)<R(41,10))
    check('quarter_complex_D_path',R(45,10)/3<R(45,10))
    sd=R(675,10000)
    check('weighted_sigma_log_constant',R(41,10)+R(675,100)/(1-sd)<12
          and R(3045,1000)/(1-sd)<12)
    Hder=(R(135,1000)+2*sd+sd**2)/(1-sd)**2+R(1,10**24)
    check('uniformizing_inverse_derivative',Hder<R(32,100) and 1/(1-Hder)<2,Hder)
    check('inverse_image_margin',R(68,100)/4>R(201,10**6))
    check('inner_shift',2*R(201,10**6)==R(402,10**6))
    check('inner_Cauchy_derivative',12/(R(1,4)-R(402,10**6))<49)
    check('inner_complex_amplitude',R(21,1000)+49*R(402,10**6)/100<R(22,1000))
    check('seam_s_displacement',R(2,100)*R(105,10000)==R(21,100000))
    seam=R(2,100)+R(21,100000)*R(1,10)/(R(75,10000)-R(21,100000))
    check('seam_seed_amplitude',seam<R(23,1000),seam)
    check('complex_chi_displacement',2*shift<R(5,10000))
    check('complex_step_branch',n*(R(1,4)+2*R(5,100)**2)**10<R(46,10))
    check('complex_cap_amplitude',R(23,1000)*(1+R(5,100)*R(46,10))<R(3,100))
    check('real_cap_first_derivative',R(2,100)*R(1,10)/(R(75,10000)-R(21,100000))+100*R(38,10)*R(23,1000)<10)
    check('cold_log_weight_margin',3*R(15,1000)<1)
    check('taper_weight_margin',R(3,100)<100*R(1,100))
    check('finite_jet_regularities',10>=8+1)
    check('full_parameter_tensor_factor',15000**2>2/rho**2)
    check('unweighted_coefficient_margin',1/(1-R(3,100))-1<R(1,10))
    check('weighted_coefficient_margin',100/(1-R(3,100))<200)
    check('base_endpoint_margin',1/(1-R(3,200))<R(51,50))
    check('weighted_inverse_margin',1/(1-R(3,100)-R(2,18))<R(6,5))
    FIXTURES.append({'fixture':'true_C10_cap_and_weighted_coefficient_jets','parameter_radius':'1e-4',
                     'k':15000,'u':'1/10','w':200,'inverse_similarity':'6/5',
                     'global_spatial_cutoff_analyticity_assumed':False,'h_terms_discarded':False})


def compositions(n):
    if n==0:
        yield ()
    else:
        for first in range(1,n+1):
            for rest in compositions(n-first):yield (first,)+rest


def ordered_ledger():
    d=R(6,5);u=R(1,10);w=S.Integer(200);e0=R(51,50);v=1+d*u
    endpoint=[e0]+[u]*8;inverse=[d]+[d*d*u*v**(j-1) for j in range(1,9)]
    U=[];out=[];counts=[]
    t=S.symbols('t')
    generating=d*(e0+(u-e0)*t)**2/((1-t)*(1-v*t))
    series=S.Poly(S.series(generating,t,0,9).removeO(),t)
    for j in range(1,9):
        comps=list(compositions(j));counts.append(len(comps))
        check('ordered_compositions_'+str(j),len(comps)==2**(j-1)
              and all(sum(len(c)==ell for c in comps)==S.binomial(j-1,ell-1) for ell in range(1,j+1)))
        zero('inverse_majorant_'+str(j),sum(d**(len(c)+1)*u**len(c) for c in comps)-inverse[j])
        value=sum(endpoint[a]*inverse[b]*endpoint[j-a-b] for a in range(j+1) for b in range(j-a+1))
        zero('endpoint_generating_function_'+str(j),value-series.nth(j))
        bound=w/u*value;normalized=bound*150**j
        check('weighted_inverse_bound_'+str(j),bound<2000 and normalized<10**21)
        U.append(value);out.append({'order':j,'B_exact':str(bound),'B_decimal':str(S.N(bound,18)),
                                   'normalized_exact':str(normalized),'normalized_decimal':str(S.N(normalized,18))})
    check('eighth_bound',w/u*U[-1]<1987 and w/u*U[-1]*150**8<R(5093,1000)*10**20)
    zero('zeroth_difference_budget',w*(d*(e0+1)+d*d)-R(7728,10))
    check('zeroth_difference_reserve',w*(d*(e0+1)+d*d)<1000)
    FIXTURES.append({'fixture':'all_ordered_insertions_and_endpoint_ledger','composition_counts':counts,
                     'weighted_entries':out,'one_weight_reserved_in_every_positive_order_term':True,
                     'both_endpoint_factors_retained':True})


def noncommuting_matrix_inverse():
    t=S.symbols('t');I=S.eye(2);D=S.Matrix([[0,3],[-3,0]]);mass=S.Integer(18)
    E=S.diag(1+t/20+t*t/100,1-t/25+t*t/200)
    T=D+mass*E*E;T0=T.subs(t,0);R0=T0.inv()
    coeff=lambda matrix,j:matrix.applyfunc(lambda x:S.expand(x).coeff(t,j))
    Tjets=[coeff(T,j) for j in range(9)]
    Ejets=[coeff(E,j) for j in range(9)]
    check('matrix_vertices_noncommute',Tjets[1]*R0!=R0*Tjets[1])
    inverse=[R0]
    for j in range(1,9):
        value=S.zeros(2)
        for comp in compositions(j):
            term=R0
            for order in comp:term=term*Tjets[order]*R0
            value+=(-1)**len(comp)*term
        inverse.append(value)
    direct=T.inv()
    direct_coeff=[S.zeros(2) for _ in range(9)]
    for a in range(2):
        for b in range(2):
            poly=S.Poly(S.series(direct[a,b],t,0,9).removeO(),t)
            for j in range(9):direct_coeff[j][a,b]=poly.nth(j)
    complete=[]
    for j in range(9):
        check('independent_noncommuting_inverse_'+str(j),inverse[j]==direct_coeff[j])
        total=sum((Ejets[a]*inverse[b]*Ejets[j-a-b] for a in range(j+1) for b in range(j-a+1)),S.zeros(2))
        complete.append(total)
    physical=E*direct*E
    first_direct=physical.diff(t).subs(t,0).applyfunc(S.simplify)
    check('canonical_endpoint_derivative',first_direct==complete[1])
    omitted=complete[1]-inverse[1]
    check('omitted_endpoint_negative_control',omitted!=S.zeros(2))
    coeffs=[complete[j][0,0] for j in range(9)]
    for j in [4,8]:
        direct_j=S.Poly(S.series(physical[0,0],t,0,j+1).removeO(),t).nth(j)
        zero('independent_canonical_matrix_'+str(j),direct_j-coeffs[j])
    FIXTURES.append({'fixture':'independent_noncommuting_canonical_inverse','order':8,
                     'omitted_endpoint_error':str(omitted),'coefficients_00':[str(v) for v in coeffs],
                     'actual_charged_source_evaluated':False})
    return coeffs


def independent_contour_integrals(coeffs):
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            radius=mp.mpf('.025')
            def value(z):
                e1=1+z/20+z*z/100;e2=1-z/25+z*z/200
                return 18*e1*e1*e2*e2/(324*e1*e1*e2*e2+9)
            for j in [1,4,8]:
                integrand=lambda angle:value(radius*mp.exp(1j*angle))*mp.exp(-1j*j*angle)
                integral=mp.quad(integrand,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi*radius**j)
                exact=mp.mpf(str(S.numer(coeffs[j])))/mp.mpf(str(S.denom(coeffs[j])))
                error=abs(integral-exact)
                check('independent_contour_'+str(dps)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,10))
                records.append({'dps':dps,'order':j,'absolute_error':mp.nstr(error,16)})
    FIXTURES.append({'fixture':'independent_complex_contour_coefficients','records':records,
                     'actual_charged_source_evaluated':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold weighted inverse jets','input_count':len(ps)}),flush=True)
    try:
        provenance();lorentz_chart();branch_and_weight_budgets();ordered_ledger()
        coefficients=noncommuting_matrix_inverse();independent_contour_integrals(coefficients)
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold weighted inverse jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_auxiliary_cold_inverse_C8_bounded':True,'cold_weight_preserved':True,
         'eighth_normalized_inverse_reserve':'5.093e20','zeroth_difference_reserve':1000,
         'additional_normalizer_transport_bounded':False,'complete_cold_source_C8_bounded':False,
         'actual_metric_tenth_derivative_bounded':False,'actual_finite_metric_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'original_cold_C3_contract_changed':False,
         'endpoint_C5_contract_changed':False,'actual_signed_compatibility_zero_proved':False,
         'nonlinear_solution_constructed':False,'physical_stability_proved':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'ordinary_proof_machine_verified':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
