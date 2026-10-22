"""Supporting checks for the complete outer seed source through eight."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
from sympy.functions.combinatorial.numbers import stirling
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-outer-eight-jets';PREVIOUS='seed-middle-eight-jets'
OWN=['checks/seed_outer_eight_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'parameter_radius':'1e-5',
          'cap_regularity':11,'actual_cap_distance':'24','positive_defect_distance':'3.9',
          'ordinary_continuum_proof_separate':True,'actual_signed_source_evaluated':False}


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
    check('prior_middle_eight_passed',p['status']=='passed' and len(p['checks'])==103 and all(row['passed'] for row in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('outer_scope',review()['complete_outer_source_order']==8 and review()['complete_outer_source_reserve']=='1e44'
          and not review()['cold_C8_bounded'] and not review()['actual_finite_metric_source_H6_bounded'] and not review()['goal_complete'])


def finite_cap_family():
    rho=R(1,100000);c=S.factorial(21)/S.factorial(10)**2
    check('real_step_first',c/R(4)**10<R(38,10))
    third=c*(20/R(4)**9+90/R(4)**8)
    check('width_half_third',8*third<45000,8*third)
    check('negative_blend_third',1000*third<6*10**6)
    capjets=[]
    for name,w,w1,w2 in [('positive',R(101,100),R(2021,1000),S.Integer(13)),
                         ('negative',S.Integer(1),S.Integer(3),S.Integer(9))]:
        cw1=R(76,10)*w+w1;cw2=160*w+2*R(76,10)*w1+w2
        f1=S.Integer(2);f2=cw1+f1**2;f3=cw2+3*f1*cw1+f1**3
        check(name+'_real_low_cap_jets',f2<20 and f3<400,[cw1,cw2,f2,f3])
        capjets.append({'cap':name,'second':str(f2),'third':str(f3)})
    check('complex_step_branch',c*(R(1,4)+2*R(2,100)**2)**10<4)
    check('complex_taper_value',1+4*R(2,100)==R(108,100))
    check('complex_branch_f_ratio',1/(1-R(432,10000))<2)
    check('cap_and_jet_regularities',10+1>=8+3)
    check('positive_cap_support',3*(24+7)+R(1,2)<94<100)
    check('dilation_support_displacement',100*rho/(1-rho)<R(11,10000)<R(1,100))
    base=[S.Integer(1),S.Integer(2),S.Integer(20),S.Integer(400)]
    perturb=[]
    for a in range(4):
        total=S.Integer(0)
        for j in range(1,9):
            coeff=S.Integer(0)
            for k in range(1,j+1):
                for ell in range(min(a,k)+1):
                    m=k+a-ell
                    coeff+=stirling(j,k,kind=2)*S.binomial(a,ell)*S.factorial(k)/S.factorial(k-ell)*2*S.factorial(m)*100**m
            total+=rho**j*coeff/S.factorial(j)
        perturb.append(total)
    check('negative_finite_value',perturb[0]<R(3,1000),perturb[0])
    for a,target in [(1,R(5,2)),(2,S.Integer(150)),(3,S.Integer(50000))]:
        check('negative_finite_spatial_'+str(a),base[a]+perturb[a]<target,S.N(base[a]+perturb[a],18))
    tail=10**6*R(1,1000)**9/(1-R(1,1000))
    check('negative_blend_tail',tail*10**8<R(1,10**6))
    check('negative_blend_product_reserve',1+3*38+3*4000+6*10**6<10**8)
    check('negative_global_low_guards',(R(5,2)+R(1,10**6))/R(997,1000)<5
          and (150+R(1,10**6))/R(997,1000)<1000 and (50000+R(1,10**6))/R(997,1000)<100000)
    FIXTURES.append({'fixture':'C11_cap_and_eight_parameter_jets','real_cap_budgets':capjets,
                     'negative_perturbations':[str(v) for v in perturb],
                     'global_spatial_analyticity_assumed':False})


def positive_blend_and_heat_guards():
    rho=R(1,100000);d=R(1,20000);base=[1,2,20,400];diff=[]
    check('positive_analytic_displacement',R(9,2)*rho/(1-rho)<d)
    for a in range(4):
        ex=1/(1-a*rho)
        value=(ex-1)*base[a]+ex*2*S.factorial(a)*40**a*((1-40*d)**(-a-1)-1)
        diff.append(value)
    targets=[R(1,100),R(4,10),S.Integer(50),S.Integer(7000)]
    check('positive_analytic_difference',all(a<b for a,b in zip(diff,targets)),[str(S.N(v,15)) for v in diff])
    taper=[1,R(76,10),160,45000]
    full=[base[j]+sum(S.binomial(j,k)*taper[k]*targets[j-k] for k in range(j+1)) for j in range(4)]
    check('positive_global_first',full[1]/R(99,100)<5)
    check('positive_global_second',full[2]/R(99,100)<1000)
    check('positive_global_third',full[3]/R(99,100)<100000)
    check('complex_coefficient_phase',R(1,100)/R(99,100)<R(5,100))
    check('spatial_strip_margin',R(24,10)*rho/(1-rho)+1/((1-rho)*40)<R(5,100))
    check('scaled_mass_real_part',1-2*rho/(1-2*rho)>R(9999,10000))
    eff=(R(9999,10000)*324-100)/(324-100)
    check('heat_growth_rescaling',eff>R(9998,10000))
    check('uniform_heat_budget_enlargement',R(9998,10000)**(-6)<R(1002,1000))
    check('finite_part_scaling_margin',1/(1-4*rho)<R(1001,1000))
    check('auxiliary_budget_factor',2*10**19*R(1002,1000)*R(1001,1000)<4*10**19)
    FIXTURES.append({'fixture':'positive_blend_and_scaled_heat_guards',
                     'spatial_difference_bounds':[str(v) for v in diff],
                     'retained_complex_low_guards':[5,1000,100000],
                     'auxiliary_bound':'4e19 q^2 L(Y)'})


def exact_scaled_heat():
    tau,y,z,m,fp=S.symbols('tau y z m fp',positive=True)
    h=lambda x:1+x+x*x
    v=lambda x:2+x*x
    k2=lambda x:3-x+x**3
    # A subtracted heat remainder with a nonconstant local coefficient and spectral damping.
    old_heat=lambda t,x:1/t+(1+x)+t*k2(x)+t*t*h(x)*S.exp(-v(x)*t)
    scaled=S.expand(z*z*old_heat(z*z*tau,z*y))
    remainder=scaled-(1/tau+z*z*(1+z*y)+tau*z**4*k2(z*y))
    zero('full_subtracted_heat_scaling',remainder-z**6*tau*tau*h(z*y)*S.exp(-z*z*v(z*y)*tau))
    original=2*m*m*h(z*y)/(m*m+v(z*y))**3
    transformed=2*m*m*z**6*h(z*y)/(z*z*(m*m+v(z*y)))**3
    zero('proper_time_and_angular_scaling',original-transformed)
    zero('finite_k2_restoration',z**(-4)*(z**4*k2(z*y))*fp-fp*k2(z*y))
    zero('boundary_cubic_and_angular_profile',z**3*z**(-4)-1/z)
    beta,w=S.symbols('beta w',positive=True)
    zero('nonzero_winding_circle_scale',-(w*beta)**2/(4*z*z*tau)+(w*beta/z)**2/(4*tau))
    wrong=2*m*m*z**6*h(z*y)/(m*m+z*z*v(z*y))**3
    missing_angular=S.factor(S.diff(wrong-original,z).subs(z,1))
    check('omitted_angular_scaling_fails',missing_angular!=0,missing_angular)
    missing_local=S.factor(S.diff(z**4*fp*k2(z*y)-fp*k2(z*y),z).subs(z,1))
    check('omitted_finite_factor_fails',missing_local!=0,missing_local)
    # For a finite mode, differentiation of the scaled weight reproduces the Mellin term.
    weight=m*m*S.exp(-m*m*z*z*tau)
    zero('first_scaled_angular_Mellin_term',S.diff(weight,z).subs(z,1)-2*tau*S.diff(weight.subs(z,1),tau))
    FIXTURES.append({'fixture':'independent_scaled_heat_and_negative_controls',
                     'omitted_angular_scaling_error':str(missing_angular),
                     'omitted_finite_restoration_error':str(missing_local),
                     'fixture_is_actual_charged_kernel':False})


def remote_and_aggregate():
    lower=(R(65,24))**2*R(25,18)
    check('rational_exponential_certificate',lower>10)
    check('far_exponential_threshold',R(864,5)>74*R(7,3))
    check('near_exponential_threshold',R(351,5)>30*R(7,3))
    convolution=2*sum(S.factorial(10)/S.factorial(10-r)*(R(40,10**5*18))**r for r in range(1,11))
    check('extended_row_coefficient_convolution',convolution<R(1,1000),convolution)
    check('extended_row_cutoff_reserve',R(480,100000)<R(5,1000))
    check('conservative_row_absorption',R(6,5)*(convolution+R(480,100000))<1)
    near=1944*R(1,10**30)*(1+2/(R(351,5)-5))
    check('positive_defect_mass_squared_moment',near<R(21,10)*R(1,10**27),near)
    for j in range(1,9):
        bound=4*S.factorial(j)*100000**j*S.Integer(21)**(j-2)*(20*j+1)
        check('telescoped_inverse_jet_'+str(j),bound<S.factorial(j)*10**(7*j))
    H={p:6*18**p*R(1,10**74)*(1+2/(R(864,5)-p-3)) for p in range(2,11)}
    check('far_moments_are_decreasing',all(p+3<R(864,5) for p in H))
    A=[S.Integer(4)*10**19*1000**j for j in range(9)]
    G=[R(352,100000)*10**74*sum(stirling(j,m,kind=2)*(m+1)*S.factorial(m)*10**(4*m)*H[m+2]
        for m in range(j+1))/(S.factorial(j)*100**j) for j in range(9)]
    Q=[S.Integer(0)]+[S.Integer(10)**16*10**(5*j)*R(352,100000)*R(21,10**28) for j in range(1,9)]
    E=[A[j]+G[j]+Q[j] for j in range(9)]
    O=[sum(stirling(j,k,kind=1)*S.factorial(k)*100**k*E[k] for k in range(j+1))/(S.factorial(j)*100**j) for j in range(9)]
    check('complete_outer_Euler_eight',max(E)<10**44,[str(S.N(v,18)) for v in E])
    check('complete_outer_ordinary_eight',max(O)<10**44,[str(S.N(v,18)) for v in O])
    check('two_defects_remain_separate',G[8]>0 and Q[8]>0 and Q[0]==0)
    x=S.symbols('x');f=S.Function('f')(x);D=lambda value:x*S.diff(value,x)
    for j in range(1,9):
        falling=f
        for n in range(j):falling=D(falling)-n*falling
        zero('ordinary_Euler_identity_'+str(j),S.expand(falling-x**j*S.diff(f,x,j)))
    FIXTURES.append({'fixture':'complete_defect_and_outer_source_ledger',
                     'far_moments':{str(p):str(v) for p,v in H.items()},
                     'auxiliary_bounds':[str(v) for v in A],
                     'actual_cap_bounds':[str(v) for v in G],
                     'positive_dilation_bounds':[str(v) for v in Q],
                     'Euler_bounds':[str(v) for v in E],'ordinary_bounds':[str(v) for v in O],
                     'Euler_bounds_decimal':[str(S.N(v,18)) for v in E],
                     'ordinary_bounds_decimal':[str(S.N(v,18)) for v in O],
                     'full_angular_tower_retained':True,'complete_outer_reserve':'1e44'})


def independent_proper_time_integrals():
    records=[]
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            for rr,ii in [('0.000005','0'),('0','0.000005'),('-0.000005','0.000003')]:
                lam=mp.mpc(rr,ii);z=mp.exp(lam);y=mp.mpf('.08');m=mp.mpf(19);fp=mp.mpf(7)
                x=z*y;h=1+x+x*x;v=2+x*x;k2=3-x+x**3
                integrand=lambda t:m*m*mp.exp(-z*z*m*m*t)*z**6*t*t*h*mp.exp(-z*z*v*t)
                integral=mp.quad(integrand,[0,mp.mpf('.001'),mp.mpf('.01'),1,mp.inf])+z**(-4)*z**4*fp*k2
                exact=2*m*m*h/(m*m+v)**3+fp*k2
                error=abs(integral-exact)
                check('independent_scaled_integral_'+str(dps)+'_'+rr+'_'+ii,error<mp.mpf('1e-42'),mp.nstr(error,9))
                records.append({'dps':dps,'lambda_real':rr,'lambda_imag':ii,'absolute_error':mp.nstr(error,15)})
    FIXTURES.append({'fixture':'independent_complex_proper_time_integrals','records':records,
                     'actual_charged_source_evaluated':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 outer eight jets','input_count':len(ps)}),flush=True)
    try:
        provenance();finite_cap_family();positive_blend_and_heat_guards();exact_scaled_heat();remote_and_aggregate();independent_proper_time_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 outer eight jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_outer_Euler_and_ordinary_C8_bounded':True,'complete_outer_C8_reserve':'1e44',
         'cold_C8_bounded':False,'actual_metric_tenth_derivative_bounded':False,
         'actual_finite_metric_source_H6_bounded':False,'actual_transfer_certificate_passed':False,
         'original_outer_C3_contract_changed':False,'endpoint_C5_contract_changed':False,
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
