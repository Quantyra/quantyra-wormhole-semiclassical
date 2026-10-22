"""Geometric budgets and source thresholds; the cold quantum-source proof remains open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-cold-admission-and-response'
PREVIOUS='finite-cold-neighborhood-transfer'
OWN=['checks/finite_cold_admission_and_response.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'digits':[60,90],'matrix_series_order':24,'identity_tolerance':'1e-40','finite_difference_steps':['1e-5','5e-6'],'actual_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)


def geometric_budgets():
    old=saved(PREVIOUS)
    check('prior_conditional_result_preserved',old['status']=='passed' and not old['general_cold_family_admitted'])
    check('prior_hashes_preserved',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    r=F(7,32);a=F('.00331');t=F('5e-10')*a
    check('tiny_endpoint_jet',t==F('1.655e-12'))
    check('normalizer_center_value',101*r+11000*r*r<600)
    check('normalizer_center_gradient',101+22000*r<5000)
    check('composed_map_value',F('1.052')*6<F('6.4'))
    check('composed_map_first',F('2.7')*6+F('1.052')*46<66)
    check('composed_map_second',310*6+2*F('2.7')*46+F('1.052')*200<2400)
    check('composed_center_first',121000*a*6+F('2.7')*600+1401*a*46+F('1.052')*5000<10000)
    check('cutoff_first',66+30*F('6.4')<420)
    check('cutoff_second',2400+60*66+1500*F('6.4')<30000)
    check('cutoff_center_first',10000+30*670<50000)
    check('inverse_gradient',9*23<210)
    check('coordinate_conformal_gradient',500*F('1.5')*11+5*420<21000)
    check('coordinate_conformal_value',5*11<100)
    check('shifted_seed_margin',(F('.018')/(F('.018')-F('2.1e-7')))**3<F('1.0001'))
    kval=9*1000;kg=2*210*1000*3+9*F('2e6')
    b=F('1.1')*(kval+100*3)
    bg=F('1.1')*(kg+21000*3+100*210+8*(kval+100*3))
    bx=F('1.1')*(2*1000*3+9*100000+10000*3+100+10*(kval+100*3))
    check('transported_coefficient_algebra',b==10230 and bg==21360240 and bx==1132010)
    FIXTURES.append({'fixture':'actual_map_budget','B':str(b),'gradient_B':str(bg),'center_B':str(bx),'complete_source_transfer_proved':False})

def cap_polynomial():
    t=S.symbols('t');v0,v1,v2=S.symbols('v0 v1 v2')
    eta=1-10*t**3+15*t**4-6*t**5
    poly=eta*(v0+v1*t+v2*t*t/2)
    for j,target in enumerate([v0,v1,v2]):
        check('cap_inner_join_'+str(j),S.diff(poly,t,j).subs(t,0)==target)
        check('cap_outer_join_'+str(j),S.diff(poly,t,j).subs(t,1)==0)
    check('cutoff_derivative_sign',S.factor(S.diff(eta,t))==-30*t*t*(t-1)**2)
    # Exact Bernstein coefficients on 16 equal intervals bound the whole polynomial.
    p=poly.subs({v0:S.Rational('3e-9'),v1:S.Rational('4e-8'),v2:S.Rational('2e-7')})
    z=S.symbols('z');maximum=S.Rational(0)
    for i in range(16):
        q=S.Poly(p.subs(t,(i+z)/16),z);n=q.degree()
        bb=[sum(q.nth(k)*S.binomial(j,k)/S.binomial(n,k) for k in range(j+1)) for j in range(n+1)]
        maximum=max(maximum,*bb)
    check('whole_cap_value_Bernstein',maximum<S.Rational('3e-8'))
    deriv=F(15,8)*(F('3e-9')+F('4e-8')+F('1e-7'))+F('4e-8')+F('2e-7')
    check('whole_cap_gradient',deriv<F('6e-7'))
    check('cap_value_cone',F('.042')+F('3e-8')+F('5e-10')<F('.043')<F('.05'))
    check('cap_gradient_cone',F('.4')+F('6e-7')+F('5e-10')<F('.401')<F('.5'))
    FIXTURES.append({'fixture':'explicit_cap','value_Bernstein_upper':str(maximum),'gradient_upper':str(deriv),'global_domain_identification_claimed':False})

def source_thresholds():
    check('transport_target',(F('2e-6')-F('1.000001e-6'))/F('5e-10')==F('1999.998'))
    check('candidate_1000_fits',F('1.000001e-6')+1000*F('5e-10')<F('2e-6'))
    check('interior_complex_tube',F('1.3e-4')+F('5e-10')<F('.001')<F('.011'))
    b0=F('1e4')*F('2e-6')*F('.00331');b1=F('2e7')*F('2e-6')*F('.00331')
    ratio=13*b0+b1*F('1.5')/36+F(3,4)*b0*F('1.5')
    check('conditional_graph_ratio',ratio<F('.007'))
    check('far_envelope_M_monotonicity',F(3,18)-F('1.9')/2<0)
    check('far_envelope_D_monotonicity',18/2-1>0)
    check('angular_inverse_cubic',12*F(1,18)==F(2,3))
    check('Newton_design_self_map',F('5e-9')+F('.25')*F('1e-8')+F('5e7')*F('1e-8')**2/2<=F('1e-8'))
    check('Newton_design_contraction',F('.25')+F('5e7')*F('1e-8')==F('.75'))
    check('cold_embedding_threshold',F('2.5e-10')/F('1e-8')==F('.025'))
    check('interior_embedding_threshold',(F('.001')-F('1.3e-4'))/F('1e-8')==87000)
    check('outer_embedding_threshold',F('4.9e-8')/F('1e-8')==F('4.9'))
    check('conditional_lapse_remainder',F('8e14')*F('1e-10')**2/2==F('4e-6'))
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            target0=mp.mpf('1e20')*mp.exp(mp.mpf('15.2'))/18**3
            target1=mp.mpf('1e20')*mp.exp(480)/18**3
            product=mp.mpf('1e20')*mp.exp(3)/48**3
            row=mp.sqrt(product/mp.mpf('1e10'))
            check('suggested_far_value_prefactor_fits_'+str(dps),mp.mpf('5e20')<target0)
            check('suggested_far_gradient_prefactor_fits_'+str(dps),mp.mpf('1e30')<target1)
            check('suggested_local_row_fits_'+str(dps),800<row)
            check('local_envelope_extremum_'+str(dps),abs(mp.diff(lambda x:x**3*mp.exp(-x/16),48))<mp.mpf('1e-40') and 48**3*mp.exp(-3)<5600)
            check('cold_distance_guard_'+str(dps),mp.mpf('.946')*mp.log(10)-mp.mpf('.25')>mp.mpf('1.928'))
            check('near_bolt_distance_guard_'+str(dps),mp.acosh(mp.mpf('.01')/mp.mpf('1e-30'))-mp.mpf('1.25')>60)
            check('local_potential_cone_'+str(dps),mp.exp(mp.mpf('.025')+mp.mpf('2e-6')*mp.mpf('.00331'))-1<mp.mpf('.026'))
            first=2*mp.e*mp.mpf('3.4e-15')*91/mp.mpf('2.5e-10')
            mixed=16*mp.e*mp.mpf('3.4e-15')*91**2/mp.mpf('2.5e-10')**2
            check('conditional_response_'+str(dps),first<mp.mpf('.007') and mixed<mp.mpf('2e10'))
            records.append({'dps':dps,'T0_threshold':mp.nstr(target0,35),'T1_threshold':mp.nstr(target1,35),'local_product_threshold':mp.nstr(product,35),'local_row_threshold_at_1e10':mp.nstr(row,35),'conditional_first':mp.nstr(first,25),'conditional_mixed':mp.nstr(mixed,25)})
    FIXTURES.append({'fixture':'missing_source_thresholds','records':records,'conditional_graph_ratio':str(ratio),'T0_known_upper':None,'T1_known_upper':None,'new_chord_row_upper':None,'full_cold_admission_proved':False})

def fixed_disk_and_frame():
    z,b=S.symbols('z b');T=(z+b)/(1+b*z)
    check('disk_transvection_generator',S.diff(T,b).subs(b,0)==1-z*z)
    check('positive_axis_boundary_not_fixed',S.factor(T-z)==-b*(z-1)*(z+1)/(b*z+1))
    gamma1=S.Matrix([[0,1],[1,0]]);gamma2=S.Matrix([[0,-S.I],[S.I,0]]);spin=gamma1*gamma2
    check('spin_frame_first_anticommutator',gamma1*spin+spin*gamma1==S.zeros(2))
    check('spin_frame_second_anticommutator',gamma2*spin+spin*gamma2==S.zeros(2))
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            for h in [mp.exp(-79),mp.mpf('1e-30')]:
                R=mp.exp(-h);theta=mp.mpf('.37');zz=R*mp.exp(1j*theta)
                deriv=mp.diff(lambda bb:(zz+bb)/(1+bb*zz),0)
                actual=mp.re(mp.exp(-1j*theta)*deriv)
                target=-mp.expm1(-2*h)*mp.cos(theta)
                error=abs(actual-target)/abs(target)
                # Cancellation in 1-R^2 costs about 35 digits; use precision-aware tolerance.
                tol=mp.power(10,-dps+37)
                check('finite_boundary_normal_velocity_'+str(dps)+'_'+mp.nstr(h,8),target>0 and error<tol)
                records.append({'dps':dps,'h':mp.nstr(h,30),'normal_velocity':mp.nstr(target,30),'relative_error':mp.nstr(error,12),'relative_tolerance':mp.nstr(tol,8)})
    FIXTURES.append({'fixture':'finite_disk_domain_diagnostic','records':records,'automatic_H2_domain_identification_rejected':True,'physical_nonexistence_claimed':False})

def matrix_chord():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            A=mp.matrix([[2,mp.mpf('.2')],[0,3]])
            E=mp.matrix([[mp.mpf('.0002'),mp.mpf('.0003')],[-mp.mpf('.0001'),mp.mpf('.0004')]])
            I=mp.eye(2);L=mp.logm(I+E)
            # Independent convergent log series, retaining matrix ordering.
            series=mp.zeros(2);power=I
            for k in range(1,25):
                power=power*E;series+=(-1)**(k+1)*power/k
            err=mp.norm(L-series)
            check('noncommuting_log_series_'+str(dps),err<mp.mpf('1e-40'))
            chord=lambda t:A*mp.expm(t*L)
            exact=chord(mp.mpf('.4'))*L
            fd=[]
            for dt in map(mp.mpf,PROTOCOL['finite_difference_steps']):
                derivative=(chord(mp.mpf('.4')+dt)-chord(mp.mpf('.4')-dt))/(2*dt)
                fd.append(mp.norm(derivative-exact)/mp.norm(exact))
            check('chord_tangent_refinement_'+str(dps),fd[0]<mp.mpf('1e-15') and fd[1]<mp.mpf('.26')*fd[0])
            check('chord_order_matters_'+str(dps),mp.norm(A*L-L*A)>mp.mpf('1e-5'))
            records.append({'dps':dps,'log_series_error':mp.nstr(err,15),'normalized_tangent_errors':[mp.nstr(e,15) for e in fd]})
    FIXTURES.append({'fixture':'noncommuting_chord_mechanism','records':records,'actual_operator_row_estimates_verified':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite cold source feasibility audit','input_count':len(ps)}),flush=True)
    try:geometric_budgets();cap_polynomial();source_thresholds();fixed_disk_and_frame();matrix_chord()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 finite cold source feasibility audit','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if isinstance(v,bool) and k not in run:run[k]=v
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
