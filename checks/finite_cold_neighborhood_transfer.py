"""Local coefficient checks and conditional budgets; cold admission is unproved."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-cold-neighborhood-transfer'
PREVIOUS='finite-outer-admission-and-response'
OWN=['checks/finite_cold_neighborhood_transfer.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'digits':[60,90],'Cauchy_nodes':[64,128],'normalized_identity_tolerance':'1e-30',
          'finite_difference_steps':['1e-5','5e-6'],'finite_difference_tolerance':'1e-7',
          'refinement_ratio':'.35','tiny_error_floor':'1e-35','actual_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)

def local_majorants():
    old=saved(PREVIOUS)
    check('prior_outer_result_preserved',old['status']=='passed' and old['outer_canonical_admission_argument'])
    check('prior_hashes_preserved',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('raw_C5_comparison_is_only_allowance',F(120)*100**5*F('1.3e-4')==F('1.56e8'))
    T=S.symbols('T');tt=S.Rational('1.655e-12');J=(1+100*tt)*T+200*tt*T*T
    D=6*tt+100*tt*T+200*tt*T*T
    comp=S.Poly(S.series(sum((100*J)**k for k in range(6)),T,0,6).removeO(),T)
    seed=S.Poly(S.series(S.Rational('.027')*((1-32*(T+D))**-2-(1-32*T)**-2),T,0,6).removeO(),T)
    log=S.Poly(100*tt+2*sum((200*tt*T/(1-100*tt))**k/k for k in range(1,6)),T)
    records=[]
    for j in range(6):
        c=comp.nth(j)/100**j;ss=seed.nth(j)/(tt*100**j);ll=log.nth(j)/(tt*100**j)
        check('local_composition_'+str(j),c<S.Rational('1.000001'))
        check('seed_shift_full_tail_'+str(j),ss<20)
        check('Jacobian_'+str(j),ll<=100)
        check('local_new_endpoint_122_'+str(j),c+ss+ll<122)
        records.append({'order':j,'composition':str(c),'seed_per_t':str(ss),'log_per_t':str(ll)})
    check('chart_image_margin',F(7,32)+6*F('1.655e-12')<F(1,4))
    check('quadratic_map_injectivity',1-2*F('1.655e-12')-200*F('1.655e-12')*F(7,32)>F('.999999'))
    check('explicit_transfer_threshold',(F('2e-6')-F('1.000001e-6'))/F('5e-10')==F('1999.998'))
    check('hypothetical_1000_has_margin',F('1.000001e-6')+1000*F('5e-10')<F('2e-6'))
    check('complex_divided_cone_slack',F('1.3e-4')+F('5e-10')<F('.001')<F('.011'))
    FIXTURES.append({'fixture':'local_endpoint_majorants','records':records,'local_allowance':122,
                    'global_transport_target':'1999.998','global_transport_upper_proved':False})

def global_algebra():
    x,t=S.symbols('x t');nu=S.Function('nu')(x,t)
    # Noncommuting inverse variation identity checked with an independent 2x2 field.
    G=S.Matrix([[2+x+t,x*t],[t,3+x*x]])
    K=G.inv();A=G.diff(t)
    check('noncommuting_inverse_variation',S.simplify(K.diff(t)+K*A*K)==S.zeros(2))
    B=S.exp(-nu)*K
    identity=S.exp(-nu)*(K.diff(t)-S.diff(nu,t)*K)
    check('transported_B_variation',S.simplify(B.diff(t)-identity)==S.zeros(2))
    gradk=-(K.diff(x)*A*K+K*A.diff(x)*K+K*A*K.diff(x))
    check('connection_variation_all_terms',S.simplify(K.diff(t,x)-gradk)==S.zeros(2))
    gradB=S.exp(-nu)*(gradk-S.diff(nu,t,x)*K-S.diff(nu,t)*K.diff(x)-S.diff(nu,x)*(K.diff(t)-S.diff(nu,t)*K))
    check('transported_gradient_B_variation',S.simplify(B.diff(t,x)-gradB)==S.zeros(2))
    kval=3**2*1000;kg=2*210*1000*3+3**2*F('2e6')
    b=F('1.1')*(kval+100*3)
    bg=F('1.1')*(kg+21000*3+100*210+8*(kval+100*3))
    kx=2*1*1000*3+9*100000
    bx=F('1.1')*(kx+10000*3+100+10*(kval+100*3))
    check('global_operator_budget_values',kval==9000 and kg==19260000 and b==10230 and bg==21360240 and bx==1132010)
    check('global_operator_algebra_under_factor_two',b<2*10000 and bg<2*F('2e7') and bx<2*F('1e6'))
    check('conditional_transported_energy',F('1.5')*F('4.2')+(3*23)*F('1.3')/18<12<20)
    FIXTURES.append({'fixture':'global_transport_algebra','B':str(b),'gradient_B':str(bg),'center_B':str(bx),
                    'input_arrays_established_for_actual_family':False,'actual_far_seam_verified':False})

def conditional_budgets():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            a=mp.mpf('3.4e-15');rho=mp.mpf('2.5e-10');L=mp.mpf(91)
            c1=2*mp.e*a*L/rho;c2=16*mp.e*a*L**2/rho**2
            check('conditional_first_response_'+str(dps),c1<mp.mpf('.007'))
            check('conditional_mixed_response_'+str(dps),c2<mp.mpf('2e10'))
            check('weighted_logarithm_'+str(dps),80+mp.log(mp.mpf('5e-11')/a)<91)
            check('local_doubled_potential_cone_'+str(dps),mp.exp(mp.mpf('.025')+mp.mpf('2e-6')*mp.mpf('.00331'))-1<mp.mpf('.026'))
            inv=1/(mp.exp(-mp.mpf('.05'))-mp.mpf(2)/18)
            check('conditional_mass_inverse_margin_'+str(dps),inv<mp.mpf('1.191') and mp.exp(mp.mpf('.05'))*inv<mp.mpf('1.3'))
            records.append({'dps':dps,'conditional_first':mp.nstr(c1,25),'conditional_mixed':mp.nstr(c2,25),'conditional_mass_inverse':mp.nstr(inv,25)})
    check('local_doubled_curvature_cone',2*100**2*F('2e-6')*F('.00331')<F('.000133'))
    a=F('.00331');bolt=F('2e-28')
    check('Ward_first_row',F('3.4')*F('.007')+1400*a*F('3e-7')<F('.024'))
    check('Ward_mixed_row',F('3.4')*F('2e10')+2800*a*F('.007')+11000*a*a*F('3e-7')<F('6.9e10'))
    check('Ward_bolt_first',F('3.4')*F('.07')+1400*bolt*F('5e-4')<F('.24'))
    check('Ward_bolt_mixed',F('3.4')*F('2e11')+2800*bolt*F('.07')+11000*bolt*bolt*F('5e-4')<F('6.9e11'))
    m=F('1.003')*F('6.9e10')+F('4.02')*F('.024')+F('5.54')*F('1.1e-6')
    c=F('1.006')*F('6.9e10')+F('10.08')*F('.024')+F('26.5')*F('1.1e-6')
    b=F('1.006')*F('6.9e11')+F('10.08')*F('.24')+F('26.5')*F('.0017')
    check('conditional_mass_primitive',23*m<F('1.6e12'))
    check('conditional_absolute_mass',F('1e-6')*m<F('7e4'))
    check('conditional_absolute_lapse',(13122*c+972*b)/2<F('8e14'))
    check('conditional_lapse_remainder_not_Newton',F('8e14')*F('1e-10')**2/2==F('4e-6'))
    scale,rr,L1,L2=S.symbols('scale rr L1 L2',positive=True)
    check('normalization_invariant_response_products',S.simplify((L1/scale)*(scale*rr)-L1*rr)==0 and S.simplify((L2/scale**2)*(scale*rr)**2-L2*rr**2)==0)
    tt=S.symbols('tt',nonnegative=True)
    for j in [1,2]:check('quarter_power_log_derivative_'+str(j),S.simplify(-S.Rational(1,4)+S.Rational(j,2)/(91+tt/2)).subs(tt,0)<0)
    FIXTURES.append({'fixture':'conditional_response_and_constraint_budgets','records':records,
                    'mass_primitive':str(23*m),'mass_density':str(F('1e-6')*m),'lapse':str((13122*c+972*b)/2),
                    'cold_neighborhood_proved':False,'full_Newton_constant_proved':False})

def real_normalization():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            eps=mp.mpf('1e-7');dt=mp.mpf('1e-10');x0=mp.mpf('.04')
            seed=lambda x:mp.mpf('.02')*x*x/(1-x)
            raw=lambda x:abs(x-x0)**mp.mpf('5.5')
            rem=lambda x:eps*(raw(x)-raw(0)-mp.diff(raw,0)*x)
            perturb=lambda x:dt*(1+2*x+3*x*x)
            av=mp.exp(-perturb(0));cv=-av**2*mp.diff(perturb,0)/2
            H=lambda x:av*x+cv*x*x
            new=lambda x:seed(H(x))+rem(H(x))+perturb(H(x))+mp.log(av+2*cv*x)-seed(x)
            check('finite_normalizer_zero_jets_'+str(dps),abs(new(0))<mp.mpf('1e-40') and abs(mp.diff(new,0))<mp.mpf('1e-40'))
            errors=[]
            for x in [mp.mpf('-.08'),mp.mpf('.01'),mp.mpf('.1')]:
                # Independently recover the old coordinate and factor after transformation.
                oldx=H(x);inverse=2*oldx/(av+mp.sqrt(av*av+4*cv*oldx))
                recovered=new(inverse)+seed(inverse)-mp.log(av+2*cv*inverse)-perturb(oldx)
                err=abs(recovered-seed(oldx)-rem(oldx));errors.append(err)
                check('finite_endpoint_inverse_'+str(dps)+'_'+str(x),err<mp.mpf('1e-40'))
            # Exact fifth derivative translation of |x|^(11/2): sqrt(delta) modulus.
            coefficient=mp.fprod([mp.mpf('5.5')-j for j in range(5)])
            ratio=(coefficient*mp.sqrt(mp.mpf('1e-12'))/mp.mpf('1e-12'))/(coefficient*mp.sqrt(mp.mpf('1e-8'))/mp.mpf('1e-8'))
            check('false_C5_translation_Lipschitz_rejected_'+str(dps),abs(ratio-100)<mp.mpf('1e-40'))
            records.append({'dps':dps,'inverse_errors':[mp.nstr(e,15) for e in errors],'translation_Lipschitz_ratio_growth':mp.nstr(ratio,15)})
    FIXTURES.append({'fixture':'real_C5_non_C6_endpoint_normalization','records':records,'actual_metric_sampled':False})

def analytic_and_Ward():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            for hh in [mp.exp(-20),mp.exp(-79)]:
                amplitude=hh*(1+mp.log(mp.mpf('.01')/hh));freq=mp.log(mp.mpf('5e-11')/(mp.mpf('3.4e-14')*amplitude))/16
                center=mp.mpf('.17');fun=lambda z:mp.cos(freq*z)
                target=-freq**2*mp.cos(freq*center)
                contour=[]
                for count in PROTOCOL['Cauchy_nodes']:
                    radius=mp.mpf('.01');val=2*sum(fun(center+radius*mp.exp(2j*mp.pi*j/count))*mp.exp(-4j*mp.pi*j/count) for j in range(count))/(count*radius**2)
                    er=abs(val-target)/(1+abs(target));contour.append(er)
                    check('analytic_Cauchy_'+str(dps)+'_'+str(hh)+'_'+str(count),er<mp.mpf('1e-30'))
                gamma=mp.mpf('.03');cut=mp.mpf('.3')
                # Divide out the cold amplitude BEFORE the comparison.
                Q=lambda lam,y:mp.sqrt(y)*(1+mp.mpf('.2')*lam+mp.mpf('.1')*lam**2)
                W=lambda lam,y:y*y*mp.exp(2*lam*gamma*y*y)
                def pressure(lam):
                    return mp.quad(lambda z:mp.diff(lambda y:W(lam,y),z*z)*Q(lam,z*z)*2*z,[0,mp.sqrt(cut)])/(2*W(lam,cut))
                exact=mp.diff(pressure,center,2);errs=[]
                for step in map(mp.mpf,PROTOCOL['finite_difference_steps']):
                    val=(pressure(center+step)-2*pressure(center)+pressure(center-step))/step**2
                    errs.append(abs(val-exact)/(1+abs(exact)))
                check('Ward_second_difference_'+str(dps)+'_'+str(hh),errs[0]<mp.mpf('1e-7') and (errs[1]<mp.mpf('.35')*errs[0] or max(errs)<mp.mpf('1e-35')))
                frozen=mp.quad(lambda z:mp.diff(lambda y:W(center,y),z*z)*mp.mpf('.2')*z*2*z,[0,mp.sqrt(cut)])/(2*W(center,cut))
                check('frozen_Ward_kernel_rejected_'+str(dps)+'_'+str(hh),abs(exact-frozen)>mp.mpf('1e-5'))
                records.append({'dps':dps,'h':mp.nstr(hh,12),'Cauchy_errors':[mp.nstr(e,12) for e in contour],'Ward_errors':[mp.nstr(e,12) for e in errs]})
    FIXTURES.append({'fixture':'normalized_analytic_and_Ward_mechanisms','records':records,'actual_source_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite cold neighborhood transfer','input_count':len(ps)}),flush=True)
    try:local_majorants();global_algebra();conditional_budgets();real_normalization();analytic_and_Ward()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 finite cold neighborhood transfer','created_utc':datetime.now(timezone.utc).isoformat(),
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
