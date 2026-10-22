"""Local cold gluing and finite-centered canonical admission; not a full metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-cold-local-gluing-and-admission'
PREVIOUS='fixed-domain-cold-far-seam'
OWN=['checks/finite_cold_local_gluing_and_admission.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'digits':[60,90],'identity_tolerance':'1e-35','finite_difference_steps':['1e-4','5e-5'],'finite_difference_tolerance':'1e-7','refinement_ratio':'.3','actual_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)




def kernel_and_graph():
    old=saved(PREVIOUS)
    check('far_result_preserved',old['status']=='passed' and old['fixed_domain_distant_seam_argument'])
    check('prior_hashes_preserved',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    delta=F('.026');m=18
    check('weighted_kernel_L1',F(36,23)+F(22,7)*F(1296,23*33)<7)
    check('weighted_Bessel_decay',F(7,10)-F(13,36)>F(1,3))
    check('weighted_L43_squared',12/(F(17,10)*F(1,3))<22<56)
    check('quadratic_kernel_L2_squared',24/(3*F(2,3)**4)<100)
    check('Neumann_inverse',7/(1-7*delta)<9)
    check('complete_nonfree_row',56*delta/(1-7*delta)<F('1.8'))
    check('normal_cap_row',2*F('1.02')*33*(F('1.5')+F('1.8'))<400)
    check('quadratic_insertion',10+100*F(9,16)<70)
    c=2*F('1.02')*70*9*33/(400*18**2)
    e=F(5,32)**2/(2*(1-delta));total=c+e
    check('first_row_Frechet',total<F('.4'))
    second=2*c*F(9,16)+2*c*e+F(3,4)*F(5,32)**4/(1-delta)**2
    check('second_row_Frechet',second<4)
    check('normalizer_second_ledger',20001+400*2+2*16*27+200+1<24000)
    check('normalizer_mixed_ledger',400*200+2*16*3000+22000==198000)
    check('normalizer_inner_map',101*F(1,8)+11000*F(1,8)**2<200 and 101+22000*F(1,8)<3000)
    check('proper_weight_gradient',F('1.05')/F('.55')<2<3)
    check('transported_energy',F('1.5')*F('4.2')+3*23*F(2,18)<20)
    check('cap_second_derivative',6*(2+15000)+F('3.75')*(1+30000)+30000<300000)
    check('cap_second_increment',300000*F('5e-10')==F('.00015'))
    check('middle_first_coordinate',2*20*20000==800000<1000000)
    mixed=2*F(1,2)*(40*2+10*20)*20000+20*2*2000000+1000
    check('middle_mixed_coordinate',mixed<100000000)
    aa=F('1e-27')
    check('common_frame_row_first',F('.4')*20000+100<10000)
    check('common_frame_row_center',F('.4')*10000000+2000<10000000)
    check('common_frame_row_mixed',F('.4')*3000000+100000+4*20000*10000000*aa+F('1e15')*aa<2000000)
    cglue=4*(2*2000000*2+2*10000*10+100000000+(2*10000000*1000000+2*10000*10000000*2)*aa)
    check('complete_mixed_gluing_coefficient',cglue<F('1e9'))
    check('flat_remainder_budget',F('1.2e21')+F('1e20')+F('1e20')<F('1e22'))
    FIXTURES.append({'fixture':'kernel_row_and_graph_ledgers','row_first_constant_per_quadratic_tag':str(total),'row_second_constant_per_product':str(second),'mixed_middle_coefficient':str(mixed),'mixed_gluing_coefficient':str(cglue),'row_R':400,'C_glue':1000000000})

def mixed_endpoint_majorants():
    T=S.symbols('T');tt=S.Rational('1.655e-12');J=(1+100*tt)*T+200*tt*T*T;D=6*tt+100*tt*T+200*tt*T*T
    HX=6*tt+100*tt*T+300*tt*T*T
    base=sum((100*T)**k for k in range(6));comp=S.Poly(S.series(base.subs(T,J),T,0,5).removeO(),T)
    extra=S.Poly(S.series(S.diff(base,T).subs(T,J)*HX,T,0,5).removeO(),T)
    seedx=S.Poly(S.series(S.Rational('.05')*((1-32*(T+D))**-2-(1-32*T)**-2),T,0,5).removeO(),T)
    seedz=S.Poly(S.series(S.Rational('.027')*64*(1-32*(T+D))**-3*HX,T,0,5).removeO(),T)
    log=S.Poly(S.series((50*tt+220*tt*T)/(1-100*tt-200*tt*T),T,0,5).removeO(),T)
    records=[]
    for j in range(5):
        c=comp.nth(j)/100**j;ex=extra.nth(j)/(tt*100**j);sx=seedx.nth(j)/(tt*100**j);sz=seedz.nth(j)/(tt*100**j);lg=log.nth(j)/(tt*100**j)
        check('mixed_composition_'+str(j),c<S.Rational('1.000001'))
        check('mixed_finite_extra_'+str(j),ex<10000 and ex*S.Rational('3.31e-9')<S.Rational('.001'))
        check('mixed_seed_shift_and_chain_'+str(j),sx+sz<60)
        check('mixed_log_'+str(j),lg<100)
        check('complete_mixed_allowance_'+str(j),c+ex*S.Rational('3.31e-9')+ex*tt+sx+sz+lg<500)
        records.append({'order':j,'composition':str(c),'extra_multiplier':str(ex),'seed_center_shift':str(sx),'seed_chain':str(sz),'log':str(lg)})
    for j,bound in [(1,1),(2,16),(3,400),(4,20000)]:
        seed=S.factorial(j)*(j+1)*S.Rational('.027')*8**j
        finite=S.Rational('1e-6')*S.Rational('.00331')*S.factorial(j)*100**j
        check('inner_background_jet_'+str(j),seed+finite<bound)
    check('endpoint_tube',F('1.000001e-6')+1000*F('5e-10')==F('1.500001e-6')<F('2e-6'))
    check('preserved_background_plus_disk',2*F('2.5e-10')==F('5e-10'))
    check('independent_complex_cone',F('1.3e-4')+F('5e-10')<F('.001')<F('.011'))
    check('solver_embedding_still_required',F('2.5e-10')/F('1e-8')==F('.025'))
    FIXTURES.append({'fixture':'finite_mixed_endpoint_contract','records':records,'mixed_allowance':500,'source_transfer_allowance':1000,'real_endpoint_radius':'5e-10','cold_background_radius':'2.5e-10','full_metric_ball_admitted':False})

def numerical_kernel_and_response():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            a=mp.mpf(13)/36;A=1-a;B=1+a
            transformed=mp.quad(lambda y:4/(A+B*y*y)**2,[0,mp.mpf('.5'),1])
            closed=1/A+2*mp.atan(mp.sqrt(B/A))/(A*mp.sqrt(A*B))
            check('kernel_L1_independent_integral_'+str(dps),abs(transformed-closed)<mp.mpf('1e-40') and closed<7)
            bessel=[]
            for s in [mp.mpf(1),mp.mpf(3)]:
                for order in [0,1]:
                    val=mp.quad(lambda t:mp.exp(-s*mp.cosh(t))*mp.cosh(order*t),[0,1,3,8])
                    err=abs(val-mp.besselk(order,s))/abs(mp.besselk(order,s));bessel.append(err)
                    check('positive_Bessel_representation_'+str(dps)+'_'+str(s)+'_'+str(order),err<mp.mpf('1e-35'))
            inv=mp.exp(mp.mpf('.043'))/(mp.exp(-mp.mpf('.043'))-mp.mpf(6)/18)
            check('transported_weight_inverse_'+str(dps),inv<2)
            check('quadratic_tag_exponentiation_'+str(dps),mp.exp(mp.mpf('.026'))*12000<20000)
            new=mp.mpf('1e9')*400**2*48**3*mp.exp(-3)
            old=mp.mpf('1e10')*400**2*48**3*mp.exp(-3)
            check('new_local_budget_'+str(dps),new<mp.mpf('9e17')<mp.mpf('1e20'))
            check('retained_path_budget_'+str(dps),old<mp.mpf('9e18')<mp.mpf('1e20'))
            first=2*mp.e*mp.mpf('3.4e-15')*91/mp.mpf('2.5e-10');mixed=16*mp.e*mp.mpf('3.4e-15')*91**2/mp.mpf('2.5e-10')**2
            check('admitted_first_response_'+str(dps),first<mp.mpf('.007'))
            check('admitted_mixed_response_'+str(dps),mixed<mp.mpf('2e10'))
            check('weighted_log_'+str(dps),80+mp.log(mp.mpf('5e-11')/mp.mpf('3.4e-15'))<91)
            records.append({'dps':dps,'weighted_kernel_L1':mp.nstr(closed,30),'Bessel_relative_errors':[mp.nstr(e,12) for e in bessel],'transported_weight_inverse':mp.nstr(inv,30),'new_local_M3_coefficient':mp.nstr(new,30),'retained_path_M3_coefficient':mp.nstr(old,30),'first_response':mp.nstr(first,30),'mixed_response':mp.nstr(mixed,30)})
    FIXTURES.append({'fixture':'kernel_comparison_and_admitted_response','records':records,'angular_cutoff_used':False,'actual_correction_norm_proved':False})

def real_C5_normalizer():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            d=mp.mpf('1e-10');join=mp.mpf('.04')
            raw=lambda z:abs(z-join)**mp.mpf('5.5')
            rem=lambda z:raw(z)-raw(0)-mp.diff(raw,0)*z
            background=lambda x,z:mp.mpf('.02')*z*z/(1-z)+mp.mpf('.001')*x*z*z+mp.mpf('1e-7')*(1+mp.mpf('.2')*x)*rem(z)
            direction=lambda x,z:d*(1+mp.mpf('.1')*x+2*z+3*z*z+raw(z))
            def normal(t,x,z):
                av=mp.exp(-t*direction(x,0));cv=-av*av*t*mp.diff(lambda zz:direction(x,zz),0)/2
                H=av*z+cv*z*z
                return background(x,H)+t*direction(x,H)+mp.log(av+2*cv*z)
            tt=mp.mpf('.31');xx=mp.mpf('.12');zz=mp.mpf('.07')
            check('normalized_value_'+str(dps),abs(normal(tt,xx,0))/d<mp.mpf('1e-35'))
            check('normalized_gradient_'+str(dps),abs(mp.diff(lambda z:normal(tt,xx,z),0))/d<mp.mpf('1e-35'))
            exact=mp.diff(lambda t:mp.diff(lambda x:normal(t,x,zz),xx),tt)/d
            errors=[]
            for step in map(mp.mpf,PROTOCOL['finite_difference_steps']):
                approx=(normal(tt+step,xx+step,zz)-normal(tt+step,xx-step,zz)-normal(tt-step,xx+step,zz)+normal(tt-step,xx-step,zz))/(4*step*step*d)
                errors.append(abs(approx-exact)/(1+abs(exact)))
            check('mixed_normalizer_real_refinement_'+str(dps),errors[0]<mp.mpf('1e-7') and (errors[1]<mp.mpf('.3')*errors[0] or max(errors)<mp.mpf('1e-35')))
            records.append({'dps':dps,'normalized_mixed_value':mp.nstr(exact,30),'normalized_errors':[mp.nstr(e,15) for e in errors],'unknown_field_is_spatially_analytic':False})
    FIXTURES.append({'fixture':'real_C5_non_C6_mixed_normalizer','records':records,'actual_metric_sampled':False})

def noncommuting_mixed_product():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            p1=mp.matrix([[0,1],[1,0]]);p2=mp.matrix([[0,-1j],[1j,0]]);p3=mp.matrix([[1,0],[0,-1]])
            A=3*mp.eye(2)+mp.mpf('.2')*1j*p1+mp.mpf('.3')*1j*p3;B=mp.mpf('.04')*1j*p1;C=mp.mpf('.03')*1j*p2;E=mp.mpf('.02')*1j*p3
            S0=A**-1;St=-S0*B*S0;Sx=-S0*C*S0;Stx=S0*C*S0*B*S0+S0*B*S0*C*S0-S0*E*S0
            left=lambda t,x:mp.matrix([[1+t+mp.mpf('.2')*x,mp.mpf('.3')+x*t]])
            right=lambda t,x:mp.matrix([[mp.mpf('.7')+mp.mpf('.1')*t],[1-mp.mpf('.2')*x+mp.mpf('.05')*x*t]])
            l=left(0,0);r=right(0,0);lt=mp.matrix([[1,0]]);lx=mp.matrix([[mp.mpf('.2'),0]]);ltx=mp.matrix([[0,1]])
            rt=mp.matrix([[mp.mpf('.1')],[0]]);rx=mp.matrix([[0],[-mp.mpf('.2')]]);rtx=mp.matrix([[0],[mp.mpf('.05')]])
            terms=[ltx*S0*r,lt*Sx*r,lt*S0*rx,lx*St*r,l*Stx*r,l*St*rx,lx*S0*rt,l*Sx*rt,l*S0*rtx]
            formula=sum(z[0] for z in terms)
            fun=lambda t,x:(left(t,x)*(A+t*B+x*C+t*x*E)**-1*right(t,x))[0]
            exact=mp.diff(lambda t:mp.diff(lambda x:fun(t,x),0),0)
            err=abs(formula-exact)/(1+abs(exact))
            check('noncommuting_nine_term_identity_'+str(dps),err<mp.mpf('1e-35'))
            check('direct_mixed_insertion_not_zero_'+str(dps),abs((l*S0*E*S0*r)[0])>mp.mpf('1e-5'))
            errors=[]
            for step in map(mp.mpf,PROTOCOL['finite_difference_steps']):
                value=(fun(step,step)-fun(step,-step)-fun(-step,step)+fun(-step,-step))/(4*step*step)
                errors.append(abs(value-exact)/(1+abs(exact)))
            check('nine_term_difference_refinement_'+str(dps),errors[0]<mp.mpf('1e-7') and errors[1]<mp.mpf('.3')*errors[0])
            records.append({'dps':dps,'nine_term_error':mp.nstr(err,15),'difference_errors':[mp.nstr(e,15) for e in errors]})
    FIXTURES.append({'fixture':'noncommuting_mixed_resolvent_and_rows','records':records,'actual_quantum_state_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite cold local gluing and admission','input_count':len(ps)}),flush=True)
    try:kernel_and_graph();mixed_endpoint_majorants();numerical_kernel_and_response();real_C5_normalizer();noncommuting_mixed_product()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 finite cold local gluing and admission','created_utc':datetime.now(timezone.utc).isoformat(),
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
