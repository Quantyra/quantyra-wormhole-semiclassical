"""Cold gap Hessian and tau-retaining products in the actual scalar residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-scalar-cold-gap-product';PREVIOUS='finite-scalar-residual-six-refined'
OWN=['checks/finite_scalar_cold_gap_product.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_refined_scalar',ROOT/'checks/finite_scalar_residual_six_refined.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R);L=R.L
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[]
PROTOCOL={'scalar_source_order':6,'gap_Hessian_order':2,'tau_gap_jet_order':8,'numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','actual_charged_source_sampled':False,'ordinary_continuum_proof_separate':True}
def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def zero(name,v):
    v=S.simplify(v);check(name,v==0,v)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def strings(a):return [str(x) for x in a]
def decimals(a):return [L.decimal(x) for x in a]
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
    old=saved(PREVIOUS)
    check('previous_refinement_passed',old['status']=='passed' and len(old['checks'])==43 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('scope_is_scalar_estimate',review()['actual_scalar_residual_H6_bounded'] and not review()['actual_transfer_certificate_passed'])
    FIXTURES.append({'fixture':'retained_sources_and_claims','previous_checks':43,'previous_inputs':577,'old_sources_edited':False})

def gap_equation_and_cover():
    x=S.symbols('rho');tau=S.Function('tau')(x);v=S.Function('V')(x);beta,ell,cs,src=S.symbols('beta L c source')
    nu=tau*v
    original=S.diff(nu,x,2)+(beta+2*ell)*S.diff(nu,x)-cs*nu-src
    divided=S.diff(v,x,2)+beta*S.diff(v,x)+(2*ell+2*S.diff(tau,x)/tau)*S.diff(v,x)
    divided+=(S.diff(tau,x,2)/tau+(beta+2*ell)*S.diff(tau,x)/tau-cs)*v-src/tau
    zero('exact_divided_radius_equation',original-tau*divided)
    check('cold_curvature_lower',1-3*F('1.334')*F('.001')-6*F('.501')*F('.001')**2-F('1e-30')>F('.99'))
    check('cold_curvature_upper',1+3*F('1.334')*F('.001')+6*F('.501')*F('.001')**2+F('1e-30')<F('1.1'))
    check('cold_bz_relative',2+3*F('1.334')*F('.001')+4*F('.501')*F('.001')**2+F('1e-28')<F('2.1'))
    check('global_gap_speed',F(6)/F('5.9')+F('1e-10')+F('1e-24')<F('1.01')**2)
    check('cold_L_and_S',F('1.01')*F('.001')<F('.0011') and F('2.1')*F('1.001')+2*F('1.01')**2*F('.001')*F('1.001')**2<3)
    check('cold_radius_potential',F('8.1')/6+4*F('1.001')**2/F('5.9')+F('1e-20')<10)
    check('cold_relative_Ltau',F('2.1')+2*F('.0011')*F('1.01')<3)
    check('cold_V_radial_gradient',F('1.001')**2*F('8e-7')**2<F('5.9')*F('3.4e-7')**2)
    source=F('6.8e-28')*82*(F('2e20')+F('1e22'))
    force=(2*F('.0011')+2*F('1.01'))*F('3.4e-7')+13*F('1.2e-8')+6*F('1.2e-8')+F('82e-10')+source
    check('cold_complete_gap_forcing',force<F('.0006'),L.decimal(force))
    check('regular_integral_factor',F('2.1')/F('.99')<F('2.13'))
    check('cold_full_Hessian',F('5.26')*F('.0006')<F('.0032')<F('.004'))
    t=S.symbols('t',positive=True)
    zero('middle_weight_extrema',S.diff(t**3/(t-1),t)-t*t*(2*t-3)/(t-1)**2)
    check('middle_weight_endpoint_values',all(a**3/(a-1)<2600 for a in [F('1.001'),F('50.1')]))
    product=F('3.432e-6')+11*F('1.2e-8');check('middle_angular_slot',product==F('3.564e-6'))
    mid=(2600*(F('1.1e-5')+product)+F('50.1')**2*2*F('2.24')*F('8e-7'))/F('5.9')
    check('middle_weighted_Hessian',mid<F('.008'),L.decimal(mid))
    check('middle_unweighted_Hessian',(1001*(F('1.1e-5')+product)+2*F('2.24')*F('8e-7'))/F('5.9')<F('.0025'))
    check('outer_weighted_Hessian',(F('2.1e-6')+(2+2*F('2.24'))*F('4.3e-8'))/F('5.9')<F('1e-6'))
    check('global_weighted_values',16*F('1.2e-8')<F('2e-7') and F('8e-9')<F('2e-7'))
    check('u_weighted_gradient',(40*F('3.2e-7'))**2<F('5.9')*F('6e-6')**2 and F('3e-6')<F('6e-6'))
    check('V_nu_core_weighted_gradient',(16*F('8e-7'))**2<F('5.9')*F('6e-6')**2)
    check('nu_outer_weighted_gradient',F('4.3e-8')+F('2.24')/4*F('2e-9')<F('5e-8'))
    cold_lap=F('6.1')*(6*F('1.2e-8')+10*F('1.2e-8')+F('82e-10')+source)
    check('cold_weighted_Lnu_over_tau',cold_lap<F('.004'),L.decimal(cold_lap))
    check('middle_weighted_Lnu_over_tau',2600*F('6.3e-7')<F('.005'))
    check('outer_nu_second',F('2.1e-6')+2*F('2.24')/40*F('4.3e-8')+21*F('2e-9')<F('2.2e-6'))
    check('outer_weighted_Lnu_over_tau',(F('2.2e-6')+F('6.48')*F('5e-8'))/F('.975')<F('3e-6'))
    FIXTURES.append({'fixture':'global_gap_Hessian_and_low_cover','cold_forcing_exact':str(force),
                     'cold_weighted_Lnu_over_tau':str(cold_lap),'middle_weighted_Hessian_exact':str(mid),
                     'unweighted_Hessian':'.004','weighted_Hessian':'.01','weighted_values':'2e-7',
                     'weighted_gradients':'6e-6','weighted_Lnu_over_tau':'.005','inverse_thermal_gap_used':False})

def low_and_tau_arrays():
    pf=fixture(PREVIOUS,'separate_regional_tangent_orders')
    U,K,V,KV=[poly([F(x) for x in pf[name]]) for name in ['U','K','Vbar','KV']]
    U[0]=min(U[0],F('1.2e-8'));V[0]=min(V[0],F('1.2e-8'))
    K[0]=min(K[0],F('2e-7'));KV[0]=min(KV[0],F('2e-7'))
    K[1]=min(K[1],F('6e-6'));KV[1]=min(KV[1],F('6e-6'))
    V[2]=min(V[2],F('.002'));KV[2]=min(KV[2],F('.005'))
    mid=fixture('middle-metric-ten-jets','radial_two_derivative_reconstruction')
    cold=fixture('cold-metric-ten-jets','weighted_radial_reconstruction')
    A=[F(x) for x in fixture('cold-metric-ten-jets','original_finite_completion_and_tensor_conversion')['Cartesian_amplifications'][:9]]
    rm=[F('1.2')*sum(F('.2')**i*F(mid['N'][j-i]) for i in range(j+1)) for j in range(9)]
    rc=[F('1.2')*sum(F('.2')**i*F(cold['N'][j-i]) for i in range(j+1)) for j in range(9)]
    one=[F(1)]*9
    mm=[F('50.1')**2*x for x in R.tensor(rm,one,6,200)]
    cn=[F('.001')*F('1.001')**2*x for x in R.tensor([82*x for x in rc],A,100,20)]
    ca=[F('.00331')*F('1.001')**2*x for x in R.tensor(rc,one,15,200)]
    outer=[F(x) for x in pf['regions']['V']['outer_weighted']]
    hh=poly([min(KV[j],max(mm[j],cn[j],ca[j],outer[j])*100**j) for j in range(9)])
    check('tau_high_order_gain',all(hh[j]<KV[j] for j in [6,7,8]))
    check('low_Hessian_replacement',V[2]==F('.002') and KV[2]==F('.005'))
    check('arrays_only_strengthened',all(a<=F(b) for arr,name in [(U,'U'),(K,'K'),(V,'Vbar'),(KV,'KV')] for a,b in zip(arr,pf[name])))
    FIXTURES.append({'fixture':'low_and_tau_retaining_gap_arrays','U':strings(U[:9]),'K':strings(K[:9]),
                     'Vbar':strings(V[:9]),'KV':strings(KV[:9]),'tau_KV':strings(hh[:9]),
                     'tau_KV_decimals':decimals(hh[:9]),'regional_tau_V':{'middle':strings(mm),'cold_near':strings(cn),'cold_away':strings(ca),'outer':strings(outer)},
                     'cold_near_ratio_monotonicity_assumed':False})
    return U,K,V,KV,hh

def gap_product(v,h,q,order=6):
    out=poly()
    for n in range(order+1):
        for k in range(n+1):
            for i in range(n-k+1):
                j=n-k-i;a=(i+1)*v[i+1];b=(j+1)*v[j+1]
                if k<2:
                    pair=min(a*(j+1)*h[j+1],(i+1)*h[i+1]*b)
                    out[n]+=pair*(1 if k==0 else F('1.01'))
                else:out[n]+=q[k]*a*b
    return out

def scalar_ledger(U,K,V,KV,h):
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    original=fixture('finite-scalar-residual-six','actual_scalar_residual_H6_ledger');previous=fixture(PREVIOUS,'refined_scalar_H6_ledger')
    q=add(poly([1]),poly([F(x) for x in pf['seed_jets']]));ll=poly([F(x) for x in pf['log_radius_jets']]);C0=poly([F(x) for x in pf['connection']])
    Tp=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    Vp=V.copy();Vp[0]=0;Up=U.copy();Up[0]=0;E=exps(add(Vp,Tp));O=scale('1.001',exps(scale(4,Tp)))
    S1=scale('1.001',mul(E,K));S1u=scale('1.001',mul(E,U));S2=scale('1.001',mul(mul(q,exps(add(Vp,scale(2,Tp)))),mul(K,KV)))
    W1=add(K,S1);J2=mul(O,add(mul(add(scale(4,power(W1,2)),scale(2,S2)),add(poly([1]),scale(4,U))),scale(16,mul(K,W1))))
    B=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')]);zabs=compose(za,q);Zd=poly([F(x) for x in original['linear_trace_W2']])
    LN=scale('6.1',mul(power(q,2),add(scale(2,der(der(K))),scale(2,mul(der(ll),der(K))))));LN[0]=min(LN[0],F('.005'))
    LT=scale('6.1',add(scale(2,mul(q,compose(B,q))),mul(power(q,2),compose(der(B),q))))
    H=mul(KV,LN)
    for n in range(7):H[n]-=KV[n]*LN[0]-min(KV[n]*LN[0],F('.005')*h[n])
    G=gap_product(KV,h,q);plain=mul(q,power(der(KV),2))
    G=poly([min(a,b) for a,b in zip(G,plain)])
    check('tau_gradient_product_improves_sixth',G[6]<plain[6],{'old_low_improved':L.decimal(plain[6]),'tau_retained':L.decimal(G[6])})
    check('tau_Lnu_product_improves_sixth',H[6]<mul(KV,LN)[6])
    QQ=scale('.501',mul(E,power(KV,2)))
    terms=[scale(F(1,2),mul(add(poly([12]),zabs),J2)),scale(2,mul(mul(O,add(U,S1u)),Zd)),
           scale(6*F('1.001'),mul(mul(mul(O,q),E),H)),scale(6,mul(mul(O,QQ),LT)),
           scale(6*F('6.1')*F('1.001'),mul(mul(mul(O,E),power(q,2)),G)),
           scale(F('4.004')*F('1.3e-10'),mul(mul(power(q,2),exps(add(scale(4,Tp),scale(2,Up)))),power(K,2)))]
    new=add(*terms);selected=poly([min(new[j],F(previous['selected_curvature_remainder'][j])) for j in range(7)])
    origP=poly([F(x) for x in original['residual_reference_W2']]);origN=scale(F(1,2),poly([F(x) for x in original['curvature_remainder_second_W2']]))
    P=poly([origP[j]-origN[j]+selected[j] for j in range(7)])
    Hphys=L.covariant(P,C0,0,6);constant=F(1001,24000)*sum(Hphys);old=F(previous['C_tau_exact']);exponent=L.decimal_power(constant)
    check('six_remainder_terms_positive',len(terms)==6 and all(min(t[:7])>0 for t in terms))
    check('previous_covariant_constant_reproduced',F(1001,24000)*sum(L.covariant(poly([F(x) for x in previous['residual_reference_W2']]),C0,0,6))==old)
    check('complete_scalar_bound_improves',constant<old,{'old':L.decimal(old),'new':L.decimal(constant),'improvement':L.decimal(old/constant)})
    check('new_decimal_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent)
    check('previous_physical_volume_retained',saved(PREVIOUS)['physical_noncompact_volume_retained'])
    budget=F(10**57)/(800*F(10)**exponent)
    check('rounded_scalar_only_budget',800*F(10)**exponent*budget==10**57)
    FIXTURES.append({'fixture':'cold_product_scalar_H6_ledger','gradient_square_product':strings(G[:7]),'V_Lnu_product':strings(H[:7]),
                     'remainder_terms':[strings(t[:7]) for t in terms],'remainder_terms_decimals':[decimals(t[:7]) for t in terms],
                     'selected_curvature_remainder':strings(selected[:7]),'residual_reference_W2':strings(P[:7]),
                     'physical_covariant_bounds':strings(Hphys),'physical_covariant_decimals':decimals(Hphys),
                     'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),'C_tau_power_of_ten':exponent,
                     'old_C_tau_exact':str(old),'improvement_ratio':L.decimal(old/constant),
                     'scalar_only_rounded_transfer_budget_over_g_epsilon_c':str(budget),'scalar_only_budget_decimal':L.decimal(budget),
                     'old_anomaly_reference_and_Lie_arrays_retained':True,'full_Schur_source_bounded':False,'target_passed':False})

def independent_regular_gap():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            a=mp.mpf('.0003');b=mp.mpf('.0001');M=a+b
            def direct(r):
                if r==0:return mp.mpf(0)
                cm1=2*mp.sinh(r/2)**2;co1=-2*mp.sin(r/2)**2
                return (a*cm1+b*(cm1*mp.cos(r)+co1+mp.sinh(r)*mp.sin(r))/2)/mp.sinh(r)
            for val in ['0','.01','.5','1','5']:
                r=mp.mpf(val)
                integ=mp.mpf(0) if r==0 else mp.quad(lambda s:mp.sinh(s)*(a+b*mp.cos(s)),[0,r])/mp.sinh(r)
                angular=M/2 if r==0 else mp.coth(r)*direct(r)
                second=a+b*mp.cos(r)-angular
                error=max(abs(integ-direct(r)),abs(mp.diff(direct,r)-second))
                check('regular_gap_integral_'+val+'_'+str(dps),error<mp.mpf('1e-42') and abs(angular)<mp.mpf('2.13')*M and abs(second)+abs(angular)<mp.mpf('5.26')*M,mp.nstr(error,18))
                records.append({'rho':val,'dps':dps,'error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_regular_bolt_gap','records':records,'includes_exact_bolt':True,'actual_charged_source_sampled':False})

def independent_cartesian_product():
    x,y=S.symbols('x y');rr=x*x+y*y
    vv=(1+2*x*x+3*y*y+rr**2+x**6+2*y**6+rr**4)/1000
    gradient=S.diff(vv,x)**2+S.diff(vv,y)**2;records=[]
    def jets(f,point,order):
        return [sum(S.binomial(n,k)*abs(S.diff(f,x,n-k,y,k).subs(point)) for k in range(n+1))/S.factorial(n) for n in range(order+1)]
    for tau0 in [S.Rational(1,10**6),S.Rational(1,5)]:
        tau=tau0*(1+rr/10000)
        for point in [{x:0,y:0},{x:S.Rational(1,5),y:S.Rational(3,10)}]:
            tc=tau.subs(point);v=[F(str(z)) for z in jets(vv,point,7)];h=[F(str(tc))*z for z in v]
            tj=[F(str(z)) for z in jets(tau,point,6)];direct=[F(str(z)) for z in jets(S.expand(tau*gradient),point,6)]
            check('Cartesian_tau_relative_'+str(tau0)+'_'+str(point),tj[1]<F('1.01')*F(str(tc)))
            bound=gap_product(poly(v),poly(h),poly(tj))
            for n in range(7):check('Cartesian_full_product_'+str(tau0)+'_'+str(point)+'_'+str(n),direct[n]<=bound[n])
            omitted=S.diff(tau*gradient,x,6)-tau*S.diff(gradient,x,6)
            check('omitted_tau_derivatives_detected_'+str(tau0)+'_'+str(point),omitted.subs(point)!=0)
            records.append({'tau0':str(tau0),'point':[str(point[x]),str(point[y])],'sixth_L1_tensor':str(direct[6]),'sixth_product_bound':str(bound[6])})
    FIXTURES.append({'fixture':'independent_full_Cartesian_gap_products','records':records,'tensor_order':6,'tau_derivatives_retained':True,'actual_charged_geometry':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold gap scalar products','input_count':len(ps)}),flush=True)
    try:
        provenance();gap_equation_and_cover();args=low_and_tau_arrays();scalar_ledger(*args);independent_regular_gap();independent_cartesian_product()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 cold gap scalar products','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
