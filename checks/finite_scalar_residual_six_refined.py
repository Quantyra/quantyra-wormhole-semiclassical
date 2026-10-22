"""Refine scalar H6 using retained per-order jets and exact curvature remainder."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial, comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-scalar-residual-six-refined';PREVIOUS='finite-scalar-residual-six'
OWN=['checks/finite_scalar_residual_six_refined.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_covariant_calculus',ROOT/'checks/finite_higher_local_coefficients.py')
L=importlib.util.module_from_spec(spec);spec.loader.exec_module(L)
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[]
PROTOCOL={'scalar_source_order':6,'new_tangent_order':8,'old_anomaly_metric_order':10,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42','actual_charged_source_sampled':False,
          'ordinary_continuum_proof_separate':True,'all_six_remainder_terms_retained':True}

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
    check('previous_scalar_run_passed',old['status']=='passed' and len(old['checks'])==48 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    for stem in ['middle-metric-ten-jets','cold-metric-ten-jets','outer-metric-ten-jets']:
        r=saved(stem);check(stem+'_passed',r['status']=='passed' and all(x['passed'] for x in r['checks']))
    check('scope_retained',review()['actual_scalar_residual_H6_bounded'] and not review()['actual_transfer_certificate_passed'])
    FIXTURES.append({'fixture':'retained_sources_and_scope','previous_checks':48,'previous_inputs':573,'old_scientific_sources_edited':False})

def tensor(radial,initial,c,s):
    n=len(radial)-1;B={(0,k):initial[k]*radial[k]*factorial(k)*100**k for k in range(n+1)}
    for r in range(n):
        for k in range(n-r):B[r+1,k]=B[r,k+1]+r*sum(F(comb(k,i))*c*factorial(i)*s**i*B[r,k-i] for i in range(k+1))
    return [B[j,0]/(factorial(j)*100**j) for j in range(n+1)]

def regional_arrays():
    mf=fixture('middle-metric-ten-jets','radial_two_derivative_reconstruction')
    cf=fixture('cold-metric-ten-jets','weighted_radial_reconstruction')
    ca=fixture('cold-metric-ten-jets','original_finite_completion_and_tensor_conversion')
    of=fixture('outer-metric-ten-jets','outer_two_derivative_reconstruction')
    oa=fixture('outer-metric-ten-jets','finite_completion_full_tensors_and_forward_map')
    A=[F(x) for x in ca['Cartesian_amplifications'][:9]];one=[F(1)]*9
    out={};records={}
    for field in ['u','nu','V']:
        key='U' if field=='u' else 'N';mid=[F(x) for x in mf[key][:9]];cold=[F(x) for x in cf[key][:9]]
        if field=='V':
            mid=[F(1200)*sum(F('.2')**i*F(mf['N'][j-i]) for i in range(j+1)) for j in range(9)]
            cold=[F('1.2')*82*sum(F('.2')**i*F(cf['N'][j-i]) for i in range(j+1)) for j in range(9)]
        mm=tensor(mid,one,6,200);cn=tensor(cold,A,100,20);cw=tensor(cold,one,15,200)
        oo=[F(x) for x in of['common_'+field][:9]]
        outer=[9*F('.0331')*F(oa['tensor_amplifications'][j])*max(oo[:j+1])+F('18e-11')*F('.1') for j in range(9)]
        cold_factor=F(1) if field=='V' else F('.00331')
        cc=[max(a,b)*cold_factor for a,b in zip(cn,cw)]
        un=[max(mm[j],cc[j],outer[j]) for j in range(9)]
        we=[max(F('50.1')**2*mm[j],F('1.001')**2*cc[j],outer[j]) for j in range(9)]
        out[field]=(un,we);records[field]={'middle':strings(mm),'cold':strings(cc),'outer_weighted':strings(outer)}
    oldU=[F('3e-8')]+[F('.01')]*5+[F('4e21')]*3
    oldK=[F('1e-4')]+[F(10)]*5+[F('4e21')]*3
    oldV=[F('1.2e-8')]+[F(100)]*5+[F('2e26')]*3
    oldKV=[F('4e-5')]+[F(4000)]*5+[F('2e26')]*3
    U=[min(oldU[j],max(out['u'][0][j],out['nu'][0][j]))*100**j for j in range(9)]
    K=[min(oldK[j],max(out['u'][1][j],out['nu'][1][j]))*100**j for j in range(9)]
    V=[min(oldV[j],out['V'][0][j])*100**j for j in range(9)]
    KV=[min(oldKV[j],out['V'][1][j])*100**j for j in range(9)]
    U[1]=min(U[1],F('4e-7'));V[1]=min(V[1],F('4e-7'))
    K[1]=min(K[1],F('.001'));KV[1]=min(KV[1],F('.001'))
    check('unweighted_gradient_conversion',F('8e-7')**2<F('5.9')*F('4e-7')**2)
    check('core_weighted_gradient',F('50.1')**4*F('8e-7')**2<F('5.9')*F('.001')**2)
    check('outer_weighted_gradient',9*F('.0331')*F('1e-7')*100+F('18e-11')*F('.1')*100<F('3e-6')<F('.001'))
    check('outer_V_gradient',F('4.3e-8')**2<F('5.9')*F('2e-8')**2)
    check('weighted_L4nu_value_core',F('50.1')**2*F('2.2e-6')<F('.006'))
    check('weighted_L4nu_value_outer',9*F('.0331')*F('1e-7')*2*100**2+F('18e-11')*F('.1')*2*100**2<F('.0006') and F('6.1')*(2*F('.0006')+2*F('3e-6'))<F('.008')<F('.1'))
    check('per_order_tangent_improved',all(K[j]<oldK[j]*100**j for j in range(6,9)))
    check('per_order_V_improved',all(KV[j]<oldKV[j]*100**j for j in range(6,9)))
    # Independent reconstruction with constant arrays reproduces frozen tensor factors.
    check('middle_tensor_reproduced',tensor(one,one,6,200)==[F(x) for x in fixture('middle-metric-ten-jets','finite_completion_and_physical_tensor_ledger')['tensor_amplifications'][:9]])
    check('cold_near_tensor_reproduced',tensor(one,A,100,20)==[F(x) for x in ca['near_tensor_amplifications'][:9]])
    check('cold_away_tensor_reproduced',tensor(one,one,15,200)==[F(x) for x in ca['away_tensor_amplifications'][:9]])
    # The underlying tensor map is linear in its full nonconstant input array.
    rr=[F(j+1,j+2) for j in range(9)];ss=[F(2*j+1,3*j+2) for j in range(9)]
    check('nonconstant_tensor_linearity',tensor([a+b for a,b in zip(rr,ss)],A,100,20)==[a+b for a,b in zip(tensor(rr,A,100,20),tensor(ss,A,100,20))])
    FIXTURES.append({'fixture':'separate_regional_tangent_orders','regions':records,'U':strings(U),'K':strings(K),'Vbar':strings(V),'KV':strings(KV),'K_decimals':decimals(K),'KV_decimals':decimals(KV),'cold_inverse_gap_used':False})
    return poly(U),poly(K),poly(V),poly(KV)

def exact_identities():
    e=S.symbols('ell',real=True);n=S.symbols('n');w=S.Function('w')(e)
    J=S.exp(-2*w)*(1+2*e*n)-1
    rhs=S.exp(-2*w)*(4*S.diff(w,e)**2*(1+2*e*n)-8*n*S.diff(w,e)-2*S.diff(w,e,2)*(1+2*e*n))
    zero('second_exponential_remainder',S.diff(J,e,2)-rhs)
    zero('first_zero',S.diff(J,e).subs({S.diff(w,e):n,w:0,e:0}))
    z=S.symbols('z');b=S.Function('b')(z);f=S.Function('f')(z);cx=S.symbols('cx')
    lap=cx*(z*z*b*S.diff(f,z,2)+(z*z*S.diff(b,z)-2*z*b)*S.diff(f,z))
    zero('seed_L4_tau',lap.subs(f,1-z).doit()-cx*(2*z*b-z*z*S.diff(b,z)))
    v,tau=S.symbols('v tau');T=1-tau+tau*S.exp(e*v);k=S.exp(e*v)
    zero('small_radius_factor',k/T-1-(1-tau)*(k-1)/T)
    zero('B_integral_derivative',S.diff(k-1-e*v*k,e)+e*v*v*k)
    x=S.symbols('x');tt=S.Function('tau')(x);vv=S.Function('V')(x);aa=S.Function('a')(x)
    LL=lambda f:S.diff(f,x,2)+aa*S.diff(f,x)
    nu=tt*vv;kk=S.exp(e*vv);TT=1+tt*(kk-1);BB=kk-1-e*vv*kk
    zero('exact_radius_operator',LL(TT)-e*kk*LL(nu)-BB*LL(tt)-e*e*kk*tt*S.diff(vv,x)**2)
    zero('log_gradient_cancellation',LL(S.log(TT))+S.diff(S.log(TT),x)**2-LL(TT)/TT)
    check('gradient_square_omission_detected',S.simplify(LL(TT)-e*kk*LL(nu)-BB*LL(tt))!=0)
    FIXTURES.append({'fixture':'exact_scalar_remainder_identities','parameter_zero_order':2,'gradient_square_retained':True,'radius_acceleration_retained':True})

def refined_ledger(U,K,V,KV):
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    old=fixture(PREVIOUS,'actual_scalar_residual_H6_ledger')
    Z=poly([F(x) for x in pf['seed_jets']]);q=add(poly([1]),Z);ll=poly([F(x) for x in pf['log_radius_jets']]);C0=poly([F(x) for x in pf['connection']])
    Tp=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    Vp=V.copy();Vp[0]=0;Up=U.copy();Up[0]=0
    E=exps(add(Vp,Tp));O=scale('1.001',exps(scale(4,Tp)))
    S1=scale('1.001',mul(E,K));S1u=scale('1.001',mul(E,U));S2=scale('1.001',mul(mul(q,exps(add(Vp,scale(2,Tp)))),mul(K,KV)))
    W1=add(K,S1)
    J2=mul(O,add(mul(add(scale(4,power(W1,2)),scale(2,S2)),add(poly([1]),scale(4,U))),scale(16,mul(K,W1))))
    B=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')]);zabs=compose(za,q)
    Zd=poly([F(x) for x in old['linear_trace_W2']])
    LN=scale('6.1',mul(power(q,2),add(scale(2,der(der(K))),scale(2,mul(der(ll),der(K))))));LN[0]=min(LN[0],F('.1'))
    LT=scale('6.1',add(scale(2,mul(q,compose(B,q))),mul(power(q,2),compose(der(B),q))))
    AA=scale('1.001',mul(mul(q,E),KV));QQ=scale('.501',mul(E,power(KV,2)))
    terms=[scale(F(1,2),mul(add(poly([12]),zabs),J2)),scale(2,mul(mul(O,add(U,S1u)),Zd)),
           scale(6,mul(mul(O,AA),LN)),scale(6,mul(mul(O,QQ),LT)),
           scale(6*F('6.1')*F('1.001'),mul(mul(mul(O,E),power(q,3)),power(der(KV),2))),
           scale(F('4.004')*F('1.3e-10'),mul(mul(power(q,2),exps(add(scale(4,Tp),scale(2,Up)))),power(K,2)))]
    new=add(*terms);oldN=scale(F(1,2),poly([F(x) for x in old['curvature_remainder_second_W2']]))
    NN=poly([min(a,b) for a,b in zip(new,oldN)])
    # Keep every non-curvature term from the original scalar residual verbatim.
    oldP=poly([F(x) for x in old['residual_reference_W2']])
    P=poly([oldP[j]-oldN[j]+NN[j] for j in range(7)])
    check('old_other_terms_positive',all(oldP[j]>oldN[j] for j in range(7)))
    check('six_exact_remainder_terms',len(terms)==6 and all(min(t[:7])>0 for t in terms))
    check('new_curvature_component_improves',new[6]<oldN[6],{'old':L.decimal(oldN[6]),'new':L.decimal(new[6])})
    check('pointwise_min_preserves_old_bounds',all(0<P[j]<=oldP[j] for j in range(7)))
    H=L.covariant(P,C0,0,6);constant=F(1001,24000)*sum(H);oldconst=F(old['C_tau_exact'])
    check('physical_conversion_old_reproduced',L.covariant(oldP,C0,0,6)==[F(x) for x in old['physical_covariant_bounds']])
    check('scalar_constant_strictly_improves',constant<oldconst,{'new':L.decimal(constant),'old':L.decimal(oldconst),'ratio':L.decimal(oldconst/constant)})
    exponent=L.decimal_power(constant)
    check('upward_decimal_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent)
    check('same_physical_volume',saved(PREVIOUS)['physical_noncompact_volume_retained'] and F('1.001')**2>F('1.001'))
    FIXTURES.append({'fixture':'refined_scalar_H6_ledger','remainder_terms':[strings(t[:7]) for t in terms],
                     'remainder_terms_decimals':[decimals(t[:7]) for t in terms],'new_curvature_remainder':strings(new[:7]),
                     'selected_curvature_remainder':strings(NN[:7]),'residual_reference_W2':strings(P[:7]),
                     'physical_covariant_bounds':strings(H),'physical_covariant_decimals':decimals(H),
                     'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),'C_tau_power_of_ten':exponent,
                     'old_C_tau_exact':str(oldconst),'improvement_ratio':L.decimal(oldconst/constant),
                     'old_anomaly_reference_and_Lie_arrays_retained':True,'full_Schur_source_bounded':False})
    budget=F(10**57)/(800*F(10)**exponent)
    check('rounded_scalar_budget',800*F(10)**exponent*budget==10**57)
    FIXTURES.append({'fixture':'quantitative_transfer_reassessment','scalar_only_rounded_transfer_budget_over_g_epsilon_c':str(budget),
                     'scalar_only_budget_decimal':L.decimal(budget),'actual_transfer_factors_quantified':False,
                     'upper_allowance_is_not_a_lower_bound':True,'target_passed':False})

def independent_finite_remainders():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for t0 in [mp.mpf('.2'),mp.mpf('.8')]:
                for ee in [mp.mpf('.3'),mp.mpf('.9')]:
                    def quantities(x):
                        ta=t0+x/100;z=1-ta;phi=-mp.log(z);pp=1/(100*z);ppp=1/(10000*z*z)
                        u=(1+2*x+x*x+x**4)/1000;up=(2+2*x+4*x**3)/1000;upp=(2+12*x*x)/1000
                        v=(2-x+3*x*x+x**3)/1000;vp=(-1+6*x+3*x*x)/1000;vpp=(6+6*x)/1000
                        nu=ta*v;nup=v/100+ta*vp;nupp=2*vp/100+ta*vpp
                        k=mp.exp(ee*v);T=1+ta*(k-1);Tx=(k-1)/100+ta*k*ee*vp
                        Txx=2*k*ee*vp/100+ta*k*(ee*vpp+ee**2*vp**2)
                        ss=mp.log(T);sp=Tx/T;spp=Txx/T-sp*sp;w=ee*u+ss
                        base=z*z*(-6*ppp-6*pp*pp+2)
                        LN=z*z*(nupp+2*pp*nup);LT=z*z*2*pp/100
                        dot=-2*(u+nu)*base-2*z*z*upp-6*LN+4*u*z*z
                        direct=mp.exp(-2*(phi+w))*(-2*(ppp+ee*upp+spp)-4*(ppp+spp)-6*(pp+sp)**2)+2*mp.exp(-2*(phi+ss))
                        B=k-1-ee*v*k
                        remainder=base*(mp.exp(-2*w)*(1+2*ee*(u+nu))-1)+ee*mp.expm1(-2*w)*dot
                        remainder-=6*mp.exp(-2*w)*(ee*(k/T-1)*LN+B/T*LT+ee**2*k/T*ta*z*z*vp*vp)
                        remainder+=2*mp.exp(-2*w)*(mp.expm1(2*ee*u)-2*ee*u)*z*z
                        return direct-base-ee*dot,remainder
                    radius=mp.mpf('.05');count=64
                    contour=sum(quantities(radius*mp.exp(2j*mp.pi*j/count))[0]*mp.exp(-12j*mp.pi*j/count) for j in range(count))/(count*radius**6)
                    real=mp.diff(lambda x:quantities(x)[1],mp.mpf(0),6)/mp.factorial(6)
                    error=abs(contour-real)
                    check('independent_finite_sixth_'+str(t0)+'_'+str(ee)+'_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    records.append({'tau0':str(t0),'ell':str(ee),'dps':dps,'sixth_coefficient_error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_finite_warped_remainders','records':records,'actual_charged_source_sampled':False,'path_parameter_nonzero':True})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 refined scalar residual H6','input_count':len(ps)}),flush=True)
    try:
        provenance();inputs=regional_arrays();exact_identities();refined_ledger(*inputs);independent_finite_remainders()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 refined scalar residual H6','created_utc':datetime.now(timezone.utc).isoformat(),
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
