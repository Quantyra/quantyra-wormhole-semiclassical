"""Low radial reconstruction, regular divided gap, and complete scalar bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-low-gap-reconstruction';PREVIOUS='finite-anomaly-tangent-geometry'
OWN=['checks/finite_low_gap_reconstruction.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_complete_scalar',ROOT/'checks/finite_anomaly_tangent_geometry.py')
D=importlib.util.module_from_spec(spec);spec.loader.exec_module(D);L=D.L;R=D.B.R
poly=L.poly;ROWS=[];FIXTURES=[];g=F('6.8e-28');A=F('.00331')
PROTOCOL={'source_order':8,'tangent_order':10,'scalar_remainder_order':8,'scalar_source_order':6,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42','actual_charged_source_sampled':False,
          'ordinary_continuum_proof_separate':True}
def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def strings(a):return [str(x) for x in a]
def decimals(a):return [L.decimal(x) for x in a]
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def cc(c,a,n):return sum(F(c)*2**i*a[n-i] for i in range(n+1))
def dd(c,a,n):return sum(F(c)*2**i*100*(n-i+1)*a[n-i+1] for i in range(n+1))
def rawc(c,a,n):return sum(F(comb(n,i))*F(c)*factorial(i)*200**i*factorial(n-i)*100**(n-i)*a[n-i] for i in range(n+1))
def rawd(c,a,n):return sum(F(comb(n,i))*F(c)*factorial(i)*200**i*factorial(n-i+1)*100**(n-i+1)*a[n-i+1] for i in range(n+1))
def invsum(a,j):return sum(F('.2')**i*a[j-i] for i in range(j+1))
def tau_product(a,j):return F('.001')*a[j]+sum(F('.004')*2**i*a[j-i] for i in range(1,j+1))

def inputs_and_coefficients():
    old=saved(PREVIOUS)
    check('previous_complete_scalar_passed',old['status']=='passed' and len(old['checks'])==71 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('middle_canonical_gradient',F('4e-7')*F('6e-6')<F('1.6e-6')**2)
    check('cold_canonical_gradient',F('4e-7')*F('1.001')<F('4.1e-7'))
    bz=2+3*F('1.334')*F('.004')+4*F('.501')*F('.004')**2+F('1e-20')+3*F('1.002')*F('1.2')*F('.001')*F('1.004')**2
    sb=F('2.1')/F('.996')+2*9*F('.004')/F('.996')**2
    check('complex_bz_over_tau',bz<F('2.1'))
    check('complex_S_over_tau',sb<4 and 2*sb<8 and 10+F('1.004')*sb<16)
    check('divided_mass_coefficients',1+3<=4 and 3*F('1.004')<4)
    check('coupled_u_V_coefficient',10*F('.004')<=F('.04'))
    check('local_Eu_before_weight_enlargement',1+3*F('1.004')**2<F('4.1'))
    check('local_Enu_division_before_enlargement',F('1.004')**2*F('1.2')*F('1.001')<F('1.3'))
    check('local_Erho_uniform',F('1e22')*F('.001')==F('1e19'))
    check('cold_regular_drift',12+2+6==20)
    t=S.symbols('t',positive=True)
    for j,bound in [(0,2600),(1,1003),(2,1001)]:
        f=t**(3-j)/(t-1);dp=S.factor(S.diff(f,t))
        check('middle_gap_weight_endpoints_'+str(j),max(F('1.001')**(3-j)/F('.001'),F('50.1')**(3-j)/F('49.1'))<=bound,str(dp))
    check('middle_high_gap_weight_monotone',S.factor(S.diff(1/(t-1),t))==-1/(t-1)**2)
    check('fixed_observation_weight_power',S.simplify(t**2*(1/t)**8-t**-6)==0)
    src=fixture('finite-source-order-scalar','complete_per_order_sources')
    rad=fixture('finite-source-order-scalar','per_order_original_reconstruction')
    FIXTURES.append({'fixture':'retained_sources_and_uniform_coefficients','source_inputs_unchanged':True,
                     'Eu_uniform':'4.1e19','Enu_over_tau_uniform':'1.3e22','Erho_uniform':'1e19',
                     'middle_initial_value':'1.2e-8','middle_initial_derivative_divided':'1.6e-8','inverse_thermal_gap_used':False})
    return [F(x) for x in src['middle']],[F(x) for x in src['cold']],rad

def middle(t,old):
    f=[F(x) for x in old['F']];u=[F('1.2e-8'),F('1.6e-8')];v=u.copy();ur=[];vr=[]
    for n in range(9):
        den=10000*(n+1)*(n+2)
        ru=(dd(2,u,n)+cc(12,u,n)+cc(10,v,n)+3*f[n]+3*g*cc(2,t,n)+g*10**18*2**n)/den
        rv=(dd(100,v,n)+cc(20000,u,n)+cc(20000,v,n)+f[n]+g*cc(2,t,n)+g*10**18*2**n)/den
        ur.append(ru);vr.append(rv);u.append(min(ru,F(old['U'][n+2])));v.append(min(rv,F(old['N'][n+2])))
    for n in [0,3,8]:
        fac=factorial(n)*100**n
        rawu=rawd(2,u,n)+rawc(12,u,n)+rawc(10,v,n)+3*fac*f[n]+3*g*rawc(2,t,n)+g*10**18*2**n*fac
        rawv=rawd(100,v,n)+rawc(20000,u,n)+rawc(20000,v,n)+fac*f[n]+g*rawc(2,t,n)+g*10**18*2**n*fac
        check('independent_middle_raw_'+str(n),rawu==ur[n]*factorial(n+2)*100**(n+2) and rawv==vr[n]*factorial(n+2)*100**(n+2))
    check('middle_low_radial_improves',all(u[j]<F(old['U'][j]) and v[j]<F(old['N'][j]) for j in range(6)))
    return {'U':u,'N':v}

def cold(t,old):
    q=[max(F(x),F(y)) for x,y in zip(old['Ward_near'],old['Ward_away'])]
    tu=[A*x for x in t];qu=[A*x for x in q];tv=[F('1.2')*82*invsum(t,j) for j in range(9)]
    mass_source=[qu[j]+tu[j]+10**19*2**j for j in range(9)];f=[F('82e-10')]
    for n in range(8):f.append((cc(4,f,n)+2*g*cc(4,mass_source,n))/(100*(n+1)))
    fullf=[min(A*F(old['F'][j]),tau_product(f,j)) for j in range(9)]
    ucap=[A*F(x) for x in old['U']];vcap=[F('1.2')*82*invsum([F(x) for x in old['N']],j) for j in range(11)]
    vcap[2]=min(vcap[2],F('3.13')*F('.0006')/(2*100**2))
    u=[F('1.2e-8'),F('4.1e-9')];v=[F('1.2e-8'),F('3.4e-9')];ur=[];vr=[]
    for n in range(9):
        den=10000*(n+1)*(n+2)
        ru=2*(dd(12,u,n)+cc(12,u,n)+cc('.04',v,n)+3*fullf[n]+3*g*cc(2,tu,n)+g*F('4.1e19')*2**n)/den
        rv=2*(dd(20,v,n)+cc(8,u,n)+cc(16,v,n)+f[n]+g*cc(2,tv,n)+g*F('1.3e22')*2**n)/den
        ur.append(ru);vr.append(rv);u.append(min(ru,ucap[n+2]));v.append(min(rv,vcap[n+2]))
    nu=[min(tau_product(v,j),A*F(old['N'][j])) for j in range(11)]
    for n in [0,3,8]:
        fac=factorial(n)*100**n
        rawu=2*(rawd(12,u,n)+rawc(12,u,n)+rawc('.04',v,n)+3*fullf[n]*fac+3*g*rawc(2,tu,n)+g*F('4.1e19')*2**n*fac)
        rawv=2*(rawd(20,v,n)+rawc(8,u,n)+rawc(16,v,n)+f[n]*fac+g*rawc(2,tv,n)+g*F('1.3e22')*2**n*fac)
        check('independent_cold_raw_'+str(n),rawu==ur[n]*factorial(n+2)*100**(n+2) and rawv==vr[n]*factorial(n+2)*100**(n+2))
    for n in [0,3,7]:
        check('independent_divided_mass_raw_'+str(n),rawc(4,f,n)+2*g*rawc(4,mass_source,n)==f[n+1]*factorial(n+1)*100**(n+1))
    check('cold_radial_caps_retained',all(x<=y for x,y in zip(u,ucap)) and all(x<=y for x,y in zip(v,vcap)))
    check('cold_gap_low_orders_improve',all(v[j]<vcap[j] for j in [3,4,5]))
    check('divided_mass_improves',all(f[j]<F('1.2')*82*invsum([F(x) for x in old['F']],j) for j in range(9)))
    return {'U':u,'N':nu,'V':v,'F_over_tau':f,'F':fullf,'T_over_tau':tv}

def tensors(mid,cold):
    prior=fixture('finite-source-order-scalar','full_tangent_and_tau_arrays')
    initial=[F(x) for x in fixture('cold-metric-ten-jets','original_finite_completion_and_tensor_conversion')['Cartesian_amplifications']]
    one=[F(1)]*11;weight=[F('50.1')**2,F('50.1')]+[F(1)]*9
    middle_t={field:R.tensor(mid[k],one,6,200) for field,k in [('u','U'),('nu','N')]}
    cold_t={};regions={}
    for field,k in [('u','U'),('nu','N'),('V','V')]:
        near=R.tensor(cold[k],initial,100,20);away=R.tensor(cold[k],one,15,200)
        cold_t[field]=[max(a,b) for a,b in zip(near,away)]
        regions[field]={'cold':strings(cold_t[field]),'outer_weighted':prior['regions'][field]['outer_weighted']}
    gap=R.tensor([invsum(mid['N'],j) for j in range(11)],one,6,200)
    arrays={}
    for field in ['u','nu','V']:
        outer=[F(x) for x in prior['regions'][field]['outer_weighted']]
        if field=='V':
            mu=[1200*x for x in gap];mw=[F('1.2')*([2600,1003][j] if j<2 else 1001)*gap[j] for j in range(11)]
        else:mu=middle_t[field];mw=[weight[j]*mu[j] for j in range(11)]
        un=[max(mu[j],cold_t[field][j],outer[j])*100**j for j in range(11)]
        we=[max(mw[j],F('1.001')**2*cold_t[field][j],outer[j])*100**j for j in range(11)]
        arrays[field]=(un,we);regions[field].update(middle_unweighted=strings(mu),middle_weighted=strings(mw))
    result=[]
    for name,field,w in [('U','u',0),('K','u',1),('Vbar','V',0),('KV','V',1)]:
        raw=arrays[field][w]
        if field=='u':raw=[max(a,b) for a,b in zip(raw,arrays['nu'][w])]
        result.append(poly([min(x,F(y)) for x,y in zip(raw,prior[name])]))
    h=poly([min(result[3][j],F(prior['tau_KV'][j]),
                max(F('1.2')*weight[j]*gap[j],F('.001')*F('1.001')**2*cold_t['V'][j],F(prior['regions']['V']['outer_weighted'][j]))*100**j) for j in range(11)])
    result.append(h)
    for name,arr in zip(['U','K','Vbar','KV','tau_KV'],result):
        check(name+'_all_old_bounds_retained',all(arr[j]<=F(prior[name][j]) for j in range(11)))
        check(name+'_middle_low_orders_improve',all(arr[j]<F(prior[name][j]) for j in [3,4,5]))
    FIXTURES.append({'fixture':'low_and_divided_radial_reconstruction','middle':{k:strings(v) for k,v in mid.items()},
                     'cold':{k:strings(v) for k,v in cold.items()},'middle_U_decimals':decimals(mid['U']),
                     'middle_N_decimals':decimals(mid['N']),'cold_V_decimals':decimals(cold['V'])})
    FIXTURES.append({'fixture':'complete_improved_tangent_arrays','regions':regions,
                     **{k:strings(v[:11]) for k,v in zip(['U','K','Vbar','KV','tau_KV'],result)},
                     'K_decimals':decimals(result[1][:11]),'KV_decimals':decimals(result[3][:11]),
                     'outer_mass_remainder_and_cubic_retained':True,'fixed_observation_z_powers_retained':True})
    return result

def complete_scalar(args):
    # Keep original old-contact recovery, but use the latest scalar caps/comparator.
    latest=fixture(PREVIOUS,'complete_refined_scalar_and_anomaly')
    original_fixture=D.fixture
    comparison={**fixture('finite-source-order-scalar','per_order_scalar_and_anomaly_ledger'),**latest}
    def retained_fixture(stem,name):
        if stem==D.PREVIOUS and name=='per_order_scalar_and_anomaly_ledger':return comparison
        return original_fixture(stem,name)
    D.fixture=retained_fixture;D.check=check;D.FIXTURES=FIXTURES
    U,K,V,KV,h=args;t=D.finite_bounds(U,V);D.scalar_ledger(t,U,K,V,KV,h)
    f=next(x for x in FIXTURES if x['fixture']=='complete_refined_scalar_and_anomaly')
    check('comparison_is_latest_complete_scalar',f['old_C_tau_exact']==latest['C_tau_exact'])
    check('latest_all_scalar_bounds_improve',all(F(a)<=F(b) for a,b in zip(f['physical_covariant_bounds'],latest['physical_covariant_bounds'])))

def independent_regular_inverse():
    x,t=S.symbols('x t',real=True);v=(1+2*x*x+3*x**4+x**8+x**10)/10**6
    source=S.expand(S.diff(v,x,2)+S.diff(v,x)/x);records=[]
    for n in range(9):
        sn=S.diff(source,x,n);rhs=sn-S.integrate(t**(n+1)*sn.subs(x,t*x),(t,0,1))
        residual=S.expand(S.diff(v,x,n+2)-rhs)
        check('independent_regular_inverse_order_'+str(n+2),residual==0)
        records.append({'field_order':n+2,'source_order':n,'residual':str(residual)})
    check('omitted_regular_integral_fails',S.expand(S.diff(v,x,2)-source)!=0)
    FIXTURES.append({'fixture':'independent_regular_radial_inverse','records':records,'exact_bolt_included':True})

def profiles(r,h):
    tau=h*mp.cosh(r);p=h*mp.sinh(r);z=1-tau;alpha=p/tau;ll=p/z
    v=(2-r*r+3*r**4+r**10)/10**6;vp=(-2*r+12*r**3+10*r**9)/10**6
    vpp=(-2+36*r*r+90*r**8)/10**6;vpr=(-2+12*r*r+10*r**8)/10**6
    rc=mp.mpf(1) if r==0 else r/mp.tanh(r)
    nu=lambda x:h*mp.cosh(x)*(2-x*x+3*x**4+x**10)/10**6
    beta_nup=tau*v+tau*rc*vpr
    direct=(mp.diff(nu,r,2)+beta_nup+2*ll*mp.diff(nu,r))/tau
    ss=2*tau/z+2*p*p/(z*z)
    divided=vpp+rc*vpr+(2*ll+2*alpha)*vp+z*ss/tau*v
    fb=lambda x:(1+x*x+x**6)/10**7
    ff=lambda x:h*mp.cosh(x)*fb(x)
    mass=mp.diff(ff,r)/tau-alpha*ff(r)/tau
    return direct,divided,mass,mp.diff(fb,r)

def independent_divided_equations():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for h in [mp.mpf('1e-8'),mp.mpf('1e-30')]:
                for at in [mp.mpf(0),mp.mpf('.3'),mp.mpf(1)]:
                    errors=[]
                    for n in [0,3,8]:
                        vals=[mp.diff(lambda r:profiles(r,h)[k],at,n)/mp.factorial(n) for k in range(4)]
                        error=max(abs(vals[0]-vals[1]),abs(vals[2]-vals[3]));errors.append(mp.nstr(error,18))
                        check('independent_division_'+str(h)+'_'+str(at)+'_'+str(dps)+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    records.append({'h':str(h),'rho':str(at),'dps':dps,'orders':[0,3,8],'errors':errors})
    FIXTURES.append({'fixture':'independent_regular_divided_equations','records':records,'exact_bolt_included':True,'actual_charged_source_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 low tangent and direct gap reconstruction','input_count':len(ps)}),flush=True)
    try:
        tm,tc,rad=inputs_and_coefficients();mid=middle(tm,rad['middle']);cold0=cold(tc,rad['cold']);args=tensors(mid,cold0)
        complete_scalar(args);independent_regular_inverse();independent_divided_equations()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 low tangent and direct gap reconstruction','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,'unchanged_scientific_algorithm_reused':PREVIOUS}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
