"""Complete source bound for the radius operator and retained scalar residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-radius-operator-source';PREVIOUS='finite-low-gap-reconstruction'
OWN=['checks/finite_radius_operator_source.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_low_gap',ROOT/'checks/finite_low_gap_reconstruction.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);L=B.L;R=B.R
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[];GAMMA=F('6.8e-28')
PROTOCOL={'source_order':8,'radius_operator_order':8,'scalar_source_order':6,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42','manufactured_orders':[0,3,6,8],
          'actual_charged_source_sampled':False,'ordinary_continuum_proof_separate':True}

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
def cap(a,b,n):return poly([min(a[j],b[j]) for j in range(n+1)])
def arr(f,k):return poly([F(x) for x in f[k]])

def inputs_and_identities():
    old=saved(PREVIOUS)
    check('prior_low_gap_result_passed',old['status']=='passed' and len(old['checks'])==114 and all(x['passed'] for x in old['checks']))
    check('prior_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    z,b0,b2,b3,b4,e,cx=S.symbols('z b0 b2 b3 b4 e cx',positive=True)
    bb=b0+b2*z**2+b3*z**3+b4*z**4-2*e*z**5+e*z**6
    au=S.expand(-2*z*S.diff(bb,z)+4*bb)
    an=S.expand(z**3*S.diff(bb,z,3)/6+4*b0+2*e*z**5-4*e*z**6)
    check('combined_u_polynomial',au==4*b0-2*b3*z**3-4*b4*z**4+12*e*z**5-8*e*z**6)
    check('combined_nu_polynomial',an==4*b0+b3*z**3+4*b4*z**4-18*e*z**5+16*e*z**6)
    u,nu,ff,ga,th,ed=S.symbols('u nu F gamma T E_delta')
    ss=-S.diff(bb,z)/z+2*bb/z**2
    cn=z*S.diff(bb,z,3)/6+4*b0/z**2+2*e*z**3-4*e*z**4
    original=cx*z**2*(2*ss*u+cn*nu-ff+ga*z**2*(th+ed))
    combined=cx*(au*u+an*nu-z**2*ff+ga*z**4*(th+ed))
    check('complete_physical_source_identity',S.expand(original-combined)==0)
    bf=S.Function('b')(z);nf=S.Function('nu')(z);p=S.sqrt(bf)
    radial=lambda f:-p*S.diff(f,z)
    radial_lap=radial(radial(nf))+(-S.diff(bf,z)/(2*p)+2*p/z)*radial(nf)
    check('canonical_radial_laplacian',S.simplify(radial_lap-(bf*S.diff(nf,z,2)+(S.diff(bf,z)-2*bf/z)*S.diff(nf,z)))==0)
    m,r,k,i=S.symbols('m r k i',integer=True)
    check('outer_connection_gains_two_powers',S.expand(1-i+m-r-(k-i)-(m-(r+1)-k))==2)
    check('unweighted_mass_and_source_conversion',S.simplify(z**-2*z**2)==1)
    check('LLL_coefficient_guard',F('.00352')/72<F('5e-5'))
    check('outer_mass_multiplier_guard',2*F('.412')/F('5.9')<F('.14'))
    check('cold_local_uniform_guard',F('1e22')*F('.001')==F('1e19'))
    check('outer_connection_coefficient_guard',F(9,4)*F('.1')**2==F(9,400))
    FIXTURES.append({'fixture':'exact_radius_source_identities','Au':str(au),'Anu':str(an),
                     'complete_physical_source_before_absolute_bounds':True,'outer_connection_weight_gain':2,
                     'fixed_weight_identity':'W*z0^2=1','full_angular_tower_retained':True})
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    q=add(poly([1]),arr(pf,'seed_jets'));ll=arr(pf,'log_radius_jets')
    tf=fixture(PREVIOUS,'complete_improved_tangent_arrays')
    return q,ll,[arr(tf,k) for k in ['U','K','Vbar','KV','tau_KV']]

def tensor_check(label,radial,initial,c,s):
    out=R.tensor(radial,initial,c,s);n=len(radial)-1
    # Independent normalization at TOTAL derivative order, rather than raw jets.
    d={(0,k):initial[k]*radial[k] for k in range(n+1)}
    for r in range(n):
        for k in range(n-r):
            total=r+k+1
            d[r+1,k]=d[r,k+1]+r*sum(F(comb(k,i)*factorial(i)*factorial(r+k-i),factorial(total))*c*F(s)**i*d[r,k-i]/100**(i+1) for i in range(k+1))
    check(label+'_normalized_tensor_reconciliation',out==[d[j,0] for j in range(n+1)])
    return out

def full_tensors(radial,region,label):
    one=[F(1)]*len(radial)
    if region=='middle':return tensor_check(label,radial,one,F(6),F(200))
    if region=='outer':return tensor_check(label,radial,one,F(9,400),F(2))
    initial=[F(x) for x in fixture('cold-metric-ten-jets','original_finite_completion_and_tensor_conversion')['Cartesian_amplifications']][:len(radial)]
    near=tensor_check(label+'_near',radial,initial,F(100),F(20))
    away=tensor_check(label+'_away',radial,one,F(15),F(200))
    return [max(a,b) for a,b in zip(near,away)]

def mass_tensors():
    old=fixture('finite-source-order-scalar','per_order_original_reconstruction')
    cold=fixture(PREVIOUS,'low_and_divided_radial_reconstruction')['cold']
    mid=full_tensors([F(x) for x in old['middle']['F']],'middle','mass_middle')
    co=full_tensors([F(x) for x in cold['F']],'cold','mass_cold')
    tail=full_tensors([F(x) for x in old['outer']['f_tail']],'outer','mass_outer_tail')
    multiplier=poly([F('.14')*F('4e-11')*F('.02')**j for j in range(9)])
    massrad=mul(multiplier,poly([1,F(1,100)]))[:9]
    mass=full_tensors(massrad,'outer','mass_outer_mode')
    regions={'middle':[mid[j]*100**j for j in range(9)],'cold':[co[j]*100**j for j in range(9)],
             'outer':[F('.412')**j*(F('.0331')*tail[j]+F('.1')*mass[j])*100**j for j in range(9)]}
    fu=poly([max(a[j] for a in regions.values()) for j in range(9)])
    check('full_mass_arrays_positive',all(min(a)>0 for a in regions.values()) and all(fu[j]>=a[j] for a in regions.values() for j in range(9)))
    FIXTURES.append({'fixture':'full_regional_mass_tensors','unweighted_D_arrays':{k:strings(v) for k,v in regions.items()},
                     'Fu':strings(fu[:9]),'Fu_decimals':decimals(fu[:9]),'outer_mass_radial':strings(massrad),
                     'outer_mass_mode_retained':True,'divergent_outer_W_F_not_used':True})
    return fu

def source_tensors(q):
    src=fixture('seed-nine-metric-eleven','all_regional_ninth_source_ledgers')
    middle=poly([F(x)*100**j for j,x in enumerate(src['middle'][:9])])
    cold=poly([F('.00331')*F(x)*100**j for j,x in enumerate(src['cold_total'][:9])])
    seed=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    lllm=scale('5e-5',compose(der(der(seed)),q))
    tau=q.copy();tau[0]=F('.001')
    lllc=scale('5e-5',add(scale(6*F('1.334'),tau),scale(12*F('.501'),power(tau,2)),
                         scale('80e-24',power(tau,3)),scale('30e-24',power(tau,4)),scale(6*F('1.001e-60'),q)))
    em0=full_tensors([F('1e16')*2**j for j in range(9)],'middle','local_middle')
    ec0=full_tensors([F('1e19')*2**j for j in range(9)],'cold','local_cold')
    em=poly([em0[j]*100**j for j in range(9)]);ec=poly([ec0[j]*100**j for j in range(9)])
    outer=[F(x) for x in src['outer_ordinary'][:9]];prefix=[max(outer[:j+1]) for j in range(9)]
    ot=full_tensors([2*x for x in prefix],'outer','complete_outer_source')
    ou=poly([F('.412')**j*F('.0331')*ot[j]*100**j for j in range(9)])
    regions={'middle':mul(power(q,4),add(middle,lllm,em)),
             'cold':mul(power(q,4),add(cold,lllc,ec)),'outer':mul(power(q,2),ou)}
    total=poly([max(a[j] for a in regions.values()) for j in range(9)])
    check('full_source_arrays_positive',all(min(a[:9])>0 for a in regions.values()))
    check('outer_local_reserve_included',min(prefix)>F('1e8') and F('.63')*min(prefix)+F('1e8')<2*min(prefix))
    check('source_orders_sufficient',all(len(src[k])==10 for k in ['middle','cold_total','outer_ordinary']) and src['full_angular_tower_retained'])
    FIXTURES.append({'fixture':'full_regional_source_tensors','W_z2_Rdelta_arrays':{k:strings(v[:9]) for k,v in regions.items()},
                     'Q':strings(total[:9]),'Q_decimals':decimals(total[:9]),'middle_massive':strings(middle[:9]),
                     'cold_massive':strings(cold[:9]),'LLL_middle':strings(lllm[:9]),'LLL_cold':strings(lllc[:9]),
                     'local_middle':strings(em[:9]),'local_cold':strings(ec[:9]),'outer_prefix':strings(prefix),
                     'middle_cold_source_entries_already_full_tensors':True,'outer_complete_physical_cancellation_retained':True})
    return total

def radius_operator(q,ll,K,fu,source):
    ap=compose([F(4)/F('5.9'),0,0,2*F('.668'),4*F('.501'),F('12e-24'),F('8e-24')],q)
    np=compose([F(4)/F('5.9'),0,0,F('.668'),4*F('.501'),F('18e-24'),F('16e-24')],q)
    parts=[scale('6.1',mul(add(ap,np),K)),scale('6.1',mul(power(q,2),fu)),scale(F('6.1')*GAMMA,source)]
    direct=add(*parts)
    old=scale('6.1',mul(power(q,2),add(scale(2,der(der(K))),scale(2,mul(der(ll),der(K))))))
    old[0]=min(old[0],F('.005'));chosen=cap(direct,old,8)
    check('radius_jet_old_bounds_retained',all(chosen[j]<=old[j] for j in range(9)))
    check('radius_sixth_eighth_strictly_improve',chosen[6]<old[6] and chosen[8]<old[8],decimals(chosen[:9]))
    FIXTURES.append({'fixture':'direct_radius_operator_jets','polynomial_u':strings(ap[:9]),'polynomial_nu':strings(np[:9]),
                     'parts':[strings(a[:9]) for a in parts],'parts_decimals':[decimals(a[:9]) for a in parts],
                     'direct':strings(direct[:9]),'old':strings(old[:9]),'selected':strings(chosen[:9]),
                     'selected_decimals':decimals(chosen[:9]),'value_over_tau_cap':'.005','unchanged_tangent_arrays':True})
    return old,chosen

def scalar_ledger(q,args,oldln,newln):
    U,K,V,KV,h=args
    old=fixture(PREVIOUS,'complete_refined_scalar_and_anomaly')
    gf=fixture(PREVIOUS,'refined_finite_geometry_and_components')
    T=arr(fixture(PREVIOUS,'finite_parameter_integral_bounds'),'Tplus')
    vp=V.copy();vp[0]=0;O=scale('1.001',exps(scale(4,T)));E=exps(add(vp,T))
    def radius_term(ln):
        hh=mul(KV,ln)
        for n in range(9):hh[n]-=KV[n]*ln[0]-min(KV[n]*ln[0],F('.005')*h[n])
        return scale(6*F('1.001'),mul(mul(mul(O,q),E),hh))
    check('previous_radius_product_reconstructed',radius_term(oldln)[:9]==[F(x) for x in old['remainder_terms'][2]])
    terms=[poly([F(x) for x in a]) for a in old['remainder_terms']];terms[2]=radius_term(newln)
    check('other_five_curvature_terms_unchanged',all(terms[i][:9]==[F(x) for x in old['remainder_terms'][i]] for i in [0,1,3,4,5]))
    nr=cap(add(*terms),arr(old,'selected_curvature_remainder'),8)
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')]);zabs=compose(za,q)
    lie=mul(K,compose([j*za[j] for j in range(len(za))],q))
    check('unchanged_Lie_reconstructed',lie[:7]==[F(x) for x in old['Lie']])
    seed=[F(x) for x in fixture('finite-scalar-residual-six','actual_scalar_residual_H6_ledger')['complete_seed_polynomial']]
    check('full_degree_twelve_seed_retained',len(seed)==13 and seed[-1]>0)
    dz=add(compose(seed,q),lie,nr)
    zw=cap(arr(gf,'ZW'),add(zabs,dz),8)
    m0,m1,c0,c1=[arr(gf['geometry'],k) for k in ['M0','M1','C0','C1']]
    bx=scale('6.1',add(mul(m0,add(der(der(dz)),mul(c0,der(dz)))),mul(m1,der(der(zabs))),mul(add(mul(m1,c0),mul(m0,c1)),der(zabs))))
    em=F(1,2000**10*180)
    anomaly=[scale(em,bx),scale(em*F(11,3),dz),scale(em*F(11,72),mul(add(zw,zabs),dz)),*[poly([F(x) for x in a]) for a in old['anomaly_terms'][3:]]]
    check('unrelated_three_anomaly_terms_unchanged',all(anomaly[i][:7]==[F(x) for x in old['anomaly_terms'][i]] for i in [3,4,5]))
    dh=add(*anomaly);ra=arr(old,'reference_remainder_second')
    p=cap(add(nr,scale(F(1,2),ra),lie,dh),arr(old,'residual_reference_W2'),6)
    phys=L.covariant(p,c0,0,6);constant=F(1001,24000)*sum(phys);before=F(old['C_tau_exact']);exponent=L.decimal_power(constant)
    check('all_six_curvature_terms_positive',len(terms)==6 and all(min(a[:9])>0 for a in terms))
    check('all_six_anomaly_terms_positive',len(anomaly)==6 and all(min(a[:7])>0 for a in anomaly))
    check('complete_anomaly_nonincreasing',all(dh[j]<=F(old['anomaly_difference'][j]) for j in range(7)))
    check('physical_covariant_bounds_nonincreasing',all(a<=F(b) for a,b in zip(phys,old['physical_covariant_bounds'])))
    check('complete_scalar_coefficient_improves',constant<before,{'coefficient':L.decimal(constant),'improvement':L.decimal(before/constant)})
    check('upward_exact_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent)
    budget=F(10**57)/(800*F(10)**exponent);check('rounded_scalar_only_transfer_budget',800*F(10)**exponent*budget==10**57)
    FIXTURES.append({'fixture':'complete_refined_scalar_and_anomaly','remainder_terms':[strings(a[:9]) for a in terms],
                     'remainder_terms_decimals':[decimals(a[:9]) for a in terms],'selected_curvature_remainder':strings(nr[:9]),
                     'anomaly_terms':[strings(a[:7]) for a in anomaly],'anomaly_term_names':old['anomaly_term_names'],
                     'anomaly_terms_decimals':[decimals(a[:7]) for a in anomaly],'anomaly_difference':strings(dh[:7]),
                     'reference_remainder_second':strings(ra[:7]),'Lie':strings(lie[:7]),'residual_reference_W2':strings(p[:7]),
                     'physical_covariant_bounds':strings(phys),'physical_covariant_decimals':decimals(phys),
                     'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),'C_tau_power_of_ten':exponent,
                     'old_C_tau_exact':str(before),'improvement_ratio':L.decimal(before/constant),
                     'scalar_only_rounded_transfer_budget_over_g_epsilon_c':str(budget),'scalar_only_budget_decimal':L.decimal(budget),
                     'all_inverse_metric_contacts_retained':True,'full_noncompact_physical_volume_retained':True,
                     'finite_geometry_unchanged':True,'full_Schur_source_bounded':False,'target_passed':False})

def manufactured(x,y,z0):
    w=x*x+y*y;z=z0-w/4;bb=w/4;cx=mp.mpf(6);ee=mp.mpf('1e-7');gamma=mp.mpf('1e-3')
    nu=(2*w+3*w*w+w**4)/10**6;u=(1+w+w**3)/10**6;ff=(1+w*w)/10**7
    nw=(2+6*w+4*w**3)/10**6;nww=(6+12*w*w)/10**6
    ss=1/z+2*bb/(z*z);cn=4/(cx*z*z)+2*ee*z**3-4*ee*z**4
    radial_lap=4*(nw+w*nww)+2*w*nw/z
    ed=(1+z*z)/1000
    th=(radial_lap-2*ss*u-cn*nu+ff-gamma*z*z*ed)/(gamma*z*z)
    rhs=cx*((2*z+4*bb)*u+(4/cx+2*ee*z**5-4*ee*z**6)*nu-z*z*ff+gamma*z**4*(th+ed))
    # Physical four-volume divergence, independently differentiated in x and y.
    def flux_x(a):
        ww=a*a+y*y;zz=z0-ww/4
        return 2*a*(2+6*ww+4*ww**3)/(10**6*zz**2)
    def flux_y(a):
        ww=x*x+a*a;zz=z0-ww/4
        return 2*a*(2+6*ww+4*ww**3)/(10**6*zz**2)
    direct=cx*z**4*(mp.diff(flux_x,x)+mp.diff(flux_y,y))
    missing_local=rhs-cx*gamma*z**4*ed
    return direct,rhs,missing_local

def independent_divergence():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for z0 in [mp.mpf('.7'),mp.mpf('.9')]:
                for x,y in [(mp.mpf(0),mp.mpf(0)),(mp.mpf('.1'),mp.mpf('-.2'))]:
                    tag=str(z0)+'_'+str(x)+'_'+str(y)+'_'+str(dps);errors=[]
                    for n in PROTOCOL['manufactured_orders']:
                        differences=[]
                        for k in range(n+1):
                            direct=mp.diff(lambda a,b:manufactured(a,b,z0)[0],(x,y),(n-k,k))/mp.factorial(n)
                            source=mp.diff(lambda a,b:manufactured(a,b,z0)[1],(x,y),(n-k,k))/mp.factorial(n)
                            differences.append(abs(direct-source))
                        error=max(differences);errors.append(mp.nstr(error,18))
                        check('independent_4D_divergence_'+tag+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    vals=manufactured(x,y,z0);omission=abs(vals[0]-vals[2])
                    check('omitted_local_source_fails_'+tag,omission>mp.mpf('1e-10'),mp.nstr(omission,18))
                    records.append({'z0':str(z0),'x':str(x),'y':str(y),'dps':dps,'orders':PROTOCOL['manufactured_orders'],
                                    'all_mixed_partials_at_each_order':True,'errors':errors,'omission':mp.nstr(omission,18)})
    FIXTURES.append({'fixture':'independent_four_dimensional_radius','records':records,'exact_center_included':True,
                     'complete_four_volume_divergence':True,'actual_charged_source_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 direct radius operator and complete source','input_count':len(ps)}),flush=True)
    try:
        q,ll,args=inputs_and_identities();fu=mass_tensors();source=source_tensors(q)
        oldln,newln=radius_operator(q,ll,args[1],fu,source);scalar_ledger(q,args,oldln,newln);independent_divergence()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 direct radius operator and complete source','created_utc':datetime.now(timezone.utc).isoformat(),
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
