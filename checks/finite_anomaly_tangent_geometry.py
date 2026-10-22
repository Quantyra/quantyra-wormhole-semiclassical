"""Finite path jets from sharper tangents and complete scalar/anomaly ledger."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from itertools import combinations
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-anomaly-tangent-geometry';PREVIOUS='finite-source-order-scalar'
OWN=['checks/finite_anomaly_tangent_geometry.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_per_order_scalar',ROOT/'checks/finite_source_order_scalar.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);L=B.L
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[]
PROTOCOL={'metric_order':10,'scalar_remainder_order':8,'scalar_source_order':6,'finite_updates':6,
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
def cap(a,b,n):return poly([min(a[j],b[j]) for j in range(n+1)])

def inputs():
    old=saved(PREVIOUS)
    check('prior_per_order_result_passed',old['status']=='passed' and len(old['checks'])==42 and all(x['passed'] for x in old['checks']))
    check('prior_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    f=fixture(PREVIOUS,'full_tangent_and_tau_arrays')
    arrays=[poly([F(x) for x in f[k]]) for k in ['U','K','Vbar','KV','tau_KV']]
    check('small_exponential_value',1/(1-F('1.2e-8')-F('3e-8'))<F('1.001'))
    check('ten_original_tangent_orders',all(len(f[k])==11 for k in ['U','K','Vbar','KV','tau_KV']) and f['analytic_mass_remainder_and_cubic_retained'])
    return arrays

def improve(T,U,V,order=10):
    Vp=V.copy();Vp[0]=0
    candidate=scale('1.001',mul(exps(add(Vp,T)),U))
    return poly([0]+[min(T[j],max(U[j],candidate[j])) for j in range(1,order+1)])

def finite_bounds(U,V):
    t=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    initial=t.copy();history=[strings(t[:11])]
    for j in range(6):
        nt=improve(t,U,V)
        check('finite_bound_step_'+str(j+1),all(U[k]<=nt[k]<=t[k] for k in range(1,11)))
        t=nt;history.append(strings(t[:11]))
    check('finite_bounds_strictly_improve',all(0<t[j]<initial[j] for j in range(1,11)))
    FIXTURES.append({'fixture':'finite_parameter_integral_bounds','initial':strings(initial[:11]),'Tplus':strings(t[:11]),
                     'Tplus_decimals':decimals(t[:11]),'history':history,'finite_updates':6,'convergence_used':False})
    return t

def geometry(T,K,V,q,ll):
    vp=V.copy();vp[0]=0
    s1=scale('1.001',mul(exps(add(vp,T)),K))
    m0=scale('3.003',exps(add(scale(2,ll),scale(4,T))));m1=scale(2,mul(m0,add(K,s1)))
    c0=scale(F(3,2),mul(m0,der(m0)));c1=scale(F(3,2),add(mul(m1,der(m0)),mul(m0,der(m1))))
    seedB=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    ref=add(compose(der(der(seedB)),q),poly(['2e-10']))
    r0=add(ref,scale(2,der(c0)),scale(2,power(c0,2)));r1=add(scale(2,der(c1)),scale(4,mul(c0,c1)))
    ri0=scale(2,r0);ri1=scale(2,r1);sc0=mul(m0,ri0);sc1=add(mul(m1,ri0),mul(m0,ri1))
    j1=scale('6.1',add(ri1,scale(F(1,4),add(mul(sc1,m0),mul(sc0,m1)))))
    w1=scale('6.1',add(mul(m1,r0),mul(m0,r1),scale(2,add(mul(m1,ri0),mul(m0,ri1))),scale(F(2,3),add(mul(sc1,power(m0,2)),scale(2,mul(mul(sc0,m0),m1))))))
    return {'M0':m0,'M1':m1,'C0':c0,'C1':c1,'J1':j1,'CW1':w1,'Sc1':scale('6.1',sc1),'S1':s1}

def contacts(g,jw,cw,T,q):
    m0,m1=g['M0'],g['M1']
    dc=add(scale(2,mul(power(m0,4),mul(cw,g['CW1']))),scale(4,mul(mul(power(m0,3),m1),power(cw,2))))
    dj=add(scale(2,mul(power(m0,2),mul(jw,g['J1']))),scale(2,mul(mul(m0,m1),power(jw,2))))
    mu=scale(F('4e-27')*F('1.001')/18,mul(exps(scale(4,T)),power(q,4)))
    dm=scale(4,mul(mu,g['S1']))
    return scale(F(7,12*2000**10*180),dc),scale(F(11,6*2000**10*180),dj),dm

def restore_old(pf,q,ll):
    t=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    k=poly([F('1e-4')]+[F(10)*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    v=poly([F('1.2e-8')]+[F(100)*100**j for j in range(1,6)]+[F('2e26')*100**j for j in range(6,11)])
    g=geometry(t,k,v,q,ll);old=fixture('finite-scalar-residual-six','weighted_first_and_second_path_geometry')
    for key,src in [('M0',pf['metric']),('C0',pf['connection']),('M1',old['metric_first']),('C1',old['connection_first'])]:
        n=10 if key.startswith('M') else 9
        check('old_'+key+'_recovered',g[key][:n+1]==[F(x) for x in src[:n+1]])
    jw=poly([F(x) for x in pf['weighted_J_initial']]);cw=poly([F(x) for x in pf['weighted_C_initial']])
    parts=contacts(g,jw,cw,t,q)
    oldrest=fixture(PREVIOUS,'per_order_scalar_and_anomaly_ledger')['retained_anomaly_rest']
    check('old_anomaly_contacts_reconciled',add(*parts)[:7]==[F(x) for x in oldrest])
    return g

def scalar_ledger(T,U,K,V,KV,h):
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    old=fixture(PREVIOUS,'per_order_scalar_and_anomaly_ledger')
    original=fixture('finite-scalar-residual-six','actual_scalar_residual_H6_ledger')
    zeta=poly([F(x) for x in pf['seed_jets']]);q=add(poly([1]),zeta);ll=poly([F(x) for x in pf['log_radius_jets']])
    oldg=restore_old(pf,q,ll);g=geometry(T,K,V,q,ll)
    for key in ['M0','M1','C0','C1','J1','CW1']:
        check('new_'+key+'_improves',all(g[key][j]<oldg[key][j] for j in range(1,7)))
    vp=V.copy();vp[0]=0;up=U.copy();up[0]=0;O=scale('1.001',exps(scale(4,T)));E=exps(add(vp,T))
    s1=g['S1'];s1u=scale('1.001',mul(E,U));s2=scale('1.001',mul(mul(q,exps(add(vp,scale(2,T)))),mul(K,KV)))
    w1=add(K,s1);j2=mul(O,add(mul(add(scale(4,power(w1,2)),scale(2,s2)),add(poly([1]),scale(4,U))),scale(16,mul(K,w1))))
    seedB=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')]);zabs=compose(za,q)
    lie=mul(K,compose([j*za[j] for j in range(len(za))],q));seed=[F(x) for x in original['complete_seed_polynomial']]
    check('full_seed_polynomial_retained',len(seed)==13 and seed[-1]>0)
    zd=add(compose(seed,q),lie)
    ln=scale('6.1',mul(power(q,2),add(scale(2,der(der(K))),scale(2,mul(der(ll),der(K))))));ln[0]=min(ln[0],F('.005'))
    lt=scale('6.1',add(scale(2,mul(q,compose(seedB,q))),mul(power(q,2),compose(der(seedB),q))))
    hh=mul(KV,ln)
    for n in range(9):hh[n]-=KV[n]*ln[0]-min(KV[n]*ln[0],F('.005')*h[n])
    gg=B.B.gap_product(KV,h,q,order=8);gg=cap(gg,mul(q,power(der(KV),2)),8)
    qq=scale('.501',mul(E,power(KV,2)))
    terms=[scale(F(1,2),mul(add(poly([12]),zabs),j2)),scale(2,mul(mul(O,add(U,s1u)),zd)),
           scale(6*F('1.001'),mul(mul(mul(O,q),E),hh)),scale(6,mul(mul(O,qq),lt)),
           scale(6*F('6.1')*F('1.001'),mul(mul(mul(O,E),power(q,2)),gg)),
           scale(F('4.004')*F('1.3e-10'),mul(mul(power(q,2),exps(add(scale(4,T),scale(2,up)))),power(K,2)))]
    nr=cap(add(*terms),poly([F(x) for x in old['selected_curvature_remainder']]),8);dz=add(zd,nr)
    ip=poly([1]);zp=poly([1])
    for j in range(1,12):zp=mul(zp,zeta);ip=add(ip,zp)
    pj=add(scale('.501',power(q,2)),scale('5e-24',power(q,3)),scale('4.5e-24',power(q,4)))
    pc=add(scale(3*F('.668'),ip),poly([6*F('.501')]),scale('20e-24',q),scale('15e-24',power(q,2)))
    jw=cap(add(g['J1'],scale(13,pj)),poly([F(x) for x in pf['weighted_J_initial']]),8)
    cw=cap(add(g['CW1'],scale(25,pc)),poly([F(x) for x in pf['weighted_C_initial']]),8)
    zw=cap(add(g['Sc1'],zabs),poly([F(x) for x in pf['weighted_Z_initial']]),8);zw=cap(zw,add(zabs,dz),8)
    m0,m1,c0,c1=[g[k] for k in ['M0','M1','C0','C1']]
    bx=scale('6.1',add(mul(m0,add(der(der(dz)),mul(c0,der(dz)))),mul(m1,der(der(zabs))),mul(add(mul(m1,c0),mul(m0,c1)),der(zabs))))
    em=F(1,2000**10*180)
    anomaly=[scale(em,bx),scale(em*F(11,3),dz),scale(em*F(11,72),mul(add(zw,zabs),dz)),*contacts(g,jw,cw,T,q)]
    dh=add(*anomaly);ra=poly()
    for n,a in [(5,4),(6,6)]:ra=add(ra,scale(F('1.001')*F('1.3e-22')*a,mul(mul(power(q,n),exps(scale(n,T))),add(scale(n*n,power(s1,2)),scale(n,s2)))))
    raw=add(nr,scale(F(1,2),ra),lie,dh);p=cap(raw,poly([F(x) for x in old['residual_reference_W2']]),6)
    phys=L.covariant(p,c0,0,6);constant=F(1001,24000)*sum(phys);before=F(old['C_tau_exact']);exponent=L.decimal_power(constant)
    check('complete_anomaly_improves',all(dh[j]<=F(old['anomaly_difference'][j]) for j in range(7)),decimals(dh[:7]))
    check('all_anomaly_terms_positive',all(min(a[:7])>0 for a in anomaly))
    check('all_curvature_remainder_terms_retained',len(terms)==6 and all(min(a[:9])>0 for a in terms))
    check('actual_covariant_bounds_improve',all(a<=F(b) for a,b in zip(phys,old['physical_covariant_bounds'])))
    check('complete_scalar_coefficient_improves',constant<before,{'coefficient':L.decimal(constant),'improvement':L.decimal(before/constant)})
    check('upward_exact_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent)
    budget=F(10**57)/(800*F(10)**exponent);check('rounded_transfer_budget',800*F(10)**exponent*budget==10**57)
    FIXTURES.append({'fixture':'refined_finite_geometry_and_components','geometry':{k:strings(v[:9]) for k,v in g.items()},
                     'JW':strings(jw[:9]),'CW':strings(cw[:9]),'ZW':strings(zw[:9]),'weighted_component_seed_terms_retained':True,
                     'metric_order':10,'curvature_order':8,'endpoint_map_differentiated':False})
    FIXTURES.append({'fixture':'complete_refined_scalar_and_anomaly','remainder_terms':[strings(a[:9]) for a in terms],
                     'remainder_terms_decimals':[decimals(a[:9]) for a in terms],'selected_curvature_remainder':strings(nr[:9]),
                     'anomaly_terms':[strings(a[:7]) for a in anomaly],'anomaly_term_names':['Laplacian','linear_trace','trace_square','Weyl_square','Ricci_square','magnetic'],
                     'anomaly_terms_decimals':[decimals(a[:7]) for a in anomaly],'anomaly_difference':strings(dh[:7]),
                     'reference_remainder_second':strings(ra[:7]),'Lie':strings(lie[:7]),'residual_reference_W2':strings(p[:7]),
                     'physical_covariant_bounds':strings(phys),'physical_covariant_decimals':decimals(phys),
                     'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),'C_tau_power_of_ten':exponent,
                     'old_C_tau_exact':str(before),'improvement_ratio':L.decimal(before/constant),
                     'scalar_only_rounded_transfer_budget_over_g_epsilon_c':str(budget),'scalar_only_budget_decimal':L.decimal(budget),
                     'all_inverse_metric_contacts_retained':True,'full_noncompact_physical_volume_retained':True,
                     'full_Schur_source_bounded':False,'target_passed':False})

def mfrac(x):return mp.mpf(x.numerator)/x.denominator

def independent_finite_log():
    vb=F('3.375e-5');wb=F('.85')*vb/(1-vb)
    check('manufactured_complex_log_initial_bound',wb/(1-wb)<F('.001') and 1/(1-F('.00004'))<F('1.001'))
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for tau0 in [F('.2'),F('.8')]:
                v=poly([F(2,100000),F(-1,100000),F(3,100000),F(1,100000)])
                nu=mul(poly([tau0,F('.1')]),v);u=poly([F(1,100000),F(2,100000),F(1,100000),0,F(1,100000)])
                U=poly([max(abs(a),abs(b)) for a,b in zip(u,nu)]);V=poly([abs(x) for x in v])
                t=poly([0]+[F('.001')*2**j for j in range(1,9)])
                for _ in range(6):t=improve(t,U,V,8)
                for ee in [mp.mpf('.3'),mp.mpf('.9')]:
                    def direct(x):return mp.log(1+(mfrac(tau0)+x/10)*mp.expm1(ee*(2-x+3*x*x+x**3)/100000))
                    ratios=[abs(mp.diff(direct,mp.mpf(0),j))/mp.factorial(j)/mfrac(t[j]) for j in range(1,9)]
                    maximum=max(ratios);check('independent_finite_log_'+str(tau0)+'_'+str(ee)+'_'+str(dps),maximum<1,mp.nstr(maximum,18))
                    records.append({'tau0':str(tau0),'ell':str(ee),'dps':dps,'order':8,'max_bound_ratio':mp.nstr(maximum,18)})
    FIXTURES.append({'fixture':'independent_finite_log_derivatives','records':records,'actual_charged_source_sampled':False})

PAIRS=list(combinations(range(4),2))
def warped(x,e,t0):
    ta=t0+x/100;z=1-ta;phi=-mp.log(z);phip=1/(100*z);phipp=1/(10000*z*z)
    u=(1+2*x+x*x+x**4)/1000;up=(2+2*x+4*x**3)/1000;upp=(2+12*x*x)/1000
    v=(2-x+3*x*x+x**3)/1000;vp=(-1+6*x+3*x*x)/1000;vpp=(6+6*x)/1000
    k=mp.exp(e*v);T=1+ta*(k-1);Tx=(k-1)/100+ta*k*e*vp;Txx=2*k*e*vp/100+ta*k*(e*vpp+e*e*vp*vp)
    ss=mp.log(T);sp=Tx/T;spp=Txx/T-sp*sp
    a=phi+e*u+ss;b=phi+ss;ap=phip+e*up+sp;app=phipp+e*upp+spp;bp=phip+sp;bpp=phipp+spp
    ia=mp.exp(-2*a);ib=mp.exp(-2*b);inverse=[ia,ia,ib,ib]
    section=[-ia*app,-ia*(bpp+bp*bp-ap*bp),-ia*(bpp+bp*bp-ap*bp),-ia*ap*bp,-ia*ap*bp,ib-ia*bp*bp]
    ric=[sum(section[k] for k,pair in enumerate(PAIRS) if i in pair) for i in range(4)];sc=sum(ric)
    jt=[r-sc/4 for r in ric];ct=[section[k]-(ric[i]+ric[j])/2+sc/6 for k,(i,j) in enumerate(PAIRS)]
    jcov=[jt[i]/inverse[i] for i in range(4)];ccov=[ct[k]/(inverse[i]*inverse[j]) for k,(i,j) in enumerate(PAIRS)]
    return inverse,jcov,ccov,4*sum(c*c for c in ct),sum(j*j for j in jt)

def contraction_derivatives(x,e,t0):
    inv,jc,cc,_,_=warped(x,e,t0)
    di=[mp.diff(lambda a:warped(x,a,t0)[0][i],e) for i in range(4)]
    dj=[mp.diff(lambda a:warped(x,a,t0)[1][i],e) for i in range(4)]
    dc=[mp.diff(lambda a:warped(x,a,t0)[2][i],e) for i in range(6)]
    cf=mp.mpf(0);cm=mp.mpf(0);jf=mp.mpf(0);jm=mp.mpf(0)
    for k,(i,j) in enumerate(PAIRS):
        inds=[i,j,i,j];weight=mp.fprod(inv[a] for a in inds)
        cf+=8*cc[k]*dc[k]*weight
        cm+=4*cc[k]**2*sum(di[inds[s]]*mp.fprod(inv[inds[r]] for r in range(4) if r!=s) for s in range(4))
    for i in range(4):
        jf+=2*jc[i]*dj[i]*inv[i]**2
        jm+=jc[i]**2*sum(di[i]*mp.fprod([inv[i]]*(2-1)) for _ in range(2))
    return cf+cm,jf+jm,cm,jm

def independent_curvature_contacts():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for t0 in [mp.mpf('.2'),mp.mpf('.8')]:
                for ee in [mp.mpf('.3'),mp.mpf('.9')]:
                    tag=str(t0)+'_'+str(ee)+'_'+str(dps);errors=[]
                    for n in [0,3,6]:
                        direct=[mp.diff(lambda x:mp.diff(lambda a:warped(x,a,t0)[k],ee),mp.mpf(0),n)/mp.factorial(n) for k in [3,4]]
                        expanded=[mp.diff(lambda x:contraction_derivatives(x,ee,t0)[k],mp.mpf(0),n)/mp.factorial(n) for k in [0,1]]
                        error=max(abs(a-b) for a,b in zip(direct,expanded));errors.append(mp.nstr(error,18))
                        check('independent_curvature_contacts_'+tag+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    omission=contraction_derivatives(mp.mpf(0),ee,t0)[2:]
                    check('omitted_inverse_contacts_fail_'+tag,min(abs(x) for x in omission)>mp.mpf('1e-15'))
                    records.append({'tau0':str(t0),'ell':str(ee),'dps':dps,'orders':[0,3,6],'errors':errors,
                                    'Weyl_omission':mp.nstr(abs(omission[0]),18),'Ricci_omission':mp.nstr(abs(omission[1]),18)})
    FIXTURES.append({'fixture':'independent_warped_curvature_contacts','records':records,'Weyl_inverse_slots':4,'Ricci_inverse_slots':2,'actual_charged_geometry':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite tangent geometry and anomaly','input_count':len(ps)}),flush=True)
    try:
        U,K,V,KV,h=inputs();t=finite_bounds(U,V);scalar_ledger(t,U,K,V,KV,h);independent_finite_log();independent_curvature_contacts()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 finite tangent geometry and anomaly','created_utc':datetime.now(timezone.utc).isoformat(),
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
