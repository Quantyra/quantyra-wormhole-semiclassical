"""Per-order complete source reconstruction and scalar anomaly refinement."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-source-order-scalar';PREVIOUS='finite-scalar-cold-gap-product'
OWN=['checks/finite_source_order_scalar.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_gap_products',ROOT/'checks/finite_scalar_cold_gap_product.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);R=B.R;L=B.L
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[]
PROTOCOL={'scalar_source_order':6,'scalar_remainder_reference_order':8,'source_order':8,'tangent_order':10,
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
def conv(a,b,n):return [sum(a[i]*b[j-i] for i in range(j+1) if i<len(a) and j-i<len(b)) for j in range(n)]
def coef(c,n=11):return [F(c)*F('.02')**j for j in range(n)]
def monomial(m,n=11):return [F(str(abs(S.binomial(m,j))))/100**j for j in range(n)]
def primitive(a,m,n):return [a[0]*F(m+2,(m+1)**2)]+[a[j-1]/(100*j) for j in range(1,n)]

def complete_sources():
    old=saved(PREVIOUS)
    check('previous_gap_result_passed',old['status']=='passed' and len(old['checks'])==86 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    f=fixture('seed-nine-metric-eleven','all_regional_ninth_source_ledgers')
    tm=[F(10**20)]+[min(F(2*10**44),F(f['middle'][j])+2**j) for j in range(1,9)]
    tc=[min(F(2*10**44),F(f['cold_total'][j])+2**j) for j in range(9)]
    ko=[max(F(x) for x in f['outer_ordinary'][:j+1]) for j in range(9)]
    check('cold_LLL_complex_polynomial',6*F('1.334')*4+12*F('.501')*F('.004')*4+F('1e-20')<33)
    check('cold_LLL_weighted_value',F('.00352')*33/72<1)
    check('middle_sources_below_old',all(x<2*10**44 for x in tm[1:]))
    check('cold_sources_below_old',all(x<2*10**44 for x in tc))
    check('outer_prefix_sources_below_old',all(F('1e8')<x<10**44 for x in ko) and all(a<=b for a,b in zip(ko,ko[1:])))
    check('LLL_not_dropped',all(tm[j]>F(f['middle'][j]) and tc[j]>F(f['cold_total'][j]) for j in range(1,9)))
    FIXTURES.append({'fixture':'complete_per_order_sources','middle':strings(tm),'cold':strings(tc),'outer_prefix':strings(ko),
                     'middle_decimals':decimals(tm),'cold_decimals':decimals(tc),'outer_prefix_decimals':decimals(ko),
                     'LLL_and_local_allocations_retained':True,'new_cutoff_used':False})
    return tm,tc,ko

def middle(t):
    g=F('6.8e-28');p=[F(10**20)];f=[F('1e-8')]
    for n in range(8):
        p.append(sum(4*2**i*(p[n-i]+t[n-i]/2) for i in range(n+1))/(100*(n+1)))
        f.append((sum(40*2**i*f[n-i] for i in range(n+1))+2*g*sum(2**i*(p[n-i]+t[n-i]+10**16*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    u=[F('1e-4')]*6;v=u.copy()
    for n in range(4,9):
        c=lambda k,a:sum(k*2**i*a[n-i] for i in range(n+1))
        d=lambda k,a:sum(k*2**i*100*(n-i+1)*a[n-i+1] for i in range(n+1))
        u.append((d(2,u)+c(12,u)+c(10,v)+3*f[n]+3*g*c(2,t)+g*10**18*2**n)/(10000*(n+1)*(n+2)))
        v.append((d(100,v)+c(20000,u)+c(20000,v)+f[n]+g*c(2,t)+g*10**18*2**n)/(10000*(n+1)*(n+2)))
    return {'P':p,'F':f,'U':u,'N':v}

def cold(t):
    g=F('6.8e-28');qn=[5*sum(F('.1')**i*t[j-i] for i in range(j+1)) for j in range(9)];qa=[t[0]/2]
    for n in range(8):qa.append(sum(5*2**i*(t[n-i]+2*qa[n-i]) for i in range(n+1))/(100*(n+1)))
    q=[max(a,b) for a,b in zip(qn,qa)];f=[F('1e-10')]
    for n in range(8):f.append((sum(2**i*f[n-i] for i in range(n+1))+2*g*sum(2**i*(q[n-i]+t[n-i]+10**22*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    u=[F('1e-4')]*6;v=u.copy()
    for n in range(4,9):
        c=lambda k,a:sum(k*2**i*a[n-i] for i in range(n+1))
        d=lambda k,a:sum(k*2**i*100*(n-i+1)*a[n-i+1] for i in range(n+1))
        u.append(2*(d(12,u)+c(12,u)+c(10,v)+3*f[n]+3*g*c(2,t)+g*10**24*2**n)/(10000*(n+1)*(n+2)))
        v.append(2*(d(12,v)+c(2,u)+c(10,v)+f[n]+g*c(2,t)+g*10**24*2**n)/(10000*(n+1)*(n+2)))
    return {'Ward_near':qn,'Ward_away':qa,'F':f,'U':u,'N':v}

def outer(k):
    gx=F('4e-27');g=F('6.8e-28');mu=[F('7e-11')]+[8*gx*k[j-1]/(100*j) for j in range(1,10)]
    f=conv(conv(coef('.14'),monomial(1),9),mu,9);small=conv(conv(coef('.15'),monomial(3),10),mu,10)
    jj=primitive(small,4,11);w=conv(conv(coef('.51'),monomial(-1),11),jj,11)
    w=[min(v,F('1e-8')) if j<=5 else v for j,v in enumerate(w)];v=[F('1e-8')]*6
    for n in range(4,9):
        dv=[100*(j+1)*v[j+1] for j in range(n+1)]
        drift=conv(conv(coef(1),monomial(1),n+1),dv,n+1)[n]/100
        lie=conv(conv(coef(4),monomial(1),n+1),w,n+1)[n]/1000
        v.append((3*f[n]+6*g*k[n]+g*10**8*F('.02')**n+lie+drift)/(10000*(n+1)*(n+2)))
    ut=conv(coef('1.001'),primitive(conv(coef('1.001'),v,11),4,11),11)
    term=conv(conv(coef('.5'),monomial(1),11),ut,11);u=[v[j]+term[j]/100 for j in range(11)]
    term=conv(conv(coef('1.01'),monomial(-1),11),ut,11);nu=[w[j]+term[j] for j in range(11)]
    gap=conv(coef('1.07'),nu,11)
    return {'mu':mu,'f_tail':f,'proper_v':v,'proper_w':w,'common_u':u,'common_nu':nu,'common_V':gap}

def reconstruction(tm,tc,ko):
    specs=[('middle-metric-ten-jets','radial_two_derivative_reconstruction',middle,[F(10**20)]+[F(2*10**44)]*8,tm),
           ('cold-metric-ten-jets','weighted_radial_reconstruction',cold,[F(2*10**44)]*9,tc),
           ('outer-metric-ten-jets','outer_two_derivative_reconstruction',outer,[F(10**44)]*9,ko)]
    data=[]
    for stem,fn,func,old_input,new_input in specs:
        old=fixture(stem,fn);recovered=func(old_input);new=func(new_input)
        check(stem+'_uniform_input_recovered',all(a==[F(x) for x in old[k]] for k,a in recovered.items()))
        check(stem+'_all_entries_nonincreasing',all(all(a<=b for a,b in zip(new[k],recovered[k])) for k in recovered))
        check(stem+'_higher_entries_improve',new['U' if 'U' in new else 'common_u'][8]<recovered['U' if 'U' in new else 'common_u'][8])
        data.append(new)
    for n in [0,3,7]:
        raw=sum(F(comb(n,i))*4*factorial(i)*200**i*factorial(n-i)*100**(n-i)*(data[0]['P'][n-i]+tm[n-i]/2) for i in range(n+1))
        check('independent_raw_middle_Ward_'+str(n),raw==data[0]['P'][n+1]*factorial(n+1)*100**(n+1))
        raw=5*sum(F(comb(n,i))*factorial(i)*10**i*factorial(n-i)*100**(n-i)*tc[n-i] for i in range(n+1))
        check('independent_raw_cold_kernel_'+str(n),raw==data[1]['Ward_near'][n]*factorial(n)*100**n)
    FIXTURES.append({'fixture':'per_order_original_reconstruction','middle':{k:strings(v) for k,v in data[0].items()},
                     'cold':{k:strings(v) for k,v in data[1].items()},'outer':{k:strings(v) for k,v in data[2].items()},
                     'middle_U_decimals':decimals(data[0]['U']),'cold_U_decimals':decimals(data[1]['U']),
                     'outer_U_decimals':decimals(data[2]['common_u']),'old_low_contracts_retained':True})
    return data

def tensors(data):
    mid,cold0,out=data;A=[F(x) for x in fixture('cold-metric-ten-jets','original_finite_completion_and_tensor_conversion')['Cartesian_amplifications']]
    oa=[F(x) for x in fixture('outer-metric-ten-jets','finite_completion_full_tensors_and_forward_map')['tensor_amplifications']];one=[F(1)]*11
    prior=fixture(PREVIOUS,'low_and_tau_retaining_gap_arrays');glob={};regions={}
    for field,key in [('u','U'),('nu','N'),('V','N')]:
        mr=mid[key];cr=cold0[key]
        if field=='V':
            mr=[1200*sum(F('.2')**i*mid['N'][j-i] for i in range(j+1)) for j in range(11)]
            cr=[F('1.2')*82*sum(F('.2')**i*cold0['N'][j-i] for i in range(j+1)) for j in range(11)]
        mm=R.tensor(mr,one,6,200);cn=R.tensor(cr,A,100,20);ca=R.tensor(cr,one,15,200)
        cf=F(1) if field=='V' else F('.00331');cc=[cf*max(a,b) for a,b in zip(cn,ca)]
        oo=[9*F('.0331')*oa[j]*(max(out['common_'+field][:j+1])+F('1e-9'))+F('18e-11')*F('.1') for j in range(11)]
        glob[field]=([max(mm[j],cc[j],oo[j])*100**j for j in range(11)],
                     [max(F('50.1')**2*mm[j],F('1.001')**2*cc[j],oo[j])*100**j for j in range(11)])
        regions[field]={'middle':strings(mm),'cold':strings(cc),'outer_weighted':strings(oo)}
    arrays=[]
    for name,field,weighted in [('U','u',0),('K','u',1),('Vbar','V',0),('KV','V',1)]:
        base=[F(x) for x in prior[name]]+[F('2e26' if field=='V' else '4e21')*100**j for j in [9,10]]
        raw=glob[field][weighted]
        if field=='u':raw=[max(a,b) for a,b in zip(raw,glob['nu'][weighted])]
        arrays.append(poly([min(a,b) for a,b in zip(base,raw)]))
    U,K,V,KV=arrays
    rm=[F('1.2')*sum(F('.2')**i*mid['N'][j-i] for i in range(j+1)) for j in range(11)]
    rc=[F('1.2')*sum(F('.2')**i*cold0['N'][j-i] for i in range(j+1)) for j in range(11)]
    mm=R.tensor(rm,one,6,200);cn=R.tensor([82*x for x in rc],A,100,20);ca=R.tensor(rc,one,15,200);oo=[F(x) for x in regions['V']['outer_weighted']]
    hh=poly([min(KV[j],max(F('50.1')**2*mm[j],F('.001')*F('1.001')**2*cn[j],F('.00331')*F('1.001')**2*ca[j],oo[j])*100**j) for j in range(11)])
    for j in range(9):hh[j]=min(hh[j],F(prior['tau_KV'][j]))
    check('higher_weighted_tangent_inputs_improve',all(K[j]<F(prior['K'][j]) for j in [6,7,8]))
    check('higher_weighted_gap_inputs_improve',all(KV[j]<F(prior['KV'][j]) for j in [6,7,8]))
    check('mass_remainder_explicit',all(F(x)>9*F('.0331')*oa[j]*F('1e-9') for field in regions.values() for j,x in enumerate(field['outer_weighted'])))
    FIXTURES.append({'fixture':'full_tangent_and_tau_arrays','regions':regions,'U':strings(U[:11]),'K':strings(K[:11]),
                     'Vbar':strings(V[:11]),'KV':strings(KV[:11]),'tau_KV':strings(hh[:11]),
                     'K_decimals':decimals(K[:11]),'KV_decimals':decimals(KV[:11]),'analytic_mass_remainder_and_cubic_retained':True})
    return U,K,V,KV,hh

def old_anomaly_parts(pf,original,q,Tp,zabs):
    geo=fixture('finite-scalar-residual-six','weighted_first_and_second_path_geometry')
    M0=poly([F(x) for x in pf['metric']]);M1=poly([F(x) for x in geo['metric_first']]);C0=poly([F(x) for x in pf['connection']]);C1=poly([F(x) for x in geo['connection_first']])
    zw=poly([F(x) for x in pf['weighted_Z_initial']]);jw=poly([F(x) for x in pf['weighted_J_initial']]);cw=poly([F(x) for x in pf['weighted_C_initial']])
    e=F(1,2000**10*180)
    def linear(dz):
        bx=scale('6.1',add(mul(M0,add(der(der(dz)),mul(C0,der(dz)))),mul(M1,der(der(zabs))),mul(add(mul(M1,C0),mul(M0,C1)),der(zabs))))
        return scale(e,add(bx,scale(F(11,3),dz),scale(F(11,72),mul(add(zw,zabs),dz))))
    seedB=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    ref=add(compose(der(der(seedB)),q),poly(['2e-10']))
    rr0=add(ref,scale(2,der(C0)),scale(2,power(C0,2)));rr1=add(scale(2,der(C1)),scale(4,mul(C0,C1)))
    ri0=scale(2,rr0);ri1=scale(2,rr1);sc0=mul(M0,ri0);sc1=add(mul(M1,ri0),mul(M0,ri1))
    j1=scale('6.1',add(ri1,scale(F(1,4),add(mul(sc1,M0),mul(sc0,M1)))))
    w1=scale('6.1',add(mul(M1,rr0),mul(M0,rr1),scale(2,add(mul(M1,ri0),mul(M0,ri1))),scale(F(2,3),add(mul(sc1,power(M0,2)),scale(2,mul(mul(sc0,M0),M1))))))
    dc=add(scale(2,mul(power(M0,4),mul(cw,w1))),scale(4,mul(mul(power(M0,3),M1),power(cw,2))))
    dj=add(scale(2,mul(power(M0,2),mul(jw,j1))),scale(2,mul(mul(M0,M1),power(jw,2))))
    mu=scale(F('4e-27')*F('1.001')/18,mul(exps(scale(4,Tp)),power(q,4)))
    dmu=scale(4,mul(mu,poly([F(x) for x in geo['S1']])))
    rest=add(scale(e,add(scale(F(7,12),dc),scale(F(11,6),dj))),dmu)
    olddz=add(poly([F(x) for x in original['linear_trace_W2']]),scale(F(1,2),poly([F(x) for x in original['curvature_remainder_second_W2']])))
    check('complete_old_anomaly_reconciled',add(linear(olddz),rest)[:7]==[F(x) for x in original['anomaly_difference_W2']])
    check('curvature_and_magnetic_contacts_retained',min(rest[:7]+dc[:7]+dj[:7]+dmu[:7])>0)
    return linear,rest,C0

def scalar_ledger(U,K,V,KV,hh):
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    original=fixture('finite-scalar-residual-six','actual_scalar_residual_H6_ledger');prior=fixture(PREVIOUS,'cold_product_scalar_H6_ledger')
    q=add(poly([1]),poly([F(x) for x in pf['seed_jets']]));ll=poly([F(x) for x in pf['log_radius_jets']])
    Tp=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    Up=U.copy();Up[0]=0;Vp=V.copy();Vp[0]=0;O=scale('1.001',exps(scale(4,Tp)));E=exps(add(Vp,Tp))
    S1=scale('1.001',mul(E,K));S1u=scale('1.001',mul(E,U));S2=scale('1.001',mul(mul(q,exps(add(Vp,scale(2,Tp)))),mul(K,KV)))
    W1=add(K,S1);J2=mul(O,add(mul(add(scale(4,power(W1,2)),scale(2,S2)),add(poly([1]),scale(4,U))),scale(16,mul(K,W1))))
    seedB=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')]);zabs=compose(za,q)
    zder=[j*za[j] for j in range(len(za))];lie=mul(K,compose(zder,q))
    seed=[F(x) for x in original['complete_seed_polynomial']];Zd=add(compose(seed,q),lie)
    check('full_degree_twelve_seed_retained',len(seed)==13 and seed[-1]>0)
    LN=scale('6.1',mul(power(q,2),add(scale(2,der(der(K))),scale(2,mul(der(ll),der(K))))));LN[0]=min(LN[0],F('.005'))
    LT=scale('6.1',add(scale(2,mul(q,compose(seedB,q))),mul(power(q,2),compose(der(seedB),q))))
    H=mul(KV,LN)
    for n in range(9):H[n]-=KV[n]*LN[0]-min(KV[n]*LN[0],F('.005')*hh[n])
    G=B.gap_product(KV,hh,q,order=8);plain=mul(q,power(der(KV),2));G=poly([min(a,b) for a,b in zip(G,plain)])
    QQ=scale('.501',mul(E,power(KV,2)))
    terms=[scale(F(1,2),mul(add(poly([12]),zabs),J2)),scale(2,mul(mul(O,add(U,S1u)),Zd)),
           scale(6*F('1.001'),mul(mul(mul(O,q),E),H)),scale(6,mul(mul(O,QQ),LT)),
           scale(6*F('6.1')*F('1.001'),mul(mul(mul(O,E),power(q,2)),G)),
           scale(F('4.004')*F('1.3e-10'),mul(mul(power(q,2),exps(add(scale(4,Tp),scale(2,Up)))),power(K,2)))]
    new=add(*terms);nr=poly([min(new[j],F(original['curvature_remainder_second_W2'][j])/2) for j in range(9)])
    for j in range(7):nr[j]=min(nr[j],F(prior['selected_curvature_remainder'][j]))
    linear,rest,C0=old_anomaly_parts(pf,original,q,Tp,zabs);dh=add(linear(add(Zd,nr)),rest)
    ra=poly()
    for n,a in [(5,4),(6,6)]:ra=add(ra,scale(F('1.001')*F('1.3e-22')*a,mul(mul(power(q,n),exps(scale(n,Tp))),add(scale(n*n,power(S1,2)),scale(n,S2)))))
    p=add(nr,scale(F(1,2),ra),lie,dh)
    p=poly([min(p[j],F(prior['residual_reference_W2'][j])) for j in range(7)])
    phys=L.covariant(p,C0,0,6);constant=F(1001,24000)*sum(phys);old=F(prior['C_tau_exact']);exponent=L.decimal_power(constant)
    check('curvature_remainder_eight_available',len(nr[:9])==9 and min(nr[:9])>0 and 8+2==10)
    check('new_anomaly_improves',dh[6]<F(original['anomaly_difference_W2'][6]),{'old':L.decimal(F(original['anomaly_difference_W2'][6])),'new':L.decimal(dh[6])})
    check('complete_scalar_bound_improves',constant<old,{'constant':L.decimal(constant),'improvement':L.decimal(old/constant)})
    check('upward_exact_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent)
    budget=F(10**57)/(800*F(10)**exponent);check('rounded_transfer_budget',800*F(10)**exponent*budget==10**57)
    FIXTURES.append({'fixture':'per_order_scalar_and_anomaly_ledger','remainder_terms':[strings(t[:9]) for t in terms],
                     'remainder_terms_decimals':[decimals(t[:9]) for t in terms],'selected_curvature_remainder':strings(nr[:9]),
                     'anomaly_difference':strings(dh[:7]),'anomaly_difference_decimals':decimals(dh[:7]),'retained_anomaly_rest':strings(rest[:7]),
                     'residual_reference_W2':strings(p[:7]),'physical_covariant_bounds':strings(phys),'physical_covariant_decimals':decimals(phys),
                     'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),'C_tau_power_of_ten':exponent,
                     'old_C_tau_exact':str(old),'improvement_ratio':L.decimal(old/constant),
                     'scalar_only_rounded_transfer_budget_over_g_epsilon_c':str(budget),'scalar_only_budget_decimal':L.decimal(budget),
                     'full_noncompact_physical_volume_retained':True,'full_Schur_source_bounded':False,'target_passed':False})

def independent_eighth_remainders():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for t0 in [mp.mpf('.2'),mp.mpf('.8')]:
                for ee in [mp.mpf('.3'),mp.mpf('.9')]:
                    def quantities(x):
                        ta=t0+x/100;z=1-ta;phi=-mp.log(z);pp=1/(100*z);ppp=1/(10000*z*z)
                        u=(1+2*x+x*x+x**4)/1000;upp=(2+12*x*x)/1000
                        v=(2-x+3*x*x+x**3)/1000;vp=(-1+6*x+3*x*x)/1000;vpp=(6+6*x)/1000
                        nu=ta*v;nup=v/100+ta*vp;nupp=2*vp/100+ta*vpp
                        k=mp.exp(ee*v);T=1+ta*(k-1);Tx=(k-1)/100+ta*k*ee*vp
                        Txx=2*k*ee*vp/100+ta*k*(ee*vpp+ee**2*vp**2)
                        ss=mp.log(T);sp=Tx/T;spp=Txx/T-sp*sp;w=ee*u+ss
                        base=z*z*(-6*ppp-6*pp*pp+2);LN=z*z*(nupp+2*pp*nup);LT=z*z*2*pp/100
                        dot=-2*(u+nu)*base-2*z*z*upp-6*LN+4*u*z*z
                        direct=mp.exp(-2*(phi+w))*(-2*(ppp+ee*upp+spp)-4*(ppp+spp)-6*(pp+sp)**2)+2*mp.exp(-2*(phi+ss))
                        bb=k-1-ee*v*k
                        rem=base*(mp.exp(-2*w)*(1+2*ee*(u+nu))-1)+ee*mp.expm1(-2*w)*dot
                        rem-=6*mp.exp(-2*w)*(ee*(k/T-1)*LN+bb/T*LT+ee**2*k/T*ta*z*z*vp*vp)
                        rem+=2*mp.exp(-2*w)*(mp.expm1(2*ee*u)-2*ee*u)*z*z
                        return direct-base-ee*dot,rem
                    radius=mp.mpf('.05');count=64
                    contour=sum(quantities(radius*mp.exp(2j*mp.pi*j/count))[0]*mp.exp(-16j*mp.pi*j/count) for j in range(count))/(count*radius**8)
                    real=mp.diff(lambda x:quantities(x)[1],mp.mpf(0),8)/mp.factorial(8);error=abs(contour-real)
                    check('independent_eighth_'+str(t0)+'_'+str(ee)+'_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    records.append({'tau0':str(t0),'ell':str(ee),'dps':dps,'order':8,'error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_eighth_finite_remainders','records':records,'actual_charged_source_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 per-order source scalar refinement','input_count':len(ps)}),flush=True)
    try:
        sources=complete_sources();data=reconstruction(*sources);args=tensors(data);scalar_ledger(*args);independent_eighth_remainders()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 per-order source scalar refinement','created_utc':datetime.now(timezone.utc).isoformat(),
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
