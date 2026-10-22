"""Refreshed finite local coefficients with explicit metric order eleven."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-refreshed-local-mixed';PREVIOUS='finite-per-order-mixed-row'
OWN=['checks/finite_refreshed_local_mixed.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_mixed_row',ROOT/'checks/finite_per_order_mixed_row.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);L=B.L
poly,add,mul,der,scale,power,compose,exps=[getattr(L,k) for k in ['poly','add','mul','der','scale','power','compose','exps']]
ROWS=[];FIXTURES=[]
PROTOCOL={'finite_metric_order':11,'parameter_tangent_order':10,'component_curvature_order':9,
          'integrated_curvature_cap_order':8,'full_trace_order':7,'local_coefficient_order':6,
          'mixed_output_order':6,'TT_input_order':8,'numerical_dps':[60,90],
          'manufactured_finite_field_order':10,'actual_charged_source_sampled':False,
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
def arr(d,k):return poly([F(x) for x in d[k]])
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})

def inputs():
    old=saved(PREVIOUS)
    check('previous_complete_mixed_row_passed',old['status']=='passed' and len(old['checks'])==73 and all(x['passed'] for x in old['checks']))
    check('all_previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    for stem in ['finite-low-gap-reconstruction','finite-radius-operator-source','finite-higher-local-coefficients','finite-higher-tt-lift']:
        r=saved(stem)
        check(stem+'_retained_pass',r['status']=='passed' and all(x['passed'] for x in r['checks']))
    check('previous_claim_boundaries_retained',not old['actual_transfer_certificate_passed'] and not old['higher_TT_inverse_quantified'] and not old['actual_TT_source_H4_bounded'])
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    z=arr(pf,'seed_jets');ll=arr(pf,'log_radius_jets');q=add(poly([1]),z)
    return pf,q,ll

def geometry(t,k,q,ll):
    m=scale('3.003',exps(add(scale(2,ll),scale(4,t))));dm=scale(4,mul(m,k))
    c=scale(F(3,2),mul(m,der(m)));dc=scale(F(3,2),add(mul(dm,der(m)),mul(m,der(dm))))
    seed=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    ref=add(compose(der(der(seed)),q),poly(['2e-10']))
    rr=add(ref,scale(2,der(c)),scale(2,mul(c,c)));dr=add(scale(2,der(dc)),scale(4,mul(c,dc)))
    ri=scale(2,rr);dri=scale(2,dr);sc=mul(m,ri);ds=add(mul(dm,ri),mul(m,dri))
    dj=add(dri,scale(F(1,4),add(mul(ds,m),mul(sc,dm))))
    dw=add(mul(dm,rr),mul(m,dr),scale(2,add(mul(dm,ri),mul(m,dri))),
           scale(F(2,3),add(mul(ds,power(m,2)),scale(2,mul(mul(sc,m),dm)))))
    ip=poly([1]);z=q.copy();z[0]=0
    for j in range(1,12):ip=add(ip,power(z,j))
    pj=add(scale('.501',power(q,2)),scale('5e-24',power(q,3)),scale('4.5e-24',power(q,4)))
    pc=add(scale(3*F('.668'),ip),poly([6*F('.501')]),scale('20e-24',q),scale('15e-24',power(q,2)))
    pz=scale('1.3e-22',add(scale(4,power(q,5)),scale(6,power(q,6))))
    return {'M':m,'DeltaM':dm,'C':c,'DeltaC':dc,'JW':add(scale('6.1',dj),scale(13,pj)),
            'CW':add(scale('6.1',dw),scale(25,pc)),'ZW':add(scale('6.1',ds),pz),
            'Rm':scale('6.1',mul(m,rr)),'mu':scale(F('4e-27')*F('1.001')/18,mul(exps(scale(4,t)),power(q,4))),
            'pZ':pz}

def finite_fields(pf,q,ll):
    ot=poly([0]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)]+[F('4e24')*100**11])
    ok=ot.copy();ok[0]=F('1e-4')
    for j in range(1,6):ok[j]=10*100**j
    restored=geometry(ot,ok,q,ll)
    for name,key,count in [('M','metric',12),('DeltaM','weighted_metric_difference',12),('C','connection',11),
                           ('DeltaC','weighted_connection_difference',11),('JW','weighted_J_initial',10),
                           ('CW','weighted_C_initial',10),('ZW','weighted_Z_initial',10)]:
        check('old_geometry_exact_'+name,restored[name][:count]==[F(x) for x in pf[key][:count]])
    tf=fixture('finite-low-gap-reconstruction','complete_improved_tangent_arrays')
    ff=fixture('finite-low-gap-reconstruction','finite_parameter_integral_bounds')
    check('only_ten_parameter_derivatives_available',all(len(tf[k])==11 for k in ['U','K','Vbar']) and len(ff['Tplus'])==11)
    t=arr(ff,'Tplus');k=arr(tf,'K');v=arr(tf,'Vbar');v[0]=0
    # Discard the padded eleventh coefficient before it can enter a field bound.
    s1=scale('1.001',mul(exps(add(v,t)),k))[:11]
    kf=[min(ok[j],max(k[j],s1[j])) for j in range(11)]+[ok[11]];t[11]=ot[11]
    gf=fixture('finite-low-gap-reconstruction','refined_finite_geometry_and_components')
    check('weighted_parameter_array_recovered',s1[:9]==[F(x) for x in gf['geometry']['S1']])
    check('finite_field_bounds_nonincreasing',all(t[j]<=ot[j] and kf[j]<=ok[j] for j in range(12)))
    check('finite_order_eleven_explicitly_retained',t[11]==kf[11]==F('4e24')*100**11 and len(s1)==11)
    check('finite_low_fields_strictly_improved',all(kf[j]<ok[j] for j in range(11)))
    FIXTURES.append({'fixture':'integrated_finite_field_bounds','Tplus':strings(t),'Kfield':strings(kf),
                     'Kfield_decimals':decimals(kf),'weighted_parameter_through_ten':strings(s1),
                     'finite_metric_order':11,'parameter_tangent_order':10,'unknown_parameter_order_eleven_used':False})
    return t,kf,restored,gf,k

def refreshed_components(pf,q,ll,t,kf,old,gf,k):
    g=geometry(t,kf,q,ll)
    check('new_metric_recovered_at_available_orders',g['M'][:9]==[F(x) for x in gf['geometry']['M0']])
    check('new_connection_recovered_at_available_orders',g['C'][:9]==[F(x) for x in gf['geometry']['C0']])
    for name,order in [('M',11),('DeltaM',11),('C',10),('DeltaC',10),('JW',9),('CW',9),('ZW',9),('Rm',7),('mu',7)]:
        check('uncapped_geometry_nonincreasing_'+name,all(0<=g[name][j]<=old[name][j] for j in range(order+1)))
    chosen={name:g[name].copy() for name in ['JW','CW','ZW']}
    for name in chosen:
        for j in range(9):chosen[name][j]=min(chosen[name][j],F(gf[name][j]))
    sf=fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly')
    seed=[F(x) for x in fixture('finite-scalar-residual-six','actual_scalar_residual_H6_ledger')['complete_seed_polynomial']]
    check('degree_twelve_seed_polynomial_retained',len(seed)==13 and seed[-1]>0)
    za=poly([0]*5+[4*F('1.3e-22'),6*F('1.3e-22')])
    lie=mul(k,compose([j*za[j] for j in range(len(za))],q))
    check('radius_source_Lie_recovered',lie[:7]==[F(x) for x in sf['Lie']])
    dz=add(compose(seed,q),lie,arr(sf,'selected_curvature_remainder'))
    for j in range(9):chosen['ZW'][j]=min(chosen['ZW'][j],g['pZ'][j]+dz[j])
    check('caps_only_through_eight',all(chosen[name][9]==g[name][9] for name in chosen))
    check('scalar_residual_contract_unchanged',sf['C_tau_power_of_ten']==20 and F(sf['C_tau_exact'])<10**20)
    tt=t.copy();kk=kf.copy();tt[11]=kk[11]=0;omitted=geometry(tt,kk,q,ll)
    full_z9=L.covariant(chosen['ZW'],g['C'],0,9)[9]
    omitted_zw=chosen['ZW'].copy();omitted_zw[9]=omitted['ZW'][9]
    omitted_z9=L.covariant(omitted_zw,g['C'],0,9)[9]
    check('order_eleven_omission_changes_ninth_scalar_derivative',full_z9>omitted_z9>0)
    check('order_eleven_omission_leaves_lower_components',all(omitted[name][:9]==g[name][:9] for name in chosen))
    FIXTURES.append({'fixture':'refreshed_finite_component_geometry',
                     'geometry':{name:strings(g[name][:order+1]) for name,order in [('M',11),('DeltaM',11),('C',10),('DeltaC',10),('Rm',7),('mu',7)]},
                     'uncapped_curvature':{name:strings(g[name][:10]) for name in chosen},
                     'weighted_components':{name:strings(v[:10]) for name,v in chosen.items()},
                     'scalar_difference_DZ':strings(dz[:9]),'cap_order':8,
                     'omitted_order_eleven_Z9_drop':str(full_z9-omitted_z9),
                     'omitted_order_eleven_E7_drop':str(2*F(1,2000**10*180)*(full_z9-omitted_z9)),
                     'omitted_array_accepted_as_bound':False,'scalar_C_tau_exact':sf['C_tau_exact']})
    return g,chosen

def full_mixed(g,chosen):
    L.check=check;L.FIXTURES=FIXTURES
    L.actual_coefficients(g['C'],chosen['JW'],chosen['CW'],chosen['ZW'],g['Rm'],g['mu'])
    cf=FIXTURES[-1];cf['fixture']='refreshed_covariant_curvature_and_local_coefficients'
    cf['higher_TT_lift_supplied']=True
    old=fixture('finite-higher-local-coefficients','actual_covariant_curvature_and_local_coefficients')
    for name in ['weighted_J','weighted_C','weighted_Z','unweighted_Rm','weighted_magnetic','trace_coefficients','local_A0','local_A1','local_A2']:
        check('actual_coefficients_nonincreasing_'+name,all(F(a)<=F(b) for a,b in zip(cf[name],old[name])) and len(cf[name])==len(old[name]))
    lf=fixture('finite-higher-tt-lift','per_input_order_lift_and_return')
    aa=[[F(x) for x in cf['local_A'+str(j)]] for j in range(3)];ee=[F(x) for x in cf['trace_coefficients']]
    hh=[[F(x) for x in row] for row in lf['tensor']];yy=[[F(x) for x in row] for row in lf['Y'][:8]]
    B.check=check;B.FIXTURES=FIXTURES
    columns=B.matrix_ledger(aa,ee,hh,yy,lf);mf=FIXTURES[-1]
    previous=fixture(PREVIOUS,'complete_per_order_mixed_matrix');before=F(previous['Cmix_exact'])
    check('every_complete_mixed_column_nonincreasing',all(v<=F(w) for v,w in zip(columns,previous['column_sums'])))
    check('complete_mixed_coefficient_strictly_improves',max(columns)<before)
    mf.update(previous_rounded_coefficient='1e55',previous_exact_coefficient=str(before),
              improvement_over_previous_rounded=L.decimal(F(10)**55/max(columns)),
              improvement_over_previous_exact=L.decimal(before/max(columns)))
    B.inverse_composition(columns);iv=FIXTURES[-1]
    oldiv=fixture(PREVIOUS,'low_TT_inverse_composition')
    for name in ['Cstart','Klow_bar','mhi']:
        check('complete_inverse_'+name+'_improves',F(iv[name+'_exact'])<F(oldiv[name+'_exact']))
        iv[name+'_previous_exact']=oldiv[name+'_exact'];iv[name+'_improvement']=L.decimal(F(oldiv[name+'_exact'])/F(iv[name+'_exact']))
    entries={}
    for l in [0,2]:
        terms=[]
        for n in range(7):
            for j in range(3):
                for k in range(n+1):
                    terms.append(('A'+str(j),n,k,F(comb(n,k),24)*aa[j][k]*hh[n-k+j][l]))
            for k in range(n+1):
                terms.append(('trace_gradient',n,k+1,F(comb(n,k),24)*ee[k+1]*yy[n-k][l]))
                terms.append(('trace_divergence',n,k,F(comb(n,k),12)*ee[k]*yy[n-k+1][l]))
        check('diagnostic_exact_column_'+str(l),sum(v[3] for v in terms)==columns[l])
        entries[str(l)]=[{'channel':a,'output_order':n,'coefficient_order':k,'value':str(v),'decimal':L.decimal(v),
                          'fraction':L.decimal(v/columns[l])} for a,n,k,v in sorted(terms,key=lambda x:x[3],reverse=True)[:8]]
    FIXTURES.append({'fixture':'complete_budget_diagnostics','largest_column_terms':entries,
                     'low_inverse_part_fractions':[L.decimal(F(v)/F(iv['Klow_bar_exact'])) for v in iv['low_coefficient_parts']],
                     'all_higher_source_and_inverse_unknowns_retained':True})

def independent_scalar_segment():
    f=poly([0]+[F((-1)**j*(j+1),10**(j+2)) for j in range(1,12)])
    h=poly([0]+[F(j+2,20**(j+2)) for j in range(1,12)]);r=add(f,scale(-1,h))
    def series(a):return add(*[scale(F(1,factorial(j)),power(a,j)) for j in range(12)])
    direct=add(series(f),scale(-1,series(h)));integral=poly()
    for n in range(11):
        for k in range(n+1):integral=add(integral,scale(F(comb(n,k),factorial(n)*(k+1)),mul(power(h,n-k),power(r,k))))
    integrated=mul(r,integral)
    check('independent_exponential_series_through_eleven',series(f)==exps(f) and series(h)==exps(h))
    check('independent_integrated_scalar_segment_through_eleven',direct==integrated)
    omitted=f.copy();omitted[11]=0
    check('independent_eleventh_spatial_coefficient_matters',series(omitted)[11]!=series(f)[11] and series(omitted)[:11]==series(f)[:11])
    FIXTURES.append({'fixture':'independent_finite_scalar_segment','order':11,'direct':strings(direct),
                     'integrated':strings(integrated),'parameter_integral_exact':True})

def independent_weighted_fields():
    x=S.symbols('x');vv=(2-x+3*x*x+x**3)/10**5;uu=(1+x*x)/10**5;records=[]
    for tau0 in [F(1,5),F(4,5)]:
        tau=S.Rational(tau0.numerator,tau0.denominator)+x/10;nu=S.expand(tau*vv)
        vc=poly([abs(F(str(S.expand(vv).coeff(x,j)))) for j in range(11)])
        uc=poly([max(abs(F(str(S.expand(uu).coeff(x,j)))),abs(F(str(nu.coeff(x,j))))) for j in range(11)])
        vp=vc.copy();vp[0]=0;t=poly([0]+[F('.001')*2**j for j in range(1,11)])
        for _ in range(6):
            su=scale('1.001',mul(exps(add(vp,t)),uc));t=poly([0]+[min(t[j],max(uc[j],su[j])) for j in range(1,11)])
        weight=F(26,5);kw=scale(weight,uc);s1=scale('1.001',mul(exps(add(vp,t)),kw))[:11]
        finite=[min(weight*F('.001')*2**j,max(kw[j],s1[j])) for j in range(11)]
        for dps in [60,90]:
            with mp.workdps(dps):
                w=mp.mpf(weight.numerator)/weight.denominator;t0=mp.mpf(tau0.numerator)/tau0.denominator
                for a in [mp.mpf('.3'),mp.mpf('.9')]:
                    fn=lambda z:mp.log1p((t0+z/10)*mp.expm1(a*(2-z+3*z*z+z**3)/10**5))
                    vals=mp.taylor(fn,0,10);uvals=mp.taylor(lambda z:a*(1+z*z)/10**5,0,10)
                    ratios=[w*max(abs(vals[j]),abs(uvals[j]))/(mp.mpf(finite[j].numerator)/finite[j].denominator) for j in range(11)]
                    check('independent_weighted_finite_field_'+str(tau0)+'_'+str(dps)+'_'+str(a),max(ratios)<1,mp.nstr(max(ratios),18))
                    records.append({'tau0':str(tau0),'a':str(a),'dps':dps,'weight':str(weight),'order':10,
                                    'maximum_ratio':mp.nstr(max(ratios),18),'ratios':[mp.nstr(v,18) for v in ratios]})
    FIXTURES.append({'fixture':'independent_integrated_weighted_finite_fields','records':records,
                     'actual_charged_geometry_sampled':False,'unknown_eleventh_tangent_used':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 refreshed local mixed coefficients','input_count':len(ps)}),flush=True)
    try:
        pf,q,ll=inputs();t,kf,old,gf,k=finite_fields(pf,q,ll)
        g,chosen=refreshed_components(pf,q,ll,t,kf,old,gf,k);full_mixed(g,chosen)
        independent_scalar_segment();independent_weighted_fields()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 refreshed local mixed coefficients','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for key,value in review().items():
        if key not in run and isinstance(value,bool):run[key]=value
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
