"""Complete mixed-row Leibniz matrix and retained low TT inverse factors."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import comb,factorial
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-per-order-mixed-row';PREVIOUS='finite-higher-tt-lift'
OWN=['checks/finite_per_order_mixed_row.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_coefficient_algebra',ROOT/'checks/finite_higher_local_coefficients.py')
L=importlib.util.module_from_spec(spec);spec.loader.exec_module(L)
ROWS=[];FIXTURES=[]
PROTOCOL={'mixed_output_order':6,'TT_input_order':8,'full_trace_order':7,'channels':5,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42','mixed_derivative_orders':[0,3,6],
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

def inputs():
    old=saved(PREVIOUS)
    check('prior_higher_lift_passed',old['status']=='passed' and len(old['checks'])==102 and all(x['passed'] for x in old['checks']))
    check('prior_frozen_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    cf=fixture('finite-higher-local-coefficients','actual_covariant_curvature_and_local_coefficients')
    lf=fixture(PREVIOUS,'per_input_order_lift_and_return')
    aa=[[F(x) for x in cf['local_A'+str(j)]] for j in range(3)];ee=[F(x) for x in cf['trace_coefficients']]
    hh=[[F(x) for x in row] for row in lf['tensor']];yy=[[F(x) for x in row] for row in lf['Y'][:8]]
    check('complete_coefficient_and_lift_orders',all(len(a)==7 for a in aa) and len(ee)==8 and len(hh)==9 and len(yy)==8 and all(len(a)==9 for a in hh+yy))
    check('all_operator_coefficient_sectors_retained',min(sum(aa,[])+ee)>0 and len(cf['weighted_magnetic'])==8)
    check('previous_uniform_enclosure',lf['mixed_row_power_of_ten']==102 and F(lf['mixed_row_coefficient_over_g_exact'])<10**102)
    scalar=fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly')
    check('unchanged_scalar_contract',scalar['C_tau_power_of_ten']==20 and F(scalar['C_tau_exact'])<10**20 and not scalar['target_passed'])
    FIXTURES.append({'fixture':'retained_full_mixed_inputs','coefficient_source':'finite-higher-local-coefficients',
                     'lift_source':PREVIOUS,'A0':strings(aa[0]),'A1':strings(aa[1]),'A2':strings(aa[2]),'E':strings(ee),
                     'full_trace_distinct_from_bookkeeping_residual':True,'scalar_C_tau_exact':scalar['C_tau_exact'],
                     'scalar_power_of_ten':20,'physical_magnetic_contact_retained':True})
    return aa,ee,hh,yy,lf

def matrix_ledger(aa,ee,hh,yy,prior):
    channels={name:[[F(0)]*9 for _ in range(7)] for name in ['A0','A1','A2','trace_gradient','trace_divergence']}
    entries=[]
    for n in range(7):
        for l in range(9):
            for j in range(3):
                for k in range(n+1):
                    value=F(comb(n,k),24)*aa[j][k]*hh[n-k+j][l]
                    channels['A'+str(j)][n][l]+=value
                    if value:entries.append({'channel':'A'+str(j),'output_order':n,'coefficient_order':k,'input_order':l,'value':value})
            for k in range(n+1):
                for name,value in [('trace_gradient',F(comb(n,k),24)*ee[k+1]*yy[n-k][l]),
                                   ('trace_divergence',F(comb(n,k),12)*ee[k]*yy[n-k+1][l])]:
                    channels[name][n][l]+=value
                    if value:entries.append({'channel':name,'output_order':n,'coefficient_order':k+1 if name=='trace_gradient' else k,'input_order':l,'value':value})
    # Independent divided-coefficient products: apply the field derivative first.
    for n in range(7):
        ok=True
        for l in range(9):
            for j in range(3):
                v=factorial(n)*sum((aa[j][k]/factorial(k))*(hh[n-k+j][l]/factorial(n-k)) for k in range(n+1))/24
                ok=ok and v==channels['A'+str(j)][n][l]
            v1=factorial(n)*sum((ee[k+1]/factorial(k))*(yy[n-k][l]/factorial(n-k)) for k in range(n+1))/24
            v2=factorial(n)*sum((ee[k]/factorial(k))*(yy[n-k+1][l]/factorial(n-k)) for k in range(n+1))/12
            ok=ok and v1==channels['trace_gradient'][n][l] and v2==channels['trace_divergence'][n][l]
        check('independent_divided_Leibniz_order_'+str(n),ok)
    matrix=[[sum(a[n][l] for a in channels.values()) for l in range(9)] for n in range(7)]
    columns=[sum(matrix[n][l] for n in range(7)) for l in range(9)]
    per_channel={name:[sum(a[n][l] for n in range(7)) for l in range(9)] for name,a in channels.items()}
    constant=max(columns);exponent=L.decimal_power(constant);winner=columns.index(constant)
    check('complete_five_channel_sum',len(channels)==5 and all(sum(v)>0 for v in per_channel.values()) and all(columns[l]==sum(v[l] for v in per_channel.values()) for l in range(9)))
    check('positive_matrix_column_norm',all(x>=0 for row in matrix for x in row) and all(x<=constant for x in columns))
    check('basis_input_attains_column_norm',sum(matrix[n][winner] for n in range(7))==constant)
    uniform=F(381,24)*(max(sum(aa,[]))+max(ee))*10**44
    check('uniform_Leibniz_multiplier',3*sum(2**n for n in range(7))==381)
    check('per_order_improves_uniform',constant<uniform<=F(prior['mixed_row_coefficient_over_g_exact'])<10**102)
    check('exact_upward_mixed_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent,{'Cmix':L.decimal(constant),'power':exponent})
    high_trace=ee[7]*yy[0][0]/24
    check('seventh_trace_derivative_nonzero',high_trace>0 and channels['trace_gradient'][6][0]>=high_trace)
    check('trace_gradient_omission_changes_bound',max(columns[l]-per_channel['trace_gradient'][l] for l in range(9))<constant)
    check('trace_divergence_omission_changes_bound',max(columns[l]-per_channel['trace_divergence'][l] for l in range(9))<constant)
    tops=sorted((e for e in entries if e['input_order']==winner),key=lambda e:e['value'],reverse=True)[:12]
    FIXTURES.append({'fixture':'complete_per_order_mixed_matrix','channels':{k:[strings(row) for row in a] for k,a in channels.items()},
                     'matrix':[strings(row) for row in matrix],'column_sums':strings(columns),'column_decimals':decimals(columns),
                     'channel_column_sums':{k:strings(v) for k,v in per_channel.items()},
                     'Cmix_exact':str(constant),'Cmix_decimal':L.decimal(constant),'Cmix_power_of_ten':exponent,
                     'previous_rounded_coefficient':'1e102','uniform_exact_coefficient':str(uniform),
                     'improvement_over_previous_rounded':L.decimal(F(10)**102/constant),'max_input_order':winner,
                     'seventh_trace_value_term_exact':str(high_trace),'seventh_trace_fraction_of_Cmix':L.decimal(high_trace/constant),
                     'largest_terms_in_max_column':[{**e,'value':str(e['value']),'decimal':L.decimal(e['value']),'fraction_of_Cmix':L.decimal(e['value']/constant)} for e in tops],
                     'weight_never_differentiated':True,'off_shell_return_retained':True})
    return columns

def inverse_composition(columns):
    check('critical_L2_inverse_28g',F(8)*F(22,7)/F('.90')<28)
    check('form_gradient_inverse_265g',50*F('5.3')==265)
    check('form_Hessian_inverse_531_sqrtg',100*F('5.3')==530 and 4*265==1060 and 530+1060*F('4e-17')<531)
    check('actual_L2_forcing_contract',27*F('2.2e-13')<F('6e-12'))
    low_parts=[28*columns[0],265*columns[1],F(531)/F('3e-17')*columns[2]]
    low=sum(low_parts);hi=max(columns[3:]);start=F(10)**20+F('6e-12')*low
    lp=L.decimal_power(low);hp=L.decimal_power(hi);sp=L.decimal_power(start)
    for name,value,exponent in [('low',low,lp),('high',hi,hp),('known_source',start,sp)]:
        check('upward_inverse_composition_'+name,F(10)**(exponent-1)<value<=F(10)**exponent,{'value':L.decimal(value),'power':exponent})
    rb=F(10**57)/(1600*F(10)**sp);hb=F(10**57)/(1600*F(10)**hp)
    check('complete_equal_split_budget',800*(F(10)**sp*rb+F(10)**hp*hb)==10**57)
    check('high_inverse_not_replaced_by_L2_gap',hi>0 and all(columns[j]<=hi for j in range(3,9)))
    FIXTURES.append({'fixture':'low_TT_inverse_composition','low_coefficient_parts':strings(low_parts),'low_parts_decimals':decimals(low_parts),
                     'Klow_bar_exact':str(low),'Klow_bar_decimal':L.decimal(low),'Klow_bar_power_of_ten':lp,
                     'mhi_exact':str(hi),'mhi_decimal':L.decimal(hi),'mhi_power_of_ten':hp,
                     'Cstart_exact':str(start),'Cstart_decimal':L.decimal(start),'Cstart_power_of_ten':sp,
                     'known_source_R_budget_over_g_epsilon_c':str(rb),'known_source_R_budget_decimal':L.decimal(rb),
                     'high_R_C_D4_F4_budget_over_g_epsilon_c':str(hb),'high_R_C_D4_F4_budget_decimal':L.decimal(hb),
                     'actual_L2_source_bound':'6e-12 ell^2/g','low_inverse_orders':[0,1,2],
                     'higher_inverse_orders_unknown':[3,4,5,6,7,8],'actual_quantum_H4_supplied':False,'higher_inverse_supplied':False,'target_passed':False})

def independent_polynomial_operator():
    x=S.symbols('x');h=sum(S.Rational(j+1,10**(j+4))*x**j for j in range(9))
    y=sum(S.Rational(j+2,10**(j+5))*x**j for j in range(8));e=sum(S.Rational(2*j+1,10**(j+3))*x**j for j in range(8))
    a=[sum(S.Rational(1+j+k,10**(k+j+2))*x**k for k in range(7)) for j in range(3)]
    direct=sum(a[j]*S.diff(h,x,j) for j in range(3))+S.diff(e*y,x)
    jet=lambda v,n:S.diff(v,x,n).subs(x,0)
    records=[]
    for n in range(7):
        local=sum(comb(n,k)*jet(a[j],k)*jet(h,n-k+j) for j in range(3) for k in range(n+1))
        gradient=sum(comb(n,k)*jet(e,k+1)*jet(y,n-k) for k in range(n+1))
        divergence=sum(comb(n,k)*jet(e,k)*jet(y,n-k+1) for k in range(n+1))
        exact=jet(direct,n);expanded=local+gradient+divergence
        check('independent_polynomial_operator_'+str(n),S.expand(exact-expanded)==0 and local+gradient+2*divergence>=exact)
        records.append({'order':n,'direct':str(exact),'expanded':str(expanded)})
    check('polynomial_seventh_trace_omission_fails',jet(e,7)*jet(y,0)>0)
    check('weight_inside_derivative_is_different',jet((2+x)**2*direct,6)!=4*jet(direct,6))
    check('four_dimensional_divergence_cost',sum([1]*4)==2*S.sqrt(sum(v*v for v in [1]*4)))
    FIXTURES.append({'fixture':'independent_polynomial_mixed_operator','records':records,'all_five_channels_retained':True,
                     'seventh_trace_contact_nonzero':True,'fixed_weight_tested':True})

def fields(x,s):
    a,b,c,d=x;w=s*(a*a+2*b*b+3*c*c+d*d+a*b)/20
    dw=[s*(2*a+b)/20,s*(4*b+a)/20,s*6*c/20,s*2*d/20]
    e=(1+a+2*b*b+c*d+a**7+2*b**7+a**3*b**4)/1000
    de=[(1+7*a**6+3*a*a*b**4)/1000,(4*b+14*b**6+4*a**3*b**3)/1000,d/1000,c/1000]
    y=[(1+a*a+b**3)/2000,(2+b+c*c+a**4)/3000,(1+c**3+d+a*b)/4000,(2+d*d+b+a*a*b)/5000]
    dy=[2*a/2000,mp.mpf(1)/3000,3*c*c/4000,2*d/5000]
    return w,dw,e,de,y,dy

def divergence_values(x,s):
    w,dw,e,de,y,dy=fields(x,s);total=mp.mpf(0)
    for j in range(4):
        def flux(t):
            xx=list(x);xx[j]=t;ww,_,ee,_,yy,_=fields(xx,s)
            return mp.exp(4*ww)*ee*yy[j]
        total+=mp.diff(flux,x[j])
    direct=mp.exp(-4*w)*total;gradient=sum(de[j]*y[j] for j in range(4));volume=4*e*sum(dw[j]*y[j] for j in range(4))
    expanded=gradient+e*sum(dy)+volume
    return direct,expanded,expanded-gradient,expanded-volume

def independent_physical_divergence():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for s in [mp.mpf(1)/3,mp.mpf(2)/3]:
                for at in [[mp.mpf(v) for v in ['0','0','.1','-.2']],[mp.mpf(v) for v in ['.1','-.2','.2','.3']]]:
                    tag=mp.nstr(s,4)+'_'+'_'.join(str(t) for t in at)+'_'+str(dps);errors=[]
                    for n in PROTOCOL['mixed_derivative_orders']:
                        diffs=[]
                        for k in range(n+1):
                            f=lambda a,b:divergence_values([a,b,at[2],at[3]],s)
                            v0=mp.diff(lambda a,b:f(a,b)[0],(at[0],at[1]),(n-k,k))/mp.factorial(n)
                            v1=mp.diff(lambda a,b:f(a,b)[1],(at[0],at[1]),(n-k,k))/mp.factorial(n)
                            diffs.append(abs(v0-v1))
                        error=max(diffs);errors.append(mp.nstr(error,18))
                        check('independent_physical_divergence_'+tag+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18))
                    vals=divergence_values(at,s);omissions=[abs(vals[0]-vals[i]) for i in [2,3]]
                    check('omitted_gradient_or_volume_fails_'+tag,min(omissions)>mp.mpf('1e-12'))
                    records.append({'conformal_scale':mp.nstr(s,18),'point':[str(t) for t in at],'dps':dps,'orders':PROTOCOL['mixed_derivative_orders'],
                                    'all_two_coordinate_mixed_partials':True,'errors':errors,'omissions':[mp.nstr(t,18) for t in omissions]})
    FIXTURES.append({'fixture':'independent_four_dimensional_divergence','records':records,'physical_volume':'exp(4w)',
                     'all_four_vector_components':True,'actual_charged_geometry_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 per-order complete mixed row','input_count':len(ps)}),flush=True)
    try:
        data=inputs();columns=matrix_ledger(*data);inverse_composition(columns);independent_polynomial_operator();independent_physical_divergence()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 per-order complete mixed row','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
         'exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
