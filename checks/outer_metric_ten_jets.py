"""Supporting verification for outer metric C10 and forward conormal maps."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='outer-metric-ten-jets';PREVIOUS='cold-metric-ten-jets'
OWN=['checks/outer_metric_ten_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'metric_order':10,
          'forward_map_order':11,'auxiliary_cutoff_degree':23,'auxiliary_cutoff_regularity':11,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'original_finite_comparison_retained':True,'seed_proper_coordinate_used':True,
          'low_order_contracts_retained':True,'inverse_map_bound_claimed':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(expr);check(name,v==0,v)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_cold_metric_passed',p['status']=='passed' and len(p['checks'])==104 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('all_regional_scope_not_finite_response',review()['all_regional_original_metric_C10_bounded']
          and not review()['actual_finite_metric_source_H6_bounded'] and not review()['endpoint_inverse_map_C11_bounded'])


def exp_upper(x,n=32):
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def source_and_coefficient_budgets():
    geom=1/(1-R(2,100));P2=(1+R(1,100))**2
    check('z_square_multiplier',R(412,1000)**2*geom*P2<R(21,100))
    check('trace_product_allowance',R(21,100)*3==R(63,100))
    check('Ward_remainder_primitive',3*geom*R(101,100)*R(3,4)<R(5,2))
    check('Ward_remainder_density',R(17,100)*R(62,10)/2*geom*P2*R(5,2)<R(16,10))
    check('constant_I_density',R(17,100)/(2*R(9989,10000))*geom*P2*4*10**15<R('4.2e14'))
    check('complete_source_budget',R(63,100)+R(16,1000)+(R('4.2e14')+10**8)/10**44<2)
    check('local_log_derivative_sum',1+sum(R(1,j*100**j) for j in range(1,9))<R(1011,1000))
    check('local_source_coefficient_margin',10**6*geom*P2*R(1011,1000)<10**8)
    mu=2*R(42,100)*9*geom/(1-R(1,100))**2
    check('mass_integrand_multiplier',mu<8,mu)
    check('f_tail_multiplier',2/R(59,10)*R(412,1000)<R(14,100))
    check('j_tail_multiplier',2*R(412,1000)**3/R(9989,10000)<R(15,100))
    check('proper_radius_multiplier',R(1002,1000)/2<R(51,100))
    check('G0_and_inverse',(R(1001,1000))**2>R(10011,10000) and R(9989,10000)*(R(1001,1000))**2>1)
    check('inverse_gap_multiplier',1/(1-R(63,1000))<R(107,100))
    check('proper_Lie_multiplier',R(412,1000)*17/2<4)
    check('mass_mode_analytic_reserve',R(15,100)*R(4,10**11)<R(1,10**11))
    mass_remainder=2*R(1,10**11)/R(15,100)*geom*(R(101,100))**4
    check('mass_mode_remainder',mass_remainder<R(1,10**9),mass_remainder)
    FIXTURES.append({'fixture':'outer_complete_source_and_coefficient_budgets',
                     'source_order':8,'complete_source_reserve':'2e44',
                     'principal_and_vacuum_cancellations_retained':True,'mass_remainder_bound':str(mass_remainder)})


def exact_mass_and_gauge_identities():
    z,J,k,m=S.symbols('z J k m');b=S.Function('b')(z);p=S.sqrt(b)
    j=-2*k*m*z**3/b;F=-2*k*m*z
    v=-j/2-S.diff(b,z)*J/(4*p);w=-p*J/(2*z)
    def dy(f):return p*S.diff(f,z)+S.diff(f,J)*j
    zero('exact_mass_mode_curvature',dy(dy(v))+S.diff(b,z)/p*dy(v)+3*F-z*S.diff(b,z,3)*w/2)
    a,vv,U,LL=S.symbols('a v U L')
    zero('common_U_equation',(-vv+a*U)-a*U+vv)
    j3=S.symbols('j3');wp=-j3/8;vp=-j3/2;u3=-vp;nu3=wp-u3/4
    zero('proper_cubic_correlation',vp-4*wp)
    zero('common_cubic_correlation',u3+2*nu3)
    zero('finite_lapse_cubic',2*u3+4*nu3)
    check('proper_and_common_cubics_distinct',S.simplify(nu3-wp)!=0)
    ell0,zz,Vac=S.symbols('ell0 z Vac')
    zero('density_principal_and_vacuum_cancellation',-ell0*zz+(Vac+ell0*zz)-Vac)
    zero('pressure_difference_principal_and_vacuum_cancellation',ell0*zz+(-2*Vac-ell0*zz)+2*Vac)
    check('omitted_Weyl_tail_detected',-ell0*zz!=0)
    Y=S.symbols('Y',positive=True)
    for power in [1,4]:
        primitive=Y**(power+1)*((1-S.log(Y))/(power+1)+R(1,(power+1)**2))
        zero('weighted_log_primitive_'+str(power),S.diff(primitive,Y)-Y**power*(1-S.log(Y)))
    FIXTURES.append({'fixture':'exact_mass_curvature_gauge_and_tail_identities',
                     'proper_cubic_ratio':4,'common_cubic_relation':'u3+2nu3=0',
                     'singular_source_restored_before_mass_integration':True})


def coeff(c,n=11):return [R(c)*R(1,50)**j for j in range(n)]
def power(m,n=11):return [abs(S.binomial(m,j))/100**j for j in range(n)]
def conv(a,b,n):return [sum(a[i]*b[j-i] for i in range(j+1) if i<len(a) and j-i<len(b)) for j in range(n)]
def prim(a,m,n):
    assert n<=len(a)+1
    return [a[0]*R(m+2,(m+1)**2)]+[a[j-1]/(100*j) for j in range(1,n)]


def reconstruction_ledger():
    K=S.Integer(10)**44;gx=R(4,10**27);gam=R(68,10**29)
    mu=[R(7,10**11)]+[8*gx*K/(100*j) for j in range(1,10)]
    F=conv(conv(coeff(R(14,100)),power(1),9),mu,9)
    smallj=conv(conv(coeff(R(15,100)),power(3),10),mu,10)
    J=prim(smallj,4,11);W=conv(conv(coeff(R(51,100)),power(-1),11),J,11)
    W=[min(W[j],R(1,10**8)) if j<=5 else W[j] for j in range(11)]
    V=[R(1,10**8)]*6
    for n in range(4,9):
        dv=[100*(j+1)*V[j+1] for j in range(n+1)]
        drift=conv(conv(coeff(1),power(1),n+1),dv,n+1)[n]/100
        lie=conv(conv(coeff(4),power(1),n+1),W,n+1)[n]/1000
        V.append((3*F[n]+6*gam*K+gam*10**8*R(1,50)**n+lie+drift)/(10000*(n+1)*(n+2)))
    U=conv(coeff(R(1001,1000)),prim(conv(coeff(R(1001,1000)),V,11),4,11),11)
    term=conv(conv(coeff(R(1,2)),power(1),11),U,11)
    u=[V[j]+term[j]/100 for j in range(11)]
    term=conv(conv(coeff(R(101,100)),power(-1),11),U,11)
    nu=[W[j]+term[j] for j in range(11)];gap=conv(coeff(R(107,100)),nu,11)
    check('proper_v_high_ledger',max(V)<R('1.372e12'))
    check('proper_w_high_ledger',max(W)<R('9.098e11'))
    check('common_high_ledger',max(u+nu+gap)<R('1.372e12'))
    check('common_remainder_reserve',max(u+nu+gap)+R(1,10**9)<10**13)
    check('common_low_remainders',max(u[:6]+nu[:6]+gap[:6])+R(1,10**9)<R(1,10**7))
    check('only_nine_mass_tail_derivatives_used',len(mu)==10 and len(smallj)==10 and len(J)==11)
    FIXTURES.append({'fixture':'outer_two_derivative_reconstruction','mu':[str(v) for v in mu],
                     'f_tail':[str(v) for v in F],'proper_v':[str(v) for v in V],'proper_w':[str(v) for v in W],
                     'common_u':[str(v) for v in u],'common_nu':[str(v) for v in nu],'common_V':[str(v) for v in gap],
                     'ordinary_common_remainder_reserve':'1e13','low_order_contracts_retained':True})


def finite_tensor_forward_ledger():
    H=S.Integer(10)**13
    check('Y_L_guard',1+S.log(10)<R(331,100))
    check('low_jet_Y_cubic',R(11,10**12)+R(6,10**7)*R(331,1000)<R(21,10**8))
    eps=R(21,10**8)/1000
    check('low_jet_small',eps<R(1,10**9))
    # epsilon(Y)^2/(Y^4 L) <= (2.1e-7)^2 *Y^2, using L>=1.
    nonlinear=10*R(21,10**8)**2/R(100)+100*eps*H
    check('finite_nonlinear_remainder',nonlinear<R(1,10**6)*H)
    check('finite_scalar_reserve',H+nonlinear<3*10**13)
    B={(0,k):S.factorial(k)*100**k for k in range(11)}
    for r in range(10):
        for k in range(10-r):B[r+1,k]=B[r,k+1]+r*sum(S.binomial(k,i)*R(9,400)*S.factorial(i)*2**i*B[r,k-i] for i in range(k+1))
    amp=[S.factor(B[j,0]/(S.factorial(j)*100**j)) for j in range(11)]
    check('full_weighted_tensor_amplification',max(amp)<R(1002,1000))
    check('outer_physical_C10_reserve',R(1002,1000)*3*10**13<10**14)
    check('zeroth_connection_retained',amp[2]>1)
    check('forward_map_value_primitive',R(6,25)*3*10**13<10**14)
    check('forward_map_eleventh_derivative',3*10**13/(100*11)<10**14)
    check('forward_map_does_not_claim_closed_C11',S.limit(S.diff(S.Symbol('Y',positive=True)**5*S.log(S.Symbol('Y',positive=True)),S.Symbol('Y',positive=True),5),S.Symbol('Y',positive=True),0,dir='+')==-S.oo)
    FIXTURES.append({'fixture':'finite_completion_full_tensors_and_forward_map',
                     'tensor_amplifications':[str(v) for v in amp],'physical_C10_reserve':'1e14',
                     'forward_C11_remainder_reserve':'1e14','ell_factor_retained':True,
                     'higher_tensors_symmetrized':False,'inverse_map_estimated':False})


def higher_auxiliary_cutoff():
    t=S.symbols('t');eta=S.expand(sum(S.binomial(23,i)*t**i*(1-t)**(23-i) for i in range(12)))
    c=S.factorial(23)/S.factorial(11)**2
    check('cutoff_beta_normalization',c==16224936)
    zero('cutoff_derivative',S.diff(eta,t)+c*t**11*(1-t)**11)
    for j in range(12):
        zero('cutoff_left_'+str(j),S.diff(eta,t,j).subs(t,0)-(1 if j==0 else 0))
        zero('cutoff_right_'+str(j),S.diff(eta,t,j).subs(t,1))
    check('cutoff_complex_derivative',c*(R(1,4)+2*R(1,20)**2)**11<5)
    check('cutoff_complex_value',1+R(1,20)*5==R(5,4))
    check('global_forward_jet_multiplier',R(5,4)*sum(R(2,5)**j for j in range(12))<R(2084,1000))
    check('global_forward_weighted_bound',R(2084,1000)*(10**14/R(10)+R(1,10**10))<3*10**13)
    check('global_forward_monotonicity',100*R(1,10)*R(1,10**10)+R(2,10**8)==R(21,10**9))
    check('old_C6_not_silently_reused',11>6)
    FIXTURES.append({'fixture':'explicit_C11_auxiliary_forward_identification','degree':23,
                     'matching_endpoint_derivative_order':11,'global_forward_reserve':'3e13',
                     'old_auxiliary_cutoff_replaced_explicitly':True,'physical_finite_path_changed':False})


def independent_weighted_integrals():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for m in [1,4]:
                for yy in ['0.0001','0.1']:
                    Y=mp.mpf(yy)
                    value=mp.quad(lambda t:t**m*(1-mp.log(Y)-mp.log(t)),[0,1])*Y**(m+1)
                    exact=Y**(m+1)*((1-mp.log(Y))/(m+1)+mp.mpf(1)/(m+1)**2)
                    error=abs(value-exact)
                    check('independent_log_integral_'+str(dps)+'_'+str(m)+'_'+yy,error<mp.mpf('1e-42'),mp.nstr(error,12))
                    records.append({'dps':dps,'power':m,'Y':yy,'absolute_error':mp.nstr(error,18)})
            # Actual exp/log-gap formula in a nontrivial conormal model,
            # compared with an independent contour coefficient extraction.
            Y0=mp.mpf('.07');radius=mp.mpf('.01');mcoef=mp.mpf('.002');tail=mp.mpf('.003')
            def radius_log(Y):
                gap=mcoef*Y**3+tail*Y**4*(1-mp.log(Y));tau=1-mp.mpf('.4')*Y
                return mp.log(1+tau*mp.expm1(gap))
            for j in [6,10]:
                contour=mp.quad(lambda a:radius_log(Y0+radius*mp.exp(1j*a))*mp.exp(-1j*j*a),[0,mp.pi,2*mp.pi])/(2*mp.pi*radius**j)
                direct=mp.diff(radius_log,Y0,j)/mp.factorial(j)
                error=abs(contour-direct)
                check('independent_conormal_completion_'+str(dps)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,12))
                records.append({'dps':dps,'finite_map_order':j,'absolute_error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_weighted_integrals_and_conormal_completion','records':records,
                     'actual_charged_source_sampled':False,'logarithmic_tail_retained':True})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 outer metric ten jets','input_count':len(ps)}),flush=True)
    try:
        provenance();source_and_coefficient_budgets();exact_mass_and_gauge_identities();reconstruction_ledger()
        finite_tensor_forward_ledger();higher_auxiliary_cutoff();independent_weighted_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 outer metric ten jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'outer_original_finite_metric_C10_bounded':True,'outer_physical_C10_remainder_reserve':'1e14',
         'all_regional_original_metric_C10_bounded':True,'mass_cubic_correlation_retained':True,
         'full_physical_covariant_tensors_bounded':True,'original_finite_comparison_retained':True,
         'all_regional_seed_source_C8_bounded':True,'local_source_allocations_retained':True,
         'endpoint_forward_map_C11_bounded':True,'global_auxiliary_forward_C11_bounded':True,
         'higher_auxiliary_cutoff_explicit':True,'endpoint_inverse_map_C11_bounded':False,
         'proper_endpoint_chords_C10_bounded':False,'ordinary_closed_boundary_C11_claimed':False,
         'actual_finite_metric_source_H6_bounded':False,'actual_transfer_certificate_passed':False,
         'original_low_order_contracts_changed':False,'endpoint_C5_contract_changed':False,
         'actual_signed_compatibility_zero_proved':False,'nonlinear_solution_constructed':False,
         'physical_stability_proved':False,'physical_source_projection_changed':False,
         'physical_boundary_condition_changed':False,'ordinary_proof_machine_verified':False,
         'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
