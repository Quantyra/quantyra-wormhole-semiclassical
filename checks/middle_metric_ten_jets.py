"""Supporting verification for middle finite-comparison metric jets."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from itertools import product
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='middle-metric-ten-jets';PREVIOUS='seed-cold-actual-domain-eight-jets'
OWN=['checks/middle_metric_ten_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'metric_order':10,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'original_finite_comparison_retained':True,'coefficient_disk_radius':'0.005',
          'inverse_tau_disk_radius':'0.05','low_order_contracts_retained':True}


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
    check('prior_complete_cold_passed',p['status']=='passed' and len(p['checks'])==170 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('all_regional_seed_source_inputs',p['all_regional_seed_source_C8_bounded'])
    check('regional_scope_not_global',review()['middle_original_finite_metric_C10_bounded']
          and not review()['global_metric_C10_bounded'] and not review()['actual_finite_metric_source_H6_bounded'])
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))


def exp_upper(x,n=36):
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def seed_coefficients():
    z=S.symbols('z');c=R(1,10**10);e=R(1,10**24)
    abs_b=R(1,1)/R(59,10)+c*z**2+R(668,1000)*z**3+R(501,1000)*z**4+2*e*z**5+e*z**6
    for j,limit in enumerate([5,12,17,13],1):
        value=S.diff(abs_b,z,j).subs(z,R(102,100))
        check('polynomial_derivative_'+str(j),value<limit,value)
    check('real_p_initial',R(1,1)/R(59,10)+c<R(42,100)**2)
    check('p_bootstrap',R(42,100)+R(5,2)*R(5,1000)<R(45,100))
    check('z_lower_bootstrap',1/R(501,10)-R(45,100)*R(5,1000)>R(17,1000))
    check('z_upper_bootstrap',R(1000,1001)+R(45,100)*R(5,1000)<R(102,100))
    check('a_Riccati_bootstrap',R(202,100)+(12+4**2/2)*R(5,1000)<4)
    check('L_complex',R(45,100)/R(17,1000)<40)
    check('pz_complex',R(45,100)*R(102,100)<1)
    check('z_square_complex',R(102,100)**2<2)
    check('Lie_complex',R(102,100)*17/2<10)
    Scoef=5/R(17,1000)+2*R(45,100)**2/R(17,1000)**2
    check('S_complex',Scoef<10000,Scoef)
    cn=R(102,100)*17/6+4/R(59,10)/R(17,1000)**2+e*(2*R(102,100)**3+4*R(102,100)**4)
    check('radius_potential_complex',cn<20000,cn)
    check('radius_drift_complex',2+2*40<100)
    check('LLL_trace_complex',12*R(352,100000)/(24*3)<1)
    check('old_trace_value',10**11+R(501,10)*10**15+1<10**20)
    check('Ward_initial_value',R(3,2)*(10**11+R(501,10)*10**15+1)<10**20)
    check('mass_initial_value',2*R(68,10**29)*10**16+R(1,10**60)<R(1,10**8))
    check('Weyl_weighted_density_budget',1000*6*2*60**4<16*10**10)
    check('Weyl_pressure_and_reference_budget',16*10**10+10**12+10**10<10**14)
    check('Box_R_coefficient_budget',2*2562+12*366==9516)
    check('Weyl_square_coefficient_budget',R(61**2,3)<1241)
    check('Euler_coefficient_budget',4*(8**2+3*40)<1000)
    check('anomaly_and_reference_budget',2*10**6*60**2<10**10)
    check('combined_local_forces',10**16+3*2*10**16<10**18)
    ex=exp_upper(R(15,100))
    check('inverse_tau_complex',1001*ex<1200)
    check('tau_complex',ex<R(12,10))
    check('z_large_disk',ex<2)
    check('g0_derivative',(2+6*R(102,100))/R(59,10)+2*c*R(102,100)+4*e*R(102,100)**3<2)
    check('thermal_alpha_budget',R(1001,10**63)*(1200**2+R(3,2)*1200)<R(1,10**53))
    check('alpha_bootstrap',R(101,100)+R(13,10)*R(5,100)<3)
    invnorm=1200*sum(R(1,5)**j for j in range(11))
    check('inverse_tau_jet_norm',invnorm<1500,invnorm)
    FIXTURES.append({'fixture':'seed_and_local_coefficient_budgets','source_analyticity_assumed':False,
                     'inverse_tau_jet_norm':str(invnorm),'local_allocations_retained':True})


def exact_equations():
    z,x,e=S.symbols('z x e');b=S.Function('b')(z);rho=S.symbols('rho')
    zz=S.Function('z')(rho);pp=S.Function('p')(rho)
    aa=-S.diff(b,z)/S.Symbol('p')
    deriv=-S.Symbol('p')*S.diff(aa,z)-S.diff(b,z)/2*S.diff(aa,S.Symbol('p'))
    zero('Riccati_identity',deriv-(S.diff(b,z,2)-aa**2/2))
    F,g,T,Ed,EA,U,N=S.symbols('F gamma theta Ed EA u nu')
    dRh=-2*g*EA-6*F+6*g*z**2*(T+Ed)+N*z*S.diff(b,z,3)
    dR4=-2*g*EA+N*(-20*e*z**3+36*e*z**4)
    R40=-12*x/z**2+4*e*z**3-6*e*z**4
    cn=z*S.diff(b,z,3)/6+4*x/z**2+2*e*z**3-4*e*z**4
    zero('radius_anomaly_and_Lie_identity',(dRh-dR4)/6-N*R40/3-(-F+g*z**2*(T+Ed)+N*cn))
    check('anomaly_still_in_conformal_equation',S.diff(-dRh/2,EA)==g)
    check('omitted_Lie_term_detected',S.simplify(N*(-20*e*z**3+36*e*z**4)/6)!=0)
    zero('canonical_pressure_cancellation',(S.Symbol('P')-(S.Symbol('P')-T))-T)
    v,tau,ell=S.symbols('v tau ell')
    zero('finite_radius_cancellation',1+tau*(S.exp(ell*v)-1)-(1+ell*tau*v*(S.exp(ell*v)-1)/(ell*v)))
    FIXTURES.append({'fixture':'exact_reconstruction_identities','common_coordinate_Lie_terms_retained':True,
                     'anomaly_removed_only_by_identity_in_radius_equation':True})


def radial_ledger():
    g=R(68,10**29);T=[S.Integer(10)**20]+[S.Integer(2)*10**44]*8
    P=[S.Integer(10)**20];F=[R(1,10**8)]
    for n in range(8):
        P.append(sum(4*2**i*(P[n-i]+T[n-i]/2) for i in range(n+1))/(100*(n+1)))
        F.append((sum(40*2**i*F[n-i] for i in range(n+1))+2*g*sum(2**i*(P[n-i]+T[n-i]+10**16*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    U=[R(1,10000)]*6;N=U.copy()
    for n in range(4,9):
        conv=lambda c,X:sum(c*2**i*X[n-i] for i in range(n+1))
        drift=lambda c,X:sum(c*2**i*100*(n-i+1)*X[n-i+1] for i in range(n+1))
        U.append((drift(2,U)+conv(12,U)+conv(10,N)+3*F[n]+3*g*conv(2,T)+g*10**18*2**n)/(10000*(n+1)*(n+2)))
        N.append((drift(100,N)+conv(20000,U)+conv(20000,N)+F[n]+g*conv(2,T)+g*10**18*2**n)/(10000*(n+1)*(n+2)))
    check('all_radial_tangents_fit_H',max(U+N)<10**15,[str(S.N(v,18)) for v in U+N])
    for j,ub,nb in [(6,'4.086e13','1.362e13'),(7,'6.042e13','2.205e13'),(8,'9.222e13','3.825e13'),(9,'1.448e14','6.642e13'),(10,'2.329e14','1.161e14')]:
        check('printed_radial_bounds_'+str(j),U[j]<R(ub) and N[j]<R(nb))
    check('low_jets_preserved',U[:6]==N[:6]==[R(1,10000)]*6)
    # Independently reconstruct differentiated ODE coefficients with binomial
    # raw-derivative products, rather than the normalized convolution code.
    for n in [0,3,7]:
        raw=sum(S.binomial(n,i)*4*S.factorial(i)*200**i*S.factorial(n-i)*100**(n-i)*(P[n-i]+T[n-i]/2) for i in range(n+1))
        zero('raw_Ward_derivative_'+str(n),raw-P[n+1]*S.factorial(n+1)*100**(n+1))
    FIXTURES.append({'fixture':'radial_two_derivative_reconstruction','P':[str(v) for v in P],'F':[str(v) for v in F],
                     'U':[str(v) for v in U],'N':[str(v) for v in N],'high_radial_reserve':'1e15'})


def finite_and_tensor_ledger():
    H=S.Integer(10)**15;low=R(6,10000)
    check('high_jet_ideal_square_zero',6+6>10)
    check('V_low_norm',1500*low==R(9,10))
    check('V_high_norm',1500*(5*H+low)<7501*H)
    check('exp_low_majorant',exp_upper(R(9,10))<R(5,2))
    high=5*H*R(5,2)+low*R(5,2)*7501*H
    check('finite_high_product_budget',high==R('23.7515')*H)
    delta=R(15,10000)
    total=delta/(1-delta)+(delta+high)/(1-delta)
    check('finite_log_radial_norm',total<25*H,total)
    check('radial_finite_reserve',25*H<3*10**16)
    B={(0,k):S.factorial(k)*100**k for k in range(11)}
    for r in range(10):
        for k in range(10-r):
            B[r+1,k]=B[r,k+1]+r*sum(S.binomial(k,i)*6*S.factorial(i)*200**i*B[r,k-i] for i in range(k+1))
    amp=[S.factor(B[j,0]/(S.factorial(j)*100**j)) for j in range(11)]
    check('physical_tensor_amplification',max(amp)<R(312,100),amp)
    check('middle_original_finite_metric_C10_reserve',R(312,100)*3*10**16<10**17)
    check('connection_zero_order_retained',amp[2]>1)
    FIXTURES.append({'fixture':'finite_completion_and_physical_tensor_ledger','radial_finite_reserve':'3e16',
                     'tensor_amplifications':[str(v) for v in amp],
                     'tensor_amplifications_decimal':[str(S.N(v,18)) for v in amp],
                     'physical_reserve':'1e17','ell_factor_retained':True,'higher_tensors_symmetrized':False})


def independent_tensors():
    r,y=S.symbols('r y');coords=[r,y];f=r+r**3;warp=S.exp(r)
    metric=S.diag(1,warp**2);inverse=metric.inv()
    G={(a,i,j):S.simplify(sum(inverse[a,b]*(S.diff(metric[b,j],coords[i])+S.diff(metric[b,i],coords[j])-S.diff(metric[i,j],coords[b]))/2 for b in range(2))) for a,i,j in product(range(2),repeat=3)}
    tensors={1:{(i,):S.diff(f,coords[i]) for i in range(2)}}
    for order in range(1,4):
        tensors[order+1]={}
        for i in range(2):
            for inds,val in tensors[order].items():
                correction=sum(G[a,i,inds[p]]*tensors[order][inds[:p]+(a,)+inds[p+1:]] for p in range(order) for a in range(2))
                tensors[order+1][(i,)+inds]=S.simplify(S.diff(val,coords[i])-correction)
    zero('independent_radial_Hessian',tensors[2][0,0]-S.diff(f,r,2))
    zero('independent_angular_Hessian',tensors[2][1,1]-warp**2*S.diff(f,r))
    difference=S.simplify((tensors[3][0,1,1]-tensors[3][1,0,1]).subs(r,0))
    check('independent_third_tensor_not_symmetric',difference!=0,difference)
    raw=[abs(S.diff(f,r,k).subs(r,0)) for k in range(5)]
    B={(0,k):raw[k] for k in range(5)}
    for order in range(4):
        for k in range(4-order):B[order+1,k]=B[order,k+1]+3*order*B[order,k]
    norms=[]
    for order in range(1,5):
        squared=sum(v.subs(r,0)**2 for v in tensors[order].values())
        check('independent_full_tensor_bound_'+str(order),squared<=B[order,0]**2,squared)
        norms.append(str(squared))
    FIXTURES.append({'fixture':'independent_warped_coordinate_tensors','third_permutation_difference':str(difference),
                     'squared_norms':norms,'actual_seed_sampled':False})


def mul(a,b,n=10):return [sum(a[k]*b[j-k] for k in range(j+1)) for j in range(n+1)]


def exp_series(a,n=10):
    out=[mp.exp(a[0])]+[mp.mpf(0)]*n
    for j in range(1,n+1):out[j]=sum(k*a[k]*out[j-k] for k in range(1,j+1))/j
    return out


def log_series(a,n=10):
    inv=[1/a[0]]+[mp.mpf(0)]*n
    for j in range(1,n+1):inv[j]=-sum(a[k]*inv[j-k] for k in range(1,j+1))/a[0]
    return [mp.log(a[0])]+[sum(k*a[k]*inv[j-k] for k in range(1,j+1))/j for j in range(1,n+1)]


def independent_finite_contours():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            cache={}
            for high in [0,3,-3]:
                nu=[mp.mpf('0.00001')]*2+[mp.mpf(0)]*9;nu[6]=mp.mpf(high)
                tau=[mp.mpf('.01')*mp.mpf('.002')**j/mp.factorial(j) for j in range(11)]
                invtau=[100*mp.mpf('-.002')**j/mp.factorial(j) for j in range(11)]
                V=mul(nu,invtau);ev=exp_series(V);ev[0]-=1
                X=mul(tau,ev);X[0]+=1;coefs=log_series(X);cache[high]=coefs
                def exact(z):
                    tv=mp.mpf('.01')*mp.exp(z/500)
                    nv=mp.mpf('.00001')*(1+z)+high*z**6
                    return mp.log(1+tv*mp.expm1(nv/tv))
                for j in [6,10]:
                    radius=mp.mpf('.1')
                    integral=mp.quad(lambda a:exact(radius*mp.exp(1j*a))*mp.exp(-1j*j*a),[0,mp.pi,2*mp.pi])/(2*mp.pi*radius**j)
                    error=abs(integral-coefs[j])
                    check('independent_radius_contour_'+str(dps)+'_'+str(high)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,12))
                    records.append({'dps':dps,'high_coefficient':high,'order':j,'absolute_error':mp.nstr(error,18)})
            error=max(abs(cache[3][j]+cache[-3][j]-2*cache[0][j]) for j in range(11))
            check('independent_high_ideal_linearity_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,12))
            check('nonlinear_completion_not_replaced_'+str(dps),abs(cache[0][0]-mp.mpf('.00001'))>mp.mpf('1e-10'))
    FIXTURES.append({'fixture':'independent_finite_radius_complex_contours','records':records,
                     'actual_seed_sampled':False,'linear_radius_interpolation_used':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 middle metric ten jets','input_count':len(ps)}),flush=True)
    try:
        provenance();seed_coefficients();exact_equations();radial_ledger();finite_and_tensor_ledger()
        independent_tensors();independent_finite_contours()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 middle metric ten jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'middle_original_finite_metric_C10_bounded':True,'middle_physical_C10_reserve':'1e17',
         'middle_radial_tangent_C10_reserve':'1e15','global_metric_C10_bounded':False,
         'full_physical_covariant_tensors_bounded':True,'original_finite_comparison_retained':True,
         'all_regional_seed_source_C8_bounded':True,'local_source_allocations_retained':True,
         'cold_metric_C10_bounded':False,'outer_metric_C10_bounded':False,
         'endpoint_map_C11_bounded':False,'actual_finite_metric_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'original_low_order_contracts_changed':False,
         'endpoint_C5_contract_changed':False,'actual_signed_compatibility_zero_proved':False,
         'nonlinear_solution_constructed':False,'physical_stability_proved':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'ordinary_proof_machine_verified':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
