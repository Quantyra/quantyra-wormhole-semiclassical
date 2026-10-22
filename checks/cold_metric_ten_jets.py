"""Supporting checks for the weighted cold finite metric C10 proof."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from itertools import product
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='cold-metric-ten-jets';PREVIOUS='middle-metric-ten-jets'
OWN=['checks/cold_metric_ten_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'metric_order':10,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'original_finite_comparison_retained':True,'bolt_weight_retained':True,
          'near_region':'0<=rho<=1','away_region':'rho>=0.5','low_order_contracts_retained':True}


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
    check('prior_middle_metric_passed',p['status']=='passed' and len(p['checks'])==89 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('cold_scope_not_global',review()['cold_original_finite_metric_C10_bounded'] and not review()['global_metric_C10_bounded'] and not review()['actual_finite_metric_source_H6_bounded'])


def exp_upper(x,n=40):
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def series_major(kind,w,d=0):
    def term(k):
        c=1/S.factorial(2*k+1) if kind=='n' else S.Integer(2)**(2*k+1)/S.factorial(2*k+2)
        if kind=='K':c=S.Integer(2)**(2*k+3)/S.factorial(2*k+4)
        return c*S.factorial(k)/S.factorial(k-d)*w**(k-d)
    return sum(term(k) for k in range(d,30))+2*term(30)


def seed_and_bolt_coefficients():
    h=R(1,10**30);s=R(1,1000)
    check('local_cold_Gronwall',exp_upper(R(1,100))*(2+R(1,100)*R(1001,1000)*s)<R(21,10))
    check('local_tau_bootstrap',4*s<R(5,1000))
    check('L_and_pz',R(4,1000)/(1-R(4,1000))<1 and R(4,1000)*(1+R(4,1000))<1)
    check('S_small',2*(R(11,1000)/R(995,1000)+2*R(13,100000)/R(995,1000)**2)<2)
    check('radius_potential',R(1004,1000)*R(81,10)/6+4/R(59,10)/R(996,1000)**2+R(1,10**20)<10)
    check('real_cold_curvature_positive',3*R(1334,1000)*s+6*R(501,1000)*s**2+R(1,10**20)<R(1,100))
    check('weighted_local_subtraction',2*10**18*4+2*10**20*R(1001,1000)*s<10**22)
    check('weighted_local_forces',10**22+3*2*10**22<10**24)
    check('away_beta_bootstrap',4+R(5,1000)*(6+25)<5)
    check('near_far_overlap',R(1,2)<1)
    check('centered_residual_budget',10*(100*h)**2+2*R(1001,1000)*h**2<110000*h**2)
    check('centered_complex_sharpening',4*h+4*110000*h**2<5*h)
    check('normal_metric_perturbation',51*10*4<10000)
    check('normal_metric_Cauchy_error',10*10000*h*2/(4-R(9,4))**3<R(1,10**24))
    eps=R(1,10**24)
    entries=[]
    for rad,limits in [(R(6,5),[R(13,10),1,1,1]),(R(3,2),[R(3,2),2,1,1])]:
        values=[series_major('n',rad**2),series_major('psi',rad**2,1),series_major('K',rad**2),series_major('K',rad**2,1)]
        check('normal_series_'+str(rad),all(v+eps<lim for v,lim in zip(values,limits)),[str(S.N(v,15)) for v in values])
        entries.append({'radius':str(rad),'upper_bounds':[str(v+eps) for v in values]})
    check('series_tail_ratio',R(4*9,4)*(31)/(30*(63*64))<R(1,2))
    check('near_regular_drift',2*R(6,5)/R(7,10)<10)
    check('Ward_kernel_complex_bound',(R(13,10)**2+R(6,5)**2)/R(7,10)**2<10)
    check('Cartesian_polydisk',1+R(3,2)*R(1,10)<R(6,5))
    check('Christoffel_complex_norm',2*R(3,2)*8*R(6,5)+R(3,2)*3*R(6,5)+17*R(6,5)**3<100)
    check('Christoffel_tensor_Cauchy',R(3,2)/R(1,10)<20)
    ex=exp_upper(R(15,100));ex2=exp_upper(R(3,10))
    check('inverse_tau_uniform',ex<R(6,5))
    check('thermal_inverse_tau_ratio',ex2<R(135,100))
    check('alpha_derivative_at_bolt',R(2,1000)+R(1002,1000)*R(101,100)*ex2+R(1,10**20)<2)
    check('alpha_bootstrap',R(101,100)+2*R(5,100)<3)
    check('weighted_inverse_tau_norm',R(6,5)*sum(R(1,5)**j for j in range(11))<R(3,2))
    check('weight_to_tau_ratio',R(1001,1000)*81<82)
    FIXTURES.append({'fixture':'uniform_seed_bolt_and_connection_coefficients','series_entries':entries,
                     'inverse_temperature_loss':False,'quantum_source_analyticity_assumed':False})


def regular_identities():
    r,t=S.symbols('r t',positive=True)
    for k in range(6):
        force=r**(2*k);v=r**(2*k+2)/(2*k+2)**2
        zero('flat_regular_inverse_'+str(k),S.diff(v,r,2)+S.diff(v,r)/r-force)
        for j in [0,3,8]:
            rhs=S.diff(force,r,j)-S.integrate(t**(j+1)*S.diff(force,r,j).subs(r,t*r),(t,0,1))
            zero('two_derivative_gain_'+str(k)+'_'+str(j),S.diff(v,r,j+2)-rhs)
    trace=2*S.sinh(r)**2;Q=S.sinh(r)**2/2
    zero('hyperbolic_Ward_solution',S.diff(Q,r)+2*S.coth(r)*Q-S.coth(r)*trace)
    q,theta0,T,Er,Ed=S.symbols('Q theta0 T Er Ed')
    zero('product_density_subtraction',(theta0/2+q)-(theta0+T)+(theta0/2+Er)-(q-T+Er))
    zero('product_trace_subtraction',theta0+T+(-theta0+Ed)-(T+Ed))
    # The first-order flat inverse without its regular integration condition
    # admits log(r); this must remain excluded at the bolt.
    zero('singular_homogeneous_mode',S.diff(S.log(r),r,2)+S.diff(S.log(r),r)/r)
    check('singular_homogeneous_mode_excluded',S.limit(S.log(r),r,0,dir='+')==-S.oo)
    FIXTURES.append({'fixture':'regular_inverse_and_product_subtraction','monomials':6,
                     'tested_source_derivative_orders':[0,3,8],'bolt_integration_constant_fixed':True})


def reconstruction_ledger():
    K=S.Integer(2)*10**44;g=R(68,10**29)
    near=[5*K*sum(R(1,10)**i for i in range(j+1)) for j in range(9)]
    away=[K/2]
    for n in range(8):away.append(sum(5*2**i*(K+2*away[n-i]) for i in range(n+1))/(100*(n+1)))
    Q=[max(a,b) for a,b in zip(near,away)];F=[R(1,10**10)]
    check('Ward_near_dominates',all(a>b for a,b in zip(near,away)))
    for n in range(8):
        F.append((sum(2**i*F[n-i] for i in range(n+1))+2*g*sum(2**i*(Q[n-i]+K+10**22*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    U=[R(1,10000)]*6;N=U.copy()
    for n in range(4,9):
        C=lambda c,X:sum(c*2**i*X[n-i] for i in range(n+1))
        D=lambda c,X:sum(c*2**i*100*(n-i+1)*X[n-i+1] for i in range(n+1))
        U.append(2*(D(12,U)+C(12,U)+C(10,N)+3*F[n]+3*g*C(2,[K]*9)+g*10**24*2**n)/(10000*(n+1)*(n+2)))
        N.append(2*(D(12,N)+C(2,U)+C(10,N)+F[n]+g*C(2,[K]*9)+g*10**24*2**n)/(10000*(n+1)*(n+2)))
    check('cold_radial_reserve',max(U+N)<2*10**15)
    check('tenth_radial_entries',U[10]<R('1.003e15') and N[10]<R('3.341e14'))
    check('old_low_radial_jets_kept',U[:6]==N[:6]==[R(1,10000)]*6)
    FIXTURES.append({'fixture':'weighted_radial_reconstruction','Ward_near':[str(v) for v in near],
                     'Ward_away':[str(v) for v in away],'F':[str(v) for v in F],
                     'U':[str(v) for v in U],'N':[str(v) for v in N],'radial_tangent_reserve':'2e15'})


def finite_and_tensor_ledgers():
    H=S.Integer(2)*10**15;aw=R(331,100000)
    check('finite_V_low',R(12,10000)*82<R(1,10))
    check('finite_V_high',(10*H+R(12,10000))*82<821*H)
    check('finite_h_majorants',exp_upper(R(1,10))<R(6,5))
    low=R(72,100000);delta=low*aw
    high=low+(6+low*821)*H
    check('finite_X_low',delta<R(24,10**7))
    check('finite_X_high',high<7*H)
    check('finite_log_radial_reserve',(low+high)/(1-delta)<10*H<3*10**16)
    C={(1,1):S.Integer(1)};A=[S.Integer(1),S.Integer(1)]
    for m in range(1,10):
        for j in range(1,m+2):C[m+1,j]=C.get((m,j-1),0)+(4*m-j)*C.get((m,j),0)
        A.append(sum(C[m+1,j]/S.factorial(m+1-j) for j in range(1,m+2)))
    check('flat_radial_to_Cartesian',A[10]<R('4.868e6'))
    def amps(initial,constant,scale):
        B={(0,k):initial[k]*S.factorial(k)*100**k for k in range(11)}
        for r in range(10):
            for k in range(10-r):B[r+1,k]=B[r,k+1]+r*sum(S.binomial(k,i)*constant*S.factorial(i)*scale**i*B[r,k-i] for i in range(k+1))
        return [S.factor(B[j,0]/(S.factorial(j)*100**j)) for j in range(11)]
    near=amps(A,100,20);away=amps([S.Integer(1)]*11,15,200)
    check('near_full_tensor_reserve',max(near)<R('1.02e7') and near[10]<10186500)
    check('away_full_tensor_reserve',max(away)<R('8.64'))
    check('connection_zero_order_retained',near[2]>A[2] and away[2]>1)
    check('complete_cold_C10_reserve',R('1.02e7')*3*10**16<10**24)
    FIXTURES.append({'fixture':'original_finite_completion_and_tensor_conversion','radial_finite_reserve':'3e16',
                     'Cartesian_amplifications':[str(v) for v in A],
                     'near_tensor_amplifications':[str(v) for v in near],'away_tensor_amplifications':[str(v) for v in away],
                     'physical_C10_reserve':'1e24','cold_weight_retained':True,'ell_factor_retained':True,
                     'higher_tensors_symmetrized':False})


def independent_Cartesian_connection():
    x,y=S.symbols('x y');coords=[x,y];w=x*x+y*y
    # An independent smooth nonconstant radial metric; the exact identity
    # must hold for every analytic psi with psi(0)=1.
    psi=1+w/3+w*w/7;K=(psi-1)/w
    metric=S.Matrix(2,2,lambda i,j:psi*(1 if i==j else 0)-K*coords[i]*coords[j]);inv=metric.inv()
    Pw=R(1,3)+2*w/7;Kw=R(1,7)
    A=Pw/psi;B=-K-Pw;D=2*Pw*K/psi-Kw
    errors=[];radial=[]
    for a,i,j in product(range(2),repeat=3):
        direct=sum(inv[a,b]*(S.diff(metric[b,j],coords[i])+S.diff(metric[b,i],coords[j])-S.diff(metric[i,j],coords[b]))/2 for b in range(2))
        formula=A*(coords[i]*(1 if a==j else 0)+coords[j]*(1 if a==i else 0))+B*coords[a]*(1 if i==j else 0)+D*coords[a]*coords[i]*coords[j]
        errors.append(S.factor(direct-formula))
    check('independent_Cartesian_Christoffel',all(v==0 for v in errors))
    zero('Cartesian_radial_geodesics',2*A+B+D*w)
    for m in [2,4,6]:
        f=(x*x+y*y)**(m//2)
        tensor=[S.diff(f,x,m-k,y,k).subs({x:0,y:0}) for k in range(m+1)]
        squared=sum(S.binomial(m,k)*tensor[k]**2 for k in range(m+1))
        check('even_radial_polynomial_tensor_'+str(m),squared<=S.factorial(m)**2*10**m)
        radial.append({'order':m,'tensor_squared_norm':str(squared)})
    FIXTURES.append({'fixture':'independent_Cartesian_seed_geometry','Christoffel_components':8,
                     'radial_polynomial_tensors':radial,'actual_seed_sampled':False})


def independent_regular_integrals():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            v=lambda r:r*r*mp.hyper([1,1],[2,2],r*r)/4
            for radial in ['0','0.2','0.9']:
                r=mp.mpf(radial)
                for j in [0,4,8]:
                    # d^j exp(r^2) from an independent finite Hermite sum.
                    force=lambda x:mp.exp(x*x)*mp.factorial(j)*sum((2*x)**(j-2*k)/(mp.factorial(k)*mp.factorial(j-2*k)) for k in range(j//2+1))
                    integral=force(r)-mp.quad(lambda t:t**(j+1)*force(t*r),[0,1])
                    direct=mp.diff(v,r,j+2)
                    error=abs(integral-direct)
                    check('independent_regular_inverse_'+str(dps)+'_'+radial+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,12))
                    records.append({'dps':dps,'rho':radial,'source_order':j,'absolute_error':mp.nstr(error,18)})
            for rr,ii in [('0.1','0.03'),('0.9','-0.03')]:
                r=mp.mpc(rr,ii)
                kernel=lambda t:(mp.sinh(t*r)*mp.cosh(t*r)*r)/(mp.sinh(r)**2)
                value=mp.quad(lambda t:kernel(t)*(2*mp.sinh(t*r)**2),[0,1])
                error=abs(value-mp.sinh(r)**2/2)
                check('independent_Ward_integral_'+str(dps)+'_'+rr,error<mp.mpf('1e-42'),mp.nstr(error,12))
    FIXTURES.append({'fixture':'independent_regular_inverse_and_Ward_integrals','records':records,
                     'includes_exact_bolt':True,'actual_charged_source_sampled':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold metric ten jets','input_count':len(ps)}),flush=True)
    try:
        provenance();seed_and_bolt_coefficients();regular_identities();reconstruction_ledger()
        finite_and_tensor_ledgers();independent_Cartesian_connection();independent_regular_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold metric ten jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'cold_original_finite_metric_C10_bounded':True,'cold_physical_C10_reserve':'1e24',
         'cold_radial_tangent_C10_reserve':'2e15','middle_original_finite_metric_C10_bounded':True,
         'full_physical_covariant_tensors_bounded':True,'original_finite_comparison_retained':True,
         'all_regional_seed_source_C8_bounded':True,'local_source_allocations_retained':True,
         'cold_weight_retained':True,'bolt_regular_inverse_used':True,
         'global_metric_C10_bounded':False,'outer_metric_C10_bounded':False,
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
