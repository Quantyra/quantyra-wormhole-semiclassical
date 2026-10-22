"""Supporting checks for the moving-normal-cap cold local source through eight."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-cold-local-eight-jets';PREVIOUS='seed-cold-weighted-inverse-jets'
OWN=['checks/seed_cold_local_eight_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','local_source_order':8,'observation_radius':'1e-2',
          'u_radius':'5e-3','thermal_coefficient_radius':'1e-4',
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'gluing_included':False,'actual_domain_error_included':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    value=S.simplify(expr);check(name,value==0,value)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_weighted_inverse_passed',p['status']=='passed' and len(p['checks'])==99 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('moving_local_scope_only',review()['local_cold_source_order']==8
          and review()['local_cold_source_reserve']=='1e8'
          and not review()['complete_cold_source_C8_bounded'] and not review()['cold_gluing_C8_bounded']
          and not review()['actual_domain_C8_bounded'] and not review()['goal_complete'])


def polynomial_and_invariant_domains():
    z,b,c,e,delta=S.symbols('z b c e delta')
    poly=z*z+b*z**3+c*z**4+4*e*z**5+e*z**6-delta*(1+z)**3
    x,d,C=S.symbols('x d C')
    u=1+z
    original=((1-u)**2*(x*(1+2*u+3*u*u)+u*u+d*u**4))/C-delta*u**3
    exact=poly.subs({b:(8*x+2+4*d)/C,c:(3*x+1+6*d)/C,e:d/C})
    zero('centered_polynomial_identity',(original-exact).subs(C,1+6*x+d))
    bmax=R(1334,1000);cmax=R(501,1000);emax=R(1,10**24);dmax=R(1,10000)
    real_d=R(1,10**12)
    check('centered_cubic_coefficient',(8000+2+4*real_d)/(6001+real_d)<bmax)
    check('centered_quartic_coefficient',(3000+1+6*real_d)/(6001+real_d)<cmax)
    r=S.symbols('r');arg=16*r-3
    step=sum(S.binomial(21,i)*arg**i*(1-arg)**(21-i) for i in range(11))
    check('fixed_normal_cap_radii',step.subs(r,R(3,16))==1 and step.subs(r,R(1,4))==0)
    major=lambda order,r:S.diff(poly,z,order).as_poly(z,b,c,e,delta)
    def budget(order,r):
        expression=S.diff(poly,z,order).expand()
        return sum(abs(coef)*r**powers[0]*bmax**powers[1]*cmax**powers[2]*emax**powers[3]*dmax**powers[4]
                   for powers,coef in S.Poly(expression,z,b,c,e,delta).terms())
    values={}
    for order,target in [(1,R(7,100)),(2,R(231,100)),(3,S.Integer(10))]:
        val=budget(order,R(3,100));check('Goursat_guard_'+str(order),val<target,S.N(val,18))
    for order,target in [(0,R(13,100000)),(1,R(11,1000)),(2,R(205,100)),(3,R(81,10)),(4,S.Integer(13))]:
        val=budget(order,R(5,1000));check('local_polynomial_guard_'+str(order),val<target,S.N(val,18));values[str(order)]=str(val)
    check('beta_disk',R(13,100000)/4<R(1,10000))
    check('fixed_p_section',R(8,1000)<R(1,100) and R(13,100000)/4/R(8,1000)<R(1,100))
    check('homogeneous_Goursat_data',all(S.diff(poly,z,j).subs({z:0,delta:0})==target for j,target in [(0,0),(1,0),(2,2)]))
    p,pbar,uu=S.symbols('p pbar uu')
    invariant=1/(1-2*(uu-1)-100*p*pbar)
    zero('rotation_invariant_parameter_reduction',invariant.subs({p:S.I*p,pbar:-S.I*pbar},simultaneous=True)-invariant)
    check('unbalanced_parameter_negative_control',S.simplify((invariant+p).subs({p:S.I*p,pbar:-S.I*pbar},simultaneous=True)-(invariant+p))!=0)
    FIXTURES.append({'fixture':'joint_polynomial_and_invariant_domains','local_polynomial_budgets':values,
                     'beta_radius':'1e-4','square_root_branch_used':False,'global_HH_state_varied':False})


def local_and_observation_budgets():
    n,q=S.symbols('n q',positive=True)
    zero('exact_angular_telescoping',(q+2*n)/(n*n*(n+q)**2)-(1/n**2-1/(n+q)**2)/q)
    ratio=R(352,100000)
    Pbound=R(205,100)**2+2*R(13,100000)*13+2*R(11,1000)*R(81,10)
    check('finite_polynomial_budget',Pbound<R(44,10),Pbound)
    born=4*ratio**3*14*10**10;finite=R(1000,480*3)*ratio**2*R(44,10)
    check('complete_local_analytic_budget',born+finite<25000,[born,finite])
    h=R(1,10**30);thermal=25000*2*h*h/(R(1,10000)-2*h*h)
    check('thermal_value_difference_retained',0<thermal/h<R(1,10**20))
    check('zeroth_local_difference',25000/R(4,1000)+R(1,10**20)<10**7)
    rho=R(1,100);hshift=2*rho/(1-2*rho)
    check('local_complex_H_shift',hshift<R(20409,10**6))
    check('uniformizing_inverse_image',R(20409,10**6)<R(17,100))
    sshift=R(40818,10**6);smax=(1+sshift)/1000
    check('local_u_displacement',sshift/(1-smax)<R(5,100))
    check('local_u_domain',smax/(1-smax)<R(11,10000))
    check('positive_local_Cauchy_budget',25000*R(5,100)/(R(5,1000)-R(11,10000))<400000)
    check('full_parameter_tensor_factor',150**2>2/rho**2)
    sd=R(1,1000)*(R(3045,1000)*R(2303,1000)+R(675,100))
    check('uniform_seed_complex_amplitude',R(41,10000)+sd/(1-sd)+R(1,10**24)<R(2,100))
    FIXTURES.append({'fixture':'full_angular_finite_and_observation_budgets','Born_angular_bound':str(born),
                     'finite_term_bound':str(finite),'joint_local_reserve':25000,
                     'zeroth_cold_reserve':10**7,'positive_parameter_reserve':400000,
                     'full_angular_tower_retained':True,'common_finite_part_retained':True})


def connection_and_independent_controls():
    w,v,z=S.symbols('w v z',real=True)
    metric=S.Matrix([[1+w*w,-v*w],[-v*w,1]]);inv=metric.inv();det=metric.det()
    deriv=lambda expr,k:0 if k==0 else S.diff(expr,w)*z+S.diff(expr,z)*(1-v*v)*w
    gamma=[[[S.factor(sum(inv[i,l]*(deriv(metric[l,k],j)+deriv(metric[l,j],k)-deriv(metric[j,k],l))/2 for l in range(2)))
             for k in range(2)] for j in range(2)] for i in range(2)]
    expected=[[[-v*w*w/z,w/z],[w/z,-v/z]],
              [[-(1+w*w)*w/z,v*w*w/z],[v*w*w/z,-v*v*w/z]]]
    for i in range(2):
        for j in range(2):
            for k in range(2):
                numerator=S.factor((gamma[i][j][k]-expected[i][j][k])*det*z)
                zero('Christoffel_'+str(i)+str(j)+str(k),S.rem(numerator,z*z-det,z))
    at=lambda expr:S.simplify(expr.subs({w:0,z:1,v:R(4,5)}))
    hessian=at(-expected[0][1][1])
    check('nonnormal_Hessian_control',hessian==R(4,5) and hessian!=0)
    t=S.symbols('t',real=True);vv=R(4,5);cc=R(3,5)
    geodesic_b1=S.atanh(vv)-S.atanh(vv/S.sqrt(1+cc*cc*S.sinh(t)**2))
    zero('independent_geodesic_Hessian',S.diff(geodesic_b1,t,2).subs(t,0)-hessian)
    third={}
    for k in range(2):
        for i in range(2):
            for j in range(2):
                third[k,i,j]=at(-deriv(expected[0][i][j],k)
                    +sum(expected[l][k][i]*expected[0][l][j]+expected[l][k][j]*expected[0][i][l] for l in range(2)))
    check('nonsymmetric_third_tensor',third[1,0,1]==-1 and third[0,1,1]==0)
    s1,s2=S.symbols('s1 s2',real=True);grad=[s1,s2]
    norm=sum(((1 if i==j else 0)*grad[k]+(1 if i==k else 0)*grad[j]-(1 if j==k else 0)*grad[i])**2
             for i in range(2) for j in range(2) for k in range(2))
    zero('Euclidean_conformal_connection_norm',norm-4*(s1*s1+s2*s2))
    rad=R(5,1000);wb=R(501,100000)
    check('connection_sinh_bound',rad/(1-rad*rad)<wb)
    check('connection_cosh_bound',rad*rad/(2*(1-rad*rad))<R(13,10**6))
    check('metric_and_inverse_bounds',wb+wb*wb<R(6,1000) and 1/(1-R(6,1000))<R(1007,1000))
    check('bare_connection_bound',(1+4*wb**2+5*wb**4+wb**6)/(1-R(13,10**6))**2<4)
    check('conformal_connection_factor',2+2*R(1415,1000)*R(6,1000)/(1-R(6,1000))<R(21,10))
    check('full_connection_bound',2+R(21,10)*R(1415,1000)*R(2,100)/rad<20)
    check('connection_tensor_Cauchy',300**2>2/rad**2)
    check('physical_tensor_conversion',1/(1-8*R(2,100))<R(6,5))
    FIXTURES.append({'fixture':'explicit_nonnormal_connection_and_independent_controls',
                     'Hessian_22_of_b1':str(hessian),'omitted_Gamma0_error':str(hessian),
                     'third_212':str(third[1,0,1]),'third_122':str(third[0,1,1]),
                     'full_connection_bound':20,'connection_derivative_scale':300})


def full_tensor_ledger():
    B={(0,k):(10**7 if k==0 else 400000)*S.factorial(k)*150**k for k in range(9)}
    gamma=lambda a:20*S.factorial(a)*300**a
    for r in range(8):
        for k in range(8-r):
            B[r+1,k]=B[r,k+1]+r*sum(S.binomial(k,a)*gamma(a)*B[r,k-a] for a in range(k+1))
    out=[]
    for j in range(9):
        value=R(6,5)*B[j,0]/(S.factorial(j)*100**j)
        check('complete_physical_local_bound_'+str(j),value<10**8,value)
        out.append({'order':j,'exact':str(value),'decimal':str(S.N(value,18))})
    zero('eighth_exact_entry',R(6,5)*B[8,0]/(S.factorial(8)*100**8)-R(5386822704,125))
    FIXTURES.append({'fixture':'full_covariant_local_source_ledger','entries':out,
                     'Gamma_zero_order_included':True,'higher_tensors_symmetrized':False})


def independent_invariant_contours():
    count=9
    def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(count)]
    def inv(a):
        out=[1/a[0]]
        for j in range(1,count):out.append(-sum(a[i]*out[j-i] for i in range(1,j+1))/a[0])
        return out
    one=[S.Integer(1)]+[S.Integer(0)]*8
    s=[(R(1,1000) if j%2==0 else R(8,10000))/S.factorial(j) for j in range(count)]
    denom=s.copy();denom[0]+=1;u=inv(denom);zz=u.copy();zz[0]-=1
    z2=mul(zz,zz);z3=mul(z2,zz);z4=mul(z3,zz);u3=mul(mul(u,u),u)
    polynomial=[z2[j]+R(4,3)*z3[j]+R(1,2)*z4[j]-R(1,10**8)*u3[j] for j in range(count)]
    coeffs=inv([one[j]-2*zz[j]-25*polynomial[j] for j in range(count)])
    uu=S.symbols('uu');poly=(uu-1)**2+R(4,3)*(uu-1)**3+R(1,2)*(uu-1)**4-R(1,10**8)*uu**3
    direct=1/(1-2*(uu-1)-25*poly);u0=R(1000,1001)
    u1=-R(8,10000)/(1+R(1,1000))**2
    u2=-R(1,1000)/(1+R(1,1000))**2+2*R(8,10000)**2/(1+R(1,1000))**3
    zero('independent_invariant_chain_first',S.diff(direct,uu).subs(uu,u0)*u1-coeffs[1])
    zero('independent_invariant_chain_second',(S.diff(direct,uu,2).subs(uu,u0)*u1*u1+S.diff(direct,uu).subs(uu,u0)*u2)/2-coeffs[2])
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            radius=mp.mpf('.025')
            def value(z):
                sval=mp.mpf('.001')*mp.cosh(z)+mp.mpf('.0008')*mp.sinh(z)
                uval=1/(1+sval);zz=uval-1
                pol=zz**2+mp.mpf(4)/3*zz**3+zz**4/2-mp.mpf('1e-8')*uval**3
                return 1/(1-2*zz-25*pol)
            for j in [1,4,8]:
                f=lambda angle:value(radius*mp.exp(1j*angle))*mp.exp(-1j*j*angle)
                integral=mp.quad(f,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi*radius**j)
                exact=mp.mpf(str(S.numer(coeffs[j])))/mp.mpf(str(S.denom(coeffs[j])))
                error=abs(integral-exact)
                check('independent_invariant_contour_'+str(dps)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,10))
                records.append({'dps':dps,'order':j,'absolute_error':mp.nstr(error,16)})
    FIXTURES.append({'fixture':'independent_invariant_composition_contours','records':records,
                     'coefficients':[str(v) for v in coeffs],'actual_Born_source_evaluated':False,
                     'fixture_parameters_are_physical_seed':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold local eight jets','input_count':len(ps)}),flush=True)
    try:
        provenance();polynomial_and_invariant_domains();local_and_observation_budgets()
        connection_and_independent_controls();full_tensor_ledger();independent_invariant_contours()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold local eight jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_moving_local_cold_C8_bounded':True,'local_cold_C8_reserve':'1e8',
         'full_physical_covariant_tensors_bounded':True,'full_angular_tower_retained':True,
         'common_finite_part_retained':True,'nonnormal_connection_retained':True,
         'complete_cold_source_C8_bounded':False,'cold_gluing_C8_bounded':False,'actual_domain_C8_bounded':False,
         'actual_metric_tenth_derivative_bounded':False,'actual_finite_metric_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'original_cold_C3_contract_changed':False,
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
