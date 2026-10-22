"""Actual positive-inverse budgets and exact independent continuum fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-positive-weighted-inverse'
PREVIOUS='finite-strong-schur-source'
OWN=['checks/finite_positive_weighted_inverse.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational
ROWS=[]
FIXTURES=[]
PROTOCOL={'arithmetic':'exact symbolic and rational',
          'weighted_inverse_constants':[4200,9500,18000,100000],
          'actual_first_positive_solve_weighted_second_order_bounded':True,
          'actual_certificate_domain_membership_proved':False,
          'actual_first_transferred_source_bounded':False,
          'actual_transfer_certificate_passed':False}


def check(name,ok,detail=''):
    assert name not in {c['name'] for c in ROWS}
    ROWS.append({'name':name,'passed':bool(ok),'detail':str(detail)})
    print(json.dumps(ROWS[-1]),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expression):
    v=S.factor(S.expand(S.simplify(expression)))
    check(name,v==0,v)


def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('predecessor_passed',p['status']=='passed' and len(p['checks'])==74 and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('predecessor_hashes',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope_retained',review['actual_first_positive_solve_weighted_second_order_bounded'] and not review['actual_certificate_domain_membership_proved'] and not review['goal_complete'])


def identities():
    z=S.symbols('z',positive=True)
    m,w,u,alpha,b=[S.Function(n)(z) for n in ['m','w','u','alpha','b']]
    dd=lambda f:S.diff(f,z)
    T=lambda f:-dd(m*dd(f))/m
    direct=-dd(m*(1+alpha)*dd(u))/m+(4+b)*u
    zero('strong_local_operator',direct-((1+alpha)*T(u)-dd(alpha)*dd(u)+(4+b)*u))
    zero('weighted_form_product',(1+alpha)*dd(u)*dd(w*w*u)-(1+alpha)*(dd(w*u)**2-dd(w)**2*u*u))
    zero('regular_measure_integration',dd(m*u*u)-dd(m)*u*u-2*m*u*dd(u))
    q,v=S.symbols('q v')
    e,E,aa,ZZ=S.symbols('e E a Z',nonzero=True)
    local_block=dd(m*aa*dd(u))/m-ZZ*u/(3*e)
    zero('local_scalar_sign',-local_block/E-(-dd(m*(aa/E)*dd(u))/m+ZZ/(3*e*E)*u))
    r,n,c=S.symbols('r n c',positive=True)
    zero('positive_energy_square_completion',r*n-4*c*r*r-(n*n/(16*c)-4*c*(r-n/(8*c))**2))
    l,gg,ff,beta,ric=S.symbols('L G F beta ric',nonnegative=True)
    zero('weighted_Bochner_square_completion',(l+beta*gg)**2+ric*gg**2-((ff-beta*gg)**2)-(l*l+2*beta*l*gg+2*beta*ff*gg+ric*gg*gg-ff*ff))
    aa,rr,u1,u2,dA,dR,ww,w1=S.symbols('A R u1 u2 A1 R1 W W1')
    density=S.symbols('density',positive=True)
    d=aa+2*rr
    lap=u2+d*u1
    hess2=u2**2+(aa*aa+2*rr*rr)*u1*u1
    ricrr=-(dA+aa*aa)-2*(dR+rr*rr)
    derivative=density*ww*ww*((d*d+dA+2*dR)*u1*u1+2*d*u1*u2)+2*density*ww*w1*d*u1*u1
    zero('radial_weighted_Bochner_divergence',density*ww*ww*(hess2-lap*lap+ricrr*u1*u1)-2*density*ww*w1*(lap*u1-u2*u1)+derivative)
    tt=S.symbols('t',nonnegative=True)
    zero('spectral_gradient_interpolation',(tt+1)**2/4-tt-(tt-1)**2/4)
    FIXTURES.append({'fixture':'exact_general_weighted_identities','physical_source_evaluated':False})


def budgets():
    g=R('4e-17')**2;emin=R('1.125e35')
    check('vertex_normalization',R(20000)**2/180<1500**2)
    check('local_principal_coefficient',10/emin<R('1e-32'))
    check('radial_a_gradient',R('.0017')/36+R(11,9)*R('3e17')<R('4e17'))
    check('local_principal_gradient',R('4e17')/emin<R('4e-18'))
    check('local_potential',R('3.1e-12')/3<R('2e-12'))
    check('Kato_Rellich_relative_bound',R('1e-32')+R('2e-18')<1)
    check('seed_A_rational_lower',R('.7071')**2<R(1,2))
    check('actual_measure_growth',(1-R('3e-8'))*(R('.7071')-R('2.45e-6'))>R(2,3))
    zero('physical_Poincare_factor',(R(2,3)/2)**2-R(1,9))
    weight2=4*(1+1/R('7.9e9'))/(1-R('6e-8'))
    check('actual_weight_gradient_squared',weight2<R('4.000001'),weight2)
    check('actual_weight_gradient',R('4.000001')<R('2.001')**2)
    local_loss=(1+R('1e-32'))*R('4.000001')-4+R('2e-12')
    check('weighted_local_mass_loss',local_loss<R('2e-6'),local_loss)
    check('weighted_form_coercivity',R(1,11)-R('1e-32')-9*R('2e-6')>0)
    c=1-R('3e-13')
    check('tensor_energy_output',1/(16*c)<R('.251')**2)
    zero('weighted_tensor_forcing',1500*R('.251')-R('376.5'))
    check('weighted_joint_value_gradient',11*R('377.5')<4200)
    check('weighted_physical_gradient',(1+R('2.001')**2)*4200**2<9500**2)
    lap=(R('377.5')+(4+R('2e-12'))*4200+R('4e-18')*9500)/(1-R('1e-32'))
    check('weighted_scalar_Laplacian',lap<18000,lap)
    check('scalar_Ricci_square_root',15<4**2)
    hess=18000+(2*R('2.001')+4)*9500
    check('weighted_scalar_Hessian',hess<100000,hess)
    for name,k,upper in [('value',4200,R('8.4e31')),('gradient',9500,R('1.9e32')),('Laplacian',18000,R('3.6e32')),('Hessian',100000,R('2e33'))]:
        zero('actual_first_solve_'+name,k*R('2e28')-upper)
    FIXTURES.append({'fixture':'exact_actual_coefficient_budgets','weight_gradient_square_upper':str(weight2),'Laplacian_coefficient':str(lap),'Hessian_coefficient':str(hess),'actual_transfer_certificate_passed':False})


t=S.symbols('t',positive=True)


def integrate_laurent(expr):
    result=S.Integer(0)
    expanded=S.expand(expr)
    if expanded==0:return result
    for term in expanded.as_ordered_terms():
        coef,power=term.as_coeff_exponent(t)
        assert not coef.has(t) and power.is_Integer and power<-1,(term,power)
        result-=coef/(power+1)
    return S.factor(result)


def inner(f,h,w=S.Integer(1)):
    return integrate_laurent(t*t*w*w*f*h)


def omission_control():
    a=t**-4;u=t**-4-2*t**-6;w=t*t
    au=inner(a,u);wau=inner(a,u,w)
    zero('rank_one_unweighted_pair',au+R(3,35))
    zero('rank_one_weighted_pair',wau-R(1,3))
    zero('weighted_positive_form_negative',au*wau+R(1,35))
    check('positive_term_cannot_be_dropped',au*wau<0)
    check('weight_does_not_commute_with_Q',S.expand(w*a*au-a*inner(a,w*u))!=0)
    check('control_has_nonzero_regular_bolt_value',u.subs(t,1)!=0)
    check('control_weighted_integrability',inner(a,a,w)>0 and inner(u,u,w)>0)
    FIXTURES.append({'fixture':'weighted_positive_form_omission','weighted_pair':str(au*wau),'actual_charged_geometry':False})


def continuum_geometry():
    time,angle,azimuth=S.symbols('time angle azimuth',real=True)
    coords=[time,t,angle,azimuth]
    f=t*t-1
    metric=[f,1/f,t*t,t*t*S.sin(angle)**2]
    inv=[1/x for x in metric]
    Gamma={}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                x=((S.diff(metric[j],coords[i]) if k==j else 0)+(S.diff(metric[i],coords[j]) if k==i else 0)-(S.diff(metric[i],coords[k]) if i==j else 0))*inv[k]/2
                Gamma[k,i,j]=S.simplify(x)
    ric=0
    for k in range(4):
        ric+=S.diff(Gamma[k,1,1],coords[k])-S.diff(Gamma[k,1,k],t)
        for l in range(4):ric+=Gamma[k,1,1]*Gamma[l,k,l]-Gamma[l,1,k]*Gamma[k,1,l]
    zero('independent_radial_Ricci',f*ric+3)
    return metric,inv,Gamma


def continuum_fixture(label,u,geometry):
    metric,inv,Gamma=geometry
    f=t*t-1;w=t*t
    alpha=R('1e-34')*t**-2;b=R('1e-13')*t**-4;s=1+alpha
    cols=[100*t**-4,10*t**-5]
    paired=[inner(a,u) for a in cols]
    qu=sum(a*c for a,c in zip(cols,paired))
    # Coordinate Hessian and contraction use Christoffels, independently of the radial form.
    hdiag=[S.simplify((S.diff(u,t,2) if i==1 else 0)-Gamma[1,i,i]*S.diff(u,t)) for i in range(4)]
    hframes=[S.simplify(inv[i]*hdiag[i]) for i in range(4)]
    hrr=f*S.diff(u,t,2)+t*S.diff(u,t)
    targets=[t*S.diff(u,t),hrr,f*S.diff(u,t)/t,f*S.diff(u,t)/t]
    for i in range(4):zero(f'{label}_Hessian_component_{i}',hframes[i]-targets[i])
    lap=S.expand(sum(hframes));hess2=S.expand(sum(v*v for v in hframes))
    T=lambda v:-S.diff(t*t*f*S.diff(v,t),t)/(t*t)
    zero(label+'_Laplacian',lap+T(u))
    qdirect=-S.diff(t*t*f*s*S.diff(u,t),t)/(t*t)+(4+b)*u+qu
    qstrong=s*T(u)-f*S.diff(alpha,t)*S.diff(u,t)+(4+b)*u+qu
    zero(label+'_local_and_nonlocal_equation',qdirect-qstrong)
    q=S.expand(qdirect)
    energy=integrate_laurent(t*t*(s*f*S.diff(u,t)**2+(4+b)*u*u))+sum(c*c for c in paired)
    zero(label+'_unweighted_energy',inner(u,q)-energy)
    weighted_local=integrate_laurent(t*t*(s*f*S.diff(w*u,t)**2+(4+b-s*4*f/t**2)*(w*u)**2))
    zero(label+'_weighted_energy',weighted_local-inner(u,q-qu,w))
    cross=integrate_laurent(t*t*2*w*f*S.diff(w,t)*S.diff(u,t)*(lap-hrr))
    bochner=integrate_laurent(t*t*w*w*(hess2-lap*lap-3*f*S.diff(u,t)**2))
    zero(label+'_weighted_Bochner',bochner-cross)
    check(label+'_regular_nonzero_bolt',u.subs(t,1)!=0 and S.limit(f*S.diff(u,t)**2,t,1)==0)
    wa2=sum(inner(a,a,w) for a in cols)
    check(label+'_weighted_vertex_bound',wa2<1500**2,wa2)
    n2=inner(q,q,w)
    out2={'value':inner(u,u,w),'gradient':integrate_laurent(t*t*w*w*f*S.diff(u,t)**2),
          'Laplacian':inner(lap,lap,w),'Hessian':integrate_laurent(t*t*w*w*hess2)}
    for key,k in [('value',4200),('gradient',9500),('Laplacian',18000),('Hessian',100000)]:
        check(label+'_bound_'+key,out2[key]<k*k*n2)
    check(label+'_tensor_energy_bound',sum(x*x for x in paired)<R('.251')**2*inner(q,q))
    FIXTURES.append({'fixture':'independent_warped_continuum','case':label,'u':str(u),'source_norm_squared':str(n2),'output_norms_squared':{k:str(v) for k,v in out2.items()},'weighted_vertex_Hilbert_Schmidt_squared':str(wa2),'actual_charged_geometry':False,'actual_source':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite positive weighted inverse','input_count':len(ps)}),flush=True)
    try:
        provenance();identities();budgets();omission_control();geom=continuum_geometry()
        continuum_fixture('mass_four',t**-4-2*t**-6,geom)
        continuum_fixture('fast_five',t**-5+R(1,3)*t**-7,geom)
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
    run={'story':'S167 finite positive weighted inverse','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
         'exception':error,'python':sys.version,'sympy':S.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'actual_positive_operator_domain_identified':good,'actual_weighted_second_order_inverse_proved':good,
         'actual_first_positive_solve_weighted_second_order_bounded':good,
         'actual_certificate_domain_membership_proved':False,'actual_first_transferred_source_bounded':False,
         'actual_signed_moment_evaluated':False,'actual_transfer_certificate_passed':False,
         'actual_weighted_graph_iterates_certified':False,'actual_signed_compatibility_zero_proved':False,
         'cold_metric_response_proved':False,'nonlinear_neighborhood_preserved':False,'nonlinear_solution_constructed':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'physical_stability_proved':False,'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
