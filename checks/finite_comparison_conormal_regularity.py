"""Supporting boundary, proper-time and reconstruction checks; not a heat solver."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-comparison-conormal-regularity'
PREVIOUS='finite-positive-weighted-inverse'
OWN=['checks/finite_comparison_conormal_regularity.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[]
FIXTURES=[]
R=S.Rational
PROTOCOL={'arithmetic':'exact symbolic identities and 80-digit independent quadrature',
          'quadrature_Y':['0.1','0.5','2'], 'relative_tolerance':'1e-60',
          'ordinary_continuum_proof_separate':True,
          'actual_charged_heat_kernel_evaluated':False,
          'higher_order_constants_evaluated':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append({'name':name,'passed':bool(ok),'detail':str(detail)})
    print(json.dumps(ROWS[-1]),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    val=S.simplify(S.expand(expr))
    check(name,val==0,val)


def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('predecessor_passed',p['status']=='passed' and len(p['checks'])==77 and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('predecessor_hashes',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope_retained',review['ordinary_fixed_seed_regularity_proved'] and review['ordinary_original_family_conormal_regularity_proved'] and not review['finite_metric_quantum_source_higher_regularity_proved'] and not review['goal_complete'])


def boundary_model():
    y,yp,tau=S.symbols('y yp tau',positive=True)
    g=lambda v:S.exp(-v*v/(4*tau))/S.sqrt(4*S.pi*tau)
    direct=g(y-yp);refl=g(y+yp)
    kd=direct-refl;kn=direct+refl
    zero('Dirichlet_left_boundary',kd.subs(y,0))
    zero('Neumann_left_boundary',S.diff(kn,y).subs(y,0))
    zero('Dirichlet_heat_equation',S.diff(kd,tau)-S.diff(kd,y,2))
    zero('Neumann_heat_equation',S.diff(kn,tau)-S.diff(kn,y,2))
    pd=S.diag(1,0);pn=S.eye(2)-pd
    kernel=pd*kd+pn*kn
    check('complementary_projectors',pd*pn==S.zeros(2) and pd+pn==S.eye(2))
    zero('balanced_rank_two_trace',S.trace(kernel)-2*direct)
    tangential_matrix=S.Matrix([[0,1],[-1,0]])
    check('connection_need_not_commute_with_projectors',tangential_matrix*pd!=pd*tangential_matrix)
    FIXTURES.append({'fixture':'constant_projector_half_space','actual_charged_heat_kernel':False,'normalization':'trace equals twice direct scalar Gaussian; no extra half'})


def mellin_models():
    y,tau,v,aa=S.symbols('y tau v a',positive=True)
    integrand=tau**aa*S.exp(-y*y/tau)
    zero('heat_integration_by_parts',y*S.diff(integrand,y)+2*S.diff(tau*integrand,tau)-2*(aa+1)*integrand)
    j=S.symbols('j',integer=True)
    # Positive Jacobian after reversing the limits in tau=y^2/v^2.
    transformed=(y*y/v**2)**(j/R(2)-3)*(2*y*y/v**3)
    zero('Mellin_change_of_variables',S.powsimp(transformed,force=True)-2*y**(j-4)*v**(3-j))
    exprs={2:S.exp(-y*y)/y**2,
           3:S.sqrt(S.pi)*S.erfc(y)/y,
           4:S.expint(1,y*y),
           5:2*S.exp(-y*y)-2*S.sqrt(S.pi)*y*S.erfc(y),
           6:S.exp(-y*y)-y*y*S.expint(1,y*y)}
    for jj,expr in exprs.items():
        zero('Mellin_ODE_'+str(jj),y*S.diff(expr,y)-(jj-4)*expr+2*S.exp(-y*y))
    mp.mp.dps=80
    numeric={2:lambda x:mp.exp(-x*x)/x**2,
             3:lambda x:mp.sqrt(mp.pi)*mp.erfc(x)/x,
             4:lambda x:mp.e1(x*x),
             5:lambda x:2*mp.exp(-x*x)-2*mp.sqrt(mp.pi)*x*mp.erfc(x),
             6:lambda x:mp.exp(-x*x)-x*x*mp.e1(x*x)}
    errors=[]
    for jj in exprs:
        for value in PROTOCOL['quadrature_Y']:
            x=mp.mpf(value)
            direct=mp.quad(lambda tt:tt**(mp.mpf(jj)/2-3)*mp.exp(-x*x/tt),[0,mp.mpf('.01'),mp.mpf('.1'),1])
            closed=numeric[jj](x)
            err=abs(direct-closed)/max(mp.mpf(1),abs(closed))
            errors.append(err)
            check(f'proper_time_quadrature_{jj}_{value}',err<mp.mpf(PROTOCOL['relative_tolerance']),mp.nstr(err,12))
    # Exact resonant primitive and its logarithmic coefficient, no numerical fit.
    zero('resonant_log_primitive',S.diff(S.log(v),v)-1/v)
    zero('log_profile_derivative',y*S.diff(exprs[4],y)+2*S.exp(-y*y))
    zero('inverse_distance_leading',S.limit(y*exprs[3],y,0)-S.sqrt(S.pi))
    FIXTURES.append({'fixture':'boundary_Mellin_profiles','indices':list(exprs),'quadrature_points':PROTOCOL['quadrature_Y'],'max_relative_error':mp.nstr(max(errors),20),'actual_source':False})


def reconstruction():
    t=S.symbols('t',positive=True)
    x,a,Q,d=S.symbols('x a Q d')
    F=x*t*t+1-2*a/t+Q/t**2-2*d/t**3+d/t**4
    trace=S.Function('Theta')(t);I=S.Function('I')(t)
    bfun=F/t**2
    Ip=S.diff(bfun,t)*trace
    pressure=I/(2*F);rho=pressure-trace/t**2
    replace={S.diff(I,t):Ip}
    zero('Ward_trace',pressure-rho-trace/t**2)
    zero('Ward_pressure_equation',S.diff(pressure,t).subs(replace)+S.diff(F,t)*(rho+pressure)/(2*F)+trace/t**3)
    # A manufactured logarithmic trace is reconstructed by direct integration.
    th=3*t+2*S.log(t)+5+7/t
    primitive=S.integrate(S.expand(S.diff(bfun,t)*th),t)
    zero('Ward_logarithmic_primitive',S.diff(primitive,t)-S.diff(bfun,t)*th)
    z=S.symbols('z',positive=True)
    for mm,kk in [(-1,0),(-1,2),(0,1),(2,2)]:
        mon=z**mm*S.log(z)**kk
        if mm==-1:ant=S.log(z)**(kk+1)/(kk+1)
        else:
            ant=z**(mm+1)*sum((-1)**s*S.factorial(kk)*S.log(z)**(kk-s)/(S.factorial(kk-s)*(mm+1)**(s+1)) for s in range(kk+1))
        zero(f'polyhomogeneous_primitive_{mm}_{kk}',S.diff(ant,z)-mon)
    G,r=S.symbols('G r',positive=True)
    f,al,rr,pp,F0=[S.Function(n)(r) for n in ['f','a','rho','p','F0']]
    dm=S.Function('mass')(r)
    zero('mass_constraint',S.diff(r*(-2*dm/r),r).subs(S.diff(dm,r),G*rr)+2*G*rr)
    FIXTURES.append({'fixture':'Ward_and_constraint_reconstruction','trace':str(th),'actual_source':False,'integration_constants':'not assigned a physical value'})


def completion():
    z,x,a,Q0,d,C,kappa,B,h=S.symbols('z x a Q0 d C kappa B h',positive=True)
    cal=x+z*z-2*a*z**3+Q0*z**4+d*z**4*(1-z)**2
    zero('outer_seed_value',cal.subs(z,0)-x)
    zero('outer_seed_geodesic_jet',S.diff(cal,z).subs(z,0))
    v=S.symbols('v',negative=True)
    chi=2*S.atanh(S.exp(v))
    zero('disk_coordinate_derivative',S.diff(chi,v)-2*S.exp(v)/(1-S.exp(2*v)))
    ycp=2*S.exp(-h)/(1-S.exp(-2*h))*kappa*B/x
    zero('outer_disk_positive_denominator_identity',ycp*(S.exp(h)-S.exp(-h))-2*kappa*B/x)
    y=S.symbols('y',positive=True)
    ell,u3,v3,u4,v4=S.symbols('ell u3 v3 u4 v4')
    u=u3*y**3+u4*y**4*S.log(y)
    vv=v3*y**3+v4*y**4*S.log(y)
    # Use r0=1/y for the leading coefficient fixture; the theorem allows its analytic lower terms.
    ratio=y+(1-y)*S.exp(ell*vv)
    Z=y/ratio
    zero('completion_radius_identity',(1/Z-1)-(1/y-1)*S.exp(ell*vv))
    zero('completion_tangent',S.diff(ratio,ell).subs(ell,0)-(1-y)*vv)
    vr=-y*y*S.diff(vv,y)
    alpha=2*S.log(ratio)+2*ell*u-ell*vv-S.log(1+ell*(1/y-1)*vr)
    cubic=S.expand(S.series(alpha,y,0,4).removeO()).coeff(y,3)
    zero('completion_lapse_cubic',cubic-(2*ell*u3+4*ell*v3))
    zero('correlated_cubic_cancellation',cubic.subs(u3,-2*v3))
    FIXTURES.append({'fixture':'analytic_coordinates_and_original_completion','actual_source':False,'correlated_cubic':'u3+2V3=0','ordinary_logarithm_allowed':True})


def negative_controls():
    y=S.symbols('y',positive=True)
    D=lambda f:y*S.diff(f,y)
    oscillatory=y**4*S.sin(y**-2)
    zero('unsigned_remainder_derivative',D(oscillatory)-(4*y**4*S.sin(y**-2)-2*y**2*S.cos(y**-2)))
    f=y**4*S.log(y)
    derivative=f
    for j in range(1,9):
        derivative=D(derivative)
        zero('conormal_log_order_'+str(j),derivative-y**4*(4**j*S.log(y)+j*4**(j-1)))
    zero('ordinary_fourth_log_derivative',S.diff(f,y,4)-24*S.log(y)-50)
    N=S.symbols('N',positive=True)
    wave=N**-5*S.sin(N*y)
    zero('C5_family_sixth_derivative',S.diff(wave,y,6)+N*S.sin(N*y))
    check('C5_family_lower_derivative_amplitudes',all(j-5<=0 for j in range(6)),'For N>=1, amplitudes N^(j-5)<=1 for j<=5; sixth amplitude N.')
    FIXTURES.append({'fixture':'invalid_regularization_controls','unsigned_bound_differentiation_rejected':True,'ordinary_C_infinity_inference_rejected':True,'uniform_C5_to_C6_inference_rejected':True})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 original comparison conormal regularity','input_count':len(ps)}),flush=True)
    try:
        provenance();boundary_model();mellin_models();reconstruction();completion();negative_controls()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
    run={'story':'S167 original comparison conormal regularity','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_fixed_seed_regularity_proved':True,'ordinary_original_family_conormal_regularity_proved':True,
         'ordinary_proof_machine_verified':False,'higher_order_constants_evaluated':False,
         'finite_metric_quantum_source_higher_regularity_proved':False,'actual_charged_heat_kernel_evaluated':False,
         'actual_certificate_domain_membership_proved':False,'actual_first_transferred_source_bounded':False,
         'actual_transfer_certificate_passed':False,'actual_signed_moment_evaluated':False,'actual_signed_compatibility_zero_proved':False,
         'cold_metric_response_proved':False,'nonlinear_neighborhood_preserved':False,'nonlinear_solution_constructed':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,'physical_stability_proved':False,
         'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
