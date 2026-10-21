"""Focused algebra and independent fixtures for S168 Ward recovery.

The continuum/domain theorem is an ordinary proof, not formal certification.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'spherical-retarded-constraint-recovery'
OWN = ['checks/spherical_retarded_constraint_recovery.py', 'notes/'+STEM+'.md',
       'evidence/spherical-constraint-source-review.json']
PREVIOUS = 'curvature-free-all-frequency-bound'
ROWS, NUMERICS = [], []
PROTOCOL = {
    'scope': 'Exact sourced Ward block, conditional constraint recovery and off-shell contacts',
    'frequencies': ['0.5+0.7i', '1+2i', '3+0.25i'],
    'boundary_reflections': ['0', '1', '-1', '0.6+0.3i'],
    'quadrature_precision': 40,
    'quadrature_relative_tolerance': '1e-24',
    'points': ['-0.2', '-1', '-3'],
    'sources': 'f_plus=(1+0.4i)exp(2x), f_minus=(0.7-0.2i)exp(3x)',
    'fixtures_are_metric_spectral_data': False,
    'physical_boundary_selected': False,
    'current_response_evaluated': False,
    'formal_verification': False,
}

def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    row = {'name': name, 'passed': bool(ok), 'detail': str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)

def zero(name, expr):
    result = S.simplify(S.sympify(expr).doit())
    check(name, result == 0, result)

def prior():
    return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]

def paths():
    return sorted(set(OWN) | set(prior()['input_sha256_after']) |
                  {'evidence/'+PREVIOUS+'.json'})

def provenance():
    p = prior()
    check('predecessor_passed', p['status']=='passed' and len(p['checks'])==370
          and all(r['passed'] for r in p['checks']))
    check('predecessor_immutable_inputs', p['input_sha256_before']==p['input_sha256_after'])
    mismatches = [f for f,h in p['input_sha256_after'].items()
                  if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_predecessor_input_hashes', not mismatches, mismatches)
    check('charged_Ward_derivation_pinned',
          'notes/charged-conformal-response-and-constraints.md' in p['input_sha256_after'])
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope', len(review['primary_sources'])==1
          and not review['actual_metric_causal_inverse_proved']
          and not review['actual_physical_boundary_selected']
          and not review['actual_current_response_evaluated'])

def metric_checks():
    t,x,th,ph,z=S.symbols('t x theta phi z', real=True)
    coords=[t,x,th,ph]
    n,l=[S.Function(c)(x) for c in ['n','l']]
    p,v=[S.Function(c)(t,x) for c in ['p','v']]
    et,er,eo,k=[S.Function(c)(t,x) for c in ['et','er','eo','k']]
    bt,br,bo=[S.Function(c)(x) for c in ['bt','br','bo']]
    U=S.exp(2*n+2*l)
    metric=S.diag(-S.exp(2*n),S.exp(2*n),S.exp(2*l),S.exp(2*l)*S.sin(th)**2)
    def christoffel(g):
        inv=g.inv()
        return [[[S.simplify(sum(inv[a,d]*(S.diff(g[d,c],coords[b])+
                     S.diff(g[d,b],coords[c])-S.diff(g[b,c],coords[d]))
                     for d in range(4))/2) for c in range(4)]
                     for b in range(4)] for a in range(4)]
    gamma=christoffel(metric)
    def div(m,gam,b):
        return S.simplify(sum(S.diff(m[a,b],coords[a])+
                   sum(gam[a][a][c]*m[c,b]-gam[c][a][b]*m[a,c]
                       for c in range(4)) for a in range(4)))
    em=S.diag(et,er,eo,eo);em[0,1]=-k;em[1,0]=k
    Dt=S.diff(et,t)+S.diff(k,x)+2*S.diff(n+l,x)*k
    Dx=-S.diff(k,t)+S.diff(er,x)+S.diff(n,x)*(er-et)+2*S.diff(l,x)*(er-eo)
    zero('four_dimensional_time_divergence',div(em,gamma,0)-Dt)
    zero('four_dimensional_radial_divergence',div(em,gamma,1)-Dx)
    zero('four_dimensional_theta_divergence',div(em,gamma,2))
    zero('four_dimensional_phi_divergence',div(em,gamma,3))

    gm=S.diag(metric[0,0]*S.exp(2*z*(p+v)),metric[1,1]*S.exp(2*z*(p+v)),
              metric[2,2]*S.exp(2*z*p),metric[3,3]*S.exp(2*z*p))
    gammaz=christoffel(gm)
    Wt=(br-bt)*S.diff(p+v,t)+2*(bo-bt)*S.diff(p,t)
    Wx=(bt-br)*S.diff(p+v,x)+2*(bo-br)*S.diff(p,x)
    bm=S.diag(bt,br,bo,bo)
    for i,(name,W) in enumerate([('time',Wt),('radial',Wx)]):
        direct=S.diff(div(bm+z*em,gammaz,i),z).subs(z,0)
        zero('off_shell_connection_'+name,direct-div(em,gamma,i)+W)
        bvals=[bt,br,bo,bo]
        connection=sum((bvals[a]-bvals[i])*S.diff(gammaz[a][a][i],z).subs(z,0)
                       for a in range(4))
        zero('connection_variation_sum_'+name,connection-W)
    jets={bt:1,br:2,bo:4,S.diff(p,t):3,S.diff(v,t):5}
    check('off_shell_forcing_is_not_zero',S.simplify(Wt.subs(jets))==26)

    tau,o,Ec,Jc=[S.Function(c)(t,x) for c in ['tau','o','Ec','Jc']]
    replace={et:o-Ec/U,er:Ec/U,eo:(tau-o)/2,k:Jc/U}
    Fc=U*(S.diff(o,t)-Wt)
    Gc=-U*(Wx+S.diff(n-l,x)*o+S.diff(l,x)*tau)
    zero('weighted_time_transport',U*(Dt-Wt).subs(replace).doit()
         -(-S.diff(Ec,t)+S.diff(Jc,x)+Fc))
    zero('weighted_radial_transport',U*(Dx-Wx).subs(replace).doit()
         -(-S.diff(Jc,t)+S.diff(Ec,x)+Gc))
    Ft,Gt=S.symbols('Fc Gc')
    evolution={S.diff(Ec,t):S.diff(Jc,x)+Ft,S.diff(Jc,t):S.diff(Ec,x)+Gt}
    zero('incoming_characteristic', (S.diff(Ec+Jc,t)-S.diff(Ec+Jc,x)).subs(evolution)-Ft-Gt)
    zero('outgoing_characteristic', (S.diff(Ec-Jc,t)+S.diff(Ec-Jc,x)).subs(evolution)-Ft+Gt)

    d=4*p+2*v;tau0=bt+br+2*bo;o0=bt+br
    zero('volume_contact',S.diff(S.exp(2*(n+z*(p+v))+2*(l+z*p)),z).subs(z,0)-U*d)
    zero('trace_gradient_contact',S.diff(U*S.exp(z*d)*(tau0+z*tau),z).subs(z,0)/U-tau-d*tau0)
    zero('orbit_gradient_contact',S.diff(U*S.exp(z*d)*(o0+z*o),z).subs(z,0)/U-o-d*o0)
    rowa,rowb=[S.Function(c)(t,x) for c in ['rowa','rowb']]
    contact={tau:rowa-d*tau0,o:rowb-d*o0}
    zero('coordinate_Hessian_time_source',Fc.subs(contact).doit()-U*(S.diff(rowb,t)-o0*S.diff(d,t)-Wt))
    zero('coordinate_Hessian_radial_source',Gc.subs(contact).doit()
         +U*(Wx+S.diff(n-l,x)*(rowb-d*o0)+S.diff(l,x)*(rowa-d*tau0)))

    pulse=S.Function('pulse')(t+x)
    residual={et:-pulse/(2*U),er:pulse/(2*U),eo:0,k:pulse/(2*U)}
    zero('boundary_injection_time_Ward',div(em,gamma,0).subs(residual).doit())
    zero('boundary_injection_radial_Ward',div(em,gamma,1).subs(residual).doit())
    inj=metric*em.subs(residual)
    target=S.zeros(4)
    for a in [0,1]:
        for b in [0,1]:target[a,b]=pulse*S.exp(-2*l)/2
    check('ingoing_covariant_tensor',all(S.simplify(z)==0 for z in inj-target))
    zero('boundary_injection_trace',sum(em[i,i] for i in range(4)).subs(residual))
    zero('boundary_injection_orbit_trace',(et+er).subs(residual))

    q,gc=S.symbols('q gc', nonzero=True, real=True)
    f=S.Function('f')(t,x)
    form=S.zeros(4);form[0,1]=z*f;form[1,0]=-z*f
    form[2,3]=q*S.sin(th)/2;form[3,2]=-form[2,3]
    inv=metric.inv();upper=inv*form*inv;density=U*S.sin(th)
    Q=S.exp(2*l-2*n)*f/gc**2
    maxt=sum(S.diff(density*upper[a,0]/gc**2,coords[a]) for a in range(4))
    maxx=sum(S.diff(density*upper[a,1]/gc**2,coords[a]) for a in range(4))
    zero('electric_Gauss_density',S.diff(maxt,z)/S.sin(th)-S.diff(Q,x))
    zero('electric_flux_evolution',S.diff(maxx,z)/S.sin(th)+S.diff(Q,t))
    zero('electric_current_integrability',S.diff(S.diff(Q,x),t)-S.diff(S.diff(Q,t),x))
    invz=gm.inv();magdensity=invz*form.subs(z,0)*invz
    check('metric_variation_has_no_orbit_magnetic_Maxwell_density',
          all(magdensity[a,b]==0 for a in range(4) for b in [0,1]))
    square=sum(form[a,b]*upper[a,b] for a in range(4) for b in range(4))
    stress=form*inv*form.T-metric*square/4
    check('electric_magnetic_stress_cross_term_zero',
          all(S.simplify(S.diff(y,z).subs(z,0))==0 for y in stress))
    jt,jx=S.symbols('jt jx')
    force=form.subs(z,0)*S.Matrix([jt,jx,0,0])
    check('spherical_magnetic_orbit_force_zero',force==S.zeros(4,1))

def spectral_algebra():
    x=S.symbols('x',real=True)
    s,a,b,rho,A,B,defect=S.symbols('s a b rho A B defect')
    cm=B*S.exp(b*x)/(s+b)
    cp=rho*B*S.exp(s*x)/(s+b)+A*(S.exp(a*x)-S.exp(s*x))/(s-a)
    zero('outgoing_resolvent_equation',s*cm+S.diff(cm,x)-B*S.exp(b*x))
    zero('incoming_resolvent_equation',s*cp-S.diff(cp,x)-A*S.exp(a*x))
    zero('homogeneous_boundary_relation',cp.subs(x,0)-rho*cm.subs(x,0))
    zero('boundary_defect_homogeneous_equation',s*defect*S.exp(s*x)-S.diff(defect*S.exp(s*x),x))
    cpr,cpi,cmr,cmi=S.symbols('cpr cpi cmr cmi',real=True)
    Cplus=cpr+S.I*cpi;Cminus=cmr+S.I*cmi
    Ec=(Cplus+Cminus)/2;Jc=(Cplus-Cminus)/2
    zero('characteristic_norm',S.expand_complex(S.conjugate(Ec)*Ec+S.conjugate(Jc)*Jc
         -(S.conjugate(Cplus)*Cplus+S.conjugate(Cminus)*Cminus)/2))
    zero('boundary_energy_flux',S.re(S.conjugate(Ec)*Jc)
         -(cpr**2+cpi**2-cmr**2-cmi**2)/4)
    gamma=S.symbols('gamma',positive=True,real=True)
    zero('boundary_defect_norm_constant',S.integrate(S.exp(2*gamma*x)/2,(x,-S.oo,0))-1/(4*gamma))

def spectral_quadrature():
    mp.mp.dps=PROTOCOL['quadrature_precision']
    A=mp.mpc(1,'.4');B=mp.mpc('.7','-.2');a=mp.mpf(2);b=mp.mpf(3)
    cases=[mp.mpc('.5','.7'),mp.mpc(1,2),mp.mpc(3,'.25')]
    rhos=[mp.mpc(0),mp.mpc(1),mp.mpc(-1),mp.mpc('.6','.3')]
    def integral_norm(coefficients):
        return mp.re(sum(c*mp.conj(d)/(lam+mp.conj(mu))
                  for c,lam in coefficients for d,mu in coefficients))
    source_norm2=(abs(A)**2/(2*a)+abs(B)**2/(2*b))/2
    for i,s in enumerate(cases):
        for j,rho in enumerate(rhos):
            def cm(x):return B*mp.exp(b*x)/(s+b)
            def cp(x):return rho*B*mp.exp(s*x)/(s+b)+A*(mp.exp(a*x)-mp.exp(s*x))/(s-a)
            errs=[]
            for px in PROTOCOL['points']:
                x=mp.mpf(px)
                minus=mp.quad(lambda y:mp.exp(-s*(x-y))*B*mp.exp(b*y),[-mp.inf,x])
                minus0=mp.quad(lambda y:mp.exp(s*y)*B*mp.exp(b*y),[-mp.inf,0])
                plus=rho*mp.exp(s*x)*minus0+mp.quad(lambda y:mp.exp(s*(x-y))*A*mp.exp(a*y),[x,0])
                errs.extend([abs(plus-cp(x)),abs(minus-cm(x))])
            check(f'complex_characteristic_quadrature_{i}_{j}',max(errs)<mp.mpf('1e-24'),mp.nstr(max(errs),8))
            norm2=mp.quad(lambda x:(abs(cp(x))**2+abs(cm(x))**2)/2,[-mp.inf,0])
            analytic=(integral_norm([(rho*B/(s+b)-A/(s-a),s),(A/(s-a),a)])
                      +integral_norm([(B/(s+b),b)]))/2
            check(f'independent_L2_integral_{i}_{j}',abs(norm2-analytic)<mp.mpf('1e-24'),mp.nstr(abs(norm2-analytic),8))
            ratio=mp.re(s)*mp.sqrt(norm2/source_norm2)
            check(f'Ward_resolvent_bound_{i}_{j}',ratio<=1,mp.nstr(ratio,16))
            defect=mp.mpc(1,'.3')
            defect_norm2=mp.quad(lambda x:abs(defect*mp.exp(s*x))**2/2,[-mp.inf,0])
            check(f'boundary_defect_norm_{i}_{j}',abs(defect_norm2-abs(defect)**2/(4*mp.re(s)))<mp.mpf('1e-24'))
            NUMERICS.append({'s':str(s),'rho':str(rho),'maximum_kernel_error':mp.nstr(max(errs),25),
                             'norm_squared':mp.nstr(norm2,25),'resolvent_ratio':mp.nstr(ratio,25)})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n')
    sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 spherical constraint recovery','input_count':len(ps)}),flush=True)
    try:
        provenance();metric_checks();spectral_algebra();spectral_quadrature()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S168 spherical constraint recovery','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'checks':ROWS,'exception':error,'numerics':NUMERICS,
         'conditional_Ward_recovery_derived':good,'Ward_right_half_plane_resolvent_derived':good,
         'off_shell_connection_and_density_contacts_retained':good,'electric_equations_derived':good,
         'actual_physical_boundary_selected':False,'actual_current_response_evaluated':False,
         'actual_metric_causal_inverse_proved':False,'actual_metric_poles_excluded':False,
         'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
         'retained_fast_factor_removed':False,'order_reduction_adopted':False,
         'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
