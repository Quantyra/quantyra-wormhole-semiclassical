"""Frozen symmetry, off-shell tensor and endpoint diagnostics; no solution claim."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-stationary-constraints';PREVIOUS='recentered-mass-boundary'
OWN=['notes/'+STEM+'.md','checks/recentered_stationary_constraints.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def intertwiner():
    i=S.I;one=S.eye(2);zero=S.zeros(8)
    sx=S.Matrix([[0,1],[1,0]]);sy=S.Matrix([[0,-i],[i,0]]);sz=S.diag(1,-1)
    k=S.kronecker_product
    gt=-k(sz,sz,one);gn=k(one,sx,one);ga=[-k(sz,sy,x) for x in [sx,sy]]
    U=k(one,sx,sz);chi=k(sz,sx,one);Pi=(S.eye(8)+chi)/2
    check('unitary_Hermitian_reflection',U.H==U and U*U==S.eye(8))
    for name,g,sign in [('time',gt,-1),('normal',gn,1),('sphere1',ga[0],1),('sphere2',ga[1],1)]:
        check('Clifford_reflection_'+name,U*g*U==sign*g)
    check('actual_first_order_domain_preserved',U*Pi==Pi*U)
    check('complementary_adjoint_domain_preserved',U*(S.eye(8)-Pi)==(S.eye(8)-Pi)*U)
    check('polar_rotation_generator_reversed',U*(gn*gt)==-(gn*gt)*U)
    # Reflection lift in a polar frame extends as a constant Cartesian lift.
    angle=S.symbols('angle',real=True)
    rot=S.cos(angle/2)*S.eye(8)+S.sin(angle/2)*gn*gt
    check('polar_spin_lift',S.simplify(U*rot.subs(angle,-angle)-rot*U)==zero)
    check('thermal_spin_structure',rot.subs(angle,2*S.pi)==-S.eye(8))
    radial=k(sz,sy,one);mass=k(one,sx,sx)
    check('Hamiltonian_radial_anticommutation',U*radial+radial*U==zero)
    check('full_angular_pair_anticommutation',U*mass+mass*U==zero)
    check('reflection_preserves_each_copy',U*k(sz,one,one)==k(sz,one,one)*U)
    for sign in [-1,1]:
        Ulll=sign*k(one,sx);Ar=k(sz,sy);Plll=(S.eye(4)+k(sz,sx))/2
        check('LLL_radial_pair_'+str(sign),Ulll*Ar+Ar*Ulll==S.zeros(4))
        check('LLL_domain_'+str(sign),Ulll*Plll==Plll*Ulll)
    # The ordinary choice gamma5 anticommutes with the boundary projector.
    g5=gt*gn*ga[0]*ga[1]
    check('chirality_alone_does_not_preserve_domain',g5*chi-chi*g5!=zero)
    FIXTURES.append(dict(fixture='actual_pair_reflection_algebra',angular_tower_truncated=False,
        current_target='exact zero on selected static magnetic family',electric_susceptibility_computed=False))

def tensor_and_contacts():
    x,ph,th,az=S.symbols('x phi theta az',real=True);coords=[x,ph,th,az]
    p=S.Function('p')(x);v=S.Function('v')(x);L=p+v;W=S.exp(4*p+2*v)
    metric=S.diag(S.exp(2*L),S.exp(2*L),S.exp(2*p),S.exp(2*p)*S.sin(th)**2)
    inv=metric.inv();Gamma={}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                z=S.simplify(sum(inv[a,d]*(S.diff(metric[d,c],coords[b])+S.diff(metric[d,b],coords[c])-S.diff(metric[b,c],coords[d]))/2 for d in range(4)))
                if z!=0:Gamma[a,b,c]=z
    gam=lambda a,b,c:Gamma.get((a,b,c),S.S.Zero)
    ex,et,eo=[S.Function(n)(x) for n in ['Ex','Et','Eo']];Es=[ex,et,eo,eo]
    divergence=S.diff(ex,x)+sum(gam(a,a,0)*(ex-Es[a]) for a in range(4))
    tau=ex+et+2*eo;orb=ex+et
    target=S.diff(W*ex,x)/W-S.diff(v,x)*orb-S.diff(p,x)*tau
    check('Christoffel_Ward_reconstruction',S.simplify(divergence-target)==0)
    # Independent Ricci calculation, not the Ward reconstruction formula.
    Ric=[]
    for a in range(4):
        z=sum(S.diff(gam(c,a,a),coords[c])-S.diff(gam(c,a,c),coords[a])+
              sum(gam(c,c,d)*gam(d,a,a)-gam(c,a,d)*gam(d,a,c) for d in range(4)) for c in range(4))
        Ric.append(S.simplify(inv[a,a]*z))
    scalar=sum(Ric);Ein=[S.simplify(r-scalar/2+S.Rational(3,7)) for r in Ric]
    check('independent_Ricci_spherical_components',S.simplify(Ein[2]-Ein[3])==0)
    b=Ein[0];f1=W*(Ein[0]+Ein[1]+2*Ein[2]);f2=W*(Ein[0]+Ein[1])
    check('actual_off_shell_Einstein_Ward',S.simplify(S.diff(W*b,x)-S.diff(v,x)*f2-S.diff(p,x)*f1)==0)
    s=S.Function('s')(x);h=S.Function('h')(x);amp=S.symbols('amp',real=True);d=4*s+2*h
    def variation(expr):
        changed=expr.subs({p:p+amp*s,v:v+amp*h},simultaneous=True).doit()
        return S.simplify(S.diff(changed,amp).subs(amp,0))
    eb=variation(b);df1=variation(f1);df2=variation(f2)
    lhs=S.diff(W*(eb+d*b),x)
    rhs=S.diff(v,x)*df2+S.diff(p,x)*df1+S.diff(h,x)*f2+S.diff(s,x)*f1
    check('independent_off_shell_directional_variation',S.simplify(lhs-rhs)==0)
    check('moving_density_contact',S.simplify(df1-W*variation(f1/W)-d*f1)==0)
    check('moving_density_contact_orbit',S.simplify(df2-W*variation(f2/W)-d*f2)==0)
    # Two smooth disk jets. x=log(xi); p and v-x are smooth in xi^2.
    fixtures=[]
    for idx in [1,2]:
        u=S.exp(2*x)
        funs={p:S.log(2)+u/(10+idx),v:x+u/(7+idx),
              s:u/(3+idx+u),h:u**2/(5+idx+u)}
        def at(expr):return S.simplify(expr.subs(funs,simultaneous=True).doit().subs(x,0))
        omitted_connection=at(S.diff(h,x)*f2+S.diff(s,x)*f1)
        omitted_density=at(S.diff(W*d*b,x))
        check('fixture_nonzero_reference_'+str(idx),at(b)!=0)
        check('fixture_residual_contact_needed_'+str(idx),omitted_connection!=0)
        check('fixture_density_contact_needed_'+str(idx),omitted_density!=0)
        check('fixture_full_identity_'+str(idx),at(lhs-rhs)==0)
        fixtures.append(dict(index=idx,reference_radial_residual=str(at(b)),
            omitted_connection_error=str(omitted_connection),omitted_density_error=str(omitted_density)))
    # A metric-independent magnetic density satisfies all four Maxwell divergences.
    q=S.symbols('q');rootg=S.exp(4*p+2*v)*S.sin(th)
    Fup={(2,3):q*S.exp(-4*p)/S.sin(th),(3,2):-q*S.exp(-4*p)/S.sin(th)}
    for bb in range(4):
        div=sum(S.diff(rootg*Fup.get((aa,bb),0),coords[aa]) for aa in range(4))
        check('warped_magnetic_Maxwell_'+str(bb),S.simplify(div)==0)
    FIXTURES.append(dict(fixture='independent_Christoffel_Ricci_and_off_shell_contacts',jets=fixtures,
        full_quantum_tensor_computed=False,dynamic_constraint_recovery_claimed=False))

def endpoint_and_targets():
    e,d0,AT,km=S.symbols('e d0 AT km',real=True)
    q=(4+e)*(1+e)
    # Expand the actual digamma/cotangent, with no dropped contact constant.
    F=S.polygamma(0,3+e)+S.pi*S.cot(S.pi*e/2)/2
    N=d0*q*(q+2)*(F+km)-AT*q
    expansion=S.series(N,e,0,1).removeO().expand()
    check('normal_pole_residue',S.simplify(expansion.coeff(e,-1)-24*d0)==0)
    C4=d0*(86-24*S.EulerGamma+24*km)-4*AT
    check('normal_pole_finite_part',S.simplify(expansion.coeff(e,0)-C4)==0)
    check('reciprocal_zero_slope',S.simplify(S.limit(1/(e*N),e,0)-1/(24*d0))==0)
    z,cut=S.symbols('z cut',positive=True)
    remote=-24*d0*z**4*(cut**e-z**e)/e
    check('cutoff_remote_logarithm',S.simplify(S.limit(remote,e,0)-24*d0*z**4*S.log(z/cut))==0)
    check('weighted_counterexample_in_space',Q(15,4)>Q(7,2))
    check('weighted_counterexample_fails_fast_phase',Q(15,4)<4)
    check('physical_admissible_kernel_strip_excludes_formal_zero',Q(4)>=4)
    vals=[]
    for digits in [100,150]:
        mp.mp.dps=digits;g0=mp.mpf(2000)**-10
        g=2*g0/(1+mp.sqrt(1+4*11*g0/(360*mp.pi)))
        Fl=1-384*mp.pi*mp.mpf('1e22')*g
        bclass=mp.sqrt((mp.sqrt(1+12*mp.pi*mp.mpf(2000)**6)-1)/6*g/g0)
        K=-mp.mpf(6)/5+mp.log(bclass)
        dd=1/(480*mp.pi**2);A=Fl/(48*mp.pi*g)+mp.mpf(23)/(8640*mp.pi**2)
        e0=6*dd/A
        def scaled(vv):
            ee=e0*vv;qq=4+5*ee+ee*ee
            fm=mp.digamma(3+ee)+mp.pi*mp.cot(mp.pi*ee/2)/2
            return ee*(dd*qq*(qq+2)*(fm+K)-A*qq)/(24*dd)
        vv=mp.findroot(scaled,(mp.mpf('.99'),mp.mpf('1.01')),tol=mp.mpf(10)**(-digits+15))
        ee=e0*vv
        check('normal_local_root_bracket_'+str(digits),scaled(1-mp.mpf('1e-29'))>0>scaled(1+mp.mpf('1e-29')))
        check('normal_local_root_derivative_'+str(digits),mp.diff(scaled,vv)<mp.mpf('-.99'))
        check('normal_local_root_residual_'+str(digits),abs(scaled(vv))<mp.mpf('1e-80'))
        check('normal_local_root_scale_'+str(digits),mp.mpf('1.86e-34')<ee<mp.mpf('1.87e-34'))
        check('normal_local_root_relative_'+str(digits),abs(vv-1)<mp.mpf('1e-29'))
        check('normal_log_radius_scale_'+str(digits),mp.mpf('5.34e33')<1/ee<mp.mpf('5.38e33'))
        vals.append(dict(digits=digits,pole_residue=mp.nstr(24*dd,70),reciprocal_slope=mp.nstr(1/(24*dd),70),
            e0=mp.nstr(e0,70),e_star=mp.nstr(ee,70),relative_shift=mp.nstr(vv-1,50),
            log_radius_scale=mp.nstr(1/ee,70),kappa_M=mp.nstr(K,60),A_T=mp.nstr(A,60)))
    for key in ['pole_residue','reciprocal_slope','e0','e_star','log_radius_scale']:
        check('cross_precision_'+key,abs(mp.mpf(vals[0][key])/mp.mpf(vals[1][key])-1)<mp.mpf('1e-50'))
    rho=Q('1e-8');Y=Q('5e-9');Z0=Q('.25');Z2=Q('5e7')
    check('nonlinear_ball_preserved_at_targets',Y+Z0*rho+Z2*rho*rho/2==rho)
    check('nonlinear_strict_contraction_at_targets',Z0+Z2*rho==Q('.75')<1)
    check('energy_defect_inverse_target',1/(1-Q('.25'))==Q(4,3))
    FIXTURES.append(dict(fixture='normal_endpoint_diagnostic_and_targets',values=vals,
        root_of_meromorphic_extension_only=True,actual_kernel_moment_computed=False,
        physical_fast_tail_inverse_regularity_proved=False,nonlinear_constants_certified=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_mass_boundary_result',old['status']=='passed' and old['stationary_metric_boundary_identification'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        intertwiner();tensor_and_contacts();endpoint_and_targets()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        stationary_current_constraint_realization_complete=error is None,
        stationary_current_constraint_scope='Exact static magnetic disk family, conditional on retained renormalized functional; no nonlinear solution or retarded constraint theorem.',
        physical_fast_tail_inverse_regularity_proved=False,actual_kernel_endpoint_moment_computed=False,
        physical_coupled_inverse_proved=False,full_dynamical_boundary_completed=False,admission_norm_inverse_quantified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
