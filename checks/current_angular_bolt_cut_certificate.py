"""Full-tower endpoint error, including finite-cut action-adjoint terms."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_finite_local_action_assembly as local
york=local.york;endpoint=local.endpoint;prior=local.prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-angular-bolt-cut-certificate';PREVIOUS='current-finite-local-action-assembly'
OWN=['notes/'+STEM+'.md','checks/current_angular_bolt_cut_certificate.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    z,b,u,e,t=sp.symbols('z beta u eps tau',real=True)
    Q=(1-z)**2*(b*(1+2*z+3*z*z)+u*z*z+e*z**4)-t*z**3
    factor=(1-z)*(12*b*(1+3*z)+12*u*z+e*(2+2*z-10*z*z+30*z**3))+6*t*z
    check('global_curvature_factorization',sp.expand((2-sp.diff(Q,z,2)-factor).subs(u,1-6*b-e))==0)
    qf=sp.Function('Q')(z)
    check('warp_energy_derivative',sp.simplify(sp.diff(sp.diff(qf,z)**2/4-qf,z)-sp.diff(qf,z)*(sp.diff(qf,z,2)/2-1))==0)
    v=sp.symbols('v',positive=True)
    for name,f,critical in [('first',v/(1+v*v)**sp.Rational(3,2),1/sp.sqrt(2)),
                           ('second',v*v/(1+v*v)**sp.Rational(3,2),sp.sqrt(2))]:
        check('c_over_K_extremum_'+name,sp.simplify(sp.diff(f,v).subs(v,critical))==0 and sp.simplify(f.subs(v,critical)-2/(3*sp.sqrt(3)))==0)
    gg=2*v/(1+v*v)**2
    check('exact_gprime_maximum',sp.simplify(gg.subs(v,1/sp.sqrt(3))-9/(8*sp.sqrt(3)))==0)
    n,q,m=sp.symbols('n q m',positive=True,integer=True)
    check('exact_angular_count',sp.summation(2*(q+2*n),(n,1,m))==2*m*q+2*m*m+2*m)
    check('angular_count_majorant',sp.expand(4*m*(m+q)-2*m*(q+m+1)-2*m*(q+m-1))==0)
    w=sp.symbols('w',positive=True)
    check('all_thermal_boundary_sum_telescopes',sp.cancel(1/(w*(1+w))-(1/w-1/(w+1)))==0)
    E,L0,m0=sp.symbols('E L0 m0',positive=True)
    f=E**-4*(L0+2*sp.log(E/m0))
    primitive=-(2*L0+4*sp.log(E/m0)+1)/E**2
    check('log_weighted_angular_primitive',sp.simplify(sp.diff(primitive,E)+E**2*sp.diff(f,E))==0)
    check('log_weighted_angular_endpoint',sp.simplify(-primitive.subs(E,m0)-(2*L0+1)/m0**2)==0)
    k=sp.Rational(99,100)
    check('boundary_angular_exponent',2*(1+k)==sp.Rational(199,50))
    check('boundary_summability_strict',2*(1+k)>2)
    check('all_thermal_endpoint_coefficient',sp.Rational(41,10)*14*2==sp.Rational(574,5))
    check('global_to_local_adjoint_coefficient',3+11*sp.Rational(1,1000)/(1-k)==sp.Rational(41,10))

def actual_geometry(p,h,zh):
    beta,u,e=p['beta'],p['u'],p['eps'];v=1+h
    check('global_polynomial_Qsecond_upper',12*beta-44*e.abs_upper()>arb('1.9') and u>0 and p['tau']>0 and zh<1)
    chi=prior.hull(arb(0),arb('.001'))
    Kv=(v+1)*(u+3*beta*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    yy=2*(chi/2).sinh()**2;yy=prior.hull(max(arb(0),yy.lower()),yy.upper());t=1+h*(1+yy)
    AA=3*beta+u+e
    H=beta*(t+2*v)+AA/(t*v*v)-2*e*(2*t+v)/(t*t*v**3)+e*(3*t*t+2*t*v+v*v)/(t**3*v**4)
    Ht=beta-AA/(t*t*v*v)+4*e/(t*t*v**3)+4*e/(t**3*v*v)-3*e/(t*t*v**4)-4*e/(t**3*v**3)-3*e/(t**4*v*v)
    B=2*Kv*v+yy*H;g=((yy+2)/(t*B)).sqrt()
    apchi=(chi/2).cosh()*B.sqrt()/(arb(2).sqrt()*Kv*t*t.sqrt())
    Bprime=chi.sinh()*(H+yy*h*Ht);tprime=h*chi.sinh()
    apchi+=arb(2).sqrt()*(chi/2).sinh()*(Bprime/(2*Kv*t*t.sqrt()*B.sqrt())-3*B.sqrt()*tprime/(2*Kv*t*t*t.sqrt()))
    aproper=apchi/g
    sinch=prior.hull(arb(1),(chi.abs_upper()/2).cosh().upper())
    achi=sinch*B.sqrt()/(arb(2).sqrt()*Kv*t*t.sqrt());arho=achi/g
    check('actual_sharp_cap_geometry',arho>arb('.999') and arho<arb('1.001') and aproper.abs_upper()<2)
    check('every_mass_inner_collar_covered',arb('.001')/18<arb('.001')*g.lower())
    c_over_K=(2+arb(1)/18)/(3*arb(3).sqrt())
    gpmax=9/(8*arb(3).sqrt())
    check('global_actual_adjoint_damping',c_over_K<arb('.4') and gpmax<arb('.65'))
    check('global_auxiliary_adjoint_bound',4*arb('.65')/arb('1.2')<3)
    delta=arb('.001');w0=arb('.5')
    check('sharp_inner_K_upper',(1/arb('.999')**2+(delta/w0)**2).sqrt()<arb('1.002'))
    check('sharp_inner_damping',2/arb('1.001')-4*delta**2/(w0**3*(1+w0))>arb('1.99'))
    check('sharp_inner_forcing',16*arb('1.002')/(1+w0)<11)
    NUMERICS.append(dict(kind='actual_geometry_and_auxiliary_adjoint',h=encoded(h),proper_over_chi=encoded(g),
        a_over_rho=encoded(arho),a_proper_derivative=encoded(aproper),c_over_K_upper=encoded(c_over_K),
        global_ell_upper='3',inner_ell_upper='4.1*(M*rho/.001)^.99',physical_paired_seam_retained=True,
        gravitational_adjoint_bounded=False))

def tower_budgets(C,eta):
    m0=arb(18);delta=arb('.001');k=arb(99)/100;r=2*(1+k)
    ss=6-arb.pi()**2/2;L0=(2/(eta*m0)).log()+arb(7)/2
    diff=12400*ss*C*eta**2*m0**4
    thermal=16*C*eta**2*m0**4*(2*L0+1)
    moment=4*C*r/(r-2)*m0**(2-r)
    boundary=arb('114.8')*delta**(-k)*(eta*m0**3)**(1+k)*moment
    return dict(difference=diff,thermal=thermal,boundary=boundary,total=diff+thermal+boundary,
        log_factor=L0,angular_boundary_moment_upper=moment)

def threshold(C,factor):
    lo=arb('1e-20');hi=arb('1e-15');target=arb('1e-10')
    for _ in range(100):
        mid=(lo+hi)/2;value=factor*tower_budgets(C,mid)['total']
        if value<target:lo=mid
        elif value>target:hi=mid
        else:break
    check('threshold_lower_sign_'+str(factor),factor*tower_budgets(C,lo)['total']<target)
    check('threshold_upper_sign_'+str(factor),factor*tower_budgets(C,hi)['total']>target)
    check('threshold_resolved_'+str(factor),hi-lo<arb('1e-32'))
    return prior.hull(lo,hi)

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    actual_geometry(p,h,zh)
    Q=local.Poly({(power,0):env[sym] for power,sym in zip(york.powers,york.qs)})
    qabs=sum((v.abs_upper() for v in Q.c.values()),arb(0))
    rho0=arb('.001')/18;zlo=zh-rho0*qabs.sqrt()
    check('global_trials_cap_coordinate_cover',zlo>arb('.999') and qabs<2)
    band=prior.hull(arb('.999'),arb(1))
    for trial in local.physical_trials(Q)[2:]:
        value=trial['u'].value(band).abs_upper()
        scaled=rho0*qabs.sqrt()*trial['u'].derivative().value(band).abs_upper()
        check('actual_TT_inner_norm_'+trial['name'],value<2 and scaled<2)
        NUMERICS.append(dict(kind='actual_physical_TT_inner_norm',name=trial['name'],
            value_upper=encoded(value),scaled_first_derivative_upper=encoded(scaled),inner_norm_upper='2',
            only_inner_collar=True,unknown_correction_or_adjoint_norm=False))
    eta=arb('2e-17');mass1=(arb(p['q']+1)/p['C']).sqrt();eps1=eta*(18/mass1)**3
    check('minimum_actual_mass',mass1>18)
    check('adaptive_cuts_meet_every_mass_cap',eta*18<arb('.001'))
    budgets=tower_budgets(p['C'],eta);target=arb('1e-10')
    final=8*budgets['total']
    NUMERICS.append(dict(kind='complete_angular_thermal_endpoint_budget',eta=encoded(eta),
        C=encoded(p['C']),
        first_mass=encoded(mass1),largest_actual_cut=encoded(eps1),
        per_half_unit_norm={key:encoded(val) for key,val in budgets.items()},
        maximum_exteriors=2,maximum_direction_norm=2,maximum_factor=8,
        full_budget=encoded(final),target='1e-10',units='unweighted paired-frequency action, complete massive angular tower',
        all_angular_levels=True,all_thermal_frequencies=True,finite_cut_adjoint_terms_included=True,
        zero_temperature_common_restoration_included=False,gravitational_residual_norm=False))
    check('complete_tower_endpoint_target',final<target,encoded(final))
    for factor in [1,2,8]:
        upper=threshold(p['C'],factor)
        check('selected_eta_meets_threshold_'+str(factor),eta<upper.lower())
        NUMERICS.append(dict(kind='permitted_adaptive_cut_threshold',factor=factor,eta_threshold=encoded(upper),
            maximum_first_cut_for_mass_lower_bound_18=True,target='1e-10',threshold_applies_only_to_endpoint_budget=True))
    FIXTURES.append(dict(kind='actual_full_tower_endpoint_certificate',fixed_cut_divergence_not_nonexistence=True,
        actual_parameter_box_retained=True,three_global_TT_direction_norms_checked=True,
        full_physical_feasibility_test_executed=False))

def main():
    ctx.prec=256;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_local_action_passed',old['status']=='passed' and old['actual_finite_local_action_assembled'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),complete_angular_thermal_endpoint_error_bounded=status=='passed',
        actual_auxiliary_adjoint_bound_proved=status=='passed',finite_cut_adjoint_terms_included=True,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
        complete_renormalized_operator=False,full_angular_curved_remainder_certified=False,
        metric_admission_proved=False,nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
