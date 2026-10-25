"""Global orbit geometry and complete seam-tail verification; no coupled inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-global-geodesic-return'
OWN = ['notes/current-global-geodesic-return.md', 'operators/radial_geodesic_return.py', 'checks/current_global_geodesic_return.py', 'evidence/current-global-geodesic-return-review.json']
INPUTS = ['notes/current-global-geodesic-return.md', 'operators/radial_geodesic_return.py', 'checks/current_global_geodesic_return.py', 'evidence/current-global-geodesic-return-review.json', 'evidence/current-global-geodesic-return-exploration.json', 'notes/current-geodesic-two-sided-remainder.md', 'notes/current-geodesic-flow-multiplier.md', 'notes/stationary-curved-quantum-tt.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-global-mixed-schur-preconditioner.md', 'notes/current-coupled-global-return-attempt.md', 'notes/current-implicit-phase-graph-domain.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'evidence/current-geodesic-two-sided-remainder.json', 'evidence/current-geodesic-flow-multiplier.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
CHECKS = []


def check(name, condition, detail=""):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def encoded(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(45), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def decode(value):
    p = value['decimal_enclosure']
    return arb(int(p['mid']), int(p['rad'])) * arb(10) ** int(p['exponent'])


from math import comb, factorial
sys.path.insert(0, str(ROOT/'checks'))
# A legacy parameter helper parses argv on import. Keep this wrapper's
# acceptance mode local while importing it, then restore the exact argv.
original_argv = sys.argv[:]
sys.argv = [v for v in sys.argv if v != '--acceptance']
import current_physical_york_dual_assembly as york
sys.argv = original_argv
york.check = check
york.prior.check = check



sys.path.insert(0,str(ROOT))
from operators.radial_geodesic_return import RadialGeodesicReturn

def execute():
    ctx.prec=512;a=arb
    p,h,zh,env=york.actual_parameters()
    flow=RadialGeodesicReturn.from_correlated_seed(p['beta'],p['u'],p['eps'],h)
    check('original_temperature_preserved',flow.gamma.overlaps(p['gamma']) and flow.tau.overlaps(p['tau']))
    z,b,u,e,t,v=sp.symbols('z b u e t v',positive=True)
    P=b*(1+2*z+3*z*z)+u*z*z+e*z**4
    Q=(1-z)**2*P-t*z**3
    aa=4*b+2*u+t;bb=3*b+u+e
    Qexpanded=b+u*z*z-aa*z**3+bb*z**4-2*e*z**5+e*z**6
    check('complete_Q_polynomial_identity',sp.expand(Q-Qexpanded)==0)
    bracket=u-(6*b+2*u)*z+e*z*z*(2-3*z)-3*t*z/(2*(1-z))
    check('whole_interval_unique_peak_identity',sp.factor(sp.diff(Q,z)-2*z*(1-z)*bracket)==0)
    tau=v*(v-1)**2*P.subs(z,1/v)
    kval=(v+1)*(u+3*b*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    check('exact_correlated_bolt_identities',sp.factor(Q.subs({t:tau,z:1/v}))==0 and sp.factor(sp.diff(Q,z).subs({t:tau,z:1/v})+2*(v-1)*kval)==0)
    x=sp.symbols('x')
    scaled=x*x-aa*x**3+bb*u*x**4-2*e*u*u*x**5+e*u**3*x**6
    check('correlated_height_identity',sp.expand(Qexpanded.subs(z,u*x)-b-u**3*scaled)==0)
    check('whole_interval_monotonicity_guard',6*p['beta']+2*p['u']-5*abs(p['eps'])>a('.99'))
    check('implicit_peak_and_separatrix_guards',flow.peak_x>a('.99') and flow.peak_x<a('1.01') and flow.separatrix_outer_x>a('1.49') and flow.separatrix_outer_x<a('1.51') and flow.peak_height_scaled>a('.3333333333') and flow.peak_height_scaled<a('.3333333334'))
    qppguard=sum((abs(flow.coefficients[k])*k*(k-1) for k in range(2,7)),a(0))
    check('global_period_lower_guard',qppguard<12 and 10/(2*a(6).sqrt())>2)
    periods=dict(radial=flow.radial_period(),crossing_half_beta=flow.crossing_period(a(1)/2),trapped_one_sixth=flow.trapped_period(a(1)/6))
    check('regular_bolt_period_enclosure',periods['radial']>a('156.704452148') and periods['radial']<a('156.704452150'))
    check('crossing_period_enclosure',periods['crossing_half_beta']>a('6.43647293852') and periods['crossing_half_beta']<a('6.43647293856'))
    check('trapped_period_enclosure',periods['trapped_one_sixth']>a('1529608.55545') and periods['trapped_one_sixth']<a('1529608.55547'))
    measure=flow.trapped_measure_enclosure()
    check('entire_trapped_measure_enclosure',measure>a('1.0993804115e-26') and measure<a('1.0993804117e-26'))
    B4=x**4-2*x**3+x*x-sp.Rational(1,30)
    check('periodic_Bernoulli_jump_and_mean',sp.diff(B4,x,3).subs(x,0)-sp.diff(B4,x,3).subs(x,1)==-24 and sp.integrate(B4,(x,0,1))==0)
    check('full_crossing_momentum_integral',sp.integrate((1-x*x)**3,(x,0,1))==sp.Rational(16,35))
    n,q=sp.symbols('n q',positive=True)
    check('complete_angular_sixth_moment_telescope',sp.factor(2*(q+2*n)/(n*n*(n+q)**2)-2/q*(1/n**2-1/(n+q)**2))==0)
    zz,j=sp.symbols('zz j');qq=sp.Function('Q')(zz);f=sp.Function('f')(zz);A=qq-j*j
    L2=lambda g:A*sp.diff(g,zz,2)+sp.diff(qq,zz)*sp.diff(g,zz)/2
    expr=A*A*sp.diff(f,zz,4)+3*A*sp.diff(qq,zz)*sp.diff(f,zz,3)+(2*A*sp.diff(qq,zz,2)+3*sp.diff(qq,zz)**2/4)*sp.diff(f,zz,2)+(A*sp.diff(qq,zz,3)/2+sp.diff(qq,zz)*sp.diff(qq,zz,2)/4)*sp.diff(f,zz)
    check('global_regular_fourth_flow_derivative_identity',sp.expand(L2(L2(f))-expr)==0)
    check('entire_tail_prefactor',16*144*sp.Rational(16,35)==sp.Rational(36864,35))
    k=a(2000);C0=(1+12*a.pi()*k**6).sqrt();b02=k**10*(C0-1)/6;L=(p['B2']/b02).log()/2
    tails={}
    for om,guard in [(100,'790.478560'),(10000,'.001004145567'),(1000000,'1.219301126e-9')]:
        bound,chi=flow.complete_mass_seam_tail(om,p['C'],p['q'],L)
        check('complete_seam_tail_'+str(om),bound>0 and bound<a(guard) and chi>2)
        tails[str(om)]=dict(form_tail_upper_expression=encoded(bound),whole_angular_chi_upper=encoded(chi),units='quadratic comparison form divided by original time period, per unit mass squared; not kappa_m')
    low=flow.mass_seam_coefficient(0,periods['radial'],1)
    check('subtraction_low_frequency_cancellation_cost',low>2000 and low<2030)
    check('unchanged_coupled_budget_arithmetic',sp.Rational('7.49e-12')-sp.Rational('5e-5')*sp.Rational('4.2953e-8')==sp.Rational('5.34235e-12') and sp.Rational('.00368')-sp.Rational('5e-5')*sp.Rational('49.462')==sp.Rational('.0012069'))
    result=dict(u=encoded(p['u']),peak_z=encoded(flow.peak_z),peak_scaled=encoded(flow.peak_x),height_over_u_cubed=encoded(flow.peak_height_scaled),separatrix_outer_z=encoded(p['u']*flow.separatrix_outer_x),trapped_measure=encoded(measure),periods={k:encoded(v) for k,v in periods.items()},mass_seam_low_frequency_coefficient=encoded(low),seam_tails=tails)
    return dict(numerics=result,global_orbit_geometry_implemented=True,all_trajectory_families_retained=True,complete_seam_tail_certified=True,smooth_periodic_H4_estimator_failed=True,full_profile_J4_certified=False,complete_global_operator_evaluated=False,complete_coupled_estimate='unresolved',new_feasibility_requirement_closed=False,approximate_graph_inverse_constructed=False,graph_defect_delta_computed=False,complete_mass_column_certified=False,other_coupled_columns_certified=False,actual_mass_computed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,closed_audit_preserved=True,next_calculation_executed=False)

if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    import flint
    flint_root = Path(flint.__file__).resolve().parent
    runtime = {'python': sys.version, 'python_executable_sha256': sha256(Path(sys.executable).read_bytes()).hexdigest(),
               'flint_version': flint.__version__, 'sympy_version': sp.__version__,
               'flint_file_sha256': {str(p.relative_to(flint_root)): sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(flint_root.rglob('*'))
                                    if p.is_file() and p.suffix in ['.py', '.pyd', '.dll']}}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(), 'mode': mode,
           'input_sha256_before': hashes, 'precision_bits':512, 'runtime':runtime,
           'frozen_sources_utf8': {p: (ROOT/p).read_text(encoding='utf-8') for p in OWN}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('all_frozen_inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs': []}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger, indent=2)+'\n').encode())
    log = ROOT/'evidence'/(STEM+'.log')
    text = json.dumps({'mode': mode, 'status': run['status'], 'checks': len(CHECKS),
                       'result': run.get('result'), 'traceback': run.get('traceback')}, indent=2)+'\n'
    with log.open('ab') as f:
        f.write(text.encode())
    print(text)
    if run['status'] != 'passed':
        raise SystemExit(1)
