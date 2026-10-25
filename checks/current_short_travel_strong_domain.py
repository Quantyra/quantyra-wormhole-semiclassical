"""Complete long-travel transfer bound; not a strong coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-short-travel-strong-domain'
OWN=['notes/'+STEM+'.md','operators/short_travel_comparison.py',
     'checks/current_short_travel_strong_domain.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
    'notes/current-coupled-return-tail-attempt.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-geodesic-flow-multiplier.md',
    'notes/current-global-geodesic-return.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/strong-interior-quantum-remainder.md',
    'notes/outer-strong-curved-quantum-remainder.md',
    'notes/finite-graph-admissibility.md',
    'notes/finite-comparison-conormal-regularity.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-global-boundary-moment-return.md',
    'notes/current-implicit-phase-graph-domain.md',
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-geodesic-flow-multiplier.json']
CHECKS=[]


def check(name, condition, detail=''):
    CHECKS.append(dict(name=name,passed=bool(condition)))
    if not condition:
        raise AssertionError(name)


def encoded(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(42),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))


sys.path.insert(0,str(ROOT/'checks'))
argv=sys.argv[:]
sys.argv=[v for v in sys.argv if v!='--acceptance']
import current_physical_york_dual_assembly as york
sys.argv=argv
york.check=check;york.prior.check=check
sys.path.insert(0,str(ROOT))
from operators.short_travel_comparison import ShortTravelComparison

# Pin the actual imported parameter implementation, including transitive
# repository dependencies, rather than only its top-level helper.
for mod in list(sys.modules.values()):
    path=getattr(mod,'__file__',None)
    if path:
        path=Path(path).resolve()
        if path.is_relative_to(ROOT) and path.suffix=='.py':
            INPUTS.append(path.relative_to(ROOT).as_posix())
INPUTS=sorted(set(INPUTS))


def execute():
    ctx.prec=512
    a=arb
    p,h,zh,env=york.actual_parameters()
    comp=ShortTravelComparison()
    x,zeta=sp.symbols('x zeta',real=True)
    full=1-sp.cos(x)-x*x/2+x**4/24
    short=zeta*(1-sp.cos(x))-x*x/2+x**4/24
    transfer=(1-zeta)*(1-sp.cos(x))
    check('exact_full_short_transfer_identity',sp.expand(full-short-transfer)==0)
    check('constant_direction_retained',short.subs(x,0)==transfer.subs(x,0)==0)
    r,radius=sp.symbols('r a',positive=True)
    integral=sp.integrate(r**-5,(r,radius,sp.oo))
    check('entire_time_tail_integral',integral==1/(4*radius**4))
    check('full_angular_spectral_and_half_normalization',sp.simplify(72*2*2*sp.pi*integral-72*sp.pi/radius**4)==0)
    Y,Z,t=sp.symbols('Y Z t',positive=True)
    crossing=Y/(Y+Z)
    cubic=sp.integrate(((1-t)*Y-t*Z)**3,(t,0,crossing))+sp.integrate((t*Z-(1-t)*Y)**3,(t,crossing,1))
    expected=(Y**4+Z**4)/(4*(Y+Z))
    check('exact_crossing_cubic_kept',sp.cancel(cubic-expected)==0)
    check('resolved_cubic_Euler_degree',sp.cancel(Y*sp.diff(expected,Y)+Z*sp.diff(expected,Z)-3*expected)==0)
    j=sp.symbols('j')
    check('single_operator_resolvent_identity',sp.cancel(j/(1-j)-(1/(1-j)-1))==0)
    orders=[sp.Rational(2)+k*sp.Rational(7,4) for k in range(9)]
    check('four_and_eight_power_regularity_accounting',orders[4]==9 and orders[8]==16)
    check('physical_bolt_trace_not_inferred_from_H9',orders[4]<13 and orders[8]>13+1)
    # These finite point tests check code against its scalar definition only.
    for factor in [a('.5'),a(1),a('1.5'),a(2),a(3)]:
        rr=comp.radius*factor
        bb=comp.brackets(rr,a(7))
        check('bracket_implementation_'+str(factor), (bb['full']-bb['short']-bb['transferred']).contains(0))
    check('smooth_cutoff_endpoint_values',comp.cutoff(comp.radius)==1 and comp.cutoff(2*comp.radius)==0)
    check('smooth_cutoff_midpoint_value',comp.cutoff(a('1.5')*comp.radius).contains(a('.5')))
    j0=1000*p['C']**2/(p['g']*p['B2']**2)
    delta=comp.physical_form_transfer_upper(p['C'],p['B2'],p['g'],p['Fl'])
    check('original_physical_conversion_inputs',j0<a('5.2e-29') and p['Fl']>a('.99999998'))
    coarse=72*a.pi()*a(128)**4*a('5.2e-29')/(a('.99')*a('.99999998'))
    check('full_transfer_certified_with_original_units',delta<coarse and coarse<a('3.2e-18'))
    eta=a('.0022000000000000032');q=eta/(1-eta)
    check('entire_actual_remainder_after_transfer',a('.0022')+delta<eta)
    check('comparison_inverse_and_energy_convergence',eta<1 and q<a('.002205'))
    uaff=a('1.1e-10');daff=a('2.671175e-12');dfb=a('.00060345')
    check('retained_graph_affine_allowance_still_sufficient',a('6.0612e157')*q**60*uaff<daff)
    check('retained_graph_feedback_product_still_sufficient',a('1.5062e156')*q**60<dfb)
    check('extra_four_powers_cost_paid',q**56/q**60>a('4e10') and q**56/q**60<a('4.3e10'))
    return dict(outcome='unresolved_quantitative_coupled_return',
        comparison_change='H_s=H0-D_a; R_s=R+D_a; original physical H unchanged',
        travel_radius=encoded(comp.radius),local_travel_support_upper=encoded(2*comp.radius),
        j0_original_physical_conversion=encoded(j0),
        complete_transfer_base_L2_coefficient=encoded(comp.full_transfer_L2_upper(p['C'])),
        complete_transfer_physical_form_upper=encoded(delta),coarse_transfer_upper=encoded(coarse),
        certified_extra_form_allowance='3.2e-18',eta_s=encoded(eta),q_s=encoded(q),
        graph_64_iteration_tail_factor=encoded(q**60),
        physical_C5_64_iteration_tail_factor=encoded(q**56),
        physical_C5_strong_product_allowance_per_unit_original_error=encoded(q**-56),
        four_extra_power_tail_cost=encoded(q**-4),
        energy_resolvent_conversion_four=encoded((1-eta)**-4),
        energy_resolvent_conversion_eight=encoded((1-eta)**-8),
        fixed_seed_domain_argument='Ordinary argument in sections 3-4; supporting arithmetic does not prove that continuum argument.',
        graph_phase_bolt_value_power=4,physical_C5_power=8,
        high_order_kernel_constants=None,actual_strong_inverse_constants=None,M_graph_4=None,M_X_8=None,
        actual_physical_operator_applications=0,all_coupled_centers_computed=False,kappa_mass=None,actual_mass=None,
        nonlinear_neighborhood_proved=False,nonlinear_admission=False,exact_compatibility=False,
        physical_matching=False,validity_noise=False,all_sector_stability=False,
        actual_controlled_solution=False,construction_obstruction=False,
        numerical_strong_return_executed=False,closed_audit_preserved=True)


if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    hashes={p:digest(ROOT/p) for p in INPUTS}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),
                 flint_version=flint.__version__,sympy_version=sp.__version__,
                 flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*'))
                                  if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,
             input_sha256_before=hashes,precision_bits=512,runtime=runtime,
             frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback'))
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:
        f.write((json.dumps(summary,indent=2)+'\n').encode())
    print(json.dumps(summary,indent=2))
    if run['status']!='passed':raise SystemExit(1)
