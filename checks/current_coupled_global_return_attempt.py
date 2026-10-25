"""Arithmetic and exact identities supporting the coupled global-return attempt; no quantum global solve."""
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
STEM = 'current-coupled-global-return-attempt'
OWN = ['notes/current-coupled-global-return-attempt.md', 'checks/current_coupled_global_return_attempt.py', 'evidence/current-coupled-global-return-attempt-review.json']
INPUTS = ['notes/current-coupled-global-return-attempt.md', 'checks/current_coupled_global_return_attempt.py', 'evidence/current-coupled-global-return-attempt-review.json', 'evidence/current-coupled-global-return-attempt-exploration.json', 'notes/current-lifted-boundary-graph-attempt.md', 'notes/current-global-return-energy-attempt.md', 'notes/current-global-boundary-moment-return.md', 'notes/current-coupled-outer-return-test.md', 'notes/current-complete-strong-response-attempt.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-curved-paired-heat-bridge.md', 'notes/current-controlled-solution-audit-closeout.md', 'notes/fixed-clock-seed-bolt-feasibility.md', 'evidence/current-lifted-boundary-graph-attempt.json', 'evidence/current-global-boundary-moment-return.json', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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

def execute():
    ctx.prec=512
    a=arb
    enc=encoded; dec=decode
    p,h,z,env=york.actual_parameters();v=1+h;be=p['beta'];ga=p['gamma'];u=p['u'];eps=p['eps'];tau=p['tau']
    coeff=[be,a(0),u,-4*be-2*u-tau,3*be+u+eps,-2*eps,eps]
    qpp=sum((k*(k-1)*coeff[k]*z**(k-2) for k in range(2,7)),a(0))
    DA=z*z*qpp/2+3*ga*z
    thermal=-DA*be/(ga*ga*v)
    cold=-be/(ga*v*h)
    Fs=p['Fl']+23*p['g']/(180*a.pi())
    scalar=20*be*Fs
    prior=json.loads((ROOT/'evidence/current-global-boundary-moment-return.json').read_bytes())['runs'][-1]['result']
    e=dec(prior['exact_slow_exponent_offset']);bp=dec(prior['Bprime_at_slow_root'])
    md=(a(3)/5).sqrt()*a('.05');cf=1/(2*(a('.5')+e)).sqrt();dm=cf*a('6.50e-39')/bp*a('1.02')
    den=a('.999')-md*dm
    far_energy=a('1.02e-46')/a('1.1e-10')
    ce=(a(3).sqrt()*a('.05')**a('-1.5')*a('2.518')+md*cf*far_energy/bp)/den
    ma=json.loads((ROOT/'evidence/current-lifted-boundary-graph-attempt.json').read_bytes())['runs'][-1]['result']['lifted_boundary_block']
    m1=dec(ma['mass_slope']);k=a('5e-5');reserve=a('.00368')-k*m1
    out={'h':enc(h),'gamma':enc(ga),'D_A_squared':enc(DA),'reference_mass_relative_temperature_derivative':enc(thermal),'reference_mass_cold_displacement_derivative':enc(cold),'physical_bolt_value_per_canonical_conformal_amplitude':enc(v/h),'canonical_conformal_amplitude_for_cold_2e_minus9':enc(a('2e-9')*h/v),'scalar_quartic_original_row_coefficient':enc(scalar),'mass_energy_coefficient':enc(ce),'mass_graph_coefficient':enc(m1),'mass_target_energy_gain':enc(k*ce),'mass_target_graph_gain':enc(k*m1),'separate_global_feedback_remaining':enc(reserve),'sufficient_energy_to_graph_difference_ratio_if_no_other_loss':enc(reserve/(k*ce))}


    # Instantiate the previously unspecified localization at its existing .05/.1 radii.
    # The beta(11,11) log cutoff has the same C10/H9 bounds as the prior .02/.04 test.
    old=json.loads((ROOT/'evidence/current-lifted-boundary-graph-attempt.json').read_bytes())['runs'][-1]['result']
    cost=old['lift_original_graph_cost'];wexp=a(5)-a('4.005')
    H9new=2*(a('.1')/a('.04'))**wexp*dec(cost['profile_H9_upper'])
    K=p['g']*be/(15*a.pi());quantum_new=6000*K*H9new
    Lnew=dec(cost['classical_leading_raw_jets'])-dec(cost['other_series_raw_jets_upper'])-dec(cost['analytic_tail_upper'])-quantum_new
    mm=dec(old['lifted_boundary_block']['mass_at_target'])
    out.update(global_cutoff_profile_H9_upper=enc(H9new),global_cutoff_complete_normal_quantum_upper=enc(quantum_new),global_cutoff_conversion_lower_per_mass=enc(Lnew),
     necessary_mass_ceiling_if_original_and_lifted_graph_radii_retained=enc(2*a('1e-11')/Lnew),conditional_mass_allowance_to_original_radius_ceiling=enc(mm*Lnew/(2*a('1e-11'))),
     unit_energy_graph_joint_mass_gain_at_proposed_target=enc(k*(ce+m1)),unit_energy_graph_mass_only_target_for_global_allowance=enc(a('.00368')/(ce+m1)),
     refined_energy_far_return_coefficient=enc(far_energy))

    # Independent symbolic checks of the ordinary reference/matching identities.
    zz=sp.symbols('z',positive=True);beta,C=sp.symbols('beta C',positive=True)
    Q=sp.Function('Q')(zz);J=sp.Function('J')(zz)
    nm=-2*beta*Q*J/zz;bm=beta*zz**3/Q+beta*sp.diff(Q,zz)*J;FF=C*Q/zz**2
    def reduceJ(ex):
        for j in [3,2,1]:ex=ex.subs(sp.diff(J,zz,j),sp.diff(zz**3/Q**2,zz,j-1))
        return sp.cancel(ex)
    check('exact_mass_lift_zero_reference_lapse',reduceJ(nm+2*bm+zz*sp.diff(nm,zz))==0)
    check('exact_mass_lift_reference_mass_tangent',reduceJ(2*FF*(-zz*sp.diff(nm,zz)-bm)+zz*sp.diff(FF,zz)*nm-2*C*beta*zz)==0)
    nn=sp.Function('n')(zz);bb=sp.Function('b')(zz);cc=sp.Function('chi')(zz)
    f=lambda n,b:2*FF*(-zz*sp.diff(n,zz)-b)+zz*sp.diff(FF,zz)*n
    al=lambda n,b:n+2*b+zz*sp.diff(n,zz)
    check('full_cutoff_transition_F',sp.expand(f(cc*nn,cc*bb)-cc*f(nn,bb)+2*FF*nn*zz*sp.diff(cc,zz))==0)
    check('full_cutoff_transition_lapse',sp.expand(al(cc*nn,cc*bb)-cc*al(nn,bb)-nn*zz*sp.diff(cc,zz))==0)
    check('scalar_quartic_reference_lapse',al(zz**4,sp.Integer(0))==5*zz**4)
    ss=sp.symbols('s');aa=(ss+3)/(3*(ss+1))
    check('entire_slow_family_reference_lapse',sp.cancel(2-(1+ss)*aa-(3-ss)/3)==0)
    CV,rr=sp.symbols('Cv r');jc=-6*beta*(CV-rr*ss*(ss-3))*(ss-4)*(ss+1)
    mixed=sp.cancel((6-ss)*jc/(3*(4-ss)))
    check('original_scalar_quartic_row',sp.expand(mixed.subs(ss,4)-20*beta*(CV-4*rr))==0)
    check('differentiated_row_loses_quartic',sp.expand(zz*sp.diff(zz**4,zz)-4*zz**4)==0)
    BB,vv,gg,DD=sp.symbols('B v gamma DA',positive=True)
    kap=gg*C/BB;Dh=DD*C/BB**2;dR=-beta*BB/(gg*vv)
    dk=(-C*beta/(BB*vv**2)+(Dh-kap/(BB*vv))*dR)/kap
    check('uncut_reference_temperature_derivative',sp.cancel(dk+DD*beta/(gg**2*vv))==0)
    check('reference_matching_denominator',DA>a('.999999') and DA<a('1.000001'))
    check('uncut_reference_temperature_scale',abs(thermal)>a('6.515e65') and abs(thermal)<a('6.516e65'))
    check('actual_bolt_norm_conversion',v/h>a('1.977e33') and v/h<a('1.978e33'))
    check('scalar_original_row_size',scalar>a('3.33333329') and scalar<a('3.33333330'))
    check('actual_localization_H9',H9new<a('1.091e22'))
    check('whole_normal_quantum_localization_cost',quantum_new<a('2.260e-10'))
    check('actual_localization_conversion_lower',Lnew>a('.00906808164'))
    check('original_and_lifted_radius_mass_ceiling',2*a('1e-11')/Lnew<a('2.205539e-9'))
    check('conditional_allowance_not_sufficient_for_both_radii',mm*Lnew/(2*a('1e-11'))>a('19.69'))
    check('mass_difference_energy_coefficient',ce<a('390.478') and ce>a('390.477'))
    check('mass_difference_graph_coefficient',m1<a('49.462'))
    check('extra_energy_gain_exceeds_remaining_budget',k*ce>reserve)
    check('restricted_family_mass_only_safe_ratio',k*(m1+ce*a('.0618179'))<a('.00368'))
    check('augmented_equal_input_mass_target_fails',k*(ce+m1)>a('.0219969'))
    check('augmented_mass_only_safe_column_target',a('8.3647e-6')*(ce+m1)<a('.00368'))
    check('conditional_fast_screen_preserved',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    out.update(outcome='unresolved_complete_coupled_global_return',outer_data_only_physical_return_norm='infinite: actual-domain kernel, ordinary argument section 1',
        exact_reference_transition_and_matching_executed=True,original_scalar_phase_row_retained=True,
        global_quantum_phase_columns_computed=False,actual_interior_strong_inverse_computed=False,
        actual_mass_coefficient_computed=False,original_physical_norm_comparisons_closed=False,
        full_nonlinear_admission=False,actual_controlled_solution=False,construction_obstruction=False,
        full_objective_complete=False,closed_audit_preserved=True,next_calculation_executed=False)
    return out

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
