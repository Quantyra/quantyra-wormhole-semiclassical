"""Arithmetic supporting the ordinary global-return attempt, not a PDE solver."""
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
STEM = 'current-global-return-energy-attempt'
OWN = ['notes/' + STEM + '.md', 'checks/current_global_return_energy_attempt.py',
       'evidence/' + STEM + '-review.json']
INPUTS = OWN + [
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-physical-york-dual-assembly.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/current-seed-boundary-phase-and-trace-defect.md',
    'notes/current-global-boundary-moment-return.md',
    'notes/current-coupled-local-graph-bound.md',
    'notes/current-coupled-outer-return-test.md',
    'notes/current-complete-strong-response-attempt.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-global-boundary-moment-return.json',
    'evidence/current-complete-strong-response-attempt.json',
]
CHECKS = []


def check(name, condition):
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


def execute():
    ctx.prec = 512
    a = arb
    gl, rq = a(2000)**-10, a('1.001') * a(2000)**-8
    fl = 1-a('1.179e-8')
    ct = a('5e15')
    trace_coeff = a(7)/12*3072+a(11)/6*484+10+a(11)/72*340+10*3148
    check('current_anomaly_polynomial_coefficient', trace_coeff < 50000)
    check('current_anomaly_trace_relative_to_U', a(50000)/(120*a.pi()) < a('1e6'))
    check('complete_quantum_tensor_component_budget',
          2*a('2e15')**2 + 2*(a('2e15')+a('5e5'))**2 < ct**2)
    check('full_radial_and_time_source_budget',
          a('1e15')+a('1e11')/2+1+2+a('1e4')+1+10/a.pi() < a('2e15'))
    egamma = ct*rq*(22/(a('.99')*fl)+16/a('23.99')).sqrt()
    cdual = 8*a.pi()*(22/(a('.99')*fl)+16/a('23.99')).sqrt()/a(5).sqrt()
    csq = a('1e22')*gl*(200*a('1.1e-27')+a('1.1e-27')**2+10*a('4.3e-25'))
    elocal = 1000*(a('5.5e-23')+csq)
    check('local_gradient_conversion', cdual < 1000)
    check('current_quantum_energy', egamma < a('9.354e-11'))
    check('current_local_energy', elocal < a('5.501e-20'))
    check('actual_current_seed_energy', egamma+elocal < a('1.1e-10'))
    c = a('.1')
    cy = a('2.5')+3*c*c/(a('.99')*a(20).sqrt()*a('8e9'))+21*c**3/(a('.99')*a(30).sqrt())
    york = cy/(a('.99')*fl).sqrt()
    # F/(lambda r^2)=1+z^2/x-2(a_mass/x)z^3+(Qc/x)z^4-2(d/x)z^5+(d/x)z^6.
    check('actual_outer_lapse_lower', 1-2*a('2.002')*c**3-3*a('8.6e-19')/a('8e9') > a('.99'))
    check('outer_transport_seven_guard',
          3*a('2.002')+2*a('3.001')*c+8*a('8.6e-19')/a('8e9') < 7)
    check('physical_York_L2_conversion', york < a('2.518'))
    k = gl/(90*a.pi())
    amplitude_factor = 120*k/(1-c)**6
    iv = c**a('2.5')/a(5).sqrt()*a('2.518')*a('1.1e-10')
    amp = amplitude_factor*iv
    far0 = amp/3
    farm = amplitude_factor*c**4/12
    check('whole_unsubtracted_tail_amplitude', amp < a('3.06e-46'))
    check('whole_unsubtracted_auxiliary_far_norm', far0 < a('1.02e-46'))
    check('mass_column_far_coefficient', farm < a('6.50e-39'))
    check('original_far_target_given_weak_mass_ceiling',
          a('1.02e-46')+a('6.50e-39')*a('1.5e28') < a('1e-10'))
    check('far_target_on_previously_admitted_family',
          a('1.02e-46')+a('6.50e-39')*a('.001')/6 < a('1.084e-42'))
    t, s = sp.symbols('t s')
    # Five ordinary derivatives of t^5/(1-t) sum the entire rising-factorial series.
    all_terms = sp.expand_func(sp.diff(t**5/(1-t), t, 5))*t
    check('infinite_rising_factorial_generating_function',
          sp.cancel(all_terms-120*t/(1-t)**6) == 0)
    odd_terms = (all_terms-all_terms.subs(t, -t))/2
    check('entire_odd_residue_kernel',
          sp.cancel(odd_terms-60*t*((1-t)**-6+(1+t)**-6)) == 0)
    n = sp.symbols('n')
    residue = ((s-4)*s*(s-3)*(s*(s-3)+2)).subs(s, -n)
    check('complete_B_residue_polynomial',
          sp.expand(residue+n*(n+1)*(n+2)*(n+3)*(n+4)) == 0)
    prior = json.loads((ROOT/'evidence/current-global-boundary-moment-return.json').read_bytes())['runs'][-1]
    check('retained_global_moment_evidence_passed', prior['status'] == 'passed')
    g0 = decode(prior['result']['g0'])
    euler = a(1).exp()
    integrals = [euler-1]
    for j in range(1, 7):
        integrals.append(euler-j*integrals[-1])
    chi_integral = c/euler*140*(integrals[3]-3*integrals[4]+3*integrals[5]-integrals[6])
    check('original_chi0_integral', chi_integral > c/euler and chi_integral < c)
    moment_rhs = -g0*chi_integral
    check('nonzero_exact_mass_compatibility_rhs',
          moment_rhs > a('1e-30') and moment_rhs < a('4.4e-30'))
    # The power statement is algebraic support, not a numerical substitute for (16).
    check('energy_normalized_transition_exponent', 4+3-2 == 5)
    eta, source, known = a('.25'), a('7.49e-12'), a('1e-14')
    joint = (source+eta*known)/(1-eta)
    check('conditional_joint_global_target', joint < a('1e-11'))
    check('sufficient_global_source_ceiling',
          (1-eta)*a('1e-11')-eta*known > a('7.49e-12'))
    check('retained_quantum_screen', a('.028913')+a('1e-11')*a('1.0706e10') < a('.1599164'))
    previous = json.loads((ROOT/'evidence/current-complete-strong-response-attempt.json').read_bytes())['runs'][-1]
    check('preceding_conditional_screen_evidence_passed',
          previous['status'] == 'passed' and previous['result']['conditional_screen_certified'])
    return {
        'outcome': 'unresolved_actual_global_graph_return',
        'current_quantum_energy_upper': encoded(egamma),
        'current_local_energy_upper': encoded(elocal),
        'actual_current_seed_energy_upper': '1.1e-10',
        'canonical_York_conversion_upper': encoded(york),
        'K_upper': encoded(k), 'unsubtracted_far_amplitude_upper': encoded(amp),
        'unsubtracted_far_norm_upper': encoded(far0), 'mass_column_far_coefficient': encoded(farm),
        'original_Efar_bound': '1.02e-46+6.50e-39 abs(m)',
        'sufficient_mass_ceiling_for_Efar': '1.5e28', 'actual_mass_ceiling_proved': False,
        'admitted_family_Efar_upper': '1.084e-42',
        'chi0_integral': encoded(chi_integral), 'required_mass_moment_rhs': encoded(moment_rhs),
        'energy_only_transition_estimator_norm': 'infinite',
        'quantum_fast_screen_upper': '.135973', 'conditional_quantum_screen_retained': True,
        'joint_data_conditional_upper': encoded(joint), 'S_coupled_computed': False,
        'actual_global_graph_return_computed': False, 'actual_controlled_solution': False,
        'full_goal_complete': False, 'normalization_changed': False,
        'auxiliary_far_kernel_is_physical_operator_replacement': False,
        'complete_angular_thermal_sums_and_contacts_retained': True,
        'research_decision': 'retain construction; reject energy-only global return; use full-equation graph parametrix',
        'next_calculation_executed': False,
    }


if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(), 'mode': mode,
           'input_sha256_before': hashes,
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
