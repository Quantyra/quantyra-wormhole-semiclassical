"""Supporting arithmetic for exact boundary representation; complete coupled estimate unresolved."""
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
STEM = 'current-implicit-phase-graph-domain'
OWN = ['notes/current-implicit-phase-graph-domain.md', 'operators/exact_slow_phase.py', 'checks/current_implicit_phase_graph_domain.py', 'evidence/current-implicit-phase-graph-domain-review.json']
INPUTS = ['notes/current-implicit-phase-graph-domain.md', 'operators/exact_slow_phase.py', 'checks/current_implicit_phase_graph_domain.py', 'evidence/current-implicit-phase-graph-domain-review.json', 'evidence/current-implicit-phase-graph-domain-exploration.json', 'notes/current-coupled-outer-return-test.md', 'evidence/current-coupled-outer-return-test.json', 'notes/current-global-boundary-moment-return.md', 'notes/current-coupled-global-return-attempt.md', 'notes/current-global-mixed-schur-preconditioner.md', 'notes/current-coupled-schur-graph-attempt.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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
from operators.exact_slow_phase import ExactSlowPhase

def execute():
    ctx.prec=512; a=arb
    p,h,zh,env=york.actual_parameters()
    r=json.loads((ROOT/'evidence/current-coupled-outer-return-test.json').read_bytes())['runs'][-1]['result']
    K=p['C']*p['beta']**2/(15*a.pi()*p['B2'])
    phase=ExactSlowPhase(decode(r['L']),decode(r['kappa_M']),K,decode(r['full_normal_shape_particular']))
    e=phase.root_bracket
    check('unique_root_endpoint_and_derivative_guards',phase.endpoint_signs[0]>0 and phase.endpoint_signs[1]<0 and phase.derivative_interval<0 and e>0 and e<a('1e-30'))
    emax=a('1e-30');den=1-emax**2/4
    cot_rem=a.pi()**4/(720*den)
    cot_drem=a.pi()**4/720*(3/den+emax**2/(2*den**2))
    check('whole_cotangent_series_remainder_bounds',cot_rem<1 and cot_drem<3)
    hp=a('.5')+a.pi()**2/12+3*emax**2
    hu=a(3).digamma()+a('.5')*emax+a('10.8')+emax**3
    check('uniform_h_and_R_derivative_bounds',hu<12 and hp<2 and 5+2*emax+(6+10*emax+3*emax**2)*12+(6*emax+5*emax**2+emax**3)*2<120)
    B4=-24*K;Beta=phase.leading_approximation_residual();diff=phase.root_minus_eta_enclosure()
    check('binary64_loses_nonzero_slow_offset',4+float(e.mid())==4.0 and B4<0)
    check('analytic_eta_is_not_exact_root',Beta<a('-1.85045370e-67') and Beta>a('-1.85045371e-67') and diff>a('4.16352087e-67') and diff<a('4.16352090e-67'))
    z,q,L,kh,kk=sp.symbols('z q L kh kk')
    s=4+z;q=s*(s-3)
    R=6-L*z+5*z+z*z+(6+5*z+z*z)*z*kh
    check('full_normal_pole_clearing_identity',sp.expand(-kk*z*q*((q+2)*(1/z+kh)-L)+kk*(4+z)*(1+z)*R)==0)
    aa=lambda t:(t+3)/(3*(t+1))
    check('canonical_trace_pair_identity',sp.factor(aa(4+z)-sp.Rational(7,15)+2*z/(15*(5+z)))==0)
    check('all_implemented_jet_polynomials_exact',all(sp.expand((4+z)**j-4**j-z*sum(sp.binomial(j,k)*4**(j-k)*z**(k-1) for k in range(1,j+1)))==0 for j in range(9)))
    jet_results={}
    for x in [a(-5),a(-100),-a('1e35')]:
        stable=phase.scaled_pair_jets(x,8)
        independent=phase.scaled_independent_slow_jets(x,8)
        ee=(e*x).exp();aslow=(7+e)/(3*(5+e))
        direct=[(phase.B0*(aslow*(4+e)**j*ee-a(7)/15*a(4)**j),phase.B0*(a(4)**j-(4+e)**j*ee)) for j in range(9)]
        check('stable_jets_match_direct_'+x.str(10),all((u-v).contains(0) for st,dr in zip(stable,direct) for u,v in zip(st,dr)))
        check('independent_slow_trace_'+x.str(10),all((nu+aslow*b).contains(0) for nu,b in independent))
        jet_results[x.str(10)]=[encoded(v) for v in stable[0]]
    rows=phase.leading_normal_rows();g0=rows['differentiated_shape_pair']
    check('original_uncomputed_rows_not_zeroed',rows['original_shape_finite_part'] is None and rows['cutoff_and_global_return'] is None)
    delta=a('.005');xa=a('.05').log()
    def tail(n):
        xmin=-a(n)*a(10).log()
        return abs(g0)*(1+delta**2)**a('1.5')*(((-2*delta*xmin).exp()-(-2*delta*xa).exp())/(2*delta)).sqrt()
    threshold=(((a('1e-11')/abs(g0))**2*2*delta/(1+delta**2)**3+(-2*delta*xa).exp()).log())/(2*delta*a(10).log())
    check('leading_term_only_finite_cut_diagnostic',tail(300)<a('1.3900e-26') and tail(3300)>a('1.3906e-11') and threshold>a('3271.3534334') and threshold<a('3271.3534336'))
    check('slow_residual_exponents_below_original_fast_weight',4+e<a('4.005') and 4+phase.eta<a('4.005'))
    aff=sp.Rational('7.49e-12')-sp.Rational('5e-5')*sp.Rational('4.2953e-8')
    fb=sp.Rational('.00368')-sp.Rational('5e-5')*sp.Rational('49.462')
    check('unchanged_complete_budget_arithmetic',aff==sp.Rational('5.34235e-12') and fb==sp.Rational('.0012069') and sp.Rational('.028913')+sp.Rational('1e-11')*sp.Rational('1.0706e10')==sp.Rational('.135973'))
    numerics=dict(exact_root=encoded(e),B4=encoded(B4),B_eta=encoded(Beta),root_minus_eta=encoded(diff),g0=encoded(g0),leading_term_H3_at_1e_minus_300=encoded(tail(300)),leading_term_H3_at_1e_minus_3300=encoded(tail(3300)),leading_term_threshold_decimal_exponent=encoded(threshold),scaled_pair_values=jet_results)
    return dict(numerics=numerics,boundary_primitive_implemented=True,rounded_phase_estimator_failure_scope='Fixed-g0 original fast source graph; actual seed asymptotic argument in ordinary note.',complete_coupled_estimate='unresolved',new_feasibility_requirement_closed=False,approximate_graph_inverse_constructed=False,graph_defect_delta_computed=False,complete_mass_column_certified=False,other_coupled_columns_certified=False,actual_mass_computed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,closed_audit_preserved=True,next_calculation_executed=False)

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
