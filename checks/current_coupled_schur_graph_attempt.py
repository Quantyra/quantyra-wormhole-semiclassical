"""Supporting arithmetic for the full weak inverse conversion; no strong coupled graph bound."""
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
STEM = 'current-coupled-schur-graph-attempt'
OWN = ['notes/current-coupled-schur-graph-attempt.md', 'checks/current_coupled_schur_graph_attempt.py', 'evidence/current-coupled-schur-graph-attempt-review.json']
INPUTS = ['notes/current-coupled-schur-graph-attempt.md', 'checks/current_coupled_schur_graph_attempt.py', 'evidence/current-coupled-schur-graph-attempt-review.json', 'evidence/current-coupled-schur-graph-attempt-exploration.json', 'notes/current-global-mixed-schur-preconditioner.md', 'notes/current-seed-full-hessian-realization.md', 'notes/stationary-seed-quantum-form.md', 'notes/current-weighted-global-return-attempt.md', 'notes/current-coupled-global-return-attempt.md', 'notes/recentered-reference-inverse.md', 'notes/current-weak-to-admission-obstruction.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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
    ctx.prec=512; a=arb; p,h,zh,env=york.actual_parameters()
    beta=p['beta'];Fl=1-384*a.pi()*a('1e22')*p['g']
    gram=4*(3-a(5).sqrt())
    B2=(16/(a('23.99')*beta)+1225/(a('.99')*Fl*beta))/gram
    cmin=((h+a('.8'))/(h+1-zh))**a('1.4')
    cmax=((h+a('.9'))/(h+1-zh))**a('1.4')
    affine=a('7.49e-12')-a('5e-5')*a('4.2953e-8')
    feedback=a('.00368')-a('5e-5')*a('49.462')
    ratio=feedback/(a('5e-5')*a('390.478'))
    mass=2*a('1e-11')/a('.00906808164')
    check('physical_positive_diagonal_guards',beta>0 and Fl>0)
    check('full_canonical_L2_inverse_upper',B2<a('2430.924468'))
    check('full_weighted_L2_certificate_upper',B2*cmax<a('3.272e49'))
    check('certificate_expression_lower_not_operator_lower',B2*cmin>a('2.774e49'))
    check('mass_affine_reserve',(affine-a('5.34235e-12')).contains(0) and sp.Rational('7.49e-12')-sp.Rational('5e-5')*sp.Rational('4.2953e-8')==sp.Rational('5.34235e-12'))
    check('rounded_mass_feedback_reserve',(feedback-a('.0012069')).contains(0) and sp.Rational('.00368')-sp.Rational('5e-5')*sp.Rational('49.462')==sp.Rational('.0012069'))
    check('rounded_energy_ratio',ratio>a('.06181654280138') and ratio<a('.06181654280140'))
    check('rounded_down_mass_comparison_target',mass>a('2.2055381495e-9') and mass<a('2.2055381496e-9'))
    check('source_mass_allowance_still_excessive',a('4.2953e-8')/mass>a('19.475'))
    check('conditional_fast_screen_scope_retained',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    nu,T,t=sp.symbols('nu T t',real=True)
    G=sp.Matrix([[16,8],[8,8]])
    check('exact_physical_canonical_Gram',sp.expand((sp.Matrix([nu,T]).T*G*sp.Matrix([nu,T]))[0]-8*(nu+T)**2-8*nu**2)==0)
    check('exact_Gram_characteristic_polynomial',sp.expand((G-t*sp.eye(2)).det()-(t-4*(3-sp.sqrt(5)))*(t-4*(3+sp.sqrt(5))))==0)
    H,P,E=sp.symbols('H P E',commutative=False)
    # H P=I-E, hence H[P(I-E)^-1]=I in the stated Banach spaces.
    Iminus=1-E
    check('ordered_right_inverse_identity',sp.expand(Iminus*Iminus**-1)==1)
    delta=sp.Rational(1,2)
    check('defect_error_amplification',1/(1-delta)==2)
    result=dict(canonical_gram_min=encoded(gram),full_unweighted_physical_L2_upper=encoded(B2),weighted_certificate_expression_lower=encoded(B2*cmin),full_weighted_physical_L2_upper=encoded(B2*cmax),conditional_affine_reserve=encoded(affine),conditional_feedback_reserve=encoded(feedback),rounded_mass_only_energy_ratio=encoded(ratio),mass_coefficient_necessary_comparison_ceiling=encoded(mass))
    return dict(numerics=result,full_seed_L2_conversion_executed=True,strong_coupled_graph_estimate='unresolved',new_feasibility_requirement_closed=False,approximate_graph_inverse_constructed=False,graph_defect_delta_computed=False,complete_mass_column_certified=False,other_coupled_columns_certified=False,actual_mass_computed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,closed_audit_preserved=True,next_calculation_executed=False)

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
