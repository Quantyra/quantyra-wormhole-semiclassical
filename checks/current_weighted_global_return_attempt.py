"""Supporting arithmetic for the constrained weighted global-return attempt; no complete Hessian bound."""
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
STEM = 'current-weighted-global-return-attempt'
OWN = ['notes/current-weighted-global-return-attempt.md', 'checks/current_weighted_global_return_attempt.py', 'evidence/current-weighted-global-return-attempt-review.json']
INPUTS = ['notes/current-weighted-global-return-attempt.md', 'checks/current_weighted_global_return_attempt.py', 'evidence/current-weighted-global-return-attempt-review.json', 'evidence/current-weighted-global-return-attempt-exploration.json', 'notes/current-coupled-global-return-attempt.md', 'notes/stationary-curved-quantum-tt.md', 'notes/strong-interior-quantum-remainder.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-controlled-solution-audit-closeout.md', 'notes/current-uniform-bolt-certificate-obstruction.md', 'evidence/current-bolt-adjoint-domain.json', 'evidence/current-coupled-global-return-attempt.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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
    a=arb; enc=encoded
    p,h,zh,env=york.actual_parameters()
    s=a('1.4');L=s*(1+2*abs(p['eps'])).sqrt();den=h+1-zh
    kwlo=((h+a('.8'))/den)**s;kwhi=((h+a('.9'))/den)**s
    infl=18/(18-L);tilt=(1-L/36)**(-8)
    moment=192/(35*a.pi()**2)*a('.0035')*tilt
    naive=6*a('1e24')*a('5.2e-29')*kwlo
    out=dict(h=enc(h),log_weight_gradient_upper=enc(L),weight_condition_lower=enc(kwlo),weight_condition_upper=enc(kwhi),actual_massive_resolvent_inflation_upper=enc(infl),four_resolvent_inflation_upper=enc(infl**4),complete_tilted_sixth_moment_inflation_upper=enc(tilt),complete_tilted_sixth_moment_over_C_upper=enc(moment),blind_gluing_upper_bound_expression_lower=enc(naive),gluing_row_ceiling_for_001_over_C=enc(a('.001')/a('5.2e-29')),gluing_inflation_ceiling_for_001=enc(a('.001')/(6*a('1e24')*a('5.2e-29'))),proposed_quantum_fourth_commutator_coefficient=enc(1/(200*a.pi())),old_quantum_budget_margin_below_01=enc(a('.01')-a('.003601000003')),conditional_mass_affine_reserve=enc(a('7.49e-12')-a('5e-5')*a('4.2953e-8')),conditional_mass_graph_reserve_rounded=enc(a('.00368')-a('5e-5')*a('49.462')))
    u0,e0=sp.symbols('u0 eps0')
    check('exact_Q_parameter_definition_identity',sp.expand(6*((1-u0-e0)/6)+u0+e0-1)==0)
    # The parameter balls share u,eps: independent interval arithmetic
    # deliberately loses that correlation. The equality is proved above.
    check('Q_parameter_interval_consistency',(6*p['beta']+p['u']+p['eps']-1).contains(0))
    check('actual_Q_sign_hypotheses',p['eps']<0 and p['tau']>0 and p['beta']>0 and p['u']>0 and zh>0 and zh<1)
    check('weight_gradient_gap',L<a('1.400000000001') and L<18)
    check('actual_massive_resolvent_inflation',infl<a('1.084338'))
    check('four_resolvent_inflation_only_not_Hessian_claim',infl**4<a('1.382477'))
    check('complete_tilted_moment',tilt<a('1.373445') and moment<a('.002671854'))
    check('weight_condition_bounds',kwlo>a('1.14113e46') and kwhi<a('1.34571e46'))
    check('blind_similarity_certificate_exceeds_budget',naive>a('3.56034e42'))
    check('rounded_weighted_gluing_sufficient_target',a('1.923076923e25')*a('5.2e-29')<a('.001'))
    check('provisional_lower_quantum_losses_fit_earlier_margin',a('1e23')*a('1e-26')+a('1e25')*a('2e-28')<a('.01')-a('.003601000003'))
    check('conditional_fast_screen_unchanged',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    r=sp.symbols('r',positive=True);F=sp.Function('F')(r);w=sp.Function('w')(r);b=sp.Function('b')(r);xi=sp.Function('xi')(r)
    constraint=lambda bb,xx:sp.diff(r*r*F*bb,r)+(r*r*sp.diff(F,r)-2*r*F)*xx
    check('exact_weighted_TT_constraint_defect',sp.simplify(constraint(w*b,w*xi)-w*constraint(b,xi)-sp.diff(w,r)*r*r*F*b)==0)
    J=sp.Function('J')(r);db=-J/(r*r*F)
    check('exact_TT_correction',sp.simplify((sp.diff(r*r*F*db,r)+sp.diff(w,r)*r*r*F*b).subs(sp.diff(J,r),sp.diff(w,r)*r*r*F*b))==0)
    n,N,q,C=sp.symbols('n N q C',positive=True,integer=True)
    check('exact_second_contact_partial_sum',sp.simplify(sp.summation(2*(q+2*n),(n,1,N))-(2*q*N+2*N*(N+1)))==0)
    # For n>=q, 2C(q+2n)/(n(n+q)) - 2C/n = 2C/(n+q)>0.
    check('fourth_contact_harmonic_lower',sp.simplify(2*C*(q+2*n)/(n*(n+q))-2*C/n-2*C/(n+q))==0)
    z=sp.symbols('z');be,uu,ee,tt=sp.symbols('beta u eps tau')
    Q=(1-z)**2*(be*(1+2*z+3*z*z)+uu*z*z+ee*z**4)-tt*z**3
    check('actual_thermal_polynomial_retained',sp.expand(Q-(be+uu*z*z-(4*be+2*uu+tt)*z**3+(3*be+uu+ee)*z**4-2*ee*z**5+ee*z**6))==0)
    out.update(outcome='unresolved_complete_coupled_return',actual_weight_domain_repair_executed=True,exact_TT_correction_executed=True,TT_correction_norm_certified=False,full_common_restored_commutator_certified=False,absolute_contact_estimator='divergent; not a divergent physical functional',blind_gluing_estimator='excessive known upper-bound expression; actual operator unknown',complete_post_elimination_columns_computed=False,actual_mass_coefficient_computed=False,original_norm_comparisons_closed=False,actual_controlled_solution=False,full_goal_complete=False,construction_obstruction=False,closed_audit_preserved=True,next_calculation_executed=False)
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
