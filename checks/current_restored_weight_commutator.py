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
STEM = 'current-restored-weight-commutator'
OWN = ['notes/current-restored-weight-commutator.md', 'checks/current_restored_weight_commutator.py', 'evidence/current-restored-weight-commutator-review.json']
INPUTS = ['notes/current-restored-weight-commutator.md', 'checks/current_restored_weight_commutator.py', 'evidence/current-restored-weight-commutator-review.json', 'evidence/current-restored-weight-commutator-exploration.json', 'notes/current-weighted-global-return-attempt.md', 'notes/current-coupled-global-return-attempt.md', 'notes/stationary-curved-quantum-tt.md', 'notes/stationary-global-tt-reduction.md', 'notes/stationary-seed-quantum-form.md', 'notes/current-seed-full-hessian-realization.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'evidence/current-weighted-global-return-attempt.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
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

def bell(xs,n,k):
    if n==0:return arb(int(k==0))
    if k==0:return arb(0)
    return sum((comb(n-1,j-1)*xs[j-1]*bell(xs,n-j,k-1) for j in range(1,n-k+2)),arb(0))

def execute():
    ctx.prec=512;a=arb;enc=encoded
    p,h,zh,env=york.actual_parameters()
    s=a('1.4');L=a('1.400000000001');dm=a('1.001')
    check('seed_Q_and_temperature_guards',h>0 and h<a('1e-30') and p['tau']/h<a('1e-30') and p['eps']<0)
    be=p['beta'];u=p['u'];ee=abs(p['eps'])
    check('P_first_second_jet_guards',8*be+2*u+4*ee<a('1.334') and 6*be+2*u+12*ee<a('1.001'))
    check('Q_high_derivative_guards',720*ee<1 and 1440*ee<1)
    theta=[a(1)]+[a(factorial(21)//factorial(21-k)*2**(k-1)*10**k) for k in range(1,6)]
    oz=[s*sum((comb(j-1,k)*theta[k]*dm**k*factorial(j-1-k) for k in range(j)),a(0)) for j in range(1,7)]
    zj=list(map(a,['1.001','2','4.51','16.6','72','400']))
    oj=[sum((oz[k-1]*bell(zj,j,k) for k in range(1,j+1)),a(0)) for j in range(1,7)];oj[0]=L
    ej=[sum((bell(oj,j,k) for k in range(1,j+1)),a(0)) for j in range(1,7)]
    bounds=list(map(a,['1.401','302','121824','4.774e7','1.803e10','6.623e12']))
    check('all_relative_weight_jets',all(x<y for x,y in zip(ej,bounds)))
    R=a(1)/400;coef=2*a.pi()*a('.002671854')/720
    p0=coef*ej[5]+a('.5')*ej[3]
    q6=ej[5]*R*R/4+6*ej[4]*R+30*ej[3]
    p2=coef*q6+2*a.pi()*(1152*L*L+576*L*L*((2/(9*R)).log()+1))+6*L*L
    alpha=1/a(265)-1/a(512)
    angular=4*(-18*alpha).exp()+8*((1/(18*alpha)).log()+1)
    far=2*a('9.8304e8')**2*(a(8)/7)**2*16*L*L*angular
    local=a('1e24')*(L/128)**2*(L/128).exp()/2
    ep=a('1e-4');K=a('1e6');b=(ep/(4*K)).sqrt();gradient=K*(a('.25')+2*(R/b).log())
    check('positive_far_exponent_and_amplified_weight',alpha>0 and 18*alpha<1 and 18/(4*L)>1)
    check('infinite_angular_weighted_moment',angular<a('39.218'))
    check('actual_local_gluing_commutator',local<a('6.048e19'))
    check('actual_far_gluing_commutator',far<a('3.105e21'))
    check('entire_gluing_commutator',far+local<a('3.166e21'))
    check('old_gluing_allowance_closed',(far+local)*a('5.2e-29')<a('.001'))
    check('geodesic_constant_contact',p0<a('1.783e8'))
    check('geodesic_second_probe',p2<a(93062))
    check('actual_local_error_constant_contact',a('1e15')*L*L+a('1e6')*ej[2]<a('3e15'))
    check('gain_two_split_and_gradient',b<R and gradient<a('1.268e7'))
    check('fourth_reserve_even_after_localization',a.pi()*ep<1/(200*a.pi()))
    check('LLL_commutator',400*L*L/(12*a.pi())<21)
    check('Weyl_weight_second_jet',2*oj[1]<600 and a(2).sqrt()*oj[1]<424)
    check('aggregate_zeroth_reserve',far+local+a('1e17')+a('1e10')+a('1e8')+21<a('4e21'))
    check('aggregate_gradient_reserve',8*a('1e17')+gradient+a('1e8')<a('1e19'))
    physical=a('1e19')*a('2.55e-27')+a('4e21')*a('6.24e-29')
    check('canonical_physical_subbudget',physical<a('2.752e-7') and physical<a('.006398999997'))
    eta=3*a(2).sqrt()/(a('1.5')-s)
    check('weak_weighted_gauge_return',eta<a('42.427'))
    z,hh,beta,uu,eps,tau=sp.symbols('z h beta u eps tau')
    PP=beta*(1+2*z+3*z*z)+uu*z*z+eps*z**4;Q=(1-z)**2*PP-tau*z**3
    lhs=(hh+1-z)*(2*Q-z*sp.diff(Q,z))-2*z*Q
    rhs=(1-z)**3*(2*beta*(1+z)-2*eps*z**4)+tau*z**3*(1+z)+hh*(2*Q-z*sp.diff(Q,z))
    check('exact_weight_to_lapse_inequality_identity',sp.expand(lhs-rhs)==0)
    # Independent pointwise Weyl product identities, before the integration by parts.
    D,c,bb,g,u0,v0,G=sp.symbols('D c b grad u v G')
    check('fourth_form_weight_identity',sp.expand((D+c*u0)**2-(2*g+bb*u0)**2-D**2-(2*c*u0*D+c*c*u0*u0-4*g*g-4*bb*u0*g-bb*bb*u0*u0))==0)
    check('mixed_third_form_weight_identity',sp.expand(((D+2*g+(bb+c)*u0)*(G-v0)+(D-2*g+(-bb+c)*u0)*(G+v0))/2-D*G-(c*u0*G-2*g*v0-bb*u0*v0))==0)
    check('conditional_fast_screen_preserved',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    return dict(relative_weight_jets=[enc(x) for x in ej],flow_log_weight_jets=[enc(x) for x in oj],constant_geodesic_contact_over_C=enc(p0),second_geodesic_probe_over_C=enc(p2),far_angular_sum_over_C=enc(angular),local_gluing_commutator_over_C=enc(local),far_gluing_commutator_over_C=enc(far),full_gluing_commutator_over_C=enc(far+local),old_gluing_E0_cost=enc((far+local)*a('5.2e-29')),gain_two_gradient=enc(gradient),quantum_commutator_fourth_coefficient=enc(a.pi()/20000),canonical_TT_input_lower_order_E0_cost=enc(physical),weighted_eta_constant=enc(eta),seed_canonical_quantum_form_subbudget_closed=True,weak_weighted_constraint_and_gauge_bounds_closed=True,strong_constrained_return_certified=False,post_elimination_mass_column_certified=False,actual_mass_computed=False,original_physical_strong_comparisons_closed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,full_objective_complete=False,closed_audit_preserved=True,next_calculation_executed=False)

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
