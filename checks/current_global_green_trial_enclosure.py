"""Global original-equation error certificate for an explicit Green INPUT.

The retained signed quantum inverses and Z_r are not evaluated by this file.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-global-green-trial-enclosure'
sys.path.insert(0, str(ROOT / 'checks'))
from current_quartic_ritz_rank_test import arb, ctx, dec, enc, receipt, pick
sys.path.insert(0, str(ROOT))
from operators.mixed_green_trial import MixedGreenTrial, legendre_bounds

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_global_green_trial_enclosure.py',
    'operators/mixed_green_trial.py', 'checks/current_quartic_ritz_rank_test.py',
    'notes/current-filtered-mixed-green-residual.md',
    'notes/current-global-mixed-schur-preconditioner.md',
    'notes/current-finite-local-action-assembly.md',
    'notes/current-signed-observation-regularization-test.md',
    'evidence/current-filtered-mixed-green-residual.json',
    'evidence/current-filtered-mixed-green-residual-verification.json',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-joint-radius-composition-test.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-signed-observation-regularization-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def last(stem):
    obj = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())
    return obj.get('runs', [obj])[-1]


def symbolic():
    z, delta, u, eps = sp.symbols('z delta u eps')
    beta = (1-u-eps)/6
    L = 1-delta
    P = beta*(1+2*z+3*z*z)+u*z*z+eps*z**4
    PL = P.subs(z, L)
    tau = delta**2*PL/L**3
    Q = (1-z)**2*P-tau*z**3
    p1, p2 = 8*beta+2*u+4*eps, 6*beta+2*u+12*eps
    ac = -3*p1*delta+3*p2*delta**2-40*eps*delta**3+15*eps*delta**4-3*delta**2*PL/L**2
    gamma = delta*PL-delta**2*sp.diff(P,z).subs(z,L)/2+3*tau*L**2/2
    assert sp.cancel(Q.subs(z,L)) == 0
    assert sp.cancel(sp.diff(Q,z).subs(z,L)+2*gamma) == 0
    assert sp.cancel(sp.diff(Q,z,2).subs(z,L)/2-1-ac) == 0
    cs, d, v, gam, ellL = sp.symbols('c d v gamma L')
    Ah = cs*ellL**2/(d*gam)
    At = cs*ellL**2/(gam*(1-d*v))
    assert sp.cancel(At*(1/d-v)-Ah) == 0
    assert sp.cancel(d*gam*Ah/ellL**2-cs) == 0


def execute():
    ctx.prec = 512
    runtime = last('current-coupled-mass-adjoint-runtime')
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected, name
    verification = last('current-filtered-mixed-green-residual-verification')
    assert verification['status'] == 'ordinary_filtered_map_argument_and_frozen_arithmetic_reviewed'
    for name, expected in verification['artifact_sha256'].items():
        assert digest(ROOT / name) == expected, name
    symbolic()

    endpoint = pick(receipt('current-finite-local-action-assembly'), 'actual_seed_correlated_endpoint')
    c = receipt('current-joint-radius-composition-test')['result']['constants']
    e = receipt('current-coupled-truncation-energy-error')['result']['parameters']
    previous = last('current-filtered-mixed-green-residual')
    assert previous['status'] == 'filtered_input_bound_passed_complete_signed_response_unresolved'
    f = previous['result']['constants']
    C, B2 = dec(e['C']), dec(e['B2'])
    u, eps = 1/C, dec(endpoint['Q_coefficients'][-1])
    h = dec(endpoint['h'])
    delta, L = h/(1+h), 1/(1+h)
    beta = (1-u-eps)/6
    assert beta.overlaps(dec(c['beta'])) and eps < 0
    p1, p2 = 8*beta+2*u+4*eps, 6*beta+2*u+12*eps
    PL = 1-p1*delta+p2*delta**2/2-4*eps*delta**3+eps*delta**4
    ac = -3*p1*delta+3*p2*delta**2-40*eps*delta**3+15*eps*delta**4-3*delta**2*PL/L**2
    a = 1+ac
    tau = delta**2*PL/L**3
    gamma = delta*PL-delta**2*(p1-p2*delta+12*eps*delta**2-4*eps*delta**3)/2+3*tau*L**2/2
    assert gamma.overlaps(dec(endpoint['gamma']))
    d = dec(c['d'])
    lam = C*beta/B2
    anom = 11*lam/(360*arb.pi())
    v0 = 2-2*eps-2*u*delta*(2-delta)+12*beta*anom
    vgap = 2*delta*(2-delta)*(1-u)-2*L**2*ac-2*eps+12*beta*anom
    xx = vgap/(L**2*a)
    ell_shift = 2*xx/(3+(9+4*xx).sqrt())
    ell = 1+ell_shift
    nu = (-1+(1+4/(d*L**2*a)).sqrt())/2
    assert ell_shift > 0 and ell > 1 and nu > 7000
    assert arb('.99') < a < arb('1.01') and 0 < d*v0 < arb('.001')
    Q = [beta, arb(0), u, -4*beta-2*u-tau, 3*beta+u+eps, -2*eps, eps]
    # The imported list is sparse. Verify every correspondence explicitly.
    for power, coefficient in zip([0,2,3,4,5,6], endpoint['Q_coefficients']):
        assert Q[power].overlaps(dec(coefficient)), power
    Mactual = sum((abs(Q[j])*j*(j-1)*(j-2) for j in range(3,7)), arb(0))
    M = arb(1000)
    assert Mactual < M
    p = dict(L=L, gamma=gamma, d=d, beta=beta, C=C, B2=B2, a=a, v0=v0,
             ell=ell, nu=nu, alpha=dec(f['alpha_original_radius_cost']),
             h=h, u=u, epsilon_seed=eps, lambda_phys=lam)
    trial = MixedGreenTrial(p)
    SD, SP, Ah, At = trial.SD, trial.SP, trial.Ah, trial.At
    assert 2*SD < SP and L-2*SP > arb('.99')
    assert 1+M*2*SP/(6*a) < 2
    near_moment = (arb(3).log()**2/3+2*arb(3).log()/9+arb(2)/27)/4+arb(1)/9
    assert near_moment < 1
    A, B, g, W = 2+M/(6*a), 5*M/(6*a), gamma/a, arb(2)
    Dgeom = W*Ah*(A/L**2*g**arb('1.5')+d*(B+4)*gamma*(2*SD).sqrt())
    Pgeom = W*At*((A/L**2*(v0+1/d)+8/C)*g**arb('1.5')+2*B*gamma*(2*SP).sqrt())
    def qp(k,S):
        return legendre_bounds(k, 1+a*S/gamma, gamma)
    def qmax(S):
        return 4*gamma*S+4*a*S*S+4*M*S**3/3, 2*gamma+4*a*S+2*M*S*S
    UnD, VnD = qp(nu, SD)
    qm, qs = qmax(SD)
    Dcut = W*SD.sqrt()*d*Ah*((15*qm/SD**2+2*(qs+2*qm)/SD)*UnD+8*VnD/SD)
    UlP, VlP = qp(ell, SP)
    UnP, VnP = qp(nu, SP)
    qm, qs = qmax(SP)
    Pcut = W*SP.sqrt()*At*((15*qm/SP**2+2*qs/SP)*(UlP+UnP)+8*(VlP+VnP)/SP)
    Ptail = W*(2*SP).sqrt()*Ah*UnD
    eD, eP = Dgeom+Dcut, Pgeom+Pcut+Ptail
    K = dec(f['preconditioned_residual_to_complete_physical_multiplier'])
    rm = dec(f['r_mix'])
    additional = K*(1+rm)*(eD+eP)
    full_input = dec(f['exact_two_resolvent_input_charge'])+additional
    assert eD < arb('1e27') and eP < arb('1e29') and eD+eP < arb('4e30')
    assert full_input < arb('8e-12')
    observation = last('current-signed-observation-regularization-test')
    assert observation['status'] == 'unresolved'  # The complete signed response is unresolved.
    assert observation['input_sha256_before'] == observation['input_sha256_after']
    assert observation['result']['scope']['observation_difference_bounded_for_complete_structured_map']
    obs = observation['result']['constants']
    known_subtotal = dec(obs['known_partial_total_with_observation_error'])+full_input
    print('Global original-operator input residual and certificate debit bounded', flush=True)

    values = []
    for label, s in [('heavy_core', d*gamma), ('heavy_outer', 100*d*gamma),
                     ('thermal_scale', gamma), ('slow_tail', arb('1e-20')),
                     ('outer_cutoff', 2*SP)]:
        point = trial.at(s)
        assert point['h'].is_finite() and point['t'].is_finite()
        values.append({'label': label, 's': enc(s), 'h_hat': enc(point['h']),
                       't_hat': enc(point['t']), 'evaluation': point['evaluation']})
        print('Validated explicit Green trial evaluation: '+label, flush=True)
    bolt = trial.at_bolt()
    assert all(x.is_finite() for x in bolt.values())
    baseline = json.loads(BASELINE.read_bytes())
    assert all(digest(Path(name)) == expected for name, expected in baseline.items())
    constants = dict(a=a, ell_minus_one=ell_shift, nu=nu, gamma=gamma, L=L,
                     A_h=Ah, A_t=At, original_Q_third_derivative_upper=Mactual,
                     D_geometry=Dgeom, D_cutoff=Dcut, P_geometry=Pgeom,
                     P_cutoff=Pcut, P_tail=Ptail, e_D=eD, e_P=eP,
                     sum_of_input_errors=eD+eP, additional_original_bolt_debit=additional,
                     complete_input_debit=full_input,
                     known_partial_subtotal_not_free_margin=known_subtotal)
    return {
        'precision_bits': ctx.prec, 'constants': {k:enc(v) for k,v in constants.items()},
        'trial_parameters': {k:enc(v) for k,v in p.items()},
        'parameter_use': 'retain the correlated exact definitions; these enclosures are not independent midpoint choices',
        'bolt_expansion': {k:enc(v) for k,v in bolt.items()}, 'trial_values': values,
        'sources': ['https://dlmf.nist.gov/14.12.E6', 'https://dlmf.nist.gov/14.3.E10',
                    'https://dlmf.nist.gov/14.8.E9'],
        'source_hypotheses': 'real x>1, order zero, degrees at least one; standard Q normalization retained',
        'unknowns': ['complete_T_t_and_F_t_quantum_applications', 'signed_Zr',
                     'subsequent_quantum_quadrature_and_application_errors', 'e_other',
                     'all_other_coupled_comparisons_and_margins'],
        'scope': {'explicit_global_computational_input':True,
                  'original_operator_global_input_errors_bounded':True,
                  'full_signed_quantum_response_evaluated':False,
                  'full_bolt_certificate_closed':False, 'original_source_defect_closed':False,
                  'graph_other_physical_and_mass_columns_closed':False,
                  'nonlinear_admission_contraction':False, 'exact_compatibility':False,
                  'physical_matching':False, 'semiclassical_validity_noise':False,
                  'all_sector_stability':False, 'construction_obstruction':False,
                  'full_objective_achieved':False},
        'protected_paths_unchanged':len(baseline),
    }


if __name__ == '__main__':
    hashes = {name:digest(ROOT/name) for name in INPUTS}
    result = {'started_utc':datetime.now(timezone.utc).isoformat(),
              'input_sha256_before':hashes,
              'frozen_sources_utf8':{name:(ROOT/name).read_text(encoding='utf-8') for name in INPUTS[:3]},
              'runtime_python_sha256':digest(Path(sys.executable))}
    try:
        result['result'] = execute()
        result['input_sha256_after'] = {name:digest(ROOT/name) for name in INPUTS}
        assert result['input_sha256_after'] == hashes
        result['status'] = 'explicit_global_input_enclosed_signed_quantum_response_unresolved'
    except Exception:
        result['status'] = 'failed'
        result['traceback'] = traceback.format_exc()
    result['finished_utc'] = datetime.now(timezone.utc).isoformat()
    out = ROOT/'evidence'/(STEM+'.json')
    data = json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    data['runs'].append(result)
    out.write_bytes((json.dumps(data,indent=2)+'\n').encode('utf-8'))
    constants = result.get('result',{}).get('constants',{})
    keys = ['e_D','e_P','sum_of_input_errors','complete_input_debit','known_partial_subtotal_not_free_margin']
    print(json.dumps({'status':result['status'],'traceback':result.get('traceback'),
                      'constants':{k:constants[k]['ball'] for k in keys if k in constants}},indent=2))
    if result['status']=='failed':raise SystemExit(1)
