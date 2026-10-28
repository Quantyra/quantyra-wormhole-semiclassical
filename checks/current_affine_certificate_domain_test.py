"""Outward arithmetic for an actual-source phase/estimator test, not a solver."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-affine-certificate-domain-test'
RUNTIME = ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_affine_certificate_domain_test.py',
    'notes/current-seed-boundary-phase-and-trace-defect.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-finite-local-action-assembly.md',
    'notes/current-coupled-application-residual-bound.md',
    'notes/current-admission-bottleneck-estimator-test.md',
    'notes/current-original-admission-map-budget-test.md',
    'notes/current-native-output-structured-source-test.md',
    'notes/current-schur-polynomial-return-test.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(v):
    d = v['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10)**int(d['exponent'])


def enc(x):
    m, r, e = x.mid_rad_10exp(65)
    return {'ball': x.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def read_run(stem):
    obj = json.loads((ROOT / 'evidence' / (stem+'.json')).read_bytes())
    run = obj.get('runs', [obj])[-1]
    assert run['status'] == 'passed'
    return run


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root']) / name) == expected

    local = read_run('current-finite-local-action-assembly')
    restored = read_run('current-canonical-heat-restoration-assembly')
    endpoint = next(x for x in local['numerics'] if x['kind'] == 'actual_seed_correlated_endpoint')
    assert endpoint['exact_Q_bolt'] == '0'
    L = dec(endpoint['z_h'])
    beta = dec(endpoint['Q_coefficients'][0])
    co = next(x for x in restored['numerics'] if x['kind'] == 'common_restoration_coefficients')
    C, B2 = dec(co['C']), dec(co['B2'])
    q = arb(2000)**4
    lam = arb(2000)**-10
    eps = 1/(60*arb.pi())
    cr = arb('1e22')
    Fl = 1-384*arb.pi()*cr*lam
    chi = Fl+arb(23)/3*eps*lam
    aref = 10*eps*q**2
    cnu = (chi-1)*aref/(C*B2)
    trace_coefficient = (chi-1)*aref/(lam*B2**2)
    MR = ((arb(16)/arb('23.99')+1225/(arb('.99')*Fl)) /
          (4*(3-arb(5).sqrt()))).sqrt()
    phase_Hmu = abs(cnu)*L**2*L.sqrt()/arb(5).sqrt()
    phase_dual_upper_expression = MR/beta*phase_Hmu
    Knew = arb('1.960636035e70')
    bolt = arb('1.667e-9')
    phase_bolt_upper_expression = Knew*phase_dual_upper_expression
    ratio = phase_bolt_upper_expression/bolt
    assert 0 < L < 1 and arb('.16') < beta < arb('.17')
    assert arb('1e-8') < 1-chi < arb('1.179e-8')
    assert -arb('4e-37') < cnu < -arb('3e-37')
    assert phase_bolt_upper_expression > arb('1e35')
    assert ratio > arb('1e44')
    # Ordinary proof in the note: nonzero cnu with exponent 4-sigma<0
    # makes the ACTUAL first residual row fail L2(d log z), hence H3.
    assert Fraction(801, 200)-4 == Fraction(1, 200)
    constants = dict(lambda_fixed=lam, C=C, B2=B2, beta=beta, z_h=L,
                     chi=chi, full_trace_coefficient_lambda_B4=trace_coefficient,
                     signed_first_row_quartic_coefficient=cnu, MR=MR,
                     phase_Hmu_norm=phase_Hmu,
                     phase_dual_upper_expression=phase_dual_upper_expression,
                     phase_bolt_upper_expression=phase_bolt_upper_expression,
                     upper_expression_over_entire_bolt_allowance=ratio,
                     optimistic_full_dual_residual_threshold=bolt/Knew)
    return {'precision_bits': ctx.prec,
            'constants': {k: enc(v) for k, v in constants.items()},
            'scope': {
                'actual_reference_residual_fast_norm_infinite': True,
                'complete_separated_estimator_fails_from_known_term': True,
                'lower_bound_on_actual_residual_norm': False,
                'lower_bound_on_actual_physical_error': False,
                'new_signed_quantum_shape_enclosure': False,
                'validated_full_reference_trial': False,
                'original_physical_admission_closed': False,
                'construction_obstruction': False,
                'full_goal_complete': False},
            'unknown': {
                'complete_signed_shape_residual': None,
                'strong_common_regulator_remainder': None,
                'actual_output_pairings_and_application_errors': None,
                'original_X_completion_and_embeddings': None,
                'remaining_bolt_allowance': None,
                'all_other_coupled_entries': None}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args = parser.parse_args()
    output = ROOT / args.output
    if output.exists():
        raise RuntimeError('Immutable receipt exists; select a fresh output path.')
    protected = json.loads(BASELINE.read_bytes())
    verify = lambda: all(digest(Path(p)) == h for p, h in protected.items())
    assert verify()
    before = {p: digest(ROOT / p) for p in INPUTS}
    data = {'started_utc': datetime.now(timezone.utc).isoformat(),
            'goal_run': 'first turn of new execute-the-plan goal',
            'input_sha256_before': before,
            'ordinary_argument_utf8': (ROOT / INPUTS[0]).read_text(encoding='utf-8'),
            'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8'),
            'python_sha256': runtime['python_sha256']}
    try:
        data['result'] = execute()
        data['status'] = 'passed'
    except Exception as exc:
        data['status'] = 'failed'
        data['error'] = repr(exc)
        raise
    finally:
        data['input_sha256_after'] = {p: digest(ROOT / p) for p in INPUTS}
        data['inputs_unchanged'] = before == data['input_sha256_after']
        data['protected_paths_unchanged'] = verify()
        data['protected_path_count'] = len(protected)
        data['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    assert data['inputs_unchanged'] and data['protected_paths_unchanged']
    for name, value in data['result']['constants'].items():
        print(name+': '+value['ball'])
    print('Receipt: '+str(output))


if __name__ == '__main__':
    main()
