"""Enclose the ordinary fixed-column estimate; not the full Ritz residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-quartic-column-energy-dual'
RUNTIME = ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/' + STEM + '.md', 'checks/current_quartic_column_energy_dual.py',
    'notes/current-disk-heat-propagation.md',
    'notes/current-local-boundary-heat-remainder.md',
    'notes/current-smeared-heat-projection-domain.md',
    'notes/current-structured-weak-response-attempt.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/current-galerkin-range-audit.md',
    'evidence/current-disk-heat-propagation.json',
    'evidence/current-local-boundary-heat-remainder.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-global-mass-conversion.json',
    'evidence/current-complete-polynomial-weak-transfer.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-quartic-column-energy-dual.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def dec(value):
    d = value['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad'])) * arb(10) ** int(d['exponent'])


def enc(x):
    m, r, e = x.mid_rad_10exp(60)
    return {'ball': x.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def prior(stem):
    r = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())['runs'][-1]
    assert r['status'] == 'passed'
    return r


def padd(*polys):
    out = [arb(0)] * max(map(len, polys))
    for p in polys:
        for j, c in enumerate(p):
            out[j] += c
    return out


def scale(p, c):
    return [a*c for a in p]


def mul(a, b):
    c = [arb(0)] * (len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i+j] += ai*bj
    return c


def derivative(p):
    return [j*p[j] for j in range(1, len(p))] or [arb(0)]


def l2(p):
    # The positive integral over [0,L] is <= this integral over [0,1].
    integral = sum((c/(j+1) for j, c in enumerate(mul(p, p))), arb(0))
    assert integral > 0
    return integral.sqrt()


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, expected in package['files'].items():
            assert digest(Path(package['root']) / name) == expected
    energy = prior('current-coupled-truncation-energy-error')['result']
    beta, C, B2 = [dec(energy['parameters'][k]) for k in ['beta', 'C', 'B2']]
    CL, CI = [dec(energy[k]) for k in ['global_shape_L2_constant',
                                      'global_shape_supremum_constant']]
    epsilon = dec(prior('current-complete-polynomial-weak-transfer')['result']['constants']['epsilon'])
    assert 0 < epsilon < arb('1.227e-74')
    CLt, CIt = CL/(1-epsilon).sqrt(), CI/(1-epsilon).sqrt()
    qt = arb('.338')/(1-epsilon)
    alpha1 = 3*qt*CLt/(arb('.99')**(arb(9)/2))
    p = prior('current-global-mass-conversion')['result']['parameters']
    u, ep, tau = [dec(p[k]) for k in ['u', 'eps', 'tau']]
    Q = [beta, arb(0), u, -4*beta-2*u-tau,
         3*beta+u+ep, -2*ep, ep]
    R = scale(derivative(derivative(Q)), -1)
    v = [arb(0)]*4+[arb(-1)]
    delta = lambda a: derivative(mul(Q, derivative(a)))
    H = lambda a: padd(delta(a), mul(R, a))
    dv = delta(v)
    O6 = padd(scale(padd(H(delta(H(v))), mul(delta(R), H(v)),
                        delta(mul(v, delta(R)))), -arb(1)/105),
              scale(padd(scale(mul(mul(mul(R, R), R), v), 16),
                         scale(mul(mul(R, R), dv), 24),
                         scale(delta(mul(mul(R, R), v)), 24),
                         scale(delta(mul(R, dv)), 24)), -arb(1)/1008))
    Oprod = padd(O6, scale(v, arb(4)/63))
    A6 = l2(Oprod)/(4*arb.pi())
    A6can = l2(O6)/(4*arb.pi())
    P4 = padd(scale(mul(mul(R, R), v), 4), scale(mul(R, dv), 8),
              scale(delta(mul(R, v)), 8), scale(delta(dv), 8))
    P4norm = l2(P4)
    boundary = prior('current-local-boundary-heat-remainder')
    whole = next(x for x in boundary['numerics'] if x['kind'] == 'whole_local_path_bound')
    disk = prior('current-disk-heat-propagation')
    assert disk['complete_actual_domain_remainder_bounded']
    pick = lambda name: next(x for x in disk['numerics'] if x['kind'] == name)
    cap = pick('conditional_complete_cap_kernel_table')['table']
    transport = pick('full_uncut_transport_norms')['tables']
    propagation = pick('complete_actual_domain_propagation_budget')
    tails = pick('complete_commutator_and_thermal_tails')
    thermal = dec(tails['full_thermal_image_guard'])
    t0 = arb('.001')
    Cbin = arb(0)
    for item in whole['terms']:
        np, nq = item['V_power'], item['t_power']
        numerator, denominator = item['coefficient'].split('/')
        coefficient = arb(numerator)/arb(denominator)
        assert np+2*nq >= 7 and coefficient > 0
        Cbin += coefficient*2**np*2*(1+arb(np)/2).gamma()*t0**(nq+arb(np-7)/2)
    Cbin *= arb('1.06')/arb.pi()
    Cbexit = 2816/arb.pi()*t0.exp()*t0**(-arb(7)/2)*(-3/(64*t0)).exp()
    assert 4*t0**(-arb(7)/2)*thermal < arb('1e-100')
    assert dec(propagation['new_boundary_D']) < 149
    bulk = 4*(dec(propagation['first_actual']['mixed'])+
              dec(propagation['higher_Volterra'])+
              dec(propagation['smooth_cap_comparison']['mixed']))+1/(120*arb.pi())
    assert bulk < 953
    Cb = Cbin+Cbexit+arb('1e-100')
    Cb2 = (149*Cb).sqrt()
    cap_at = lambda eta: 4*C*CLt/(B2*beta)*(A6*eta+arb(4)/7*Cb2*eta**(arb(7)/4)+arb(953)/2*eta**2)
    capbound = cap_at(t0)
    # The existing L-infinity route is finite but pays the original CI.
    A6L1 = A6  # L<=1, hence L1<=L2 on the enlarged unit interval.
    oldcap = 4*C*CIt/(B2*beta)*(A6L1*t0+arb(1101)/2*t0**2)
    A_cap = 2*(dec(cap['base']['inf'])+thermal)
    A_int = dec(transport['P']['base']['inf'])+dec(propagation['cutoff_jet_norm']['inf'])*t0**3
    A = arb(max(A_cap.upper(), A_int.upper()))
    ts = t0/2
    M = 2*(A/ts+(dec(propagation['first_actual']['base'])+
                  dec(propagation['higher_Volterra']))*ts**3)
    assert M < 734
    assert (arb(2000)**4+1)/C > 324
    CH, C0 = 64*M*CLt/3, 2*CLt/(3*arb.pi())
    C2, Cp = CLt*P4norm/(480*arb.pi()), CLt/(189*arb.pi())
    common = 4*C/(B2*beta*t0*324)
    raw, area, fourth = common*CH/t0, common*C0/t0**2, common*C2
    product = common*Cp*(t0+arb(1)/324)
    positive = raw+area+fourth+product
    massive = capbound+positive
    capcan = 4*C*CLt/(B2*beta)*(A6can*t0+arb(4)/7*Cb2*t0**(arb(7)/4)+arb(953)/2*t0**2)
    positive_can = raw+area+fourth
    massive_can = capcan+positive_can
    assembly = prior('current-canonical-heat-restoration-assembly')
    diagonal = next(x for x in assembly['numerics'] if x['kind']=='assembled_known_TT_entry' and x['i']==2 and x['j']==2)
    known_energy = -dec(diagonal['value'])/beta
    assert known_energy > 0
    etail = epsilon/(1-epsilon)
    c_lower, c_upper = arb(known_energy.lower()), arb(known_energy.upper())
    x_lower = (-massive_can+(massive_can**2+4*(1+etail)*c_lower).sqrt())/(2*(1+etail))
    x_upper = (massive_can+(massive_can**2+4*(1-etail)*c_upper).sqrt())/(2*(1-etail))
    assert x_lower > 0 and x_upper < 2
    alpha_exact = qt/x_lower
    tail_debit = etail*x_upper
    d_cut = alpha_exact*(capcan+tail_debit)
    d_unsampled = alpha_exact*(massive_can+tail_debit)
    inverse = json.loads((ROOT/'evidence/current-spatial-inverse-output-test.json').read_bytes())['result']['constants']
    W0 = dec(inverse['W0_complete_scalar_physical_upper'])
    W1 = dec(inverse['W1_complete_mixed_physical_upper'])
    ft = qt/(1+qt)
    def physical(debit):
        z = debit/ft
        return W0*(3*z+3*z*z+z*z*z)+W1*(2*z+z*z)
    oldpositive = common*(64*1521*CIt/t0+C0/t0**2+C2+Cp*(t0+arb(1)/324))
    target = arb('1e-30')
    constants = {'CLt': CLt, 'CIt': CIt, 'q_t': qt,
        'one_dimensional_Ritz_coefficient_envelope': alpha1,
        'A6_L2': A6, 'P4_L2': P4norm, 'boundary_Linfinity_collar_constant': Cbin,
        'boundary_Linfinity_exit_constant': Cbexit, 'boundary_L2_constant': Cb2,
        'bulk_pointwise_expression': bulk, 'heat_diagonal_bound_M': M,
        'cap_at_001': capbound, 'old_supremum_cap_expression': oldcap,
        'cap_estimator_improvement_factor': oldcap/capbound,
        'positive_raw': raw, 'positive_area': area, 'positive_fourth': fourth,
        'positive_product': product, 'positive_time_bound': positive,
        'complete_massive_remainder_column_bound': massive,
        'old_fixed_column_positive_expression': oldpositive,
        'positive_estimator_improvement_factor': oldpositive/positive,
        'one_dimensional_cap_residual_debit': alpha1*capbound,
        'one_dimensional_massive_residual_upper_expression': alpha1*massive,
        'one_dimensional_cap_over_sufficient_target': alpha1*capbound/target,
        'one_dimensional_massive_over_sufficient_target': alpha1*massive/target,
        'conditional_positive_column_error_target': (target-alpha1*capbound)/alpha1}
    constants.update({'canonical_A6_L2': A6can, 'canonical_cap_at_001':capcan,
        'canonical_positive_time_bound':positive_can,
        'canonical_complete_massive_column_bound':massive_can,
        'known_local_diagonal_energy':known_energy,
        'complete_psi_energy_lower':x_lower,'complete_psi_energy_upper':x_upper,
        'sharper_exact_one_dimensional_coefficient_envelope':alpha_exact,
        'retained_transport_column_error':tail_debit,
        'sharper_one_dimensional_cut_residual_debit':d_cut,
        'sharper_one_dimensional_unsampled_residual_expression':d_unsampled,
        'sharper_unsampled_residual_over_target':d_unsampled/target,
        'sharper_conditional_positive_column_error_target':(target-d_cut)/alpha_exact,
        'positive_error_enclosure_improvement_required':positive_can/((target-d_cut)/alpha_exact),
        'best_case_cut_bolt_error_expression':physical(d_cut),
        'best_case_unsampled_bolt_error_expression':physical(d_unsampled),
        'best_case_unsampled_bolt_over_entire_allowance':physical(d_unsampled)/arb('1.667e-9')})
    return {'constants': {k:enc(x) for k,x in constants.items()},
        'Q_coefficients': [enc(x) for x in Q],
        'O6prod_coefficients': [enc(x) for x in Oprod],
        'P4_coefficients': [enc(x) for x in P4],
        'actual_F_n': None, 'full_residual_bound': None, 'signed_quantum_center': None,
        'multicolumn_coefficient_norms': None, 'actual_numeric_application_errors': None,
        'scope': {'full_energy_dual_fixed_quartic_column': True,
            'fixed_column_common_regulator_removal': True,
            'all_angular_and_thermal_sums_retained': True,
            'signed_column_evaluated': False, 'complete_quantum_action_evaluated': False,
            'complete_Ritz_residual_certified': False, 'physical_budget_closed': False,
            'nonlinear_neighborhood_transfer': False, 'construction_obstruction': False,
            'existence_or_stability_proved': False}}


if __name__ == '__main__':
    path = ROOT / 'evidence' / (STEM+'-assembly.json')
    assert not path.exists(), 'Keep receipts immutable.'
    before = {p:digest(ROOT/p) for p in INPUTS}
    receipt = {'started_utc': datetime.now(timezone.utc).isoformat(),
        'input_sha256_before': before, 'python': sys.executable,
        'ordinary_argument_utf8': (ROOT/INPUTS[0]).read_text(encoding='utf-8'),
        'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result'] = execute()
        receipt['status'] = 'column_bound_arithmetic_passed_full_residual_unresolved'
    except BaseException as exc:
        receipt['status'] = 'failed'
        receipt['failure'] = repr(exc)
        raise
    finally:
        receipt['input_sha256_after'] = {p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        path.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    assert before == receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'], 'constants':{
        k:v['ball'] for k,v in receipt['result']['constants'].items()}}, indent=2))
