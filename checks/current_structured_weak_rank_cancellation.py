"""Full structured weak residual lower bounds; not an inverse implementation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import argparse
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-structured-weak-rank-cancellation'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = ['notes/'+STEM+'.md', 'checks/current_structured_weak_rank_cancellation.py',
          'notes/current-strong-residual-rank-test.md',
          'checks/current_strong_residual_rank_test.py',
          'evidence/current-strong-residual-rank-test-final.json',
          'evidence/current-structured-source-graph-test.json',
          'evidence/current-geodesic-two-sided-remainder.json',
          'evidence/current-canonical-heat-restoration-assembly.json',
          'evidence/current-spatial-inverse-output-test.json',
          'evidence/current-coupled-mass-adjoint-runtime.json']


def sha(p):
    return sha256(Path(p).read_bytes()).hexdigest()


def execute():
    # The imported file only defines helpers after loading its pinned runtime;
    # its scientific execution is guarded by __main__ and is not invoked.
    spec = importlib.util.spec_from_file_location(
        'accepted_rank_helpers', ROOT/'checks/current_strong_residual_rank_test.py')
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    arb, ctx, dec, enc = old.arb, old.ctx, old.dec, old.enc
    ctx.prec = 512
    assert sha(sys.executable) == old.runtime['python_sha256']
    for pkg in old.runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert sha(Path(pkg['root'])/name) == expected
    accepted = json.loads((ROOT/INPUTS[4]).read_bytes())
    assert accepted['status'] == 'finite_output_estimator_test_executed_full_application_unresolved'
    assert all(sha(ROOT/p) == h for p, h in accepted['input_sha256_after'].items())
    geo = old.receipt('current-geodesic-two-sided-remainder')
    ep = old.pick(geo, 'actual_seed_correlated_endpoint')
    gp = old.pick(geo, 'actual_finite_seed_two_sided_parameters')
    rc = old.pick(old.receipt('current-canonical-heat-restoration-assembly'),
                  'common_restoration_coefficients')
    sc = old.receipt('current-structured-source-graph-test',
                    'source_graph_bound_checked_complete_application_unresolved')['constants']
    sp = old.receipt('current-spatial-inverse-output-test',
                    'error_estimator_test_passed_signed_action_unresolved')['result']['constants']
    add, mul, diff, shift, scale = old.add, old.mul, old.diff, old.shift, old.scale
    fa, fm, fd, fs = old.fadd, old.fmul, old.fdiff, old.fscale
    Q = {k: dec(v) for k, v in zip([0, 2, 3, 4, 5, 6], ep['Q_coefficients'])}
    C, B2, lam, q = (dec(gp[k]) for k in ['C', 'B2', 'lambda_value', 'q'])
    beta, eps = Q[0], dec(sp['epsilon'])
    assert beta.overlaps(lam*B2/C)
    d, Cv, Krest = (dec(sc[k]) for k in ['d', 'C_v', 'mixed_remainders_total'])
    qt = arb('.338')/(1-eps)
    ds, xx = Q[6]*C, lam*B2
    tau = {4: (1+3*xx+ds)/xx, 5: -5*ds/xx, 6: arb(9)/2*ds/xx}
    drift = add(diff(Q), scale(shift(Q, -1), -2))
    dtau = scale(shift(add(mul(Q, diff(diff(tau))), mul(drift, diff(tau))), 2), 1/beta)
    m = diff(shift(diff(Q), -4))
    J = [scale(diff(m), 2), m]
    K = fm(J, Q)
    wfield = [shift(mul(m, m), 4)]
    xi = fa(fd(K), fs(wfield, -1))
    nu = fa(fs(xi, arb(1)/2), [shift(p, -1) for p in K])
    T = fa(fs(xi, -1), fs(fd(K), -arb(1)/2))
    coefficient = add(scale(shift(Q, -2), 24), scale(shift(diff(Q), -1), -6))
    dfield = fa(fm(fd(fd(wfield)), Q),
                fm(fd(wfield), add(diff(Q), scale(shift(Q, -1), -6))),
                fm(fa(fd(K), fs([shift(p, -1) for p in K], -1)), coefficient))
    assert len(dfield) == 3
    # Ordinary argument (5) proves the exact identity; this supports its algebra.
    sumfield = fa(T, fs(nu, 3))
    direct = fa(fs(fa(fm(fd(fd(sumfield)), Q), fm(fd(sumfield), diff(Q))), -2),
                fm(fd(nu), scale(shift(Q, -1), 12)))
    difference = fa(direct, fs(dfield, -1))
    assert all(v.contains(0) for p in difference for v in p.values())
    Dseed = add(scale(diff(diff(Q)), -1), scale(shift(diff(Q), -1), 6),
                scale(shift(Q, -2), -12))
    f = []
    for rational in accepted['result']['exact_bump_derivative_squared_norms']:
        v = Fraction(rational)
        f.append((arb(v.numerator)/v.denominator).sqrt())
    matrices = {name: [[dec(v) for v in row] for row in mat]
                for name, mat in accepted['result']['local_coefficient_absolute_matrices'].items()}
    fields = {'xi': xi, 'nu': nu, 'T': T, 'dfield': dfield}
    derivatives = {}
    for name, co in fields.items():
        derivatives[name] = []
        for p in co:
            row = [p]
            for j in range(3):
                row.append(diff(row[-1]))
            derivatives[name].append(row)
    xiprime = fd(xi)
    alphaR = 8*arb.pi()*arb('1e22')/(C*B2)
    s2 = 4*C*C/q
    assert 2*C*C*(q+2)/(q+1)**2+2*C*C/(q+1) < s2
    ranks = [1000000, 10000000, 100000000]
    results = []
    for rank in ranks:
        M = rank+1
        a, b = arb('.4'), arb('.6')
        h = (b-a)/M
        omega = 2/h
        b0s, energy, oldenergy, ledgers = [], [], [], []
        for j in range(128):
            # Expand a cover of cell centers to cover each entire physical cell.
            lo = a if j == 0 else a+(b-a)*j/128-h/2
            hi = b if j == 127 else a+(b-a)*(j+1)/128+h/2
            interval = lo.union(hi)
            sup = lambda p: old.value(p, interval).abs_upper()
            inf = lambda p: old.value(p, interval).lower()
            cb = {name: [[sup(p) for p in row] for row in rows]
                  for name, rows in derivatives.items()}
            def nk(name, k):
                return sum((comb(k, jj)*cb[name][i][jj]*omega**(i+k-jj)*f[i+k-jj]
                            for i in range(3) for jj in range(k+1)), arb(0))
            D0 = inf(xi[2])*omega**2*f[2]-sup(xi[1])*omega*f[1]-sup(xi[0])*f[0]
            D1 = inf(xi[2])*omega**3*f[3]-sum(
                (sup(xiprime[i])*omega**i*f[i] for i in range(3)), arb(0))
            assert D0 > 0 and D1 > 0
            jets = {name: [hi**2*nk(name, k)/D0 for k in range(4)]
                    for name in ['xi', 'nu', 'T']}
            qmax, qmin, q1, q2 = sup(Q), inf(Q), sup(diff(Q)), sup(diff(diff(Q)))
            L1 = qmax.sqrt()*jets['xi'][1]/(beta.sqrt()*lo)
            L2 = (qmax*jets['xi'][2]+sup(drift)*jets['xi'][1])/beta
            glow = qmin*lo**4/(beta*hi**2)*(D1/nk('xi', 0))**2
            Aup = 6*d*L2**2+arb('6.01')*L1**2+arb('24.01')
            pm = 4*(Cv*inf(tau)+d*inf(tau)*glow-d*sup(dtau)/2)
            assert pm > 0
            b0s.append((pm**2/Aup).lower())
            v = jets['nu'][:3]+jets['T'][:3]
            local = {name: sum((mat[ii][jj]*v[ii]*v[jj] for ii in range(6)
                               for jj in range(6)), arb(0))/beta
                     for name, mat in matrices.items()}
            lt0, lt1, lt2, lt3 = jets['T']
            ld = hi**2*nk('dfield', 0)/D0
            zz = sup(Dseed)
            newR2 = alphaR/beta*(16*lt0**2+(4*C+2*C*C*zz)*12*qmax*jets['nu'][1]**2
                +C*C*(4*zz**2*lt0**2+8*zz*lt0*ld+2*ld**2))
            lap = qmax*lt2+q1*lt1
            ea2 = dec(rc['beta2']).abs_upper()*(4*q2*q2*lt0**2+16*q2*lt0*lap+8*lap*lap)/beta
            elll = dec(rc['LLL_Hessian_coefficient']).abs_upper()*qmax*lt1**2/beta
            evac = dec(rc['vacuum_coefficient_upper'])*lt0**2/beta
            aa, dd = qmax*lt3, arb('1.5')*q1*lt2+q2*lt1/2
            egeo = 2*s2/(525*arb.pi()*B2*beta)*qmax*(arb(5)/8*aa*aa+arb('1.5')*aa*dd+dd*dd)
            rest = local['Einstein']+local['magnetic']+local['Weyl']+ea2+elll+evac+egeo
            pref = (1+eps)/(1-arb('.0022'))
            energy.append((pref*(rest+newR2)).upper())
            oldenergy.append((pref*(rest+local['R2'])).upper())
            ledgers.append({'D0_lower': D0.lower(), 'D1_lower': D1.lower(),
                            'R2_before': local['R2'].upper(), 'R2_after': newR2.upper(),
                            'local_Einstein': local['Einstein'].upper(),
                            'complete_nonlocal': egeo.upper(), 'complete_other':
                            (local['magnetic']+local['Weyl']+ea2+elll+evac).upper()})
        b0 = min(b0s)
        assert b0.sqrt() > Krest
        bstar = (b0.sqrt()-Krest)**2
        eup, oldup = max(energy), max(oldenergy)
        mu = bstar/eup
        fullscreen, innerscreen = arb('1e-30'), arb('4.9887427328e-32')
        vals = {'rank': arb(rank), 'cell_width': h, 'main_B_lower': b0,
                'complete_B_lower': bstar, 'complete_Et_energy_upper': eup,
                'independent_curvature_jet_energy_upper_same_new_source': oldup,
                'energy_upper_expression_improvement': oldup/eup,
                'complete_original_residual_lower': mu,
                'complete_original_residual_lower_before_curvature_cancellation': bstar/oldup,
                'ratio_to_full_sufficient_screen': mu/fullscreen,
                'ratio_to_inner_sufficient_screen': mu/innerscreen,
                'complement_norm_necessary_if_positive': (mu-fullscreen)/(1+qt)}
        vals.update({'maximum_'+key: max(row[key] for row in ledgers) for key in
                     ['R2_before', 'R2_after', 'local_Einstein', 'complete_nonlocal', 'complete_other']})
        result = {'fixed_output_dimension': rank, 'cover_intervals': 128,
                  'constants': {k: enc(v) for k, v in vals.items()},
                  'full_original_residual_screen_rejected': bool(mu > fullscreen),
                  'structured_inner_weak_screen_rejected': bool(mu > innerscreen)}
        results.append(result)
        print(json.dumps({'rank': rank, 'mu_lower': mu.str(18),
                          'full_screen_rejected': result['full_original_residual_screen_rejected'],
                          'inner_screen_rejected': result['structured_inner_weak_screen_rejected']}), flush=True)
    return {'precision_bits': 512, 'rank_tests': results,
            'complete_quantum_energy_retained': True,
            'actual_domain_original_weak_residual_test': True,
            'numerical_S_K_implemented': False, 'actual_upper_residual_certificate': None,
            'physical_gate_closed': False, 'goal_complete': False,
            'construction_obstruction': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args = parser.parse_args()
    output = ROOT/args.output
    if output.exists():
        raise SystemExit('Refusing to replace an immutable receipt.')
    protected = json.loads(BASELINE.read_bytes())
    before = {p: sha(ROOT/p) for p in INPUTS}
    receipt = {'started_utc': datetime.now(timezone.utc).isoformat(),
               'input_sha256_before': before,
               'ordinary_argument_utf8': (ROOT/INPUTS[0]).read_text(encoding='utf-8'),
               'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8'),
               'previous_turn': 'no progress; resumed after implementation impasse'}
    try:
        assert all(sha(p) == h for p, h in protected.items())
        receipt['result'] = execute()
        receipt['status'] = 'original_weak_rank_test_executed_full_solver_unresolved'
    except Exception as exc:
        receipt['status'] = 'failed'
        receipt['error'] = repr(exc)
        raise
    finally:
        receipt['input_sha256_after'] = {p: sha(ROOT/p) for p in INPUTS}
        receipt['inputs_unchanged'] = before == receipt['input_sha256_after']
        receipt['protected_paths_unchanged'] = all(sha(p) == h for p, h in protected.items())
        receipt['protected_path_count'] = len(protected)
        receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8', newline='\n')
    assert receipt['inputs_unchanged'] and receipt['protected_paths_unchanged']
