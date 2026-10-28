"""Interval substitutions for the ordinary full-source/graph argument, not an inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-structured-source-graph-test'
INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_structured_source_graph_test.py',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-finite-local-action-assembly.md',
    'notes/strong-effective-tt-source.md',
    'notes/weighted-york-source-conversion.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-preconditioned-source-attempt.md',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-preconditioned-source-attempt.json',
    'evidence/current-quartic-ritz-rank-test-final.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
runtime = json.loads((ROOT/INPUTS[-1]).read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def dec(x):
    d = x['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad']))*arb(10)**int(d['exponent'])


def enc(x):
    m, r, e = x.mid_rad_10exp(60)
    return {'ball': x.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def main():
    ctx.prec = 512
    out = ROOT/'evidence'/(STEM+'.json')
    assert not out.exists(), 'Preserve immutable execution evidence.'
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for package in runtime['packages'].values():
        for name, h in package['files'].items():
            assert digest(Path(package['root'])/name) == h
    protected = json.loads(Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json').read_bytes())
    assert all(digest(Path(p)) == h for p, h in protected.items())
    before = {p: digest(ROOT/p) for p in INPUTS}
    started = datetime.now(timezone.utc).isoformat()
    result = {'started_utc': started, 'status': 'started',
              'input_sha256_before': before,
              'ordinary_argument_utf8': (ROOT/INPUTS[0]).read_text(encoding='utf-8'),
              'arithmetic_source_utf8': Path(__file__).read_text(encoding='utf-8')}
    try:
        geo = json.loads((ROOT/'evidence/current-geodesic-two-sided-remainder.json').read_bytes())['runs'][-1]
        prev = json.loads((ROOT/'evidence/current-preconditioned-source-attempt.json').read_bytes())
        rank = json.loads((ROOT/'evidence/current-quartic-ritz-rank-test-final.json').read_bytes())
        assert geo['status'] == rank['status'] == 'passed'
        assert prev['status'] == 'accounting_checked_complete_source_application_unresolved'
        p = next(x for x in geo['numerics'] if x['kind'] == 'actual_finite_seed_two_sided_parameters')
        rc = rank['result']['constants']
        lam, Fl, B2, C = [dec(p[k]) for k in ['lambda_value', 'F_lambda', 'B2', 'C']]
        beta = dec(rc['beta'])
        assert (lam*B2/C).overlaps(beta)
        assert lam*B2 > arb('8e9')
        assert abs(dec(p['guards']['finite_seed_d'])) < arb('8.6e-19')
        thermal = dec(p['guards']['thermal_delta'])
        assert thermal > 0 and thermal < 1
        e = lam/(60*arb.pi())
        d = 96*arb.pi()*arb('1e22')*lam-e
        Cv = 1+11*e/3
        assert d > 0 and d < arb('2.946e-9')
        assert Fl > 1-arb('1.179e-8')
        q, eps = [dec(rc[k]) for k in ['q_t', 'epsilon']]
        eq = dec(prev['constants']['combined_Et_reference_error'])
        assert eq < arb('.0022')
        h0, h1, h2 = [1/x.sqrt() for x in [arb('23.99'), arb('5.99'), 6*d]]
        hH = h2+4*h1
        # Covariant tensor/projector and transpose coefficient arithmetic.
        proj2 = (4*arb(16)**2+4*arb(11)**2+24*arb(5)**2).sqrt()
        j_hess = 3000+2*30*9+arb('3.5')*600
        c_hess = 7*3000+2*43*108+arb('5.1')*7680
        c2 = arb(7)/12*6*40+arb(11)/6*14
        c1 = arb(7)/12*12*10000+arb(11)/6*200
        c0 = arb(7)/12*(600000+12320)+arb(11)/6*(20000+(120+26+arb('6.01'))*7+98)
        assert proj2 < 300 and j_hess < 10000 and c_hess < 100000
        assert c2 < 200 and c1 < 80000 and c0 < 1000000
        v = dec(p['guards']['Gq2_over_lambdaB4'])/(60*arb.pi())
        main = arb('12.04')*(Cv**2/arb('23.99')+d/6).sqrt()
        curv = arb(22).sqrt()*e*(200*hH+80000*h1+1000000*h0)
        w = arb('12.04')*arb(11)/36*e*arb('1.1e-27')*h0
        hw = arb(22).sqrt()*d*arb('4.3e-25')*h0
        mag = 20*v*h0
        ret = arb(18).sqrt()*arb('9e-28')*h1
        rest = curv+w+hw+mag+ret
        kb = main+rest
        cb = kb*q.sqrt()
        ch = (1-eps)*arb('.99')*Fl
        target = arb('5e-32')
        cutoff = cb**2/((1-eq)*target**2)
        rho_to_a = 1/(1-eq)
        strong_to_a = 1/((1-eq)*ch.sqrt())
        beta_target = target/strong_to_a
        # A12 is the already accepted ordered application multiplier.
        a12 = dec(rc['source_application_error_multiplier'])
        debit_multiplier = (1+q)*a12*strong_to_a
        safe_beta = arb('4.964e-32')
        safe_a = strong_to_a*safe_beta
        debit = debit_multiplier*safe_beta
        assert rest < arb('2e-27')
        assert kb < arb('2.459') and cb < arb('1.43')
        assert kb**2 < arb('6.047')
        assert cutoff < arb('8.2e62') and cutoff > arb('8e62')
        assert safe_beta < beta_target and safe_a < target
        assert debit < arb('6.705e-32')
        values = dict(lambda_value=lam, d=d, C_v=Cv,
            physical_L2_equals_existing_Hmu_factor=lam*B2/(C*beta),
            scalar_Hessian_constant=hH,
            mixed_main=main, mixed_curvature=curv, mixed_w=w,
            mixed_Hess_w=hw, mixed_magnetic=mag, mixed_return=ret,
            mixed_remainders_total=rest, complete_mixed_A_to_H_constant=kb,
            complete_B_H_to_H_constant=kb**2,
            complete_B_Et_to_H_and_S_operator_graph_constant=cb,
            original_Et_to_H_constant=1/ch.sqrt(), c_H=ch,
            actual_reference_comparison=eq,
            positive_A_Q_H_gap=(1-eq)*ch,
            generic_graph_tail_prefactor=cb/(1-eq).sqrt(),
            generic_sufficient_spectral_cutoff_NOT_executed=cutoff,
            energy_dual_residual_to_application_multiplier=rho_to_a,
            strong_residual_to_application_multiplier=strong_to_a,
            sufficient_strong_residual_target_NOT_attained=beta_target,
            safe_strong_residual_target_NOT_attained=safe_beta,
            strong_residual_to_COMPLETE_residual_multiplier=debit_multiplier,
            conditional_source_application_error=safe_a,
            conditional_source_debit_in_complete_residual=debit,
            unchanged_conditional_complete_residual=dec(rc['conditional_complete_residual_upper']))
        result.update(status='source_graph_bound_checked_complete_application_unresolved',
            constants={k: enc(v) for k, v in values.items()},
            complete_source_strong_L2_bound=True,
            structured_associated_operator_graph_bound=True,
            spatial_fourth_flow_derivative_constant=None,
            numerical_spectral_projector=None,
            actual_source_strong_residual=None,
            actual_source_application_error=None,
            actual_K_application_error=None,
            complete_numerical_residual=None,
            physical_margin_closed=False, goal_complete=False,
            formal_verification=False)
        for k, val in values.items():
            print(k+': '+val.str(28))
    except Exception as exc:
        result.update(status='failed', exception_type=type(exc).__name__)
        raise
    finally:
        after = {p: digest(ROOT/p) for p in INPUTS}
        result.update(input_sha256_after=after,
            inputs_unchanged=before == after,
            protected_paths_unchanged=all(digest(Path(p)) == h for p, h in protected.items()),
            protected_path_count=len(protected),
            python_sha256=runtime['python_sha256'],
            finished_utc=datetime.now(timezone.utc).isoformat())
        out.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    assert result['inputs_unchanged'] and result['protected_paths_unchanged']
    print('Full strong source and associated graph bounded; implemented residual remains unknown.')


if __name__ == '__main__':
    main()
