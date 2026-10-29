"""Actual coupled mass-functional bound; test, not certify, signed residuals."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_filtered_test_mass_trace as m
arb, ctx, dec, enc, sp = m.arb, m.ctx, m.dec, m.enc, m.sp
ROOT, BASELINE, runtime = m.ROOT, m.BASELINE, m.runtime
STEM = 'current-filtered-mass-adjoint-bound'
digest = lambda p: sha256(p.read_bytes()).hexdigest()


def reference_identities():
    sys.path.insert(0, str(ROOT))
    from operators.york_mixed_return import YorkMixedReturn
    mixed = YorkMixedReturn()
    c, z = mixed.core, mixed.z
    Q, Qp = c.Q[:2]
    directions = [[mixed.xi[0]/2+Q*mixed.J/z],
                  [-(mixed.br+3*mixed.xi[0]+Qp*mixed.J)/2]]
    for row in directions:
        for _ in range(2):
            row.append(mixed.D(row[-1]))
    H = c.build('reference')['strong']
    seed = c.seed_substitution()
    whole = sp.cancel((-z**4*sum(H[0][j][k]*directions[j][k]
                               for j in range(2) for k in range(3))).subs(seed))
    bref = 2*c.Fl*z**4*(6*c.beta+18*c.eps*z**2-15*c.eps*z+2*c.eps+2/c.C)
    assert sp.cancel(whole-bref*mixed.xi[0]) == 0
    Qseed = seed[Q]
    V = 24+12*c.eps*z**5/c.beta-24*c.eps*z**6/c.beta
    expected = [c.Fl*V, -6*c.Fl*(z*sp.diff(Qseed,z)-3*Qseed)/c.beta,
                -6*c.Fl*Qseed/c.beta, sp.Integer(0), sp.Integer(0)]
    assert all(sp.cancel((z**4*H[0][0][k]/c.beta).subs(seed)-expected[k]) == 0
               for k in range(5))
    # Exact abstract block algebra; no physical finite-mode matrix is used.
    t = sp.Symbol('t', nonnegative=True)
    assert sp.simplify((t/(1+t))*(t/(1+t)**2)-t**2/(1+t)**3) == 0
    assert sp.simplify(sp.diff(t/(1+t)**2,t)-(1-t)/(1+t)**3) == 0
    return dict(whole_reference_mixed=sp.factor(whole).__str__(),
                normalized_scalar_euler_coefficients=[str(sp.factor(v)) for v in expected],
                mixed_derivatives_and_York_returns_cancel=True,
                scalar_all_coefficients_agree=True)


def execute():
    ctx.prec = 768
    identities = reference_identities()
    params = {k: dec(v) for k,v in m.result('current-global-green-trial-enclosure')['trial_parameters'].items()}
    massp = {k: dec(v) for k,v in m.result('current-global-mass-conversion')['parameters'].items()}
    oldmass = m.result('current-complete-weak-mass-residual')
    prior = m.result('current-filtered-j-regularity-test-attempt2')
    transfer = m.constants('current-filtered-test-mass-trace-attempt2')
    green = m.constants('current-filtered-mixed-green-residual')
    beta, Cgeom, B2 = params['beta'], params['C'], params['B2']
    eps, Fl = massp['eps'], massp['Fl']
    Fmin = 1-arb('1.179e-8')
    assert Fl > Fmin and Fl < 1
    assert arb('.16') < beta < arb('.17')
    assert params['L'] < 1
    Braw = 12+70*abs(eps)/beta+4/(Cgeom*beta)
    Vdelta = 36*abs(eps)/beta
    assert Braw < arb('12.001') and 24-Vdelta > arb('23.99') and 24+Vdelta < arb('24.001')
    B0 = arb('1.000001')
    # B0 is the exact declared guard in the frozen ordinary proof/checker.
    # Decimal receipt re-encoding encloses that constant with a small radius;
    # comparing its upper endpoint to the same exact constant is invalid.
    assert dec(oldmass['joint_supremum']).overlaps(B0)
    assert "B0=arb('1.000001')" in (ROOT/'checks/current_complete_weak_mass_residual.py').read_text()
    R = arb('2.1e-5')
    assert dec(oldmass['complete_adjoint_residual_upper']) < R
    Jref = arb('12.001')*8*B0
    scalar_L2 = Jref/(Fmin*arb('23.99'))
    reference_energy = Jref/(Fmin*arb('23.99')).sqrt()
    laplacian = (arb('24.001')*scalar_L2+Jref/Fmin)/6
    ratio1 = arb('6.0001')/(6*Fmin)
    ratio2 = arb('24.001')/(arb('23.99')*Fmin)
    assert ratio2 > ratio1
    Sref = (ratio2*reference_energy**2+6*arb('2.946e-9')*laplacian**2).sqrt()
    q = arb('.338')
    Klo = arb('1.9999999764')
    Kf = 2*q.sqrt()*(Sref+R)/(3*Klo)
    ef = 2*q.sqrt()*R/(3*Klo)
    cf = q/(1+q)**2
    Kg, eg = cf*Kf, cf*ef
    assert Sref < 20 and Kg < 1 and eg < arb('1e-6')

    # Physical lift is unchanged. Only its fixed inner part is needed to
    # lower-test the positive expression of this specified upper estimator.
    coeff = [beta, arb(0), 1/Cgeom, -(4*beta+2/Cgeom+massp['tau']),
             3*beta+1/Cgeom+eps, -2*eps, eps]
    zc = arb('.1')
    Qt = sum((abs(v)*zc**k for k,v in enumerate(coeff) if k), arb(0))
    assert beta-Qt > arb('.15') and beta+Qt < arb('.17')
    Qpt = sum((k*abs(v)*zc**(k-1) for k,v in enumerate(coeff) if k), arb(0))
    assert Qpt < arb('.3')
    bm_lower = arb('.16')/arb('.17')-arb('.17')*arb('.3')*zc/(4*arb('.15')**2)
    assert bm_lower > arb('.88')
    Hmin = arb('1.32')*(zc**7/7).sqrt()
    Srec = dec(prior['constants']['actual_zero_derivative_reconstruction_S0'])
    L0, L1 = [dec(p['physical_profile_L2_upper']) for p in prior['whole_spectrum_profile_bounds'][:2]]
    profile_combination = Srec*(q*L0+L1)
    alpha_q = 2*Cgeom/(B2*beta)  # G=1, same as the accepted physical action.
    glue = arb('6e24')
    floor = Kg*alpha_q*glue*Hmin*profile_combination
    target = arb('1.49e34')
    assert floor/target > arb('1e5')
    ceiling = target/(Kg*alpha_q*Hmin*profile_combination)
    eps_tail = transfer['retained_full_form_error']
    bridge = transfer['full_to_retained_j_energy_error']
    full_ceiling = (target-bridge)/(1+eps_tail).sqrt()
    vals = dict(reference_mixed_coefficient_upper=Braw,
        reference_potential_deviation_upper=Vdelta,
        reference_scalar_source_Hmu_upper=Jref,
        reference_scalar_Hmu_upper=scalar_L2,
        reference_scalar_energy_upper=reference_energy,
        reference_scalar_laplacian_Hmu_upper=laplacian,
        actual_scalar_reference_energy_upper=Sref,
        complete_coupled_adjoint_residual_upper_used=R,
        mass_F_functional_norm_upper=Kf,
        mass_F_functional_reference_error_upper=ef,
        actual_filtered_mass_functional_norm_upper=Kg,
        filtered_mass_functional_reference_error_upper=eg,
        original_Q_inner_lower=beta-Qt, original_Q_inner_upper=beta+Qt,
        original_Qprime_inner_upper=Qpt, mass_profile_inner_lower_coefficient=bm_lower,
        original_normalized_mass_shape_L2_lower=Hmin,
        quantum_action_conversion_alpha_q=alpha_q,
        complete_structured_profile_combination=profile_combination,
        separated_gluing_estimator_expression_lower=floor,
        separated_gluing_expression_over_complete_target=floor/target,
        separated_gluing_expression_original_bolt_lower=green['pref_t']*floor,
        optimistic_gluing_coefficient_ceiling_ALL_other_errors_zero=ceiling,
        current_gluing_coefficient_over_optimistic_ceiling=glue/ceiling,
        unchanged_full_response_bridge=bridge,
        conditional_full_E_error_ceiling=full_ceiling,
        conditional_mass_residual_coefficient_ceiling_ALL_other_errors_zero=full_ceiling/Kg)
    assert all(v.is_finite() and v > 0 for v in vals.values())
    return dict(constants={k:enc(v) for k,v in vals.items()}, exact_identities=identities,
        signed_r0_mass_pairing=None, signed_r1_mass_pairing=None,
        actual_signed_mass_vector_norm=None, actual_Dj=None, actual_signed_Zr=None,
        scope=dict(complete_adjoint_transferred_to_filtered_mass_functional=True,
            original_mass_normalization_and_York_return_retained=True,
            specified_additive_global_gluing_estimator_rejected=True,
            actual_gluing_or_mass_term_lower_bounded=False,
            full_signed_mass_term_certified=False, actual_signed_response_certified=False,
            retained_boundary_trace_transferred=False, available_physical_margin=False,
            original_P32_defect=False, actual_kappa_m=False, complete_original_admission=False,
            nonlinear_existence=False, exact_compatibility=False, physical_matching=False,
            semiclassical_validity_noise=False, all_sector_stability=False,
            construction_obstruction=False, full_goal_complete=False))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--receipt', default=STEM)
    args = ap.parse_args()
    out = ROOT/'evidence'/(args.receipt+'.json')
    assert not out.exists(), 'Preserve immutable receipts'
    old = json.loads((ROOT/'evidence/current-filtered-test-mass-trace-attempt2.json').read_bytes())
    inputs = list(dict.fromkeys([
        'notes/'+STEM+'.md', 'checks/current_filtered_mass_adjoint_bound.py',
        'evidence/current-filtered-test-mass-trace-attempt2.json',
        'evidence/current-filtered-test-mass-trace-attempt2-verification.json',
        'notes/current-coupled-input-energy-bound.md',
        'notes/current-reference-clock-adjoint-trial.md',
        'notes/current-complete-weak-mass-residual.md',
        'notes/current-coupled-mass-adjoint.md',
        'notes/stationary-curved-quantum-tt.md',
        'evidence/current-complete-weak-mass-residual.json',
        'evidence/current-complete-weak-mass-residual-review.json',
        'evidence/current-complete-weak-mass-residual-runtime.json',
        *old['input_sha256_after']]))
    before = {p:digest(ROOT/p) for p in inputs}
    protected = json.loads(BASELINE.read_bytes())
    rec = dict(started_utc=datetime.now(timezone.utc).isoformat(), input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        assert all(digest(ROOT/p)==h for p,h in old['input_sha256_after'].items())
        assert all(digest(Path(p))==h for p,h in protected.items())
        assert digest(Path(sys.executable)) == runtime['python_sha256']
        for pkg in runtime['packages'].values():
            assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
        oldres = json.loads((ROOT/'evidence/current-complete-weak-mass-residual.json').read_bytes())['runs'][-1]
        assert oldres['status']=='passed'
        assert all(digest(ROOT/p)==h for p,h in oldres['input_sha256_after'].items())
        rec['result'] = execute()
        after = {p:digest(ROOT/p) for p in inputs}
        assert before == after
        assert all(digest(Path(p))==h for p,h in protected.items())
        rec.update(status='filtered_mass_functional_bounded_additive_gluing_estimator_rejected',
            input_sha256_after=after, protected_paths_unchanged=len(protected),
            python_sha256=runtime['python_sha256'],
            pinned_packages={k:len(v['files']) for k,v in runtime['packages'].items()})
    except Exception:
        rec.update(status='failed', traceback=traceback.format_exc())
    rec['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'], error=rec.get('traceback'),
        constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
    if rec['status']=='failed':
        raise SystemExit(1)
