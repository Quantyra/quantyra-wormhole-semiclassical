"""Actual structured L1 return and its contribution to the open mass budget."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import json
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_filtered_test_mass_trace as m
import current_source_adapted_fractional_return as fr
arb, ctx, dec, enc = m.arb, m.ctx, m.dec, m.enc
ROOT, BASELINE, runtime = m.ROOT, m.BASELINE, m.runtime
STEM = 'current-mass-pairing-integrable-return'
digest = lambda p: sha256(p.read_bytes()).hexdigest()


def rational(x):
    return Fraction(int(x['numerator']), int(x['denominator']))


def ball(x):
    return arb(x.numerator)/arb(x.denominator)


def execute():
    ctx.prec = 768
    p, geo, bs = fr.load()
    v = {k:ball(rational(x)) for k,x in
         m.result('current-correlated-signed-trial-test-attempt3')['exact_nominal_parameters'].items()}
    adj = m.constants('current-filtered-mass-adjoint-bound-attempt2')
    bridge = m.constants('current-filtered-test-mass-trace-attempt2')
    green = m.constants('current-filtered-mixed-green-residual')
    detail = m.result('current-source-adapted-fractional-return')['derivative_accounting']
    profiles = m.result('current-filtered-j-regularity-test-attempt2')['whole_spectrum_profile_bounds']
    mp = {k:dec(x) for k,x in m.result('current-global-mass-conversion')['parameters'].items()}
    L, ga, a, beta, B2, C = [p[k] for k in ['L','gamma','a','beta','B2','C']]
    muY = 4*a/3
    assert muY.overlaps(geo['muY']) and muY > 1
    assert arb('.99') < L < 1 and arb('.16') < beta < arb('.17')
    assert all(v[k] > 0 for k in ['alpha','gamma','A','twice_a','muS','muT','muY','scalar_c','kappa','delta'])

    A0 = 3*v['alpha']/(4*v['gamma']*v['A'])
    Lf = (1+abs(v['muY']-v['twice_a'])/v['twice_a'])*(1+2*v['scalar_c']/v['muS'])/v['muT']
    h = v['d1']/2
    Delta = h*h-v['d0']
    u = abs(Delta)/(h*h)
    assert h > 0 and u < 1 and v['Cv']-v['delta']*h > 0
    LF = v['kappa']*((v['Cv']-v['delta']*h)/(h*h)+v['delta']/h)/(1-u)
    assert LF < arb('.251')

    R1 = sum(((j-1)*abs(bs[j])*L**(j-3)/a for j in range(3,7)),arb(0))
    Am = arb(2)
    Bm = (2*(L*arb('2.1')+arb('1.2'))/3+2*L*muY)/L
    Cm = (10+L*R1)/L
    Em = 2/(3*L)
    Dm = arb('6.3')+arb('8.1')*L
    Krec = Am+Bm/muY+2*L*Cm+Em*L*Dm/muY
    U0 = ga*A0*Lf*Krec
    U1 = LF*U0

    # Original polynomial geometry on the WHOLE mass support, including
    # its transition. Negative coefficients are not replaced by positive
    # ones in the upper test when that would lose a required guard.
    coeff = [beta, arb(0), mp['u'], -4*beta-2*mp['u']-mp['tau'],
             3*beta+mp['u']+mp['eps'], -2*mp['eps'], mp['eps']]
    zc = arb('.2')
    Qt = sum((abs(x)*zc**k for k,x in enumerate(coeff) if k),arb(0))
    Qupper = beta+sum((abs(x)*zc**k for k,x in enumerate(coeff) if k and not x < 0),arb(0))
    Qp = sum((k*abs(x)*zc**(k-1) for k,x in enumerate(coeff) if k),arb(0))
    assert beta-Qt > arb('.15') and Qupper < arb('.17') and Qp < arb('.3')
    H3 = arb('1.5')*(arb('.17')/arb('.15')+arb('.17')*arb('.3')*zc/(4*arb('.15')**2))
    Hsup = H3*zc**3
    Hcontact = H3*(zc**3+zc**2/arb(5).sqrt())
    alpha_q = 2*C/(B2*beta)
    Lq = arb(2000)**8/(C*beta*B2)
    assert alpha_q.overlaps(adj['quantum_action_conversion_alpha_q'])
    Kg = adj['actual_filtered_mass_functional_norm_upper']
    q = arb('.338')
    weighted = Kg*(q*U0+U1)
    Eglue = weighted*alpha_q*arb('6e24')*Hsup
    Ep0 = weighted*2*Lq*arb('1e20')*Hcontact
    Ezero = Eglue+Ep0
    target = arb('1.49e34')
    assert Ezero < arb('1e28') and Ezero/target < arb('1e-6')

    # The accepted Bernstein inequalities hold on the entire auxiliary
    # spectrum. Keep the prescribed averaging; the point source has no
    # order-one comparison for f0.
    K0, K1 = [ball(rational(x['coefficient_sum'])/rational(x['Bernstein_guard'])) for x in profiles[:2]]
    radius = arb('1e-110')
    P0 = ga/radius
    Wmax = 2*(arb(4)/3).exp()
    X10 = ga.sqrt()*K0*((1+P0).log()/2+Wmax*P0/(1+P0)).sqrt()
    X11 = ga.sqrt()*K1/2
    S1 = dec(detail['reconstructed_reference_derivative_constants'][1])
    C1 = arb(10).sqrt()*S1
    Ageo = 2*alpha_q*Kg*(2*arb.pi()).sqrt()*C1*(q*X10+X11)
    norm_conversion = bridge['full_E_to_retained_Et_norm_comparison']
    Ebridge = bridge['full_to_retained_j_energy_error']
    conditional_remainder = (target-Ebridge)/norm_conversion-Ezero
    Mh_ceiling = conditional_remainder/Ageo
    oldfloor = adj['separated_gluing_estimator_expression_lower']
    assert conditional_remainder > 0 and oldfloor/Ezero > arb('1e12')
    vals = dict(actual_gamma=ga, nominal_gamma=v['gamma'], source_amplitude_A0=A0,
        complete_fj_L1_factor=Lf, complete_Fref_Neumann_ratio=u, complete_Fref_L1_bound=LF,
        actual_Rq_derivative_bound=R1, actual_reconstruction_L1_constant=Krec,
        complete_trial0_shape_L1_upper=U0, complete_trial1_shape_L1_upper=U1,
        original_mass_support_Q_lower=beta-Qt, original_mass_support_Q_upper=Qupper,
        original_mass_support_Qprime_upper=Qp, actual_mass_shape_cubic_bound=H3,
        actual_mass_shape_Linfinity_upper=Hsup, actual_weighted_contact_factor=Hcontact,
        quantum_action_alpha_q=alpha_q, constant_contact_Lq=Lq,
        complete_gluing_block_upper=Eglue, complete_constant_contact_block_upper=Ep0,
        complete_zero_order_mass_block_upper=Ezero,
        zero_order_block_fraction_of_conditional_target=norm_conversion*Ezero/target,
        zero_order_block_original_bolt_charge=green['pref_t']*norm_conversion*Ezero,
        rejected_old_estimator_lower_over_new_zero_order_upper=oldfloor/Ezero,
        actual_averaging_radius=radius, complete_spectral_split_P0=P0,
        trial0_order_one_auxiliary_norm_upper=X10, trial1_order_one_auxiliary_norm_upper=X11,
        actual_order_one_reconstruction_constant=C1,
        principal_pairing_coefficient_Ageo=Ageo,
        conditional_mass_shape_ceiling_ALL_other_errors_zero=Mh_ceiling,
        unchanged_response_bridge=Ebridge, full_to_retained_norm_conversion=norm_conversion,
        conditional_full_E_account_after_bridge_and_zero_order=conditional_remainder)
    assert all(x.is_finite() and x>0 for x in vals.values())
    return dict(constants={k:enc(x) for k,x in vals.items()},
        mass_shape_fractional_norm=None, principal_mass_block=None,
        complete_derivative_and_scalar_mass_remainder=None,
        signed_r0_mass_pairing=None, signed_r1_mass_pairing=None,
        actual_signed_mass_vector_norm=None, mass_free_residual=None,
        complete_application_error=None, actual_Dj=None, actual_signed_Gt=None,
        scope=dict(actual_structured_L1_return_bounded=True,
            complete_zero_order_mass_block_bounded=True,
            actual_first_flow_derivative_bounded=True,
            principal_pairing_well_defined_on_stated_fixed_inputs=True,
            full_signed_mass_term_certified=False, signed_response_certified=False,
            complete_bolt_allowance=False, available_physical_margin=False,
            original_P32_defect=False, actual_kappa_m=False, complete_original_admission=False,
            nonlinear_neighborhood_transfer=False, nonlinear_existence=False,
            exact_compatibility=False, physical_matching=False, semiclassical_validity_noise=False,
            all_sector_stability=False, construction_obstruction=False, full_goal_complete=False))


if __name__ == '__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--receipt',default=STEM)
    args=ap.parse_args()
    out=ROOT/'evidence'/(args.receipt+'.json')
    assert not out.exists(), 'Preserve immutable receipts'
    old=json.loads((ROOT/'evidence/current-filtered-mass-adjoint-bound-attempt2.json').read_bytes())
    frac=json.loads((ROOT/'evidence/current-source-adapted-fractional-return.json').read_bytes())
    inputs=list(dict.fromkeys([
        'notes/'+STEM+'.md','checks/current_mass_pairing_integrable_return.py',
        'evidence/current-filtered-mass-adjoint-bound-attempt2.json',
        'evidence/current-filtered-mass-adjoint-bound-attempt2-verification.json',
        'evidence/current-source-adapted-fractional-return.json',
        'evidence/current-source-adapted-fractional-return-verification.json',
        'checks/current_source_adapted_fractional_return.py',
        *old['input_sha256_after'],*frac['input_sha256_after']]))
    before={p:digest(ROOT/p) for p in inputs}
    protected=json.loads(BASELINE.read_bytes())
    rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
        arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
    try:
        for prior in [old,frac]:
            assert all(digest(ROOT/p)==h for p,h in prior['input_sha256_after'].items())
        assert all(digest(Path(p))==h for p,h in protected.items())
        assert digest(Path(sys.executable))==runtime['python_sha256']
        for pkg in runtime['packages'].values():
            assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
        rec['result']=execute()
        after={p:digest(ROOT/p) for p in inputs}
        assert before==after
        assert all(digest(Path(p))==h for p,h in protected.items())
        rec.update(status='structured_L1_return_and_zero_order_mass_block_bounded_full_response_open',
            input_sha256_after=after,protected_paths_unchanged=len(protected),
            python_sha256=runtime['python_sha256'],
            pinned_packages={k:len(v['files']) for k,v in runtime['packages'].items()})
    except Exception:
        rec.update(status='failed',traceback=traceback.format_exc())
    rec['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
        constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
    if rec['status']=='failed':
        raise SystemExit(1)
