"""Frozen return/threshold audit; does not solve the actual metric equation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json
import traceback
import sympy as S
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'auxiliary-physical-return-feasibility'
PREVIOUS = 'anomaly-complete-auxiliary-response'
OWN = ['notes/'+STEM+'.md', 'checks/auxiliary_physical_return_feasibility.py',
       'evidence/'+STEM+'-review.json']
ROWS, FIXTURES = [], []


def digest(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()


def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name, passed=bool(ok), detail=str(detail)))
    if not ok:
        raise AssertionError(name)


def algebra():
    t, B, r, u, e, v, a, c, z = S.symbols('t B r u e v a c z', positive=True)
    # Exact sphere return with r_E,new=r_E,0 exp(ell_E).
    rnew = S.exp(-u-t*v)*S.exp(u)*r*S.exp(t*c/2)
    Vnew = S.log((rnew-B)/(r-B))
    check('linear_radius_return', S.simplify(S.diff(Vnew,t).subs(t,0)
          -(c/2-v)*r/(r-B)) == 0)
    orbit = S.exp(-2*(u+t*v))*S.exp(2*u)*S.exp(t*a)
    sphere = (rnew/r)**2
    check('linear_conformal_cancellation', S.simplify(
          S.diff(S.log(orbit/sphere)/2,t).subs(t,0)-(a-c)/2) == 0)
    check('finite_conformal_cancellation', S.simplify(orbit/sphere-S.exp(t*(a-c))) == 0)
    check('exact_finite_radius_return', S.simplify(
          (r*S.exp(z)-B)/(r-B)-(1+(S.exp(z)-1)/(1-B/r))) == 0)
    s, h = S.symbols('s h', positive=True)
    lp = S.Rational(3,5)/(1+s)+S.Rational(2,5)/(s+h)-1/s
    check('radius_weight_monotonicity', S.simplify(lp*s*(1+s)*(s+h)
          +(3*s+2*h*s+5*h)/5) == 0)
    n, y = S.symbols('n y', positive=True)
    wave = S.sin(n*y)/n**2
    check('plateau_second_derivative', S.simplify(S.diff(wave,y,2)+S.sin(n*y)) == 0)
    check('plateau_fifth_derivative', S.simplify(S.diff(wave,y,5)-n**3*S.cos(n*y)) == 0)
    check('HS_b_projection', S.simplify(2*a*a+2*c*c-(a-c)**2-(a+c)**2) == 0)
    FIXTURES.append(dict(fixture='exact_return_and_local_frequency_witness',
                         source_specific_metric_solution=False,
                         cutoff_and_connection_lower_orders='ordinary compact-annulus proof in note'))


def budgets():
    h = Q(1,8*2000**10)
    Wr, Wc = Q('3e20'), Q('4.7e18')
    check('radius_weight_rational_enclosure', 4*(1+1/h)**3 < Wr**5)
    check('cold_weight_rational_enclosure', Q(4,64**5)/h**3 < Wc**5)
    a0, b0 = Q('1.51e-28')/Q('3.2e8'), Q('1e-27')/Q('3.2e8')
    check('scalar_source_substitution', a0 == Q('4.71875e-37'))
    check('scalar_mixed_substitution', b0 == Q('3.125e-36'))
    target = Q('5e-9')
    M = Q('3.3e-29')
    check('rounded_radius_target', Wr*((Q(1,2)+b0)*M+a0) < target)
    check('rounded_cold_b_target', Wc*M/2 < target)
    check('scalar_radius_with_metric_input_one', Wr*(a0+b0) < Q('1.1e-15'))
    check('scalar_radius_without_metric_input', Wr*a0 < Q('1.42e-16'))
    check('scalar_cold_logradius', Wc*(a0+b0) < Q('2e-17'))
    rho, Y0, Z0, Z2 = Q('1e-8'), target, Q('.25'), Q('5e7')
    check('ball_preservation_allocation', Y0+Z0*rho+Z2*rho*rho/2 == rho)
    check('ball_contraction_allocation', Z0+Z2*rho == Q('.75'))
    check('interior_admission_threshold', (Q('.001')-Q('1.3e-4'))/rho == 87000)
    check('cold_background_target', Q('5e-7')-rho == Q('4.9e-7'))
    check('outer_background_target', Q('1e-7')-rho == Q('9e-8'))
    check('old_outer_upper_has_no_slack', Q('1.01e-7') > Q('1e-7'))
    check('reference_energy_embedding_threshold', target/Q('1.44e-10') == Q(625,18))
    check('transferred_energy_embedding_threshold', target/Q('1.92e-10') == Q(625,24))
    check('known_cold_embedding_exceeds_targets', Q('2.6e15') > Q(625,18))
    check('block_response_allocation', (Q('3e-9')+Q('7.5e-10'))/(1-Q('.25')) == target)
    # Published decimal allowances are diagnostics; their exact prior
    # rational provenance remains in the unchanged retained receipt.
    g = Q(1,2000**10)
    Cf, Ci = Q('1.069626977892751e201'), Q('4.141944868545865e118')
    start, mhi, R, CD = Q('3.123143830205975e28'), Q('2.577600186401461e12'), Q('.97'), Q(2)
    records = []
    for ec in [Q(1),Q('1e-4')]:
        budget = Q('1e57')*g*ec/800
        rt = budget/start
        prod = budget*g/(mhi*Cf)
        cft = budget*g/(R*mhi*CD)
        check('legacy_low_budget_negative_'+str(ec), R*start > budget)
        check('legacy_interior_failure_'+str(ec), 800*R/g*mhi*CD*Ci/g/(Q('1e57')*ec) > Q('1.7374e143')/ec)
        check('legacy_global_failure_'+str(ec), 800*R/g*mhi*CD*Cf/g/(Q('1e57')*ec) > Q('4.4868e225')/ec)
        records.append(dict(epsilon_c=str(ec), R_target=str(rt), R_CD_target=str(prod), CF_target=str(cft)))
    FIXTURES.append(dict(fixture='exact_rational_allowance_tests', legacy_records=records,
                         old_action_bounds_transferred_to_new_action=False,
                         radius_N2_threshold_rounded=str((target/Wr-a0)/(Q(1,2)+b0))))


def numerical():
    records = []
    for digits in [90,140]:
        mp.mp.dps = digits
        h = 1/(8*mp.mpf(2000)**10)
        Wr = mp.power(2,mp.mpf('.4'))*mp.power(1+1/h,mp.mpf('.6'))
        Wc = mp.power(2,mp.mpf('.4'))/(64*mp.power(h,mp.mpf('.6')))
        a0, b0, y = mp.mpf('4.71875e-37'),mp.mpf('3.125e-36'),mp.mpf('5e-9')
        linear = (y/Wr-a0)/(mp.mpf('.5')+b0)
        finite = (-mp.expm1(-y)/Wr-a0)/(mp.mpf('.5')+b0)
        check('N2_threshold_range_'+str(digits), mp.mpf('3.3e-29') < linear < mp.mpf('3.5e-29'))
        check('finite_radius_threshold_smaller_'+str(digits), 0 < finite < linear)
        worst = mp.mpf(0)
        for tau in [h/(1+h),mp.mpf('1e-20'),mp.mpf('1e-6'),mp.mpf('.1'),mp.mpf('.99')]:
            for sign in [-1,1]:
                zz = sign*tau*(-mp.expm1(-y))
                # Independent exact physical radius and stable formula.
                B = 1-tau
                direct = mp.log((mp.exp(zz)-B)/tau)
                stable = mp.log1p(mp.expm1(zz)/tau)
                err = abs(direct-stable)
                check('finite_radius_direct_'+str((digits,str(tau),sign)), err < mp.mpf('1e-50'))
                check('finite_radius_inclusion_'+str((digits,str(tau),sign)), abs(stable) <= y)
                worst = max(worst,err)
        records.append(dict(digits=digits,h_min=mp.nstr(h,65),Wr=mp.nstr(Wr,65),Wc=mp.nstr(Wc,65),
                            linear_N2_threshold=mp.nstr(linear,65),finite_N2_threshold=mp.nstr(finite,65),
                            max_direct_radius_error=mp.nstr(worst,20)))
    for key in ['h_min','Wr','Wc','linear_N2_threshold','finite_N2_threshold']:
        check('cross_precision_'+key, abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1) < mp.mpf('1e-55'))
    FIXTURES.append(dict(fixture='independent_finite_return_and_thresholds',records=records,
                         continuum_solution_computed=False))


def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
                 {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
             frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
             previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_auxiliary_result',old['status']=='passed' and old['anomaly_complete_auxiliary_row_equivalence_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();budgets();numerical()
    except Exception:
        error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,
               input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
               finished_utc=datetime.now(timezone.utc).isoformat(),
               exact_physical_radius_return_derived=error is None,
               scalar_product_conformal_cancellation_derived=error is None,
               N2_to_raw_C5_embedding_rejected=error is None,
               Einstein_frame_metric_response_bounded=False,actual_coupled_cold_response_bounded=False,
               full_transformed_boundary_domain_quantified=False,admission_norm_inverse_quantified=False,
               nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
               physical_stability_proved=False,component_complete=False,goal_complete=False,
               formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json')
    receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=run['status'],exception=error)]:
            log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':
        raise SystemExit(1)


if __name__=='__main__':
    main()
