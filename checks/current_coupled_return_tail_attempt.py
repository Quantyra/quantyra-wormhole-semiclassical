"""Exact conditional tail arithmetic. Does NOT apply H, H0, R or T to fields."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from datetime import datetime, timezone
import fractions
import json
import sys
import traceback

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-coupled-return-tail-attempt'
OWN = ['notes/'+STEM+'.md', 'checks/current_coupled_return_tail_attempt.py']
INPUTS = OWN + [
    'evidence/current-coupled-return-tail-attempt-startup-failure.json',
    'notes/current-geodesic-two-sided-remainder.md',
    'notes/current-geodesic-flow-multiplier.md',
    'notes/current-global-geodesic-return.md',
    'notes/current-coupled-schur-graph-attempt.md',
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-global-return-energy-attempt.md',
    'notes/current-implicit-phase-graph-domain.md',
    'notes/current-complete-strong-response-attempt.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-geodesic-flow-multiplier.json']
CHECKS = []


def check(name, condition):
    CHECKS.append(dict(name=name, passed=bool(condition)))
    if not condition:
        raise AssertionError(name)


def enclosure(value):
    ends = []
    for rounding in [ROUND_FLOOR, ROUND_CEILING]:
        with localcontext() as ctx:
            ctx.prec = 24
            ctx.rounding = rounding
            ends.append(str(Decimal(value.numerator)/Decimal(value.denominator)))
    if not F(ends[0]) <= value <= F(ends[1]):
        raise ArithmeticError('Outward decimal enclosure failed')
    return dict(lower=ends[0], upper=ends[1])


def execute():
    eta = F('.0022')
    q = eta/(1-eta)
    check('exact_full_form_ratio', q == F(11, 4989) and q < F('.002205'))
    # Formal polynomial identity sum (-k)^j(1+k)=1-(-k)^N,
    # checked coefficientwise, not by a finite scalar/matrix sample.
    for N in [16, 32, 64, 128, 256]:
        coeff = [0]*(N+1)
        for j in range(N):
            coeff[j] += (-1)**j
            coeff[j+1] += (-1)**j
        check('formal_telescoping_'+str(N), coeff == [1]+[0]*(N-1)+[-(-1)**N])
    km, am, cg, ce = map(F, ['5e-5', '4.2953e-8', '49.462', '390.478'])
    saff, sfb = F('7.49e-12'), F('.00368')
    raff, rfb = saff-km*am, sfb-km*cg
    check('conditional_reserves', raff == F('5.34235e-12') and rfb == F('.0012069'))
    check('independent_energy_debit', km*ce == F('.0195239'))
    ratio = rfb/(km*ce)
    check('restricted_family_ratio', F('.06181654280138') < ratio < F('.06181654280139'))
    daff, dfb, Uaff = raff/2, rfb/2, F('1.1e-10')
    rows = {}
    for N in [16, 32, 64, 128, 256]:
        tail = q**(N-4)
        rows[str(N)] = dict(q_N_minus_4=enclosure(tail),
            affine_M_graph_4_ceiling=enclosure(daff/(Uaff*tail)),
            feedback_M_graph_4_times_U_feedback_ceiling=enclosure(dfb/tail),
            exact_operator_affine_energy_tail=enclosure(Uaff*q**N))
    check('64_step_rounded_down_affine_allowance', F('6.0612e157')*q**60*Uaff < daff)
    check('64_step_rounded_down_feedback_product', F('1.5062e156')*q**60 < dfb)
    check('32_step_allowances', F('5.8978e72')*q**28*Uaff < daff and F('1.4656e71')*q**28 < dfb)
    check('128_step_allowances', F('6.4016e327')*q**124*Uaff < daff and F('1.5908e326')*q**124 < dfb)
    screen = F('.028913')+F('1e-11')*F('1.0706e10')
    check('existing_fast_screen_unchanged', screen == F('.135973') and F('.1599164')-screen == F('.0239434'))
    return dict(outcome='unresolved_complete_coupled_return',
        arithmetic='exact rational; outward decimal endpoints checked against exact fractions',
        q=enclosure(q), candidate_power=4,
        source_budget=enclosure(saff), feedback_budget=enclosure(sfb),
        conditional_mass_debit=enclosure(km*am), conditional_mass_feedback=enclosure(km*cg),
        residual_affine_budget=enclosure(raff), residual_feedback_budget=enclosure(rfb),
        independent_energy_debit=enclosure(km*ce), restricted_energy_ratio=enclosure(ratio),
        proposed_error_affine_allocation=enclosure(daff), proposed_error_feedback_allocation=enclosure(dfb),
        conditional_tables=rows,
        physical_operator_applications_executed=0,
        strong_power_domain_inclusion_proved=False, M_graph_4=None, U_feedback=None,
        complete_signed_center=None, kappa_mass=None, other_coupled_columns=None,
        actual_mass=None, original_physical_norm_comparisons_closed=False,
        nonlinear_admission=False, exact_compatibility=False, physical_matching=False,
        validity_noise=False, all_sector_stability=False,
        new_feasibility_requirement_closed=False, actual_controlled_solution=False,
        construction_obstruction=False, decisive_next_calculation_executed=False,
        provenance_limit='One prior interactive exact-rational threshold exploration; its stdout is retained in the conversation, not a separate immutable source artifact. No exploratory scientific failure.')


if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    digest = lambda p: sha256(p.read_bytes()).hexdigest()
    hashes = {p: digest(ROOT/p) for p in INPUTS}
    module_path = Path(fractions.__file__)
    # This Python distribution stores its stdlib inside python313.zip.
    module_artifact = (module_path if module_path.is_file() else
                       next(p for p in module_path.parents if p.is_file()))
    run = dict(started_utc=datetime.now(timezone.utc).isoformat(), mode=mode,
        input_sha256_before=hashes,
        frozen_sources_utf8={p: (ROOT/p).read_text(encoding='utf-8') for p in OWN},
        runtime=dict(python=sys.version, python_executable_sha256=digest(Path(sys.executable)),
                     fractions_artifact=str(module_artifact),
                     fractions_artifact_sha256=digest(module_artifact)))
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p: digest(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs': []}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger, indent=2)+'\n').encode())
    log = dict(mode=mode, status=run['status'], checks=len(CHECKS),
               result=run.get('result'), traceback=run.get('traceback'))
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:
        f.write((json.dumps(log, indent=2)+'\n').encode())
    print(json.dumps(log, indent=2))
    if run['status'] != 'passed':
        raise SystemExit(1)
